import itertools as it
from math import comb

from boostree import RTree
import numpy as np

from pyfr.mpiutil import AlltoallMixin, Sorter, get_comm_rank_root, mpi
from pyfr.nputil import morton_encode, zip_chunks
from pyfr.shapes import BaseShape
from pyfr.util import subclass_where


def get_interpolator(name, ndims, is_admissible, opts, order=None):
    cls = subclass_where(BaseInterpolator, name=name)

    return cls.from_opts(ndims, is_admissible, opts, order=order)


def _ploc_op(etype, nspts, cfg):
    # Shape and associated basis
    shapecls = subclass_where(BaseShape, name=etype)
    basis = shapecls(nspts, cfg)

    # Interpolate the shape points to the solution points
    return basis.sbasis.nodal_basis_at(basis.upts)


def str_group(gpts, nsplit):
    def recursive_tile(pts, depth=0):
        if depth == pts.shape[1]:
            return [pts]

        sidx = np.argsort(pts[:, depth])
        groups = []

        # Recurse in the next dimension
        for chunk in np.array_split(pts[sidx], nsplit):
            groups.extend(recursive_tile(chunk, depth + 1))

        return groups

    return recursive_tile(gpts)


class BaseInterpolator:
    name = None

    int_opts = {}
    float_opts = {}
    enum_opts = {}
    dflt_opts = {}

    @classmethod
    def from_opts(cls, ndims, is_admissible, opts={}, *, order=None):
        # Normalise hyphens to underscores in user-supplied keys
        opts = {k.replace('-', '_'): v for k, v in opts.items()}

        kwargs = {}
        for k, v in dict(cls.dflt_opts, **opts).items():
            if k in cls.int_opts:
                kwargs[k] = int(v)
            elif k in cls.float_opts:
                kwargs[k] = float(v)
            elif k in cls.enum_opts:
                kwargs[k] = cls.enum_opts[k][v]
            else:
                raise ValueError('Invalid option')

        if order is not None:
            kwargs['order'] = order

        return cls(ndims, is_admissible, **kwargs)


class IDWInterpolator(BaseInterpolator):
    name = 'idw'

    # Integer options
    int_opts = {'n'}

    # Floating point options
    float_opts = {'rho'}

    # Default options
    dflt_opts = {'n': 0, 'rho': 0}

    def __init__(self, ndims, is_admissible, *, order=None, n=0, rho=0):
        self.is_admissible = is_admissible
        self.n = n or 2**ndims
        self.rho = rho or ndims + 1
        self.eps = 10*np.finfo(float).eps

    def __call__(self, pts, spts, svals):
        dists = np.linalg.norm(spts - pts[:, None], axis=2)

        swts = np.maximum(dists, self.eps)**-self.rho
        swts /= swts.sum(axis=1, keepdims=True)
        out = (swts[:, None] @ svals)[:, 0]

        # Exact hits take the value of the coincident source point
        amin = dists.argmin(axis=1)
        dmin = np.take_along_axis(dists, amin[:, None], 1)[:, 0]
        if (ex := dmin < self.eps).any():
            out[ex] = svals[ex, amin[ex]]

        return out


class WENOInterpolator(BaseInterpolator):
    name = 'weno'

    int_opts = {'n', 'degree', 'sub_degree', 'nsub', 'sub_n'}
    float_opts = {'rho', 'q', 'ct', 'cond', 'dir_bias', 'gamma0'}
    enum_opts = {'mode': {'wenoz': 'wenoz', 'teno': 'teno'}}
    dflt_opts = {'n': 0, 'degree': 0, 'sub_degree': 0, 'nsub': 0,
                 'sub_n': 0, 'rho': 0, 'q': 4, 'ct': 1.0e-3,
                 'cond': 1.0e8, 'dir_bias': 2.5, 'gamma0': 0.85,
                 'mode': 'wenoz'}

    def __init__(self, ndims, is_admissible, *, order=None, n=0, degree=0,
                 sub_degree=0, nsub=0, sub_n=0, rho=0, q=4, ct=1.0e-3,
                 cond=1.0e8, dir_bias=2.5, gamma0=0.85, mode='wenoz'):
        # Auto-derive degree from the source polynomial order
        if not degree:
            degree = min(order, 3) if order else 2

        if not sub_degree:
            sub_degree = max(degree - 1, 1)

        # Auto-derive the stencil size to be ~3x overdetermined
        nterms = comb(ndims + degree, degree)
        if not n:
            n = max(3*nterms, 24 if ndims == 2 else 48)

        self.is_admissible = is_admissible
        self.n = n
        self.q = q
        self.ct = ct
        self.cond = cond
        self.dir_bias = dir_bias
        self.gamma0 = gamma0
        self.eps = 10*np.finfo(float).eps

        # Keep an IDW fallback for exact hits and rejected stencils
        self._idw = IDWInterpolator(ndims, is_admissible, n=self.n,
                                    rho=rho or ndims + 1)

        # Precompute the polynomial bases used by the central and side fits
        cpowers = self._monomial_powers(ndims, degree)
        cexp = np.array(cpowers, dtype=int)
        self._central_power_data = (cpowers, cexp, cexp.sum(axis=1))

        spowers = self._monomial_powers(ndims, min(sub_degree, degree))
        sexp = np.array(spowers, dtype=int)
        self._sub_power_data = (spowers, sexp, sexp.sum(axis=1))

        nsub = nsub or 2**ndims
        self._sub_nterms = len(self._sub_power_data[0])

        # Use an overdetermined directional fit when enough points are present
        self.sub_n = sub_n or max(2*self._sub_nterms, self._sub_nterms + 3)
        self._directions = self._direction_vectors(ndims, nsub)

        # Stencil sizes determine which candidate fits are possible
        self._has_central = self.n >= len(cpowers)
        self._has_dirs = self.n >= self._sub_nterms

        # Select the appropriate nonlinear weighting
        if mode == 'teno':
            self._weights = self._teno_weights
        else:
            self._weights = self._wenoz_weights

    def __call__(self, pts, spts, svals):
        rel = spts - pts[:, None]
        dists = np.linalg.norm(rel, axis=2)

        # Sort each stencil by its distance from the query point
        order = np.argsort(dists, axis=1)
        dists = np.take_along_axis(dists, order, 1)
        rel = np.take_along_axis(rel, order[:, :, None], 1)
        svals = np.take_along_axis(svals, order[:, :, None], 1)

        # Work in a local scaled coordinate system to improve conditioning
        scale = np.maximum(dists[:, -1:], 1.0e-14)
        x = rel / scale[:, :, None]

        cands = []

        # First attempt to build a central fit
        if self._has_central:
            cands.append(self._fit_stencils(x, svals,
                                            self._central_power_data))

        # Then attempt to build directional substencils
        if self._has_dirs:
            ncore = min(max(self._sub_nterms - 1, 2), self.n)
            sub_n = min(self.sub_n, self.n)

            unit = rel / np.maximum(dists, 1e-300)[:, :, None]
            for direction in self._directions:
                align = unit @ direction

                # Penalise misaligned points by up to (1 + 2*dir_bias)
                score = dists*(1 + self.dir_bias*(1 - align))

                # Stencils are sorted so the core is the leading block
                score[:, :ncore] = -np.inf
                idx = np.argpartition(score, sub_n - 1, axis=1)[:, :sub_n]

                xs = np.take_along_axis(x, idx[:, :, None], 1)
                vs = np.take_along_axis(svals, idx[:, :, None], 1)
                cands.append(self._fit_stencils(xs, vs,
                                                self._sub_power_data))

        # Combine the candidate fits with WENO-Z or TENO weights
        if cands:
            vals = np.stack([c[0] for c in cands], axis=1)
            betas = np.stack([c[1] for c in cands], axis=1)
            oks = np.stack([c[2] for c in cands], axis=1)

            weights = self._nonlinear_weights(betas, oks)
            val = (weights[:, None] @ vals)[:, 0]
            nc = oks.any(axis=1)
        else:
            val = np.zeros_like(svals[:, 0])
            nc = np.zeros(len(pts), dtype=bool)

        # Rejected and inadmissible fits revert to the convex IDW fallback
        nc &= self.is_admissible(val.T)

        idw = self._idw(np.zeros_like(pts), rel, svals)
        out = np.where(nc[:, None], val, idw)

        # Exact hits take the value of the coincident source point
        if (ex := dists[:, 0] < self.eps).any():
            out[ex] = svals[ex, 0]

        return out

    @staticmethod
    def _monomial_powers(ndims, degree):
        powers = [p for p in it.product(range(degree + 1), repeat=ndims)
                  if sum(p) <= degree]
        powers.sort(key=lambda p: (sum(p), p))

        return powers

    @staticmethod
    def _monomial_matrix(pts, power_data):
        _, exponents, _ = power_data

        # Each column is one monomial evaluated at every sample point
        vdm = np.ones((*pts.shape[:-1], len(exponents)))
        for i, exp in enumerate(exponents.T):
            vdm *= pts[..., i:i + 1]**exp

        return vdm

    @staticmethod
    def _wendland_c2(r):
        # Compactly supported C2 kernel: psi(r) = (1-r)_+^4 (4r + 1)
        return np.clip(1 - r, 0.0, None)**4*(4*r + 1)

    @staticmethod
    def _direction_vectors(ndims, nsub):
        # In 2D use evenly spaced angular substencils
        if ndims == 2:
            ang = np.linspace(0.0, 2*np.pi, nsub, endpoint=False)
            return np.column_stack([np.cos(ang), np.sin(ang)])
        # In higher dimensions use axis-balanced corner directions
        else:
            dirs = np.array(list(it.product((-1.0, 1.0), repeat=ndims)))
            dirs /= np.linalg.norm(dirs, axis=1)[:, None]

            return dirs[:nsub]

    def _fit_stencils(self, x, svals, power_data):
        _, _, orders = power_data

        radii = np.linalg.norm(x, axis=2)
        rmax = np.maximum(radii.max(axis=1, keepdims=True), 1e-14)

        weights = self._wendland_c2(radii / rmax)
        weights = np.maximum(weights, 1e-12)

        # Normalised up front; a uniform row scaling leaves the fit alone
        weights /= weights.sum(axis=1, keepdims=True)

        # The compact weights taper distant points without a hard cutoff
        vdm = self._monomial_matrix(x, power_data)
        sqrtw = np.sqrt(weights)
        amat = sqrtw[:, :, None]*vdm
        rhs = sqrtw[:, :, None]*svals

        # A single SVD gives both the condition estimate and the solve
        u, sval, vh = np.linalg.svd(amat, full_matrices=False)
        ok = sval[:, 0] <= self.cond*sval[:, -1]

        svs = np.where(sval > 0, sval, 1)
        urhs = (u.mT @ rhs) / svs[:, :, None]
        coeffs = vh.mT @ urhs
        fit = vdm @ coeffs

        # Weighted variance and residual, aggregated across all nvars
        mean = (weights[:, None] @ svals)[:, 0]
        var = (weights*((svals - mean[:, None])**2).sum(axis=2)).sum(axis=1)
        resid = (weights*((fit - svals)**2).sum(axis=2)).sum(axis=1)

        # Higher-order coefficient energy, weighted by total order
        cvar = ((orders[1:, None]*coeffs[:, 1:])**2).sum(axis=(1, 2))

        # Blend residual and higher-order energy into one smoothness beta
        blend = resid + 1e-2*cvar
        norm = var + 1e-14*(mean**2).sum(axis=1)

        nz = norm > 0
        beta = np.where(nz, blend / np.where(nz, norm, 1), 0)

        return coeffs[:, 0], np.where(ok, beta, np.inf), ok

    def _wenoz_weights(self, betas, gamma, tau, eps=1.0e-14):
        alpha = gamma*(1 + (tau[:, None] / (betas + eps))**self.q)
        return alpha / np.maximum(alpha.sum(axis=1, keepdims=True), eps)

    def _teno_weights(self, betas, gamma, tau, eps=1.0e-14):
        alpha = gamma*(1 + (tau[:, None] / (betas + eps))**self.q)
        chi = alpha / np.maximum(alpha.sum(axis=1, keepdims=True), eps)

        # TENO hard-rejects candidates whose normalized weight is too small
        mask = chi > self.ct
        if (none := ~mask.any(axis=1)).any():
            mask[none, betas[none].argmin(axis=1)] = True

        weights = gamma*mask
        return weights / np.maximum(weights.sum(axis=1, keepdims=True), eps)

    def _nonlinear_weights(self, betas, oks):
        ncands = np.maximum(oks.sum(axis=1), 1)
        rcpn = 1 / ncands

        # Bias the linear weights toward an accepted central stencil
        if self._has_central:
            okc = oks[:, 0]
        else:
            okc = np.zeros(len(oks), dtype=bool)

        goth = np.where(okc, (1 - self.gamma0)/np.maximum(ncands - 1, 1), rcpn)
        gamma = np.where(oks, goth[:, None], 0)
        gamma[:, 0] = np.where(okc, self.gamma0, gamma[:, 0])

        # Tau measures how different the candidate smoothness values are
        bmean = np.where(oks, betas, 0).sum(axis=1)*rcpn
        tau = np.where(oks, np.abs(betas - bmean[:, None]), 0).max(axis=1)

        return self._weights(betas, gamma, tau)


class BaseCloudResampler(AlltoallMixin):
    # Interpolation block size, in points
    BLK_SZ = 2**12

    def __init__(self, pts, solns, interp, progress):
        self.ndims = pts.shape[1]
        self.nvars = solns.shape[1]
        self.interp = interp
        self.progress = progress

        with progress.start('Distribute source points/solutions'):
            self.pts, self.solns = self._distribute_src_points(pts, solns)

        with progress.start('Index source points'):
            # Index the points
            self.pts_tree = self._compute_pts_tree(self.pts)

            # Index the bounding boxes (inclusives and exclusive of ourself)
            self.ibbox_tree, self.ebbox_tree = self._compute_bbox_trees(
                self.pts
            )

        # Track what points we have sent to other ranks
        comm, rank, root = get_comm_rank_root()
        self._fringe_idxs = [set() for i in range(comm.size)]

    def _distribute_src_points(self, pts, solns):
        comm, rank, root = get_comm_rank_root()

        # Obtain the bounding box for our ranks points
        pmin, pmax = pts.min(axis=0), pts.max(axis=0)

        # Reduce over all ranks
        comm.Allreduce(mpi.IN_PLACE, pmin, op=mpi.MIN)
        comm.Allreduce(mpi.IN_PLACE, pmax, op=mpi.MAX)

        # Normalise our points and compute their Morton codes
        fact = 2**21 if len(pmin) == 3 else 2**32
        ipts = (fact*((pts - pmin) / (pmax - pmin).max())).astype(np.uint64)
        keys = morton_encode(ipts, [fact]*len(pmin))

        # Use these codes to sort our points
        sorter = Sorter(comm, keys)

        # With this redistribute our points and solutions
        return sorter(pts), sorter(solns)

    def _compute_pts_tree(self, pts):
        return RTree.from_points(pts, storage='point')

    def _compute_bbox_trees(self, pts):
        comm, rank, root = get_comm_rank_root()

        # Number of splits along each dimension
        ns = 14

        # Allocate global bounding box array
        bboxes = np.full((comm.size, ns**self.ndims, 2, self.ndims), np.nan)

        # Fill out our portion of the array
        for i, spts in enumerate(str_group(pts, ns)):
            if len(spts):
                bboxes[rank, i] = spts.min(axis=0), spts.max(axis=0)

        # Exchange
        comm.Allgather(mpi.IN_PLACE, bboxes)

        # Flatten the array, tagging each box with the rank which owns it
        branks = np.repeat(np.arange(comm.size), bboxes.shape[1])
        bmins, bmaxs = np.moveaxis(bboxes, 2, 0).reshape(2, -1, self.ndims)

        # Discard the entries associated with empty groups
        valid = ~np.isnan(bmins[:, 0])

        # Construct bounding box trees inclusive and exclusive of our rank
        trees = []
        for i in (-1, rank):
            sel = valid & (branks != i)
            trees.append(RTree.from_boxes(bmins[sel], bmaxs[sel], branks[sel]))

        return trees

    def sample_with_mesh_config(self, mesh, cfg):
        pts, sidxs = [], []

        # Numeric data type
        dtype = np.dtype(cfg.get('backend', 'precision', 'double')).type

        # Interpolate the shape points to the solution points
        for etype in mesh.eidxs:
            op = _ploc_op(etype, len(mesh.spts[etype]), cfg)
            plocs = op @ mesh.spts[etype].reshape(op.shape[1], -1)
            plocs = plocs.reshape(-1, self.ndims)

            pts.append(plocs)
            sidxs.append(len(plocs) + (sidxs[-1] if sidxs else 0))

        # Sample the solution at these solution points
        solns = self.sample_with_pts(np.vstack(pts), dtype)
        solns = np.split(solns, sidxs[:-1])

        # Reshape these solutions into their canonical forms
        esolns = {}
        for (etype, eidxs), soln in zip(mesh.eidxs.items(), solns):
            esoln = soln.reshape(-1, len(eidxs), self.nvars)
            esolns[etype] = esoln.transpose(1, 2, 0)

        return esolns

    def sample_with_pts(self, tpts, dtype):
        with self.progress.start('Distribute target points'):
            tpts, tcountdisp, tidx = self._distribute_tgt_pts(tpts)

        with self.progress.start('Sample target points'):
            tsolns = self._sample_tgt_points(tpts, dtype)

        with self.progress.start('Distribute target samples'):
            return self._distribute_tgt_samples(tsolns, tcountdisp, tidx)

    def _distribute_tgt_pts(self, pts):
        comm, rank, root = get_comm_rank_root()

        if comm.size > 1:
            # Attempt to assign points to ranks
            tranks, resid = self._compute_initial_tgt_dist(pts)

            # Exchange residual information between ranks
            rranks = self._compute_resid_ranks(resid)

            # Assign residual points to the rank with the nearest point
            for i, (r, _) in rranks.items():
                tranks[i] = r
        else:
            tranks = np.zeros(len(pts), dtype=int)

        # Sort the points according to their rank
        tidx = np.argsort(tranks)

        # Distribute the points
        tdisps = np.searchsorted(tranks[tidx], np.arange(comm.size))
        tcount = self._disp_to_count(tdisps, len(tranks))
        rbuf = self._alltoallcv(comm, pts[tidx], tcount, tdisps)

        return *rbuf, np.argsort(tidx)

    def _compute_initial_tgt_dist(self, pts):
        comm, rank, root = get_comm_rank_root()

        tranks = np.full(len(pts), -1, dtype=int)

        # See which ranks, if any, have a claim to each point
        iranks, icounts = self.ibbox_tree.intersect(pts, pts)
        icounts = icounts.astype(int)
        ioffs = np.cumsum(icounts) - icounts

        # Points with no intersecting box go to the rank which is nearest
        if len(nidx := np.flatnonzero(icounts == 0)):
            tranks[nidx] = self.ibbox_tree.knn(pts[nidx], 1)[0][:, 0]

        # Residual list for target points claimed by multiple ranks
        resid = [[] for i in range(comm.size)]

        for i in np.flatnonzero(icounts):
            off, c = ioffs[i], icounts[i]
            uranks = set(iranks[off:off + c].tolist())

            # All boxes belong to the same rank so assign
            if len(uranks) == 1:
                tranks[i] = iranks[off]
            # Multiple distinct ranks; defer
            else:
                for r in uranks:
                    resid[r].append((i, pts[i].tolist()))

        return tranks, resid

    def _compute_resid_ranks(self, resid):
        comm, rank, root = get_comm_rank_root()
        rranks, rdists = {}, []

        # Exchange residual points between ranks
        for rresid in comm.alltoall([[p for i, p in r] for r in resid]):
            if rresid:
                # Identify our closest points
                rpts = np.array(rresid)
                dists = self.pts_tree.knn(rpts, 1)[1]

                rdists.append(dists[:, 0].tolist())
            else:
                rdists.append([])

        # Exchange distance information and reduce
        for r, (re, rd) in enumerate(zip(resid, comm.alltoall(rdists))):
            for (i, _), d in zip(re, rd):
                if (rd := rranks.get(i, None)) is None or d < rd[1]:
                    rranks[i] = (r, d)

        return rranks

    def _interp_blocks(self, pts, idxs):
        out = np.empty((len(pts), self.nvars))

        # Chunked as the gathered stencils are much larger than the points
        for ochunk, pchunk, ichunk in zip_chunks(self.BLK_SZ, out, pts, idxs):
            ochunk[:] = self.interp(pchunk, self.pts[ichunk],
                                    self.solns[ichunk])

        return out

    def _sample_tgt_points(self, pts, dtype):
        comm, rank, root = get_comm_rank_root()

        # Allocate the interpolated solution array
        solns = np.empty((len(pts), self.nvars), dtype=dtype)

        nn = self.interp.n
        deferred, off = [], 0
        fpreqs = [[] for i in range(comm.size)]

        # Determine the nearest points to each sample point
        nidxs, ndists = self.pts_tree.knn(pts, nn)
        ndists = ndists[:, -1]

        # Compute the associated bounding boxes
        pmins, pmaxs = pts - ndists[:, None], pts + ndists[:, None]

        # Determine which ranks intersect with these bounding boxes
        iranks, icounts = self.ebbox_tree.intersect(pmins, pmaxs)

        for i, (p, count) in enumerate(zip(pts, icounts)):
            # If any other ranks intersect then defer
            if count:
                deferred.append(i)
                for j in np.unique(iranks[off:off + count]):
                    fpreqs[j].append([*p, ndists[i]])

                off += count

        # Points which no peer rank can contribute to are done locally
        if len(local := np.flatnonzero(icounts == 0)):
            solns[local] = self._interp_blocks(pts[local], nidxs[local])

        self._exchange_fringe_pts(fpreqs, nn)

        # Process the deferred points
        if len(dpts := pts[deferred]):
            didxs = self.pts_tree.knn(dpts, nn)[0]
            solns[deferred] = self._interp_blocks(dpts, didxs)

        comm.barrier()

        # Return the interpolated solution values
        return solns

    def _exchange_fringe_pts(self, frboxes, nn):
        comm, rank, root = get_comm_rank_root()

        # Tally up how many requests we have for each rank
        scount = np.array([len(fr) for fr in frboxes])
        if scount.any():
            sbuf = np.vstack([fr for fr in frboxes if fr])
        else:
            sbuf = np.empty((0, self.ndims + 1))

        # Exchange the requests
        rdata, (_, rdisps) = self._alltoallcv(comm, sbuf, scount)

        fidxs, fcount = [], []

        # Split these up on a per-rank basis
        for i, ibboxes in enumerate(np.split(rdata, rdisps[1:])):
            fcnt = len(fidxs)

            # Identify nearby points on our rank
            pts, ndists = ibboxes[:, :self.ndims], ibboxes[:, self.ndims]
            idxs = self.pts_tree.knn(pts, nn, radius=ndists)[0]
            idxs = set(np.unique(idxs[idxs >= 0]).tolist())

            # Exclude points that have already been sent over
            fidxs.extend(idxs - self._fringe_idxs[i])
            self._fringe_idxs[i].update(idxs)

            # Tally up how many points we are sending to this rank
            fcount.append(len(fidxs) - fcnt)

        # Convert these lists to arrays
        fidxs = np.array(fidxs, dtype=int)
        fcount = np.array(fcount, dtype=int)
        fdisps = self._count_to_disp(fcount)

        # Exchange fringe points and solutions
        pts = self._alltoallcv(comm, self.pts[fidxs], fcount, fdisps)[0]
        solns = self._alltoallcv(comm, self.solns[fidxs], fcount, fdisps)[0]

        # Add the received fringe points to our tree
        self.pts_tree.add(np.arange(len(pts)) + len(self.pts), pts)

        # Incorporate this fringe data into our point and solution arrays
        if len(pts):
            self.pts = np.vstack([self.pts, pts])
            self.solns = np.vstack([self.solns, solns])

    def _distribute_tgt_samples(self, tsolns, tcountdisp, tidx):
        comm, rank, root = get_comm_rank_root()

        # Send the samples back to the ranks they came from
        tsolns = self._alltoallcv(comm, tsolns, *tcountdisp)[0]
        return tsolns[tidx]


class NativeCloudResampler(BaseCloudResampler):
    def __init__(self, mesh, soln, interp, progress):
        cpts, csolns = self._concat_pts_solns(mesh, soln)

        super().__init__(cpts, csolns, interp, progress)

    def _concat_pts_solns(self, mesh, soln):
        pts, solns = [], []

        for etype in mesh.eidxs:
            # Interpolate the shape points to the solution points
            op = _ploc_op(etype, len(mesh.spts[etype]), soln.config)
            ploc = op @ mesh.spts[etype].reshape(op.shape[1], -1)
            pts.append(ploc.reshape(-1, mesh.ndims))

            # Extract the solution at the solution points
            supts = soln.data[etype].swapaxes(1, 2)
            solns.append(supts.reshape(-1, supts.shape[2]))

        return np.vstack(pts), np.vstack(solns)
