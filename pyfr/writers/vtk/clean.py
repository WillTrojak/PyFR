from collections import namedtuple

import numpy as np

from pyfr.mpiutil import DistributedDirectory, get_comm_rank_root, mpi
from pyfr.writers.vtk.shapes import get_vtk_shape


KindGlobal = namedtuple('KindGlobal',
                        'keys defidx nint emap count reducer')


def _unique_rows(keys):
    n, k = keys.shape
    if k == 1:
        uniq, inv = np.unique(keys.ravel(), return_inverse=True)
        return uniq.reshape(-1, 1), inv

    # lexsort + first-of-group is faster than struct-view for k>1
    order = np.lexsort(keys.T)
    srt = keys[order]
    fmask = np.empty(n, dtype=bool)
    fmask[:1] = True
    fmask[1:] = (srt[1:] != srt[:-1]).any(axis=1)
    inv = np.empty(n, dtype=int)
    inv[order] = np.cumsum(fmask) - 1
    return srt[fmask], inv


class Reducer:
    # FxHash multiplicative-XOR constants (golden ratio and sign-bit mask)
    FX_SEED = np.uint64(0x9E3779B97F4A7C15)
    FX_MASK = np.uint64(0x7FFFFFFFFFFFFFFF)

    def __init__(self, comm, fkeys):
        self.dd = DistributedDirectory(comm, self._hash_rows(fkeys))
        self.uniq, self.inv = _unique_rows(self.dd.scatter(fkeys))

    def __call__(self, values):
        rv = self.dd.scatter(values)
        gv = np.zeros((len(self.uniq), *values.shape[1:]), dtype=values.dtype)
        np.add.at(gv, self.inv, rv)
        return self.dd.gather(gv[self.inv])

    @classmethod
    def _hash_rows(cls, keys):
        h = np.zeros(len(keys), dtype=np.uint64)
        for col in keys.view(np.uint64).T:
            h = (h*cls.FX_SEED) ^ col

        return h & cls.FX_MASK


class CleanToGrid:
    @staticmethod
    def _class_rows(classes, nint):
        return (classes[:, None]*nint + np.arange(nint)).ravel()

    def __init__(self, cnodemap, divmap, nsvptsmap, shared):
        self.comm, _, _ = get_comm_rank_root()
        self.cnodemap = cnodemap

        topos, edef, ktypes = {}, {}, set()
        for etype, cnodes in cnodemap.items():
            shape = get_vtk_shape(etype, divmap[etype])
            topos[etype] = shape.topology(nsvptsmap[etype])
            edef[etype] = np.isin(cnodes, shared).any(axis=1)
            ktypes.update(topos[etype].kinds)

        # Reduce to obtain the global set of entities
        ktypes = set().union(*self.comm.allgather(ktypes))

        # Per-kind cross-etype dedup -> self.kinds (with emap populated)
        self.kinds = {kind: self._init_kind(kind, topos, cnodemap, edef)
                      for kind in sorted(ktypes)}

        self.layouts = {etype: self._build_layout(etype, topos[etype])
                        for etype in cnodemap}

    def _init_kind(self, kind, topos, cnodemap, edef):
        evdofs, klist, dlist, nint, ncols = {}, [], [], 0, 0
        for etype, cn in cnodemap.items():
            if kind not in topos[etype].kinds:
                continue

            keys, fsrc, cpos, nint = topos[etype].canonical_vdofs(kind, cn)
            evdofs[etype] = (keys, fsrc, cpos)
            klist.append(keys)
            dlist.append(np.repeat(edef[etype], len(keys) // len(cn)))
            ncols = keys.shape[1]

        # Agree on nint and key width across ranks
        nint = self.comm.allreduce(nint, op=mpi.MAX)
        ncols = self.comm.allreduce(ncols, op=mpi.MAX)

        # Deduplicate canonical keys across all contributing etypes
        akeys = np.concatenate(klist or [np.empty((0, ncols), dtype=int)])
        adefs = np.concatenate(dlist or [np.empty(0, dtype=bool)])

        uniq, inv = _unique_rows(akeys)
        defidx = np.unique(inv[adefs])

        # Map each etype's vdofs to their global class indices
        emap, off = {}, 0
        for etype, (keys, fsrc, cpos) in evdofs.items():
            cls = inv[off:off + len(keys)]
            emap[etype] = (np.unique(cls), fsrc, cls*nint + cpos)
            off += len(keys)

        # Determine counts
        count = np.bincount(inv, minlength=len(uniq)) // max(nint, 1)
        reducer = Reducer(self.comm, uniq[defidx])
        count[defidx] = reducer(count[defidx])

        return KindGlobal(uniq, defidx, nint, emap, count, reducer)

    def _build_layout(self, etype, topo):
        neles, npts = len(self.cnodemap[etype]), topo.npts
        remap = np.empty((neles, npts), dtype=int)
        kept, off = [], 0

        # Merge coincident vdofs into unique kept positions
        for g in self.kinds.values():
            if (kp := g.emap.get(etype)) is None:
                continue

            # Local vdof positions within this etype's classes
            classes, fsrc, fdst = kp
            cls, cpos = divmod(fdst, g.nint)
            lcls = np.searchsorted(classes, cls)
            within = lcls*g.nint + cpos

            # First-occurrence source indices for each kept vdof
            ukept, fidx = np.unique(within, return_index=True)
            ksrc = np.empty(len(classes)*g.nint, dtype=int)
            ksrc[ukept] = fsrc[fidx]

            # Update the remap table and accumulate kept sources
            remap[divmod(fsrc, npts)] = off + within
            kept.append(ksrc)
            off += len(classes)*g.nint

        # Body vdofs: unique by construction (sequential kept ids)
        bidxs = topo.body_idxs
        bpos = ((np.arange(neles)*npts)[:, None] + bidxs).ravel()
        seq = np.arange(bpos.size).reshape(neles, len(bidxs))
        remap[:, bidxs] = off + seq
        kept.append(bpos)

        return remap, np.concatenate(kept), bpos

    def select(self, etype, values):
        _, kept, _ = self.layouts[etype]
        flat = values.swapaxes(0, 1).reshape(-1, *values.shape[2:])
        return flat[kept]

    def average(self, fields, ncomp, dtype):
        # Per-kind: accumulate sums -> exchange deficient sums -> divide by
        # precomputed topological count
        vdofvals = {etype: vals.swapaxes(0, 1).reshape(-1, ncomp)
                    for etype, vals in fields.items()}

        out = {etype: [] for etype in fields}
        for g in self.kinds.values():
            psum = np.zeros((len(g.keys)*g.nint, ncomp), dtype=dtype)
            for etype, (_, fsrc, fdst) in g.emap.items():
                np.add.at(psum, fdst, vdofvals[etype][fsrc])

            rows = self._class_rows(g.defidx, g.nint)
            vals = psum[rows].reshape(-1, g.nint, ncomp)
            psum[rows] = g.reducer(vals).reshape(-1, ncomp)

            psum /= np.repeat(g.count, g.nint)[:, None]
            for etype, (classes, _, _) in g.emap.items():
                # When this etype's classes cover every global class for
                # this kind, the fancy-index slice is the whole array
                if len(classes) == len(g.keys):
                    out[etype].append(psum)
                else:
                    rows = self._class_rows(classes, g.nint)
                    out[etype].append(psum[rows])

        # Append body slice and concat per-etype
        for etype in fields:
            _, _, body = self.layouts[etype]
            out[etype].append(vdofvals[etype][body])
            out[etype] = np.concatenate(out[etype])

        return out
