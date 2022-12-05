<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
   // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr, gl, gr;
    fpdtype_t nf_sub, nf_fl, nf_fr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'rhol', 'pl', 'vl', 'gl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'rhor', 'pr', 'vr', 'gr')};

    // Get the normal left and right velocities
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    fpdtype_t cl = sqrt(gl*pl/rhol);
    fpdtype_t cr = sqrt(gr*pr/rhor);

    // Estimate the left and right wave speed, sl and sr
    fpdtype_t sl = min(nvl - cl, nvr - cr);
    fpdtype_t sr = max(nvl + cl, nvr + cr);
    fpdtype_t rcpsrsl = 1/(sr - sl);

    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join('n[{j}]*fl[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
    nf_fr = ${' + '.join('n[{j}]*fr[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
    nf_sub = (sr*nf_fl - sl*nf_fr + sl*sr*(ur[${i}] - ul[${i}]))*rcpsrsl;
    nf[${i}] = (0 <= sl) ? nf_fl : (0 >= sr) ? nf_fr : nf_sub;
% endfor
</%pyfr:macro>