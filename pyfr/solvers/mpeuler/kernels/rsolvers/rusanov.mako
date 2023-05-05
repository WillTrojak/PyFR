<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr, gl, gr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'rhol', 'pl', 'vl', 'gl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'rhor', 'pr', 'vr', 'gr')};

    // Sum the left and right velocities and take the normal
    //fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    //fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};

    // wave speeds
    //pdtype_t cl = sqrt(gl*pl/rhol);
    //fpdtype_t cr = sqrt(gr*pr/rhor);
    //fpdtype_t c = max(fabs(nvl) + cl, fabs(nvr) + cr);
    // Sum the left and right velocities and take the normal

    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed
    fpdtype_t a = sqrt(0.5*(gl + gr)*(pl + pr)/(rhol + rhor)) + 0.5*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
