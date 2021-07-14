# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.edeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    fpdtype_t vl[${ndims}] = ${pyfr.array('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t vr[${ndims}] = ${pyfr.array('ur[{i}]', i=(1, ndims + 1))};

    // Normal of the average interface velocity
    fpdtype_t nv = 0.5*${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Mean Pressure
    fpdtype_t pavg = 0.5*(ul[0] + ur[0]);

    // Estimate the wave speed
    fpdtype_t a = 1.5*fabs(nv) + sqrt(0.5*nv*nv + ${c['ed-zeta']} + pavg);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
