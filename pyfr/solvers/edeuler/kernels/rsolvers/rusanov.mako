# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.edeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    // Normal of the average interface velocity
    fpdtype_t unl = ${pyfr.dot('n[{i}]', 'ul[1+{i}]', i=ndims)};
    fpdtype_t unr = ${pyfr.dot('n[{i}]', 'ur[1+{i}]', i=ndims)};

    // Estimate the wave speed
    fpdtype_t a = max(1.5*fabs(unl) + sqrt(0.25*unl*unl + ${c['ed-zeta']} + ul[0]), 
                      1.5*fabs(unr) + sqrt(0.25*unr*unr + ${c['ed-zeta']} + ur[0])); 

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
