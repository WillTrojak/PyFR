# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.hypenavstokes.kernels.flux'/>

<% nurtr = c['nu']/c['tr'] %>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t a[${nvars}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    fpdtype_t vl[${ndims}] = ${pyfr.array('ul[{i}]', i=(1, ndims + 1))};
    fpdtype_t vr[${ndims}] = ${pyfr.array('ur[{i}]', i=(1, ndims + 1))};

    // Normal of the average interface velocity
    fpdtype_t nv = 0.5*${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the wave speed
    // fpdtype_t s = fabs(nv) + sqrt(nv*nv + ${c['ac-zeta'] + nurtr});
    fpdtype_t s = fabs(nv) + sqrt(nv*nv + ${c['ac-zeta']});
% for i in range(nvars):
% if i < ndims + 1:
    a[${i}] = s;
% else:
    a[${i}] = 0;
% endif
% endfor

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a[${i}]*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
