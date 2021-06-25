# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};
    fpdtype_t srl = sqrt(ul[0]), srr = sqrt(ur[0]);

    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}]*srl + vr[{i}]*srr', i=ndims)};

    // Estimate the maximum wave speed
    fpdtype_t h = ((ul[${nvars-1}] + pl)*srr + (ur[${nvars-1}] + pr)*srl) / 
                  ((srl + srr)*srl*srr);
    fpdtype_t a = sqrt(${gamma-1}*(h - 0.5*nv*nv));

    fpdtype_t l = max(abs(u - a), abs(u + a));

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*l*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
