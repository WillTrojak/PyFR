# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='grad_array',params='s,q'>

% for i,j in pyfr.ndrange(ndims, ndims):
    q[${i}][${j}] = s[${ndims + 1 + i + (j*ndims)}];
% endfor

</%pyfr:macro>

<% rtr = 1./c['tr'] %>

<%pyfr:macro name='inviscid_flux' params='s, f'>
    // Velocity in the indices 1 to ndims+1 of the conservative variable array
    fpdtype_t v[] = ${pyfr.array('s[{i}]', i=(1, ndims + 1))};
    fpdtype_t q[${ndims}][${ndims}];

    ${pyfr.expand('grad_array','s','q')};

    // Pressure in the conservative variable array index 0
    fpdtype_t p = s[0];

    // Mass flux
% for i in range(ndims):
    f[${i}][0] = ${c['ac-zeta']}*v[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    // Momentum fluxes
    f[${i}][${j + 1}] = -${c['nu']}q[${i}][${j}] + v[${i}]*v[${j}]${' + p' if i == j else ''};
    
    // Gradient fluxes
    f[${i}][${i + 1 + ndims + j*ndims}] = ${rtr}*v[${j}];
% endfor

   // Gradient fluxes

</%pyfr:macro>
