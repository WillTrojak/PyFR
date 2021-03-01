# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% rtr = 1./c['tr'] %>

<%pyfr:macro name='inviscid_flux' params='s, f'>
    // Velocity in the indices 1 to ndims+1 of the conservative variable array
    fpdtype_t v[] = ${pyfr.array('s[{i}]', i=(1, ndims + 1))};

    // Pressure in the conservative variable array index 0
    fpdtype_t p = s[0];

    // Mass flux
% for i in range(ndims):
    f[${i}][0] = ${c['ac-zeta']}*v[${i}];
% endfor

% for i, j in pyfr.ndrange(ndims, ndims):

    // Momentum fluxes
    f[${i}][${j + 1}] = -${c['nu']}*s[${1 + ndims + i + j*ndims}]
                        + v[${i}]*v[${j}] ${' + p' if i == j else ''};

    // Gradient fluxes
% for k in range(ndims):
% if k == i:
    f[${i}][${1 + ndims + k + j*ndims}] = -${rtr}*v[${j}];
% else:
    f[${i}][${1 + ndims + k + j*ndims}] = 0;
% endif
% endfor

% endfor
</%pyfr:macro>

<%pyfr:macro name='inviscid_1dflux' params='s, f'>
    fpdtype_t v[] = ${pyfr.array('s[{i}]', i=(1, ndims + 1))};

    // Mass flux
    f[0] = ${c['ac-zeta']}*v[0];

    // Momentum fluxes
    f[1] = v[0]*v[0] + s[0] - ${c['nu']}*s[${1 + ndims}];
% for j in range(2, ndims):
    f[${j + 2}] = v[0]*v[${j}] - ${c['nu']}*s[${1 + ndims + j*ndims}];
% endfor

    // Gradient fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
% if i == 0:
    f[${1 + ndims + i + j*ndims}] = -${rtr}*v[${j}];
% else:
    f[${1 + ndims + i + j*ndims}] = 0;
% endif
% endfor
</%pyfr:macro>
