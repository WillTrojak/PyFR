# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, f, d, p, v, g'>
    d = ${' + '.join('s[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t invrho = 1/d, E = s[${ndims + nspec}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${nspec + i}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};
    g = (${' + '.join('{cp}*s[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*s[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    p = rhoe*(g - 1);

    // Mass flux
% for i, j in pyfr.ndrange(ndims, nspec):
    f[${i}][${j}] = v[${i}]*s[${j}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + nspec}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor

    // Energy fluxes
% for i in range(ndims):
    f[${i}][${ndims + nspec}] = (E + p)*v[${i}];
% endfor

</%pyfr:macro>
