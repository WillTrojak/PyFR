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

    // Fraction
% for i, j in pyfr.ndrange(ndims, nspec-1):
    f[${i}][${j + nspec + ndims + 1}] = v[${i}]*s[${j + nspec + ndims + 1}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='s, f, d, p, v, g'>
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
% for j in range(nspec):
    f[${j}] = v[0]*s[${j}];
% endfor

    // Momentum fluxes
    f[${0 + nspec}] = rhov[0]*v[0] + p;
% for j in range(1, ndims):
    f[${j + nspec}] = rhov[0]*v[${j}];
% endfor

    // Energy fluxes
    f[${ndims + nspec}] = (E + p)*v[0];

    // Fraction flux
% for i in range(nspec - 1):
    f[${i + nspec + ndims + 1}] = v[0]*s[${i + nspec + ndims + 1}];
% endfor

</%pyfr:macro>
