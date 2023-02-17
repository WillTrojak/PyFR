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

    // Compute d(G)/dairhoi
    fpdtype_t dgdar[${nspec}];
    fpdtype_t arr_inv = ${' + '.join(f"{c[f'cp{j}'] - c[f'cv{j}']}*s[{j}]" for j in range(nspec))};
    arr_inv = 1 / (arr_inv*arr_inv);
% for i in range(nspec):
    dgdar[${i}] = (${' + '.join('{b}*s[{j}]'.format(
                     b=(c[f'cv{i}']*(c[f'cp{j}'] - c[f'cv{j}']) -
                        c[f'cv{j}']*(c[f'cp{i}'] - c[f'cv{i}'])), j=j)
                     for j in range(nspec) if j != i)})*arr_inv;
% endfor

    // Compute the pressure
    fpdtype_t ek = 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};
    //fpdtype_t rhoe = E - ek;
    g = (${' + '.join('{cp}*s[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*s[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    p = (E - ek)*(g - 1);

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
    f[${i}][${nvars - 1}] = (ek + p)*v[${i}];
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
    fpdtype_t ek = 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)};
    g = (${' + '.join('{cp}*s[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*s[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    p = (E - ek)*(g - 1);

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
    f[${nvars - 1}] = (ek + p)*v[0];
</%pyfr:macro>
