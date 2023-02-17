<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='eflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='inout fpdtype_t[${str(nvars)}]'>

    // Compute the density
    fpdtype_t d = ${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t invrho = 1/d;

    // Compute d(G)/dairhoi
    fpdtype_t dgdar[${nspec}];
    fpdtype_t arr_inv = ${' + '.join(f"{c[f'cp{j}'] - c[f'cv{j}']}*u[{j}]" for j in range(nspec))};
    arr_inv = 1 / (arr_inv*arr_inv);
% for i in range(nspec):
    dgdar[${i}] = (${' + '.join('{b}*u[{j}]'.format(
                     b=(c[f'cv{i}']*(c[f'cp{j}'] - c[f'cv{j}']) -
                        c[f'cv{j}']*(c[f'cp{i}'] - c[f'cv{i}'])), j=j)
                     for j in range(nspec) if j != i)})*arr_inv;
% endfor

    // Compute the pressure
    fpdtype_t p, g, E = u[${ndims}];
    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('u[{i}]', i=ndims)};
    g = (${' + '.join('{cp}*u[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*u[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    p = rhoe*(g - 1);

    // Add the additional internal energy term to the energy flux
    f[${nvars - 1}] += ${' + '.join(f'p*dgdar[{j}]*f[{j}]' for j in range(nspec))};

</%pyfr:kernel>