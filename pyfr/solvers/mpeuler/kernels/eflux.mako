<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:kernel name='etflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'>

    // Compute the density
    fpdtype_t d = ${' + '.join('s[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t invrho = 1/d;

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
    fpdtype_t p, E = u[${ndims}];
    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('u[{i}]', i=ndims)};
    g = (${' + '.join('{cp}*s[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*s[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    p = rhoe*(g - 1);

    // Add the additional internal energy term to the energy flux
% for i in range(ndims):
    f[${i}][${ndims}] += ${' + '.join(f'dgdar[{j}]*f[{i}][{j}]' for j in range(nspec))};
% endfor

</%pyfr:kernel>