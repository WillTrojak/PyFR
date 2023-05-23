<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_entropy' params='u, d, ad, p, e'>
    d = ${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rcpd = 1.0/d, E = u[${nspec + ndims}];

    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${nspec + i}];
% endfor

    ad = u[0];
% for i in range(1, nspec):
    ad = fmin(ad, u[${i}]);
% endfor

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*rcpd*${pyfr.dot('rhov[{i}]', i=ndims)};
    fpdtype_t cv_inv = 1/(${' + '.join('{cv}*u[{i}]'.format(i=i, cv=c[f'cv{i}'])
                            for i in range(nspec))});
    fpdtype_t cp = ${' + '.join('{cp}*u[{i}]'.format(i=i, cp=c[f'cp{i}'])
                     for i in range(nspec))};
    p = rhoe*(cp*cv_inv - 1);
    fpdtype_t T = rhoe*cv_inv;

    // Compute combined species entropies
    e = ${' + '.join('{cv}*u[{i}]*log(T*pow(fmax(u[{i}], {d_min}), {gam}))'.format(i=i,
          cv=c[f'cv{i}'], gam=1-c[f'cp{i}']/c[f'cv{i}'], d_min=d_min) for i in range(nspec))};
    e = e*rcpd;

    // Compute specific physical entropy
    e = ((T > 0) && (d > 0)) ? e : ${fpdtype_max};
</%pyfr:macro>
