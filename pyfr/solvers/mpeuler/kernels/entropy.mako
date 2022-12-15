<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% inf = 1e20 %>
<%pyfr:macro name='compute_entropy' params='u, ad, T, e'>
    fpdtype_t d = ${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rcpd = 1.0/d, E = u[${nspec + ndims}];

    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${nspec + i}];
% endfor

    ad = u[0];
% for i in range(nspec-1):
    ad = min(ad, u[${i+1}]);
% endfor

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*rcpd*${pyfr.dot('rhov[{i}]', i=ndims)};
    fpdtype_t cv = ${' + '.join('{cv}*u[{i}]'.format(i=i, cv=c[f'cv{i}']) 
                     for i in range(nspec))};
    T = rhoe/cv;

    // Compute combined species entropies
    e = ${' + '.join('{cv}*u[{i}]*log(T*pow(u[{i}], {gam}))'.format(i=i, 
          cv=c[f'cv{i}'], gam=1-c[f'cp{i}']/c[f'cv{i}']) for i in range(nspec))};
    e = e*rcpd;

    // Compute specific physical entropy
    e = ((T > 0) && (ad > 0)) ? e : ${inf};
</%pyfr:macro>
