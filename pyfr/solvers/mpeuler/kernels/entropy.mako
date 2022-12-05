# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% inf = 1e20 %>
<%pyfr:macro name='compute_entropy' params='u, ad, amin, amax, T, e'>
    fpdtype_t d = ${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rcpd = 1.0/d, E = u[${ndims + nspec}];

    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = u[${nspec + i}];
% endfor

    amin = ${inf};
    amax = ${-inf};
    ad = ${inf};
    fpdtype_t a_acc = 1;
% for i in range(nspec-1):
    amin = min(amin, u[${nspec + ndims + 1 + i}]);
    amax = max(amax, u[${nspec + ndims + 1 + i}]);
    ad = min(ad, u[${i}]);
    a_acc -= u[${nspec + ndims + 1 + i}];
% endfor
    amin = min(amin, a_acc);
    amax = max(amax, a_acc);
    amax = 1 - amax;

    // Compute the pressure
    fpdtype_t rhoe = E - 0.5*rcpd*${pyfr.dot('rhov[{i}]', i=ndims)};
    //fpdtype_t rcp = (${' + '.join('{cp}*u[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))});
    fpdtype_t rcv_inv = 1/(${' + '.join('{cv}*u[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    //fpdtype_t gamma = rcp*rcv_inv;
    T = rhoe*rcv_inv;

    // Compute combined species entropies
    e = ${' + '.join('{cv}*u[{i}]*log(T*pow(u[{i}], {gam}))'.format(i=i, cv=c[f'cv{i}'], gam=1-c[f'cp{i}']/c[f'cv{i}']) for i in range(nspec))};
    e = e*rcpd;

    // Compute specific physical entropy
    // e = ((T > 0) && (amin > 0) && (amax > 0) && (ad > 0)) ? exp(e) : ${inf};
    e = ((T > 0) && (amin > 0) && (amax > 0) && (ad > 0)) ? e : ${inf};
</%pyfr:macro>
