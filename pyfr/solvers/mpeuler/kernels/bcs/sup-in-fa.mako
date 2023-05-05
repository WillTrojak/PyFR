<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nspec):
    ur[${i}] = ${c[f'a{i}rho{i}']};
% endfor

    fpdtype_t rho = ${'+'.join('{x}'.format(x=c[f'a{i}rho{i}']) for i in range(nspec))};
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + nspec}] = rho*(${c[v]});
% endfor
    fpdtype_t gamma = (${'+'.join('{x}*{y}'.format(x=c[f'a{i}rho{i}'], y=c[f'cp{i}']) for i in range(nspec))}) /
                      (${'+'.join('{x}*{y}'.format(x=c[f'a{i}rho{i}'], y=c[f'cv{i}']) for i in range(nspec))});
    ur[${nvars - 1}] = ${c['p']}/(gamma - 1) +
                       0.5*(1.0/rho)*${pyfr.dot('ur[{i}]', i=(nspec, ndims + nspec))};
</%pyfr:macro>
