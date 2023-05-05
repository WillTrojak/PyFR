<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nspec):
    ur[${i}] = ${c[f'a{i}rho{i}']};
% endfor

    fpdtype_t rho_r = ${'+'.join(f"{c[f'a{i}rho{i}']}" for i in range(nspec))};
% for i, v in enumerate('uvw'[:ndims]):
    ur[${i + nspec}] = rho_r * (${c[v]});
% endfor
    fpdtype_t invrho, g;

    invrho = 1/${' + '.join('ul[{i}]'.format(i=i) for i in range(nspec))};

    ur[${nvars - 1}] = ul[${nvars - 1}]
                     - 0.5*invrho*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))}
                     + 0.5*(1 / rho_r)*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
