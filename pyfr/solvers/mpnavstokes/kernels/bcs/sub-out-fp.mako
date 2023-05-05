<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nvars - 1):
    ur[${i}] = ul[${i}];
% endfor
    fpdtype_t g, invrho;

    invrho = 1/${' + '.join('ul[{i}]'.format(i=i) for i in range(nspec))};
    g = (${' + '.join('{cp}*ul[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
        (${' + '.join('{cv}*ul[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    ur[${nvars - 1}] = ${c['p']}/(g - 1);
                     + 0.5*invrho*${pyfr.dot('ul[{i}]', i=(spec, ndims + nspec))};
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
