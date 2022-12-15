<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.bcs.common'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.flux'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur'>
    fpdtype_t nor = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + nspec)
                                 for i in range(ndims))};
    // Mass
% for i in range(nspec):
    ur[${i}] = ul[${i}];
% endfor
    // Momentum
% for i in range(ndims):
    ur[${i + nspec}] = ul[${i + nspec}] - 2*nor*nl[${i}];
% endfor
    // Energy
    ur[${nspec + ndims}] = ul[${nspec + ndims}];
</%pyfr:macro>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>
    // Ghost state r
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

% for i in range(nvars):
    ul[${i}] = magnl*(ficomm[${i}]);
% endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
