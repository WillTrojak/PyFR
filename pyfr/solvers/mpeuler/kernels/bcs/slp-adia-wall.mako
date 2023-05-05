<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};
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
