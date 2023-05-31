<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

% for i in range(nspec):
    utl[${i}] = ul[${i}];
    utr[${i}] = ur[${i}];
% endfor
    utl[${nvars - 1}] = ul[${nvars - 1}];
    utr[${nvars - 1}] = ur[${nvars - 1}];

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', str(nspec))};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', str(nspec))};

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ntf')};

% for i in range(nspec):
    nf[${i}] = ntf[${i}];
% endfor
    nf[${nvars - 1}] = ntf[${nvars - 1}];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', str(nspec))};
</%pyfr:macro>
