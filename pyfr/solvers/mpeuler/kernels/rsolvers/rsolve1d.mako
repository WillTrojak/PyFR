<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.transform'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

% for i in range(ndims, ndims + nspec + 1):
    utl[${i}] = ul[${i}];
    utr[${i}] = ur[${i}];
% endfor
    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', '0')};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', '0')};

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ntf')};

% for i in range(ndims, ndims + nspec + 1):
    nf[${i}] = ntf[${i}];
% endfor
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', '0')};
</%pyfr:macro>
