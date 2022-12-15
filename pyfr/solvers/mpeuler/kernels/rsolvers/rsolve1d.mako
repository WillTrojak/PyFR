<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.backends.base.mako.transform'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

% for i in range(nspec):
    utl[${i}] = ul[${i}];
    utr[${i}] = ur[${i}];
% endfor
    utl[${nspec + ndims}] = ul[${nspec + ndims}];
    utr[${nspec + ndims}] = ur[${nspec + ndims}];

    ${pyfr.expand('transform_to', 'n', 'ul', 'utl', str(nspec))};
    ${pyfr.expand('transform_to', 'n', 'ur', 'utr', str(nspec))};

    ${pyfr.expand('rsolve_1d', 'utl', 'utr', 'ntf')};

% for i in range(nspec):
    nf[${i}] = ntf[${i}];
% endfor
    nf[${nspec + ndims}] = ntf[${nspec + ndims}];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', str(nspec))};
</%pyfr:macro>
