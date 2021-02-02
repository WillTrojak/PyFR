# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.backends.base.makocommon.transform'/>
<%include file='pyfr.backends.base.makocommon.transform_grad'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];

    utl[0] = ul[0]; utr[0] = ur[0];
    ${pyfr.expand('transform_to','n', 'ul', 'utl', '1')};
    ${pyfr.expand('transform_to','n', 'ur', 'utr', '1')};

% if ndims == 2:
    ${pyfr.expand('transform_grad_to', 'n', 'ul', 'utl', '4')};
    ${pyfr.expand('transform_grad_to', 'n', 'ur', 'utr', '4')};
% elif ndims == 3:
    ${pyfr.expand('transform_grad_to', 'n', 'ul', 'utl', '5')};
    ${pyfr.expand('transform_grad_to', 'n', 'ur', 'utr', '5')};
% endif

    ${pyfr.expand('rsolve_t1d', 'utl', 'utr', 'ntf')};

    nf[0] = ntf[0];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', '1')};
% if ndims == 2:
    ${pyfr.expand('transform_grad_from', 'n', 'ntf', 'nf', '4')};
% elif ndims == 3:
    ${pyfr.expand('transform_grad_from', 'n', 'ntf', 'nf', '5')};
% endif
</%pyfr:macro>
