# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.hypenavstokes.kernels.flux'/>

<% nurtr = c['nu']/c['tr'] %>
<% rtr = 1/c['tr'] %>
<% znurtr = c['ac-zeta'] + c['nu']/c['tr'] %>

<%pyfr:macro name='rsolve_t1d' params='ql, qr, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];

    ${pyfr.expand('inviscid_1dflux', 'ql', 'fl')};
    ${pyfr.expand('inviscid_1dflux', 'qr', 'fr')};

    fpdtype_t s = max(fabs(ql[1]) + sqrt(ql[1]*ql[1] + ${c['ac-zeta'] + nurtr}),
                      fabs(qr[1]) + sqrt(qr[1]*qr[1] + ${c['ac-zeta'] + nurtr}));
    
    // Output
% for i in range(ndims + 1):
    nf[${i}] = 0.5*(fl[${i}] + fr[${i}]) + 0.5*s*(ql[${i}] - qr[${i}]);
% endfor

% for i in range(ndims + 1, nvars):
nf[${i}] = 0.5*(fl[${i}] + fr[${i}]);
% endfor
    
</%pyfr:macro>

<%include file='pyfr.solvers.hypenavstokes.kernels.rsolvers.rsolve_trans'/>
