# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>

    ur[0] = ul[0];
% for i, v in enumerate(c['v']):
    ur[${i + 1}] = -ul[${i + 1}] + ${2*v};
% endfor

% for i in range(1+ndims, nvars):
	ur[${i}] = ul[${i}];
% endfor

</%pyfr:macro>