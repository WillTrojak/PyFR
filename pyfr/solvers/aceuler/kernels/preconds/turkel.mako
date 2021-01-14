# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='localprecond' params='u, fi, fo'>

	fo[0] = fi[0];
% for i in range(ndims):
	fo[${i + 1}] = fi[${i + 1}] - ${c['ac-alpha']}*u[${i + 1}]*fi[0];
% endfor

</%pyfr:macro>