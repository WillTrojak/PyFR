# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% rtr = 1/c['tr'] %>

<%pyfr:macro name='sources' params='df, u, x'>
	
% for i in range(ndims*ndims):
	df[${i + 1 + ndims}] -= u[${i + 1 + ndims}]*${rtr};
% endfor
</%pyfr:macro>