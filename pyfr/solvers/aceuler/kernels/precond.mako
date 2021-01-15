# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.aceuler.kernels.preconds.${precond}'/>

<%pyfr:kernel name='precond' ndim='2'
			  u='in fpdtype_t[${str(nvars)}]'
              f='inout fpdtype_t[${str(nvars)}]'>
    fpdtype_t f_c[${nvars}];

% for i in range(nvars):
	f_c[${i}] = f[${i}];
% endfor
	${pyfr.expand('localprecond', 'u', 'f_c', 'f')};	
</%pyfr:kernel>
