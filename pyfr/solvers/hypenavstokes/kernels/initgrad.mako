# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='initgrad' ndim='2'
              gradu='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
           	  u='inout fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'>
    fpdtype_t tmpgradu[${ndims}];

% for j in range(ndims):
% for i in range(ndims):
    tmpgradu[${i}] = gradu[${i}][${j+1}];
% endfor
% for i in range(ndims):
	u[${1 + i + ndims + j*ndims}] = rcpdjac*(${' + '.join(f'smats[{k}][{i}]*tmpgradu[{k}]'
                                             for k in range(ndims))});
% endfor
% endfor
</%pyfr:kernel>