# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.${source_sys}.kernels.${source}'/>

<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

    fpdtype_t rcpdjac_c = ${(32./3.14159265359)**3};

% for i, ex in enumerate(srcex):
    tdivtconf[${i}] = -rcpdjac_c*tdivtconf[${i}] + ${ex};
% endfor
    ${pyfr.expand('sources', 'tdivtconf', 'u', 'ploc')};
</%pyfr:kernel>
