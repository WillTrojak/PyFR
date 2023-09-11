<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.transform_grad'/>

<%pyfr:kernel name='gradcoru' ndim='2'
              gradu='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'>
    ${pyfr.expand('transform_grad', 'gradu', 'smats', 'rcpdjac')};
</%pyfr:kernel>
