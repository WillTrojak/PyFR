<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='ugrad_copy' ndim='2'
              grad='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              ugrad='out fpdtype_t[${str(ndims)}][${str(nspec + ndims)}]'>
% for i, j in pyfr.ndrange(ndims, nspec + ndims):
    ugrad[${i}][${j}] = grad[${i}][${j}];
% endfor
</%pyfr:kernel>
