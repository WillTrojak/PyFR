<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='fgrad_copy' ndim='2'
              grad='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              fgrad='out fpdtype_t[${str(ndims)}][${str(nspec-1)}]'>
% for i, j in pyfr.ndrange(ndims, nspec-1):
    fgrad[${i}][${j}] = grad[${i}][${j + nspec + ndims + 1}];
% endfor
</%pyfr:kernel>
