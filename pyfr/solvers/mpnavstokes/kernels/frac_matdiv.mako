<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='frac_matdiv' ndim='2'
              fgrad='in fpdtype_t[${str(ndims)}][${str(nspec-1)}]'
              u='in fpdtype_t[${str(nvars)}]'
              mat='out fpdtype_t[${str(nspec-1)}]'>
    fpdtype_t inv_rho = 1/(${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))});
% for i in range(nspec-1):
    mat[${i}] = inv_rho*(${' + '.join('u[{k}]*fgrad[{j}][{i}]'.format(k=nspec+j, j=j, i=i) 
                           for j in range(ndims))});
% endfor
</%pyfr:kernel>
