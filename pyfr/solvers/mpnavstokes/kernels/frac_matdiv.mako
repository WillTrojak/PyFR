<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='frac_matdiv' ndim='2'
              ugrad='in fpdtype_t[${str(ndims)}][${str(nspec+ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              mat='out fpdtype_t[${str(nspec-1)}]'>
    fpdtype_t inv_rho = 1/(${' + '.join('u[{i}]'.format(i=i) 
                             for i in range(nspec))});
    fpdtype_t divrhou = ${' + '.join('ugrad[{i}][{j}]'.format(i=i, j=i + nspec) 
                          for i in range(ndims))};
    
    fpdtype_t ugradrho = 0;
% for i in range(ndims):
    ugradrho += u[${nspec + i}]*(${' + '.join('ugrad[{i}][{j}]'.format(i=i, j=j) 
                                  for j in range(nspec))});
% endfor
    ugradrho *= inv_rho;
    fpdtype_t divu = (divrhou - ugradrho)*inv_rho;

% for i in range(nspec-1):
    mat[${nspec - 1 - i}] = u[${nvars - i}]*divu;
% endfor
</%pyfr:kernel>
