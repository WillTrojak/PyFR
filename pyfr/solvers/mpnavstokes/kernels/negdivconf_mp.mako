<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf_mp' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              grad='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
    fpdtype_t inv_rho = 1/(${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))});
    fpdtype_t divrhou = ${' + '.join('grad[{i}][{j}]'.format(i=i, j=i + nspec) 
                          for i in range(ndims))};
% for i in range(ndims):
% if i == 0:
    fpdtype_t ugradrho = u[${nspec}]*(${' + '.join('grad[0][{j}]'.format(j=j) 
                                        for j in range(nspec))});
% else:
    ugradrho += u[${nspec + i}]*(${' + '.join('grad[{i}][{j}]'.format(i=i, j=j) 
                                   for j in range(nspec))});
% endif
% endfor
    fpdtype_t divu = (divrhou - ugradrho*inv_rho)*inv_rho;

% for i, ex in enumerate(srcex):
% if i <= nspec + ndims:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
% else:
    tdivtconf[${i}] = -(rcpdjac*tdivtconf[${i}] - inv_rho*u[${i}]*divu) + ${ex};
% endif
% endfor
</%pyfr:kernel>
