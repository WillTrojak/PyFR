<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconf_mp_aa' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              mat='in fpdtype_t[${str(nspec-1)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
    fpdtype_t inv_rho = 1/(${' + '.join('u[{i}]'.format(i=i) for i in range(nspec))});
    
% for i, ex in enumerate(srcex):
% if i <= nspec + ndims:
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex};
% else:
    tdivtconf[${i}] = -mat[${i-nspec-ndims-1}] + ${ex};
% endif
% endfor
</%pyfr:kernel>
