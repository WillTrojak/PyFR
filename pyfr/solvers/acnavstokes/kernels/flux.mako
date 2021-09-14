# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
% if c['ac-alpha'] > 0.:
% for i in range(ndims):
    fout[${i}][0] = -${c['nu']*c['ac-alpha']}*(grad_uin[${i}][0] - 
                     ${' + '.join('uin[{j}]*grad_uin[{i}][{j}]'.format(i=i, j=j+1) for j in range(ndims))});
% endfor
% endif

% for i, j in pyfr.ndrange(ndims, ndims):
    fout[${i}][${j+1}] += -${c['nu']}*grad_uin[${i}][${j+1}];
% endfor
</%pyfr:macro>
