# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.hypenavstokes.kernels.flux'/>

<%pyfr:kernel name='divsparse' ndim='1'
              df='inout fpdtype_t[${str(nvars)}]'>
    // Compute the flux
    fpdtype_t ftemp[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'smats[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
