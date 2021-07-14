# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.aceuler.kernels.flux'/>
<%include file='pyfr.solvers.acnavstokes.kernels.flux'/>

<%pyfr:kernel name='tflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'>
    // Compute the flux (F = Fi + Fv)
    fpdtype_t ftemp[${ndims}][${nvars}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp')};
    ${pyfr.expand('viscous_flux_add', 'u', 'f', 'ftemp')};

    fpdtype_t smats_c[${ndims}][${ndims}];
% for i, j in pyfr.ndrange(ndims, ndims):
    smats_c[${i}][${j}] = 0.;
% endfor
    smats_c[0][0] = ${(3.1415926/32.)**2};
    smats_c[1][1] = ${(3.1415926/32.)**2};
    smats_c[2][2] = ${(3.1415926/32.)**2};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'smats_c[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
