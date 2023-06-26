<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.flux'/>

<%pyfr:kernel name='tflux_vis' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              artvisc='in broadcast-col fpdtype_t'
              f='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'>
    // Compute the flux (F = Fi + Fv)
    fpdtype_t ftemp[${ndims}][${nvars}];
    fpdtype_t d, p, v[${ndims}], g;
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'd', 'p', 'v', 'g')};
    ${pyfr.expand('viscous_flux_add', 'u', 'f', 'ftemp')};
    ${pyfr.expand('artificial_viscosity_add', 'f', 'ftemp', 'artvisc')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, nvars):
    f[${i}][${j}] = ${' + '.join(f'smats[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};
% endfor
</%pyfr:kernel>
