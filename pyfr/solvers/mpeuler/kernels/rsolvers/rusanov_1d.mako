<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve_1d' params='ul, ur, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr, gl, gr;

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'rhol', 'pl', 'vl', 'gl')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'rhor', 'pr', 'vr', 'gr')};

    // wave speeds
    fpdtype_t cl = sqrt(gl*pl/rhol);
    fpdtype_t cr = sqrt(gr*pr/rhor);
    fpdtype_t a = max(fabs(vl[0]) + cl, fabs(vr[0]) + cr);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(fl[${i}] + fr[${i}]) + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.rsolve1d'/>
