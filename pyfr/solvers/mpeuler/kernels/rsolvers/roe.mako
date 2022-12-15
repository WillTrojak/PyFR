<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve_1d' params='ul, ur, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${nvars}], fr[${nvars}], Adu[${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t El = ul[${nvars - 1}], Er = ur[${nvars - 1}];
    fpdtype_t rhol, rhor, pl, pr, gl, gr;
    fpdtype_t v_roe[${ndims}], arho_roe[${nspec}];

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'rhol', 'pl', 'vl', 'gl')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'rhor', 'pr', 'vr', 'gr')};

    fpdtype_t sqr_rhol = sqrt(rhol);
    fpdtype_t sqr_rhor = sqrt(rhor);

    fpdtype_t inv_rho_roe = 1 / (sqr_rhol + sqr_rhor);
% for i in range(ndims):
    v_roe[${i}] = (sqr_rhol*vl[${i}] + sqr_rhor*vr[${i}]) * inv_rho_roe;
% endfor
    fpdtype_t u_abs = fabs(v_roe[0]);
    
% for i in range(nspec):
    arho_roe[${i}] = (sqr_rhol*ul[${i}] + sqr_rhor*ur[${i}]) * inv_rho_roe;
% endfor
    fpdtype_t g_roe = (${' + '.join('arho_roe[{i}]*{cp}'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) / 
                      (${' + '.join('arho_roe[{i}]*{cv}'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});

    fpdtype_t h_roe = (sqr_rhol*(El/rhol + pl) + sqr_rhor*(Er/rhor + pr)) 
                       * inv_rho_roe;
    //fpdtype_t p_roe = ((g_roe - 1)/g_roe)*(h_roe - 0.5*${pyfr.dot('v_roe[{i}]', i=ndims)})*sqrt(rhol*rhor);

    fpdtype_t p_roe = (sqr_rhol*pl + sqr_rhor*pr) * inv_rho_roe;

    // Fraction mass
% for i in range(nspec):
    Adu[${i}] = u_abs*(ur[${i}] - ul[${i}]);
% endfor

    // Momentum
    fpdtype_t drho = rhor - rhol;
% for i in range(ndims):
    Adu[${i + nspec}] = u_abs*drho*v_roe[${i}];
% endfor

    // Energy
    Adu[${nspec + ndims}] = 0.5*u_abs*v_roe[0]*v_roe[0]*drho + p_roe*(1/(gr - 1) - 1/(gl - 1));

% for i in range(nvars):
    nf[${i}] = 0.5*(fl[${i}] + fr[${i}]) - 0.5*Adu[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.rsolve1d'/>
