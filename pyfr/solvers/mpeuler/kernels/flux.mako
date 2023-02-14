<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% t_kmax = 3 %>

<%pyfr:macro name='temperature' params='q, rho, rhoe, T, p'>
    fpdtype_t rhoe_dash, rhoh, de;

    p = ${s['Rconst']}*T*(${' + '.join(f'q[{ndims + 1 + i}]' for i in range(nspec))});

    // Newton solve for temperature
% for k in range(t_kmax):

    rhoh = ${s['R0']*s['M0']}*q[${ndims + 1}]*${pyfr.polyval(s['Hr0'], 'T')};
% for i in range(1, nspec):
    rhoh += ${s[f'R{i}']*s[f'M{i}']}*q[${ndims + 1 + i}]*${pyfr.polyval(s[f'Hr{i}'], 'T')};
% endfor
    rhoe_dash = rhoh - p;

    de = ${s[f'R0']*s[f'M0']}*q[${ndims + 1}]*${pyfr.polyval(s[f'Cv0'], 'T')};
% for i in range(1, nspec):
    de += ${s[f'R{i}']*s[f'M{i}']}*q[${ndims + 1 + i}]*${pyfr.polyval(s[f'Cv{i}'], 'T')};
% endfor

    T += (rhoe - rhoe_dash) / de;
    p = ${s['Rconst']}*T*(${' + '.join(f'q[{ndims + 1 + i}]' for i in range(nspec))});
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux' params='q, f, p, T, rho, v'>
    rho = ${' + '.join(f"({s[f'M{i}']}*q[{ndims + 1 + i}])" for i in range(nspec))};
    fpdtype_t invrho = 1/rho, E = q[${ndims}];

    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('q[{i}]', i=ndims)};

    // Linear initial guess for temperature based T = e / cv(298.15K)
    T = rhoe / (${' + '.join(f"{s[f'Cv_ref{i}']*s[f'M{i}']}*q[{ndims + 1 + i}]" for i in range(nspec))});
    ${pyfr.expand('temperature', 'q', 'rho', 'rhoe', 'T', 'p')};

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = q[${i}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j}] = rhov[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor

    // Energy flux
% for i in range(ndims):
    f[${i}][${ndims}] = (E + p)*v[${i}];
% endfor

    // Concentration flux
% for i, j in pyfr.ndrange(ndims, nspec):
    f[${i}][${ndims + 1 + j}] = v[${i}]*q[${ndims + 1 + j}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='q, f, p, T, rho, v'>
    rho = ${' + '.join(f"({s[f'M{i}']}*q[{ndims + 1 + i}])" for i in range(nspec))};
    fpdtype_t invrho = 1/rho, E = q[${ndims}];

    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('q[{i}]', i=ndims)};

    // Linear initial guess for temperature based T = e / cv(298.15K)
    T = rhoe / (${' + '.join(f"{s[f'Cv_ref{i}']*s[f'M{i}']}*q[{ndims + 1 + i}]" for i in range(nspec))});

    ${pyfr.expand('temperature', 'q', 'rho', 'rhoe', 'T', 'p')};

    // Compute the velocities
% for i in range(ndims):
    v[${i}] = invrho*q[${i}];
% endfor

    // Momentum fluxes
    f[0] = q[0]*v[0] + p;
% for i in range(1, ndims):
    f[${i}] = q[0]*v[${i}];
% endfor

    // Energy flux
    f[${ndims}] = v[0]*(E + p);

    // Concentration flux
% for i in range(ndims + 1, ndims + 1 + nspec):
    f[${i}] = v[0]*q[${i}];
% endfor
</%pyfr:macro>
