<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='enthalpy' params='q, T, h'>

</%pyfr:macro>

<%pyfr:macro name='temperature' params='q, rho, rhoe, T, p, t_kmax'>
    fpdtype_t rhoe_dash, rhoh, de;

    p = ${s['Rconst']}*T*(${' + '.join(f'q[{ndims + 1 + i}]' for i in range(nspec))});

    // Newton solve for temperature
% for k in range(t_kmax):

    rhoh = ${s['R0']*s['M0']}*q[${ndims + 1}]*${pyfr.polyval(s['Hr0'], 'T')};
% for i in range(1, nspec):
    rhoh += ${s[f'R{i}']*s[f'M{i}']}*q[{ndims + 1 + i}]*${pyfr.polyval(s[f'Hr{i}'], 'T')};
% endfor
    rhoe_dash = rhoh - p

    de = ${s[f'R0']*s[f'M0']}*q[${ndims + 1}]*${pyfr.polyval(s[f'Cv0'], 'T')};
% for i in range(1, nspec):
    de += ${s[f'R{i}']*s[f'M{i}']}*q[${ndims + 1 + i}]*${pyfr.polyval(s[f'Cv{i}'], 'T')};
% endfor

    T += (rhoe - rhoe_dash) / de;
    p = ${s['Rconst']}*T*(${' + '.join(f'q[{ndims + 1 + i}]' for i in range(nspec))});
% endfor
</%pyfr:macro>

<%pyfr:macro name='inviscid_flux' params='s, f, p, T, rho, v'>
    rho = ${' + '.join(f"({s[f'M{i}']}*s[{ndims + 1 + i}])" for i in range(nspec))};
    fpdtype_t invrho = 1.0/rho, E = s[${nvars - 1}];

    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('s[{i}]', i=ndims)};

    // Linear initial guess for temperature based T = e / cv(298.15K)
    T = rhoe / (${' + '.join(f"{s[f'Cv{i}_ref']*s[f'M{i}']}*s[{ndims + 1 + i}]") for i in range(nspec)});

    ${pyfr.expand('temperature', 'rho', 'rhoe', 'T', 'p', t_kmax=3)};

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
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
    f[${i}][${ndims + 1 + j}] = v[${i}]*s[${ndims + 1 + j}];
% endfor

</%pyfr:macro>

<%pyfr:macro name='inviscid_flux_1d' params='s, f, p, T, rho, v'>
    rho = ${' + '.join(f"({s[f'M{i}']}*s[{ndims + 1 + i}])" for i in range(nspec))};
    fpdtype_t invrho = 1.0/rho, E = s[${nvars - 1}];

    fpdtype_t rhoe = E - 0.5*invrho*${pyfr.dot('s[{i}]', i=ndims)};

    // Linear initial guess for temperature based T = e / cv(298.15K)
    T = rhoe / (${' + '.join(f"{s[f'Cv{i}_ref']*s[f'M{i}']}*s[{ndims + 1 + i}]") for i in range(nspec)});

    ${pyfr.expand('temperature', 'rho', 'rhoe', 'T', 'p', t_kmax=3)};

    // Compute the velocities
% for i in range(ndims):
    v[${i}] = invrho*s[${i}];
% endfor

    // Momentum fluxes
    f[0] = s[0]*v[0] + p;
% for i in range(1, ndims):
    f[${i}] = s[0]*v[${i}];
% endfor

    // Energy flux
    f[${ndims}] = v[0]*(E + p);

    // Concentration flux
% for i in range(ndims + 1, ndims + 1 + nspec):
    f[${i}] = v[0]*s[${i}];
% endfor

</%pyfr:macro>
