<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve_1d' params='ul, ur, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t usl[${nvars}], usr[${nvars}], fsl, fsr;
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr, Tl, Tr;

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'pl', 'Tl', 'rhol', 'vl')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'pr', 'Tr', 'rhor', 'vr')};

    // Calculate cpl, cv, and gamma
    fpdtype_t cpl_s = ${pyfr.polyval(s['Cp0'], 'Tl')};
    fpdtype_t cpr_s = ${pyfr.polyval(s['Cp0'], 'Tr')};
    fpdtype_t cpl = ul[${ndims + 1}]*${s['R0']*s['M0']}*cpl_s;
    fpdtype_t cpr = ur[${ndims + 1}]*${s['R0']*s['M0']}*cpr_s;
    fpdtype_t cvl = ul[${ndims + 1}]*${s['R0']*s['M0']}*(cpl_s - 1);
    fpdtype_t cvr = ur[${ndims + 1}]*${s['R0']*s['M0']}*(cpr_s - 1);
% for i in range(1, nspec):
    cpl_s = ${pyfr.polyval(s[f'Cp{i}'], 'Tl')};
    cpr_s = ${pyfr.polyval(s[f'Cp{i}'], 'Tr')};
    cpl += ${s[f'R{i}']*s[f'M{i}']}*ul[${ndims + 1 + i}]*cpl_s;
    cpr += ${s[f'R{i}']*s[f'M{i}']}*ur[${ndims + 1 + i}]*cpr_s;
    cvl += ${s[f'R{i}']*s[f'M{i}']}*ul[${ndims + 1 + i}]*(cpl_s - 1);
    cvr += ${s[f'R{i}']*s[f'M{i}']}*ur[${ndims + 1 + i}]*(cpr_s - 1);
%endfor
    fpdtype_t gl = cpl / cvl;
    fpdtype_t gr = cpr / cvr;

    // Wave speeds
    fpdtype_t cl = sqrt(gl*pl/rhol);
    fpdtype_t cr = sqrt(gr*pr/rhor);
    fpdtype_t sl = min(vl[0] - cl, vr[0] - cr);
    fpdtype_t sr = min(vl[0] + cl, vr[0] + cr);
    fpdtype_t sstar = (pr - pl + ul[0]*(sl - vl[0]) - ur[0]*(sr - vr[0])) /
                      (rhol*(sl - vl[0]) - rhor*(sr - vr[0]));

    fpdtype_t ul_com = (sl - vl[0]) / (sl - sstar);
    fpdtype_t ur_com = (sr - vr[0]) / (sr - sstar);

    // Momentum
    usl[0] = ul_com*rhol*sstar;
    usr[0] = ur_com*rhor*sstar;
% for i in range(1, ndims):
    usl[${i}] = ul_com*ul[${i}];
    usr[${i}] = ur_com*ur[${i}];
% endfor

    // Energy
    usl[${ndims}] = ul_com*(ul[${ndims}] + rhol*(sstar - vl[0])*(sstar + pl/(rhol*(sl - vl[0]))));
    usr[${ndims}] = ur_com*(ur[${ndims}] + rhor*(sstar - vr[0])*(sstar + pr/(rhor*(sr - vr[0]))));

    // Concentration
% for i in range(ndims + 1, ndims + 1 + nspec):
    usl[${i}] = ul_com*ul[${i}];
    usr[${i}] = ur_com*ur[${i}];
% endfor

    // Output
% for i in range(nvars):
    fsl = fl[${i}] + sl*(usl[${i}] - ul[${i}]);
    fsr = fr[${i}] + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? fl[${i}] : (sl <= 0 && 0 <= sstar) ? fsl : (sstar <= 0 && 0 <= sr) ? fsr : fr[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.rsolve1d'/>
