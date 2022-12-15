<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpeuler.kernels.flux'/>

<%pyfr:macro name='rsolve_1d' params='ul, ur, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t usl[${nvars}], usr[${nvars}], fsl, fsr;
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t rhol, rhor, pl, pr, gl, gr;

    ${pyfr.expand('inviscid_flux_1d', 'ul', 'fl', 'rhol', 'pl', 'vl', 'gl')};
    ${pyfr.expand('inviscid_flux_1d', 'ur', 'fr', 'rhor', 'pr', 'vr', 'gr')};

    // Wave speeds
    fpdtype_t cl = sqrt(gl*pl/rhol);
    fpdtype_t cr = sqrt(gr*pr/rhor);
    fpdtype_t sl = min(vl[0] - cl, vr[0] - cr);
    fpdtype_t sr = min(vl[0] + cl, vr[0] + cr);
    fpdtype_t sstar = (pr - pl + ul[${nspec}]*(sl - vl[0]) - ur[${nspec}]*(sr - vr[0])) / 
                      (rhol*(sl - vl[0]) - rhor*(sr - vr[0]));
    
    fpdtype_t ul_com = (sl - vl[0]) / (sl - sstar);
    fpdtype_t ur_com = (sr - vr[0]) / (sr - sstar);

% for i in range(nspec):
    usl[${i}] = ul_com*ul[${i}];
    usr[${i}] = ur_com*ur[${i}];
% endfor

    usl[${nspec}] = ul_com*rhol*sstar;
    usr[${nspec}] = ur_com*rhor*sstar;
% for i in range(1, ndims):
    usl[${nspec + i}] = ul_com*ul[${nspec + i}];
    usr[${nspec + i}] = ur_com*ur[${nspec + i}];
% endfor

    usl[${nspec + ndims}] = ul_com*(ul[${nspec + ndims}] + rhol*(sstar - vl[0])*(sstar + pl/(rhol*(sl - vl[0]))));
    usr[${nspec + ndims}] = ur_com*(ur[${nspec + ndims}] + rhor*(sstar - vr[0])*(sstar + pr/(rhor*(sr - vr[0]))));

    // Output
% for i in range(nvars):
    fsl = fl[${i}] + sl*(usl[${i}] - ul[${i}]);
    fsr = fr[${i}] + sr*(usr[${i}] - ur[${i}]);
    nf[${i}] = (0 <= sl) ? fl[${i}] : (sl <= 0 && 0 <= sstar) ? fsl : (sstar <= 0 && 0 <= sr) ? fsr : fr[${i}];
% endfor
</%pyfr:macro>

<%include file='pyfr.solvers.mpeuler.kernels.rsolvers.rsolve1d'/>
