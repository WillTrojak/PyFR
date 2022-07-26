# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t V_e = ${' + '.join('{0}*nl[{1}]'.format(c['uvw'[i]], i)
                                 for i in range(ndims))};
    fpdtype_t V_i = ${' + '.join('ul[{0}]*nl[{1}]'.format(i + 1, i)
                                 for i in range(ndims))};

    fpdtype_t C_i = sqrt(V_i*V_i + ul[0] + ${c['bc-ed-zeta']});
    fpdtype_t C_e = sqrt(V_e*V_e + ${c['p']} + ${c['bc-ed-zeta']});

    fpdtype_t mVC_i = 0.5*V_i + C_i;
    fpdtype_t mVC_e = 0.5*V_e - C_e;
    fpdtype_t dp = ul[0] - ${c['p']};

    fpdtype_t V_b = min(V_i, V_e);
    fpdtype_t C_bl = sqrt(0.25*V_b*V_b + ${c['bc-ed-zeta']} + min(ul[0], ${c['p']}));
    fpdtype_t C_br = C_bl;

    fpdtype_t dpl, dpr, gdpl, gdpr;
% for i in range(c['niters']):    
    if (1.5*V_i - C_i < 1.5*V_b - C_bl)
    {   // Rarefaction
        dpl = 0.25*(V_i*V_i - V_b*V_b);
        gdpl = -0.5*V_b;
    }
    else
    {   // Shock
        dpl  = -(V_b - V_i)*mVC_i;
        gdpl = -mVC_i;
    }

    if ( (1.5*V_b + C_br < 1.5*V_e + C_e) & (1.5*V_i - C_i >= 1.5*V_b - C_bl) )
    {   // Rarefaction
        dpr = -0.25*(V_e*V_e - V_b*V_b);
        gdpr = 0.5*V_b;
    }
    else
    {   // Shock
        dpr = -(V_e - V_b)*mVC_e;
        gdpr = mVC_e;
    }

    C_bl = sqrt(0.25*V_b*V_b + ${c['bc-ed-zeta']} + dpl + ul[0]);
    C_br = sqrt(0.25*V_b*V_b + ${c['bc-ed-zeta']} + ${c['p']} - dpr);
    V_b -= (dp + dpr + dpl)/(gdpl + gdpr);
% endfor

% for i in range(ndims):
    ur[${i + 1}] = (V_i >= 0)
                 ? (ul[${i + 1}]+ (V_b - V_i)*nl[${i}])
                 : (${c['uvw'[i]]} + (V_b - V_e)*nl[${i}]);
% endfor
    ur[0] = 0.5*(dpl + ul[0] + ${c['p']} - dpr);
</%pyfr:macro>
