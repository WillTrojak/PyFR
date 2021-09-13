# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.edeuler.kernels.flux'/>

<% t_tol = 0.99 %>
<%pyfr:macro name='transform_to' params='n, u, t, offset'>
% if ndims == 2:
    t[offset + 0] =  n[0]*u[offset + 0] + n[1]*u[offset + 1];
    t[offset + 1] = -n[1]*u[offset + 0] + n[0]*u[offset + 1];
% elif ndims == 3:
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);
        t[offset + 0] =  n[0]*u[offset + 0] + n[1]*u[offset + 1] + n[2]*u[offset + 2];
        t[offset + 1] = -n[1]*u[offset + 0] + (n[0] + h*n[2]*n[2])*u[offset + 1] - h*n[1]*n[2]*u[offset + 2];
        t[offset + 2] = -n[2]*u[offset + 0] - h*n[1]*n[2]*u[offset + 1] + (n[0] + h*n[1]*n[1])*u[offset + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);
        t[offset + 0] = n[0]*u[offset + 0] + n[1]*u[offset + 1] + n[2]*u[offset + 2];
        t[offset + 1] =  (1 - h*n[0]*n[0])*u[offset + 0] + n[0]*u[offset + 1] - h*n[0]*n[2]*u[offset + 2];
        t[offset + 2] = -h*n[0]*n[2]*u[offset + 0] + n[2]*u[offset + 1] + (1 - h*n[2]*n[2])*u[offset + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);
        t[offset + 0] = n[0]*u[offset + 0] + n[1]*u[offset + 1] + n[2]*u[offset + 2];
        t[offset + 1] = -h*n[0]*n[1]*u[offset + 0] + (1 - h*n[1]*n[1])*u[offset + 1] + n[1]*u[offset + 2];
        t[offset + 2] =  (1 - h*n[0]*n[0])*u[offset + 0] - h*n[0]*n[1]*u[offset + 1] + n[0]*u[offset + 2];
    }
% endif
</%pyfr:macro>

// Transforms from m=[1,0,0]^T
<%pyfr:macro name='transform_from' params='n, t, u, offset'>
% if ndims == 2:
    u[offset + 0] = n[0]*t[offset + 0] - n[1]*t[offset + 1];
    u[offset + 1] = n[1]*t[offset + 0] + n[0]*t[offset + 1];
% elif ndims == 3:
    u[0] = t[0];
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);
        u[offset + 0] =  n[0]*t[offset + 0] - n[1]*t[offset + 1] - n[2]*t[offset + 2];
        u[offset + 1] =  n[1]*t[offset + 0] + (n[0] + h*n[2]*n[2])*t[offset + 1] - h*n[1]*n[2]*t[offset + 2];
        u[offset + 2] =  n[2]*t[offset + 0] - h*n[1]*n[2]*t[offset + 1] + (n[0] + h*n[1]*n[1])*t[offset + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);
        u[offset + 0] = n[0]*t[offset + 0] +  (1 - h*n[0]*n[0])*t[offset + 1] - h*n[0]*n[2]*t[offset + 2];
        u[offset + 1] = n[1]*t[offset + 0] + n[0]*t[offset + 1] + n[2]*t[offset + 2];
        u[offset + 2] = n[2]*t[offset + 0] - h*n[0]*n[2]*t[offset + 1] + (1 - h*n[2]*n[2])*t[offset + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);
        u[offset + 0] = n[0]*t[offset + 0] - h*n[0]*n[1]*t[offset + 1] + (1 - h*n[0]*n[0])*t[offset + 2];
        u[offset + 1] = n[1]*t[offset + 0] + (1 - h*n[1]*n[1])*t[offset + 1] - h*n[0]*n[1]*t[offset + 2];
        u[offset + 2] = n[2]*t[offset + 0] + n[1]*t[offset + 1] + n[0]*t[offset + 2];
    }
% endif
</%pyfr:macro>

<%pyfr:macro name='rsolve_t1d' params='ul, ur, nf'>
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t nf_flstar, nf_frstar;

    fpdtype_t pl = ul[0];
    fpdtype_t pr = ur[0];

    fpdtype_t um = 0.5*(ul[1] + ur[1]);
    fpdtype_t pm = 0.5*(pl + pr);

    // Estimate the wave speeds
    fpdtype_t sl = 1.5*um - sqrt(0.25*um*um + ${c['ed-zeta']} + pm);
    fpdtype_t sr = 1.5*um + sqrt(0.25*um*um + ${c['ed-zeta']} + pm);
    fpdtype_t ss = (sl*sr*(pl - pr) - sr*ul[1]*(pl + ${c['ed-zeta']})
                                    + sl*ur[1]*(pr + ${c['ed-zeta']})) /
                   ((sl - ul[1])*(pl + ${c['ed-zeta']}) -
                    (sr - ur[1])*(pr + ${c['ed-zeta']}));

    // Precompute some inverses
    fpdtype_t slmss_inv = 1 / (sl - ss);
    fpdtype_t srmss_inv = 1 / (sr - ss);

    // Estimate star pressure
    fpdtype_t ps = ((sl - ul[1])*pl + ${c['ed-zeta']}*(ss - ul[1])) * slmss_inv;
    
    // Set left and right flux
    fl[0] = ul[1]*(pl + ${c['ed-zeta']});
    fr[0] = ur[1]*(pr + ${c['ed-zeta']});
    fl[1] = ul[1]*ul[1] + pl;
    fr[1] = ur[1]*ur[1] + pr;
% for i in range(2, nvars):
    fl[${i}] = ul[${i}]*ul[1];
    fr[${i}] = ur[${i}]*ur[1];
% endfor

    // Output
% for i in range(nvars):
% if i == 0:
    nf_flstar = ss*(ps + ${c['ed-zeta']});
    nf_frstar = nf_flstar;
% elif i == 1:
    nf_flstar = ss*ss + ps;
    nf_frstar = nf_flstar;
% else:
    nf_flstar = ss*ul[${i}]*(sl - ul[1]) * slmss_inv;
    nf_frstar = ss*ur[${i}]*(sr - ur[1]) * srmss_inv;
% endif

    nf[${i}] = (sl >= 0) ? fl[${i}] : (sl <= 0 && ss >= 0) ? nf_flstar :
               (ss <= 0 && sr >= 0) ? nf_frstar : fr[${i}];
% endfor
</%pyfr:macro>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t utl[${nvars}], utr[${nvars}], ntf[${nvars}];
    utl[0] = ul[0]; utr[0] = ur[0];
    ${pyfr.expand('transform_to','n', 'ul', 'utl', '1')};
    ${pyfr.expand('transform_to','n', 'ur', 'utr', '1')};
    ${pyfr.expand('rsolve_t1d', 'utl', 'utr', 'ntf')};
    nf[0] = ntf[0];
    ${pyfr.expand('transform_from', 'n', 'ntf', 'nf', '1')};
</%pyfr:macro>
