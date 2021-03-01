# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.hypenavstokes.kernels.flux'/>

<% nurtr = c['nu']/c['tr'] %>
<% rtr = 1/c['tr'] %>
<% znurtr = c['ac-zeta'] + c['nu']/c['tr'] %>

<%pyfr:macro name='rsolve_t1d' params='ql, qr, nf'>
% if ndims == 2:
    fpdtype_t fl[${nvars}], fr[${nvars}];
    fpdtype_t dq[${nvars}], a[4], l[4];

    ${pyfr.expand('inviscid_1dflux', 'ql', 'fl')};
    ${pyfr.expand('inviscid_1dflux', 'qr', 'fr')};

    // Preconditioned jump array
    dq[0] = qr[6] - ql[6];
    dq[1] = qr[4] - ql[4];
    dq[2] = qr[3] - ql[3];
    dq[3] = qr[5] - ql[5];
    dq[4] = qr[2] - ql[2];
    dq[5] = qr[1] - ql[1];
    dq[6] = qr[0] - ql[0];

    fpdtype_t uh = 0.5*(ql[1] + qr[1]);
    fpdtype_t vh = 0.5*(ql[2] + qr[2]);
    fpdtype_t ah = sqrt(uh*uh + ${znurtr});
    fpdtype_t bh = sqrt(0.25*uh*uh + ${nurtr});
    fpdtype_t ch = sqrt(uh*uh + ${c['ac-zeta']});

    fpdtype_t rcpau = 1/(ch*ch + ah*uh);
    fpdtype_t rcmau = 1/(ch*ch - ah*uh);

    l[0] = fabs(uh + bh);
    l[1] = fabs(uh - bh);
    l[2] = fabs(uh + ah);
    l[3] = fabs(uh - ah);
      
    a[3] = -((dq[6] - ${c['nu']}*dq[2])*(uh + ah)*${0.5/znurtr} - 0.5*dq[5])/ah;
    a[2] = -(dq[5] + (uh - ah)*a[3])*rau; 
    a[1] =  (dq[4] + ${c['tr']}*(0.5*uh + bh)*dq[3] + 
                  + a[2]*vh*(uh + ah)*(0.5*uh - bh + ah)*rcpau 
                  + a[3]*vh*(uh - ah)*(0.5*uh - bh - ah)*rcmau)*${0.5*rtr}/bh;
    a[0] = (dq[3] - a[1] - a[2]*vh*(uh + ah)*${rtr}*rcpau 
                         - a[3]*vh*(uh - ah)*${rtr}*rcmau);
    
    nf[0] = 0.5*(fr[0] + fl[0]) + ${0.5*c['ac-zeta']}*(l[2]*a[2] + l[3]*a[3]);
    nf[1] = 0.5*(fr[1] + fl[1] + l[3]*a[3]*(uh - ah) + l[2]*a[2]*(uh + ah));
    nf[2] = 0.5*(fr[2] + fl[2] + 
                 ${c['tr']}*(l[0]*a[0]*(0.5*uh + bh) + l[1]*a[1]*(0.5*uh - bh)) + 
                 vh*(l[2]*a[2]*rcpau*(uh + ah)*(uh + ah) + l[3]*a[3]*rcmau*(uh - ah)*(uh - ah)));
    nf[3] = 0.5*(fr[3] + fl[3]) - ${0.5*rtr}*(l[2]*a[2] + l[3]*a[3]);
    nf[4] = 0.5*(fr[4] + fl[4]);
    nf[5] = 0.5*(fr[5] + fl[5] - l[0]*a[0] - l[1]*a[1] -
                 vh*${rtr}*((uh + ah)*rcpau*l[5]*a[2] + (uh - ah)*rcmau*l[3]*a[3]));
    nf[6] = 0.5*(fr[6] + fl[6]);
% elif ndims == 3:

%endif
</%pyfr:macro>

<%include file='pyfr.solvers.hypenavstokes.kernels.rsolvers.rsolve_trans'/>
