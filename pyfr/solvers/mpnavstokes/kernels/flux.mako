<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:
<% uidx = nspec + 0 %>
<% vidx = nspec + 1 %>
<% Eidx = nspec + ndims %>
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho = ${' + '.join('uin[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rhou = uin[${uidx}], rhov = uin[${vidx}], E = uin[${Eidx}];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;

    fpdtype_t rho_x = ${' + '.join('grad_uin[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_uin[1][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v])
    fpdtype_t u_x = grad_uin[0][${uidx}] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][${uidx}] - u*rho_y;
    fpdtype_t v_x = grad_uin[0][${vidx}] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][${vidx}] - v*rho_y;

    fpdtype_t E_x = grad_uin[0][${Eidx}];
    fpdtype_t E_y = grad_uin[1][${Eidx}];

    fpdtype_t mu_c = (${' + '.join('{mu}*uin[{i}]'.format(i=i, mu=c[f'mu{i}'])
                       for i in range(nspec))})*rcprho;
    fpdtype_t Pr_c = (${' + '.join('{Pr}*uin[{i}]'.format(i=i, Pr=c[f'Pr{i}'])
                        for i in range(nspec))})*rcprho;
    fpdtype_t g = (${' + '.join('{cp}*uin[{i}]'.format(i=i, cp=c[f'cp{i}'])
                     for i in range(nspec))}) /
                  (${' + '.join('{cv}*uin[{i}]'.format(i=i, cv=c[f'cv{i}'])
                     for i in range(nspec))});
    fpdtype_t kappa = -mu_c*g/Pr_c;

    // Compute temperature derivatives (c_v*dT/d[x,y])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);

    fout[0][${uidx}] += t_xx; fout[1][${uidx}] += t_xy;
    fout[0][${vidx}] += t_xy; fout[1][${vidx}] += t_yy;

    fout[0][${Eidx}] += u*t_xx + v*t_xy + kappa*T_x;
    fout[1][${Eidx}] += u*t_xy + v*t_yy + kappa*T_y;
</%pyfr:macro>
% elif ndims == 3:
<% uidx = nspec + 0 %>
<% vidx = nspec + 1 %>
<% widx = nspec + 2 %>
<% Eidx = nspec + ndims %>
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t rho  = ${' + '.join('uin[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rhou = uin[${uidx}], rhov = uin[${vidx}], rhow = uin[${widx}];
    fpdtype_t E    = uin[${Eidx}];

    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = ${' + '.join('grad_uin[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_uin[1][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_z = ${' + '.join('grad_uin[2][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_uin[0][${uidx}] - u*rho_x;
    fpdtype_t u_y = grad_uin[1][${uidx}] - u*rho_y;
    fpdtype_t u_z = grad_uin[2][${uidx}] - u*rho_z;
    fpdtype_t v_x = grad_uin[0][${vidx}] - v*rho_x;
    fpdtype_t v_y = grad_uin[1][${vidx}] - v*rho_y;
    fpdtype_t v_z = grad_uin[2][${vidx}] - v*rho_z;
    fpdtype_t w_x = grad_uin[0][${widx}] - w*rho_x;
    fpdtype_t w_y = grad_uin[1][${widx}] - w*rho_y;
    fpdtype_t w_z = grad_uin[2][${widx}] - w*rho_z;

    fpdtype_t E_x = grad_uin[0][${Eidx}];
    fpdtype_t E_y = grad_uin[1][${Eidx}];
    fpdtype_t E_z = grad_uin[2][${Eidx}];

    fpdtype_t mu_c = (${' + '.join('{mu}*uin[{i}]'.format(i=i, mu=c[f'mu{i}'])
                       for i in range(nspec))})*rcprho;
    fpdtype_t Pr_c = (${' + '.join('{Pr}*uin[{i}]'.format(i=i, Pr=c[f'Pr{i}'])
                        for i in range(nspec))})*rcprho;
    fpdtype_t g = (${' + '.join('{cp}*uin[{i}]'.format(i=i, cp=c[f'cp{i}'])
                     for i in range(nspec))}) /
                  (${' + '.join('{cv}*uin[{i}]'.format(i=i, cv=c[f'cv{i}'])
                     for i in range(nspec))});
    fpdtype_t kappa = -mu_c*g/Pr_c;

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t T_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x + w*w_x));
    fpdtype_t T_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y + w*w_y));
    fpdtype_t T_z = rcprho*(E_z - (rcprho*rho_z*E + u*u_z + v*v_z + w*w_z));

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*mu_c*rcprho*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*mu_c*rcprho*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*mu_c*rcprho*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -mu_c*rcprho*(v_x + u_y);
    fpdtype_t t_xz = -mu_c*rcprho*(u_z + w_x);
    fpdtype_t t_yz = -mu_c*rcprho*(w_y + v_z);

    fout[0][${uidx}] += t_xx; fout[1][${uidx}] += t_xy; fout[2][${uidx}] += t_xz;
    fout[0][${vidx}] += t_xy; fout[1][${vidx}] += t_yy; fout[2][${vidx}] += t_yz;
    fout[0][${widx}] += t_xz; fout[1][${widx}] += t_yz; fout[2][${widx}] += t_zz;

    fout[0][${Eidx}] += u*t_xx + v*t_xy + w*t_xz + kappa*T_x;
    fout[1][${Eidx}] += u*t_xy + v*t_yy + w*t_yz + kappa*T_y;
    fout[2][${Eidx}] += u*t_xz + v*t_yz + w*t_zz + kappa*T_z;
</%pyfr:macro>
% endif
