<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mpnavstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
% for i in range(nspec):
    ur[${i}] = ul[${i}];
% endfor
% for i in range(nspec, ndims + nspec):
    ur[${i}] = -ul[${i}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}]
                     - (0.5/ul[0])*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))};
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_grad_state' params='ul, nl, grad_ul, grad_ur'>
    fpdtype_t rho = ${' + '.join('ul[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rcprho = 1.0/rho;

% if ndims == 2:
<% uidx = nspec + 0 %>
<% vidx = nspec + 1 %>
<% Eidx = nspec + ndims %>
    fpdtype_t rhou = ul[${uidx}], rhov = ul[${vidx}], E = ul[${Eidx}];

    fpdtype_t u = rcprho*rhou, v = rcprho*rhov;

    fpdtype_t rho_x = ${' + '.join('grad_ul[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_ul[1][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v])
    fpdtype_t u_x = grad_ul[0][${uidx}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${uidx}] - u*rho_y;
    fpdtype_t v_x = grad_ul[0][${vidx}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vidx}] - v*rho_y;

    fpdtype_t E_x = grad_ul[0][${Eidx}];
    fpdtype_t E_y = grad_ul[1][${Eidx}];

    // Compute temperature derivatives (c_v*dT/d[x,y])
    fpdtype_t Tl_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x));
    fpdtype_t Tl_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y));

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Correct copied across in-fluid temp gradients to in-wall gradients
    grad_ur[0][3] -= nl[0]*nl[0]*Tl_x + nl[0]*nl[1]*Tl_y;
    grad_ur[1][3] -= nl[1]*nl[0]*Tl_x + nl[1]*nl[1]*Tl_y;

% elif ndims == 3:
<% uidx = nspec + 0 %>
<% vidx = nspec + 1 %>
<% widx = nspec + 2 %>
<% Eidx = nspec + ndims %>
    fpdtype_t rhou = ul[${uidx}], rhov = ul[${vidx}], rhow = ul[${widx}];
    fpdtype_t E    = ul[${Eidx}];

    fpdtype_t u = rcprho*rhou, v = rcprho*rhov, w = rcprho*rhow;

    fpdtype_t rho_x = ${' + '.join('grad_ul[0][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_y = ${' + '.join('grad_ul[1][{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rho_z = ${' + '.join('grad_ul[2][{i}]'.format(i=i) for i in range(nspec))};

    // Velocity derivatives (rho*grad[u,v,w])
    fpdtype_t u_x = grad_ul[0][${uidx}] - u*rho_x;
    fpdtype_t u_y = grad_ul[1][${uidx}] - u*rho_y;
    fpdtype_t u_z = grad_ul[2][${uidx}] - u*rho_z;
    fpdtype_t v_x = grad_ul[0][${vidx}] - v*rho_x;
    fpdtype_t v_y = grad_ul[1][${vidx}] - v*rho_y;
    fpdtype_t v_z = grad_ul[2][${vidx}] - v*rho_z;
    fpdtype_t w_x = grad_ul[0][${widx}] - w*rho_x;
    fpdtype_t w_y = grad_ul[1][${widx}] - w*rho_y;
    fpdtype_t w_z = grad_ul[2][${widx}] - w*rho_z;

    fpdtype_t E_x = grad_ul[0][${Eidx}];
    fpdtype_t E_y = grad_ul[1][${Eidx}];
    fpdtype_t E_z = grad_ul[2][${Eidx}];

    // Compute temperature derivatives (c_v*dT/d[x,y,z])
    fpdtype_t Tl_x = rcprho*(E_x - (rcprho*rho_x*E + u*u_x + v*v_x + w*w_x));
    fpdtype_t Tl_y = rcprho*(E_y - (rcprho*rho_y*E + u*u_y + v*v_y + w*w_y));
    fpdtype_t Tl_z = rcprho*(E_z - (rcprho*rho_z*E + u*u_z + v*v_z + w*w_z));

    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ul', 'nl', 'grad_ul', 'grad_ur')};

    // Correct copied across in-fluid temp gradients to in-wall gradients
    grad_ur[0][4] -= nl[0]*nl[0]*Tl_x + nl[0]*nl[1]*Tl_y + nl[0]*nl[2]*Tl_z;
    grad_ur[1][4] -= nl[1]*nl[0]*Tl_x + nl[1]*nl[1]*Tl_y + nl[1]*nl[2]*Tl_z;
    grad_ur[2][4] -= nl[2]*nl[0]*Tl_x + nl[2]*nl[1]*Tl_y + nl[2]*nl[2]*Tl_z;
% endif
</%pyfr:macro>
