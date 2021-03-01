# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% t_tol = 0.99 %>

// Transforms to m=[1,0,0]^T
// See Moler and Hughes 1999
<%pyfr:macro name='transform_to' params='n, u, t, o'>

% if ndims == 2:
    t[o + 0] =  n[0]*u[o + 0] + n[1]*u[o + 1];
    t[o + 1] = -n[1]*u[o + 0] + n[0]*u[o + 1];
% elif ndims == 3:
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);
        t[o + 0] =  n[0]*u[o + 0] +                 n[1]*u[o + 1] +                 n[2]*u[o + 2];
        t[o + 1] = -n[1]*u[o + 0] + (n[0] + h*n[2]*n[2])*u[o + 1] -          h*n[1]*n[2]*u[o + 2];
        t[o + 2] = -n[2]*u[o + 0] -          h*n[1]*n[2]*u[o + 1] + (n[0] + h*n[1]*n[1])*u[o + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);

        t[o + 0] =              n[0]*u[o + 0] + n[1]*u[o + 1] +              n[2]*u[o + 2];
        t[o + 1] = (1 - h*n[0]*n[0])*u[o + 0] + n[0]*u[o + 1] -       h*n[0]*n[2]*u[o + 2];
        t[o + 2] =      -h*n[0]*n[2]*u[o + 0] + n[2]*u[o + 1] + (1 - h*n[2]*n[2])*u[o + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);

        t[o + 0] =              n[0]*u[o + 0] +              n[1]*u[o + 1] + n[2]*u[o + 2];
        t[o + 1] =      -h*n[0]*n[1]*u[o + 0] + (1 - h*n[1]*n[1])*u[o + 1] + n[1]*u[o + 2];
        t[o + 2] = (1 - h*n[0]*n[0])*u[o + 0] -       h*n[0]*n[1]*u[o + 1] + n[0]*u[o + 2];
    }
% endif
</%pyfr:macro>

// Transforms from m=[1,0,0]^T
<%pyfr:macro name='transform_from' params='n, t, u, o'>

% if ndims == 2:
    u[o + 0] = n[0]*t[o + 0] - n[1]*t[o + 1];
    u[o + 1] = n[1]*t[o + 0] + n[0]*t[o + 1];
% elif ndims == 3:
    if (fabs(n[0]) < ${t_tol})
    {
        fpdtype_t h = 1/(1 + n[0]);

        u[o + 0] = n[0]*t[o + 0] -                 n[1]*t[o + 1] -                 n[2]*t[o + 2];
        u[o + 1] = n[1]*t[o + 0] + (n[0] + h*n[2]*n[2])*t[o + 1] -          h*n[1]*n[2]*t[o + 2];
        u[o + 2] = n[2]*t[o + 0] -          h*n[1]*n[2]*t[o + 1] + (n[0] + h*n[1]*n[1])*t[o + 2];
    }
    else if (fabs(n[1]) < fabs(n[2]))
    {
        fpdtype_t h = 1/(1 - n[1]);

        u[o + 0] = n[0]*t[o + 0] + (1 - h*n[0]*n[0])*t[o + 1] -       h*n[0]*n[2]*t[o + 2];
        u[o + 1] = n[1]*t[o + 0] +              n[0]*t[o + 1] +              n[2]*t[o + 2];
        u[o + 2] = n[2]*t[o + 0] -       h*n[0]*n[2]*t[o + 1] + (1 - h*n[2]*n[2])*t[o + 2];
    }
    else
    {
        fpdtype_t h = 1/(1 - n[2]);

        u[o + 0] = n[0]*t[o + 0] -       h*n[0]*n[1]*t[o + 1] + (1 - h*n[0]*n[0])*t[o + 2];
        u[o + 1] = n[1]*t[o + 0] + (1 - h*n[1]*n[1])*t[o + 1] -       h*n[0]*n[1]*t[o + 2];
        u[o + 2] = n[2]*t[o + 0] +              n[1]*t[o + 1] +              n[0]*t[o + 2];
    }
% endif
</%pyfr:macro>
