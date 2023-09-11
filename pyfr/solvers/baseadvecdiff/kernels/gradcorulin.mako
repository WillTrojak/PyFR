<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.smats'/>
<%include file='pyfr.solvers.baseadvecdiff.kernels.transform_grad'/>

<%pyfr:kernel name='gradcorulin' ndim='2'
              gradu='inout fpdtype_t[${str(ndims)}][${str(nvars)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'>
    // Compute the S matrices
    fpdtype_t smats[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', 'smats', 'djac')};

    fpdtype_t rcpdjac = 1 / djac;

    ${pyfr.expand('transform_grad', 'gradu', 'smats', 'rcpdjac')};
</%pyfr:kernel>
