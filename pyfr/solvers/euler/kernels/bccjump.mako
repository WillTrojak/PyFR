<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[3]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t pl = ${c['gamma'] - 1}*(ul[${nvars - 1}] - 0.5/ul[0]*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))});

    // Write out the jumps
    jumpl[0] = 0;
    jumpl[1] = fabs(pl);
    jumpl[2] = mag_nl;
</%pyfr:kernel>
