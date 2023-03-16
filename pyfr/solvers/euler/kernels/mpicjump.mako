<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='mpicjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              ur='in mpi fpdtype_t[${str(nvars)}]'
              jumpl='out fpdtype_t[1]'>
    fpdtype_t pl = ${c['gamma'] - 1}*(ul[${nvars - 1}] - 0.5/ul[0]*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))});
    fpdtype_t pr = ${c['gamma'] - 1}*(ur[${nvars - 1}] - 0.5/ur[0]*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))});

    // Write out the jumps
    jumpl[${i}] = pl - pr;
</%pyfr:kernel>