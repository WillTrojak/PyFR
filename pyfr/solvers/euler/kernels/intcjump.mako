<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              ur='in view fpdtype_t[${str(nvars)}]'
              jumpl='out view fpdtype_t[1]'
              jumpr='out view fpdtype_t[1]'>
    fpdtype_t pl = ${c['gamma'] - 1}*(ul[${nvars - 1}] - 0.5/ul[0]*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))});
    fpdtype_t pr = ${c['gamma'] - 1}*(ur[${nvars - 1}] - 0.5/ur[0]*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))});

    // Write out the jumps
    jumpl[${i}] = pl - pr;
    jumpr[${i}] = pr - pl;
</%pyfr:kernel>
