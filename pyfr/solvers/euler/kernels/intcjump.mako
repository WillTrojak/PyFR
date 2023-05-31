<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              ur='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[3]'
              jumpr='out view fpdtype_t[3]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t pl = ${c['gamma'] - 1}*(ul[${nvars - 1}] - 0.5/ul[0]*${pyfr.dot('ul[{i}]', i=(1, ndims + 1))});
    fpdtype_t pr = ${c['gamma'] - 1}*(ur[${nvars - 1}] - 0.5/ur[0]*${pyfr.dot('ur[{i}]', i=(1, ndims + 1))});

    // Write out the jumps
    jumpl[0] = mag_nl*(pl - pr);
    jumpl[1] = fabs(pl);
    jumpl[2] = mag_nl;

    jumpr[0] = mag_nl*(pr - pl);
    jumpr[1] = fabs(pr);
    jumpr[2] = mag_nl;
</%pyfr:kernel>
