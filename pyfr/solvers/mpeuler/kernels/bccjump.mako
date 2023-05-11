<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.mpeuler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[3]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});

    // Calculate pressure
    fpdtype_t invdl = 1/${' + '.join('ul[{i}]'.format(i=i) for i in range(nspec))};
    fpdtype_t rhoel = ul[${ndims + nspec}] - 0.5*invdl*${pyfr.dot('ul[{i}]', i=(nspec, nspec + ndims))};
    fpdtype_t gl = (${' + '.join('{cp}*ul[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
                   (${' + '.join('{cv}*ul[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    fpdtype_t pl = rhoel*(gl - 1);

    // Write out the jumps
    jumpl[0] = 0;
    jumpl[1] = abs(pl);
    jumpl[2] = mag_nl;
</%pyfr:kernel>
