<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='mpicjump' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              ur='in mpi fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              jumpl='out view fpdtype_t[3]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});

    // Calculate pressure
    fpdtype_t invdl = 1/(${' + '.join('ul[{i}]'.format(i=i) for i in range(nspec))});
    fpdtype_t rhoel = ul[${ndims + nspec}] - 0.5*invdl*${pyfr.dot('ul[{i}]', i=(nspec, nspec + ndims))};
    fpdtype_t gl = (${' + '.join('{cp}*ul[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
                   (${' + '.join('{cv}*ul[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    fpdtype_t pl = rhoel*(gl - 1);

    fpdtype_t invdr = 1/(${' + '.join('ur[{i}]'.format(i=i) for i in range(nspec))});
    fpdtype_t rhoer = ul[${ndims + nspec}] - 0.5*invdr*${pyfr.dot('ur[{i}]', i=(nspec, nspec + ndims))};
    fpdtype_t gr = (${' + '.join('{cp}*ur[{i}]'.format(i=i, cp=c[f'cp{i}']) for i in range(nspec))}) /
                   (${' + '.join('{cv}*ur[{i}]'.format(i=i, cv=c[f'cv{i}']) for i in range(nspec))});
    fpdtype_t pr = rhoer*(gr - 1);

    // Write out the jumps
    jumpl[0] = mag_nl*(pl - pr);
    jumpl[1] = fabs(pl);
    jumpl[2] = mag_nl;
</%pyfr:kernel>
