<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='kxrcf' ndim='1'
              jump='in fpdtype_t[${str(nfpts)}][3]'
              sensor='out fpdtype_t[1][3]'
              mass='in broadcast fpdtype_t[1][${str(nfpts)}]'>

    // Integrate jump around element
    sensor[0][0] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][0]', k=nfpts)};
    sensor[0][1] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][1]', k=nfpts)};
    sensor[0][2] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][2]', k=nfpts)};

    fpdtype_t sarea = sensor[0][2];
    fpdtype_t s_norm = sqrt(sensor[0][1]);
% if ndims == 2:
    fpdtype_t h = sensor[0][2]*${0.5/pi};
% elif ndims == 3:
    fpdtype_t h = sqrt(sensor[0][2]*${0.25/pi});
% endif

    sensor[0][0] = abs(sensor[0][0])/(pow(h, ${0.5*(order + 1)})*sarea*s_norm);
</%pyfr:kernel>
