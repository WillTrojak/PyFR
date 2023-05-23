<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='kxrcf' ndim='1'
              jump='in fpdtype_t[${str(nfpts)}][3]'
              sensor='out fpdtype_t[1][3]'
              mass='in broadcast fpdtype_t[1][${str(nfpts)}]'>

    // Integrate jump around element
    sensor[0][0] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][0]', k=nfpts)};
    sensor[0][1] = jump[0][1];
% for i in range(1, nfpts):
    sensor[0][1] = fmax(sensor[0][1], jump[${i}][1]);
% endfor
    sensor[0][2] = ${pyfr.dot('mass[0][{k}]', 'jump[{k}][2]', k=nfpts)};

    fpdtype_t s_norm = sensor[0][1];
    fpdtype_t sarea = sensor[0][2];
% if ndims == 2:
    fpdtype_t h = sarea*${1/pi};
% elif ndims == 3:
    fpdtype_t h = sqrt(sarea*${1/pi});
% endif

    sensor[0][0] = fabs(sensor[0][0])/(pow(h, ${order + 1.0})*sarea*s_norm);
</%pyfr:kernel>
