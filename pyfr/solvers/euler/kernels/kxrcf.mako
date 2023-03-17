<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='kxrcf' ndim='1'
              sensor='inout fpdtype_t[1][3]'>
    fpdtype_t sarea = sensor[0][2];
    fpdtype_t s_norm = sqrt(sensor[0][1]);
    fpdtype_t h = pow(sensor[0][2], ${1/(ndims - 1)})*${1/pi};

    sensor[0][0] = abs(sensor[0][0])/(pow(h, ${0.5*(order + 1)})*sarea*s_norm);
</%pyfr:kernel>
