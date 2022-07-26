# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.edeuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t nf_fl, nf_fr, nf_sub;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

    // Normal of the average interface velocity
    fpdtype_t unl = ${pyfr.dot('n[{i}]', 'ul[1+{i}]', i=ndims)};
    fpdtype_t unr = ${pyfr.dot('n[{i}]', 'ur[1+{i}]', i=ndims)};

    // Estimate the wave speed
    fpdtype_t sl = min(1.5*unl - sqrt(0.25*unl*unl + ${c['ed-zeta']} + ul[0]),
                       1.5*unr - sqrt(0.25*unr*unr + ${c['ed-zeta']} + ur[0])); 
    fpdtype_t sr = max(1.5*unl + sqrt(0.25*unl*unl + ${c['ed-zeta']} + ul[0]),
                       1.5*unr + sqrt(0.25*unr*unr + ${c['ed-zeta']} + ur[0])); 
    fpdtype_t srmsl_inv = 1 / (sr - sl);
    
    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join('n[{j}]*fl[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
    nf_fr = ${' + '.join('n[{j}]*fr[{j}][{i}]'.format(i=i, j=j)
                         for j in range(ndims))};
    nf_sub = (sr*nf_fl - sl*nf_fr + sl*sr*(ur[${i}] - ul[${i}]))*srmsl_inv;
    nf[${i}] = (0 <= sl) ? nf_fl : (0 >= sr) ? nf_fr : nf_sub;
% endfor
</%pyfr:macro>