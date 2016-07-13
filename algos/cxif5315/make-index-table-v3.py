#!/usr/bin/env python
#------------------------------
"""
Created on Nov 3, 2015

@author: Mikhail Dubrovin
"""
#------------------------------

#from algos.diffraction.FiberIndexing import *
from algos.diffraction.FiberIndexing import *

#------------------------------

def make_index_table(prefix='./v03-') :

    from algos.core.GlobalUtils import str_tstamp
    fname = '%slut-cxif5315-r0169-%s.txt' % (prefix, str_tstamp())
    fout = open(fname,'w')
    fout.write('# file name: %s\n' % fname)

    #------------------------------
    # Photon energy
    Egamma_eV  = 6003.1936                               # eV SIOC:SYS0:ML00:AO541
    wavelen_nm = wavelength_nm_from_energy_ev(Egamma_eV) # nm
    evald_rad  = wave_vector_value(Egamma_eV)            # 1/A
    #-------
    sigma_ql   = 0.003 * evald_rad
    sigma_qt   = 0.001 * evald_rad
    #-------
    rec  = '\n# photon energy = %.4f eV' % (Egamma_eV)\
         + '\n# wavelength = %.4f A' % (wavelen_nm*10)\
         + '\n# wave number/Evald radius k = 1/lambda = %.6f 1/A' % (evald_rad)\
         + '\n# sigma_ql = %.6f 1/A (approximately = k * <pixel size>/' % (sigma_ql)\
         + '\n# sigma_qt = %.6f 1/A (approximately = k * <pixel size>/' % (sigma_qt)\
         + '<sample-to-detector distance> = k*100um/100mm)'\
         + '\n# 3*sigma_ql = %.6f 1/A\n' % (3*sigma_ql)\
         + '\n# 3*sigma_qt = %.6f 1/A\n' % (3*sigma_qt)
    print rec
    fout.write(rec)

    #------------------------------
    # Lattice parameters
    # from previous analysis note:
    #a, b, c = 18.36, 26.65, 4.81        # Angstrom
    #alpha, beta, gamma = 90, 90, 77.17  # 180 - 102.83 degree
    a= 18.55 # Angstrom
    b, c = 1.466*a, 0.262*a              # Angstrom
    alpha, beta, gamma = 90, 90, 78.47   # 180 - 101.53 degree
    hmax, kmax, lmax = 4, 6, 0           # size of lattice to consider

    a1, a2, a3 = triclinic_primitive_vectors(a, b, c, alpha, beta, gamma)
    b1, b2, b3 = reciprocal_from_bravias(a1, a2, a3)

    msg1 = '\n# Triclinic crystal cell parameters:'\
         + '\n#   a = %.2f A\n#   b = %.2f A\n#   c = %.2f A' % (a, b, c)\
         + '\n#   alpha = %.2f deg\n#   beta  = %.2f deg\n#   gamma = %.2f deg' % (alpha, beta, gamma)

    fmt = '%10.6f'
    msg2 = '\n# 3-d space primitive vectors:\n#   a1 = (%s)\n#   a2 = (%s)\n#   a3 = (%s)' %\
           (', '.join([fmt % v for v in a1]),\
            ', '.join([fmt % v for v in a2]),\
            ', '.join([fmt % v for v in a3]))
           #(str(a1), str(a2), str(a3))

    msg3 = '\n# reciprocal space primitive vectors:\n#   b1 = (%s)\n#   b2 = (%s)\n#   b3 = (%s)' %\
           (', '.join([fmt % v for v in b1]),\
            ', '.join([fmt % v for v in b2]),\
            ', '.join([fmt % v for v in b3]))
           #(str(b1), str(b2), str(b3))

    rec = '%s\n%s\n%s\n' % (msg1, msg2, msg3)
    print rec
    fout.write(rec)

    fout.write('\n# %s\n\n' % (89*'_'))

    #for line in triclinic_primitive_vectors.__doc__.split('\n') : fout.write('\n# %s' % line)

    # 2-d lattice test
    test_lattice       (b1, b2, b3, hmax, kmax, lmax, cdtype=np.float32)
    plot_lattice       (b1, b2, b3, hmax, kmax, lmax, cdtype=np.float32, evald_rad=evald_rad, qtol=3*sigma_ql, prefix=prefix) #, do_movie=True, delay=1000)
    lattice_node_radius(b1, b2, b3, hmax, kmax, lmax, cdtype=np.float32)
    # 3-d lattice test
    #test_lattice       (b1, b2, b3, hmax, kmax, 1,    cdtype=np.float32)
    #lattice_node_radius(b1, b2, b3, hmax, kmax, 1,    cdtype=np.float32)

    #------------------------------
    #return
    #------------------------------

    # binning for look-up table and plots

    # bin parameters for q in units of k = Evald's sphere radius [1/A]
    bpq = BinPars((-0.25, 0.25), 1000, vtype=np.float32, endpoint=False)

    # bin parameters for omega [degree] - fiber rotation angle around axis
    bpomega = BinPars((0.,  180.), 360, vtype=np.float32, endpoint=False)
    
    # bin parameters for beta [degree] - fiber axis tilt angle
    #bpbeta = BinPars((15.,  195.),  2, vtype=np.float32, endpoint=True)
    #bpbeta = BinPars((15.,   15.),  1, vtype=np.float32, endpoint=False)
    #bpbeta = BinPars((5.,    25.),  2, vtype=np.float32, endpoint=True)
    #bpbeta  = BinPars((0.,    180.), 2, vtype=np.float32, endpoint=True)
    bpbeta  = BinPars((0.,      0.), 1, vtype=np.float32, endpoint=True)
    bpbeta2 = BinPars((180.,  180.), 1, vtype=np.float32, endpoint=True)
    #bpbeta  = BinPars((0.,    50.), 11, vtype=np.float32, endpoint=True)
    #bpbeta2 = BinPars((180., 230.), 11, vtype=np.float32, endpoint=True)
    str_beta = 'for-beta:%s' % (bpbeta.strrange)
     
    print '\n%s\nIndexing lookup table\n' % (91*'_')
    #lut  = make_lookup_table(b1, b2, b3, hmax, kmax, lmax, np.float32, evald_rad, sigma_ql, fout, bpq, bpomega, bpbeta)
    lut  = make_lookup_table_v2(b1, b2, b3, hmax, kmax, lmax, np.float32, evald_rad, sigma_ql, sigma_qt, fout, bpq, bpomega, bpbeta)
    lut2 = make_lookup_table_v2(b1, b2, b3, hmax, kmax, lmax, np.float32, evald_rad, sigma_ql, sigma_qt, fout, bpq, bpomega, bpbeta2)

    fout.close()
    print '\nIndexing lookup table is saved in the file: %s' % fname

    #------------------------------
    # produce and save plots
    import algos.graph.GlobalGraphics as gg

    #img = lut2 # or lut2
    img = lut + lut2

    img_range = (bpq.vmin, bpq.vmax, bpomega.vmax, bpomega.vmin) 
    axim = gg.plotImageLarge(lut, img_range=img_range, amp_range=None, figsize=(15,13),\
                      title='Non-symmetrized for beta', origin='upper', window=(0.05,  0.06, 0.94, 0.92), cmap='gray_r')
    axim.set_xlabel('$q_{H}$ ($1/\AA$)', fontsize=18)
    axim.set_ylabel('$\omega$ (degree)', fontsize=18)
    gg.save('%splot-img-prob-omega-vs-qh-%s.png' % (prefix, str_beta), pbits=1)

    axim = gg.plotImageLarge(img, img_range=img_range, amp_range=None, figsize=(15,13),\
                      title='Symmetrized for beta (beta, beta+pi)', origin='upper', window=(0.05,  0.06, 0.94, 0.92), cmap='gray_r') # 'Greys')
    axim.set_xlabel('$q_{H}$ ($1/\AA$)', fontsize=18)
    axim.set_ylabel('$\omega$ (degree)', fontsize=18)
    gg.save('%splot-img-prob-omega-vs-qh-sym-%s.png' % (prefix, str_beta), pbits=1)

    arrhi = np.sum(img,0)    
    fighi, axhi, hi = gg.hist1d(bpq.binedges, bins=bpq.nbins-1, amp_range=(bpq.vmin, bpq.vmax), weights=arrhi,\
                                color='b', show_stat=True, log=False,\
                                figsize=(15,5), axwin=(0.05, 0.12, 0.85, 0.80),\
                                title=None, xlabel='$q_{H}$ ($1/\AA$)', ylabel='Intensity', titwin=None)
    gg.show()

    gg.save_fig(fighi, '%splot-his-prob-vs-qh-%s.png' % (prefix, str_beta), pbits=1)

    qh_weight = zip(bpq.bincenters, arrhi)
    fname = '%sarr-qh-weight-%s.npy' % (prefix, str_beta)
    print 'Save qh:weigt array in file %s' % fname
    np.save(fname, qh_weight)

#------------------------------

if __name__ == "__main__" :
    make_index_table()

#------------------------------
