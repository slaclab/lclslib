#!/usr/bin/env python

#------------------------------
"""Example of analysis script to check look-up table and qh_to_xy transformation
"""
#------------------------------

import sys
import os
import math
import numpy as np
from Detector.GlobalUtils import print_ndarr
import pyimgalgos.GlobalGraphics as gg    
#from pyimgalgos.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, fig_axes, plot_img, plot_peaks_on_img

#------------------------------

R_EVALD = 0.484187 # [1/A]
sigma_qh = 0.003 * R_EVALD

#from pyimgalgos.FiberIndexing import BinPars
from pyimgalgos.HBins import HBins
bpq     = HBins((-0.25, 0.25), 1500)
#bpq     = HBins((-0.25, 0.25), 500)
bpomega = HBins((0., 180.), 360)
 
#------------------------------
#------------------------------

def list_omega_qhrow() :
    """Returns a test list of parameters [(omega, <1-d-array-of-intensities-for-omega>)]
    """
    from time import time
    t0_sec = time()

    lst_omega_qhrow = []

    from pyimgalgos.TDFileContainer import TDFileContainer
    from pyimgalgos.TDNodeRecord import TDNodeRecord
 
    #--- load table of nodes from lookup table

    WDIR = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work'
    fname = os.path.join(WDIR, 'lut-cxif5315-r0169-2016-02-03T15:10:48.txt')
    print 'Load nodes from file %s' % fname

    fc = TDFileContainer(fname, indhdr='index', objtype=TDNodeRecord) #, pbits=256)
    fc.print_content(nlines=50)

    # loop over look-up table of crystal orientations
    for orinum in fc.group_num_iterator() :
        origrp = fc.next()
        # FILLS FOR beta=0 ONLY!
        if origrp()[0].beta > 0 : continue
        omega, row = omega_qhrow(origrp, sigma_qh, bpq)

        lst_omega_qhrow.append((omega, row))

    print 'Consumed time to generate list = %7.3f sec' % (time()-t0_sec)
    print 'Test list of [(omega, <1-d-array-of-intensities-for-omega>)] is generated from lookup table.'

    return lst_omega_qhrow
 
#------------------------------

def omega_qhrow(origrp, sigma_q, bpq) :
    """For each orientation group of records
       returns omega, and histogram array (row) for horizontal q component.
    """
    qhrow = np.zeros((bpq.nbins(),), dtype=np.float32)

    range_q = 3 * sigma_q
    factor_q = -1./(2.*sigma_q*sigma_q)

    beta, omega = 0, 0
    
    # loop over lattice nodes
    for rec in origrp() :
        beta, omega, qh, prob = rec.beta, rec.omega, rec.qh, rec.P
        #print 'omega =%6.1f,  qh(1/A) = %8.4f,  Prob. = %8.4f' % (omega, qh, prob)

        qcenters = bpq.bincenters()
        iqmin, iqmax = bpq.bin_indexes((qh-range_q, qh+range_q))
        for iq in range(iqmin, iqmax+1) :
            dq = qcenters[iq] - qh
            qhrow[iq] = prob * math.exp(factor_q*dq*dq)

    #print 'group omega =%6.1f' % omega

    return omega, qhrow

#------------------------------

def lut_as_image(list_oq) :
    """Returns look-up table as a 2-d image
    """
    img = np.zeros((bpomega.nbins(), bpq.nbins()), dtype=np.float32)
    binw = bpomega.binwidth()
    for omega, row in list_oq :
        #print omega, row[0:10]
        iomega = math.floor(omega/binw)
        if iomega > 359 : continue
        img[iomega,:] += row
    return img

#------------------------------

def plot_lut_as_omega_vs_qh(list_oq) :
    """Plots content of the lookup table as an image of intensities for omega(deg) vs. hq(1/A)
    """
    img = lut_as_image(list_oq)
    print_ndarr(img, 'img')

    img_range = (bpq.vmin(), bpq.vmax(), bpomega.vmax(), bpomega.vmin()) 
    axim = gg.plotImageLarge(img, img_range=img_range, amp_range=None, figsize=(15,13),\
                      title='Plot reconstructed from look-up table', origin='upper',\
                             window=(0.06,  0.06, 0.94, 0.92), cmap='gray_r')
    axim.set_xlabel('$q_{H}$ ($1/\AA$)', fontsize=18)
    axim.set_ylabel('$\omega$ (degree)', fontsize=18)
    gg.save('img-lut-prob-omega-vs-qh.png', pbits=1)
    gg.show('do not block')

    #import matplotlib.pyplot  as plt
    #plt.imshow(img)
    #plt.show()
    
#------------------------------
#------------------------------
#------------------------------
#------------------------------

def qh_to_xy(qh, R) :
    """Returns reciprocal (xe,ye) coordinates of the qh projection on Evald sphere.
       qh - (numpy array) horizontal component of q values (1/A)
       R  - (float scalar) Evald sphere radius (1/A)
       Assuming that center of the Evald sphere is in (-R,0); qh is oriented along y.
       NOTE: qh, L, sina, cosa, xe, ye - are the numpy arrays of the same shape as qh
    """
    sina = qh/R
    cosa = np.sqrt(1.-sina*sina)
    xe = R * (cosa-1.)
    return  xe, qh

    #L = np.sqrt(R*R + qh*qh) 
    #sina = qh/L
    #cosa = R/L
    #xe = R * (cosa-1.)
    #ye = R * sina
    #return xe, ye

#------------------------------

def xy_lattice_image(list_oq) :
    """Returns 2-d image of the crystal lattice in x-y reciprocal (x,y) space.
    """

    from pyimgalgos.FiberAngles import rotation

    img = np.zeros((bpq.nbins(), bpq.nbins()), dtype=np.float32)

    qh = bpq.bincenters()
    xe, ye = qh_to_xy(qh, R_EVALD)
    print_ndarr(qh, 'qh')

    for omega, row in list_oq :

        #if omega%10 > 0 : continue 
        #print omega, #, row[0:10]

        xrot, yrot = rotation(xe, ye, -omega)

        iX = bpq.bin_indexes(xrot)
        iY = bpq.bin_indexes(yrot)

        #print_ndarr(xrot, 'xrot')
        #print_ndarr(iX, 'iX')
        #print_ndarr(iY, 'iY')
        img[iX,iY] += row

    return img

#------------------------------

def arr_2d_gauss(rank=2, sigma=1.5) :
    """returns 2-d Gaussian distribution centred in the center of square matrix of shape=(2*rank+1,2*rank+1)
    """    
    rank1 = rank+1 if rank>1 else 2
    arrq1 = np.zeros(shape=(rank1,rank1), dtype=np.float32)    
    f = -0.5/(sigma*sigma)
    for r in range(rank1) :
        for c in range(rank1) :
            arrq1[r,c] = math.exp(f*(r*r+c*c))
    arrbot = np.hstack([arrq1 [:,:0:-1],arrq1])
    arr    = np.vstack([arrbot[:0:-1,:],arrbot])
    print '2-d Gaussian array of shape %s\n' % (str(arr.shape)), arr
    return arr

#------------------------------

def plot_xy_lattice(list_oq) :
    """Plots image of the crystal lattice, using list of [(omega,<1-d-array-of-intensities-for-omega>)]
    """

    img = xy_lattice_image(list_oq)
    print_ndarr(img, 'img')

    #--- Convolution of image 
    from scipy.signal import convolve2d
    g2d = arr_2d_gauss(2, 1.5)
    img = convolve2d(img, g2d, mode='same', boundary='fill', fillvalue=0)
    #---

    img_range = (bpq.vmin(), bpq.vmax(), bpq.vmin(), bpq.vmax()) 
    axim = gg.plotImageLarge(img, img_range=img_range, amp_range=None, figsize=(15,13),\
                      title='Lattice', origin='upper',\
                             window=(0.08,  0.06, 0.94, 0.92)) # , cmap='gray_r')
    axim.set_xlabel('Reciprocal x ($1/\AA$)', fontsize=18)
    axim.set_ylabel('Reciprocal y ($1/\AA$)', fontsize=18)
    gg.save('img-lut-lattice-xy.png', pbits=1)
    gg.show()

#------------------------------

if __name__ == "__main__" :
    
    list_oq = list_omega_qhrow()
    plot_lut_as_omega_vs_qh(list_oq)
    plot_xy_lattice(list_oq)
    sys.exit('Done')

#------------------------------
# EOF
#------------------------------


    


