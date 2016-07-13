#!/usr/bin/env python

##-----------------------------
"""Example of analysis script:
   processes cxif5315-r0169 indexed data file;
   loop over all entries of the indexed data file and creates a couple of plots for 
   1) intensity in h,k space (also print_crystal_in_hk_space),
   2) q_y vs q_x intensity in reciprocal space.
"""
##-----------------------------

import sys
import os
import math
import numpy as np
from time import strftime, localtime #, gmtime, time

from pyimgalgos.TDFileContainer import TDFileContainer
from pyimgalgos.TDMatchRecord   import TDMatchRecord

from pyimgalgos.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, fig_axes, plot_img, plot_peaks_on_img

#from pyimgalgos.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, fig_axes, plot_img, plot_peaks_on_img
#from pyimgalgos.FiberAngles import fraser, rotation_phi_beta, calc_phi, calc_beta, funcy#, fraser_xyz

#------------------------------

class Storage :
    """Store for shared parameters."""

    DO_HIST = True            # save histograms
    DO_PLOT = True            # save plots

    def __init__(self) :
        self.exp  = None
        self.run  = None
        self.tsec = None
        self.nda  = None

        self.h1d_omega_ind = []
        self.h1d_beta_ind  = []
        self.h1d_nmatched  = []
        self.h1d_npeak_nm  = []
        self.h1d_nnode_nm  = []

sp = Storage()

#------------------------------

def store_exp_run_info(rec) :
    """Preserve a few record parameters to compose the output file name
    """
    if sp.exp is None :
        sp.exp  = rec.exp
        sp.run  = rec.run
        sp.tsec = rec.tsec
        sp.tstamp  = strftime('%Y-%m-%dT%H:%M:%S', localtime())
        sp.exp_run = '%s-r%04d' % (sp.exp, sp.run)

#------------------------------

def peak_dist(w=5) :
    """ returns np.array with shape=(2*w+1,2*w+1) with 2-d Gaussian distribution 
    """
    if sp.nda is None :
        sigma = 0.5*w
        f = 1./(sigma*sigma)    
        sp.nda = np.zeros((2*w+1, 2*w+1), dtype=np.float32)
        for x in range(-w,w+1) :
           for y in range(-w,w+1) :
               sp.nda[x+w,y+w] = math.exp(-f*(x*x+y*y)) 
    return sp.nda

#------------------------------

def add_peak_to_img(img, bpq, r, sina, cosa, amp, w=10) :
    """Add peak (r, sina, cosa, amp) to the 2-d numpy array img with bin parameters bpq 
    """
    ix = (r*cosa - bpq.vmin)/bpq.binwidth
    iy = (r*sina - bpq.vmin)/bpq.binwidth
    img[ix-w:ix+w+1, iy-w:iy+w+1] += amp*peak_dist(w)

#------------------------------

def proc_file(fname, hmax=6, kmax=6, bpq=None) :
    """Process file with indexing records (matching of data peaks to crystal nodes) 
    """
    sp.fc = TDFileContainer(fname, indhdr='Evnum', objtype=TDMatchRecord) #, pbits=256)
    sp.fc.print_content(nlines=20)

    img_space = np.zeros((2*hmax+1, 2*kmax+1), dtype=np.float32)
    img_recip = np.zeros((bpq.nbins, bpq.nbins), dtype=np.float32)

    for evnum in sp.fc.group_num_iterator() :
        event = sp.fc.next()
        lst_recs = event.get_objs()
        rec0 = lst_recs[0]
        store_exp_run_info(rec0)
        print '%s\nEvent# %6d exp=%s:run=%d %s %s' % (80*'_', evnum, rec0.exp, rec0.run, rec0.date, rec0.time)
        print 'HDR:', sp.fc.hdr
        omega_rad = math.radians(rec0.omega)
        sina = math.sin(omega_rad)
        cosa = math.cos(omega_rad)

        # accumulate data on plots
        for rec in lst_recs :

            print 'REC:', rec.line[:-1]
            h, k, amp = rec.h, rec.k, rec.peak_signal()
            #print rec.h, rec.k, rec.peak_signal(), rec.phi_fit, rec.beta_fit, rec.beta, rec.omega
            print 'state=%s  omega=%.2f  R=%.6f  qh=%.6f  qh_fit=%.6f' % (rec.state, rec.omega, rec.R, rec.qh, rec.qh_fit)

            if not rec.state=='MATCHED' : continue
            #if not rec.state=='PEAK-NM' : continue
            ##if not rec.state=='NODE-NM' : continue

            img_space[h+hmax, k+kmax] += amp
            r = rec.qh_fit
            if not math.fabs(r)<bpq.vmax : continue

            add_peak_to_img(img_recip, bpq, r, sina, cosa, amp, w=10)

        # fill histogram arrays for event
        sp.h1d_omega_ind.append(rec.omega)
        sp.h1d_beta_ind.append(rec.beta)

        nmatched = np.array([1 for rec in lst_recs if rec.state=='MATCHED'], dtype=np.int16).sum()
        npeak_nm = np.array([1 for rec in lst_recs if rec.state=='PEAK-NM'], dtype=np.int16).sum()
        nnode_nm = np.array([1 for rec in lst_recs if rec.state=='NODE-NM'], dtype=np.int16).sum()

        sp.h1d_nmatched.append(nmatched)
        sp.h1d_npeak_nm.append(npeak_nm)
        sp.h1d_nnode_nm.append(nnode_nm)

    img_space /= img_space.max()
    img_recip /= img_recip.max()
    return img_space, img_recip

#------------------------------

def print_crystal_in_hk_space(img, hmax, kmax) :
    """ Prints img of unit-normalized intensities as a 2-d table for h and k indexes
    """
    print 'Crystal in h,k space:'
    print 'k   :',
    for k in range(-kmax, kmax+1) : print '%8d' % k,
    for h in range(-hmax, hmax+1) :        
        print '\nh=%2d:' % h,
        for k in range(-kmax, kmax+1) : print '%8.4f' % (img[h+hmax,k+kmax]),
    print ''

#------------------------------

def h1d(hlst, bins=None, amp_range=None, weights=None, color=None, show_stat=True, log=False,\
        figsize=(6,5), axwin=(0.15, 0.12, 0.78, 0.80), title='Title', xlabel='x', ylabel='y', titwin=None, fnm='fnm.png') :
    """Wrapper for hist1d, move, and save methods, using common store parameters
    """
    fig, axhi, hi = hist1d(np.array(hlst), bins, amp_range, weights, color, show_stat,\
                           log, figsize, axwin, title, xlabel, ylabel, titwin)

    #move(sp.hwin_x0y0[0], sp.hwin_x0y0[1])
    save('%s-%s' % (sp.prefix, fnm), sp.DO_HIST)
    return fig, axhi, hi

#------------------------------

def plot_histograms() :

    #---------
    h1d(np.array(sp.h1d_omega_ind), bins=360, amp_range=(0.,180.), \
        title ='Equ: after indexing', xlabel='$\omega$ (deg)', ylabel='Events', fnm='equ-omega.png')
    #---------
    h1d(np.array(sp.h1d_beta_ind), bins=130, amp_range=(0.,260.), \
        title ='Equ: after indexing', xlabel='beta (deg)', ylabel='Events', fnm='equ-beta.png')
    #---------
    h1d(np.array(sp.h1d_nmatched), bins=10, amp_range=(-0.5, 9.5), \
        title ='N matched peaks', xlabel='N matched peaks', ylabel='Events', fnm='nmatched.png')
    #---------
    h1d(np.array(sp.h1d_npeak_nm), bins=10, amp_range=(-0.5, 9.5), \
        title ='N peak non-matched', xlabel='N peak non-matched', ylabel='Events', fnm='npeak-nm.png')
    #---------
    h1d(np.array(sp.h1d_nnode_nm), bins=10, amp_range=(-0.5, 9.5), \
        title ='N node non-matched', xlabel='N node non-matched', ylabel='Events', fnm='nnode-nm.png')
    #---------

#------------------------------

def do_main() :
    """ Main method to do work
    """
    from pyimgalgos.FiberIndexing import BinPars
    from pyimgalgos.GlobalUtils import create_directory

    # h-k space image parameters
    hmax = 4
    kmax = 6

    # recipical space image parameters, bins for 2-d image 
    bpq = BinPars((-0.25, 0.25), 1200, vtype=np.float32, endpoint=True)

    #fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/peak-idx-cxif5315-r0169-2015-11-13T17:04:37.txt'
    #fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/peak-idx-cxif5315-r0169-2015-12-01T15:44:49.txt'
    fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/peak-idx-cxif5315-r0169-2016-05-12T18:17:57.txt'

    rdir = './results-idx'
    create_directory(rdir)

    sp.prefix  = '%s/2016-05-13-v01-idx-res-matched' % rdir # file name prefix for histograms
    #sp.prefix  = '%s/2016-05-13-v01-idx-res-peak-nm' % rdir # file name prefix for histograms

    img_space, img_recip = proc_file(fname, hmax, kmax, bpq)
    print_crystal_in_hk_space(img_space, hmax, kmax)

    print 'img_recip.shape=', img_recip.shape


    if sp.DO_HIST : plot_histograms()


    if sp.DO_PLOT :

        import pyimgalgos.GlobalGraphics as gg

        img = img_space
        img_range=(-kmax-0.5, kmax+0.5, -hmax-0.5, hmax+0.5)
        axim = gg.plotImageLarge(img, img_range=img_range, amp_range=(0,1), figsize=(8,6),\
                      title='Crystal structure in h-k space', origin='upper', window=(0.1, 0.1, 0.9, 0.86))
        axim.set_xlabel('k index', fontsize=18)
        axim.set_ylabel('h index', fontsize=18)

        gg.savefig('%s-%s-crystal-in-hk-space.png' % (sp.prefix, sp.exp_run)) # sp.tstamp


        img = img_recip
        ave, rms = img.mean(), img.std()
        amin, amax = 0, ave+5*rms
        img_range=(bpq.vmin, bpq.vmax, bpq.vmin, bpq.vmax)
        axim = gg.plotImageLarge(img, img_range=img_range, amp_range=(amin, amax), figsize=(10,8),\
                  title='Crystal structure in reciprocal space', origin='upper', window=(0.1, 0.1, 0.86, 0.86))
        axim.set_xlabel('$q_x$ ($1/\AA$)', fontsize=18)
        axim.set_ylabel('$q_y$ ($1/\AA$)', fontsize=18)

        gg.savefig('%s-%s-crystal-in-recip-space.png' % (sp.prefix, sp.exp_run))

        gg.show()

        #np.save(fname, qh)

#------------------------------
if __name__ == "__main__" :
    do_main()
    sys.exit('Processing is completed')

#------------------------------
