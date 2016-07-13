#!/usr/bin/env python

##-----------------------------
"""Example of analysis script:
   processes cxif5315-r0169 peak data from file, select good peaks, good events,
   find fiber phi and beta angles, apply image phi rotation and fraser transformation to q-space,
   plot image (if and which requested), accumulate/plot/save histograms (if all events are processed).
   NOTE: in v3 fraser transformation to q-space is applied for beta=0 !!!

Usage::
   # This script uses peak data from file in directory wdir, rdir,
   # so run it from release directory, which have users package cxif5315:
   python ./cxif5315/proc-cxif5315-r0169-peaks-from-file.py
    
   # Parameters to edit:
   EVTMAX, SKIP, DO_PLOT, DO_PHI_ROTATION, DO_FRASER, wdir, rdir, file names etc.
"""
##-----------------------------

import sys
import os
import math
import numpy as np
from time import time
from scipy.optimize import curve_fit
from scipy.stats import chi2

from algos.core.GlobalUtils import print_ndarr, subtract_bkgd
from algos.graph.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, fig_axes, plot_img, plot_peaks_on_img
from algos.core.PeakStore import PeakStore
from algos.diffraction.FiberAngles import fraser, rotation_phi_beta, calc_phi, calc_beta, funcy, funcy_l0, funcy_l1_v0, funcy_l1_v1#, fraser_xyz
from algos.histo.RadialBkgd import RadialBkgd, polarization_factor

#------------------------------
#------------------------------
#------------------------------
#------------------------------
class Storage :
    """Store for shared parameters. Assumed to be used as a singleton
    """
    PIXEL_SIZE  = 109.92            # pixel size [um]
    DETD        = 913.27            # sample to detector distance in number of pixels
    DETD_um     = DETD*PIXEL_SIZE   # sample to detector distance [um]

    WAVE_NUMBER = 0.484187          # [1/A]
    SCALE       = WAVE_NUMBER       # 
    #-------
    SIGMA_Q     = 0.002*WAVE_NUMBER # [1/A] resolution in q
    #-------

    DQHMAX      = 3*SIGMA_Q         # [1/A] dqh_max for evaluation of probability
    FACTOR  = -1./(2*SIGMA_Q*SIGMA_Q)
    QVMAX       = 0.02              # [1/A] qv_max - event is descarded from indexing if any peak has |qv| > QVMAX

    PROBTHR     = 0.001             # threshold on probability to consider peak indexed

    USE_PEAK_SELECTION = True

    BETA_FROM_ARC = True
    BETA_FROM_ARC = False
    
    ###--- For plots only:
    PLOT_HIST = False
    PLOT_HIST = True

    #DO_PLOT = True
    DO_PLOT = False

    SKIP    = 0 # should be >=0
    #EVTMAX  = 100 + SKIP
    EVTMAX  = 100000 + SKIP

    DO_PHI_ROTATION = False
    DO_PHI_ROTATION = True

    DO_FRASER  = False
    DO_FRASER  = True

    ###------------------- 

    winds_bkgd = [(s, 10, 100, 270, 370) for s in (4,12,20,28)] # use part of segments 4 and 20 to subtr bkgd

    #xoffset, yoffset = 0, 0
    #xsize,   ysize   = 1850, 1850
    #xsize,   ysize   = 2000, 2000

    xy0_off_pix = (1000,1000)       # offset of the detector center in pixels
    #xoffset, yoffset = 400, 400    # offset of the image lowest corner
    #xsize,   ysize   = 1150, 1150  # image size

    # Image trim parameters
    xoffset, yoffset = 200, 200    # offset of the image lowest corner
    xsize,   ysize   = 1600, 1600  # image size

    # parameters of the fraser transformation
    img_center = (xy0_off_pix[0] - xoffset,\
                  xy0_off_pix[1] - yoffset)

    fras_shape = (1200,1200)


    def __init__(sp) :
        """Init parameters for histograms and counters
        """
        sp.pstore_idx    = None
        sp.pstore_fit    = None

        sp.count_plt = 0
        sp.psana_and_graphics_is_initialized = False
        sp.fit_phi0  = 0
        sp.img_sum   = None
        sp.img_num   = 0

        sp.count_evt_sel = 0
        sp.count_evt_idxed = 0
    
        sp.lst_plt_evt_peaks  = []
    
        # Arc region
        sp.count_arc_pks_tot  = 0
        sp.count_arc_pks_sel  = 0
        sp.count_arc_2peaks   = 0
    
        sp.lst_arc_evt_peaks  = []
        sp.lst_arc_1peak_r    = []
        sp.lst_arc_1peak_phi  = []
        sp.lst_arc_2peak_dist = []
        sp.lst_arc_2peak_dr   = []
        sp.lst_arc_2peak_phi  = []
        sp.lst_arc_npktot     = []
        sp.lst_arc_npksel     = []
        sp.lst_arc_amax       = []
        sp.lst_arc_atot       = []
        sp.lst_arc_npix       = []
        sp.lst_arc_r          = []
        sp.lst_arc_phi        = []
        sp.lst_arc_bkgd       = []
        sp.lst_arc_rms        = []
        sp.lst_arc_son        = []
        sp.lst_arc_sonc       = []

        # Arc phi-beta fit parameters
        sp.lst_arc_sel_2dist  = []
        sp.lst_arc_sel_2dcen  = []
        sp.lst_arc_ev_p0      = []
        sp.lst_arc_ev_p1      = []
        sp.lst_arc_ev_dp0     = []
        sp.lst_arc_ev_dp1     = []

        sp.lst_arc_phi_err    = []
        sp.lst_arc_beta_err   = []
        sp.lst_arc_prob       = []

    
        # Equator region
        sp.count_equ_2peaks   = 0
        sp.count_equ_pks_tot  = 0
        sp.count_equ_pks_sel  = 0
    
        sp.lst_equ_evt_peaks  = []
        sp.lst_equ_npktot     = []
        sp.lst_equ_npksel     = []
        sp.lst_equ_amax       = []
        sp.lst_equ_atot       = []
        sp.lst_equ_npix       = []
        sp.lst_equ_r          = []
        sp.lst_equ_r_raw      = []
        sp.lst_equ_phi        = []
        sp.lst_equ_bkgd       = []
        sp.lst_equ_rms        = []
        sp.lst_equ_son        = []
        sp.lst_equ_sonc       = []
        sp.lst_equ_r1         = []
        sp.lst_equ_r2         = []
        sp.lst_equ_phi1       = []
        sp.lst_equ_phi2       = []
        sp.lst_equ_phi1c      = []
        sp.lst_equ_phi2c      = []
        sp.lst_equ_npks_2inarc= []
        sp.lst_equ_dphicmb    = []
        sp.lst_equ_dphi000    = []
        sp.lst_equ_dphi180    = []
    
        sp.lst_equ_ev_phi     = []
        sp.lst_equ_ev_beta    = []

        sp.lst_equ_phi_err    = []
        sp.lst_equ_beta_err   = []
        sp.lst_equ_prob       = []
    
        sp.lst_equ_ev_p0        = []
        sp.lst_equ_ev_p1        = []
        sp.lst_equ_ev_dp0       = []
        sp.lst_equ_ev_dp1       = []
        sp.lst_equ_r_fit_selec  = []
        sp.lst_equ_qh_fit_selec = []
        sp.lst_equ_qv_fit_selec = []
        sp.lst_equ_omega        = []
        sp.lst_equ_beta         = []
        sp.lst_equ_probsel      = []
        sp.lst_equ_npk_idxed    = []

#------------------------------
#------------------------------
#------------------------------

# create singleton object with shared parameters
sp = Storage()

#------------------------------
#------------------------------
#------------------------------

def h1d(hlst, bins=None, amp_range=None, weights=None, color=None, show_stat=True, log=False,\
        figsize=(6,5), axwin=(0.15, 0.12, 0.78, 0.80), title='Title', xlabel='x', ylabel='y', titwin=None, fnm='fnm.png') :
    """Wrapper for hist1d, move, and save methods, using common store parameters
    """
    fig, axhi, hi = hist1d(np.array(hlst), bins, amp_range, weights, color, show_stat,\
                           log, figsize, axwin, title, xlabel, ylabel, titwin)

    move(sp.hwin_x0y0[0], sp.hwin_x0y0[1])
    save('%s-%s' % (sp.prefix, fnm), sp.do_save)
    return fig, axhi, hi

#------------------------------
def plotResultsArc(prefix='plot', do_save=True, hwin_x0y0 = (0,0)) :
    """Plots/saves histograms for raw parameters in Arc region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0
    #---------
    h1d(np.array(sp.lst_arc_1peak_r), bins=60, amp_range=(425,455), \
           title ='Arc: 1-peak radius', xlabel='r1 [pixel]', ylabel='Peaks',\
           fnm='arc-1peak-r.png')
    #---------
    #h1d(np.array(sp.lst_arc_1peak_phi), bins=360, amp_range=(-360,360), \
    #       title ='Arc: angle of 1 peak', xlabel='phi1[deg]', ylabel='Events',\
    #       fnm='arc-1peak-phi.png')
    #---------
    h1d(np.array(sp.lst_arc_1peak_phi), bins=160, amp_range=(60,140), \
           title ='Arc: angle of 1 peak', xlabel='phi1[deg]', ylabel='Events',\
           fnm='arc-1peak-phi-zoomed.png')
    #---------
    h1d(np.array(sp.lst_arc_2peak_dist), bins=200, amp_range=(1,401), \
           figsize=(12,5), axwin=(0.08, 0.12, 0.88, 0.80), \
           title ='Arc: Distance between 2 peaks', xlabel='Distance[pixel]', ylabel='Events',\
           fnm='arc-2peak-dist.png')
    #---------
    h1d(np.array(sp.lst_arc_2peak_dr), bins=40, amp_range=(-20,20), \
           title ='Arc: dr between 2 peaks', xlabel='dr[pixel]', ylabel='Events',\
           fnm='arc-2peak-dr.png')
    #---------
    #h1d(np.array(sp.lst_arc_2peak_phi), bins=360, amp_range=(-360,360), \
    #       title ='Arc: Orientation angle from 2 peaks', xlabel='<phi>[deg]', ylabel='Events',\
    #       fnm='arc-2peak-phi.png')
    #---------
    h1d(np.array(sp.lst_arc_2peak_phi), bins=80, amp_range=(60,140), \
           title ='Arc: Orientation angle from 2 peaks', xlabel='<phi>[deg]', ylabel='Events',\
           fnm='arc-2peak-phi-zoomed.png')
    #---------
    #return
    #---------
    h1d(np.array(sp.lst_arc_npktot), bins=25, amp_range=(-0.5,24.5), \
           title ='Arc: Number of peaks after peakfinder', xlabel='N peaks in event', ylabel='Events',\
           fnm='arc-npkraw.png')
    #---------
    h1d(np.array(sp.lst_arc_amax), bins=200, amp_range=(0,1000), \
           title ='Arc: Amax', xlabel='Amax [ADU]', ylabel='Peaks',\
           fnm='arc-amax.png')
    #---------
    h1d(np.array(sp.lst_arc_atot), bins=200, amp_range=(0,10000), \
           title ='Arc: Atot', xlabel='Atot [ADU]', ylabel='Peaks',\
           fnm='arc-atot.png')
    #---------
    h1d(np.array(sp.lst_arc_npix), bins=200, amp_range=(0,400), \
           title ='Arc: Number of pixels/peak', xlabel='Npix', ylabel='Peaks',\
           fnm='arc-npix.png')
    #---------
    #h1d(np.array(sp.lst_arc_npix), bins=150, amp_range=(200,500), \
    #       title ='Arc: Number of pixels (ext.)', xlabel='Npix', ylabel='Peaks',\
    #       fnm='arc-npix-ext.png')
    #---------
    h1d(np.array(sp.lst_arc_r), bins=60, amp_range=(425,455), \
           title ='Arc: Peak radius', xlabel='r [pixel]', ylabel='Peaks',\
           fnm='arc-rpeak.png')
    #---------
    h1d(np.array(sp.lst_arc_phi), bins=360, amp_range=(-180,180), \
           figsize=(12,5), axwin=(0.10, 0.12, 0.83, 0.80), \
           title ='Arc: Peak azimuthal angle', xlabel='phi [rad]', ylabel='Peaks',\
           fnm='arc-phi.png')
    #---------
    h1d(np.array(sp.lst_arc_npksel), bins=15, amp_range=(-0.5,14.5), \
           title ='Arc: Number of peaks selected', xlabel='N peaks selected', ylabel='Events',\
           fnm='arc-npksel.png')
    #---------
    #---------
    h1d(np.array(sp.lst_arc_bkgd), bins=100, amp_range=(-100,100), \
           title ='Arc: bkgd', xlabel='bkgd', ylabel='Events',\
           fnm='arc-bkgd.png')
    #---------
    h1d(np.array(sp.lst_arc_rms), bins=80, amp_range=(0,80), \
           title ='Arc: rms', xlabel='rms', ylabel='Events',\
           fnm='arc-rms.png')
    #---------
    h1d(np.array(sp.lst_arc_son), bins=100, amp_range=(0,100), \
           title ='Arc: S/N', xlabel='S/N', ylabel='Events',\
           fnm='arc-son.png')
    #---------

#------------------------------
def plotResultsEqu(prefix='plot', do_save=True, hwin_x0y0 = (0,600)) :
    """Plots/saves histograms for raw parameters in Equ region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0
    #---------
    h1d(np.array(sp.lst_equ_amax), bins=200, amp_range=(0,1000), \
           title ='Equ: Amax', xlabel='Amax [ADU]', ylabel='Peaks',\
           fnm='equ-amax.png')
    #---------
    h1d(np.array(sp.lst_equ_atot), bins=200, amp_range=(0,10000), \
           title ='Equ: Atot', xlabel='Atot [ADU]', ylabel='Peaks',\
           fnm='equ-atot.png')
    #---------
    h1d(np.array(sp.lst_equ_npix), bins=200, amp_range=(0,400), \
           title ='Equ: Number of pixels/peak', xlabel='Npix', ylabel='Peaks',\
           fnm='equ-npix.png')
    #---------
    h1d(np.array(sp.lst_equ_r), bins=300, amp_range=(0,600), \
           title ='Equ: Peak radius', xlabel='r [pixel]', ylabel='Peaks',\
           fnm='equ-rpeak.png')
    #---------
    h1d(np.array(sp.lst_equ_r_raw), bins=300, amp_range=(0,600), \
           title ='Equ: Peak radius', xlabel='r [pixel]', ylabel='Peaks',\
           fnm='equ-rpeak-raw.png')
    #---------
    h1d(np.array(sp.lst_equ_phi), bins=360, amp_range=(-180,180), \
           figsize=(12,5), axwin=(0.10, 0.12, 0.83, 0.80), \
           title ='Equ: Peak azimuthal angle', xlabel='phi [deg]', ylabel='Peaks',\
           fnm='equ-phi.png')
    #---------
    h1d(np.array(sp.lst_equ_dphicmb), bins=90, amp_range=(-45,45), \
           #figsize=(12,5), axwin=(0.10, 0.12, 0.83, 0.80), \
           title ='Equ: Peak azimuthal angle combined diff.', xlabel='dphi comb.[deg]', ylabel='Peaks',\
           fnm='equ-dphicmb.png')
    #---------
    h1d(np.array(sp.lst_equ_dphi000), bins=90, amp_range=(-45,45), \
           #figsize=(12,5), axwin=(0.10, 0.12, 0.83, 0.80), \
           title ='Equ: Peak azimuthal angle relative 0', xlabel='dphi [deg]', ylabel='Peaks',\
           fnm='equ-dphi000.png')
    #---------
    h1d(np.array(sp.lst_equ_dphi180), bins=90, amp_range=(-45,45), \
           #figsize=(12,5), axwin=(0.10, 0.12, 0.83, 0.80), \
           title ='Equ: Peak azimuthal angle relative 180', xlabel='dphi [deg]', ylabel='Peaks',\
           fnm='equ-dphi180.png')
    #---------

    #---------
    h1d(np.array(sp.lst_equ_bkgd), bins=100, amp_range=(-100,100), \
           title ='Equ: bkgd', xlabel='bkgd', ylabel='Events',\
           fnm='equ-bkgd.png')
    #---------
    h1d(np.array(sp.lst_equ_rms), bins=80, amp_range=(0,80), \
           title ='Equ: rms', xlabel='rms', ylabel='Events',\
           fnm='equ-rms.png')
    #---------
    h1d(np.array(sp.lst_equ_son), bins=100, amp_range=(0,100), \
           title ='Equ: S/N', xlabel='S/N', ylabel='Events',\
           fnm='equ-son.png')
    #---------
    #---------
    h1d(np.array(sp.lst_equ_npktot), bins=20, amp_range=(-0.5,19.5), \
           title ='Equ: Number of peaks after peakfinder', xlabel='N peaks in event', ylabel='Events',\
           fnm='equ-npkraw.png')
    #---------
    h1d(np.array(sp.lst_equ_npksel), bins=10, amp_range=(-0.5,9.5), \
           title ='Equ: Number of peaks selected', xlabel='N peaks selected', ylabel='Events',\
           fnm='equ-npksel.png')
    #---------
    h1d(np.array(sp.lst_equ_npks_2inarc), bins=15, amp_range=(-0.5,14.5), \
           title ='Equ: N peaks in Equ  @ 2 in Arc', xlabel='N peaks in Equ @ 2 in Arc', ylabel='Events',\
           fnm='equ-npks-2inarc.png')
    #---------

    return

    #---------
    h1d(np.array(sp.lst_equ_phi1c), bins=100, amp_range=(-150,-50), \
           title ='Equ: Azimuthal angle at phi<0 corrected on orientation', xlabel='phi [deg]', ylabel='Peaks',\
           fnm='equ-phi1c.png')
    #---------
    h1d(np.array(sp.lst_equ_phi2c), bins=100, amp_range=(50,150), \
           title ='Equ: Azimuthal angle at phi>0 corrected on orientation', xlabel='phi [deg]', ylabel='Peaks',\
           fnm='equ-phi2c.png')
    #---------
    h1d(np.array(sp.lst_equ_phi1), bins=100, amp_range=(-150,-50), \
           title ='Equ: Azimuthal angle at phi<0', xlabel='phi [deg]', ylabel='Peaks',\
           fnm='equ-phi1.png')
    #---------
    h1d(np.array(sp.lst_equ_phi2), bins=100, amp_range=(50,150), \
           title ='Equ: Azimuthal angle at phi>0', xlabel='phi [deg]', ylabel='Peaks',\
           fnm='equ-phi2.png')
    #---------
    h1d(np.array(sp.lst_equ_r1), bins=200, amp_range=(0,600), \
           title ='Equ: Peak radius at phi<0', xlabel='r1 [pixel]', ylabel='Peaks',\
           fnm='equ-r1.png')
    #---------
    h1d(np.array(sp.lst_equ_r2), bins=200, amp_range=(0,600), \
           title ='Equ: Peak radius at phi>0', xlabel='r2 [pixel]', ylabel='Peaks',\
           fnm='equ-r2.png')
    #---------

#------------------------------
def plotResultsEquSel(prefix='plot', do_save=True, hwin_x0y0 = (0,600)) :
    """Plots/saves histograms for selected events in EQU region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0
    #---------
    h1d(np.array(sp.lst_equ_ev_phi), bins=80, amp_range=(-40,40), \
           title ='Equ: event angle phi for 2-peak events', xlabel='$\phi$ (deg)', ylabel='Events',\
           fnm='equ-ev-phi.png')
    #---------
    h1d(np.array(sp.lst_equ_ev_beta), bins=100, amp_range=(-60,40), \
           title ='Equ: event angle beta for 2-peak events', xlabel='beta (deg)', ylabel='Events',\
           fnm='equ-ev-beta.png')
    #---------
    #---------
    h1d(np.array(sp.lst_equ_r_fit_selec), bins=500, amp_range=(0,500), figsize=(12,5), \
           title ='Equ: selected after fit', xlabel='r (pix)', ylabel='Events',\
           fnm='equ-ev-fit-selec-r.png')
    #---------
    h1d(np.fabs(sp.lst_equ_qh_fit_selec), bins=500, amp_range=(0.,0.25), figsize=(12,5), \
           title ='Equ: selected after fit', xlabel='$|q_{H}|$ ($1/\AA$)', ylabel='Events',\
           fnm='equ-ev-fit-selec-abs-qh.png')
    #---------
    arr_lut_qh_prob = lut_qh_hist_array()
    fig3, axhi3, hi3 = h1d(np.fabs(arr_lut_qh_prob[:,0]), bins=500, amp_range=(0.,0.25), weights=arr_lut_qh_prob[:,1], figsize=(12,5), \
           title ='|qh| in look-up table', xlabel='$|q_{H}|$ ($1/\AA$)', ylabel='$\propto$ Probability',\
           fnm='equ-lut-abs-qh.png')
    arr, bins, patches = hi3
    bincenters = bins + 0.5*(bins[1]-bins[0])
    #---------
    fig4, axhi4, hi4 = h1d(np.fabs(sp.lst_equ_qh_fit_selec), bins=500, amp_range=(0.,0.25), figsize=(12,5), \
           title ='Comparison for |qh|', xlabel='$|q_{H}|$ ($1/\AA$)', ylabel='Events',\
           fnm='equ-ev-fit-selec-abs-qh.png')
    f = hi4[0].max()/arr.max()
    axhi4.plot(bincenters[:-1], f*arr,'r-')
    save('%s-%s' % (sp.prefix, 'equ-ev-fit-selec-abs-qh-comp.png'), sp.do_save)
    #---------
    fig0, axhi0, hi0 = h1d(np.array(sp.lst_equ_qh_fit_selec), bins=1000, amp_range=(-0.25,0.25), figsize=(24,5), \
           title ='Equ: selected after fit', xlabel='$q_{H}$ ($1/\AA$)', ylabel='Events',\
           fnm='equ-ev-fit-selec-qh.png')
    #---------
    fig2, axhi2, hi2 = h1d(arr_lut_qh_prob[:,0], bins=1000, amp_range=(-0.25,0.25), weights=arr_lut_qh_prob[:,1], figsize=(24,5), \
           title ='qh in look-up table', xlabel='$q_{H}$ ($1/\AA$)', ylabel='$\propto$ Probability',\
           fnm='equ-lut-qh.png')
    arr, bins, patches = hi2
    bincenters = bins + 0.5*(bins[1]-bins[0])
    #---------
    fig1, axhi1, hi1 = h1d(np.array(sp.lst_equ_qh_fit_selec), bins=1000, amp_range=(-0.25,0.25), figsize=(24,5), \
           title ='Comparison of data with prediction', xlabel='$q_{H}$ ($1/\AA$)', ylabel='Events',\
           fnm='equ-ev-fit-selec-qh-comp.png')
    f = hi1[0].max()/arr.max()
    axhi1.plot(bincenters[:-1], f*arr,'r-')
    save('%s-%s' % (sp.prefix, 'equ-ev-fit-selec-qh-comp.png'), sp.do_save)
    #---------
    h1d(np.array(sp.lst_equ_qv_fit_selec), bins=100, amp_range=(-0.04,0.01), \
           title ='Equ: selected after fit', xlabel='$q_{V}$ ($1/\AA$)', ylabel='Events',\
           fnm='equ-ev-fit-selec-qv.png')
    #---------
    h1d(np.array(sp.lst_equ_omega), bins=360, amp_range=(0.,180.), \
           title ='Equ: after indexing', xlabel='$\omega$ (deg)', ylabel='Events',\
           fnm='equ-omega.png')
    #---------
    h1d(np.array(sp.lst_equ_beta), bins=130, amp_range=(0.,260.), \
           title ='Equ: after indexing', xlabel='beta (deg)', ylabel='Events',\
           fnm='equ-beta.png')
    #---------
    h1d(np.array(sp.lst_equ_probsel), bins=100, amp_range=(0.,1.), \
           title ='Equ: after indexing', xlabel='$\propto$ Probability', ylabel='Events',\
           fnm='equ-probsel.png')
    #---------
    #---------
    fig, axhi, hi = h1d(np.array(sp.lst_equ_npk_idxed), bins=5, amp_range=(-0.5,4.5), \
           title ='Equ: after indexing', xlabel='Number of indexed peaks', ylabel='Events',\
           fnm='equ-nidxed-peaks.png')

    print 'Number of indexed peaks:', hi[0], ' for bins:', hi[1]
    #---------

#------------------------------
def plotResultsEquFit(prefix='plot', do_save=True, hwin_x0y0=(0,400)) :
    """Plots/saves histograms for selected events in EQU region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0
    #---------
    h1d(np.array(sp.lst_equ_ev_dp0), bins=80, amp_range=(-40,40), \
           title ='Equ: fit parameter phi-phi0', xlabel='d$\phi$ (deg)', ylabel='Events',\
           fnm='equ-ev-fit-dphi.png')
    #---------
    h1d(np.array(sp.lst_equ_ev_dp1), bins=120, amp_range=(-60,60), \
           title ='Equ: fit parameter beta-beta0', xlabel='d(beta) (deg)', ylabel='Events',\
           fnm='equ-ev-fit-dbeta.png')
    #---------
    h1d(np.array(sp.lst_equ_ev_p0), bins=80, amp_range=(-40,40), \
           title ='Equ: fit parameter p0=phi', xlabel='$\phi$ (deg)', ylabel='Events',\
           fnm='equ-ev-fit-phi.png')
    #---------
    h1d(np.array(sp.lst_equ_ev_p1), bins=120, amp_range=(-80,40), \
           title ='Equ: fit parameter p1=beta', xlabel='beta (deg)', ylabel='Events',\
           fnm='equ-ev-fit-beta.png')
    #---------
    h1d(np.array(sp.lst_equ_phi_err), bins=100, amp_range=(0,10), \
           title ='Equ: event fit angle phi error', xlabel='Error on $\phi$ (deg)', ylabel='Events',\
           fnm='equ-ev-fit-phi-err.png')
    #---------
    h1d(np.array(sp.lst_equ_beta_err), bins=100, amp_range=(0,50), \
           title ='Equ: event fit angle beta error', xlabel='Error on beta (deg)', ylabel='Events',\
           fnm='equ-ev-fit-beta-err.png')
    #---------
    h1d(np.array(sp.lst_equ_prob), bins=100, amp_range=(0,1), \
           title ='Equ: event fit probability', xlabel='Fit probability', ylabel='Events',\
           fnm='equ-ev-fit-prob.png')
    #---------

#------------------------------
def plotResultsArcFit(prefix='plot', do_save=True, hwin_x0y0 = (0,0)) :
    """Plots/saves histograms for phi-beta fit in ARC region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0

    #---------
    h1d(np.array(sp.lst_arc_ev_dp0), bins=80, amp_range=(-40,40), \
           title ='ARC: fit parameter phi-phi0', xlabel='d$\phi$ (deg)', ylabel='Events',\
           fnm='arc-ev-fit-dphi.png')
    #---------
    h1d(np.array(sp.lst_arc_ev_dp1), bins=100, amp_range=(-25,25), \
           title ='ARC: fit parameter beta-beta0', xlabel='d(beta) (deg)', ylabel='Events',\
           fnm='arc-ev-fit-dbeta.png')
    #---------
    h1d(np.array(sp.lst_arc_ev_p0), bins=80, amp_range=(-40,40), \
           title ='ARC: fit angle phi', xlabel='$\phi$ (deg)', ylabel='Events',\
           fnm='arc-ev-fit-phi.png')
    #---------
    h1d(np.array(sp.lst_arc_ev_p1), bins=100, amp_range=(-40,10), \
           title ='ARC: fit angle beta', xlabel='beta (deg)', ylabel='Events',\
           fnm='arc-ev-fit-beta.png')
    #---------
    h1d(np.array(sp.lst_arc_phi_err), bins=100, amp_range=(0,5), \
           title ='ARC: fit angle phi error', xlabel='Error on $\phi$ (deg)', ylabel='Events',\
           fnm='arc-ev-fit-phi-err.png')
    #---------
    h1d(np.array(sp.lst_arc_beta_err), bins=100, amp_range=(0,50), \
           title ='ARC: fit angle beta error', xlabel='Error on beta (deg)', ylabel='Events',\
           fnm='arc-ev-fit-beta-err.png')
    #---------
    h1d(np.array(sp.lst_arc_prob), bins=100, amp_range=(0,1), \
           title ='ARC: phi-beta fit probability', xlabel='Fit probability', ylabel='Events',\
           fnm='arc-ev-fit-prob.png')
    #---------
    #---------
    return
    #---------
    #---------
    h1d(np.array(sp.lst_arc_2peak_dist), bins=200, amp_range=(1,401), \
           figsize=(12,5), axwin=(0.08, 0.12, 0.88, 0.80), \
           title ='Arc: Distance between 2 peaks', xlabel='Distance[pixel]', ylabel='Events',\
           fnm='arc-2peak-dist.png')
    #---------
    h1d(np.array(sp.lst_arc_sel_2dist), bins=200, amp_range=(1,401), \
           figsize=(12,5), axwin=(0.08, 0.12, 0.88, 0.80), \
           title ='Arc: 2-peak distance sel for fit', xlabel='Distance[pixel]', ylabel='Events',\
           fnm='arc-2peak-dist-sel.png')
    #---------
    #h1d(np.array(sp.lst_arc_sel_2dcen), bins=120, amp_range=(380,440), \
    #       title ='Arc: 2-peak center radius', xlabel='Distance[pixel]', ylabel='Events',\
    #       fnm='arc-2peak-dcen-sel.png')
    #---------

#------------------------------
def init_psana_and_graphics() :
    """ Initialize psana and event imaging stuff
        Access to psana is needed in case if we need more data than it is available in peak-file.
    """
    if not sp.DO_PLOT : return
    if sp.psana_and_graphics_is_initialized : return

    #-----------
    import psana
    sp.psana = psana

    from Detector.AreaDetector import AreaDetector

    # Non-standard calib directory
    #sp.psana.setOption('psana.calib-dir', './empty/calib')
    #sp.psana.setOption('psana.calib-dir', '/reg/d/psdm/CXI/cxif5315/calib')
    #-----------

    dsname = 'exp=cxif5315:run=169:idx'
    src = psana.Source('DetInfo(CxiDs2.0:Cspad.0)')
    print '%s\nExample for\n  dataset: %s\n  source : %s' % (85*'_',dsname, src)

    ds = psana.DataSource(dsname)
    sp.env = ds.env()
    sp.runobj = ds.runs().next()
    #sp.evt = ds.events().next()
    #sp.runnum = sp.evt.run()

    sp.det = AreaDetector(src, sp.env, pbits=0, iface='P')
    #sp.det = psana.Detector(src, sp.env)

    sp.fig, sp.axim, sp.axcb = fig_axes() # if sp.DO_PLOT else (None, None, None)

    #sys.exit('init_psana_and_graphics is done')

    sp.psana_and_graphics_is_initialized = True

#------------------------------
def saveImageAveraged(prefix='plot', fnm='img-averaged.npy') :
    if sp.img_num :
        fname = '%s-%s' % (prefix, fnm)
        print 'Save averaged over %d events image in file %s' % (sp.img_num, fname)
        np.save(fname, sp.img_sum/sp.img_num)
        sp.img_num = 0

#------------------------------
def plotEvent(event) :

    if not sp.DO_PLOT : return

    sp.count_plt += 1
    if sp.count_plt <= sp.SKIP   : return
    if sp.count_plt > sp.EVTMAX :
        print 'Number of events exceeds requested %d' % sp.EVTMAX
        show() # hold image untill it is closed
        saveImageAveraged()
        sys.exit()

    init_psana_and_graphics()

    #tsec, tnsec, fid = 1434301977, 514786085, 44835
    #pk = sp.peak
    #pk = event.get_objs()[0]
    pk = event()[0] # get the 1st peak from the list for time stamp

    et = sp.psana.EventTime(int((pk.tsec<<32)|pk.tnsec), pk.fid)
    evt = sp.runobj.event(et)
    #for key in evt.keys() : print key
    runnum = evt.run()

    if sp.count_plt == 1 \
    or sp.count_plt == sp.SKIP+1 :

        sp.nda_peds  = sp.det.pedestals(runnum)
        sp.nda_bkgd  = sp.det.bkgd(runnum)
        sp.nda_smask = sp.det.mask(evt, calib=False, status=True, edges=True, central=True, unbond=True, unbondnbrs=True)

        print_ndarr(sp.nda_peds,  'nda_peds')
        print_ndarr(sp.nda_bkgd,  'nda_bkgd')
        print_ndarr(sp.nda_smask, 'nda_smask')

        # Pixel image indexes
        #sp.det.tilt_geo(evt,0,0,20)

        iX,iY = sp.det.indexes_xy(evt, xy0_off_pix=sp.xy0_off_pix)
        sp.iX = np.array(iX, dtype=np.int64)
        sp.iY = np.array(iY, dtype=np.int64)

        # Protect indexes (should be POSITIVE after offset subtraction)
        sp.imRow = np.select([sp.iX<sp.xoffset], [0], default=sp.iX-sp.xoffset)
        sp.imCol = np.select([sp.iY<sp.yoffset], [0], default=sp.iY-sp.yoffset)

        print_ndarr(sp.iX,    'iX')
        print_ndarr(sp.iY,    'iY')
        print_ndarr(sp.imRow, 'imRow')
        print_ndarr(sp.imCol, 'imCol')

        Xarr =  sp.det.coords_x(runnum)
        Yarr =  sp.det.coords_y(runnum)

        sp.rb = RadialBkgd(Xarr, Yarr, mask=sp.nda_smask, radedges=(5200, 80000), nradbins=200,\
                                                          phiedges=(40, 320), nphibins=1)
        sp.pf = polarization_factor(sp.rb.pixel_rad(), sp.rb.pixel_phi()+90, sp.DETD_um)

        print 80*'_'
    #sys.exit('plotEvent is done')

    nda_raw = sp.det.raw(evt)

    if nda_raw is not None :

        nda =  np.array(nda_raw, dtype=np.float32, copy=True)
        nda -= sp.nda_peds

        #sp.det.common_mode_apply(runnum, nda, cmpars=(5,50))

        # Subtract background shape averaged for pure water
        nda = subtract_bkgd(nda, sp.nda_bkgd, mask=sp.nda_smask, winds=sp.winds_bkgd, pbits=0)

        # Subtract dynamically evaluater redial background
        #nda = sp.rb.subtract_bkgd(nda.flatten() * sp.pf)
        #nda.shape = sp.nda_peds.shape

        nda *= sp.nda_smask
        
        img = getImage(evt, nda)

        ave, rms = img.mean(), img.std()
        amin, amax = ave-1*rms, ave+8*rms
        plot_img(img, mode='do not hold', amin=amin, amax=amax)
        #plot_peaks_on_img(peaks_arc, sp.axim, sp.imRow, sp.imCol, color='w') #, pbits=3)

        # drow peaks on the top of image. For now it is not implemented for Fraser transformed image.
        if not sp.DO_FRASER :        
            lst_peaks = np.array(sp.lst_plt_evt_peaks)
            plot_peaks_on_img(lst_peaks, sp.axim, sp.imRow, sp.imCol, color='w') #, pbits=3)

        sp.fig.canvas.set_window_title('Event: %d selected: %d ' % (pk.evnum, sp.count_plt)) 
        sp.fig.canvas.draw() # re-draw figure content

        if sp.img_sum is None :
            sp.img_sum  = np.array(img, dtype=np.float64)
        else :
            sp.img_sum += img
        sp.img_num += 1

#------------------------------
def getImage(evt, nda) :
    """Returns image to plot 
    """

    if sp.DO_PHI_ROTATION :

        #apply tilt angle to current geometry object
        sp.det.tilt_geo(evt,0,0,sp.fit_phi - sp.fit_phi0)
        sp.fit_phi0 = sp.fit_phi

        #get indexes for tilted geometry
        iX,iY = sp.det.indexes_xy(evt, xy0_off_pix=sp.xy0_off_pix, do_update=True)
        sp.iX = np.array(iX, dtype=np.int64)
        sp.iY = np.array(iY, dtype=np.int64)

        # Protect indexes (should be POSITIVE after offset subtraction)
        sp.imRow = np.select([sp.iX<sp.xoffset], [0], default=sp.iX-sp.xoffset)
        sp.imCol = np.select([sp.iY<sp.yoffset], [0], default=sp.iY-sp.yoffset)

    #img = sp.det.image(evt, nda)
    img = sp.det.image(evt, nda, xy0_off_pix=sp.xy0_off_pix)[sp.xoffset:sp.xoffset+sp.xsize, sp.yoffset:sp.yoffset+sp.ysize]


    if sp.DO_FRASER :
        s12, s3, recimg = fraser(img, sp.fit_beta, sp.DETD, center=sp.img_center, oshape=sp.fras_shape)
        return recimg

    return img

#------------------------------
def funcy_equ(x, phi_deg, bet_deg) :
    """Wrapper for funcy_l0
    """    
    return funcy_l0(x, phi_deg, bet_deg)
    #return funcy(x, phi_deg, bet_deg) 

#------------------------------
def funcy_arc(x, phi_deg, bet_deg) :
    """Wrapper for funcy_l1 with setting of parameter DoR and sgnrt
    """    
    #return funcy_l1_v0(x, phi_deg, bet_deg, DoR=433/sp.DETD, sgnrt=-1.) # 433/sp.DETD
    return funcy_l1_v1(x, phi_deg, bet_deg, DoR=390/sp.DETD, sgnrt=1.) # 391/sp.DETD - GOOD SOLUTION
    #return funcy_l1_v1(x, phi_deg, bet_deg, DoR=391/sp.DETD, sgnrt=-1) # 391/sp.DETD
 
#------------------------------
def fitArcPhiBeta(e) :
    """ Fit 2-peak event in ARC region
    """
    if len(sp.lst_arc_evt_peaks) != 2 : return
    #if sp.d_2peak < 125 : return # v1
    if sp.d_2peak < 220 : return # v2
    if sp.d_2peak > 260 : return # v2

    L = sp.DETD_um

    # Assume regular Cortesian orientation for x,y 
    xn = np.array([pk.x/L for pk in sp.lst_arc_evt_peaks], dtype=np.double)
    yn = np.array([pk.y/L for pk in sp.lst_arc_evt_peaks], dtype=np.double)
    en = np.array([pk.csigma*sp.PIXEL_SIZE/L for pk in sp.lst_arc_evt_peaks], dtype=np.double)
    #en = np.array([max(pk.rsigma,pk.csigma)*sp.PIXEL_SIZE/L for pk in sp.lst_arc_evt_peaks], dtype=np.double)
    #en *= 0.5

    print 'Fit ARC beta-phi xn: ', xn
    print 'Fit ARC beta-phi yn: ', yn
    print 'Fit ARC beta-phi en: ', en

    #p0 = [-3.2,-9] # phi, beta angle central values
    p0 = [-3.2,-17] # phi, beta angle central values

    #if True :
    try :
        p01, pcov1 = curve_fit(funcy_arc, xn, yn, p0)
        print 'Fit ARC beta-phi: fit results popt1: %s' % (str(p01))

        popt, pcov = curve_fit(funcy_arc, xn, yn, p01, en, absolute_sigma=True) # absolute_sigma is missing in old vers
        #popt, pcov = curve_fit(funcy_arc, xn, yn, p01, en)
        print 'Fit ARC beta-phi: fit results popt : %s' % (str(popt))
        ###print 'Fit ARC beta-phi: pcov:\n', str(pcov)

        sp.fia_phi,  sp.fia_beta  = popt[0:2]
        sp.fia_dphi, sp.fia_dbeta = popt[0:2]-p0[0:2]

        dy = yn - funcy_arc(xn, *popt)
        print 'Fit ARC beta-phi dy: ', dy, '  <===  yn - f(xn)'

        sp.fia_phi_err, sp.fia_beta_err = np.sqrt(np.diag(pcov))
        print 'Fit ARC beta-phi: phi_err=%.3f  beta_err=%.3f' % (sp.fia_phi_err, sp.fia_beta_err)

        sp.fia_chi2 = np.sum(((funcy_arc(xn, *popt) - yn) / en)**2)


        sp.fia_ndof = len(xn) - 1
        sp.fia_prob = chi2.sf(sp.fia_chi2, sp.fia_ndof, loc=0, scale=1)
        print 'Fit ARC beta-phi: chi2=%.2f  npks=%d  prob=%.6f' % (sp.fia_chi2, sp.fia_ndof, sp.fia_prob)

        #if eventArcFittedIsSelected():

        sp.lst_arc_sel_2dist.append(sp.d_2peak)
        sp.lst_arc_sel_2dcen.append(sp.dc_2peak)

        sp.lst_arc_ev_p0.append(sp.fia_phi)
        if True : # sp.fia_beta_err > 0.1 and sp.fia_beta_err < 100 :
            sp.lst_arc_ev_p1.append(sp.fia_beta)
        sp.lst_arc_ev_dp0.append(sp.fia_dphi)
        sp.lst_arc_ev_dp1.append(sp.fia_dbeta)

        sp.lst_arc_phi_err.append(sp.fia_phi_err)
        sp.lst_arc_beta_err.append(sp.fia_beta_err)
        sp.lst_arc_prob.append(sp.fia_prob)

        sp.fit_arc_status = True

    #except OptimizeWarning:
    except (RuntimeError), e:
        print 'RuntimeError in curve_fit in %s' % sys._getframe().f_code.co_name
        print 'reason:', e

#------------------------------
def procEventArc(e) :
    """ Process event for ARC region
    """
    # summary for event
    sp.lst_arc_npksel.append(sp.count_arc_pks_sel)
    sp.lst_arc_npktot.append(sp.count_arc_pks_tot)

    sp.d_2peak       = None
    sp.dc_2peak      = None
    sp.phi_2peak_ori = None

    if sp.count_arc_pks_sel == 1 :
        #r1, phi1 = sp.lst_arc_evt_peaks[0][6:8]
        pk1 = sp.lst_arc_evt_peaks[0]
        r1, phi1 = pk1.r, pk1.phi

        sp.lst_arc_1peak_r  .append(r1)
        sp.lst_arc_1peak_phi.append(phi1)

    elif sp.count_arc_pks_sel == 2 :
    # process 2-peak event
        sp.count_arc_2peaks += 1
        pk1, pk2 = sp.lst_arc_evt_peaks[0:2]
        x1,y1,r1,phi1 = pk1.x, pk1.y, pk1.r, pk1.phi 
        x2,y2,r2,phi2 = pk2.x, pk2.y, pk2.r, pk2.phi  
        dx, dy = x2-x1, y2-y1
        d = math.sqrt(dx*dx + dy*dy) / sp.PIXEL_SIZE

        xc, yc = (x2+x1)/2, (y2+y1)/2
        dc = math.sqrt(xc*xc + yc*yc) / sp.PIXEL_SIZE
        #print '  Event %6d, 2 peak coordinates:  %6d %6d %6d %6d d=%8.1f[pix]' % (e()[0].evnum, x1, y1, x2, y2, d)

        sp.d_2peak       = d
        sp.dc_2peak      = dc
        sp.phi_2peak_ori = 0.5*(phi1+phi2)

        sp.lst_arc_2peak_dist.append(d)
        sp.lst_arc_2peak_dr  .append(r2-r1)
        sp.lst_arc_2peak_phi .append(sp.phi_2peak_ori)

        #phi-beta fit to two peaks in ARC region
        fitArcPhiBeta(e)

#------------------------------
def procPeakDataArc(pk) :
    """ Process peak for ARC region; accumulate peak statistics in histogram arrays.
    """
    #===================
    # discard from all histograms except its own
    sp.lst_arc_atot.append(pk.atot)
    if pk.atot<1600 : return
    #===================
    sp.lst_arc_amax.append(pk.amax)
    sp.lst_arc_npix.append(pk.npix)
    sp.lst_arc_r   .append(pk.r)
    sp.lst_arc_phi .append(pk.phi)
    sp.lst_arc_bkgd.append(pk.bkgd)
    sp.lst_arc_rms .append(pk.rms)
    sp.lst_arc_son .append(pk.son)
    sp.lst_arc_sonc.append(pk.sonc)

    sp.count_arc_pks_tot += 1
    if not peakIsSelectedArc(pk) : return
    sp.count_arc_pks_sel += 1

    sp.lst_arc_evt_peaks.append(pk)
    if sp.DO_PLOT : sp.lst_plt_evt_peaks.append((pk.seg, pk.row, pk.col, pk.amax, pk.atot, pk.npix))

#------------------------------
def procPeakDataEqu(pk) :
    """ Process peak for EQU region; accumulate peak data
    """
    # discard small radius peaks from all histograms
    #===================
    sp.lst_equ_atot.append(pk.atot)
    if pk.atot<1600 : return
    #sp.lst_equ_r_raw.append(pk.r)
    if pk.r<100 : return
    #===================

    sp.lst_equ_amax.append(pk.amax)
    sp.lst_equ_npix.append(pk.npix)
    sp.lst_equ_r   .append(pk.r)
    sp.lst_equ_phi .append(pk.phi)
    sp.lst_equ_bkgd.append(pk.bkgd)
    sp.lst_equ_rms .append(pk.rms)
    sp.lst_equ_son .append(pk.son)
    sp.lst_equ_sonc.append(pk.sonc)

    sp.count_equ_pks_tot += 1
    if not peakIsSelectedEqu(pk) : return
    sp.count_equ_pks_sel += 1

    pk.dphicmb = pk.dphi000-1.5 if -90<pk.phi and pk.phi<90 else pk.dphi180-4.8

    sp.lst_equ_dphicmb.append(pk.dphicmb)
    sp.lst_equ_dphi000.append(pk.dphi000)
    sp.lst_equ_dphi180.append(pk.dphi180)

    sp.lst_equ_evt_peaks.append(pk)
    if sp.DO_PLOT : sp.lst_plt_evt_peaks.append((pk.seg, pk.row, pk.col, pk.amax, pk.atot, pk.npix))

#------------------------------
def peakIsSelectedArc(pk) :
    """Apply peak selection criteria to each peak from file
    """
    if sp.USE_PEAK_SELECTION :
        #if pk.son<9     : return False
        #if pk.amax<150  : return False
        if pk.atot<1600 : return False
        if pk.npix>200  : return False
        if pk.r<435     : return False
        if pk.r>443     : return False
        if pk.rms<10    : return False
        if pk.rms>60    : return False
        if math.fabs(pk.bkgd)>20: return False

    return True

#------------------------------
def peakIsSelectedEqu(pk) :
    """Apply peak selection criteria to each peak from file
    """
    if sp.USE_PEAK_SELECTION :
        #if pk.son<9     : return False
        #if pk.amax<150  : return False
        if pk.atot<1600  : return False
        if pk.npix>200  : return False
        if pk.rms>60    : return False
        if math.fabs(pk.bkgd)>20 : return False

    if pk.r<100     : return False
    if pk.r>454     : return False

    return True

#------------------------------
def fitEquPhiBeta(e) :
    """ Fit event peaks in EQU region
    """
    if not sp.event_is_selected : return
    if sp.count_equ_pks_sel < 2 : return

    L = sp.DETD_um

    # Assume regular Cortesian orientation for x,y 
    xn = np.array([pk.x/L for pk in sp.lst_equ_evt_peaks], dtype=np.double)
    yn = np.array([pk.y/L for pk in sp.lst_equ_evt_peaks], dtype=np.double)
    en = np.array([pk.csigma*sp.PIXEL_SIZE/L for pk in sp.lst_equ_evt_peaks], dtype=np.double)
    #en = np.array([max(pk.rsigma,pk.csigma)*sp.PIXEL_SIZE/L for pk in sp.lst_equ_evt_peaks], dtype=np.double)
    #en = np.ones_like(xn)

    #print 'Fit beta-phi xn: ', xn
    #print 'Fit beta-phi yn: ', yn
    #print 'Fit beta-phi en: ', en

    p0 = [-3.2,-16.0] # phi, beta angle central values

    #if True :
    try :
        popt, pcov = curve_fit(funcy_equ, xn, yn, p0, en, absolute_sigma=True)
        #print 'Fit beta-phi: fit results popt: ', str(popt)
        #print 'Fit beta-phi: pcov:\n', str(pcov)

        sp.fib_phi,  sp.fib_beta  = popt[0:2]
        sp.fib_dphi, sp.fib_dbeta = popt[0:2]-p0[0:2]

        sp.fib_phi_err, sp.fib_beta_err = np.sqrt(np.diag(pcov))
        #print 'Fit beta-phi: phi_err=%.3f  beta_err=%.3f' % (sp.fib_phi_err, sp.fib_beta_err)

        sp.fib_chi2 = np.sum(((funcy_equ(xn, *popt) - yn) / en)**2)
        sp.fib_ndof = len(xn) - 1
        sp.fib_prob = chi2.sf(sp.fib_chi2, sp.fib_ndof, loc=0, scale=1)
        #print 'Fit beta-phi: chi2=%.2f  npks=%d  prob=%.6f' % (sp.fib_chi2, sp.fib_ndof, sp.fib_prob)

        #if eventEquFittedIsSelected():
        sp.lst_equ_ev_p0.append(sp.fib_phi)
        sp.lst_equ_ev_p1.append(sp.fib_beta)
        sp.lst_equ_ev_dp0.append(sp.fib_dphi)
        sp.lst_equ_ev_dp1.append(sp.fib_dbeta)

        sp.lst_equ_phi_err.append(sp.fib_phi_err)
        sp.lst_equ_beta_err.append(sp.fib_beta_err)
        sp.lst_equ_prob.append(sp.fib_prob)

        sp.fit_equ_status = True

    except :
        print 'RuntimeError in curve_fit in %s' % sys._getframe().f_code.co_name
        #print 'RuntimeError in curve_fit'

#------------------------------
def procEventEqu(e) :
    """Process event for EQU region
    """
    # summary for event
    sp.lst_equ_npksel.append(sp.count_equ_pks_sel)
    sp.lst_equ_npktot.append(sp.count_equ_pks_tot)
    if sp.count_arc_pks_sel == 2 :
        sp.lst_equ_npks_2inarc.append(sp.count_equ_pks_sel)

    L = sp.DETD_um
    orient = sp.phi_2peak_ori # evaluated from 2-peak event in Arc region

    # Histogram for orientation from 2-peak ARC 
    #===========================================
    if  sp.event_is_selected \
        and (orient is not None) \
        and abs(orient) < 20 :

        for peak in sp.lst_equ_evt_peaks :
            r, phi = peak.r, peak.phi
            if(phi<0) :
                sp.lst_equ_r1.append(r)
                sp.lst_equ_phi1.append(phi)
                sp.lst_equ_phi1c.append(phi-orient)
            else :
                sp.lst_equ_r2.append(r)
                sp.lst_equ_phi2.append(phi)
                sp.lst_equ_phi2c.append(phi-orient)

    # Evaluated angles for 2-peak in EQU
    #====================================
    if  sp.event_is_selected \
        and sp.count_equ_pks_sel == 2 :

        sp.count_equ_2peaks += 1
        pk1,pk2 = sp.lst_equ_evt_peaks[0:2]
        evnum1, amax1, atot1, npix1, x1, y1, r1, phi1 = pk1.evnum, pk1.amax, pk1.atot, pk1.npix, pk1.x, pk1.y, pk1.r, pk1.phi 
        evnum2, amax2, atot2, npix2, x2, y2, r2, phi2 = pk2.evnum, pk2.amax, pk2.atot, pk2.npix, pk2.x, pk2.y, pk2.r, pk2.phi 
        phi  = calc_phi (x1, y1, x2, y2, L) # sp.DETD)
        beta = calc_beta(x1, y1, phi, L) # sp.DETD)

        phi_deg  = math.degrees(phi)
        beta_deg = math.degrees(beta)

        #print 'Found 2-peak event in equ, r1=%f, r2=%f, phi1=%f, phi2=%f, Event phi=%f, beta=%f' % (r1, r2, phi1, phi2, phi_deg, beta_deg)

        #if math.fabs(phi_deg) > 70 :
        sp.lst_equ_ev_phi .append(phi_deg)
        sp.lst_equ_ev_beta.append(beta_deg)

    # Fit angles for > 1,2,etc.-peak in EQU
    #======================================
    fitEquPhiBeta(e)

#------------------------------
def savePeaksWithFitResults(evgrp) :
    """ save file with peak + fit info in records
    """
    lst_peaks = sp.lst_arc_evt_peaks + sp.lst_equ_evt_peaks
    
    rec0 = evgrp.get_objs()[0]
    pk0  = lst_peaks[0]

    # create output store object
    if  sp.pstore_fit is None :
        header_fit = ' results for exp=%s:run=%d' % (pk0.exp, pk0.run)
        sp.pstore_fit = PeakStore(pk0.exp, pk0.run, sp.ofn_prefix_fit, header_fit, add_header='', pbits=1)

    # add event header
    print 40*'=', ' save in file:'

    #rec_hdr = '# %s  phi-fit beta-fit   qh-fit    qv-fit  dqh[1/A]' % (sp.fc.hdr)
    rec_hdr = '# %s  fit-phi fit-beta  phi-err  beta-err  fit-chi2  ndof  fit-prob' % (sp.fc.hdr[1:])
    print rec_hdr
    sp.pstore_fit.save_comment(rec_hdr)

    fmt = '%s  %7.2f %7.2f %9.6f %9.6f %9.6f %5d %9.6f'

    # add records about non-matched peaks
    for pk in lst_peaks :
        msg = fmt % (pk.line[:-1], sp.fit_phi, sp.fit_beta, sp.fit_phi_err, sp.fit_beta_err, sp.fit_chi2, sp.fit_ndof, sp.fit_prob)
        print msg
        sp.pstore_fit.save_peak(peak_rec=msg)

    # add space between events for readability
    sp.pstore_fit.save_peak(peak_rec='')

#------------------------------
def lut_qh_hist_array() :
    """Loops over lookup table, makes array of (qh,P) tuples for histogram comparison 
    """
    fc = sp.fc_lut
    lst_lut_qh_prob=[]
    # loop over look-up table of crystal orientations
    for orinum in fc.group_num_iterator() :
        origrp = fc.next()
        for rec in origrp() :
            lst_lut_qh_prob.append((rec.qh, rec.P))

    return np.array(lst_lut_qh_prob, dtype=np.float32)

#------------------------------
def findEventIndexes(e) :
    """ 1. apply phi and beta(fraser) transformations
           to (sample-to-detector distance) normalized peak coordinate arrays xn, yn.
        2. make qv histogram
        3. discard events with |qv| > sp.QVMAX
        4. make qh histogram
        5. find most probable orientation using look-up table
        6. fill istograms
        7. save peaks with index info in file
    """
    print 'procEventIndexing: Begin indexing with phi=%.1f  beta=%.1f' % (sp.fit_phi, sp.fit_beta)

    # apply phi and beta(fraser) transformations
    xn = [pk.x for pk in sp.lst_equ_evt_peaks]
    yn = [pk.y for pk in sp.lst_equ_evt_peaks]
    beta = 0 # sp.fib_beta
    qh_arr, qv_arr = rotation_phi_beta(xn, yn, sp.DETD_um, sp.fit_phi, beta, sp.SCALE) 

    fc = sp.fc_lut

    # fill qv histogram list, save evaluated qv, qh in peak
    for qh, qv, pk in zip(qh_arr, qv_arr, sp.lst_equ_evt_peaks) :            
        sp.lst_equ_qv_fit_selec.append(qv)
        pk.qh, pk.qv = qh, qv        
        #print 'peak qh=%8.4f   qv=%8.4f' % (qh, qv)

    # discard events with |qv| > sp.QVMAX
    for qv in qv_arr :
        if math.fabs(qv) > sp.QVMAX : return

    # fill qh histogram list
    for qh in qh_arr :            
        sp.lst_equ_qh_fit_selec.append(qh)

    npk_idxed_max = 0
    origrpsel     = None
    probgrpsel    = 0
    sp.probgrp    = 0

    # loop over look-up table of crystal orientation groups
    for orinum in fc.group_num_iterator() :
        origrp = fc.next()
        npk_idxed = 0


        #if orinum>100 : break

        # loop over records in crystal orientation group
        for rec in origrp() :

            rec.match_peak = None
            rec.match_qh   = None
            rec.match_qv   = None
            rec.match_dqh  = None

            # loop over event peaks and try to match them with group record
            for qh, qv, pk in zip(qh_arr, qv_arr, sp.lst_equ_evt_peaks) :
                dqh = math.fabs(qh - rec.qh)
                if dqh > sp.DQHMAX : continue
                prob = rec.P * math.exp(sp.FACTOR * dqh * dqh)
                #print 'rec beta=%6.1f  omega=%6.1f  qh=%7.4f  |dqh|=%8.6f  prob=%8.6f' % (rec.beta, rec.omega, qh, dqh, prob)

                if prob > sp.PROBTHR :
                    npk_idxed += 1
                    sp.probgrp = prob if npk_idxed == 1 else sp.probgrp * prob 

                    # add in oreintation record info about matched peak
                    rec.match_peak  = pk
                    rec.match_qh    = qh
                    rec.match_qv    = qv
                    rec.match_dqh   = dqh

        if  npk_idxed  > npk_idxed_max\
        or (npk_idxed == npk_idxed_max and sp.probgrp > probgrpsel) :
            npk_idxed_max = npk_idxed
            origrpsel  = origrp
            probgrpsel = sp.probgrp
            origrpsel.match_prob = probgrpsel
            #print 'rec beta=%6.1f  omega=%6.1f' % (rec.beta, rec.omega)

    if origrpsel is None : return

    sp.lst_equ_npk_idxed.append(npk_idxed_max)  
    if npk_idxed_max < 2 : return # select >=2 indexed peaks
    #if npk_idxed_max < 1 : return

    sp.count_evt_idxed += 1

    # for found orientation fill histograms, save peaks in file
    rec = origrpsel()[0]

    print 'Most probable orientation: beta=%6.1f  omega=%6.1f  prob=%8.6f' % (rec.beta, rec.omega, probgrpsel)
    # fill histogram lists
    sp.lst_equ_omega.append(rec.omega)    
    sp.lst_equ_beta.append(rec.beta)    
    sp.lst_equ_probsel.append(probgrpsel)  

    saveEventIndexingResults(origrpsel) 

#------------------------------
def saveEventIndexingResults(origrp) : # and sp.lst_equ_evt_peaks
    """ save peaks with index info in file
    """
    # save in peak object reference to the matched orientation record 
    for pk in sp.lst_equ_evt_peaks : pk.match_orirec = None
    for rec in origrp() :
        pk = rec.match_peak
        if pk is not None :
            pk.match_orirec = rec

    # save event indexing results in file
    # -----------------------------------
    rec0 = origrp.get_objs()[0]
    pk0  = sp.lst_equ_evt_peaks[0]
    z = 0

    # create output store object
    if  sp.pstore_idx is None :
        header_idx = 'Indexing results for exp=%s:run=%d' % (pk0.exp, pk0.run)
        sp.pstore_idx = PeakStore(pk0.exp, pk0.run, sp.ofn_prefix_idx, header_idx, add_header='', pbits=1)

    # add event header
    print 40*'=', ' save in file:'
    #evt_hdr = '#  %s' % sp.fc.hdr[:73]
    #evt_rec = '# %s' % pk0.line[:73]
    #print evt_hdr
    #print evt_rec
    #sp.pstore_idx.save_comment(evt_hdr)
    #sp.pstore_idx.save_comment(evt_rec)

    rec_hdr = '# STATE phi-fit beta-fit   qh-fit    qv-fit  dqh[1/A]  %s  %s  %s' % (sp.fc.hdr[:12], sp.fc.hdr[35:], sp.fc_lut.hdr)
    print rec_hdr
    sp.pstore_idx.save_comment(rec_hdr)

    fmt = '%7s %7.2f %7.2f %9.6f %9.6f %9.6f %s %s %7d  %s %s'

    # add records about non-matched lattice nodes
    for rec in origrp() :
        if rec.match_peak is None :
            msg = fmt % ('NODE-NM', sp.fit_phi, sp.fib_beta, z,z,z, pk0.line[:13], pk0.line[35:65], pk0.evnum, pk0.empty[75:], rec.line[:-1])
            print msg
            sp.pstore_idx.save_peak(peak_rec=msg)
    
    # add records about non-matched peaks
    for pk in sp.lst_equ_evt_peaks :
        if pk.match_orirec is None :
            msg = fmt % ('PEAK-NM', sp.fit_phi, sp.fib_beta, pk.qh, pk.qv, z, pk.line[:13], pk.line[35:65], pk.evnum, pk.line[75:-1], rec0.empty)
            print msg
            sp.pstore_idx.save_peak(peak_rec=msg)

    # add records about matched peaks with lattice nodes 
    for rec in origrp() :
        pk = rec.match_peak
        if pk is not None :
            msg = fmt %\
                  ('MATCHED', sp.fit_phi, sp.fib_beta, rec.match_qh, rec.match_qv, rec.match_dqh, pk.line[:13], pk.line[35:65], pk.evnum, pk.line[75:-1], rec.line[:-1])
            print msg
            sp.pstore_idx.save_peak(peak_rec=msg)

    # add space between events for readability
    sp.pstore_idx.save_peak(peak_rec='')

#------------------------------
def printEventEquPeaks() :
    """ print equatorial peaks
    """
    if sp.lst_equ_evt_peaks == [] : return
    print ' %s\n%s\n %s' % (sp.fc.hdr[:72], sp.lst_equ_evt_peaks[0].line[:73], sp.fc.hdr[73:])
    for peak in sp.lst_equ_evt_peaks :
        print peak.line[73:-1]

#------------------------------
def initEvent() :
    """ Initialization for the event processing
    """
    # Arc region
    sp.count_arc_pks_tot = 0
    sp.count_arc_pks_sel = 0
    sp.lst_arc_evt_peaks = []
    sp.fit_arc_status    = False

    # Equ region
    sp.count_equ_pks_tot = 0
    sp.count_equ_pks_sel = 0
    sp.lst_equ_evt_peaks = []
    sp.fit_equ_status    = False

    sp.lst_plt_evt_peaks = []

#------------------------------
def eventIsSelected() :
    """Apply selection criteria to entire event, based on list of peaks 
    """
    sp.event_is_selected = False

    if sp.count_arc_pks_sel > 2 : return False
    if sp.count_equ_pks_sel > 5 : return False
    #if sp.count_equ_pks_sel  < 2 : return False

    # Require all peaks in the narrow range of dphi[deg] 
    #for pk in sp.lst_equ_evt_peaks :
    #    if math.fabs(pk.dphicmb) > 6 : return False

    sp.event_is_selected = True
    return True

#------------------------------
def eventFittedIsSelected() :
    """Selection is true if fitted angles phi and beta are within [-3*sigma,+3*sigma] of their central value
    """
    if sp.BETA_FROM_ARC :
        if not sp.fit_arc_status : return False 
        if math.fabs(sp.fia_dphi)  > 5 : return False # remove tails
        if math.fabs(sp.fia_dbeta) > 10 : return False # remove tails

        sp.fit_phi      =  sp.fia_phi
        sp.fit_beta     = -sp.fia_beta
        sp.fit_beta_err =  sp.fia_beta_err
        sp.fit_phi_err  =  sp.fia_phi_err
        sp.fit_chi2, sp.fit_ndof, sp.fit_prob = sp.fia_chi2, sp.fia_ndof, sp.fia_prob

    else :
        if not sp.fit_equ_status : return False 
        if math.fabs(sp.fib_dphi) > 5 : return False # remove tails
        if math.fabs(sp.fib_dbeta) > 10 : return False # remove tails

        sp.fit_phi      =  sp.fib_phi
        sp.fit_beta     = -sp.fib_beta
        sp.fit_beta_err =  sp.fib_beta_err
        sp.fit_phi_err  =  sp.fib_phi_err
        sp.fit_chi2, sp.fit_ndof, sp.fit_prob = sp.fib_chi2, sp.fib_ndof, sp.fib_prob

    if sp.fit_beta_err > 15 : return False
    if sp.fit_phi_err  > 3  : return False

    return True
 
#------------------------------
def procEvent(e) :
    """Process event when all its peaks are loaded.
    """
    if eventIsSelected() : sp.count_evt_sel += 1

    procEventArc(e) 
    procEventEqu(e)

    if sp.event_is_selected and eventFittedIsSelected() :

        savePeaksWithFitResults(e)

        plotEvent(e)
        for pk in sp.lst_equ_evt_peaks : sp.lst_equ_r_fit_selec.append(pk.r) # fill hist list

        printEventEquPeaks()
        findEventIndexes(e)

#------------------------------

def procPeaksFromFile(fname, fname_lut, ofn_prefix='plot') :
    """Main method for processing of peaks from files;
       1. load peak data from file in TDFileContainer object,
       2. initialize analysis parameters, histogram arrays,
       3. loop over file events/peaks, get and process peak data,
       4. plot and show results.
    """
    from algos.core.TDFileContainer import TDFileContainer
    from algos.core.TDNodeRecord    import TDNodeRecord
    from algos.core.TDPeakRecord    import TDPeakRecord

    sp.ofn_prefix_idx = os.path.join(os.path.dirname(ofn_prefix),'peak-idx')
    sp.ofn_prefix_fit = os.path.join(os.path.dirname(ofn_prefix),'peak-fit')

    #--- load look-up table for crystal orientation indexing
    sp.fc_lut = TDFileContainer(fname_lut, indhdr='index', objtype=TDNodeRecord) #, pbits=256)
    sp.fc_lut.print_content(nlines=50)

    #--- load table of peaks
    t0_sec = time()
    sp.fc = TDFileContainer(fname, indhdr='Evnum', objtype=TDPeakRecord) #, pbits=256)
    sp.fc.print_content(nlines=20)

    #----------
    # sys.exit('Test exit') ######### TEST EXIT
    #----------

    for evnum in sp.fc.group_num_iterator() :
        event = sp.fc.next()
        lst_peaks = event.get_objs()

        initEvent()

        print '%s Event# %6d %s' % (4*'_', evnum, 4*'_')
        #print '%s\n %s\n%s\n%s' % (71*'_', sp.fc.hdr[:70], lst_peaks[0].line[:71], sp.fc.hdr[72:])
        for peak in lst_peaks :
            #print peak.line.rstrip('\n')[73:]
            if peak.reg == 'ARC' : procPeakDataArc(peak) ############
            if peak.reg == 'EQU' : procPeakDataEqu(peak) ############

        procEvent(event)                                 ############

    #----------
    #    if evnum > 2000 : break # sys.exit('TEST EXIT') # break ######### TEST EXIT
    #    if evnum > 575 : break ######### TEST EXIT
    #----------

    print 'Consumed time = %7.3f sec' % (time()-t0_sec)
    print 'Number of selected events: %d or %7.3f of total number of events: %d' % \
          (sp.count_evt_sel, float(sp.count_evt_sel)/evnum, evnum)
    print 'Number of 2-peak in arc events: %d, of total number of processed events: %d' % \
          (sp.count_arc_2peaks, evnum)
    print 'Number of 2-peak in equ events: %d, of total number of processed events: %d' % \
          (sp.count_equ_2peaks, evnum)
    print 'Number of indexed events: %d' % (sp.count_evt_idxed)

    if sp.PLOT_HIST :
        plotResultsArc   (ofn_prefix, do_save=True, hwin_x0y0=(0,0))
        plotResultsEqu   (ofn_prefix, do_save=True, hwin_x0y0=(0,600))
        plotResultsEquSel(ofn_prefix, do_save=True, hwin_x0y0=(0,300)) ############
        plotResultsEquFit(ofn_prefix, do_save=True, hwin_x0y0=(0,400))
        plotResultsArcFit(ofn_prefix, do_save=True, hwin_x0y0=(0,100))
        show()
        saveImageAveraged(ofn_prefix, fnm='img-averaged.npy')          ############

#------------------------------
def do_main() :
    """ Main method to do work
    """
    # Extend a permitted number of windows for histograms
    import os
    import matplotlib as mpl
    from algos.core.GlobalUtils import create_directory

    #print mpl.rcParams.keys()
    #mpl.rcParams['figure.max_open_warning'] = 60 

    wdir = '/reg/neh/home/dubrovin/LCLS/rel-mengning/work'
    rdir = './results'
    create_directory(wdir)
    create_directory(rdir)

    fname_lut = '%s/lut-cxif5315-r0169-2016-02-03T15:10:48.txt' % wdir

    # Data processing since 2016-03-30
    #fname_pft = '%s/pfv2-cxif5315-r0169-2016-04-01.txt' % wdir
    #fname_pft = '%s/pfv3-cxif5315-r0169-2016-04-01.txt' % wdir
    #fname_pft = '%s/pfv4-cxif5315-r0169-2016-04-01.txt' % wdir
    #fname_pft = '%s/2016-04-12-meng-r169_filter1.txt' % wdir
    #fname_pft = '%s/2016-04-20-meng-r169_filter1.txt' % wdir
    #pfvers = 's-bkgd'

    # Data processing since 2016-04-20
    #pfvers = 'xpfv2r1'
    pfvers = 'xpfv3r1'
    #pfvers = 'xpfv4r1'
    #fname_pft = '%s/%s-cxif5315-r0169-2016-04-25.txt' % (wdir, pfvers)
    #fname_pft = '%s/%s-cxif5315-r0169-2016-04-22.txt' % (wdir, pfvers)
    #fname_pft = '%s/%s-cxif5315-r0169-2016-05-03.txt' % (wdir, pfvers) # bkgd w/o cm
    #fname_pft = '%s/%s-cxif5315-r0169-2016-05-05.txt' % (wdir, pfvers) # bkgd w/o loc max
    #fname_pft = '%s/%s-cxif5315-r0169-2016-05-09.txt' % (wdir, pfvers) # bkgd w/o loc max
    fname_pft = '%s/%s-cxif5315-r0169-2016-05-10.txt' % (wdir, pfvers) # cmod #5

    oprefix   = '%s/plot-cxif5315-r0169-%s-2016-05-10-v01' % (rdir, pfvers)

    procPeaksFromFile(fname_pft, fname_lut, oprefix)

#------------------------------
if __name__ == "__main__" :
    do_main()
    sys.exit('Processing is completed')

#------------------------------
