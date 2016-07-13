#!/usr/bin/env python

#------------------------------

import sys
import math
import numpy as np
from time import time
from pyimgalgos.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, fig_axes, plot_img, plot_peaks_on_img

#------------------------------

class Storage :
    """Store for shared parameters"""
    def __init__(self) :
       self.lst_npeaks_ch = []
       self.lst_npeaks_pf = []
       self.lst_matched   = []

#------------------------------

sp = Storage() # make a singleton for Storage

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

def plot_retults(prefix='plot', do_save=True, hwin_x0y0 = (0,0)) :
    """Plots/saves histograms for raw parameters in Arc region
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0

    #---------
    h1d(np.array(sp.lst_npeaks_ch), bins=25, amp_range=(-0.5, 24.5), \
           title ='Cheetah: Number of peaks', xlabel='Number of peaks', ylabel='Frames',\
           fnm='his-npeaks-cheetah.png')
    #---------
    h1d(np.array(sp.lst_npeaks_pf), bins=25, amp_range=(-0.5, 24.5), \
           title ='Pkfinder: Number of peaks', xlabel='Number of peaks', ylabel='Frames',\
           fnm='his-npeaks-pkfdr.png')
    #---------
    fig, axhi, hi = h1d(np.array(sp.lst_matched), bins=2, amp_range=(-0.5,1.5), \
           title ='Cheetah-Pkfinder peak matching', xlabel='Non-matched : Matched', ylabel='Cheetah peaks',\
           fnm='his-peaks-matched.png')

    nnotm, nmatch= hi[0]
    print '#peaks matched: %d  not-matched: %d  of total in Cheetah: %d  frac matched: %.3f' %\
          (nmatch, nnotm, nmatch+nnotm, float(nmatch)/(nmatch+nnotm))
    #---------

    show()

#------------------------------

def peak_matched(pk_ch, pk_pf) :
    """Returns True if Cheetah peak is matched to the peakfinder peak, False othervise
    """
    #pk_pf.seg, pk_pf.row, pk_pf.col
    s, r, c = pk_ch.seg_row_col()
    if pk_pf.seg != s              : return False
    if math.fabs(pk_pf.row-r) > 5. : return False
    if math.fabs(pk_pf.col-c) > 5. : return False
    return True
  
#------------------------------
#------------------------------
#------------------------------
#------------------------------

def comparePeakLists(fnpeaks, fncheet, fnmask, ofn_prefix='plot') :
    """Main method for processing of peaks from files;
       1. load peak data from files using TDFileContainer object,
       2. compare peaks in the event
       3. plot and show results.
    """
    from pyimgalgos.TDFileContainer     import TDFileContainer
    from pyimgalgos.TDPeakRecord        import TDPeakRecord
    from pyimgalgos.TDCheetahPeakRecord import TDCheetahPeakRecord

    t0_sec = time()

    #--- load mask
    mask = np.loadtxt(fnmask)
    mask.shape = (32,185,388)

    #--- load Cheetah table of peaks
    fc_ch = TDFileContainer(fncheet, indhdr='frameNumber', objtype=TDCheetahPeakRecord) # , pbits=256)
    #fc_ch.print_content(nlines=5)

    #--- load PeakFinder table of peaks
    fc_pf = TDFileContainer(fnpeaks, indhdr='Evnum', objtype=TDPeakRecord) #, pbits=256)
    #fc_pf.print_content(nlines=5)

    #----------
    #sys.exit('Test exit') ######### TEST EXIT
    #----------
    #sp.ofn_prefix = os.path.join(os.path.dirname(ofn_prefix),'peak-idx')

    for evnum in fc_ch.group_num_iterator() :
        gr_ch = fc_ch.next()

        print '%s\nEvent# %d' % (12*'_', evnum)
        #print fc_ch.header())
    
        peaks_ch =  gr_ch.get_objs()
        sp.lst_npeaks_ch.append(len(peaks_ch))

        gr_pf = fc_pf.group(evnum-1) # accounts for offset in event numeration in Cheetah and pf.

        if gr_pf is None :
            print 'WARNING: %d peaks found in Cheetah and 0 in peakfinder...' % (len(peaks_ch))
            sp.lst_matched.append(0)
            continue

        peaks_pf = gr_pf.get_objs()
        sp.lst_npeaks_pf.append(len(peaks_pf))

        if peaks_pf[0].fid != peaks_ch[0].fid :
            print 'WARNING: Mismatch in fidelity between Cheetah and peakfinder frames %s : %s'%\
                  (peaks_ch[0].fid, peaks_pf[0].fid) 
            continue


        # Too large number of peaks
        if len(peaks_ch)>5 : continue
   
        for pk_ch in peaks_ch :
            #print pk_ch.line.rstrip('\n')
            #pk_ch.print_peak_data_short()
            #pk_ch.print_seg_row_col()

            # Discard peaks with small radius
            #if pk_ch.peak_r_assembled < 100 : continue

            # Discard peaks with large radius
            #if pk_ch.peak_r_assembled > 430 : continue

            # Cheetah peak is in our ROI ?
            s, r, c = pk_ch.seg_row_col()
            if not mask[s, int(r), int(c)] : continue

            matched = 0
            for pk_pf in peaks_pf :
                #pk_pf.print_peak_data()
                if peak_matched(pk_ch, pk_pf) :
                    matched = 1
                    break

            sp.lst_matched.append(matched)
            if not matched : 
                #pk_ch.print_peak_data_short()
                #pk_ch.print_seg_row_col()
                print 'seg:%2d  row:%3d  col:%3d  R: %7.1f  npix:%4d  totInt:%8.1f  maxInt:%7.1f  sigmaBG:%7.3f  SNR:%6.3f  resA:%.3f' %\
                    (s, int(r), int(c), pk_ch.peak_r_assembled, pk_ch.nPixels, pk_ch.totalIntensity, pk_ch.maxIntensity, pk_ch.sigmaBG, pk_ch.SNR, pk_ch.peak_resA)
                
                print 'pf peaks:'
                for pk in peaks_pf :
                    #pk.print_peak_data()
                    print 'seg:%2d  row:%3d  col:%3d  R: %7.1f  npix:%4d  totInt:%8.1f  maxInt:%7.1f  sigmaBG:%7.3f  SNR:%6.3f  phi:%.3f' %\
                        (pk.seg, pk.row, pk.col, pk.r, pk.npix, pk.atot, pk.amax, pk.rms, pk.sonc, pk.phi)

    #----------
    #    if evnum > 1000 : break ######### TEST EXIT
    #----------

    print 'Consumed time = %7.3f sec' % (time()-t0_sec)

    plot_retults(ofn_prefix, do_save=True, hwin_x0y0 = (10,10))

#------------------------------

def do_main() :
    """ Main method to do work
    """
    # Extend a permitted number of windows for histograms
    import os

    wdir = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/'
    fnpeaks = os.path.join(wdir,'pfv2-cxif5315-r0169-2016-03-28T15:02:47.txt')
    fncheet = os.path.join(wdir,'r0169-cheetah-peaks.txt')
    fnmask  = os.path.join(wdir,'../cxif5315/masks/2016-03-28-roi-mask-nda-equ.txt')
    comparePeakLists(fnpeaks, fncheet, fnmask, '2016-03-30-figs/plot-cxif5315-r0169-pfv2-vs-cheetah-v01')

#------------------------------

if __name__ == "__main__" :
    do_main()
    sys.exit('Processing is completed')

#------------------------------
