#!/usr/bin/env python

##-----------------------------
"""Example of analysis script:
   processes cxif5315-r0169 peak list with phi-beta fit info;
   loop over all entries of the peak-fit file and create histograms for fit parameters
"""
##-----------------------------

import sys
import os
import math
import numpy as np
from time import time # strftime, localtime #, gmtime, time

from pyimgalgos.TDFileContainer import TDFileContainer
from pyimgalgos.TDPeakRecord    import TDPeakRecord
from pyimgalgos.GlobalGraphics import hist1d, show, move, save #, move_fig, save_fig, fig_axes, plot_img, plot_peaks_on_img

#------------------------------

class Storage :
    """Store for shared parameters."""

    prefix  = './his-equ-fit' # file name prefix for histograms
    DO_HIST = True        # save histograms

    def __init__(self) :

        self.h1d_phi      = []
        self.h1d_beta     = []
        self.h1d_phi_err  = []
        self.h1d_beta_err = []
        self.h1d_prob     = []

sp = Storage()

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
    h1d(np.array(sp.h1d_phi), bins=80, amp_range=(-40.,40.), \
        title ='Equ: fit phi angle', xlabel='$\phi$ (deg)', ylabel='Events', fnm='phi.png')
    #---------
    h1d(np.array(sp.h1d_beta), bins=100, amp_range=(-60.,40.), \
        title ='Equ: fit beta angle', xlabel='beta (deg)', ylabel='Events', fnm='beta.png')
    #---------
    h1d(np.array(sp.h1d_phi_err), bins=100, amp_range=(0.,10.), \
        title ='Equ: fit phi angle', xlabel='error of $\phi$ (deg)', ylabel='Events', fnm='phi-err.png')
    #---------
    h1d(np.array(sp.h1d_beta_err), bins=100, amp_range=(0.,50.), \
        title ='Equ: fit beta angle', xlabel='error of beta (deg)', ylabel='Events', fnm='beta-err.png')
    #---------
    h1d(np.array(sp.h1d_prob), bins=100, amp_range=(0.,1.), \
        title ='Equ: fit probability', xlabel='probability', ylabel='Events', fnm='prob.png')
    #---------

    show()

#------------------------------

def proc_file(fname) :
    """Process file with peak records extended by the phi-beta fit info
    """

    t0_sec = time()
    sp.fc = TDFileContainer(fname, indhdr='Evnum', objtype=TDPeakRecord) #, pbits=256)
    print 'File load time = %7.3f sec' % (time()-t0_sec)
    #sp.fc.print_content(nlines=20)

    for evnum in sp.fc.group_num_iterator() :
        event = sp.fc.next()
        lst_peaks = event.get_objs()
        pk0 = lst_peaks[0]
        print '%s\nEvent# %6d exp=%s:run=%d %s %s' % (80*'_', evnum, pk0.exp, pk0.run, pk0.date, pk0.time)
        print ' %s\n%s' % (sp.fc.hdr, pk0.line)

        #for pk in lst_peaks : print 'x=%.6f, y=%.6f' % (px.x, pk.y)

        # fill histogram arrays for event
        sp.h1d_phi     .append(pk0.fit_phi)
        sp.h1d_beta    .append(pk0.fit_beta)
        sp.h1d_phi_err .append(pk0.fit_phi_err)
        sp.h1d_beta_err.append(pk0.fit_beta_err)
        sp.h1d_prob    .append(pk0.fit_prob)

    if sp.DO_HIST : plot_histograms()

#------------------------------

if __name__ == "__main__" :
    proc_file('/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/pfv2-fit-cxif5315-r0169-2016-04-06T08:32:23.txt')
    sys.exit('Processing is completed')

#------------------------------
