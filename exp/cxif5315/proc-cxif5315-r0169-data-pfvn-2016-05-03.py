#!/usr/bin/env python

##-----------------------------
"""Example: processes cxif5315-r0169 cspad data, finds peaks,
   makes plot (if requested) and saves peaks in (small data) file.

Usage::
   # This script uses masks explicitly from ./cxif5315/masks/
   # so run it from release directory, which have users package cxif5315:
   python ./cxif5315/proc-cxif5315-r0169-data-pfvn-2016-05-03.py

   bsub -q psfehq -o log-r169-pfvnr1.log python ./cxif5315/proc-cxif5315-r0169-data-pfvn-2016-05-03.py
    
   # Parameters to edit:
   EVTMAX, DO_PLOT, file names etc.
"""
##-----------------------------

import os
import sys
import psana
import math
import numpy as np
from time import time
from Detector.AreaDetector  import AreaDetector
from ImgAlgos.PyAlgos       import PyAlgos, print_arr, print_arr_attr

from pyimgalgos.PeakStore   import PeakStore
from pyimgalgos.GlobalUtils import subtract_bkgd
from pyimgalgos.RadialBkgd  import RadialBkgd, polarization_factor

##-----------------------------
# Initialization of graphics
#from pyimgalgos.GlobalGraphics import store as sp
import pyimgalgos.GlobalGraphics as gg
import pyimgalgos.Graphics       as gr

from CalibManager.PlotImgSpeWidget import add_stat_text

##-----------------------------

ntest = int(sys.argv[1]) if len(sys.argv)>1 else 1
print 'Example # %d' % ntest

##-----------------------------
SKIP        = 0

EVTMAX      = 100000 + SKIP
#EVTMAX      = 10 + SKIP

DO_SPEC     = False
#DO_SPEC     = True
DO_PLOT     = False
#DO_PLOT     = True
EVTPLOT     = 1 

DIST_STOD   = 94000 # Sample-to-detector distance in um for  polarization_factor, etc

#WDIR = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work'
WDIR = 'cxif5315/masks'

##-----------------------------

def do_print(i) :
    """Returns True for events which need to be printed"""
    return True
    #return False
    #if i==1 : return True
    #return not i%10

##-----------------------------

runnum = 169
dsname = 'exp=cxif5315:run=%d' % (runnum)
src    = 'CxiDs2.0:Cspad.0'
print '%s\nExample for\n  dataset: %s\n  source : %s' % (85*'_',dsname, src)

# Non-standard calib directory
#psana.setOption('psana.calib-dir', './empty/calib')
psana.setOption('psana.calib-dir', '/reg/d/psdm/CXI/cxif5315/calib')

ds  = psana.DataSource(dsname)
env = ds.env()

#runnum = evt.run()
#evt = ds.events().next()
#run = ds.runs().next()
#runnum = run.run()

##-----------------------------

det = AreaDetector(src, env)
print 85*'_', '\nInstrument: %s  run number: %d' % (det.instrument(), runnum)

nda_peds  = det.pedestals(runnum)
nda_bkgd  = det.bkgd(runnum)
nda_smask = det.mask(runnum, calib=False, status=True, edges=True, central=True, unbond=True, unbondnbrs=True)

#print_arr_attr(nda_peds,  'nda_peds')
#print_arr_attr(nda_bkgd,  'nda_bkgd')
#print_arr_attr(nda_smask, 'nda_smask')
##-----------------------------

shape_cspad = (32,185,388)

mask_arc = np.loadtxt(os.path.join(WDIR, '2016-03-28-roi-mask-nda-arc.txt'))
mask_equ = np.loadtxt(os.path.join(WDIR, '2016-03-28-roi-mask-nda-equ.txt'))
mask_tot = np.loadtxt(os.path.join(WDIR, '2016-03-28-roi-mask-nda-equ-arc.txt'))

mask_arc.shape = mask_equ.shape = mask_tot.shape = shape_cspad
print_arr_attr(mask_arc, 'mask_arc')

seg1 = np.ones((185,388))
regs_bkgd_norm = np.zeros(shape_cspad, dtype=np.int16)
for s in (4,12,20,28) : regs_bkgd_norm[s,10:100,270:370] = 20*seg1[10:100,270:370]

winds_bkgd = [(s, 10, 100, 270, 370) for s in (4,12,20,28)] # use part of segments 4 and 20 to subtr bkgd

winds_arc = [(s, 0, 185, 0, 388) for s in (0,7,8,15)]
winds_equ = [(s, 0, 185, 0, 388) for s in (0,1,9,15,16,17,25,31)]
winds_tot = [(s, 0, 185, 0, 388) for s in (0,1,7,8,9,15,16,17,23,24,25,31)]

mask_winds_tot = np.zeros(shape_cspad, dtype=np.int16)
mask_winds_tot[(0,1,7,8,9,15,16,17,23,24,25,31),:,:] = seg1
#mask_winds_equ[(0,1,9,15,16,17,25,31),:,:] = seg1
#mask_winds_arc[(0,7,8,15),:,:] = seg1

print_arr(winds_arc, 'winds_arc')
print_arr_attr(winds_arc, 'winds_arc')

alg_arc = PyAlgos(windows=winds_arc, mask=mask_arc, pbits=2)
alg_arc.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=10)
#alg_arc.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=500, son_min=6) # for v2r1

alg_equ = PyAlgos(windows=winds_equ, mask=mask_equ, pbits=0)
alg_equ.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=10)
#alg_equ.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=500, son_min=6) # for v2r1

#alg_equ.print_attributes()
#alg_equ.print_input_pars()

##-----------------------------

xoffset, yoffset = 300, 300
xsize,   ysize   = 1150, 1150

# Pixel image indexes
iX = np.array(det.indexes_x(runnum), dtype=np.int64) #- xoffset
iY = np.array(det.indexes_y(runnum), dtype=np.int64) #- yoffset

ix_cent, iy_cent = det.point_indexes(runnum) 
print 'Center indexes for CSPAD geometry: ix, iy = %.1f, %.1f' % (ix_cent, iy_cent)

# Protect indexes (should be POSITIVE after offset subtraction)
imRow = np.select([iX<xoffset], [0], default=iX-xoffset)
imCol = np.select([iY<yoffset], [0], default=iY-yoffset)

Xarr =  det.coords_x(runnum)
Yarr =  det.coords_y(runnum)

# Pixel coordinates [um] (transformed as needed)
Xum =  Yarr
Yum = -Xarr

mask_bkgd = nda_smask # * mask_winds_tot

#rb = RadialBkgd(Xarr, Yarr, mask=mask_bkgd, radedges=(5200, 80000), nradbins=200, nphibins=1)
#rb = RadialBkgd(Xarr, Yarr, mask=mask_bkgd, radedges=None, nradbins=300, phiedges=(-45,315), nphibins=4)
#pf = polarization_factor(rb.pixel_rad(), rb.pixel_phi()+90, DIST_STOD)
#rb.print_attrs()
#rb.print_ndarrs()
print 80*'_'

# Derived pixel raduius in [um] and angle phi[degree]
Rum = np.sqrt(Xum*Xum + Yum*Yum)
Phi = np.arctan2(Yum,Xum) * 180 / np.pi

imRow.shape  = imCol.shape  = \
Xum.shape    = Yum.shape    = \
Rum.shape    = Phi.shape    = shape_cspad

#addhdr = '  Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent '+\
#         'rsigma  csigma rmin rmax cmin cmax    bkgd     rms     son  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]'
addhdr = '  Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent '+\
         'rsigma  csigma    bkgd     rms     son  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]'

fmt = '%7d  %3s  %3d %4d %4d  %4d  %8.1f  %8.1f  %6.1f  %6.1f %6.2f  %6.2f'+\
      '  %6.2f  %6.2f  %6.2f'+\
      ' %6d  %6d  %8.0f  %8.0f  %8.0f  %8.2f'

pstore = PeakStore(env, runnum, prefix='peaks', add_header=addhdr, pbits=1)
pstore.print_attrs()

##-----------------------------
fig, axim, axcb = gg.fig_axes() if DO_PLOT else (None, None, None)

fighi, axhi = None, None
if DO_SPEC :
    fighi = gr.figure(figsize=(6,5), title='Spectrum of pixel intensities', move=(800,0))
    axhi  = gr.add_axes(fighi, (0.17, 0.12, 0.78, 0.78))

##-----------------------------

def geo_pars(s,r,c) :
    inds = (s,r,c)
    return imRow[inds], imCol[inds], Xum[inds], Yum[inds], Rum[inds], Phi[inds]

##-----------------------------

t0_sec_evloop = time()
npeaks_tot = 0
nda = None
peaks = None

# loop over events in data set
for i, evt in enumerate(ds.events()) :

    if do_print(i) and i%100==0 : print 'Event %d' % (i)

    #for key in evt.keys() : print key

    if i<SKIP    : continue
    if i>=EVTMAX : break

    # get calibrated data ndarray and proccess it if it is available
    t1_sec = time()

    #nda_data = det.calib(evt, cmpars=(1,50,50,100)) # in calib (1 25 25 100)
    nda_data = det.raw(evt)

    if nda_data is not None :

        nda =  np.array(nda_data, dtype=np.float32, copy=True)
        nda -= nda_peds

        #det.common_mode_apply(runnum, nda, cmpars=(1,50,50,1000))
        det.common_mode_apply(runnum, nda, cmpars=(5,50))

        nda = subtract_bkgd(nda, nda_bkgd, mask=nda_smask, winds=winds_bkgd, pbits=0)

        ##nda = rb.subtract_bkgd_interpol(nda.flatten() * pf)
        #nda = rb.subtract_bkgd(nda.flatten() * pf)
        #nda.shape = shape_cspad

        nda *= nda_smask

        print '  ----> calibration dt = %f sec' % (time()-t1_sec)

        #print_arr_attr(nda, 'calibrated data')
        t0_sec = time()

        # run peakfinders and get list of peak records for each region
        #peaks_arc = alg_arc.peak_finder_v2r1(nda, thr=30, r0=7, dr=2)
        peaks_arc = alg_arc.peak_finder_v3r1(nda, rank=5, r0=7, dr=2, nsigm=0) # 1.64 (5%)
        #peaks_arc = alg_arc.peak_finder_v4r1(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)

        #peaks_equ = alg_equ.peak_finder_v2r1(nda, thr=30, r0=7, dr=2)
        peaks_equ = alg_equ.peak_finder_v3r1(nda, rank=5, r0=7, dr=2, nsigm=0) # 1.64 (5%)
        #peaks_equ = alg_equ.peak_finder_v4r1(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)

        # available after v2r1 ONLY!
        #maps_of_conpix_arc = alg_arc.maps_of_connected_pixels()
        #maps_of_conpix_equ = alg_equ.maps_of_connected_pixels()

        # available after v3r1 ONLY!
        #maps_of_locmax_arc = alg_arc.maps_of_local_maximums()
        #maps_of_locmax_equ = alg_equ.maps_of_local_maximums()

        ###===================
        if do_print(i) : print '%s\n%s\n%s\n%s' % (85*'_', pstore.header[0:66], pstore.rec_evtid(evt), addhdr)
        ###===================
        npeaks = 0

        peak_reg_lists = zip(('ARC','EQU'), (peaks_arc, peaks_equ)) 

        # loop over ARC and EQU regions
        for reg, peak_list in peak_reg_lists :

            # loop over peaks found in the region
            for peak in peak_list :

                # get peak parameters
                seg,row,col,npix,amax,atot,rcent,ccent,rsigma,csigma,\
                rmin,rmax,cmin,cmax,bkgd,rms,son = peak[0:17]

                # get pixel coordinates
                imrow, imcol, xum, yum, rum, phi = geo_pars(seg, row, col)
                
                # make peak-record and save it in the file
                rec = fmt % (i, reg, seg, row, col, npix, amax, atot, rcent, ccent, rsigma, csigma,\
                      bkgd, rms, son,\
                      imrow, imcol, xum, yum, rum, phi)
            
                pstore.save_peak(evt, rec)
                npeaks += 1
                ###===================
                # if do_print(i) : print '%s' % rec
                ###===================

        npeaks_tot += npeaks

        ###===================
        if do_print(i) : print 'Event %d, npeaks=%3d of total=%6d --- dt/evt = %f sec' % (i, npeaks, npeaks_tot, time()-t0_sec)
        ###===================

        if DO_PLOT and i%EVTPLOT==0 :
            #nda = bkg
            #nda = pf
            #nda = nda_bkgd
            #nda = nda_bkgd + regs_bkgd_norm      
            #img = det.image(evt, nda)

            #img = det.image(evt, maps_of_conpix_equ)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, mask_tot*nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, mask_tot*nda + nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            img = det.image(evt, nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, bkg)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]

            #ave, rms = img.mean(), img.std()
            #print 'ave, rms', ave, rms
            #amin, amax = ave-1*rms, ave+8*rms
            #amin, amax = -20, 20
            amin, amax = -50, 100
            gg.plot_img(img, mode='do not hold', amin=amin, amax=amax)
            gg.plot_peaks_on_img(peaks_arc, axim, imRow, imCol, color='w') #, pbits=3)
            gg.plot_peaks_on_img(peaks_equ, axim, imRow, imCol, color='w') #, pbits=3)

            fig.canvas.set_window_title('Event: %d' % i)    
            fig.canvas.draw() # re-draw figure content

            #gg.plotHistogram(nda, amp_range=(-100,100), bins=200, title='Event %d' % i)
            
        if DO_SPEC :
            nda4 = np.array(nda)
            nda4.shape = (4,8,185,388)
            seg=0
            hiarr = nda4[:,seg,:,:].flatten()
            rangehi = (-50,100)
            hi = gr.hist(axhi, hiarr, bins=150, amp_range=rangehi, weights=None, color=None, log=False)
            gr.add_title_labels_to_axes(axhi, title='Spectrum seg-%d, event: %d'%(seg,i), xlabel='Intensity (ADU)', ylabel='Pixels')
            fighi.canvas.draw() # re-draw figure content
            #gr.save_fig(fighi, fname='img-spe-seg%d-ev%06d-sub-bkgd.png'%(seg,i))

##-----------------------------

print ' ----> Event loop time = %f sec, npeaks = %d' % (time()-t0_sec_evloop, npeaks_tot)
#pstore.close_file()
gg.show() # hold image untill it is closed

##-----------------------------

sys.exit('Processing is completed')

##-----------------------------
