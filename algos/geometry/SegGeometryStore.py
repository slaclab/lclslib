#--------------------------------------------------------------------------
# File and Version Information:
#  $Id: SegGeometryStore.py 11400 2016-02-23 00:07:50Z dubrovin@SLAC.STANFORD.EDU $
#
# Description:
#  Module SegGeometryStore...
#
#------------------------------------------------------------------------

"""
:py:class:`algos.geometry.SegGeometryStore` - is a factory class/method to switch between different device-dependent
segments/sensors to access their pixel geometry uling :py:class:`algos.geometry.SegGeometry` interface.

Usage::

    from algos.geometry.SegGeometryStore import sgs

    sg = sgs.Create('SENS2X1:V1', pbits=0377)
    sg2= sgs.Create('EPIX100:V1', pbits=0377)
    sg3= sgs.Create('PNCCD:V1',   pbits=0377)
    sg4= sgs.Create('ANDOR3D:V1', pbits=0377)

    sg.print_seg_info(pbits=0377)
    size_arr = sg.size()
    rows     = sg.rows()
    cols     = sg.cols()
    shape    = sg.shape()
    pix_size = sg.pixel_scale_size()
    area     = sg.pixel_area_array()
    mask     = sg.pixel_mask(mbits=0377)    
    sizeX    = sg.pixel_size_array('X')
    sizeX, sizeY, sizeZ = sg.pixel_size_array()
    X        = sg.pixel_coord_array('X')
    X,Y,Z    = sg.pixel_coord_array()
    xmin = sg.pixel_coord_min('X')
    ymax = sg.pixel_coord_max('Y')
    xmin, ymin, zmin = sg.pixel_coord_min()
    xmax, ymax, zmax = sg.pixel_coord_mas()
    ...

@see other interface methods in :py:class:`algos.geometry.SegGeometry`, :py:class:`algos.geometry.SegGeometryCspad2x1V1`

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@version $Id: 2013-03-08$

@author Mikhail S. Dubrovin
"""

#--------------------------------
__version__ = "$Revision: 11400 $"
#--------------------------------

from algos.geometry.SegGeometryCspad2x1V1 import cspad2x1_one
from algos.geometry.SegGeometryEpix100V1  import epix2x2_one
from algos.geometry.SegGeometryMatrixV1   import SegGeometryMatrixV1, segment_one, matrix_pars

#------------------------------

class SegGeometryStore() :
    """Factory class for SegGeometry-base objects of different detectors"""

#------------------------------

    def __init__(sp) :
        pass

#------------------------------

    def Create(sp, segname='SENS2X1:V1', pbits=0 ) :
        """ Factory method returns device dependent SINGLETON object with interface implementation  
        """        
        if segname=='SENS2X1:V1' : return cspad2x1_one # SegGeometryCspad2x1V1(use_wide_pix_center=False)
        if segname=='EPIX100:V1' : return epix2x2_one  # SegGeometryEpix100V1(use_wide_pix_center=False)
        if segname=='PNCCD:V1'   : return segment_one  # SegGeometryMatrixV1()
        #if segname=='ANDOR3D:V1' : return seg_andor3d  # SegGeometryMatrixV1(rows=2048, cols=2048, ...)
        if segname[:4]=='MTRX'   :
            rows, cols, psize_row, psize_col = matrix_pars(segname)
            return SegGeometryMatrixV1(rows, cols, psize_row, psize_col,\
                                       pix_size_depth=100,\
                                       pix_scale_size=min(psize_row, psize_col))
        return None

#------------------------------

sgs = SegGeometryStore()

#------------------------------
#----------- TEST -------------
#------------------------------

def test_seggeom() :

    import sys

    from time import time
    t0_sec = time()

    if len(sys.argv)==1   : print 'For test(s) use command: python', sys.argv[0], '<test-number=1-3>'

    elif(sys.argv[1]=='1') :
        sg = sgs.Create('SENS2X1:V1', pbits=0377)
        sg.print_seg_info(pbits=0377)
        
    elif(sys.argv[1]=='2') :
        sg = sgs.Create('EPIX100:V1', pbits=0377)
        sg.print_seg_info(pbits=0377)

    elif(sys.argv[1]=='3') :
        sg = sgs.Create('PNCCD:V1', pbits=0377)
        sg.print_seg_info(pbits=0377)

    elif(sys.argv[1]=='4') :
        sg = sgs.Create('MTRX:512:512:54:54', pbits=0377)
        print 'Consumed time for MTRX:512:512:54:54 (sec) =', time()-t0_sec
        sg.print_seg_info(pbits=0377)
  
    else : print 'Non-expected arguments: sys.argv=', sys.argv, ' use 0,1,2,...'

#------------------------------

if __name__ == "__main__" :
    test_seggeom()
    print 'End of test.'

#------------------------------
