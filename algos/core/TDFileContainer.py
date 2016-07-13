#!/usr/bin/env python

#--------------------------------------------------------------------------
# File and Version Information:
#  $Id: TDFileContainer.py 11634 2016-04-05 18:47:50Z dubrovin@SLAC.STANFORD.EDU $
#
# Description:
#  class TDFileContainer
#
#------------------------------------------------------------------------

"""TDFileContainer - text/table data file container - load/hold/provide access to data from text file.

It is assumed that text data file contains records of the same format and occasionally record-header
beginning with character # (hash in [0] position).
Example of the file content::

    # Exp     Run  Date       Time      time(sec)   time(nsec) fiduc  Evnum  Reg  Seg  Row  Col  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  486382070  104421     0  EQU   17  153   48  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  494719789  104424     1  EQU    1  161   32  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  494719789  104424     1  EQU   17  170   51  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  503058551  104427     2  EQU   25  170  310  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  503058551  104427     2  EQU   25  180  292  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  511393301  104430     3  EQU    1  162   27  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  536405573  104439     6  ARC    8   11   41  ...
    cxif5315  169  2015-02-22 02:20:47  1424600447  536405573  104439     6  ARC    8   10   20  ...
    ...

Header (without #) should have the same as data number of literal fields separated by spaces.
Records in the file should be grupped by unique group-id parameter,
for example a group of records may have the same group number or some unique index.


Originaly it is designed to work with text file containing record data generated by peak-finder.
It is adopted to work with any other object type beside peak data.

Usage::

    # !!! NOTE: None is returned whenever requested information is missing.

    # Import
    from algos.core.TDFileContainer     import TDFileContainer
    from algos.core.TDPeakRecord        import TDPeakRecord # use it by default in TDFileContainer
    from algos.core.TDNodeRecord        import TDNodeRecord
    from algos.core.TDCheetahPeakRecord import TDCheetahPeakRecord

    # Initialization
    # for peakfinder records
    fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/pfv2-cxif5315-r0169-2015-09-14T14:28:04.txt'
    fc = TDFileContainer(fname, indhdr='Evnum', objtype=TDPeakRecord, pbits=0)

    # for index table:
    fc = TDFileContainer(fname, indhdr='index', objtype=TDNodeRecord)

    # for Cheetah file with peaks:
    fc = TDFileContainer(fname, indhdr='frameNumber', objtype=TDCheetahPeakRecord)

    gr_nums = fc.group_numbers()
    ngrps   = fc.number_of_groups()
    grnum   = fc.current_group_number()

    gr_curr = fc.group(grpnum) # returns current or specified group
    gr_next = fc.next()        # returns next group
    gr_prev = fc.previous()    # returns previous group
    hdr     = fc.header()        

    # Print
    fc.print_attrs()
    fc.print_content(nlines=None) # prints nline (or all by default) lines from file conteiner
 
    # ____________________________________

    # Example of iterations over groups
    for grnum in fc.group_num_iterator() :
        group = fc.next()
        group.print_attrs()
        peaks = group.get_objs()

        for pk in peaks :
            pk.print_short()

            # Information available through the TDPeakRecord object pk
            # ________________________________________________________
            # pk.exp, pk.run, pk.evnum, pk.reg
            # pk.date, pk.time, pk.tsec, pk.tnsec, pk.fid
            # pk.seg, pk.row, pk.col, pk.amax, pk.atot, pk.npix
            # pk.rcent, pk.ccent, pk.rsigma, pk.csigma
            # pk.rmin, pk.rmax, pk.cmin, pk.cmax
            # pk.bkgd, pk.rms, pk.son
            # pk.imrow, pk.imcol
            # pk.x, pk.y, pk.r, pk.phi
            # pk.sonc
            # pk.dphi000
            # pk.dphi180
            # pk.line

    # Example of direct access to group by its number
    grpnum = 8 # but grpnum is not necessaraly conecutive number, it should be in fc.group_num_iterator() ...
    group = fc.group(grpnum) # returns current or specified group
    group.print_attrs()

This software was developed for the LCLS project.
If you use all or part of it, please give an appropriate acknowledgment.

@see classes
\n  :py:class:`pyimgalgos.TDFileContainer` - file records container. 
\n  :py:class:`pyimgalgos.TDGroup` - holds a list of records associated with a single group.
\n  :py:class:`pyimgalgos.TDPeakRecord` - provides access to the peak record.
\n  :py:class:`pyimgalgos.TDNodeRecord` - provides access to the look-up table with crystal orientation record.
\n  :py:class:`pyimgalgos.TDCheetahPeakRecord` - provides access to the Cheetah peak record.

@version $Id: TDFileContainer.py 11634 2016-04-05 18:47:50Z dubrovin@SLAC.STANFORD.EDU $

@author Mikhail S. Dubrovin
"""
#------------------------------
__version__ = "$Revision: 11634 $"
# $Source$
##-----------------------------

import os
#import sys
from time import time

from algos.core.TDGroup import TDGroup
from algos.core.TDPeakRecord  import TDPeakRecord

##-----------------------------
##-----------------------------

class TDFileContainer :
    """ Load and hold record list from file and provide access by group index
    """
    def __init__(self, fname, indhdr='Evnum', objtype=TDPeakRecord, pbits=0) :
        """Constructor
           Args:
           fname   (str) - text table data file name 
           indhdr  (str) - header of the field used for group indexing
           objtype (TD*Recor) - object type used for data record processing/access
           pbits   (int) - print control bit-word; pbits & 256 - tracking
        """
        if pbits & 256 : print 'c-tor of class %s' % self.__class__.__name__
        self.indhdr = indhdr
        self.objtype = objtype
        self.pbits = pbits
        self.hdr = None
        self.grnum = -1
        #self.lst_of_recs = [] # list of recs loaded from record data file
        self.lst_grnum    = [] # list of group numbers in the data file
        self.lst_begin    = [] # list of record indexes in the lst_of_recs
        self.lst_nrecords = [] # list of numbor of records in group 
        
        self._load_recs_from_file(fname)
        self._group_indexing()
        self._reset_indexes()

##-----------------------------

    def _reset_indexes(self) :
        """ resets indexes for iterator
        """        
        self.first_iteration = True
        self.grnum_curr = self.lst_grnum[0] # reset current group after indexing
        self.indlst_curr = 0                # reset current index of internal lists
         
##-----------------------------

    def __del__(self) :
        """d-tor
        """
        if self.pbits & 256 : print 'd-tor of class %s' % self.__class__.__name__
        pass

##-----------------------------

    def __call__(self) :
        """ Alias to group_num_iterator()
        """
        self.group_num_iterator()

##-----------------------------

    def print_content(self, nlines=None) :
        """ Prints content of the file-container; by default-entire file.
        """
        if self.pbits & 256 : print """default method of class %s""" % self.__class__.__name__

        print '\n', 120*'_', '\n%s holds data from file:\n  %s\n' % (self.__class__.__name__, self.fname)
        for i,rec in enumerate(self.lst_of_recs) :
            if nlines is not None and i>nlines : break
            print rec,
        print 'etc.' if nlines is not None else 'End of file'

##-----------------------------

    def print_attrs(self) :
        print 'Attributes of the class %s object' % self.__class__.__name__
        print '  fname : %s' % self.fname,\
              '\n  pbits : %d' % self.pbits,\
              '\n  hdr   : %s' % self.hdr,\
              '\n  nrecs : %d' % len(self.lst_of_recs),\
              '\n  Auto-defined grnum index in the record data : %d' % self.igrnum

##-----------------------------

    def _load_recs_from_file(self, fname) :
        if not os.path.lexists(fname) : raise IOError('File %s is not found' % fname)
        self.fname = fname
        t0_sec = time()
        f=open(fname,'r')
        self.lst_of_recs = []
        for rec in f : self.lst_of_recs.append(rec.replace(',',' '))
        f.close()
        if self.pbits & 256 : print 'File loading time %.3f sec' % (time()-t0_sec)

##-----------------------------

    def _load_recs_from_file_v0(self, fname) :
        if not os.path.lexists(fname) : raise IOError('File %s is not found' % fname)
        self.fname = fname
        t0_sec = time()
        f=open(fname,'r')
        self.lst_of_recs = f.readlines()
        f.close()
        if self.pbits & 256 : print 'File loading time %.3f sec' % (time()-t0_sec)

##-----------------------------

    def _part_rec_parser(self, rec) :
        """ 1. saves the 1st header in self.hdr, return None for header
            2. defines index of the field self.indhdr (='Evnum')
            3. returns None for empty recs (if any)
            4. returns group number found in the record data
        """
        if len(rec)==1 : return None # ignore empty records

        if rec[0]=='#' : # rec is header or comment
            if  self.hdr is None :
                if not (self.indhdr in rec) : return None

                self.hdr = rec.lstrip('#').rstrip('\n')
                self.igrnum = self.hdr.split().index(self.indhdr)
                if self.pbits & 256 : print 'self.igrnum', self.igrnum
            return None

        # partly split data fields and return group number
        fields = rec.split(None,self.igrnum+1)    
        return int(fields[self.igrnum])

##-----------------------------

    def _group_indexing(self) :
        """loops over list of records, makes lists for indexing
        """
        if self.pbits & 256 : print '_group_indexing'
        t0_sec = time()

        self.count = 0

        for ind, rec in enumerate(self.lst_of_recs) :

            grnum = self._part_rec_parser(rec)
            if grnum is None : continue # in case of comments and empty recs

            # check if record is from the next group and add it to the list
            if grnum != self.grnum :                
                if not (self.grnum < 0) : # skip 1st record
                    self.lst_nrecords.append(self.count) 
                self.count = 1
                self.grnum = grnum
                self.lst_grnum.append(grnum)
                self.lst_begin.append(ind)
                 
                #print 'New group number: %d' % grnum
            else :
                self.count += 1

            #==== TEST ======
            #if ind>100 : break
            #print rec
            #================
 
        self.lst_nrecords.append(self.count) # add for last record

        if self.pbits & 256 :
            print 'Last group %d contains %d records' % (self.grnum, self.count)
            print 'Group indexing time %.3f sec' % (time()-t0_sec)

##-----------------------------
# This is time consuming operation
#    def list_of_groups(self) :
#        """returns list of group objects
#        """
#        self._reset_indexes()
#        return [self.next() for grnum in self.lst_grnum]
#
##-----------------------------

    def group_numbers(self) :
        """returns list of group numbers in the file
        """
        return self.lst_grnum

##-----------------------------

    def group_num_iterator(self) :
        """resets indexes to the beginning of arrays and returns list of group numbers
        """
        self._reset_indexes()
        return self.lst_grnum

##-----------------------------

    def number_of_groups(self) :
        """returns number of groups in file
        """
        return len(self.lst_grnum)

##-----------------------------

    def current_group_number(self) :
        """returns current group number
        """
        return self.grnum_curr

##-----------------------------

    def header(self) :
        """returns string header
        """
        return self.hdr

##-----------------------------

    def _group_for_index(self) :
        """returns group for specified range of indexes
        """
        self.grnum_curr = self.lst_grnum [self.indlst_curr]
        begin           = self.lst_begin [self.indlst_curr]
        nrecords          = self.lst_nrecords[self.indlst_curr]

        if self.pbits & 256 : 
            print 'grnum_curr=%d  indlst_curr=%d  begin=%d  nrecords=%d' %\
                  (self.grnum_curr, self.indlst_curr, begin, nrecords)

        evt_recs = self.lst_of_recs[begin:begin+nrecords]
        #print '%s\nList of records for group %d' % (80*'_', self.grnum)
        #for rec in recs : print rec

        return TDGroup(evt_recs, self.objtype, pbits=self.pbits)

##-----------------------------

    def group(self, grnum=None) :
        """returns current or specified group
        """
        if self.pbits & 256 : print 'group(evnum=%s)' % str(grnum)
        if grnum is not None :
            if not (grnum in self.lst_grnum) : return None
            self.indlst_curr = self.lst_grnum.index(grnum)        
        return self._group_for_index()

##-----------------------------

    def next(self) :
        """returns next group
        """
        if self.pbits & 256 : print 'next group'

        if  self.first_iteration :
            self.first_iteration = False
            return self._group_for_index() # do not increment indexes on first iteration

        if  self.indlst_curr < len(self.lst_grnum)-1 :
            self.indlst_curr += 1
            return self._group_for_index()
        else :
            if self.pbits : print 'WARNING: %s.next() reached the end of the list, return None'%\
               self.__class__.__name__
            return None

##-----------------------------

    def previous(self) :
        """returns previous group
        """
        if self.pbits & 256 : print 'previous group'

        if  self.first_iteration :
            self.first_iteration = False
            return self._group_for_index() # do not decrement indexes on first iteration

        if  self.indlst_curr > 0 :
            self.indlst_curr -= 1
            return self._group_for_index()
        else :
            if self.pbits : print 'WARNING: %s.previous() reached the beginning of the list, return None'%\
               self.__class__.__name__
            return None


##-----------------------------
##-----------------------------
## Aliases for depricated names
##-----------------------------
##-----------------------------

    def event_numbers(self) :
        """Depricated, see group_numbers()"""
        return self.group_numbers()

##-----------------------------

    def evnum_iterator(self) :
        """Depricated, see group_num_iterator()"""
        return self.group_num_iterator()

##-----------------------------

    def number_of_events(self) :
        """Depricated, see number_of_groups()"""
        return self.number_of_groups()

##-----------------------------

    def current_event_number(self) :
        """Depricated, see current_group_number()"""
        return self.current_group_number()

##-----------------------------

    def event(self, evnum=None) :
        """Depricated, see group(evnum)"""
        return self.group(grnum=evnum)

##-----------------------------
##-----------------------------
##-----------------------------
##-----------------------------
##-----------------------------

def do_work() :
    """ Test
    """
    fname = '/reg/neh/home1/dubrovin/LCLS/rel-mengning/work/pfv2-cxif5315-r0169-2015-09-14T14:28:04.txt'
    fc = TDFileContainer(fname, indhdr='Evnum', objtype=TDPeakRecord, pbits=0)
    fc.print_attrs()
    fc()

    # Direct access to TDGroup object
    group = fc.group(8)
    group.print_attrs()

    t0_sec = time()

    for grpnum in fc.group_num_iterator() :
        group = fc.next()
        print '%s\nGroup %d  ' % (80*'_', grpnum)
        for record in group() :
            print ' ',
            record.print_short()

        for i, peak in enumerate(group()) :
            print '  peak#%2d  bkgd=%5.1f  rms=%5.1f  S/N=%5.1f' % (i, peak.bkgd, peak.rms, peak.son)
        
    print '\nTime to iterate using next() %.3f sec' % (time()-t0_sec)

    #t0_sec = time()
    #groups = fc.list_of_groups()  
    #print 'Time to generate list of group objects %.3f sec' % (time()-t0_sec)

##-----------------------------

if __name__ == "__main__" :
    do_work()
    print('Test is completed')
    #sys.exit('Processing is completed')

##-----------------------------