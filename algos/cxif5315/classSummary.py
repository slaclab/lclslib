#bsub -a mympi -n 36 -q psanaq -o %J.log python classSummary.py

from psana import *
from Detector.PyDetector import PyDetector
import numpy as np
import sys

class Stats:
    def __init__(self,detarr,exp,run,detname):
        self.sum=detarr.astype(np.float64)
        self.sumsq=detarr.astype(np.float64)*detarr.astype(np.float64)
        self.maximum=detarr.astype(np.float64)
        self.nevent=1
        self.exp = exp
        self.run = run
        self.detname = detname
    def update(self,detarr):
        self.sum+=detarr
        self.sumsq+=detarr*detarr
        self.maximum=np.maximum(self.maximum,detarr)
        self.nevent+=1
    def store(self,tag=''):
        self.totevent = comm.reduce(self.nevent)
        if rank==0:
            comm.Reduce(MPI.IN_PLACE,self.sum)
            comm.Reduce(MPI.IN_PLACE,self.sumsq)
            comm.Reduce(MPI.IN_PLACE,self.maximum,op=MPI.MAX)
            # Accumulating floating-point numbers introduces errors,
            # which may cause negative variances.  Since a two-pass
            # approach is unacceptable, the standard deviation is
            # clamped at zero.
            self.mean = self.sum / float(self.totevent)
            self.variance = (self.sumsq / float(self.totevent)) - (self.mean**2)
            self.variance[self.variance < 0] = 0
            self.stddev = np.sqrt(self.variance)
            file = '%s_%4.4d_%s_%s'%(self.exp,self.run,self.detname,tag)
            print 'writing file',file
            np.savez(file,mean=self.mean,stddev=self.stddev,max=self.maximum)
        else:
            comm.Reduce(self.sum,self.sum)
            comm.Reduce(self.sumsq,self.sumsq)
            comm.Reduce(self.maximum,self.maximum,op=MPI.MAX)

ds = DataSource('exp=cxif5315:run=169:idx')
env = ds.env()

from smallData import getSmallData
cheetahFile,smallData = getSmallData('peaks-cxif5315-r0169-2015-09-16T10:52:30.txt')
hitTimes = [EventTime((sd.timesec<<32)|sd.timensec,sd.fiduc) for sd in smallData]

# set this to sys.maxint to analyze all events
maxevents = sys.maxint

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

detname = ['CxiDs2.0:Cspad.0']
srclist = [Source('DetInfo(%s)'%d) for d in detname]
detlist = [PyDetector(s, env) for s in srclist]
for d,n in zip(detlist,detname):
    d.detname = n

nevent = np.array([0])

for run in ds.runs():
    runnumber = run.run()
    # list of all events
    times = run.times()
    nevents = min(len(times),maxevents)
    mylength = nevents/size # easy but sloppy. lose few events at end of run.

    # chop the list into pieces, depending on rank
    mytimes= times[rank*mylength:(rank+1)*mylength]

    hits = []
    print len(mytimes),len(hitTimes)
    for t in mytimes:
        hit = False
        for ht in hitTimes:
            if t.seconds()==ht.seconds() and t.fiducial()==ht.fiducial():
                hit = True
                break
        hits.append(hit)
    
    for i in range(mylength):
        if i%10==0: print 'Rank',rank,'processing event',rank*mylength+i
        evt = run.event(mytimes[i])
        # very useful for seeing what data is in the event
        #print evt.keys()
        if evt is None:
            print '*** event fetch failed'
            continue
        for d in detlist:
            try:
                detarr = d.calib_data(evt)
            except ValueError:
                id = evt.get(EventId)
                print 'Value Error!'
                print id
                print id.time(),id.fiducials()
                continue
            if detarr is None:
                print '*** failed to get detarr'
                continue
            if not hasattr(d,'class0'):
                d.class0 = Stats(detarr,env.experiment(),evt.run(),d.detname)
                d.class1 = Stats(detarr,env.experiment(),evt.run(),d.detname)
            else:
                if hits[i]:
                    print 'hit'
                    d.class1.update(detarr)
                else:
                    d.class0.update(detarr)
        nevent+=1

for d in detlist:
    d.class0.store('class0')
    d.class1.store('class1')
