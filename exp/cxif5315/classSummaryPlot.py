import numpy as np
from psana import *
from Detector.PyDetector import PyDetector

import sys
filename=sys.argv[1]
plottype=sys.argv[2]

print filename,plottype

from psmon import publish
from psmon.plots import MultiPlot,XYPlot,Image

ds = DataSource('exp=cxif5315:run=169')
env = ds.env()
evt = ds.events().next()

data = np.load(filename)
src = Source('DetInfo(CxiDs2.0:Cspad.0)')
det = PyDetector(src,env)

img = det.image(evt,data[plottype])

publish.local = True
plt = Image(0,filename,img) 
publish.send('image', plt)
