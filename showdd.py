#! /usr/bin/env python

from ppclass import pp
import numpy as np
from scipy import ndimage

dx = 50. ; filefile = "/home/aymeric/Big_Data/LES_dd/psfc_f18.nc"
lim = 4.
#lim = -3.0
#lim = -3.2

#dx = 10. ; filefile = "/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.nc"

## getplot pressure field
## ... generic settings
psfc = pp()
psfc.file = filefile
psfc.var = "PSFC"
psfc.xcoeff = dx/1000.
psfc.ycoeff = dx/1000.
psfc.xlabel = "x distance (km)"
psfc.ylabel = "y distance (km)"
psfc.title = "Surface pressure $P_s$"
#psfc.vmin = 482.5
#psfc.vmax = 483.5
psfc.div = 20
#psfc.out = "png" ; psfc.includedate = False ; psfc.folder = "dd/"

## calculate std
std = np.std(psfc.getf())

## calculate mean pressure
## ... same results whether or not time average
psfcmean = pp()
psfcmean << psfc
psfcmean.x="0,10000"
psfcmean.y="0,10000"
psfcmean.t="0,10000"
psfcmean.getplot()

#for t in [300]:
#for t in range(0,500,50):
#for t in [390,395,400,405,410,415]:
#for t in [178,180,182]:

for t in range(265,275,1):

   ## plot pressure field
   psfc.t = t
   psfc.filename = "dd"+str(int(t))
   psfc.getplot(extraplot=1)

   ## calculate anomaly normalized by std
   anomaly = (psfc - psfcmean) / std
      
   ## plot anomaly
   anomaly.vmax = lim
   anomaly.vmin = -4.0
   anomaly.div = 32
   anomaly.colorb = "gist_earth_r"
   anomaly.colorb = "gist_stern"
   anomaly.plotin = psfc
   anomaly.title = "Anomaly $P_s - <P_s>^{xy}$"
   anomaly.units = "$\sigma$"
   #anomaly.plot()

   ## sobel transform
   field = anomaly.request[0][0][0][0][0][0].field
   sx = ndimage.sobel(field, axis=0, mode='constant') 
   sy = ndimage.sobel(field, axis=1, mode='constant')
   tmp = np.hypot(sx, sy)

   ## more prominent circles if we mutiply by field
   tmp = tmp*field
   #lim = -12.
   #tmp[tmp>lim]=0.
   #tmp[tmp<lim]=1.

   #tmp = np.exp(-tmp)   

   anomaly.request[0][0][0][0][0][0].field = tmp
   anomaly.title = "Sobel transform" 
   anomaly.units = "dimless" ; anomaly.vmax = None ; anomaly.vmin = None ; anomaly.div = 10
   anomaly.colorb = "binary_r"
   #anomaly.vmax = np.max(tmp)
   anomaly.plot()

   ### laplace transform
   #anomaly.request[0][0][0][0][0][0].field=ndimage.laplace(anomaly.request[0][0][0][0][0][0].field)
   #anomaly.title = "Sobel + Laplace transform"
   #anomaly.units = "dimless" ; anomaly.vmax = None ; anomaly.vmin = None ; anomaly.div = 10
   #anomaly.plot()
