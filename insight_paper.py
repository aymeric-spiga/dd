#! /usr/bin/env python

from finddd import finddd
from ddstat import statdd

plotplot=False
timelist=None

ltstart=11.
dt_out=50.
ffftab = ["25M_LS120_W10.nc", \
          "25M_LS120_W20.nc", \
          "25M_LS30_W10.nc", \
          "25M_LS30_W20.nc", \
          "25M_LS300_W10.nc", \
          "25M_LS300_W20.nc"]

ffftab = ["25M_LS300_W10.nc", \
          "25M_LS300_W20.nc"]

ffftab = ["25M_LS300_W20.nc"]

ffftab = ["25M_LS0_W10.nc", \
          "25M_LS0_W20.nc"]

ffftab = ["25M_LS0_W20.nc"]



## for checks
#ffftab = ["25M_LS300_W10_extract.nc"]
#timelist = range(1)
#plotplot = True


for fff in ffftab:
    finddd(fff,dt_out=dt_out,lt_start=ltstart,method=1,save=True,plotplot=plotplot,timelist=timelist,prescribe=True)
    statdd(fff+"m1_2")


