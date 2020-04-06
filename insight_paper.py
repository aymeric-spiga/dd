#! /usr/bin/env python

from finddd import finddd
from ddstat import statdd

ffftab = ["25M_LS120_W10.nc", \
          "25M_LS120_W20.nc", \
          "25M_LS30_W10.nc", \
          "25M_LS30_W20.nc", \
          "25M_LS300_W10.nc", \
          "25M_LS300_W20.nc"]

for fff in ffftab:
    finddd(fff,dt_out=20.,lt_start=14.,method=1,save=True,plotplot=False)
    statdd(fff+"m1_2")

