#! /bin/bash
fd=/scratch/cnt0027/lmd1167/aspigaplaneto/2018_spiga_LES4InSight/
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS120_W10/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS120_W10.nc
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS120_W20/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS120_W20.nc
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS30_W10/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS30_W10.nc
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS30_W20/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS30_W20.nc
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS300_W10/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS300_W10.nc
#ncrcat -v PSFC -O -d Time,,,50 \
#  ${fd}/25M_LS300_W20/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
#  25M_LS300_W20.nc
ncrcat -v PSFC -O -d Time,,,50 \
  ${fd}/25M_LS300_W10/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
  25M_LS0_W10.nc
ncrcat -v PSFC -O -d Time,,,50 \
  ${fd}/25M_LS300_W20/wrfout_d01_9999-01-01_0[4-7]*:??:?? \
  25M_LS0_W20.nc

