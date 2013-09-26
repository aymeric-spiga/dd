#! /usr/bin/env python
from finddd import finddd

####################
####################
#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",range(0,591,1))
#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",range(0,591,1),method=2)
####################
####################


#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",range(0,591,10),method=2)
#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",range(0,591,10),method=2,plotplot=True)
#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",[300.,350.,400.],method=2,plotplot=True)
#### checks
#plotplot = True ; save = False
#timelist = [270.,302.,305.,352.,358.,412.,423.] #joli
#timelist = [300.,310.,320.]
#timelist = [210.,250.,280.,400.]
#finddd("/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",range(350,450,5),method=2,plotplot=True,save=False)

####################
####################
finddd("/home/aymeric/Big_Data/LES_dd/press_ustm_exomars.nc",range(0,54,1),dx=12.,halolim=100,method=2)
#timelist = [28.,48.,49.]
#plotplot = True
#save = False
####################
####################

####################
####################
#finddd("/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.nc",\
#           range(0,336,1),dx=10.,lt_start=9.,halolim=150,method=2)
###halolim 50 idem
####################
####################



###############################################################################
#dx = 100. # horizontal resolution
#lt_out = 10. # frequency of outputs (s) --- wrong?????
#lt_start = 6. # starting local time 
#filefile = "/home/aymeric/Big_Data/LES_dd/psfc_oldinsight100m.nc"
#halolim = 20 # limit for halos
#timelist = range(0,4440,1) # times in file
#timelist = range(0,4440,100) # times in file
#method = 1 # method used to detect vortices
#plotplot = False # plot or not
#save = True # save results or not
###############################################################################


##namefile = "/home/aymeric/Big_Data/LES_dd/press_ustm_exomars.nc2.txt"
##dt_out = 50.
##lt_start = 8.


