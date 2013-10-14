#! /usr/bin/env python
from finddd import finddd

finddd("/planeto/aslmd/LESdata/case_A_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)

exit()

# LES_DUST_DEVIL_OLDPHYSok
finddd("/planeto/aslmd/LESdata/188324p.nc",dt_out=50.,lt_start=9.)
# LES
finddd("/planeto/aslmd/LESdata/160564p.nc",dt_out=50.,lt_start=9.)
finddd("/planeto/aslmd/LESdata/156487.nc",dt_out=50.,lt_start=8.)
# LES_INSIGHT
finddd("/planeto/aslmd/LESdata/2007p.nc",dt_out=10.,lt_start=6.)
finddd("/planeto/aslmd/LESdata/13526p.nc",dt_out=10.,lt_start=6.)
# EXOMARSshear
finddd("/planeto/aslmd/LESdata/172097.nc",dt_out=50.,lt_start=13.)
# QJ cases
finddd("/planeto/aslmd/LESdata/case_A_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)
finddd("/planeto/aslmd/LESdata/case_B_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)
finddd("/planeto/aslmd/LESdata/case_C_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)
finddd("/planeto/aslmd/LESdata/case_I_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)
finddd("/planeto/aslmd/LESdata/case_HIGH_50m_145_145_201_12km.nc",dt_out=100.,lt_start=8.)

#####################################################################
#####################################################################
#####################################################################


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
#finddd("/planeto/aslmd/LESdata/press_ustm_exomars.nc",range(0,54,1),halolim=100,method=2)
#timelist = [28.,48.,49.]
#plotplot = True
#save = False
####################
####################

####################
####################
#finddd("/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.nc",\
#           range(0,336,1),lt_start=9.,halolim=150,method=2)
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


