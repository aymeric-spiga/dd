#! /usr/bin/env python
import matplotlib.pyplot as mpl
import ppplot
import numpy as np

##########################################################
namefile = "/home/aymeric/Big_Data/LES_dd/psfc_f18.ncm1_1.txt"
namefile = "/home/aymeric/Big_Data/LES_dd/psfc_f18.ncm2_1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.ncm1_1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.ncm2_1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/press_ustm_exomars.nc1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/psfc_oldinsight100m.ncm1_1.txt"
##########################################################

#case = "188324p"
#case = "191798"
#case = "160564p"
#case = "156487"
#case = "2007p"
#case = "13526p"
#case = "172097"


namefile = "/planeto/aslmd/LESdata/"+case+".ncm1_1.txt"

# load data
data = np.loadtxt(namefile,delimiter=";")
t = data[:,0] ; s = data[:,1] ; d = data[:,2]
i = data[:,3] ; j = data[:,4] # int?

drop = ppplot.plot1d() ; drop.f = d ; drop.x = s
drop.linestyle = '' ; drop.marker = '.' ; drop.color = 'r' ; drop.fmt = "%.1f"
drop.xlabel = "Vortex size (m)" ; drop.ylabel = "Pressure drop (Pa)"
drop.makeshow()

