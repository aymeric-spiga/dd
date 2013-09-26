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

# read data np.loadtxt?
## import pylab as plb data = plb.loadtxt('data.dat') x = data[:,0] y= data[:,1]
f = open(namefile,'r')
t = [] ; s = [] ; d = [] ; i = [] ; j = []
for line in f:
    ct,cs,cd,ci,cj = line.strip().split(';')
    t.append(float(ct))
    s.append(float(cs)) ; d.append(float(cd))
    i.append(int(ci)) ; j.append(int(cj))
f.close()
# ensure numpy array
t = np.array(t) ; s = np.array(s) ; d = np.array(d) ; i = np.array(i) ; j = np.array(j)

drop = ppplot.plot1d() ; drop.field = d ; drop.absc = s
drop.lstyle = '' ; drop.marker = '.' ; drop.color = 'r' ; drop.fmt = "%.1f"
drop.xlabel = "Vortex size (m)" ; drop.ylabel = "Pressure drop (Pa)"
drop.make() ; ppplot.save()

