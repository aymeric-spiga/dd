#! /usr/bin/env python
import matplotlib.pyplot as mpl
import ppplot
import numpy as np
import scipy.optimize as sciopt

### SAVED
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc_f18.ncm2_1.txt" ; drop = False ; typefit=1
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc_f18.ncm2_1.txt" ; drop = True ; typefit=2
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc_f18.ncm1_1.txt" ; drop = False ; typefit=1
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc_f18.ncm1_1.txt" ; drop = True ; typefit=1
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc.LMD_LES_MARS.160564.ncm2_1.txt" ; drop = False ; typefit=1 #bof. tous les 10 est mieux.
#
###########################################################
##namefile = "/home/aymeric/Big_Data/LES_dd/psfc_f18.ncm1_1.txt"
##namefile = "/home/aymeric/Big_Data/LES_dd/psfc_f18.ncm2_1.txt"
##namefile = "/home/aymeric/Big_Data/LES_dd/psfc.LMD_LES_MARS.160564.ncm1_1.txt"
##namefile = "/home/aymeric/Big_Data/LES_dd/sav/psfc.LMD_LES_MARS.160564.ncm2_1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/press_ustm_exomars.nc1.txt"
#namefile = "/home/aymeric/Big_Data/LES_dd/sav/press_ustm_exomars.ncm2_1.txt"
##namefile = "/home/aymeric/Big_Data/LES_dd/psfc_oldinsight100m.ncm1_1.txt"




case = "188324p"
case = "191798"
case = "160564p"
case = "156487"
case = "2007p"
case = "13526p"
case = "172097"
namefile = "/planeto/aslmd/LESdata/"+case+".ncm1_1.txt"



drop = False ; typefit=1
#drop = True ; typefit=1
#drop = True ; typefit=2
##########################################################
## define bins (we expect fit do not depend too much on this -- to be checked)
nbins = 7 ; www = 3.
#nbins = 10 ; www = 2.5
#nbins = 15 ; www = 2. 
##nbins = 20 ; www = 1.75
##nbins = 30 ; www = 1.5
##nbins = 50 ; www = 1.2
##########################################################
limrest = 4. # restrict limit (multiple of dx)
#limrest = 2.
#limrest = 3.
#limrest = 0.
##########################################################

# -- functions for fitting
# -- http://wiki.scipy.org/Cookbook/FittingData
# power law
def fitfunc(x,a,b): return a*(x**(-b))
# power law + constant
def fitfunc3(x,a,b,c): return a*(x**(-b)) + c * (x**(-0.5))
# exponential law
def fitfunc2(x,a,b): return a*np.exp(-b*x)

# load data
data = np.loadtxt(namefile,delimiter=";")
t = data[:,0] ; s = data[:,1] ; d = data[:,2]
i = data[:,3] ; j = data[:,4] # int?

# a way to guess resolution
dx = np.min(s)/2. ; print "resolution: ", dx

# choose a variable
if drop: var = d
else: var = s

# restrictions
restrict = (s > 0) # initialization (True everywhere)
restrict = restrict*(s >= limrest*dx) # remove lowest sizes (detection limit) 
#restrict = restrict*(np.abs(i-j) <= 6.*dx) # condition sur i,j (width,height)
#restrict = restrict*(d > 0.9) # limit on drop for size (casse la power law? seulement si keep smaller devils)
#restrict = restrict*(d > 0.5)
var = var[restrict]

## define bins
zebins = [np.min(var)]
for i in range(0,nbins):  zebins.append(zebins[i]*(www**0.5))
zebins = np.array(zebins)
middle = 0.5*(zebins[1:] + zebins[:-1])
binwidth = zebins[1:] - zebins[:-1]

# plot histogram
yeah = mpl.hist(var,log=True,bins=zebins,normed=True,color='white') ; mpl.xscale('log')
#yeah = mpl.hist(var,bins=zebins)#,normed=True)

# print info
print "min %5.2e // max %5.2e" % (np.min(var),np.max(var))
for iii in range(len(zebins)-1):
    print "%5.2e in [%5.0f %5.0f] %5.0f" % (yeah[0][iii],zebins[iii],zebins[iii+1],np.abs(zebins[iii+1]-zebins[iii]))

# fitting function to superimpose
if typefit == 1: xx = sciopt.curve_fit(fitfunc, middle, yeah[0])
elif typefit == 2: xx = sciopt.curve_fit(fitfunc2, middle, yeah[0])
elif typefit == 3: xx = sciopt.curve_fit(fitfunc3, middle, yeah[0])
print "exponent",xx[0][1]
print "variance %",100.*xx[1][1][1]/xx[0][1]

# plot obtained fit along with actual points
if typefit == 1: func = fitfunc(middle,xx[0][0],xx[0][1])
elif typefit == 2: func = fitfunc2(middle,xx[0][0],xx[0][1])
elif typefit == 3: func = fitfunc3(middle,xx[0][0],xx[0][1],xx[0][2]) 
func[func<0]=np.nan
ind = yeah[0]>0
mpl.plot(middle[ind],yeah[0][ind],'k.')
mpl.plot(middle[ind],func[ind],'r-')
mpl.plot(middle[ind],func[ind],'r.')

# print fit vs. actual populations
total = var.shape[0]
fit = func*binwidth*total
pop = yeah[0]*binwidth*total
for iii in range(len(fit)):
    if fit[iii] > 0.99 or pop[iii] != 0:
        print "fit %5.0f actual %5.0f" % (fit[iii],pop[iii])

# plot settings
ax = mpl.gca() ; ax.set_xbound(lower=np.min(zebins),upper=np.max(zebins))
if drop: mpl.xlabel("Pressure drop (Pa)") ; ax.set_xbound(lower=1.e-1,upper=10.)
else: mpl.xlabel("Vortex size (m)") ; ax.set_xbound(lower=10.,upper=5000.)
mpl.ylabel('Population density $n / N w_{bin}$')
mpl.title("Statistics on "+str(total)+" detected vortices")

# show plot and end
mpl.show()
exit()

################################################################################
################################################################################
################################################################################

## UNCOMMENT
##import plfit,plpva,plplot

#[alpha, xmin, L] = plfit.plfit(var,'xmin',0.3)
[alpha, xmin, L] = plfit.plfit(var)
print alpha,xmin

#a = plpva.plpva(var,0.75,'xmin',0.75)
#print a


h = plplot.plplot(var,xmin,alpha)

ppplot.save(mode="png")

#mpl.loglog(h[0], h[1], 'k--',linewidth=2)

#mpl.show()



