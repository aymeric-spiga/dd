#! /usr/bin/env python
import matplotlib.pyplot as mpl
import ppplot
import numpy as np
import scipy.optimize as sciopt

# -- functions for fitting
# -- http://wiki.scipy.org/Cookbook/FittingData
# power law
def fitfunc(x,a,b): return a*(x**(-b))
# power law + constant
def fitfunc3(x,a,b,c): return a*(x**(-b)) + c * (x**(-0.5))
# exponential law
def fitfunc2(x,a,b): return a*np.exp(-b*x)

################################################################################
### HISTODD
### namefile: text file to be analyzed
### drop: look at size (False) or pressure drop (True)
### typefit: use fitfunc above 1, 2, or 3
### nbins: bins (we expect fit do not depend too much on this -- to be checked)
### limrest: restrict limit (multiple of dx) -- greater equal
### limtime: use only data earlier than this local time 
### limdrop: use only data with deeper drop than value
### limwind: use only data with friction velocity larger than this value
### addtitle: add a name for the examined case
################################################################################
def histodd(namefile,drop=False,typefit=1,nbins=12,limrest=4,limtime=None,limdrop=0.3,addtitle="",limwind=None):

    # width adapted to number of bins
    widthbin = {7:3.,10:2.5,12:2.2,15:2.0,20:1.75,30:1.5,50:1.2}
    www = widthbin[nbins]
 
    # load data
    data = np.loadtxt(namefile,delimiter=";")
    t = data[:,0] ; s = data[:,1] ; d = data[:,2]
    i = data[:,3] ; j = data[:,4]
    
    # a way to guess resolution [smallest radius is sqrt(3*dx*dx)]
    dx = np.ceil(np.min(s)/np.sqrt(3))
    
    # choose a variable
    if drop: var = d
    else: var = s
    
    # restrictions
    restrict = (s > 0) # initialization (True everywhere)
    restrict = restrict*(s >= limrest*dx) # remove lowest sizes (detection limit) 
    restrict = restrict*(np.abs(i-j) <= 6.*dx) # condition sur i,j (width,height) pour ~round shape
    #restrict = restrict*(d > 0.9) # limit on drop for size 
                                   # (casse la power law? seulement si keep smaller devils)
    #restrict = restrict*(d > 0.5)
    #restrict = restrict*(d > 0.3)
    poum=False
    if poum:
      out = var[np.logical_not(restrict)]
      outi = i[np.logical_not(restrict)]
      outj = j[np.logical_not(restrict)]
      print out, outi, outj
      print out.size
    var = var[restrict]
    
    ## define bins
    zebins = [np.min(var)]
    for i in range(0,nbins):  zebins.append(zebins[i]*(www**0.5))
    zebins = np.array(zebins)
    middle = 0.5*(zebins[1:] + zebins[:-1])
    binwidth = zebins[1:] - zebins[:-1]
    total = var.shape[0]
    
    # plot histogram
    yeah = mpl.hist(var,log=True,bins=zebins,normed=True,color='white') ; mpl.xscale('log')
    #yeah = mpl.hist(var,bins=zebins,normed=True,color='white')
    
    # add error bars
    incertitude_nvortex = 10.
    if not drop:
      err = incertitude_nvortex/(yeah[0]*binwidth*total+0.0001)
      ind = err < 1000. ; err = err[ind] ; y = yeah[0][ind] ; x = middle[ind]
      err = y*err
      mpl.errorbar(x,y,yerr=[err,err],fmt='k.')
    
    ## print info
    #print "min %5.2e // max %5.2e" % (np.min(var),np.max(var))
    #for iii in range(len(zebins)-1):
    #    this = zebins[iii] ; next = zebins[iii+1]
    #    print "%5.2e in [%5.0f %5.0f] %5.0f" % (yeah[0][iii],this,next,np.abs(next-this))
    
    # fitting function to superimpose
    if typefit == 1: xx = sciopt.curve_fit(fitfunc, middle, yeah[0])
    elif typefit == 2: xx = sciopt.curve_fit(fitfunc2, middle, yeah[0])
    elif typefit == 3: xx = sciopt.curve_fit(fitfunc3, middle, yeah[0])
    print "exponent",xx[0][1],"variance %",100.*xx[1][1][1]/xx[0][1]
    
    # label
    lablab = r"$\alpha=$%4.1f"%(xx[0][1])
    titi = addtitle+"LES "+str(int(dx))+"m. "+str(total)+" detected vortices"
    
    # plot obtained fit along with actual points
    if typefit == 1: func = fitfunc(middle,xx[0][0],xx[0][1])
    elif typefit == 2: func = fitfunc2(middle,xx[0][0],xx[0][1])
    elif typefit == 3: func = fitfunc3(middle,xx[0][0],xx[0][1],xx[0][2]) 
    func[func<0]=np.nan
    ind = yeah[0]>0
    mpl.plot(middle[ind],yeah[0][ind],'k.')
    mpl.plot(middle[ind],func[ind],'r-',label=lablab)
    mpl.plot(middle[ind],func[ind],'r.')
    
    # display bin limits
    tata = '[bins:'
    yorgl = ind[np.where(ind)].size + 1 # to add the last limit
    for el in zebins[0:yorgl]:
       tata = tata + "%0.f "%(el)
    tata=tata+']'
    
    # print fit vs. actual populations
    fit = func*binwidth*total
    pop = yeah[0]*binwidth*total
    for iii in range(len(fit)):
        if fit[iii] > 0.99 or pop[iii] != 0:
            yorgl = 100.*np.abs(fit[iii]-pop[iii])/fit[iii]
            yargl = 100.*incertitude_nvortex/fit[iii]
            print "fit %4.0f real %4.0f pc %3.0f ipc %3.0f" % (fit[iii],pop[iii],yorgl,yargl)
    
    # plot settings
    ax = mpl.gca() ; ax.set_xbound(lower=np.min(zebins),upper=np.max(zebins))
    if drop: 
        mpl.xlabel("Pressure drop (Pa)")
        ax.set_xbound(lower=1.e-1,upper=10.)
    else: 
        mpl.xlabel("Vortex size (m)")
        #ax.set_xbound(lower=30.,upper=3000.)
        ax.set_xbound(lower=30.,upper=1000.)
        ax.set_ybound(lower=5.e-6,upper=5.e-2)
    mpl.ylabel('Population density $n / N w_{bin}$')
    #mpl.title(titi+'\n '+tata)
    mpl.title(titi)
    mpl.legend(loc="upper right",fancybox=True)
        
    # show plot and end
    mpl.show()

#################################################################################
#################################################################################
#################################################################################
#
### UNCOMMENT
###import plfit,plpva,plplot
#
##[alpha, xmin, L] = plfit.plfit(var,'xmin',0.3)
#[alpha, xmin, L] = plfit.plfit(var)
#print alpha,xmin
#
##a = plpva.plpva(var,0.75,'xmin',0.75)
##print a
#
#
#h = plplot.plplot(var,xmin,alpha)
#
#ppplot.save(mode="png")
#
##mpl.loglog(h[0], h[1], 'k--',linewidth=2)
#
##mpl.show()
#
#

