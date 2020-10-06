#! /usr/bin/env python
import matplotlib.pyplot as mpl
import ppplot
import numpy as np
import scipy.optimize as sciopt
from matplotlib.ticker import FormatStrFormatter

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
def histodd(namefile,folder='./',drop=False,typefit=1,nbins=12,limrest=4,limtime=None,limdrop=0.3,addtitle="",limwind=None):

    mpl.rcParams['lines.markersize'] = 10

    # width adapted to number of bins
    widthbin = {7:3.,10:2.5,12:2.2,15:2.0,20:1.75,30:1.5,50:1.2,100:1.1,200:1.05}
    www = widthbin[nbins]
 
    # load data
    data = np.loadtxt(namefile+"_1.txt",delimiter=";")
    t = data[:,0] ; s = data[:,1] ; d = data[:,2]
    i = data[:,3] ; j = data[:,4] ; v = data[:,5]
    
    # a way to guess resolution 
    # [smallest radius is sqrt((limnmesh)*dx*dx) with limnmest defined in gethalo]
    dx = np.ceil(np.min(s)/np.sqrt(4))
    
    # choose a variable
    if drop: var = d
    else: var = s
    
    # restrictions
    restrict = (s > 0) # initialization (True everywhere)
    #restrict = restrict*(np.abs(i-j) <= 6.*dx) # condition sur i,j (width,height) pour ~round shape
                                                # -- does not really change the results anyway
    if limrest is not None: restrict = restrict*(s >= limrest*dx) # remove lowest sizes (detection limit) 
    if limdrop is not None: restrict = restrict*(d >= limdrop) # remove lowest drop
    if limtime is not None: restrict = restrict*(t <= limtime) # remove later local times
    if limwind is not None: restrict = restrict*(v > limwind) # remove lowest velocity (often false positives)
    poum=False
    #poum=True
    if poum:
      out = var[np.logical_not(restrict)]
      outi = i[np.logical_not(restrict)]
      outj = j[np.logical_not(restrict)]
      ##print out, outi, outj
      #print np.min(out),np.max(out),np.median(out)
      #print 100.*out.size/var.size
    ###
    var2 = var[restrict]
    total = var2.shape[0]
    #titi = addtitle+"LES "+str(int(dx))+"m. N="+str(total)+" detected vortices"
    titi = addtitle+"N="+str(total)+" detected drops"

    ## define bins
    zebins = [np.min(var2)]
    for i in range(0,nbins):  zebins.append(zebins[i]*(www**0.5))
    zebins = np.array(zebins)
    middle = 0.5*(zebins[1:] + zebins[:-1])
    binwidth = zebins[1:] - zebins[:-1]
    #if (not drop) and (np.round(binwidth[0]) < dx): 
    #  print "too much bins make it binwidth < dx. decrease number of bins."
    #  print dx, binwidth
    #  exit()
    minfunc = 1./(binwidth*total) # to show in histo the minimum population: 1

    #print zebins
 
    # plot histogram
    mpl.figure(figsize=(18,6)) ; mpl.xscale('log')
    yeah = mpl.hist(var2,log=True,bins=zebins,normed=True,color='white')
    #yeah = mpl.hist(var2,bins=zebins,normed=True,color='white')
    #yeah = mpl.hist(var2,bins=zebins,color='white')

    # add error bars
    incertitude_nvortex = 5. #10.
    if not drop:
      err = incertitude_nvortex/(yeah[0]*binwidth*total+0.0001)
      ind = err < 1000. ; err = err[ind] ; y = yeah[0][ind] ; x = middle[ind]
      err = y*err
      mpl.errorbar(x,y,yerr=[err,err],fmt='k.')
  
    ## print info
    #print "min %5.2e // max %5.2e" % (np.min(var2),np.max(var2))
    #for iii in range(len(zebins)-1):
    #    this = zebins[iii] ; next = zebins[iii+1]
    #    print "%5.2e in [%5.0f %5.0f] %5.0f" % (yeah[0][iii],this,next,np.abs(next-this))
    
    ### fitting function to superimpose
    ydata = yeah[0] ; xdata = middle
    if typefit == 1: xx = sciopt.curve_fit(fitfunc, xdata, ydata)
    elif typefit == 2: xx = sciopt.curve_fit(fitfunc2, xdata, ydata)
    elif typefit == 3: xx = sciopt.curve_fit(fitfunc3, xdata, ydata)
    #print "exponent",xx[0][1] #,"variance %",100.*xx[1][1][1]/xx[0][1]
    pcov = xx[1]
    perr = np.sqrt(np.diag(pcov))
    #print "one-sigma error on parameter", perr[1]
    #print "------> interval for parameter", np.round(xx[0][1]-3.*perr[1],2), np.round(xx[0][1]+3.*perr[1],2)
    
    ## label
    lablab = r"exponent %3.1f $\pm$ %3.1f"%(xx[0][1],3.*perr[1])
    
    # plot obtained fit along with actual points
    if typefit == 1: func = fitfunc(middle,xx[0][0],xx[0][1])
    elif typefit == 2: func = fitfunc2(middle,xx[0][0],xx[0][1])
    elif typefit == 3: func = fitfunc3(middle,xx[0][0],xx[0][1],xx[0][2]) 
    func[func<0]=np.nan
    ind = yeah[0]>0
    mpl.bar(yeah[1][:-1],yeah[0],width=binwidth,color="w",align="edge",linewidth=1,edgecolor="k")
    mpl.plot(middle[ind],yeah[0][ind],'ks')           # data
    mpl.plot(middle[ind],func[ind],'r-',label=lablab) # fit
    #mpl.plot(middle[ind],func[ind],'r.')              # fit points
    mpl.plot(middle[ind],10*minfunc[ind],'g:',label="10-element limit")        # 10 limit
    mpl.plot(middle[ind],minfunc[ind],'r:',label="1-element limit")           # 1 limit
    divbound = 1.35 # larger for larger spaces around bins. 1.5 pretty good.
    minbin = np.min(middle[ind])/divbound # for plotting window
    maxbin = np.max(middle[ind])*divbound # for plotting window
    minval = np.min(minfunc[ind])/divbound
 
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
            #print "fit %4.0f real %4.0f pc %3.0f ipc %3.0f" % (fit[iii],pop[iii],yorgl,yargl)

    ## additional conditional histogram in blue
    ##var3 = var[restrict*(d > 0.5)]
    #var3 = var2 ; zebins = zebins[yeah[0]*binwidth*total >= incertitude_nvortex]
    #mpl.hist(var3,log=True,bins=zebins,normed=True,color='blue')
    
    # plot settings
    ax = mpl.gca()
    ax.set_xbound(lower=minbin,upper=maxbin)
    ax.set_ybound(lower=minval)
    if drop: 
        mpl.xlabel("Pressure drop (Pa)")
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1g"))
        ##ax.set_xbound(lower=1.e-1,upper=10.)
    else: 
        mpl.xlabel("Vortex size (m)")
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
        ##ax.set_xbound(lower=30.,upper=600.)
        ##ax.set_ybound(lower=5.e-5,upper=5.e-2)

    mpl.ylabel('Population density $n / N w_{bin}$')
    #mpl.title(titi+'\n '+tata)
    mpl.title(titi)
    mpl.legend(loc="upper right",fancybox=True)

    ## show plot and end
    #mpl.show()
    if drop:
      mpl.savefig(folder+"/"+"ddhistdrop_"+namefile+"_bin"+str(nbins)+"_limdrop%i"%(limdrop*10.)+".pdf")
    else:
      mpl.savefig(folder+"/"+"ddhist_"+namefile+"_bin"+str(nbins)+".pdf")
    mpl.close()

################################################################################
### FDROPSIZE
################################################################################
def fdropsize(namefile):
    # load data
    data = np.loadtxt(namefile+"_1.txt",delimiter=";")
    s = data[:,1] ; d = data[:,2]
    # plot drop = f(size)
    restrict = (s > 0)
    dropl = ppplot.plot1d() ; dropl.f = d[restrict] ; dropl.x = s[restrict]
    dropl.linestyle,dropl.marker,dropl.color,dropl.fmt = '','.','r',"%.1f"
    dropl.xlabel = "Vortex size (m)" ; dropl.ylabel = "Pressure drop (Pa)"
    #dropl.logx= True ; dropl.logy = True
    dropl.makeshow()
################################################################################

################################################################################
### FDROPWIND
################################################################################
def fdropwind(namefile):
    # load data
    data = np.loadtxt(namefile+"_1.txt",delimiter=";")
    d = data[:,2] ; v = data[:,5]
    # plot drop = f(size)
    restrict = (d > 0)
    dropl = ppplot.plot1d() ; dropl.f = d[restrict] ; dropl.x = v[restrict]
    dropl.linestyle,dropl.marker,dropl.color,dropl.fmt = '','.','r',"%.1f"
    dropl.xlabel = r"Friction velocity (m s$^{-1}$)" ; dropl.ylabel = "Pressure drop (Pa)"
    #dropl.logx= True ; dropl.logy = True
    dropl.makeshow()
################################################################################

################################################################################
### FSIZEWIND
################################################################################
def fsizewind(namefile):
    # load data
    data = np.loadtxt(namefile+"_1.txt",delimiter=";")
    s = data[:,1] ; v = data[:,5]
    # plot drop = f(size)
    restrict = (s > 0)
    dropl = ppplot.plot1d() ; dropl.f = v[restrict] ; dropl.x = s[restrict]
    dropl.linestyle,dropl.marker,dropl.color,dropl.fmt = '','.','r',"%.1f"
    dropl.xlabel = "Vortex size (m)" ; dropl.ylabel = r"Friction velocity (m s$^{-1}$)"
    #dropl.logx= True ; dropl.logy = True
    dropl.makeshow()
################################################################################

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

