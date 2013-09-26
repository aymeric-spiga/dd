#! /usr/bin/env python

def detsize( xx, res=1, thres=3, loga=False ):
    import numpy as np
    import math
    size = []
    sizecalc = 1
    diff = np.asarray( np.roll(xx,-1) - xx )
    for i in diff:
        if abs(i) > 1:
            if sizecalc >= thres: 
                if loga: addthis = math.log(sizecalc*res)
                else:    addthis = sizecalc*res
                size.append(addthis)
            sizecalc = 1
        else:
            sizecalc += 1
    return size

def getsize(filename):

    import numpy as np
    from scipy.ndimage.measurements import minimum_position
    from scipy import ndimage
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    import myplot as myp
   
    ### LOAD NETCDF DATA
    nc = Dataset(filename) 
    psfc = nc.variables["PSFC"]
    print "yeah"
    
    ### LOOP on TIME
    ### NB: a same event could be counted several times...
    shape = np.array(psfc).shape
    allsizesx = []
    allsizesy = []
    depression = []
    stride = 1 #5
    stride = 20
    #stride = 50
    stride = 100
    start = 0
    start = stride
    for i in range(start,shape[0],stride):

        psfc2d = np.array ( psfc [ i, : , : ] )
    
        ############### CRITERION
        ave = np.mean(psfc2d,dtype=np.float64)  ## dtype otherwise inaccuracy

        #limdp = -0.2 ## on en loupe pas mal
        #where = np.where(psfc2d - ave < limdp)  ## comme le papier Phoenix
    
        std = np.std(psfc2d,dtype=np.float64)   ## dtype otherwise inaccuracy
        fac = 4.   ## how many sigmas. not too low, otherwise vortices are not caught. 4 good choice.
                   ## 2.5 clearly too low, 3.5 not too bad, 4 probably good
        fac = 3.5
        fac = 3.2
        ##fac = 2.5
        #fac = 3. ## final choice
        #fac = 2.5
        lim = ave - fac*std
        where = np.where(psfc2d < lim)
        ############### END CRITERION

        depression = np.append(depression,np.ravel(psfc2d[where])-ave)
   
        ## lab is 0 or 1
        lab = np.zeros(np.array(psfc2d).shape) ## points to be treated by the minimum_position routine
        lab[where] = 1.  ## do not treat points close to 'mean' (background) pressure
  
        xx = []
        yy = []
        while 1 in lab:
            p = minimum_position(psfc2d,labels=lab)
            lab[p] = 0 ## once a minimum has been found in a grid point, do not search here again.
            if p[0] not in xx: xx.append(p[0]) ## if x coordinate not yet in the list add it
            if p[1] not in yy: yy.append(p[1]) ## if y coordinate not yet in the list add it
        xx.sort()
        yy.sort()
        ### now xx and yy are sorted arrays containing grid points with pressure minimum
       
        ######## DETERMINE SIZE OF STRUCTURES
        ######## this is rather brute-force...
        sizex = detsize( xx, res = 10, loga=False, thres=2 )
        sizey = detsize( yy, res = 10, loga=False, thres=2 )
        #sizex = detsize( xx, res = 10, loga=False, thres=3 )
        #sizey = detsize( yy, res = 10, loga=False, thres=3 )
        sizex = detsize( xx, res = 15, loga=False, thres=2 )
        sizey = detsize( yy, res = 15, loga=False, thres=2 )
        ###
        print sizex, sizey
        #if ( mym.max(sizex) > mym.max(sizey) ): sizey = sizex  ### un peu limite dans certains cas
        if (len(sizex) > len(sizey))     : sizey = sizex        ### plus fidele mais petit souci lorsque PBC
        elif (len(sizex) == len(sizey))  : 
          if ( mym.max(sizex) > mym.max(sizey) ): sizey = sizex
          else                                  : sizex = sizey
        else                            : sizex = sizey
        allsizesx = np.append(allsizesx,sizex)
        allsizesy = np.append(allsizesy,sizey)
        print i, ' on ', shape[0], ' caught ', len(sizex), ' vortices ', sizex
        ########
  
    allsizesx.sort()
    allsizesy.sort()
    
    return allsizesx, allsizesy, depression

#########################################################################
#########################################################################

import matplotlib.pyplot as plt
import pickle
import numpy as np
import matplotlib.mlab as mlab
import mymath as mym
import myplot as myp

import plfit
import plplot
import randht
import plpva

save = True
#save = False
pression = False
pression = True

filename = "/home/aymeric/Big_Data/psfc_f18.nc"

if save:
    ### getsize
    allsizesx, allsizesy, depression = getsize(filename)
    ### sauvegarde texte pour inspection
    mym.writeascii(allsizesx,'allsizex.txt')
    mym.writeascii(allsizesy,'allsizey.txt')
    mym.writeascii(depression,'alldepression.txt')
    ### sauvegarde binaire pour utilisation python
    myfile = open('allsizex.bin', 'wb') ; pickle.dump(allsizesx, myfile) ; myfile.close()
    myfile = open('allsizey.bin', 'wb') ; pickle.dump(allsizesy, myfile) ; myfile.close()
    myfile = open('alldepression.bin', 'wb') ; pickle.dump(depression, myfile) ; myfile.close()

### load files
myfile = open('allsizex.bin', 'r')
allsizesx = pickle.load(myfile)
myfile = open('allsizey.bin', 'r')
allsizesy = pickle.load(myfile)
myfile = open('alldepression.bin', 'r')
depression = pickle.load(myfile)
depression = np.array(abs(depression))#*1000.

### sizes
#plothist = np.append(allsizesx,allsizesy)
plothist = allsizesx
if pression: plothist = depression
plothist.sort()
print 'mean ', np.mean(plothist,dtype=np.float64)  
print 'std ', np.std(plothist,dtype=np.float64)
print 'max ', np.max(plothist)
print 'min ', np.min(plothist)
print 'len ', len(plothist)


### MAKE BINS
nbins = 100
zebins = [2.0]
#nbins = 8
#zebins = [19.0]
#nbins = 15
#zebins = [11.] ##12 non mais donne un peu la meme chose
#zebins = [20.] ##20 non car trop pres du premier
nbins = 100
zebins = [2./np.sqrt(2.)]  ## ne pas tomber sur une dizaine ronde
nbins = 200

if pression: zebins = [0.3]

for i in range(0,nbins):  zebins.append(zebins[i]*np.sqrt(2))
zebins = np.array(zebins)
#### select reasonable bins for DD
if not pression:
    zebins = zebins [ zebins > 15. ] 
    #zebins = zebins [ zebins > 20. ]
    zebins = zebins [ zebins > 25. ]
    zebins = zebins [ zebins < 1000. ]
else:
    zebins = zebins [ zebins < 10. ]
print 'corrected bins ',zebins

#### HISTOGRAM
plt.figure(1)
plt.hist(    plothist,\
             log=True,\
             bins=zebins,\
#             cumulative=-1,\
             normed=True,\
        )
plt.xscale('log')
if pression:    plt.xlabel('Pressure (Pa)')
else:           plt.xlabel('Diameter (m)')
plt.ylabel('Population (normalized)')
if pression: prefix="p"
else:        prefix=""
myp.makeplotres(prefix+"histogram",res=200,disp=False)
plt.close(1)

### COMPARED HISTOGRAMS
### --- FIT WITH POWER LAW
if pression: [alpha, xmin, L] = plfit.plfit(plothist,'xmin',0.3)
else:        [alpha, xmin, L] = plfit.plfit(plothist,'limit',20.)
print alpha,xmin

#a = plpva.plpva(plothist,0.75,'xmin',0.75)
#print a

#### DEUXIEME ROUTINE
####IL FAUT UTILISER LE DISCRET POUR LA TAILLE !!!
#if pression:   myplfit = plfit.plfit(plothist,verbose=True,xmin=0.75)
#else:          myplfit = plfit.plfit(plothist,verbose=True,xmin=20.)
#myplfit.plotppf()
#plt.show()
#exit()


#plt.figure(1)
#h = plplot.plplot(plothist,xmin,alpha)
#myp.makeplotres(prefix+"fit",res=200,disp=False)
plt.figure(2)
### --- POWER LAW (factor does not really matter)
power = (xmin/2.2)*np.array(randht.randht(10000,'powerlaw',alpha))
#power = (xmin/2.2)*np.array(randht.randht(10000,'cutoff',alpha,10.)) ##marche pas si trop grand
print 'mean ', np.mean(power,dtype=np.float64)
### --- EXPONENTIAL LAW
expo = randht.randht(10000,'exponential',1./(np.mean(power,dtype=np.float64)*1.00))
print 'mean ', np.mean(expo,dtype=np.float64)
### --- PLOT
plt.hist(    [plothist,power,expo],\
             label=['LES vortices','Power law '+'{:.1f}'.format(alpha),'Exponential law'],\
             log=True,\
             bins=zebins,\
#            cumulative=-1,\
             normed=True,\
        )
plt.legend()
plt.xscale('log')
if pression:    plt.xlabel('Pressure (Pa)')
else:           plt.xlabel('Diameter (m)')
plt.ylabel('Population (normalized)')
myp.makeplotres(prefix+"comparison",res=200,disp=False)
plt.close(2)

########################
########################
zebins = [30.,42.,60.,84.,120.,170.,240.,340.]
plothist = []
plothist = np.append(plothist,30 *np.ones(306))
plothist = np.append(plothist,42 *np.ones(58) )
plothist = np.append(plothist,60 *np.ones(66) )
plothist = np.append(plothist,84 *np.ones(41) )
plothist = np.append(plothist,120*np.ones(19) )
plothist = np.append(plothist,170*np.ones(9)  )
plothist = np.append(plothist,240*np.ones(2)  )
plothist = np.append(plothist,340*np.ones(1)  )

#zebins = [50.,71.,100.,141.,200.,282.,400.]
#plothist = []
#plothist = np.append(plothist,50. *np.ones(36))
#plothist = np.append(plothist,71. *np.ones(18) )
#plothist = np.append(plothist,100. *np.ones(12) )
#plothist = np.append(plothist,141. *np.ones(6) )
#plothist = np.append(plothist,200.*np.ones(4) )
#plothist = np.append(plothist,282.*np.ones(1)  )
#plothist = np.append(plothist,400.*np.ones(2)  )

exit()

plt.figure(3)
[alpha, xmin, L] = plfit.plfit(plothist,'xmin',30)#50.)
print alpha,xmin
#a = plpva.plpva(plothist,30,'xmin',30)
h = plplot.plplot(plothist,xmin,alpha)
plt.loglog(h[0], h[1], 'k--',linewidth=2)
plt.hist(    plothist,\
             log=True,\
             bins=zebins,\
#             cumulative=-1,\
             normed=True,\
        )
plt.xscale('log')
plt.xlabel('Pressure (micro Pa)')
plt.ylabel('Population (normalized)')
myp.makeplotres("data",res=200,disp=False)

#plt.figure(4)
#[alpha, xmin, L] = plfit.plfit(plothist,'xmin',50.) #,'xmin',30.)
#print alpha,xmin
#h = plplot.plplot(plothist,xmin,alpha)
#myp.makeplotres("datafit",res=200,disp=False)
