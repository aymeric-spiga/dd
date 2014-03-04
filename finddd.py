#! /usr/bin/env python
from ppclass import pp,ncattr
from ppplot import plot2d,save,writeascii
import numpy as np
from scipy.ndimage.measurements import minimum_position
import matplotlib.pyplot as mpl
import ppcompute

## method 3
#from skimage import filter,transform,feature
#import matplotlib.patches as mpatches 
###########

###############################################################################
# FIND DUST DEVILS
# filefile --> file
# dx --> horizontal resolution
# dt_out --> frequency of outputs (s) 
# lt_start --> starting local time
# halolim --> limit for halos
# method --> method used to detect vortices
# plotplot --> plot or not
# save --> save results or not
###############################################################################
def finddd(filefile,\
           timelist=None,\
           dt_out=50.,\
           lt_start=8.,\
           halolim=None,\
           method=1,\
           plotplot=False,\
           save=True):

    ###############################################################################
    ########################## FOR METHOD 1 FOR METHOD 1 ##########################
    ## FACLIST : how many sigmas below mean we start to consider this could be a vortex
    ## ... < 3.0 dubious low-intensity minima are caught
    ## ... 3 ~ 3.1 is a bit too low but helps? because size of pressure drop could be an underestimate of actual dust devil
    ## ... 3.2 ~ 3.3 is probably right, this is the one we choose for exploration [3.25]
    ## ... NB: values 3.2 to 3.6 yields the same number of devils but measured sizes could vary a bit
    ## ... > 3.6 makes strong minima disappear... especially >3.8
    ## ... EVENTUALLY 3.5-4 is better for problematic cases 
    ###############################################################################
    faclist = [3.75]

    ## (see below) NEIGHBOR_FAC is the multiple of std used to evaluate size
    ## --> 1: limit not discriminative enough. plus does not separate neighbouring vortices.
    ##        ... but interesting: gives an exponential law (because vortices are artificially merged?)
    ## --> 2.7: very good for method 1. corresponds usually to ~0.3
    ## --> 2: so-so. do not know what to think.
    neighbor_fac = 2.7

    ###############################################################################
    if save:
        myfile1 = open(filefile+'m'+str(method)+'_'+'1.txt', 'w')
        myfile2 = open(filefile+'m'+str(method)+'_'+'2.txt', 'w')
    ###############################################################################

    ## get the resolution within the file
    dx = ncattr(filefile,'DX') ; print "resolution in meters is: ",dx
    ## if no halolim is given, guess it from resolution
    if halolim is None:
        extentlim = 2000. # the putative maximum extent of a vortex in m
        halolim = extentlim / dx
        print "maximum halo size is: ",halolim

    ## mean and std calculations
    ## -- std is used in both methods for limits
    ## -- mean is only used in method 1
    print "calculate mean and std, please wait."
    ## -- get time series of 2D surface pressure
    psfc = pp(file=filefile,var="PSFC").getf()
    ## -- calculate mean and standard deviation
    ## -- ... calculating std at all time is not right!
    ## -- ... for mean value though, similar results with both methods
    mean = np.mean(psfc,dtype=np.float64)
    std = np.std(psfc,dtype=np.float64)
    damax = np.max(psfc)
    damin = np.min(psfc)
    ## some information about inferred limits
    print "**************************************************************"
    print "MEAN",mean
    print "STD",std
    print "LIMIT FOR PRESSURE MINIMUM:",-np.array(faclist)*std
    print "LIMIT FOR EVALUATING SIZE OF A GIVEN LOCAL MINIMUM",-neighbor_fac*std
    print "**************************************************************"
    
    # if no timelist is given, take them all
    if timelist is None:
        sizet = psfc.shape[0]
        print "treat all time values: ",sizet
        timelist = range(0,sizet-1,1)

    ## LOOP ON TIME
    for time in timelist:
    
     ## get 2D surface pressure at a given time
     ## (this is actually so quick we don't use psfc above)
     psfc2d = pp(file=filefile,var="PSFC",t=time).getf()
     #ustm = pp(file=filefile,var="USTM",t=time).getf()

     ## MAIN ANALYSIS. LOOP ON FAC. OR METHOD.
     for fac in faclist:
     #fac = 3.75
     #for method in [2,1]:
  
      ## initialize arrays
      tabij = [] ; tabsize = [] ; tabdrop = []
      tabijcenter = [] ; tabijvortex = [] ; tabijnotconv = [] ; tabdim = []

      ################ FIND RELEVANT POINTS TO BE ANALYZED
      ## lab is 1 for points to be treated by minimum_position routine
      ## we set elements at 1 where pressure is under mean-fac*std
      ## otherwise we set to 0 because this means we are close enough to mean pressure (background)
      lab = np.zeros(psfc2d.shape)
      if method == 1:
          # method 1: standard deviation
          lab[np.where(psfc2d < mean-fac*std)] = 1
      elif method == 2:
          # method 2: polynomial fit
          # ... tried smooth but too difficult, not accurate and too expensive
          deg = 5 #plutot bien (deg 10 ~pareil) mais loupe les gros (~conv cell)
          #deg = 2 #pas mal mais false positive (pareil 3-4 mm si un peu mieux)
          #        #(OK now with fixing the false positive bug)
          nx = psfc2d.shape[1] ; ny = psfc2d.shape[0]
          xxx = np.array(range(nx)) ; yyy = np.array(range(ny))
          anopsfc2d = psfc2d*0. ; polypsfc2d = psfc2d*0.
          for iii in range(0,nx,1):
             poly = np.poly1d(np.polyfit(yyy,psfc2d[iii,:],deg))
             polypsfc2d[iii,:] = poly(yyy)
          for jjj in range(0,ny,1):
             poly = np.poly1d(np.polyfit(xxx,psfc2d[:,jjj],deg))
             polypsfc2d[:,jjj] = 0.5*polypsfc2d[:,jjj] + 0.5*poly(xxx)
          ## smooth a little to avoid 'crosses' (plus, this removes PBC problems)
          polypsfc2d = ppcompute.smooth2diter(polypsfc2d,n=deg)
          # compute anomaly and find points to be explored
          anopsfc2d = psfc2d - polypsfc2d
          limlim = fac*std ## same as method 1
          lab[np.where(anopsfc2d < -limlim)] = 1
      elif method == 3:
          # method 3 : find centers of circle features using a Hough transform
          # initialize the array containing point to be further analyzed
          lab = np.zeros(psfc2d.shape)
          ## prepare the field to be analyzed by the Hough transform
          ## normalize it in an interval [-1,1]
          field = psfc2d
          # ... test 1. local max / min used for normalization.
          mmax = np.max(field) ; mmin = np.min(field) ; dasigma = 2.5
          ## ... test 2. global max / min used for normalization.
          #mmax = damax ; mmin = damin ; dasigma = 1.0 #1.5 trop restrictif
          spec = 2.*((field-mmin)/(mmax-mmin) - 0.5)
          # perform an edge detection on the field
          # ... returns an array with True on edges and False outside
          # http://sciunto.wordpress.com/2013/03/01/detection-de-cercles-par-une-transformation-de-hough-dans-scikit-image/    
          edges = filter.canny(filter.sobel(spec),sigma=dasigma)
          # initialize plot for checks
          if plotplot:
            fig, ax = mpl.subplots(ncols=1, nrows=1, figsize=(10,10))
            ax.imshow(field, cmap=mpl.cm.gray)
          # detect circle with radius 3dx. works well.
          # use an Hough circle transform
          radii = np.array([3]) 
          hough_res = transform.hough_circle(edges, radii)
          # analyze results of the Hough transform
          nnn = 0 
          sigselec = neighbor_fac
          #sigselec = 3.
          for radius, h in zip(radii, hough_res):
            # number of circle features to keep
            # ... quite large. but we want to be sure not to miss anything.
            nup = 30 
            maxima = feature.peak_local_max(h, num_peaks=nup)
            # loop on detected circle features
            for maximum in maxima:
              center_x, center_y = maximum - radii.max()
              # nup is quite high so there are false positives.
              # ... but those are easy to detect
              # ... if pressure drop is unclear (or inexistent)
              # ... we do not take the point into account for further analysis
              # ... NB: for inspection give red vs. green color to displayed circles
              diag = field[center_x,center_y] - (mean-sigselec*std)
              if diag < 0:  
                  col = 'green'
                  nnn = nnn + 1
                  lab[center_x,center_y] = 1
              else:
                  col = 'red'
              # draw circles
              if plotplot:
                circ = mpatches.Circle((center_y, center_x), radius,fill=False, edgecolor=col, linewidth=2)
                ax.add_patch(circ)
          if plotplot:
            mpl.title(str(nnn)+" vortices") ; mpl.show()
            mpl.close()

      ## while there are still points to be analyzed...
      while 1 in lab:
        ## ... get the point with the minimum field values
        if method == 1 or method == 3:
            ij = minimum_position(psfc2d,labels=lab)
        elif method == 2:
            ij = minimum_position(anopsfc2d,labels=lab)
        ## ... store the indexes of the point in tabij
        tabij.append(ij)
        ## ... remove the point from labels to be further explored by minimum_position
        lab[ij] = 0
    
      ################ GET SIZES BASED ON THOSE FOUND POINTS
      ## reslab is the same as lab
      ## except for scanning purpose we keep the information
      ## about how went the detection
      ## --> and we set a lower fac
      ## --> above a high fac is good to catch only strong vortices
      ## --> but here a casual fac=3 is better to get accurate sizes
      ## --> or even lower as shown by plotting reslab 
      reslab = np.zeros(psfc2d.shape)
      if method == 1 or method == 3:
          reslab[np.where(psfc2d < mean-neighbor_fac*std)] = 1
      elif method == 2:
          reslab[np.where(anopsfc2d < -neighbor_fac*std)] = 1
     
      ## initialize halomax and while loop
      ## HALOMAX : maximum halo defined around a minima to evaluate the size
      ## ... halomax must be large enough to encompass a vortex
      ## ... but not too large otherwise neighboring vortex are caught
      halomax = 3 ; notconv = 9999 ; yorgl = 9999
    
      ## WHILE LOOP on HALOMAX exploration
      while ( notconv > 0 and halomax < halolim and yorgl != 0 ):
       # now browse through all points caught in previous loop
       for ij in tabij:
        ## ... OK. take each indexes found before with minimum_position
        i,j = ij[0],ij[1]
        ## EITHER
        ## ... if reslab is already 0, we do not have to do anything
        ## ... because this means point is already part of detected vortex
        ## OR
        ## ... if the ij couple is already in a previously detected vortex
        ## ... we don't actually need to consider it again
        ## ... this is necessary otherwise (sometimes a lot) of false positives
        if reslab[i,j] <= 0 or ij in tabijvortex:
          pass
        else:
          ## ... then define a growing halo around this point
          ## ... we make the halo grow until convergence of 
          ## ... the number of under-limit points (i.e. with lab=1) within halo
          ## ... which means the whole vortex is encompassed
          ## ... we start with a halo of 1 point around minimum point
          ## ... and we end when we reached halomax which likely means dubious case
          halo = 1 ; nmesh = -9999. ; prevmesh = 9999. ; notconverged = False
          while nmesh != prevmesh and halo <= halomax:
              ## ... save the number of vortex points calculated at previous iteration
              prevmesh = nmesh
              ## ... define a halo around the minimum point
              minx,maxx,miny,maxy = i-halo,i+halo+1,j-halo,j+halo+1
              ## ... treat the boundary case (TBD: periodic boundary conditions)
              if minx < 0: minx = 0
              if miny < 0: miny = 0
              if maxx > psfc2d.shape[0]: maxx = psfc2d.shape[0]
              if maxy > psfc2d.shape[1]: maxy = psfc2d.shape[1]
              ## ... define the patch, made of the halo of points
              patch = reslab[minx:maxx,miny:maxy]       
              ## ... count how many 1 are inside this patch
              ## ... these are points for which value is >0
              ## ... because close to mean is 0 and already caught is <0
              ## ... and not converged are >1 so still to be scanned
              nmesh = len(np.where(patch >= 1)[0])
              ## ... in case widening halo adds 1-2 points only
              ## ... we consider vortex is encompassed
              ## ... this helps to separate neighboring vortices
              ## ... and only yields marginal inaccuracies
              if halo > 3 and abs(nmesh-prevmesh)<=2: prevmesh = nmesh
              ## ... if with halo=1 we caught a 1-point vortex, end here because spurious
              if halo == 1 and nmesh == 1: prevmesh = nmesh
              ## ... with the last halo=halomax test, we can detect not-converged cases
              if halo == halomax and nmesh != prevmesh: notconverged = True
              ## ... increment halo size before the next iteration
              halo = halo + 1
          ## ... now we got the grid points encompassed by the vortex in nmesh
          ## ... then to continue scanning we store results in reslab
          if notconverged:
            ## multiply reslab by 2. 
            ## --> if it is >1, point could be part of another vortex found later
            ## --> if it is <0, point is already within a caught vortex, so reslab should remain <0
            ## --> if it is =0, this should remain 0 because close to average
            reslab[i,j] = reslab[i,j]*2
            tabijnotconv.append(ij)
          else:
            ## OK. this is most likely an actual vortex. we get the drop.
            ## we multiply by mesh area, then square to get approx. size of vortex
            if method == 1 or method ==3:
                drop = -psfc2d[i,j]+mean
            else:
                drop = -anopsfc2d[i,j]
            size = int(np.sqrt(nmesh*dx*dx))
            ## if vortex is too small (i.e. too close to mesh grid resolution) we don't store information
            ## ... we just remove the patch from labels to be further explored
            facdx = 2.
            if size < facdx*dx:
                reslab[minx:maxx,miny:maxy] = 0
            ## otherwise it is a VORTEX! we store info in arrays
            else:
                ## we put the vortex points into tabijvortex so that 
                ## a (i,j) couple into a vortex is not considered in further explorations
                ## -- also we evaluate the x-size and y-size of the vortex (max distance to center)
                ## -- this can be useful to detect not-so-round fake vortices (convective gusts?)
                maxw = -9999 ; maxh = -9999 ; maxu = -9999
                for iii in range(minx,maxx):
                 for jjj in range(miny,maxy):
                  if reslab[iii,jjj] != 0:
                     tabijvortex.append((iii,jjj))
                     width = np.abs(iii-ij[0])
                     if width > maxw: maxw = width
                     height = np.abs(jjj-ij[1])
                     if height > maxh: maxh = height
                     #ustar = ustm[iii,jjj]
                     #if ustar > maxu: maxu = ustar
                #print size,drop,ustar
                ## store info in dedicated arrays
                tabdim.append((maxw*dx,maxh*dx))
                tabsize.append(size)
                tabdrop.append(drop)
                tabijcenter.append(ij)
                #print "... VORTEX!!!! size %.0f drop %.1f coord %.0f %.0f" % (size,drop,i,j)
                ## we remove the patch from labels to be further explored
                ## we multiply reslab by -1 to plot detection maps
                reslab[minx:maxx,miny:maxy] = patch*-1
     
       ## count how many points are not converged and left to be analyzed
       notconv = len(np.where(reslab > 1)[0])
       yorgl = len(np.where(reslab == 1)[0])
    
       ## increment halomax
       ## to speed-up the increment is slightly increasing with considered halomax
       halomax = halomax + halomax / 2
    
      ## just for simpler plots.
      reslab[reslab > 2] = 4
      reslab[reslab < -2] = -4

      ## give some info to the user
      if len(tabsize) > 0:
        nvortex = len(tabsize)
        maxsize = np.max(tabsize)
        maxdrop = np.max(tabdrop)
      else:
        nvortex = 0
        maxsize = 0
        maxdrop = 0.
      notconv = len(np.where(reslab > 1)[0])
      print "t=%3.0f / n=%2.0f / s_max=%4.0f / d_max=%4.1f / halo_out=%3.0f / notconvp=%3.1f" \
            % (time,nvortex,maxsize,maxdrop,halomax,100.*notconv/float(reslab.size))        
    
      ## save results in a text file
      if save:
          # convert t in local time
          ttt = lt_start + time*dt_out/3700.      
          # write files
          myfile2.write( "%5.2f ; %5.0f ; %5.0f ; %7.2f\n" % (ttt,nvortex,maxsize,maxdrop) )
          for iii in range(len(tabsize)):
              myfile1.write( "%5.2f ; %5.0f ; %7.2f ; %5.0f ; %5.0f\n" \
              % (ttt,tabsize[iii],tabdrop[iii],tabdim[iii][0],tabdim[iii][1]) )
    
      #### PLOT PLOT PLOT PLOT
      if (nvortex>0 and plotplot) or (nvortex>0 and maxsize > 800.):
       mpl.figure(figsize=(12,8))
       myplot = plot2d()
       myplot.x = np.array(range(psfc2d.shape[1]))*dx/1000.
       myplot.y = np.array(range(psfc2d.shape[0]))*dx/1000.
       myplot.title = str(nvortex)+" vortices found (indicated diameter / pressure drop)"
       myplot.xlabel = "x distance (km)"
       myplot.ylabel = "y distance (km)"
       if method > 0:
       #if method == 1:
           #myplot.field = ustm 
           myplot.f = psfc2d
           #myplot.field = polypsfc2d
           myplot.vmin = -2.*std 
           myplot.vmax = +2.*std
           myplot.vmin = mean - 6.*std
           myplot.vmax = mean + 6.*std
       else:
           myplot.field = anopsfc2d
           myplot.vmin = -1.5
           myplot.vmax = 0.5
       myplot.fmt = "%.1f"
       myplot.div = 20
       myplot.colorbar = "spectral"
       myplot.make()
      
       ### ANNOTATIONS
       for iii in range(len(tabsize)):
        ij = tabijcenter[iii]
        coord1 = ij[1]*dx/1000.
        coord2 = ij[0]*dx/1000.
        txt = "%.0f/%.1f" % (tabsize[iii],tabdrop[iii])
        mpl.annotate(txt,xy=(coord1,coord2),
             xytext=(-10,-30),textcoords='offset points',ha='center',va='bottom',\
             bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),\
             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3',color='red'),\
             size='small')
    
       ###show detection
       #lev = [-4,-2,-1,0,1,2,4]
       #lev = [-1,0]
       #mpl.contourf(myplot.absc,myplot.ordi,reslab,alpha=0.9,cmap=mpl.cm.get_cmap("binary_r"),levels=lev)
    
       ### SHOW OR SAVE IN FILE
       mpl.show()
       #save(mode="png",filename="detectm"+"_"+str(time)+"_"+str(method),folder="detect/",includedate=False)
    
    ## close data files
    if save:
        myfile1.close()
        myfile2.close()
