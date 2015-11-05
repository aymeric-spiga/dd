#! /usr/bin/env python
from ppclass import pp,ncattr
from ppplot import plot2d,save,writeascii
import numpy as np
from scipy.ndimage.measurements import minimum_position
import matplotlib.pyplot as mpl
import ppcompute

from gethalo import gethalo


###############################################################################
# FIND DUST DEVILS
# filefile --> file
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
           filewind=None,\
           save=True):

    if method == 3:
       print "importing additional scikit-image packages"
       from skimage import filter,transform,feature
       import matplotlib.patches as mpatches
       from scipy import ndimage


    ###############################################################################
    ########################## FOR METHOD 1 FOR METHOD 2 ##########################
    ## FACLIST : how many sigmas below mean we start to consider this could be a vortex
    ## ... < 3.0 dubious low-intensity minima are caught
    ## ... 3 ~ 3.1 is a bit too low but helps? because size of pressure drop could be an underestimate of actual dust devil
    ## ... 3.2 ~ 3.3 is probably right, this is the one we choose for exploration [3.25]
    ## ... NB: values 3.2 to 3.6 yields the same number of devils but measured sizes could vary a bit
    ## ... > 3.6 makes strong minima disappear... especially >3.8
    ## ... EVENTUALLY 3.5-4 is better for problematic cases 
    ###############################################################################
    faclist = [3.75]

    ################################ FOR ALL METHODS
    ## (see below) NEIGHBOR_FAC is the multiple of std used to evaluate size
    ## --> 1: limit not discriminative enough. plus does not separate neighbouring vortices.
    ##        ... but interesting: gives an exponential law (because vortices are artificially merged?)
    ## --> 2.7: very good for method 1. corresponds usually to ~0.3
    ## --> 2: so-so. do not know what to think. but usually too low.
    neighbor_fac = 2.7
    #### METHOD 3 --> optimizing neighbor_fac with visual checks and superimposing wind friction (max must be at boundaries)
    ##neighbor_fac = 1.5 # too low --> vortices too large + false positives
    ##neighbor_fac = 2.7 # optimal --> good for separation, only a slight underestimation of size
    #neighbor_fac = 3.0 # too high --> excellent for separation, but size quite underestimated
    ###############################################################################
    ###############################################################################

    ###############################################################################
    ###############################################################################
    if save:
        myfile1 = open(filefile+'m'+str(method)+'_'+'1.txt', 'w')
        myfile2 = open(filefile+'m'+str(method)+'_'+'2.txt', 'w')
        if filewind is not None:
           myfile3 = open(filewind+'m'+str(method)+'_'+'1.txt', 'w')
    ###############################################################################
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
    psfc = pp(file=filefile,var="PSFC",verbose=True).getf()

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
     if filewind is not None:
       ustm = pp(file=filewind,var="USTM",t=time).getf()

     ## MAIN ANALYSIS. LOOP ON FAC. OR METHOD.
     for fac in faclist:
     #fac = 3.75
     #for method in [2,1]:
  
      ## initialize arrays
      tabij = [] ; tabsize = [] ; tabdrop = []
      tabijcenter = [] ; tabijvortex = [] ; tabdim = [] 
      tabwind = []

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
          # method 3 : find centers of circle features using image processing techniques

          # initialize the array containing point to be further analyzed
          lab = np.zeros(psfc2d.shape)
          ### field to analyze: pressure
          ### --- apply a Laplace transform to highlight drops
          field = ndimage.laplace(psfc2d)

          ### prepare the field to be analyzed 
          ### by the Hough transform or Blob detection
          ### --> normalize it in an interval [-1,1]
          ### --> NB: dasigma serves later for Hough transform
          ### --> NB: polynomial de-trending does not seem to help
          # ... test 1. local max / min used for normalization.
          mmax = np.max(field) ; mmin = np.min(field) ; dasigma = 2.5
          ## ... test 2. global max / min used for normalization. bof.
          #mmax = damax ; mmin = damin ; dasigma = 1.0 #1.5 trop restrictif
          spec = 2.*((field-mmin)/(mmax-mmin) - 0.5)

          #### **** BLOB DETECTION ****
          #### Better than Hough transform for multiple adjacent vortices
          #### log: best / dog or doh: miss small vortices, hence the majority
          #### SITE: http://scikit-image.org/docs/dev/auto_examples/plot_blob.html
          #### PUBLISHED: https://peerj.com/articles/453/
          ### --------------------------------------------        
          ### the parameters below are aimed for efficiency
          ### ... because anyway the actual size is not detected
          ### ... so setting max_sigma to a high value is not needed
          ### --------------------------------------------
          blobs = feature.blob_log(spec, max_sigma=3, num_sigma=3, threshold=0.05)
          ### a plot to check detection
          if plotplot:
            fig, ax = mpl.subplots(1, 1)
            what_I_plot = psfc2d #spec #field
            ax.imshow(what_I_plot, cmap=mpl.cm.gray)
          ### store the detected points in lab
          for blob in blobs:
            center_x, center_y, r = blob
            lab[center_x,center_y] = 1
            if plotplot:
              circ = mpatches.Circle((center_y, center_x), r*np.sqrt(2), fill=False, edgecolor='green', linewidth=2)
              ax.add_patch(circ)
          if plotplot: mpl.show()

#################################### BEGIN TEST HOUGH TRANSFORM
#          # perform an edge detection on the field
#          # ... returns an array with True on edges and False outside
#          # http://sciunto.wordpress.com/2013/03/01/detection-de-cercles-par-une-transformation-de-hough-dans-scikit-image/    
#          edges = filter.canny(filter.sobel(spec),sigma=dasigma)
#          # initialize plot for checks
#          if plotplot:
#            fig, ax = mpl.subplots(ncols=1, nrows=1, figsize=(10,8))
#            ax.imshow(field, cmap=mpl.cm.gray)
#          ## detect circle with radius 3dx. works well. 5dx detection pretty similar.
#          ## use an Hough circle transform
#          radii = np.array([2,3])
#          hough_res = transform.hough_circle(edges, radii)
#          # analyze results of the Hough transform
#          nnn = 0 
#          sigselec = neighbor_fac
#          #sigselec = 3.
#          for radius, h in zip(radii, hough_res):
#            # number of circle features to keep
#            # ... quite large. but we want to be sure not to miss anything.
#            nup = 30 
#            maxima = feature.peak_local_max(h, num_peaks=nup)
#            # loop on detected circle features
#            for maximum in maxima:
#              center_x, center_y = maximum #- radii.max()
#              # nup is quite high so there are false positives.
#              # ... but those are easy to detect
#              # ... if pressure drop is unclear (or inexistent)
#              # ... we do not take the point into account for further analysis
#              # ... NB: for inspection give red vs. green color to displayed circles
#              diag = field[center_x,center_y] - (mean-sigselec*std)
#              ## uncomment below to keep all detections
#              #diag = -1
#              if diag < 0:  
#                  col = 'green'
#                  nnn = nnn + 1
#                  lab[center_x,center_y] = 1
#              else:
#                  col = 'red'
#              # draw circles
#              if plotplot:
#                circ = mpatches.Circle((center_y, center_x), radius,fill=False, edgecolor=col, linewidth=2)
#                ax.add_patch(circ)
#          if plotplot:
#            mpl.title(str(nnn)+" vortices")
#            if nnn>0: mpl.show()
#            mpl.close()
#################################### END TEST HOUGH TRANSFORM

      ## while there are still points to be analyzed...
      while 1 in lab:
        ## ... get the point with the minimum field values
        if method == 1: # or method == 3:
            ij = minimum_position(psfc2d,labels=lab)
        elif method == 2:
            ij = minimum_position(anopsfc2d,labels=lab)
        elif method == 3:
            ij = minimum_position(1-lab)
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
      #reslabf = np.zeros(psfc2d.shape) # TESTS
      if method == 1 or method == 3:
          reslab[np.where(psfc2d < mean-neighbor_fac*std)] = 1
          #reslabf[np.where(psfc2d < mean-neighbor_fac_fine*std)] = 1 # TESTS
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
          ## GET HALOS. SEE FUNCTION ABOVE.
          nmesh,maxw,maxh,reslab,tabijvortex=gethalo(ij,reslab,halomax,tabijvortex)

          ## OK. check this is actually a vortex.
          ## get the size. get the drop.
          ## store results in file
          if nmesh is not None:

            ## calculate size
            ## we multiply by mesh area, then square to get approx. size of vortex
            ## [if one wants to obtain equivalent diameter, multiply by 2.*np.sqrt(np.pi)]
            size = np.sqrt(nmesh*dx*dx)

            ## check size. if not OK recompute halo with more stringent zone around pressure minimum.
            ## -- NB: reslab and tabijvortex do not need to be changed again, was done just before
            ## --     however, we could have been a little bit more subtle to disentangle twin vortices        
            # if (np.abs(maxw-maxh)*dx/size > 0.33):
            # #if (np.sqrt(maxw*maxh*dx*dx) > size):
            #    #print "asymmetry!",np.abs(maxw-maxh)*dx,size
            #    nmesh,maxw,maxh,dummy,dummy=gethalo(ij,reslabf,halomax,tabijvortex)
            #    if nmesh is not None: size = int(np.sqrt(nmesh*dx*dx))
            #    #print "new values",np.abs(maxw-maxh)*dx,size

            ## calculate drop.
            if method == 1 or method == 3: drop = -psfc2d[i,j]+mean
            else: drop = -anopsfc2d[i,j]

            #############################################################
            ##### Check this is the actual minimum (only tested so far with method=3)
            #if method == 1 or method ==3:
            #  ## ... define a halo around the minimum point
            #  ix,ax,iy,ay = i-maxw,i+maxw+1,j-maxh,j+maxh+1
            #  ## ... treat the boundary case (TBD: periodic boundary conditions)
            #  nx = reslab.shape[1] ; ny = reslab.shape[0]
            #  if ix < 0: ix = 0
            #  if iy < 0: iy = 0
            #  if ax > nx: ax = nx
            #  if ay > ny: ay = ny
            #  ## ... keep real minimal value
            #  ## DOMAINMIN --> does not change a lot results (not worth it)
            #  domainmin = np.max(-psfc2d[ix:ax,iy:ay])+mean
            #  if drop < domainmin:
            #     print "corrected drop",drop,domainmin
            #     drop = domainmin
            #  ### DOMAINDROP --> leads to underestimate drops in most cases
            #  #domaindrop = np.max(psfc2d[ix:ax,iy:ay])-np.min(psfc2d[ix:ax,iy:ay])
            #  #drop = domaindrop
            #############################################################

            ## if available get info on friction velocity
            if filewind is not None:
              ## ... define a halo around the minimum point
              ix,ax,iy,ay = i-maxw,i+maxw+1,j-maxh,j+maxh+1
              ## ... treat the boundary case (TBD: periodic boundary conditions)
              nx = reslab.shape[1] ; ny = reslab.shape[0]
              if ix < 0: ix = 0
              if iy < 0: iy = 0
              if ax > nx: ax = nx
              if ay > ny: ay = ny
              ## WINDMAX
              windmax = np.max(ustm[ix:ax,iy:ay])
              tabwind.append(windmax)
            else:
              tabwind.append(0.) 

            ## store info in dedicated arrays
            tabdim.append((maxw*dx,maxh*dx))
            tabsize.append(size)
            tabdrop.append(drop)
            tabijcenter.append(ij)
            #print "... VORTEX!!!! size %.0f drop %.1f coord %.0f %.0f" % (size,drop,i,j)

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
        maxwind = np.max(tabwind)
      else:
        nvortex = 0
        maxsize = 0
        maxdrop = 0.
        maxwind = 0.
      notconv = len(np.where(reslab > 1)[0])
      print "t=%3.0f / n=%2.0f / s_max=%4.0f / d_max=%4.1f / halo_out=%3.0f / notconvp=%3.1f / wind=%4.1f" \
            % (time,nvortex,maxsize,maxdrop,halomax,100.*notconv/float(reslab.size),maxwind)        
    
      ## save results in a text file
      if save:
          # convert t in local time
          ttt = lt_start + time*dt_out/3700.      
          # write files
          myfile2.write( "%5.2f ; %5.0f ; %6.1f ; %8.3f ; %8.3f\n" % (ttt,nvortex,maxsize,maxdrop,maxwind) )
          for iii in range(len(tabsize)):
              myfile1.write( "%5.2f ; %6.1f ; %8.3f ; %5.0f ; %5.0f ; %8.3f\n" \
              % (ttt,tabsize[iii],tabdrop[iii],tabdim[iii][0],tabdim[iii][1],tabwind[iii]) )

      #### PLOT PLOT PLOT PLOT
      damaxsize = 10000.
      #damaxsize = 400.
      if (nvortex>0 and plotplot) or (nvortex>0 and maxsize > damaxsize):
      #if nvortex > 200:
       mpl.figure(figsize=(12,8))
       myplot = plot2d()
       myplot.x = np.array(range(psfc2d.shape[1]))*dx/1000.
       myplot.y = np.array(range(psfc2d.shape[0]))*dx/1000.
       myplot.title = str(nvortex)+" vortices found (indicated diameter / pressure drop)"
       myplot.xlabel = "x distance (km)"
       myplot.ylabel = "y distance (km)"
       if method != 2:
           #myplot.f = ustm 
           myplot.f = psfc2d
           #myplot.vmin = damin
           #myplot.vmax = damax
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
        txt = "%.0f/%.2f/%.0f" % (tabsize[iii],tabdrop[iii],100*np.abs(tabdim[iii][0]-tabdim[iii][1])/tabsize[iii])
        txt = "%.0f m / %.2f Pa" % (tabsize[iii],tabdrop[iii])
        mpl.annotate(txt,xy=(coord1,coord2),
             xytext=(-10,-30),textcoords='offset points',ha='center',va='bottom',\
             bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),\
             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3',color='red'),\
             size='small')
    
       ###show detection
       lev = [-4,-2,-1,0,1,2,4] ## all colours for detection cases
       lev = [-1,0] # show dubious areas as detection areas
       lev = [-1,1] # show dubious areas as no-detection areas
       mpl.contourf(myplot.x,myplot.y,reslab,alpha=0.9,cmap=mpl.cm.get_cmap("binary_r"),levels=lev)
    
       ### SHOW OR SAVE IN FILE
       mpl.show()
       #save(mode="png",filename="detectm"+"_"+str(time)+"_"+str(method),folder="detect/",includedate=False)
    
    ## close data files
    if save:
        myfile1.close()
        myfile2.close()
