import numpy as np

## nmesh,maxw,maxh,reslab,tabijvortex=gethalo(ij,reslab,halomax,tabijvortex)
def gethalo(ij,reslab,halomax,tabijvortex):
          i,j = ij[0],ij[1]
          ## ... then define a growing halo around this point
          ## ... we make the halo grow until convergence of 
          ## ... the number of under-limit points (i.e. with lab=1) within halo
          ## ... which means the whole vortex is encompassed
          ## ... we start with a halo of 1 point around minimum point
          ## ... and we end when we reached halomax which likely means dubious case
          halo = 1 ; nmesh = -9999. ; prevmesh = 9999. ; notconverged = False
          nx = reslab.shape[1] ; ny = reslab.shape[0]
          while nmesh != prevmesh and halo <= halomax:
              ## ... save the number of vortex points calculated at previous iteration
              prevmesh = nmesh
              ## ... define a halo around the minimum point 
              ## ... (the +1 in upper bounds is necessary for correct patch)
              minx,maxx,miny,maxy = i-halo,i+halo+1,j-halo,j+halo+1
              ## ... treat the boundary case (TBD: periodic boundary conditions)
              if minx < 0: minx = 0
              if miny < 0: miny = 0
              if maxx > nx: maxx = nx
              if maxy > ny: maxy = ny
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
            nmesh,maxw,maxh=None,None,None
          else:
            ## if vortex is too small (i.e. too close to mesh grid resolution) we don't store information
            ## ... we just remove the patch from labels to be further explored
            ## ... not exactly size < facdx*dx given that size is a circle, but should be enough
            if nmesh <= 2:
                reslab[minx:maxx,miny:maxy] = 0
                nmesh,maxw,maxh=None,None,None
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
                ## we remove the patch from labels to be further explored
                ## we multiply reslab by -1 to plot detection maps
                reslab[minx:maxx,miny:maxy] = patch*-1
          return nmesh,maxw,maxh,reslab,tabijvortex

