#! /usr/bin/env python
import ppplot
import numpy as np

### STATDD
def statdd(namefile,limtime=None):

    # load data
    data = np.loadtxt(namefile+".txt",delimiter=";")
    t = data[:,0] ; n = data[:,1] ; s = data[:,2] ; d = data[:,3] 
    
    # remove size and drop point when no vortex detected
    d[np.where(n==0)] = np.nan ; s[np.where(n==0)] = np.nan
    
    ## PLOTS
    number = ppplot.plot1d()
    number.f = n
    number.x = t
    number.linestyle = ''
    number.marker = '.'
    number.color = 'b'
    number.xlabel = "Local time (hour)"
    number.ylabel = "Detected vortices"
    number.xmax = limtime
    #number.makeshow()
    number.make()
    ppplot.save(mode="pdf",filename=namefile+"_detect")
    
    drop = ppplot.plot1d()
    drop.f = d
    drop.x = t
    drop.linestyle = ''
    drop.marker = '.'
    drop.color = 'r'
    drop.fmt = "%.1f"
    drop.xlabel = "Local time (hour)"
    drop.ylabel = "Maximum drop of detected vortices (Pa)"
    drop.xmax = limtime
    #drop.makeshow()
    drop.make()
    ppplot.save(mode="pdf",filename=namefile+"_maxdrop")
    
    size = ppplot.plot1d()
    size.f = s
    size.x = t
    size.linestyle = ''
    size.marker = '.'
    size.color = 'g'
    size.xlabel = "Local time (hour)"
    size.ylabel = "Maximum size of detected vortices (m)"
    size.xmax = limtime
    #size.makeshow()
    ppplot.save(mode="pdf",filename=namefile+"_maxsize")

#    try:
#      w = data[:,4]
#      w[np.where(n==0)] = np.nan
#      fric = ppplot.plot1d()
#      fric.f = w
#      fric.x = t
#      fric.fmt = "%.2f"
#      fric.linestyle = ''
#      fric.marker = '.'
#      fric.color = 'm'
#      fric.xlabel = "Local time (hour)"
#      fric.ylabel = "Maximum friction velocity within a vortex (m/s)"
#      fric.xmax = limtime
#      fric.makeshow()
#    except:
#      pass
