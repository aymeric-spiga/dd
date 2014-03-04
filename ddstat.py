#! /usr/bin/env python
import ppplot
import numpy as np

### STATDD
def statdd(namefile):

    # load data
    data = np.loadtxt(namefile,delimiter=";")
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
    number.makeshow()
    
    drop = ppplot.plot1d()
    drop.f = d
    drop.x = t
    drop.linestyle = ''
    drop.marker = '.'
    drop.color = 'r'
    drop.fmt = "%.1f"
    drop.xlabel = "Local time (hour)"
    drop.ylabel = "Maximum drop of detected vortices (Pa)"
    drop.makeshow()
    
    size = ppplot.plot1d()
    size.f = s
    size.x = t
    size.linestyle = ''
    size.marker = '.'
    size.color = 'g'
    size.xlabel = "Local time (hour)"
    size.ylabel = "Maximum size of detected vortices (m)"
    size.makeshow()
