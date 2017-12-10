#! /usr/bin/env python
from ddhist import histodd
from ddstat import statdd

vlimwind = None

################################################################################
### HISTODD
### namefile: text file to be analyzed
### drop: look at size (False) or pressure drop (True)
### typefit: use fitfunc above 1, 2, or 3
### nbins: bins (we expect fit do not depend too much on this -- to be checked)
### limrest: restrict limit (multiple of dx) -- greater equal 
### limtime: use only data earlier than this local time 
### limdrop: use only data with deeper drop than value
### addtitle: add a name for the examined case
################################################################################
#def histodd(namefile,drop=False,typefit=1,nbins=12,limrest=4,limtime=None,limdrop=0.3,addtitle=""):
### bins = {7,10,12,15,20,30,50}
################################################################################

###------------------------------------------------------------------------------------------------
### number -- comment   -- out -- ztop -- dx -- nx  -- nz  -- wind -- nphys -- lt
### 160564 -- phoenix   -- 50s -- 6km  -- 10 -- 369 -- 101 -- 0    -- n     -- 9-14 (16 for r)
### 188324 -- phoenix   -- 50s -- 6km  -- 10 -- 329 -- 201 -- 0    -- n     -- 9-13
### 172097 -- exm_shear -- 50s -- 5km  -- 12 -- 249 -- 401 -- 5    -- y     -- 13h~13h45
###------------------------------------------------------------------------------------------------

## power law: sizes OK, drops not OK (exponential-ish)
vlimtime = 12 #None egalement OK en general
vlimdrop = 0.3 #0.0 moins bon pressure drop
vlimrest = 5 #3 moins bon pressure drop
vbins = 10 #15-20 relativement idem. 30-50 pour verifier.
mmm = "m1" #m2 not enough. anyway m2 is weird (see finddd)

#fff = "188324p.nc" #10m cool. PL beautiful for sizes. EXP bof for drops.
#fff = "r160564p.nc" # 10m cool. may be better with more vortices.
#mmm = "m3"
   ## m3 includes a bit more larger vortices than m1;
   ## -- and seems closer to exponential to m3 for drops
#fff = "r160564p.ncsigma15" # peu de differences en m1 et m3 avec r160564p
#fff = "160564p.nc" # 10m bof. not enough vortices.
#fff = "13526p.nc" # 25m cool. PL for sizes AND drops. PL un peu limite sur les large devils? because merging?
#fff = "156487.nc" # 25m less cool
#fff = "172097.nc" # does not work?
#fff = "191798.nc" # bad
#fff = "2007p.nc" # 100m bad.

mmm = "m3"
fff = "test160564"
vbins = 15
vlimdrop = 0.25 ##0.21 limit for a pressure minimum
vlimrest = 3
vlimtime = 13

fff = "188324p.nc"
vlimtime = 11
#vlimtime = 12
##vlimdrop = 0.3 #interessant aussi avec 0.25 PL 2.4 au lieu de 2.8 bonne idee du range
##vlimwind = 0.3 #ne change pas grand chose
mmm = "m1"
mmm = "m3"

#fff = "test160564"
#vlimtime = 12

#from ddhist import fdropsize,fdropwind,fsizewind
#fdropsize("BIGLES10m_wind5_PSFC_9-11_stride6.ncm3")
#fdropwind("BIGLES10m_wind5_PSFC_9-11_stride6.ncm3")
#fsizewind("BIGLES10m_wind5_PSFC_9-11_stride6.ncm3")
#exit()
#
#vlimwind = 0.5 ; vlimtime = None
#vlimdrop = 0.20 ; vlimrest = 4 
#histodd("BIGLES10m_wind5_PSFC_9-11_stride6.ncm3",nbins=vbins,typefit=1,limrest=vlimrest,limtime=vlimtime,limdrop=vlimdrop,limwind=vlimwind)
#histodd("BIGLES10m_wind5_PSFC_9-11_stride6.ncm3",drop=True,nbins=vbins,typefit=2,limrest=vlimrest,limdrop=vlimdrop,limtime=vlimtime,limwind=vlimwind)
#exit()



##################################################
##################################################
##################################################
ffftab = []
#ffftab.append("BIGLES10m_wind10_PSFC_9-11_stride20.nc")
#ffftab.append("BIGLES10m_wind5_PSFC_9-11_stride20.nc")
#ffftab.append("BIGLES10m_wind10_PSFC_9-11.nc")
#ffftab.append("BIGLES10m_wind5_PSFC_9-11.nc")
###
bintab = []
bintab.append(12)
bintab.append(30)
bintab.append(100)
###
mmm = "m3"
###
vlimdrop = None ; vlimtime = None
###vlimdrop = 0.20 ; vlimrest = 3 # 3 points de grille trop faible. population 30>40m systematiquement dessous
vlimdrop = 0.20 ; vlimrest = 4 
###vlimdrop = 0.30 ; vlimrest = 4 
###vlimdrop = 0.35 ; vlimrest = 5
###vlimdrop = 0.25 ; vlimrest = 4
###vlimdrop = 0.50 ; vlimrest = 8
###


### test avec limitation aux wind les plus forts
ffftab = ["BIGLES10m_wind5_PSFC_9-11_stride6.nc"] ; vlimwind = 0.8 ; vlimdrop = None ; vlimtime = None
vlimwind = 0.5 # 0.5 équivalent // 0.8 n'importe quoi

for fff in ffftab:
 for vbins in bintab:
  ##################################################
  typefit_size = 1
  typefit_drop = 2
  ##################################################
  histodd(fff+mmm,          nbins=vbins,typefit=typefit_size,limrest=vlimrest,limtime=vlimtime,limdrop=vlimdrop,limwind=vlimwind)
  histodd(fff+mmm,drop=True,nbins=vbins,typefit=typefit_drop,limrest=vlimrest,limdrop=vlimdrop,limtime=vlimtime,limwind=vlimwind)
  #fdropsize(fff+mmm)
  #statdd(fff+mmm+"_2.txt",limtime=vlimtime)
  ##################################################
  ##################################################
