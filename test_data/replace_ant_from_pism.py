#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

ice5g_file = '/home/albrecht/vilma/data/ice5g.nc'
pism2vilma_file = '/home/albrecht/vilma/data/pism2vilma.nc' 

secperyear=3600.0*24.0*365.0

ncf = nc.Dataset(ice5g_file,'r')
ice            = ncf.variables['Ice'][:]
epoch          = ncf.variables['epoch'][:]
ncf.close()
ntemp,nlat,nlon=np.shape(ice)
print(np.shape(ice))

ncf2 = nc.Dataset(pism2vilma_file,'a')
thk            = ncf2.variables['thk'][:]
time           = ncf2.variables['time'][:]/secperyear/1000.0
#ncf2.close()

#FIXME: Actually PISM extra time is defined at 124.5 ... 0.5 kyr BP


mintj=np.min(epoch)
print(mintj)
print(time)

notinice5g=[]
for i,ti in enumerate(time):
    if ti not in epoch:
      notinice5g.append(ti) 
#print notinice5g

for i,ti in enumerate(time):
  for j,tj in enumerate(epoch):
    if tj>=mintj: # and tj<-110.0:
      if ti==tj:
        print(ti,'match')
        for ilat in range(nlat):
          for ilon in range(nlon):
            if thk[i,ilat,ilon]<0.0: #-9.e+33:
              thk[i,ilat,ilon]=ice[j,ilat,ilon]
      elif ti in notinice5g:
        if (ti-tj)==1.0:
          print(ti,'interpolate',tj,epoch[j+1])
          for ilat in range(nlat):
            for ilon in range(nlon):
              if thk[i,ilat,ilon]<0.0: #-9.e+33:
                thk[i,ilat,ilon]=(ice[j,ilat,ilon]+ice[j+1,ilat,ilon])*0.5

#ncf.variables['Ice'][:]=ice[:]
ncf2.variables['thk'][:,:,:]=thk[:,:,:]
ncf2.variables['time'][:]=time[:]
ncf2.close()

