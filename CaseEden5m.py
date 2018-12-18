#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 00:09:11 2018

@author: Xiaodong Ming

"""
import os
import numpy as np
import time
from ArcGridDataProcessing import arcgridread,MaskInterp
import HiPIMS_IO as HPIO
rootPath = '/home/b1055010/GeoClasses/release/bin/Eden5mMGPUs/'
os.chdir(rootPath)

start = time.clock()
rainMaskName = 'Eden_rainfall_mask_50m.asc';
maskMat,maskHead,maskExtent = arcgridread(rainMaskName)
demName = 'eden_5m_withriver.asc'
demMat,demHead,demExtent = arcgridread(demName)
rainMask = MaskInterp(maskMat,maskHead,demMat,demHead)
rainMask[np.isnan(rainMask)]=0
# read gauge data
gauges_pos  = np.loadtxt('gauges_pos_20m.dat')
rain_source = np.loadtxt('rainRadar2015_120300.txt')
import matplotlib.pyplot as mplPy
#create original rain mask map
fig, ax = mplPy.subplots(1, figsize=(6, 10))
img = mplPy.imshow(maskMat,extent=demExtent)
mplPy.colorbar(img,fraction=0.05, pad=0.01,shrink=0.92,ax=ax)
ax.axes.grid(linestyle='--',linewidth=1)
fig.savefig('EdenMaskoriginal.png', dpi=300)

#create new rain mask map
fig, ax = mplPy.subplots(1, figsize=(6, 10))
img = mplPy.imshow(rainMask,extent=demExtent)
mplPy.colorbar(img,fraction=0.05, pad=0.01,shrink=0.92,ax=ax)
ax.axes.grid(linestyle='--',linewidth=1)
fig.savefig('EdenMask.png', dpi=300)

end = time.clock()
print('total time elapse: '+str(end-start))
#%%
fileToCreate = ['precipitation_mask','precipitation_source','boundary','gauges_pos']
start = time.clock()
HPIO.InputSetup_MG(rootPath,demMat,demHead,numSection=8,h0=0,
                        fileToCreate=fileToCreate,
                        rain_mask=rainMask,
                        rain_source = rain_source,
                        gauges_pos=gauges_pos)
end = time.clock()
print('total time elapse: '+str(end-start))