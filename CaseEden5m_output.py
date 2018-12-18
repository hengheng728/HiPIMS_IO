#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:36:25 2018

@author: Xiaodong Ming
"""
import os
import numpy as np
import OutputManage as OM
import matplotlib.pyplot as mplPy
rootPath = '/home/b1055010/GeoClasses/release/bin/Eden6_5mMGPUs/'
os.chdir(rootPath)
numSection=6
fileTag = 'h_21600'
grid,head,extent=OM.CombineGridFile(rootPath,numSection,fileTag)
fileName = fileTag+'.asc'
print(fileName+' is combined')
OM.ArcgridwriteGZip(fileName,grid,head)

np.warnings.filterwarnings('ignore')
grid[grid<0.1]=np.nan
fig, ax = mplPy.subplots(1, figsize=(6, 10))
img = mplPy.imshow(grid,extent=extent)
mplPy.colorbar(img,fraction=0.05, pad=0.01,shrink=0.92,ax=ax)
ax.axes.grid(linestyle='--',linewidth=0.5)
fig.savefig(fileTag+'.png', dpi=800)
print('figure '+fileTag+'.png is drawn')

Time=21600
OM.Backup2Initial(rootPath,numSection,Time,remove=True)
