#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:08:32 2018

@author: Xiaodong Ming
"""
import numpy as np
import matplotlib.pyplot as plt
import InputSetupFuncs as IS
import ArcGridDataProcessing as AP
#%% draw dem map with boundary cell scatters
def DEMwithBoundary(demMat,demHead,
                    boundList=[],
                    figureName=[],
                    figsize=(10, 6),dpi=300):
    # fig,ax = DEMwithBoundary(demMat,demHead,boundList,figureName)
    # plot DEM map with user-defined boundaries
    # if figureName is empty, then the figure will not be saved
    _, bnMat_outline = IS.Get_ID_BNoutline_Mat(demMat)
    bnMat = IS.BoundaryClassify(bnMat_outline,demHead,boundList)
    demExtent = AP.demHead2Extent(demHead)   
    np.warnings.filterwarnings('ignore')    
    boundVec = bnMat[bnMat>0]
    boundValues = list(np.unique(boundVec)) # value for each boundary type 0,1,2,...
    boundSubs = np.where(bnMat>0)
    boundX,boundY = AP.Sub2Map(boundSubs[0],boundSubs[1],demHead)
#***************************draw map    
    fig, ax = plt.subplots(1, figsize=figsize)
    # draw DEM
    img=plt.imshow(demMat,extent=demExtent)
    # colorbar
    fig.colorbar(img,ax=ax,fraction=0.1,shrink=0.9)
    # draw boundary cells
    cStr='mbkgcmbkgc'
    n=0
    for value in boundValues:
        pointcolor = cStr[n]
        n=n+1
        ind = boundVec==value
        ax.scatter(boundX[ind],boundY[ind],s=0.5,facecolor=pointcolor)
    # save figure
    if len(figureName)>0:
        fig.savefig(figureName, dpi=dpi)
    
    return fig,ax
#%% draw inundation map with domain outline
def InundationMap(zMat,zHead,depth=0.2,figureName=[],
                    figsize=(10, 6),dpi=300):
    # fig,ax = InundationMap(demMat,demHead,boundList,figureName)
    # plot inundation map with domain outline
    # depth: the minimum water depth to show in the map
    # if figureName is empty, then the figure will not be saved
    # outline cells:0
    _, bnMat_outline = IS.Get_ID_BNoutline_Mat(zMat)# outline cells:0,non-outline cells:-2
    zExtent = AP.demHead2Extent(zHead)   
    np.warnings.filterwarnings('ignore')    
    boundSubs = np.where(bnMat_outline==0)
    boundX,boundY = AP.Sub2Map(boundSubs[0],boundSubs[1],zHead)
#***************************draw map    
    fig, ax = plt.subplots(1, figsize=figsize)
    # draw inundation
    zMat = zMat+0
    zMat[zMat<depth]=np.nan
    img=plt.imshow(zMat,extent=zExtent)
    # colorbar
    fig.colorbar(img,ax=ax,fraction=0.1,shrink=0.8)
    # draw domain outline cells
    ax.scatter(boundX,boundY,s=0.2,facecolor='r')
    ax.axes.grid(linestyle='-.',linewidth=0.2)
    # save figure
    if len(figureName)>0:
        fig.savefig(figureName, dpi=dpi)
        
    return fig,ax
