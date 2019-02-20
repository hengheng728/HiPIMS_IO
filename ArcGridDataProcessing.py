#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:08:42 2018

@author: b4042552
"""
' a RasterDataProcess module '

__author__ = 'Xiaodong Ming'

#%read ArcGIS ascii file
#%reset -f reset the console and delete all variables
#import linecache
import numpy as np
import glob
import gzip
import math
#import os
def arcgridread(fileName,headrows = 6):
    """
    read ArcGrid format raster file and return the gridded data array, 
    coordinates and cellsize information and the extent of the grid
    """
    try:
        fh = open(fileName, 'r')
        fh.close()
    # Store configuration file values
    except FileNotFoundError:
        # Keep preset values
        print('Error: '+fileName+' does not appear to exist')
        return
# read head
    head = {} # store head information including ncols, nrows,...
    numheadrows = 6
    n=1
    with open(fileName, 'rt') as f:
    # read head
        for line in f:
            if n<=numheadrows:
                line = line.split(" ",1)
                head[line[0]] = float(line[1])
            else:
                break
            n = n+1
# read value array
    gridArray  = np.loadtxt(fileName, skiprows=numheadrows,dtype='float64')
    gridArray[gridArray == head['NODATA_value']] = float('nan')
    head['ncols']=int(head['ncols'])
    head['nrows']=int(head['nrows'])
    left = head['xllcorner']
    right = head['xllcorner']+head['ncols']*head['cellsize']
    bottom = head['yllcorner']
    top = head['yllcorner']+head['nrows']*head['cellsize']
    extent = (left,right,bottom,top)
    #gridArray = float(gridArray)
    return gridArray,head,extent
#%%
def arcgridwrite(fileName,Z,head):    
    Z = Z+0
    Z[np.isnan(Z)]= head['NODATA_value']
    if not isinstance(head,dict):
        raise TypeError('bad argument: head')
    with open(fileName, 'wb') as f:       
        f.write(b"ncols    %d\n" % head['ncols'])
        f.write(b"nrows    %d\n" % head['nrows'])
        f.write(b"xllcorner    %g\n" % head['xllcorner'])
        f.write(b"yllcorner    %g\n" % head['yllcorner'])
        f.write(b"cellsize    %g\n" % head['cellsize'])
        f.write(b"NODATA_value    %g\n" % head['NODATA_value'])    
        np.savetxt(f,Z,fmt='%g', delimiter=' ')
    return None
#%%
def demHead2Extent(demHead):
    # convert dem head file (dict) to a spatial extent of the DEM
    R = demHead
    left = R['xllcorner']
    right = R['xllcorner']+R['ncols']*R['cellsize']
    bottom = R['yllcorner']
    top = R['yllcorner']+R['nrows']*R['cellsize']
    extent = (left,right,bottom,top)
    return extent

#%% shapePoints= makeDiagonalShape(extent)
def makeDiagonalShape(extent):
    #extent = (left,right,bottom,top)
    shapePoints = np.array([[extent[0],extent[2]],
                           [extent[1],extent[2]],
                           [extent[1],extent[3]],
                           [extent[0],extent[3]]])
    return shapePoints
#%% rows,cols = Map2Sub(X,Y,zHead)
def Map2Sub(X,Y,zHead):
    # convert map points to subscripts of a matrix with geo reference zHead
    # x and y coordinate of the centre of the first cell in the matrix
    x0 = zHead['xllcorner']+0.5*zHead['cellsize']
    y0 = zHead['yllcorner']+(zHead['nrows']-0.5)*zHead['cellsize']
    rows = (y0-Y)/zHead['cellsize'] # row and col number starts from 0
    cols = (X-x0)/zHead['cellsize']
    if isinstance(rows,np.ndarray):
        rows = rows.astype('int64')
        cols = cols.astype('int64') #.astype('int64')
    else:
        rows = int(rows)
        cols = int(cols)
    return rows,cols
#%%
def Sub2Map(rows,cols,zHead):
    # convert row and cols in a matrix with geo reference zHead to X and Y coordinates
    # x and y coordinate of the centre of the first cell in the matrix
    x11 = zHead['xllcorner']+0.5*zHead['cellsize']
    y11 = zHead['yllcorner']+(zHead['nrows']+0.5)*zHead['cellsize']
    X = x11+cols*zHead['cellsize']
    Y = y11-rows*zHead['cellsize']
    return X,Y
#%%
def MaskInterp(maskMat,maskHead,zMat,zHead):
    # zMask = MaskInterp(maskMat,maskHead,zMat,zHead)
    rows_Z,cols_Z = np.where(~np.isnan(zMat))
    X,Y = Sub2Map(rows_Z,cols_Z,zHead)
    rows_Mask,cols_mask = Map2Sub(X,Y,maskHead)
    rows_Mask[rows_Mask<0]=0
    cols_mask[cols_mask<0]=0
    rows_Mask[rows_Mask>maskHead['nrows']-1]=maskHead['nrows']-1
    cols_mask[cols_mask>maskHead['ncols']-1]=maskHead['ncols']-1
    values = maskMat[rows_Mask,cols_mask]
    zMask = zMat+0
    zMask[rows_Z,cols_Z]=values
    return zMask
#%% resample grid data to a new resolution via nearest interpolation
def GridResample(zMat,head,newsize):
    """
    resample a grid to a new grid resolution
    """
    if isinstance(newsize, dict):
        head_new = newsize.copy()
    else:            
        head_new = head.copy()
        head_new['cellsize'] = newsize
        ncols = math.floor(head['cellsize']*head['ncols']/newsize)
        nrows = math.floor(head['cellsize']*head['nrows']/newsize)
        head_new['ncols']=ncols
        head_new['nrows']=nrows
    #centre of the first cell in zMat
    x11 = head_new['xllcorner']+0.5*head_new['cellsize']
    y11 = head_new['yllcorner']+(head_new['nrows']-0.5)*head_new['cellsize']
    xAll = np.linspace(x11,x11+(head_new['ncols']-1)*head_new['cellsize'],head_new['ncols'])
    yAll = np.linspace(y11,y11-(head_new['nrows']-1)*head_new['cellsize'],head_new['nrows'])
    rowAll,colAll = Map2Sub(xAll,yAll,head)
    rows_Z,cols_Z = np.meshgrid(rowAll,colAll) # nrows*ncols array
    zNew = zMat[rows_Z,cols_Z]
    zNew = zNew.transpose()
    extent_new = demHead2Extent(head_new)
    return zNew, head_new, extent_new
#%% zMatClip,headClip=ArrayClip(zMat,head,clipExtent)
def ArraySquareClip(zMat,head,clipExtent):
    """
    clip array to a smaller one according to mask and its geoinformation
    extent = (left,right,bottom,top)
    """
    
    X = np.array([clipExtent[0],clipExtent[1]]) # left to right
    Y = np.array([clipExtent[3],clipExtent[2]]) # top to bottom
    rows,cols = Map2Sub(X,Y,head)
    zMatClip = zMat[rows[0]:rows[1]+1,cols[0]:cols[1]+1]
    xllcorner = head['xllcorner']+cols[0]* head['cellsize']
    yllcorner = head['yllcorner'] + (zMat.shape[0]-rows[1]) * head['cellsize']
    headClip = head.copy()
    headClip['xllcorner']= xllcorner
    headClip['yllcorner']= yllcorner
    headClip['nrows'] = zMatClip.shape[0]
    headClip['ncols'] = zMatClip.shape[1]
    extentClip = demHead2Extent(headClip)
    return zMatClip,headClip,extentClip
    
#%% zMask = MaskExtraction(maskMat,maskHead,zHead)
def MaskExtraction(maskMat,maskHead,zHead,maskValue=False):
    """
    extract rainfall mask to model domian with its size and resolution
    # zMask = MaskExtraction(maskMat,maskHead,zHead)
    # maskValue=False:mask value is not given in maskMat, so a mask value mask is to be created
    """
    if ~maskValue:
        maskMat = np.arange(np.size(maskMat)).reshape((maskHead['nrows'],maskHead['ncols']),order='F')
    zMask = np.zeros((zHead['nrows'],zHead['ncols']))
    rows_Z,cols_Z = np.where(~np.isnan(zMask))
    X,Y = Sub2Map(rows_Z,cols_Z,zHead)
    rowsInMask,colsInMask = Map2Sub(X,Y,maskHead)
    
    # make sure rows and cols of domain scells are inside mask 
    rowsInMask[rowsInMask<0]=0
    colsInMask[colsInMask<0]=0
    rowsInMask[rowsInMask>maskHead['nrows']-1]=maskHead['nrows']-1
    colsInMask[colsInMask>maskHead['ncols']-1]=maskHead['ncols']-1
    
    values = maskMat[rowsInMask,colsInMask] # mask values
    zMask[rows_Z,cols_Z]=values
    return zMask
#%% Write and compress asc file
def ArcgridwriteGZip(fileName,Z,head):
    # fileName: compressed file name (automatically added by a suffix '.gz') 
    # example：
    # ArcgridwriteGZip('file.txt',zMat,zHead)
    Z = Z+0
    if not isinstance(head,dict):
        raise TypeError('bad argument: head')
    if fileName[-3:]!='.gz':
        fileName = fileName+'.gz'  
    the_ascfile=gzip.open(fileName, 'wb')        
    the_ascfile.write(b"ncols    %d\n" % head['ncols'])
    the_ascfile.write(b"nrows    %d\n" % head['nrows'])
    the_ascfile.write(b"xllcorner    %g\n" % head['xllcorner'])
    the_ascfile.write(b"yllcorner    %g\n" % head['yllcorner'])
    the_ascfile.write(b"cellsize    %g\n" % head['cellsize'])
    the_ascfile.write(b"NODATA_value    %g\n" % head['NODATA_value'])
    Z[np.isnan(Z)]= head['NODATA_value']
    np.savetxt(the_ascfile,Z,fmt='%.3f', delimiter=' ')
    the_ascfile.close()
    print(fileName)
    return None
#%%
def ArcgridreadGZip(fileName):
    """
    read ArcGrid format raster file and return the gridded data array, 
    coordinates and cellsize information and the extent of the grid
    """
# read head
    head = {} # store head information including ncols, nrows,...
    numheadrows = 6
    n=1
    with gzip.open(fileName, 'rt') as f:
    # read head
        for line in f:
            if n<=numheadrows:
                line = line.split(" ",1)
                head[line[0]] = float(line[1])
            else:
                break
            n = n+1
    gridArray  = np.loadtxt(fileName, skiprows=numheadrows,dtype='float64')
    gridArray[gridArray == head['NODATA_value']] = float('nan')
    head['ncols']=int(head['ncols'])
    head['nrows']=int(head['nrows'])
    left = head['xllcorner']
    right = head['xllcorner']+head['ncols']*head['cellsize']
    bottom = head['yllcorner']
    top = head['yllcorner']+head['nrows']*head['cellsize']
    extent = (left,right,bottom,top)
    #gridArray = float(gridArray)
    return gridArray,head,extent
#%% zG,headG = CombineRaster(inputFolder,outputName=[],fileTag='*.asc')
def CombineRaster(inputFolder,outputName=[],fileTag='*.asc'):
    fileList = glob.glob(inputFolder+'/'+fileTag)
    _,_,extent = arcgridread(fileList[0])
    left = extent[0]
    right = extent[1]
    bottom = extent[2]
    top = extent[3]
    zList = []
    headList = []
    for file in fileList:
        print(file[-17:])
        z,head,extent = arcgridread(file)
        zList.append(z)
        headList.append(head)
        left = min(left,extent[0])
        right = max(right,extent[1])
        bottom = min(bottom,extent[2])
        top = max(top,extent[3])
    # global raster
    headG = head.copy()
    headG['xllcorner'] = left
    headG['yllcorner'] = bottom
    headG['ncols'] = int((right-left)/headG['cellsize'])
    headG['nrows'] = int((top-bottom)/headG['cellsize'])
    zG = np.zeros((headG['nrows'],headG['ncols']))+np.nan
    for i in range(len(fileList)):
        z = zList[i]
        head = headList[i]        
        # centre coordinats of the first element in array z
        x = head['xllcorner']+head['cellsize']/2
        y = head['yllcorner']+head['cellsize']*(head['nrows']-0.5)        
        r0,c0 = Map2Sub(x,y,headG)
        #print(head['nrows'],head['ncols'])        
        zG[r0:r0+int(head['nrows']),c0:c0+int(head['ncols'])] = z
        print(i)
    if len(outputName)>0:
        if outputName[-2:]=='gz':
            ArcgridwriteGZip(outputName,zG,headG)
        else:
            arcgridwrite(outputName,zG,headG)    
    extentG = (left,right,bottom,top)
    return zG,headG,extentG
#%% zNew,headNew=ArcgridReplace(z0,head0,zRe,headRe)
def ArcgridReplace(z0,head0,zRe,headRe):
    """ To replace part of values in z0 with values in zRe
    # the cellsize in head0 and headRe must be the same
    """        
    # centre coordinats of the first element in array z
    if head0['cellsize']!=headRe['cellsize']:
        raise TypeError('the cellsize in head0 and headRe must be the same')
    x0 = headRe['xllcorner']+headRe['cellsize']/2
    y0 = headRe['yllcorner']+headRe['cellsize']*(headRe['nrows']-0.5)
    r0,c0 = Map2Sub(x0,y0,head0)
    rowMin = min(r0,0)
    rowMax = max(r0+headRe['nrows']-1,head0['nrows']-1)
    colMin = min(c0,0)
    colMax = max(c0+headRe['ncols']-1,head0['ncols']-1)

    if rowMin<0:
        Vtop = -rowMin
    else:
        Vtop = 0
    if colMin<0:
        Vleft = -colMin
    else:
        Vleft = 0
    if rowMax>head0['nrows']-1:
        Vbottom = rowMax-(head0['nrows']-1)
    else:
        Vbottom = 0
    if colMax>head0['ncols']-1:
        Vright = colMax-(head0['ncols']-1)
    else:
        Vright = 0
        
    padWidth = [(Vtop,Vbottom),(Vleft,Vright)]
    zNew = np.pad(z0,padWidth,mode='constant',constant_values=head0['NODATA_value'])
    zNew[zNew==head0['NODATA_value']]=np.nan
    headNew = head0.copy()
    headNew['nrows']=zNew.shape[0]
    headNew['ncols']=zNew.shape[1]
    headNew['xllcorner'] = head0['xllcorner']+head0['cellsize']*colMin
    headNew['yllcorner'] = head0['yllcorner']+head0['cellsize']*(head0['ncols']-rowMax)
    r0,c0 = Map2Sub(x0,y0,headNew)
    zRe1 = zNew[r0:r0+headRe['nrows'],c0:c0+headRe['ncols']]
    zRe1[~np.isnan(zRe)]=zRe[~np.isnan(zRe)]
    zNew[r0:r0+headRe['nrows'],c0:c0+headRe['ncols']]=zRe1
    
    return zNew,headNew