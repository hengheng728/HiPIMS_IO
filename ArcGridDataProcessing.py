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
#%%
def Map2Sub(X,Y,zHead):
    # convert map points to subscripts of a matrix with geo reference zHead
    # x and y coordinate of the centre of the first cell in the matrix
    x11 = zHead['xllcorner']+0.5*zHead['cellsize']
    y11 = zHead['yllcorner']+(zHead['nrows']+0.5)*zHead['cellsize']
    rows = -(Y-y11)/zHead['cellsize']
    rows = int(rows)-1
    cols = (X-x11)/zHead['cellsize']+1
    cols = int(cols)-1#.astype('int64')
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
        arcgridwrite(outputName,zG,headG)
    extentG = (left,right,bottom,top)
    return zG,headG,extentG