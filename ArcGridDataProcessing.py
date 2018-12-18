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
    rows = rows.astype('int64')
    cols = (X-x11)/zHead['cellsize']+1
    cols = cols.astype('int64')
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