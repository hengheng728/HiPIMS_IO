#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 16:35:40 2018

@author: Xiaodong Ming
"""
import numpy as np
import pandas as pd
import shutil
import glob
import gzip
import os    
from ArcGridDataProcessing import arcgridread

#%% Combine Grid files from Multiple GPU outputs:
def CombineGridFile(rootPath,numSection,fileTag):
    """
    demMatGlobal,demHeadGlobal,extentGlobal=CombineGridFile(rootPath,numSection,fileTag)
    Input:
    rootPath: directory of case input files containing domain ID 0,1,... folder
    numSection: number of subsection domains
    fileTag: DEM, h_TTTT, h_max_TTTT_max
    Return:
    gauges_ind_pos_All: 3-col array representing gauge index, X, and Y coordinates
    gaugeArray_All: gauge value array with recorded time in the first column
    """
    if 'DEM' in fileTag:
        fileTail = '/input/mesh/DEM.txt'
    else:
        fileTail = '/output/' + fileTag + '.asc'
    if rootPath[-1]!='/':
        rootPath = rootPath+'/'
    demMatGlobal = []
    demHeadGlobal = []
    for i in range(numSection):
        fileName = rootPath+str(i)+fileTail
        if i==0:
            demMatGlobal,demHeadGlobal,_ = arcgridread(fileName)
        else:
            demMat,demHead,_ = arcgridread(fileName)
            demMatGlobal = np.r_[demMat[0:-2,:],demMatGlobal]
    row,col = demMatGlobal.shape
    demHeadGlobal['ncols'] = col
    demHeadGlobal['nrows'] = row
    left = demHeadGlobal['xllcorner']
    right = demHeadGlobal['xllcorner']+demHeadGlobal['ncols']*demHeadGlobal['cellsize']
    bottom = demHeadGlobal['yllcorner']
    top = demHeadGlobal['yllcorner']+demHeadGlobal['nrows']*demHeadGlobal['cellsize']
    extentGlobal = (left,right,bottom,top)
    return demMatGlobal,demHeadGlobal,extentGlobal


#%% combine gauge files from Multiple GPU outputs:
def CombineGaugeFile(rootPath,numSection,fileTag,writeFolder=False):
    """
    Input:
    rootPath: directory of case input files containing domain ID 0,1,... folder
    numSection: number of subsection domains
    fileTag: h_gauges, hU_gauges, eta_gauges
    Return:
    gauges_ind_pos_All: 3-col array representing gauge index, X, and Y coordinates
    gaugeArray_All: gauge value array with recorded time in the first column
    or
    gaugeGeoTable
    gaugeTimeSeries
    """
    if rootPath[-1]!='/':
            rootPath = rootPath+'/'    
    gaugeFileTail = '/output/' + fileTag + '.dat'
    for i in range(numSection):
        fileName = rootPath+str(i)+'/input/field/gauges_ind.dat'
        gauges_ind = np.loadtxt(fileName,dtype='float64')
        fileName = rootPath+str(i)+'/input/field/gauges_pos.dat'
        gauges_pos = np.loadtxt(fileName,dtype='float64',ndmin=2)
        gauges_ind_pos = np.c_[gauges_ind,gauges_pos]
        fileName = rootPath+str(i)+gaugeFileTail
        gaugeArray  = np.loadtxt(fileName,dtype='float64')
        numGauge = gauges_ind.size
        if i==0:
            gauges_ind_pos_All = gauges_ind_pos+0
            if 'hU' in fileTag:
                gaugeArray_All = gaugeArray[:,0:numGauge*2+1]
            else:
                gaugeArray_All = gaugeArray[:,0:numGauge+1]
        else:
            gauges_ind_pos_All = np.r_[gauges_ind_pos_All,gauges_ind_pos]
            if 'hU' in fileTag:
                gaugeArray_All = np.c_[gaugeArray_All,gaugeArray[:,1:numGauge*2+1]]
            else:
                gaugeArray_All = np.c_[gaugeArray_All,gaugeArray[:,1:numGauge+1]]
    gaugeGeoTable = pd.DataFrame(gauges_ind_pos_All[:,1:],\
                                     index = gauges_ind_pos_All[:,0].astype('int'),\
                                     columns=['GaugeX','GaugeY'])
    gaugeIndex = gaugeGeoTable.index.values
    columnsName = ['Time'] 
    for i in gaugeIndex:
        if 'hU' in fileTag:
            columnsName.append('Gauge_' + str(i)+'_u')
            columnsName.append('Gauge_' + str(i)+'_v')
        else:
            columnsName.append('Gauge_' + str(i))
    gaugeTimeSeries = pd.DataFrame(gaugeArray_All,columns=columnsName)
    if writeFolder==True:
        gaugeGeoTable.to_csv(rootPath+'gaugeGeo.dat',sep=' ')
        gaugeTimeSeries.to_csv(rootPath+fileTag+'TimeSeries.dat',sep=' ',index=False)
    return gaugeGeoTable,gaugeTimeSeries

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
#%% delete output file or files
def DeleteMultiOutputFiles(rootPath,numSec,fileStr):
    """
    example:
        DeleteMultiOutputFiles(rootPath,numSec,'h_0.dat')
        DeleteMultiOutputFiles(rootPath,numSec,'*.txt')
        DeleteMultiOutputFiles(rootPath,numSec,'*_backup_*')
        DeleteMultiOutputFiles(rootPath,numSec,'*') delete all files
    """
    if rootPath[-1]!='/':
            rootPath = rootPath+'/'
    for i in range(numSec):
        fileToRemove = rootPath+str(i)+'/output/'+fileStr
        fileToRemove = glob.glob(fileToRemove)
        for file in fileToRemove:
            os.remove(file)
    return None

#%% copy backup files as initial condition files
def Backup2Initial(rootPath,numSection,Time,remove=False):
    if rootPath[-1]!='/':
            rootPath = rootPath+'/'
    for i in range(numSection):
        outputDir = rootPath+str(i)+'/output/'
        inputDir = rootPath+str(i)+'/input/field/'
        fileBackup = outputDir+'h_backup__'+str(Time)+'.dat'
        fileInitial = inputDir+'h.dat'
        print(fileBackup)
        print(fileInitial)
        shutil.copy(fileBackup,fileInitial)
        if remove==True:
            os.remove(fileBackup)
        fileBackup = outputDir+'hU_backup__'+str(Time)+'.dat'
        fileInitial = inputDir+'hU.dat'
        shutil.copy(fileBackup,fileInitial)
        if remove==True:
            os.remove(fileBackup)
    return None