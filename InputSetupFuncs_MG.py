#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HiPIMS input setup functions for Multiple GPUs
Created on Thu Dec  6 12:31:31 2018

@author: Xiaodong Ming
"""
import numpy as np
import time
import os
import shutil
import matplotlib.patches as mplP
from ArcGridDataProcessing import arcgridwrite,demHead2Extent,makeDiagonalShape
from InputSetupFuncs import CreateIOFolders,Get_ID_BNoutline_Mat,\
    WriteRainSource,Write_ZtypeFile,\
    Create_ID_BoundCode_Array,Create_ID_zValue_Array,\
    WriteRainMask,WriteBoundSource,WriteGaugePos,Ztype2Grid
#%%############################################################################
#-------------------------Data processing functions----------------------------
############################################################################### 
#%% Global operation: divide Grid into sections and return rows on border
    # return sectionRowInd
def DivideGrid(demMat,demHeader,numSection):
    Z_valid = ~np.isnan(demMat)
    Z_valid_rowsum = np.sum(Z_valid,axis=1)
    Z_valid_rowCumsum = np.cumsum(Z_valid_rowsum)
    sectionRowInd = np.zeros((numSection,2))
    # divide Z matrix by rows
    rowStart = 0
    exchangeRows = 2
    demSecHeadList = []
    for i in range(numSection):
        numSectionCells=Z_valid_rowCumsum[-1]*(i+1)/numSection
        rowEnd = -1+np.sum(Z_valid_rowCumsum-numSectionCells<=0)
        sectionRowInd[i,0] = rowStart
        if i==numSection: #to cover the last row
            rowEnd = demMat.size-1
        sectionRowInd[i,1] = rowEnd        
        subDemHead = demHeader.copy()
        subDemHead['yllcorner']=demHeader['yllcorner']+(
                                demMat.shape[0]-1-rowEnd
                                )*demHeader['cellsize']
        subDemHead['nrows'] = rowEnd-rowStart+1
        demSecHeadList.append(subDemHead)
        rowStart = rowEnd-(exchangeRows-1)
    sectionRowInd = sectionRowInd.astype('int64')
    return sectionRowInd,demSecHeadList
#%% Global operation: generate section Grid ID and bound type including shared border
#   generate CellIDList, BoundTypeList, and shareBorderList
#   call two sub functions: Get_ID_BNoutline_Mat, genSharedBound     
def Get_ID_BN_mat_sec(sectionRowInd,BoundTypeMat):
    numSection = sectionRowInd.shape[0]
    CellIDList=[] # cell ID based on each sub section
    BoundTypeList= [] # bound type with share-bound added
    shareBorderList=[] # ID of shared bound cells in each section
    # divide the global grid into small sections
    for i in range(numSection):
        n = numSection-(i+1)
        B_sec = BoundTypeMat[sectionRowInd[i,0]:sectionRowInd[i,1]+1,:]#BoundType
        C_sec,_ = Get_ID_BNoutline_Mat(B_sec)#Cell ID for one section
        CellIDList.append(C_sec)        
        if n==0:
            shareTag='top'
        elif n==numSection-1:
            shareTag='bottom'
        else:
            shareTag='both'
        BoundType_sec,shareBorder = GenSharedBound(B_sec,C_sec,shareTag)
        BoundTypeList.append(BoundType_sec)
        shareBorderList.append(shareBorder)
    return CellIDList,BoundTypeList,shareBorderList

#%% Gobal operation
#   convert global zValue to Local zValue
    # called in writeZTypeFile_Sec
def GlobalZvalue2Local(zValue,rowInd):
    # global Zvalue can be a scalar, a matrix or a list of two Zvalues(hU)
    # return a local zValue for a section and keep its original format
    if isinstance(zValue,np.ndarray):
       zV = zValue[rowInd[0]:rowInd[1]+1,:]
    elif isinstance(zValue, list) and len(zValue)==2:
        if np.isscalar(zValue[0]):
            zV0 = zValue[0]
        else:
            zV0 = zValue[0][rowInd[0]:rowInd[1]+1,:]
        if np.isscalar(zValue[1]):
            zV1 = zValue[1]
        else:
            zV1 = zValue[1][rowInd[0]:rowInd[1]+1,:]   
        zV = [zV0,zV1]
    else: #scalar, don't change
        zV = zValue+0
    return zV
#%% Local operation: GenSharedBound for each section
    # called in function DivideGrid
def GenSharedBound(B_sec,C_sec,shareTag='both'):
    #B_sec: Bound type matrix in one section
    #C_sec: Cell ID martix in one section
    #BoundTagArray: add -1 to each domain-shared bound cell
    #shareBorder: topH,topL are the top two lines
    #             bottomH,bottomL are the bottom two lines
    BoundTagArray = B_sec+0
    if shareTag=='bottom':
        line0=[]
        line1=[]
    else:
        ind0 = np.where(~np.isnan(B_sec[0,:]))
        line0 = C_sec[0,ind0]
        ind0 = np.where(B_sec[0,:]==-2)
        BoundTagArray[0,ind0]=-1
        
        ind1 = np.where(~np.isnan(B_sec[1,:]))
        line1 = C_sec[1,ind1]
        ind1 = np.where(B_sec[1,:]==-2)    
        BoundTagArray[1,ind1]=-1
    # domain-shared boundary    
    if shareTag=='top':
        line_2=[]
        line_1=[]
    else:
        ind_2 = np.where(~np.isnan(B_sec[-2,:]))
        line_2 = C_sec[-2,ind_2]
        ind_2 = np.where(B_sec[-2,:]==-2) 
        BoundTagArray[-2,ind_2]=-1
        
        ind_1 = np.where(~np.isnan(B_sec[-1,:]))
        line_1 = C_sec[-1,ind_1]
        ind_1 = np.where(B_sec[-1,:]==-2) 
        BoundTagArray[-1,ind_1]=-1
    shareBorder = {'topH':line0,
                   'topL':line1,
                   'bottomH':line_2,
                   'bottomL':line_1}
    return BoundTagArray,shareBorder


#%%############################################################################
#-----------------------------Writing functions--------------------------------
############################################################################### 
#%% CreateFolders(rootPath,numSection)
def CreateIOFolders_MG(rootPath,numSection):
    sectionPathList=[]
    if rootPath[-1]!='/':
        rootPath = rootPath+'/'
    for i in range(numSection):
        sectionRoot = rootPath+str(i)
        if not os.path.exists(sectionRoot):
            os.makedirs(sectionRoot)
        dirInput,dirOutput,dirMesh,dirField = CreateIOFolders(sectionRoot)
        sectionPath = {'field':dirField,'mesh':dirMesh,
                 'input':dirInput,
                 'output':dirOutput}
        sectionPathList.append(sectionPath)
    return sectionPathList    
    
#%% WriteSecDEM      
def WriteSecDEM(sectionPathList,demMat,demHeaderList,sectionRowInd):
    # write DEM for each section
    start = time.perf_counter()
    numSection = sectionRowInd.shape[0]
    fileTail = 'DEM.txt'    
    for i in range(numSection):
        n = numSection-(i+1)
        Z_sec = demMat[sectionRowInd[i,0]:sectionRowInd[i,1]+1,:]#DEM
        subDemHead = demHeaderList[i]
        fileName = sectionPathList[n]['mesh']+fileTail
        arcgridwrite(fileName,Z_sec,subDemHead)
    #print(fileTail+' is created')
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse
#%% WriteHaloFile       
def WriteHaloFile(rootPath,shareBorderList):
    numSection = len(shareBorderList)    
    if rootPath[-1]!='/':
        rootPath = rootPath+'/'
    fileName = rootPath+'halo.dat'
    with open(fileName,'w') as fileTowrite:
        fileTowrite.write("No. of Domains\n")
        fileTowrite.write("%d\n" % numSection)
        for domainID in range(numSection):
            i = numSection-(domainID+1)            
            fileTowrite.write("#%d\n" % domainID)
            shareBorder = shareBorderList[i]
            for line in ['bottomL','bottomH','topH','topL']:
                if len(shareBorder[line])==0:
                    fileTowrite.write(' \n')
                else:
                    np.savetxt(fileTowrite,shareBorder[line],fmt='%d', delimiter=' ')
    return None
#%% Write zType file to each section
def WriteZTypeFile_Sec(fileTag,zValue,CellIDList,BoundTypeList,sectionPathList,sectionRowInd,hCode,hUCode):
    # write Z-type file for each section
    # fileTag: 'h', 'hU', 'manning',...
    # four function called:
    #   globalZvalue2Local, create_ID_zValue_Arr, 
    #   create_ID_BoundCode_Array, write_ZtypeFile
    start = time.perf_counter()
    numSection = len(sectionPathList)    
    for domainID in range(numSection):
        i = numSection-(domainID+1)
        CellID = CellIDList[i]
        BoundType = BoundTypeList[i]
        rowInd =  sectionRowInd[i,:] 
        zValue_sec = GlobalZvalue2Local(zValue,rowInd)                                
        # generate boudary array for different z type files
        id_bCodes = Create_ID_BoundCode_Array(CellID,BoundType,hCode,hUCode)
        if fileTag == 'h':
            id_bCode = id_bCodes['hArray']
        elif fileTag == 'hU':
            id_bCode = id_bCodes['hUArray']
            if np.isscalar(zValue_sec):
                zValue_sec = [zValue_sec,zValue_sec]
        else:
            id_bCode = id_bCodes['otherArray']
        # generate cell ID-value array for different z type files
        id_zv = Create_ID_zValue_Array(CellID,zValue_sec)
        fileTail = fileTag+'.dat'
        writeName = sectionPathList[domainID]['field']+fileTail
        Write_ZtypeFile(writeName,id_zv,id_bCode)
    #print(fileTail+' is created')
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse

#%% writeRainMask_Sec
def WriteRainMask_Sec(sectionPathList,CellIDList,rain_mask,sectionRowInd):
    # rain_mask can be a scalar or a matrix with the same size of DEM
    start = time.perf_counter()
    numSection = len(CellIDList)
    for domainID in range(numSection):
        i = numSection-(domainID+1)
        CellID = CellIDList[i]
        rowInd =  sectionRowInd[i,:] 
        rainMask_sec = GlobalZvalue2Local(rain_mask,rowInd)
        id_zV = Create_ID_zValue_Array(CellID,zValue=rainMask_sec)
        WriteRainMask(sectionPathList[domainID]['field'],id_zV)
    #print('precipitation_mask.dat is created')
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse
#%% write Rain Source in each section
def WriteRainSource_Sec(sectionPathList,rain_source):
    # create one file in the first folder and copy to other folders
    start = time.perf_counter()
    numSection = len(sectionPathList)
    fileFolder = sectionPathList[0]['field']
    fileName = WriteRainSource(fileFolder,rain_source)
    if isinstance(fileName, list):
        for i in range(1,numSection):
            for file in fileName:
                shutil.copy(file,sectionPathList[i]['field'])
    else:
       for i in range(1,numSection):
           shutil.copy(fileName,sectionPathList[i]['field'])
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse
#%% write Boundary source
def WriteBoundSource_Sec(sectionPathList,h_BC_source,hU_BC_source):
    start = time.perf_counter()
    numSection = len(sectionPathList)
    fileFolder = sectionPathList[0]['field']
    fileName_h,fileName_hU=WriteBoundSource(fileFolder,h_BC_source,hU_BC_source)
    for file in fileName_h:
        for i in range(1,numSection):
            shutil.copy(file,sectionPathList[i]['field'])  
    for file in fileName_hU:
        for i in range(1,numSection):
            shutil.copy(file,sectionPathList[i]['field'])
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse
#%% write gauges postion and index
def WriteGaugesPos_Ind(sectionPathList,demHeadList,gauges_pos):
    # write gauges postion and index to each section 
    #   two function called:
    #       demHead2Extent,makeDiagonalShape
    start = time.perf_counter()
    if len(gauges_pos)==0:
        gauges_pos = np.array([0,0],ndmin=2)
    numSection = len(sectionPathList)    
    for domainID in range(numSection):
        i = numSection-(domainID+1)
        extentSec = demHead2Extent(demHeadList[i])
        shapePoints= makeDiagonalShape(extentSec)
        poly = mplP.Polygon(shapePoints, closed=True)
        ind = np.where(poly.contains_points(gauges_pos))  
        ind = ind[0]
        if ind.size==0:
            gauges_pos_Sec = np.array([0,0],ndmin=2)
            ind = np.array([-1]) # if no points are located in the sub-section
        else:
            gauges_pos_Sec = gauges_pos[ind,:]
            ind = ind.astype('int64')
        # gauges_pos.dat    
        WriteGaugePos(sectionPathList[domainID]['field'],gauges_pos_Sec)
        fileTail2 = 'gauges_ind.dat'
        fileName = sectionPathList[domainID]['field']+fileTail2
        with open(fileName,'w') as file2write:
            np.savetxt(file2write,ind,fmt='%d',delimiter=' ')
    end = time.perf_counter()
    timeElapse = end-start
    return timeElapse
#%% Convert ZtypeFiles to Grid file in multiGPU folder
def Ztype2Grid_MG(caseFolder,numSection,fileTag):
    #Convert ZtypeFiles to Grid file in multiGPU folder
    #fileTag: 'z','h',
    if caseFolder[-1]!='/':
        caseFolder=caseFolder+'/' 
    zHeadGlobal = []
    for i in range(numSection):
        rootPath = caseFolder+str(i)+'/'
        zMat,zHead,_ = Ztype2Grid(rootPath,fileTag)
        if fileTag=='hU':
            if i==0:
                zMatGlobal0=zMat[0]
                zMatGlobal1=zMat[1]
                zHeadGlobal = zHead
            else:
                zMatGlobal0 = np.r_[zMat[0][0:-2,:],zMatGlobal0]
                zMatGlobal1 = np.r_[zMat[1][0:-2,:],zMatGlobal1]
            zMatGlobal = [zMatGlobal0,zMatGlobal1]
        else:
            if i==0:
                zMatGlobal=zMat
                zHeadGlobal = zHead
            else:
                zMatGlobal = np.r_[zMat[0:-2,:],zMatGlobal]
    if isinstance(zMatGlobal,list):
        row,col = zMatGlobal[0].shape
    else:
        row,col = zMatGlobal.shape
    zHeadGlobal['ncols'] = col
    zHeadGlobal['nrows'] = row
    left = zHeadGlobal['xllcorner']
    right = zHeadGlobal['xllcorner']+zHeadGlobal['ncols']*zHeadGlobal['cellsize']
    bottom = zHeadGlobal['yllcorner']
    top = zHeadGlobal['yllcorner']+zHeadGlobal['nrows']*zHeadGlobal['cellsize']
    zExtentGlobal = (left,right,bottom,top)
    return zMatGlobal,zHeadGlobal,zExtentGlobal
    