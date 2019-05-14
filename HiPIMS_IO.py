#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 12:22:28 2018

@author: Xiaodong Ming
"""
import numpy as np
import sys
import time
from ArcGridDataProcessing import arcgridwrite
import InputSetupFuncs as ISF
# import InputSetupFuncs
# import InputSetupFuncs_MG
import InputSetupFuncs_MG as ISF_MG

def HiPIMS_setup(folderName, demMat, demHead, numSection=1, boundList=[],
               fileToCreate='all', h0=0, hU0=[0,0], manning=0.035,
               rain_mask=0, rain_source=np.array([[0,0],[60,0]]),
               sewer_sink=0, cumulative_depth=0, hydraulic_conductivity=0,
               capillary_head=0, water_content_diff=0,               
               gauges_pos=[], timesValue=[0,3600,600,3600]):
    """
    folderName: a path to store input folder
    fileToCreate: 'all'|'z','h','hU','manning','sewer_sink',
                        'cumulative_depth','hydraulic_conductivity',
                        'capillary_head','water_content_diff'
                        'precipitation_mask','precipitation_source',
                        'boundary_condition','gauges_pos'
                        
    """
    if numSection==1: # call function for single GPU
        InputSetup(folderName, demMat, demHead, boundList, fileToCreate,
                   h0, hU0, manning, rain_mask,rain_source,
                   sewer_sink, cumulative_depth, hydraulic_conductivity,
                   capillary_head, water_content_diff, gauges_pos)
    else: # call function for multi-GPU
        InputSetup_MG(folderName, demMat, demHead, numSection, boundList,
                   fileToCreate, h0, hU0, manning, rain_mask, rain_source,
                   sewer_sink, cumulative_depth, hydraulic_conductivity,
                   capillary_head, water_content_diff, gauges_pos)
    if fileToCreate=='all' or 'device_setup' in fileToCreate:    
        GenDeviceFile(folderName,numSection)
    # start, end, output interval, backup interval, cooldown interval
    if fileToCreate=='all' or 'times_setup' in fileToCreate:
        GenTimeFile(folderName,numSection,Values=timesValue) 
    return None

#%%****************************************************************************

#%================================Single GPU Input=============================


# InputSetup for single GPU
def InputSetup(folderName,demMat,demHead,boundList=[],
               fileToCreate='all', h0=0, hU0=[0,0], manning=0.035,
               rain_mask=0,rain_source = np.array([[0,0],[60,0]]),
               sewer_sink=0, cumulative_depth=0, hydraulic_conductivity=0,
               capillary_head=0, water_content_diff=0,               
               gauges_pos=[]):
    """
    folderName: a path to store input folder
    fileToCreate: 'all'|'z','h','hU','manning',
                        'precipitation_mask','precipitation_source',
                        'boundary_condition','gauges_pos'                        
    """
    start=time.perf_counter()    
    #create input and output folders
    _,_,dirMesh,dirField = ISF.CreateIOFolders(folderName)
#---------- initiate field names and values-------------------    
    fieldFiles = {'z':demMat,'h':h0,'hU':hU0,'precipitation':0, # first ztype file
                      'manning':manning,'sewer_sink':sewer_sink,
                      'cumulative_depth':cumulative_depth,
                      'hydraulic_conductivity':hydraulic_conductivity,
                      'capillary_head':capillary_head,
                      'water_content_diff':water_content_diff, # last ztype file
                      'precipitation_mask':rain_mask, 
                      'precipitation_source':rain_source,
                      'boundary_condition':[],
                      'gauges_pos':gauges_pos}
      
#%------------ generate computing cell ID and define boundary cells-------------    
    # outline cells are classified according to user-defined boundList: 1,2,...
    # and boundary code for h and hU are created
    CellIDMat,OutlineBoundMat = ISF.Get_ID_BNoutline_Mat(demMat)
    BoundTypeMat = ISF.BoundaryClassify(OutlineBoundMat,demHead,boundList) #oultine cells 0,1,2...
    
    # create boundary array for writing files
    hCode,hUCode = ISF.Get3ElementBoundCode(boundList)
    # boundary arrays 
    id_bCodes = ISF.Create_ID_BoundCode_Array(CellIDMat,BoundTypeMat,hCode,hUCode)
    h_BC_source, hU_BC_source = ISF.BoundSourceProcessing(BoundTypeMat,boundList,demHead)
#%******************  create files ********************************
    if not isinstance(fileToCreate,list):        
        if fileToCreate=='all':
            fileToCreate=list(fieldFiles.keys())
            fileToCreate.insert(0,'DEM')
        else:
            fileToCreate = [fileToCreate]
    totalFileNum = len(fileToCreate)
    print('Files to create:')
    print(fileToCreate)
    progress = 1
    fileLeft = totalFileNum
    if fileLeft>13:
        fileLeft=13
    end=time.perf_counter()
    timeElapse = end-start #seconds
#------------ create DEM.txt in mesh-------------------------------------------
    if 'DEM' in fileToCreate:
        fileNameToRite = dirMesh+'DEM.txt'
        updt(totalFileNum, progress, 'DEM',timeElapse*fileLeft*10)
        start = time.perf_counter()
        arcgridwrite(fileNameToRite,demMat,demHead)        
        progress = progress+1
        #print(fileNameToRite+' is created')
        end = time.perf_counter()
        timeElapse = end-start
        fileLeft = fileLeft-1
    for key,value in fieldFiles.items():        
        if key in fileToCreate:
            start = time.perf_counter()
            if len(key)>13:
                showTag = key[0:5]+'..'+key[-5:]
            else:
                showTag = key
            updt(totalFileNum, progress, showTag,timeElapse*fileLeft)
            if key=='precipitation_mask': # time 1
                id_zV = ISF.Create_ID_zValue_Array(CellIDMat,value)
                ISF.WriteRainMask(dirField,id_zV)
            elif key=='precipitation_source': # time 0.5
                ISF.WriteRainSource(dirField,value)
            elif key=='boundary_condition':
                ISF.WriteBoundSource(dirField,h_BC_source,hU_BC_source)
            elif key=='gauges_pos':
                ISF.WriteGaugePos(dirField,value)
            else:
                # create valid cell array for writing files
                id_zV = ISF.Create_ID_zValue_Array(CellIDMat,value)
                if key == 'h':
                    id_bCode = id_bCodes['hArray']
                elif key == 'hU':
                    id_bCode = id_bCodes['hUArray']
                else:
                    id_bCode = id_bCodes['otherArray']
                fileNameToRite = dirField+key+'.dat'                
                ISF.Write_ZtypeFile(fileNameToRite,id_zV,id_BoundCode=id_bCode)                                
            fileLeft=fileLeft-1
            if fileLeft<0:
                fileLeft=0
            end = time.perf_counter()
            timeElapse = end-start            
            #updt(totalFileNum, progress, showTag,timeElapse*fileLeft)
            progress = progress+1
            #print('write_File: '+str(end-start))

    return None
#%%****************************************************************************

#%===============================Multiple GPU Input============================


# InputSetup for multiple GPU
def InputSetup_MG(folderName,demMat,demHead,numSection=1,boundList=[],
               fileToCreate='all', h0=0, hU0=[0,0], manning=0.035,
               rain_mask=0,rain_source = np.array([[0,0],[60,0]]),
               sewer_sink=0, cumulative_depth=0, hydraulic_conductivity=0,
               capillary_head=0, water_content_diff=0,               
               gauges_pos=[]):
    """
    folderName: a path to store input folder
    fileToCreate: 'all'|'z','h','hU','manning',
                        'precipitation_mask','precipitation_source',
                        'boundary_condition','gauges_pos'
                        
    """
    #create input and output folders
    start = time.clock()
    sectionPathList = ISF_MG.CreateIOFolders_MG(folderName,numSection)
#---------- initiate field names and values-------------------    
    fieldFiles = {'z':demMat,'h':h0,'hU':hU0,'precipitation':0, # first ztype file
                      'manning':manning,'sewer_sink':sewer_sink,
                      'cumulative_depth':cumulative_depth,
                      'hydraulic_conductivity':hydraulic_conductivity,
                      'capillary_head':capillary_head,
                      'water_content_diff':water_content_diff, # last ztype file
                      'precipitation_mask':rain_mask, 
                      'precipitation_source':rain_source,
                      'boundary_condition':[],
                      'gauges_pos':gauges_pos}
    # divide grid into small sections
    sectionRowInd,demHeadList = ISF_MG.DivideGrid(demMat,demHead,numSection)

    # define outline cell and inner cell, 0:outline cell; -2: inner cell
    _,OutlineBoundMat = ISF_MG.Get_ID_BNoutline_Mat(demMat)

    # outline cells are classified according to user-defined boundList: 1,2,...
    # and boundary code for h and hU are created 
    BoundTypeMat = ISF.BoundaryClassify(OutlineBoundMat,demHead,boundList) #oultine cells 0,1,2...
    hCode,hUCode = ISF.Get3ElementBoundCode(boundList) #n*3 array, n: number of user-defined boundary
    h_BC_source, hU_BC_source  = ISF.BoundSourceProcessing(BoundTypeMat,boundList,demHead)
    # generate cell ID, boundary type, cell IDs on shared-border for each section
    # border cells are added to the boundary type as -1
    CellIDList,BoundTypeList,shareBorderList = ISF_MG.Get_ID_BN_mat_sec(sectionRowInd,BoundTypeMat)
    ISF_MG.WriteHaloFile(folderName,shareBorderList)
    end = time.clock()
#%*****************************  create files ********************************
    if not isinstance(fileToCreate,list):        
        if fileToCreate=='all':
            fileToCreate=list(fieldFiles.keys())
            fileToCreate.insert(0,'DEM')
        else:
            fileToCreate = [fileToCreate]
    totalFileNum = len(fileToCreate)
    print('Files to create:')
    print(fileToCreate)
    progress = 1
    timeElapse = end-start #seconds
    fileLeft = totalFileNum
    if fileLeft>13:
        fileLeft=13
#------------ create DEM.txt in mesh-------------------------------------------
    if 'DEM' in fileToCreate:
        updt(totalFileNum, progress-1, 'DEM',timeElapse*fileLeft) 
        ISF_MG.WriteSecDEM(sectionPathList,demMat,demHeadList,sectionRowInd)
        #print('DEM '+str(timeElapse))
        fileLeft = totalFileNum-1
        updt(totalFileNum, progress, 'DEM',timeElapse*fileLeft)        
        progress = progress+1
#------------ create ztype, rainfall, and boundary files-----------------------
    for key,value in fieldFiles.items():        
        if key in fileToCreate:
            if len(key)>13:
                showTag = key[0:5]+'..'+key[-5:]
            else:
                showTag = key
            updt(totalFileNum, progress-1, showTag,timeElapse*fileLeft)            
            if key=='precipitation_mask':
                ISF_MG.WriteRainMask_Sec(sectionPathList,CellIDList,value,sectionRowInd)
            elif key=='precipitation_source':
                ISF_MG.WriteRainSource_Sec(sectionPathList,value)
            elif key=='boundary_condition':
                ISF_MG.WriteBoundSource_Sec(sectionPathList,h_BC_source,hU_BC_source)
            elif key=='gauges_pos':
                ISF_MG.WriteGaugesPos_Ind(sectionPathList,demHeadList,value)
            else:
                timeElapse=ISF_MG.WriteZTypeFile_Sec(key,value,CellIDList,BoundTypeList,sectionPathList,sectionRowInd,hCode,hUCode)                        
            fileLeft=fileLeft-1
            #print(key+' '+str(timeElapse))
            if fileLeft<0:
                fileLeft=0
            updt(totalFileNum, progress, showTag,timeElapse*fileLeft)
            progress = progress+1
    return None
#%%**********************************Independent functions*********************
    #functions do not based on DEM data
#%% WriteBoundSource
def WriteRainSource(rootPath,rain_source,numSection):
    if numSection==1: # single GPU
        _,_,_,dirField = ISF.CreateIOFolders(rootPath)
        ISF.WriteRainSource(dirField,rain_source)
    else: # multi-GPU
        sectionPathList = ISF_MG.CreateIOFolders_MG(rootPath,numSection)
        ISF_MG.WriteRainSource_Sec(sectionPathList,rain_source)
    return None
#%% device_setup.dat
def GenDeviceFile(rootPath,numGPU,Values=[]):
    if len(Values)==0:
        Values=np.array(range(numGPU))
    Values = Values.reshape((1,Values.size))
    if numGPU==1:
        np.savetxt(rootPath+'/input/device_setup.dat',Values,fmt='%g')
    else:
        np.savetxt(rootPath+'/device_setup.dat',Values,fmt='%g')
    return None
#%% times_setup.dat
def GenTimeFile(rootPath,numGPU,Values=[0,3600,1800,3600]):
    Values=np.array(Values)
    Values = Values.reshape((1,Values.size))
    if numGPU==1:
        np.savetxt(rootPath+'/input/times_setup.dat',Values,fmt='%g')
    else:
        np.savetxt(rootPath+'/times_setup.dat',Values,fmt='%g')
    return None
#%% Displays or updates a console progress bar
def updt(total, progress, fileTag, timeLeft):
    """
    Displays or updates a console progress bar.
    """
    if total==progress:
        fileTag = 'finished'
    else:
        fileTag = fileTag+'...'    
    barLength, status = 50, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r|{}| {:.0f}% {:<16} time left: {:.0f}s {}".format(
        chr(9608) * block + "-" * (barLength - block), round(progress * 100, 0),
        fileTag,timeLeft,status)
    sys.stdout.write(text)
    sys.stdout.flush()