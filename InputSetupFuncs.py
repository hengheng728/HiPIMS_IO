#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:55:15 2018

@author: b4042552
"""
import numpy as np
import os
import pandas as pd
#import matplotlib as mpl
#import time
import matplotlib.patches as mplP
import scipy.signal
from ArcGridDataProcessing import arcgridread,makeDiagonalShape,demHead2Extent
"""
Key varaible names:

    demMat: a matrix of the DEM grid with nan
    demHead: a dictionary with spatial reference infromation of DEM grid
    
    boundList: a list of dict, each dict contain keys (polyPoints,type,h,hU) to 
        define a boundary including position, type, and IO sources 
        1.polyPoints is a numpy array giving X(1st col) and Y(2nd col) of points 
          to define the position of a boundary. An empty polyPoints means outline boundary.
        2.type: 'open'|'rigid'
        3.h: a two-col numpy array, 1st col is time(s), 2nd col is water depth(m)
        4.hU: a two-col numpy array, 1st col is time(s), 2nd col is discharge(m3/s)
             or a three-col numpy array, 2nd col and 3rd col are velocities(m/s) 
             in x and y direction, respectively
    
    ID: (0,1,2,...) ID of each valid cell in DEM grid excluding nan
    BN: (-1,0,1,2,...) Boundary Number, an intermediate variable to classify
        different boundary types: -1:section-shared cell; 0: outline cell;
        1,~: user-defined IO bound cell on the outline
    idMat: a martix of ID with the same size of demMat
    bnMat: a martix of BN with the same size of demMat,-2 means non-bound cell
    bnMat_outline: a martix of BN with the same size of demMat,-2 means non-bound cell
    

"""
#%% Global and local operation: generate cell ID and outline boundary
def Get_ID_BNoutline_Mat(demMat):
    # cellID, bnMat_outline = Get_ID_BNoutline_Mat(demMat)
    Z = demMat*0+1
    Z_flip = np.flipud(Z)
    D1 = np.nancumsum(Z_flip)
    Z_flip1D = np.reshape(Z_flip,np.shape(D1))
    D1[np.isnan(Z_flip1D)]=np.nan
    D1 = np.reshape(D1,np.shape(Z_flip))
    # Series number of valid cells: 0 to number of cells-1
    # from bottom left corner towards right and top
    idMat = np.flipud(D1)-1
    del Z,Z_flip,D1,Z_flip1D
    D = idMat*0
    D[np.isnan(idMat)]=-1
    h_hv   = np.array([[0,1,0], [1,0,1], [0,1,0]])
    D = scipy.signal.convolve2d(D, h_hv,mode='same')
    D[D<0]=np.nan    
    D[0,:] = np.nan
    D[-1,:] = np.nan
    D[:,0] = np.nan
    D[:,-1] = np.nan
    # boundary cells with valid cell ID are extracted
    bnMat_outline = idMat*0-2 # non-outline cells:-2
    Outline_Cell_index = np.isnan(D)&~np.isnan(idMat)
        # outline boundary cell
    bnMat_outline[Outline_Cell_index] = 0 # outline cells:0
    return idMat,bnMat_outline
#%% Global operation: classify outline boundary cells into different types
    # shared border is not considered yet in this step
def BoundaryClassify(bnMat_outline,demHead,boundList=[]):
    #bnMat = BoundaryClassify(bnMat_outline,demHead,boundList)
    #%get rows and columns of outline bound cells
    # bnMat: nan: invalida cell; -2: non-bound cell; 0: outline cell;
    #               1,~: user-defined IO bound cell on the outline
    bnMat = bnMat_outline
    R = demHead
    Extent = demHead2Extent(R)
    BoundSubs = np.where(bnMat_outline==0)
    Bound_Cell_X = R['xllcorner']+(BoundSubs[1]+0.5)*R['cellsize']
    Bound_Cell_Y = R['yllcorner']+(R['nrows']-BoundSubs[0]-0.5)*R['cellsize']    
    BoundSubs = np.array([BoundSubs[0],BoundSubs[1]])
    BoundSubs = np.transpose(BoundSubs)
    n=1 # sequence number of boundaries
    for dBound in boundList:        
        if len(dBound['polyPoints'])==0: #outline boundary
            polyPoints = makeDiagonalShape(Extent)
        elif len(dBound['polyPoints'])==2:
            xyv = dBound['polyPoints']
            polyPoints = makeDiagonalShape([np.min(xyv[:,0]),
                                            np.max(xyv[:,0]),
                                            np.min(xyv[:,1]),
                                            np.max(xyv[:,1])])
        else:
            polyPoints = dBound['polyPoints']
        
        poly = mplP.Polygon(polyPoints, closed=True)
        Bound_Cell_XY = np.array([Bound_Cell_X,Bound_Cell_Y])
        Bound_Cell_XY = np.transpose(Bound_Cell_XY)
        ind1 = poly.contains_points(Bound_Cell_XY)
        row = BoundSubs[ind1,0]
        col = BoundSubs[ind1,1]
        bnMat[row,col]=n
        #BoundNum[ind1,0] = n+1
        n=n+1
    return bnMat
#%% Global operation: create the 3-element bound code for IO boundary
    # return hCode and hUCode    
def Get3ElementBoundCode(boundList):
    #create the 3-element bound code for each bound type(user-defined boundary)
    # create bound code for each IO boundary according to BounList
    # Input: boundList
    # Output: hCode,hUCode, n*3 array, n is the number of boundaries 
    if len(boundList)==0:
        boundList= [{'polyPoints': [],'type': 'open'}] #default outline boundary    
    hCode = np.zeros((len(boundList),3))
    hCode[:,0] = 2
    hUCode = hCode+0    
    n=0 #sequence of boundary
    m_h=0 #sequence of boundary with IO source files  
    m_hU=0
    for dBound in boundList:        
        if dBound['type']=='rigid':
            hUCode[n,1] = 2 #[2 2 0]
        elif dBound['type']=='open':
            hUCode[n,1] = 1
            if 'h' in dBound:
                hCode[n,:] = [3,0,m_h] #[3 0 m]
                m_h=m_h+1
            if 'hU' in dBound:
                hUCode[n,:] = [3,0,m_hU] #[3 0 m]
                m_hU=m_hU+1
        n=n+1
    return hCode,hUCode
#%% Local operation: Create_ID_BoundCode_Array for each section
    # called in writeZTypeFile_Sec
def Create_ID_BoundCode_Array(CellID,BoundType,hCode,hUCode):
    """
    id_bCode = Create_ID_BoundCode_Array(CellID,BoundType,hCode,hUCode)
    CellID: Cell ID matrix 0,1,2,...
    BoundType: Bound type number matrix -1[share],0[outline],1,...[user-defined]
    # hCode,hUCode: 3-col arrays to store boundary code for each boundary type
    # id_bCode: return a dict consisting of three numpy arrays:
    # hArray: Bound cell ID and boundary codes for h
    # hUArray: Bound cell ID and boundary codes for hU
    # otherArray: Bound cell ID and boundary codes for other files
    """
    np.warnings.filterwarnings('ignore')
    ind = BoundType>=-1
    #id_BN: 1.cellID-2.BoundTypeNum
    #    location and IDs of each boundary cell and its boundary numbers 
    #               corresponding to boundList
    id_BN = np.c_[CellID[ind],BoundType[ind]]
    id_BN = id_BN[id_BN[:,0].argsort()]
    BoundNum = id_BN[:,1]
    id_b_hCode = np.zeros((len(BoundNum),3))    
    id_b_hCode[:,0] = 2 # default open h [2, 0, 0]
    id_b_hCode[BoundNum==-1,0] = 4 # share bound [4,0,0] 
    id_b_hUCode = id_b_hCode+0
    id_b_hUCode[:,1] = 1 # default open hU [2, 1, 0]
    id_b_hUCode[BoundNum==-1,0] = 4 # share bound [4,0,0]
    id_b_hUCode[BoundNum==-1,1] = 0
    id_b_oCode = id_b_hCode+0
    id_b_oCode[:,0] = 2 # default for other files [2, 0, 0]
    id_b_oCode[BoundNum==-1,0] = 4 # share bound [4,0,0] 
    
    for k in range(0,hCode.shape[0]):        
        id_b_hCode[BoundNum==k+1,:] = hCode[k,:]
    for k in range(0,hUCode.shape[0]):        
        id_b_hUCode[BoundNum==k+1,:] = hUCode[k,:]
    id_b_hCode  = np.c_[id_BN[:,0],id_b_hCode]
    id_b_hUCode = np.c_[id_BN[:,0],id_b_hUCode]
    id_b_oCode  = np.c_[id_BN[:,0],id_b_oCode]    
    id_bCode = {'hArray':id_b_hCode,'hUArray':id_b_hUCode,'otherArray':id_b_oCode}
    return id_bCode
#%% Global or local operation
    # called in writeZTypeFile_Sec
def Create_ID_zValue_Array(Cell_ID_mat,zValue=0):
    """
    Cell_ID_mat: Cell ID matrix
    zValue: a scalar, a numpy array with the same shape of DEM
    or a list of two elements representing u and v
    """
    ind = np.where(~np.isnan(Cell_ID_mat))
    #Cell_ID_Subs: 1st col is valid cell ID, 
    #              2nd & 3rd cols are rows and cols in DEM
    Cell_ID_Subs = np.c_[Cell_ID_mat[ind],ind[0],ind[1]]
    Cell_ID_Subs = Cell_ID_Subs[Cell_ID_Subs[:,0].argsort()]
    ids = Cell_ID_Subs[:,0]
    subs = Cell_ID_Subs[:,[1,2]]
    subs = subs.astype('int64')
    if np.isscalar(zValue):
        zV = np.zeros(ids.shape)+zValue
    elif isinstance(zValue, list) and len(zValue)==2: # [u,v]
        if np.isscalar(zValue[0]):
            zV0 = np.zeros(ids.shape)+zValue[0]
        else:
            zV0 = zValue[0][subs[:,0],subs[:,1]]
        if np.isscalar(zValue[1]):
            zV1 = np.zeros(ids.shape)+zValue[1]
        else:
            zV1 = zValue[1][subs[:,0],subs[:,1]]    
        zV = np.c_[zV0,zV1]
    else: #numpy array
        zV = zValue[subs[:,0],subs[:,1]]
    id_zV = np.c_[ids,zV]
    return id_zV
#%% Global Operation: prepare boundary source    
def BoundSourceProcessing(bnMat,boundList,demHead):
    #h_BC_source, hU_BC_source  = BoundSourceProcessing(bnMat,boundList,demHead)
    # prepare h and hU boundary sources to the writable format
    # h_BC_source: two column-array
    # hU_BC_source: three column-array
    np.warnings.filterwarnings('ignore')
    subs = np.where(bnMat>0) # user-defined boundary
    id_BN = np.c_[subs[0]*0,bnMat[bnMat>0],subs[0],subs[1]]
    cellsize = demHead['cellsize']
    h_BC_source = []
    hU_BC_source = []
    n = 1 # sequence of boundary
    for dBound in boundList:        
        if 'h' in dBound:
            if len(dBound['h'])==0:
                h_BC_source.append(np.array([[0,0],[3600,0]]))
            else:
                h_BC_source.append(np.array(dBound['h'],ndmin=2))
        if 'hU' in dBound:
            if len(dBound['hU'])==0:
                hU_BC_array = np.array([[0,0,0],[3600,0,0]])
            else:
                sourceArray = np.array(dBound['hU'],ndmin=2)
                if np.shape(sourceArray)[1]==2:
                    Q = sourceArray[:,1] # discharge given, need to transfered to velocity
                    indBound = id_BN[:,1]==n                    
                    rows = id_BN[indBound,2]
                    cols = id_BN[indBound,3]
                    boundLength = cellsize*len(rows)
                    rowsRange = np.ptp(rows)
                    colsRange = np.ptp(cols)
                    hypotenRange = np.sqrt(colsRange**2+rowsRange**2)
                    cosV = colsRange/hypotenRange
                    sinV = rowsRange/hypotenRange
                    hUx = Q*sinV/boundLength
                    hUy = Q*cosV/boundLength
                    hU_BC_array = np.c_[sourceArray[:,0],hUx,hUy] #velocity                    
                else:
                    hU_BC_array = np.array(dBound['hU'],ndmin=2)
            hU_BC_source.append(hU_BC_array)
        n=n+1
    return h_BC_source, hU_BC_source                                        
#%%==============================write file====================================
    #%% create IO Folders for single GPU
def CreateIOFolders(folderName):
    #dirInput,dirOutput,dirMesh,dirField = CreateIOFolders(folderName)
    if folderName[-1]!='/':
        folderName = folderName+'/'        
    dirInput = folderName+'input/'
    dirOutput = folderName+'output/'
    if not os.path.exists(dirOutput):
        os.makedirs(dirOutput)
    if not os.path.exists(dirInput):
        os.makedirs(dirInput)
    dirMesh = dirInput+'mesh/'
    if not os.path.exists(dirMesh):
        os.makedirs(dirMesh)
    dirField = dirInput+'field/'
    if not os.path.exists(dirField):
        os.makedirs(dirField)
    return dirInput,dirOutput,dirMesh,dirField
#%% write Rain Mask
def WriteRainMask(fileFolder,id_zV):
    fileName = fileFolder+'precipitation_mask.dat'
    fmt = ['%12d %d']
    fmt = '\n'.join(fmt*id_zV.shape[0])
    id_zV_Str = fmt % tuple(id_zV.ravel()) 
    with open(fileName, 'w') as file2write:
        file2write.write("$Element Number\n")
        file2write.write("%d\n" % id_zV.shape[0])
        file2write.write("$Element_id  Value\n")
        file2write.write(id_zV_Str)
#        np.savetxt(file2write,id_zV,fmt=fmt,delimiter=' ')
    #print(fileName+' is created')
    return None
#%% write Rainfall Source
def WriteRainSource(fileFolder,rain_source):
    #fileName = WriteRainSource(fileFolder,rain_source)
    fmt1 = '%g'
    fmt2 = '%.8e'
    if isinstance(rain_source, np.ndarray): # write precipitation_source_all
        numSource = rain_source.shape[1]-1
        formatSpec = [fmt2]*numSource
        formatSpec.insert(0,fmt1)
        fileTail = 'precipitation_source_all.dat'
        fileName = fileFolder+fileTail
        with open(fileName,'w') as file2write:
            file2write.write("%d\n" % numSource)
            np.savetxt(file2write,rain_source,fmt=formatSpec,delimiter=' ')
        #print(fileName+' created')
    elif isinstance(rain_source, list):
        fileName = []
        N=0
        for tRain in rain_source:
            fileTail = 'precipitation_source_'+str(N)+'.dat'
            fileName1 = fileFolder+fileTail
            N = N+1
            with open(fileName1,'w') as file2write:
                np.savetxt(file2write,tRain,fmt=[fmt1,fmt2],delimiter=' ')
            fileName.append(fileName1)
        #print(fileName[0]+'~'+fileName1+' created')
    return fileName
#%% Universal write Write_ZtypeFile
#   Secondary function, called in writeZTypeFile_Sec                    
def Write_ZtypeFile(fileName,id_zV,id_BoundCode):
    # write z type files
    # add a 0.00001 to boundary with IO file
    fmt = ['%d %g']
    if fileName[-4:]!='.dat':
        fileName = fileName+'.dat'
    if fileName[-5:]=='h.dat':        
        ind=id_BoundCode[id_BoundCode[:,1]==3,0]
        ind = ind.astype('uint32')
        if len(ind)>0:
            if np.max(id_zV[ind,1])==0:
                id_zV[ind,1] = id_zV[ind,1]+0.0001
    if fileName[-6:]=='hU.dat':
        fmt = ['%d %.8e %.8e']
        ind = id_BoundCode[id_BoundCode[:,1]==3,0]   
        ind = ind.astype('uint32')
        if len(ind)>0:
            if np.max(id_zV[ind,1])==0:
                id_zV[ind,1] = id_zV[ind,1]+0.0001

    fmt = '\n'.join(fmt*id_zV.shape[0])
    id_zV_Str = fmt % tuple(id_zV.ravel())
    fmt=['%-12d %2d %2d %2d']
    fmt = '\n'.join(fmt*id_BoundCode.shape[0])
    id_BoundCode_Str = fmt % tuple(id_BoundCode.ravel()) 
    with open(fileName, 'w') as file2write:
        file2write.write("$Element Number\n")
        file2write.write("%d\n" % id_zV.shape[0])
        file2write.write("$Element_id  Value\n")
        file2write.write(id_zV_Str)
        file2write.write("\n$Boundary Numbers\n")
        file2write.write("%d\n" % id_BoundCode.shape[0])
        file2write.write("$Element_id  Value\n") 
        file2write.write(id_BoundCode_Str)
    return None
#%% write Boundary source
def WriteBoundSource(fileFolder,h_BC_source,hU_BC_source):
    #fileName_h,fileName_hU=WriteBoundSource(fileFolder,h_BC_source,hU_BC_source)
    fmt_h  = ['%g','%g']
    fmt_hU = ['%g','%.8e','%.8e']
    m_h = 0
    m_hU = 0
    fileName_h = []
    fileName_hU = []
    for sourceValue in h_BC_source:            
        fileName = fileFolder+'h_BC_'+str(m_h)+'.dat'
        np.savetxt(fileName,sourceValue,fmt=fmt_h,delimiter=' ')
        fileName_h.append(fileName)
        #print(fileName+' is created')
        m_h=m_h+1
    for sourceValue in hU_BC_source:            
        fileName = fileFolder+'hU_BC_'+str(m_hU)+'.dat'
        np.savetxt(fileName,sourceValue,fmt=fmt_hU,delimiter=' ')
        fileName_hU.append(fileName)
        #print(fileName+' is created')
        m_hU=m_hU+1
    return fileName_h,fileName_hU
#%% write Gauge Position
def WriteGaugePos(fileFolder,gauges_pos):
    if len(gauges_pos)==0:
        gauges_pos = np.array([0,0],ndmin=2)
    fileName = fileFolder+'gauges_pos.dat'
    fmt=['%g %g']
    fmt = '\n'.join(fmt*gauges_pos.shape[0])
    gauges_pos_Str = fmt % tuple(gauges_pos.ravel()) 
    with open(fileName,'w') as file2write:
        file2write.write(gauges_pos_Str)
#        np.savetxt(file2write,gauges_pos,fmt=fmt,delimiter=' ')
    #print(fileName+' is created')
#%%    
def Ztype2Grid(rootPath,fileTag):
    # convert ztype field file to a grid file
    # zMat,zHead,zExtent = Ztype2Grid(rootPath,fileTag)
    # require the existance of DEM.txt file in input/mesh
    # example: zMat,zHead,zExtent = Ztype2Grid('CaseP/input','h')
    if rootPath[-1]!='/':
        rootPath=rootPath+'/'
    inputFolder = rootPath+'input/'    
    demFile = inputFolder+'mesh/DEM.txt'
    zMat,zHead,zExtent = arcgridread(demFile)
    zTypeFile = inputFolder+'field/' +fileTag+'.dat'
    with open(zTypeFile,'r') as f:    
        headlines = [next(f) for x in range(2)]
        N = int(headlines[1])
        vecArray=pd.read_csv(f,delim_whitespace=True ,nrows=N)
    cellID,_ = Get_ID_BNoutline_Mat(zMat)
    subs = np.where(~np.isnan(cellID))
    ind = np.where(~np.isnan(cellID))
    #Cell_ID_Subs: 1st col is valid cell ID, 
    #              2nd & 3rd cols are rows and cols in DEM
    Cell_ID_Subs = np.c_[cellID[ind],ind[0],ind[1]]
    if Cell_ID_Subs.shape[0]!=N:
        raise TypeError('The size of DEM not consistent with '+fileTag+'.dat')
    Cell_ID_Subs = Cell_ID_Subs[Cell_ID_Subs[:,0].argsort()]
    subs = Cell_ID_Subs[:,[1,2]]
    subs = subs.astype('int64')
    if fileTag=='hU':
        zMat0=zMat+0
        zMat1=zMat+0
        zMat0[subs[:,0],subs[:,1]]=np.array(vecArray.iloc[:,0])
        zMat1[subs[:,0],subs[:,1]]=np.array(vecArray.iloc[:,1])
        zMat = [zMat0,zMat1]
    else:
        zMat[subs[:,0],subs[:,1]]=np.array(vecArray.Value)
    return zMat,zHead,zExtent