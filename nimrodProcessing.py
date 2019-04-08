# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 12:30:12 2019

@author: Xiaodong Ming
"""

#%% read downloaded NIMROD tar file
import os
import tarfile
import gzip
import nimrod
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from datetime import datetime
import matplotlib.pyplot as plt
import ArcGridDataProcessing as AP
#%% download tar data from internet
def downloadNimrodTar(year,dateStr='0131',localDir=[]):
    if len(localDir)==0:
        localDir = os.getcwd()+'/'
    from ftplib import FTP
    
    ftp = FTP('ftp.ceda.ac.uk')
    ftp.login('mxiaodong','mxd88719')
    remoteDir = 'badc/ukmo-nimrod/data/composite/uk-1km/'+str(year)+'/'
    ftp.cwd(remoteDir)
    #files = ftp.nlst()# Get All Files
    fileString = '*'+str(year)+dateStr+'*'
    files = ftp.nlst(fileString)
    # Print out the files
    for file in files:
        print("Downloading..." + file)
        ftp.retrbinary("RETR " + file ,open(localDir + file, 'wb').write)
    ftp.close()
#%% read one grid from NIMROD tar file
def extractOneGridFromNimrodTar(tarfileName,datetimeStr):
    """
    input:
        tarfileName: NIMROD tar file name
        datetimeStr: yyyyMMddHHmm('201012310830') the date and time string
    output:
        a nimrod object
    """
    gridObj = []
    tar = tarfile.open(tarfileName)
    tar.getmembers()
    for member in tar.getmembers():
        fgz = tar.extractfile(member)
        #print(member.name)
        if datetimeStr in member.name:
            print(member.name)
            f=gzip.open(fgz,'rb')
        # using nimrod package
            gridObj=nimrod.Nimrod(f)
            f.close()
    tar.close()   
    return gridObj

    
#%% extract data from NIMROD tar file
def getzMatFromNimrodTar(fileName,clipExtent=[]):
    #  clipExtent = (left,right,bottom,top)
    zMatList=[]
    dtStrList = []
    tar = tarfile.open(fileName)
    tar.getmembers()
    for member in tar.getmembers():
        fgz = tar.extractfile(member)
        #print(member.name)
        datetimeStr = member.name[31:31+12]
        f=gzip.open(fgz,'rb')
        # using nimrod package
        gridObj=nimrod.Nimrod(f)
        f.close()
        #del f,fgz,member
        if len(clipExtent)==4:
            gridObj = gridObj.apply_bbox(clipExtent[0],clipExtent[1],clipExtent[2],clipExtent[3])
        zMat,head,zExtent = gridObj.extract_asc()
        zMat[zMat==0]=np.nan #unit: mm/h
        zMatList.append(zMat)
        #del zMat
        dtStrList.append(datetimeStr)        
    tar.close()
    #del tar
    #import gc
    #gc.collect()
    zMatList = [x for _,x in sorted(zip(dtStrList,zMatList))]
    dtStrList.sort()
    return dtStrList,zMatList,head,zExtent
#%% create rainfall mask and rainfall source array for HiPIMS from 
def getRainMaskSource(tarList,demHead,datetimeStart,datetimeEnd=[]):
    
    """
    INPUTS
    tarList: a list of NIMROD tart files downloaded from badc
    demHead: a dictionary of dem file information
    datetimeStart: datetime object, the start date and time for rainfall source array
    datetimeEnd: datetime object, the end date and time for rainfall source array
    """
    #%% define clip extent
    R = demHead
    left = R['xllcorner']
    right = R['xllcorner']+R['ncols']*R['cellsize']
    bottom = R['yllcorner']
    top = R['yllcorner']+R['nrows']*R['cellsize']
    clipExtent = (left,right,bottom,top)
    #%% create rainfall mask grid according to the size of dem
    zMatList = []
    dtStrList = []
    for tarName in tarList:
        dtStr1,zMat1,head,zExtent = getzMatFromNimrodTar(tarName,clipExtent)
        zMatList = zMatList+zMat1
        dtStrList = dtStrList+dtStr1
#    zMatClip,headClip,extentClip = AP.ArraySquareClip(zMatList[0],head,zExtent)
     #mask value start from 0 and increase colum-wise
    maskValue = np.arange(np.size(zMat1)).reshape((head['nrows'],head['ncols']),order='F')
    rainMask= AP.MaskExtraction(maskValue,head,demHead)
    #%% create rainfall source array
    zMatList = [x for _,x in sorted(zip(dtStrList,zMatList))]
    dtStrList.sort()
    dates_list = [datetime.strptime(oneDate,'%Y%m%d%H%M')-datetimeStart for oneDate in dtStrList]
    timeSeries = [oneInterval.total_seconds() for oneInterval in dates_list]
    zMatSelected = [a.flatten(order='F') for a, b in zip(zMatList, timeSeries) if b>=0]
    timeSeries = np.array(timeSeries)
    timeSeries = timeSeries[timeSeries>=0]
    rainArray = np.array(zMatSelected)
    rainArray[np.isnan(rainArray)]=0
    rainSource = np.c_[timeSeries,rainArray/3600/1000]
    if not isinstance(datetimeEnd,list): 
        endTime = datetimeEnd-datetimeStart
        endTime = endTime.total_seconds()
        rainSource = rainSource[timeSeries<=endTime,:]
        
    return rainMask,rainSource
    
#%% ============================Visualization==================================
#%% plot the initial map, then renew raster files only
def initialMap(zMat,zExtent,poly_df=[],mapExtent=[],figsize=(6,8),vmin=0,vmax=10):        
    # create figure
    fig,ax = plt.subplots(1, figsize=figsize)
    # plot shapefile outline
    if len(poly_df)!=0:
        poly_df.plot(ax=ax,facecolor='none',edgecolor='r',linewidth=0.5, animated=True)
    else:
        mapExtent=zExtent 
    # create raster map    
    img = ax.imshow(zMat,extent=zExtent,vmin=vmin,vmax=vmax)    
    plt.axis('equal')
    if len(mapExtent)==0:
        mapExtent = (min(poly_df.bounds.minx),max(poly_df.bounds.maxx),
                     min(poly_df.bounds.miny),max(poly_df.bounds.maxy))
    ax.set_xlim(mapExtent[0], mapExtent[1])
    ax.set_ylim(mapExtent[2], mapExtent[3])
    
    # deal with x and y tick labels
    if mapExtent[1]-mapExtent[0]>10000:
        labels = [str(int(value/1000)) for value in ax.get_xticks()]
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(labels)
        labels = [str(int(value/1000)) for value in ax.get_yticks()]
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(labels,rotation=90)
        ax.set_xlabel('km towards east')
        ax.set_ylabel('km towards north')
    else:
        ax.set_xlabel('metre towards east')
        ax.set_ylabel('metre towards north')
        plt.yticks(rotation=90)
    
    ax.axes.grid(linestyle='--',linewidth=0.5)
    
    axins = inset_axes(ax,
               width="5%",  # width = 5% of parent_bbox width
               height="100%",  # height : 50%
               loc='lower right',
               bbox_to_anchor=(0.06, 0., 1, 1),
               bbox_transform=ax.transAxes,
               borderpad=0,
               )
    cbar=plt.colorbar(img,pad=0.05,cax=axins)
#        labels = cbar.ax.get_yticklabels()
#        cbar.ax.set_yticklabels(labels, rotation='vertical')
    cbar.ax.set_title('mm/h',loc='left')        
    return fig,ax,img

#%% plot all maps with data from zMat List
# need inputs from func getzMatFromNimrodTar and func initialMap
def plotRainRateFromList(zMatList,zExtent,dtStrList,
                         poly_df=[],mapExtent=[],
                         figsize=(6,8),vmin=0,vmax=10):
    """
    zMatList: List of rainfall rate matrix
    zExtent: extent of the rainfall matrix
    dtStrList: List of datetime strings
    poly_df: background polygon
    vmin,vmax: the range of rain value to be shown in figure
    """
    # plot the initial map, then renew raster files only
    fig,ax,img=initialMap(zMatList[0],zExtent,
                          poly_df,mapExtent,figsize,vmin,vmax)
    # renew raster file
    for i in range(len(zMatList)):
        titleStr = datetime.strptime(dtStrList[i], '%Y%m%d%H%M')
        titleStr=titleStr.strftime("%Y-%m-%d %H:%M")
        ax.set_title(titleStr)
        img.set_data(zMatList[i])
        fig.savefig(dtStrList[i]+'.jpg', dpi=100)
        print(titleStr)
    plt.close(fig)
#%% make an animation for rainfall rate    
def animateRainRate(videoName,zMatList=[],zExtent=[],dtStrList=[],
                         poly_df=[],mapExtent=[],imagesList=[],
                         figsize=(6,8),vmin=0,vmax=10,fps=30):
    # plot all the rainfall files
    import cv2

    if len(imagesList)==0:
        plotRainRateFromList(zMatList,zExtent,dtStrList,
                         poly_df,mapExtent,
                         figsize,vmin,vmax)
        images = [s + '.jpg' for s in dtStrList]
        images.sort()
    else:
        images = imagesList
    frame = cv2.imread(images[0])
    height, width, layers = frame.shape
    video = cv2.VideoWriter(videoName, -1, fps, (width,height))
    for image in images:
        video.write(cv2.imread(image))
        os.remove(image)
    cv2.destroyAllWindows()
    video.release()
#%%
def animateRainRateFromNimrodTar(videoName,tarList,clipExtent=[],poly_df=[],mapExtent=[],
                         figsize=(6,8),vmin=0,vmax=10,fps=30):
    imagesList = []
    for tarName in tarList:
        dtStrList,zMatList,head,zExtent = getzMatFromNimrodTar(tarName,clipExtent)
        plotRainRateFromList(zMatList,zExtent,dtStrList,
                         poly_df,mapExtent,
                         figsize,vmin,vmax)
        imagesList = imagesList+dtStrList
    imagesList = [s + '.jpg' for s in imagesList]
    imagesList.sort()
    animateRainRate(videoName,imagesList=imagesList,fps=fps)
    