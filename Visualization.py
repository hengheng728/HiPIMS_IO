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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gzip
from osgeo import osr, gdal 
import shutil
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
	# create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
	
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
    
    return bnMat
#%% draw inundation map with domain outline
def InundationMap(zMat,zHead,depth=0.2,figureName=[],
                    figsize=(10, 6),dpi=300,vmin=0.2,vmax=False):
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
    if vmax!=False:
        img=plt.imshow(zMat,extent=zExtent,vmin=vmin,vmax=vmax)
    else:
        img=plt.imshow(zMat,extent=zExtent)
    # colorbar
	# create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    # draw domain outline cells
    ax.scatter(boundX,boundY,s=0.2,facecolor='r')
    ax.axes.grid(linestyle='-.',linewidth=0.2)
    # save figure
    if len(figureName)>0:
        fig.savefig(figureName, dpi=dpi)
        
    return fig,ax

#%% convert ascii file to geotif
def Asc2GeoTiff(ascfileName,tiffFileName,srcEPSG=27700,destEPSG=4326):
    # convert a projected asc file to a tiff with geographical coordinates
    # srcEPSG = 27700 #BNG
    # destEPSG = 4326 #WGS84
    source = ascfileName
    target = tiffFileName
    if ~target.endswith('.tif'):
        target=target+'.tif'
    # open CSV source file
    if source.endswith('.gz'):
        with gzip.open(source, 'rb') as f_in:
            with open(source[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        source = source[:-3]

    cvs = gdal.Open(source)
    #cvs = gdal.Open(source)
    if cvs is None:
        print( 'ERROR: Unable to open %s' % source)
        return

    # get GeoTIFF driver
    geotiff = gdal.GetDriverByName("GTiff")
    if geotiff is None:
        print( 'ERROR: GeoTIFF driver not available.')
        return None

    # set source coordinate system of coordinates in CSV file
    src_crs = osr.SpatialReference()
    src_crs.ImportFromEPSG(srcEPSG)

    # set destination projection parameters
    dest_crs = osr.SpatialReference()
    dest_crs.ImportFromEPSG(destEPSG)

    # set coordinate transformation
    tx = osr.CoordinateTransformation(src_crs, dest_crs)

    # get raster dimension related parameters of source dataset
    xo, xs, xr, yo, yr, ys = cvs.GetGeoTransform()
    xsize = cvs.RasterXSize
    ysize = cvs.RasterYSize

    # convert corner coordinates from old to new coordinate system
    (ulx, uly, ulz) = tx.TransformPoint(xo, yo)
    (lrx, lry, lrz) = tx.TransformPoint(xo + xs * xsize + xr * ysize,\
                                        yo + yr * xsize + ys * ysize)

    # create blank in-memory raster file with same dimension as CSV raster
    mem = gdal.GetDriverByName('MEM')
    dest_ds = mem.Create('', xsize, ysize, 1, gdal.GDT_Float32)

    # get new transformation
    dest_geo = (ulx, (lrx - ulx) / xsize, xr,\
                uly, yr, (lry - uly) / ysize)

    # set the geotransformation
    dest_ds.SetGeoTransform(dest_geo)
    dest_ds.SetProjection(dest_crs.ExportToWkt())

    # project the source raster to destination coordinate system
    gdal.ReprojectImage(cvs, dest_ds, \
                        src_crs.ExportToWkt(), dest_crs.ExportToWkt(),\
                        gdal.GRA_Bilinear, 0.0, 0.0)

    # save projected in-memory raster to disk
    geotiff.CreateCopy(target, dest_ds, 0 )
    return None

