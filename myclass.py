#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 13:03:25 2019

@author: b4042552
"""
import os
import numpy as np
import gzip
import matplotlib.pyplot as plt
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable

class Test(object):
    def __init__(self,name,score):
        self.name = name
        self.score = score
        self.__add_attributes()
        
    def __add_attributes(self):
        self.name_upper = self.name.upper()
        self.score_rate = self.score/100
#%% To deal with raster data
class raster(object):
    """
    Created on Tue Apr 7 2019
    @author: Xiaodong Ming 
    
    To deal with raster data with a ESRI ASCII or GTiff format  
    
    Properties:
        sourceFile: file name to read grid data
        outputFile: file name to write in grid data
        array: a numpy array storing grid cell values
        header: a dict storing reference information of the grid
        extent: a tuple storing outline limits of the raster (left,right,bottom,top)
        extent_dict: a dictionary storing outline limits of the raster
        projection: (string) the Well-Known_Text (wkt) projection information
    
    Methods(public):
        Write_asc: write grid data into an asc file compressed(.gz) or decompressed 
        To_osgeo_raster: convert this object to an osgeo raster object
        RectClip: clip raster according to a rectangle extent
        Clip: clip raster according to a polygon
        Rasterize: rasterize the shapefile to the raster object and return a bool array
            with Ture value in and on the polygon/polyline
        Mapshow:    
            
    
    Methods(private):
        __header2extent: convert header to extent
        __read_asc
        __map2sub
        __sub2map
        __read_asc
        __read_tif
        
        
        
    """
#%%======================== initialization function ===========================   
    def __init__(self,sourceFile=None,array=None,header=None,epsg=None):
        """
        sourceFile: name of a asc/tif file if a file read is needed
        
        """
        self.sourceFile = sourceFile
        if epsg is not None:
            self.projection = self.__SetWktProjection(epsg)
        else:
            self.projection = None
        if sourceFile is None:
            self.array = array
            self.header = header
            self.sourceFile = 'sourceFile.asc'
        else:
            if sourceFile.endswith('.tif'):
                self.__read_tif() # only read the first band
            else:
                self.__read_asc()

        if isinstance(self.header,dict)==0:
            raise ValueError('header is not a dictionary')
        else:
            # create self.extent and self.extent_dict 
            self.__header2extent()
            
#%%============================= Spatial analyst ==============================
   
    def RectClip(self,clipExtent):
        """
        clipExtent: left,right,bottom,top
        clip raster according to a rectangle extent
        return:
           a new raster object
        """
        new_obj = copy.deepcopy(self)
        X = clipExtent[0:2]
        Y = clipExtent[2:4]
        rows,cols = self.__map2sub(X,Y)
        Xcentre,Ycentre = self.__sub2map(rows,cols)
        xllcorner = min(Xcentre)-0.5*self.header['cellsize']
        yllcorner = min(Ycentre)-0.5*self.header['cellsize']
        # new array
        new_obj.array = self.array[min(rows):max(rows),min(cols):max(cols)]
        # new header
        new_obj.header['nrows'] = new_obj.array.shape[0]
        new_obj.header['ncols'] = new_obj.array.shape[1]
        new_obj.header['xllcorner'] = xllcorner
        new_obj.header['yllcorner'] = yllcorner
        # new extent
        new_obj.__header2extent()
        new_obj.sourceFile = None       
        return new_obj
    
    def Clip(self,mask=None):
        """
        clip raster according to a mask
        mask: 
            1. string name of a shapefile
            2. numpy vector giving X and Y coords of the mask points
        
        return:
            a new raster object
        """
        from osgeo import ogr
        if isinstance(mask, str):
            shpName =  mask
        # Open shapefile datasets        
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDataset = shpDriver.Open(shpName, 0) # 0=Read-only, 1=Read-Write
        layer = shpDataset.GetLayer()
        shpExtent = np.array(layer.GetExtent()) #(minX,maxY,maxX,minY)           
        # 1. rectangle clip raster
        new_obj = self.RectClip(shpExtent)
        new_raster = copy.deepcopy(new_obj)                
        indexArray = new_raster.Rasterize(shpDataset)
        arrayClip = new_raster.array
        arrayClip[indexArray==0]=new_raster.header['NODATA_value']
        new_raster.array = arrayClip        
        shpDataset.Destroy()
        
        return new_raster
    
    def Rasterize(self,shpDSName,rasterDS=None):
        """
        rasterize the shapefile to the raster object and return a bool array
            with Ture value in and on the polygon/polyline
        shpDSName: string for shapefilename, dataset for ogr('ESRI Shapefile') object
        
        return numpy array
        """
        from osgeo import gdal,ogr
        if isinstance(shpDSName,str):
            shpDataset = ogr.Open(shpDSName)
        else:
            shpDataset = shpDSName
        layer = shpDataset.GetLayer()
        if rasterDS is None:
            target_ds = self.To_osgeo_raster()
        else:
            target_ds = rasterDS
        gdal.RasterizeLayer(target_ds, [1], layer,burn_values=[-3333])
        indexArray = target_ds.ReadAsArray()
        indexArray[indexArray!=-3333]=0
        indexArray[indexArray==-3333]=1
        target_ds=None
        return indexArray
         
#%%=============================Visualization==================================
    #%% draw inundation map with domain outline
    def Mapshow(self,figureName=None,figsize=None,dpi=300,vmin=None,vmax=None,
                cax=True):
        """
        Display raster data without projection
        figureName: the file name to export map,if figureName is empty, then
            the figure will not be saved
        figsize: the size of map
        dpi: The resolution in dots per inch
        vmin and vmax define the data range that the colormap covers
        """
        np.warnings.filterwarnings('ignore')    
        fig, ax = plt.subplots(1, figsize=figsize)
        # draw inundation
        zMat = self.array
        zMat[zMat==self.header['NODATA_value']]=np.nan
        img=plt.imshow(zMat,extent=self.extent,vmin=vmin,vmax=vmax)
        # colorbar
    	# create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        if cax==True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(img, cax=cax)
        ax.axes.grid(linestyle='-.',linewidth=0.2)
        # save figure
        if figureName is not None:
            fig.savefig(figureName, dpi=dpi)
            
        return fig,ax
#%%===============================output=======================================   
    def Write_asc(self,outputFile,EPSG=None,compression=False):
        
        """
        write raster as asc format file 
        outputFile: output file name
        EPSG: epsg code, if it is given, a .prj file will be written
        compression: logic, whether compress write the asc file as gz
        """
        if compression:
            if not outputFile.endswith('.gz'):
                outputFile=outputFile+'.gz'        
        self.outputFile = outputFile
        Z = self.array
        header = self.header
        Z[np.isnan(Z)]= header['NODATA_value']
        if not isinstance(header,dict):
            raise TypeError('bad argument: header')
                     
        if outputFile.endswith('.gz'):
            f = gzip.open(outputFile, 'wb') # write compressed file
        else:
            f = open(outputFile, 'wb')
        f.write(b"ncols    %d\n" % header['ncols'])
        f.write(b"nrows    %d\n" % header['nrows'])
        f.write(b"xllcorner    %g\n" % header['xllcorner'])
        f.write(b"yllcorner    %g\n" % header['yllcorner'])
        f.write(b"cellsize    %g\n" % header['cellsize'])
        f.write(b"NODATA_value    %g\n" % header['NODATA_value'])
        np.savetxt(f,Z,fmt='%g', delimiter=' ')
        f.close()
        if EPSG is not None:
            self.__SetWktProjection(EPSG)
        # if projection is defined, write .prj file    
        if self.projection is not None:    
            prjFileName = outputFile
            wkt = self.projection
            if outputFile.endswith('.asc'):
                prjFileName=prjFileName[0:-4]+'.prj'
            elif outputFile.endswith('.asc.gz'):
                prjFileName=prjFileName[0:-7]+'.prj'
                
            prj = open(prjFileName, "w")            
            prj.write(wkt)
            prj.close()
            
        return None
    
    # convert this object to an osgeo raster object
    def To_osgeo_raster(self,filename=None,fileformat = 'GTiff',destEPSG=27700):        
        """
        convert this object to an osgeo raster object, write a tif file if 
            necessary
        filename: the output file name, if it is given, a tif file will be produced
        fileformat: GTiff or AAIGrid
        destEPSG: the EPSG projection code default: British National Grid
        
        return:
            an osgeo raster dataset
            or a tif filename if it is written
        """
        from osgeo import gdal,osr
        if filename is None:
            dst_filename = 'temp.tif'
        else:
            dst_filename = filename
        if not dst_filename.endswith('.tif'):
            dst_filename = dst_filename+'.tif'
    
        # You need to get those values like you did.
        PIXEL_SIZE = self.header['cellsize']  # size of the pixel...        
        x_min = self.extent[0] # left  
        y_max = self.extent[3] # top
        dest_crs = osr.SpatialReference()
        dest_crs.ImportFromEPSG(destEPSG)
        driver = gdal.GetDriverByName(fileformat)
    
        dataset = driver.Create(dst_filename,
            xsize=self.header['ncols'],
            ysize=self.header['nrows'],
            bands=1,
            eType=gdal.GDT_Float32)
    
        dataset.SetGeoTransform((
            x_min,    # 0
            PIXEL_SIZE,  # 1
            0,                      # 2
            y_max,    # 3
            0,                      # 4
            -PIXEL_SIZE))  
    
        dataset.SetProjection(dest_crs.ExportToWkt())
        array = self.array
        array[array==self.header['NODATA_value']]=np.nan
        dataset.GetRasterBand(1).WriteArray(self.array)
        if filename is not None:
            dataset.FlushCache()  # Write to disk.
            dataset = None
            return dst_filename
        else:
            os.remove(dst_filename)
            return dataset

#%%=========================== private functions ==============================        
    def __header2extent(self):
        """
        To convert header (dict) to a spatial extent of the DEM
        extent: (left,right,bottom,top)
        """
        R = self.header
        left = R['xllcorner']
        right = R['xllcorner']+R['ncols']*R['cellsize']
        bottom = R['yllcorner']
        top = R['yllcorner']+R['nrows']*R['cellsize']
        self.extent = (left,right,bottom,top)
        self.extent_dict = {'left':left, 'right':right, 'bottom':bottom, 'top':top}

    def __map2sub(self,X,Y):
        """
        convert map points to subscripts of a matrix with geo reference header
        X,Y: coordinates in map units
        return
            rows,cols: (int) subscripts of the data matrix
        """
        #x and y coordinate of the centre of the first cell in the matrix
        X = np.array(X)
        Y = np.array(Y)
        header = self.header
        
        x0 = header['xllcorner']+0.5*header['cellsize']
        y0 = header['yllcorner']+(header['nrows']-0.5)*header['cellsize']
        rows = (y0-Y)/header['cellsize'] # row and col number starts from 0
        cols = (X-x0)/header['cellsize']
        if isinstance(rows,np.ndarray):
            rows = rows.astype('int64')
            cols = cols.astype('int64') #.astype('int64')
        else:
            rows = int(rows)
            cols = int(cols)
        return rows,cols

    def __sub2map(self,rows,cols):
        """
        convert subscripts of a matrix to map coordinates 
        rows,cols: subscripts of the data matrix, starting from 0
        return
            X,Y: coordinates in map units
        """
        #x and y coordinate of the centre of the first cell in the matrix
        if not isinstance(rows,np.ndarray):
            rows = np.array(rows)
            cols = np.array(cols)        
        
        header = self.header
        left = self.extent[0] #(left,right,bottom,top)
        top = self.extent[3]
        X = left + (cols+0.5)*header['cellsize']
        Y = top  - (rows+0.5)*header['cellsize']
         
        return X,Y

# read ascii file        
    def __read_asc(self):
        """
        read asc file and return array,header
        if self.sourceFile ends with '.gz', then read the compressed file
        """
        fileName = self.sourceFile
        try:
            fh = open(fileName, 'r')
            fh.close()
        # Store configuration file values
        except FileNotFoundError:
            # Keep preset values
            print('Error: '+fileName+' does not appear to exist')
            return
        # read header
        header = {} # store header information including ncols, nrows,...
        numheaderrows = 6
        n=1
        if fileName.endswith('.gz'):
            # read header
            with gzip.open(fileName, 'rt') as f:                
                for line in f:
                    if n<=numheaderrows:
                        line = line.split(" ",1)
                        header[line[0]] = float(line[1])
                    else:
                        break
                    n = n+1
        else:
            # read header
            with open(fileName, 'rt') as f:            
                for line in f:
                    if n<=numheaderrows:
                        line = line.split(" ",1)
                        header[line[0]] = float(line[1])
                    else:
                        break
                    n = n+1
    # read value array
        array  = np.loadtxt(fileName, skiprows=numheaderrows,dtype='float64')
        array[array == header['NODATA_value']] = float('nan')
        header['ncols']=int(header['ncols'])
        header['nrows']=int(header['nrows'])
        self.array = array
        self.header = header
        prjFile = self.sourceFile[:-4]+'.prj'
        if os.path.isfile(prjFile):
            with open(prjFile, 'r') as file:
                projection = file.read()
            self.projection = projection
        return None

# read GTiff file
    def __read_tif(self):
        """
        read tif file and return array,header
        only read the first band
        """
        from osgeo import gdal
        tifName = self.sourceFile
        ds = gdal.Open(tifName)
        
        ncols = ds.RasterXSize
        nrows = ds.RasterYSize
        geoTransform = ds.GetGeoTransform()
        x_min = geoTransform[0]
        cellsize = geoTransform[1]
        y_max = geoTransform[3]
        xllcorner = x_min
        yllcorner = y_max - nrows*cellsize
        rasterBand = ds.GetRasterBand(1)
        NODATA_value = rasterBand.GetNoDataValue()
        array = rasterBand.ReadAsArray()
        header = {'ncols':ncols, 'nrows':nrows, 
                  'xllcorner':xllcorner, 'yllcorner':yllcorner,
                  'cellsize':cellsize,'NODATA_value':NODATA_value}        
        self.header = header
        self.array = array
        self.projection = ds.GetProjection()
        rasterBand = None
        ds = None
        return None
    def __SetWktProjection(self,epsg_code):
        """
        get coordinate reference system (crs) as Well Known Text (WKT) 
            from https://epsg.io
        epsg_code: the epsg code of a crs, e.g. BNG:27700, WGS84:4326
        return wkt text
        """
        import requests
        # access projection information
        wkt = requests.get('https://epsg.io/{0}.prettywkt/'.format(epsg_code))
        # remove spaces between charachters
        remove_spaces = wkt.text.replace(" ","")
        # place all the text on one line
        output = remove_spaces.replace("\n", "")
        self.projection = output
        return output         

        
