#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 21:17:28 2019

@author: b4042552
"""
#%%
import sys
sys.path.insert(0,'/Users/b4042552/Dropbox/Python/HiPIMS')
from myclass import raster

objDEM = raster('Data/UpperLee50m.asc')
#Shapefile = '/Data/ExamplePolygon.shp'
#raster_ds = objDEM.to_osgeo_raster()
#objDEM_clip = objDEM.Clip(Shapefile)
#objDEM_clip.Write_asc('example_clip.asc',EPSG=27700)
#%%
#Shapefile = 'RoadsLine.shp'
#data = gdal.Open('DEM2m_NewcastleFrame2.tif')