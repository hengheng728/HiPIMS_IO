#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 21:17:28 2019

@author: b4042552
"""
#%%
from myclass import raster

objDEM = raster('ExampleDEM.asc')
objDEM.mapshow()