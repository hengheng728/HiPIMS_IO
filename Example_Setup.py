#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To demonstrate how to generate input data for HiPIMS for both single and multiple GPUs
Created on Sun Dec  2 21:40:54 2018

@author: Xiaodong MIng
"""
import os
import sys
import numpy as np
import time
scriptsPath = '/Users/b4042552/Dropbox/Python/HiPIMS' # position storing HiPIMS_IO.py and ArcGridDataProcessing.py
sys.path.insert(0,scriptsPath)

from HiPIMS_IO import HiPIMS_setup
import ArcGridDataProcessing as AP

# define root path for the example case
rootPath='/Users/b4042552/Dropbox/Python/CaseP'

# read example DEM data
demMat,demHead,demExtent = AP.arcgridread('ExampleDEM.asc') # stored in the same dir with HiPIMS_IO.py

# define boundary condition
bound1Points = np.array([[535, 206], [545, 206], [545, 210], [535, 210]])*1000
bound2Points = np.array([[520, 230], [530, 230], [530, 235], [520, 235]])*1000
dBound0 = {'polyPoints': [],'type': 'open','h': [],'hU': []}
dBound1 = {'polyPoints': bound1Points,'type': 'open','h': [[0,10],[60,10]]}
dBound2 = {'polyPoints': bound2Points,'type': 'open','hU': [[0,50000],[60,30000]]}
boundList = [dBound0,dBound1,dBound2]
del dBound0,dBound1,dBound2,bound1Points,bound2Points

# define rainfall source, a same rainfall source for the whole model domain
# the rainfall mask is default defined as 0 for all the domain cells
rain_source = np.array([[0,100/1000/3600/24],
                        [86400,100/1000/3600/24],
                        [86401,0]])
# define monitor positions
gauges_pos = np.array([[534.5,231.3],
                       [510.2,224.5],
                       [542.5,225.0],
                       [538.2,212.5],
                       [530.3,219.4]])*1000

# define the number of GPUs to use in HiPIMS
numSection=1
os.chdir(rootPath)

# generate input files for HiPIMS
start = time.perf_counter()
HiPIMS_setup(rootPath,demMat,demHead,numSection=numSection,h0=0,
                        boundList=boundList,fileToCreate='all',
                        rain_source = rain_source,
                        gauges_pos=gauges_pos)
end = time.perf_counter()
print('total time elapse: '+str(end-start))
