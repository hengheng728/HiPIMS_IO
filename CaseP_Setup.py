#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To decompose domain into sub-sections for multi-GPU computing
Created on Sun Dec  2 21:40:54 2018

@author: Xiaodong MIng
"""
import os
import numpy as np
import time
#os.chdir('/Users/b4042552/Dropbox/Python')
import HiPIMS as HI
import HiPIMS.HiPIMS_IO as HPIO
start=time.perf_counter()
demMat,demHead,demExtent = HI.ArcGridDataProcessing.arcgridread('/Users/b4042552/Dropbox/Python/Data/UpperLee500m.asc')
end=time.perf_counter()
print(end-start)
#%%
bound1Points = np.array([[535, 206], [545, 206], [545, 210], [535, 210]])*1000
bound2Points = np.array([[520, 230], [530, 230], [530, 235], [520, 235]])*1000
dBound0 = {'polyPoints': [],'type': 'open','h': [],'hU': []}
dBound1 = {'polyPoints': bound1Points,'type': 'open','h': [[0,10],[60,10]]}
dBound2 = {'polyPoints': bound2Points,'type': 'open','hU': [[0,50000],[60,30000]]}
boundList = [dBound0,dBound1,dBound2]
del dBound0,dBound1,dBound2,bound1Points,bound2Points
rain_source = np.array([[0,100/1000/3600/24],
                        [86400,100/1000/3600/24],
                        [86401,0]])
gauges_pos = np.array([[534.5,231.3],
                       [510.2,224.5],
                       [542.5,225.0],
                       [538.2,212.5],
                       [530.3,219.4]])*1000
rootPath='/Users/b4042552/Dropbox/Python/CaseP'
numSection=4
#%% Single GPU
start = time.perf_counter()
HPIO.InputSetup(rootPath,demMat,demHead,h0=0,
                        boundList=boundList,fileToCreate='all',
                        rain_source = rain_source,
                        gauges_pos=gauges_pos)
end = time.perf_counter()
print('total time elapse: '+str(end-start))
#%% Multiple GPU
start = time.perf_counter()
HPIO.InputSetup_MG(rootPath,demMat,demHead,numSection=numSection,h0=0,
                        boundList=boundList,fileToCreate='all',
                        rain_source = rain_source,
                        gauges_pos=gauges_pos)
end = time.perf_counter()
print('total time elapse: '+str(end-start))
#%%
from InputSetupFuncs_MG import Ztype2Grid_MG
A,B,C=Ztype2Grid_MG(rootPath,numSection,'precipitation_mask')