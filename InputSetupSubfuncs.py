#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:24:53 2018

@author: Xiaodong MIng
"""
import numpy as np
#import matplotlib as mpl
#import time

#%% write Boundary source
def writeBoundSource(fileFolder,h_BC_source,hU_BC_source):
    m_h = 0
    m_hU = 0
    fileName_h = []
    fileName_hU = []
    for sourceValue in h_BC_source:            
        fileName = fileFolder+'h_BC_'+str(m_h)+'.dat'
        np.savetxt(fileName,sourceValue,fmt=['%16.2f','%16.4f'],delimiter=' ')
        fileName_h.append(fileName)
        print(fileName+' is created')
        m_h=m_h+1
    for sourceValue in hU_BC_source:            
        fileName = fileFolder+'hU_BC_'+str(m_hU)+'.dat'
        np.savetxt(fileName,sourceValue,fmt=['%16.2f','%16.6f','%16.6f'],delimiter=' ')
        fileName_hU.append(fileName)
        print(fileName+' is created')
        m_hU=m_hU+1
    return fileName_h,fileName_hU
  
#%% write Gauge Position
def writeGaugePos(fileFolder,gauges_pos):
    if len(gauges_pos)==0:
        gauges_pos = np.array([0,0],ndmin=2)
    fileName = fileFolder+'gauges_pos.dat'
    with open(fileName,'w') as file2write:
        np.savetxt(file2write,gauges_pos,fmt=['%16.4f','%16.4f'],delimiter=' ')
    print(fileName+' is created')
