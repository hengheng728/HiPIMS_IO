#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:36:25 2018

@author: Xiaodong Ming
"""
import os
import OutputManage as OM
import Visualization as VL
rootPath = '/home/b1055010/GeoClasses/release/bin/Eden6_5mMGPUs/'
os.chdir(rootPath)
numSection=6
fileTag = 'h_21600'
grid,head,extent=OM.CombineGridFile(rootPath,numSection,fileTag)
fileName = fileTag+'.asc'
print(fileName+' is combined')
OM.ArcgridwriteGZip(fileName,grid,head)
figureName=fileTag+'.png'
VL.InundationMap(grid,head,depth=0.1,figureName=figureName,
                    figsize=(10, 6),dpi=800)

print(figureName+'.png is drawn')

#Time=21600
#OM.Backup2Initial(rootPath,numSection,Time,remove=True)
