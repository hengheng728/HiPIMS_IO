#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 14:43:41 2018

@author: Xiaodong Ming
"""
#%%
from OutputManage import CombineGridFile,ArcgridwriteGZip


rootPath = '/Users/b4042552/Dropbox/Python/CaseP'
fileTag = 'h_7200'
numSection = 4
zMat,zHead=CombineGridFile(rootPath,numSection,fileTag)

fileName = rootPath+'/'+fileTag+'.asc'
ArcgridwriteGZip(fileName,zMat,zHead)