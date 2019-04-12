#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 13:03:25 2019

@author: b4042552
"""
import numpy as np
import gzip
import matplotlib.pyplot as plt
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
    
    To deal with raster data
    
    Properties:
        sourceFile: file name to read grid data
        outputFile: file name to write in grid data
        dataArray: numpy array storing grid cell value
        dataHead: dict storing reference information of the grid
        extent: 
    
    Methods(public):
        write_asc: write grid data into an asc file compressed(.gz) or decompressed 
    
    Methods(private):
        __head2extent: convert head to extent
        __read_asc
        
    """
    def __init__(self,sourceFile=None,z=None,head=None):
        """
        sourceFile: name of a asc file if it is needed
        
        """
        self.sourceFile = sourceFile
        if sourceFile is not None:
            z,head = self.__read_asc()
        self.dataArray = z
        self.dataHead = head
        if isinstance(self.dataHead,dict)==0:
            raise ValueError('dataHead is not a dictionary')
        else:
            self.__head2extent()
        
    def __head2extent(self):
        """
        To convert head (dict) to a spatial extent of the DEM
        extent: (left,right,bottom,top)
        """
        R = self.dataHead
        left = R['xllcorner']
        right = R['xllcorner']+R['ncols']*R['cellsize']
        bottom = R['yllcorner']
        top = R['yllcorner']+R['nrows']*R['cellsize']
        self.extent = (left,right,bottom,top)
        
    def __read_asc(self):
        """
        read asc file and return dataArray,head,extent
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
        # read head
        head = {} # store head information including ncols, nrows,...
        numheadrows = 6
        n=1
        if fileName.endswith('.gz'):
            # read head
            with gzip.open(fileName, 'rt') as f:                
                for line in f:
                    if n<=numheadrows:
                        line = line.split(" ",1)
                        head[line[0]] = float(line[1])
                    else:
                        break
                    n = n+1
        else:
            # read head
            with open(fileName, 'rt') as f:            
                for line in f:
                    if n<=numheadrows:
                        line = line.split(" ",1)
                        head[line[0]] = float(line[1])
                    else:
                        break
                    n = n+1
    # read value array
        dataArray  = np.loadtxt(fileName, skiprows=numheadrows,dtype='float64')
        #dataArray[dataArray == head['NODATA_value']] = float('nan')
        head['ncols']=int(head['ncols'])
        head['nrows']=int(head['nrows'])
        return dataArray,head
#%%===============================output=======================================   
    def write_asc(self,outputFile,compression=False):
        """
        write raster as asc format file
        outputFile: output file name
        compression: logic, whether compress write the file as gz
        """
        if compression:
            if not outputFile.endswith('.gz'):
                outputFile=outputFile+'.gz'        
        self.outputFile = outputFile
        Z = self.dataArray
        head = self.dataHead
        Z[np.isnan(Z)]= head['NODATA_value']
        if not isinstance(head,dict):
            raise TypeError('bad argument: head')
                     
        if outputFile.endswith('.gz'):
            f = gzip.open(outputFile, 'wb') # write compressed file
        else:
            f = open(outputFile, 'wb')
        f.write(b"ncols    %d\n" % head['ncols'])
        f.write(b"nrows    %d\n" % head['nrows'])
        f.write(b"xllcorner    %g\n" % head['xllcorner'])
        f.write(b"yllcorner    %g\n" % head['yllcorner'])
        f.write(b"cellsize    %g\n" % head['cellsize'])
        f.write(b"NODATA_value    %g\n" % head['NODATA_value'])
        np.savetxt(f,Z,fmt='%g', delimiter=' ')
        f.close()        
#%%=============================Visualization==================================
    #%% draw inundation map with domain outline
    def mapshow(self,figureName=None,figsize=None,dpi=300,vmin=None,vmax=None,
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
        zMat = self.dataArray
        zMat[zMat==self.dataHead['NODATA_value']]=np.nan
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

        
