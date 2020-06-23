
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from clawpack.geoclaw.topotools import Topography
from numpy import *

#from pyclaw.data import Data
#probdata = Data('setprob.data')

grav = 9.81

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints=12893
    nypoints=7
    xupper=25782.e0
    yupper=5.e0
    xlower = 0.e0
    ylower = -5.e0
    outfile= "BMvsGC_80.topotype2"
    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")

def topo(x,y):
    """
    Run-up Channel
    """
    x2=25782.
    h0=-80.
    from numpy import where

    
    z=where(x<x2,h0,h0)
    

    for i in range(len(x)):
        for  j in range(len(x[i])):
            if (x[i,j]<x2):
                z[i,j]=h0
                #else:
                #    z[i,j]=(1.0/5.0)*(x[i,j]-149164.0)
    
        
    return z

if __name__=='__main__':
    maketopo()

