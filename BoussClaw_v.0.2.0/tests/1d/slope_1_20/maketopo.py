
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 1001
    nypoints = 21
    xlower = -50.e0
    xupper =  50.e0
    yupper =  1.e0
    ylower = -1.e0
    outfile= "tank.tt2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 101
    nypoints = 101
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    x1= 19.85
    y1=-1.
    from numpy import where
    #z = where( x<19.85, -1. , -1.)
    z = where( x<19.85, -x/19.85 , -1.)	
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    ze = -((x+0e0)**2 + (y+0e0)**2)/10.
    z = where(ze>-10., 40.e0*exp(ze), 0.)
    return z

if __name__=='__main__':
    maketopo()
    #makeqinit()
