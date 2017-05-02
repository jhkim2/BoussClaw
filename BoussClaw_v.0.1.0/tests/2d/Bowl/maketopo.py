
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 200+1
    nypoints = 200+1
    xlower = -1000.e0
    xupper =  1000.e0
    ylower = -1000.e0
    yupper =  1000.e0
    outfile= "bowl.tt3"     
    topotools.topo3writer(outfile,topo2,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    s1 =-1./150. # trench angle
    s2 =-3.*pi/180. # continental slope angle
    s3 =-.2*pi/180. # continental shelf angle
    s4 = -1./10.     # slope at beachmek
    y1 = -7000.
    x1 = 0.
    y2 = -2000.
    x2 = (y2-y1)/s1
    y3 = -150.
    x3 = x2 + (y3-y2)/s2
    y4 = -20.
    x4 = x3 + (y4-y3)/s3
    from numpy import where
    z = where( x>x1 , y1 , where(x>x2, s1*x+y1, \
        where(x>x3, s2*(x-x2)+y2, where(x>x4, s3*(x-x3)+y3, s4*(x-x4)+y4 )))) 
    print x1,x2,x3,x4
    return z

def topo1(x,y):
    s1 =-5./75. # trench angle
    s2 =-2./150. # continental shelf
    y1 = -7000.
    x1 = 0.
    y2 = -2000.
    x2 = (y2-y1)/s1
    from numpy import where
    z = where( x>x1 , y1 , where(x>x2, s1*x+y1, s2*(x-x2)+y2)) *0.-1.
    return z

def topo2(x,y):
    from numpy import where,sqrt
    r = sqrt(x**2+y**2)
    z = -4.*(1.-(r/80.)**2)
    return z

if __name__=='__main__':
    maketopo()

