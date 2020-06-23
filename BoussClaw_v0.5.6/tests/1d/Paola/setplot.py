"""
Set up the plot figures, axes, and items to be done for each frame.
This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
from six.moves import range

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'ascii'    # 'ascii' or 'binary' to match setrun.py


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Slice', figno=1)
    #plotfigure.kwargs = {'figsize': (9,3)}
    #plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        eta = q[3,:,m]
        return x,eta

    def y_slice(current_data):
        y = current_data.y[0,:]
        q = current_data.q
        eta = q[3,1,:]
        return y,eta

    def B_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        h = q[0,:,m]
        eta = q[3,:,m]
        B = eta - h
        return x,B

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('water surface')
    plotaxes.title = 'eta'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.5,4]

    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    # Topography
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = B_slice
    plotitem.color = 'k'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':2}

    def fixup(current_data):
        import pylab
        import numpy as np
        t = current_data.t
        pylab.title('Surface at t = %4.2f' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []  # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.format = 'ascii'                # Format of output
    # plotdata.format = 'netcdf'

    return plotdata
