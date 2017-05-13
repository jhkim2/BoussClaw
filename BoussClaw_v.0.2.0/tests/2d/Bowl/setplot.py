""" 
Set up the plot figures, axes, and items to be done for each frame.
This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters. 
""" 

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
    
import os

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)


    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    #plotfigure.kwargs = {'figsize': (7,14)}
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    drytol = 1e-5
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def surface_or_depth(current_data):
        from numpy import ma, where
        q = current_data.q
        h = q[0,:,:]
        eta = q[3,:,:]
        topo = eta - h
        surface = ma.masked_where(h<=drytol, eta)
        depth = ma.masked_where(h<=drytol, h)
	surface_or_depth =  where(topo<0, surface, depth)
	return surface_or_depth

    def fixup(current_data):
        import pylab
        t = current_data.t
        addgauges(current_data)
        pylab.title('Surface at %4.2f seconds' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
        #pylab.xticks(rotation=20)
    plotaxes.afteraxes = fixup
	
    #plotitem.plot_var = surface_or_depth
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin =  -.2
    plotitem.pcolor_cmax =   .2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='x-slice', figno=1)
    plotfigure.kwargs = {'figsize': (12,5)}
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
    plotaxes.ylimits = [-.2,.5]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    # Water Surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.outdir = '../bous_cfl9_dx05'
    #plotitem.map_2d_to_1d = x_slice
    #plotitem.color = 'r'
    #plotitem.plotstyle = '--'
    #plotitem.kwargs = {'linewidth':1}
    
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
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='zoom', figno=2)
    #plotfigure.kwargs = {'figsize': (10,20)}
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    drytol = 1e-5

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        pylab.title('Surface at %4.2f seconds' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def surface_or_depth(current_data):
        from numpy import ma, where
        q = current_data.q
        h = q[0,:,:]
        eta = q[3,:,:]
        topo = eta - h
        surface = ma.masked_where(h<=drytol, eta)
        depth = ma.masked_where(h<=drytol, h)
	surface_or_depth =  where(topo<0, surface, depth)
	return surface_or_depth

    def aa(current_data):
        from numpy import where,sqrt
        q = current_data.q
        h = q[0,:,:]
        u = where(h>0.,q[1,:,:]/h,0.)
        t = current_data.t
        t = (t-2.)*sqrt(9.81)
	return
	
    #plotitem.plot_var = surface_or_depth
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin =  -.01
    plotitem.pcolor_cmax =   .01
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotaxes.xlimits = [-1., 5.]
    plotaxes.ylimits = [15., 25.]

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.4,.4]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    #plotitem.plot_var = gaugetopo
    #plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel
        t = current_data.t 
        gaugeno = current_data.gaugeno
        #plot(t, 0*t, 'k')

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    # Figures for fgmax - max values on fixed grids
    #-----------------------------------------
    #otherfigure = plotdata.new_otherfigure(name='max flow depth', 
    #                fname='flowdepth.png')

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
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

    
