""" 
Set up the plot figures, axes, and items to be done for each frame.
This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters. 
""" 
try:
    from setplotfg import setplotfg
except:
    print "Did not find setplotfg.py"
    setplotfg = None
    
import os

infile='../energy_estimate.txt'
fname='energy_estimate.txt'
if os.path.exists(fname):
    os.remove(fname)
fout=open(fname,'w')

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


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from pyclaw.plotters import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    #plotaxes.scaled = True
    #plotaxes.afteraxes = addgauges

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

    def aa(current_data):
        from numpy import where,sqrt
        q = current_data.q
        h = q[0,:,:]
        u = where(h>0.,q[1,:,:]/h,0.)
        t = current_data.t
        t = (t-2.)*sqrt(9.81)
        dx= current_data.dx
        eta = q[3,:,:]
        b = h-eta
        energy=where(h>0., .5*h*u**2 + .5*9.81*eta**2 , 0.)
        t_energy=sum(dx*energy)
        fout=open(infile,'a')
        fout.write("%20.10e %20.10e \n" % (t, t_energy[0]))
	return

    def aa_bous(current_data):
        import pylab
        from numpy import where,sqrt,zeros,max
        x = current_data.x[:,0]
        q = current_data.q
        h = q[0,:,:]
        u = where(h>0.,q[1,:,:]/h,0.)
        t = current_data.t
        t = -6.+t*sqrt(9.81)
        dx= current_data.dx
        eta = q[3,:,:]
        eta = where(h>0., eta, 0.)
        m = len(q[0,:,0])
        b = h-eta
        dudx= 0.*q[0,:,:]
        dbdx= 0.*q[0,:,:]
        for i in range(1,m-1):
            #dudx[i,:]=(u[i+1,:]-u[i-1,:])/2./dx
            #dbdx[i,:]=(b[i+1,:]-b[i-1,:])/2./dx
            dudx[i,:]=(u[i+1,:]-u[i,:])/dx
            dbdx[i,:]=(b[i+1,:]-b[i,:])/dx
        e_0 = where(h>0., .5*h*u**2 + .5*9.81*eta**2 , 0.)
        #e_1 = where(h>0., b**3/6.*dudx**2 + b**2/2.*u*dudx*dbdx , 0.)
        b_param = 0. #1./15.
        e_1 = where(h>0., (1./6.+b_param)*h**3*dudx**2 + (.5+0.)*h**2*u*dudx*dbdx + (.5+0.)*dbdx**2*u**2*h, 0.)
        e_t = e_0 + e_1
        te_0=sum(dx*e_0)
        te_1=sum(dx*e_1)
        te_t=sum(dx*e_t)
        eta_slice = eta[:,1]
        eta_slice = map(lambda x: (x), eta_slice)
        hmax=max(eta_slice)
        xloc= x[eta_slice.index(hmax)]
        eta_b=max(eta_slice)/(xloc/19.95)
        fout=open(infile,'a')
        fout.write("%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e \n" % (t, te_0[1], te_1[1],  te_t[1], eta_b, xloc))
        pylab.title('Surface at t = %4.2f' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
	return
	
    #plotitem.plot_var = surface_or_depth
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin =  -.1
    plotitem.pcolor_cmax =   .1
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

    plotaxes.afteraxes = aa_bous

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Slice', figno=1)
    plotfigure.kwargs = {'figsize': (9,3)}
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
    plotaxes.xlimits = [-10, 30]
    plotaxes.ylimits = [-.1, .5]
    
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
        t = -6.+t*np.sqrt(9.81) # Dimensionless
        pylab.title('Surface at t = %4.2f' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Energy', figno=21)
    plotfigure.kwargs = {'figsize': (12,5)}
    plotfigure.show = False

    def x_energy(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        dx= current_data.dx
        q = current_data.q
        m = len(y)
        m = int(m/2)
        h = q[0,:,m]
        hu= q[1,:,m]
        from numpy import where
        u = where(h>0., hu/h, 0.)
        eta = q[3,:,m]
        g = 9.81
        from numpy import where
        eta = q[3,:,m]
        k = len(h)
        b = h-eta
        dudx= 0.*q[0,:,m]
        dbdx= 0.*q[0,:,m]
        for i in range(k-2):
            dudx[i+1]=(u[i+2]-u[i])/2./dx
            dbdx[i+1]=(b[i+2]-b[i])/2./dx
        E=where(h>0., where(u<0., .5*h*u**2 + .5*9.81*eta**2+b**3/6.*dudx**2+b**2/2.*u*dudx*dbdx, 0.) , 0.)
        #E = where(h>0.,.5*hu**2/h + .5*g*eta**2 ,0.)
        return x,E
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('water surface')
    plotaxes.title = 'eta'
    plotaxes.xlimits = [-5.,40.]
    plotaxes.ylimits = [-.1,1.]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_energy
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    def fixup(current_data):
        import pylab
        import numpy as np
        t = current_data.t
        t = (t-2.)*np.sqrt(9.81)  # dimensionless
        pylab.title('Surface at t = %4.2f' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Bouss vs SW', figno=3)
    plotfigure.kwargs = {'figsize': (12,4)}
    plotfigure.show = False

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
    plotaxes.xlimits = [-5.,10.]
    plotaxes.ylimits = [-.2, .6]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_sw_dx005'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'r'
    plotitem.plotstyle = '--'
    plotitem.kwargs = {'linewidth':1}
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_bous_dx005'
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
    
    def add_legend(current_data):
        from pylab import plot, legend
        legend(('SWE','Boussinesq'),loc='lower left')
        
    plotaxes.afteraxes = add_legend


    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Bouss', figno=4)
    plotfigure.kwargs = {'figsize': (12,4)}
    plotfigure.show = False

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
    plotaxes.xlimits = [-5.,5]
    plotaxes.ylimits = [-.1, .6]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_bouss/dx_10'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'c'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_bouss/dx_05'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'g'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_bouss/dx_02'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'r'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_bouss/dx_01'
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
        from pylab import plot, legend
        legend(('dx = 0.1','0.05','0.02','0.01'),loc='lower left')
        import numpy as np
        t = current_data.t
        t = (t-2.)*np.sqrt(9.81)  # dimensionless
        pylab.title('Surface at t = %4.2f' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)

    plotaxes.afteraxes = fixup


    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Convergence SWE', figno=5)
    plotfigure.kwargs = {'figsize': (12,4)}
    plotfigure.show = False

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
    plotaxes.title = 'line'
    plotaxes.xlimits = [-10.,10]
    plotaxes.ylimits = [-.2, .5]
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_sw_dx02'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'r'
    plotitem.plotstyle = ':'
    plotitem.kwargs = {'linewidth':1}
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_sw_dx01'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'g'
    plotitem.plotstyle = '--'
    plotitem.kwargs = {'linewidth':1}

    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_sw_dx005'
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
    
    def add_legend(current_data):
        from pylab import plot, legend
        legend(('dx = 0.2','0.1','0.05'),loc='lower left')
        
    plotaxes.afteraxes = add_legend


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
                    type='each_gauge')
                    
    plotfigure.show = False
    plotfigure.kwargs = {'figsize': (10,4)}

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.01, .2]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.outdir='../_output_uni'
    #plotitem.plot_var = 3
    #plotitem.plotstyle = 'b-'
    #plotitem.kwargs = {'linewidth':1}
    
    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir='../_output'
    plotitem.plot_var = 3
    plotitem.plotstyle = 'r-'
    plotitem.kwargs = {'linewidth':1}
    
    # Plot surface as blue curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.outdir='../_output_amr'
    #plotitem.plot_var = 3
    #plotitem.plotstyle = 'b-'
    #plotitem.kwargs = {'linewidth':1}

    # Plot topo as green curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    #plotitem.plot_var = gaugetopo
    #plotitem.plotstyle = 'g-'
    
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        plot(t, 0*t, 'k')
        
    def add_amrlegend(current_data):
        from pylab import plot, legend
        legend(('Uniform','AMR'),loc='upper right')#prop={'size':6})
        #plot(t, 0*t, 'k')

    #plotaxes.afteraxes = add_amrlegend

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

    
