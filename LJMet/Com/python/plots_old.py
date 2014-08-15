#!/usr/bin/env python
#########################################################################
#
# plots_old.py
#
# plots_old module
#
# Usage: ./analyze.py --old-plot <number>
#
#
# Authors:
#          Gena Kukartsev, July 2012
#
#
#########################################################################



import sys



setlogscale = False



legend = '[plots]:'


def MakePlot(options):
    #
    # old-style plotting, for backward compatibility
    #
    
    legend = '[plots.MakePlot]:'

    if options.plot_number==None:
        print legend, "no plot specified, exiting"
        sys.exit(-1)


    # HT
    if int(options.plot_number)==1:
        import DataLoader as dl
        import Plotter as pl
        ds = dl.LoadData()
        pl.DrawPlot('ht', ds, 20, 200, 1200, './',
                    "H_{T} [GeV]", "arbitrary units",
                    xLegend = 0.68, yLegend = 0.66,
                    format=options.format,
                    legendWidth = 0.25, legendHeight = 0.15,
                    xStat=0.60, yStat=0.87,
                    cut='',
                    ymin=0, ymax=0.2,
                    log_scale = setlogscale
                    )
        
    # mlb1
    if int(options.plot_number)==2:
        import DataLoader as dl
        import Plotter as pl
        ds = dl.LoadData()
        pl.DrawPlot('mlb1', ds, 20, 0, 300, './',
                    "m_{lb} [GeV]", "arbitrary units",
                    xLegend = 0.68, yLegend = 0.66,
                    format=options.format,
                    legendWidth = 0.25, legendHeight = 0.15,
                    xStat=0.60, yStat=0.87,
                    cut='',
                    ymin=0, ymax=0.2,
                    log_scale = setlogscale
                    )
        
    # mlb2
    if int(options.plot_number)==3:
        import DataLoader as dl
        import Plotter as pl
        ds = dl.LoadData()
        pl.DrawPlot('mlb2', ds, 20, 0, 300, './',
                    "m_{lb} [GeV]", "arbitrary units",
                    xLegend = 0.68, yLegend = 0.66,
                    format=options.format,
                    legendWidth = 0.25, legendHeight = 0.15,
                    xStat=0.60, yStat=0.87,
                    cut='',
                    ymin=0, ymax=0.2,
                    log_scale = setlogscale
                    )
        
    # bjet1 pt
    if int(options.plot_number)==4:
        import DataLoader as dl
        import Plotter as pl
        ds = dl.LoadData()
        pl.DrawPlot('bjet_1_pt', ds, 20, 0, 300, './',
                    "p_{T}^{b} [GeV]", "arbitrary units",
                    xLegend = 0.68, yLegend = 0.66,
                    format=options.format,
                    legendWidth = 0.25, legendHeight = 0.15,
                    xStat=0.60, yStat=0.87,
                    cut='',
                    ymin=0, ymax=0.3,
                    log_scale = setlogscale
                    )
        
    # bjet2 pt
    if int(options.plot_number)==5:
        import DataLoader as dl
        import Plotter as pl
        ds = dl.LoadData()
        pl.DrawPlot('bjet_2_pt', ds, 20, 0, 300, './',
                    "p_{T}^{b} [GeV]", "arbitrary units",
                    xLegend = 0.68, yLegend = 0.66,
                    format=options.format,
                    legendWidth = 0.25, legendHeight = 0.15,
                    xStat=0.60, yStat=0.87,
                    cut='',
                    ymin=0, ymax=0.5,
                    log_scale = setlogscale
                    )
        


