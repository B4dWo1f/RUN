#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import matplotlib as mpl
mpl.use('Agg')

import log_help
import logging
LG = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.units import units
import matplotlib.lines as mlines
import metpy.calc as mpcalc
from matplotlib import gridspec
from colormaps import HEIGHTS


# import os
# here = os.path.dirname(os.path.realpath(__file__))
# import rasterio
# from rasterio.merge import merge
# from cartopy.feature import NaturalEarthFeature
# import cartopy.crs as ccrs
# import matplotlib.patheffects as PathEffects
# from matplotlib.colors import LightSource, BoundaryNorm

## Dark Theme ##################################################################
# #  COLOR = 'black'
# #  ROLOC = '#e0e0e0'
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 15.0
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['axes.labelsize'] = 'large' # fontsize of the x any y labels

# mpl.rcParams['axes.facecolor'] = (1,1,1,0)
# mpl.rcParams['figure.facecolor'] = (1,1,1,0)
# mpl.rcParams["savefig.facecolor"] = (1,1,1,0)

#  mpl.rcParams['text.color'] = ROLOC #COLOR
#  mpl.rcParams['axes.labelcolor'] = COLOR
#  mpl.rcParams['axes.facecolor'] = COLOR #'black'
#  mpl.rcParams['savefig.facecolor'] = COLOR #'black'
#  mpl.rcParams['xtick.color'] = COLOR
#  mpl.rcParams['ytick.color'] = COLOR
#  mpl.rcParams['axes.edgecolor'] = COLOR
#  mpl.rcParams['axes.titlesize'] = 20
################################################################################


@log_help.timer(LG)
def skewt_plot(p,tc,tdc,t0,date,u=None,v=None,fout='sounding.png',latlon='',title='',show=False):
   """
   h: heights
   p: pressure
   tc: Temperature [C]
   tdc: Dew point [C]
   date: date of the forecast
   u,v: u,v wind components
   adapted from:
   https://geocat-examples.readthedocs.io/en/latest/gallery/Skew-T/NCL_skewt_3_2.html#sphx-glr-gallery-skew-t-ncl-skewt-3-2-py
   """
   print('Checking units')
   if p.attrs['units'] != 'hPa':
      print('P wrong units')
      exit()
   if tc.attrs['units'] != 'degC':
      print('Tc wrong units')
      exit()
   if tdc.attrs['units'] != 'degC':
      print('Tdc wrong units')
      exit()
   if t0.attrs['units'] != 'degC':
      print('T0 wrong units')
      exit()
   if type(u) != type(None) and type(v) != type(None):
      if u.attrs['units'] != 'm s-1':
         print('Wind wrong units')
         exit()
   LG.info('Inside skewt plot')
   p = p.values
   tc = tc.values
   tdc = tdc.values
   t0 = t0.mean().values
   u = u.values * 3.6  # km/h
   v = v.values * 3.6  # km/h
   # Grid plot
   LG.info('creating figure')
   fig = plt.figure(figsize=(11, 12))
   LG.info('created figure')
   LG.info('creating axis')
   gs = gridspec.GridSpec(1, 4)
   fig.subplots_adjust(wspace=0.,hspace=0.)
   # ax1 = plt.subplot(gs[1:-1,0])
   # Adding the "rotation" kwarg will over-ride the default MetPy rotation of
   # 30 degrees for the 45 degree default found in NCL Skew-T plots
   LG.info('created axis')
   LG.info('Creatin SkewT')
   skew = SkewT(fig, rotation=45, subplot=gs[0,0:-1])
   ax = skew.ax
   LG.info('Created SkewT')

   if len(latlon) > 0:
       ax.text(0,1, latlon, va='top', ha='left', color='k',
                     fontsize=12, bbox=dict(boxstyle="round",
                                            ec=None, fc=(1., 1., 1., 0.9)),
                     zorder=100, transform=ax.transAxes)
   # Plot the data, T and Tdew vs pressure
   skew.plot(p, tc, 'C3')
   skew.plot(p, tdc, 'C0')
   LG.info('plotted dew point and sounding')

   # LCL
   lcl_pressure, lcl_temperature = mpcalc.lcl(p[0]*units.hPa,
                                              t0*units.degC,
                                              tdc[0]*units.degC)
   skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

   # Calculate the parcel profile  #XXX units workaround
   parcel_prof = mpcalc.parcel_profile(p* units.hPa,
                                       t0 * units.degC,
                                       tdc[0]* units.degC).to('degC')
   # Plot the parcel profile as a black line
   skew.plot(p, parcel_prof, 'k', linewidth=2)
   LG.info('plotted parcel profile')

   # shade CAPE and CIN
   skew.shade_cape(p* units.hPa,
                   tc * units.degC, parcel_prof)
   skew.shade_cin(p * units.hPa,
                  tc * units.degC,
                  parcel_prof,
                  tdc * units.degC)
   LG.info('plotted CAPE and CIN')

   if type(u) != type(None) and type(v) != type(None):
      LG.info('Plotting wind')
      ax2 = plt.subplot(gs[0,-1], sharey=ax)
      # ax2.yaxis.set_visible(False)
      ax2.yaxis.tick_right()
      ax2.xaxis.tick_top()
      wspd = np.sqrt(u*u + v*v)
      ax2.scatter(wspd, p, c=p, cmap=HEIGHTS) #'viridis_r')
      gnd = mpcalc.pressure_to_height_std(np.array(p[0])*units.hPa)
      gnd = gnd.to('m')
      ax2.axhline(p[0],c='k',ls='--')
      ax2.text(56,p[0],f'{int(gnd.magnitude)}m',horizontalalignment='right')
      ax2.set_xlim(0,56)
      ax2.set_xlabel('Wspeed (km/h)')
      ax2.grid()
      def p2h(x):
         """
         x in hPa
         """
         y = mpcalc.pressure_to_height_std(np.array(x)*units.hPa)
         # y = y.metpy.convert_units('m')
         y = y.to('m')
         return y.magnitude
      def h2p(x):
         """
         x in m
         """
         # x = x.values
         y = mpcalc.height_to_pressure_std(np.array(x)*units.m)
         # y = y.metpy.convert_units('hPa')
         y = y.to('hPa')
         return y.magnitude
      # print('------')
      # print(p[0])
      # print('---')
      # print(p2h(p[0]))
      # print('---')
      # print(h2p(p2h(p[0])))
      # print('------')
      # exit()
      ax2y = ax2.secondary_yaxis(1.,functions=(p2h,h2p))
      ax2y.set_ylabel('height (m)')
      ax2y.set_yticks(np.array([1000,2000,3000,4000,6000,6000,10000]))
      # ax2y.set_yticks([1,2,3,4])
      ax2.set_xticks([0,5,10,15,20,30,40,50])
      ax2.set_xticklabels(['0','','10','','20','30','40','50'])
      # from matplotlib.ticker import ScalarFormatter
      # ax2y.set_major_formatter(ScalarFormatter())
      plt.setp(ax2.get_yticklabels(), visible=False)
      # Hodograph
      ax_hod = inset_axes(ax2, '110%', '30%', loc=1)
      ax_hod.set_yticklabels([])
      L = 60
      ax_hod.text(  0, L-5,'N', horizontalalignment='center',
                               verticalalignment='center')
      ax_hod.text(L-5,  0,'E', horizontalalignment='center',
                               verticalalignment='center')
      ax_hod.text(-(L-5),0 ,'W', horizontalalignment='center',
                               verticalalignment='center')
      ax_hod.text(  0,-(L-5),'S', horizontalalignment='center',
                               verticalalignment='center')
      h = Hodograph(ax_hod, component_range=L)
      h.add_grid(increment=20)
      h.plot_colormapped(-u, -v, p, cmap=HEIGHTS)  #'viridis_r')  # Plot a line colored by pressure (altitude)
      LG.info('Plotted wind')



      # # print('=-=-=-=-=-=-=')
      # # print(ax2.get_yscale())
      # # print(ax2y.get_yscale())
      # # print('=-=-=-=-=-=-=')
      # # ax2y.yaxis.tick_right()
      # # # ax2y.set_ylim(*reversed(ax.get_ylim()))
      # # ax2y.set_ylim(*ax.get_ylim())
      # # # calc pressure to height
      # # locs = ax2.get_yticks()
      # # labels = [mpcalc.pressure_to_height_std(h*units.hPa) for h in locs]
      # # labels = [round(l.to('m'),1) for l in labels]
      # # for xp,xh in zip(locs,labels):
      # #    print(xp,xh)
      # # ax2y.set_yticks(locs)
      # # ax2y.set_yticklabels([f'{int(l.magnitude)}' for l in labels])
   # #    # Plot only every n windbarb
   # #    n = 3
   # #    skew.plot_barbs(pressure=p[::n].values * units.hPa,
   # #                    u=u[::n], v=v[::n],
   # #                    xloc=1.05, fill_empty=True,
   # #                    sizes=dict(emptybarb=0.075, width=0.1, height=0.2))
   # #    # Draw line underneath wind barbs
   # #    line = mlines.Line2D([1.05, 1.05], [0, 1], color='gray', linewidth=0.5,
   # #                         transform=ax.transAxes,
   # #                         dash_joinstyle='round',
   # #                         clip_on=False,
   # #                         zorder=0)
   # #    ax.add_line(line)

   # Add relevant special lines
   # Choose starting temperatures in Kelvin for the dry adiabats
   LG.info('Plot adiabats, and other grid lines')
   t0 = units.K * np.arange(243.15, 473.15, 10)
   skew.plot_dry_adiabats(t0=t0, linestyles='solid', colors='gray', linewidth=1)

   # Choose temperatures for moist adiabats
   t0 = units.K * np.arange(281.15, 306.15, 4)
   msa = skew.plot_moist_adiabats(t0=t0,
                                  linestyles='solid',
                                  colors='lime',
                                  linewidths=.75)

   # Choose mixing ratios
   w = np.array([0.001, 0.002, 0.003, 0.005, 0.008, 0.012, 0.020]).reshape(-1, 1)

   # Choose the range of pressures that the mixing ratio lines are drawn over
   p_levs = units.hPa * np.linspace(1000, 400, 7)
   skew.plot_mixing_lines(mixing_ratio=w, pressure=p_levs,
                          linestyle='dotted',linewidths=1, colors='lime')

   LG.info('Plotted adiabats, and other grid lines')

   skew.ax.set_ylim(1000, 150)
   skew.ax.set_xlim(-30, 30)
   # skew.ax.set_xlim(-30, 40)
   # XXX gvutil not installed
   # gvutil.set_titles_and_labels(ax, maintitle="ATS Rawinsonde: degC + Thin wind")
   # Set axes limits and ticks
   # gvutil.set_axes_limits_and_ticks(ax=ax, xlim=[-30, 50],
   #                 yticks=[1000, 850, 700, 500, 400, 300, 250, 200, 150, 100])

   # Change the style of the gridlines
   ax.grid(True, which='major', axis='both',
            color='tan', linewidth=1.5, alpha=0.5)

   ax.set_xlabel("Temperature (C)")
   ax.set_ylabel("P (hPa)")
   if len(title) == 0:
       title = f"{(date).strftime('%d/%m/%Y-%H:%M')} (local time)"
       ax.set_title(title,fontsize=20)
   else:
       ax.set_title(title, fontsize=20)
   if show: plt.show()
   LG.info('saving')
   fig.savefig(fout, bbox_inches='tight', pad_inches=0.1, dpi=150, quality=90)
   LG.info('saved')
   plt.close('all')
