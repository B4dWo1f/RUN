#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
here = os.path.dirname(os.path.realpath(__file__))
import numpy as np
import rasterio
from rasterio.merge import merge
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
# try: plt.style.use('rasp')
# except: pass
import matplotlib.patheffects as PathEffects
from matplotlib.colors import LightSource, BoundaryNorm

## Dark Theme ##################################################################
import matplotlib as mpl
#  COLOR = 'black'
#  ROLOC = '#e0e0e0'
mpl.rcParams['axes.facecolor'] = (1,1,1,0)
mpl.rcParams['figure.facecolor'] = (1,1,1,0)
mpl.rcParams["savefig.facecolor"] = (1,1,1,0)

#  mpl.rcParams['text.color'] = ROLOC #COLOR
#  mpl.rcParams['axes.labelcolor'] = COLOR
#  mpl.rcParams['axes.facecolor'] = COLOR #'black'
#  mpl.rcParams['savefig.facecolor'] = COLOR #'black'
#  mpl.rcParams['xtick.color'] = COLOR
#  mpl.rcParams['ytick.color'] = COLOR
#  mpl.rcParams['axes.edgecolor'] = COLOR
#  mpl.rcParams['axes.titlesize'] = 20
################################################################################


def setup_plot(ref_lat,ref_lon,left,right,bottom,top,transparent=True):
   orto = ccrs.PlateCarree()
   projection = ccrs.LambertConformal(ref_lon,ref_lat)
   # projection = ccrs.RotatedPole(pole_longitude=-150, pole_latitude=37.5)

   extent = left, right, bottom, top
   fig = plt.figure(figsize=(11,9)) #, frameon=False)
   # ax = plt.axes(projection=projection)
   ax = fig.add_axes([0,0,0.99,1],projection=projection)
   ax.set_extent(extent, crs=orto)
   if transparent:
       # ax.outline_patch.set_visible(False)
       ax.background_patch.set_visible(False)
   return fig,ax,orto

def save_figure(fig,fname,dpi=150, quality=90):
   fig.savefig(fname, transparent=True, bbox_inches='tight', pad_inches=0,
                      dpi=dpi, quality=quality)

def terrain_plot(reflat,reflon,left,right,bottom,top):
   fig, ax, orto = setup_plot(reflat,reflon,left,right,bottom,top)
### RASTER ###################################################################
   files = os.popen('ls geb*').read().strip().splitlines()
   srcs = [rasterio.open(fname, 'r') for fname in files]
   # fname0 = files[0]
   # fname1 = files[1]
   # src0 = rasterio.open(fname0, 'r')
   # src1 = rasterio.open(fname1, 'r')
   D = 2
   mosaic, out_trans = merge(srcs, (left-D, bottom-D, right+D, top+D))
   terrain = mosaic[0,:,:]
   ls = LightSource(azdeg=315, altdeg=45)
   ve = 0.7
   terrain = ls.hillshade(terrain, vert_exag=ve)
   ax.imshow(terrain, extent=(left-D, right+D, bottom-D, top+D),
                      origin='upper', cmap='gray',
                      aspect='equal', interpolation='lanczos',
                      zorder=0, transform=orto)
   return fig,ax,orto

def parallel_and_meridian(fig,ax,orto,left,right,bottom,top,nx=1,ny=1):
   lcs = 'k--'
   D = 1
   # Plotting meridian
   for x in range(int(left-D), int(right+D)):
      if x%nx ==0:
         ax.plot([x,x],[bottom-D,top+D], lcs, transform=orto)
   # Plotting parallels
   for y in range(int(bottom-D), int(top+D)):
      if y%ny == 0:
         ax.plot([left-D,right+D],[y,y], lcs, transform=orto)
   return fig,ax,orto

def rivers_plot(fig,ax,orto):
   rivers = NaturalEarthFeature('physical',
                                'rivers_lake_centerlines_scale_rank',
                                '10m', facecolor='none')
   ax.add_feature(rivers, lw=2 ,edgecolor='C0',zorder=50)
   rivers = NaturalEarthFeature('physical', 'rivers_europe',
                                '10m', facecolor='none')
   ax.add_feature(rivers, lw=2 ,edgecolor='C0',zorder=50)
   lakes = NaturalEarthFeature('physical', 'lakes', '10m') 
   ax.add_feature(lakes, lw=2 ,edgecolor='C0',zorder=50)
   lakes = NaturalEarthFeature('physical', 'lakes_historic', '10m')
   ax.add_feature(lakes, lw=2 ,edgecolor='C0',zorder=50)
   lakes = NaturalEarthFeature('physical', 'lakes_pluvial', '10m')
   ax.add_feature(lakes, lw=2 ,edgecolor='C0',zorder=50)
   lakes = NaturalEarthFeature('physical', 'lakes_europe', '10m')
   ax.add_feature(lakes, lw=2 ,edgecolor='C0',zorder=50)
   return fig,ax,orto

def sea_plot(fig,ax,orto):
   sea = NaturalEarthFeature('physical', 'bathymetry_all', '10m') #, facecolor='none')
   ax.add_feature(sea, lw=2) # ,edgecolor='C0',zorder=50)
   return fig,ax,orto

def ccaa_plot(fig,ax,orto):
   provin = NaturalEarthFeature('cultural', 'admin_1_states_provinces_lines',
                                            '10m', facecolor='none')
   country = NaturalEarthFeature('cultural', 'admin_0_countries', '10m',
                                                   facecolor='none')
   ax.add_feature(provin, lw=2 ,edgecolor='k',zorder=50)
   ax.add_feature(country,lw=2.3, edgecolor='k', zorder=51)
   return fig,ax,orto

def csv_plot(fig,ax,orto, fname):
   Yt,Xt = np.loadtxt(fname,usecols=(0,1),delimiter=',',unpack=True)
   names = np.loadtxt(fname,usecols=(2,),delimiter=',',dtype=str)
   ax.scatter(Xt,Yt,s=40,c='r',marker='x',transform=orto,zorder=53)
   return fig,ax,orto

def csv_names_plot(fig,ax,orto, fname):
   # Cities
   Yt,Xt = np.loadtxt(fname,usecols=(0,1),delimiter=',',unpack=True)
   names = np.loadtxt(fname,usecols=(2,),delimiter=',',dtype=str)
   for x,y,name in zip(Xt,Yt,names):
      txt = ax.text(x,y,name, horizontalalignment='center',
                              verticalalignment='center',
                              color='k',fontsize=13,
                              transform=orto,zorder=52)
      txt.set_path_effects([PathEffects.withStroke(linewidth=5,
                                                   foreground='w')])
   return fig,ax,orto



def scalar_plot(fig,ax,orto, lons,lats,prop, delta,vmin,vmax,cmap,
                                                 levels=[], creation_date=''):
   if len(levels) > 0:
      norm = BoundaryNorm(levels,len(levels))
   else:
      levels = np.arange(vmin,vmax,delta)
      norm = None
   try:
      C = ax.contourf(lons,lats,prop, levels=levels, extend='max',
                                      antialiased=True, norm=norm,
                                      cmap=cmap, vmin=vmin, vmax=vmax,
                                      zorder=10, transform=orto)
   except:
       # LG.warning('NaN values found, unable to plot')
       C = None
   if len(creation_date) > 0:
       ax.text(1,0., creation_date, va='bottom', ha='right', color='k',
                     fontsize=12, bbox=dict(boxstyle="round",
                                            ec=None, fc=(1., 1., 1., 0.9)),
                     zorder=100, transform=ax.transAxes)
   # cbaxes = fig.add_axes([0., -0.05, 0.999, 0.03]) 
   # fig.colorbar(C, cax=cbaxes,orientation="horizontal")
   # cbaxes.set_xlabel('test')
   # return C, cbaxes
   return C

def vector_plot(fig,ax,orto,lons,lats,U,V, dens=1.5,color=(0,0,0,0.75)):
   x = lons[0,:].values
   y = lats[:,0].values
   # try:  #XXX workaround for DrJack's wind calculations
   #     U = U.values
   #     V = V.values
   # except: pass
   ax.streamplot(x,y, U,V, color=color, linewidth=1, density=dens,
                           arrowstyle='->',arrowsize=2.5,
                           zorder=11,
                           transform=orto)



from mpl_toolkits.axes_grid1 import make_axes_locatable
def plot_colorbar(cmap,delta=4,vmin=0,vmax=60,levels=None,name='cbar',
                                        units='',fs=18,norm=None,extend='max'):
   fig, ax = plt.subplots()
   img = np.random.uniform(vmin,vmax,size=(40,40))
   if len(levels) == 0:
      levels=np.arange(vmin,vmax,delta)
   img = ax.contourf(img, levels=levels,
                          extend=extend,
                          antialiased=True,
                          cmap=cmap,
                          norm=norm,
                          vmin=vmin, vmax=vmax)
   plt.gca().set_visible(False)
   divider = make_axes_locatable(ax)
   cax = divider.new_vertical(size="2.95%", pad=0.25, pack_start=True)
   fig.add_axes(cax)
   cbar = fig.colorbar(img, cax=cax, orientation="horizontal")
   cbar.ax.set_xlabel(units,fontsize=fs)
   fig.savefig(f'{name}.png', transparent=True,
                              bbox_inches='tight', pad_inches=0)


def skewt_plot(p,tc,tdc,t0,date,u=None,v=None,show=False):
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
   p = p.values
   tc = tc.values
   tdc = tdc.values
   t0 = t0.values
   u = u.values * 3.6  # km/h
   v = v.values * 3.6  # km/h
   from metpy.plots import SkewT
   from metpy.units import units
   import matplotlib.lines as mlines
   import metpy.calc as mpcalc
   from matplotlib import gridspec
   fig = plt.figure(figsize=(11, 12))
   gs = gridspec.GridSpec(1, 4)
   fig.subplots_adjust(wspace=0.,hspace=0.)
   # ax1 = plt.subplot(gs[1:-1,0])
   # Adding the "rotation" kwarg will over-ride the default MetPy rotation of
   # 30 degrees for the 45 degree default found in NCL Skew-T plots
   skew = SkewT(fig, rotation=45, subplot=gs[0,0:-1])
   ax = skew.ax

   # Plot the data, T and Tdew vs pressure
   skew.plot(p, tc, 'C3')
   skew.plot(p, tdc, 'C0')

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

   # shade CAPE and CIN
   skew.shade_cape(p* units.hPa,
                   tc * units.degC, parcel_prof)
   skew.shade_cin(p * units.hPa,
                  tc * units.degC,
                  parcel_prof,
                  tdc * units.degC)

   if type(u) != type(None) and type(v) != type(None):
      ax2 = plt.subplot(gs[0,-1], sharey=ax)
      # ax2.yaxis.set_visible(False)
      ax2.yaxis.tick_right()
      ax2.xaxis.tick_top()
      wspd = np.sqrt(u*u + v*v)
      ax2.plot(wspd, p)
      ax2.set_xlim(0,40)
      ax2.set_xlabel('Wspeed (km/h)')
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
      ax2y.set_yticks([10**3,5*10**3,10**4])
      ax2y.set_yticks([1,2,3,4])
      ax2.set_xticks([0,5,10,15,20,30,40])
      # from matplotlib.ticker import ScalarFormatter
      # ax2y.set_major_formatter(ScalarFormatter())
      plt.setp(ax2.get_yticklabels(), visible=False)


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
   ax.set_title(f"{(date).strftime('%d/%m/%Y-%H:%M')} (local time)")
   if show: plt.show()
   return fig,ax

# Backup
# def skewt_plot(p,tc,tdc,date,u=None,v=None,show=False):
#    # adapted from https://geocat-examples.readthedocs.io/en/latest/gallery/Skew-T/NCL_skewt_3_2.html#sphx-glr-gallery-skew-t-ncl-skewt-3-2-py
#    from metpy.plots import SkewT
#    from metpy.units import units
#    import matplotlib.lines as mlines
#    import metpy.calc as mpcalc
#    fig = plt.figure(figsize=(9, 12))
#    # Adding the "rotation" kwarg will over-ride the default MetPy rotation of
#    # 30 degrees for the 45 degree default found in NCL Skew-T plots
#    skew = SkewT(fig, rotation=45)
#    ax = skew.ax

#    # Plot the data, T and Tdew vs pressure
#    skew.plot(p, tc, 'r')
#    skew.plot(p, tdc, 'b')

#    # LCL
#    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], tc[0], tdc[0])
#    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

#    # Calculate the parcel profile  #XXX units workaround
#    parcel_prof = mpcalc.parcel_profile(p.values * units.hPa,
#                                        tc[0] * units.degC,
#                                        tdc[0].values * units.degC).to('degC')
#    # Plot the parcel profile as a black line
#    skew.plot(p, parcel_prof, 'k', linewidth=2)

#    # shade CAPE and CIN
#    skew.shade_cape(p.values * units.hPa,
#                    tc.values * units.degC, parcel_prof)
#    skew.shade_cin(p.values * units.hPa,
#                   tc.values * units.degC,
#                   parcel_prof,
#                   tdc.values * units.degC)

#    if type(u) != type(None) and type(v) != type(None):
#       # Plot only every n windbarb
#       n = 3
#       skew.plot_barbs(pressure=p[::n].values * units.hPa,
#                       u=u[::n], v=v[::n],
#                       xloc=1.05, fill_empty=True,
#                       sizes=dict(emptybarb=0.075, width=0.1, height=0.2))
#       # Draw line underneath wind barbs
#       line = mlines.Line2D([1.05, 1.05], [0, 1], color='gray', linewidth=0.5,
#                            transform=ax.transAxes,
#                            dash_joinstyle='round',
#                            clip_on=False,
#                            zorder=0)
#       ax.add_line(line)

#    # Add relevant special lines
#    # Choose starting temperatures in Kelvin for the dry adiabats
#    t0 = units.K * np.arange(243.15, 473.15, 10)
#    skew.plot_dry_adiabats(t0=t0, linestyles='solid', colors='gray', linewidth=1)

#    # Choose temperatures for moist adiabats
#    t0 = units.K * np.arange(281.15, 306.15, 4)
#    msa = skew.plot_moist_adiabats(t0=t0,
#                                   linestyles='solid',
#                                   colors='lime',
#                                   linewidths=.75)

#    # Choose mixing ratios
#    w = np.array([0.001, 0.002, 0.003, 0.005, 0.008, 0.012, 0.020]).reshape(-1, 1)

#    # Choose the range of pressures that the mixing ratio lines are drawn over
#    p_levs = units.hPa * np.linspace(1000, 400, 7)
#    skew.plot_mixing_lines(mixing_ratio=w, pressure=p_levs,
#                           linestyle='dotted',linewidths=1, colors='lime')

#    skew.ax.set_ylim(1000, 150)
#    skew.ax.set_xlim(-30, 30)
#    # skew.ax.set_xlim(-30, 40)
#    # XXX gvutil not installed
#    # gvutil.set_titles_and_labels(ax, maintitle="ATS Rawinsonde: degC + Thin wind")
#    # Set axes limits and ticks
#    # gvutil.set_axes_limits_and_ticks(ax=ax, xlim=[-30, 50],
#    #                 yticks=[1000, 850, 700, 500, 400, 300, 250, 200, 150, 100])

#    # Change the style of the gridlines
#    ax.grid(True, which='major', axis='both',
#             color='tan', linewidth=1.5, alpha=0.5)

#    ax.set_xlabel("Temperature (C)")
#    ax.set_ylabel("P (hPa)")
#    ax.set_title(f"{(date).strftime('%d/%m/%Y-%H:%M')} (local time)")
#    if show: plt.show()
#    return fig,ax
