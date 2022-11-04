#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
is_cron = False
is_cron = bool( os.getenv('RUN_BY_CRON') )
from time import time
import datetime as dt
import numpy as np
from netCDF4 import Dataset
from metpy.units import units
import metpy.calc as mpcalc
import pygrib # import pygrib interface to grib_api
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.ndimage.filters import gaussian_filter
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
try: plt.style.use('mystyle')
except: pass

################################# LOGGING ####################################
import logging
import log_help
log_file = '.'.join( __file__.split('/')[-1].split('.')[:-1] ) + '.log'
# log_file = here+'/'+'.'.join( __file__.split('/')[-1].split('.')[:-1] ) + '.log'
lv = logging.DEBUG
logging.basicConfig(level=lv,
                 format='%(asctime)s %(name)s:%(levelname)s-%(message)s',
                 datefmt='%Y/%m/%d-%H:%M',
                 filename = log_file, filemode='w')
LG = logging.getLogger('main')
if not is_cron: log_help.screen_handler(LG, lv=lv)
LG.info(f'Starting: {__file__}')
##############################################################################

#### Colormaps ##################################################################
# Frentes
red = np.array([255,0,0,255])/255
white = np.array([255,255,255,0])/255
blue = np.array([0,0,255,255])/255
cols  = [(1-x)*blue+x*white for x in np.linspace(0,1,10)]
cols += [(1-x)*white+x*red  for x in np.linspace(0,1,10)]
cm_frentes = LinearSegmentedColormap.from_list(name='frentes',colors=cols)

# Clouds
col0 = np.array([1,1,1,0])
col1 = np.array([.1,.1,.1,.7])
cols = [(1-x)*col0+x*col1 for x in np.linspace(0,1,5)]
cm_clouds = LinearSegmentedColormap.from_list(name='clouds',colors=cols)
# Rain
col0 = np.array([0,0,0,0])
col1 = np.array([82,117,118,50])
col2 = np.array([163,233,235,120])
cols = [col0, col1, col2]
cm_rain = ListedColormap( [C/255 for C in cols] )


# def WTFF(fname,left,right,bottom,top):
def WTFF(inps):
   """
   fname: [str] GFS data file to plot
   left,right,bottom,top: [float] min,max longitude and min,max latitude of the
   data
   """
   fname,left,right,bottom,top = inps
   extent = left, right, bottom, top

   # Read GFS' grib2 file
   grbs = pygrib.open(fname)

   # Get dates
   t,T = str(grbs.message(1)).split(':')[-2:]
   T = int(T.split()[-1])
   dataDate = dt.datetime.strptime(str(T), '%Y%m%d%H%M')
   t = int(t.split()[-2])
   t = dt.timedelta(hours=t)
   fcstDate = dataDate + t 

   # U & V
   index = 0    # surface
   index = 36   # 850hPa
   Us = grbs.select(name='U component of wind')
   U  = Us[index]
   Vs = grbs.select(name='V component of wind')
   V  = Vs[index]
   # Temperature
   temps = grbs.select(name='Temperature')
   # index = 0    # surfeace
   # index = 35   # 850hPa
   temp = temps[index-1]
   # Rain
   R = grbs.select(name='Categorical rain')[0]
   lats_r, lons = R.latlons()
   lons_r = np.unwrap(lons, discont=180)
   # Pressure
   press = grbs.select(name='Pressure reduced to MSL')
   press = press[0]   # in Pa
   press.values = gaussian_filter(press.values, 3)
   # Clouds
   clouds = grbs.select(name='Total Cloud Cover')
   clouds = clouds[22]
   # Lat, Lon
   lats, lons = press.latlons()
   lons = np.unwrap(lons, discont=180)
   dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)


   # Fix units
   assert U['level'] == V['level'] == temp['level']
   isopress = temp['level'] * units(temp['pressureUnits'])

   temp = temp.values * units(temp['parameterUnits'])
   potential_temperature = mpcalc.potential_temperature(isopress, temp)

   F = mpcalc.frontogenesis(potential_temperature,
                            U.values * units(U['parameterUnits']),
                            V.values * units(V['parameterUnits']),
                            dx=dx, dy=dy)
   convert_to_per_100km_3h = 1000*100*3600*3
   fronto =  F.magnitude * convert_to_per_100km_3h

   ref_lon = -3.7
   ref_lat = 40.4
   # left,right,bottom,top = -18, 8, 29, 49
   # left,right,bottom,top = -40, 30, 20, 70

   orto = ccrs.PlateCarree()
   mercator = True
   if mercator:
      proj = ccrs.Mercator(ref_lon, bottom, top)
      fig = plt.figure(figsize=(14,9.15), frameon=False)
   else:
      proj = ccrs.Orthographic(ref_lon, ref_lat)
      fig = plt.figure(figsize=(14,8.6), frameon=False)
   # ax = plt.axes(projection=projection)
   # ax = fig.add_axes([0,0,0.99,1], projection=proj)
   # ax = fig.add_axes([.125, .11, .775, .77], projection=proj)
   ax = fig.add_axes([0,0,1,1], projection=proj)
   ax.set_extent(extent, crs=orto)

   color = np.array([130,42,27])/255
   # ax.stock_img()
   # ax.coastlines(resolution='50m', color='w', linewidth=2.3, zorder=100)
   ax.coastlines(resolution='50m', color=color, linewidth=1.7, zorder=101)
   C = ax.contour(lons, lats, press.values/100, levels=range(804,1204,4),
                                    colors='k', linewidths=.75,
                                    transform=orto)
   ax.clabel(C, inline=True, fontsize=10, fmt='%d')
   C = ax.contour(lons, lats, press.values/100, levels=range(804,1204,16),
                                    colors='k', linewidths=1.75, zorder=102,
                                    transform=orto)
   ax.clabel(C, inline=True, fontsize=10, fmt='%d')
   # print('Plotting frontogenesis')
   # print(np.nanmin(F))
   # print(np.nanmax(F))

   l = 1
   print('Plotting fronts')
   C = ax.contourf(lons[::l,::l],
                   lats[::l,::l],
                 fronto[::l,::l],
                   # vmin=-8, vmax=8,
                   np.arange(-6, 6.5, .5),
                   cmap=cm_frentes,
                   # cmap=plt.cm.bwr,
                   extend='both', transform=orto)
   print('Plotting clouds')
   ax.contourf(lons[::l,::l],
               lats[::l,::l],
      clouds.values[::l,::l],
               cmap=cm_clouds,
               vmin=0, vmax=100, transform=orto)
   print('Plotting wind')
   n = 5
   f = 2
   ax.barbs(lons[::n,::n],lats[::n,::n],
            U.values[::n,::n], V.values[::n,::n],
            color='k', length=4, #pivot='middle',
            sizes=dict(emptybarb=0.25/f, spacing=0.5/f, height=0.9/f),
            linewidth=0.75, transform=orto)
   print('Plotting rain')
   rain = ax.contourf(lons_r[::l,::l],
                      lats_r[::l,::l],
                    R.values[::l,::l], cmap=cm_rain,
                      transform=orto)
   # plt.colorbar(rain, orientation='horizontal', pad=0, aspect=50, extendrect=True)
   print('Ploted')
   ax.text(0.5, .99, f"{(fcstDate+UTCshift).strftime('%a %d %b')}",
           transform = ax.transAxes, va='top', ha='center', fontsize=30,
           bbox=dict(facecolor=(1,1,1,.8), ec='k', pad=10), zorder=1000)
   msg  = [f"valid: {(fcstDate+UTCshift).strftime('%Y-%m-%d %H:%M')}"]
   msg += [f"GFS: {(dataDate).strftime('%Y-%m-%d %H:%M')}"]
   msg += [f"plot: {dt.datetime.now().strftime('%Y-%m-%d %H:%M')}"]
   ax.text(0, .008, f'{isopress:~}',
           transform = ax.transAxes,
           bbox=dict(facecolor=(1,1,1,.9), ec='k', pad=5.0), zorder=1000)
   ax.text(1, .008, '\n'.join(msg),
           transform = ax.transAxes, ha='right',
           bbox=dict(facecolor=(1,1,1,.9), ec='k', pad=5.0), zorder=1000)
           # backgroundcolor=(1,1,1,.5), edgecolor='k')
   # ax.text(.008,.008,fcstDate+UTCshift,transform = ax.transAxes,
   #         backgroundcolor=(1,1,1,.5), edgecolor='k')
# cb = plt.colorbar(C, orientation='horizontal', pad=0, aspect=50, extendrect=True)

   auxx = press.values/100 - 1013
   auxx = gaussian_filter(press.values, 3)/100 - 1013
   aux = np.where(auxx>0, auxx, 0)
   th = 30
   for ii,p in enumerate( peak_local_max(aux, min_distance=th) ):
      i,j = p
      # print(ii,i,j,lats[i,j], lons[i,j], f'H{ii}')
      if press.values[i,j]/100 > 1022: text = 'A'
      else: text = 'a'
      ax.text(lons[i,j], lats[i,j], text,
              color='r', ha='center', va='center', weight='bold',
              fontsize=30,
              # bbox=dict(facecolor=(1,1,1,.9), ec='k', pad=5.0),
              transform=orto, zorder=999)
   th = 7
   aux = np.where(auxx<0, -auxx, 0)
   for ii,p in enumerate( peak_local_max(aux, min_distance=th) ):
      i,j = p
      if press.values[i,j]/100 < 980: text = 'B'
      else: text = 'b'
      # print(ii,i,j,lats[i,j], lons[i,j], f'L{ii}')
      ax.text(lons[i,j], lats[i,j], text,
              color='b', ha='center', va='center', weight='bold',
              # bbox=dict(facecolor=(1,1,1,.9), ec='k', pad=5.0),
              fontsize=30,
              transform=orto, zorder=999)

   # plt.show()
   print(f"Saving: {fcstDate.strftime('%Y%m%d_%H%M')}.png")
   fig.savefig(f"{fcstDate.strftime('%Y%m%d_%H%M')}.png")
   plt.close('all')


import generate_config as gen
import download_gfs_data as download
import gfs
folder = '.'
domain = 'Spain6_1'
days = 14
UTCshift = dt.datetime.now() - dt.datetime.utcnow()
UTCshift = dt.timedelta(hours = round(UTCshift.total_seconds()/3600))
start_date = gfs.get_GFS_calc(dt.datetime.now()) + UTCshift
end_date   = start_date + dt.timedelta(days=days)
# left = -75
# right = 45
# bottom = 15
# top = 70
x_ = 40
y_ = 20
lon0 = -9
lat0 = 43
left, right, bottom, top = lon0-1.3*x_, lon0+x_, lat0-y_, lat0+1.05*y_
# left, right, bottom, top = -180, 180, 0, 90
gen.main(folder, domain, start_date, end_date, timedelta=3,
         left=left, right=right,
         bottom=bottom, top=top, folder_out='~/Documents/storage')
print('Starting Download')
download.main()
print('Finished Download')

files = os.popen('ls dataGFS/gfs.*').read().strip().split()

from multiprocessing import Pool

parallel = True
if parallel:
   inputs = []
   for f in files:
      inputs.append((f,left,right,bottom,top))
   n = 5
   with Pool(n) as pool:
      Res = pool.map(WTFF, inputs)
   pool.close()
   pool.join()
   del pool
else:
   for fname in files:
      told = time()
      print(fname)
      WTFF((fname, left,right,bottom,top))
      print(time()-told)
      exit()

com = 'rm *.mp4 *.webm'
os.system(com)
com = 'ls -1 20*.png > files.txt'
os.system(com)
com = 'mencoder -nosound -ovc lavc -lavcopts vcodec=mpeg4 -o fronto.mp4 -mf type=png:fps=3 mf://@files.txt > /dev/null 2> /dev/null'
os.system(com)
com = 'rm files.txt'
os.system(com)
com = 'ffmpeg -i fronto.mp4 frontogenesis.webm'
os.system(com)
com = 'ffmpeg -i fronto.mp4 -an -vcodec libx264 -crf 23 frontogenesis.mp4'
os.system(com)

com = 'scp frontogenesis.mp4 olympus:CODES/web/RUNweb/assets/images/'
os.system(com)
com = 'scp frontogenesis.webm olympus:CODES/web/RUNweb/assets/images/'
os.system(com)

com = 'rm 20*.png'
os.system(com)
