#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
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
import matplotlib.pyplot as plt
try: plt.style.use('mystyle')
except: pass


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

   # F = mpcalc.frontogenesis(potential_temperature,
   #                          U.values * units(U['parameterUnits']),
   #                          V.values * units(V['parameterUnits']),
   #                          dx=dx, dy=dy)
   # convert_to_per_100km_3h = 1000*100*3600*3
   # fronto =  F.magnitude * convert_to_per_100km_3h

   ref_lon = -3.7
   ref_lat = 40.4
   # left,right,bottom,top = -18, 8, 29, 49
   # left,right,bottom,top = -40, 30, 20, 70

   orto = ccrs.PlateCarree()
   proj = ccrs.Mercator(ref_lon, bottom, top)
   # proj = ccrs.Orthographic(ref_lon, ref_lat)

   fig = plt.figure(figsize=(13,8.5), frameon=False)
   # ax = plt.axes(projection=projection)
   # ax = fig.add_axes([0,0,0.99,1], projection=proj)
   # ax = fig.add_axes([.125, .11, .775, .77], projection=proj)
   ax = fig.add_axes([0,0,1,1], projection=proj)
   ax.set_extent(extent, crs=orto)

   color = np.array([102,42,27])/255
   # ax.stock_img()
   # ax.coastlines(resolution='50m', color='w', linewidth=2.3, zorder=100)
   ax.coastlines(resolution='50m', color=color, linewidth=1.5, zorder=101)
   C = ax.contour(lons, lats, press.values/100, levels=range(804,1204,4),
                                    colors='k', linewidths=.5,
                                    transform=orto)
   ax.clabel(C, inline=True, fontsize=10, fmt='%d')
   C = ax.contour(lons, lats, press.values/100, levels=range(804,1204,16),
                                    colors='k', linewidths=1.5, zorder=102,
                                    transform=orto)
   ax.clabel(C, inline=True, fontsize=10, fmt='%d')
   # print('Plotting frontogenesis')
   # print(np.nanmin(F))
   # print(np.nanmax(F))

   print('Plotting fronts')
   l = 1
   # C = ax.contourf(lons[::l,::l],
   #                 lats[::l,::l],
   #               fronto[::l,::l],
   #                 # vmin=-8, vmax=8,
   #                 np.arange(-6, 6.5, .5),
   #                 cmap=cm_frentes,
   #                 # cmap=plt.cm.bwr,
   #                 extend='both', transform=orto)
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
   msg  = [f"valid: {(fcstDate+UTCshift).strftime('%Y-%m-%d %H:%M')}"]
   msg += [f"GFS: {(dataDate).strftime('%Y-%m-%d %H:%M')}"]
   msg += [f"plot: {dt.datetime.now().strftime('%Y-%m-%d %H:%M')}"]
   ax.text(0, .008, f'{isopress:~}',
           transform = ax.transAxes,
           bbox=dict(facecolor=(1,1,1,.9), edgecolor='k', pad=5.0))
   ax.text(1, .008, '\n'.join(msg),
           transform = ax.transAxes, ha='right',
           bbox=dict(facecolor=(1,1,1,.9), edgecolor='k', pad=5.0))
           # backgroundcolor=(1,1,1,.5), edgecolor='k')
   # ax.text(.008,.008,fcstDate+UTCshift,transform = ax.transAxes,
   #         backgroundcolor=(1,1,1,.5), edgecolor='k')
# cb = plt.colorbar(C, orientation='horizontal', pad=0, aspect=50, extendrect=True)

   # from skimage.feature import peak_local_max
   # aux = press.values/100 - 1013
   # XX = peak_local_max(aux, threshold_rel=.57)
   # labels,X,Y = [],[],[]
   # for p in XX:
   #    i,j = p
   #    print(lats[i,j], lons[i,j], 'H')
   #    ax.text(lons[i,j], lats[i,j], 'H', transform=orto)   

   # plt.show()
   print(f"Saving: {fcstDate.strftime('%Y%m%d_%H%M')}.png")
   fig.savefig(f"{fcstDate.strftime('%Y%m%d_%H%M')}.png")
   plt.close('all')


import generate_config as gen
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

files = os.popen('ls dataGFS/gfs.*').read().strip().split()

from multiprocessing import Pool

inputs = []
for f in files:
   inputs.append((f,left,right,bottom,top))
n = 5
with Pool(n) as pool:
   Res = pool.map(WTFF, inputs)
pool.close()
pool.join()
del pool

# for fname in files:
#    told = time()
#    print(fname)
#    WTFF((fname, left,right,bottom,top))
#    print(time()-told)
#    exit()

com = 'rm *.mp4 *.webm'
os.system(com)
com = 'ls -1 20*.png > files.txt'
os.system(com)
com = 'mencoder -nosound -ovc lavc -lavcopts vcodec=mpeg4 -o fronto.mp4 -mf type=png:fps=2 mf://@files.txt > /dev/null 2> /dev/null'
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
