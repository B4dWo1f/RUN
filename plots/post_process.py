#!/usr/bin/python3
# -*- coding: UTF-8 -*-

# Recompile DrJack's Fortran if necessary
import os
diff = os.popen('diff drjack.f90 .drjack.f90').read().strip()
if len(diff) > 0:
   os.system('f2py3 -c -m drjack drjack.f90 && cp drjack.f90 .drjack.f90')
else: print('Already compiled')
import drjack
import plot_functions as PF   # My plotting functions

import numpy as np
import matplotlib.pyplot as plt
try: plt.style.use('rasp')
except: pass
import datetime as dt
import wrf
from netCDF4 import Dataset
import metpy
import rasterio
from rasterio.merge import merge
from colormaps import WindSpeed, Convergencias, greys

# from random import choice
# files = os.popen('ls wrfout_d01_2021-04-10*').read().strip().splitlines()
# fname = choice(files)

import sys
try: fname = sys.argv[1]
except IndexError:
   print('File not specified')
   exit()



# Get UTCshift automatically
UTCshift = dt.datetime.now() - dt.datetime.utcnow()
UTCshift = dt.timedelta(hours = round(UTCshift.total_seconds()/3600))


print('File:',fname)
ncfile = Dataset(fname)
## Calculation parameters
# for x in ncfile.ncattrs():
#    print(x)
#    print(ncfile.getncattr(x))
#    print('')
reflat = ncfile.getncattr('CEN_LAT')
reflon = ncfile.getncattr('CEN_LON')

## all WRF variables
for v,k in ncfile.variables.items():
   print(v)
   # print(k.ncattrs())
   try: print(k.getncattr('description'))
   except: print('Description:')
   try: print(k.getncattr('units'))
   except: print('units: None')
   # print(k.dimensions)
   print(k.shape)
   print('')
print('*******')

# prefix to save files
prefix = fname.split('/')[-1].replace('wrfout_','').split(':')[0] + '_'

# Date in UTC
date = str(wrf.getvar(ncfile, 'times').values)
date = dt.datetime.strptime(date[:-3], '%Y-%m-%dT%H:%M:%S.%f')
print('Forecat for:',date)


### Read region data
# The domain is a rectangular (regular) grid in Lambert projection
# Latitude, longitude  [degrees]
lats = wrf.getvar(ncfile, "lat")
lons = wrf.getvar(ncfile, "lon")
# bounds contain the bottom-left and upper-right corners of the domain
# Notice that bounds will not be the left/right/top/bottom-most
# latitudes/longitudes since the grid is only regular in Lambert
bounds = wrf.geo_bounds(wrfin=ncfile)
print('lat/lon',lats.shape)

# useful to setup the extent of the maps
left   = bounds.bottom_left.lon
right  = bounds.top_right.lon
bottom = bounds.bottom_left.lat
top    = bounds.top_right.lat
# left   = np.min(wrf.to_np(lons))
# right  = np.max(wrf.to_np(lons))
# bottom = np.min(wrf.to_np(lats))
# top    = np.max(wrf.to_np(lats))

# WRF-Terrain
# Topography in metres used in the calculations_______________________[m] (ny,nx)
terrain = wrf.getvar(ncfile, "ter", units='m') # = HGT
# Vertical levels of the grid______________________________________[m] (nz,ny,nx)
# Also called Geopotential Heights. heights[0,:,:] is the first level, 15m
# above ground 
# XXX why 15 and not 10?
heights = wrf.getvar(ncfile, "height", units='m') # = z


# Pressure_______________________________________________________[hPa] (nz,ny,nx)
pressure = wrf.getvar(ncfile, "pressure")
print('Pressure:',pressure.shape)
# Perturbation pressure___________________________________________[Pa] (nz,ny,nx)
p = wrf.getvar(ncfile, "P")
print('P:',p.shape)
# Base state pressure_____________________________________________[Pa] (nz,ny,nx)
pb = wrf.getvar(ncfile, "PB")
print('PB:', pb.shape)
# Sea Level Pressure_________________________________________________[mb] (ny,nx)
slp = wrf.getvar(ncfile, "slp", units='mb')
print('SLP:',slp.shape)


# Sounding properties
# Temperature_____________________________________________________[°C] (nz,ny,nx)
tc = wrf.getvar(ncfile, "tc")
print('tc',tc.shape)
# Temperature Dew Point___________________________________________[°C] (nz,ny,nx)
td = wrf.getvar(ncfile, "td", units='degC')
print('td',td.shape)
# Temperature 2m above ground____________________________________[K-->°C] (ny,nx)
t2m = wrf.getvar(ncfile, "T2").metpy.convert_units('degC')
t2m.attrs['units'] = 'degC'
# SOIL TEMPERATURE AT LOWER BOUNDARY_____________________________[K-->°C] (ny,nx)
tmn = wrf.getvar(ncfile, "TMN").metpy.convert_units('degC')
tmn.attrs['units'] = 'degC'
# SOIL TEMPERATURE AT LOWER BOUNDARY_____________________________[K-->°C] (ny,nx)
tsk = wrf.getvar(ncfile, "TSK").metpy.convert_units('degC')
tsk.attrs['units'] = 'degC'

print(t2m)
print(tsk)
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
ax.contourf(lons,lats,tsk-t2m, vmin=-5,vmax=5,transform=orto,alpha=0.5)
fig.tight_layout()
plt.show()
exit()

# Planetary Boundary Layer Height_____________________________________[m] (ny,nx)
# Atmospheric Boundary layer thickness above ground
bldepth = wrf.getvar(ncfile, "PBLH")
print('PBLH',bldepth.shape)

# Surface sensible heat flux in____________________________________[W/m²] (ny,nx)
hfx = wrf.getvar(ncfile, "HFX") 
print('HFX:',bldepth.shape)

# Wind___________________________________________________________[m/s] (nz,ny,nx)
ua = wrf.getvar(ncfile, "ua")  # U wind component
va = wrf.getvar(ncfile, "va")  # V wind component
wa = wrf.getvar(ncfile, "wa")  # W wind component
wspd_wdir = wrf.getvar(ncfile, "wspd_wdir")
wspd = np.sqrt(ua*ua + va*va ) # Wind instensity is only 2D
print('Wind:',ua.shape)

# Cloud water mixing ratio____________________________________[Kg/kg?] (nz,ny,nx)
qcloud = wrf.getvar(ncfile, "QCLOUD")#"Cloud water mixing ratio"
print('Qcloud:',qcloud.shape)

# Water vapor mixing ratio__________________________________________[] (nz,ny,nx)
qvapor = wrf.getvar(ncfile, "QVAPOR")
print('Qvapor:',qvapor.shape)


# Rain
rain = wrf.getvar(ncfile, "RAINC") + wrf.getvar(ncfile, "RAINNC")
print('rain:', rain.shape)

# Clouds
low_cloudfrac  = wrf.getvar(ncfile, "low_cloudfrac")
mid_cloudfrac  = wrf.getvar(ncfile, "mid_cloudfrac")
high_cloudfrac = wrf.getvar(ncfile, "high_cloudfrac")
print('low cloud:', low_cloudfrac.shape)
print('mid cloud:', mid_cloudfrac.shape)
print('high cloud:', high_cloudfrac.shape)

# Derived Quantities by DrJack
# BL Max. Up/Down Motion (BL Convergence)__________________________[cm/s] (ny,nx)
wblmaxmin = drjack.calc_wblmaxmin(0, wa.transpose(),
                                     heights.transpose(),
                                     terrain.transpose(),
                                     bldepth.transpose()).transpose()
print('WBLmaxmin:',wblmaxmin.shape)

# Thermal Updraft Velocity (W*)_____________________________________[m/s] (ny,nx)
wstar = drjack.calc_wstar( hfx, bldepth )
print('W*:',wstar.shape)

# BLcwbase____________________________________________________________[m] (ny,nx)
laglcwbase = 0 # lagl = 0 --> height above sea level
               # lagl = 1 --> height above ground level
# criteriondegc = 1.0
maxcwbasem = 5486.40
cwbasecriteria = 0.000010
# all the transposes are necessary to agree with DrJack's notation
blcwbase = drjack.calc_blcloudbase( qcloud.transpose(),  heights.transpose(),
                                    terrain.transpose(), bldepth.transpose(),
                                    cwbasecriteria, maxcwbasem, laglcwbase)
blcwbase = blcwbase.transpose()
print('BLcwbase:',blcwbase.shape)

# Height of Critical Updraft Strength (hcrit)_________________________[m] (ny,nx)
hcrit = drjack.calc_hcrit( wstar, terrain, bldepth)
print('Hcrit:',hcrit.shape)

# Height of SFC.LCL___________________________________________________[m] (ny,nx)
# Cu Cloudbase ~I~where Cu Potential > 0~P~
zsfclcl = drjack.calc_sfclclheight( pressure.transpose(),
                                    tc.transpose(), td.transpose(),
                                    heights.transpose(),
                                    terrain.transpose(),
                                    bldepth.transpose() )
zsfclcl = zsfclcl.transpose()
print('zsfclcl:',zsfclcl.shape)

# OvercastDevelopment Cloudbase___________________________________[m?] (nz,ny,nx)
qvaporblavg = drjack.calc_blavg( qvapor.transpose(),
                                 heights.transpose(),
                                 terrain.transpose(),
                                 bldepth.transpose() )
qvaporblavg = qvaporblavg.transpose()
pmb=var = 0.01*(p.values+pb.values) # press is vertical coordinate in mb
zblcl = drjack.calc_blclheight( pmb.transpose(),
                                tc.transpose(),
                                qvaporblavg.transpose(),
                                heights.transpose(),
                                terrain.transpose(),
                                bldepth.transpose() )
zblcl = zblcl.transpose()
print('zblcl:',zblcl.shape)


# Thermalling Height__________________________________________________[m] (ny,nx)
# Oriol
# hglider_aux = np.minimum(hwcrit,zsfclcl)
# hglider=np.minimum(hglider_aux,zblcl)
hglider = np.minimum(np.minimum(hcrit,zsfclcl), zblcl)
print('hglider:',hglider.shape)

# from xarray.core.dataarray import DataArray
# hglider = DataArray(hglider)
# hglider.attrs["units"] = "metres"


# BL Avg Wind______________________________________________________[m/s?] (ny,nx)
# uv NOT rotated to grid in m/s
uv = wrf.getvar(ncfile, "uvmet")
uEW = uv[0,:,:,:]
vNS = uv[1,:,:,:]
ublavgwind = drjack.calc_blavg(uEW.transpose(),
                               heights.transpose(),
                               terrain.transpose(),
                               bldepth.transpose())
vblavgwind = drjack.calc_blavg(vNS.transpose(),
                               heights.transpose(),
                               terrain.transpose(),
                               bldepth.transpose())
ublavgwind = ublavgwind.transpose()
vblavgwind = vblavgwind.transpose()
print('uBLavg:',ublavgwind.shape)
print('vBLavg:',vblavgwind.shape)

# BL Top Wind______________________________________________________[m/s?] (ny,nx)
utop,vtop = drjack.calc_bltopwind(uEW.transpose(),
                                  vNS.transpose(),
                                  heights.transpose(),
                                  terrain.transpose(),
                                  bldepth.transpose())
utop = utop.transpose()
vtop = vtop.transpose()
print('utop:',utop.shape)
print('vtop:',vtop.shape)








################   Skew T
lat,lon = 41.078854,-3.707029
def sounding(lat,lon,date,ncfile,tc,td,t0,ua,uv,fout='sounding.png'):
               # ,UTCshift=dt.timedelta(days=0),fout='sounding.png'):
   """
   lat,lon: spatial coordinates for the sounding
   date: UTC date-time for the sounding
   ncfile: ntcd4 Dataset from the WRF output
   tc: Model temperature in celsius
   tdc: Model dew temperature in celsius
   t0: Model temperature 2m above ground
   ua: Model X wind (m/s)
   va: Model Y wind (m/s)
   fout: save fig name
   """
   j,i = wrf.ll_to_xy(ncfile, lat, lon)
   # Get sounding data for specific location
   h = heights[:,i,j]
   p = pressure[:,i,j]
   tc = tc[:,i,j]
   tdc = td[:,i,j]
   t0 = t0[i,j]
   u = ua[:,i,j]
   v = va[:,i,j]

   fig,ax = PF.skewt_plot(p,tc,tdc,t0,date,u,v,show=False)
   fig.savefig(fout)


def strip_plot(fig,ax,lims,aspect,fname,dpi=65):
   ax.set_aspect(aspect)
   ax.set_xlim(lims[0:2])
   ax.set_ylim(lims[2:4])
   ax.get_xaxis().set_visible(False)
   ax.get_yaxis().set_visible(False)
   plt.axis('off')
   fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
   fig.savefig(fname, transparent=True, bbox_inches='tight', pad_inches=0,
                      dpi=dpi, quality=90)   # compression

soundings = [('bustarviejo', (40.87575, -3.68661)),
             ('abantos', (40.611774, -4.154882)),
             ('canencia', (40.93716440325912, -3.7299387579645904)),
             ('arcones', (41.078854, -3.707029))]
for name,point in soundings:
   lat,lon = point
   name = prefix + 'sounding_' + name + '.png'
   sounding(lat,lon,date,ncfile,tc,td,t2m,ua,uv,fout=name)


def cross_path(start,end):
   start = wrf.ll_to_xy(ncfile, start[0], start[1])
   end = wrf.ll_to_xy(ncfile, end[0], end[1])

   start_point = wrf.CoordPair(x=start[0], y=start[1])
   end_point   = wrf.CoordPair(x=end[0], y=end[1])

   levels = np.linspace(900,3000,100)
   p_vert = wrf.vertcross(wspd, heights, start_point=start_point,
                                         end_point=end_point,
                                         levels=levels,
                                         latlon=True)
   fig, ax = plt.subplots()
   dx = 0,p_vert.shape[1]
   ax.imshow(p_vert,origin='lower',extent=[*dx,900,3000])
   ax.set_aspect(1/100)
   ax.set_title('Wind speed')
   fig.tight_layout()
   plt.show()


# start = 41.260096661378306, -3.6981749563104716
# end = 40.78691893349439, -3.6903966736445306
# cross_path(start,end)

################








#################################################################################
#                                     Plots                                     #
#################################################################################
from time import time

# BLdepth #######################################################################
told = time()
vmin = None
vmax = None
delta= None
cmap = 'viridis'  #Convergencias
# levels = [-3,-2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3]
title = 'BLdepth ' + date.strftime('%d/%m/%Y-%H:%M')
told1 = time()
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
print('Base plot:',time()-told1)
ax.set_title(title)
PF.scalar_plot(fig,ax,orto,lons,lats,bldepth,delta,vmin,vmax,cmap)
# extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
# fig.savefig(prefix + f'bldepth.png', bbox_inches = extent)
fname = prefix + f'bldepth.png'

fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
fig.savefig(fname, transparent=False, bbox_inches='tight', pad_inches=0,
                   dpi=90, quality=90)
print('bldepth plot:',time()-told)





# extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
# fig.savefig(prefix + f'bldepth.png', bbox_inches = extent)
# plt.show()
# exit()




# Wind ##########################################################################
told = time()
vmin = 0
vmax = 60
delta = 4
cmap = WindSpeed
title = 'Wind ' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M')
# scalar
spd = wspd[0,:,:].values * 3.6  #km/h
clabel = 'km/h'
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
PF.scalar_plot(fig,ax,orto,lons,lats,spd,delta,vmin,vmax,cmap)
U = ua[0,:,:].values
V = va[0,:,:].values
# vector
PF.vector_plot(fig,ax,orto,lons,lats,U,V,dens=2)
# n=10
# ax.barbs(lons[::n,::n].values, lats[::n,::n].values,
#           U[::n, ::n], V[::n, ::n], color='r',
#           length=6, pivot='middle',
#           transform=orto)
ax.set_title(title)
fig.savefig(prefix + 'sfcwind.png')
print('sfcwind plot:',time()-told)


# WBLmaxmin #####################################################################
told = time()
vmin = -3
vmax = 3
delta= 0.1
cmap = Convergencias
levels = [-3,-2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3]
title = 'Convergencia ' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M')
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
PF.scalar_plot(fig,ax,orto,lons,lats,wblmaxmin/100,delta,vmin,vmax,cmap,levels=levels)
ax.set_title(title)
fig.savefig(prefix + 'wblmaxmin.png')
print('WBLmaxmin plot:',time()-told)


# Clouds ########################################################################
told = time()
# vmin = -3
# vmax = 3
# delta= 0.1
cmap = greys
# levels = [-3,-2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3]
# title = 'Convergencia ' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M')
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
C = ax.contourf(lons,lats,low_cloudfrac+mid_cloudfrac+high_cloudfrac,
                          cmap=cmap, transform=orto,zorder=100)
ax.set_title('Cloudfrac' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M'))
fig.colorbar(C)
fig.tight_layout()
fig.savefig(prefix + 'cloudfrac.png')


# Clouds ########################################################################
told = time()
vmin = 200
vmax = 3800
cmap = WindSpeed
delta = 240
# levels = [-3,-2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3]
# title = 'Convergencia ' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M')
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
PF.scalar_plot(fig,ax,orto,lons,lats,zblcl,delta,vmin,vmax,cmap)
ax.set_title('zblcl' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M'))
fig.savefig(prefix + 'overcast.png')



# Clouds ########################################################################
told = time()
vmin = 200
vmax = 3800
cmap = WindSpeed
delta = 240
# levels = [-3,-2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3]
# title = 'Convergencia ' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M')
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
PF.scalar_plot(fig,ax,orto,lons,lats,zsfclcl,delta,vmin,vmax,cmap)
ax.set_title('zsfclcl' + (date+UTCshift).strftime('%d/%m/%Y-%H:%M'))
fig.savefig(prefix + 'overcast2.png')

exit()










# print(np.min(low_cloudfrac).values, np.max(low_cloudfrac).values)
# print(np.min(mid_cloudfrac).values, np.max(mid_cloudfrac).values)
# print(np.min(high_cloudfrac).values, np.max(high_cloudfrac).values)
# ax.contourf(lats,lons,low_cloudfrac,  c='C1', transform=orto)
# ax.contourf(lats,lons,mid_cloudfrac,  c='C2', transform=orto)
# ax.contourf(lats,lons,high_cloudfrac, c='C3', transform=orto)
   # C = ax.contourf(lons,lats,prop, levels=levels, extend='max',
   #                                 antialiased=True, norm=norm,
   #                                 cmap=cmap, vmin=vmin, vmax=vmax,
   #                                 alpha=0.5, transform=orto)
# PF.scalar_plot(fig,ax,orto,lons,lats,wblmaxmin/100,delta,vmin,vmax,cmap,levels=levels)
# ax.set_title(title)
fig.savefig(prefix + 'clouds.png')
# print('WBLmaxmin plot:',time()-told)

exit()



################################### OBSOLETE ###################################
fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
ax.scatter([lons[0,0],lons[0,0],lons[-1,-1],lons[-1,-1]],
           [lats[0,0],lats[-1,-1],lats[0,0],lats[-1,-1]],
           s=90,c='g', label='lat/lon', transform=orto)
ax.scatter([left,left,right,right],
           [bottom,top,bottom,top],s=30,c='r', label='lat/lon max/min',
           transform=orto)
left   = bounds.bottom_left.lon
right  = bounds.top_right.lon
bottom = bounds.bottom_left.lat
top    = bounds.top_right.lat
ax.scatter([left,left,right,right],
           [bottom,top,bottom,top],s=20,c='b', label='bound', transform=orto)

ax.contourf(lons,lats,terrain, alpha=0.5, transform=orto)
plt.show()

exit()










# Wind
slp = wrf.getvar(ncfile, "slp", units='mb')   # Sea Level Pressure
pressure = wrf.getvar(ncfile, "pressure")     # Pressure (hPa)
tc = wrf.getvar(ncfile, "tc")                 # temperature (Celsius)
ua = wrf.getvar(ncfile, "ua",units='km h-1')  # U wind component
va = wrf.getvar(ncfile, "va",units='km h-1')  # V wind component
wa = wrf.getvar(ncfile, "wa",units='km h-1')  # W wind component
#alternative
# wspd,wdir = wrf.getvar(ncfile, "wspd_wdir",units='km h-1')  # W wind component
wspd = np.sqrt(ua*ua + va*va ) # + wa*wa)

# Clouds
low_cloudfrac  = wrf.getvar(ncfile, "low_cloudfrac")
mid_cloudfrac  = wrf.getvar(ncfile, "mid_cloudfrac")
high_cloudfrac = wrf.getvar(ncfile, "high_cloudfrac")

# CAPE
cape2d = wrf.getvar(ncfile, "cape_2d")
MCAPE = cape2d[0,:,:]  # CAPE
MCIN = cape2d[1,:,:]   #CIN
LCL =cape2d[2,:,:]   # Cloud base when forced lifting occurs

LFC =cape2d[3,:,:]   # Level of Free Convection. Altitude in the atmosphere
                     # where the temperature of the environment decreases faster
                     # than the moist adiabatic lapse rate of a saturated air
                     # parcel at the same level. 
cape3d = wrf.getvar(ncfile, "cape_3d")
CAPE3d = cape3d[0,:,:,:]
CIN3d = cape3d[1,:,:,:]

print(MCAPE.shape)
print(CAPE3d.shape)


fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
# ax.contourf(lons,lats,MCAPE, alpha=0.5, transform=orto)
# ax.contourf(lons,lats,wspd[0,:,:], vmin=0,vmax=60, extend='max',
#                                   cmap = WindSpeed,
#                                   alpha=0.5, zorder=10,
#                                   transform=orto)
ax.scatter(lons,lats,c=wspd[0,:,:], vmin=0,vmax=60, #extend='max',
                                  cmap = WindSpeed,
                                  alpha=0.5, zorder=10,
                                  transform=orto)
ax.set_title('CAPE')




plt.show()
exit()

fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
ax.contourf(lons,lats,CAPE3d[0,:,:], alpha=0.5, transform=orto)
ax.set_title('CAPE3d (0)')




fig, ax, orto = PF.base_plot(reflat,reflon,left,right,bottom,top)
ax.contourf(lons,lats,wspd[0,:,:], vmin=0,vmax=60, extend='max',
                                  cmap = WindSpeed,
                                  alpha=0.5, zorder=10,
                                  transform=orto)
ax.set_title('Wind')

plt.show()
exit()
x = lons[0,:].values
y = lats[:,0].values
U = ua[0,:,:].values
V = va[0,:,:].values
dens = 1.5
ax.streamplot(x,y, U,V, color=(0,0,0,0.75), linewidth=1, density=dens,
                        arrowstyle='->',arrowsize=2.5,
                        zorder=11,
                        transform=orto)

n=10
ax.barbs(lons[::n,::n].values, lats[::n,::n].values,
          U[::n, ::n], V[::n, ::n], color='r',
          length=6, pivot='middle',
          transform=orto)

plt.show()




exit()

from matplotlib import gridspec
fig = plt.figure()
gs = gridspec.GridSpec(2, 1)
fig.subplots_adjust() #wspace=0.1,hspace=0.1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax1.scatter(lons,lats,c=heights[:,:,0],vmin=0,vmax=2300)
ax2.scatter(lons,lats,c=terrain,vmin=0,vmax=2300)
ax1.scatter(-5.798180335840546, 40.19957065875734,s=100,c='r')
ax2.scatter(-5.798180335840546, 40.19957065875734,s=100,c='r')
ax1.set_aspect('equal')
ax2.set_aspect('equal')
fig.tight_layout()
plt.show()


# exit()

slp = wrf.getvar(ncfile, "slp", units='mb')   # Sea Level Pressure
z = wrf.getvar(ncfile, "z", units='m')        # Geopotential Height
pressure = wrf.getvar(ncfile, "pressure")     # Pressure (hPa)
tc = wrf.getvar(ncfile, "tc")                 # temperature (Celsius)
ua = wrf.getvar(ncfile, "ua",units='km h-1')  # U wind component
va = wrf.getvar(ncfile, "va",units='km h-1')  # V wind component
wa = wrf.getvar(ncfile, "wa",units='km h-1')  # W wind component
wspd_wdir = wrf.getvar(ncfile, "wspd_wdir", units='km h-1')
print(ua.shape)
print(wspd_wdir.shape)
WW = np.sqrt(ua*ua+va*va)
print(np.max(np.abs(wspd_wdir[0,0,:,:]-WW)))
from matplotlib import gridspec
fig = plt.figure()
gs = gridspec.GridSpec(2, 1)
fig.subplots_adjust() #wspace=0.1,hspace=0.1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
ax1.imshow(WW[0,:,:])
ax2.imshow(wspd_wdir[0,0,:,:])
fig.tight_layout()
plt.show()

# exit()
uv10        = getvar(ncfile, "uvmet10")
wspd_wdir10 = getvar(ncfile, "wspd_wdir10")
cape_2d = getvar(ncfile, "cape_2d")
cloudfrac = getvar(ncfile, "cloudfrac") # clouds
ctt = getvar(ncfile, "ctt") # Cloud Top Temperature
pw = getvar(ncfile, "pw")# Precipitable water
bldepth = getvar(ncfile, "PBLH")#bl depth (metres)
ter = getvar(ncfile, "HGT")#Terrain Height (metres)
qcloud = getvar(ncfile, "QCLOUD")#"Cloud water mixing ratio"
t2 = getvar(ncfile, "T2")# TEMP at 2 M
uv = getvar(ncfile, "uvmet") # u,v NOT rotated to grid   in m/s
p = getvar(ncfile, "P") #"perturbation pressure"
pb = getvar(ncfile, "PB") #"BASE STATE PRESSURE"
hfx = getvar(ncfile, "HFX") #for sfc. sensible heat flux in w/m2        142
thetac = getvar(ncfile, "T")+ 26.85 #"perturbation potential tempera    ture theta-t0"
qvapor = getvar(ncfile, "QVAPOR") #"Water vapor mixing ratio"
vhf = getvar(ncfile, "LH") #"LATENT HEAT FLUX AT THE SURFACE"
td = getvar(ncfile, "td") #dew point temperature (C)
td2 = getvar(ncfile, "td2") #dew point temperature at 2 M (C)
snow = getvar(ncfile, "SNOW")#"SNOW WATER EQUIVALENT"
rainc = getvar(ncfile, "RAINC")#"ACCUMULATED TOTAL CUMULUS PRECIPIT    ATION"
rainnc = getvar(ncfile, "RAINNC")#"ACCUMULATED TOTAL GRID SCALE PRE    CIPITATION"
snowh = getvar(ncfile, "SNOWH")#"PHYSICAL SNOW DEPTH"

pblh = bldepth
#lats, lons = latlon_coords(slp)
lats = getvar(ncfile, "XLAT")
lons = getvar(ncfile, "XLONG")
dx = slp.projection.dx
dy = slp.projection.dy



U = wrf.to_np(ua[0,:,:]) * 3.6  #km/h
V = wrf.to_np(va[0,:,:]) * 3.6  #km/h
windspd = np.sqrt(U*U + V*V)

left   = np.min(wrf.to_np(lons))
right  = np.max(wrf.to_np(lons))
bottom = np.min(wrf.to_np(lats))
top    = np.max(wrf.to_np(lats))



def plot_scalar(lats,lons,prop,title=''):
   fig, ax = plt.subplots()
   ax.contourf(lons,lats,prop)
   if len(title)>0:
      ax.set_title(title)
   fig.tight_layout()
   plt.show()




# BL Max Up/Down motion
# wa   = wa.transpose()
# z    = z.transpose()
# ter  = ter.transpose()
# pblh = pblh.transpose()
wblMxMn = drjack.calc_wblmaxmin(0,wa,z,ter,pblh)


laglcwbase = 1
# criteriondegc = 1.0
maxcwbasem = 5486.40
cwbasecriteria = 0.000010

blcwbase = drjack.calc_blcloudbase( qcloud, z, ter, pblh, cwbasecriteria, maxcwbasem, laglcwbase)
print(blcwbase)



zsfclcl = drjack.calc_sfclclheight( p,t,td,z,ter,pblh )
print(zsfclcl)


uEW = uv[:,:,:,0]
vNS = uv[:,:,:,1]
ublavgwind = drjack.calc_blavg(uEW,z,ter,pblh)
vblavgwind = drjack.calc_blavg(vNS,z,ter,pblh)
print(ublavgwind)
print(vblavgwind)




wstar = drjack.calc_wstar( hfx, pblh )
print(wstar)
# plot_scalar(lats,lons,wstar,title='wstar')



# XXX NOT DEFINED????
# wstar = drjack.calc_blmax( a,z,ter,pblh, bparam )
# print(wstar)
# plot_scalar(lats,lons,wstar,title='blmax')


hcrit = drjack.calc_hcrit( wstar, ter, pblh)
print(hcrit)
# plot_scalar(lats,lons,wstar,title='hcrit')


utop,vtop = drjack.calc_bltopwind(uEW,vNS,z,ter,pblh  )
# plot_scalar(lats,lons,np.sqrt(utop*utop+vtop*vtop),title='Wind top')



pmb=var = 0.01*(p.values+pb.values) # press is vertical coordinate in mb

blcldpct = drjack.calc_subgrid_blcloudpct_grads( qvapor, qcloud, t,pmb, z, ter, pblh, cwbasecriteria  )
plot_scalar(lats,lons,blcldpct,title='blcldpct')


exit()











########################################
exit()
print(np.min(wblMxMn), np.max(wblMxMn))

fig, ax = plt.subplots()
A = ax.contourf(lons.transpose(),lats.transpose(),wblMxMn/100,
            extend='max',
            levels=[-3,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,3],
            vmin=-3,vmax=3,
            cmap=Convergencias)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.new_vertical(size="2.95%", pad=0.25, pack_start=True)
fig.add_axes(cax)
cbar = fig.colorbar(A, cax=cax, orientation="horizontal")
# cbar.ax.set_xlabel(units,fontsize=fs)

fig.tight_layout()
plt.show()
