#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
here = os.path.dirname(os.path.realpath(__file__))
import util as ut

import pathlib
import numpy as np
import matplotlib.pyplot as plt
# try: plt.style.use('rasp')
# except: pass
import datetime as dt
import wrf
from netCDF4 import Dataset
import metpy
import rasterio
from rasterio.merge import merge
from colormaps import WindSpeed, Convergencias, greys, CAPE, Rain

# from random import choice
# files = os.popen('ls wrfout_d01_2021-04-10*').read().strip().splitlines()
# fname = choice(files)

import sys
try: INfname = sys.argv[1]
except IndexError:
   print('File not specified')
   exit()

wrfout_folder, OUT_folder = ut.get_outfolder('../config.ini')
DOMAIN = INfname.split('/')[-1].replace('wrfout_','').split('_')[0]


creation_date = pathlib.Path(INfname).stat().st_mtime
creation_date = dt.datetime.fromtimestamp(creation_date)
# Report here GFS batch and calculation time
gfs_batch = open('batch.txt','r').read().strip()
gfs_batch = dt.datetime.strptime(gfs_batch, '%d/%m/%Y-%H:%M')
date_label =  'GFS: ' + gfs_batch.strftime('%d/%m/%Y-%H:%Mz') + '\n'
date_label += 'plot: ' + creation_date.strftime('%d/%m/%Y-%H:%M ')

# Get UTCshift automatically
UTCshift = dt.datetime.now() - dt.datetime.utcnow()
UTCshift = dt.timedelta(hours = round(UTCshift.total_seconds()/3600))


print('File:',INfname)
ncfile = Dataset(INfname)

# Date in UTC
# prefix to save files
date = str(wrf.getvar(ncfile, 'times').values)
date = dt.datetime.strptime(date[:-3], '%Y-%m-%dT%H:%M:%S.%f')
print('Forecat for:',date)

# Variables for saving outputs
OUT_folder = '/'.join([OUT_folder,DOMAIN,date.strftime('%Y/%m/%d')])
com = f'mkdir -p {OUT_folder}'
os.system(com)
HH = date.strftime('%H%M')



## Calculation parameters
# for x in ncfile.ncattrs():
#    print(x)
#    print(ncfile.getncattr(x))
#    print('')
reflat = ncfile.getncattr('CEN_LAT')
reflon = ncfile.getncattr('CEN_LON')

## all WRF variables
# for v,k in ncfile.variables.items():
#    print(v)
#    # print(k.ncattrs())
#    try: print(k.getncattr('description'))
#    except: print('Description:')
#    try: print(k.getncattr('units'))
#    except: print('units: None')
#    # print(k.dimensions)
#    print(k.shape)
#    print('')
# print('*******')




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
# Topography in metres used in the calculations______________________[m] (ny,nx)
terrain = wrf.getvar(ncfile, "ter", units='m') # = HGT
# Vertical levels of the grid_____________________________________[m] (nz,ny,nx)
# Also called Geopotential Heights. heights[0,:,:] is the first level, 15m
# above ground 
# XXX why 15 and not 10?
heights = wrf.getvar(ncfile, "height", units='m') # = z


# Pressure______________________________________________________[hPa] (nz,ny,nx)
pressure = wrf.getvar(ncfile, "pressure")
print('Pressure:',pressure.shape)
# Perturbation pressure__________________________________________[Pa] (nz,ny,nx)
p = wrf.getvar(ncfile, "P")
print('P:',p.shape)
# Base state pressure____________________________________________[Pa] (nz,ny,nx)
pb = wrf.getvar(ncfile, "PB")
print('PB:', pb.shape)
# Sea Level Pressure________________________________________________[mb] (ny,nx)
slp = wrf.getvar(ncfile, "slp", units='mb')
print('SLP:',slp.shape)


# Sounding properties
# Temperature____________________________________________________[°C] (nz,ny,nx)
tc = wrf.getvar(ncfile, "tc")
print('tc',tc.shape)
# Temperature Dew Point__________________________________________[°C] (nz,ny,nx)
td = wrf.getvar(ncfile, "td", units='degC')
print('td',td.shape)
# Temperature 2m above ground___________________________________[K-->°C] (ny,nx)
t2m = wrf.getvar(ncfile, "T2").metpy.convert_units('degC')
t2m.attrs['units'] = 'degC'
print('t2m',t2m.shape)
# SOIL TEMPERATURE AT LOWER BOUNDARY____________________________[K-->°C] (ny,nx)
tmn = wrf.getvar(ncfile, "TMN").metpy.convert_units('degC')
tmn.attrs['units'] = 'degC'
print('tmn',tmn.shape)
# SOIL TEMPERATURE AT LOWER BOUNDARY____________________________[K-->°C] (ny,nx)
tsk = wrf.getvar(ncfile, "TSK").metpy.convert_units('degC')
tsk.attrs['units'] = 'degC'
print('tsk',tsk.shape)

# lat,lon = 40.87575,-3.68661
# j,i = wrf.ll_to_xy(ncfile, lat, lon)
# print('tc:',tc[0,j,i].values)
# print('t2m:',t2m[j,i].values)
# print('tmn:',tmn[j,i].values)
# print('tsk:',tsk[j,i].values)
# exit()


# Planetary Boundary Layer Height____________________________________[m] (ny,nx)
# Atmospheric Boundary layer thickness above ground
bldepth = wrf.getvar(ncfile, "PBLH")
print('PBLH',bldepth.shape)

# Surface sensible heat flux in___________________________________[W/m²] (ny,nx)
hfx = wrf.getvar(ncfile, "HFX") 
print('HFX:',bldepth.shape)

# Wind__________________________________________________________[m/s] (nz,ny,nx)
ua = wrf.getvar(ncfile, "ua")  # U wind component
va = wrf.getvar(ncfile, "va")  # V wind component
wa = wrf.getvar(ncfile, "wa")  # W wind component
wspd_wdir = wrf.getvar(ncfile, "wspd_wdir")
wspd = np.sqrt(ua*ua + va*va ) # Wind instensity is only 2D
sfcwind = wspd[0,:,:]
print('Wind:',ua.shape)

# Cloud water mixing ratio___________________________________[Kg/kg?] (nz,ny,nx)
qcloud = wrf.getvar(ncfile, "QCLOUD")#"Cloud water mixing ratio"
print('Qcloud:',qcloud.shape)

# Water vapor mixing ratio_________________________________________[] (nz,ny,nx)
qvapor = wrf.getvar(ncfile, "QVAPOR")
print('Qvapor:',qvapor.shape)


# Rain______________________________________________________________[mm] (ny,nx)
rain = wrf.getvar(ncfile, "RAINC") + wrf.getvar(ncfile, "RAINNC")
print('rain:', rain.shape)

# Clouds_____________________________________________________________[%] (ny,nx)
low_cloudfrac  = wrf.getvar(ncfile, "low_cloudfrac")
mid_cloudfrac  = wrf.getvar(ncfile, "mid_cloudfrac")
high_cloudfrac = wrf.getvar(ncfile, "high_cloudfrac")
print('low cloud:', low_cloudfrac.shape)
print('mid cloud:', mid_cloudfrac.shape)
print('high cloud:', high_cloudfrac.shape)
blcloudpct = low_cloudfrac+mid_cloudfrac+high_cloudfrac
blcloudpct = np.clip(blcloudpct, None, 100)


# CAPE____________________________________________________________[J/kg] (ny,nx)
cape2d = wrf.getvar(ncfile, "cape_2d")
MCAPE = cape2d[0,:,:]  # CAPE
MCIN = cape2d[1,:,:]   #CIN
LCL =cape2d[2,:,:]   # Cloud base when forced lifting occurs



## Derived Quantities by DrJack ################################################
# Using utils wrappers to hide the transpose of every variable XXX Inefficient
# BL Max. Up/Down Motion (BL Convergence)_________________________[cm/s] (ny,nx)
wblmaxmin = ut.calc_wblmaxmin(0, wa, heights, terrain, bldepth)
print('WBLmaxmin:',wblmaxmin.shape)

# Thermal Updraft Velocity (W*)____________________________________[m/s] (ny,nx)
wstar = ut.calc_wstar( hfx, bldepth )
print('W*:',wstar.shape)

# BLcwbase___________________________________________________________[m] (ny,nx)
laglcwbase = 0 # lagl = 0 --> height above sea level
               # lagl = 1 --> height above ground level
# criteriondegc = 1.0
maxcwbasem = 5486.40
cwbasecriteria = 0.000010
blcwbase = ut.calc_blcloudbase( qcloud,  heights, terrain, bldepth,
                                cwbasecriteria, maxcwbasem, laglcwbase)
print('BLcwbase:',blcwbase.shape)

# Height of Critical Updraft Strength (hcrit)________________________[m] (ny,nx)
hcrit = ut.calc_hcrit( wstar, terrain, bldepth)
print('Hcrit:',hcrit.shape)

# Height of SFC.LCL__________________________________________________[m] (ny,nx)
# Cu Cloudbase ~I~where Cu Potential > 0~P~
zsfclcl = ut.calc_sfclclheight( pressure, tc, td, heights, terrain, bldepth )
print('zsfclcl:',zsfclcl.shape)
# Mask Cu Pot > 0
zsfclcldif = bldepth + terrain - zsfclcl
null = 0. * zsfclcl
# cu_base_pote = np.where(zsfclcldif>0, zsfclcl, null)
zsfclcl = np.where(zsfclcldif>0, zsfclcl, null)

# OvercastDevelopment Cloudbase__________________________________[m?] (nz,ny,nx)
pmb = 0.01*(p.values+pb.values) # press is vertical coordinate in mb
zblcl = ut.calc_blclheight(qvapor,heights,terrain,bldepth,pmb,tc)
print('zblcl:',zblcl.shape)
# Mask Overcast dev Pot > 0
zblcldif = bldepth + terrain - zblcl
null = 0. * zblcl
# over_base_pote = np.where(zblcldif>0, zblcl, null)
zblcl = np.where(zblcldif>0, zblcl, null)


# Thermalling Height_________________________________________________[m] (ny,nx)
# From Oriol Cervello's
hglider = np.minimum(np.minimum(hcrit,zsfclcl), zblcl)
print('hglider:',hglider.shape)
# XXX warning!! hglider becomes numpy.array
# from xarray.core.dataarray import DataArray
# hglider = DataArray(hglider)
# hglider.attrs["units"] = "metres"


# BL Avg Wind_____________________________________________________[m/s?] (ny,nx)
# uv NOT rotated to grid in m/s
uv = wrf.getvar(ncfile, "uvmet")
uEW = uv[0,:,:,:]
vNS = uv[1,:,:,:]
ublavgwind = ut.calc_blavg(uEW, heights, terrain, bldepth)
vblavgwind = ut.calc_blavg(vNS, heights, terrain, bldepth)
blwind = np.sqrt(ublavgwind*ublavgwind + vblavgwind*vblavgwind)
print('uBLavg:',ublavgwind.shape)
print('vBLavg:',vblavgwind.shape)
print('BLwind:',blwind.shape)

# BL Top Wind_____________________________________________________[m/s?] (ny,nx)
utop,vtop = ut.calc_bltopwind(uEW, vNS, heights,terrain,bldepth)
bltopwind = np.sqrt(utop*utop + vtop*vtop)
print('utop:',utop.shape)
print('vtop:',vtop.shape)
print('BLtopwind:',bltopwind.shape)

















################################################################################
#                                     Plots                                     #
################################################################################

## Soundings ###################################################################
f_cities = f'{here}/soundings.csv'
Yt,Xt = np.loadtxt(f_cities,usecols=(0,1),delimiter=',',unpack=True)
names = np.loadtxt(f_cities,usecols=(2,),delimiter=',',dtype=str)
soundings = [(n,(la,lo))for n,la,lo in zip(names,Yt,Xt)]
for place,point in soundings:
   lat,lon = point
   name = f'{OUT_folder}/{HH}_sounding_{place}.png'
   title = f"{place.capitalize()}"
   title += f" {(date+UTCshift).strftime('%d/%m/%Y-%H:%M')}"
   ut.sounding(lat,lon,lats,lons,date,ncfile,pressure,tc,td,t2m,ua,va,title,fout=name)



## Scalar properties ###########################################################
import plot_functions as PF   # My plotting functions
import matplotlib as mpl
#  COLOR = 'black'
#  ROLOC = '#e0e0e0'
mpl.rcParams['axes.facecolor'] = (1,1,1,0)
mpl.rcParams['figure.facecolor'] = (1,1,1,0)
mpl.rcParams["savefig.facecolor"] = (1,1,1,0)
mpl.rcParams["figure.dpi"] = 150
dpi = 150
from configparser import ConfigParser, ExtendedInterpolation
def get_properties(fname,section):
   """
   Load the config options and return it as a class
   """
   # LG.info(f'Loading config file: {fname}')
   # if not os.path.isfile(fname): return None
   config = ConfigParser(inline_comment_prefixes='#')
   config._interpolation = ExtendedInterpolation()
   config.read(fname)
   # LG.debug('Trying to read start/end dates')
   #XXX We shouldn't use eval
   factor = float(eval(config[section]['factor']))
   vmin   = float(eval(config[section]['vmin']))
   vmax   = float(eval(config[section]['vmax']))
   delta  = float(eval(config[section]['delta']))
   try: 
       levels = config[section]['levels'].replace(']','').replace('[','')
       levels = list(map(float,levels.split(',')))
   except KeyError: levels = []
   levels = [float(l) for l in levels]
   cmap = config[section]['cmap']
   units = config[section]['units']
   return factor,vmin,vmax,delta,levels,cmap,units




fontsize_title = 30
from time import time
# Background plots #############################################################
## Terrain 
fig,ax,orto = PF.terrain_plot(reflat,reflon,left,right,bottom,top)
fname = f'{OUT_folder}/terrain.png'
PF.save_figure(fig,fname,dpi=dpi)

# ## Ocean
# fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
# PF.sea_plot(fig,ax,orto)
# fig.savefig(f'{OUT_folder}/ocean.png', transparent=True, bbox_inches='tight',
#                    pad_inches=0, #dpi=90,
#                    quality=90)
## Parallel and meridian
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.parallel_and_meridian(fig,ax,orto,left,right,bottom,top)
fname = f'{OUT_folder}/meridian.png'
PF.save_figure(fig,fname,dpi=dpi)

## Rivers
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.rivers_plot(fig,ax,orto)
fname = f'{OUT_folder}/rivers.png'
PF.save_figure(fig,fname,dpi=dpi)

## CCAA
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.ccaa_plot(fig,ax,orto)
fname = f'{OUT_folder}/ccaa.png'
PF.save_figure(fig,fname,dpi=dpi)

## Cities
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.csv_plot(fig,ax,orto,f'{here}/cities.csv')
fname = f'{OUT_folder}/cities.png'
PF.save_figure(fig,fname,dpi=dpi)

## Citiy Names
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.csv_names_plot(fig,ax,orto,f'{here}/cities.csv')
fname = f'{OUT_folder}/cities_names.png'
PF.save_figure(fig,fname,dpi=dpi)

## Takeoffs 
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.csv_plot(fig,ax,orto,f'{here}/takeoffs.csv')
fname = f'{OUT_folder}/takeoffs.png'
PF.save_figure(fig,fname,dpi=dpi)

## Takeoffs Names
fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
PF.csv_names_plot(fig,ax,orto,f'{here}/takeoffs.csv')
fname = f'{OUT_folder}/takeoffs_names.png'
PF.save_figure(fig,fname,dpi=dpi)


# Properties ###################################################################
wrf_properties = {'sfcwind':sfcwind, 'blwind':blwind, 'bltopwind':bltopwind,
                  'hglider':hglider, 'wstar':wstar, 'zsfclcl':zsfclcl,
                  'zblcl':zblcl, 'cape':MCAPE, 'wblmaxmin':wblmaxmin,
                  'bldepth': bldepth,  #'bsratio':bsratio,
                  'rain':rain, 'blcloudpct':blcloudpct}

colormaps = {'WindSpeed': WindSpeed, 'Convergencias':Convergencias,
             'greys':greys, 'CAPE': CAPE, 'Rain': Rain, 'None':None}

titles = {'sfcwind':'Viento Superficie', 'blwind':'Viento Promedio',
          'bltopwind':'Viento Altura', 'hglider':'Techo (azul)',
          'wstar':'Térmica', 'zsfclcl':'Base nube', 'zblcl':'Cielo cubierto',
          'cape':'CAPE', 'wblmaxmin':'Convergencias',
          'bldepth': 'Altura Capa Convectiva', 'bsratio': 'B/S ratio',
          'rain': 'Lluvia', 'blcloudpct':'Nubosidad (%)'}

# prop_units = {'sfcwind':'km/h', 'blwind':'km/h',
#               'bltopwind':'km/h', 'hglider':'m',
#               'wstar':'m/s', 'zsfclcl':'m', 'zblcl':'m',
#               'cape':'J/Kg', 'wblmaxmin':'m/s',
#               'bldepth': 'm', #'bsratio': 'B/S ratio',
#               'rain': 'mm', 'blcloudpct':'%'}

# plot scalars #################################################################
ftitles = open(f'{OUT_folder}/titles.txt','w')
for prop in ['sfcwind', 'blwind', 'bltopwind', 'hglider', 'wstar', 'zsfclcl',
             'zblcl', 'cape', 'wblmaxmin', 'bldepth', # 'bsratio',
             'rain', 'blcloudpct']:
   factor,vmin,vmax,delta,levels,cmap,units = get_properties('plots.ini', prop)
   cmap = colormaps[cmap]   #Convergencias
   title = titles[prop]
   title = f"{title} {(date+UTCshift).strftime('%d/%m/%Y-%H:%M')}"
   M = wrf_properties[prop]
   fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
   C = PF.scalar_plot(fig,ax,orto, lons,lats,wrf_properties[prop]*factor,
                      delta,vmin,vmax,cmap,
                      levels=levels,
                      creation_date=date_label)
   fname = f'{OUT_folder}/{HH}_{prop}.png'
   PF.save_figure(fig,fname,dpi=dpi)

   ftitles.write(f"{fname} ; {title}\n")
   PF.plot_colorbar(cmap,delta,vmin,vmax,levels,name=f'{OUT_folder}/{prop}',
                    units=units,fs=15,norm=None,extend='max')
   plt.close('all')
ftitles.close()


## Vector properties ###########################################################
names = ['sfcwind','blwind','bltopwind']
winds = [[ua[0,:,:].values, va[0,:,:].values],
         [ublavgwind, vblavgwind],
         [utop, vtop]]

for wind,name in zip(winds,names):
   fig,ax,orto = PF.setup_plot(reflat,reflon,left,right,bottom,top)
   U = wind[0]
   V = wind[1]
   PF.vector_plot(fig,ax,orto,lons,lats,U,V, dens=1.5,color=(0,0,0))
   # fname = OUT_folder +'/'+ prefix + name + '_vec.png'
   fname = f'{OUT_folder}/{HH}_{name}_vec.png'
   PF.save_figure(fig,fname,dpi=dpi)

#XXX shouldn't do this here
wrfout_folder += gfs_batch.strftime('/%Y/%m/%d/%H')
com = f'mkdir -p {wrfout_folder}'
print('****')
print(com)
os.system(com)
com = f'mv {INfname} {wrfout_folder}'
print('****')
print(com)
os.system(com)
