#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from os.path import expanduser
import datetime as dt
import logging
try: LG
except NameError: LG = logging.getLogger(__name__)

def main(folder, domain, start_date, end_date, timedelta=1,
         folder_out='/storage', lglv='debug',
         left=-17, right=8, bottom=30, top=48,
         Ncores=14, wait4batch=40):
   """
   folder: [str] folder keeping the RUN codes (could be a copy)
   domain: [str] WRF Domain name to run. It should be one of the domains
                 in RUN/Domains.
   start_date: [datetime.datetime] Start time of the calculation
   end_date: [datetime.datetime] End time of the calculation
   left, right, bottom, top: [float] Limits of the bounding box of the WRF
     calculation
     - left  : smallest longitude
     - right : biggest longitude
     - bottom: smallest latitude
     - top   : biggest latitude
   Ncores: [int] Number of cores to run WRF
   wait4batch: [float] Minutes to keep trying for the lastest GFS batch
   """
   LG.info('Parameters recived:')
   LG.info(f'RUN folder: {folder}')
   LG.info(f'Domain: {domain}')
   LG.info(f'Dates: {start_date} - {end_date}')
   LG.info(f'Time delta: {timedelta}')
   LG.info(f'OUT folder: {folder_out}')
   fmt = '%d/%m/%Y-%H:%M'
   # daily_hours = [start_date.hour,end_date.hour]
   # Folders
   # folder = '~/METEO'  # now its an argument
   run_folder = f'{folder}/RUN'
   # domain = 'Spain6_1'  # now its an argument
   domains = 1,2   #XXX ignored since 18/6/2021 if it works, it should be removed
   gfs_folder = f'{folder}/dataGFS'

   config = configparser.ConfigParser()
   config['run'] = {}
   config['run']['start_date']  = start_date.strftime(fmt)
   config['run']['end_date']    = end_date.strftime(fmt)
   # config['run']['daily_hours'] = ','.join([str(x) for x in daily_hours])
   # Domain folder
   domain_folder = f'{run_folder}/Domains/{domain}'
   if not os.path.isdir(domain_folder):
      print(f'DOMAIN does not exist: {domain_folder}')
      # LG.warning(f'Folder {domain_folder} does not exsit. Creating it')
      # os.system(f'mkdir -p {domain_folder}')
   config['run']['domain_folder'] = domain_folder
   config['run']['domains'] = ','.join([str(x) for x in domains])
   config['run']['GFS_data_folder'] = gfs_folder
   config['run']['GFS_timedelta'] = str(timedelta)
   config['run']['output_folder'] = f'{folder_out}/WRFOUT/{domain}'  #XXX
   config['run']['plots_folder'] =  f'{folder_out}/PLOTS/{domain}'   #XXX
   config['run']['leftlon']   = str(left)
   config['run']['rightlon']  = str(right)
   config['run']['toplat']    = str(top)
   config['run']['bottomlat'] = str(bottom)
   config['run']['Ncores']    = str(Ncores)
   config['run']['wait4batch']    = str(wait4batch)
   config['run']['log_level'] = lglv
   with open('config.ini', 'w') as configfile:
       config.write(configfile)


if __name__ == '__main__':
   is_cron = False
   is_cron = bool( os.getenv('RUN_BY_CRON') )
################################# LOGGING ####################################
   import log_help
   log_file = '.'.join( __file__.split('/')[-1].split('.')[:-1] ) + '.log'
   lv = logging.DEBUG
   fmt='%(asctime)s:%(name)s:%(levelname)s: %(message)s'
   logging.basicConfig(level=lv, format=fmt, datefmt='%Y/%m/%d-%H:%M:%S',
                                 filename = log_file, filemode='w')
   LG = logging.getLogger('main')
   if not is_cron: log_help.screen_handler(LG, lv=lv)
##############################################################################
   LG.info(f'Starting: {__file__}')
   import sys
   try:
      args = sys.argv[1:]
      folder = args[0]
      domain = args[1]
      days = int(args[2])
      start,end = args[3].split(',')
      start = int(start)
      end = int(end)
   except:
      print('Error reading input, please use the following format:')
      print('  $ python generate_config.py <day> <start>,<end>')
      print('\nEx: Today since 12:00 until 20:00')
      print('  $ python generate_config.py 0 12,20')
      print('\nEx: Tomorrow since 8:00 until 20:00')
      print('  $ python generate_config.py 1 8,20')
      exit()
   # Time
   start_date = dt.datetime.now().replace(hour=start,minute=0,second=0)
   start_date = start_date.replace(microsecond=0)
   start_date = start_date + dt.timedelta(days = days)
   end_date = start_date.replace(hour=end)
   timedelta = 1   # hourly data
   main(folder, domain, start_date, end_date, timedelta, lglv='debug',
        left=-17, right=8, bottom=30, top=48,
        Ncores=13, wait4batch=40)
