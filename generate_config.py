#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from os.path import expanduser
import datetime as dt

import sys
try:
   args = sys.argv[1:]
   days = int(args[0])
   start,end = args[1].split(',')
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


fmt = '%d/%m/%Y-%H:%M'

# Time
start_date = dt.datetime.now().replace(hour=start,minute=0,second=0)
start_date = start_date.replace(microsecond=0)
start_date = start_date + dt.timedelta(days = days)
end_date = start_date.replace(hour=end)
daily_hours = [start,end]

# Folders
folder = '~/METEO'
domain = 'Spain6_1'
domains = 1,2   #XXX this should be automatic
Ncores = 26
left,right = -17,8
bottom,top = 30,48
wait4batch = 40  # Minutes to keep trying for the last GFS batch

config = configparser.ConfigParser()
config['run'] = {}
config['run']['start_date']  = start_date.strftime(fmt)
config['run']['end_date']    = end_date.strftime(fmt)
config['run']['daily_hours'] = '8,20'   #','.join([str(x) for x in daily_hours])
config['run']['domain_folder'] = f'{folder}/RUN/Domains/{domain}'
config['run']['domains'] = ','.join([str(x) for x in domains])
config['run']['GFS_data_folder'] = f'{folder}/dataGFS'
config['run']['output_folder'] = f'/storage/WRFOUT/{domain}'
config['run']['plots_folder'] = f'/storage/PLOTS/{domain}'
config['run']['leftlon']   = str(left)
config['run']['rightlon']  = str(right)
config['run']['toplat']    = str(top)
config['run']['bottomlat'] = str(bottom)
config['run']['Ncores']    = str(Ncores)
config['run']['wait4batch']    = str(wait4batch)
config['run']['log_level'] = 'debug'

with open('config.ini', 'w') as configfile:
    config.write(configfile)
