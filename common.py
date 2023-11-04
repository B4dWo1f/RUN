#!/usr/bin/python3
# -*- coding: UTF-8 -*-

from configparser import ConfigParser, ExtendedInterpolation
from os.path import expanduser
from os import listdir
from os.path import isfile, join
import datetime as dt
import os
here = os.path.dirname(os.path.realpath(__file__))
import logging
import log_help
LG = logging.getLogger(__name__)
import json
fmt = '%d/%m/%Y-%H:%M'


class RunParams(object):
   def __init__(self, start_date, end_date, GFS_timedelta, domains,
                      domain_folder, GFS_data_folder,
                      output_folder, plots_folder,
                      leftlon, rightlon, toplat, bottomlat, Ncores,wait4batch):
      UTCshift = dt.datetime.now() - dt.datetime.utcnow()
      UTCshift = dt.timedelta(hours = round(UTCshift.total_seconds()/3600))
      self.UTCshift = UTCshift
      LG.info(f'Got UTCshift: {self.UTCshift}')

      # In this house we obey the laws of UTC time!!!
      self.start_date = min(start_date, end_date) - UTCshift
      self.end_date   = max(start_date, end_date) - UTCshift
      msg = f'User wants: {self.start_date.strftime(fmt)} - '
      msg += f'{self.end_date.strftime(fmt)}'
      LG.info(msg)
      # aux = []
      # for h in daily_hours:
      #    aux.append(dt.datetime.combine(dt.date(2007,10,12),h) - UTCshift)
      # self.daily_hours = [h.time()  for h in aux]
      # LG.info(f'Daily hour mask: {[x.hour for x in self.daily_hours]}')
      self.domains = domains
      self.domain_folder = domain_folder
      self.GFS_timedelta = GFS_timedelta
      self.GFS_data_folder = GFS_data_folder
      com = f'mkdir -p {self.GFS_data_folder}'
      LG.info(f'creating folder {self.GFS_data_folder}')
      self.output_folder = output_folder
      os.system(com)
      com = f'mkdir -p {self.output_folder}'
      LG.info(f'creating folder {self.output_folder}')
      os.system(com)
      self.plots_folder = plots_folder
      com = f'mkdir -p {self.plots_folder}'
      LG.info(f'creating folder {self.plots_folder}')
      os.system(com)
      self.leftlon   = leftlon
      self.rightlon  = rightlon
      self.toplat    = toplat
      self.bottomlat = bottomlat
      self.Ncores = Ncores
      self.wait4batch = wait4batch
   def __str__(self):
      msg = f'Forecast for: {self.start_date} - {self.end_date}\n'
      # msg += 'Hours mask: '
      # hs = [h.strftime('%H:%M') for h in self.daily_hours]
      # msg += ' - '.join(hs) + '\n'
      msg += f'GFS data: {self.GFS_data_folder}\n'
      msg += f'GFS timedelta: {self.GFS_timedelta}\n'
      msg += f'domain: {self.domain_folder}\n'
      msg += f'regions: {self.domains}\n'
      msg += f'Longitudes: {self.leftlon} - {self.rightlon}\n'
      msg += f'Latitudes : {self.toplat} - {self.bottomlat}\n'
      msg += f'Ncores: {self.Ncores}'
      return msg


class GFS_Data(object):
   """
   This class wraps all the information about the FTP server
   """
   def __init__(self, server, folder_root, Nparallel=2, folder_data=''):
      self.server = server
      self.folder_root = folder_root
      self.folder_data = folder_data
      self.Nparallel = Nparallel
      if len(self.folder_data) > 0:
         self.folder = f'{self.folder_root}/{self.folder_data}'
   def add_data_folder(self,folder_data):
      self.folder_data = folder_data
      if len(self.folder_data) > 0:
         self.folder = f'{self.folder_root}/{self.folder_data}'
   def __str__(self):
      msg =  f'  Nparallel: {self.Nparallel}\n'
      msg += f'     server: {self.server}\n'
      msg += f'root folder: {self.folder_root}\n'
      msg += f'date folder: {self.folder_data}'
      try: msg += f'\n     folder: {self.folder}'
      except AttributeError: pass
      return msg


class DataNotAvailableError(Exception):
   def __init__(self, data, message='Data is not available'):
      # Call the base class constructor with the parameters it needs
      self.message = message
      self.data = data
      super().__init__(message)
      # self.errors = errors
   def __str__(self):
      return f"{self.data.strftime(fmt)} -> {self.message}"


def load(fname='config.ini'):
   """
   Load the config options and return it as a class
   """
   LG.info(f'Loading config file: {fname}')
   if not os.path.isfile(fname): return None
   config = ConfigParser(inline_comment_prefixes='#')
   config._interpolation = ExtendedInterpolation()
   config.read(fname)

   # Date and time options
   try:
      LG.debug('Trying to read start/end dates')
      start_date = config['run']['start_date']
      start_date = dt.datetime.strptime(start_date, '%d/%m/%Y-%H:%M')
      end_date = config['run']['end_date']
      end_date = dt.datetime.strptime(end_date, '%d/%m/%Y-%H:%M')
      LG.info(f'Run for: {start_date} - {end_date}')
   except KeyError:
      LG.debug('Start/end dates not specified, going for now+days_run')
      days_run = int(config['run']['days_run'])
      start_date = dt.datetime.now()
      start_date = start_date.replace(hour=8,minute=0,second=0,microsecond=0)
      end_date = start_date + dt.timedelta(days=days_run)
      LG.info(f'Run for: {start_date} - {end_date}')

   # try:
   #    hours_mask = config['run']['daily_hours'].split(',')    #XXX TODO
   #    hours_mask = [dt.time(*tuple(map(int,t.split(':')))) for t in hours_mask]
   # except KeyError:
   #    LG.warning('Hours mask (daily_hours) not provided! Using 0-24')
   #    hours_mask = [dt.time(0),dt.time(23)]

   # Domain
   domains_run = tuple(map(int,config['run']['domains'].split(',')))
   domain_folder = expanduser(config['run']['domain_folder'])
   GFS_data_folder = expanduser(config['run']['GFS_data_folder'])
   GFS_timedelta = int(config['run']['GFS_timedelta'])
   output_folder = expanduser(config['run']['output_folder'])
   plots_folder = expanduser(config['run']['plots_folder'])
   leftlon   = float(config['run']['leftlon'])
   rightlon  = float(config['run']['rightlon'])
   toplat    = float(config['run']['toplat'])
   bottomlat = float(config['run']['bottomlat'])

   # System
   Ncores = int(config['run']['Ncores'])
   wait4batch = float(config['run']['wait4batch'])

   R = RunParams(start_date, end_date, GFS_timedelta, domains_run,
                 domain_folder, GFS_data_folder,output_folder,plots_folder,
                 leftlon, rightlon, toplat, bottomlat, Ncores, wait4batch)

   return R


