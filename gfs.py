#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
here = os.path.dirname(os.path.realpath(__file__))
import datetime as dt
import numpy as np
from ftplib import FTP
from multiprocessing import Pool
import requests
from bs4 import BeautifulSoup
from time import sleep
from random import random
import log_help
import logging
LG = logging.getLogger(__name__)
log_help.screen_handler(LG, lv=logging.DEBUG)

fmt = '%d/%m/%Y-%H:%M'


def get_GFS_calc(date):
   """
   This function returns the expected datetime of the latest-calculated data
   for the date provided.
   date: datetime.datetime object (at least it needs 'hour' attribute)
   """
   hours = [0,6,12,18,date.hour]
   inds = list(np.argsort(hours))
   ind = inds.index(4)-1
   return date.replace( hour = hours[ind], minute=0, second=0,microsecond=0)


def checker(folder,files,mode='ftp'):
   """
   Check if the remote `folder` is hosting all the `files`
   USE mode=ftp!!!!
   """
   LG.debug(f'Checking files in {mode} server')
   folder = f'pub/data/nccf/com/gfs/prod/{folder}/atmos'
   # Get remote list of files
   if mode == 'ftp':
      base = 'ftp.ncep.noaa.gov'
      ftp_connection = make_login(base)
      try:
         server_files = change_directory(ftp_connection, folder)
      except myFTPError:
         LG.warning('Folder does not exist')
         return False
   elif mode == 'http':
      LG.critical('SUUUPER buggy. Use mode=ftp')
      url = f'https://nomads.ncep.noaa.gov/{folder}'
      server_files = check_files_http(url,files)
   # Check if all the necessary files are present
   if set(files) <= set(server_files):
      LG.info('All files present')
      return True
   else:
      missing_files = set(files) - (set(files) & set(server_files))
      LG.warning(f'Missing {len(missing_files)} files')
      return False


def get_files(Params, dates, data_folder='data'):
   """
   Downloads the GFS data for each of the dates provided in dates
   """
   data_folder = Params.GFS_data_folder
   base_url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?'
### Ensure all data from same GFS run XXX check!
   now = dt.datetime.utcnow()
   batch = min([min(dates), now])
   batch = get_GFS_calc(batch)
   LG.info(f"Tentative GFS batch: {batch.strftime('%Y%m%d/%H')}")
   delta = dt.timedelta(hours=6)
   for i_batch in range(4):   # loop for previous GFS batchs
      batch = batch - i_batch*delta
      LG.info(f'Trying calculation: {batch.strftime(fmt)}')
      # Generate all the urls
      folders,files = [],[]
      for date in dates:
         folder,fname = ncep_folder_fname(batch,date)
         files.append(fname)
      is_data_present = checker(folder,files)

      if is_data_present:
         # If all data is present break the GFS batchs loop
         LG.info('all the data is in the server, starting download')
         break

### Download
   LG.info(f"Correct GFS batch: {batch.strftime('%Y%m%d/%H')}")
   LG.info(f"Starting download of GFS data")
   with open(f'{here}/plots/batch.txt','w') as f_batch:
       f_batch.write(batch.strftime(fmt))
   LG.critical(f'{Params.output_folder}/batch.txt')
   with open(f'{Params.output_folder}/batch.txt','w') as f_batch:
       f_batch.write(batch.strftime(fmt))
   # Prepare all inputs for parallel download
   urls,files = [],[]
   for date in dates:
      url,fname = URL_subregion(batch,date)
      urls.append(url)
      files.append(f'{data_folder}/{fname}')
   all_inps = [(u,f) for u,f in zip(urls,files)]
   ## for inp in all_inps:
   ##    status = download_file(*inp)
   ##    # sleep(1)
   # Parallel option. The server seems to penalize this option
   npool = 5  # int(min([len(all_inps)/2,10])) # XXX  check depending on ISP
   LG.debug(f'Starting download with {npool} requests in parallel')
   pool = Pool(npool)
   remaining = all_inps
   Ntries = 5
   cont = 0
   # Download until all files are downloaded
   while len(remaining) > 0 and cont < Ntries:
      status = pool.starmap(download_file, remaining)
      remaining = [inp for inp,stat in zip(remaining, status) if not stat]
      if len(remaining) > 0:
         LG.warning(f'Missing {len(remaining)}. Retry')
      else: break
      #TODO add sleep here
      cont += 1
   LG.info('All GFS data downloaded')
   return True


### Harcoded server data
def ncep_folder_fname(batch,date):
   """
   returns the correct folder in the server and file name of the file with the
   GFS forecast for day `date` made on day `batch`
   """
   remote_foler =  f"gfs.{batch.strftime('%Y%m%d/%H')}"
   ind = int((date-batch).total_seconds()/3600)
   fname = f"gfs.t{batch.strftime('%H')}z.pgrb2.0p25.f{ind:03d}"
   return remote_foler,fname


def URL_subregion(batch,date):
   """
   returns the url of the file with the GFS forecast for day `date` made on day
   `batch`, cropped by subregion
   TODO: limits are hardcoded for my domain, they should be an input
       leftlon: -17
       rightlon: 8
       toplat: 48
       bottomlat: 30
   For GFS subregion, this is the url
   batch: date of the GFS data
   date: date of the forecast
   """
   folder, fname = ncep_folder_fname(batch,date)
   url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?'
   url += f'file={fname}'
   url += f'&all_lev=on&all_var=on'
   url += f'&subregion='
   url += f'&leftlon=-17'
   url += f'&rightlon=8'
   url += f'&toplat=48'
   url += f'&bottomlat=30'
   url += f'&dir=/{folder}/atmos'
   return url, fname
########################


#################################################################################
#                                    HTTP(s)                                    #
#################################################################################
def check_files_http(url,filenames):
    """
    Checks if a remote file exists
    files are too heavy, better to download the html and look for the
    files in the text
    """
    #XXX hotfix
    # filenames shoud be the name of the file in the server. NOT local
    filenames = [f.split('/')[-1] for f in filenames]
    LG.debug(f'Checking {url} for {len(filenames)} files')
    headers = {'User-Agent': 'Mozilla/5.0'}
    req = requests.get(url, headers=headers, stream=False)
    LG.info(f'HTTP response: {req.status_code}')
    # req = requests.get(url, headers=headers, timeout=(20,30),stream=False)
    if req.status_code == 404:
       LG.debug('Error 404')
       return False
    with open('gfs.html', 'w') as f:
       f.write(req.text)
    soup = BeautifulSoup(req.text, 'lxml')
    server_files = [f.text for f in soup.findAll('a')]
    return server_files


def download_file(url,fout):
   """
   Download file from url and store it in fout
   """
   sleep(5*random())  # to avoid banning
   LG.info(f'Downloading {url}')
   req = requests.get(url)
   fid = open(fout, 'wb')
   fid.write(req.content)
   fid.close()
   #TODO add check
   if not os.path.isfile(fout) or os.path.getsize(fout) == 0:
      return False
   LG.info(f'saved in {fout}')
   return True


#################################################################################
#                                      FTP                                      #
#################################################################################

class myFTPError(Exception):
   pass

def make_login(server,credentials=None,Ntries=10):
   """
   Make a python-FTP connection to the FPT server
   """
   LG.debug(f'Connecting to {server}')
   out = False
   cont = 0
   while not out and cont <= Ntries:
      try:
         LG.info(f'FTP attempt {cont}/{Ntries}')
         ftp = FTP(server, timeout=5)
         out = True
      except:
         LG.warning(f'Unable to connect to {server} ({cont}/{Ntries})')
         sleep(1)
      cont += 1
   if not out: raise myFTPError('The server cannot be reached')
   LG.debug('Connection successful')

   if credentials == None:
      ftp.login()
      LG.info('FTP Login sucessful')
   else: raise myFTPError(f'FTP login with credentials is not implemented yet')
   return ftp

def change_directory(ftp, folder):
   """
   change folder to `folder` destination. One step at a time, for easier
   debugging
   Return the list of files in the latest folder
   """
   LG.debug(f'cd to {folder}')
   for f in list(filter(None, folder.split('/'))):
      try:
         ftp.cwd(f)
         LG.debug(f'cd {f}')
      except: raise myFTPError(f'Unable to cd to folder {f}')
   LG.info(f'Changed directory to {folder}')
   return ftp.nlst()

