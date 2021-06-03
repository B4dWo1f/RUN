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


def get_GFS_calc(date, shift = 3+40/60):
   """
   This function returns the expected datetime of the latest GFS data for the
   provided date.
   date: [datetime.datetime] target date
   shift: [float] hours since batch name. ej: data from batch 06 is not usually
          available before 9:40
   Returns:
   dateGFS: [datetime.datetime]
   """
   # hours = [0,6,12,18,date.hour]
   # inds = list(np.argsort(hours))
   # ind = inds.index(4)-1
   # dateGFS = date.replace( hour=hours[ind], minute=0, second=0,microsecond=0)
   GFS = list(range(0,24,6))
   hours = [h+shift for h in GFS] + [date.hour+date.minute/60]
   inds = list(np.argsort(hours))
   ind = (inds.index(4)-1) % len(GFS)
   dateGFS = date.replace(hour=GFS[ind],minute=0,second=0,microsecond=0)
   if dateGFS>date: dateGFS -= dt.timedelta(days=1)
   return dateGFS


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


def get_files(Params, dates, data_folder='data',wait4batch=40):
   """
   Downloads the GFS data for each of the dates provided in dates
   wait4batch: [float] Minutes to keep trying for the last GFS batch
   """
   data_folder = Params.GFS_data_folder
   leftlon = Params.leftlon
   rightlon = Params.rightlon
   toplat = Params.toplat
   bottomlat = Params.bottomlat
   base_url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?'
### Ensure all data from same GFS run XXX check!
   now = dt.datetime.utcnow()
   batch = min([min(dates), now])
   batch = get_GFS_calc(batch)
   LG.info(f"Tentative GFS batch: {batch.strftime('%Y%m%d/%H')}")
   delta = dt.timedelta(hours=6)
   is_data_present = False
   for i_batch in range(4):   # loop for previous GFS batchs
      batch = batch - i_batch*delta

      start = dt.datetime.now()
      waiting = (dt.datetime.now() - start).total_seconds()/60
      LG.info(f'Is data present: {is_data_present} trying since {waiting:.0f}s ago')
      while not is_data_present and waiting < wait4batch:
         waiting = (dt.datetime.now() - start).total_seconds()/60
         LG.info(f'Trying batch: {batch.strftime(fmt)} ({waiting:.0f}/{wait4batch}m)')
         # Generate all the urls
         folders,files = [],[]
         for date in dates:
            folder,fname = ncep_folder_fname(batch,date)
            files.append(fname)
         LG.info(f'Files: {files[0]} - {files[-1]}')
         is_data_present = checker(folder,files)
         if not is_data_present:
             LG.info('Waiting 60 seconds')
             sleep(60)

         # if is_data_present:
         #    # If all data is present break the GFS batchs loop
         #    LG.info('all the data is in the server, starting download')
         #    exit = True
         #    break
         # else: exit = False
      if is_data_present: break
      LG.critical('FAILED. Falling back to previous batch')


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
      url,fname = URL_subregion(batch,date,leftlon,rightlon,toplat,bottomlat)
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


def URL_subregion(batch,date,leftlon,rightlon,toplat,bottomlat):
   """
   returns the url of the file with the GFS forecast for day `date` made on day
   `batch`, cropped by subregion
   For GFS subregion, this is the url
   batch: date of the GFS data
   date: date of the forecast
   """
   folder, fname = ncep_folder_fname(batch,date)
   url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?'
   url += f'file={fname}'
   url += f'&all_lev=on&all_var=on'
   url += f'&subregion='
   url += f'&leftlon={leftlon}'
   url += f'&rightlon={rightlon}'
   url += f'&toplat={toplat}'
   url += f'&bottomlat={bottomlat}'
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


if __name__ == '__main__':
    date = dt.datetime.now()
    dateGFS = get_GFS_calc(date)
    print('Time:',date)
    print('Corresponding GFS:',dateGFS)
