#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
here = os.path.dirname(os.path.realpath(__file__))
is_cron = False
fmt = '%d/%m/%Y-%H:%M'
################################# LOGGING ####################################
import logging
import log_help
log_file = here+'/'+'.'.join( __file__.split('/')[-1].split('.')[:-1] ) + '.log'
lv = logging.DEBUG
logging.basicConfig(level=lv,
                 format='%(asctime)s %(name)s:%(levelname)s - %(message)s',
                 datefmt='%Y/%m/%d-%H:%M',
                 filename = log_file, filemode='w')
LG = logging.getLogger('main')
if not is_cron: log_help.screen_handler(LG, lv=lv)
LG.info(f'Starting: {__file__}')
##############################################################################
import common
import datetime as dt

fini = 'config.ini'

R,FTP = common.load(fini)


################
# namelist.wps #
################
fname = f'{R.domain_folder}/namelist.wps'
with open(fname,'r') as f:
   all_text = []
   for line in f.read().strip().splitlines():
      if line.strip().startswith('max_dom'):
         line = f' max_dom = {len(R.domains)},'
      elif line.strip().startswith('start_date'):
         line = f"start_date = "
         for _ in R.domains:
            line += f"\'{R.start_date.strftime('%Y-%m-%d_%H:%M:%S')}\', "
         line = ' ' + line.strip()
      elif line.strip().startswith('end_date'):
         line = f"end_date   = "
         for _ in R.domains:
            line += f"\'{R.end_date.strftime('%Y-%m-%d_%H:%M:%S')}\', "
         line = ' ' + line.strip()
      all_text.append(line)

with open(f'{here}/namelist.wps','w') as f:
   f.write('\n'.join(all_text)+'\n\n')


##################  XXX nooooo no no no no no no no
# namelist.input #  XXX for fuck's sake!!! no!!!!
##################  XXX rewrite this
duration = R.end_date - R.start_date
duration = duration.total_seconds()
# time_control
run_days = int(duration/3600/24)
duration -= run_days * 3600 * 24
run_hours = int(duration/3600)
duration -= run_hours * 3600
run_minutes = int(duration/60)
duration -= run_minutes * 60
run_seconds = int(duration)
start_year   = R.start_date.year
start_month  = R.start_date.month
start_day    = R.start_date.day
start_hour   = R.start_date.hour
start_minute = R.start_date.minute
start_second = R.start_date.second
end_year   = R.end_date.year
end_month  = R.end_date.month
end_day    = R.end_date.day
end_hour   = R.end_date.hour
end_minute = R.end_date.minute
end_second = R.end_date.second

fname = f'{R.domain_folder}/namelist.input'
with open(fname,'r') as f:
   all_text = []
   for line in f.read().strip().splitlines():
      # Ignored for real cases?
      # https://www2.mmm.ucar.edu/wrf/users/namelist_best_prac_wrf.html#run_days
      if line.strip().startswith('run_days'):
         line = f'run_days                 = {run_days},'
      elif line.strip().startswith('run_hours'):
         line = f'run_hours                = {run_hours},'
      elif line.strip().startswith('run_minutes'):
         line = f'run_minutes              = {run_minutes},'
      elif line.strip().startswith('run_seconds'):
         line = f'run_seconds              = {run_seconds},'
      if line.strip().startswith('start_year'):
         line = f'start_year               = '
         for _ in R.domains:
            line += f'{start_year},     '
         line = line.strip()
      elif line.strip().startswith('start_month'):
         line = f'start_month              = '
         for _ in R.domains:
            line += f'{start_month:02d},       '
         line = line.strip()
      elif line.strip().startswith('start_day'):
         line = f'start_day                = '
         for _ in R.domains:
            line += f'{start_day:02d},       '
         line = line.strip()
      elif line.strip().startswith('start_hour'):
         line = f'start_hour               = '
         for _ in R.domains:
            line += f'{start_hour:02d},       '
         line = line.strip()
      elif line.strip().startswith('start_minute'):
         line = f'start_minute             = '
         for _ in R.domains:
            line += f'{start_minute:02d},       '
         line = line.strip()
      elif line.strip().startswith('start_second'):
         line = f'start_second             = '
         for _ in R.domains:
            line += f'{start_second:02d},       '
         line = line.strip()
      elif line.strip().startswith('end_year'):
         line = f'end_year                 = '
         for _ in R.domains:
            line += f'{end_year},     '
         line = line.strip()
      elif line.strip().startswith('end_month'):
         line = f'end_month                = '
         for _ in R.domains:
            line += f'{end_month:02d},       '
         line = line.strip()
      elif line.strip().startswith('end_day'):
         line = f'end_day                  = '
         for _ in R.domains:
            line += f'{end_day:02d},       '
         line = line.strip()
      elif line.strip().startswith('end_hour'):
         line = f'end_hour                 = '
         for _ in R.domains:
            line += f'{end_hour:02d},       '
         line = line.strip()
      elif line.strip().startswith('end_minute'):
         line = f'end_minute               = '
         for _ in R.domains:
            line += f'{end_minute:02d},       '
         line = line.strip()
      elif line.strip().startswith('end_second'):
         line = f'end_second               = '
         for _ in R.domains:
            line += f'{end_second:02d},       '
         line = line.strip()
      all_text.append(line)


with open(f'{here}/namelist.input','w') as f:
   f.write('\n'.join(all_text)+'\n\n')
