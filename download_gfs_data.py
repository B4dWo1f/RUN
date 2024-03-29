#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
# here = os.path.dirname(os.path.realpath(__file__))
is_cron = False
is_cron = bool( os.getenv('RUN_BY_CRON') )
fmt = '%d/%m/%Y-%H:%M'
import common
import datetime as dt
import gfs
import logging
LG = logging.getLogger('main')


def main(fini = 'config.ini'):
   R = common.load(fini)

   LG.info(f'{R.start_date} - {R.end_date}')
   current_date = R.start_date
   step = dt.timedelta(hours=R.GFS_timedelta)  #XXX should be in config.ini
   max_tries = 5
   dates_calc = []
   while current_date <= R.end_date:
      # if R.daily_hours[0] <= current_date.time() <= R.daily_hours[1]:
      dates_calc.append((current_date))
      current_date += step

   cont = 0
   while cont < max_tries:
      try:
         got_all_files = gfs.get_files(R, dates_calc, R.GFS_data_folder,
                                       wait4batch=R.wait4batch)
         if got_all_files: cont = 2*max_tries  # XXX dumb ways to exit...
      except:
         cont += 1

if __name__ == '__main__':
################################# LOGGING ####################################
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
   main()
