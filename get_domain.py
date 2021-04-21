#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""
crappy workaround to extract the domain from run_rasp.sh
"""

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
# if not is_cron: log_help.screen_handler(LG, lv=lv)
LG.info(f'Starting: {__file__}')
##############################################################################
import common
import datetime as dt

fini = 'config.ini'

R = common.load(fini)


print(R.domain_folder)
print(R.GFS_data_folder)
print(R.output_folder)
print(R.Ncores)
