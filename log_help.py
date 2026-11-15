#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import logging
import datetime as dt
import utils as ut
from pathlib import Path
from functools import wraps
from time import perf_counter


def batch_logger(script_path, batch, domain='', is_cron=True, log_dir=None):
   """
   Sets up a logger that writes to a log file named according to the current
   date and batch
   Optionally adds screen logging if not run by cron
   """
   log_dir = Path(log_dir)
   if not log_dir.is_absolute():
      # Resolve relative to the file where this function lives
      log_dir = Path(__file__).resolve().parent / log_dir
   script_path = Path(script_path).resolve()
   base_name = script_path.stem
   # log_name = f"{base_name}_{domain}_GFS{batch.strftime('%H')}"
   names = [base_name]
   if len(domain) > 0: names.append(domain)
   names.append(f"GFS{batch.strftime('%H')}")
   log_name = '_'.join(names)

   log_dir = Path(log_dir) if log_dir else script_path.parent
   ut.check_directory(log_dir)

   log_path  = log_dir / f"{log_name}.log"
   perf_path = log_dir / f"{log_name}.perform"

   # Create and configure loggers
   logger = logging.getLogger('main')
   logger.setLevel(logging.INFO)

   perform = logging.getLogger('perform')
   perform.setLevel(logging.DEBUG)

   formatter = logging.Formatter(
       '%(asctime)s %(name)s:%(levelname)s - %(message)s',
       datefmt='%Y/%m/%d-%H:%M:%S'
   )

   # Main logger file handler
   file_handler = logging.FileHandler(log_path, mode='w')
   file_handler.setFormatter(formatter)
   logger.addHandler(file_handler)

   # Performance logger file handler
   perf_handler = logging.FileHandler(perf_path, mode='w')
   perf_handler.setFormatter(formatter)
   perform.addHandler(perf_handler)

   if not is_cron:
       screen_handler(logger, lv=logging.DEBUG)

   msg = f"Started script: {script_path} | domain: {domain} | batch: {batch}"
   logger.info(msg)
   return logger, perform



def screen_handler(lg=None,lv='debug',fmt='%(name)s -%(levelname)s- %(message)s'):
   """
     This function adds a screenHandler to a given logger. If no logger is
     specified, just return the handler
   """
   if lv == 'debug': lv = logging.DEBUG
   elif lv == 'info': lv = logging.INFO
   elif lv == 'warning': lv = logging.WARNING
   elif lv == 'error': lv = logging.ERROR
   elif lv == 'critical': lv = logging.CRITICAL
   sh = logging.StreamHandler()
   sh.setLevel(lv)
   fmt = logging.Formatter(fmt)
   sh.setFormatter(fmt)
   if lg != None: lg.addHandler(sh)
   return sh

## Timer Decorator
def timer(lg_progression,lg_performance):
   """Logs the execution time of a certain function to the provided logger"""
   def real_timer(wrapped):
      @wraps(wrapped)   # preserve metadata about the function being decorated
      def inner(*args, **kwargs):
         t0 = perf_counter()
         lg_progression.debug(f'Entering {wrapped.__name__}')
         ret = wrapped(*args, **kwargs)
         lg_progression.debug(f'Finished {wrapped.__name__}')
         t1 = perf_counter()
         lg_performance.debug(f'Time for {wrapped.__name__}: {t1-t0:.4f}s')
         return ret
      return inner
   return real_timer


def log2screen(lg,lv='info'):
   """ This decorator adds *temporarily* a screenHandler to a given logger """
   def do_it(wrapped):
      def inner(*args, **kwargs):
         sh = screen_handler(lg,lv=lv)  # add ScreenHandler
         ret = wrapped(*args, **kwargs)
         lg.removeHandler(sh)     # remove ScreenHandler
         return ret
      return inner
   return do_it


def disable(lg):
   """Temporarily raise the log level to CRITICAL to avoid over logging"""
   def do_it(wrapped):
      def inner(*args, **kwargs):
         lv = lg.getEffectiveLevel()
         lg.setLevel(logging.CRITICAL)
         ret = wrapped(*args, **kwargs)
         lg.setLevel(lv)
         return ret
      return inner
   return do_it

def disable2(lg):
   """Temporarily raise the log level to INFO to avoid over logging"""
   def do_it(wrapped):
      def inner(*args, **kwargs):
         lv = lg.getEffectiveLevel()
         lg.setLevel(logging.INFO)
         ret = wrapped(*args, **kwargs)
         lg.setLevel(lv)
         return ret
      return inner
   return do_it
