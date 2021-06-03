#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import gfs
import datetime as dt

date = dt.datetime.now()
dateGFS = gfs.get_GFS_calc(date)
print('Time:',date)
print('Corresponding GFS:',dateGFS)
