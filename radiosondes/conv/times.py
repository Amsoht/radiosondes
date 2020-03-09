# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 09:25:30 2020

@author: Thomas Wagenhäuser

module name: times
"""
import numpy as np
import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
import pathlib
from metpy.calc import height_to_geopotential as h_t_g
from metpy.units import units
import random
from random import randrange

from ..plot_sc.plot_building import make_one_spine_visible



def datetime_to_secofday(datetimedata=[], refdatetime=None): #datetimes or Timestamps, supports overnight timeseries
    """Convert datetime objects to seconds of the day

    arguments:
        - datetimedata, list or pandas.series of datetimes, now also handels single Timestamps
    keywordarguments:
        - refdatetime, seconds will be calculated relative to midnight of refdatetime
            by default, refdatetime=datetimedata[0]"""
    try:
        if refdatetime:
            midnight = refdatetime.normalize()
        else:
            midnight = datetimedata[0].normalize()
        secofday=[]
        for val in datetimedata:
            secofday.append(timedelta.total_seconds(val-midnight))
        return secofday
    except TypeError:
        midnight = datetimedata.normalize()
        secofday = timedelta.total_seconds(datetimedata-midnight)
        return secofday





class Datetime_Creator():
    """Created on Tue Nov 19 16:05:23 2019

        @author: Thomas Wagenhäuser
        purpose: convert seconds of the day to hh:mm:ss
        First create an instance of Datetime_Creator with the start Date.
        Later feed secofday single data into instance.create_datetime()
    """
    def __init__(self, Dates=[], secofday=[]):
        self.dates = Dates
        self.secofday = secofday

    def create_date(self, secofday):
        y = self.dates[0]
        m = self.dates[1]

        if secofday < 86400:
            d = self.dates[2]
        elif secofday >= 86400:
            d = self.dates[2] + 1
        return dt.date(y,m,d)

    def create_time(self, secofday):

        m, s = divmod(secofday, 60)

        h, m = divmod(m, 60)

        h = int(h)
        m = int(m)
        s = int(s)

        return dt.time(h,m,s)

    def create_datetime(self, secofday):
        return datetime.combine(self.create_date(secofday), self.create_time(secofday))
