# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:57:47 2020

@author: Thomas Wagenhäuser
"""

#TODO:
# check dont_mirror
# in interpol_data_gaps 'combined' only use reasonable columns!
# - decorate loop over datasets with plot functions
# put Deltas in a dict
# put P in a dict
# - take Thomas data
#   try to recreate his data with different E_sat -> this seems to be it!
#   try with RH not overwritten
#   try with geometric altitude
#   compare P_LSCE to P_old
#   add old Pressure column
#   add diff to all plots
# Deal with mirror and Arduino data,
# place axes and subplot adjustment


import numpy as np
import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
import tkinter as tk
from tkinter import filedialog
import pathlib
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
# import clevercsv as csv

from radiosondes.radiosondes.radiosondes import RS_Trainou2019, RS_Lindenberg2018, RS_Sodankyla2018, comparison_RS
from radiosondes.manip.randomsample import randomblox

plt.ion()

#read in data and provide base values for P calculation

#Lindenberg 2018:
RS_paths = [
    # r"C:\Sodankyla_2018\Sodankyla_20180621\Radiosonde\FLEDT20180621085245.tsv",
    r"C:\Sodankyla_2018\Sodankyla_20180625\Radiosonde\FLEDT20180625074941.tsv",
    # r"C:\Sodankyla_2018\Sodankyla_20180625\Radiosonde\FLEDT20180625084548.tsv",
    # r"C:\Sodankyla_2018\Sodankyla_20180625\Radiosonde\FLEDT20180625155553.tsv",
    # r"C:\Sodankyla_2018\Sodankyla_20180626_AC3\Radiosonde\FLEDT20180626154917.tsv"
    ]

titles = [
           # 'GUF003 June 21\nRS92',
           'GUF003 June 25a\nRS92',
          # 'GUF003 June 25b\nRS92',
          # 'GUF003 June 25c\nRS92',
           # 'GUF003 June 26\nRS92'
          ]

# RS = RS_Sodankyla2018(path=RS_paths[0], P_Calc=True, ipol=False, geometric=True)




#prepare plotting different methods for each flightday
RS_methods = {'mirror':[True, False], 'ipol':['combined', True, False]}
rndm_missing = True

def get_RSmethods(RS_dict):
    ls=[]
    for mir in RS_dict['mirror']:
        for ip in RS_dict['ipol']:
            ls.append((mir, ip))
    return ls

method_list = get_RSmethods(RS_methods)

# plt.figure()
RS_daylist = []
for path in RS_paths:
    RS_list = []
    flag=True
    for mirror, ipol in method_list:
        RS = RS_Sodankyla2018(path=path, mirror=mirror, P_Calc=False, ipol=ipol, geometric=True, dirksenTv=True, freq=timedelta(seconds=1))
        if rndm_missing:
            if flag:
                rndm_blx = randomblox(RS.df, nstarts=15, blocklim=200, random_state=9)  # 9 for 25a
                flag = False
            RS.df.loc[RS.df.index.isin(rndm_blx), :] = np.nan

        RS.calculation()

        # plt.plot(RS.df.index, RS.df[RS.cols_out['P']['col']], 'o', label="mirror: {}, ipol: {}".format(mirror, ipol))


        if rndm_missing:
            RS_old = RS_Sodankyla2018(path=path, mirror=mirror, P_Calc=True, ipol=ipol, geometric=True, dirksenTv=True)
            RS.df = RS.df.merge(RS_old.df['P_old'], left_index=True, right_index=True, how='right', suffixes=('_gaps',''))
        RS_list.append(RS)
    RS_daylist.append(RS_list)
# plt.legend()
# plt.show()
#%%
#############################################################################
# Save plots?
saveplots = pathlib.Path(r"C:\Trainou_2019\radiosonde\Joram_Mail_20200401")
saveplots = False

# Plots relative delta P instead of absolute delta P?
relP = True
T_RH_z_plot=True
dP_plot=True
P_plot=True



if saveplots:
    #Set special string variable for plotname
    special = 'Std_z0Pstart_Geometric_iterAuto_ModemTv'

    if relP:
        saveplots = saveplots.parent / (saveplots.name + '_relP')
else:
    special = ''


#%% add Arduino data
Arduino_paths = [
        # pathlib.Path(r"C:\Sodankyla_2018\Sodankyla_20180621\Arduino\20180621_noerror.300"),
        pathlib.Path(r"C:\Sodankyla_2018\Sodankyla_20180625\Arduino\20180625_table.304"),
        # pathlib.Path(r"C:\Sodankyla_2018\Sodankyla_20180625\Arduino\20180625_table.304"),
        # pathlib.Path(r"C:\Sodankyla_2018\Sodankyla_20180625\Arduino\20180625_table.304"),
        # pathlib.Path(r"C:\Sodankyla_2018\Sodankyla_20180626_AC3\Arduino\20180626_table.301")
]

parse_dates = ['Date/Time']

Arduino_dfs = []
for path in Arduino_paths:
    try:
        df = pd.read_csv(path, delimiter=' *, *', parse_dates=parse_dates, index_col=parse_dates, engine='python', skiprows=[1])
    except Exception:
        df = pd.read_csv(path, delimiter=' *, *', engine='python', skiprows=[1])
        df.index = pd.to_datetime(df['Date'], format='%y/%m/%d %H:%M:%S')
    Arduino_dfs.append(df)


if relP:
    yname = r'$\frac{\Delta P}{P_{static}}$ [%]'
else:
    yname = r'$\Delta$P [hPa]'
ynameP = 'P [hPa]'
xname = 'time [dd hh:mm]'


P_list = ['GUF', 'old_height'] #, 'old_sensor', 'Ard'

P_dict = {'GUF':{'col':'P',
                 'name':r'$P_{GUF}$ [hPa]'},
          'old_height':{'col':'P_old',
                  'name':r'$P_{old}$ [hPa]'},
          'Ard':{'col':'prsP', #'Pressure' for June21, else: 'prsP'
                 'name':r'$P_{Arduino}$ [hPa]'}
          }

P_dict = {key: P_dict[key] for key in P_list}


Alt_list = ['RS_GPS', 'RS_GEOPOT'] #, 'Ard_GPS' add Ard_geopot?

Alt_dict = {'RS_GPS':{'col':'Height',
                      'name':'r$z_{RS}$ [m]'},
            'RS_GEOPOT':{'col':'geopotential altitude',
                         'name':r'$z_{geopot, RS}$ [m]'},
            'Ard_GPS':{'col':'Altitude',
                       'name':r'$z_{Arduino}$ [m]'}
            }

Alt_dict = {key: Alt_dict[key] for key in Alt_list}


T_list = ['RS_T']

T_dict = {'RS_T':{'col':'T', #unit is °C
                  'name': r'$T_{RS}$ [°C]'},
          'Tv':{'col':'Tv',
                'name':r'$T_v$ [K]'}
          }


RH_list = ['RS_RH']

RH_dict = {'RS_RH':{'col':'RH',
                    'name':r'$RH_{RS}$ [%]'}
           }


dict_list = [P_dict, Alt_dict, T_dict, RH_dict]
dict_dict = {'P_dict':P_dict,
             'Alt_dict':Alt_dict,
             'T_dict':T_dict,
             'RH_dict':RH_dict}

tidy_cols_base = [dic[key]['col'] for dic in dict_list for key in dic] #get relevant column names from dictionaries

if 'Ard' in P_list:
    set_dont_mirror = [{'P_dict':'Ard'}] #use dict name and key as a new dictionary
else:
    set_dont_mirror = [{}]

save_df = comparison_RS(RS_daylist, method_list, Arduino_dfs,
                  titles, yname, ynameP, xname,
                  P_list, P_dict, Alt_list, Alt_dict, T_list, T_dict, RH_list, RH_dict,
                  dict_list, dict_dict, tidy_cols_base,
                  set_dont_mirror, saveplots, relP, special,
                  T_RH_z_plot, dP_plot, P_plot)


    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # axes[0].legend()
    # axes[3].legend()

save_df['m_dfs'][0][3].head()

save_df['m_dfs'][0][3]['n_iterations_Tv'].describe()
