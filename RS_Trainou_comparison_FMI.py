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

from radiosondes.radiosondes.radiosondes import RS_Trainou2019, comparison_RS

plt.ion()

RS_old_paths = [
            r"C:\Trainou_2019\radiosonde\FMI_20200409\AC_TRN_RINGO_20190616_006.mto",
            r"C:\Trainou_2019\radiosonde\FMI_20200409\AC_TRN_RINGO_20190617_014.mto"
            ]

P0s = [
        1006.1,
       1000.5
       ]

z0s = [
        134,
       134
       ]

starts = [
           datetime(2019,6,16,9,25),
          datetime(2019,6,17,12,30)
          ]

# only take the second day, so it doesn't take too long
# RS_paths = [RS_paths[1]]
# P0s = [P0s[1]]
# z0s = [z0s[1]]
# RS_old_paths = [RS_old_paths[1]]


# prepare plotting different methods for each flightday
RS_methods = {'mirror':[True, False], 'ipol':['combined', True, False]}

def get_RSmethods(RS_dict):
    ls=[]
    for mir in RS_dict['mirror']:
        for ip in RS_dict['ipol']:
            ls.append((mir, ip))
    return ls

method_list = get_RSmethods(RS_methods)


RS_daylist = []
for P0, z0, start, path_old in zip(P0s, z0s, starts, RS_old_paths):
    RS_list = []
    for mirror, ipol in method_list:
        RS = RS_Trainou2019(path=path_old, delim=' ', P0=P0, z0=z0, mirror=mirror, P_Calc=False, ipol=ipol, geometric=False, dirksenTv=True, start=start)
        # RS.df.rename(columns={'P_old':'P_LSCE'}, inplace=True)
        # RS_old = RS_Trainou2019(path_old, delim=' ', init_P_calc=False)
        # RS.df = RS.df.merge(RS_old.df['P_old'], left_index=True, right_index=True, how='left')
        RS.calculation(drop='RH')  # start pressure calculation here
        RS.df = RS.df[~RS.df.index.duplicated(keep='first')]  # "~" bitwise operator, turns False into True and vice versa
        RS_list.append(RS)
    RS_daylist.append(RS_list)

#%%
#############################################################################
# Save plots?
saveplots = pathlib.Path(r"C:\Trainou_2019\radiosonde\FMI_20200409")
saveplots = False  # if you want to save the plots: comment this out! dont overwrite the path with a True value!

# Plots relative delta P instead of absolute delta P?
relP = False
T_RH_z_plot=True
dP_plot=True
P_plot=True

if saveplots:
    #Set special string variable for plotname
    special = 'Std_z0Pstart_Geopotential_iterAuto_dirksenTv'

    if relP:
        saveplots = saveplots.parent / (saveplots.name + '_relP')

else:
    special = ''


#%% add Arduino data
Arduino_paths = [
        pathlib.Path(r"C:\Trainou_2019\radiosonde\FMI_20200409\Datalogger\20190616.B00"),
        pathlib.Path(r"C:\Trainou_2019\radiosonde\FMI_20200409\Datalogger\20190617.B00")
        ]

col_names = ['year', 'month', 'day', 'hour', 'minute', 'second', 'T1', 'prsP',
             'T2', 'flag1', 'flag2', 'r1', 'r2', 'r3', 'r4', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'flag3', 'flag4']
Arduino_dfs = []
for path in Arduino_paths:
    df = pd.read_csv(path, delimiter=',', names=col_names)
    df.set_index(pd.to_datetime(df.iloc[:, 0:6]), inplace=True)
    df = df.iloc[:, 6:]
    df = df[~df.index.duplicated(keep='first')]  # "~" bitwise operator, turns False into True and vice versa

    Arduino_dfs.append(df)

titles = [
           'GUF003 June 12\nRS006',
          'GUF003 June 17\nRS014'
          ]
if relP:
    yname = r'$\frac{\Delta P}{P_{static}}$ [%]'
else:
    yname = r'$\Delta$P [hPa]'
ynameP = 'P [hPa]'
xname = 'time [dd hh:mm]'


P_list = ['GUF', 'old_height', 'Ard']  #, 'Ard', 'LSCE',

P_dict = {'GUF':{'col':'P',
                 'name':r'$P_{GUF}$ [hPa]'},
          # 'LSCE':{'col':'P_LSCE',
          #         'name':r'$P_{LSCE}$ [hPa]'},
          'old_height':{'col':'P_old',
                 'name':r'$P_{U.S. stand. atm.}$ [hPa]'},
          'Ard':{'col':'prsP',
                 'name':r'$P_{Arduino}$ [hPa]'}
          }

P_dict = {key: P_dict[key] for key in P_list}


Alt_list = ['RS_GPS', 'RS_GEOPOT']  # , 'Ard_GPS'

Alt_dict = {'RS_GPS':{'col':'Alt',
                      'name':'r$z_{RS}$ [m]'},
            'RS_GEOPOT':{'col':'geopotential altitude',
                         'name':r'$z_{geopot, RS}$ [m]'},
            'Ard_GPS':{'col':'Altitude',
                       'name':r'$z_{Arduino}$ [m]'}
            }

Alt_dict = {key: Alt_dict[key] for key in Alt_list}

T_list = ['RS_T']

T_dict = {'RS_T':{'col':'T',
                  'name': r'$T_{RS}$ [K]'},
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


set_dont_mirror = [{'P_dict':'Ard'}] #use dict name and key as a new dictionary

save_df = comparison_RS(RS_daylist, method_list, Arduino_dfs,
                  titles, yname, ynameP, xname,
                  P_list, P_dict, Alt_list, Alt_dict, T_list, T_dict, RH_list, RH_dict,
                  dict_list, dict_dict, tidy_cols_base,
                  set_dont_mirror, saveplots, relP, special=special,
                  T_RH_z_plot=T_RH_z_plot, dP_plot=dP_plot, P_plot=P_plot)


    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # axes[0].legend()
    # axes[3].legend()

# %% export data
save_df.keys()
dir(save_df['m_dfs'])
dir(save_df['tidy_dfs'])

len(save_df['m_dfs'])
len(save_df['tidy_dfs'])

type(Arduino_dfs[1].index[0])
len(Arduino_dfs[1].columns)
len(col_names) - 6

# only use the seconds df (June 17, RS014)
m_df = save_df['m_dfs'][1]
tidy_df = save_df['tidy_dfs'][1]

m_df[4].columns
tidy_df[4].columns

df = RS_daylist[1][4].df  # this should be the not mirrored - ipol True DataFrame; see method_list[4]
print(df.columns)
