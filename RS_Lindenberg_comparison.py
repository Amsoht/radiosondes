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
import matplotlib.pyplot as plt


from radiosondes.radiosondes.radiosondes import RS_Trainou2019, RS_Lindenberg2018, comparison_RS
from radiosondes.manip.randomsample import randomblox





#read in data and provide base values for P calculation

#Lindenberg 2018:
RS_paths = [{'GPS':r"C:\Lindenberg_2018\20181128_GUF003\Radiosonde\Lindenberg-Forschung-2_20181128_083057_GpsResults.csv",
            'PTU':r"C:\Lindenberg_2018\20181128_GUF003\Radiosonde\Lindenberg-Forschung-2_20181128_083057_PtuResults.csv"
            }]

# RS = RS_Lindenberg2018(paths=RS_paths[0], P_Calc=True, ipol=False, geometric=True)

#%%


#prepare plotting different methods for each flightday
RS_methods = {'mirror':[True, False], 'ipol':['combined', True, False]}
rndm_missing = False

def get_RSmethods(RS_dict):
    ls=[]
    for mir in RS_dict['mirror']:
        for ip in RS_dict['ipol']:
            ls.append((mir, ip))
    return ls

method_list = get_RSmethods(RS_methods)


RS_daylist = []
for path in RS_paths:
    RS_list = []
    flag = True
    for mirror, ipol in method_list:
        RS = RS_Lindenberg2018(paths=path, mirror=mirror, P_Calc=False, ipol=ipol, geometric=True, dirksenTv=True)
        if rndm_missing:
            if flag:
                rndm_blx = randomblox(RS.df, nstarts=15, blocklim=400)
                flag = False
                
            RS.df = RS.df[~RS.df.index.isin(rndm_blx)]
            
        RS.calculation()
        RS_list.append(RS)
        
    RS_daylist.append(RS_list)

#%%
#############################################################################
# Save plots?
saveplots = pathlib.Path(r"C:\Trainou_2019\radiosonde\Vergleich_Druckberechnung\Lindenberg2018")
saveplots = False #commt out, if you dont want to save the plots

# Plots relative delta P instead of absolute delta P?
relP = False
T_RH_z_plot=False
dP_plot=True
P_plot=False

special = '' #probably dont need this

if saveplots:
    #Set special string variable for plotname
    special = 'Std_Geometric_iterAuto_dirksenTv'
    
    if relP:
        saveplots = saveplots.parent / (saveplots.name + '_relP')
    
    
#%% add Arduino data
Arduino_paths = [
        pathlib.Path(r"C:\Lindenberg_2018\20181128_GUF003\Arduino\20181128_table.302")
        ]
parse_dates = ['Date/Time']

Arduino_dfs = []
for path in Arduino_paths:
    df = pd.read_csv(path, delimiter=' *, *', parse_dates=parse_dates, index_col=parse_dates, engine='python')
    Arduino_dfs.append(df)

titles = ['GUF003 November 28\nRS41']
if relP:
    yname = r'$\frac{\Delta P}{P_{static}}$ [%]'
else:
    yname = r'$\Delta$P [hPa]'
ynameP = 'P [hPa]'
xname = 'time [dd hh:mm]'


P_list = ['GUF', 'old_height', 'old_sensor'] #, 'Ard'

P_dict = {'GUF':{'col':'P',
                 'name':r'$P_{GUF}$ [hPa]'},
          'old_height':{'col':'PressureFromHeight [hPa]',
                  'name':r'$P_{RS from height}$ [hPa]'},
          'old_sensor':{'col':'SensorPressure [hPa]',
                 'name':r'$P_{RS from sensor}$ [hPa]'},
          'Ard':{'col':'prsP',
                 'name':r'$P_{Arduino}$ [hPa]'}
          }

P_dict = {key: P_dict[key] for key in P_list}


Alt_list = ['RS_GPS', 'RS_GEOPOT'] #, 'Ard_GPS' #add Ard_geopot?

Alt_dict = {'RS_GPS':{'col':'Height [m]',
                      'name':'r$z_{RS}$ [m]'},
            'RS_GEOPOT':{'col':'geopotential altitude',
                         'name':r'$z_{geopot, RS}$ [m]'},
            'Ard_GPS':{'col':'Altitude',
                       'name':r'$z_{Arduino}$ [m]'}
            }

Alt_dict = {key: Alt_dict[key] for key in Alt_list}


T_list = ['RS_T']

T_dict = {'RS_T':{'col':'Temperature [K]', #unit is °C
                  'name': r'$T_{RS}$ [°C]'},
          'Tv':{'col':'Tv',
                'name':r'$T_v$ [K]'}
          }


RH_list = ['RS_RH']

RH_dict = {'RS_RH':{'col':'Humidity [%]',
                    'name':r'$RH_{RS}$ [%]'}
           }


dict_list = [P_dict, Alt_dict, T_dict, RH_dict]
dict_dict = {'P_dict':P_dict,
             'Alt_dict':Alt_dict,
             'T_dict':T_dict,
             'RH_dict':RH_dict}

tidy_cols_base = [dic[key]['col'] for dic in dict_list for key in dic] #get relevant column names from dictionaries


set_dont_mirror = [{'P_dict' : 'Ard'}] #use dict name and key as a new dictionary


save_df = comparison_RS(RS_daylist, method_list, Arduino_dfs, 
                  titles, yname, ynameP, xname, 
                  P_list, P_dict, Alt_list, Alt_dict, T_list, T_dict, RH_list, RH_dict, 
                  dict_list, dict_dict, tidy_cols_base, 
                  set_dont_mirror, saveplots, relP, special,
                  T_RH_z_plot, dP_plot, P_plot)




    
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # axes[0].legend()
    # axes[3].legend()


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    