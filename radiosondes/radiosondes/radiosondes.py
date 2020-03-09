# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 16:58:29 2020

@author: Thomas Wagenhäuser
"""

#TODO:
# - Composite Class RS_calcP and RS_read_data Classes
# - use dicts for column names


import sys
import re
import numpy as np
import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
from dateutil import parser
import tkinter as tk
from tkinter import filedialog
import pathlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from metpy.calc import height_to_geopotential as h_t_g

from ..conv.times import Datetime_Creator, datetime_to_secofday
from ..calc.p_from_z_T_RH import P_From_GPS_T_RH
from ..plot_sc.plot_building import make_one_spine_visible



def get_idx_col_from_parse_dates(parse_dates):
    if type(parse_dates) != bool:
        if type(parse_dates[0]) == list:
            return '_'.join(parse_dates[0])
        return parse_dates
    else:
        return



def get_kwarg_from_file(path, kwarg):
    try:
        with open(path, 'r') as f:
            for l in f:
                if kwarg in l:
                    key, arg = re.split(': \s*', l)
                    return arg
            return
    except Exception as excp:
        print(f"\n\nError: {excp}\nCouldn't get keyword argument...")



class RS_calcP(object):

    def __init__(self, P0=1013.15, z0=0.0, mirror=False, P_Calc=True, ipol='combined', freq='S', geometric=True, start=False, cols_in={}, cols_out={}, drop=False, dirksenTv=True,**kwargs): #geometric: treat as is: True, calc geopotential: False
        self.mirror = mirror
        self.P_Calc = P_Calc
        self.ipol = ipol
        self.freq = freq
        self.P0 = P0
        self.z0 = z0
        self.start = start
        self.geometric = geometric
        self.dirksenTv = dirksenTv
        # self.drop_outliers(['RH']) #hard-coded so far
        if cols_in != {}:
            self.cols_in = cols_in
        else:
            self.cols_in = {'z':{'col':'GPS_ALT',
                              'unit':'[m]'},
                         'RH':{'col':'RH',
                               'unit':'[%]'},
                         'T':{'col':'T',
                              'unit':'[°C]'}
                         }
        if cols_out != {}:
            self.cols_out = cols_out
        else:
            self.cols_out = {'P':       {'col':'P', #same name needed for initialization in radiosondes.py
                                         'unit':'[hPa]'},
                             'z_geopot':{'col':'geopotential altitude',
                                         'unit':'[m]'},
                             'Tv':      {'col':'Tv',
                                         'unit':'[K]'}}
        if self.cols_in['T']['unit'] == '[K]':
            self.df[cols_in['T']['col']] -= 273.15 #in order to get temperature in °C
            self.cols_in['T']['unit'] = '[°C]'
        if P_Calc == True:
            self.calculation(**kwargs)


    def calculation(self, drop=False, **kwargs):
        #prepare dataframe for P update
        #self.df.rename(columns={'P':'P_old'}, inplace=True) #get this from a dict
        # self.df['P'] = np.full(self.df.shape[0], np.NaN) #initialize empty P column
        # self.df.loc[self.df.index[0],'P'] = self.P0 #enter pressure start value
        # self.df.loc[self.df.index[0],'GPS_ALT'] = self.z0 #enter altitude start value
        if drop:
            try:
                if type(drop) != list:
                    drop = [drop]
                for var in drop:
                    self.drop_outliers(var)
            except Exception as excp:
                print(f"\n\nError: {excp}\nAn error occured while trying to drop outliers. Check 'drop' keyword.")
        #maybe add a Temperature start value (so far taken from radiosonde)

        #update P
        self.P_Calc_Obj = P_From_GPS_T_RH(self.df, P0=self.P0, z0=self.z0, start=self.start, mirror=self.mirror, cols_in=self.cols_in, cols_out=self.cols_out, **kwargs) #create
        self.P_Calc_Obj.update_P(ipol=self.ipol, freq=self.freq, geometric=self.geometric, dirksenTv=self.dirksenTv) #perform
        self.df = self.P_Calc_Obj.df #extend DataFrame (!no copy!)


    def drop_outliers(self, cols=[], wins=[], diffs=[]):
        wins_default = 10
        diffs_default = 10
        if cols == []:
            return
        if len(wins) != len(cols):
            for num in range(len(cols)):
                wins.append(wins_default)
        if len(diffs) != len(cols):
            for num in range(len(cols)):
                diffs.append(diffs_default)

        print('Dropping outliers',
              f'Variables: {cols}',
              f'Window size left and right each: {wins}',
              f'Sigma factor: {diffs}',
              '...',
              sep='\n')

        for col, win, diff in zip(cols, wins, diffs):
            self.df['mean_r'] = self.df[col].rolling(win).mean().shift(1)
            self.df['std_r'] = self.df[col].rolling(win).std().shift(1)
            #flip df
            inverse_df = self.df.sort_index(ascending=False)
            inverse_df['mean'] = inverse_df[col].rolling(win).mean()
            inverse_df['std'] = inverse_df[col].rolling(win).std()
            inverse_df.sort_index(ascending=True, inplace=True)
            self.df['mean_l'] = inverse_df['mean'].shift(-1)
            self.df['std_l'] = inverse_df['std'].shift(-1)

            self.df['diff_r'] = self.df[col] - self.df['mean_r']
            self.df['diff_l'] = self.df[col] - self.df['mean_l']
            crit_r = self.df['diff_r'] > diff * self.df['std_r']
            crit_l = self.df['diff_l'] > diff * self.df['std_l']

            # self.df['outlier'] = crit_r & crit_l
            # self.df.loc[self.df['outlier'],col] = np.nan
            outlier = crit_r & crit_l
            self.df.loc[outlier,col] = np.nan
        print('Done dropping outliers!\n')



class RS_Sodankyla2018(RS_calcP):
    """RS object with RS_DataFrame. Set columnnames here"""

    def __init__(self, path='', delim='\t', parse_dates=False, init_P_Calc=True, P_Calc=False, cols_in={}, cols_out={}, P0=False, z0=False, **kwargs):
        """create RS_Sodankyla2018 object with """
        self.path = pathlib.Path(path)
        self.parse_dates=parse_dates
        self.df = pd.DataFrame()
        self.read_data(delim, skiprows=47) #skip header comments
        self.make_utc_idx()

        if cols_in != {}:
            self.cols_in = cols_in
        else:
            self.cols_in = {'z':{'col':'Height',
                              'unit':'[m]'},
                         'RH':{'col':'RH',
                               'unit':'[%]'},
                         'T':{'col':'T',
                              'unit':'[K]'}
                         }
        if cols_out != {}:
            self.cols_out = cols_out
        else:
            self.cols_out = {'P':       {'col':'P', #same name needed for initialization in radiosondes.py
                                         'unit':'[hPa]'},
                             'z_geopot':{'col':'geopotential altitude',
                                         'unit':'[m]'},
                             'P_old_combined': {'col':'P',
                                                'unit':'[hPa]'},
                             'Tv':      {'col':'Tv',
                                         'unit':'[K]'}
                             }

        #rename old Pressure column
        #TODO: write a function in order to create ICARTT variables dictionary
        if self.cols_out['P_old_combined']['col'] == self.cols_out['P']['col']:
            new_col_name = self.cols_out['P_old_combined']['col'] + '_old'
            self.df.rename(columns={self.cols_out['P_old_combined']['col']:new_col_name}, inplace=True)
            self.cols_out['P_old_combined']['col'] = new_col_name


        if init_P_Calc:
            if P0:
                self.P0 = P0
            else:
                self.P0 = self.df.loc[self.df.index[0],self.cols_out['P_old_combined']['col']]
            if z0:
                self.z0 = z0
            else:
                self.z0 = self.df.loc[self.df.index[0],self.cols_in['z']['col']]
            super().__init__(cols_in=self.cols_in, cols_out=self.cols_out, P0=self.P0, z0=self.z0, P_Calc=P_Calc, **kwargs)

        #self.df.rename(columns={'Esat':'e_s', 'Alt':'GPS_ALT'}, inplace=True) #needs to be changed in conversions.py


    def read_data(self, delim='\t', skiprows=False):
        """read initial data"""

        print('\nRead in RS data...')
        if self.path == '':
            self.path = pathlib.Path(filedialog.askopenfilename())
        try:
            idx_col = get_idx_col_from_parse_dates(self.parse_dates)
            self.df = pd.read_csv(self.path, skiprows=skiprows, delimiter=delim, parse_dates=self.parse_dates, index_col=idx_col, skipinitialspace=True, na_values=['-32768', '-33041'])
        except Exception as excp: #maybe auto try with other delimiter?
            print(f'\n\nError: {excp}\nAn error occured while reading in RS data. Check delimiter in radiosondes.py!')
            sys.exit()
        print('Done!')


    def make_utc_idx(self, sec_data_col=''):
        if sec_data_col == '':
            sec_data_col = 'time'
        df = self.df
        start = parser.parse(get_kwarg_from_file(path=self.path, kwarg='Launch time'))
        start = pd.Timestamp(start)
        fix_time = 'Date_Time'
        start_sec = datetime_to_secofday(start)
        secofday = df[sec_data_col] - df[sec_data_col][0] + start_sec
        timeconverter = Datetime_Creator(Dates=[start.year,start.month,start.day])
        df[fix_time] = secofday.apply(timeconverter.create_datetime)
        df.set_index(fix_time, inplace=True)



class RS_Lindenberg2018(RS_calcP):

    def __init__(self,paths={'PTU':'','GPS':''}, delim=';', parse_dates=['DataSrvTime'], init_P_Calc=True, P_Calc=False, cols_in={}, cols_out={}, P0=False, z0=False, **kwargs):
        self.paths = {key: pathlib.Path(path) for key,path in paths.items()}
        self.parse_dates=parse_dates
        self.dfs = {'PTU':pd.DataFrame(),
                    'GPS':pd.DataFrame()}
        self.read_data(delim)

        if cols_in != {}:
            self.cols_in = cols_in
        else:
            self.cols_in = {'z':{'col':'Height [m]',
                              'unit':'[m]'},
                         'RH':{'col':'Humidity [%]',
                               'unit':'[%]'},
                         'T':{'col':'Temperature [K]',
                              'unit':'[K]'}
                         }
        if cols_out != {}:
            self.cols_out = cols_out
        else:
            self.cols_out = {'P':       {'col':'P', #same name needed for initialization in radiosondes.py
                                         'unit':'[hPa]'},
                             'z_geopot':{'col':'geopotential altitude',
                                         'unit':'[m]'},
                             'P_old_from_height': {'col':'PressureFromHeight [hPa]',
                                                   'unit':'[hPa]'},
                             'P_old_from_sensor': {'col':'SensorPressure [hPa]',
                                                   'unit':'[hPa]'},
                             'Tv':      {'col':'Tv',
                                         'unit':'[K]'}}
        self.fix_timestamp()

        if init_P_Calc:
            self.df = self.dfs['PTU']
            if P0:
                self.P0 = P0
            else:
                self.P0 = self.df.loc[self.df.index[0],self.cols_out['P_old_from_sensor']['col']]
            if z0:
                self.z0 = z0
            else:
                self.z0 = self.df.loc[self.df.index[0],self.cols_in['z']['col']]
            super().__init__(cols_in=self.cols_in, cols_out=self.cols_out, P0=self.P0, z0=self.z0, P_Calc=P_Calc, **kwargs)


    def read_data(self, delim=';'):
        """read initial data"""

        print('\nRead in RS data...')
        for key, path in self.paths.items():
            if path.is_dir():
                path = pathlib.Path(filedialog.askopenfilename(title=f'Please select {key} file'))
            try:
                idx_col = get_idx_col_from_parse_dates(self.parse_dates)
                df = pd.read_csv(path, delimiter=delim, parse_dates=self.parse_dates, index_col=idx_col)
                self.dfs[key] = df
            except Exception as excp: #maybe auto try with other delimiter?
                print(f'\n\nError: {excp}\nAn error occured while reading in RS data. Check delimiter in radiosondes.py!')
                sys.exit()
        print('Done!')

    def fix_timestamp(self, start=False, sec_data_col='', dfs={}, merge_cols={}):
        if start == False:
            start = self.dfs['PTU'].index[0] - timedelta(minutes=16, seconds=20)

        if sec_data_col == '':
            sec_data_col = 'RadioRxTimePk [s]'

        df = self.dfs['PTU']
        fix_time = 'Date_Time'
        start_sec = datetime_to_secofday(start)
        secofday = df[sec_data_col] - df[sec_data_col][0] + start_sec
        timeconverter = Datetime_Creator(Dates=[start.year,start.month,start.day])
        df[fix_time] = secofday.apply(timeconverter.create_datetime)
        df.set_index(fix_time, inplace=True)



class RS_Trainou2019(RS_calcP):
    """RS object with RS_DataFrame. Set columnnames here"""

    def __init__(self, path='', delim='\t', parse_dates=[['Date','Time']], init_P_Calc=True, P_Calc=False, cols_in={}, cols_out={}, P0=False, z0=False, drop='RH', **kwargs):
        """create RS_Trainou2019 object with """
        self.path = pathlib.Path(path)
        self.parse_dates=parse_dates
        self.df = pd.DataFrame()
        self.read_data(delim)

        if cols_in != {}:
            self.cols_in = cols_in
        else:
            self.cols_in = {'z':{'col':'Alt',
                              'unit':'[m]'},
                         'RH':{'col':'RH',
                               'unit':'[%]'},
                         'T':{'col':'T',
                              'unit':'[°C]'}
                         }
        if cols_out != {}:
            self.cols_out = cols_out
        else:
            self.cols_out = {'P':       {'col':'P', #same name needed for initialization in radiosondes.py
                                         'unit':'[hPa]'},
                             'z_geopot':{'col':'geopotential altitude',
                                         'unit':'[m]'},
                             'P_old_from_height': {'col':'P',
                                                   'unit':'[hPa]'},
                             'Tv':      {'col':'Tv',
                                         'unit':'[K]'}
                             }

        #rename old Pressure column
        #TODO: write a function in order to create ICARTT variables dictionary
        if self.cols_out['P_old_from_height']['col'] == self.cols_out['P']['col']:
            new_col_name = self.cols_out['P_old_from_height']['col'] + '_old'
            self.df.rename(columns={self.cols_out['P_old_from_height']['col']:new_col_name}, inplace=True)
            self.cols_out['P_old_from_height']['col'] = new_col_name


        if init_P_Calc:
            if P0:
                self.P0 = P0
            else:
                self.P0 = self.df.loc[self.df.index[0],self.cols_out['P_old_from_height']['col']]
            if z0:
                self.z0 = z0
            else:
                self.z0 = self.df.loc[self.df.index[0],self.cols_in['z']['col']]
            super().__init__(cols_in=self.cols_in, cols_out=self.cols_out, P0=self.P0, z0=self.z0, P_Calc=P_Calc,drop=drop, **kwargs)

        #self.df.rename(columns={'Esat':'e_s', 'Alt':'GPS_ALT'}, inplace=True) #needs to be changed in conversions.py


    def read_data(self, delim='\t'):
        """read initial data"""

        print('\nRead in RS data...')
        if self.path == '':
            self.path = pathlib.Path(filedialog.askopenfilename())
        try:
            idx_col = get_idx_col_from_parse_dates(self.parse_dates)
            self.df = pd.read_csv(self.path, delimiter=delim, parse_dates=self.parse_dates, index_col=idx_col)
        except Exception as excp: #maybe auto try with other delimiter?
            print(f'\n\nError: {excp}\nAn error occured while reading in RS data. Check delimiter in radiosondes.py!')
            sys.exit()
        print('Done!')



def dz_from_dp(row, Tv_name, p1_colname, p2_colname):
    """calculate difference in height from pressure difference, based on the hypsometric equation

    This i just a round about function!"""

    Rd = 287.05
    g = 9.80616
    return Rd/g*row[Tv_name]*np.log(row[p1_colname]/row[p2_colname])



def comparison_RS(RS_daylist, method_list, Arduino_dfs,
                  titles, yname, ynameP, xname,
                  P_list, P_dict, Alt_list, Alt_dict, T_list, T_dict, RH_list, RH_dict,
                  dict_list, dict_dict, tidy_cols_base,
                  set_dont_mirror, saveplots, relP, special,
                  T_RH_z_plot=True, dP_plot=True, P_plot=True):

    dont_mirror = []
    for entry in set_dont_mirror:
        for dic,key in entry.items():
            try:
                val = dict_dict[dic][key]
            except KeyError:
                continue
            dont_mirror += [key, val['col'], val['name']]

    save_df = {'m_dfs':[],
               'tidy_dfs':[]}
    for RS_list, a_df, title in zip(RS_daylist, Arduino_dfs, titles): #get all (6) methods for each day

        figs = []
        axes_set = []
        names = []

        if T_RH_z_plot:
            #plot ascent and descent T, RH profile
            z_fig, z_axs = plt.subplots(nrows=2, ncols=1, figsize=[6, 9], sharex=True)
            z_fig.suptitle(title, size='x-large', weight='bold')

            asc_df = RS_list[-1].df[:RS_list[-1].df[P_dict['GUF']['col']].idxmin()]
            des_df = RS_list[-1].df[RS_list[-1].df[P_dict['GUF']['col']].idxmin():]

            key = {'alt':Alt_dict['RS_GPS']['col'],
                   'T':T_dict['RS_T']['col'],
                   'RH':RH_dict['RS_RH']['col']}

            for i, z_ax in enumerate(z_axs.flat):
                RH_twinax = z_ax.twinx()
                if i == 0:
                    z_ax.plot(asc_df[key['alt']], asc_df[key['T']], label = 'T ascent')
                    z_ax.plot(des_df[key['alt']], des_df[key['T']], label = 'T descent')
                    z_ax.set_ylabel('temperature')

                    RH_twinax._get_lines.prop_cycler = z_ax._get_lines.prop_cycler #for better autocolors
                    RH_twinax.plot(asc_df[key['alt']], asc_df[key['RH']], label = 'RH ascent')
                    RH_twinax.plot(des_df[key['alt']], des_df[key['RH']], label = 'RH descent')
                    RH_twinax.set_ylabel('RH')

                    handles, labels = z_ax.get_legend_handles_labels()
                    RH_handles, RH_lables = RH_twinax.get_legend_handles_labels()
                    handles += RH_handles
                    labels += RH_lables
                    z_ax.legend(handles, labels, loc='upper right', markerscale=10)

                if i == 1:
                    asc_df = asc_df[[col for k, col in key.items()]].groupby(key['alt']).mean() #ascent groupby altitude
                    des_df = des_df[[col for k, col in key.items()]].groupby(key['alt']).mean() #descent groupby altitude

                    des_df.rename(columns = {key['T']:key['T']+'_descent',
                                             key['RH']:key['RH']+'_descent'},
                                  inplace=True)
                    comb_df = pd.concat([asc_df,des_df], axis=1, sort=False).interpolate()
                    comb_df['T_diff'] = comb_df[key['T']] - comb_df[key['T']+'_descent']
                    comb_df['RH_diff'] = comb_df[key['RH']] - comb_df[key['RH']+'_descent']

                    z_ax.plot(comb_df.index, comb_df['T_diff'], label = r'$\Delta$T')
                    z_ax.set_ylabel('temperature')
                    z_ax.set_xlabel('altitude')

                    RH_twinax._get_lines.prop_cycler = z_ax._get_lines.prop_cycler #for better autocolors
                    RH_twinax.plot(comb_df.index, comb_df['RH_diff'], label = r'$\Delta$RH')
                    RH_twinax.set_ylabel('RH')
                    z_ax.set_xlabel('altitude')

                    handles, labels = z_ax.get_legend_handles_labels()
                    RH_handles, RH_lables = RH_twinax.get_legend_handles_labels()
                    handles += RH_handles
                    labels += RH_lables
                    z_ax.legend(handles, labels, loc='upper right', markerscale=10)

            figs.append(z_fig)
            axes_set.append([z_axs])
            names.append('RH_T_z')

        if dP_plot or P_plot:
            dfs = [RS.df.copy() for RS in RS_list]
            # for df in dfs:
            #     df[r'Delta_P'] = df['P']-df['P_old']

            #merge dataframes, be careful with overlapping columnnames?
            if 'Ard' in P_list:
                m_dfs = []
                for df in dfs:
                    if (df.index.tzinfo is not None):
                        df.index = df.index.tz_convert(None) #overwrite timezone data, because Pandas can't deal with it
                    m_df = a_df.join(df, how='outer') #sort=False by default
                    m_dfs.append(m_df)
            else:
                m_dfs = dfs

            save_df['m_dfs'].append(m_dfs)

            tidy_dfs = []
            for df in m_dfs:

                if 'Ard' in P_list:
                    ndf = df.copy()[~df[P_dict['Ard']['col']].isna().values & ~df['P'].isna().values] #TODO: use a key from a dict
                else:
                    ndf = df.copy()
                # ndf = df.copy()[~df['Lat'].isna().values]
                # ndf = df.copy()[~df['secofday'].isna().values]

                #calculate pressure differences
                d_names = [] #actually it would be sufficent to produce them only once. However, this way, I can use one loop for two things
                dalt_names = []
                for num, P in enumerate(P_list[:2]): #only compare GUF and LSCE to everything
                    for i in range(num+1, len(P_list)): #skip already calculated differences
                        d_name = ' '.join(P_dict[P]['name'].split(' ')[:-1]) +' - '+ ' '.join(P_dict[P_list[i]]['name'].split(' ')[:-1])
                        dalt_name = f'$\Delta$z for ' + d_name
                        ndf[d_name] = ndf[P_dict[P]['col']] - ndf[P_dict[P_list[i]]['col']]
                        d_names.append(d_name)
                        Tv_name = T_dict['Tv']['col']
                        if P_plot:
                            ndf[dalt_name] = df.apply(dz_from_dp, args=(Tv_name, P_dict[P]['col'], P_dict[P_list[i]]['col']), axis=1)
                            dalt_names.append(dalt_name)

                        if (P in dont_mirror) or (P_list[i] in dont_mirror):
                            dont_mirror.append(d_name)
                            dont_mirror.append(dalt_name)

                tidy_cols = tidy_cols_base.copy()
                tidy_cols += d_names
                tidy_cols += dalt_names
                ndf = ndf[tidy_cols]
                # ndf.rename(columns={'Altitude':'alt_Arduino', 'prsP':'P_Arduino','P_old':'P_rs','P':'P_calc', 'GPS_ALT':'alt_rs', 'Delta_P':'Delta_P_rs'}, inplace=True)
                # ndf['Delta_P_Arduino'] = ndf['P_calc'] - ndf['P_Arduino']
                # ndf['Delta_P_LSCE'] = ndf['P_calc'] - ndf['P_LSCE']
                tidy_dfs.append(ndf)

            save_df['tidy_dfs'].append(tidy_dfs)
            #tidy_dfs[0].head()


            if dP_plot:
                #plot pressure differences
                fig, axs = plt.subplots(nrows=3, ncols=2, figsize=[22, 10], sharex='col', sharey=True)
                fig.suptitle(title, size='x-large', weight='bold')

                axes = axs.T.reshape(-1)
                twinax=[] #list of twinaxes made from the array of axes
                rhax=[]
                for ax in axes:
                    twinax.append(ax.twinx())
                    rhax.append(ax.twinx())

                figs.append(fig)
                axes_set.append([axes, twinax, rhax])
                names.append('deltaP')


            if P_plot:
                #plot pressure on log scale
                figP, axsP = plt.subplots(nrows=3, ncols=2, figsize=[22, 10], sharex='col', sharey=True)
                figP.suptitle(title, size='x-large', weight='bold')
                axesP = axsP.T.reshape(-1)

                daltax = []
                for ax in axesP:
                    daltax.append(ax.twinx())

                daltax[0].get_shared_y_axes().join(*[ax for ax in daltax])
                figs.append(figP)
                axes_set.append([axesP, daltax])
                # axes_set.append(daltax)
                names.append('P')



            for i, df in enumerate(tidy_dfs):

                if dP_plot:

                    #plot RS_T
                    col, name = T_dict['RS_T'].items()
                    a = twinax[i].plot(df.index, df[col[1]], 'o', ms=0.125)
                    twinax[i].set_ylabel(name[1], color=a[0].get_color())

                    #plot RH
                    # Offset the right spine of par2.  The ticks and label have already been
                    # placed on the right by twinx above.
                    rhax[i].spines["right"].set_position(("axes", 1.1))
                    make_one_spine_visible(rhax[i], side='right')


                    for key, val in RH_dict.items():
                        rhax[i]._get_lines.prop_cycler = twinax[i]._get_lines.prop_cycler #for better autocolors
                        a = rhax[i].plot(df.index, df[val['col']], 'o', ms=0.125)
                        rhax[i].set_ylabel(val['name'], color=a[0].get_color())

                    axes[i]._get_lines.prop_cycler = twinax[i]._get_lines.prop_cycler #for better autocolors
                    #plot deltaP
                    for d_name in d_names:
                        if not((d_name in dont_mirror) and (method_list[i][0] == True)): #e.g. don't use Arduino data for mirrored data, since it won't be mirrored
                            if relP:
                                y = df[d_name]/df[P_dict['old_height']['col']]*100
                            else:
                                y = df[d_name]
                            axes[i].plot(df.index, y, 'o', ms=0.25, label=d_name)
                            axes[i].set_zorder(twinax[i].get_zorder()+1)
                            axes[i].patch.set_visible(False)
                            # print(f'plotting {d_name}...')
                        else:
                            axes[i]._get_lines.prop_cycler.__next__()




                    if i >= axs.shape[0]: #axs number of cols: 3, the loop is iterating over axs.T
                        axes[i].yaxis.set_tick_params(labelbottom=True)

                    #add text and labels
                    axes[i].set_title(method_list[i])
                    axes[i].set_ylabel(yname)
                    axes[i].axis(xmin=df.index[1], xmax=df.index[-1])
                    #nice time axes
                    axes[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                    plt.setp(axes[i].get_xticklabels(), ha="right", rotation=45)
                    axes[i].grid(True)


                if P_plot:
                    #plot P
                    for key, val in P_dict.items():
                        if not((key in dont_mirror) and (method_list[i][0] == True)) and (key in P_list):
                            axesP[i].plot(df.index, df[val['col']], 'o', fillstyle='full', ms=0.125, label=' '.join(val['name'].split(' ')[:-1]))
                            # print(f"plotting {val['col']}...")
                        else:
                            axesP[i]._get_lines.prop_cycler.__next__()

                    for dalt_name in dalt_names:
                        # print('\n\n\n\n\n--------------------------\n',df[dalt_name], sep='\n\n', end='\n\n\n\n')
                        if not((dalt_name in dont_mirror) and (method_list[i][0] == True)):
                            y=df[dalt_name]
                            daltax[i]._get_lines.prop_cycler = axesP[i]._get_lines.prop_cycler
                            daltax[i].plot(df.index, y, 'o', ms=0.25, label=dalt_name)
                        else:
                            daltax[i]._get_lines.prop_cycler.__next__()
                    daltax[i].set_ylabel(f'$\Delta$z [m]')
                    axesP[i].set_ylabel(ynameP)
                    axesP[i].set_yscale('log')
                    axesP[i].grid(True)
                    axesP[i].axis(xmin=df.index[1], xmax=df.index[-1])
                    axesP[i].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                    plt.setp(axesP[i].get_xticklabels(), ha="right", rotation=45)

                    if i >= axsP.shape[0]: #axs number of cols: 3, the loop is iterating over axs.T
                        axesP[i].yaxis.set_tick_params(labelbottom=True)


            if saveplots:
                names = ['_'.join(title.split('\n')[0].split(' ')+[special]+[name]) for name in names] #really depends on titles and special!


            for fi, axx, name in zip(figs, axes_set, names):
                if 'RH_T_z' not in name:
                    handles = []
                    labels = []
                    for ax in axx:
                        n_handles, n_labels = ax[-1].get_legend_handles_labels() #ax is a list of axes (i.e. for all subplots), only use the last axis
                        handles += n_handles
                        labels += n_labels
                        fi.delaxes(ax[0]) #remove mirrored/'combined', since it is the same as mirrored/True
                    emp = fi.add_subplot(321)
                    emp.axis('off')
                    emp.legend(handles, labels, loc='lower right', markerscale=30)
                    fi.subplots_adjust(top=0.91,
                                        bottom=0.105,
                                        left=0.055,
                                        right=0.905,
                                        hspace=0.11,
                                        wspace=0.3)
                if saveplots:
                    path = saveplots
                    fi.savefig(path/(name+'.png'), dpi=300)

    return save_df
