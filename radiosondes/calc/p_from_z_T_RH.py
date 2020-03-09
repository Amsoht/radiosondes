# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 09:25:30 2020

@author: Thomas Wagenhäuser

module name: p_from_z_T_RH
purpose: use the barometric equation and in situ measurements of P0, z0, z, T, RH to calculate pressure profile

comments:
-currently uses timeseries data from radiosondes

ToDo:
-move asc_mirror() to manip folder?
-support df without time information (only z profile)
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


#%%
#------------------------------------------------------------------------------
#                      -- Declare some constants --
#
#  Parameter         Description                 Units
#  ---------   ---------------------------     -----------
#     Rd       Gas constant for dry air        J/degree/kg
#
#     g        Acceleration due to gravity       m/s^2
#
#    eps       Ratio of the molec. weights        None
#              of water and dry air
#
#------------------------------------------------------------------------------
#  Rd  = 286.9968933
#  g   = 9.80616
#  eps = 0.621970585
#  ratio  = ( 1.0 - eps ) / eps / 1000. ; DB: factor 1000. to account for units g/kg
#
#------------------------------------------------------------------------------
#            -- Check if levels start at top or bottom of atmosphere  --
#------------------------------------------------------------------------------


class P_From_GPS_T_RH(object):
    """Integrational Pressure Calculations from Height, Temperature and Water Vapour

    make sure, GPS data is offset corrected with station height
    first value should be groundstation pressure, height, Temperature and Water Vapour
    """

    def __init__(self, df, P0=1013.15, z0=0.0, start=False, mirror=False, cols_in={}, cols_out={}): #df with certain format!
        """Input pandas.DataFrame with date_time index, 'GPS_ALT' [m], 'RH' [%], 'T' [°C], 'P' [hPa]

        column names may be specified in a dict of the following type:
            cols_in = {'z': {'col':'GPS_ALT',
                          'unit':'[m]'},
                    'RH':{'col':'RH',
                          'unit':'[%]'},
                    'T': {'col':'T',
                          'unit':'[°C]'},
                    'P': {'col':'P',
                          'unit':'[hPa]'},
                    'z_geopot':{'col':'geopotential altitude',
                                'unit':'[m]'}
                    }""" #maybe create dict rom DataFrame from Excel-csv file

        self.P0 = P0
        self.z0 = z0
        self.start = start
        self.df = df.copy()
        self.Rd = 287.05 #286.9968933 where to get accurate value?
        #self.Rw = 461.5226519
        self.g = 9.80616
        self.eps = 0.62197
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
        if mirror:
            self.asc_mirror()


    def calc_esat(self, T):
        """Calculate water vapour saturation pressure from T [°C]"""
        return 6.112*np.exp((17.67*T)/(T+243.5))


    def calc_Tv_from_w(self, w, temp):
        """Calculate Tv from water vapour mixing ratio (kg/kg). Currently not used"""
        #linear Tv with ln(press), take the mean
        Tv = temp*(1+((1-self.eps)/self.eps)*w)
        return Tv


    def calc_Tv_from_e(self, rh, e_s, press, temp, mode='std'):
        """Calculate Tv from water vapour pressure, temperature [K] and pressure [hPa]"""
        e = rh*e_s
        if mode == 'std':
            Tv = temp/(1-e/press*(1-self.eps))
        elif mode == 'Modem':
            Tv = temp*(1+(1-self.eps)*e/press)
        return Tv



    def asc_mirror(self):
        """Overwrite df with Ascent-Descent data, where descent is a mirror of ascent"""
        #create ascent DataFrame
        print('Overwriting DataFrame with mirrored ascent...')
        df = self.df.copy()
        asc_df = df[:df[self.cols_in['z']['col']].idxmax()]

        #create descent DataFrame from ascent copy
        des_df = asc_df.copy()
        idx_name = des_df.index.name
        des_df.reset_index(inplace=True) #reset to be able to sort
        des_df.sort_index(ascending=False, inplace=True, ignore_index=True)
        diff = -1 * pd.to_timedelta(des_df[idx_name].diff().dt.total_seconds().cumsum(), unit='s') #get this from dict
        des_df[idx_name] = diff + des_df.loc[0,idx_name]
        des_df.drop(0, inplace=True) #remove top record (otherwise: non-unique index)
        des_df.set_index(idx_name, inplace=True)

        #merge ascent and descent DataFrames, overwrite self.df
        self.df = pd.concat([asc_df, des_df], axis=0)
        print('Done overwriting with mirrored ascent!\n')


    def interpol_data_gaps(self, mode='combined', freq='S', expol=False):
        #TODO: account for latitude, longitude, wind, etc. they should be Nans when not original?
        #is there a difference, when using mirror and combined vs standard?
        print('Interpolating ascent and descent',
              f'mode: {mode}',
              f'freq: {freq}',
              '...',
              sep='\n')

        #resample (timebased):
        if expol: #maybe give the option to provide pressure or Temperature value
                try:
                    time, alt = expol
                    self.df.append(pd.DataFrame([alt], index=[time], columns=[self.cols_in['z']['col']]))
                except Exception as excp:
                    print(f'\n\nError: {excp}\nAn error occured while trying to extrapolate.\nContinuing without extrapolating...')

        # resam_df = self.df.resample(freq).mean() #this means, that ALL data is interpolated
        resam_df = self.df[[key['col'] for var,key in self.cols_in.items()]] \
                        .resample(freq).mean()

        if mode == 'combined':
            #combine ascent and descent met profiles and interpolate profile to timeseries
            z_df = self.df[[key['col'] for var,key in self.cols_in.items()]] \
                    .groupby(self.cols_in['z']['col']).mean()

                    # asc_z_df = asc_df.reset_index().set_index('GPS_ALT')
                    # des_z_df = des_df.reset_index().set_index('GPS_ALT')
                    # comb_z_df = df.groupby('GPS_ALT').mean()
                    #reidx_asc_df = z_df.reindex(asc_df['GPS_ALT'].values)

            #create a time and GPS_ALT dataframe that is resampled
            resam_df = self.df[[self.cols_in['z']['col']]].resample(freq).mean().interpolate()
            # print(f'\n--------\nresam_df.shape now is: {resam_df.shape}\n--------\n')

            #reindex altitude-grouped dataframe by resampled altitude values (i.e. time based)
            reidx_df = z_df.reindex(resam_df[self.cols_in['z']['col']].values) #reindex to resampled z
            reidx_df = pd.concat([reidx_df, z_df], axis=0).sort_index() #merge reindexed z_df with original z_df, in order to preserve records between resampling timestamps
            reidx_df = reidx_df[~reidx_df.index.duplicated(keep='first')] #drop duplicates, that were caused by this merge

            #fill nan values
            ipol_reidx_df = reidx_df.interpolate()

            #merge combined and interpolated data to Date_Time by altitude
            idx_name = resam_df.index.name
            ipol_resam_df = resam_df.reset_index()\
                                    .merge(ipol_reidx_df, on=self.cols_in['z']['col'], how='left', suffixes=(False,False))\
                                    .set_index(idx_name)
            #this should work. maybe interpolate until ground level?
            # print(f'\n--------\nipol_resam_df.shape now is: {ipol_resam_df.shape}\n--------\n')

        else:
            ipol_resam_df = resam_df.interpolate()

        # print(f'\n--------\nipol_resam_df.shape now is: {ipol_resam_df.shape}\n--------\n')

        print('Done interpolating!\n')
        return ipol_resam_df #overwrite previous dataframe (be careful with non-profile-data like lat and lon!) #TODO



    def calc_press(self, press, z, temp, rh, dirksenTv=True): #list press, list z, list temp, list rh, list e_s

        # if z[0] > 12000:
        #     rh = [0,0]
        #Calculate virtual temperature and pressure
        press2 = [press[1]]
        if dirksenTv:
            T_avg = np.array(temp).mean()
            e_s_avg = self.calc_esat(T_avg-273.15)
            RH_avg = np.array(rh).mean()
            P_avg = np.sqrt(press[0]*press2[0])
            Tv = [self.calc_Tv_from_e(RH_avg, e_s_avg, P_avg, T_avg)]*2
        else:
            e_s = [self.calc_esat(temp_i-273.15) for temp_i in temp]
            Tv = [self.calc_Tv_from_e(rh[0], e_s_i, press_i, temp_i) for press_i, temp_i, e_s_i in zip(press,temp, e_s)] #

        i = 0 #viele Iterationen 10, 20, 50 mal rechnen und Tv evolution anschauen
#        print(f'\nTv: = {Tv}')
        iterate = True
        while iterate and i < 20:
            if dirksenTv:
                Tv_avg = Tv[i+1]
            else:
                Tv_avg = pd.Series([Tv[i], Tv[i+1]]).mean()

            press2.append(press[0]*np.exp((z[0]-z[1])*self.g/(self.Rd*Tv_avg)))

            if dirksenTv:
                P_avg = np.sqrt(press[0]*press2[i+1])
                Tv.append(self.calc_Tv_from_e(RH_avg, e_s_avg, P_avg, T_avg))
            else:
                Tv.append(self.calc_Tv_from_e(rh[1], e_s[1], press2[i+1], temp[1]))
            i += 1
            if Tv[-1] == Tv[-2]:
                iterate = False

        n_iterations = i

        return press2[-1], n_iterations, Tv_avg


    def update_P(self, ipol='combined', freq='S', geometric=False, dirksenTv=True):
        print('Calculating pressure...')

        if self.start:
            try:
                self.df.loc[:self.start,self.cols_out['P']['col']] = self.P0
                self.df.loc[:self.start,self.cols_in['z']['col']] = self.z0
            except Exception as excp:
                print(f"An Error occured: {excp}\nCouldn't find starting time index. Proceeding with only the first record being overwritten...\n")
        if ipol:
            self.df_orig = self.df.copy()
            self.df = self.interpol_data_gaps(mode=ipol, freq=freq) #overwrites self.df

        #drop NaN rows, where T and RH is NaN
        self.df.dropna(how='all', subset=[self.cols_in[var]['col'] for var in ['T','RH','z']], inplace=True)

        self.df[self.cols_out['P']['col']] = np.full(self.df.shape[0], np.NaN) #initialize empty P column
        self.df.loc[self.df.index[0],self.cols_out['P']['col']] = self.P0 #enter pressure start value
        self.df.loc[self.df.index[0],self.cols_in['z']['col']] = self.z0 #enter altitude start value



        press = self.df[self.cols_out['P']['col']].values #.interpolate().values #TODO: load column names from a dictionary
        temp = self.df[self.cols_in['T']['col']].interpolate().values + 273.15 #in K
        rh = self.df[self.cols_in['RH']['col']].interpolate().apply(lambda x: 100 if x > 100 else x) #account for wrong RH values above 100%
        rh = rh.values / 100 #to account for %
        # e_s = pd.DataFrame(temp-273.15).iloc[:,0].apply(self.calc_esat).values #write a function to calculate this
        geometric_alt = self.df[self.cols_in['z']['col']].interpolate().values
        alt = h_t_g(geometric_alt*units.meter)/self.g
        alt = [alt[i].magnitude for i in range(len(alt))]
        self.df[self.cols_out['z_geopot']['col']] = alt
        if geometric:
            alt = geometric_alt
        n_iterations_Tv = [np.NaN]*self.df.shape[0]
        Tv = [np.NaN]*self.df.shape[0]

        if (True in np.isnan(press)):
            for i, press_i in enumerate(press[:-1]):
                #print(f'press: {press_i}\ntemp: {temp[i:i+2]}\nrh: {rh[i:i+2]}\ne_s: {e_s[i:i+2]}')
                press[i+1], n_iterations_Tv[i+1], Tv[i+1] = self.calc_press([press_i]*2, alt[i:i+2], temp[i:i+2], rh[i:i+2], dirksenTv=dirksenTv)

        else:
            for i in range(len(press)-1):
                press[i+1], n_iterations_Tv[i+1], Tv[i+1] = self.calc_press(press[i:i+2], alt[i:i+2], temp[i:i+2], rh[i:i+2], dirksenTv=dirksenTv)

        self.df[self.cols_out['P']['col']] = press #rename this? in order to keep old and new P values. check for unique values?
        self.df['n_iterations_Tv'] = n_iterations_Tv
        self.df[self.cols_out['Tv']['col']] = Tv
        if ipol:
            #add old columns to new columns
            #do I introduce some Problems here?
            newkeys = self.df.columns
            oldkeys = [col for col in self.df_orig.columns if not(col in newkeys)]
            self.df = self.df.merge(self.df_orig[oldkeys], left_index=True, right_index=True, how='left')
        #print(f'\n{self.df}')
        print('Done calculating pressure!\n')

#check units (Temperature etc.)
#add Arduino pressure readings, DONE
#add legend DONE
#add humidity DONE
#leave out humidity above 12 km? DONE 2020 02 04: better keep those values?
#mirror ascent and calculate DONE

#%% test
if __name__ == '__main__':
    pass