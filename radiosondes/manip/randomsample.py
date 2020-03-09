# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 09:25:30 2020

@author: Thomas Wagenhäuser

module name: randomsample
purpose: (pseudo-)randomly select blocks of data, e.g. to remove it from a datasets

"""
import numpy as np
import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
import random
from random import randrange


def randomblox(df, nstarts=15, blocklim=400, random_state=None):
    """Function to randomly select multiple blocks of a dataframe.

    df: pandas.DataFrame
    nstarts: number of random block starting points
    blocklim: max size of randomly generated blocks
    random_state: set seed for reproducability
    returns: index"""

    ranidx = df.sample(nstarts, random_state=random_state).index
    ranblox = ranidx
    for idx in ranidx:
        random.seed(random_state)
        blocksize = randrange(blocklim)
        s = df.index[:-blocksize].searchsorted(idx)
        ranblox = ranblox.append(df.iloc[s:s+blocksize-1].index)
    return ranblox


def set_nan_df_idx(df, idx):
    """Function to set rows to numpy.nan values.

    df: pandas.DataFrames
    idx: indices in df to set to nan"""
    
    df.loc[df.index.isin(idx),:] = np.nan
    return df


def random_nan_blox(df, nstarts=15, blocklim=400, random_state=None):
    """Function to randomly set rows to numpy.nan values.

    df: pandas.DataFrame
    nstarts: number of random block starting points
    blocklim: max size of randomly generated blocks
    random_state: set seed for reproducability
    returns: pandas.DataFrame"""

    rndm_blx = randomblox(df=df, nstarts=nstarts, blocklim=blocklim, random_state=random_state)
    return set_nan_df_idx(df, rndm_blx)
