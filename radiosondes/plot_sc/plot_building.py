# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:52:41 2020

@author: Thomas Wagenhäuser
"""
import matplotlib.pyplot as plt


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
        

def make_one_spine_visible(ax, side='right'):
    make_patch_spines_invisible(ax)
    ax.spines[side].set_visible(True)