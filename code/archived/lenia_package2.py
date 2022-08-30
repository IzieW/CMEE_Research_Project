# !/usr/bin/env python3

"""Script to organise Lenia configurations into cleaner classes.
First work towards developing a package"""

# IMPORTS #
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import csv
from copy import deepcopy
from matplotlib import animation
from scipy.signal import convolve2d
import scipy as sc
from IPython.display import HTML, Image
import IPython
import multiprocessing as mp
import time
#import statsmodels.tsa.stattools as t
# Silence warnings
np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)  # Silence


############ INITIATE CLASSES  ####################
class Creature:
    """Defines a life form, their kernel and growth functions
    R: radius
    T: Time
    m: mean (growth function)
    s: standard deviation (growth function)
    b: kernel peaks (growth function)"""
    keys = ["R", "T", "m", "s", "b"]
    size = 64
    mid = size // 2

    bell = lambda x, m, s: np.exp(-((x - m) / s) ** 2 / 2)

    species_cells = {"orbium": [[0, 0, 0, 0, 0, 0, 0.1, 0.14, 0.1, 0, 0, 0.03, 0.03, 0, 0, 0.3, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0.08, 0.24, 0.3, 0.3, 0.18, 0.14, 0.15, 0.16, 0.15, 0.09, 0.2, 0, 0, 0,
                                 0],
                                [0, 0, 0, 0, 0, 0.15, 0.34, 0.44, 0.46, 0.38, 0.18, 0.14, 0.11, 0.13, 0.19, 0.18, 0.45,
                                 0, 0, 0],
                                [0, 0, 0, 0, 0.06, 0.13, 0.39, 0.5, 0.5, 0.37, 0.06, 0, 0, 0, 0.02, 0.16, 0.68, 0, 0,
                                 0],
                                [0, 0, 0, 0.11, 0.17, 0.17, 0.33, 0.4, 0.38, 0.28, 0.14, 0, 0, 0, 0, 0, 0.18, 0.42, 0,
                                 0],
                                [0, 0, 0.09, 0.18, 0.13, 0.06, 0.08, 0.26, 0.32, 0.32, 0.27, 0, 0, 0, 0, 0, 0, 0.82, 0,
                                 0],
                                [0.27, 0, 0.16, 0.12, 0, 0, 0, 0.25, 0.38, 0.44, 0.45, 0.34, 0, 0, 0, 0, 0, 0.22, 0.17,
                                 0],
                                [0, 0.07, 0.2, 0.02, 0, 0, 0, 0.31, 0.48, 0.57, 0.6, 0.57, 0, 0, 0, 0, 0, 0, 0.49, 0],
                                [0, 0.59, 0.19, 0, 0, 0, 0, 0.2, 0.57, 0.69, 0.76, 0.76, 0.49, 0, 0, 0, 0, 0, 0.36, 0],
                                [0, 0.58, 0.19, 0, 0, 0, 0, 0, 0.67, 0.83, 0.9, 0.92, 0.87, 0.12, 0, 0, 0, 0, 0.22,
                                 0.07],
                                [0, 0, 0.46, 0, 0, 0, 0, 0, 0.7, 0.93, 1, 1, 1, 0.61, 0, 0, 0, 0, 0.18, 0.11],
                                [0, 0, 0.82, 0, 0, 0, 0, 0, 0.47, 1, 1, 0.98, 1, 0.96, 0.27, 0, 0, 0, 0.19, 0.1],
                                [0, 0, 0.46, 0, 0, 0, 0, 0, 0.25, 1, 1, 0.84, 0.92, 0.97, 0.54, 0.14, 0.04, 0.1, 0.21,
                                 0.05],
                                [0, 0, 0, 0.4, 0, 0, 0, 0, 0.09, 0.8, 1, 0.82, 0.8, 0.85, 0.63, 0.31, 0.18, 0.19, 0.2,
                                 0.01],
                                [0, 0, 0, 0.36, 0.1, 0, 0, 0, 0.05, 0.54, 0.86, 0.79, 0.74, 0.72, 0.6, 0.39, 0.28, 0.24,
                                 0.13, 0],
                                [0, 0, 0, 0.01, 0.3, 0.07, 0, 0, 0.08, 0.36, 0.64, 0.7, 0.64, 0.6, 0.51, 0.39, 0.29,
                                 0.19, 0.04, 0],
                                [0, 0, 0, 0, 0.1, 0.24, 0.14, 0.1, 0.15, 0.29, 0.45, 0.53, 0.52, 0.46, 0.4, 0.31, 0.21,
                                 0.08, 0, 0],
                                [0, 0, 0, 0, 0, 0.08, 0.21, 0.21, 0.22, 0.29, 0.36, 0.39, 0.37, 0.33, 0.26, 0.18, 0.09,
                                 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0.03, 0.13, 0.19, 0.22, 0.24, 0.24, 0.23, 0.18, 0.13, 0.05, 0, 0, 0,
                                 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.06, 0.08, 0.09, 0.07, 0.05, 0.01, 0, 0, 0, 0, 0]],
                     "Hydrogeminium": [
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.03, 0.04, 0.04, 0.04, 0.03, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0.04, 0.1, 0.16, 0.2, 0.23, 0.25, 0.24, 0.21, 0.18, 0.14, 0.1, 0.07, 0.03, 0, 0,
                          0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0.01, 0.09, 0.2, 0.33, 0.44, 0.52, 0.56, 0.58, 0.55, 0.51, 0.44, 0.37, 0.3, 0.23, 0.16,
                          0.08, 0.01, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0.13, 0.29, 0.45, 0.6, 0.75, 0.85, 0.9, 0.91, 0.88, 0.82, 0.74, 0.64, 0.55, 0.46, 0.36,
                          0.25, 0.12, 0.03, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0.14, 0.38, 0.6, 0.78, 0.93, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99, 0.89, 0.78, 0.67, 0.56,
                          0.44, 0.3, 0.15, 0.04, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.08, 0.39, 0.74, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.85, 0.74,
                          0.62, 0.47, 0.3, 0.14, 0.03, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.32, 0.76, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.88, 0.75,
                          0.61, 0.45, 0.27, 0.11, 0.01, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.35, 0.83, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.88, 0.73, 0.57, 0.38, 0.19, 0.05, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 0.85, 0.67, 0.47, 0.27, 0.11, 0.01],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.55, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.93, 0.83, 0.79, 0.84, 0.88, 0.89, 0.9, 0.93, 0.98, 1.0, 1.0,
                          1.0, 1.0, 0.98, 0.79, 0.57, 0.34, 0.15, 0.03],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.47, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.72, 0.54, 0.44, 0.48, 0.6, 0.7, 0.76, 0.82, 0.91, 0.99,
                          1.0, 1.0, 1.0, 1.0, 0.91, 0.67, 0.41, 0.19, 0.05],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.27, 0.99, 1.0,
                          1.0, 1.0, 1.0, 0.9, 0.71, 0.65, 0.55, 0.38, 0.2, 0.14, 0.21, 0.36, 0.52, 0.64, 0.73, 0.84,
                          0.95, 1.0, 1.0, 1.0, 1.0, 1.0, 0.78, 0.49, 0.24, 0.07],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.14, 0.63, 0.96,
                          1.0, 1.0, 1.0, 0.84, 0.17, 0, 0, 0, 0, 0, 0, 0, 0.13, 0.35, 0.51, 0.64, 0.77, 0.91, 0.99, 1.0,
                          1.0, 1.0, 1.0, 0.88, 0.58, 0.29, 0.09],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.07, 0.38, 0.72,
                          0.95, 1.0, 1.0, 1.0, 0.22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.11, 0.33, 0.5, 0.67, 0.86, 0.99, 1.0,
                          1.0, 1.0, 1.0, 0.95, 0.64, 0.33, 0.1],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.32, 0.49, 0.71, 0.93,
                          1.0, 1.0, 1.0, 0.56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.31, 0.52, 0.79, 0.98, 1.0, 1.0,
                          1.0, 1.0, 0.98, 0.67, 0.35, 0.11],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.6, 0.83, 0.98, 1.0,
                          1.0, 0.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0.38, 0.71, 0.97, 1.0, 1.0, 1.0,
                          1.0, 0.97, 0.67, 0.35, 0.11],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.51, 0.96, 1.0, 1.0,
                          0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09, 0.34, 0.68, 0.95, 1.0, 1.0, 1.0, 1.0,
                          0.91, 0.61, 0.32, 0.1],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.13, 0.56, 0.99, 1.0, 1.0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.17, 0.45, 0.76, 0.96, 1.0, 1.0, 1.0, 1.0,
                          0.82, 0.52, 0.26, 0.07],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.33, 0.7, 0.94, 1.0, 1.0,
                          0.44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.33, 0.68, 0.91, 0.99, 1.0, 1.0, 1.0,
                          1.0, 0.71, 0.42, 0.19, 0.03],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.53, 0.89, 1.0, 1.0, 1.0, 0.8,
                          0.43, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.47, 0.86, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.95, 0.58, 0.32, 0.12, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.77, 0.99, 1.0, 0.97, 0.58, 0.41,
                          0.33, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.54, 0.95, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.8, 0.44, 0.21, 0.06, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.39, 0.83, 1.0, 1.0, 0.55, 0.11,
                          0.05, 0.15, 0.22, 0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.58, 0.99, 1.0, 1.0, 1.0,
                          1.0, 1.0, 0.59, 0.29, 0.11, 0.01, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.04, 0.55, 0.81, 0.86, 0.97, 1.0, 1.0, 0.5,
                          0, 0, 0.01, 0.09, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.26, 0.78, 1.0, 1.0, 1.0, 1.0,
                          1.0, 0.66, 0.35, 0.13, 0.03, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.33, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.93,
                          0.11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.23, 0.73, 0.95, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.62, 0.35, 0.12, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.51, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.72,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.56, 0.25,
                          0.09, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12, 0.38, 1.0, 1.0, 1.0, 0.66, 0.08, 0.55, 1.0,
                          1.0, 1.0, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.35, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          0.67, 0.12, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.49, 0, 0, 0.87, 1.0,
                          0.88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7, 0.07, 0, 0,
                          0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.04, 0.21, 0.48, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0, 0,
                          0.04, 0.42, 0.26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12, 0.21, 0.34, 0.58, 1.0, 1.0, 1.0, 0.99,
                          0.97, 0.99, 0.46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.96, 0, 0.31, 1.0, 1.0, 1.0, 0.53, 0, 0,
                          0, 0, 0, 0, 0, 0, 0.2, 0.21, 0, 0, 0, 0.27, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.87, 0.52, 0.01, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0.84, 1.0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0, 0.83, 1.0, 1.0, 0.52, 0, 0, 0,
                          0, 0, 0, 0, 0.26, 0.82, 0.59, 0.02, 0, 0, 0.46, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9, 0.55, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0.39, 0.99, 1.0, 1.0, 1.0, 1.0, 0.78, 0.04, 0, 0, 0, 0.93, 0.92, 0, 0, 0,
                          0, 0, 0, 0, 0, 0.69, 1.0, 1.0, 0.36, 0, 0, 1.0, 1.0, 0.65, 0.66, 0.97, 0.87, 0.54, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0.55, 0.75, 0.59, 0.74, 1.0, 1.0, 0, 0, 0.75, 0.71, 0.18, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0.29, 0, 0, 0.45, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.47, 0.39, 0.71, 0.25, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0.69, 0.81, 0.8, 0.92, 1.0, 0.13, 0, 0, 0.13, 0.94, 0.58, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 1.0, 1.0, 0.34, 0, 0.04, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.24, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0.63, 0.85, 0.9, 0.98, 1.0, 0.09, 0, 0, 0.02, 1.0, 0.64, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.59, 1.0, 1.0, 0.84, 0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0.64, 0.65, 0.67, 1.0, 1.0, 0.21, 0.01, 0, 0.04, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.69, 1.0, 1.0, 1.0, 0.29, 0.37, 1.0, 1.0, 0.6, 0.63, 1.0, 0.84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0.44, 0.73, 0.73, 0.85, 1.0, 0.97, 0.23, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0.06, 0, 0,
                          0, 0.97, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.33, 0.24, 0.67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0.12, 0.55, 0.9, 0.9, 1.0, 1.0, 1.0, 0.43, 0.04, 0, 0, 0, 0, 0, 0, 0, 0.31, 0.54, 0,
                          0, 0, 0.88, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0.29, 0.71, 1.0, 1.0, 1.0, 1.0, 0.79, 0.28, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0.77, 0.54,
                          0, 0, 0.87, 1.0, 1.0, 1.0, 1.0, 1.0, 0.31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0],
                         [0, 0.16, 0.27, 0.41, 0.72, 0.99, 1.0, 1.0, 0.82, 0.42, 0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0.1,
                          0.55, 0.58, 0.58, 0.77, 0.99, 1.0, 1.0, 1.0, 1.0, 0.63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0.31, 0.48, 0.45, 0.46, 0.63, 0.88, 1.0, 0.83, 0.59, 0.28, 0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0.32, 0.7, 0.95, 1.0, 1.0, 1.0, 1.0, 0.7, 0.58, 0.12, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0.23, 0.54, 0.53, 0.48, 0.57, 0.59, 0.65, 0.63, 0.55, 0.35, 0.13, 0.03, 0.02, 0.09, 0.74, 1.0,
                          0.09, 0, 0, 0, 0.32, 0.86, 1.0, 1.0, 1.0, 1.0, 0.57, 0.44, 0.31, 0.16, 0.01, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0.31, 0.45, 0.31, 0.18, 0.28, 0.39, 0.47, 0.54, 0.5, 0.35, 0.2, 0.16, 0.28, 0.75, 1.0,
                          0.42, 0.01, 0, 0, 0.6, 1.0, 1.0, 1.0, 1.0, 0.51, 0.29, 0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0.14, 0.3, 0.4, 0.54, 0.71, 0.74, 0.65, 0.49, 0.35, 0.27, 0.47, 0.6, 0.6, 0.72,
                          0.98, 1.0, 1.0, 1.0, 1.0, 0.65, 0.33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0.06, 0.33, 0.53, 0.69, 0.94, 0.99, 1.0, 0.84, 0.41, 0.16, 0.15, 0.96, 1.0,
                          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.73, 0.13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0.42, 0.86, 0.98, 0.98, 0.99, 1.0, 0.94, 0.63, 0.32, 0.62, 1.0, 1.0, 1.0,
                          1.0, 1.0, 1.0, 1.0, 0.65, 0.23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0.07, 0.62, 0.95, 1.0, 1.0, 0.99, 0.98, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          1.0, 0.98, 0.14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0.03, 0.46, 0.89, 1.0, 1.0, 0.97, 0.83, 0.75, 0.81, 0.94, 1.0, 1.0, 1.0, 1.0,
                          0.99, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0.14, 0.57, 0.88, 0.93, 0.81, 0.58, 0.45, 0.48, 0.64, 0.86, 0.97, 0.99,
                          0.99, 0.42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0.23, 0.45, 0.47, 0.39, 0.29, 0.19, 0.2, 0.46, 0.28, 0.03, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08, 0.22, 0.24, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.07, 0.22, 0.14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                     "random": np.random.random_sample([64,64]),
                     "Not": np.zeros([64,64])}

    def __init__(self, filename, dict=0, species="orbium", cluster=False, scale = None, random=False, cx=20, cy=20, dir=0, n=1, absorption=0.00003):
        """Initiate creature from parameters filename, or if file is false, load dictionary"""
        if filename:
            dict = {"organism_count": None, "absorption": None, "kernel_m": None, "correlation_length": None, "L": None, "absorption": None, "fix": None, "generations":None}
            name = deepcopy(filename)
            # Load parameters #
            if cluster:
                filename = filename + "_parameters.csv"
            else:
                filename = "../parameters/" + filename.lower() + "_parameters.csv"
            with open(filename, "r") as f:
                csvread = csv.reader(f)
                for row in csvread:
                    if row[0] == "b":  # Where b is list of values
                        dict[row[0]] = [float(i) for i in row[1].strip("[]").split(",")]
                    else:
                        dict[row[0]] = float(row[1])
            self.name = name
        else:
            self.name = dict["orbium"]
        self.cells = Creature.species_cells[species]  # Cells only
        # Each parameter
        self.R = dict["R"]
        self.T = dict["T"]
        self.m = dict["m"]
        self.s = dict["s"]
        self.b = np.asarray(dict["b"])
        self.dict = dict

        self.nutrition = 1

        ## records ##
        self.time_log = np.array([["mean", "var"]])
        self.theta_log = np.array([Creature.keys])
        self.selection_log = []

        if scale:
            self.cells = sc.ndimage.zoom(self.cells, scale, order=0); self.R *= scale
            self.scale = scale

        if not dict["absorption"]:
            dict["absorption"] = 0.0000003
        if absorption:
            dict["absorption"] = absorption
        self.absorption = dict["absorption"]

        if dict["kernel_m"]:
            self.kernel_m = float(dict["kernel_m"])
            self.kernel_s = float(dict["kernel_s"])
        else:
            self.kernel_m = 0.5
            self.kernel_s = 0.15

        self.running_mean = deepcopy(dict["m"])

        if dict["correlation_length"]:
            self.enviro_L = int(dict["correlation_length"])
            self.enviro_l = int(dict["l"])
            self.enviro_T = int(dict["correlation_time"])
            self.enviro = StochasticEnviro(l=self.enviro_l, L=self.enviro_L, T=self.enviro_T)

        self.cx = cx
        self.cy = cy

        self.dict = dict
        self.generations = dict["generations"]


        if dict["organism_count"]:
            self.n = int(dict["organism_count"])
        else:
            self.n = int(n)

        if dict["fix"]:
            self.fix = dict["fix"]

        self.mutations = 0
        self.evolved_in = 0

        if dir:
            for i in range(dir):
                self.cells = np.rot90(self.cells)
        self.dir = dir

        self.survival_mean = dict["survival_mean"]  # Mean survival time in evolved environment
        self.survival_var = dict["survival_var"]  # Survival var in evolved environment

        self.A = 0
        self.K = self.kernel()

        self.initiate(random=random)  # load A

    def figure_world(A, cmap="viridis"):
        """Set up basic graphics of unpopulated, unsized world"""
        global img  # make final image global
        fig = plt.figure()  # Initiate figure
        img = plt.imshow(A, cmap=cmap, interpolation="nearest", vmin=0)  # Set image
        plt.title = ("World A")
        plt.close()
        return fig

    def figure_asset(self, cmap="viridis", K_sum=1, bar_K=False):
        """ Chart creature's kernel and growth based on parameter solutions
        Subplot 1: Graph of Kernel in matrix form
        Subplot 2: Cross section of Kernel around center. Y: gives values of cell in row, X: gives column number
        Subplot 3: Growth function according to values of U (Y: growth value, X: values in U)
        """
        R = self.R
        K = self.kernel(fourier=False)
        growth = self.growth
        K_size = K.shape[0];
        K_mid = K_size // 2  # Get size and middle of Kernel
        fig, ax = plt.subplots(1, 3, figsize=(14, 2),
                               gridspec_kw={"width_ratios": [1, 1, 2]})  # Initiate figures with subplots

        ax[0].imshow(K, cmap=cmap, interpolation="nearest", vmin=0)
        ax[0].title.set_text("Kernel_K")

        if bar_K:
            ax[1].bar(range(K_size), K[K_mid, :], width=1)  # make bar plot
        else:
            ax[1].plot(range(K_size), K[K_mid, :])  # otherwise, plot normally
        ax[1].title.set_text("K cross-section")
        ax[1].set_xlim([K_mid - R - 3, K_mid + R + 3])

        if K_sum <= 1:
            x = np.linspace(0, K_sum, 1000)
            ax[2].plot(x, growth(x))
        else:
            x = np.arange(K_sum + 1)
            ax[2].step(x, growth(x))
        ax[2].axhline(y=0, color="grey", linestyle="dotted")
        ax[2].title.set_text("Growth G")
        return fig

    def save(self, verbose=True, cluster=False):
        """Save creature configuration to csv"""
        if cluster:
            filename = self.name.lower() + "_parameters.csv"
        else:
            filename = "../parameters/" + self.name.lower() + "_parameters.csv"
        with open(filename, "w") as f:
            csvwrite = csv.writer(f)
            for i in Creature.keys:
                csvwrite.writerow([i, self.__dict__[i]])
            csvwrite.writerow(["mutations", self.mutations])
            csvwrite.writerow(["gradient", self.evolved_in])
            csvwrite.writerow(["survival_mean", self.survival_mean])
            csvwrite.writerow(["survival_var", self.survival_var])
            csvwrite.writerow(["organism_count", self.n])
        if verbose:
            print(self.name + " configuration saved to parameters/")

    def initiate(self, size=size, random=False):
        """Initiate learning channel with creature cell configurations"""
        A = np.zeros([size, size])
        new_A = deepcopy(A)
        if random:  # if multiple orbium
            temp = deepcopy(A)
            temp[np.random.randint(64), np.random.randint(64)] = 1  # seed grid randomely for each n of orbium
            c = self.cells
            for i in range(np.random.randint(4)):  # Rotate cells randomely
                c = np.rot90(c)
            self.A = convolve2d(temp, c, mode="same", boundary="wrap")
        else:
            A[self.cx, self.cy] = 1
            A = convolve2d(A, self.cells, mode="same", boundary="wrap")  # Update grid
            self.A = A

        self.kernel()

        self.nutrition = 1  # reset nutrition
        self.running_mean = deepcopy(self.m)  # reset to starting mean value

    def kernel(self, mid=mid, fourier=True, show=False):
        """ Learning kernel for parameter solution. Default fourier transformed"""
        D = np.linalg.norm(np.ogrid[-mid:mid, -mid:mid]) / self.R * len(self.b)  # define distance matrix
        """Take all cells within distance 1 and transform along gaussian gradient. 
        Produces a smooth ring-shaped kernel"""
        K = (D < len(self.b)) * self.b[np.minimum(D.astype(int), len(self.b) - 1)] * Creature.bell(D % 1, self.kernel_m,
                                                                                                   self.kernel_s)
        K = K / np.sum(K)  # normalise
        if show:
            plt.matshow(K)
        if fourier:
            K = np.fft.fft2(np.fft.fftshift(K))  # fourier transform kernel
        return K

    def growth(self, U):
        """Defines growth of cells based input neighbourhood sum, U
        and parameter configuration. Specifically, mean and standard deviation"""
        return Creature.bell(U, self.running_mean, self.s) * 2 - 1

    def obstacle_growth(self):
        """Obstacle growth function: Obstacle creates severe negative growth in Life form"""
        return -10 * np.maximum(0, (self.enviro - 0.001))

    def update(self, i):
        """Update creature according to any number of layered environmental grids"""
        global img
        U = np.real(np.fft.ifft2(self.K * np.fft.fft2(self.A)))  # Convolve by kernel to get neighbourhood sums
        self.A = np.clip(self.A + 1 / self.T * (self.growth(U) + sum([i.growth() for i in self.enviro])), 0, 1)
        img.set_array(sum([self.A, sum([i.grid for i in self.enviro])]))
        return img,

    def update_naive(self, i):
        """Update learning channel by 1/T according to values in the learning channel"""
        #global img
        U = np.real(np.fft.ifft2(self.K * np.fft.fft2(self.A)))  # Convolve by kernel to get neighbourhood sums
        self.A = np.clip(self.A + 1 / self.T * (self.growth(U)), 0, 1)  # Update A by growth function *1/T
        #img.set_array(self.A)
        #return img

    def update_obstacle(self, i):
        """Update learning channel by 1/T according to values in the learning channel and obstacle channel"""
        global img
        U = np.real(np.fft.ifft2(self.K * np.fft.fft2(self.A)))  # Convolve by kernel to get neighbourhood sums
        self.A = np.clip(self.A + 1 / self.T * (self.growth(U) + self.obstacle_growth()), 0,
                         1)  # Update A by growth function *1/T
        img.set_array(sum([self.A, self.enviro]))
        return img,

    def render_naive(self, name):
        fig = Creature.figure_world(self.A)
        anim = animation.FuncAnimation(fig, self.update_naive, frames=200, interval = 20)
        anim.save(name, writer = "imagemagick")

    def render(self, enviro, name=None):
        """Render orbium in any number sof layered environments"""
        if name:
            name = "../results/" + name + "_anim.gif"
        else:
            name = "../results/" + self.name + "_anim.gif"
        print("Rendering animation...")
        self.initiate()
        self.enviro = enviro
        fig = Creature.figure_world(self.A + enviro.grid)
        anim = animation.FuncAnimation(fig, self.update_food, frames=200, interval=20)
        anim.save(name, writer="imagemagick")

    def render2(self, *enviro, name=None):
        """Render orbium in any number of layered environments"""
        if name:
            name = "../results/" + name + "_anim.gif"
        else:
            name = "../results/" + self.name + "_anim.gif"
        print("Rendering animation...")
        self.initiate()
        self.enviro = enviro
        fig = Creature.figure_world(self.A + sum([i.grid for i in enviro]))
        anim = animation.FuncAnimation(fig, self.update, frames=200, interval=20)
        anim.save(name, writer="imagemagick")

    def update_theta(self, muse):
        """Update parameters from parameters of input instance, muse.
        Update kernel"""
        for key in Creature.keys:
            self.__dict__[key] = muse.__dict__[key]

    def theta(self, dict=True):
        if dict:
            return {key:self.__dict__[key] for key in Creature.keys}
        else:
            return [self.__dict__[key] for key in Creature.keys]

    def show(self):
        plt.matshow(self.A)

class StochasticEnviro:
    def __init__(self, l=64, L=64, T=1000):
        self.size = l
        self.l = l  # board length
        self.L = L  # Correlation Length
        self.T = T  # Correlation time
        self.pstar = 0

        # stochastic weights
        self.xpt = 0
        self.ypt = 0

        self.Ls =[15, 20, 30, 30, 20, 15]
        self.past_L = []

        self.dt = 0  # difference in time
        self.grid = 0  # grid
        self.pseq = []  # all p from p=1 -> p*
        self.zpq = []

        self.i = 0  # counter for use in fluctating L

        self.initiate()  # initiate weights and grid

    def initiate(self):
        self.pstar = int(self.l / self.L)
        self.pseq = [[x, y] for x in range(1, self.pstar + 1) for y in range(1, self.pstar + 1)]
        self.zeta_pq()
        self.reset_weights()
        self.dt = self.T / 100  # reset dt
        self.grid = np.zeros([self.size, self.size])  # reset grid
        self.i = 0  # reset counter (for fluctating L)

    def reset_weights(self):
        """Reset stochastic weights xpt, ypt"""
        self.xpt = np.zeros(1050)
        self.ypt = np.zeros(1050)

    def zeta(self, x, p):
        return (2 * np.pi * x * p) / self.l


    def zeta_pq(self):
        """Return list of zeta matrix"""
        zpq = []
        dex = np.indices([self.size, self.size])  # array of indices dex[0]=x, dex[1]=y
        self.zpq = [self.zeta(dex[0], pq[0])+self.zeta(dex[1], pq[1]) for pq in self.pseq]


    def update_weights(self, x):
        return x * (1 - self.dt / self.T) + np.random.normal() * np.sqrt(2 * self.dt / (self.T * self.pstar))

    def get_board(self, time=1, seed=None):
        """Get entire state of grid for any amount of total time"""
        if seed:
            np.random.seed(seed)
        self.initiate()
        for t in range(0, time):
            if t % self.dt == 0:  # for every difference in time, update weights
                self.update()

    def update_L(self, random=False):
        if random:
            self.L = self.Ls[np.random.randint(0,len(self.Ls))]
        else:
            self.L = self.Ls[self.i % len(self.Ls)]
            self.i += 1
        self.pstar = int(self.l/self.L)
        self.pseq = [[x, y] for x in range(1, self.pstar + 1) for y in range(1, self.pstar + 1)]
        self.zeta_pq()
        self.reset_weights()


    def update(self):
        """Update grid by self.dt"""
        self.xpt = [self.update_weights(self.xpt[i]) for i in range(len(self.pseq))]
        self.ypt = [self.update_weights(self.ypt[i]) for i in range(len(self.pseq))]

        temp = [(self.xpt[i]+1j*self.ypt[i])*np.exp(1j*self.zpq[i]) for i in range(len(self.pseq))]

        """        for i in range(len(self.pseq)):
            temp.append((self.xpt[i]+1j*self.ypt[i])*np.exp(1j*self.zpq[i]))
            #temp.append(self.xpt[i] * np.cos(self.zpq[i]) + 1j*self.ypt[i] * np.sin(self.zpq[i]))"""

        self.grid = np.absolute(sum(temp))

    def figure_asset(self, t=1, seed=1):
        fig, ax = plt.subplots(1, 2)
        self.get_board(time=t, seed=seed)
        ax[0].imshow(self.grid)
        ax[0].title.set_text("2D")
        ax[1].plot(np.arange(0, self.l), self.grid[int(self.l / 2),], label="t0")
        for i in range(int(11 % self.dt)): self.update()
        ax[1].plot(np.arange(0, self.l), self.grid[int(self.l/2),], label = "t10")
        for i in range(int(90 % self.dt)): self.update()
        ax[1].plot(np.arange(0, self.l), self.grid[int(self.l/2),], label = "t100")
        ax[1].title.set_text("Cross section of matrix")
        plt.legend()
        plt.suptitle("L=" + str(self.L) + ", T=" + str(self.T) + ", t=" + str(t) + ", seed=" + str(seed))
        return fig

    def range_over(self, par, seed=1):
        """Plot enviro ranging over several parameters"""
        if par == "L":
            r = [8, 32, 64]
        if par == "T":
            r = [1, 10, 50, 100, 1000]
        if par == "l":
            r = [64, 100, 500]

        fig, ax = plt.subplots(3, 3)
        for x in range(5):
            self.__dict__[par] = r[x]
            self.initiate()
            self.get_board(time=1, seed=seed)
            ax[x, 0].imshow(self.grid)
            ax[x, 0].title.set_text("L=" + str(self.L) + ", T=" + str(self.T))
            ax[x, 1].plot(np.arange(0, self.size), self.grid[int(self.size / 2),])
            ax[x, 1].title.set_text("t=1")
            self.get_board(time=100, seed=seed)
            ax[x, 2].plot(np.arange(0, self.size), self.grid[int(self.size / 2),])
            ax[x, 2].title.set_text("t=1")

    def figure_world(self):
        """Configure figure for animation"""
        global img
        fig = plt.figure()
        img = plt.imshow(self.grid, cmap="viridis", interpolation="nearest", vmin=0)
        plt.close()
        return fig

    def update_anim(self, i):
        """Update for grid by t for animation"""
        global img
        if i % self.dt == 0:
            for i in range(len(self.pseq)):
                self.xpt[i] = self.update_weights(self.xpt[i])
                self.ypt[i] = self.update_weights(self.ypt[i])
            temp = []
            for i in range(len(self.pseq)):
                temp.append((self.xpt[i]+1j*self.ypt[i])*np.exp(1j*self.zpq[i]))
                #temp.append(self.xpt[i] * np.cos(self.zpq[i]) + 1j*self.ypt[i] * np.sin(self.zpq[i]))
            self.grid = np.absolute(sum(temp))

        img.set_array(self.grid)
        return img,

    def render(self, frames=200, interval=20, name=None):
        """render animation"""
        self.initiate()
        self.get_board(1)
        fig = self.figure_world()
        anim = animation.FuncAnimation(fig, self.update_anim, frames=frames, interval=interval)
        if name:
            print("rendering animation...")
            anim.save(name, writer="imagemagick")
        else:
            return anim

    def show(self):
        plt.matshow(self.grid)

class ObstacleChannel:
    def __init__(self, n=3, r=5, seed=0, dir="up", gradient=0, peak=1, wrap=True):
        """Defines obstacle environment.
        n = number of obstacles per QUARTER of grid
        r = obstacle radius
        (if moving): dir = direction of movement [up, down, left, right]
        """
        self.name = "enviro_" + str(n) + "_obstacle_radius_" + str(r)
        self.n = n
        self.r = r
        self.seed = seed
        self.gradient = gradient
        self.grid = 0
        self.peak = peak  # max value in obstacle grid
        self.wrap = wrap

        self.direction = dir
        directions = {"up": (0, 1), "down": (2, 1), "left": (1, 0), "right": (1, 2)}
        k = np.zeros([3, 3])
        k[directions[self.direction]] = 1
        self.dir_kernel = k

        mid = Creature.mid
        D = np.linalg.norm(np.ogrid[-mid:mid, -mid:mid]) / self.r

        if gradient:
            exponential = lambda x, l: l * np.exp(-l * x)
            self.kernel = ((D < 1) * exponential(D, self.gradient)) / self.gradient  # normalise to keep between 0 and 1
        else:
            self.kernel = np.ones([self.r, self.r])

        self.initiate()  # Load grid of obstacles

    def figure_asset(self):
        x = np.linspace(0, np.sum(self.kernel), 1000)
        mid = self.kernel.shape[0] // 2
        self.grid = x
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(self.kernel)
        ax[0].title.set_text("Obstacle kernel")
        ax[1].plot(self.kernel[mid, :])
        ax[1].title.set_text("Kernel cross-section")
        ax[2].plot(x, self.growth())
        ax[2].title.set_text("Obstacle growth")
        plt.show()

    def initiate(self, seed=0, size=Creature.size, mid=Creature.mid):
        """Initiate obstacle channel at random.
        Done by initiating half of grid with random obstacle configurations and
        stiching two halves together to allow more even spacing"""
        if seed:
            np.random.seed(seed)
        o = np.zeros(size * size)
        o[np.random.randint(0, len(o), self.n)] = 1
        o.shape = [size, size]
        # Convolve by kernel shape
        if self.wrap:
            self.grid = convolve2d(o, self.kernel, mode="same", boundary="wrap") * self.peak
        else:
            self.grid = convolve2d(o, self.kernel, mode="same") * self.peak

    def initiate_equal(self, seed=0):
        if seed:
            np.random.seed(seed)
        size = Creature.size // 2
        mid = size // 2
        n = self.n // 4
        o = np.zeros(size * size)
        o[np.random.randint(0, len(o), n)] = 1
        o.shape = [size, size]
        D = np.linalg.norm(np.ogrid[-mid:mid, -mid:mid]) / self.r
        k = (D < 1)
        o = convolve2d(o, k, mode="same", boundary="wrap")  # One half loaded
        o = np.vstack((o, o))
        return np.hstack((o, o))

    def growth(self):
        return -10 * np.maximum(0, (self.grid - 0.001))

    def move(self):
        self.grid = convolve2d(self.grid, self.dir_kernel, mode="same", boundary="wrap")

    def change_dir(self, direction):
        self.direction = direction
        directions = {"up": (0, 1), "down": (2, 1), "left": (1, 0), "right": (1, 2)}
        k = np.zeros([3, 3])
        k[directions[self.direction]] = 1
        self.dir_kernel = k

    def show(self):
        """Show obstacle configuration"""
        plt.matshow(self.grid)


############## INITIATE MEASUREMENTS ####################
def overlap(creature, feature):
    """Return sum of cells overlapping overlap"""
    overlapping = sum([creature.A, feature.grid]) > feature.grid
    return np.sum(overlapping * feature.grid)

global time_log, theta_log, selection_coef
time_log = pd.DataFrame(columns=["wild_mean", "wild_var"])
theta_log = pd.DataFrame(columns=Creature.keys)
selection_coef = []

def record_time(wild_mean, wild_var):
    """Record timeline"""
    global time_log
    x = pd.DataFrame([[wild_mean, wild_var]],
                     columns=["wild_mean", "wild_var"])  # record averages
    time_log = pd.concat([time_log, x])


def record_theta(creature):
    global theta_log
    theta_log = pd.concat([theta_log, pd.DataFrame(creature.theta())])


##################### EVOLUTION #############################
def mutate(p):
    """Mutate input parameter p"""
    return np.exp(np.log(p) + np.random.uniform(low=-0.2, high=0.2))


def prob_fixation(wild_time, mutant_time, N):
    """Return probability of fixation given time survived by mutant and wild type,
    and psuedo-population size N"""
    s = (mutant_time - wild_time) / wild_time  # selection coefficient

    """If s is zero, there is no selective difference between types and probability of 
    fixation is equal to 1/populdation size"""
    if s:
        return (1 - np.exp(-2 * s)) / (1 - np.exp(-2 * N * s))
    else:
        return 1 / N


def update_man(creature, enviro, moving=False, give_sums=False):
    """Update learning channel by 1/T according to values in learning channel A,
    and obstacle channel O"""
    U = np.real(np.fft.ifft2(creature.K * np.fft.fft2(creature.A)))
    if enviro:
        creature.A = np.clip(creature.A + 1 / creature.T * (creature.growth(U) + sum([e.growth() for e in enviro])), 0,
                             1)
    else:
        creature.A = np.clip(creature.A + 1 / creature.T * (creature.growth(U)), 0, 1)
    if give_sums:
        print(np.sum(creature.A))


def selection(t_wild, t_mutant):
    """Return winning solution based on survival times of wild type (t_wild) and
    mutant type (t_mutant)."""
    pfix = prob_fixation(wild_time=t_wild, mutant_time=t_mutant, N=population_size)  # get probability of fixation
    if pfix >= np.random.uniform(0, 1):
        # ACCEPT MUTATION
        return True
    else:
        # REJECT MUTATION
        return False


def run_one(creature, enviro, show_after=0, moving=False, verbose=True, give_sums=False):
    """Run creature of given parameters in given obstacle configuration until it dies.
    Show after specifies number of timesteps at when it will show what the grid looks like"""
    t = 0  # set timer
    global sums
    sums = np.zeros(10000)
    while np.sum(creature.A) and (
            t < 10000):  # While there are still cells in the learning channel, and timer is below cut off
        t += 1  # update timer by 1
        if verbose & (t % 1000 == 0):
            print(t)  # Show that it is working even after long waits
        if give_sums:
            sums[t - 1] = update_man(creature, enviro, moving=moving, give_sums=True)
        else:
            update_man(creature, enviro, moving=moving)  # Run update and show
        # if t == show_after:
        #   plt.matshow(sum([creature.A, obstacle.grid]))
    return t


def mutate_and_select(creature, enviro, moving=False, runs=100):
    """Mutate one parameter from creature and assess fitness of new solution agaisnt wild type
    in input obstacle environment. Save winning parameters to Creature.

    Method involve running wild type and mutant over ten distinct obstacle environment, and
    summing the survival time of each."""
    wild_type = creature
    mutant = deepcopy(creature)

    ## Choose parameter at random and mutate in mutant_type
    x = np.random.randint(0, len(Creature.keys) - 4)
    mutant.__dict__[Creature.keys[x]] = mutate(mutant.__dict__[Creature.keys[x]])
    mutant.K = mutant.kernel()  # update mutant kernel

    # Run mutant and wild over runs number of obstacle configurations
    t_wild = np.zeros(runs)
    t_mutant = np.zeros(runs)
    for i in range(runs):
        for e in enviro: e.initiate()  # configure environments
        wild_type.initiate()
        mutant.initiate()
        t_wild[i] = run_one(wild_type, enviro, moving=moving)
        t_mutant[i] = run_one(mutant, enviro, moving=moving)

    # Record mean and variance of survival times
    wild_mean = t_wild.mean()
    print(wild_mean)
    mutant_mean = t_mutant.mean()
    record_time(wild_mean=wild_mean,
                wild_var=t_wild.var(),
                mutant_mean=mutant_mean,
                mutant_var=t_mutant.var())

    # Select winning parameter
    if selection(wild_mean, mutant_mean):
        print("Accept mutation")
        creature.update_theta(mutant)  # Update creature parameters
        return True
    else:
        print("Reject mutation")
        return False


def optimise(creature, obstacle, N, seed=0, fixation=10, moving=False, gradient=False):
    """Mutate and select input creature in psuedo-population of size N
    until wild type becomes fixed over fixation number of generations"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    """Evolve until parameters become fixed over fixation number of generations"""
    fix = 0  # Initiate fixation count
    while fix < fixation:
        if mutate_and_select(creature, obstacle, moving=moving):  # Updates creature values
            fix = 0  # Mutation has been accepted, reset count
        else:
            fix += 1
            print(fix)

    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    creature.name = str(creature.n) + "_orbium" + "_f" + str(fixation) + "_s" + str(seed) + "N" + str(
        N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, obstacle, summary=True)
    creature.survival_mean = survival_time[0]
    creature.survival_var = survival_time[1]
    creature.evolved_in = obstacle.__dict__

    time_log.to_csv("../results/" + creature.name + "_times.csv")  # Save timelog to csv
    creature.save()  # save parameters
    return 1


def optimise_timely(creature, obstacle, N, seed=0, run_time=10, moving=False, cluster=False, name=None):
    """Mutate and select input creature in psuedo-population of size N
    until wild type becomes fixed over fixation number of generations"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    run_time = run_time * 60  # Translate to seconds
    """Evolve until parameters become fixed over fixation number of generations"""
    gen = 0  # time_count
    mutation = 0
    start = time.time()
    while (time.time() - start) < run_time:
        if mutate_and_select(creature, obstacle, moving=moving):  # Updates creature values
            mutation += 1
        gen += 1

    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    if name:
        creature.name = name
    else:
        creature.name = str(creature.n) + "_orbium_t" + str(run_time) + "s" + str(seed) + "N" + str(
            N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, obstacle, summary=True)
    creature.mutations = mutation
    creature.survival_mean = survival_time[0]
    creature.survival_var = survival_time[1]
    if cluster:
        time_log.to_csv(creature.name + "_times.csv")  # Save timelog to csv
    else:
        time_log.to_csv("../results/" + creature.name + "_times.csv")  # Save timelog to csv

    creature.save(cluster=cluster)  # save parameters
    return 1


def optimise_layered(creature, enviro, N, seed=0, run_time=10, moving=False, cluster=False, name=None):
    """Mutate and select input creature in psuedo-population of size N
    until wild type becomes fixed over fixation number of generations"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    if not enviro:
        print("WARNING: No environment specified.")

    run_time = run_time * 60  # Translate to seconds
    """Evolve until parameters become fixed over fixation number of generations"""
    gen = 0  # time_count
    mutation = 0
    start = time.time()
    while (time.time() - start) < run_time:
        if mutate_and_select(creature, enviro, moving=moving):  # Updates creature values
            mutation += 1
        gen += 1

    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    if name:
        creature.name = name
    else:
        creature.name = str(creature.n) + "_orbium_t" + str(run_time) + "s" + str(seed) + "N" + str(
            N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, enviro, summary=True)
    creature.mutations = mutation
    creature.survival_mean = survival_time[0]
    creature.survival_var = survival_time[1]
    if cluster:
        time_log.to_csv(creature.name + "_times.csv")  # Save timelog to csv
    else:
        time_log.to_csv("../results/" + creature.name + "_times.csv")  # Save timelog to csv
    creature.save(cluster=cluster)  # save parameters
    return 1


def get_survival_time(creature, obstacle=None, runs=10, summary=False, verbose=False):
    """Calculate average run time over seeded 10 configurations.
    Return mean and variance."""
    times = np.zeros(runs)
    if obstacle:
        for i in range(1, runs + 1):
            np.random.seed(i)
            creature.initiate()  # Reset grid
            for o in obstacle: o.initiate(seed=i)  # set obstacle
            times[i - 1] = run_one(creature, obstacle, verbose=verbose)
    else:
        for i in range(1, runs + 1):
            creature.initiate()
            times[i - 1] = run_one(creature, obstacle, verbose=verbose)

    if summary:
        return times.mean(), times.var()
    else:
        return times


