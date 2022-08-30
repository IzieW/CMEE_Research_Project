# !/usr/bin/env python3

"""Cluster evolve orbium in stochastic attractant distribution made up of
nutrient B.

Run over environments with T=[500, 1000, 5000, 10000]
L = [15, 20, 30]

N = [100, 1000, 10000]

ie. 270 RUNS TOTAL

repeated over 3 different absorptions

ie. 810 runs"""

from stochastic_nutrient_B import *
import os

## FUNCTIONS ##

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    # Initiate different environments
    TL = [[t, l] for t in [100, 1000, 10000] for l in [15, 20, 30]]  # all combinations of T and L
    TL = TL + TL + TL

    # Initiate orbium
    orbium = Creature("orbium", cluster=True)

    if i < 270:
        orbium.absorption = 0.00008 # 15%
        a=15
    elif i < 540:
        orbium.absorption = 0.0001 # 30%
        i -=270
        a=30
    elif i < 810:
        orbium.absorption = 0.00012 # 50%
        i -= 540
        a=50

    ind = i // 10  # Run 10 seeds for each- ie. update T and L every multiple of 10

    T = TL[ind][0]
    L = TL[ind][1]

    # Initiate Environment
    enviro = StochasticEnviro(T=T, L=L, l=64)

    # Specify run parameters
    max_time = 60*12  # max time to run simulation
    fixation = 600  # solution stays fixed over fixation number generations before optimisation succeeds

    if i < 120: # Split into three and run on different population sizes
        N=100
        name = "naive_B_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)+"a"+str(a)
        optimise(orbium, enviro, N=N, k=1, seed=i, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    elif i < 240:
        N=1000
        name = "naive_B_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)+"a"+str(a)
        optimise(orbium, enviro, N=N, k=1, seed=i, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't evolve
    else:
        N=10000
        name = "naive_B_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)+"a"+str(a)
        optimise(orbium, enviro, N=N, k=1, seed=i, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't evolve




run_simulations(iter)
