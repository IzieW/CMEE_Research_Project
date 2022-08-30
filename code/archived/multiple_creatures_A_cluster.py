# !/usr/bin/env python3

"""Cluster evolve orbium in stochastic attractant distribution made up of
nutrient A.

Run over environments with T=[500, 1000, 5000, 10000]
L = [15, 20, 30]

N = [100, 1000, 10000]

ie. 360 RUNS TOTAL

Each run contains two orbium in the same spatial and evolutionary arena which
are mutated individually each generation"""

from two_nutrient_A import *
import os

## FUNCTIONS ##

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    # Initiate different environments
    TL = [[t, l] for t in [500, 1000, 5000, 10000] for l in [15, 20, 30]]  # all combinations of T and L
    TL = TL + TL + TL
    ind = i // 10  # Run 10 seeds for each- ie. update T and L every multiple of 10

    T = TL[ind][0]
    L = TL[ind][1]

    # Initiate orbium
    orbium = Creature("orbium", cluster=True)
    orbium.absorption = 0.0000003

    # Initiate Environment
    enviro = StochasticEnviro(T=T, L=L, l=64)

    # Specify run parameters
    max_time = 60*70  # max time to run simulation
    fixation = 600  # solution stays fixed over fixation number generations before optimisation succeeds

    if i < 120: # Split into three and run on different population sizes
        N=100
        name = "multiple_creatures_A_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    elif i < 240:
        N=1000
        name = "multiple_creatures_A_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't evolve
    else:
        N=10000
        name = "multiple_creatures_A_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't evolve


run_simulations(iter)
