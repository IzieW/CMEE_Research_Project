# !/usr/bin/env python3

"""Cluster evolve orbium in stochastic attractant distribution made up of
nutrient A, where L changes every enviro.T timesteps

Run over environments with T=[500, 1000, 5000, 10000]
L = [15, 20, 30]

N = [100, 1000, 10000]

For both L that fluctates randomely and non-randomely.


ie. 240 RUNS TOTAL """

from stochastic_nutrient_T_flucL_A import *
import os

## FUNCTIONS ##

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    # Initiate different environments
    Ts = [500, 1000, 5000, 10000]  # all T
    Ts = Ts + Ts + Ts

    # For half of runs, fluctate L randomely
    if i < 120:
        random_fluc = False
        fluc = "stable"
    else:
        random_fluc = True
        fluc = "random"
        i -= 120  # reset

    ind = i // 10  # Run 10 seeds for each- ie. update T every multiple of 10

    T = Ts[ind]

    # Initiate orbium
    orbium = Creature("orbium", cluster=True)
    orbium.absorption = 0.0000003

    # Initiate Environment
    enviro = StochasticEnviro(T=T, L=15, l=64)

    # Specify run parameters
    max_time = 60*70  # max time to run simulation
    fixation = 600  # solution stays fixed over fixation number generations before optimisation succeeds

    if i < 40: # Split into three and run on different population sizes
        N=100
        name = "T_flucL_A_s"+str(i)+"l64L"+fluc+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    elif i < 80:
        N=1000
        name = "T_flucL_A_s"+str(i)+"l64L"+fluc+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    else:
        N=10000
        name = "T_flucL_A_s"+str(i)+"l64L"+fluc+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted




run_simulations(iter)
