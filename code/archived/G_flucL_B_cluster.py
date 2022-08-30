# !/usr/bin/env python3

"""Cluster evolve orbium in stochastic attractant distribution made up of
nutrient B, where L changes every g generatons.

Run over environments with T=[500, 1000, 5000, 10000]
g = [10, 100, 500], ie. generation interval before L updates

N = [100, 1000, 10000]

For both L that fluctates randomely and non-randomely.


ie. 720 RUNS TOTAL """

from stochastic_nutrient_G_flucL_B import *
import os

## FUNCTIONS ##

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    # Initiate different environments and generation inverals
    Tg = [[t, g] for t in [500, 1000, 5000, 10000] for g in [10, 100, 500]]  # all T g combos
    Tg = Tg + Tg + Tg

    # For half of runs, fluctate L randomely
    if i < 360:
        random_fluc = False
        fluc = "stable"
    else:
        random_fluc = True
        fluc = "random"
        i -= 360  # reset

    ind = i // 10  # Run 10 seeds for each- ie. update T every multiple of 10

    T = Tg[ind][0]
    g = Tg[ind][1]

    # Initiate orbium
    orbium = Creature("orbium", cluster=True)
    orbium.absorption = 0.00003

    # Initiate Environment
    enviro = StochasticEnviro(T=T, L=15, l=64)

    # Specify run parameters
    max_time = 60*70  # max time to run simulation
    fixation = 600  # solution stays fixed over fixation number generations before optimisation succeeds

    if i < 120: # Split into three and run on different population sizes
        N=100
        name = "G_flucL_B_s"+str(i)+"l64L"+fluc+"g"+str(g)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, g=g, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    elif i < 240:
        N=1000
        name = "G_flucL_B_s"+str(i)+"l64L"+fluc+"g"+str(g)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, g=g, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted
    else:
        N=10000
        name = "G_flucL_B_s"+str(i)+"l64L"+fluc+"g"+str(g)+"T"+str(T)+"fixation"+str(fixation)+"N"+str(N)
        optimise(orbium, enviro, g=g, fluc_random= random_fluc, N=N, k=1, seed=i, fixation_mark=fixation, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted



run_simulations(iter)
