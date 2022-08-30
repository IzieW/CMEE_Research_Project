# !/usr/bin/env python3

"""Cluster evolve TWO orbium in stochastic attractant distribution made up of
nutrient A.

Run over environments with

T=[500, 1000, 5000, 10000]

L = [10, 15, 20, 30]

N = [10, 100, 1000, 10000, 100000]

TWENTY SEEDS in each

ie. 1600 RUNS TOTAL  - qsub 0-1599 """

from two_depleting_nutrient_A import *
import os

## FUNCTIONS ##

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    # Initiate different environments
    TL = [[n, t, l] for n in [10, 100, 1000, 10000, 100000] for t in [500, 1000, 5000, 10000] for l in [10, 15, 20, 30]]  # all combinations of T and L
    ind = i // 20  # Run 10 seeds for each- ie. update T and L and N every multiple of 10

    N = TL[ind][0]
    T = TL[ind][1]
    L = TL[ind][2]

    # Initiate orbium
    orbium = Creature("orbium", cluster=True)
    orbium.absorption = (0.001/28.5)

    # Initiate Environment
    enviro = StochasticEnviro(T=T, L=L, l=64)

    # Specify run parameters
    max_time = 60*70  # max time to run simulation

    name = "two_depletion_nutrient_A_s"+str(i)+"l64L"+str(L)+"T"+str(T)+"N"+str(N)
    optimise(orbium, enviro, N=N, k=1, seed=i, max_time = max_time, name=name, minus_key=1, median=False)  # nutrient can't be adjusted




run_simulations(iter)
