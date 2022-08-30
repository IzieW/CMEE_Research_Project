# !/usr/bin/env python3

"""Test nutrient A naive environment to check for correct selection coefficients etc"""

import os
from stochastic_nutrient_A import *


iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]

    LN = [[L, N] for N in Ns for L in Ls]

    L = LN[i][0]
    N = LN[i][1]

    orbium = Creature("orbium", cluster=True)
    orbium.absorption = 0.0000003
    enviro = StochasticEnviro(L=L)
    optimise(orbium, enviro, N=N, max_time=60, seed=i, name="test_A_s"+str(i)+"L"+str(L)+"N"+str(N))

run_simulations(iter)
