# !/usr/bin/env python3

"""Test nutrient A naive environment to check for correct selection coefficients etc"""

import os
from two_nutrient_B import *


iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]
    absorption = [0.00008, 0.0001, 0.00012]

    LN = [[L, N, a] for N in Ns for L in Ls for a in absorption]

    L = LN[i][0]
    N = LN[i][1]
    a = LN[i][2]

    orbium = Creature("orbium", cluster=True)
    orbium.absorption = a
    enviro = StochasticEnviro(L=L)
    optimise(orbium, enviro, N=N, max_time=60, seed=i, name="test_multiple_B_s"+str(i)+"L"+str(L)+"N"+str(N)+"absorption"+str(a))

run_simulations(iter)
