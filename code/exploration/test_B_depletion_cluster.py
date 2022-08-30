# !/usr/bin/env python3

"""Test nutrient A naive environment to check for correct selection coefficients etc"""

from archived.nutrient_B_depletion import *


iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_simulations(i):
    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]
    absorption = [0.01/20, 0.01/26, 0.01/35]


    LN = [[L, N, a] for N in Ns for L in Ls for a in absorption]

    L = LN[i][0]
    N = LN[i][1]
    a = LN[i][2]

    if a == absorption[0]:
        pa = 50
    elif a == absorption[1]:
        pa = 30
    else:
        pa = 15

    orbium = Creature("orbium", cluster=True)
    orbium.absorption = a
    enviro = StochasticEnviro(L=L)
    optimise(orbium, enviro, N=N, max_time=60, seed=i, name="test_depletion_B_s"+str(i)+"L"+str(L)+"N"+str(N)+"absorption"+str(pa))

run_simulations(iter)
