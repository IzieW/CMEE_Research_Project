# !/usr/bin/env python3

"""Find best starting nutrient and depletion values for alt_nutrient2

Orbium's nutrition is depleted by 0.01 each timestep"""

from archived.alt_nutrient import *


# With depletion of 0.01 every timestep
def plot_nutrition(absorption = 0.00003):
    """Plot nutrition across different environments"""
    orbium = Creature("orbium", absorption = absorption)
    enviros = [StochasticEnviro(L=L) for L in [10, 15, 20, 30]]
    Ls = [10, 15, 20, 30]
    for i in range(len(enviros)):
        means = run_one(orbium, enviros[i], return_means=True)
        orbium.initiate()
        plt.plot(means, label="L="+str(Ls[i]))
    plt.legend()
    plt.xlabel("Timesteps")
    plt.ylabel("Nutrition")
    plt.title("Orbium nutrition with depletion of 0.01 and nutrient of 3e-5")
    plt.show()

def test_means(runs = 10, T=1000, L=20):
    """Test how survival time changes over different runs"""
    orbium = Creature("orbium")
    enviro = StochasticEnviro(L=L, T=T)
    times =[]
    all_means = []
    for i in range(runs):
        orbium.initiate()
        enviro.initiate()
        enviro.update()
        times.append(run_one(orbium, enviro))

    return times, all_means
