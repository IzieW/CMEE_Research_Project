# !/usr/bin/env python
"""Render every output from evolution of orbium in nutrient A

NAIVE EVOLUTION


4 unique renderings for each N (12 total - 0-11)


"""
from stochastic_nutrient_A import *
import os
# Group by T, group by N

def cluster_render(i):
    Ts = [500, 1000, 5000, 10000]
    Ts = Ts + Ts + Ts # one for each N
    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]

    # Assign N and T according to i
    N = Ns[i//3]
    T = Ts[i]

    # Find files with T and L
    files = os.listdir("../results/nutrient_A/naive/")
    files = [i for i in files if re.search(r"parameters.csv$", i)]
    files = [i for i in files if get_N(i)==N and get_T(i)==T]

    orbium = [Creature("../results/nutrient_A/naive/"+file.split("_parameters.csv")[0], cluster=True) for file in files]

    filler = Creature("orbium", species="Not")
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    new = []
    for L in Ls:
        temp = [i for i in orbium if i.enviro_L == L]
        while len(temp) < 10:
            temp.append(filler)
        new += temp

    orbium = deepcopy(new)


    seeds = [orbium[n].dict["seed"] for n in range(10)] # get first 10 seeds

    dim = [10, 3]

    labels = {"y": [np.arange(32, (64*dim[1])+32, 64), Ls],
              "ylab": "Correlation Length (L)",
              "xlab": "T="+str(T),
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds]
                       }

    render_grid(orbium, name="naive_A_N"+str(N)+"T"+str(T)+"rendered", dim = dim, **labels)

def do_render(T, N, path):
    Ls = [15, 20, 30]

    # Find files with T and L
    files = os.listdir(path)
    files = [i for i in files if re.search(r"parameters.csv$", i)]
    files = [i for i in files if get_N(i)==N and get_T(i)==T]

    orbium = [Creature(path+file.split("_parameters.csv")[0], cluster=True) for file in files]

    filler = Creature("orbium", species="Not")
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    new = []
    for L in Ls:
        temp = [i for i in orbium if i.enviro_L == L]
        while len(temp) < 10:
            temp.append(filler)
        new += temp

    orbium = deepcopy(new)


    seeds = [orbium[n].dict["seed"] for n in range(10)] # get first 10 seeds

    dim = [10, 3]

    labels = {"y": [np.arange(32, (64*dim[1])+32, 64), Ls],
              "ylab": "Correlation Length (L)",
              "xlab": "T="+str(T),
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds],
              "title": "N="+str(N)+", T="+str(T)
                       }

    render_grid(orbium, name=path+"/gifs/naive_A_N"+str(N)+"T"+str(T)+"rendered", dim = dim, **labels)

def render_all():
    Ts = [500, 1000, 5000, 10000]
    Ns = [100, 1000, 10000]

    [do_render(T, N, "../results/nutrient_A/naive/") for N in Ns for T in Ts]

