# !/usr/bin/env python3
# Render all results from run naive_nutrient_A_render
from stochastic_nutrient_B import *

import os
#iter = int(os.environ.get("PBS_ARRAY_INDEX"))
def do_render_A(i):
    # read in data
    df = pd.read_csv("../results/naive_nutrient_A.csv")

    TL = [[n, l, t] for n in [10, 100, 1000, 10000, 100000] for t in [100, 500, 1000, 5000, 10000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]
    T = TL[i][2]
    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]
    df = df[df.correlation_time==T]



    df = df.sort_values("seed")
    names = [i.split("_parameters.csv")[0] for i in df.files]
    df["name"] = names

    dicts = [{colname: df[colname].iloc[i] for colname in df.columns} for i in range(len(df))]

    dim = [5, 4]

    # Load orbium
    orbium = [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]

    filler = Creature("orbium", species="Not")
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    if len(orbium) != 20:
        for i in range(20-len(orbium)):
            orbium.append(filler)

    return orbium
    seeds = list(df.seed.iloc[0:5])  # get first seeds

    labels = {"xlab":"top row seeds",
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds],
              "ylab":"correlation_time (T)",
              "y": [np.arange(32, (64*dim[1])+32, 64), Ts],
              "title": "N="+str(N)+", L="+str(L)
                       }

    render_grid(orbium, dim=dim, name="../results/naive_B/gifs/naive_B_N"+str(N)+"L"+str(L), **labels)

def do_render_B(i):
    # read in data
    df = pd.read_csv("../results/naive_nutrient_B.csv")

    TL = [[n, l] for n in [10, 100, 1000, 10000, 100000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]

    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]


    df = df.sort_values("seed")
    names = [i.split("_parameters.csv")[0] for i in df.files]
    df["name"] = names

    dicts = [{colname: df[colname].iloc[i] for colname in df.columns} for i in range(len(df))]

    dim = [10, 3]

    # Load orbium
    orbium = [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]

    filler = Creature("orbium", species="Not")
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    new = []
    for T in Ts:
        temp = [i for i in orbium if i.enviro_T == T]
        if len(temp) != 10:
            for i in range(10-len(temp)): temp.append(filler)  # adjust length to make sure 10 in each row
        new += temp  # add to list

    orbium = deepcopy(new)

    seeds = [orbium[i].dict["seed"] for i in range(10)]  # get first seeds

    labels = {"xlab":"top row seeds",
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds],
              "ylab":"correlation_time (T)",
              "y": [np.arange(32, (64*dim[1])+32, 64), Ts],
              "title": "N="+str(N)+", L="+str(L)
                       }

    render_grid(orbium, dim=dim, name="../results/naive_B/gifs/naive_B_N"+str(N)+"L"+str(L), **labels)

def show_max(T, parameter="seed"):
    df = pd.read_csv("../results/naive_nutrient_A.csv")
    Ls = [10, 15, 20, 30]
    d = df[df.correlation_time==T]
    orbium =[]
    for L in Ls:
            dat = d[d.correlation_length==L]
            dat = dat.sort_values("survival_mean").iloc[-5:]  # take top 5 values for each
            dat["name"]="name"
            dicts = [{colname: dat[colname].iloc[i] for colname in dat.columns} for i in range(len(dat))]
            orbium += [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]

    orbium = [i.dict[parameter] for i in orbium]
    dim = [5,4]
    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = orbium[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0
    return grid


def render_max(seed=0):
    df = pd.read_csv("../results/naive_nutrient_A.csv")
    Ts = [100, 1000, 10000]
    Ls = [10, 15, 20, 30]
    dim = [5,4]
    for T in Ts:
        d = df[df.correlation_time==T]
        orbium = []
        for L in Ls:
            dat = d[d.correlation_length==L]
            dat = dat.sort_values("survival_mean").iloc[-5:]  # take top 5 values for each
            dat["name"]="name"
            dicts = [{colname: dat[colname].iloc[i] for colname in dat.columns} for i in range(len(dat))]
            orbium += [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]


        labels = {
              "ylab":"correlation_length (L)",
              "y": [np.arange(32, (64*dim[1])+32, 64), Ls],
              "title": "max survivors T="+str(T)
                       }

        render_grid(orbium, dim=dim, name="../results/naive_A/gifs/maximums/T"+str(T), **labels)




def overview_N(T, L):
    df = pd.read_csv("../results/naive_nutrient_A.csv")

    df = df[df.correlation_time==T]
    df = df[df.correlation_length==L]
    d = df.groupby("N").describe().survival_mean
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel("Population size")
    plt.ylabel("mean survival time")

def overview_T(T):
    df = pd.read_csv("../results/naive_nutrient_A.csv")
    df = df[df.correlation_time==T]
    d = df.groupby("correlation_length").describe().survival_mean
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel("correlation length (L)")
    plt.ylabel("mean survival time")

def overview_L(L):
    df = pd.read_csv("../results/naive_nutrient_A.csv")
    df = df[df.correlation_length==L]
    d = df.groupby("correlation_time").describe().survival_mean
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel("correlation time (L)")
    plt.ylabel("mean survival time")

def show_seeds_A(T, L, N, par="seed", dim=[5, 4], dataframe=None):
    df = pd.read_csv("../results/naive_nutrient_A.csv")
    df = df[df.correlation_time==T]
    df = df[df.correlation_length==L]
    df = df[df.N==N]
    df = df.sort_values("seed")

    creatures = list(df[par])

    if len(creatures) != 20:
        for i in range(20-len(creatures)):
            creatures.append(0)

    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = creatures[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0

    df= pd.DataFrame({"seed": list(df.seed),
                         "survival": list(df.survival_mean),
                         "survival_var": list(df.survival_var)})
    df = df.sort_values("seed")

    plt.errorbar(df.seed, df.survival, yerr=np.sqrt(df.survival_var/10), capsize=4)
    plt.xlabel("seed")
    plt.ylabel("survival mean")

    if dataframe:
        return df
    else:
        return grid

def show_seeds_B(L, N, par="seed", dim=[10, 3], dataframe=None):
    df = pd.read_csv("../results/naive_nutrient_B.csv")
    df = df[df.correlation_length==L]
    df = df[df.N==N]
    df = df.sort_values("seed")

    creatures = []
    for T in [100, 1000, 10000]:
        temp = df[df.correlation_time==T]
        temp = list(temp[par])
        if len(temp) != 10:
            for i in range(10-len(temp)):
                temp.append(0)
        creatures += temp

    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = creatures[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0

    df= pd.DataFrame({"seed": list(df.seed),
                         "survival": list(df.survival_mean),
                         "survival_var": list(df.survival_var)})
    df = df.sort_values("seed")

    plt.errorbar(df.seed, df.survival, yerr=np.sqrt(df.survival_var/10), capsize=4)
    plt.xlabel("seed")
    plt.ylabel("survival mean")

    if dataframe:
        return df
    else:
        return grid


def organise(creatures, parameter="seed"):
    c = [i.dict[parameter] for i in creatures]
    dim = [5,4]
    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = c[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0

    return grid

def animate_one(seed):
    df = pd.read_csv("../results/naive_nutrient_B.csv")
    d = df[df.seed==seed]
    dict = {colname: d[colname].iloc[0] for colname in d.columns}
    dict["name"]= "this guy"
    orbium = Creature(filename=None, dict=dict)


    render_seeds(orbium, dim=[2,2], seed_range=4, name="../results/naive_B/gifs/seed"+str(seed)+".gif")

def give_one(seed):
    df = pd.read_csv("../results/depleting_nutrient_B.csv")
    d = df[df.seed==seed]
    dict = {colname: d[colname].iloc[0] for colname in d.columns}
    dict["name"]= "../results/depleting_B/data/"+dict["files"].split("_parameters.csv")[0]
    orbium = Creature(filename=None, dict=dict, vitals=True)
    return orbium
