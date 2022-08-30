#!/usr/bin/env python3

from two_nutrient_B import *

#iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def do_render(i):
    df = pd.read_csv("two_nutrient_B.csv")

    TL = [[n, l] for n in [10, 100, 1000, 10000, 100000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]

    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]

    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names
    df["orbium"] = names
    df["kernel_m"] = np.repeat(None, len(df))
    df["fix"] = np.repeat(None, len(df))

    filler = Creature("orbium", species="Not", cluster=True)
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    orbium = []
    for T in Ts:
        d = df[df.correlation_time == T]
        d = d.sort_values("seed")
        temp_dict = [{colname: d[colname].iloc[n] for colname in d.columns} for n in range(len(d))]
        t1 = [n for n in temp_dict if n["creature"]=="c1"]
        t2 = [n for n in temp_dict if n["creature"]=="c2"]
        temp = []
        for t in range(len(t1)):
            temp.append([Creature(filename=None, dict=t1[t], cluster=True), Creature(filename=None, dict=t2[t], cluster=True)])
        if len(temp) != 10:
            for t in range(10-len(temp)):
                temp.append([filler, filler])
        orbium += temp

    dim = [10, 3]

    seeds = [orbium[t][0].dict["seed"] for t in range(10)]

    return seeds
    labels = {"xlab":"top row seeds",
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds],
              "ylab":"correlation_time (T)",
              "y": [np.arange(32, (64*dim[1])+32, 64), Ts],
              "title": "N="+str(N)+", L="+str(L)}

    render_grid(orbium, dim = [10, 3], name= "two_B_N"+str(N)+"L"+str(L), **labels)

#do_render(iter)
