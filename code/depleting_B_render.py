# !/usr/bin/env python3

from depleting_nutrient_B import *

#iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def do_render_B(i):
    # read in data
    df = pd.read_csv("depleting_nutrient_B.csv")

    TL = [[n, l] for n in [10, 100, 1000, 10000, 100000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]

    Ts = [500, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]


    df = df.sort_values("seed")
    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names

    dicts = [{colname: df[colname].iloc[n] for colname in df.columns} for n in range(len(df))]

    dim = [10, 3]

    # Load orbium
    orbium = [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]

    filler = Creature("orbium", species="Not", cluster=True)
    filler.enviro = StochasticEnviro()
    filler.dict["seed"] = 0

    new = []
    for T in Ts:
        temp = [n for n in orbium if n.enviro_T == T]
        if len(temp) != 10:
            for n in range(10-len(temp)): temp.append(filler)  # adjust length to make sure 10 in each row
        new += temp  # add to list

    orbium = deepcopy(new)

    seeds = [orbium[n].dict["seed"] for n in range(10)]  # get first seeds

    labels = {"xlab":"top row seeds",
                "x": [np.arange(32, (64*dim[0])+32, 64), seeds],
              "ylab":"correlation_time (T)",
              "y": [np.arange(32, (64*dim[1])+32, 64), Ts],
              "title": "N="+str(N)+", L="+str(L)
                       }

    render_grid(orbium, dim=dim, name="../results/depleting_B/video/depleting_B_N"+str(N)+"L"+str(L), video=True, frames=10000, **labels)
    #render_grid(orbium, dim=dim, name="depleting_B_N"+str(N)+"L"+str(L), video=False, frames=1000, **labels)


#do_render_B(iter)
