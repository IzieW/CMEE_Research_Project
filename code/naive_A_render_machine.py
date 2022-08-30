# !/usr/bin/env python3

from stochastic_nutrient_A import *

iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def do_render(i):
    # read in data
    df = pd.read_csv("naive_nutrient_A.csv")
    del df["Unnamed: 0"]
    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names

    forms = [i for i in df.qualitative.unique() if i != "nan" if i != "na"]

    form = forms[i] # assign form

    df = df[df.qualitative==form].sort_values("survival_mean").iloc[-3:]  # take top three


    if form == "slow_sensing":
        frames = 30000
        standard= True
    else:
        frames = int(df.survival_mean.max() + 100)
        standard=False

    dicts = [{colname: df[colname].iloc[n] for colname in df.columns} for n in range(len(df))]

    for dict in dicts:
        orbium = Creature(filename=None, dict=dict, cluster=True)
        s = int(dict["seed"])
        N = int(dict["N"])
        T = int(dict["correlation_time"])
        L = int(dict["correlation_length"])

        labels = {"title": form+" (seed="+str(s)+", N="+str(N)+", L="+str(L)+", T="+str(T)+")"}

        for t in range(3):
            animate(orbium, orbium.enviro, name=form+"_s"+str(s)+"N"+str(N)+"L"+str(L)+"T"+str(T)+"seed"+str(t)+".mp4",
                    seed=t, frames=frames, standard=standard, **labels)
            animate(orbium, orbium.enviro, name=form+"_s"+str(s)+"N"+str(N)+"L"+str(L)+"T"+str(T)+"seed"+str(t)+".gif",
                    seed=t, frames=int(frames/10), standard=standard, **labels)


do_render(iter)
