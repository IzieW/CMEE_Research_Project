# !/usr/bin/env python3

"""Figure out right nutrient absoprtion parameter for orbium"""

from stochastic_nutrient_A import *

def get_available_nutrition(seed=0):
    np.random.seed(0)
    orbium = Creature("orbium")
    enviros = []
    overlaps = []
    orbium.absorption = 0.01
    for L in [10, 15, 20, 30]:
        enviro = StochasticEnviro(L=L)
        for i in range(10):
            orbium.initiate()
            enviro.initiate()
            o= run_one(orbium, enviro, overlaps=True)
            print(len(o))
            enviros += list(np.repeat(L, len(o)))
            overlaps += o
    return pd.DataFrame({"enviro": enviros,
                         "overlap": overlaps})

#df = get_available_nutrition()

"""with open("../results/available_nutrient_depletion.csv", "w") as f:
    csvwrite = csv.writer(f)
    for row in df:
        csvwrite.writerow([row])"""


def load_data(depletion=False):
    if depletion:
        name = "../results/available_nutrient_depletion.csv"
    else:
        name = "../results/available_nutrient.csv"

    overlaps=[]
    with open(name, "r") as f:
        csvread = csv.reader(f)
        for row in csvread:
            overlaps.append((float(row[0])))
    return overlaps

def percents(overlaps):
    x = np.arange(int(np.min(overlaps)), int(np.max(overlaps)))
    p = []
    for i in x:
        p.append(np.sum(np.asarray(overlaps) > i)/len(overlaps))

    return pd.DataFrame({"value": x,
                         "percent": p})


def get_times(absorption):
    orbium = Creature("orbium")
    orbium.absorption = absorption
    Ls = [10, 15, 20, 30]
    Ts = [500, 1000, 5000, 10000]
    df = pd.DataFrame()
    for T in Ts:
        temp = []
        for L in Ls:
            enviro = StochasticEnviro(L=L, T=T)
            runs = []
            for i in range(10):
                orbium.initiate()
                enviro.initiate()
                runs.append(run_one(orbium, enviro))
            temp.append(np.mean(runs))
        df["T"+str(T)] = temp

    return df
