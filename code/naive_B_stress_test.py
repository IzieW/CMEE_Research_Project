#!/usr/bin/env python3

"""Script outlining orbium stress test: To be copied into
main scripts for different nutrient types"""
import csv

import matplotlib.pyplot as plt
import pandas as pd

from stochastic_nutrient_B import *

#################### RUN STRESS TEST ###########################
#iter = int(os.environ.get("PBS_ARRAY_INDEX"))

def run_one_stress_test(creature, enviro, snapshot=None):
    """Run creature of given parameters in given environmental configuration until it dies.
    Return number of timesteps creature survived for before death"""
    t = 0  # set timer
    mass =[]
    volume = []
    mean_centroid = [] # mean center
    var_centroid = [] # variance from centroid
    nutrition = []

    grid = np.zeros([64,64])

    while np.sum(creature.A) and (
            t < 10000):  # While there are still cells in the learning channel, and timer is below cut off

        if t % enviro.dt == 0:
            enviro.update()  # Update enviro

        t += 1  # update timer by 1

        creature.mean_centroid()
        mean_centroid.append(creature.centroid)
        var_centroid.append(creature.centroid_var())
        mass.append(np.sum(creature.A))
        volume.append(np.sum(creature.A > 0))
        nutrition.append(creature.nutrition)

        update_man(creature, enviro)  # Run update and show


        if t == snapshot:
            grid = [creature.A, enviro.grid]



    return (pd.DataFrame({
                        "nutrient": nutrition,
                         "centroid_x": [c[0] for c in mean_centroid],
                        "centroid_y": [c[1] for c in mean_centroid],
                        "var_x": [c[0] for c in var_centroid],
                          "var_y": [c[1] for c in var_centroid],
                         "mass": mass,
                         "volume": volume}), grid)

def save_view(grids, name, dim, seeds=None, Ts=None, title=None):
    enviro_fill = np.zeros([64*dim[1], 64*dim[0]])
    creature_fill = np.zeros([64*dim[1], 64*dim[0]])

    enviro_grids = [i[1] for i in grids]
    creature_grids = [i[0] for i in grids]

    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_fill[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            creature_fill[y_from:y_to * 64, x_from:x_to * 64] = creature_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0


    plt.matshow(np.dstack([np.zeros([64*dim[1], 64*dim[0]]), creature_fill*2, enviro_fill*0.5]))

    if Ts:
        plt.xticks = [np.arange(32, (64*dim[0])+32, 64), seeds]
        plt.yticks = [np.arange(32, (64*dim[1])+32, 64), Ts]
        plt.ylabel("Correlation time (T)")
        plt.xlabel("Top row seeds")
    if title:
        plt.suptitle(title)

    if name:
        plt.savefig(name)
        plt.close()

def show_cluster(cluster, dim, name=None, path="../results/naive_B/morpho_cluster1.csv"):
    """save png of cluster results"""
    df = pd.read_csv(path)
    d = df[df.cluster==cluster]

    names = [n.split("_parameters.csv")[0] for n in d.files]
    d["name"] = names

    dict = [{colname: d[colname].iloc[n] for colname in d.columns} for n in range(len(d))]

    orbium = [Creature(filename=None, dict=i, cluster=True) for i in dict]

    grids = []

    for o in orbium:
        o.initiate()
        o.enviro.initiate()
        for t in range(int(o.survival_mean/2)):
            if t % o.enviro.dt == 0:
                o.enviro.update()
            update_man(o, o.enviro)
        grids.append([o.A, o.enviro.grid])

    save_view(grids, name=name, dim=dim, title="cluster "+str(cluster))

def stress_test(creature):
    """Run stress test three times and return best run"""
    best = []
    for i in range(3):
        creature.initiate()
        creature.enviro.initiate()
        temp = run_one_stress_test(creature, creature.enviro, snapshot=int(creature.dict["survival_mean"]/2))
        if len(temp)>len(best):
            best = deepcopy(temp) # get best run out of them
    return best

def stress_test_set(i, path, name):
    """path to files, name to use for saved files"""
    df = pd.read_csv(path)
    TL = [[n, l] for n in [10, 100, 1000, 10000, 100000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]

    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]

    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names

    grids = []


    # Load in creatures for stress test
    for T in Ts:  # for each row of T
        d = df[df.correlation_time == T]
        d = d.sort_values("seed")
        temp_dict = [{colname: d[colname].iloc[n] for colname in d.columns} for n in range(len(d))]
        orbium = [Creature(filename=None, dict=i, cluster=True) for i in temp_dict]  # load in creatures
        for o in orbium:
            r = stress_test(o)
            grids.append(r[1]) # add snapshot to grids
            temp = pd.DataFrame(r[0]) # save results
            temp.to_csv("naive_nutrient_B_s"+str(int(o.dict["seed"]))+"l64L"+str(L)+"T"+str(T)+"N"+str(N)+"_stress_test.csv")
        if len(orbium) !=10:
            for n in range(10-len(orbium)):
                grids.append([np.zeros([64,64]), np.zeros([64,64])])

    d = df[df.correlation_time==Ts[0]]
    d = d.sort_values("seed")
    seeds = list(d.seed)
    if len(seeds) != 10:
        for i in range(10-len(seeds)): seeds.append(0)

    save_view(grids, name=name+"s"+str(int(o.dict["seed"]))+"N"+str(N)+"L"+str(L)+"T"+str(T)+"snapshot.png", dim=[10,3], seeds=seeds, Ts=Ts)

def get_morphology(i, path):
    """run until at morphology half way and save grids to csv"""
    df = pd.read_csv(path)
    TL = [[n, l] for n in [10, 100, 1000, 10000, 100000] for l in [10, 15, 20, 30]]  # all combinations of T and L

    N = TL[i][0]
    L = TL[i][1]

    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]

    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names

    grids = []

    # Load in creatures for stress test
    for T in Ts:  # for each row of T
        d = df[df.correlation_time == T]
        d = d.sort_values("seed")
        temp_dict = [{colname: d[colname].iloc[n] for colname in d.columns} for n in range(len(d))]
        orbium = [Creature(filename=None, dict=i, cluster=True) for i in temp_dict]  # load in creatures
        for o in orbium:
            temp = stress_test(o)[1]
            name = "s"+str(int(o.dict["seed"]))+"l64L"+str(L)+"T"+str(T)+"N"+str(N)
            np.savetxt("naive_B_"+name+"_cells.csv", temp[0])  # save to text
            np.savetxt("naive_B_"+name+"_enviro.csv", temp[1])


#get_morphology(iter, "naive_nutrient_B.csv")

######### VIEWING INFO IN GRIDS  #########
def organise_grid(grids, dim):
    organ_grid = np.zeros([dim[1], dim[0]])

    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            organ_grid[y_from:y_to, x_from:x_to] = grids[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0

    return organ_grid

def get_set(N, L, mkfile= False, parameter="seed", path="../results/naive_nutrient_B.csv"):
    """path to files, name to use for saved files"""
    df = pd.read_csv(path)
    Ts = [100, 1000, 10000]

    df = df[df.N==N]
    df = df[df.correlation_length==L]

    names = [n.split("_parameters.csv")[0] for n in df.files]
    df["name"] = names

    grids = []


    # Load in creatures for stress test
    for T in Ts:  # for each row of T
        d = df[df.correlation_time == T]
        d = d.sort_values("seed")
        temp_dict = [{colname: d[colname].iloc[n] for colname in d.columns} for n in range(len(d))]
        for o in temp_dict:
            grids.append(o[parameter]) # add snapshot to grids
        if len(temp_dict) !=10:
            for n in range(10-len(temp_dict)):
                grids.append(0)

    if mkfile:
        fgrids =["seed"]
        fgrids += grids
        with open("../results/naive_B/N"+str(N)+"L"+str(L)+"qualitative.csv", "w") as f:
            csvwrite = csv.writer(f)
            csvwrite.writerow(fgrids)
            csvwrite.writerow(["forms"])

    return organise_grid(grids, [10,3])

###################### LOAD IN RESULTS ############################
def load_one(seed, path="../results/naive_B/data/", file=False):
    files = os.listdir(path)
    files = [i for i in files if re.search(r"stress_test.csv$", i)]
    files = [i for i in files if get_seed(i)==seed]

    df = pd.read_csv(path+files[0])
    del df["Unnamed: 0"]

    if file:
        return df, files[0]
    else:
        return df

def trajectory(seed):
    temp = load_one(seed, file=True)

    file = temp[1]
    d = temp[0]

    plt.plot(d.centroid_x, d.centroid_y)
    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.suptitle("Centroid phase space trajectory")
    plt.title("Seed "+str(seed)+", N="+get_N(file)+", L="+get_L(file)+", T="+get_T(file))

def centroid_phase(seed):
    temp = load_one(seed, file=True)
    file = temp[1]
    d = temp[0]
    plt.errorbar(d.centroid_x, d.centroid_y, yerr=np.sqrt(d.var_y/len(d.var_y)),
                   xerr= np.sqrt(d.var_x/len(d.var_x)), capsize=2, ecolor="blue")
    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.suptitle("Centroid phase space trajectory with variance")
    plt.title("Seed "+str(seed)+", N="+str(get_N(file))+", L="+str(get_L(file))+", T="+str(get_T(file)))


def load_tests(path = "../results/naive_B/data/"):
    """Load in summary of stress tests"""
    files = os.listdir(path)
    df = load_vitals(path)
    files = [i for i in files if re.search(r"stress_test.csv$", i)]

    ## List all measures ##
    seed = np.zeros(len(files))
    mean_mass = np.zeros(len(files))
    var_mass = np.zeros(len(files))
    mean_volume = np.zeros(len(files))
    var_volume = np.zeros(len(files))
    mean_density = np.zeros(len(files))
    var_density = np.zeros(len(files))

    mean_centroid_x = np.zeros(len(files))
    var_centroid_x = np.zeros(len(files))
    mean_centroid_y = np.zeros(len(files))
    var_centroid_y = np.zeros(len(files))

    mean_nutrient = np.zeros(len(files))
    var_nutrient = np.zeros(len(files))

    snap_centroid_x = np.zeros(len(files))  # centroid value at snapshot time
    snap_var_x = np.zeros(len(files)) # centroid var at snapshot time
    snap_centroid_y = np.zeros(len(files))  # centroid value at snapshot time
    snap_var_y = np.zeros(len(files)) # centroid var at snapshot time
    snap_mass = np.zeros(len(files))  # mass at snapshot time
    snap_volume = np.zeros(len(files)) # volume at snapshot time


    final_nutrient = np.zeros(len(files)) # final value of nutrient

    # Fill measures #
    for i in range(len(files)):
        seed[i] = get_seed(files[i])
        t = pd.read_csv(path+files[i])
        mean_mass[i] = t.mass.mean()
        var_mass[i] = t.mass.var()
        mean_volume[i] = t.volume.mean()
        var_volume[i] = t.volume.var()
        mean_density[i] = (t.mass/t.volume).mean()
        var_density[i] = (t.mass/t.volume).var()
        mean_centroid_x[i] = t.centroid_x.mean()
        var_centroid_x[i] = t.var_x.mean()
        mean_centroid_y[i] = t.centroid_y.mean()
        var_centroid_y[i] = t.var_y.mean()
        mean_nutrient[i] = t.nutrient.mean()
        var_nutrient[i] = t.nutrient.var()

        d = df[df.seed==seed[i]]
        snap = int(d.survival_mean/2)

        if len(t) > snap:
            snap_centroid_x[i] = t.centroid_x.iloc[snap]
            snap_var_x[i] = t.var_x.iloc[snap]
            snap_centroid_y[i] = t.centroid_y.iloc[snap]
            snap_var_y[i] = t.var_y.iloc[snap]
            snap_mass[i] = t.mass.iloc[snap]
            snap_volume[i] = t.volume.iloc[snap]
        final_nutrient[i] = t.nutrient.iloc[-1]

    big_boy = pd.DataFrame({
        "seed2": seed,
        "mean_mass": mean_mass,
        "var_mass": var_mass,
    "mean_volume": mean_volume,
    "var_volume": var_volume,
    "mean_density": mean_density,
    "var_density": var_density,
    "mean_centroid_x": mean_centroid_x,
    "var_centroid_x": var_centroid_x,
    "mean_centroid_y": mean_centroid_y,
    "var_centroid_y": var_centroid_y,
    "mean_nutrient": mean_nutrient,
    "var_nutrient": var_nutrient,
    "snap_centroid_x": snap_centroid_x,
    "snap_var_x": snap_var_x,
    "snap_centroid_y": snap_centroid_y,
    "snap_var_y": snap_var_y,
    "snap_mass": snap_mass,
    "snap_volume":snap_volume,
    "final_nutrient": final_nutrient
    })

    big_boy = big_boy.sort_values("seed2")
    df = df.sort_values("seed")

    for colname in big_boy.columns:
        df[colname]= list(big_boy[colname])

    df["snap_density"] = df.snap_mass/df.snap_volume
    return df


def load_qualitative(main_file = "../results/naive_nutrient_B.csv", path = "../results/naive_B/"):
    """load qualitative files and save to main df"""
    files = os.listdir(path)
    files = [i for i in files if re.search(r"qualitative.csv$", i)]
    qual = []
    seed = []

    for file in files:
        temp = pd.read_csv(path+file, header=None)
        temp = temp.T
        temp = temp[1:]
        seed += [int(i) for i in list(temp[0])]
        qual += [str(i) for i in list(temp[1])]

    new_d = pd.DataFrame({"seed": seed,
                         "form": qual})

    df = pd.read_csv(main_file)
    new_d = new_d.sort_values("seed")
    df = df.sort_values("seed")
    df["qualitative"] = list(new_d.form)

    return df






