#!/usr/bin/env python3

"""Script outlining orbium stress test: To be copied into
main scripts for different nutrient types"""

from stochastic_nutrient_A import *

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
        nutrition.append(creature.running_mean)

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

def save_view(grids, name, dim):
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
    plt.savefig(name)
    plt.close()


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

    Ts = [500, 1000, 10000]

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
            temp.to_csv(name+"s"+str(int(o.dict["seed"]))+"N"+str(N)+"L"+str(L)+"T"+str(T)+"stress_test.csv")

    save_view(grids, name=name+"s"+str(int(o.dict["seed"]))+"N"+str(N)+"L"+str(L)+"T"+str(T)+"snapshot.png", dim=[10,3])
