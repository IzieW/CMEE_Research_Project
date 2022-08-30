# !/usr/bin/env python3
"""Develop stochastic nutrient environment using 2D stochastic distribution
described in Godany et. al. Nutrient modelled as follows: Overlap with the nutrient
adds to the growth mean proportional to overlap, scaled by the orbium's nutrient parameter.
Meanwhile, every timestep the environment downticks the orbium's growth mean.

To survive, orbium must selectively seek food to remain within mean."""
import matplotlib.pyplot as plt
import matplotlib

from lenia_package2 import *
import os
import re
path = "../results/stochastic_nutrient3/14_hour/"
############################## FUNCTIONS #########################
#### ANALYTICS  #####
def get_files(path):
    """Pulls list of files from input directory in results directory"""
    return os.listdir(path)
def get_seed(file):
    """Return seed from input filename"""
    return int(re.search(r"s\d+", file).group().split("s")[1])
def get_T(file):
    """Return seed from input filename"""
    return int(re.search(r"T\d+", file).group().split("T")[1])
def get_L(file):
    """Return seed from input filename"""
    return int(re.search(r"L\d+", file).group().split("L")[1])
def get_seed(file):
    """Return seed from input filename"""
    return int(re.search(r"s\d+", file).group().split("s")[1])
def get_T(file):
    """Return seed from input filename"""
    return int(re.search(r"T\d+", file).group().split("T")[1])
def get_L(file):
    """Return seed from input filename"""
    return int(re.search(r"L\d+", file).group().split("L")[1])
def load_vitals(path):
    """Load in parameter information from all files. Return in one dataframe"""
    files = os.listdir(path)
    files = [i for i in files if re.search(r"parameters.csv$", i)]
    # files = [i.lower() for i in files]
    # Kick off data frame with first value in list
    df = pd.read_csv(path + files[0], header=None)
    df = df.transpose()
    df = df.set_axis(list(df.iloc[0]), axis=1)  # Tranpose and rename header with column names
    df = df[1:]

    for i in range(1, len(files)):
        temp = pd.read_csv(path + files[i], header=None)
        temp = temp.transpose()
        temp = temp.set_axis(list(temp.iloc[0]), axis=1)  # rename headers
        temp = temp[1:]

        df = pd.concat([df, temp])

    del df["b"]
    df["b"] = 1

    df = df.astype("float64")
    df["files"] = files

    return df
def sort(files, T=None, L=None):
    new = []
    if L:
        for l in L:
            ls = [i for i in files if get_L(i) == l]
            df = pd.DataFrame(ls, [get_T(i) for i in ls])
            for i in list(df.sort_index()[0]): new.append(i)
    else:
        df = pd.DataFrame(files, [get_T(i) for i in files])
        for i in list(df.sort_index()[0]): new.append(i)

    return new
def survival_times(creatures, seed, dim=[4, 3]):
    np.random.seed(seed)
    for i in creatures: i.initiate(), i.enviro.initiate()
    survival_time = [run_one(i, i.enviro) for i in creatures]

    return organise(survival_time, dim=dim)
def get_creature(creatures, seed):
    """find creature from list of creatures"""
    i = 0
    while creatures[i].dict["seed"] != seed:
        i += 1

    return creatures[i]
def load_all(directory, sorting=True, notebook=False, T=[1000, 10000], L=[20, 32, 64]):
    """Load all orbium into list"""
    files = os.listdir(directory)
    files = [i for i in files if re.search(r"parameters.csv$", i)]

    if sorting:
        files = sort(files, T, L)

    organisms = []
    for i in files:
        orbium = Creature(directory + i.split("_parameters.csv")[0], cluster=True)
        if notebook:
            orbium.enviro = StochasticEnviro(l=100, L=get_L(i), T=get_T(i))
        organisms.append(orbium)

    return organisms

def load_creature(directory, seed):
    """Search properties in organisms"""
    files = os.listdir(directory)
    files = [i for i in files if re.search(r"parameters.csv$", i)]
    i = 0
    while get_seed(files[i]) != seed:
        i += 1

    orbium = Creature(directory + files[i].split("_parameters.csv")[0], cluster=True)

    return orbium
bell = lambda x, m, s: np.exp(-((x - m) / s) ** 2 / 2)
def show_seeds(creatures, par="seed", dim=[4, 3]):
    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    parameter = [i.dict[par] for i in creatures]
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = parameter[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0
    return grid
def plot_means(df, parameter, sample_size, *Ts, log=False, title=None):
    """Sorting by correlation lengths and time, plot the means of a input parameter.
    df = data frame of vitals
    parameter= parameter of choice to compare on y
    sample_size = number of seeds per group (used to calculate error)
    *Ts = nonkeyword argument with as many values of T (correlation time) as wanting to sort by"""
    i = 0
    cmap = matplotlib.cm.get_cmap("gist_rainbow")
    for T in Ts:
        means = df[df.correlation_time == T].groupby("correlation_length").describe()[parameter]
        if log:
            plt.plot(np.log(means["mean"]), label = "T="+str(T), c=cmap(i))
            plt.errorbar(means.index, np.log(means["mean"]), yerr = np.log(means["std"]/sample_size), elinewidth=1, capsize=4, color=cmap(i))
        else:
            plt.plot(means["mean"], label = "T="+str(T), c=cmap(i))
            plt.errorbar(means.index, means["mean"], yerr = means["std"]/sample_size, elinewidth=1, capsize=4, color=cmap(i))
        i += 10

    """mean1000 = df[df.correlation_time == 1000].groupby("correlation_length").describe()[parameter]
    mean10000 = df[df.correlation_time == 10000].groupby("correlation_length").describe()[parameter]
    mean100 = df[df.correlation_time == 100].groupby("correlation_length").describe()[parameter]
    plt.plot(mean100["mean"], label = "T=100", c="blue")
    plt.errorbar(mean100.index, mean100["mean"], yerr=mean100["std"]/2, elinewidth=1, capsize=4, color="blue")
    plt.plot(mean1000["mean"], label = "T=1000", c= "green")
    plt.errorbar(mean1000.index, mean1000["mean"], yerr=mean1000["std"]/2, elinewidth=1, capsize=4, color="green")
    plt.plot(mean10000["mean"], label = "T=10000", c="red")
    plt.errorbar(mean10000.index, mean10000["mean"], yerr=mean10000["std"]/2, elinewidth=1, capsize=4, color = "red")"""
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.xlabel("Correlation Length (L)")
    if title:
        plt.suptitle(title)
    if log:
        plt.ylabel("log mean" + parameter)
    else:
        plt.ylabel("mean "+parameter)

    plt.show()
def plot_T(df, parameter, sample_size):
    x = df.groupby("correlation_time").describe()[parameter]
    plt.plot(x["mean"])
    plt.errorbar(x.index, x["mean"], yerr=x["std"]/np.sqrt(sample_size), elinewidth=1, capsize=4)
    plt.xlabel= "Correlation Time (T)"
    plt.ylabel= "Mean "+ parameter

    plt.show()
def organise(data, dim=[4, 3]):
    x_from, y_from = 0, 0
    ind = 0
    grid = np.zeros([dim[1], dim[0]])
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to, x_from:x_to] = data[ind]
            x_from = x_to
            ind += 1
        y_from = y_to
        x_from = 0
    return grid


#### RENDERING ###
def update_neutral(i):
    """Update in neutral environment"""
    global img
    for self in creatures:
        U = np.real(np.fft.ifft2(self.K * np.fft.fft2(self.A)))  # Convolve by kernel to get neighbourhood sums
        self.A = np.clip(self.A + 1 / self.T * (self.growth(U)), 0, 1)

    all_grids = [i.A for i in creatures]

    grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to * 64, x_from:x_to * 64] = all_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    img.set_array(grid)
    return img
def set_neutral(i):
    """Update in neutral environment"""
    all_grids = [i.A for i in creatures]

    grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to * 64, x_from:x_to * 64] = all_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    return grid
def neutral_grid(creatures, dim=[4, 3], initiate=True):
    globals().update({"creatures": creatures, "dim": dim})

    if initiate:
        for i in creatures: i.initiate()

    fig = figure_world(set_neutral(0))

    anim = animation.FuncAnimation(fig, update_neutral, frames=400, interval=20)
    return anim

def figure_world(A, **labels):
    global img
    fig = plt.figure()
    if labels.get("x"):
        plt.xticks(labels.get("x")[0], labels.get("x")[1])
    if labels.get("y"):
        plt.yticks(labels.get("y")[0], labels.get("y")[1])
    if labels.get("xlab"):
        plt.xlabel(labels.get("xlab"))
    if labels.get("ylab"):
        plt.ylabel(labels.get("ylab"))

    img = plt.imshow(A, cmap="viridis", interpolation="nearest", vmin=0)
    #plt.close()
    return fig

def set_grid(i):
    # update board, and creatures
    for n in range(len(creatures)):
        if i % creatures[n].enviro.dt == 0:
            creatures[n].enviro.update()
        update_man(creatures[n], creatures[n].enviro)

    enviro_grids = [c.enviro.grid for c in creatures]
    creature_grids = [c.A for c in creatures]

    # Populate grid pieces with updates
    enviro_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    creature_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_grid[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            creature_grid[y_from:y_to * 64, x_from:x_to * 64] = creature_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    return [np.zeros([dim[1]*64, dim[0]*64]),  creature_grid*2, enviro_grid*0.5]  # need 3 to render RGB

def update_grid(i):
    global img
    # update board, and creatures
    for n in range(len(creatures)):
        if i % creatures[n].enviro.dt == 0:
            creatures[n].enviro.update()
        update_man(creatures[n], creatures[n].enviro)

    enviro_grids = [c.enviro.grid for c in creatures]
    creature_grids = [c.A for c in creatures]

    # Populate grid pieces with updates
    enviro_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    creature_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_grid[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            creature_grid[y_from:y_to * 64, x_from:x_to * 64] = creature_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    img.set_array(np.dstack([np.zeros([dim[1]*64, dim[0]*64]), creature_grid*2, enviro_grid*0.5]))  # need three to render RBG
    return img

def update_grid_std(t):
    global img
    # update board, and creatures
    t *= 25
    for i in range(t, t + 25):
        for n in range(len(creatures)):
            if i % creatures[n].enviro.T == 0:
                creatures[n].enviro.update_L()
            if i % creatures[n].enviro.dt == 0:
                creatures[n].enviro.update()
            update_man(creatures[n], creatures[n].enviro)

    all_grids = [sum([c.A, c.enviro.grid]) for c in creatures]

    # Populate grid pieces with updates
    grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to * 64, x_from:x_to * 64] = all_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    img.set_array(grid)
    return img

def render_grid(creatures, name=None, dim=[4, 3], seed=0, standard=False, initiate=True, frames=400, **labels):
    globals().update({"creatures": creatures, "dim": dim})

    np.random.seed(seed)

    if initiate:
        for i in creatures:
            i.enviro.initiate()
            i.initiate()

    grid = set_grid(0)
    fig = figure_world(np.dstack(grid), **labels)
    print("rendering animation...")

    if standard:  # Standardise to show 10,000 timesteps (jumping)
        anim = animation.FuncAnimation(fig, update_grid_std, frames=400, interval=20)
    else:
        anim = animation.FuncAnimation(fig, update_grid, frames=frames, interval=20)
    if name:
        anim.save(name + "_seed" + str(seed) + ".gif", writer="imagemagick")
    else:
        return anim

def update_seed_grids(i):
    global img
    grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to * 64, x_from:x_to * 64] = all_grids[ind][i]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    img.set_array(grid)
    return img

def set_seed_grids():
    global dim
    grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            grid[y_from:y_to * 64, x_from:x_to * 64] = all_grids[ind][0]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    return grid

def get_seeded_grids(creature, seed):
    np.random.seed(seed)
    creature.enviro.initiate()
    creature.initiate()
    grids = []
    for t in range(400):
        if t % creature.enviro.dt == 0:
            creature.enviro.update()
        update_man(creature, creature.enviro)
        grids.append(sum([creature.A, creature.enviro.grid]))
    return grids

def render_seeds(creature, seed_range=12, dim = [4,3], name=None):
    dict = {"creature": creature,
            "dim": dim,
            "all_grids": [get_seeded_grids(creature, i) for i in range(seed_range)]}
    globals().update(dict)
    fig = figure_world(set_seed_grids())
    print("Rendering animation...")
    anim = animation.FuncAnimation(fig, update_seed_grids, frames=400, interval=20)
    if name:
        anim.save(name, writer="imagemagick")
    else:
        return anim

### ANIMATION ###
def update(i):
    global img
    if i % enviro.dt == 0:
        enviro.update()

    creature.running_mean = creature.running_mean - 0.001  # minus 0.01 from growth mean for every time step
    creature.running_mean = creature.running_mean + (overlap(creature, enviro)*creature.nutrient)  # add nutrients from food
    U = np.real(np.fft.ifft2(creature.K * np.fft.fft2(creature.A)))  # Get neighborhood sums
    # Update board by 1/T * growth sum and obstacle distribution (enviro)
    creature.A = np.clip(creature.A + 1 / creature.T * creature.growth(U), 0, 1)
    img.set_array(sum([creature.A, enviro.grid]))

    return img,

def animate(creature, enviro, name=None):
    globals().update({"creature": creature,
                      "enviro": enviro})

    means = []
    enviro.initiate()
    enviro.update()  # lift off
    creature.initiate()

    fig = figure_world(sum([creature.A, enviro.grid]))
    print("rendering animation...")
    anim = animation.FuncAnimation(fig, update, frames=400, interval=20)
    if name:
        anim.save(name, writer="imagemagick")
    else:
        return anim

###   EVOLUTION   #####
def overlap(creature, feature):
    """Return sum of cells overlapping overlap"""
    overlapping = sum([creature.A, feature.grid]) > feature.grid
    return np.sum(overlapping * feature.grid)

global time_log, theta_log, selection_coef
time_log = pd.DataFrame(columns=["wild_mean", "wild_var"])
theta_log = pd.DataFrame(columns=Creature.keys)
selection_coef = []

def record_time(wild_mean, wild_var):
    """Record timeline"""
    global time_log
    x = pd.DataFrame([[wild_mean, wild_var]],
                     columns=["wild_mean", "wild_var"])  # record averages
    time_log = pd.concat([time_log, x])


def record_theta(creature):
    global theta_log
    theta_log = pd.concat([theta_log, pd.DataFrame(creature.theta())])


def mutate(p):
    """Mutate input parameter p"""
    return np.exp(np.log(p) + np.random.uniform(low=-0.2, high=0.2))


def prob_fixation(wild_time, mutant_time, N):
    """Return probability of fixation given time survived by mutant and wild type,
    and psuedo-population size N"""
    s = (mutant_time - wild_time) / wild_time  # selection coefficient
    selection_coef.append(s)
    """If s is zero, there is no selective difference between types and probability of 
    fixation is equal to 1/populdation size"""
    if s:
        return (1 - np.exp(-2 * s)) / (1 - np.exp(-2 * N * s))
    else:
        return 1 / N

def selection(t_wild, t_mutant):
    """Return winning solution based on survival times of wild type (t_wild) and
    mutant type (t_mutant)."""
    pfix = prob_fixation(wild_time=t_wild, mutant_time=t_mutant, N=population_size)  # get probability of fixation
    if pfix >= np.random.uniform(0, 1):
        # ACCEPT MUTATION
        return True
    else:
        # REJECT MUTATION
        return False


def update_man(creature, food, moving=False, give_sums=False):
    """Update learning channel by 1/T according to values in learning channel A,
    and obstacle channel O"""
    creature.running_mean = creature.running_mean - 0.001  # minus 0.01 from growth mean for every time step
    creature.running_mean = creature.running_mean + (
                overlap(creature, food) * creature.absorption)  # add nutrients from food

    U = np.real(np.fft.ifft2(creature.K * np.fft.fft2(creature.A)))
    creature.A = np.clip(creature.A + 1 / creature.T * creature.growth(U), 0, 1)


def run_one(creature, enviro, show_after=0, moving=False, verbose=False, give_sums=False, return_means=False):
    """Run creature of given parameters in given obstacle configuration until it dies.
    Show after specifies number of timesteps at when it will show what the grid looks like"""
    t = 0  # set timer
    global sums
    sums = np.zeros(10000)
    means = []
    while np.sum(creature.A) and (
            t < 10000):  # While there are still cells in the learning channel, and timer is below cut off

        if t % enviro.dt == 0:
            enviro.update()  # Update enviro

        t += 1  # update timer by 1

        if verbose & (t % 1000 == 0):
            print(t)  # Show that it is working even after long waits

        if give_sums:
            sums[t - 1] = update_man(creature, enviro, moving=moving, give_sums=True)
        else:
            update_man(creature, enviro)  # Run update and show
            means.append(creature.running_mean)
        # if t == show_after:
        #   plt.matshow(sum([creature.A, obstacle.grid]))
    if return_means:
        return means
    else:
        return t

def mutate_and_select(creature, enviro, moving=False, runs=10, minus_keys=4, median=False):
    """Mutate one parameter from creature and assess fitness of new solution agaisnt wild type
    in input obstacle environment. Save winning parameters to Creature.

    Method involve running wild type and mutant over ten distinct obstacle environment, and
    summing the survival time of each."""
    wild_type = creature
    mutant = deepcopy(creature)

    ## Choose parameter at random and mutate in mutant_type
    x = np.random.randint(0, len(Creature.keys) - minus_keys)
    mutant.__dict__[Creature.keys[x]] = mutate(mutant.__dict__[Creature.keys[x]])
    mutant.K = mutant.kernel()  # update mutant kernel

    # Run mutant and wild over runs number of obstacle configurations
    t_wild = np.zeros(runs)
    t_mutant = np.zeros(runs)
    for i in range(runs):
        enviro.initiate()  # configure environments
        wild_type.initiate()
        mutant.initiate()
        t_wild[i] = run_one(wild_type, enviro, moving=moving)
        t_mutant[i] = run_one(mutant, enviro, moving=moving)

    # Record mean and variance of survival times
    if median:
        wild_mean = np.median(t_wild)
        mutant_mean = np.median(t_mutant)
    else:
        wild_mean = t_wild.mean()
        mutant_mean = t_mutant.mean()
    print(wild_mean)
    record_time(wild_mean=wild_mean,
                wild_var=t_wild.var())
    record_theta(wild_type)
    # Select winning parameter
    if selection(wild_mean, mutant_mean):
        print("Accept mutation")
        creature.update_theta(mutant)  # Update creature parameters
        return True
    else:
        print("Reject mutation")
        return False

def optimise(creature, obstacle, N=1000, k=1, seed=0, fixation_mark=600, max_time = 100, name=None, minus_key=1, median=False):
    """Mutate and select input creature in psuedo-population of size N
    until wild type becomes fixed over fixation number of generations.

    Creature = Starting creature to mutate
    Obstacle = Environment to Creature is evolved in
    N = psuedo-population size
    fixation_mark = number of generations without mutation that specifies arrival at a global optima
    max_time = maximum time to run for before manually ending it"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    max_time = max_time * 60  # Translate to seconds
    """Evolve until parameters become fixed over fixation number of generations"""
    gen = 0  # time_count
    mutation = 0
    fix = 0
    start = time.time()
    while (fix < fixation_mark) and ((time.time() - start) < max_time):
        print("mutations: ", mutation)
        if mutate_and_select(creature, obstacle, minus_keys=minus_key, median=median):  # Updates creature values
            mutation += 1
            fix = 0  # reset everytime there is a mutation
        else:
            fix += 1
        gen += 1
        print(fix)

    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    if name:
        creature.name = name
    else:
        creature.name = str(creature.n) + "_orbium_fix" + str(fix) + "s" + str(seed) + "N" + str(
            N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, obstacle, summary=True, median=median)
    creature.mutations = mutation
    creature.generations = gen
    creature.survival_mean = survival_time[0]
    creature.survival_var = survival_time[1]
    creature.fix = fix

    with open(name+"_parameters.csv", "w") as f:
        csvwrite = csv.writer(f)
        for key in Creature.keys:
            csvwrite.writerow([key, creature.__dict__[key]])
        csvwrite.writerow(["mutations", mutation])
        csvwrite.writerow(["survival_mean", creature.survival_mean])
        csvwrite.writerow(["survival_var", creature.survival_var])
        csvwrite.writerow(["absorption", creature.absorption])
        csvwrite.writerow(["fixation", fix])
        csvwrite.writerow(["generations", gen])
        csvwrite.writerow(["organism_count", creature.n])
        csvwrite.writerow(["l", 64])
        csvwrite.writerow(["correlation_length", obstacle.L])
        csvwrite.writerow(["correlation_time", obstacle.T])
        csvwrite.writerow(["seed", seed])
        csvwrite.writerow(["N", N])
        csvwrite.writerow(["k", k])

    time_log.to_csv(creature.name + "_times.csv")  # Save timelog to csv
    theta_log.to_csv(creature.name + "_theta.csv")

    with open(creature.name + "_selection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in selection_coef:
            csvwrite.writerow([i])

    return 1


def optimise_timely(creature, obstacle, N, seed=0, run_time=10, moving=False, cluster=False, name=None, minus_key=1, median=False):
    """Mutate and select input creature in psuedo-population of size N
    until wild type becomes fixed over fixation number of generations"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    run_time = run_time * 60  # Translate to seconds
    """Evolve until parameters become fixed over fixation number of generations"""
    gen = 0  # time_count
    mutation = 0
    start = time.time()
    while (time.time() - start) < run_time:
        print("mutations: ", mutation)
        if mutate_and_select(creature, obstacle, moving=moving, minus_keys=minus_key, median=median):  # Updates creature values
            mutation += 1
        gen += 1

    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    if name:
        creature.name = name
    else:
        creature.name = str(creature.n) + "_orbium_t" + str(run_time) + "s" + str(seed) + "N" + str(
            N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, obstacle, summary=True, median=median)
    creature.mutations = mutation
    creature.survival_mean = survival_time[0]
    creature.survival_var = survival_time[1]
    if cluster:
        time_log.to_csv(creature.name + "_times.csv")  # Save timelog to csv
        theta_log.to_csv(creature.name + "_theta.csv")
    else:
        time_log.to_csv("../results/" + creature.name + "_times.csv")  # Save timelog to csv

    #creature.save(cluster=cluster)  # save parameters
    return 1


def get_survival_time(creature, obstacle=None, runs=10, summary=False, verbose=False, median=False):
    """Calculate average run time over seeded 10 configurations.
    Return mean and variance."""
    times = np.zeros(runs)
    if obstacle:
        for i in range(1, runs + 1):
            np.random.seed(i)
            creature.initiate()  # Reset grid
            obstacle.initiate()  # set obstacle
            times[i - 1] = run_one(creature, obstacle, verbose=verbose)
    else:
        for i in range(1, runs + 1):
            creature.initiate()
            times[i - 1] = run_one(creature, obstacle, verbose=verbose)

    if summary:
        if median:
            return np.median(times), times.var()
        else:
            return times.mean(), times.var()
    else:
        return times


def testing_limits():
    orbium = Creature("orbium")
    enviro = StochasticEnviro()
    for i in range(3):
        orbium.initiate()
        plt.plot(run_one(orbium, enviro)[1], label = str(orbium.nutrient))
        orbium.nutrient = orbium.nutrient * 10
    plt.legend()
    plt.hlines(0.15, -10, 300, colors="black", linestyles="dotted")

def test_equilibrium_m(creature):
    """Test survival time as growth mean varies"""

    m = np.arange(-0.1, 0.17, 0.01)
    survival_time = list(np.zeros(len(m)))

    for i in range(len(m)):
        creature.m = m[i]
        creature.initiate()
        survival_time[i] = lenia_package.run_one(creature, enviro=None, verbose=False)

    return pd.DataFrame({
        "m": m,
        "survival_time": survival_time
    })


def render_results(path, filename):
    orbium = load_all(path, T=[100, 1000, 10000], L=[10, 20, 32])
    print("Rendering all results...")
    # Render all
    #render_grid(orbium, name=filename+"all", dim=[12,3])

    # Render by group
    orbium100 = [i for i in orbium if i.enviro_T==100]
    render_grid(orbium100, name=filename+"T100", dim=[4,3])
    orbium1000 = [i for i in orbium if i.enviro_T==1000]
    render_grid(orbium1000, name=filename+"T1000", dim=[4,3])
    orbium10000 = [i for i in orbium if i.enviro_T==10000]
    render_grid(orbium10000, name=filename+"T10000", dim=[4,3])

