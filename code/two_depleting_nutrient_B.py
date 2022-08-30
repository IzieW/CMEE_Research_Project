# !/usr/bin/env python3
"""House evolution of orbium in stochastic attractant distribution made up of
nutrient source B.

In nutrient source A, the orbium has a "nutrition" parameter, which is multiplied by U
just before the growth distribution is calculated. The nutrition parameter starts at 1, and it constantly
depleted by 0.01 each timestep. The orbium is able to recover nutrition by seeking nutrients in its environment.

Crucially, the orbium cannot eat too much- its only task is constantly seeking. The max nutrient which can be stored is 2,
giving the orbium 200 timesteps before it runs out and begins depleting.

In this variation, n number of orbium are placed in the same grid-space and evolutionary arena.
The two are not killed if they overlap, but they are in competition to absorb nutrients which is depleting.

Each generation they are mutated independently. Taking turns which is mutated first.
"""

####### IMPORTS ########
import pandas as pd

from lenia_package import *
import os
import re

################################# ANALYTICS ################################
def get_T(file):
    """Return seed from input filename"""
    return int(re.search(r"T\d+", file).group().split("T")[1])
def get_N(file):
    """Return seed from input filename"""
    return int(re.search(r"N\d+", file).group().split("N")[1])

def get_a(file):
    return re.search(r"absorption\d......", file).group().split("absorption")[1].split("_")[0]

def get_L(file):
    """Return seed from input filename"""
    return int(re.search(r"L\d+", file).group().split("L")[1])

def get_seed(file):
    """Return seed from input filename"""
    return int(re.search(r"s\d+", file).group().split("s")[1])

def load_s(file):
    s = []
    with open(file, "r") as f:
        csvread = csv.reader(f)
        for row in csvread:
            s.append(float(row[0]))
    return s

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

    df["creature"] = [re.search(r"c\d", i).group() for i in df.files]


    return df

def load_all(directory, sorting=True, notebook=False, L=[20, 32, 64]):
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

def sort(files, Ns, Ls):
    """Sort files according to three lists"""
    final = []
    # get all absorption values
    As = np.unique([re.search(r"absorption\d......", i).group().split("absorption")[1].split("_")[0] for i in files])
    for a in As:
        fa = [i for i in files if get_a(i)==a]
        for N in Ns:
            f = [i for i in fa if get_N(i)==N]
            for L in Ls:
                ln = [i for i in f if get_L(i)==L]
                final += ln
    return final

def sort_pairs(path, splitby="parameters"):
    """Sort multiple orbium into their pairs, return df of files according to splitby"""
    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]

    files = os.listdir(path)
    f1 = [i for i in files if re.search(r"c1", i)]
    f2 = [i for i in files if re.search(r"c2", i)]

    f1 = sort(f1, Ns, Ls)
    f2 = sort(f2, Ns, Ls)

    f1 = [i for i in f1 if re.search(r"_"+splitby+".csv$", i)]
    f2 = [i for i in f2 if re.search(r"_"+splitby+".csv$", i)]

    return pd.DataFrame({"c1": f1,
                         "c2": f2})

def load_pairs(path):
    """load creature list with pairs of creatures for render"""
    df = sort_pairs(path)
    creatures = []

    c1 = [Creature(path+i.split("_parameters.csv")[0], cluster=True) for i in df.c1]
    c2 = [Creature(path+i.split("_parameters.csv")[0], cluster=True) for i in df.c2]

    for i in range(len(df)):
        creatures.append([c1[i], c2[i]])

    return creatures

def plot_s(creatures, dat="selection"):
    plt.plot(creatures[0].__dict__[dat], label="c1, mutations "+str(creatures[0].dict["mutations"]))
    plt.plot(creatures[1].__dict__[dat], label="c2, mutations "+str(creatures[1].dict["mutations"]))

    plt.legend()
    plt.xlabel("Generations")
    plt.ylabel(dat)
    plt.suptitle("Multiple B selection over time")
    plt.title("T=1000, N="+str(creatures[0].dict["N"])+" L="+str(creatures[0].enviro_L)+" absorption "+str(creatures[0].dict["absorption"]))

################################# RENDERING #################################

def update_grid(i):
    # update board, and creatures
    global img

    for n in range(len(creatures)):
        if i % creatures[n][0].enviro.dt == 0:
            creatures[n][0].enviro.update()
        update_man(*creatures[n], creatures[n][0].enviro)

    enviro_grids = [c[0].enviro.grid for c in creatures]
    c1_grids = [c[0].A for c in creatures]
    c2_grids = [c[1].A for c in creatures]

    # Populate grid pieces with updates
    enviro_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c1_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c2_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_grid[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            c1_grid[y_from:y_to * 64, x_from:x_to * 64] = c1_grids[ind]
            c2_grid[y_from:y_to * 64, x_from:x_to * 64] = c2_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0
    img.set_array(np.dstack([c1_grid*2, c2_grid*2, enviro_grid*0.5]))
    return img

def update_grid_std(i):
    # update board, and creatures
    global img

    i *= 25
    for i in range(t, t+25):
        for n in range(len(creatures)):
            if i % creatures[n][0].enviro.dt == 0:
                creatures[n][0].enviro.update()
            update_man(*creatures[n], creatures[n][0].enviro)

    enviro_grids = [c[0].enviro.grid for c in creatures]
    c1_grids = [c[0].A for c in creatures]
    c2_grids = [c[1].A for c in creatures]

    # Populate grid pieces with updates
    enviro_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c1_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c2_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_grid[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            c1_grid[y_from:y_to * 64, x_from:x_to * 64] = c1_grids[ind]
            c2_grid[y_from:y_to * 64, x_from:x_to * 64] = c2_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0
    img.set_array(np.dstack([c1_grid*2, c2_grid*2, enviro_grid*0.5]))
    return img


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
    if labels.get("title"):
        plt.title(labels.get("title"))

    img = plt.imshow(A, cmap="viridis", interpolation="nearest", vmin=0)
    #plt.close()
    return fig


def set_grid(i):
    # update board, and creatures
    for n in range(len(creatures)):
        if i % creatures[n][0].enviro.dt == 0:
            creatures[n][0].enviro.update()
        update_man(*creatures[n], creatures[n][0].enviro)

    enviro_grids = [c[0].enviro.grid for c in creatures]
    c1_grids = [c[0].A for c in creatures]
    c2_grids = [c[1].A for c in creatures]

    # Populate grid pieces with updates
    enviro_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c1_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    c2_grid = np.zeros([64 * dim[1], 64 * dim[0]])
    x_from, y_from = 0, 0
    ind = 0  # index from grids
    for y_to in range(1, dim[1] + 1):
        for x_to in range(1, dim[0] + 1):
            enviro_grid[y_from:y_to * 64, x_from:x_to * 64] = enviro_grids[ind]
            c1_grid[y_from:y_to * 64, x_from:x_to * 64] = c1_grids[ind]
            c2_grid[y_from:y_to * 64, x_from:x_to * 64] = c2_grids[ind]
            x_from = x_to * 64
            ind += 1
        y_from = y_to * 64
        x_from = 0

    return [c1_grid*2,  c2_grid*2, enviro_grid*0.5]  # need 3 to render RGB


def render_grid(creatures, name=None, dim=[4, 3], seed=0, standard=False, initiate=True, frames=400, **labels):
    # Creatures are tuple of creatures c1 and c2
    globals().update({"creatures": creatures, "dim": dim})

    np.random.seed(seed)

    if initiate:
        for i in creatures:
            i[0].initiate(random=True)
            i[1].initiate(random=True)
            i[0].enviro.initiate()

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

def show(c1, c2, enviro):
    plt.matshow(sum([c1.A, c2.A, enviro.grid]))

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

def update(i):
    global img
    if i % enviro.dt == 0:
        enviro.update()
    update_man(c1, c2, enviro)
    img.set_array(np.dstack([c1.A, c2.A, enviro.grid/2]))

def animate(c1, c2, enviro, name=None, seed=0):
    np.random.seed(seed)
    globals().update({"c1": c1,
                      "c2": c2,
                      "enviro": enviro})
    enviro.initiate()
    enviro.update()
    c1.initiate(random=True)
    c2.initiate(random=True)

    fig = figure_world(np.dstack([c1.A, c2.A, enviro.grid/2]))

    print("rendering animation...")
    anim = animation.FuncAnimation(fig, update, frames=400, interval=20)
    if name:
        anim.save(name, writer="imagemagick")
    else:
        return anim

def save(data, name):
    with open(name, "w") as f:
        csvwrite = csv.writer(f)
        for row in data:
            csvwrite.writerow(row)

################################### EVOLUTION ################################
def overlap(creature, feature, sums=True, two_creatures=False):
    """Return sum of cells in environment which overlap with creature"""
    if two_creatures:
        overlapping = sum([creature.A, feature.A]) > feature.A
        if sums:
            return np.sum(overlapping * feature.A)
        else:
            return (overlapping * feature.A) > 0
    else:
        overlapping = sum([creature.A, feature.grid]) > feature.grid
        if sums:
            return np.sum(overlapping * feature.grid)
        else:
            return (overlapping * feature.grid) > 0

def mutate(p):
    """Mutate input parameter p"""
    return np.exp(np.log(p) + np.random.uniform(low=-0.2, high=0.2))

def prob_fixation(c, wild_time, mutant_time, N, k=1):
    """Return probability of fixation given time survived by mutant and wild type,
    and psuedo-population size N and scaling coefficient k"""
    s = ((mutant_time - wild_time) / wild_time)*k  # selection coefficient
    c.selection_log.append(s)  # record selection coefficient
    """If s is zero, there is no selective difference between types and probability of 
    fixation is equal to 1/populdation size"""
    if s:
        return (1 - np.exp(-2 * s)) / (1 - np.exp(-2 * N * s))
    else:  # if s is zero, prob is equal to 1/N
        return 1 / N

def selection(c, t_wild, t_mutant, k):
    """Return winning solution based on survival times of wild type (t_wild) and
    mutant type (t_mutant)."""
    pfix = prob_fixation(c, wild_time=t_wild, mutant_time=t_mutant, N=population_size, k=k)  # get probability of fixation
    if pfix >= np.random.uniform(0, 1):
        # ACCEPT MUTATION
        return True
    else:
        # REJECT MUTATION
        return False

def metabolism(creature, enviro):
    """Update creature's running mean according to depletion and consumption of nutrients"""
    creature.nutrition -= 0.01  # minus 0.01 from nutrition for every timestep
    creature.nutrition += overlap(creature, enviro) * creature.absorption  # add nutrients from food
    creature.nutrition = min(creature.nutrition, 1)  # max they can have is 1- gives them 40 timesteps to death

def update_man(c1, c2, enviro):
    """Update learning channel by 1/T according to values in learning channel A of creature1 and creature2,
    and environmental nutrient gradient. Update running mean according to depletion and
    recovery of nutrients"""

    if np.random.randint(0, 10) < 5:  # draw random number to determine who goes first
        order = [c1, c2]
    else:
        order = [c2, c1]

    # Update nutrients and board
    for c in order:
        metabolism(c, enviro)  # update nutrient values
        U = np.real(np.fft.ifft2(c.K * np.fft.fft2(c.A)))
        c.A = np.clip(c.A + 1 / c.T * c.growth(U*min(c.nutrition, 1)), 0, 1)
        enviro.grid = np.clip(enviro.grid-c.A, 0, np.max(enviro.grid)) # delete from nutrient gradient


def run_one(c1, c2, enviro, verbose=False, return_means=False):
    """Run creature of given parameters in given environmental configuration until it dies.
    Return number of timesteps creature survived for before death.

    c1 = creature of time interest."""
    t = 0  # set timer
    global sums
    sums = np.zeros(10000)
    means = []
    while np.sum(c1.A) and (  # While main creature is still alive
            t < 10000):

        if t % enviro.dt == 0:
            enviro.update()  # Update enviro

        t += 1  # update timer by 1

        if verbose & (t % 1000 == 0):
            print(t)  # Show that it is working even after long waits

        update_man(c1, c2, enviro)  # Run update and show
        means.append(c1.running_mean)

    if return_means:
        return means
    else:
        return t # return time taken for main creature to die


def mutate_and_select(c1, c2, enviro, k=1, moving=False, runs=10, minus_keys=1, median=False):
    """Mutate one parameter from c1 and assess fitness of new solution agaisnt wild type
    in input environment with c2. Save winning parameters to Creature.

    Method involve running wild type and mutant over ten distinct obstacle environment, and
    summing the survival time of each."""
    wild_type = c1
    mutant = deepcopy(c1)
    ## Choose parameter at random and mutate in mutant_type
    x = np.random.randint(0, len(Creature.keys) - minus_keys)
    mutant.__dict__[Creature.keys[x]] = mutate(mutant.__dict__[Creature.keys[x]])

    # Run mutant and wild over runs number of obstacle configurations

    t_wild = np.zeros(runs)
    t_mutant = np.zeros(runs)

    for i in range(runs):
        enviro.initiate()
        wild_type.initiate(random=True)
        mutant.initiate(random=True)
        c2.initiate(random=True)
        t_wild[i] = run_one(wild_type, c2, enviro)
        t_mutant[i] = run_one(mutant, c2, enviro)

    # Record mean and variance of survival times
    if median:
        wild_mean = np.median(t_wild)
        mutant_mean = np.median(t_mutant)
    else:
        wild_mean = t_wild.mean()
        mutant_mean = t_mutant.mean()
    print(wild_mean)
    c1.time_log = np.append(c1.time_log, [[wild_mean, t_wild.var()]], axis=0)
    c1.theta_log = np.append(c1.theta_log, [c1.theta(dict=False)], axis=0)
    # Select winning parameter
    if selection(c1, wild_mean, mutant_mean, k=k):
        print("Accept mutation")
        c1.mutations += 1
        c1.update_theta(mutant)  # Update creature parameters
        return True
    else:
        print("Reject mutation")
        return False

def optimise(creature, obstacle, N=1000, k=1, seed=0, max_time = 100, name=None, minus_key=1, median=False):
    """Evolve two creatures from creatures independently- alternating which is mutated and selected first.

    Creature = Starting creature to mutate
    Obstacle = Environment to Creature is evolved in
    N = psuedo-population size
    fixation_mark = number of generations without mutation that specifies arrival at a global optima
    max_time = maximum time to run for before manually ending it"""
    global population_size
    population_size = N
    np.random.seed(seed)  # set seed

    c1 = deepcopy(creature)
    c2 = deepcopy(creature)

    max_time = max_time * 60  # Translate to seconds
    """Evolve until parameters become fixed over fixation number of generations"""
    gen = 0  # time_count

    # initiate fixation and mutation counts for all creatures
    fix = 0
    start = time.time()

    while (np.absolute(c1.dfit) > 1/N) and (np.absolute(c2.dfit) > 1/N) and ((time.time() - start) < max_time):  # both fixations need to be met
        print("c1 mutations: ", c1.mutations)
        print("c2 mutations: ", c2.mutations)
        # Mutate one and then the other in environment where they both exist
        if mutate_and_select(c1, c2, obstacle, k=k, minus_keys=minus_key, median=median) or mutate_and_select(c2, c1, obstacle, k=k, minus_keys=minus_key, median=median):
            fix = 0  # if either mutate, fix is reset
        else:
            fix += 1  # otherwise, fix is added to
        gen +=1
        if (c1.time_log[-1, 0] == 10000 and not c1.time_log[-1,1]) and (c2.time_log[-1, 0] == 10000 and not c2.time_log[-1, 1]):
            c1.dfit = 0  # if both max out, stop simulation
            c2.dfit = 0
        if gen > 600:
            c1.dfit = c1.time_log[-601:-1, 0].mean() - c1.time_log[-600:, 0].mean()
            c2.dfit = c2.time_log[-601:-1, 0].mean() - c2.time_log[-600:, 0].mean()


        # switch names
        temp = deepcopy(c1)
        c1 = deepcopy(c2)
        c2 = deepcopy(temp)



    print("Saving configuration...")
    """Save winning parameters and timelogs"""


    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")

    c1_survival_time = get_survival_time(c1, c2, obstacle, summary=False, median=median)
    c1.survival_mean = c1_survival_time.mean()
    c1.survival_var = c1_survival_time.var()

    with open(name+"_c1_survivalraw.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in c1_survival_time:
            csvwrite.writerow([i])


    c2_survival_time = get_survival_time(c2, c1, obstacle, summary=False, median=median)
    c2.survival_mean = c2_survival_time.mean()
    c2.survival_var = c2_survival_time.var()

    with open(name+"_c2_survivalraw.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in c2_survival_time:
            csvwrite.writerow([i])


    with open(name+"_c1_parameters.csv", "w") as f:
        csvwrite = csv.writer(f)
        for key in Creature.keys:
            csvwrite.writerow([key, c1.__dict__[key]])
        csvwrite.writerow(["mutations", c1.mutations])
        csvwrite.writerow(["survival_mean", c1.survival_mean])
        csvwrite.writerow(["survival_var", c1.survival_var])
        csvwrite.writerow(["absorption", c1.absorption])
        csvwrite.writerow(["fixation", fix])
        csvwrite.writerow(["generations", gen])
        csvwrite.writerow(["organism_count", 2])
        csvwrite.writerow(["l", 64])
        csvwrite.writerow(["correlation_length", obstacle.L])
        csvwrite.writerow(["correlation_time", obstacle.T])
        csvwrite.writerow(["seed", seed])
        csvwrite.writerow(["dfit", c1.dfit])
        csvwrite.writerow(["N", N])
        csvwrite.writerow(["k", k])

    with open(name+"_c2_parameters.csv", "w") as f:
        csvwrite = csv.writer(f)
        for key in Creature.keys:
            csvwrite.writerow([key, c2.__dict__[key]])
        csvwrite.writerow(["mutations", c2.mutations])
        csvwrite.writerow(["survival_mean", c2.survival_mean])
        csvwrite.writerow(["survival_var", c2.survival_var])
        csvwrite.writerow(["absorption", c2.absorption])
        csvwrite.writerow(["fixation", fix])
        csvwrite.writerow(["generations", gen])
        csvwrite.writerow(["organism_count", 2])
        csvwrite.writerow(["l", 64])
        csvwrite.writerow(["correlation_length", obstacle.L])
        csvwrite.writerow(["correlation_time", obstacle.T])
        csvwrite.writerow(["seed", seed])
        csvwrite.writerow(["dfit", c2.dfit])
        csvwrite.writerow(["N", N])
        csvwrite.writerow(["k", k])

    save(c1.time_log, name=name+"_c1_times.csv")
    save(c1.theta_log, name=name+"_c1_theta.csv")
    save(c2.time_log, name=name+"_c2_times.csv")
    save(c2.theta_log, name=name+"_c2_theta.csv")

    with open(name+"_c1_selection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for row in c1.selection_log:
            csvwrite.writerow([row])

    with open(name+"_c2_selection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for row in c2.selection_log:
            csvwrite.writerow([row])

    with open(name+"_c1_winningselection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for row in c1.s_winning:
            csvwrite.writerow([row])

    with open(name+"_c2_winningselection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for row in c2.s_winning:
            csvwrite.writerow([row])

    return 1

def get_survival_time(c1, c2, enviro=None, runs=10, summary=False, verbose=False, median=False):
    """Calculate average run time over seeded 10 configurations.
    Return mean and variance."""
    times = np.zeros(runs)
    for i in range(1, runs + 1):
            np.random.seed(i)
            c1.initiate(random=True)  # Reset grid
            c2.initiate(random=True)
            enviro.initiate()  # set obstacle
            times[i - 1] = run_one(c1, c2, enviro, verbose=verbose)

    if summary:
        if median:
            return np.median(times), times.var()
        else:
            return times.mean(), times.var()
    else:
        return times
