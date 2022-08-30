# !/usr/bin/env python3
"""House evolution of orbium in stochastic attractant distribution made up of
nutrient source B.

In nutrient source A, the orbium has a "nutrition" parameter, which is multiplied by U
just before the growth distribution is calculated. The nutrition parameter starts at 1, and it constantly
depleted by 0.01 each timestep. The orbium is able to recover nutrition by seeking nutrients in its environment.

Crucially, the orbium cannot eat too much- its only task is constantly seeking. The max nutrient which can be stored is 2,
giving the orbium 200 timesteps before it runs out and begins depleting."""
####### IMPORTS ########
from lenia_package import *
import os
import re
############################ ANALYTICS ############################
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
def get_a(file):
    """Return seed from input filename"""
    return int(re.search(r"\d_selection.csv", file).group().split("_selection.csv")[0])

def get_N(file):
    """Return seed from input filename"""
    return int(re.search(r"N\d+", file).group().split("N")[1])
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
        i += 100

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
############################ RENDERING #############################

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

    [do_render(T, N, "../results/nutrient_B/naive/") for N in Ns for T in Ts]

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
    if labels.get("title"):
        plt.title(labels.get("title"))

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

def render_grid(creatures, name=None, dim=[4, 3], seed=0, video=False, standard=False, initiate=True, frames=400, **labels):
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
        anim = animation.FuncAnimation(fig, update_grid_std, frames=frames, interval=20, save_count=frames)
    else:
        anim = animation.FuncAnimation(fig, update_grid, frames=frames, interval=20, save_count=frames)
    if name:
        writervideo = animation.FFMpegWriter(fps=60)
        if video:
            anim.save(name + "_seed" + str(seed) + ".mp4", writer=writervideo)
        else:
            anim.save(name + "_seed" + str(seed) + ".gif", writer=writervideo)
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
############################ EVOLUTION #############################
### RECORDS ###
#global time_log, theta_log, selection_coef
#time_log = pd.DataFrame(columns=["wild_mean", "wild_var"])
#theta_log = pd.DataFrame(columns=Creature.keys)
#selection_coef = []

def record_time(wild_mean, wild_var):
    """Record timeline"""
    global time_log
    x = pd.DataFrame([[wild_mean, wild_var]],
                     columns=["wild_mean", "wild_var"])  # record averages
    time_log = pd.concat([time_log, x])


def record_theta(creature):
    global theta_log
    theta_log = pd.concat([theta_log, pd.DataFrame(creature.theta())])


def run_one_stress_test(creature, enviro, snapshot=None, img_save=None):
    """Run creature of given parameters in given environmental configuration until it dies.
    Return number of timesteps creature survived for before death"""
    t = 0  # set timer
    global sums
    sums = np.zeros(10000)
    mass =[]
    volume = []
    centroid = [] # difference of centroid
    means = []
    grid = np.zeros([64,64])
    while np.sum(creature.A) and (
            t < 10000):  # While there are still cells in the learning channel, and timer is below cut off

        if t % enviro.dt == 0:
            enviro.update()  # Update enviro

        t += 1  # update timer by 1

        centroid.append(creature.centroid())
        mass.append(np.sum(creature.A))
        volume.append(np.sum(creature.A > 0))
        means.append(creature.running_mean)

        update_man(creature, enviro)  # Run update and show


        if t == snapshot:
            grid = [creature.A, enviro.grid]

    return (pd.DataFrame({"nutrient": means,
                         "centroid": centroid,
                         "mass": mass,
                         "volume": volume}), grid)

def stress_test(creature):
    mass = []
    volume = []
    nutrient = []
    density = []
    yspeed = []
    xspeed = []
    grids = []
    for i in range(3):
        creature.initiate()
        creature.enviro.initiate()
        temp = run_one_stress_test(creature, creature.enviro, snapshot=int(creature.dict["survival_mean"]/2))
        grids.append(temp[1])
        temp = temp[0]
        mass += list(temp.mass)
        volume += list(temp.volume)
        density += list(temp.mass/temp.volume)
        nutrient.append(temp.nutrient.iloc[-1])  # append last nutrient value
        xspeed.append(np.sum([i[0] for i in temp.centroid])/len(temp))  # average speed
        yspeed.append(np.sum([i[1] for i in temp.centroid])/len(temp))

    enviro_fill = np.zeros([64, 64*3])
    creature_fill = np.zeros([64, 64*3])
    enviro_grids = [i[1] for i in grids]
    creature_grids = [i[0] for i in grids]

    for i in range(3):
        n = i*64
        enviro_fill[0:64,n:n+64] = enviro_grids[i]
        creature_fill[0:64,n:n+64] = creature_grids[i]

    img_save = "../results/naive_B/snapshots/s"+str(int(creature.dict["seed"]))+"N"+str(int(creature.dict["N"]))+"T"+str(int(creature.enviro_T))+"L"+str(int(creature.enviro_L))+".png"
    plt.matshow(np.dstack([np.zeros([64, 64*3]), creature_fill*2, enviro_fill*0.5]))
    plt.suptitle("t="+str(int(creature.dict["survival_mean"]/2)))
    plt.savefig(img_save)
    plt.close()

    return {"mass": np.mean(mass),
                       "volume": [np.mean(volume)],
                       "density": [np.mean(density)],
                       "final_nutrition":[np.mean(nutrient)],
                       "xspeed": [np.mean(xspeed)],
                       "yspeed": [np.mean(yspeed)]}

def stress_test_all():

    df = pd.read_csv("../results/naive_nutrient_B.csv")
    df = df.sort_values("seed")
    names = [i.split("_parameters.csv")[0] for i in df.files]
    df["name"] = names

    dicts = [{colname: df[colname].iloc[i] for colname in df.columns} for i in range(len(df))]
    orbium = [Creature(filename=None, dict=dicts[i], cluster=True) for i in range(len(dicts))]
    df = pd.DataFrame(stress_test(orbium[0]))

    for i in range(1, len(orbium)):
        df = pd.concat([df, pd.DataFrame(stress_test(orbium[i]))])

    df.to_csv("../results/naive_A/naive_A_stress.csv")
    return df
### PROCESS ###
def overlap(creature, feature, sums=True):
    """Return sum of cells in environment which overlap with creature"""
    overlapping = sum([creature.A, feature.grid]) > feature.grid
    if sums:
        return np.sum(overlapping * feature.grid)
    else:
        return overlapping * feature.grid

def mutate(p):
    """Mutate input parameter p"""
    return np.exp(np.log(p) + np.random.uniform(low=-0.2, high=0.2))

def prob_fixation(wild_time, mutant_time, N, k=1):
    """Return probability of fixation given time survived by mutant and wild type,
    and psuedo-population size N and scaling coefficient k"""
    s = ((mutant_time - wild_time) / wild_time)*k  # selection coefficient
    selection_coef.append(s)  # record selection coefficient
    """If s is zero, there is no selective difference between types and probability of 
    fixation is equal to 1/populdation size"""
    if s:
        return (1 - np.exp(-2 * s)) / (1 - np.exp(-2 * N * s))
    else:  # if s is zero, prob is equal to 1/N
        return 1 / N

def selection(t_wild, t_mutant, k):
    """Return winning solution based on survival times of wild type (t_wild) and
    mutant type (t_mutant)."""
    pfix = prob_fixation(wild_time=t_wild, mutant_time=t_mutant, N=population_size, k=k)  # get probability of fixation
    if pfix >= np.random.uniform(0, 1):
        # ACCEPT MUTATION
        winning_s.append(selection_coef[-1])  # record winning coefficient
        return True
    else:
        # REJECT MUTATION
        return False

def update_man(creature, enviro):
    """Update learning channel by 1/T according to values in learning channel A,
    and environmental nutrient gradient. Update nutrition count according to constant decay and
    any amounts replenished through the environment"""

    creature.nutrition -= 0.0002*(np.sum(creature.A>0)**.75)  # depleted proportional to living cells
    creature.nutrition += overlap(creature, enviro) * creature.absorption  # add nutrients from food
    creature.nutrition = min(creature.nutrition, 1)  # max they can have is 1 - gives 40 timesteps to death

    # enviro.grid = np.clip(enviro.grid-creature.A, 0, np.max(enviro.grid)) # delete from enviro

    U = np.real(np.fft.ifft2(creature.K * np.fft.fft2(creature.A)))
    creature.A = np.clip(creature.A + 1 / creature.T * creature.growth(U*min(creature.nutrition, 1)), 0, 1)  # Grow based on nutrition levels

def update_naive(creature, enviro):
    enviro.grid = np.clip(enviro.grid-creature.A, 0, np.max(enviro.grid)) # delete from enviro
    U = np.real(np.fft.ifft2(creature.K * np.fft.fft2(creature.A)))
    creature.A = np.clip(creature.A + 1 / creature.T * creature.growth(U*min(creature.nutrition, 1)), 0, 1)  # Grow based on nutrition levels


def run_one(creature, enviro, show_after=0, moving=False, verbose=False, give_sums=False, overlaps = False, return_means=False):
    """Run creature of given parameters in given environmental configuration until it dies.
    Return number of timesteps creature survived for before death"""
    t = 0  # set timer
    global sums
    sums = np.zeros(10000)
    means = []
    o = []
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
            means.append(creature.nutrition)
            o.append(overlap(creature, enviro))
        # if t == show_after:
        #   plt.matshow(sum([creature.A, obstacle.grid]))
    if return_means:
        return means
    elif overlaps:
        return o
    else:
        return t

def mutate_and_select(creature, enviro, k=1, moving=False, runs=10, minus_keys=1, median=False):
    """Mutate one parameter from creature and assess fitness of new solution agaisnt wild type
    in input obstacle environment. Save winning parameters to Creature.

    Method involve running wild type and mutant over ten distinct obstacle environment, and
    summing the survival time of each."""
    wild_type = creature
    mutant = deepcopy(creature)

    ## Choose parameter at random and mutate in mutant_type
    x = np.random.randint(0, len(Creature.keys) - minus_keys)
    mutant.__dict__[Creature.keys[x]] = mutate(mutant.__dict__[Creature.keys[x]])

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
    if selection(wild_mean, mutant_mean, k=k):
        print("Accept mutation")
        creature.update_theta(mutant)  # Update creature parameters
        return True
    else:
        print("Reject mutation")
        return False

def optimise(creature, obstacle, N=1000, k=1, seed=0, max_time = 100, name=None, minus_key=1, median=False):
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
    dfit = 10 # set dfit
    while (np.absolute(dfit) > 1/N) and ((time.time() - start) < max_time):
        print("mutations: ", mutation)
        if mutate_and_select(creature, obstacle, k=k, minus_keys=minus_key, median=median):  # Updates creature values
            mutation += 1
            fix = 0  # reset everytime there is a mutation
        else:
            fix += 1
        gen += 1
        if gen >= 3:
            if (time_log.wild_mean.iloc[-3:].mean()== 10000) and not time_log.wild_mean.iloc[-3:].var():
                dfit = 0 # terminate when reach max out of 10000 over last three generations
        if gen > 600:  # after 600 generations
            dfit = time_log.wild_mean.iloc[-601:-1].mean() - time_log.wild_mean.iloc[-600:].mean()  # get difference of fitness over last 600 generations

        print(gen)


    print("Saving configuration...")
    """Save winning parameters and timelogs"""

    if name:
        creature.name = name
    else:
        creature.name = str(creature.n) + "_orbium_fix" + str(fix) + "s" + str(seed) + "N" + str(
            N)

    """Update survival time mean and variance by running over 10 configurations with seed"""
    print("Calculating survival means...")
    survival_time = get_survival_time(creature, obstacle, summary=False, median=median)

    with open(name+"_survivalraw.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in survival_time:
            csvwrite.writerow([i])

    creature.mutations = mutation
    creature.generations = gen
    creature.survival_mean = survival_time.mean()
    creature.survival_var = survival_time.var()
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
        csvwrite.writerow(["dfit", dfit])
        csvwrite.writerow(["N", N])
        csvwrite.writerow(["k", k])

    time_log.to_csv(creature.name + "_times.csv")  # Save timelog to csv
    theta_log.to_csv(creature.name + "_theta.csv")


    with open(creature.name + "_selection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in selection_coef:
            csvwrite.writerow([i])

    with open(creature.name + "_winningselection.csv", "w") as f:
        csvwrite = csv.writer(f)
        for i in winning_s:
            csvwrite.writerow([i])

    return 1

def get_survival_time(creature, enviro=None, runs=10, summary=False, verbose=False, median=False):
    """Calculate average run time over seeded 10 configurations.
    Return mean and variance."""
    times = np.zeros(runs)
    if enviro:
        for i in range(1, runs + 1):
            np.random.seed(i)
            creature.initiate()  # Reset grid
            enviro.initiate()  # set obstacle
            times[i - 1] = run_one(creature, enviro, verbose=verbose)
    else:
        for i in range(1, runs + 1):
            creature.initiate()
            times[i - 1] = run_one(creature, enviro, verbose=verbose)

    if summary:
        if median:
            return np.median(times), times.var()
        else:
            return times.mean(), times.var()
    else:
        return times


def test_nut(orbium, absorption=0.00003, means=False):
    enviro = orbium.enviro
    orbium.absorption = absorption
    orbium.initiate()
    enviro.initiate()
    return run_one(orbium, enviro, return_means=means)

def get_rings():
    df = load_vitals("../results/nutrient_B/naive/")
    df = df[df.correlation_time==1000]
    df = df[df.correlation_length==30]
    file = df.files.iloc[0]
    orbium = Creature("../results/nutrient_B/naive/"+file.split("_parameters.csv")[0], cluster=True)
    #run_one(orbium, orbium.enviro)

    return orbium

def test_n(a):
    orbium = Creature("orbium")
    orbium.absorption = a
    enviro = StochasticEnviro(L=30)
    means = []
    for i in range(10):
        orbium.initiate()
        enviro.initiate()
        means.append(run_one(orbium, enviro))
    return (np.mean(means), np.var(means))

"""orbium = load_all("../results/test_B/", sorting=False)

Ls = [15, 20, 30]
Ns =[100, 1000, 10000]
abs = [0.00008, 0.0001, 0.00012]

dim = [3,3]

for a in abs:
    sorted = []
    o = [i for i in orbium if i.absorption == a]
    for N in Ns:
        for L in Ls:
            sorted += [i for i in o if i.dict["N"] == N and i.enviro_L==L]

    labels = {"x": [np.arange(32, (64*dim[1])+32, 64), Ls],
              "xlab": "Correlation Length (L)",
              "ylab": "Population size (N)",
                "y": [np.arange(32, (64*dim[0])+32, 64), Ns],
              "title": "Nutrient B. Absorption="+str(a)
                       }
    render_grid(sorted, name="../results/test_B_absorption_"+str(a), dim=dim, **labels)
"""


"""
df = load_vitals("../results/test_B/")
Ns = [100, 1000, 10000]
Ls = [15, 20, 30]
fig, ax= plt.subplots(1,3)
for n in range(3):
    a = 5
    files = os.listdir("../results/test_B/")
    files = [i for i in files if get_L(i)==Ls[n]]
    f = [i for i in files if re.search(r"_selection.csv$", i) and get_N(i)==Ns[n]]
    f = [i for i in f if get_a(i)==a][0]
    mutations = df[df.seed==get_seed(f)].mutations.iloc[0]
    d = pd.read_csv("../results/test_B/"+f, header=None)
    ax[n].plot(d[0])
    ax[n].set_title("L="+str(Ls[n])+" mutations="+str(mutations))

plt.legend()
plt.xlabel("Generations")
plt.ylabel("Selection coefficient")
plt.suptitle("Selection coefficient over evolution of orbium in enviro T=1000, absorption=0.00008")
"""



