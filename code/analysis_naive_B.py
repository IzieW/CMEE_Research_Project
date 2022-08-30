#! usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

from stochastic_nutrient_B import *
from mpl_toolkits import mplot3d
import matplotlib

colours = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499',  '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00']

name = "../results/naive_B_solution_count.csv"
df = pd.read_csv("../results/naive_nutrient_B.csv")
times = pd.read_csv("../results/naive_B_means.csv")
del df["Unnamed: 0"]
del times["Unnamed: 0"]
df = df.replace(np.nan, "NA")

Ns = [10, 100, 1000, 10000, 100000]
Ts = [100, 1000, 10000]
Ls = [10, 15, 20, 30]


enviros = {"N": Ns,
           "correlation_length": Ls,
           "correlation_time": Ts}
"""Script takes full results - including stress test and qualitative analysis and creatures various
charts based on data"""


##### SOLUTION COUNTS #####

def counts_by_solution(df):
    df = df[df.qualitative != "NA"]
    N_list=[]
    L_list = []
    forms = []
    T_list = []
    count = []
    for form in df.qualitative.unique():
        d = df[df.qualitative==form]
        for N in Ns:
            da = d[d.N==N]
            for L in Ls:
                dat = da[da.correlation_length==L]
                for T in Ts:
                    data = dat[dat.correlation_time==T]
                    forms.append(form)
                    N_list.append(N)
                    L_list.append(L)
                    T_list.append(T)
                    count.append(len(data))

    return pd.DataFrame({
        "qualitative": forms,
        "N": N_list,
        "T": T_list,
        "L": L_list,
        "count": count
    })


def overview_forms(df):
    """produces summary chart counting the discrete morphological forms"""
    d = df.groupby("qualitative").describe()
    plt.bar(d.index, d.R["count"])
    plt.ylabel("Count")
    plt.xlabel("Discrete form")
    plt.title("Nutrient B: Overview of discrete forms after evolution")

def count_solutions(df, name=None):
    """Count unique solutions arrived at in each environment"""
    df=df[df.qualitative != "NA"]
    dict = {
        "N":[],
        "correlation_time": [],
        "correlation_length": [],
        "N_solutions": []
    }
    for N in Ns:
        d = df[df.N==N]
        for L in Ls:
            da = d[d.correlation_length==L]
            for T in Ts:
                dat = da[da.correlation_time==T]
                dict["N"].append(N)
                dict["correlation_length"].append(L)
                dict["correlation_time"].append(T)
                dict["N_solutions"].append(len(dat.qualitative.unique()))

    d = pd.DataFrame(dict)
    if name:
        d.to_csv(name)
    return d

def plot_solutions(df):
    """Plot number of solutions across different environments"""
    d = count_solutions(df)
    Ns = list(d.N.unique())
    fig, ax= plt.subplots(5, 1)
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colours)
    for i in range(len(Ns)):
        dat = d[d.N == Ns[i]]
        for T in d.correlation_time.unique():
            temp = dat[dat.correlation_time==T]
            ax[i].plot(temp.correlation_length, temp.N_solutions, label="T="+str(T))
            ax[i].set_title("Population size ="+str(Ns[i]))
    ax[0].legend()
    ax[4].set_xlabel("Correlation time (L)")
    ax[2].set_ylabel("Number of unique solutions")

def phase_counts(df):
    """L~T phase diagram of counts for each population size N"""


def plot_solutions_L(df):
    """plot number of solutions across L"""
    d = count_solutions(df)
    d = d.groupby("correlation_length").describe().N_solutions
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize = 4)
    plt.xlabel("Correlation length (L)")
    plt.ylabel("Number of diverse solutions")
    plt.title("Number of diverse solutions arrived at in environments with increasing nutrient abundance")


def plot_solutions_T(df):
    """plot number of solutions across T"""
    d = count_solutions(df)
    d = d.groupby("correlation_time").describe().N_solutions
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize = 4)
    plt.xlabel("Correlation time (T)")
    plt.ylabel("Number of diverse solutions")
    plt.title("Number of diverse solutions arrived at in environments with decreasing variation in environment")


def prepare_data(df):
    """Run all relevant functions upon data load in"""
    count_solutions(df)

########### Nutrient B analysis only ###########
# Main relationship between solutions by nutrient abundance (correlation_length)
def forms_by_L(df, forms_as_categories=True, show=True):
    dat = pd.DataFrame()
    for L in Ls:
        d = df[df.correlation_length == L]
        tot = np.sum(d.groupby("qualitative").describe().R["count"])
        percents = list((d.groupby("qualitative").describe().R["count"]/tot)*100)
        forms = list(d.groupby("qualitative").describe().index)
        if forms_as_categories:
            temp = pd.DataFrame()
            for i in range(len(forms)):
                temp[forms[i]] = [percents[i]]
            if show:
                print("L =" + str(L))
                print(temp)
        else:
            temp = pd.DataFrame({"qualitative": forms,
                                "percents": percents})
            temp["correlation_length"] = L
        dat = pd.concat([dat, temp])

    return dat

def forms_by_T(df, forms_as_categories=True, show=True):
    dat = pd.DataFrame()
    for T in Ts:
        d = df[df.correlation_time == T]
        tot = np.sum(d.groupby("qualitative").describe().R["count"])
        percents = list((d.groupby("qualitative").describe().R["count"]/tot)*100)
        forms = list(d.groupby("qualitative").describe().index)
        if forms_as_categories:
            temp = pd.DataFrame()
            for i in range(len(forms)):
                temp[forms[i]] = [percents[i]]
            if show:
                print("T=" + str(T))
                print(temp)
        else:
            temp = pd.DataFrame({"qualitative": forms,
                                "percents": percents})
            temp["correlation_time"] = T
        dat = pd.concat([dat, temp])

    return dat


def plot_forms_by_L(df):
    forms = list(df.groupby("qualitative").describe().index)
    d = forms_by_L(df, show=False)
    d = d.fillna(0)
    d["correlation_length"] = [str(i) for i in Ls]
    d.plot(x="correlation_length", kind="bar", stacked=True,
           title="Solution makeup across environments with increasing nutrient abundance",
           xlim="Correlation Length (L)", ylim="Percent of total outcomes")

def plot_forms_by_T(df):
    forms = list(df.groupby("qualitative").describe().index)
    d = forms_by_T(df, show=False)
    d = d.fillna(0)
    d["correlation_time"] = [str(i) for i in Ts]

    d.plot(x="correlation_time", kind="bar", stacked=True,
           title="Solution makeup across environments with decreasing levels of variation",
           xlim="Correlation Time (T)", ylim="Percent of total outcomes")


def filter_env(L, T, d):
    d = d[d.correlation_length == L]
    return d[d.correlation_time==T].N_solutions

def get_table_means(dat, y):
    z = np.ndarray([3,4])
    for i in range(3):
        d = dat[dat.correlation_time==Ts[i]]
        z[i,] = [d[d.correlation_length==L][y].mean() for L in Ls]

    return z

def get_table(N, y, d):
    if N:
        d = d[d.N==N]
    z = np.ndarray([3,4])
    if y == "N_solutions":
        for i in range(3):
            z[i,] = list(d[d.correlation_time==Ts[i]][y])
    else:
        for i in range(3):
            temp = d[d.correlation_time==Ts[i]]
            z[i,] = [d[d.correlation_length==L][y].mean() for L in Ls]

    return z

def plot_3D(df):
    d = count_solutions(df)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    ax.scatter3D(d.correlation_length, np.log10(d.correlation_time), d.N_solutions, c=d.N, marker="^")
    plt.legend()

def plot_counts_3D(df, z="N_solutions", means=True):
    """3D plot of L~T for a given variable z. If means is false, will plot for each population size N"""
    if z == "N_solutions":
        d = count_solutions(df)
    else:
        d=df
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    X,Y = np.meshgrid(Ls, Ts)
    if means:
        ax.plot_surface(np.log(X),  np.log10(Y), get_table_means(d, y=z), alpha = 0.75)
    else:
        for n in Ns:
            ax.plot_surface(np.log(X),  np.log10(Y), get_table(N=n, y=z, d=d), alpha = 0.75, label="N="+str(n))
        plt.legend()
    #ax.contourf(np.log(Ls),  np.log10(Ts), z, levels=500, alpha = 0.5)
    #ax.contourf(np.log(Ls),  np.log10(Ts), get_table(N=100,d=d), levels=500, alpha=0.5)
    ax.set_ylabel("Correlation Time (T)")
    ax.set_xlabel("Correlation Length (L)")
    ax.set_zlabel("Number of distinct solutions")


########## SURVIVAL MEANS ##########

# 1) Survival time by enviro
def plot_survival_by_L(times):
    i = 0
    for T in Ts:
        d = times[times.correlation_time==T]
        d = d.groupby("correlation_length").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, c=colours[i], label="T="+str(T))
        i +=1
    plt.legend()
    plt.xlabel("Correlation Length (L)")
    plt.ylabel("Mean survival time")
    plt.title("Mean survival time of orbium evolved in stochasticly distributed nutrient environments")

def plot_survival_by_T(times):
    i = 0
    for L in Ls:
        d = times[times.correlation_length==L]
        d = d.groupby("correlation_time").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, c=colours[i], label="L="+str(L))
        i +=1
    plt.legend()
    plt.xlabel("Correlation Time (T)")
    plt.ylabel("Mean survival time")
    plt.title("Mean survival time of orbium evolved in stochasticly distributed nutrient environments")

# 2) PARAMETERS BY ENVIRO
def plot_theta_by_enviro(theta, label= None, df=df, enviro=None, forms=False):
    """Plot parameters theta by enviro"""
    if not enviro and not forms:
        plt.scatter(df[theta], df.survival_mean, c=colours[0])
        if label:
            plt.xlabel(label)
        else:
            plt.xlabel(theta)
        plt.ylabel("Mean survival time")
    elif enviro:
        i=0
        for p in enviros[enviro]:
            d = df[df[enviro]==p]
            plt.scatter(d[theta], d.survival_mean, label=str(enviro)+"="+str(p), c=colours[i])
            i+=1
        plt.legend()
        if label:
            plt.xlabel(label)
        else:
            plt.xlabel(theta)
    else:
        i=0
        for form in df.qualitative.unique():
            d = df[df.qualitative==form]
            plt.scatter(d[theta], d.survival_mean, label=form, c=colours[i])
            i+=1
        plt.xlabel(label)
        plt.ylabel("Mean survival time")

def plot_theta_by_TL(theta, label=None, df=df):
    i=0
    for L in Ls:
        d = df[df.correlation_length == L]
        for T in Ts:
            d = df[df.correlation_time==T]
            plt.scatter(d[theta], d.survival_mean, label="L="+str(L)+", T="+str(T), c=colours[i])
            i +=1
    plt.legend()
    plt.xlabel(label)
    plt.ylabel("Mean survival time")


# 3) FORMS BY SURVIVAL TIME
def plot_survival_forms(df=df):
    """Plot survival mean of forms"""
    d = df.groupby("qualitative").describe().survival_mean
    plt.bar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel = "End solution"
    plt.ylabel = "Mean survival time"
    plt.title("Mean survival time across different forms")

def plot_survival_forms_byenviro(df=df):
    """Plot survival mean of forms"""
    d = df.groupby("qualitative").describe().survival_mean
    plt.bar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel = "End solution"
    plt.ylabel = "Mean survival time"
    plt.title("Mean survival time across different forms")


##### FORMS BY ENVIRO ######
def plot_forms_by_enviro(enviro, df=df):
    i=0
    forms = [i for i in list(df.qualitative.unique()) if type(i)==str]
    for form in forms:
        d = df[df.qualitative==form]
        d = d.groupby(enviro).describe().R["count"]
        plt.plot(d.index, d.values, label=form, c=colours[i])
        i+=1
    plt.legend()

def scatter_log_parameter_space(df=df, orbium=False, R=False):
    """3D plot of L~T for a given variable z. If means is false, will plot for each population size N"""
    d = df[df.qualitative != "NA"]
    cmap = matplotlib.cm.get_cmap('hsv')
    forms = d.qualitative.unique()
    fig = plt.figure(figsize = plt.figaspect(0.5))
    ax = fig.add_subplot(1,2,1, projection="3d")
    X,Y = np.meshgrid(d["T"], d.R)
    i=0
    #ax = plt.axes(projection="3d")
    ax = fig.add_subplot(1,2,1, projection="3d")
    for form in forms:
        dat = d[d.qualitative==form]
        ax.scatter3D(np.log(dat.m), np.log(dat.s),  np.log(dat.R), alpha = 0.75, label=form, c=colours[i])
        #ax.scatter3D(np.log(dat.m.mean()), np.log(dat.s.mean()),  np.log(dat["T"].mean()), alpha = 0.75, label=form, c=colours[i])
        i+=1
    if orbium:
        ax.scatter3D(np.log(0.15), np.log(0.0015), np.log(10), label="orbium", c="r")
    #ax.contourf(np.log(Ls),  np.log10(Ts), z, levels=500, alpha = 0.5)
    #ax.contourf(np.log(Ls),  np.log10(Ts), get_table(N=100,d=d), levels=500, alpha=0.5)
    ax.set_xlabel("log Growth mean (m)")
    ax.set_ylabel("log Growth std (s)")
    ax.set_zlabel("log R")

    ax=fig.add_subplot(1,2,2, projection="3d")
    i=0
    for form in forms:
        dat = d[d.qualitative==form]
        ax.scatter3D(np.log(dat.m), np.log(dat.s),  np.log(dat["T"]), alpha = 0.75, label=form, c=colours[i])
        i+=1
    if orbium:
        ax.scatter3D(np.log(0.15), np.log(0.0015), np.log(10), label="orbium", c="r")
    plt.legend(bbox_to_anchor=(1.05, 1))
    #ax.contourf(np.log(Ls),  np.log10(Ts), z, levels=500, alpha = 0.5)
    #ax.contourf(np.log(Ls),  np.log10(Ts), get_table(N=100,d=d), levels=500, alpha=0.5)
    ax.set_xlabel("log Growth mean (m)")
    ax.set_ylabel("log Growth std (s)")
    ax.set_zlabel("log T")
    plt.suptitle("Average parameter space across all solutions")
    plt.legend(bbox_to_anchor=(1.05, 1))


def plot_survival_forms_by_enviro(enviro, df=df):
    i=0
    forms = [i for i in list(df.qualitative.unique()) if type(i)==str]
    for form in forms:
        d = df[df.qualitative==form]
        d = d.groupby(enviro).describe().survival_mean
        plt.scatter(d.index, d["mean"], c=colours[i], s=5)
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), label=form, c=colours[i], capsize=4)
        i+=1
    plt.legend()

def stats_by_enviro(enviro, stat, df=df):
    """How does mass etc. vary with environment"""
    plt.scatter(df[enviro], df[stat])
    plt.xlabel(enviro)
    plt.ylabel(stat)

##### QUALITATIVE ANALYSIS OF FORMS  #####
def form_stats(stat, df=df):
    d = df.groupby("qualitative").describe()
    plt.errorbar(d.index, d[stat]["mean"], yerr=d[stat]["std"]/np.sqrt(d[stat]["count"]), capsize=4)
    plt.xlabel= "Qualitative form"
    plt.ylabel= stat
