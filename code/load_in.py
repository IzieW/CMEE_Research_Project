# !/usr/bin/env python3
import matplotlib.pyplot as plt

from stochastic_nutrient_A import *

def read_raw(file):
    with open(file, "r") as f:
        csvread = csv.reader(f)
        return [float(i[0]) for i in csvread]

def load_raw_means(path):
    files = os.listdir(path)
    files = [i for i in files if re.search(r"survivalraw.csv$", i)]

    Ns =[]
    Ts=[]
    Ls=[]
    seeds=[]

    means = []

    for file in files:
        Ns += list(np.repeat(get_N(file), 10))
        Ts += list(np.repeat(get_T(file), 10))
        Ls += list(np.repeat(get_L(file), 10))
        seeds += list(np.repeat(get_seed(file), 10))
        means += read_raw(path+file)

    df = pd.DataFrame({"N": Ns,
                         "correlation_time": Ts,
                         "correlation_length":Ls,
                         "seed": seeds,
                         "survival_time": means})

    return df

def plot_means(df, explan, title=None):
    """Plot single bar plot with means across a given explanatory variable"""
    d = df.groupby(explan).describe().survival_time
    plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4)
    plt.xlabel(explan)
    plt.ylabel("mean survival time")
    if title:
        plt.title(title)

def plot_NbyT(df):
    Ns = np.unique(df.N)
    for n in Ns:
        d = df[df.N==n].groupby("correlation_time").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, label = "N="+str(n))

    plt.legend()
    plt.ylabel("Mean survival time")
    plt.xlabel("Correlation time (T)")
    plt.title("Survival time across environments with decreasing variation")

def plot_NbyL(df):
    Ns = np.unique(df.N)
    for n in Ns:
        d = df[df.N==n].groupby("correlation_length").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, label = "N="+str(n))

    plt.legend()
    plt.ylabel("Mean survival time")
    plt.xlabel("Correlation length (L)")
    plt.title("Survival time across environments with increasing nutrient availability")


def plot_TbyL(df):
    Ns = np.unique(df.correlation_time)
    for n in Ns:
        d = df[df.correlation_time==n].groupby("correlation_length").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, label = "T="+str(n))

    plt.legend()
    plt.ylabel("Mean survival time")
    plt.xlabel("Correlation length (L)")
    plt.title("Survival time across different stochastic nutrient environments")

def plot_LbyT(df):
    Ns = np.unique(df.correlation_length)
    for n in Ns:
        d = df[df.correlation_length==n].groupby("correlation_time").describe().survival_time
        plt.errorbar(d.index, d["mean"], yerr=d["std"]/np.sqrt(d["count"]), capsize=4, label = "L="+str(n))

    plt.legend()
    plt.ylabel("Mean survival time")
    plt.xlabel("Correlation time (T)")
    plt.title("Survival time across environments with decreasing variation")
