#!/usr/bin/env python3
import matplotlib.pyplot as plt

from two_nutrient_B import *

df = pd.read_csv("two_nutrient_B.csv")

def plot_by_N(N=10, df=df):
    d = df[df.N == N]
    d1 = d[d.creature=="c1"]
    d2 = d[d.creature=="c2"]
    plt.scatter(d1.seed, d1.survival_mean, label = "c1")
    plt.scatter(d2.seed, d2.survival_mean, label="c2")
    plt.legend()
    plt.title("Overview of two orbium survival times")
    plt.xlabel("Seed")
    plt.ylabel("Mean survival time")




