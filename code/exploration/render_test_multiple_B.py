# !/usr/bin/env python3

"""Render nutrient B on cluster"""

from two_nutrient_B import *
import os


def do_render(i):
    creatures = load_pairs("../results/multiple_B/")
    As = np.unique([c[0].absorption for c in creatures])

    orbium = [c for c in creatures if c[0].absorption == As[i]]  # filter by absorption
    print(len(orbium))

    Ls = [15, 20, 30]
    Ns = [100, 1000, 10000]

    dim = [3,3]

    labels = {"y": [np.arange(32, (64*dim[1])+32, 64), Ns],
              "ylab": "Population size (N)",
              "xlab": "Correlation length (L)",
            "x": [np.arange(32, (64*dim[0])+32, 64), Ls],
              "title": "Multiple B. Absorption "+str(As[i])
                       }

    render_grid(orbium, dim=[3,3], name="../results/test_multiple_B_absorption"+str(As[i]), **labels)

do_render(0)
do_render(1)
do_render(2)
