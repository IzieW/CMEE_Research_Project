# !/usr/bin/env python3

from archived.nutrient_B_depletion import *

orbium = load_all("./", sorting=False)

Ls = [15, 20, 30]
Ns =[100, 1000, 10000]
abs = [0.01/20, 0.01/26, 0.01/35]

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
    render_grid(sorted, name="test_B_depletion_absorption_"+str(a), dim=dim, **labels)
