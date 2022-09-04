# Sensing in self-organising systems: Exploring the evolutionary conditions for the emergence of "sensing" behaviours in Lenia cellular automata. 
This repo contains all the code and data involved in my Msci Computational Methods in Ecology and Evolution research project. 

In this project I attempted to evolve lenia creatures capable of seeking food. The main finding? Evolving intelligent life is difficult. 

## What does an orbium eat? 
I fashioned a two-dimensional stochastic nutrient distribution modelled after the one used in Godany et al.'s work on the evolution of different chemotactic strategies. This produces various nutrient peaks across the board which update stochastically over time. 

Setting how the nutrient would interact with orbium to make it play the role of a nutrient source was more tricky- Orbium are already at a dynamic equilibrium without any resources from the environment. I took the parts that interact in the equilibrium and tried to distribute them throughout the nutrient gradient, such that the orbium was dependent on the right kind of exchanges with the nutrient source in order to maintain its steady state. 

This was largely shooting in the dark... and there is definitely no one way to configure this, so I developed two different nutrient types to test interactions at different levels of complexity. 

**NUTRIENT A** Orbium's growth mean is depleted each timestep and can be recovered subject to consuming nutrients in the environment. Orbium must consume just enough nutrients to maintain its growth mean within a range of viability to maintain equilibrium of cells. 


**NUTRIENT B** Orbium are endowed with a further "nutrition" parameter which is depleted each tiemstep and can be recovered subject to consuming nutrients. The nutrition parameter starts at 1 and each timestep is multiplied by the neighborhood sum of cells such that if the nutrition parameter is ever less than 1, growth of cells will slow expontentially until the orbium dies. 

## Usage
### Languages
All scripts in the code directory were written in Python 3.8.10, Bash 5.0.17, and R4.1.1.
### Structure and Usage
The **code** repository contains all scripts which encode processes for optimisation, simulation and rendering of the creatures. 

lenia_package.py is the main source script loaded in by all others. Modelled after the processes in Bert Chan's notebook, this script brings in a new architexture for dealing with the CA that defines creatures as a class of object, making it is easier to update, load in and maintain. The other python scripts then define unique evolutionary environments to run and optimise the orbium in.
