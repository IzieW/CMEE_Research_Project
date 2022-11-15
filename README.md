# Sensing in self-organising systems: Exploring the evolutionary conditions for the emergence of "sensing" behaviours in Lenia cellular automata. 
This repo contains all the code and data involved in my Msci Computational Methods in Ecology and Evolution research project. 

![](https://github.com/IzieW/CMEE_Research_Project/blob/master/gifs/weak_sensing_B_464.gif)
*orbium (green) demonstrates "weak sensing" behaviour as it turns trajectory to maximise overlap with nutrient peaks (blue)*

Cellular automata (CA) represent promising candidates for studying the origins of life and mind as autonomous, self-assembling dynamical systems. Bert Chan's "Lenia" is a variety of smooth, semi-continuous cellular automata, capable of producing a wide range of interesting behaviours and morphologies. Recently, the flower's lab has used machine learning to produce Lenia life forms capable of avoidance behaviours. After optimisation, these novel artificial life forms are able to successfully navigate a field of obstacles. 


In this project I attempted to evolve the performance of rudimentary "seeking" or "sensing" behaviours in self-assembling artificial life forms. That is, the capability to swim towards a given food source.  

## What does an orbium eat? 
This project entail developing a novel "food" source. 

I fashioned a two-dimensional stochastic nutrient distribution, which produces nutrient peaks of different sizes (L) across the board which update stochastically over set time intervals(T)

Defining how the nutrient would interact with orbium to make it play the role of a nutrient source was more tricky- Orbium are already at a dynamic equilibrium without any resources from the environment. I attempted to contrive this type of dependence by making the variables in the orbium's parameter set dynamically dependent on interaction with the nutrient source in order to maintain its steady state. Thus if the orbium eats too much, or too little, they will spin out of equilibrium and die. 

This was fashioned in two ways: 

**NUTRIENT A** Orbium's growth mean is depleted each timestep and can be recovered subject to consuming nutrients in the environment. Orbium must consume just enough nutrients to maintain its growth mean within a range of viability to maintain equilibrium of cells. 


**NUTRIENT B** Orbium are endowed with a further "nutrition" parameter which is depleted each tiemstep and can be recovered subject to consuming nutrients. The nutrition parameter starts at 1 and each timestep is multiplied by the neighborhood sum of cells such that if the nutrition parameter is ever less than 1, growth of cells will slow expontentially until the orbium dies. 

The project found both A and B were able to successfully simulate the approximate role of nutrient to a living system. 

## Usage
### Languages
All scripts in the code directory were written in Python 3.8.10, Bash 5.0.17, and R4.1.1.
### Structure and Usage
The **code** repository contains all scripts which encode processes for optimisation, simulation and rendering of the creatures. 

**lenia_package.py** includes the main script for all functions relating to the simulation and rendering of orbium, as well as their nutrient environment. 

lenia_package.py is the main source script loaded in by all others. Modelled after the processes in Bert Chan's notebook, this script brings in a new architexture for dealing with the CA that defines creatures as a class of object, making it is easier to update, load in and maintain. The other python scripts then define unique evolutionary environments to run and optimise the orbium in.
