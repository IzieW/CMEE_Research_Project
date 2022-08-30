# Sensing in self-organising systems: Exploring the evolutionary conditions for the emergence of "sensing" behaviours in Lenia cellular automata. 
This repo contains all the code and data involved in my Msci Computational Methods in Ecology and Evolution research project. 

In this project I attempted to evolve lenia creatures capable of seeking food. The main finding? Evolving intelligent life is difficult. 

## What does an orbium eat? 
I fashioned a stochastic nutrient gradient which updates as set intervals overtime. 

## Usage
### Languages
All scripts in the code directory were written in Python 3.8.10, Bash 5.0.17, and R4.1.1.
### Structure and Usage
The **code** repository contains all scripts which encode processes for optimisation, simulation and rendering of the creatures. 

lenia_package.py is the main source script loaded in by all others. Modelled after the processes in Bert Chan's notebook, this script brings in a new architexture for dealing with the CA that defines creatures as a class of object, making it is easier to update, load in and maintain. The other python scripts then define unique evolutionary environments to run and optimise the orbium in.
