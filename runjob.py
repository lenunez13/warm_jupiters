import numpy as np
from src import nbody as nb

generalRelativity = True
tides = False
epsilon = 1e-9
years = 2.*np.pi
tmax = 0.05e6*years
Nout = 1e3
sim = nb.makesim()
file = nb.runsim(sim,tmax,Nout,generalRelativity=generalRelativity,tides=tides,epsilon=epsilon)[0]

generalRelativity = True
tides = False
epsilon = 1e-6
years = 2.*np.pi
tmax = 0.05e6*years
Nout = 1e3
sim = nb.makesim()
file = nb.runsim(sim,tmax,Nout,generalRelativity=generalRelativity,tides=tides,epsilon=epsilon)[0]