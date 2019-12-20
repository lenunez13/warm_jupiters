import numpy as np
from src import nbody as nb
years = 2.*np.pi

tmax = 1e5*years
Nout = 1e3
GR = False
tides = True
epsilon = 1e-9
notes = ''

sim = nb.makesim()
file = nb.runsim(sim,tmax=tmax,Nout=Nout,GR=GR,tides=tides,epsilon=epsilon,notes=notes)
