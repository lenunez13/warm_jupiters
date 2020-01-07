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

arg1=[10,100,1000] #tmax
arg2=[10,100,1000] #Nout
arg3=[False,True] #GR,tides
arg4=[1e-9,1e-6,1e-3] #epsilon


def mytestfunc(a_float, b_float, a_bool, c_float):
    print(a_float, b_float, a_bool, c_float)
   
   
parameterlist=[]
for i in arg1:
    for j in arg2:
        for k in arg3:
            for l in arg4:
                parameterlist.append([i,j,k,l])
