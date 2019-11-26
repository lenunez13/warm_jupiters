import rebound
import reboundx
from reboundx import constants
import numpy as np
import datetime
from astropy import constants as const
import matplotlib.pyplot as plt
from src import nbody as nb
import pandas as pd

Nout = 1e3
years = 2.*np.pi
t_up = 5
times = np.logspace(0,t_up,t_up+1)*years
eps = [0,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
rt1111 = np.zeros((len(eps),len(times)))


def f(three_body,inner_mass,generalRelativity,tides,tmax,Nout,epsilon):
    runtime = nb.runsim(nb.makesim(three_body=three_body,inner_mass=inner_mass),
              tmax,Nout,generalRelativity=generalRelativity,tides=tides,epsilon=epsilon)[1] 
    return runtime


for epsilon in eps:
    for i,tmax in enumerate(times):
        rt1111[epsilon][i] = f(1,1,1,1,tmax,Nout,epsilon)
    
times = times/years

df = pd.DataFrame({'times (yr)':times,
                  'rt1111(eps=0)':rt1111[0],
                  'rt1111(eps=1e-9)':rt1111[1],
                  'rt1111(eps=1e-8)':rt1111[2],
                  'rt1111(eps=1e-7)':rt1111[3],
                  'rt1111(eps=1e-6)':rt1111[4],
                  'rt1111(eps=1e-5)':rt1111[5],
                  'rt1111(eps=1e-4)':rt1111[6],
                  'rt1111(eps=1e-3)':rt1111[7],
                  'rt1111(eps=1e-2)':rt1111[8],
                  'rt1111(eps=1e-1)':rt1111[9],
                  'rt1111(eps=1)':rt1111[10],})


start = datetime.datetime.now()
datetag = start.isoformat().replace('-','').replace(':','').replace('.','')[:-6]
savedir = 'data/'
savetag = 'rt.eps.'+datetag
savetype = '.csv'
df.to_csv(savedir+savetag+savetype,index=False)