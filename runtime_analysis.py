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
t_up = 6
times = np.logspace(0,t_up,t_up+1)*years
rt0000 = np.zeros(len(times))
rt1100 = np.zeros(len(times))
rt1110 = np.zeros(len(times))
rt1101 = np.zeros(len(times))
rt1111 = np.zeros(len(times))
epsilon=1e-6

def f(three_body,inner_mass,generalRelativity,tides,tmax,Nout,epsilon):
    runtime = nb.runsim(nb.makesim(three_body=three_body,inner_mass=inner_mass),
              tmax,Nout,generalRelativity=generalRelativity,tides=tides,epsilon=epsilon)[1] 
    return runtime
    
for i,tmax in enumerate(times):
    rt0000[i] = f(0,0,0,0,tmax,Nout,epsilon)
    rt1100[i] = f(1,1,0,0,tmax,Nout,epsilon)
    rt1110[i] = f(1,1,1,0,tmax,Nout,epsilon)
    rt1101[i] = f(1,1,0,1,tmax,Nout,epsilon)
    rt1111[i] = f(1,1,1,1,tmax,Nout,epsilon)
    
times = times/years

df = pd.DataFrame({'times (yr)':times,
                  'rt0000':rt0000,
                  'rt1100':rt1100,
                  'rt1110':rt1110,
                  'rt1101':rt1101,
                  'rt1111':rt1111})


start = datetime.datetime.now()
datetag = start.isoformat().replace('-','').replace(':','').replace('.','')[:-6]
savedir = 'data/'
savetag = 'rt.eps'+str(epsilon)+'.'+datetag
savetype = '.csv'
df.to_csv(savedir+savetag+savetype,index=False)