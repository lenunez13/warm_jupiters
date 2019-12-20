import rebound
import reboundx
from reboundx import constants
import numpy as np
import datetime
from astropy import constants as const
import pandas as pd

logpath = '/storage/home/len56/work/warm_jupiters/data/datalog.csv'
sapath = '/storage/home/len56/work/warm_jupiters/simulation_archive/'

years = 2.*np.pi
G = const.G.value
c = const.c.value
M_jup = const.M_jup.value
M_sun = const.M_sun.value
au = const.au.value
R_jup = const.R_jup.value
Q = 1.0e5 #tidal quality factor
k = 0.26 #tidal Love number

def update_log(tmax,Nout,GR,tides,epsilon,notes):
    start = datetime.datetime.now()
    tag = start.isoformat().replace('T','').replace(':','')\
                                             .replace('-','').replace('.','')[2:-5]
    df = pd.read_csv(logpath)
    df = df.append({'tag':tag,
                   'tmax': str(round(tmax/(2.*np.pi))),
                   'Nout':str(round(Nout)),
                   'GR':str(GR),
                   'tides':str(tides),
                   'epsilon':epsilon,
                   'notes':notes,
                   'runtime':'RUNNING...'},
                   ignore_index=True)
    df.to_csv(logpath,index=False)
    del df
    return tag,start

def makesim(three_body=True,inner_mass=True):
    """Makes a rebound.Simulation() of system HD 147018."""
    sim = rebound.Simulation()
    sim.add(m = 0.927)
    
    m_inner = 2.12*M_jup/M_sun
    if inner_mass == False:
        m_inner = 0.0
    else: pass
    sim.add(m = m_inner,
            a = 0.2389, #1.0,
            e = 0.4686, #0.9,
            inc = np.radians(35.614629), #np.radians(65.0),
            Omega = np.radians(0.0),
            omega = np.radians(66.0054), #np.radians(38.4),
            M = np.radians(0.698350))
    
    if three_body == True:
        sim.add(m = 6.59*M_jup/M_sun,
                a =1.923,
                e = 0.133,
                inc = np.radians(3.3853710), #0.0,
                Omega = np.radians(180.0),
                omega = np.radians(136.865), #np.radians(17.2),
                M = np.radians(-293.214))
    else: pass
    return sim

def makesim_past():
    """Makes a rebound.Simulation() of system past configuration
    of HD 147018, based on Fig. 4 of Dawson & Chiang (2014)
    https://arxiv.org/abs/1410.2604"""
    sim = rebound.Simulation()
    sim.add(m = 0.927)
    sim.add(m = 0.0, #2.12*M_jup/M_sun,
            a = 1.0,
            e = 0.9,
            inc = np.radians(65.0),
            Omega = np.radians(0.0),
            omega = np.radians(38.4),
            M = np.radians(0.698350))
    sim.add(m = 6.59*M_jup/M_sun,
            a =1.923,
            e = 0.133,
            inc = 0.0,
            Omega = np.radians(180.0),
            omega = np.radians(17.2),
            M = np.radians(-293.214))
    return sim
    
def runsim(sim,tmax,Nout,GR=False,tides=False,Q=1.0e5,epsilon=1e-9,notes='-'):
    sim.integrator='ias15'
    sim.ri_ias15.epsilon = epsilon
    sim.move_to_com()
    interval =  tmax/Nout
    ps = sim.particles
    rebx = reboundx.Extras(sim)
     
    if GR == True:
        gr = rebx.load_force("gr")
        rebx.add_force(gr)
        gr.params["c"] = constants.C  

    if tides == True:
        mod = rebx.load_operator("modify_orbits_direct2")
        rebx.add_operator(mod)
        tau = Q/(3.0*k)*((2.12*M_jup/M_sun)/ps[0].m)/(R_jup/au)**5.0
        ps[1].params["tau_e"] = tau
        
    tag,start = update_log(tmax,Nout,GR,tides,epsilon,notes)
    file = sapath+'sa'+tag+'.bin'
    sim.automateSimulationArchive(file,interval=interval,deletefile=True)
    
    sim.integrate(tmax,exact_finish_time=0)
    
    runtime = datetime.datetime.now() - start
    df = pd.read_csv(logpath)
    df.set_index('tag', inplace=True)
    df.loc[float(tag),'runtime'] = runtime.total_seconds()
    df.to_csv(logpath,index=True)
    del df
    
    return file
    
def calc_imut(inc1,inc2,Omega1,Omega2):
    deltaOmega = Omega1-Omega2
    cosimut = np.cos(inc1)*np.cos(inc2)+np.sin(inc1)*np.sin(inc2)*np.cos(deltaOmega)
    imut = np.arccos(cosimut)
    return np.array(imut)

def calc_deltapomega(pomega1,pomega2):
    deltapom = (pomega2-pomega1+2.*np.pi)%(2.*np.pi)
    return np.array(deltapom)


    
    