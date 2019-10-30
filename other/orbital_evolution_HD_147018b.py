#%matplotlib inline
import matplotlib.pyplot as plt
import os
import rebound
import reboundx
import numpy as np
import datetime
plt.switch_backend('agg')

Msun = 1988500 #1e24 kg
Mjup = 1898.19 #1e24 kg
Mj = Mjup/Msun #...Jupiter mass in Solar units
years = 2.*np.pi
AU = 149597870700 #meters per AU
Rj = 7.1492e7 #meters per Jupiter radius
Q = 1.0e5 #tidal quality factor
k = 0.26 #tidal Love number
Mp = 2.12*Mj
Mstar = 0.927

#...add HD 147018 parameters:
def simulation():
    sim = rebound.Simulation()
    sim.add(m=0.927)
    sim.add(a=0.238900,e=0.468600,m=0.0,omega=np.radians(66.0054),inc=np.radians(35.614629),Omega=0.,M=np.radians(0.698350))
    sim.add(a=1.9230,e=0.133000,m=0.0062886114,omega=np.radians(136.865),inc=np.radians(3.3853710),Omega=np.radians(180.0),M=np.radians(-293.214))
    return sim

def run_simulation(sim,times,tides=False):
    Nout = len(times)    
    e1 = np.zeros(Nout)
    
    sim.integrator='ias15'
    sim.ri_ias15.epsilon=0 #...to be used with fixed timestep in IAS
    sim.move_to_com()
    ps = sim.particles
    
    if tides == True:
        rebx = reboundx.Extras(sim)
        mod = rebx.load_operator("modify_orbits_direct2")
        rebx.add_operator(mod)
        tau = Q/(3*k)*(Mp/Mstar)/(Rj/AU)**5
        ps[1].params["tau_e"] = tau

    for i,time in enumerate(times):
        sim.integrate(time)
        e1[i] = ps[1].e
    return e1

times = np.linspace(0,1e6,1000)*years
sim = simulation()
e_tT = run_simulation(sim,times,tides=True)
sim = simulation()
e_tF = run_simulation(sim,times,tides=False)

now = datetime.datetime.now().isoformat().replace('-','').replace(':','').replace('.','')[:-6]
dir_name = '/storage/home/len56/work/warm_jupiters/'
savetag = 'hd147018_3-body.'+now
header = "HD 147018 3-body simulation with tides. Columns: time(years),e_tT,e_tF; where e = eccentricity inner planet, tT = tides on, tF = tides off"
np.savetxt(dir_name+'data/'+savetag+'.txt',np.c_[times/years,e_tT,e_tF],header=header)




