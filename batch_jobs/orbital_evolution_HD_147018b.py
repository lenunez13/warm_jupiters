import matplotlib.pyplot as plt
import os
import rebound
import reboundx
import numpy as np
import datetime
plt.switch_backend('agg')

#...Function to run REBOUND simulation and return inner planet's eccentricity
def run_simulation(times,add_tidal_migration=False):
    sim = rebound.Simulation()
    #...add HD 147018 parameters:
    sim.add(m=0.927)
    sim.add(a=0.238900,e=0.468600,m=0.0,omega=np.radians(66.0054),inc=np.radians(35.614629),Omega=0.,M=np.radians(0.698350))
    sim.add(a=1.9230,e=0.133000,m=0.0062886114,omega=np.radians(136.865),inc=np.radians(3.3853710),Omega=np.radians(180.0),M=np.radians(-293.214))
    sim.move_to_com() # Moves to the center of momentum frame
    #sim.move_to_hel()
    ps = sim.particles
    
    if add_tidal_migration == True:
        rebx = reboundx.Extras(sim)
        mod = rebx.load_operator("modify_orbits_direct2")
        rebx.add_operator(mod)
        ps[1].params["tau_e"] = 1

    Nout = len(times)    
    e1 = np.zeros(Nout)
    sim.integrator='ias15'
    sim.ri_ias15.epsilon=0 #...to be used with fixed timestep in IAS
    for i,time in enumerate(times):
        sim.integrate(time)
        e1[i] = ps[1].e
    return e1

#...Set time array and run simulation w/ and w/o tidal migration:
years = 2.0*np.pi
times = np.linspace(0,1e5,101)*years
e1_tm_f = run_simulation(times, add_tidal_migration=False)
e1_tm_t = run_simulation(times, add_tidal_migration=True)

#...Save results in txt file:
now = datetime.datetime.now().isoformat().replace('-','').replace(':','').replace('.','')[:-6]
dir_name = '/storage/home/len56/work/warm_jupiters/'
savetag = 'hd147018bc_tm.'+now
header = "Columns: time (years), eccentricity (tides off), eccentricity (tides on)"
np.savetxt(dir_name+'data/'+savetag+'.txt',np.c_[times/years,e1_tm_f,e1_tm_t],header=header)

#...Plot and save figure:
fig = plt.figure()
plt.plot(times/years,e1_tm_f,label='Tides off',color='black')
plt.plot(times/years,e1_tm_t,label='Tides on',linestyle='--',color='red')
plt.xlabel("Time (yr)" )
plt.ylabel("Eccentricity")
plt.title('HD 147018b')
plt.legend()
plt.savefig(dir_name+'plots/'+savetag+'.pdf',bbox_inches='tight')
plt.close()
