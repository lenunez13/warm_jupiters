import rebound
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import reboundx
from reboundx import constants
import math
import os
Mj = 955e-6 #...Jupiter mass in Solar units

#---HD 147018---#
def setup_system():
    sim = rebound.Simulation()
    #sim.add(m=1.0)
    sim.add(m=0.927)
    sim.add(a=0.238900,e=0.468600,m=0.0,omega=np.radians(66.0054),inc=np.radians(35.614629),Omega=0.,M=np.radians(0.698350))
    sim.add(a=1.9230,e=0.133000,m=0.0062886114,omega=np.radians(136.865),inc=np.radians(3.3853710),Omega=np.radians(180.0),M=np.radians(-293.214))
    #sim.move_to_com() # Moves to the center of momentum frame
    sim.move_to_hel()
    return sim

def sim_integrate(sim,time_stop,Nout=1000,addOrbitalEvolution=False,add_gr=False):
    sim.dt = sim.particles[1].P/60.
    #sim.integrator = "mercurius"
    sim.integrator = "ias15"
    times = np.linspace(0.,2.*np.pi*time_stop,Nout)
    pomega1 = np.zeros(Nout)
    pomega2 = np.zeros(Nout)
    e1 = np.zeros(Nout) 
    inc1 = np.zeros(Nout)
    inc2 = np.zeros(Nout)
    Omega1 = np.zeros(Nout) 
    Omega2 = np.zeros(Nout) 


        
    if add_gr == True:
        rebx = reboundx.Extras(sim)
        gr = rebx.load_force("gr")
        rebx.add_force(gr)
        gr.params["c"] = constants.C
    
        
    for i,t in enumerate(times):
        sim.move_to_hel()
        sim.integrate(t,exact_finish_time=0)
        pomega1[i] = sim.particles[1].pomega
        pomega2[i] = sim.particles[2].pomega
        e1[i] = sim.particles[1].e
        inc1[i] = sim.particles[1].inc
        inc2[i] = sim.particles[2].inc
        Omega1[i] = sim.particles[1].Omega
        Omega2[i] = sim.particles[2].Omega
        #if i%10 == 0:
        #    print(str(int((i+1)/10)) + "%")
        
    times = times/(2.*np.pi)
    
    return np.array([times,pomega1,pomega2,e1,inc1,inc2,Omega1,Omega2])
    
#...Vector-addition of inclination
def vector_addition(data):
    times,pomega1,pomega2,e1,inc1,inc2,Omega1,Omega2 = data
    deltaOmega = Omega1-Omega2
    cosimut = np.cos(inc1)*np.cos(inc2)+np.sin(inc1)*np.sin(inc2)*np.cos(deltaOmega)
    imut = np.degrees(np.arccos(cosimut))
    deltapom = np.degrees(pomega2-pomega1+2*np.pi)%360.
    
    return np.array([deltapom,imut])

tmax = 0.05*1e6
sim = setup_system()
data = sim_integrate(sim,tmax,add_gr=False)
sim = setup_system()
data_GR = sim_integrate(sim,tmax,add_gr=True)

times,pomega1,pomega2,e1,inc1,inc2,Omega1,Omega2 = data
times_GR,pomega1_GR,pomega2_GR,e1_GR,inc1_GR,inc2_GR,Omega1_GR,Omega2_GR = data_GR
deltapom,imut = vector_addition(data)
deltapom_GR,imut_GR = vector_addition(data_GR)


fig, [ax1, ax2, ax3] = plt.subplots(3, 1, sharex=True,figsize=(11,11)) 
plt.subplots_adjust(wspace=0, hspace=0.05)
ax = plt.gca()
ax.autoscale(enable=True, axis='x', tight=True)


#ax1.set_title(r'$M_{star} =$ '+str(Mstar),fontsize=20)
ax1.plot(times,imut_GR,color='black',label='GR')
ax1.plot(times,imut,label="Newtonian")
ax1.set_ylabel(r'$i_{\mathrm{mut}}$ (deg)',fontsize=20)

ax2.plot(times,e1_GR,color='black',label="GR")
ax2.plot(times,e1,label="Newtonian")
ax2.set_ylabel(r'$e$',fontsize=20)


ax3.plot(times,deltapom_GR,color='black',label="GR")
ax3.plot(times,deltapom,label="Newtonian")
ax3.set_ylabel(r'$\Delta\varpi$ (deg)',fontsize=20)
ax3.set_xlabel('Time (yr)',fontsize=20)
ax3.legend(loc="lower left",fontsize=16)


savetag = 'apsidal_evolution'
savetype = '.pdf'
plt.savefig(savetag+savetype,dpi=200,bbox_inches='tight')
plt.close()

