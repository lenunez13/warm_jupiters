import rebound
import numpy as np
import matplotlib.pyplot as plt
import reboundx
from reboundx import constants
import math
Mj = 955e-6 #...Jupiter mass in Solar units

#---HD 147018---#

def setup_system2(e0=0.468600,a0=0.238900,mp=0.00202983,mstar=0.927):
    sim = rebound.Simulation()
    sim.add(m=mstar)
    sim.add(a=a0,e=e0,m=mp,inc=np.radians(35.614629))#,omega=np.radians(66.0054),Omega=0.,M=np.radians(0.698350))
    #sim.move_to_com() # Moves to the center of momentum frame
    #sim.move_to_hel()
    return sim


force_is_velocity_dependent = False
def tidal_evolution2(reb_sim,rebx_effect,dt,timing):
    #...Get parameters from simulation
    sim = reb_sim.contents
    t = sim.t
    nyears = 2.*np.pi
    τ = 1.*nyears
    f = sim.particles[1].f
    n = sim.particles[1].n
    e = sim.particles[1].e
    a = sim.particles[1].a
    i = sim.particles[1].inc
    Ω = sim.particles[1].Omega
    ω = sim.particles[1].omega
    idot = 0.
    x = sim.particles[1].x
    y = sim.particles[1].y
    z = sim.particles[1].z
    xdot = sim.particles[1].vx
    ydot = sim.particles[1].vy
    zdot = sim.particles[1].vz
    
    #...Modularize recurrent formulas
    F0 = 1. - e**2.
    r = a*F0 / (1.+e*math.cos(f))
    rdot = n*a / math.sqrt(F0) * e*math.sin(f)
    rfdot = n*a / math.sqrt(F0) * (1. + e*math.cos(f)) 
    drde = (-2.*e*r / F0) - ((r**2. * math.cos(f)) / (a * F0))
    drdotde = rdot / (e * F0)
    drfdotde = rfdot*(e + math.cos(f))/(F0*(1. + e*math.cos(f)))
    F1 = drdotde * math.cos(ω+f) - drfdotde * math.sin(ω+f)
    F2 = drdotde * math.cos(i) * math.sin(ω+f) + drfdotde * math.cos(i) * math.cos(ω+f) - (idot*drde*z/r)
    F3 = (zdot + idot/i * z) * idot
    F4 = math.sin(i) * math.sin(ω+f)
    F5 = math.sin(i) * math.cos(ω+f)
    F6 = math.cos(i) * math.sin(ω+f)
    F7 = math.cos(i) * math.cos(ω+f)

    #...Calculate orbital evolution
    #edot = -(e0/τ)*np.exp(-t/τ)
    edot = (2.*np.pi*e0/τ)*np.cos(2.*np.pi*t/τ)
    adot = 2*a*e*edot/F0
    #adot = 0.
    xdotdot = (-xdot + 3.*z*idot*math.sin(Ω)) / (2.*a) * adot \
              + (math.cos(Ω)*F1 - math.sin(Ω)*F2)*edot \
              + math.sin(Ω)*F3
    ydotdot = (-ydot - 3.*z*idot*math.cos(Ω)) / (2.*a) * adot \
              + (math.sin(Ω)*F1 + math.cos(Ω)*F2) * edot \
              - math.cos(Ω) * F3
    zdotdot = ((-rdot/(2.*a)) * F4  + (-rfdot/(2.*a)) * F5 + (r/a*idot*F6)) * adot \
              + (drdotde*F4 + drfdotde*F5 + drde*idot*F6) * edot \
              + (rdot*F6 + rfdot*F7 + r*idot/i*F6 - idot*z) * idot
        
    sim.particles[1].ax += xdotdot
    sim.particles[1].ay += ydotdot
    sim.particles[1].az += zdotdot
    
    
#...Newtownian integration
year = 2.*np.pi # One year in units where G=1
def sim_integrate2(sim,time_stop,Nout=100,add_tidal_evolution=False):
    times = np.linspace(0.,time_stop*year,Nout)
    sim.integrator = "mercurius"
    e = np.zeros(Nout) 
    a = np.zeros(Nout)

    if add_tidal_evolution == True: 
        rebx = reboundx.Extras(sim)
        custom_effect = rebx.add_custom_force(tidal_evolution2, force_is_velocity_dependent)
        
    for i,t in enumerate(times):
        sim.integrate(t)
        e[i] = sim.particles[1].e
        a[i] = sim.particles[1].a
    times = times/(2.*np.pi)  
    return np.array([times,e,a])

def eEvol(e0,t,τ,evolType:str):
    if evolType == 'exp':
        e = e0*np.exp(-t/τ)
    elif evolType == 'trig':
        e = e0*(np.sin(2.*np.pi*t/τ) + 1)
    elif evolType == 'const':
        e = e0 * np.ones(len(t))
    return e

e0=0.468600
nyears = 2.*np.pi
τ = 1.
t_stop = 1e3*τ
sim = setup_system2(mp=0.0,e0=e0)
data = sim_integrate2(sim,t_stop,add_tidal_evolution=True)
times,e,a = data

Mstar = sim.particles[0].m
Mp = sim.particles[1].m

e0 = e[0]
e_theory = eEvol(e0,times,τ,'trig')
fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True,figsize=(11,11)) 
plt.subplots_adjust(wspace=0, hspace=0)
ax = plt.gca() 
ax.autoscale(enable=True, axis='x', tight=True)

ax1.set_title(r'$(M_{\mathrm{star}},M_{P})=($'+str(Mstar)+'$,$'+str(Mp)+'$)$',fontsize=20)
ax1.plot(times,e,label='Simulation')
ax1.set_ylabel(r'$e$ ',fontsize=20)
ax1.plot(times,e_theory,'--',label='Theory')
#ax1.set_ylim(e0*0.9,e0*1.1)
#ax1.axhline(y=1,color='black')
#ax1.axvline(x=5,color='black')
ax1.legend()

ax2.plot(times,a)
ax2.set_ylabel(r'$a$ (AU)',fontsize=20)
ax2.set_xlabel('time (yr)', fontsize=20)
#ax2.set_ylim(a0*0.9,a0*1.1);

plt.legend()

savetag = 'eEvol_exp_python'
savetype = '.png'

plt.savefig('plots/'+savetag+savetype,dpi=300,bbox_inches='tight')