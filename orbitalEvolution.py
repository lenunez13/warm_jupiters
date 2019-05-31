import rebound
import numpy as np
import matplotlib.pyplot as plt
import reboundx
from reboundx import constants
import math
import datetime
from mpl_toolkits.mplot3d import Axes3D

#...Initialize global parameters corresponding to HD 147018b
global e0, a0, Mp, Mstar, A, tau0, inc0
e0 = 0.5
a0 = 0.2
Mp0 = 0.002
Mstar0 = 1.0
tau0 = 1.0
A0 = 1.0
inc0 = np.radians(35.614629)

"A function for the theoretical expectation of the time-evolution of eccentricity."
def eFxn(init_e,time,tau,A):
    #e = A*init_e*(np.sin(2.0*np.pi*time/tau) + 1.0)
    e = init_e*np.sin(time)
    return e

def edotovereFxn(init_e,time,tau,A,e):
    #edot = (2.0*np.pi*A*init_e/tau)*np.cos(2.0*np.pi*time/tau)/e
    edotovere = init_e*np.cos(time)/e
    return edotovere

def adotoveraFxn(e,edotovere):
    adotovera = 2.0*e**2.0*edotovere/(1.0 - e**2.0)
    return adotovera
    


def setupSystem(a,e,inc,Mstar,Mp):
    sim = rebound.Simulation()
    sim.add(m=Mstar)
    sim.add(a=a,e=e,m=Mp,inc=inc)#,omega=np.radians(66.0054),Omega=0.,M=np.radians(0.698350))
    #sim.move_to_com() # Moves to the center of momentum frame
    #sim.move_to_hel()
    return sim

force_is_velocity_dependent = True
def orbitalEvolution(reb_sim,rebx_effect,dt,timing):
    sim = reb_sim.contents
    t = sim.t #...time
    nyears = 2.0*np.pi
    
    #...Orbital parameters: REBOUND
    n = sim.particles[1].n #mean motion
    semi = sim.particles[1].a #semimajor axis
    e = sim.particles[1].e #eccentricity
    inc = sim.particles[1].inc #inclination 
    node = sim.particles[1].Omega #longitude of ascending node
    ω = sim.particles[1].omega #argument of pericenter
    pomega = sim.particles[1].pomega #longitude of pericenter
    f = sim.particles[1].f #true anomaly
    Manom = sim.particles[1].M #mean anomaly
    periq = semi*(1.0 - e) #perihelion distance (see mercury6_2.for)
    msum = sim.particles[0].m + sim.particles[1].m
    
    #...Position: REBOUND
    x = sim.particles[1].x
    y = sim.particles[1].y
    z = sim.particles[1].z
    d = sim.particles[1].d #...or see r below
    
    #...Velocity (Option 1): REBOUND
    vx = sim.particles[1].vx
    vy = sim.particles[1].vy
    vz = sim.particles[1].vz
    

    #...Setup eccentricity evolution
    edotovere = edotovereFxn(e0,t,tau0*nyears,A0,e)
    adotovera = adotoveraFxn(e,edotovere)
    idotoveri = 0.0

    if edotovere != 0.0:
        if e == 0.0:
            f = Manom
        else:
            temp = (periq*(1.0 + e)/(x**2.0+y**2.0+z**2.0)**(1.0/2.0) - 1.0)/e
            temp = math.copysign(min(abs(temp),1.0),temp)
            f = math.acos(temp)
            if math.sin(Manom) < 0.0:
                f = 2.*np.pi - f
    
        #...Orbital equations WDMC (A5 - A8)
        #r = semi*(1.0-e**2.0)/(1.0+e*math.cos(f)) #...or sim.particles[1].d
        r = d
        rdot = math.sqrt(msum)*semi**(-1.0/2.0)*e*math.sin(f)*(1.0-e**2.0)**(-1.0/2.0)
        rfdot = math.sqrt(msum)*semi**(-1.0/2.0)*(1.0+e*math.cos(f))*(1.0-e**2.0)**(-1.0/2)
        drde = -2.0*e*r/(1.0-e**2.0)-r**2.0*math.cos(f)/(semi*(1.0-e**2.0))
        drdotde = rdot/(e*(1.0-e**2.0))
        drfdotde = rfdot*(e+math.cos(f))/((1.0-e**2.0)*(1.0+e*math.cos(f)))
        peri = (pomega - node + 2.0*np.pi) % 2.0*np.pi
        
    #...Equations WDMC (A2); possibly to be moved inside/outside for-loop
    xdot = math.cos(node)*(rdot*math.cos(ω+f) - rfdot*math.sin(ω+f)) \
           - math.sin(node)*(rdot*math.cos(inc)*math.sin(ω+f) + rfdot*math.cos(inc)*math.cos(ω+f))
    ydot = math.sin(node)*(rdot*math.cos(ω+f) - rfdot*math.sin(ω+f)) \
           + math.cos(node)*(rdot*math.cos(inc)*math.sin(ω+f) + rfdot*math.cos(inc)*math.cos(ω+f))
    zdot = rdot*math.sin(inc)*math.sin(ω+f) + rfdot*math.sin(inc)*math.cos(ω+f)
    
    #...Velocity (Option 2): Orbital calculation
    #vx = xdot
    #vy = ydot
    #vz = zdot
    
    #...Update velocity: equations WDMC (A9 - A11)
    vux = x*adotovera 
    vuy = y*adotovera 
    vuz = z*adotovera 
    
    #...Update acceleration: equations WDMC (A12 - A14)
    ax = -vx*adotovera/2.0
    ay = -vy*adotovera/2.0
    az = -vz*adotovera/2.0

    if edotovere != 0.0 and e < 1.0 and e > 0.0:
        #...Update velocity: equations WDMC (A9 - A11)
        vux = vux + (r/semi - (1.0+e**2.0))*x*edotovere/(1.0-e**2.0)
        vuy = vuy + (r/semi - (1.0+e**2.0))*y*edotovere/(1.0-e**2.0)
        vuz = vuz + (r/semi - (1.0+e**2.0))*z*edotovere/(1.0-e**2.0)

        #...Update acceleration: equations WDMC (A12 - A14)
        ax = ax + (math.cos(node)*(drdotde*math.cos(ω+f) - drfdotde*math.sin(ω+f)) \
                - math.sin(node)*(drdotde*math.cos(inc)*math.sin(ω+f) \
                + drfdotde*math.cos(inc)*math.cos(ω+f)))*edotovere*e 
        ay = ay + (math.sin(node)*(drdotde*math.cos(ω+f) - drfdotde*math.sin(ω+f)) \
                + math.cos(node)*(drdotde*math.cos(inc)*math.sin(ω+f) \
                + drfdotde*math.cos(inc)*math.cos(ω+f)))*edotovere*e 
        az = az + (drdotde*math.sin(inc)*math.sin(ω+f) + drfdotde*math.sin(inc)*math.cos(ω+f))*edotovere*e

        sim.particles[1].ax += ax
        sim.particles[1].ay += ay
        sim.particles[1].az += az
        sim.particles[1].vx += vux
        sim.particles[1].vy += vuy
        sim.particles[1].vz += vuz

        
    else:
        sim.particles[1].ax += ax
        sim.particles[1].ay += ay
        sim.particles[1].az += az
        sim.particles[1].vx += vux
        sim.particles[1].vy += vuy
        sim.particles[1].vz += vuz    

#---INTEGRATION---#
nyears = 2.0*np.pi
def simIntegrate(sim,time_stop,Nout=100,addOrbitalEvolution=False):
    times = np.linspace(0.,time_stop*nyears,Nout)
    sim.dt = sim.particles[1].P/60.0 # small fraction of planetary period
    sim.integrator = "mercurius"
    e = np.zeros(Nout) 

    if addOrbitalEvolution == True: 
        rebx = reboundx.Extras(sim)
        custom_effect = rebx.add_custom_force(orbitalEvolution, force_is_velocity_dependent)
        
        for i,t in enumerate(times):
            sim.move_to_com()
            sim.integrate(t)
            e[i] = sim.particles[1].e  
    return np.array([times,e])


t_stop = 100
sim = setupSystem(a0,e0,inc0,Mstar0,Mp0)
data = simIntegrate(sim,t_stop,addOrbitalEvolution=True)
times,e = data

#---PLOTTING---#
eTheory = eFxn(e0,times,tau0,A0)
s=2
fig, ax1 = plt.subplots(1, 1, sharex=True,figsize=(6*s,4*s))
ax = plt.gca() 
ax.autoscale(enable=True, axis='x', tight=True)

ax1.plot(times,e,label='Simulation')
ax1.set_ylabel(r'$e$ ',fontsize=20)
ax1.set_xlabel('time (yr)', fontsize=20)
ax1.plot(times,eTheory,'--',label='Theory')
ax1.legend(loc="upper right");



datestamp = str(datetime.datetime.now()).replace(" ","_")
savetag = 'orbitalEvolution_'
savetype = '.png'
plt.savefig('plots/'+savetag+datestamp+savetype,dpi=100,bbox_inches='tight')