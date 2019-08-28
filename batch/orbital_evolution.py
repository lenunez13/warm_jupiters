import rebound
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import reboundx
from reboundx import constants
import math
import datetime
import matplotlib.gridspec as gridspec

#...Initialize global parameters corresponding to HD 147018b
global e0, a0, Mp0, Mstar0, A0, tau0, inc0, year
e0 = 0.5
a0 = 0.2
Mp0 = 0.0
Mstar0 = 1.0
tau0 = 1.0e5
A0 = 1.0
inc0 = np.radians(45.0)
year = 2.0*np.pi

"A function for the theoretical expectation of the time-evolution of eccentricity."
def e_fxn(init_e,time,tau,A):
    e = A*init_e*(np.sin(2.0*np.pi*time/tau) + 1.0)
    #e = A*init_e*(np.cos(time) + 1)
    #e = A*init_e*(np.sin(time) + 1)
    #e = init_e/2.0*time**2.0
    #e = init_e
    #e = A*init_e*np.exp(-time/tau)
    #e = A*time/tau + init_e
    return e

def edotovere_fxn(init_e,time,tau,A,e):
    edot = (2.0*np.pi*A*init_e/tau)*np.cos(2.0*np.pi*time/tau)
    #edot = -A*init_e*np.sin(time)
    #edot = A*init_e*np.cos(time)
    #edot = init_e*time
    #edot = 0.0
    #edot = -A*init_e/tau*np.exp(-time/tau)
    #edot = A/tau
    return edot/e

def adotovera_fxn(e,edotovere):
    adotovera = 2.0*e**2.0*edotovere/(1.0 - e**2.0)
    return adotovera
    


def setup_system(a,e,inc,Mstar,Mp):
    sim = rebound.Simulation()
    sim.add(m=Mstar)
    sim.add(a=a,e=e,m=Mp,inc=inc)#,omega=np.radians(66.0054),Omega=0.,M=np.radians(0.698350))
    sim.move_to_com()
    return sim

#sim = setup_system(a0,e0,inc0,Mstar0,Mp0)
#fig = rebound.OrbitPlot(sim,slices=True,color=True,unitlabel="[AU]")

def orbitalEvolution(reb_sim,rebx_effect,dt,timing):
    sim = reb_sim.contents
    t = sim.t #...time
    ps = sim.particles

    ax = 0.0
    ay = 0.0
    az = 0.0
    vux = 0.0
    vuy = 0.0
    vuz = 0.0
    
    #...Orbital parameters: REBOUND
    n = ps[1].n #mean motion
    semi = ps[1].a #Semimajor axis
    ecc = ps[1].e #Eccentricity
    inc = ps[1].inc #inclination 
    Ω = ps[1].Omega #longitude of ascending node
    ω = ps[1].omega #argument of pericenter
    pomega = ps[1].pomega #longitude of pericenter
    tru_anom = ps[1].f #True Anomaly
    Manom = ps[1].M #mean anomaly
    periq = semi*(1.0 - ecc) #perihelion distance (see mercury6_2.for)
    msum = ps[0].m + ps[1].m
    
    #...Position and velocity: REBOUND
    x = ps[1].x
    y = ps[1].y
    z = ps[1].z
    vx = ps[1].vx
    vy = ps[1].vy
    vz = ps[1].vz
    
    #...Setup eccentricity evolution
    edotovere = edotovere_fxn(e0,t,tau0*year,A0,ecc)
    adotovera = 0.0

    if edotovere != 0.0:
        
        #[Mercury definition of true anomaly would be inserted here.]

        #...Orbital equations WDMC (A5 - A8)
        r = semi*(1.0-ecc**2.0)/(1.0+ecc*math.cos(tru_anom)) #...or ps[1].d
        rdot = math.sqrt(msum)*semi**(-1.0/2.0)*ecc*math.sin(tru_anom)*(1.0-ecc**2.0)**(-1.0/2.0)
        rfdot = math.sqrt(msum)*semi**(-1.0/2.0)*(1.0+ecc*math.cos(tru_anom))*(1.0-ecc**2.0)**(-1.0/2)
        drde = -2.0*ecc*r/(1.0-ecc**2.0)-r**2.0*math.cos(tru_anom)/(semi*(1.0-ecc**2.0))
        drdotde = rdot/(ecc*(1.0-ecc**2.0))
        drfdotde = rfdot*(ecc+math.cos(tru_anom))/((1.0-ecc**2.0)*(1.0+ecc*math.cos(tru_anom)))
        #[Mercury definition of argument of periastron would be inserted here.]
        
   #...end if
    #...Equations WDMC (A2); possibly to be moved inside/outside for-loop
    xdot = math.cos(Ω)*(rdot*math.cos(ω+tru_anom) - rfdot*math.sin(ω+tru_anom)) - math.sin(Ω)*(rdot*math.cos(inc)*math.sin(ω+tru_anom) + rfdot*math.cos(inc)*math.cos(ω+tru_anom))
    ydot = math.sin(Ω)*(rdot*math.cos(ω+tru_anom) - rfdot*math.sin(ω+tru_anom)) + math.cos(Ω)*(rdot*math.cos(inc)*math.sin(ω+tru_anom) + rfdot*math.cos(inc)*math.cos(ω+tru_anom))
    zdot = rdot*math.sin(inc)*math.sin(ω+tru_anom) + rfdot*math.sin(inc)*math.cos(ω+tru_anom)
    

    if edotovere != 0.0 and ecc < 1.0 and ecc > 0.0:
        #...Update acceleration: equations WDMC (A12 - A14)
        ax += (drdotde*math.cos(ω+tru_anom) - drfdotde*math.sin(ω+tru_anom))*edotovere*ecc  
        ay += (drdotde*math.cos(inc)*math.sin(ω+tru_anom) + drfdotde*math.cos(inc)*math.cos(ω+tru_anom))*edotovere*ecc 
        az += (drdotde*math.sin(inc)*math.sin(ω+tru_anom) + drfdotde*math.sin(inc)*math.cos(ω+tru_anom))*edotovere*ecc
        ps[1].ax += ax.real
        ps[1].ay += ay.real
        ps[1].az += az.real

        #...Update velocity: equations WDMC (A9 - A11)
        vux += (r/semi - 1.0 - ecc**2.0)*x/(1.0-ecc**2.0)*edotovere
        vuy += (r/semi - 1.0 - ecc**2.0)*y/(1.0-ecc**2.0)*edotovere
        vuz += (r/semi - 1.0 - ecc**2.0)*z/(1.0-ecc**2.0)*edotovere
        ps[1].vx += vux.real
        ps[1].vy += vuy.real
        ps[1].vz += vuz.real
        


        #sim.integrator_synchronize()



            #sim.ri_whfast.recalculate_jacobi_this_timestep = 1
        #...end if 
   
#     else:
#         #sim.integrator_synchronize()
#         ps[1].vx += vux
#         ps[1].vy += vuy
#         ps[1].vz += vuz
#         ps[1].ax += ax
#         ps[1].ay += ay
#         ps[1].az += az

#         #sim.ri_whfast.recalculate_jacobi_this_timestep = 1

#---INTEGRATION---#
def simulation_integration(sim,integration_time,Nout=100,add_orbital_evolution=False):
    times = np.linspace(0.0,integration_time*year,Nout)
    ps = sim.particles
    #sim.dt = ps[1].P/60.0
    sim.integrator = "ias15"
    ecc = np.zeros(Nout)
    edotovere = np.zeros(Nout)
    dist = np.zeros(Nout)
    posx = np.zeros(Nout)
    posy = np.zeros(Nout)
    posz = np.zeros(Nout)
    accelx = np.zeros(Nout)
    accely = np.zeros(Nout)
    accelz = np.zeros(Nout)
    velx = np.zeros(Nout)
    vely = np.zeros(Nout)
    velz = np.zeros(Nout)
    if add_orbital_evolution == True: 
        rebx = reboundx.Extras(sim)
        custom_effect = rebx.add_custom_force(orbitalEvolution, force_is_velocity_dependent=True)
        custom_effect.force_is_velocity_dependent = 1
    print("INTEGRATION PROGRESS:")
    for i,time in enumerate(times):
        if i%1 == 0: print(str(i)+"%", end=" ")
        sim.integrate(time)
        ecc[i] = ps[1].e
        edotovere[i] = edotovere_fxn(e0,time,tau0,A0,ps[1].e)
        dist[i] = ps[1].d
        posx[i] = ps[1].x
        posy[i] = ps[1].y
        posz[i] = ps[1].z
        accelx[i] = ps[1].ax
        accely[i] = ps[1].ay 
        accelz[i] = ps[1].az            
        velx[i] = ps[1].vx
        vely[i] = ps[1].vy
        velz[i] = ps[1].vz

    return np.array([times,ecc,edotovere,dist,posx,posy,posz,accelx,accely,accelz,velx,vely,velz])


sim = setup_system(a0,e0,inc0,Mstar0,Mp0)
integration_time = 1e4 #...years
results = simulation_integration(sim,integration_time,Nout=100,add_orbital_evolution=True)
times,ecc,edotovere,dist,posx,posy,posz,accelx,accely,accelz,velx,vely,velz = results
print(" ")
print("INTEGRATION COMPLETE.")

#---PLOTTING---#
ecc_theory = e_fxn(e0,times/year,tau0,A0)
speed = np.sqrt(velx**2.0+vely**2.0+velz**2.0)
error = (ecc - ecc_theory)/ecc_theory

s=2
fig, [ax1,ax2] = plt.subplots(2, 1, sharex=True,figsize=(6*s,4*s))
plt.subplots_adjust(wspace=0, hspace=0.1)
ax = plt.gca() 
ax.autoscale(enable=True, axis='x', tight=True)
linewidth=3.0
linewidth=linewidth
legend_size=16
labelsize=12
fontsize=16


ax1.plot(times/year,ecc,label=r'$e_{sim}$',linewidth=linewidth)
ax1.plot(times/year,ecc_theory,'--',label=r'$e_{\mathrm{theory}}$',color='black',linewidth=linewidth)
ax1.set_ylabel('Eccentricity',fontsize=fontsize)
ax1.legend(prop={'size': legend_size})
ax1.tick_params(axis="y", labelsize=labelsize)



ax2.plot(times/year,error,linewidth=linewidth,color='red')
ax2.set_ylabel('Error',fontsize=fontsize)
ax2.set_xlabel('Time (yr)',fontsize=fontsize)
ax2.tick_params(axis="y", labelsize=labelsize)
ax2.tick_params(axis="x", labelsize=labelsize)
ax2.set_ylim(0,.05)




datestamp = str(datetime.datetime.now()).replace(" ","_")
savetag = 'orbital_evolution_'
savetype = '.pdf'
print("SAVING PLOT ...")
plt.savefig('plots/'+savetag+datestamp+savetype,dpi=200,bbox_inches='tight')
print("DONE.")

