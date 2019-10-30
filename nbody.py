import rebound
import reboundx
from reboundx import constants
import numpy as np
import datetime
from astropy import constants as const

years = 2.*np.pi
G = const.G.value
c = const.c.value
M_jup = const.M_jup.value
M_sun = const.M_sun.value
au = const.au.value
R_jup = const.R_jup.value
Q = 1.0e5 #tidal quality factor
k = 0.26 #tidal Love number

def makesim():
    sim = rebound.Simulation()
    sim.add(m = 0.927)
    sim.add(m = 0.0,
            a = 0.2389, #1.0,
            e = 0.4686, #0.9,
            inc = np.radians(35.614629), #np.radians(65.0),
            Omega = np.radians(0.0),
            omega = np.radians(66.0054), #np.radians(38.4),
            M = np.radians(0.698350))
    sim.add(m = 6.59*M_jup/M_sun,
            a =1.923,
            e = 0.133,
            inc = np.radians(3.3853710), #0.0,
            Omega = np.radians(180.0),
            omega = np.radians(136.865), #np.radians(17.2),
            M = np.radians(-293.214))
    return sim
    
def runsim(sim,tmax,Nout,gr=False,tides=False,integrator='ias15'):
    sim_dir = 'simulation_archive/'
    start = datetime.datetime.now()
    datetag = start.isoformat().replace('-','').replace(':','').replace('.','')[:-6]
    file = sim_dir+'sa'+datetag+'.bin'
    print("Simulation start: "+start.isoformat()+"\n")
    print("INTEGRATING . . .\n")
    
    if integrator == 'ias15':
        sim.integrator='ias15'
        sim.ri_ias15.epsilon = 0

    if integrator == 'whfast':
        sim.integrator = "whfast"
        sim.ri_whfast.safe_mode = 0

    sim.move_to_com()
    interval =  tmax/Nout
    ps = sim.particles
    rebx = reboundx.Extras(sim)
    sim.automateSimulationArchive(file,interval=interval,deletefile=True)
    
    
    if gr == True:
        gr = rebx.load_force("gr")
        rebx.add_force(gr)
        gr.params["c"] = constants.C  

    if tides == True:
        mod = rebx.load_operator("modify_orbits_direct2")
        rebx.add_operator(mod)
        tau = Q/(3.0*k)*(ps[1].m/ps[0].m)/(R_jup/au)**5
        #tau = tmax*1.5
        ps[1].params["tau_e"] = tau
        
    sim.integrate(tmax,exact_finish_time=0)
    del sim
    end = datetime.datetime.now()
    runtime = end - start
    print("\n")
    print("Simulation end: "+end.isoformat())
    print("Runtime: "+str(runtime.total_seconds())+" s")
    print("Simulation archive: "+file)
    
def calc_imut(inc1,inc2,Omega1,Omega2):
    deltaOmega = Omega1-Omega2
    cosimut = np.cos(inc1)*np.cos(inc2)+np.sin(inc1)*np.sin(inc2)*np.cos(deltaOmega)
    imut = np.arccos(cosimut)
    return np.array(imut)

def calc_deltapomega(pomega1,pomega2):
    deltapom = (pomega2-pomega1+2.*np.pi)%(2.*np.pi)
    return np.array(deltapom)


    
    