import numpy as np

def calc_imut(inc1,inc2,Omega1,Omega2):
    deltaOmega = Omega1-Omega2
    cosimut = np.cos(inc1)*np.cos(inc2)+np.sin(inc1)*np.sin(inc2)*np.cos(deltaOmega)
    imut = np.arccos(cosimut)
    return np.array(imut)

def calc_deltapomega(pomega1,pomega2):
    deltapom = (pomega2-pomega1+2.*np.pi)%(2.*np.pi)
    return np.array(deltapom)