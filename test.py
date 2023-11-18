import matplotlib.pyplot as plt
import numpy as np
from labellines import labelLine, labelLines

temp=[0.01,0.05,0.1,0.5,1,10,100,1000,10000]


kb   = 1.38*10**(-23) # J.K^-1

def fermi_distrib(E,mu,T):
    return 1/(1+np.exp((E-mu)/(kb*T)))

def chemical_potential(T,Ef):
    mu=Ef*(1-np.pi**2/12*(kb*T/Ef)**2)
    return mu

Ef=3*10**(-18)

for i in range(0,len(temp)):
    T=temp[i]
    chem=chemical_potential(T,Ef)
    print(chem)
    energies=[0.01*Ef*k for k in range(0,200)]
    fermi_dirac=[fermi_distrib(E, chem, T) for E in energies]
    plt.plot(energies, fermi_dirac)

plt.legend()
plt.show()