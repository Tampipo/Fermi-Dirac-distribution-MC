from config import *
import numpy as np
import math as mp
print(N)

def init(): #initialize parameters
    Ex=(hbar*2*np.pi)**2/(2*me*Lx**2*kb*T) #dimensionless energy in x direction
    Ey=(hbar*2*np.pi)**2/(2*me*Ly**2*kb*T) #y direction
    Ez=(hbar*2*np.pi)**2/(2*me*Lz**2*kb*T) #z direction
    E_0 = min(Ex, Ey, Ez) 
    n_cut=-mp.log2(0.01)/E_0 #max number of states in a direction
    config_dict={} #initialize particle states
    n_x=0
    n_y=0
    n_z=0
    s=-1
    for i in range (0,N):
        if s==-1:
            s=1
        elif n_x>n_y:
            n_y+=1
            s=-s
        elif n_y>n_z:
            n_z+=1
            s=-s
        else :
            n_x+=1
            s=-s
        config_dict[f'{i}']=[n_x,n_y,n_z,s] #state is given by n_x,n_y,n_z,s
    return config_dict


def main():
    print(init_states())

if __name__=="__main__":
    print(36*hbar**2,me*kb)
    main()
