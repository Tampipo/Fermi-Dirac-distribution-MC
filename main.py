from config import *
import numpy as np
import math as mp


def init(): #initialize parameters
    Ex=(hbar*2*np.pi)**2/(2*me*Lx**2*kb*T) #dimensionless energy in x direction
    Ey=(hbar*2*np.pi)**2/(2*me*Ly**2*kb*T) #y direction
    Ez=(hbar*2*np.pi)**2/(2*me*Lz**2*kb*T) #z direction
    E_0 = min(Ex, Ey, Ez) 
    n_cut=-mp.log2(0.01)/E_0 #max number of states in a direction
    return Ex,Ey,Ez,E_0,n_cut

def init_states():
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

def proba(old_state, new_state, Ex,Ey,Ez, config_dict): #new_state is a list of the incoming numbers
    old_numbers=config_dict[f'{old_state}']
    return min(1,np.exp(-Ex(old_numbers[0]**2-new_state[0]**2)-Ey(old_numbers[1]**2-new_state[1]**2)-Ez(old_numbers[2]**2-new_state[2]**2)))

def main():
    print("Please select your parameters")
    N=input("Number of particles:")
    T=input("Temperature (K):")
    Lx=input("Box dimensions x (m):")
    Ly=input("Box dimensions y (m):")
    Lz=input("Box dimensions z (m):")
    print("Initializing parameters...")
    init_param=init()

if __name__=="__main__":
    print(36*hbar**2,me*kb)
    main()
