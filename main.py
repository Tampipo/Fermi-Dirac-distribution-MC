import random
import numpy as np
import math as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

### Physics Constants ###

import numpy as np
import math as mp 
import random
import matplotlib.pyplot as plt

### Physics Constants ###

hbar = 6.02*10**(-34) # J.s
kb   = 1.38*10**(-23) # J.K^-1
me   = 9.11*10**(-31) # kg

def init(T,Lx,Ly,Lz,N): #initialize parameters
    Ex=(hbar*2*np.pi)**2/(2*me*Lx**2*kb*T) #dimensionless energy in x direction
    Ey=(hbar*2*np.pi)**2/(2*me*Ly**2*kb*T) #y direction
    Ez=(hbar*2*np.pi)**2/(2*me*Lz**2*kb*T) #z direction
    E_0 = min(Ex, Ey, Ez) 
    n_cut=min(10*N,int(-mp.log2(0.1)/E_0)) #max number of states in a direction
    return Ex,Ey,Ez,E_0,n_cut
    
def create_liste(N):
    """Crée la nouvelle liste qui définis l'ordre de tirage des N états"""
    L = np.array(range(N)) 
    random.shuffle(L)
    return(L)
    
def init_states(N):
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
        config_dict[f'{i}']=np.array([n_x,n_y,n_z,s]) #state is given by n_x,n_y,n_z,s
    return config_dict

def choose_new_state(dict_config, position,n_cut):

    test = 0 
    # elec = dict_config[f'{position}'] 
    while test == 0:
        n_x =  random.randrange(0, n_cut+1)
        n_y =  random.randrange(0, n_cut+1)
        n_z =  random.randrange(0, n_cut+1)
        s   =  2*random.randrange(0, 2) - 1
        print(s)
        test = 1
        l = dict_config.keys()
        for k in l:
            if k[0] == n_x and k[1] == n_y and k[2] == n_z and k[3] == s:
                test = 0
    return(position, np.array([n_x, n_y, n_z, s]))

def proba(old_state, new_state, Ex,Ey,Ez, config_dict): #new_state is a list of the incoming numbers
    old_numbers=config_dict[f'{old_state}']
    proba=min(1,np.exp(Ex*(old_numbers[0]**2-new_state[0]**2)+Ey*(old_numbers[1]**2-new_state[1]**2)+Ez*(old_numbers[2]**2-new_state[2]**2)))
    x=random.random()
    #print(proba, 'proba')
    print(old_numbers[0]**2-new_state[0]**2)
    if x<proba:
        config_dict[f'{old_state}']=new_state


def main():
    
    ask = int(input('0 pour choisir 1 par defaut : '))
    if ask == 0:
        print(type(me))
        print("Please select your parameters")
        N=int(input("Number of particles:"))
        T=float(input("Temperature (K):"))
        Lx=float(input("Box dimensions x (m):"))
        Ly=float(input("Box dimensions y (m):"))
        Lz=float(input("Box dimensions z (m):"))
        n_step = int(input("Number of step :"))
    else:
        N = 10
        T = 1
        Lx = 10**(-6)
        Ly = 10**(-6)
        Lz = 10**(-6)
        n_step = 1000
 
    print("Initializing parameters...")
    print(Lx)
    init_param=init(T,Lx,Ly,Lz,N)
    n_cut=init_param[-1]
    Ex=init_param[0]
    Ey=init_param[1]
    Ez=init_param[2]
    print("**********Simulation parameters*********")
    print("Number of particles: ", N)
    print("Temperature: ", T, "(K)")
    print("Box dimensions: ", Lx, Ly, Lz, "(m)" )
    print("States cut : ", n_cut)
    print("Energie along x",Ex)
    
    x = np.array([0])
    energie_moy = np.array([])
    
    config_dict = init_states(N)
    
    e = 0
    for cle, valeur in config_dict.items():
        e += mp.sqrt(valeur[0]**2+valeur[1]**2+valeur[2]**2)/N
    x = [0]
    energie_moy = [e]
    
    for k in range(n_step):
        liste = create_liste(N)
        for l in liste:
            old_state = l
            new_state = choose_new_state(config_dict, l, n_cut)[1]
            proba(old_state, new_state, Ex, Ey, Ez, config_dict)
            
        e = 0
        for cle, valeur in config_dict.items():
            e += mp.sqrt(valeur[0]**2+valeur[1]**2+valeur[2]**2)/N
        
        x += [k+1]
        energie_moy += [e]
        plt.plot(x, energie_moy)
        plt.show()
            
            

if __name__=="__main__":
    print(36*hbar**2,me*kb)
    main()
