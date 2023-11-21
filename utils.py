import random
import numpy as np
import math as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D

### Physics Constants ###

hbar = 6.02*10**(-34) # J.s
kb   = 1.38*10**(-23) # J.K^-1
me   = 9.11*10**(-31) # kg

### Useful functions ###

def fermi_distrib(E,mu):
    return 1/(1+np.exp(E-mu))

def chemical_potential_low(Tnorm):
    mu=(1-np.pi**2/12*(Tnorm)**2)
    return mu

def chemical_potential_high(Tnorm):
    mu=3/2*Tnorm*np.log(4*np.pi/(6*np.pi**2)**(2/3)*(1/Tnorm))
    return mu

def init(T,Lx,Ly,Lz): #initialize parameters
    Ex=(hbar*2*np.pi)**2/(2*me*Lx**2*kb*T) #dimensionless energy in x direction
    Ey=(hbar*2*np.pi)**2/(2*me*Ly**2*kb*T) #y direction
    Ez=(hbar*2*np.pi)**2/(2*me*Lz**2*kb*T) #z direction
    return Ex,Ey,Ez
    
def create_liste(N):
    """Crée la nouvelle liste qui définis l'ordre de tirage des N états"""
    L = np.array(range(N)) 
    random.shuffle(L)
    return(L)
    
def init_states(N,Ex,Ey,Ez):
    config_dict={}
    dict_init={} 
    
    L = []
    c = 0
    n = N//2 + 1
    for k in range(-n,n):
        for j in range(-n,n):
            for i in range(-n,n):
                c += 1
                e  = Ex*k**2 + Ey*j**2 + Ez*i**2
                L+=[[e,k,j,i]]
                
    L1 = sorted(L, key=lambda x: x[0])
  
    if N%2 == 0:
        c = 0
        l = 0
        for k in range(N//2): 
            config_dict[f'{c}'] = np.array([L1[l][1],L1[l][2],L1[l][3],-1])
            c+=1
            config_dict[f'{c}'] = np.array([L1[l][1],L1[l][2],L1[l][3], 1])
            c+=1 
            l+=1
        return(config_dict)
            
    if N%2 == 1:
        c = 0
        l = 0
        for k in range(N//2): 
            config_dict[f'{c}'] = np.array([L1[l][1],L1[l][2],L1[l][3],-1])
            c+=1
            config_dict[f'{c}'] = np.array([L1[l][1],L1[l][2],L1[l][3], 1])
            c+=1 
            l+=1
        config_dict[f'{c}'] = np.array([L1[l][1],L1[l][2],L1[l][3],-1])
        return(config_dict)

def fermi_energy(N, Lx,Ly,Lz):
    return (3*np.pi**2*N/(Lx*Ly*Lz))**(2/3)*hbar**2/(2*me)

def Ncut(T,Ef,Lx,Ly,Lz):
    ET = 1/2*kb*T
    nx = int(mp.sqrt((ET + Ef*kb*T)/(hbar*2*np.pi)**2*(2*me*Lx**2)))
    ny = int(mp.sqrt((ET + Ef*kb*T)/(hbar*2*np.pi)**2*(2*me*Ly**2)))
    nz = int(mp.sqrt((ET + Ef*kb*T)/(hbar*2*np.pi)**2*(2*me*Lz**2)))
    
    return nx, ny, nz
    
#print(init_states2(10))

def choose_new_state(dict_config, position,n_cut):

    test = 0 
    # elec = dict_config[f'{position}'] 
    while test == 0:
        n_x =  random.randrange(-n_cut, n_cut+1)
        n_y =  random.randrange(-n_cut, n_cut+1)
        n_z =  random.randrange(-n_cut, n_cut+1)
        s   =  2*random.randrange(0, 2) - 1
        test = 1
        for cle, valeur in dict_config.items():
            if valeur[0] == n_x and valeur[1] == n_y and valeur[2] == n_z and valeur[3] == s:
                test = 0
    return(position, np.array([n_x, n_y, n_z, s]))

def proba(old_state, new_state, Ex,Ey,Ez, config_dict): #new_state is a list of the incoming numbers
    old_numbers=config_dict[f'{old_state}']
    if Ex*(old_numbers[0]**2-new_state[0]**2)+Ey*(old_numbers[1]**2-new_state[1]**2)+Ez*(old_numbers[2]**2-new_state[2]**2)>0:
        proba=1
    else:
        proba=np.exp(Ex*(old_numbers[0]**2-new_state[0]**2)+Ey*(old_numbers[1]**2-new_state[1]**2)+Ez*(old_numbers[2]**2-new_state[2]**2))
    x=random.random()
    #print(proba, 'proba')
    #print(old_numbers[0]**2-new_state[0]**2)
    if x<proba:
        config_dict[f'{old_state}']=new_state

def get_energy(particle, config_dict,Ex,Ey,Ez):
    return np.round(Ex*config_dict[f'{particle}'][0]**2+Ey*config_dict[f'{particle}'][1]**2+Ez*config_dict[f'{particle}'][2]**2,5)
