import random
import numpy as np
import math as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
### Physics Constants ###

import numpy as np
import math as mp 
import random
import matplotlib.pyplot as plt

### Physics Constants ###

hbar = 6.02*10**(-34) # J.s
kb   = 1.38*10**(-23) # J.K^-1
me   = 9.11*10**(-31) # kg



def fermi_distrib(E,mu,T):
    return 1/(1+np.exp(E-mu))

def chemical_potential(T,Ef):
    mu=Ef*(1-np.pi**2/12*(1/Ef)**2)
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
        N = 100
        T = 0.1
        Lx = 10**(-6)
        Ly = 10**(-6)
        Lz = 10**(-6)
        n_step = 10000
 
    print("Initializing parameters...")
    print(Lx)
    init_param=init(T,Lx,Ly,Lz)
    Ex=init_param[0]
    Ey=init_param[1]
    Ez=init_param[2]
    config_dict = init_states(N,Ex,Ey,Ez)
    Ef=get_energy(N-1, config_dict, Ex,Ey, Ez)
    mu=chemical_potential(T,Ef)
    energies_plot=np.linspace(0,100*Ex,1000)
    fermi_dirac=[fermi_distrib(E, Ef, T) for E in energies_plot]
    #print(mu, T)
    n_cut=min(Ncut(T,Ef,Lx,Ly,Lz))
    #n_cut=30
    print("**********Simulation parameters*********")
    print("Number of particles: ", N)
    print("Temperature: ", T, "(K)")
    print("Box dimensions: ", Lx, Ly, Lz, "(m)" )
    print("States cut : ", n_cut)
    print("Energie along x",Ex)
    print("Fermi energy",Ef)

    
    #print(config_dict)

    e = 0
    for cle, valeur in config_dict.items():
        e += mp.sqrt(valeur[0]**2+valeur[1]**2+valeur[2]**2)/N
    x = [0]
    energie_moy = [e]
    
    Fermi_Dirac=[[[0,0,0,-1],0]]
    for k in range(n_step):
        liste = create_liste(N)
        for l in liste:
            old_state = l
            new_state = choose_new_state(config_dict, l, n_cut)[1]
            proba(old_state, new_state, Ex, Ey, Ez, config_dict)
        E=[]
        e = 0
        for i in range (0,N):
            E.append(get_energy(i, config_dict,Ex,Ey,Ez))
        x += [k]
        E=np.array(E)
        e= np.mean(E)
        x += [k+1]
        energie_moy += [e]
        part=[]
        for cle,valeur in config_dict.items():
            part.append(list(valeur))
        #print(len(part))

        for i in range(0,len(Fermi_Dirac)):
            n=0
            states_exists=False
            index=0
            for j in range (0,len(part)):
                if Fermi_Dirac[i][0]==part[j]:
                    index=j
                    states_exists=True
            if states_exists:
                Fermi_Dirac[i][1]=k/(k+1)*Fermi_Dirac[i][1]+1/(k+1)
                part.pop(index)
        # print(len(Fermi_Dirac))
        
        #print(part)
        treated_states=[]

        for j in range (0,len(part)):
            if part[j] not in treated_states:
                n=0
                for s in range (0,len(part)):
                    if part[s]==part[j]:
                        n+=1
                Fermi_Dirac.append([part[j],n/(k+1)])
                treated_states.append(part[j])
        #plt.plot(x, energie_moy)
        #plt.show()
        if k%(n_step/10)==0:
            print(f'Step {k}')
            fig1 = plt.figure()
            ax1 = fig1.add_subplot()     
            Fermi_Dirac_energies=[]
            Fermi_Dirac_part=[]
            Fermi_Dirac_part_mean=[]
            Fermi_Dirac_part_std=[]
            for i in range (0,len(Fermi_Dirac)):
                Fermi_Dirac_energies.append(np.round(Ex*Fermi_Dirac[i][0][0]**2+Ey*Fermi_Dirac[i][0][1]**2+Ez*Fermi_Dirac[i][0][2]**2,5))
                Fermi_Dirac_part.append(Fermi_Dirac[i][1])
            distinct_energies=set(Fermi_Dirac_energies)
            #print(distinct_energies)
            for energy in distinct_energies:
                correct_energies_part=[]
                for i in range (0,len(Fermi_Dirac_energies) ):
                    if Fermi_Dirac_energies[i]==energy:
                        correct_energies_part.append(Fermi_Dirac_part[i])
                array=np.array(correct_energies_part)
                Fermi_Dirac_part_mean.append(np.mean(array))
                Fermi_Dirac_part_std.append(np.std(array))
            print(len(Fermi_Dirac_part_mean), len(list(distinct_energies)))
            ax1.scatter(list(distinct_energies),np.array(Fermi_Dirac_part_mean), color='blue', label='mean')
            ax1.plot(energies_plot, fermi_dirac, label="Fermi-Dirac", color="red")
            ax1.set_xlim(0,3*Ef)
            ax1.axvline(x = Ef, label = 'Fermi Energy',linestyle='--')
            ax1.set_ylim(0,max(Fermi_Dirac_part_mean)+0.1*max(Fermi_Dirac_part_mean))
            ax1.set_ylabel("Occupation number")
            ax1.set_xlabel("Energy (dimensionless)")
            ax1.set_title("Fermi_Dirac distribution")
            ax1.legend()

            plt.savefig(f'./img/{k}_fermi_dirac')


            ax=plt.gca()
            ax.cla()
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')     
            ax.set_title("Momentum space")
            ax.set_xlabel("Nx")
            ax.set_ylabel("Ny")
            ax.set_zlabel("Nz")
            ax.set_xlim(-3,3)
            ax.set_ylim(-3,3)
            ax.set_zlim(-3,3)   
            #ax.tick_params(axis='x', labelrotation = 20)
            xu,yu,zu=[],[],[]
            xd,yd,zd=[],[],[]
            pts_u=[]
            pts_d=[]
            #print(config_dict)
            for i in range (0, len(config_dict)):
                valeur=config_dict[f'{i}']
                if valeur[3]==1:
                    xu.append(valeur[0])
                    yu.append(valeur[1])
                    zu.append(valeur[2])
                    pts_u.append([valeur[0],valeur[1],valeur[2]])
                if valeur[3]==-1:
                    xd.append(valeur[0])
                    yd.append(valeur[1])
                    zd.append(valeur[2])
                    pts_d.append([valeur[0],valeur[1],valeur[2]])
            array_up=np.array(pts_u)
            array_down=np.array(pts_d)
            hull=ConvexHull(array_up)
            ax.plot(array_up.T[0], array_up.T[1], array_up.T[2], "o", color="blue", label="Particles")
            for s in hull.simplices:
                s = np.append(s, s[0])  # Here we cycle back to the first coordinate
                ax.plot(array_up[s, 0], array_up[s, 1], array_up[s, 2], "r-")
            #ax.scatter(xu,yu,zu, label='Up', marker='o',color='b')   
            #ax.scatter(xd,yd,zd, label='Down', marker='^',color='r') 
            plt.legend()        
            plt.savefig(f'./img/{k}_momentum_space')

    


if __name__=="__main__":
    main()
