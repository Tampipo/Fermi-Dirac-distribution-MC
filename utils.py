import random
import numpy as np
import math as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
plt.rcParams['text.usetex'] = True
### Physics Constants ###

hbar = 6.02*10**(-34) # J.s
kb   = 1.38*10**(-23) # J.K^-1
me   = 9.11*10**(-31) # kg

### Useful functions ###

def fermi_distrib(E,mu):
    return 1/(1+np.exp(E-mu))

def fermi_dirac_temp(E,mu,T):
    return 1/(1+np.exp((E-mu)/(kb*T)))

def chemical_potential_low(Tnorm):
    mu=(1-np.pi**2/12*(Tnorm)**2)
    return mu

def chemical_potential_high(Tnorm):
    mu=3/2*Tnorm*np.log(4*np.pi/(6*np.pi**2)**(2/3)*(1/Tnorm))
    return mu

def chemical_potential_low_part(N,Lx,Ly,Lz,T):
    c=hbar**2/(2*me)*(3*np.pi**2/(Lx*Ly*Lz))**(2/3)
    #print(c)
    mu=c*N**(2/3)*(1-np.pi**2/12*T**2*kb**2/c**2*N**(-4/3))
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


def UpdateFermiDiracStates(Fermi_Dirac, part, k):
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
        else:
            Fermi_Dirac[i][1]=k/(k+1)*Fermi_Dirac[i][1]

    treated_states=[]

    for j in range (0,len(part)):
        if part[j] not in treated_states:
            n=0
            for s in range (0,len(part)):
                if part[s]==part[j]:
                    n+=1
            Fermi_Dirac.append([part[j],n/(k+1)])
            treated_states.append(part[j])


class Simulate:
    def __init__(self, N, Lx, Ly, Lz, T, n_step):
        self.N = N
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.T = T
        self.n_step = n_step
        self.init_param = init(T, Lx, Ly, Lz)
        self.Ex = self.init_param[0]
        self.Ey = self.init_param[1]
        self.Ez = self.init_param[2]
        self.config_dict = init_states(N, self.Ex, self.Ey, self.Ez)
        self.Ef = fermi_energy(N, Lx, Ly, Lz)
        self.n_cut =  30
        self.current_step = 0
        self.states = [[[0,0,0,-1],0]]

    def update_states(self):
        part=[]
        for cle,valeur in self.config_dict.items():
            part.append(list(valeur))
        UpdateFermiDiracStates(self.states, part, self.current_step)

    def run(self,save_phase_space=False):
        from tqdm import tqdm
        for k in tqdm(range(self.n_step), desc = "Running MC steps..." ):
            liste = create_liste(self.N)
            for l in liste:
                old_state = l
                new_state = choose_new_state(self.config_dict, l, self.n_cut)[1]
                proba(old_state, new_state, self.Ex, self.Ey, self.Ez, self.config_dict)
            self.update_states()
            self.current_step += 1
            if k%(self.n_step/10)==0 and save_phase_space:
                _,_,dist_energies, dist_part = self.distrib()
                fig, ax = plt.subplots()
                ax.scatter(list(dist_energies),np.array(dist_part), color='blue', label='mean')
                energies_plot=np.linspace(0,100*self.Ex,1000)
                fermi_Ef = lambda x: fermi_distrib(x,self.Ef/(kb*self.T))
                ax.plot(energies_plot, fermi_Ef(energies_plot), label="Fermi-Dirac", color="red")
                ax.set_xlim(0,7*self.Ef/(kb*self.T))
                ax.axvline(x = self.Ef/(kb*self.T), label = 'Fermi Energy',linestyle='--')
                ax.set_ylim(0,max(dist_part)+0.1*max(dist_part))
                ax.set_ylabel("Occupation number")
                ax.set_xlabel("Energy (dimensionless)")
                ax.set_title("Fermi Dirac distribution")
                ax.legend()

                plt.savefig(f'./Plots/{k}_fermi_dirac')

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
                for i in range (0, len(self.config_dict)):
                    valeur=self.config_dict[f'{i}']
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
                plt.legend()        
                plt.savefig(f'./Plots/{k}_momentum_space')

    def distrib(self):
        Fermi_Dirac_energies=[]
        Fermi_Dirac_part=[]
        Fermi_Dirac_part_mean=[]
        Fermi_Dirac_part_std=[]
        for i in range (0,len(self.states)):
            Fermi_Dirac_energies.append(np.round(self.Ex*self.states[i][0][0]**2+self.Ey*self.states[i][0][1]**2+self.Ez*self.states[i][0][2]**2,5))
            Fermi_Dirac_part.append(self.states[i][1])
        distinct_energies=set(Fermi_Dirac_energies)
        
        for energy in distinct_energies:
            correct_energies_part=[]
            for i in range (0,len(Fermi_Dirac_energies) ):
                if Fermi_Dirac_energies[i]==energy:
                    correct_energies_part.append(Fermi_Dirac_part[i])
            array=np.array(correct_energies_part)
            Fermi_Dirac_part_mean.append(np.mean(array))
            Fermi_Dirac_part_std.append(np.std(array))
        popt, pcov= curve_fit(fermi_distrib, np.array(list(distinct_energies)), np.array(Fermi_Dirac_part_mean))
        return popt, self.Ef, distinct_energies, Fermi_Dirac_part_mean
