from utils import *
from labellines import labelLine, labelLines
from scipy.optimize import curve_fit

def fit_fermi_dirac(x,a):
    return 1/(1+np.exp(x-a))

def fermi_dirac_temp(E,mu,T):
    return 1/(1+np.exp((E-mu)/(kb*T)))
temp=[500,1000,2500,5000,7500,10000]
colors=["r", "b", "g", "m","c", "y"]


def main():
    c=0
    N = 50
    Lx = 10**(-8)
    Ly = 10**(-8)
    Lz = 10**(-8)
    n_step = 100000
    fig1 = plt.figure()
    ax1 = fig1.add_subplot()  
    chemical_potentials=[]
    for T in temp:
        print("Initializing parameters...")
        init_param=init(T,Lx,Ly,Lz)
        Ex=init_param[0]
        Ey=init_param[1]
        Ez=init_param[2]
        config_dict = init_states(N,Ex,Ey,Ez)
        #Ef=get_energy(N-1, config_dict, Ex,Ey, Ez)
        Ef=fermi_energy(N,Lx,Ly,Lz)/(kb*T)
        energies_plot=np.linspace(0,100*Ex,1000)
        fermi_dirac=[fermi_distrib(E, Ef) for E in energies_plot]
        n_cut=5*min(Ncut(T,Ef,Lx,Ly,Lz))
        print("**********Simulation parameters*********")
        print("Number of particles: ", N)
        print("Temperature: ", T, "(K)")
        print("Box dimensions: ", Lx, Ly, Lz, "(m)" )
        print("States cut : ", n_cut)
        print("Energie along x",Ex)
        print("Fermi energy",Ef)
        e = 0
        for cle, valeur in config_dict.items():
            e += mp.sqrt(valeur[0]**2+valeur[1]**2+valeur[2]**2)/N
        x = [0]
        energie_moy = [e]
        
        Fermi_Dirac=[[[0,0,0,-1],0]]
        for k in range(n_step):
            if k%(n_step/10)==0:
                print(f'Step {k}, Temprature {T} K')
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

            treated_states=[]

            for j in range (0,len(part)):
                if part[j] not in treated_states:
                    n=0
                    for s in range (0,len(part)):
                        if part[s]==part[j]:
                            n+=1
                    Fermi_Dirac.append([part[j],n/(k+1)])
                    treated_states.append(part[j])

        Fermi_Dirac_energies=[]
        Fermi_Dirac_part=[]
        Fermi_Dirac_part_mean=[]
        Fermi_Dirac_part_std=[]
        for i in range (0,len(Fermi_Dirac)):
            Fermi_Dirac_energies.append(np.round(Ex*Fermi_Dirac[i][0][0]**2+Ey*Fermi_Dirac[i][0][1]**2+Ez*Fermi_Dirac[i][0][2]**2,5))
            Fermi_Dirac_part.append(Fermi_Dirac[i][1])
        distinct_energies=set(Fermi_Dirac_energies)
        
        for energy in distinct_energies:
            correct_energies_part=[]
            for i in range (0,len(Fermi_Dirac_energies) ):
                if Fermi_Dirac_energies[i]==energy:
                    correct_energies_part.append(Fermi_Dirac_part[i])
            array=np.array(correct_energies_part)
            Fermi_Dirac_part_mean.append(np.mean(array))
            Fermi_Dirac_part_std.append(np.std(array))
        print(len(Fermi_Dirac_part_mean), len(list(distinct_energies)))
        popt, pcov= curve_fit(fit_fermi_dirac, np.array(list(distinct_energies)), np.array(Fermi_Dirac_part_mean))
        xdata=np.linspace(0,5*Ef,1000)
        data=kb*T*xdata
        ax1.plot(data, fermi_dirac_temp((data), popt[0]*kb*T, T),label=f'{T} K',color=colors[c])
        chemical_potentials.append(popt[0]*kb*T)
        ax1.scatter(np.array(list(distinct_energies))*kb*T,np.array(Fermi_Dirac_part_mean), color=colors[c])
        c+=1

    labelLines(ax1.get_lines(),zorder=2.5)
    #ax1.plot(energies_plot, fermi_dirac, label="Fermi-Dirac", linestyle="--", color="gray")
    ax1.set_xlim(0,1.6*10**(-18))
    #ax1.axvline(x = Ef, label = 'Fermi Energy',linestyle='--')
    #ax1.set_ylim(0,max(Fermi_Dirac_part_mean)+0.1*max(Fermi_Dirac_part_mean))
    ax1.set_ylabel("Occupation number")
    ax1.set_xlabel("Energy")
    ax1.set_title("Fermi_Dirac distribution")
    plt.savefig(f'./img_simu/chem_potential')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    ax2.scatter(temp,chemical_potentials)
    ax2.set_ylabel("Chemical potential")
    ax2.set_xlabel("Temperature (K)")
    ax2.set_title("Chemical potential")
    plt.savefig(f'./img_simu/chem_potential_all_temps')


if __name__=="__main__":
    main()
