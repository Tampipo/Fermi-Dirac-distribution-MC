from utils import *
plt.rcParams['text.usetex'] = True





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
        T = 10
        Lx = 10**(-8)
        Ly = 10**(-8)
        Lz = 10**(-8)
        n_step = 100000
 
    print("Initializing parameters...")
    print(Lx)
    init_param=init(T,Lx,Ly,Lz)
    Ex=init_param[0]
    Ey=init_param[1]
    Ez=init_param[2]
    config_dict = init_states(N,Ex,Ey,Ez)
    #Ef=get_energy(N-1, config_dict, Ex,Ey, Ez)
    Ef=fermi_energy(N,Lx,Ly,Lz)/(kb*T)
    energies_plot=np.linspace(0,100*Ex,1000)
    fermi_dirac=[fermi_distrib(E, Ef) for E in energies_plot]
    #print(mu, T)
    n_cut=5*min(Ncut(T,Ef,Lx,Ly,Lz))
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
            else:
                Fermi_Dirac[i][1]=k/(k+1)*Fermi_Dirac[i][1]
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
            ax1.set_xlim(0,7*Ef)
            ax1.axvline(x = Ef, label = 'Fermi Energy',linestyle='--')
            ax1.set_ylim(0,max(Fermi_Dirac_part_mean)+0.1*max(Fermi_Dirac_part_mean))
            ax1.set_ylabel("Occupation number")
            ax1.set_xlabel("Energy (dimensionless)")
            ax1.set_title("Fermi_Dirac distribution")
            ax1.legend()

            plt.savefig(f'./img_simu/{k}_fermi_dirac')


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
            plt.savefig(f'./img_simu/{k}_momentum_space')

    


if __name__=="__main__":
    main()
