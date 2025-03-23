from utils import *
from labellines import labelLine, labelLines


colors=["r", "b", "g", "m","c", "y", "C1", "C6", "C7", "C4", "C2"]
plt.rcParams['text.usetex'] = True


def simulate(N, Lx, Ly, Lz, T, n_step):
    simul = Simulate(N, Lx, Ly, Lz, T, n_step)
    simul.run()
    popt, Ef, distinct, Part_mean =  simul.distrib()
    return popt[0], Ef, distinct, Part_mean



def main(job_args):
    colors=["r", "b", "g", "m","c", "y", "C1", "C6", "C7", "C4", "C2"]
    N = job_args["global"]["N"]
    Lx = job_args["global"]["Lx"]
    Ly = job_args["global"]["Ly"]
    Lz = job_args["global"]["Lz"]
    parallel = job_args["global"]["parallel"]
    T = job_args["params"]["T"]
    n_step = job_args["params"]["n_toys"]
    from tqdm import tqdm
    if parallel:
        from joblib import Parallel, delayed
        from tqdm_joblib import tqdm_joblib
        import multiprocessing
        with tqdm_joblib(tqdm(desc="Simulating", total=len(T))) as progress_bar:
            res = Parallel(n_jobs=4)(delayed(simulate)(N, Lx, Ly, Lz, T[i], n_step) for i in range(len(T)))
            mus = [res[i][0] for i in range(len(T))]
            Ef = [res[i][1]/(kb*T[i]) for i in range(len(T))]
            distinct_energies = [res[i][2] for i in range(len(T))]
            Fermi_Dirac_part_means = [res[i][3] for i in range(len(T))]
    else:
        mus = []
        Ef = []
        distinct_energies = []
        Fermi_Dirac_part_means = []
        for i in tqdm(range(len(T))):
            mu, ef, dist, part_mean= simulate(N, Lx, Ly, Lz, T[i], n_step)
            mus.append(mu)
            Ef.append(ef/(kb*T[i]))
            distinct_energies.append(dist)
            Fermi_Dirac_part_means.append(part_mean)
    
    ## Plotting
    mus = np.array(mus)
    Ef = np.array(Ef)
    fig, ax = plt.subplots()
    for i in range(len(T)):
        xdata=np.linspace(0,5*Ef[i],1000)
        data=kb*T[i]*xdata
        ax.plot(data, fermi_dirac_temp((data), mus[i]*kb*T[i], T[i]),label=f'$T={T[i]}$ K',color=colors[i])
        ax.scatter(np.array(list(distinct_energies[i]))*kb*T[i],np.array(Fermi_Dirac_part_means[i]), color=colors[i])

    ax.set_xlim(0,0.8*10**(-18))
    ax.axvline(x = Ef[-1]*kb*T[-1], label = 'Fermi Energy',linestyle='--')
    ax.set_ylabel(r"$n(E)$")
    ax.set_xlabel(r"$E$ (J)")
    ax.set_title("Fermi Dirac distribution")
    plt.legend()
    plt.savefig(f'./Plots/chem_potential_temp')

    fig2, ax2 = plt.subplots()
    ax2.scatter(1/Ef,mus/Ef, label="Simulation")
    ax2.set_ylabel(r"$\mu/E_{f}$")
    low_temp=np.linspace(0,1, 2000)
    low_chem=[chemical_potential_low(T) for T in low_temp]
    ax2.plot(low_temp, low_chem, label=r"Low $T$, $\mu$", color='c', linestyle='--')
    high_temp=np.linspace(0.5,2, 2000)
    high_chem=[chemical_potential_high(T) for T in high_temp]
    ax2.plot(high_temp, high_chem, label=r"High $T$, $\mu$", color='r', linestyle='--')
    ax2.set_xlabel(r"$T/T_{f}$")
    ax2.set_title("Chemical potential")
    plt.legend()
    plt.savefig(f'./Plots/chem_potential_all_temps')
