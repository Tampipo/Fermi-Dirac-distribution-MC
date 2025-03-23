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
        with tqdm_joblib(tqdm(desc="Simulating", total=len(N))) as progress_bar:
            res = Parallel(n_jobs=4)(delayed(simulate)(N[i], Lx, Ly, Lz, T, n_step) for i in range(len(N)))
            mus = [res[i][0] for i in range(len(N))]
            Ef = [res[i][1]/(kb*T) for i in range(len(N))]
            distinct_energies = [res[i][2] for i in range(len(N))]
            Fermi_Dirac_part_means = [res[i][3] for i in range(len(N))]
    else:
        mus = []
        Ef = []
        distinct_energies = []
        Fermi_Dirac_part_means = []
        for i in tqdm(range(len(N))):
            mu, ef, dist, part_mean= simulate(N[i], Lx, Ly, Lz, T, n_step)
            mus.append(mu)
            Ef.append(ef/(kb*T))
            distinct_energies.append(dist)
            Fermi_Dirac_part_means.append(part_mean)
    
    ## Plotting
    mus = np.array(mus)
    Ef = np.array(Ef)
    fig, ax = plt.subplots()
    for i in range(len(N)):
        xdata=np.linspace(0,5*Ef[i],1000)
        data=kb*T*xdata
        ax.plot(data, fermi_dirac_temp((data), mus[i]*kb*T, T),label=f'$N={N[i]}$',color=colors[i])
        ax.scatter(np.array(list(distinct_energies[i]))*kb*T,np.array(Fermi_Dirac_part_means[i]), color=colors[i])

    ax.set_xlim(0,1.3*10**(-18))
    ax.set_ylabel(r"$n(E)$")
    ax.set_xlabel(r"$E$ (J)")
    ax.set_title("Fermi Dirac distribution")
    plt.legend()
    plt.savefig(f'./Plots/chem_potential_part')
    plt.close(fig)

    fig2, ax2 = plt.subplots()
    parts=np.linspace(1,210, 10000)
    chem_theo=[chemical_potential_low_part(part,Lx,Ly,Lz,T) for part in parts]
    ax2.scatter(N,mus*kb*T, label="Simulation")
    ax2.plot(parts, chem_theo, label=r"Low $T$ theory", linestyle='--', color='r')
    ax2.set_ylabel(r"$\mu$")
    ax2.set_xlabel(r"$N$")
    ax2.set_title("Chemical potential")
    plt.legend()
    plt.savefig(f'./Plots/chem_potential_all_parts')
