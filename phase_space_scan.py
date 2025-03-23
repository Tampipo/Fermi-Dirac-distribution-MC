from utils import *
from labellines import labelLine, labelLines

colors=["r", "b", "g", "m","c", "y", "C1", "C6", "C7", "C4", "C2"]
plt.rcParams['text.usetex'] = True


def simulate(N, Lx, Ly, Lz, T, n_step, save_phase_space=False):
    simul = Simulate(N, Lx, Ly, Lz, T, n_step)
    simul.run(save_phase_space)
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
    mus = []
    Ef = []
    distinct_energies = []
    Fermi_Dirac_part_means = []
    _ = simulate(N, Lx, Ly, Lz, T, n_step, save_phase_space=True)


    
