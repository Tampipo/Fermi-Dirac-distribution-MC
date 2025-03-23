import argparse


def SetupArgs():
    parser = argparse.ArgumentParser(description='Fermi-Dirac distribution simulation')
    parser.add_argument('--Lx', type=float, default=10**(-8), help='Box dimension in x direction (m)')
    parser.add_argument('--Ly', type=float, default=10**(-8), help='Box dimension in y direction (m)')
    parser.add_argument('--Lz', type=float, default=10**(-8), help='Box dimension in z direction (m)')
    parser.add_argument('--parallel', action='store_true', default=False, help='Run the simulation in parallel')



    subparsers = parser.add_subparsers(dest="job", help="Choose a job to run", required=True)

    part_parser = subparsers.add_parser('part_scan', help='Scan the chemical potential')
    part_parser.add_argument('--T', type=float, default=2000, help='Temperature (K)')
    part_parser.add_argument('--n_toys', type=int, default=100000, help='Number of steps')
    part_parser.add_argument('--N', type=int, nargs='+', default=[10,20,50,70,100,150,200], help='Number of particles to scan')

    temp_parser = subparsers.add_parser('temp_scan', help='Scan the temperature')
    temp_parser.add_argument('--N', type=int, default=100, help='Number of particles')
    temp_parser.add_argument('--n_toys', type=int, default=100000, help='Number of steps')
    temp_parser.add_argument('--T', type=float, nargs='+', default=[1000,2000,3000,4000, 5000,10000,20000,30000,40000, 50000, 60000], help='Temperature (K) to scan')

    phase_parser = subparsers.add_parser('phase_space', help='Plot phase space')
    phase_parser.add_argument('--T', type=float, default=10, help='Temperature (K)')
    phase_parser.add_argument('--n_toys', type=int, default=100000, help='Number of steps')
    phase_parser.add_argument('--N', type=int, default=100, help='Number of particles')

    args = parser.parse_args()
    print(args)

    return args


def OrganizeArgs(args):
    job_args = {
        "global": {
            "N": args.N,
            "Lx": args.Lx,
            "Ly": args.Ly,
            "Lz": args.Lz,
            "parallel": args.parallel
        },
        "job": args.job,
        "params": vars(args)
    }
    return job_args


