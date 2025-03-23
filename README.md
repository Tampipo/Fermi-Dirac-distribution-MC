# Fermi-Dirac Distribution Simulation

This project simulates the Fermi-Dirac distribution using Monte Carlo methods. It includes scripts to scan the chemical potential, temperature, and plot phase space.

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/Fermi-Dirac-distribution-MC.git
    cd Fermi-Dirac-distribution-MC
    ```

2. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Arguments

The main script accepts several arguments to configure the simulation. You can specify the box dimensions, whether to run in parallel, and the specific job to run.

### Jobs

There are three main jobs you can run:

1. **part_scan**: Scan the chemical potential.
    ```bash
    python main.py part_scan --T 2000 --n_toys 100000 --N 10 20 50 70 100 150 200
    ```

2. **temp_scan**: Scan the temperature.
    ```bash
    python main.py temp_scan --N 100 --n_toys 100000 --T 1000 2000 3000 4000 5000 10000 20000 30000 40000 50000 60000
    ```

3. **phase_space**: Plot phase space.
    ```bash
    python main.py phase_space --T 10 --n_toys 100000 --N 100
    ```

### Example

To run a simulation scanning the chemical potential:
```bash
python main.py part_scan --T 2000 --n_toys 100000 --N 10 20 50 70 100 150 200
```

## Output

The simulation results are saved in the `Plots` directory. The plots include the Fermi-Dirac distribution and the chemical potential.

## License

This project is licensed under the MIT License.
