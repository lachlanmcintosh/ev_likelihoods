import os
import subprocess
import argparse
import time
from typing import List, Tuple
import logging

# Setup logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class SimulationJob:
    """Class that represents a simulation job."""

    SIMULATION_SCRIPT = "run_standard_simulations.py"
    TIME_LIMIT = "30:00"
    MEMORY = "40G"
    OUTPUT_DIR = "SIMULATIONS"

    def __init__(self, lam: float, alpha: List[float], meta_only: bool = False):
        self.lam = lam
        self.alpha = alpha
        self.meta_only = meta_only
        self.filename = self._construct_filename()
        self.log_filename = self._construct_log_filename()
        self.command = self._construct_command()
        self.sbatch_command = self._construct_sbatch_command()

    def _construct_filename(self) -> str:
        lam_str = str(self.lam)
        alpha_str = '_'.join(map(str, self.alpha))
        return f"simulation_lam{lam_str}_alpha{alpha_str}"

    def _construct_log_filename(self) -> str:
        """Construct the log file name."""
        return os.path.join(self.OUTPUT_DIR, f"{self.filename}_log.txt")

    def _construct_command(self) -> str:
        """Construct the command to run the simulation."""
        meta_only = "--meta_only" if self.meta_only else ""
        return (f"python {self.SIMULATION_SCRIPT} {meta_only} --filename {self.filename} "
                f"--alpha {' '.join(map(str, self.alpha))} --time_limit {self.TIME_LIMIT}")

    def _construct_sbatch_command(self) -> List[str]:
        """Construct the sbatch command to submit the job."""
        return ["sbatch", "-p", "regular", "-c", "1", "--mem", self.MEMORY, "--wrap", 
                self.command, "-o", self.log_filename]

    def submit(self):
        """Submit the job."""
        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        logger.info(f"Submitting job with command: {self.command}")
        logger.info(f"Alpha values for this job: {', '.join(map(str, self.alpha))}")
        subprocess.run(self.sbatch_command)
        time.sleep(1)  # Delay to prevent overloading the system


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta_only", action="store_true",
                        help="If specified, skips the simulations and only performs the meta analysis.")
    return parser.parse_args()

def submit_meta_analysis_job(base_filename: str):
    """Submits a Meta-Analysis job to the SLURM scheduler."""
    command = ["python", "meta_analysis/meta_analysis.py", "--base_simulation_filename", base_filename]
    sbatch_command = ["sbatch", "-p", "regular", "-c", "1", "--mem", "20G", "--wrap", 
                ' '.join(command), "-o", f"{base_filename}_meta_log.txt"]
    subprocess.run(sbatch_command)
    time.sleep(1)

def main():
    alphas: List[List[float]] = [
        [0.5, 0.5, 0.5],
        [1, 1, 1],
        [2, 2, 2],
        [3, 3, 3],
        [4, 4, 4],
        [2, 5, 1],
        [5, 2, 1],
        [1, 2, 5],
        [1, 5, 2],
        [5, 1, 2]
    ]

    lams: List[float] = [0.5, 1, 2, 4]

    args = parse_args()

    # Iterate over all alpha values
    for lam in lams:
        for alpha in alphas:
            job = SimulationJob(lam, alpha, args.meta_only)
            job.submit()

            # Call meta analysis for this combination of lambda and alpha
            base_filename = f"simulation_lam{lam}_alpha{'_'.join(map(str, alpha))}"
            submit_meta_analysis_job(base_filename)

            #subprocess.run(["python", "meta_analysis/meta_analysis.py", "--base_simulation_filename", base_filename])
  
    # Call meta analysis for all simulations
    #subprocess.run(["python", "meta_analysis/meta_analysis.py", "--base_simulation_filename", "simulation"])
    submit_meta_analysis_job("simulation")

if __name__ == "__main__":
    main()

