import itertools
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
    TIME_LIMIT = "2:00:00"
    MEMORY = "1G"
    OUTPUT_DIR = "SIMULATIONS"

    def __init__(self, lam: float, alpha: List[float], base_filename: str, meta_only: bool = False):
        self.lam = lam
        self.alpha = alpha
        self.base_filename = base_filename
        self.meta_only = meta_only
        self.filename = self._construct_filename()
        self.log_filename = self._construct_log_filename()
        self.command = self._construct_command()
        self.sbatch_command = self._construct_sbatch_command()

    def _construct_filename(self) -> str:
        return self.base_filename

    def _construct_log_filename(self) -> str:
        """Construct the log file name."""
        return os.path.join(self.OUTPUT_DIR, f"{self.filename}_log.txt")


    def _construct_command(self) -> str:
        """Construct the command to run the simulation."""
        meta_only = "--meta_only" if self.meta_only else ""
        return (f"python {self.SIMULATION_SCRIPT} {meta_only} --filename {self.filename} "
                f"--lam {self.lam} --alpha {' '.join(map(str, self.alpha))} --time_limit {self.TIME_LIMIT}")

    def _construct_sbatch_command(self) -> List[str]:
        """Construct the sbatch command to submit the job."""
        return ["sbatch", "-p", "regular", "-c", "1", "--mem", self.MEMORY, "--wrap", 
                self.command, "-o", self.log_filename]

    def submit(self) -> str:
        """Submit the job."""
        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        logger.info(f"Submitting job with command: {self.command}")
        logger.info(f"Alpha values for this job: {', '.join(map(str, self.alpha))}")
        result = subprocess.run(self.sbatch_command, capture_output=True, text=True)
        time.sleep(1)  # Delay to prevent overloading the system
        return result.stdout.split()[-1]  # return job id

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_filename", type=str, required=True,
                        help="Base filename for the output files.")
    parser.add_argument("--meta_only", action="store_true",
                        help="If specified, skips the simulations and only performs the meta analysis.")
    return parser.parse_args()

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_filename", type=str, required=True,
                        help="Base filename for the output files.")
    parser.add_argument("--meta_only", action="store_true",
                        help="If specified, skips the simulations and only performs the meta analysis.")
    parser.add_argument("--cutoffs_best_estimate", nargs=2, type=float, required=True,
                        help="Pair of cutoff values for best cohort estimate")
    parser.add_argument("--cutoffs_arbitrary_estimate", nargs=2, type=float, required=True,
                        help="Pair of cutoff values for arbitrary cohort estimate")
    return parser.parse_args()


def submit_meta_analysis_job(base_filename: str, dependency: str):
    """Submits a Meta-Analysis job to the SLURM scheduler."""
    command = ["python", "meta_analysis/meta_analysis.py", "--base_simulation_filename", base_filename]
    if dependency:
        sbatch_command = ["sbatch", "--dependency=afterok:" + dependency, "-p", "regular", "-c", "1", "--mem", "20G", "--wrap",
                ' '.join(command), "-o", f"SIMULATIONS/{base_filename}_meta_log.txt"]
    else:
        sbatch_command = ["sbatch", "-p", "regular", "-c", "1", "--mem", "20G", "--wrap",
                ' '.join(command), "-o", f"SIMULATIONS/{base_filename}_meta_log.txt"]
    subprocess.run(sbatch_command)

def submit_meta_analysis_job(base_filename: str, cutoffs_best_estimate: tuple, cutoffs_arbitrary_estimate: tuple, dependency: str):
    """Submits a Meta-Analysis job to the SLURM scheduler."""
    command = ["python", "meta_analysis/meta_analysis.py", "--base_simulation_filename", base_filename,
               "--cutoffs_best_estimate", str(cutoffs_best_estimate[0]), str(cutoffs_best_estimate[1]),
               "--cutoffs_arbitrary_estimate", str(cutoffs_arbitrary_estimate[0]), str(cutoffs_arbitrary_estimate[1])]
    if dependency:
        sbatch_command = ["sbatch", "--dependency=afterok:" + dependency, "-p", "regular", "-c", "1", "--mem", "20G", "--wrap",
                ' '.join(command), "-o", f"SIMULATIONS/{base_filename}_meta_log.txt"]
    else:
        sbatch_command = ["sbatch", "-p", "regular", "-c", "1", "--mem", "20G", "--wrap",
                ' '.join(command), "-o", f"SIMULATIONS/{base_filename}_meta_log.txt"]
    subprocess.run(sbatch_command)


def main():
    os.makedirs("SIMULATIONS", exist_ok=True)
    args = parse_args()
    
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

    # Base alphas
    base_alphas = [
        [1, 1, 1],
        [1, 2, 5],
        [1, 5, 2],
        [2, 1, 5],
        [2, 5, 1],
        [5, 1, 2],
        [5, 2, 1]
    ]

    # Multipliers
    multipliers = [1, 10, 100]

    # Compute the actual alphas by multiplying base alphas with multipliers
    alphas = [[x*y for x in alpha] for alpha in base_alphas for y in multipliers]


    lams: List[float] = [1, 2, 4]


    # Iterate over all alpha values
    job_ids = [] # List to collect job ids
    for lam in lams:
        for alpha in alphas:
            base_filename = f"{args.base_filename}_lam{lam}_alpha{'_'.join(map(str, alpha))}"
            if not args.meta_only:
                job = SimulationJob(lam, alpha, base_filename, args.meta_only)
                job_id = job.submit()
                job_ids.append(job_id)
                # Call meta analysis for this combination of lambda and alpha
                submit_meta_analysis_job(base_filename, args.cutoffs_best_estimate, args.cutoffs_arbitrary_estimate, job_id)
            else:
                # No dependencies when running meta-only
                submit_meta_analysis_job(base_filename, args.cutoffs_best_estimate, args.cutoffs_arbitrary_estimate, "")

            print(base_filename)
            time.sleep(1)

    # Call meta analysis for all simulations
    # All job ids are joined with ':' to ensure meta analysis runs after all jobs are done
    if job_ids:  # If jobs have been submitted
        submit_meta_analysis_job(args.base_filename, args.cutoffs_best_estimate, args.cutoffs_arbitrary_estimate, ':'.join(job_ids))
    else:  # If no jobs have been submitted (i.e., --meta_only was used)
        submit_meta_analysis_job(args.base_filename, args.cutoffs_best_estimate, args.cutoffs_arbitrary_estimate, "")


    print(args.base_filename)


if __name__ == "__main__":
    main()

