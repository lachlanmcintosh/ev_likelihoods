import os
import subprocess
import glob
import argparse
import time
import logging
from typing import List
from textwrap import dedent

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# User name for submitting and managing jobs.
USER_NAME = "lmcintosh"
# Directory where all simulations are saved.
SIMULATIONS_DIR = "/vast/scratch/users/lmcintosh/CN_SV_SNV_evolutionary_likelihood/SIMULATIONS"
WAITING_TIME = 20  # seconds

SBATCH_SCRIPT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={filename}_job_{job_id}
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --requeue
#SBATCH --qos=bonus
#SBATCH --output={sim_dir}/{filename}-%j_{job_id}.log
#SBATCH --error={sim_dir}/{filename}-%j_{job_id}.log
#SBATCH --time={time_limit}

source ~/anaconda3/etc/profile.d/conda.sh
conda activate sage

python do_all.py -t {job_id} -f {filename} --run_scripts 1 2 3 4 --lam {lam} --alpha {alpha}
"""


class SimulationManager:
    """
    Class to manage simulation jobs on a high performance computing cluster.
    """
    def __init__(self, filename: str, lam: float, alpha: List[float], max_cpus: int, time_limit: str):
        """
        Constructor for the SimulationManager object.
        """
        self.filename = filename
        self.lam = lam
        self.alpha = alpha
        self.max_cpus = max_cpus
        self.time_limit = time_limit


    @property
    def script_path(self) -> str:
        return os.path.join(SIMULATIONS_DIR, f"{self.filename}_singlecore")


    def _create_sbatch_script(self, job_id: int) -> str:
        script_path = f"{self.script_path}_{job_id}.sh"
        with open(script_path, "w") as script_file:
            script_file.write(SBATCH_SCRIPT_TEMPLATE.format(
                filename=self.filename, 
                job_id=job_id, 
                sim_dir=SIMULATIONS_DIR, 
                time_limit=self.time_limit, 
                lam=self.lam, 
                alpha=' '.join(map(str, self.alpha))))
        return script_path


    def _count_running_jobs(self) -> int:
        """
        Counts the number of currently running jobs associated with the filename.
        """
        completed_process = subprocess.run(["squeue", "-u", USER_NAME], capture_output=True, text=True)
        # Subtract 1 as the first line is a header line.
        return len(completed_process.stdout.splitlines()) - 1

    def _get_jobs(self) -> List[str]:
        """
        Get the jobs related to the filename.
        """
        completed_process = subprocess.run(["squeue", "-u", USER_NAME, "-o", "%.100j"], capture_output=True, text=True)
        return [line.strip() for line in completed_process.stdout.splitlines()[1:] if line.strip().startswith(f"{self.filename}_job")]

    def submit_jobs(self, total_reps: int) -> None:
        """
        Submits the specified number of jobs to the cluster.
        """
        job_id = 0
        while job_id < total_reps:
            sbatch_script = self._create_sbatch_script(job_id)
            sbatch_command = ["sbatch", sbatch_script]
            subprocess.run(sbatch_command, check=True)
            job_id += 1

    def remove_logs(self) -> None:
        """
        Removes all logs, pickles, and sbatch scripts associated with the filename.
        """
        file_patterns = [f"{self.filename}-*.log", f"{self.filename}_*.pickle", f"{self.filename}_singlecore_*.sh"]
        for file_pattern in file_patterns:
            for file in glob.glob(os.path.join(SIMULATIONS_DIR, file_pattern)):
                os.remove(file)

    def cancel_existing_jobs(self):
        """
        Cancel existing jobs related to the filename.
        """
        jobs = self._get_jobs()
        for job_name in jobs:
            job_id = job_name.split("_")[-1]
            logging.warning(f"Found existing job {job_id} with name {job_name}. Cancelling this job.")
            subprocess.run(["scancel", job_id], check=True)
            # Assert that the job has indeed been cancelled.
            assert job_id not in self._get_jobs(), f"Failed to cancel job {job_id}."

    def are_jobs_finished(self) -> bool:
        """
        Check if all jobs related to the filename have finished.
        """
        return len(self._get_jobs()) == 0

def parse_args() -> argparse.Namespace:
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, default="simulation",
                        help="The base name for the log files. Default is 'simulation'.")
    parser.add_argument("--lam", type=float, default=2.5,
                        help="Lambda value for the simulation. Default is 1.0.")
    parser.add_argument("--alpha", nargs="+", type=float, default=[20.0, 20.0, 100.0],
                        help="Alpha values for the simulation as a list of floats.")
    parser.add_argument("--total_reps", type=int, default=50,
                        help="Total number of repetitions for the simulation. Default is 40.")
    parser.add_argument("--max_cpus", type=int, default=1000,
                        help="The maximum number of CPUs to use. Default is 300.")
    parser.add_argument("--time_limit", type=str, default='20:00',
                        help="The time limit for each job in the format 'MM:SS'. Default is '01:00' for 1 minute.")
    return parser.parse_args()

def main() -> None:
    """
    The main function that is run when the script is executed.
    """
    args = parse_args()

    manager = SimulationManager(
        filename=args.filename,
        lam=args.lam,
        alpha=args.alpha,
        max_cpus=args.max_cpus,
        time_limit=args.time_limit
    )

    # Cancel any existing jobs with the same filename
    manager.cancel_existing_jobs()
    
    manager.remove_logs()
    manager.submit_jobs(args.total_reps)

    while not manager.are_jobs_finished():
        logging.info("Jobs are still running, waiting...")
        time.sleep(WAITING_TIME)


if __name__ == "__main__":
    main()

