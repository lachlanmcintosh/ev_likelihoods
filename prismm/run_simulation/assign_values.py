import numpy as np
from scipy.stats import poisson
from typing import Optional, List, Tuple
from dataclasses import dataclass

from clonal_trees.run_simulation.simulation_priors.random_number_generator import generate_dirichlet_probability, random_integer_log_scale
from clonal_trees.run_simulation.simulation_priors.print import print_simulation_parameters


# Constants
DEFAULT_MAX_EPOCHS = 8
DEFAULT_LAM = 1
DEFAULT_PROB_GD_ROUNDS = [0.4, 0.4, 0.2]
DEFAULT_ALPHA = [2, 2, 10]
MIN_RATE = 10000
MAX_RATE = 1000000


class EpochAssigner:
    """The class responsible for assigning values for pre, mid, and post epochs."""
    def __init__(self, num_gd_rounds: int, lam: float, max_epochs: int) -> None:
        """
        Initialize the EpochAssigner with the number of genome doubling rounds,
        Poisson parameter lambda, and max epochs.

        :param num_gd_rounds: Number of genome doubling rounds.
        :param lam: Poisson distribution parameter lambda for generating the number of epochs.
        :param max_epochs: Maximum number of epochs. In each epoch a single round of genome doubling or anueploidy can occur. 
        Furthermore, SNVs can occur in each epoch prior to that rounds genome doubling or anueploidy.
        """
        self.num_gd_rounds = num_gd_rounds
        self.prob_E = self._compute_prob_E(lam, max_epochs)
        self.pre, self.mid, self.post = None, None, None # TODO: why are these initialized to 0?
        self.max_epochs = max_epochs

    @staticmethod
    def _compute_prob_E(lam: float, max_epochs: int) -> np.ndarray:
        """Compute the Poisson probabilities and normalize them.

        :param lam: Poisson distribution parameter.
        :param max_epochs: Maximum number of epochs.
        :return: Normalized Poisson probabilities.
        """
        prob_E = [poisson.pmf(i, lam) for i in range(max_epochs+1)]  # it is +1 because max_epochs is inclusive
        prob_E = prob_E / np.sum(prob_E) # normalize because the poisson distribution here is truncated and don't sum to 1.
        return prob_E

    def _get_random_choice(self) -> int:
        """Returns a random choice for the number of epochs given the prior truncated poisson distribution."""
        return np.random.choice(range(len(self.prob_E)), p=self.prob_E)

    def assign_values(self) -> None:
        """Assign values to the pre, mid, and post epochs based on the number of genome doubling rounds."""
        if self.num_gd_rounds == 0:
            self._assign_values_for_no_GD_rounds()
        elif self.num_gd_rounds == 1:
            self._assign_values_for_one_GD_round()
        elif self.num_gd_rounds == 2:
            self._assign_values_for_two_GD_rounds()

    def _get_valid_random_choice(self, number_of_periods: int) -> Tuple[int, ...]:
        """Return valid random choices that do not exceed a given limit."""
        choices = tuple(self._get_random_choice() for _ in range(number_of_periods))
        while sum(choices) >= len(self.prob_E) - number_of_periods +1:
            choices = tuple(self._get_random_choice() for _ in range(number_of_periods))
        return choices

    def _assign_values_for_no_GD_rounds(self) -> None:
        """Assign values for the case when there are no genome doubling rounds."""
        self.pre, = self._get_valid_random_choice(number_of_periods=1)
        self.mid = self.post = -1

    def _assign_values_for_one_GD_round(self) -> None:
        """Assign values for the case when there is one genome doubling round."""
        self.pre, self.mid = self._get_valid_random_choice(number_of_periods=2)
        self.post = -1

    def _assign_values_for_two_GD_rounds(self) -> None:
        """Assign values for the case when there are two genome doubling rounds."""
        self.pre, self.mid, self.post = self._get_valid_random_choice(number_of_periods=3)


class SimulationArgs:
    """Class to encapsulate simulation arguments."""
    max_epochs: Optional[int] = None
    lam: Optional[float] = None
    pre: Optional[int] = None
    mid: Optional[int] = None
    post: Optional[int] = None
    p_up: Optional[float] = None
    p_down: Optional[float] = None
    rate: Optional[int] = None
    alpha: Optional[List[int]] = None
    total_epochs: Optional[int] = None
    prob_gd_rounds: Optional[List[float]] = None

    def simulate_parameters_not_given_as_arguments(self):
        """Simulates default parameters if not provided in the arguments when running one of the entry points."""
        self.max_epochs = self.max_epochs if self.max_epochs is not None else DEFAULT_MAX_EPOCHS

        self.lam = self.lam if self.lam is not None else DEFAULT_LAM

        num_gd_rounds = np.random.choice([0, 1, 2], p=self.prob_gd_rounds) if self.prob_gd_rounds is not None else np.random.choice([0, 1, 2], p=DEFAULT_PROB_GD_ROUNDS)

        if self.mid is None or self.post is None:
            epoch_assigner = EpochAssigner(num_gd_rounds, self.lam, self.max_epochs)
            epoch_assigner.assign_values()
            self.pre, self.mid, self.post = epoch_assigner.pre, epoch_assigner.mid, epoch_assigner.post

        self.total_epochs = self.pre + self.mid + self.post + 2

        if self.p_up is None or self.p_down is None:
            if self.alpha is None:
                self.alpha = DEFAULT_ALPHA
            probabilities = generate_dirichlet_probability(self.alpha)
            self.p_up, self.p_down, _ = probabilities

        if self.rate is None:
            self.rate = random_integer_log_scale(MIN_RATE, MAX_RATE)

        print_simulation_parameters(
            pre=self.pre,
            mid=self.mid,
            post=self.post,
            p_up=self.p_up,
            p_down=self.p_down,
            rate=self.rate
        )

