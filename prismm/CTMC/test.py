import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

# Function to create markov matrix
def create_markov_matrix(n):
    G = np.zeros((n, n))
    for i in range(n):
        if 2*i+1 < n:
            G[i, 2*i+1] = 1
    return G

# Function to create rate matrix for Continuous Time Markov Chain (CTMC)
def create_rate_matrix(n, lam, mu):
    Q = np.zeros((n, n))
    for i in range(n):
        if i != 0:
            Q[i, i-1] = lam
        if i != n-1:
            Q[i, i+1] = mu
        Q[i, i] = -np.sum(Q[i])
    return Q

# Define a function that calculates the negative log likelihood
def neg_log_likelihood(params, observations, n, model=1):
    # unpack the parameters
    s, t, lam, mu = params
    # Create the matrices and calculate P(t)
    G = create_markov_matrix(n)
    Q = create_rate_matrix(n, lam, mu)
    if model == 1:
        P_t = np.dot(np.dot(expm(t*Q), G), expm(s*Q))
    elif model == 2:
        P_t = expm(t*Q)
    # Compute the likelihood based on the observations and return the negative log likelihood
    likelihood = np.prod([P_t[1, i]**observations[i] for i in range(n)])  # assuming independent observations
    return -np.log(likelihood)

# Function to calculate the AIC
def calculate_aic(nll, k):
    return 2*k - 2*(-nll)  # Minus nll because it's the neg log likelihood

# Run the main loop
def main():
    # Matrix size
    n = 16

    # Base parameters
    base_s = 1
    base_t = 1.3
    base_lam = 2
    base_mu = 2

    s=0.2
    # Generate slightly varying values for each parameter
    hidden_s = np.random.normal(base_s, s)  # mean = base_s, std dev = 0.1
    hidden_t = np.random.normal(base_t, s)
    hidden_lam = np.random.normal(base_lam, s)
    hidden_mu = np.random.normal(base_mu, s)

    # Create the hidden matrices and calculate the hidden state probabilities
    hidden_G = create_markov_matrix(n)
    hidden_Q = create_rate_matrix(n, hidden_lam, hidden_mu)
    hidden_P_t = np.dot(np.dot(expm(hidden_t*hidden_Q), hidden_G), expm(hidden_s*hidden_Q))


    # Simulate the observations as multinomially distributed random variables
    total_count = 46
    observations = np.random.multinomial(total_count, hidden_P_t[1]/np.sum(hidden_P_t[1]))

    # Initial guess for the parameters
    initial_guess = [1, 1, 1, 1]

    # Bounds for the parameters (s, t, lam, mu)
    bounds = [(0, None), (0, None), (0, None), (0, None)]

    # Model 1: maximize the likelihood using the 'L-BFGS-B' method

    # Model 2: maximize the likelihood using the 'L-BFGS-B' method
    result1 = minimize(neg_log_likelihood, initial_guess, args=(observations, n, 1), bounds=bounds, method='L-BFGS-B')

    # Then update the minimize function calls in the main function to pass n:
    result2 = minimize(neg_log_likelihood, initial_guess, args=(observations, n, 2), bounds=bounds, method='L-BFGS-B')

    # Calculate AICs
    aic1 = calculate_aic(result1.fun, len(initial_guess))
    aic2 = calculate_aic(result2.fun, len(initial_guess))

    # Compare the models
    if aic1 < aic2:
        better_model = 1
    elif aic2 < aic1:
        better_model = 2
    else:
        better_model = 0  # both models are equally good

    return better_model


# Run the main loop for 1000 iterations and count how often each model is the best
counts = np.zeros(3)
for _ in range(100):
    result = main()
    counts[result] += 1

# Print the results
print('Model 1 is better: ', counts[1], ' times')
print('Model 2 is better: ', counts[2], ' times')
print('Both models are equally good: ', counts[0], ' times')



print("\nTODO:\n"
      "1. Update it to also infer the rate matrix from the data. This will require a different "
      "approach to the one above.\n"
      "2. Work out how to add in the SNV information.\n"
      "3. Is there a way to compare DTMC and CTMC models? The CTMC model is simpler and more easy "
      "to compute a likelihood on potentially, however we know that cell division is an inherently discrete "
      "process but there may be so many cellular divisions that it is computationally infeasable "
      "to do inference on the DTMC in which case a CTMC would be a good approximation."
     )


# thoughts on the CTMC model and how to use the informaiton from SNVs.

# I don't know if I even know how to do it in the DTMC case. I think I need to do some reading of my code again.