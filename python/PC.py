### PC Algorithm for Causal Discovery
### Clark Kaminsky (clarkkam@umich.edu), Mohammed Ombadi (ombadi@umich.edu); 2023

import pandas as pd
import numpy as np


# Setup hyperparameters and import data

tau_max = 5     # maximum time lag
alpha = 0.05    # significance threshold
q_max = 5       # maximum number of combinations

data = np.genfromtxt('Datasets/det.txt', delimiter=',')

# Extract the number of variables (columns) and observations (rows)
n, v = data.shape
p_max = v * tau_max # max number of parents

print(data[:5, :])
print(f'v: {v}, n: {n}')

# Set up time lagged variables
data_time_lagged = np.empty((v, tau_max + 1, n - tau_max))

for var_index in range(v):
    lag = 0
    for j in range(tau_max, 1, -1):
        data_time_lagged[var_index, j, :] = data[(tau_max - lag):(n - lag), var_index]
        lag += 1

print(data_time_lagged)



# Loop over all variables
    # Identify parents
    # Loop over the max number of parents, check if there are enough parents to consider
        #n = 0, N(Xt) = {Xt-1 ... Xt-taumax, Yt-1 ... Yt-taumax, ...}
        # While |N(Xt)| <= n
            # Select n dimensional subset of N(Xt)
            # Calculate conditional mutual information and partial correlation
    # Remove links below significance threshold