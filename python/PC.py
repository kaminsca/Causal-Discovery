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
    for j in range(tau_max):
        lag = tau_max - j
        data_time_lagged[var_index, j, :] = data[(tau_max - lag):(n - lag), var_index]

print(data_time_lagged)
print(data_time_lagged.shape)

# Loop over all variables
for var_index, var in enumerate(data_time_lagged):
    n = 0
    print(f'current variable shape: {var.shape}')
    # Find parents of current variable
    parent_vars = []
    for y_index, y in enumerate(data_time_lagged):
        for lag, var_lagged in enumerate(y):
            if lag == 0 and y_index == var_index:
                continue
            parent_vars.append(var_lagged)
    num_parents = len(parent_vars)
    print(num_parents)
    while num_parents >= n:
        for parent in parent_vars:
            # Select n dimensional subset of N(Xt)
            j = 1
            # Calculate conditional mutual information and partial correlation
            # Remove links below significance threshold
        n += 1