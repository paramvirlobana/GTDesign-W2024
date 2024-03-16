import pandas as pd
import numpy as np

# Define the given formulas as functions

def calc_Y_pAMDC(row):
    term1 = row['Y_p_beta0'] + abs(row['beta1']/row['alpha2']) * (row['beta1']/row['alpha2']) * row['Y_p_beta_alpha2'] - row['Y_p_beta0']
    term2 = (row['t_MAX']/row['C'] / 0.2) ** (-row['beta1']/row['alpha2'])
    return term1 * term2

def calc_K_accl(row):
    return 1 - row['K_2'] * (1 - row['K_1'])

def calc_K_p(row):
    return 0.914 * (2/3 * row['Y_pAMDC'] * row['K_accl'] + row['K_sh'])

def calc_shock_dp_q1(row):
    return (row['R_H']/row['R_T']) * row['dp_q1_HUB']

def calc_Y_SHOCK(row):
    term1 = (1 - (1 + (row['gamma'] - 1)/2 * row['M1']**2) ** (-row['gamma']/(row['gamma'] - 1)))
    term2 = (1 - (1 + (row['gamma'] - 1)/2 * row['M2']**2) ** (-row['gamma']/(row['gamma'] - 1)))
    return row['dp_q1_SHOCK'] * (row['P1']/row['P2']) * (term1/term2)

# Create a DataFrame for input parameters
# Replace this with your actual input data

input_data = {
    'Y_p_beta0': [0.02],
    'beta1': [0.5],
    'alpha2': [1.0],
    'Y_p_beta_alpha2': [0.03],
    't_MAX': [0.1],
    'C': [0.02],
    'K_1': [0.05],
    'K_2': [0.08],
    'K_sh': [0.04],
    'R_H': [1.1],
    'R_T': [1.2],
    'dp_q1_HUB': [0.06],
    'P1': [101325],
    'P2': [90000],
    'gamma': [1.4],
    'M1': [0.8],
    'M2': [0.6]
}

df = pd.DataFrame(input_data)

# Calculate the parameters and add them to the DataFrame

df['Y_pAMDC'] = df.apply(calc_Y_pAMDC, axis=1)
df['K_accl'] = df.apply(calc_K_accl, axis=1)
df['dp_q1_SHOCK'] = df.apply(calc_shock_dp_q1, axis=1)
df['Y_SHOCK'] = df.apply(calc_Y_SHOCK, axis=1)
df['K_p'] = df.apply(calc_K_p, axis=1)

# Print the DataFrame with the results
print(df[['Y_pAMDC', 'K_accl', 'dp_q1_SHOCK', 'Y_SHOCK', 'K_p']])


# For iterations:

from itertools import product

# Define your range for Y_p_beta0 and any other parameters
Y_p_beta0_range = np.linspace(0.02, 1.0, num=50)  # 50 evenly spaced values from 0.02 to 1.0
# Add ranges for other parameters as needed
# ...

# Use itertools.product to create a Cartesian product of all parameter combinations
parameter_combinations = list(product(Y_p_beta0_range, ...))  # replace '...' with other parameter ranges

# Now convert this into a DataFrame
df = pd.DataFrame(parameter_combinations, columns=['Y_p_beta0', ...])  # Add the other parameter names

# Continue as before, calculating the parameters and adding them to the DataFrame
