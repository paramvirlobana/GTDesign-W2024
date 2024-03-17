import pandas as pd

#2 - TRAILING EDGE LOSSES

#t/o is trailing edge opening to throat opening ratio
#Need to add Beta2 and alpha3 from output code

""""
To find TET from graph
""""

def figure_2_10(X_value, x=True):
    """
    x = true for t/o
    """
    fig_2_10 = pd.read_csv(r'_input_database\figure_2_10.csv')
    t/o = fig_2_10['x'].values
    y1 = fig_2_10['TET'.values]

    if x == True:
        interp_func = interpld(t/o, y1, kind='cubic')
        interpolated_y = interp_func(x_value)
        return interpolated_y

#After getting the TET Energy Coeff (delta_phi), we can calculate TET loss K_TE using:
def calculate_K_TE(M_2, delta_phi_squared_TE, gamma_g):
    term1 = (gamma_g - 1) / 2 * M_2**2
    term2 = 1 / (1 - delta_phi_squared_TE) - 1
    numerator = 1 - term1 * (term2)**(-gamma_g / (gamma_g - 1)) - 1
    denominator = 1 - (1 + term1)**(-gamma_g / (gamma_g - 1))
    K_TE = numerator / denominator
    return K_TE

    #beta2 and alpha3 values
    #B2 = 
    #a3 = 

    # values for M_2, delta_phi_squared_TE, and gamma
    #M_2 = 
    #delta_phi_squared_TE = 
    #gamma_g = 

if a3 != 0:
    # Calculate Δϕ^2 TE using the provided formula
    delta_phi_squared_TE = delta_phi_squared_TE_beta1_0 + abs(B2/a3) * (B2/a3) * (delta_phi_squared_TE_beta2_a3 - delta_phi_squared_TE_beta2_0)
else:
    print("Error: a3 cannot be zero.")

# Calculate K_TE using the defined function
K_TE = calculate_K_TE(M_2, delta_phi_squared_TE, gamma_g)


# The result is stored in delta_phi_squared_TE
print("Calculated Δϕ^2 TE:", delta_phi_squared_TE)
print("Calculated K_TE:", K_TE)
