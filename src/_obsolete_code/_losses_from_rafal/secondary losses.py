import pandas as pd
import numpy as np


#3 – Secondary Losses


#- Y_s: Overall loss coefficient.
#- K_s: Loss coefficient related to secondary flow.
#- K_3: Constant related to secondary flow losses.
#- K_p: Profile loss coefficient.
#- Y_{S,AMDiC}: Profile loss coefficient specific to the axial mean diameter chord.
#- f_{(AR)}: Aspect ratio function.
#- C_l: Lift coefficient.
#- s/c: Blade spacing to chord ratio.
#- h/c: Blade height to chord ratio.



#how to calculate Cl_sc used in Y_AMDiC
def calculate_Cl_sc (alpha2, alpha3, alpha_m):
    tan_alpha1 = math.tan(alpha2)
    tan_alpha2 = math.tan(alpha3)
    cos_alpha_m = math.cos(alpha_m)
    Cl_sc = 2 * (tan_alpha2 + tan_alpha3) * cos_alpha_m
    return Cl_sc

# Function to calculate alpha_m
def calculate_alpha_m(alpha2, alpha3):
    return math.atan(0.5 * (math.tan(alpha2) - math.tan(alpha3)))

# Function to calculate f(AR)
def calculate_f_AR(h_over_c):
    if h_over_c <= 2:
        return (1 - 0.25 * np.sqrt(2 - h_over_c)) / h_over_c
    else:
        return 1 / h_over_c

# Define the function to calculate secondary loss coefficient Y_s
def calculate_Y_s(Y_sAMDiC, K_s):
    return 1.2 * Y_sAMDiC * K_s

# Define the function to calculate Y_sAMDiC
def calculate_Y_sAMDiC(f_AR, cos_alpha3, Cl_sc)**2, cos_beta2, cos_alpha_m):
    return 0.0334 * (f_AR * (cos_alpha3 / cos_beta2) * (Cl_sc)**2 * (cos_alpha3**2 / cos_alpha_m**3))

# Define the function to calculate K_s
def calculate_K_s(K_3, K_accel):
    return 1 - K_3 * (1 - K_accel)
#K_3 comes from below, from figure 2_16
def calculate_Cl_sc (alpha2, alpha3, alpha_m):
    tan_alpha2 = math.tan(alpha2)
    tan_alpha3 = math.tan(alpha3)
    cos_alpha_m = math.cos(alpha_m)
    Cl_sc = 2 * (tan_alpha2 + tan_alpha3) * cos_alpha_m
    return Cl_sc
  


""""
To find K_3 from graph
""""

def figure_2_16(X_value, x=True):

""""
x = true for b_x/h
""""
fig_2_16 = pd.read_csv(r'_input_database\figure_2_16.csv')
b_x/h = fig_2_16['x'].values
y1 = fig_2_16['K_3'.values]

if x == True:
    interp_func = interpld(b_x/h, y1, kind='cubic')
    interpolated_y = interp_func(x_value)
return interpolated_y

# Calculate Y_sAMDiC, K_s, Y_s, alpha_m, f_AR
df['Y_sAMDiC'] = df.apply(lambda row: calculate_Y_sAMDiC(row['f_AR'], np.cos(row['alpha3']), 
                                                          row['Cl_sc)'], np.cos(row['beta2']), 
                                                          np.cos(row['alpha_m'])), axis=1)
df['K_s'] = df.apply(lambda row: calculate_K_s(row['K_3'], row['K_accel']), axis=1)
df['Y_s'] = df.apply(lambda row: calculate_Y_s(row['Y_sAMDiC'], row['K_s']), axis=1)
df['alpha_m'] = df.apply(lambda row: calculate_alpha_m(row['alpha2'], row['alpha3']), axis=1)
df['f_AR'] = df.apply(lambda row: calculate_f_AR(row['h_over_c']), axis=1)

print(df)
