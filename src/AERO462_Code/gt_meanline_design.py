import math
import numpy as np

mass_flow_rate = 20                     # [kg/s]
turbine_isentropic_efficiency = 0.9     # [/]
T_01 = 1100                             # [k] Turbine inlet temperature
turbine_temperature_drop = 145          # [k]
p_ratio_t = 1.87                        # [/]
p_01 = 400                              # [kPa] Turbine inlet total pressure 
revs = 250                              # [revs/s]
U = 340                                 # [m/s]   

p_0 = 101                               # [kPa] Ambient Pressure
T_0 = 25 + 273.15                       # [K] Ambient Temperature
R = 287
#-----------------------------------------------------------------
# For turbine section
c_p_t = 1243.67                         # [J/kg.K]
gamma_t = 1.30                          # For combusted gas in the combustion chamber.
#==================================================================
# Design Choice Input
phi =  0.9

def calc_psi(c_p, delta_T_0_s, U):
    """
    This function calculates the blade loading coefficient
    ___________
    Input: c_p, delta_T_0_s, U
    """
    psi = (2 * c_p * delta_T_0_s)/(U**2)
    print("psi = ", psi)
    return psi 

def calc_lambda(alpha, phi, psi):
    """
    This function calculates the degree of reaction
    ___________
    Input: Swirl angle "alpha", flow coefficient "phi", blade loading coefficient "psi"
    ___________
    Output: Returns the value of the degree of coefficient and the tan_beta_3
    """
    alpha = math.radians(alpha)
    tan_beta_3 = math.tan(alpha) + (1/phi)
    lambda_val = (4 * tan_beta_3 * phi - psi)/4
    print("lambda = ", lambda_val)
    return lambda_val

def calc_T_static(T_total, M):
    T_static = T_total/(1 + 0.15*(M**2))
    return T_static


def calc_area(T_static, p_static, C_a, n):
    """
    This function calculate air properties and also the area from mass flow rate
    ___________
    Input: Total Temperature "T_01", Total Pressure "p_01", Flow velocity "c", axial component of the velocity "c_a", station number "n"
    ___________
    Output: A plot showing stability in real and imaginary axis.
    """
    rho_ = p_static/(R*T_static)
    A = mass_flow_rate/(rho_* C_a)/1000
    print("rho_",n, " = ", rho_ * 1000)
    print("A_",n, " = ", A)
    return A

def calc_height(A, r_m, n):
    """
    This calculated geometric values related to the cross section of the turbine
    ___________
    Input: area "A", station number "n"
    ___________
    Output: A plot showing stability in real and imaginary axis.
    """
    h = (revs * A)/(U)
    rtrm = (r_m + 0.5*h)/(r_m - 0.5*h)
    print("h_",n,"=", h)
    print("rtrm_",n,"=", rtrm)
    return h

def calc_freevortex_nozzle(r_m, r_r, r_t, alpha):
    tan_alpha_r_2 = (r_m / r_r) * math.tan(math.radians(alpha))
    tan_alpha_t_2 = (r_m / r_t) * math.tan(math.radians(alpha))

    alpha_2_r_fv = np.arctan(tan_alpha_r_2)
    alpha_2_t_fv = np.arctan(tan_alpha_t_2)
    print("nozzle alpha values: root, tip")
    print(math.degrees(alpha_2_r_fv), math.degrees(alpha_2_t_fv))

    beta_2_r_fv = math.degrees(np.arctan(tan_alpha_r_2 - ((r_m / r_r) * (1/phi))))
    beta_2_t_fv = math.degrees(np.arctan(tan_alpha_t_2 - ((r_m / r_t) * (1/phi))))
    print("nozzle beta values: root, tip")
    print(beta_2_r_fv, beta_2_t_fv)

    return math.degrees(alpha_2_r_fv), beta_2_r_fv

def calc_freevortex_rotor(r_m, r_r, r_t, alpha):
    tan_alpha_r_3 = (r_m / r_r) * math.tan(math.radians(alpha))
    tan_alpha_t_3 = (r_m / r_t) * math.tan(math.radians(alpha))

    alpha_3_r_fv = math.degrees(np.arctan(tan_alpha_r_3))
    alpha_3_t_fv = math.degrees(np.arctan(tan_alpha_t_3))
    print("rotor alpha values: root, tip")
    print(alpha_3_r_fv, alpha_3_t_fv)

    beta_3_r_fv = math.degrees(np.arctan(tan_alpha_r_3 + ((r_m / r_r) * (1/phi))))
    beta_3_t_fv = math.degrees(np.arctan(tan_alpha_t_3 + ((r_m / r_t) * (1/phi))))
    print("rotor beta values: root, tip")
    print(abs(beta_3_r_fv), abs(beta_3_t_fv))

    return alpha_3_r_fv, beta_3_r_fv


# 1st Design Choice:M_2 = 1.1 M


def gt_design(M, M_3r):

    # Turbine Nozzle Pitchline Design
    T_2 = calc_T_static(T_01, M)
    print("T_2 = ", T_2)

    # The speed of sound in station 2 is
    a_2 = math.sqrt(gamma_t * R * T_2)
    print("a_2 = ", a_2)

    # absolute gas speed at the nozzle exit
    C_2 = a_2 * M
    print("C_2 = ", C_2)

    # Calculate the constant axial velocity
    C_a = U * phi
    print("C_a = ", C_a)

    # Calculate alpha_2
    alpha_2 = math.degrees(np.arccos(C_a / C_2))
    print("alpha_2 = ",alpha_2)

    # Calculate c_w_2
    C_w_2 = C_2 * math.sin(math.radians(alpha_2))
    print("C_w_2 = ",C_w_2)

    # Calculate beta_2
    beta_2 = math.degrees(np.arctan((C_w_2 - U) / C_a))
    print("beta_2 = ",beta_2)

    V_2 = C_a / math.cos(math.radians(beta_2))
    print("V_2 = ",V_2)

    # exit mach number = 0.85 -> design choice
    T_3 = calc_T_static(T_01 - turbine_temperature_drop, M_3r)
    print("T_3 = ",T_3)

    a_3 = math.sqrt(gamma_t * R * T_3)
    print("a_3 = ", a_3)

    # absolute gas speed at the rotor exit
    V_3 = a_3 * M_3r
    print("V_3 = ", V_3)

    # calculate beta_3
    beta_3 = math.degrees(np.arccos(C_a / V_3))
    print("beta_3 = ",beta_3)
    
    print("C_w_3 = ", V_3 * math.sin(math.radians(beta_3))-U)

    alpha_3 = math.degrees(np.arctan(abs((V_3 * math.sin(math.radians(beta_3))-U)/(C_a))))
    print("alpha_3 = ",alpha_3)
    C_3 = C_a / (math.cos(math.radians(alpha_3)))
    print("C_3 = ",C_3)
    psi = calc_psi(c_p_t, turbine_temperature_drop, U)
    lambda_val = calc_lambda(alpha_3, phi, psi)

    T_1 = T_01 - ((C_a/(math.cos(math.radians(alpha_3))))**2/(2*c_p_t))
    print("T_1 = ", T_1)

    p_1 = p_01 /(T_01/T_1)**(gamma_t/(gamma_t - 1))
    print("p_1 = ",p_1)

    A_1  = calc_area(T_1, p_1, C_a, 1) 

    p_2 = (p_01 )/(T_01/T_2)**(gamma_t/(gamma_t - 1))
    #print(p_2)
    A_2  = calc_area(T_2, p_2, C_a, 2)


    p_3 = (p_01/1.873)*(T_3/955)**(gamma_t/(gamma_t - 1))
    print("p_3 = ", p_3)
    A_3  = calc_area(T_3, p_3, C_a, 3)  
    r_m = U/(2*math.pi*revs)
    h_1 = calc_height(A_1, r_m, 1)
    h_2 = calc_height(A_2, r_m, 2)
    h_3 = calc_height(A_3, r_m, 3)

    r_t2 = r_m + 0.5*h_2
    r_t3 = r_m + 0.5*h_2

    r_r2 = r_m - 0.5*h_2
    r_r3 = r_m - 0.5*h_2

    alpha_2r, beta_2r = calc_freevortex_nozzle(r_m, r_r2, r_t2, alpha_2)
    alpha_3r, beta_3r = calc_freevortex_rotor(r_m, r_r3, r_t3, alpha_3)
    lambda_cr = calc_lambda(alpha_3r, phi, psi)

    # Calculation of M_V2
    V_2r = C_a * (1 / np.cos(math.radians(beta_2r)))
    print("V_2r = ",V_2r)
    C_2r = C_a * (1 / np.cos(math.radians(alpha_2r)))
    print("C_2r = ",C_2r)

    T_2r = T_01 - ((C_2r**2)/(2 * c_p_t))
    print("T_2r = ",T_2r) 

    M_V2_r = V_2r / (math.sqrt(gamma_t * R * T_2r))
    print("M_V2_r = ",M_V2_r)
    p_03 = p_01 / p_ratio_t
    # Aero Losses
    Y_r = (p_01 - p_03)/(p_03 - p_3)
    print("Y_r = ",Y_r)

    # Stress Analysis
    r_t_stress = r_m + 0.5*h_3
    U_t = 2 * math.pi *250*r_t_stress
    sigma_max_straight = 8470 * (0.5 * U_t**2)*(1 - (r_r3/r_t3)**2) * 10**(-6)
    print("sigma_straight = ", sigma_max_straight)
    print("sigma_tapered = ", (2/3)*sigma_max_straight)



##########################################
# RUN THE DESIGN FUNCTION

gt_design(1.1, 0.875)