"""

Created by Paramvir Lobana on 25 November, 2023

Gas Turbines and Turbomachinary Projects

This class is created for functions to be used for meanline turbine design

"""
import math
import matplotlib.pyplot as plt
import numpy
#-----------------------------------------------------------------
# The following class contains all the functions to be used for meanline analysis
#-----------------------------------------------------------------

class myfunctions():
    def __init__(self):
        pass

    def degtorad(deg):
        """
        Returns angle in radians
        ___________
        Input: Angle in degrees
        """
        rads = deg * (math.pi / 180)
        return rads
    
    def radtodeg(rad):
        """
        Returns angle in degrees
        ___________
        Input: Angle in radians
        """
        degs = rad * (180/math.pi)
        return degs
    
    def calc_psi(c_p, delta_T_0_s, U):
        """
        This function calculates the blade loading coefficient
        ___________
        Input: c_p, delta_T_0_s, U
        """
        psi = (2 * c_p * delta_T_0_s)/(U**2)
        print("psi = ", psi)
        return psi 
    
    def calc_lambda(alpha_3, phi, psi):
        """
        This function calculates the degree of reaction
        ___________
        Input: Swirl angle "alpha_3", flow coefficient "phi", blade loading coefficient "psi"
        ___________
        Output: Returns the value of the degree of coefficient and the tan_beta_3
        """
        alpha_3 = myfunctions.degtorad(alpha_3)
        tan_beta_3 = math.tan(alpha_3) + (1/phi)
        lambda_val = (4 * tan_beta_3 * phi - psi)/4
        print("lambda = ", lambda_val)
        return lambda_val, tan_beta_3
    
    def calc_angles(tan_beta_3, psi, phi, lambda_val):

        """
        This function required multiple angles for the velocity triangles
        ___________
        Input: tan_beta_3, psi, phi, lambda_val
        ___________
        Output: A plot showing stability in real and imaginary axis.
        """
        beta_3 = myfunctions.radtodeg(numpy.arctan(tan_beta_3))
        print("beta_3 = ", beta_3)

        beta_2 = myfunctions.radtodeg(numpy.arctan((1/(2*phi)) * (0.5 * psi - 2*lambda_val)))
        print("beta_2 = ", beta_2)

        tan_alpha_2 = numpy.tan(math.radians(beta_2)) + (1 / phi)
        alpha_2 = myfunctions.radtodeg(numpy.arctan(tan_alpha_2))
        print("alpha_2 = ", alpha_2)

        return beta_3, beta_2, alpha_2
    
    def calc_temperature(M, T_total):
        """
        Calculates static temperature for given fluid.
        ___________
        Input: Total Temperature "T_total", Mach number "M"
        ___________
        Output: Static Temperature "T"
        """
        scratch = (gamma_t - 1)/2
        T = T_total / (1 + scratch*(M)**2)
        print("T_static = ", T)
        return T
    
    def calc_a(M, T_static):
        a = math.sqrt(gamma_t*R*T_static*1000)
        v = M*a
        print("v = ", v)
        return a, v
    
    def calc_components(alpha, C):
        # convert alpha to degrees
        alpha = math.radians(alpha)
        c_w = C*math.sin(alpha)
        c_a = C*math.cos(alpha)
        print("c_w = ", c_w)
        print("c_a = ", c_a)
        return c_a, c_w


    def calc_properties(T_t_1, p_t_1, c, c_a,n):
        """
        This function calculate air properties and also the area from mass flow rate
        ___________
        Input: Total Temperature "T_t_1", Total Pressure "p_t_1", Flow velocity "c", axial component of the velocity "c_a", station number "n"
        ___________
        Output: A plot showing stability in real and imaginary axis.
        """
        T_static =  T_t_1 - (c**2)/(2 * c_p_t)
        p_static =  p_t_1 / ((T_t_1/T_static)**(gamma_t/(gamma_t - 1)))
        rho_ = p_static/(R*T_static)
        A = mass_flow_rate/(rho_* c_a)
        print("T_",n, " = ", T_static)
        print("P_",n, " = ", p_static)
        print("rho_",n, " = ", rho_)
        print("A_",n, " = ", A)
        return T_static, p_static, rho_, A

         
    def velocity_vectors(U, phi, alpha_2):
        c_a_2 = U*phi
        print("c_a_2 =", c_a_2)
        c_2 = (c_a_2)/(math.cos(myfunctions.degtorad(alpha_2)))
        print("c_2 ", c_2)

        # Computation of temperatures
        T_02_T_2 = (c_2**2)/(2 * c_p_t)
        #print("T_02-T_2 =", T_02_T_2)
        T_2 = T_t_1 - T_02_T_2
        print("T_2 ", T_2)
        T_2_dash = T_2 - 0.05*T_02_T_2
     
        # Computation of P_2
        p_2 = p_t_1 / ((T_t_1/T_2_dash)**(gamma_t/(gamma_t - 1)))
        print("p_2 =", p_2)
        # Computation of actual pressure ratio
        p_01_p_2 = p_t_1/p_2
        #print("p_01/p_2 =", p_01_p_2)

        # Computation of Critical pressure ratio
        p_01_p_c = ((gamma_t + 1)/2)**(gamma_t/(gamma_t-1))
        #print("p_01/p_c =", p_01_p_c)
        # Computation of density
        rho_2 = p_2/(R * T_2)
        print("rho_2 =", rho_2)
        # Computation of area at section 2
        A_2 = mass_flow_rate/(rho_2*c_a_2)
        print("A_2 =", A_2)
        # Nozzle cross section 
        A_2_N = A_2 * math.cos(myfunctions.degtorad(alpha_2))
        print("area_N ", A_2_N)

        c_a_3 = c_a_2
        c_a_1 = c_a_3/math.cos(myfunctions.degtorad(alpha_3))
        print("c_a_1 = ", c_a_1)
        
        
        c_1 = c_a_1
        T_1, p_1, rho_1, A_1 = myfunctions.calc_properties(T_t_1, p_t_1, c_1, c_a_1, 1)

        c_3 = c_a_1
        T_03 = T_t_1 - turbine_temperature_drop
        p_03 = p_t_1/p_ratio_t
        T_3, p_3, rho_3, A_3 = myfunctions.calc_properties(T_03, p_03, c_1, c_a_1, 3)


        return A_1, A_2, A_3
    
    def calc_height(A, n):
        """
        This calculated geometric values related to the cross section of the turbine
        ___________
        Input: area "A", station number "n"
        ___________
        Output: A plot showing stability in real and imaginary axis.
        """
        h = (revs * A)/(U)
        print("A_",n,"=", h)
        return h

    def velocityvectors(U, C, alpha_1, alpha_2, beta_2, beta_3, plot):
        #U = 340     
        #C = 200

        alpha_1 = myfunctions.degtorad(alpha_1)             # entry angle in degrees 
        c_1 = C         # inlet velocity         
        c_a_1 = c_1 * math.cos(alpha_1) 
        c_w_1 = c_1 * math.sin(alpha_1)

        #============================================
        # evaluation of all velocities at point 2
        #============================================
        c_a_2 = c_a_1
        alpha_2 = myfunctions.degtorad(alpha_2)
        beta_2 = myfunctions.degtorad(35)
        c_2 = (c_a_2)/(math.cos(alpha_2))
        c_w_2 = c_a_2 * math.tan(alpha_2)

        V_w_2 = c_w_2 - U
        V_2 = V_w_2 / (math.sin(beta_2))

        #============================================
        # evaluation of all velocities at point 3
        #============================================
        alpha_3 = alpha_1
        beta_3 = myfunctions.degtorad(55)
        c_a_3 = c_a_1
        V_3 = c_a_3/(math.cos(beta_3))
        V_w_3 = V_3 * math.sin(beta_3)
        c_w_3 = V_w_3 - U
        c_3 = c_a_3/math.cos(alpha_3)


        print(c_1)
        print(c_2)
        print(c_3)

        print(V_2)
        print(V_3)

        #============================================
        # PLOTTING THE VECTORS
        #============================================

        # Function to plot a vector with a labeled offset from the arrow tip
        def plot_vector(ax, origin, vector, color, label, offset=25, fontsize=10):
            ax.quiver(*origin, *vector, angles='xy', scale_units='xy', scale=1, color=color)
            label_position = (origin[0] + vector[0] + offset, origin[1] + vector[1] + offset)
            ax.text(*label_position, label, color=color, ha='center', va='center', fontsize=fontsize)

        # Create subplots
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        fig.suptitle('Velocity Vectors at Different Points')

        plot_vector(axes[0], (0, 0), (c_w_1, c_a_1), 'r', r'$\vec{C}_1$', fontsize=12)
        plot_vector(axes[0], (0, 0), (c_w_1, 1), 'b', r'$\vec{c}_{w_1}$', fontsize=12)
        plot_vector(axes[0], (0, 0), (1, c_a_1), 'b', r'$\vec{c}_{a_1}$', fontsize=12)
        axes[0].set_title('Before Stator')
        axes[0].set_xlim([-50, 600])
        axes[0].set_ylim([-50, 600])
        axes[0].set_xlabel(r'$\theta \; [m/s]$')
        axes[0].set_ylabel(r'$x\; [m/s]$ ')

        plot_vector(axes[1], (0, 0), (-1*c_w_2, c_a_2), 'r', r'$\vec{C}_1$', fontsize=12)
        plot_vector(axes[1], (0, 0), (-1*c_w_2, 1), 'b', r'$\vec{c}_{w_1}$', fontsize=12)
        plot_vector(axes[1], (0, 0), (1, c_a_2), 'b', r'$\vec{c}_{a_1}$', fontsize=12)
        axes[1].set_title('Between Stator and Rotor')
        axes[2].set_xlim([-50, 600])
        axes[2].set_ylim([-50, 600])
        axes[1].set_xlabel(r'$\theta \; [m/s]$')
        axes[1].set_ylabel(r'$x\; [m/s]$ ')

        plot_vector(axes[2], (0, 0), (c_w_3, c_a_3), 'r', r'$\vec{C}_1$', fontsize=20)
        plot_vector(axes[2], (0, 0), (c_w_3, 1), 'b', r'$\vec{c}_{w_1}$', fontsize=12)
        plot_vector(axes[2], (0, 0), (1, c_a_3), 'b', r'$\vec{c}_{a_1}$', fontsize=12)
        axes[2].set_title('After rotor')
        axes[2].set_xlim([-50, 600])
        axes[2].set_ylim([-50, 600])
        axes[2].set_xlabel(r'$\theta \; [m/s]$')
        axes[2].set_ylabel(r'$x\; [m/s]$ ')
        # Shows the plots only if the user wants to output plot
        if plot==True:
            plt.show()

        else:
            pass
#-----------------------------------------------------------------
# END OF FUNCTIONS
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# INPUT VALUES
#-----------------------------------------------------------------
mass_flow_rate = 20                     # [kg/s]
turbine_isentropic_efficiency = 0.9     # [/]
T_01 = 1100                             # [k] Turbine inlet temperature
turbine_temperature_drop = 145          # [k]
p_ratio_t = 1.87                        # [/]
p_01 = 400                              # [kPa] Turbine inlet total pressure 
revs = 250                              # [revs/s]
U = 340                                 # [m/s]  
M_0 = 0                                 # Flight Mach Number
p_0 = 101                               # [kPa] Ambient Pressure
T_0 = 25 + 237.15                       # [K] Ambient Temperature
R = 0.287
# For turbine section
c_p_t = 1243.67                         # [J/kg.K]
gamma_t = 1.30                          # For combusted gas in the combustion chamber.
#-----------------------------------------------------------------
# ASSUMPTION BLOCK TO MODIFY
alpha_3 = 5                            # [deg]
phi = 0.8                               # [/]
#-----------------------------------------------------------------
print("===============================================")
print("===================STAGE 2=====================")
print("===============================================")
# SINCE WE ARE FOLLOWING CONSTANT AXIAL VELOCITY DESIGN APPROACH
# We use the same axial velocity as c_a
# alpha_1 is chosen to be 0
v_1 = phi * U
T_1 = 1500 - ((v_1**2)/(2 * c_p_t))
a_1 = math.sqrt(gamma_t*R*T_1*1000)
M_1 = v_1/a_1
print(M_1)
print(T_1)
c_a, c_w = myfunctions.calc_components(60, C_2)
print("===============================================")


print("===============================================")
print("===================STAGE 2=====================")
print("===============================================")


T_2 = myfunctions.calc_temperature(1.1, 1500) # CHANGE
# Assume M_2 = 1.1
# absolute gas velocity at this stage
a_2, C_2 = myfunctions.calc_a(1.1, T_2)

# First assumption exit angle ~60 [deg]
# Calculate the angles 
c_a, c_w = myfunctions.calc_components(60, C_2)
print("===============================================")





"""
psi = myfunctions.calc_psi(c_p_t, turbine_temperature_drop, U)
lambda_val, tanbeta_3 = myfunctions.calc_lambda(alpha_3, phi, psi)

beta_3, beta_2, alpha_2 = myfunctions.calc_angles(tanbeta_3, psi, phi, lambda_val)
A_1, A_2, A_3 = myfunctions.velocity_vectors(U, phi, alpha_2)



r_m = U/(2*math.pi*revs)
h_1 = myfunctions.calc_height(A_1, 1)
h_2 = myfunctions.calc_height(A_2, 2)
h_3 = myfunctions.calc_height(A_3, 3)
#y = myfunctions.velocityvectors(U, 500, 10, alpha_2, beta_2, beta_3, plot=True)
"""