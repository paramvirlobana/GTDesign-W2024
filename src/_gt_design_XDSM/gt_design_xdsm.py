# This file creates the XDSM diagram for the preliminary design phase that will be followed in the aero capstone class.

from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT, RIGHT

x = XDSM(use_sfmath=True)

## ADD ALL SYSTEMS VARIABLE FUNCTIONS
meanline = "meanline_design"
structural = "structural_design"
freevortex = "blade_vortex_analysis"
off_design = "off_design"
aero_losses = "aerodynamic_losses"
opt = "opt"
sol = "solver"

## ADD THE OPTIMIZER AND THE SOLVER
x.add_system(opt, OPT, r"\text{Optimizer}")
x.add_system(sol, SOLVER, r"\text{MDA}")




## ADD ALL THE SYSTEMS
x.add_system(meanline, FUNC, (r"\text{Analysis - 1}",r"\text{Meanline Design}"))
x.add_system(structural, FUNC, (r"\text{Analysis - 2}",r"\text{Structural Analysis}"))
x.add_system(freevortex, FUNC, (r"\text{Analysis - 3}",r"\text{Blade Vortex Analysis}"))
x.add_system(off_design, FUNC, (r"\text{Analysis - 4}",r"\text{Off Design}"))
x.add_system(aero_losses, FUNC, (r"\text{Analysis - 5}",r"\text{Aerodynamic Losses}"))

## ADD THE INPUTS FOR MEANLINE
x.add_input(sol, (r"\Lambda, \psi", r"\alpha_3, M_{exit}"))


x.connect(aero_losses, opt, (r"\eta_0", r"\eta_{0_{off}}"))
x.connect(meanline, structural, r"A_2, U")
x.connect(structural, freevortex, r"r_t, r_h, r_m")
x.connect(meanline, freevortex, (r"\phi_2, \phi_3, \alpha_2", r"U, C_{a_2}"))

# BACK TO THE MDA SOLVER
x.connect(meanline, sol, (r"\phi_2, \phi_3, \beta_2, \alpha_2", r"A_2, A_3, C_1, C_2, C_3"))
x.connect(structural, sol, (r"N, \omega, r_t", r"r_h, r_m, h"))
x.connect(freevortex, sol, (r"\alpha_{2h}, \alpha_{2t}, \beta_{2h}, \beta_{2t}", r"\alpha_{3h}, \alpha_{3t}, \beta_{3h}, \beta_{3t}"))




x.add_input(meanline, (r"\text{Cycle Calculations}", r"T_{01}, T_{02}, T_{03}", r"P_{01}, P_{03}"))
x.add_output(opt, r"\text{Optimized Design *}", side=LEFT)



## Write PDF
x.write("gt_design_XDSM")
