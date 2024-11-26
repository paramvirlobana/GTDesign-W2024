# This file creates the XDSM diagram for the preliminary design phase that will be followed in the aero capstone class.

from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT, RIGHT, METAMODEL

x = XDSM(use_sfmath=True)

## ADD ALL SYSTEMS VARIABLE FUNCTIONS
meanline = "meanline_design"
structural = "structural_design"
freevortex = "blade_vortex_analysis"
off_design = "off_design"
aero_losses = "aerodynamic_losses"
opt = "opt"
funcs = "funcs"

## ADD THE OPTIMIZER AND THE SOLVER
x.add_system(opt, OPT, r"\text{Optimizer}")

## ADD ALL THE SYSTEMS
x.add_system(meanline, FUNC, (r"\text{1:}",r"\text{Meanline Design}"))
x.add_system(structural, FUNC, (r"\text{2:}",r"\text{Structural Analysis}"))
x.add_system(freevortex, FUNC, (r"\text{3:}",r"\text{Blade Vortex Analysis}"))
x.add_system(off_design, FUNC, (r"\text{4:}",r"\text{Off Design}"))
x.add_system(aero_losses, FUNC, (r"\text{5:}",r"\text{Aerodynamic Losses}"))
x.add_system(funcs, METAMODEL, (r"\text{6:}",r"\text{Objective Function}"))

## ADD THE INPUTS FOR MEANLINE
x.add_input(opt, r"x^{(0)}")



x.connect(opt, off_design, r"4: i_{des}")
x.connect(opt, meanline, r"1: \psi, \alpha_3, M_3")
x.connect(opt, structural, r"2: AN^2")


x.connect(meanline, structural, r"2: y_2")
x.connect(structural, freevortex, r"3: y_3")
x.connect(meanline, freevortex, (r"3: y_2"))

x.connect(opt, aero_losses, r"5: \zeta_{rotor}, y_1")
x.connect(aero_losses, funcs, r"6: N_{rotor}, \zeta_{rotor}")
x.connect(meanline, aero_losses, r"5: y_2")
x.connect(off_design, aero_losses, r"5: y_5")


#x.connect()
x.connect(funcs, aero_losses, r"N_{rotor}")


x.connect(meanline, opt, r"\eta_{design}")
x.connect(meanline, funcs, r"\eta_{design}")

x.connect(off_design, opt, r"\eta_{off-design}")
x.connect(off_design, funcs, r"\eta_{off-design}")


x.connect(off_design, opt, r"i_{des}")
x.connect(structural, opt, r"AN^2")


x.connect(freevortex, aero_losses, r"5: y_4")
x.connect(structural, funcs, r"6: AN^2")
x.connect(aero_losses, freevortex, r"K_{T_{off-design}}")
x.connect(aero_losses, meanline, r"K_{T_{design}}")


x.add_output(off_design, r"\eta_{off-design}^{\star}")
x.add_output(meanline, r"\eta_{design}^{\star}")
x.add_output(off_design, r"i_{des}^{\star}")
x.add_output(structural, r"AN^2{\star}")

x.connect(funcs, opt, (r"f"))


# BACK TO THE MDA SOLVER
x.add_output(opt, r"\text{Optimized Design *}", side=LEFT)



## Write PDF
x.write("gt_design_XDSM")
