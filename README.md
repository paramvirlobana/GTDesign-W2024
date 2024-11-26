# Gas Turbine Design and Optimization Framework
Program written for the course: ```Gas Turbine Design``` ```Winter 2024``` ```Concordia University, Montreal```

## Description
A very basic turbine design cycle for the core of a gas turbine engine. Uses basic optimization methods to obtain the design that meets the requirements. The image below shows the turbine schematic used in the project presented via this code.

<div style="text-align: center;">
  <img src="report/figures/engine_core.png?raw=true" alt="XDSM Diagram" width="500"/>
</div>

<br />

Complete report can be found here: [Complete Report](report/Gas_Turbine_Design_2024.pdf)

## Usage
The program is used from the ```main.ipynb``` file present in the src code folder. This file contains the main iterative loop that iterates over a user specified design space and produces a design matrix. 

The program returns a large matrix with multiple possible design. The matrix is sorted in an order of preference that meets/maximizes most of the design requirements. 

In the current version of the program, the following parameters define the design space:
1. Reaction ($\Lambda$) -> 0.35 - 0.65
2. Flow Coefficient ($\phi$) -> 0.6 - 1.0
3. Structural Limit ($AN^2$) -> 10000000 - 22580600 (checking for a wide range)
4. Zweifel for vane -> 0.75 - 0.90
5. Zweifel for rotor -> 0.80 - 0.95

The current optimization is defined below:
```python
# OPTIMIZATION FUNCTION
data_efficiency['func_optimize'] = 0.4 * (data_efficiency['eta_opt_normalized']) + 0.3 * (data_meanline_losses['N_rotor_normalized']) + 0.2 * (data_efficiency['delta_eta_optimize_normalized'])
```

The weights of the optimization function can be set to meet a particular design requirement while trading off the others. NOTE: the sum of the weights should be 1. The image below shows the design and optimization framework.

![Alt text](report/figures/gt_xdsm.png?raw=true "XDSM Diagram")
