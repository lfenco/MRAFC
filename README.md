# Code for Nonlinear trajectory tracking with a 6DOF AUV using an MRAFC controller
Main code and functions required for the article "Nonlinear trajectory tracking with a 6DOF AUV using an MRAFC controller" published in IEEE Latin America Transactions. This code is composed of the following files:

1) mrafc_23112024.m
2) var_def.m
3) trackpath.m
4) trayectory_generation.m
5) initial_conditions.m
6) lmi_calculation.m
7) membership_function.m
8) update_data.m
9) lmiModel09102023
------------------------------------------------
## Description
1) (main code)
2) (definition of variables and parameters for the AUV model and reference model)
3) (function that allows to make a spiral trajectory based on 8 reference points)
4) (base function for the generation of straight line or circular arc type trajectories)
5) (definition of initial conditions for the AUV in the different controllers to be evaluated)
6) (Calculation of the linear matrix inequality for the linear subsystems of the AUV)
7) (definition of the membership functions and the working ranges)
8) (function for updating the parameters for the control law)
9) (database that saves the results obtained in the lmi_calculation.m function)

## Running the code
When running the main code, a series of input parameters will appear which will serve to indicate whether we want 
to perform a simulation when the system is in a noisy environment, or in an environment with external disturbance. 
If in both cases a value of zero is entered, the system identifies that it is evaluating an environment with ideal 
conditions. In addition, there are some additional intuitive parameters such as enabling saturation for the 
thrusters. Or even if what is desired is a simulation for regulation or trajectory tracking. So, below is a 
configuration for an ideal environment in trajectory tracking:

* input enable noise [Not noise = 0 -- noise = 1] : 0;
* input disturbance amplitude [0.8] : 0;
* input enable break [Not Break = 0 -- Break = 1]: 0;
* input Cnt Points: [8Points = 0 -- 16Points = 1] : 0;
* not saturation=0 -- saturation=1: 1;
* model reference: new calculation = 0 -- load data = 1: 1;
* input [regulation  = 0 -- tracking = 1] : 1;
