#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# carb.chem
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
A simple model for exploring the links between the sedimentary rock cycle, electron transfer, and carbonate chemistry in the oceans [DOI TEST]
#
The model solves for ocean-atmosphere carbon cycling and marine carbonate chemistry based on parameterized fluxes through the sedimentary rock cycle.  All model equations and parameters are described in detail in the Supplementary Information to C.T. Reinhard and W.W. Fischer [2020] "Mechanistic links between the sedimentary rock cycle and marine acid-base chemistry", Geochemistry, Geophysics, Geosystems, doi:.  These scripts will spin up the model, perform the time-dependent perturbation analyses, and reproduce all figures from C.T. Reinhard + W.W. Fischer [2020]. 
#
As a courtesy, we request that others who use this code please cite Reinhard + Fischer [2020].  We also request that those who use/modify the code please send publications and/or modified code to CTR [chris.reinhard@eas.gatech.edu].
#
REQUIREMENTS: Written and tested on Matlab R2018a/R2018b/R2019a [but should be fully compatible with other versions]
#
# Scripts for running the model include:
- carb_chem             --> main model script
- make_clean		--> deletes all existing spinup files, model output, and figures
- make_spinup		--> sequentially spins up the model scenarios below from cold; stored in 'spinup'
- make_run		--> sequentially runs model from spinup; results stored in 'output'
- make_plot		--> loads output and makes all figures from Reinhard + Fischer [2020]; plots stored in 'figures'
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# input scenarios
- 111-113	--> doubled pyrite weathering                                        | f_py_CaCO3  = 1.0
- 114-116       --> doubled pyrite weathering                                        | f_py_CaSiO3 = 1.0
- 121-123       --> doubled pyrite and organic carbon weathering                     | f_py_CaCO3  = 1.0
- 124-126       --> doubled pyrite and organic carbon weathering                     | f_py_CaSiO3 = 1.0
- 131-133       --> 50% increase in weathering of all crustal phases                 | f_py_CaCO3  = 1.0 
- 134-136       --> 50% increase in weathering of all crustal phases                 | f_py_CaSiO3 = 1.0 
- 211-213       --> doubled pyrite weathering, power function                        | f_py_CaCO3  = 1.0
- 214-216       --> doubled pyrite weathering, power function                        | f_py_CaSiO3 = 1.0
- 221-223       --> doubled pyrite and organic carbon weathering, power function     | f_py_CaCO3  = 1.0
- 224-226       --> doubled pyrite and organic carbon weathering, power function     | f_py_CaSiO3 = 1.0
- 231-233       --> 50% increase in weathering of all crustal phases, power function | f_py_CaCO3  = 1.0
- 234-236       --> 50% increase in weathering of all crustal phases, power function | f_py_CaSiO3 = 1.0
- 311-313       --> transient perturbations                                          | f_py_CaCO3  = 1.0    
- 314-316       --> transient perturbations                                          | f_py_CaSiO3 = 1.0  
- 411-413       --> 50% increase in volcanic CO2 outgassing                          | f_py_CaCO3  = 1.0
- 414-416       --> 50% increase in volcanic CO2 outgassing                          | f_py_CaSiO3 = 1.0
- 421-423       --> 50% increase in volcanic CO2 outgassing, +4 Tmol O2 equiv.       | f_py_CaCO3  = 1.0
- 424-426       --> 50% increase in volcanic CO2 outgassing, +4 Tmol O2 equiv.       | f_py_CaSiO3 = 1.0
#
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
