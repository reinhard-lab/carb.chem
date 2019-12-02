%% integration time (yrs):
   p.integration_time=10e6;

%% perturbation parameters:
 % form of forcing function:
   p.F_form = 'logistic'; % 'step'/'gaussian'/'logistic'/'periodic'
 % adjustment timescales (yrs):
   p.tau_pert = 2.0e6;    % gaussian/step
   p.tau_log  = 0.5e6;    % logistic
 % perturbation mid-points (yrs):
   p.t0     = 3.0e6;      % gaussian/step
   p.t0_log = 4.0e6;      % logistic
 % periodic forcing parameters (yrs);
   p.t_per_1 = 0.5;
   p.t_per_2 = 1.e6;
   p.t_per_3 = 1.0;
 % perturbation magnitudes:
   p.Delta_J_volc   = 0.0;
   p.Delta_J_MOR    = 0.0;
   p.Delta_J_w_org  = 0.0;
   p.Delta_J_w_py   = 0.0;
   p.Delta_J_w_ev   = 0.0;
   p.Delta_J_b_org  = 0.0;
   p.Delta_J_b_py   = 0.0;
   p.Delta_J_b_ev   = 0.0;
   p.Delta_J_red_v  = 0.0;
   
   p.Delta_f_w_sil  = 0.5;
   p.Delta_f_w_carb = 0.5;
   p.Delta_f_w_org  = 0.5;
   p.Delta_f_w_py   = 0.5;
   p.Delta_f_w_ev   = 0.5;
   
   p.Delta_py_CaCO3  = 0.0;
   p.Delta_py_CaSiO3 = 0.0;
   p.Delta_py_H      = 0.0;

%% dimensional/conversion parameters:
   p.ocean_mass_kg   = 1.41e21;       % mass of ocean [kg]
   p.atm_moles       = 1.773e20;      % molar mass of atmosphere [mol]
   p.tSY             = 3.1536*(1.e7); % seconds per year
   p.Salinity0       = 35;            % salinity [permil]
   p.K_sp_pressure   = 500;           % pressure [atm]
   
%% surface temperature model:
   p.t_Earth    = 4570;                 % age of the Earth [Ma]
   p.t_scenario = 750;                    % age of model scenario [positive]
   p.T0 = 288;                          % reference base temperature [K]
   p.S_relative = (1 + 0.38*(p.t_scenario/p.t_Earth))^(-1.0);   % relative solar luminosity
   p.S_0        = 1368;                 % modern solar constant [W/m^2]
   p.S_t        = p.S_0*p.S_relative;   % solar constant at model age [W/m^2]
   p.alpha_T    = 0.3;                  % planetary albedo
   p.emis_T     = 0.773;                % emissivity of atmosphere [see Mills et al., 2011]
   p.sig_T      = 5.67e-8;              % Boltzmann constant [W / m^2 * K]
   p.mu_T       = 0.5;                  % constant for derived climate sensitivity
     
   p.pCO2_0     = 278;        % baseline pCO2 [ppm]
   p.pCH4_0     = 0.715;      % baseline pCH4 (low) [ppm]
   p.pCH4_1     = 2.5;        % baseline pCH4 (high) [ppm]
   p.pN2O_0     = 0.290;      % baseline pN2O (low) [ppm]
   p.pN2O_1     = 2.5;        % baseline pN2O (high) [ppm]
  
%% deep ocean temperature:
   p.T_constant = 0;                 % set to '0' for constant deep ocean temperature
   p.T_eq0      = 4.0;               % reference equilibrium temperature [degC]
   p.t_rel      = 1500;              % relaxation timescale for deep ocean temperature [years]
   p.sig_deep   = 2;                 % climate 'sensitivity' of deep ocean [see Archer & Buffet, 2005]
  
%% surface weathering parameters:
 % carbonate/silicate cycle:
   p.f_w0_sil       = 1.0;        % baseline uplift/erodability factor, silicates
   p.f_w0_carb      = 1.0;        % baseline uplift/erodability factor, carbonates
   p.f_w0_ev        = 1.0;        % baseline uplift/erodability factor, evaporites
   p.ACT            = 0.09;       % activation energy for silicate weathering
   p.RUN            = 0.038;      % runoff dependence on temp, silicates
   p.n_CO2_sil      = 0.25;       % CO2 catalysis exponent [see GEOCARB II]
   p.n_CO2_carb     = 0.25;       % CO2 catalysis exponent [see GEOCARB II]
 
 % transport-limited weathering terms:
   p.W_max_sil    = 2.5;      % transport-limited silicate weathering rate [relative]
   p.W_max_carb   = 1.0;      % transport-limited carbonate weathering rate [relative]
   p.k_w_sil      = 100.0;    % kinetic constant for transport-limited weathering
   p.k_w_carb     = 100.0;    % kinetic constant for transport-limited weathering
 
 % factors partitioning H+ equivalents:
   p.f_py_CaCO3_0  = 1.0;       % fraction of pyrite H+ to carbonate weath
   p.f_py_CaSiO3_0 = 0.0;       % fraction of pyrite H+ to silicate weath
   p.f_py_H_0      = 0.0;       % fraction of pyrite H+ transferred to ocean
 
 % organic/pyrite cycle:
   pO2_0              = (0.2046)*p.atm_moles; % baseline pO2 [atm]
   p.M_O2_0           = pO2_0/1e12;      % baseline pO2 [mol]
   p.SO4_0            = 2000;           % baseline [SO4] [µM]
   p.M_SO4_0          = (p.SO4_0/1e18)*p.ocean_mass_kg;
   p.f_w0_org         = 1.0;             % fractional area of exposed orgC
   p.f_w0_py          = 1.0;             % fractional area of exposed pyS
   p.C_weath_exponent = 0.5;             % exponent for orgC weathering
   p.S_weath_exponent = 0.5;             % exponent for pyS weathering
   p.C_w_max          = 1.01;            % saturation param 1, orgC weathering
   p.C_w_sat          = 1.e-3;           % saturation param 2, orgC weathering
   p.S_w_max          = 1.0005;          % saturation param 1, pyS weathering
   p.S_w_sat          = 1.e-4;           % saturation param 2, pyS weathering
   p.S_b_max          = 1.0691;          % saturation param 1, pyS burial
   p.S_b_sat          = 1.8685;          % saturation param 2, pyS burial
   p.S_O              = 15/8;            % O2/FeS2 stoichiometry
   p.k_SO4            = 10;
   
%% some isotopic parameters:   
   p.d13_volc         = -6;
   p.d13_MOR          = -6;
   p.D_crust          = 1.0;
   p.D34_max          = 50;
   p.Km_D             = 10;

%% carbonate burial parameters:
   p.Omega_reference = 1.6;          % reference CaCO3 saturation state   
   p.n_carb          = 1.7;          % exponent for CaCO3 precipitation
   p.B_T0            = 400;          % total boron concentration [µmol/kg]
  
%% seafloor weathering parameters:
 % basic phys/chem constants:
   p.R        = 8.3144621*(1.0e-3);   % gas constant [kJ/K*mol]
   p.A_crust  = 2.63060374583339e+20; % surface area of weathering crust [m^2]
   p.gamma_H  = 0.65;                 % activity coefficient for H+
   p.gamma_OH = 0.18;                 % activity coefficient for OH-
  
 % mineral fractions (dimensionless):
   p.f_labradorite = 0.540;   
   p.f_diopside    = 0.320;
   p.f_glass       = 0.090;
   p.f_apatite     = 0.043;
   p.f_forsterite  = 0.013;
  
 % labradorite parameters:  
   p.kH_labradorite   = 10^(-8.28);   % rate constant (H+) [mol/m2*s]
   p.kOH_labradorite  = 10^(-9.07);   % rate constant (OH-) [mol/m2*s]
   p.kw_labradorite   = 10^(-11.5);   % rate constant (H2O) [mol/m2*s]
   p.EaH_labradorite  = 65;           % activation energy (H+) [kJ/mol]
   p.EaOH_labradorite = 50;           % activation energy (OH-) [kJ/mol]
   p.Eaw_labradorite  = 68;           % activation energy (H2O) [kJ/mol]
   p.nH_labradorite   = 0.7;          % reaction order (H+)
   p.nOH_labradorite  = 0.3;          % reaction order (OH-)
  
 % diopside parameters:  
   p.kH_diopside  = 10^(-9.85);       % rate constant (H+) [mol/m2*s]
   p.EaH_diopside = 42;               % activation energy (H+) [kJ/mol]
   p.nH_diopside  = 0.14;             % reaction order (H+)
 
 % basaltic glass parameters: 
   p.kH_glass   = 10^(-6.7);          % rate constant (H+) [mol/m2*s]
   p.kOH_glass  = 10^(-9.45);         % rate constant (OH-) [mol/m2*s] 
   p.kw_glass   = 10^(-10.35);        % rate constant (H2O) [mol/m2*s]
   p.EaH_glass  = 30;                 % activation energy (H+) [kJ/mol]
   p.EaOH_glass = 50;                 % activation energy (OH-) [kJ/mol]
   p.Eaw_glass  = 50;                 % activation energy (H2O) [kJ/mol]
   p.nH_glass   = 0.5;                % reaction order (H+)
   p.nOH_glass  = 0.175;              % reaction order (OH-)
 
 % apatite parameters: 
   p.kH_apatite  = 10^(-5.08);        % rate constant (H+) [mol/m2*s]
   p.EaH_apatite = 34.7;              % activation energy (H+) [kJ/mol]
   p.nH_apatite  = 0.87;              % reaction order (H+) 
 
 % forsterite parameters:
   p.kH_forsterite  = 10^(-6.6);      % rate constant (H+) [mol/m2*s]
   p.EaH_forsterite = 75;             % activation energy (H+) [kJ/mol]
   p.nH_forsterite  = 0.5;            % reaction order (H+)
   
%% initial conditions, state variables:
 % concentrations:
   p.pCO2_init = 500;                   % initial pCO2 [ppm]
   p.DIC_init  = 2100;                  % initial [DIC] [µmol/kg]
   p.Alk_init  = 2280;                  % initial [Alk] [µmol/kg]
   p.Ca_init   = 10280;                 % initial [Ca] [µmol/kg]
   p.SO4_init  = 20000;                  % initial [SO4] [µmol/kg]
   p.pO2_init  = (0.2046)*0.5;          % initial pO2 [atm]
  
 % masses:
   p.M_CO2_init = ((p.pCO2_init/(1e6))*p.atm_moles)/(1e12); % initial M_CO2 [Tmol]
   p.M_DIC_init = (p.DIC_init/1e18)*p.ocean_mass_kg;        % initial M_DIC [Tmol]
   p.M_Alk_init = (p.Alk_init/1e18)*p.ocean_mass_kg;        % initial M_Alk [Tmol]
   p.M_Ca_init  = (p.Ca_init/1e18)*p.ocean_mass_kg;         % initial M_Ca [Tmol]
   p.M_SO4_init = (p.SO4_init/1e18)*p.ocean_mass_kg;        % initial M_SO4 [Tmol]
   p.M_O2_init  = (p.pO2_init*p.atm_moles)/(1e12);          % initial M_O2 [Tmol] 
  
 % crustal masses:   
   p.M_org0_mol  = 1.0e21;                % initial M_org [mol]
   p.M_carb0_mol = 7.0e21;                % initial M_carb [mol]
   p.M_py0_mol   = 2.0e20;                % initial M_py [mol]
   p.M_ev0_mol   = 0.1e20;                % initial M_ev [mol]
   p.M_org0      = p.M_org0_mol/(1.e12);  % initial M_org [Tmol]
   p.M_carb0     = p.M_carb0_mol/(1.e12); % initial M_carb [Tmol]
   p.M_py0       = p.M_py0_mol/(1.e12);   % initial M_py [Tmol]
   p.M_ev0       = p.M_ev0_mol/(1.e12);   % initial M_ev [Tmol]
   
 % isotope compositions
   p.d13C_sw_init    = 1;
   p.d34S_sw_init    = 15;
   p.d13C_Morg_init  = -25;
   p.d13C_Mcc_init   = 1;
   p.d34S_Mpy_init   = -15;
   p.d34S_Mev_init   = 15;
   
%% initial conditions, fluxes:
 % input fluxes:  
   p.J_volc0   = 9.0;       % volcanic CO2 flux [TmolC/yr]
   p.J_MOR0    = 1.6;       % mid-ocean ridge CO2 flux [TmolC/yr]
   p.J_w_sil0  = 7.0;       % silicate weathering flux [TmolC/yr]
   p.J_w_carb0 = 10.0;      % carbonate weathering flux [TmolC/yr]
   p.J_w_org0  = 9.0;       % orgC weathering flux [TmolC/yr]
   p.J_w_py0   = 6.0;       % pyS weathering flux [TmolS/yr]
   p.J_w_ev0   = 1.0;       % evapS weathering flux [TmolS/yr]
   p.J_red_v0  = 0.0;
 
 % output fluxes:
   p.J_b_carb0  = 16.0;     % carbonate burial flux [TmolC/yr]
   
   p.k_cc       = p.J_b_carb0/((p.Omega_reference-1)^p.n_carb);  % rate constant for CaCO3 burial [y^-1]
   
   p.J_b_crust0 = 1.6;      % seafloor weathering flux [TmolC/yr]
   p.J_b_org0   = 9.0;      % orgC burial flux [TmolC/yr]
   p.J_b_py0    = 6.0;      % pyS burial flux [TmolS/yr]
   p.J_b_ev0    = 1.0;      % evapS burial flux [TmolS/yr]

%% scaling and offseting carbonate system:
   p.carb_scale=[1.0e-2,1.0,100,10,200,100 ...
           ,100,100];
   p.carb_offset=[0,0,300,10,2000,100 ...
           ,300,100];
   p.carb_label=[ ...
     'H    ';'OH   ';'CO2g ';'H2CO3';'HCO3 ';'CO3  ' ...
     ;'BOH3 ';'BOH4 '];
   p.carb_init_guess=[0.006, 3.6, 300, 10, 1800, 200, ...
                 300, 100];

   indx=0;
   indx=indx+1; p.index.T_deep=indx;
   indx=indx+1; p.index.M_DIC=indx;
   indx=indx+1; p.index.M_Alk=indx;
   indx=indx+1; p.index.M_Ca=indx;
   indx=indx+1; p.index.M_SO4=indx;
   indx=indx+1; p.index.M_O2=indx;
   indx=indx+1; p.index.M_org=indx;
   indx=indx+1; p.index.M_carb=indx;
   indx=indx+1; p.index.M_py=indx;
   indx=indx+1; p.index.M_ev=indx;
   indx=indx+1; p.index.d13C_sw=indx;
   indx=indx+1; p.index.d34S_sw=indx;

   p.offset(p.index.T_deep)=p.T_eq0/1.5;
   p.offset(p.index.M_DIC)=p.M_DIC_init/1.1;
   p.offset(p.index.M_Alk)=p.M_Alk_init/1.1;
   p.offset(p.index.M_Ca)=p.M_Ca_init/1.1;
   p.offset(p.index.M_SO4)=p.M_SO4_init/1.1;
   p.offset(p.index.M_O2)=p.M_O2_init/1.1;
   p.offset(p.index.M_org)=p.M_org0/1.1;
   p.offset(p.index.M_carb)=p.M_carb0/1.1;
   p.offset(p.index.M_py)=p.M_py0/1.1;
   p.offset(p.index.M_ev)=p.M_ev0/1.1;
   p.offset(p.index.d13C_sw)=p.d13C_sw_init/1.1;
   p.offset(p.index.d34S_sw)=p.d34S_sw_init/1.1;

   p.scale(p.index.T_deep)=1;
   p.scale(p.index.M_DIC)=1e6;
   p.scale(p.index.M_Alk)=1e6;
   p.scale(p.index.M_Ca)=1e7;
   p.scale(p.index.M_SO4)=1e7;
   p.scale(p.index.M_O2)=1e7;
   p.scale(p.index.M_org)=1e7;
   p.scale(p.index.M_carb)=1e7;
   p.scale(p.index.M_py)=1e7;
   p.scale(p.index.M_ev)=1e7;
   p.scale(p.index.d13C_sw)=1;
   p.scale(p.index.d34S_sw)=1;
  
   p.label(p.index.T_deep,:)     ='T_deep    ';
   p.label(p.index.M_DIC,:)      ='M_DIC     ';
   p.label(p.index.M_Alk,:)      ='M_Alk     ';
   p.label(p.index.M_Ca,:)       ='M_Ca      ';
   p.label(p.index.M_SO4,:)      ='M_SO4     ';
   p.label(p.index.M_O2,:)       ='M_O2      ';
   p.label(p.index.M_org,:)      ='M_org     ';
   p.label(p.index.M_carb,:)     ='M_carb    ';
   p.label(p.index.M_py,:)       ='M_py      ';
   p.label(p.index.M_ev,:)       ='M_ev      ';
   p.label(p.index.d13C_sw,:)    ='d13C_sw   ';
   p.label(p.index.d34S_sw,:)    ='d34S_sw   ';

 % fsolve options:
   p.optim_options=optimset('MaxIter',1000,'MaxFunEvals',1000 ...
                           ,'DerivativeCheck','off','Display','off' ...
                           ,'Jacobian','on' ...
                           );
 % ode solver options:
   p.ode_options=odeset('RelTol',1.0e-6,'AbsTol',1.e-6,'MaxStep',1e3);
  
%% logical switches:
   p.do_land_plants=0;
   p.CO2_forcing='goldblatt'; % 'mills' / 'goldblatt'
   p.do_silicate_weathering_feedback=1;
   p.do_carbonate_weathering_feedback=1;
   p.do_transport_limited_weathering=1; % transport-limited weathering
   p.do_carbonate_burial_feedback=1;
   p.do_seafloor_weathering=1;
   p.seafloor_weath='chem'; % 'chem' / 'temp'
   p.volc='forced'; % 'constant' / 'forced' / 'parameterized'
   p.organic_C_weath='saturation'; % 'power' / 'saturation' / 'forced'
   p.pyrite_S_weath='saturation'; % 'power' / 'saturation' / 'forced'
   p.evaporite_S_weath='parameterized'; % 'constant' / 'forced' / 'parameterized'
   p.organic_C_burial='parameterized'; % 'constant' / 'forced' / 'parameterized'
   p.pyrite_S_burial='parameterized'; % 'constant' / 'forced' / 'parameterized'
   p.evaporite_burial='parameterized'; % 'constant' / 'forced' / 'parameterized'
   p.do_constant_calcium=0;
   p.csys_solver='csys3';  % 'implicit' / 'gradient' / 'csys3'
   p.do_write_spinup=1;
   p.do_read_spinup=1;
   p.do_use_Ksp_aragonite=0;
   p.do_run_model=1;
