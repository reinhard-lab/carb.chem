%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function carb_chem(scenario,spin,time)
%% A simple model for exploring the links between the sedimentary
%% rock cycle, electron transfer, and carbonate chemistry
%% in the ocean
%%
%% Units: Tmol (10^12 mol), years, µmol/kg 
%% 
%% To run this program: specify scenario number and clear all first:
%% cd ~/Documents/MATLAB/carb_chem_mod/; clear all; carb_chem(...)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% close all:
fclose('all');
if nargin<1
  disp('*** must specify an argument with number of scenario to run.');
  return
end

%% declare output variable:
global o 

%% establish parameter set:
p=set_parameters(scenario);
% overwrite default integration time
p.integration_time  = time;
% overwrite forcing params if spinning up
if spin == 1
   p.Delta_J_volc   = 0.0;
   p.Delta_J_MOR    = 0.0;
   p.Delta_J_w_org  = 0.0;
   p.Delta_J_w_py   = 0.0;
   p.Delta_J_w_ev   = 0.0;
   p.Delta_J_b_org  = 0.0;
   p.Delta_J_b_py   = 0.0;
   p.Delta_J_b_ev   = 0.0;
   p.Delta_J_red_v  = 0.0;
   
   p.Delta_f_w_sil  = 0.0;
   p.Delta_f_w_carb = 0.0;
   p.Delta_f_w_org  = 0.0;
   p.Delta_f_w_py   = 0.0;
   p.Delta_f_w_ev   = 0.0;
   
   p.Delta_py_CaCO3  = 0.0;
   p.Delta_py_CaSiO3 = 0.0;
   p.Delta_py_H      = 0.0;
end

%% run model:
if p.do_run_model
  %% set (nondimensional) initial conditions from spinup or defaults:
  [initial_conditions,spinup]=set_initial_conditions(p);

  %% integrate: 
  [t,x] = ode15s(@rhs,[0:p.integration_time/1000:p.integration_time] ...
                 ,initial_conditions,p.ode_options,p,spinup);

  %% save to output file:
  save(sprintf('output/out_scenario%2.2d.mat',p.scenario),'p','o');
  
  %% dimensionalize output of ode solver:
  N=length(x(1,:)); %% number of state variables
  for i=1:N; x(:,i)=x(:,i)*p.scale(i)+p.offset(i);end

  %% write spinup:
  if p.do_write_spinup==1
    e=o.nt;
    spinup.J_volc=o.J_volc(e);
    spinup.J_w_sil=o.J_w_sil(e);
    spinup.J_w_sil_py=o.J_w_sil_py(e);
    spinup.J_w_carb=o.J_w_carb(e);
    spinup.J_w_carb_py=o.J_w_carb_py(e);
    spinup.J_H_py=o.J_H_py(e);
    spinup.J_MOR=o.J_MOR(e);
    spinup.J_w_org=o.J_w_org(e);
    spinup.J_b_carb=o.J_b_carb(e);
    spinup.J_b_crust=o.J_b_crust(e);
    spinup.J_b_org=o.J_b_org(e);
    spinup.J_w_py=o.J_w_py(e);
    spinup.J_w_ev=o.J_w_ev(e);
    spinup.J_b_py=o.J_b_py(e);
    spinup.J_b_ev=o.J_b_ev(e);
    spinup.J_red_v=o.J_red_v(e);
    spinup.M_org=o.M_org(e);
    spinup.M_carb=o.M_carb(e);
    spinup.M_py=o.M_py(e);
    spinup.M_ev=o.M_ev(e);
    spinup.M_DIC=o.M_DIC(e);
    spinup.M_Alk=o.M_Alk(e);
    spinup.M_Ca=o.M_Ca(e);
    spinup.M_SO4=o.M_SO4(e);
    spinup.M_O2=o.M_O2(e);
    spinup.M_ev=o.M_ev(e);
    spinup.Alk=o.Alk(e);
    spinup.BOH3=o.BOH3(e);
    spinup.BOH4=o.BOH4(e);
    spinup.B_T=o.B_T(e);
    spinup.CO2g=o.CO2g(e);
    spinup.CO3=o.CO3(e);
    spinup.CO3_sat=o.CO3_sat(e);
    spinup.Ohm=o.Ohm(e);
    spinup.DIC=o.DIC(e);
    spinup.Ca=o.Ca(e);
    spinup.SO4=o.SO4(e);
    spinup.pO2=o.pO2(e);
    spinup.pCH4=o.pCH4(e);
    spinup.pN2O=o.pN2O(e);
    spinup.f_org=o.f_org(e);
    spinup.f_py=o.f_py(e);
    spinup.H2CO3=o.H2CO3(e);
    spinup.HCO3=o.HCO3(e);
    spinup.H=o.H(e);
    spinup.OH=o.OH(e);
    spinup.S=o.S(e);
    spinup.T=o.T(e);
    spinup.pH=o.pH(e);
    spinup.T_surf=o.T_surf(e);
    spinup.T_deep=o.T_deep(e);
    spinup.d13_cc=o.d13_cc(e);
    spinup.D_org=o.D_org(e);
    spinup.D_py=o.D_py(e);
    spinup.d13C_sw=o.d13C_sw(e);
    spinup.d34S_sw=o.d34S_sw(e);
    last=length(t); 
    write_spinup(t(last),x(last,:),p,spin,spinup); 
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initial_conditions,spinup]=set_initial_conditions(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial conditions for state variables:
    initial_conditions=[p.T_eq0 ...
                        p.M_DIC_init p.M_Alk_init p.M_Ca_init ...
                        p.M_SO4_init p.M_O2_init ...
                        p.M_org0 p.M_carb0 ...
                        p.M_py0 p.M_ev0 ...
                        p.d13C_sw_init p.d34S_sw_init];

%% default values for these parameters in case spinup does not
%% exist and they are still needed:
    spinup.J_w_sil=p.J_w_sil0; 
    spinup.J_w_carb=p.J_w_carb0;
    spinup.J_b_carb=p.J_b_carb0; 
    spinup.J_b_crust=p.J_b_crust0;
    spinup.J_w_org=p.J_w_org0;
    spinup.J_b_org=p.J_b_org0;
    spinup.J_w_py=p.J_w_py0;
    spinup.J_b_py=p.J_b_py0;
    spinup.J_w_ev=p.J_w_ev0;
    spinup.J_b_ev=p.J_b_ev0;

%% read spinup data:
    if p.do_read_spinup==1 
      [t0,x0,spinup]=read_spinup(p);
      if isnan(t0)==0 
        t=t0; 
        initial_conditions=x0; 
        spinup=spinup;
      end
    end

%% nondimensionalize initial conditions:
    initial_conditions=(initial_conditions-p.offset)./p.scale;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = rhs(t,x,p,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solution of time-dependent equations for state variables

%% save carbonate system solution between calls:
persistent pH H OH CO2g H2CO3 HCO3 CO3 BOH3 BOH4
persistent nt % time step
%% declare global output variable:
global o

%% obtain variables from the state vector:
%% nondim=(dim-offset)/scale
%% dim=nondim*scale+offset
%% d(nondim)/dt=(d(dim)/dt)/scale
n=p.index.T_deep;
T_deep=     x(n)*p.scale(n)+p.offset(n);
n=p.index.M_DIC;
M_DIC=      x(n)*p.scale(n)+p.offset(n);
n=p.index.M_Alk;
M_Alk=      x(n)*p.scale(n)+p.offset(n);
n=p.index.M_Ca;
M_Ca=       x(n)*p.scale(n)+p.offset(n);
n=p.index.M_SO4;
M_SO4=      x(n)*p.scale(n)+p.offset(n);
n=p.index.M_O2;
M_O2=       x(n)*p.scale(n)+p.offset(n);
n=p.index.M_org;
M_org=      x(n)*p.scale(n)+p.offset(n);
n=p.index.M_carb;
M_carb=     x(n)*p.scale(n)+p.offset(n);
n=p.index.M_py;
M_py=       x(n)*p.scale(n)+p.offset(n);
n=p.index.M_ev;
M_ev=       x(n)*p.scale(n)+p.offset(n);
n=p.index.d13C_sw;
d13C_sw=    x(n)*p.scale(n)+p.offset(n);
n=p.index.d34S_sw;
d34S_sw=    x(n)*p.scale(n)+p.offset(n);

%% other physical parameters for solubility constants:
T = T_deep;
S = p.Salinity0;
P = p.K_sp_pressure;

%% calculate carbonate saturation using Ca: 
Ca      = (M_Ca*(1e18))/p.ocean_mass_kg;
CO3_sat = K_sp(p,P,T,S)/Ca;
  
%% solve for concentrations from reservoir masses:
DIC  = (M_DIC*1.e18)/p.ocean_mass_kg;
Alk  = (M_Alk*1.e18)/p.ocean_mass_kg;
B_T  = p.B_T0;
SO4  = (M_SO4*1.e18)/p.ocean_mass_kg;
pO2  = (M_O2*1.e12)/p.atm_moles;
pCH4 = p.pCH4_0;
pN2O = p.pN2O_0;

%% solve for initial carbonate system:
if isempty(H)
  init_guess=p.carb_init_guess; 
else
  init_guess=[H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4]; 
end
if     strcmp(p.csys_solver,'implicit')==1
  [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
    =csys_iterative(t,p,Alk,DIC,B_T,T,S,init_guess);
elseif strcmp(p.csys_solver,'gradient')==1
  [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
    =csys_gradient(t,p,Alk,DIC,B_T,T,S,init_guess);
elseif strcmp(p.csys_solver,'approx')==1
  [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
    =csys_approx(t,p,Alk,DIC,B_T,T,S,init_guess);
elseif strcmp(p.csys_solver,'csys3')==1
  [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
    =csys_csys3(t,p,Alk,DIC,B_T,T,S,init_guess);
else
  disp('*** unrecognized carbonate solver, pausing.')
  pause
end

%% check for unphysical values for (mass) state variables:
if min([M_DIC,M_Alk,M_Ca,M_SO4,M_O2])<0
  fprintf(1,'*** negative state variable!  pausing.  t=%gMyr\n',t/1.e6)
  for i=1:length(x)
    fprintf('%s=%g;\n',p.label(i,:),x(i)*p.scale(i)+p.offset(i)); 
  end
  pause
end

%% check for unphysical values for carbonate system:
if min([CO2g,H2CO3,OH,H,pH,HCO3,CO3,BOH3,BOH4])<0
  fprintf(1,['*** negative carbonate variable!    t=%gMyr,\n' ...
             ' Alk=%g,\n DIC=%g,\n B_T=%g,\n ' ...
             'CO2g=%g,\n H2CO3=%g,\n OH=%g,\n H=%g,\n pH=%g,\n ' ...
             'HCO3=%g,\n CO3=%g,\n BOH3=%g,\n BOH4=%g\n pausing.\n'] ...
          ,t/1.e6,Alk,DIC,B_T,CO2g,H2CO3,OH,H,pH, ...
          HCO3,CO3,BOH3,BOH4);  
  pause
end

%% calculate surface temperature:
T_surf=T_surf(t,p,CO2g,pCH4,pN2O);

%% calculate atmospheric fraction of CO2:
C_atm   = (CO2g/1.e6)*p.atm_moles;
F_atm   = C_atm/(M_DIC*(1.e12));

%% calculate isotope fractionation factors:
d13_cc=d13_cc(t,p,T_surf,F_atm,d13C_sw);
D_org=D_org(t,p,CO2g);
D_py=D_py(t,p,SO4);

%% obtain flux values for a given timestep: 
T_eq=T_eq(t,p,CO2g);
J_volc=J_volc(t,p);
J_w_sil=J_w_sil(t,p,T_surf,CO2g,spinup);
J_w_sil_py=f_py_CaSiO3(t,p)*J_w_py(t,p,M_O2,M_py,spinup);
J_w_carb=J_w_carb(t,p,T_surf,CO2g,M_carb,spinup);
J_w_carb_py=f_py_CaCO3(t,p)*J_w_py(t,p,M_O2,M_py,spinup);
J_H_py=f_py_H(t,p)*J_w_py(t,p,M_O2,M_py,spinup);
J_MOR=J_MOR(t,p);
J_w_org=J_w_org(t,p,M_O2,M_org,spinup);
J_b_carb=J_b_carb(t,p,CO3,CO3_sat,spinup);
J_b_crust=J_b_crust(t,p,H,OH,T_deep,spinup);
J_b_org=J_b_org(t,p,M_O2,spinup);
J_w_py=J_w_py(t,p,M_O2,M_py,spinup);
J_w_ev=J_w_ev(t,p,M_ev,spinup);
J_b_py=J_b_py(t,p,J_b_org,SO4,spinup);
J_b_ev=J_b_ev(t,p,Ca,SO4,spinup);
J_red_v=J_red_v(t,p);
f_py_CaCO3=f_py_CaCO3(t,p);
f_py_CaSiO3=f_py_CaSiO3(t,p);
f_py_H=f_py_H(t,p);

%% trigger for constant Ca:
if p.do_constant_calcium==1
  FACTOR_Ca=0;
else
  FACTOR_Ca=1;
end

%% initialize array:
Nprognostic=length(p.index);
answer=zeros(Nprognostic,1)*NaN;

%% ------------------------------------------------------------------------
%% Specify all equations for state variables:
%% ------------------------------------------------------------------------

%% T_deep:
answer(p.index.T_deep) = p.T_constant ...
                        *( (T_eq - T_deep) / p.t_rel );
%% M_DIC:
answer(p.index.M_DIC) = ( ( J_volc + J_w_carb ...
                           + (2*f_py_CaCO3*J_w_py) ...
                           + J_MOR + J_w_org ) ...
                         -( J_b_carb + J_b_org + J_b_crust ) );
%% M_Alk:
answer(p.index.M_Alk) = ( ( 2*J_w_sil + 2*J_w_carb ...
                           +(2*f_py_CaCO3*J_w_py) ...
                           +2*J_b_py ) ...
                         -( 2*J_b_carb ...
                           +(2*f_py_H*J_w_py) ) );
%% M_Ca:
answer(p.index.M_Ca) = FACTOR_Ca ...
                      * ( J_w_sil + J_w_carb ...
                         +(f_py_CaSiO3*J_w_py) ...
                         +(2*f_py_CaCO3*J_w_py) ...
                         +J_w_ev ) ...
                       -( J_b_carb+J_b_ev );
%% M_SO4:
answer(p.index.M_SO4) = ( ( J_w_py + J_w_ev ) ...
                         -( J_b_py + J_b_ev ) );
%% M_O2:
answer(p.index.M_O2) = ( ( J_b_org ...
                          + p.S_O*J_b_py ) ...
                        -( J_w_org ...
                          + p.S_O*J_w_py ...
                          + J_red_v ) );
%% M_org:
answer(p.index.M_org) = ( J_b_org - J_w_org );
%% M_carb:
answer(p.index.M_carb) = ( J_b_carb - J_w_carb );
%% M_py:
answer(p.index.M_py) = ( J_b_py - J_w_py );
%% M_ev:
answer(p.index.M_ev) = ( J_b_ev - J_w_ev );
                  
%% d13C_sw:
answer(p.index.d13C_sw) = ( ( (J_volc*(p.d13_volc-d13C_sw)) ...
                             +(J_w_carb*(p.d13C_Mcc_init-d13C_sw)) ...
                             +((2*f_py_CaCO3*J_w_py)*(p.d13C_Mcc_init-d13C_sw)) ...
                             +(J_MOR*(p.d13_MOR-d13C_sw)) ...
                             +(J_w_org*(p.d13C_Morg_init-d13C_sw)) ) ...
                           -( (J_b_carb*(d13C_sw-d13_cc)) ...
                             +(J_b_org*(d13C_sw-D_org)) ...
                             +(J_b_org*(d13C_sw-p.D_crust)) ) ) * (1/M_DIC);
                         
%% d34S_sw:
answer(p.index.d34S_sw) = ( ( (J_w_py*(p.d34S_Mpy_init-d34S_sw)) ...
                             +(J_w_ev*(p.d34S_Mev_init-d34S_sw)) ) ...
                           -( (J_b_py*(d34S_sw-D_py)) ) ) * (1/M_SO4);

%% rescale calculated state variables:
%% nondim=(dim-offset)/scale
%% dim=nondim*scale+offset
%% d(nondim)/dt=(d(dim)/dt)/scale
answer=answer'./p.scale';

%% save diagnostics to output variable:
% first, advance time step used for saving diagnostics.
if isempty(nt); nt=1; else; nt=nt+1; end
if nt==1; init=1; else; init=0; end
o.max_nt=1000000;
if nt>o.max_nt
  fprintf('*** nt>o.max_nt, nt=%d; o.max_nt=%d; t=%g\n',nt,o.max_nt,t);
  pause
end
%% avoid saving trial steps from ode solver:
while nt>2 && t<=o.t(nt-1)
  nt=nt-1;
end
o.nt=nt;

if init; ZZ=zeros(o.max_nt,1)*NaN; end
if init; o.Alk=ZZ; end;             		o.Alk(nt)=Alk;
if init; o.BOH3=ZZ; end;            		o.BOH3(nt)=BOH3;
if init; o.BOH4=ZZ; end;            		o.BOH4(nt)=BOH4;
if init; o.B_T=ZZ; end;             		o.B_T(nt)=B_T;
if init; o.CO2g=ZZ; end;            		o.CO2g(nt)=CO2g;
if init; o.CO3=ZZ; end;             		o.CO3(nt)=CO3;
if init; o.CO3_sat=ZZ; end;       	    	o.CO3_sat(nt)=CO3_sat;
if init; o.Ohm=ZZ; end;                     o.Ohm(nt)=CO3/CO3_sat;
if init; o.DIC=ZZ; end;             		o.DIC(nt)=DIC;
if init; o.Ca=ZZ; end;              		o.Ca(nt)=Ca;
if init; o.SO4=ZZ; end;                     o.SO4(nt)=SO4;
if init; o.J_volc=ZZ; end;                  o.J_volc(nt)=J_volc;
if init; o.J_b_carb=ZZ; end; 		        o.J_b_carb(nt)=J_b_carb;
if init; o.J_b_crust=ZZ; end;               o.J_b_crust(nt)=J_b_crust;
if init; o.pO2=ZZ; end;                     o.pO2(nt)=pO2;
if init; o.pCH4=ZZ; end;                    o.pCH4(nt)=pCH4;
if init; o.pN2O=ZZ; end;                    o.pN2O(nt)=pN2O;
if init; o.J_b_org=ZZ; end;         		o.J_b_org(nt)=J_b_org;
if init; o.J_w_sil=ZZ; end;                 o.J_w_sil(nt)=J_w_sil;
if init; o.J_w_sil_py=ZZ; end;              o.J_w_sil_py(nt)=J_w_sil_py;
if init; o.f_py_CaSiO3=ZZ; end;             o.f_py_CaSiO3(nt)=f_py_CaSiO3;
if init; o.J_w_carb=ZZ; end;                o.J_w_carb(nt)=J_w_carb;
if init; o.J_w_carb_py=ZZ; end;             o.J_w_carb_py(nt)=J_w_carb_py;
if init; o.f_py_CaCO3=ZZ; end;              o.f_py_CaCO3(nt)=f_py_CaCO3;
if init; o.J_H_py=ZZ; end;                  o.J_H_py(nt)=J_H_py;
if init; o.f_py_H=ZZ; end;                  o.f_py_H(nt)=f_py_H;
if init; o.J_MOR=ZZ; end;                   o.J_MOR(nt)=J_MOR;
if init; o.J_w_org=ZZ; end;                 o.J_w_org(nt)=J_w_org;
if init; o.J_w_py=ZZ; end;                  o.J_w_py(nt)=J_w_py;
if init; o.J_w_ev=ZZ; end;                  o.J_w_ev(nt)=J_w_ev;
if init; o.J_b_py=ZZ; end;                  o.J_b_py(nt)=J_b_py;
if init; o.J_b_ev=ZZ; end;                  o.J_b_ev(nt)=J_b_ev;
if init; o.J_red_v=ZZ; end;                 o.J_red_v(nt)=J_red_v;
if init; o.f_org=ZZ; end;                   o.f_org(nt)=J_b_org/(J_b_org+J_b_carb+J_b_crust);
if init; o.f_py=ZZ; end;                    o.f_py(nt)=J_b_py/(J_b_py+J_b_ev);
if init; o.H2CO3=ZZ; end;           		o.H2CO3(nt)=H2CO3;
if init; o.HCO3=ZZ; end;             		o.HCO3(nt)=HCO3;
if init; o.H=ZZ; end;               		o.H(nt)=H;
if init; o.M_DIC=ZZ; end;           		o.M_DIC(nt)=M_DIC;
if init; o.M_Alk=ZZ; end;                   o.M_Alk(nt)=M_Alk;
if init; o.M_Ca=ZZ; end;                    o.M_Ca(nt)=M_Ca;
if init; o.M_SO4=ZZ; end;                   o.M_SO4(nt)=M_SO4;
if init; o.M_O2=ZZ; end;                    o.M_O2(nt)=M_O2;
if init; o.M_org=ZZ; end;                   o.M_org(nt)=M_org;
if init; o.M_carb=ZZ; end;                  o.M_carb(nt)=M_carb;
if init; o.M_py=ZZ; end;                    o.M_py(nt)=M_py;
if init; o.M_ev=ZZ; end;                    o.M_ev(nt)=M_ev;
if init; o.OH=ZZ; end;              		o.OH(nt)=OH;
if init; o.S=ZZ; end;               		o.S(nt)=S;
if init; o.T=ZZ; end;               		o.T(nt)=T;
if init; o.pH=ZZ; end;              		o.pH(nt)=pH;
if init; o.T_surf=ZZ; end;                  o.T_surf(nt)=T_surf;
if init; o.T_deep=ZZ; end;                  o.T_deep(nt)=T_deep;
if init; o.F_atm=ZZ; end;                   o.F_atm(nt)=F_atm;
if init; o.d13_cc=ZZ; end;                  o.d13_cc(nt)=d13_cc;
if init; o.D_org=ZZ; end;                   o.D_org(nt)=D_org;
if init; o.D_py=ZZ; end;                    o.D_py(nt)=D_py;
if init; o.d13C_sw=ZZ; end;                 o.d13C_sw(nt)=d13C_sw;
if init; o.d34S_sw=ZZ; end;                 o.d34S_sw(nt)=d34S_sw;
if init; o.t=ZZ; end;               		o.t(nt)=t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% forcing function -------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=F(t,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specifies shape of forcing function:

  if strcmp(p.F_form,'gaussian')
      out = exp(-(t-p.t0)^2/p.tau_pert^2);
  elseif strcmp(p.F_form,'step')
      out = (sqrt(pi)/2)* ...
            0.5*(1+tanh((t-(p.t0-p.tau_pert))/0.01e6)) .*  ...
            0.5*(1+tanh(((p.t0+p.tau_pert)-t)/0.01e6));
  elseif strcmp(p.F_form,'logistic')
      out = ((1+exp(-(t-(p.t0_log-p.tau_log))/(0.5*p.tau_log)))^(-1.0)) ...
           *((1+exp(-(t-(p.t0_log+p.tau_log))/(0.5*p.tau_log)))^(-1.0));
  elseif strcmp(p.F_form,'periodic')
      out = ( 1+(p.t_per_1*sin((pi())*t/p.t_per_2)) ) / p.t_per_3;
  else
    disp('*** Error: unknown form for F(t), ^C to quit.')
    pause
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proton flux functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_py_CaCO3(t,p)
out = p.f_py_CaCO3_0*(1+p.Delta_py_CaCO3*F(t,p));

function out=f_py_CaSiO3(t,p)
out = p.f_py_CaSiO3_0+(p.Delta_py_CaSiO3*F(t,p));

function out=f_py_H(t,p)
out = p.f_py_H_0*(1+p.Delta_py_H*F(t,p));
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temperature functions --------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=T_eq(t,p,CO2g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% equilibrium deep ocean temperature at given pCO2

    out = p.T_eq0 ...
         + ( p.sig_deep*log(CO2g/p.pCO2_0)/log(2) );
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=T_surf(t,p,CO2g,pCH4,pN2O)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculates Earth surface temperature based on luminosity, pCO2, and pCH4

  % base (effective) temperature based on luminosity and albedo:
        T_base  = ( (p.S_t*(1-p.alpha_T)) ...
                 /((4*p.sig_T)*(1-(0.5*p.emis_T))) )^(0.25);
 
  % radiative forcing due to CO2 (see Byrne & Goldblatt, 2014):
        pCO2    = CO2g/1.e6;
        pCO2_0  = p.pCO2_0/1.e6;
        lam_CO2 = 5.32*log(pCO2/pCO2_0) + 0.39*log(pCO2/pCO2_0)^2.0;
  
  % radiative forcing due to CH4 (see Byrne & Goldblatt, 2014):     
    if pCH4 <= 2.5
        pCH4    = pCH4/1.e6;
        pCH4_0  = p.pCH4_0/1.e6;
        lam_CH4 = 1173*(sqrt(pCH4)-sqrt(pCH4_0)) ...
                 -71636*(sqrt(pCH4)-sqrt(pCH4_0))^2;
    elseif pCH4 > 2.5
        pCH4    = pCH4/1.e6;
        pCH4_1  = p.pCH4_1/1.e6;
        lam_CH4 = 0.824 + (4/5)*log(pCH4/pCH4_1) + (1/5)*log(pCH4/pCH4_1)^2;   
    end
    
  % radiative forcing due to N2O (see Byrne & Goldblatt, 2014):
    if pN2O <= 2.5
        pN2O    = pN2O/1.e6;
        pN2O_0  = p.pN2O_0/1.e6;
        lam_N2O = 3899*(sqrt(pN2O)-sqrt(pN2O_0)) ...
                 +38256*(sqrt(pN2O)-sqrt(pN2O_0))^2;
    elseif pN2O > 2.5
        pN2O    = pCH4/1.e6;
        pN2O_1  = p.pN2O_1/1.e6;
        lam_N2O = 4.182 + 3*log(pN2O/pN2O_1) + 0.5469*log(pN2O/pN2O_1)^2;   
    end
        del_T   = p.mu_T * (lam_CO2+lam_CH4+lam_N2O);
        out     = T_base + del_T;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions specifying input fluxes --------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_volc(t,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% volcanic CO2 flux:

  if strcmp(p.volc,'constant')
      out = p.J_volc0;
  elseif strcmp(p.volc,'forced') && strcmp(p.F_form,'gaussian')
      out = p.J_volc0*(1+p.Delta_J_volc*F(t,p));
  elseif strcmp(p.volc,'forced') && strcmp(p.F_form,'step')
      out = p.J_volc0*(1+p.Delta_J_volc*F(t,p));
  elseif strcmp(p.volc,'forced') && strcmp(p.F_form,'logistic')
      out = p.J_volc0*(1+p.Delta_J_volc*F(t,p));
  elseif strcmp(p.volc,'forced') && strcmp(p.F_form,'periodic')
      out = p.J_volc0*(p.Delta_J_volc*F(t,p));
  else
      disp('*** Error: unknown form for volcanic CO2 flux, ^C to quit.')
      pause
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_red_v(t,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% volcanic reductant flux:

  out = p.J_red_v0+(p.Delta_J_red_v*F(t,p));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_MOR(t,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mid-ocean ridge CO2 flux:

  out = p.J_MOR0*(1+p.Delta_J_MOR*F(t,p));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_w_sil(t,p,T_surf,CO2g,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% silicate weathering flux:

   if p.do_silicate_weathering_feedback==1
       if p.do_transport_limited_weathering==1
        % non-dimensional uplift/erodability factor  
           f_w_sil  = p.f_w0_sil*(1+p.Delta_f_w_sil*F(t,p));
        % temperature function
           f_T_sil   = ( exp(p.ACT*(T_surf-p.T0)) ) ...
                      *( (1+p.RUN*(T_surf-p.T0))^(0.65) );
        % CO2 function
           f_CO2_sil = (CO2g/p.pCO2_0)^p.n_CO2_sil;
        % combined kinetic silicate weathering expression
           W_K_sil   = f_w_sil*f_T_sil*f_CO2_sil;
        % total (transport+kinetic) silicate weathering
           W_tot_sil = (p.W_max_sil-W_K_sil) ...
                        *(1+exp(-p.k_w_sil*(W_K_sil-p.W_max_sil)))^(-1.0) ...
                        +W_K_sil;
           out       = W_tot_sil*p.J_w_sil0;
       else
        % non-dimensional uplift/erodability factor  
           f_w_sil   = p.f_w0_sil*(1+p.Delta_f_w_sil*F(t,p));
        % temperature function
           f_T_sil   = ( exp(p.ACT*(T_surf-p.T0)) ) ...
                      *( (1+p.RUN*(T_surf-p.T0))^(0.65) );
        % CO2 function
           f_CO2_sil = (CO2g/p.pCO2_0)^p.n_CO2_sil;
        % combined kinetic silicate weathering expression
           out   = f_w_sil*f_T_sil*f_CO2_sil*p.J_w_sil0;
       end
   else
       out       = spinup.J_w_sil;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_w_carb(t,p,T_surf,CO2g,M_carb,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% carbonate weathering flux:

   if p.do_carbonate_weathering_feedback==1
       if p.do_transport_limited_weathering==1
        % non-dimensional uplift/erodability factor
           f_w_carb   = p.f_w0_carb*(1+p.Delta_f_w_carb*F(t,p));
        % temperature function
           f_T_carb   = 1+(0.087*(T_surf-p.T0));
        % CO2 function
           f_CO2_carb = (CO2g/p.pCO2_0)^p.n_CO2_carb;
        % sedimentary reservoir function
           f_sedM     = (M_carb/p.M_carb0);
        % combined kinetic carbonate weathering expression
           W_K_carb   = f_w_carb*f_T_carb*f_CO2_carb*f_sedM;
        % total (transport+kinetic) carbonate weathering
           W_tot_carb = (p.W_max_carb-W_K_carb) ...
                       *(1+exp(-p.k_w_carb*(W_K_carb-p.W_max_carb)))^(-1.0) ...
                       +W_K_carb;
           out        = W_tot_carb*p.J_w_carb0;
       else
        % non-dimensional uplift/erodability factor
           f_w_carb   = p.f_w0_carb*(1+p.Delta_f_w_carb*F(t,p));
        % temperature function
           f_T_carb   = 1+(0.087*(T_surf-p.T0));
        % CO2 function
           f_CO2_carb = (CO2g/p.pCO2_0)^p.n_CO2_carb;
        % combined kinetic carbonate weathering expression
           out        = f_w_carb*f_T_carb*f_CO2_carb*p.J_w_carb0;
       end
   else
       out            = spinup.J_w_carb;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_w_org(t,p,M_O2,M_org,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organic carbon weathering flux:

  if strcmp(p.organic_C_weath,'power')
      f_w_org = p.f_w0_org*(1+p.Delta_f_w_org*F(t,p));
      out     = f_w_org ...
               *p.J_w_org0 ...
               *(M_org/p.M_org0) ...
               *((M_O2/p.M_O2_0)^p.C_weath_exponent);
  elseif strcmp(p.organic_C_weath,'saturation')
      pO2     = (M_O2*1e12)/p.atm_moles;
      f_w_org = p.f_w0_org*(1+p.Delta_f_w_org*F(t,p));
      r_C     = p.C_w_max*(pO2/(pO2+p.C_w_sat));
      out     = f_w_org ...
               *r_C ...
               *p.J_w_org0 ...
               *(M_org/p.M_org0);
  elseif strcmp(p.organic_C_weath,'forced')
      out=p.J_w_org0*(1+p.Delta_J_w_org*F(t,p));
  elseif strcmp(p.organic_C_weath,'constant')
      out=spinup.J_w_org;
  else
      disp('*** Error: unknown form for oxidative weathering feedback, ^C to quit.')
      pause
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_w_py(t,p,M_O2,M_py,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pyrite sulfur weathering flux:

  if strcmp(p.pyrite_S_weath,'power')
      f_w_py = p.f_w0_py*(1+p.Delta_f_w_py*F(t,p));
      out    = f_w_py ...
              *p.J_w_py0 ...
              *(M_py/p.M_py0) ...
              *((M_O2/p.M_O2_0)^p.S_weath_exponent);
  elseif strcmp(p.pyrite_S_weath,'saturation')
      pO2    = (M_O2*1e12)/p.atm_moles;
      f_w_py = p.f_w0_py*(1+p.Delta_f_w_py*F(t,p));
      r_S    = p.S_w_max*(pO2/(pO2+p.S_w_sat));
      out    = f_w_py ...
              *r_S ...
              *p.J_w_py0 ...
              *(M_py/p.M_py0);
  elseif strcmp(p.pyrite_S_weath,'forced')
      out=p.J_w_py0*(1+p.Delta_J_w_py*F(t,p));
  elseif strcmp(p.pyrite_S_weath,'constant')
      out=spinup.J_w_py;
  else
      disp('*** Error: unknown form for oxidative weathering feedback, ^C to quit.')
      pause
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_w_ev(t,p,M_ev,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% evaporite sulfur weathering flux:  

  if strcmp(p.evaporite_S_weath,'constant')
      out = spinup.J_w_ev;
  elseif strcmp(p.evaporite_S_weath,'forced')
      out=p.J_w_ev0*(1+p.Delta_J_b_ev*F(t,p));
  elseif strcmp(p.evaporite_S_weath,'parameterized')
      f_w_ev = p.f_w0_ev*(1+p.Delta_f_w_ev*F(t,p));
      out = f_w_ev ...
           *p.J_w_ev0 ...
           *(M_ev/p.M_ev0);
  else
      disp('*** Error: unknown form for evaporite weathering flux, ^C to quit.')
      pause
  end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function specifying output fluxes --------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_b_carb(t,p,CO3,CO3_sat,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% carbonate burial flux:

  if p.do_carbonate_burial_feedback==1
     Omega = CO3/CO3_sat;
     dO    = (Omega-1.0);
     out   = p.k_cc*sign(dO)*(abs(dO)^p.n_carb);
  else
     Omega = CO3/CO3_sat;
     dO    = (Omega-1.0);
     if Omega >= 1.0
       out = p.k_cc*(dO^p.n_carb);
     else
       out = 0;
     end
  end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_b_crust(t,p,H,OH,T_deep,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% seafloor weathering flux:

  if p.do_seafloor_weathering==1 && strcmp(p.seafloor_weath,'chem')
    T_perc          = 313;
    H_activity      = (H/(1.0e6))*p.gamma_H;
    OH_activity     = (OH/(1.0e6))*p.gamma_OH;
    R_labradorite   = ( p.kH_labradorite*exp(-p.EaH_labradorite/(p.R*T_perc)) ...
                       *(H_activity^p.nH_labradorite) ) ...
                     +( p.kOH_labradorite*exp(-p.EaOH_labradorite/(p.R*T_perc)) ...
                       *(OH_activity^p.nOH_labradorite) ) ...
                     +( p.kw_labradorite*exp(-p.Eaw_labradorite/(p.R*T_perc)) );
    R_diopside      = ( p.kH_diopside*exp(-p.EaH_diopside/(p.R*T_perc)) ...
                       *(H_activity^p.nH_diopside) );
    R_glass         = ( p.kH_glass*exp(-p.EaH_glass/(p.R*T_perc)) ...
                       *(H_activity^p.nH_glass) ) ...
                     +( p.kOH_glass*exp(-p.EaOH_glass/(p.R*T_perc)) ...
                       *(OH_activity^p.nOH_glass) ) ...
                     +( p.kw_glass*exp(-p.Eaw_glass/(p.R*T_perc)) );
    R_apatite       = ( p.kH_apatite*exp(-p.EaH_apatite/(p.R*T_perc)) ...
                       *(H_activity^p.nH_apatite) );
    R_forsterite    = ( p.kH_forsterite*exp(-p.EaH_forsterite/(p.R*T_perc)) ...
                       *(H_activity^p.nH_forsterite) );
    R_tot           = (p.f_labradorite*R_labradorite) ...
                     +(p.f_diopside*R_diopside) ...
                     +(p.f_glass*R_glass) ...
                     +(p.f_apatite*R_apatite) ...
                     +(p.f_forsterite*R_forsterite);
    R_tot_converted = R_tot*p.tSY;
    out             = (p.A_crust*R_tot_converted)/(1.0e12);
  elseif p.do_seafloor_weathering==1 && strcmp(p.seafloor_weath,'temp')
    f_sr = 1.0;
    out  = p.J_b_crust0*f_sr*((T_deep-p.T_eq0)/15.12);
  else
    out  = p.J_b_crust0;
  end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_b_org(t,p,M_O2,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organic carbon burial flux:
   
   if strcmp(p.organic_C_burial,'constant')
     out = spinup.J_b_org;
   elseif strcmp(p.organic_C_burial,'forced')
     out = p.J_b_org0*(1+p.Delta_J_b_org*F(t,p));
   elseif strcmp(p.organic_C_burial,'parameterized')
     s1 = 1.0;
     s2 = 7;
     s3 = 0.05;
     RO2 = M_O2/p.M_O2_0;
     out = p.J_b_org0 ...
          *(1+p.Delta_J_b_org*F(t,p)) ...
          *(s1+(s2*exp(-RO2/s3)));
   else
     disp('*** Error: unknown form for C_org burial flux, ^C to quit.')
     pause
   end
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_b_py(t,p,J_b_org,SO4,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pyrite sulfur burial flux:

   if strcmp(p.pyrite_S_burial,'constant')
      out = spinup.J_b_py;
   elseif strcmp(p.pyrite_S_burial,'forced')  
      out = p.J_b_py0*(1+p.Delta_J_b_py*F(t,p));
   elseif strcmp(p.pyrite_S_burial,'parameterized')
      out = p.J_b_py0 ...
           *(1+p.Delta_J_b_py*F(t,p)) ...
           *(J_b_org/p.J_b_org0) ...
           *(SO4/p.SO4_init);
   else
      disp('*** Error: unknown form for FeS2 burial flux, ^C to quit.')
      pause
   end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=J_b_ev(t,p,Ca,SO4,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% evaporite sulfur burial flux:  

  if strcmp(p.evaporite_burial,'constant')
      out = spinup.J_b_ev;
  elseif strcmp(p.evaporite_burial,'forced')
      out = p.J_b_ev0*(1+p.Delta_J_b_ev*F(t,p));
  elseif strcmp(p.evaporite_burial,'parameterized')
      IAP_0 = p.Ca_init * p.SO4_init;
      IAP   = Ca * SO4;
      out   = p.J_b_ev0 * ( IAP/IAP_0 );
  else
      disp('*** Error: unknown form for evaporite burial flux, ^C to quit.')
      pause
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions specifying isotope fractionations ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=d13_cc(t,p,T_surf,F_atm,d13C_sw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fractionation between surface carbon reservoir and CaCO3

    d_ocean = d13C_sw + ( F_atm*((9483/T_surf)-23.89) );
    out = d_ocean + 15.10 - (4232/T_surf);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=D_org(t,p,CO2g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% photosynthetic C-isotope fractionation

    out = abs( (80/(0.034*CO2g)) - 33 );    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=D_py(t,p,SO4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sulfate reduction S-isotope fractionation

    out = p.D34_max * (SO4/(p.Km_D+SO4));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions for calculating carbonate system parameters ------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K0=K_0(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solubility of CO2 in seawater
% Note: units are mol/kg*atm (same as µmol/kg*ppm)

  T  = T_deep+273.15;
  K0 = exp( ...
      (9345.17/T)-60.2409+23.3585*log(T/100.) ...
      +S*( 0.023517-0.00023656*T  +  ...
           0.00000047036*T*T ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K1=K_1(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first dissociation constant of carbonic acid
% Note: calculated in mol/kg, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  tmp1    = 2.83655 - 2307.1266 ./ T - 1.5529413 .* log(T);
  tmp2    =         - (0.20760841 + 4.0484 ./ T) .* sqrt(S);
  tmp3    =         + 0.08468345 .* S - 0.00654208 .* S .* sqrt(S);   
  tmp4    =         + log(1 - 0.001005 .* S);
  lnK1roy = tmp1 + tmp2 + tmp3 + tmp4;
  K1      = exp(lnK1roy);

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end;

  K1      = K1*exp(p.lnkpok0(1));

%% convert units (µmol/kg):
  K1      = K1*1.e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K2=K_2(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% second dissociation constant of carbonic acid
% Note: calculated in mol/kg, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  tmp1    = -9.226508 - 3351.6106 ./ T - 0.2005743 .* log(T);
  tmp2    = (-0.106901773 - 23.9722 ./ T) .* sqrt(S);
  tmp3    = 0.1130822 .* S - 0.00846934 .* S.^1.5 + log(1 - 0.001005 * S);
  lnK2roy = tmp1 + tmp2 + tmp3;
  K2      = exp(lnK2roy);

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end

  K2      = K2*exp(p.lnkpok0(2));

%% convert units (µmol/kg):
  K2      = K2*1.e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Kb=K_B(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dissociation constant for boric acid
% Note: calculated in mol/kg, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S.^(3./2.)-0.0996*S.*S);
  tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
  tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S).*log(T);
  lnKb = tmp1 ./ T + tmp2 + tmp3 + 0.053105*sqrt(S).*T;
  Kb   = exp(lnKb);

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end

  Kb   = Kb*exp(p.lnkpok0(3));

%% convert units (µmol/kg):
  Kb   = Kb*1.e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Kw=K_w(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ion product for water
% Note: calculated in mol^2/kg^2, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  Kw = exp( 148.96502-13847.26/T-23.6521*log(T) ...
           +(118.67/T-5.977+1.0495*log(T))*sqrt(S)-0.01615*S );

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end

  Kw = Kw*exp(p.lnkpok0(4));

%% convert units (µmol^2/kg^2):
  Kw = Kw*1.e12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=K_sp(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% choosing solubility product constant for CaCO3 (calcite/aragonite)

  if p.do_use_Ksp_aragonite==1
     K = K_sp_aragonite(p,P,T_deep,S);
  else
     K = K_sp_calcite(p,P,T_deep,S);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=K_sp_calcite(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solubility product constant for calcite
% Note: calculated in mol^2/kg^2, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10

tmp1      = -171.9065-0.077993.*T+2839.319./T+71.595.*log10(T);
tmp2      = +(-0.77712+0.0028426.*T+178.34./T).*sqrt(S);
tmp3      = -0.07711.*S+0.0041249.*S.^1.5;
log10Kspc = tmp1 + tmp2 + tmp3;
Kspc      = 10.^(log10Kspc);

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end

Kspc      = Kspc*exp(p.lnkpok0(5));

%% convert units (µmol^2/kg^2):
out       = Kspc*1.e12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Kspa=K_sp_aragonite(p,P,T_deep,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solubility product constant for aragonite
% Note: calculated in mol^2/kg^2, converted at end of function

  TC      = T_deep;
  T       = TC + 273.15;
  P       = p.K_sp_pressure;
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10

tmp1      = -171.945-0.077993.*T+2903.293./T+71.595.*log10(T);
tmp2      = +(-0.068393+0.0017276.*T+88.135./T).*sqrt(S);
tmp3      = -0.10018.*S+0.0059415.*S.^1.5;
log10Kspa = tmp1 + tmp2 + tmp3;
Kspa      = 10.^(log10Kspa);

%% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
  end

Kspa      = Kspa*exp(p.lnkpok0(6));

%% convert units (µmol^2/kg^2):
Kspa      = Kspa*1.e12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% carbonate system solvers -----------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
          =csys_iterative(t,p,Alk,DIC,B_T,T_deep,S,init_guess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent n_error_messages

%% calculate dissolution coefficients to save operations below:
  P  = p.K_sp_pressure; % atm
  K0 = K_0(p,P,T_deep,S);
  K1 = K_1(p,P,T_deep,S);
  Kb = K_B(p,P,T_deep,S);
  K2 = K_2(p,P,T_deep,S);
  Kw = K_w(p,P,T_deep,S);

%% initialize a=1/[H^+] and ALK_C from previous time step:
  a                 = 1/init_guess(1);
  max_csys_iterates = 50;
  min_csys_iterates = 5;
  if isempty(n_error_messages); n_error_messages=0; end
    iterate=1;
    epsilon_a=1.0;
  while iterate<min_csys_iterates || ( iterate<max_csys_iterates && epsilon_a > 1.0e-6 )
    iterate=iterate+1;
%% calculate carbonate alkalinity:
    ALK_C     = Alk-( (a*Kb*B_T)/(1+a*Kb) + a*Kw - 1/a );
%% evaluate gamma:
    gamma     = ALK_C/DIC;
%% solve for new a=1/[H^+]:
    a_new     = (-(1-gamma)*K1...
               +sqrt( (1-gamma)^2*K1^2 ...
                     + 4*(2-gamma)*K1*K2*gamma )...
                ) / ( 2*(2-gamma)*K1*K2 );
%% evaluate Del-a in this iteration:
    epsilon_a = abs(a_new-a)/a;
    a         = a_new;
  end

%% track iterations:  
  if iterate == max_csys_iterates && n_error_messages<6
    n_error_messages=n_error_messages+1;
    fprintf(1,['reached max carbonate system iterates,' ...
               ' t=%g, iterate=%d, Alk=%8.2f, DIC=%8.2f\n'],t,iterate,Alk,DIC);
   elseif n_error_messages==6
    n_error_messages=n_error_messages+1;
     fprintf(1,'not printing any more maximum number of iterates warnings.\n');
  end

%% Now solve for all other carbonate system variables:
  H       = 1/a;
  pH      = -log10(H)+6;  
  OH      = Kw*a;
  H2CO3   = (1/(a*K1+a^2*K1*K2+1))*DIC;
  HCO3    = ((a*K1)/(a*K1+a^2*K1*K2+1))*DIC;
  CO3     = ( (a^2*K1*K2)...
             /(a*K1+a^2*K1*K2+1) )*DIC;
  CO2g    = H2CO3/K0;
  BOH4    = ((a*Kb)/(1+a*Kb))*B_T;
  BOH3    = BOH4/(a*Kb);
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
          =csys_gradient(t,p,Alk,DIC,B_T,T_deep,S,init_guess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  p_fsolve.DIC = DIC;
  p_fsolve.B_T = B_T;
  p_fsolve.Alk = Alk;
  p_fsolve.T   = T_deep;
  p_fsolve.S   = S;

  P = p.K_sp_pressure;
  T = T_deep;
  p = pressure_effect_on_K(t,p);
  S = p.Salinity0;
  
  K.K_w = K_w(p,P,T_deep,S);
  K.K_0 = K_0(p,P,T_deep,S);
  K.K_1 = K_1(p,P,T_deep,S);
  K.K_2 = K_2(p,P,T_deep,S);
  K.K_B = K_B(p,P,T_deep,S);

%% nondimensionalize initial guess:
  init_guess = (init_guess-p.carb_offset)./p.carb_scale;
  [x,fval,exitflag,output,J] = fsolve(@(X)csys_eqns(X,t,p,p_fsolve,K),init_guess,p.optim_options);

%% dimensionalize the fsolve solution:
  X = x.*p.carb_scale+p.carb_offset;

%% specify output arguments:
  H       = X(1);
  OH      = X(2);
  CO2g    = X(3);
  H2CO3   = X(4);
  HCO3    = X(5);
  CO3     = X(6);
  BOH3    = X(7);
  BOH4    = X(8);
  pH      = -log10(H)+6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,J]=csys_eqns(X,t,p,p_fsolve,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% returns carbonate system equations and Jacobian

%% get variables from unknown vector X:  X is nondimensional, so
%% for example X(1)=(Hp-p.carb_sulf_offset(1))/p.carb_sulf_scale(1):
  H       = X(1) * p.carb_scale(1) + p.carb_offset(1);
  OH      = X(2) * p.carb_scale(2) + p.carb_offset(2);
  CO2g    = X(3) * p.carb_scale(3) + p.carb_offset(3);
  H2CO3   = X(4) * p.carb_scale(4) + p.carb_offset(4);
  HCO3    = X(5) * p.carb_scale(5) + p.carb_offset(5);
  CO3     = X(6) * p.carb_scale(6) + p.carb_offset(6);
  BOH3    = X(7) * p.carb_scale(7) + p.carb_offset(7);
  BOH4    = X(8) * p.carb_scale(8) + p.carb_offset(8);

%% get some additional quantities:
  DIC   = p_fsolve.DIC;
  B_T   = p_fsolve.B_T;
  Alk   = p_fsolve.Alk;

%% Calculate equations:
  F = [ ...
      K.K_w-H*OH,
      K.K_0*CO2g-H2CO3,
      K.K_1*H2CO3-H*HCO3,
      K.K_2*HCO3-H*CO3,
      K.K_B*BOH3-H*BOH4,
      H2CO3+HCO3+CO3-DIC,
      BOH4+BOH3-B_T,
      HCO3+2*CO3+BOH4+OH-H-Alk ...
      ]';
  
%% Calculate Jacobian:
  
  J=[  -OH,   -H,    0,	     0,	     0,	    0,	   0,     0  ;
         0,	   0,  K.K_0,   -1,	     0,	    0,	   0,     0  ;
      -HCO3,   0,    0,	   K.K_1,   -H,	    0,	   0,	  0  ;
       -CO3,   0,    0,	     0,	   K.K_2,  -H,     0,	  0  ;
      -BOH4,   0,    0,	     0,	     0,	    0,	 K.K_B,  -H  ;
         0,    0,    0,      1,      1,     1,     0,     0  ;
         0,    0,    0,      0,      0,     0,     1,     1  ;
        -1,    1,    0,      0,      1,     2,     0,     1 ];

  if 0
    %% Add inequality penalties for negative carbonate variables:
    %% dimensional variables:
      Xdim = X.*p.carb_scale+p.carb_offset;
    %% find indices of negative variables:
      I_penalty = find(Xdim<0);
    %% add penalty:
      F(I_penalty) = F(I_penalty)+Xdim(I_penalty).^2;
    %% add to Jacobian:
      for i=I_penalty
        J(i,i)=J(i,i)+2*Xdim(i);
      end
  end
  
  %% scale the Jacobian:
  N=length(F);
  for i=1:N
    J(i,:)=J(i,:).*p.carb_scale(:)';
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
          =csys_approx(t,p,Alk,DIC,B_T,T_deep,S,init_guess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  P     = 1;
  HCO3  = 2*DIC-Alk;
  CO3   = Alk-DIC;
  H     = K_2(p,P,T_deep,S)*(2*DIC-Alk)/(Alk-DIC);
  H2CO3 = (K_2(p,P,T_deep,S)/K_1(p,P,T_deep,S))*(2*DIC-Alk)^2/(Alk-DIC);
  CO2g  = H2CO3/K_0(p,P,T_deep,S);
  pH    = -log10(H)+6;
  OH    = NaN;
  HCO3  = NaN;
  CO3   = NaN;
  BOH3  = NaN;
  BOH4  = NaN;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pH,H,OH,CO2g,H2CO3,HCO3,CO3,BOH3,BOH4] ...
          =csys_csys3(t,p,Alk,DIC,B_T,T_deep,S,init_guess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate coefficients:
  P   = p.K_sp_pressure; % atm
  K0  = K_0(p,P,T_deep,S);
  K1  = K_1(p,P,T_deep,S);
  Kb  = K_B(p,P,T_deep,S);
  K2  = K_2(p,P,T_deep,S);
  Kw  = K_w(p,P,T_deep,S);
  B_T = p.B_T0;
  
%% csys3 routine for solving carbonate system:
  p5  = -1.;
  p4  = -Alk-Kb-K1;
  p3  = DIC*K1-Alk*(Kb+K1)+Kb*B_T+Kw-Kb*K1-K1*K2;
  tmp = DIC*(Kb*K1+2.*K1*K2)-Alk*(Kb*K1+K1*K2)+Kb*B_T*K1;
  p2  = tmp+(Kw*Kb+Kw*K1-Kb*K1*K2);
  tmp = 2.*DIC*Kb*K1*K2-Alk*Kb*K1*K2+Kb*B_T*K1*K2;
  p1  = tmp+(+Kw*Kb*K1+Kw*K1*K2);
  p0  = Kw*Kb*K1*K2;
  p   = [p5 p4 p3 p2 p1 p0];
  r   = roots(p);
  
%% calculate other carbonate system variables:
  H     = max(real(r));
  pH    = -log10(H)+6;
  OH    = Kw/H;
  H2CO3 = DIC / (1.+K1/H+K1*K2/H/H);
  HCO3  = DIC / (1+H/K1+K2/H);
  CO3   = DIC / (1+H/K2+H*H/K1/K2);
  CO2g  = H2CO3/K0;
  BOH4  = (B_T*Kb)/(Kb+H);
  BOH3  = (H*BOH4)/Kb;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_spinup(t,x,p,spin,spinup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if spin == 1
  save(sprintf('spinup/spinup-%2.2d.mat',p.scenario),'t','x','spinup');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,spinup]=read_spinup(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename=sprintf('spinup/spinup-%2.2d.mat',p.scenario);
  if exist(filename,'file')==0
    disp(' ');
    fprintf(1,'\n *** spinup file does not exist: %s \n\n',filename);
    x=NaN; t=NaN;
    spinup.J_b_carb=NaN; 
    spinup.J_b_crust=NaN;
  else
    load(filename);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=set_parameters(scenario)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  t=0;
  p.scenario=scenario;
  fprintf(1,'running carb_chem scenario %2.2d\n',p.scenario);

%% read scenario parameters:
  addpath input;  % Could also put this line in my ~/matlab/startup.m
  eval(sprintf('scenario%2.2d',p.scenario));
