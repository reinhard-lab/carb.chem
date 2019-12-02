%
tic
%
% navigate to home directory
cd ~/Documents/MATLAB/carb.chem/
%
% make directory for figures [if none exists]
if ~exist('~/Documents/MATLAB/carb.chem/figures/','dir')
    mkdir('~/Documents/MATLAB/carb.chem/figures/')
end
%
%% Main Text Figure 2, Supplementary Figures S1+S2
% load and parse data, including spinup data for forcing normalization
cd ~/Documents/MATLAB/carb.chem/spinup
%
load('spinup-111');
J_py0_111  = spinup.J_w_py;
%
load('spinup-121');
J_py0_121  = spinup.J_w_py;
J_org0_121 = spinup.J_w_org;
%
load('spinup-131');
J_py0_131  = spinup.J_w_py;
J_org0_131 = spinup.J_w_org;
J_sil0_131 = spinup.J_w_sil;
J_ev0_131  = spinup.J_w_ev;
%
cd ~/Documents/MATLAB/carb.chem/output
%
load('out_scenario111');
out_111    = o;
params_111 = p;
%
load('out_scenario112');
out_112    = o;
params_112 = p;
%
load('out_scenario113');
out_113    = o;
params_113 = p;
%
load('out_scenario114');
out_114    = o;
params_114 = p;
%
load('out_scenario115');
out_115    = o;
params_115 = p;
%
load('out_scenario116');
out_116    = o;
params_116 = p;
%
load('out_scenario121');
out_121    = o;
params_121 = p;
%
load('out_scenario122');
out_122    = o;
params_122 = p;
%
load('out_scenario123');
out_123    = o;
params_123 = p;
%
load('out_scenario124');
out_124    = o;
params_124 = p;
%
load('out_scenario125');
out_125    = o;
params_125 = p;
%
load('out_scenario126');
out_126    = o;
params_126 = p;
%
load('out_scenario131');
out_131    = o;
params_131 = p;
%
load('out_scenario132');
out_132    = o;
params_132 = p;
%
load('out_scenario133');
out_133    = o;
params_133 = p;
%
load('out_scenario134');
out_134    = o;
params_134 = p;
%
load('out_scenario135');
out_135    = o;
params_135 = p;
%
load('out_scenario136');
out_136    = o;
params_136 = p;
%
% make and save figures 
% Main Text Figure 2
figure(1);
%
subplot(5,3,1);
plot(out_111.t./1.e6,out_111.J_w_py./J_py0_111,'g','linewidth',1.5);
xlim([2 10]);
ylim([1 2]);
xlabel('time [Myr]');
ylabel('J_i/J_i^0');
box on;
grid on;
%
subplot(5,3,2);
plot(out_121.t./1.e6,out_121.J_w_py./J_py0_121,'g','linewidth',1.5);
hold on;
plot(out_121.t./1.e6,out_121.J_w_org./J_org0_121,'b--','linewidth',1.5);
xlim([2 10]);
ylim([1 2]);
xlabel('time [Myr]');
ylabel('J_i/J_i^0');
box on;
grid on;
%
subplot(5,3,3);
plot(out_131.t./1.e6,out_131.J_w_py./J_py0_131,'g','linewidth',1.5);
hold on;
plot(out_131.t./1.e6,out_131.J_w_org./J_org0_131,'b--','linewidth',1.5);
plot(out_131.t./1.e6,out_131.J_w_sil./J_sil0_131,'k-','linewidth',1.5);
plot(out_131.t./1.e6,out_131.J_w_ev./J_ev0_131,'y','linewidth',1.5);
xlim([2 10]);
ylim([1 2]);
xlabel('time [Myr]');
ylabel('J_i/J_i^0');
box on;
grid on;
%
subplot(5,3,[4,7]);
plot(out_111.t./1.e6,out_111.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_112.t./1.e6,out_112.CO2g,'b--','linewidth',1.5);
plot(out_113.t./1.e6,out_113.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(5,3,[5,8]);
plot(out_121.t./1.e6,out_121.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_122.t./1.e6,out_122.CO2g,'b--','linewidth',1.5);
plot(out_123.t./1.e6,out_123.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(5,3,[6,9]);
plot(out_131.t./1.e6,out_131.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_132.t./1.e6,out_132.CO2g,'b--','linewidth',1.5);
plot(out_133.t./1.e6,out_133.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(5,3,[10,13]);
plot(out_114.t./1.e6,out_114.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_115.t./1.e6,out_115.CO2g,'r--','linewidth',1.5);
plot(out_116.t./1.e6,out_116.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(5,3,[11,14]);
plot(out_124.t./1.e6,out_124.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_125.t./1.e6,out_125.CO2g,'r--','linewidth',1.5);
plot(out_126.t./1.e6,out_126.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(5,3,[12,15]);
plot(out_134.t./1.e6,out_134.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_135.t./1.e6,out_135.CO2g,'r--','linewidth',1.5);
plot(out_136.t./1.e6,out_136.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Main_Text_Figure_2.pdf');
%
% Supplementary Figure S1
figure(2);
%
subplot(2,3,1);
plot(out_111.t./1.e6,out_111.d13C_sw,'b-','linewidth',1.5);
hold on;
plot(out_112.t./1.e6,out_112.d13C_sw,'b--','linewidth',1.5);
plot(out_113.t./1.e6,out_113.d13C_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,2);
plot(out_121.t./1.e6,out_121.d13C_sw,'b-','linewidth',1.5);
hold on;
plot(out_122.t./1.e6,out_122.d13C_sw,'b--','linewidth',1.5);
plot(out_123.t./1.e6,out_123.d13C_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,3);
plot(out_131.t./1.e6,out_131.d13C_sw,'b-','linewidth',1.5);
hold on;
plot(out_132.t./1.e6,out_132.d13C_sw,'b--','linewidth',1.5);
plot(out_133.t./1.e6,out_133.d13C_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,4);
plot(out_114.t./1.e6,out_114.d13C_sw,'r-','linewidth',1.5);
hold on;
plot(out_115.t./1.e6,out_115.d13C_sw,'r--','linewidth',1.5);
plot(out_116.t./1.e6,out_116.d13C_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,5);
plot(out_124.t./1.e6,out_124.d13C_sw,'r-','linewidth',1.5);
hold on;
plot(out_125.t./1.e6,out_125.d13C_sw,'r--','linewidth',1.5);
plot(out_126.t./1.e6,out_126.d13C_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,6);
plot(out_134.t./1.e6,out_134.d13C_sw,'r-','linewidth',1.5);
hold on;
plot(out_135.t./1.e6,out_135.d13C_sw,'r--','linewidth',1.5);
plot(out_136.t./1.e6,out_136.d13C_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([-4.0 3.0]);
xlabel('time [Myr]');
ylabel('d^{13}C_{sw} [o/oo]');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Supplementary_Figure_S1.pdf');
%
% Supplementary Figure S2
figure(3);
%
subplot(2,3,1);
plot(out_111.t./1.e6,out_111.d34S_sw,'b-','linewidth',1.5);
hold on;
plot(out_112.t./1.e6,out_112.d34S_sw,'b--','linewidth',1.5);
plot(out_113.t./1.e6,out_113.d34S_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,2);
plot(out_121.t./1.e6,out_121.d34S_sw,'b-','linewidth',1.5);
hold on;
plot(out_122.t./1.e6,out_122.d34S_sw,'b--','linewidth',1.5);
plot(out_123.t./1.e6,out_123.d34S_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,3);
plot(out_131.t./1.e6,out_131.d34S_sw,'b-','linewidth',1.5);
hold on;
plot(out_132.t./1.e6,out_132.d34S_sw,'b--','linewidth',1.5);
plot(out_133.t./1.e6,out_133.d34S_sw,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,4);
plot(out_114.t./1.e6,out_114.d34S_sw,'r-','linewidth',1.5);
hold on;
plot(out_115.t./1.e6,out_115.d34S_sw,'r--','linewidth',1.5);
plot(out_116.t./1.e6,out_116.d34S_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,5);
plot(out_124.t./1.e6,out_124.d34S_sw,'r-','linewidth',1.5);
hold on;
plot(out_125.t./1.e6,out_125.d34S_sw,'r--','linewidth',1.5);
plot(out_126.t./1.e6,out_126.d34S_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
subplot(2,3,6);
plot(out_134.t./1.e6,out_134.d34S_sw,'r-','linewidth',1.5);
hold on;
plot(out_135.t./1.e6,out_135.d34S_sw,'r--','linewidth',1.5);
plot(out_136.t./1.e6,out_136.d34S_sw,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([13 19]);
xlabel('time [Myr]');
ylabel('d^{34}S_{sw} [o/oo]');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Supplementary_Figure_S2.pdf');
%
% clear output
clear;
%
% directory control
cd ~/Documents/MATLAB/carb.chem/output/
%
%% Main Text Figure 3
% load and parse data
load('out_scenario311');
out_311    = o;
params_311 = p;
%
load('out_scenario312');
out_312    = o;
params_312 = p;
%
load('out_scenario313');
out_313    = o;
params_313 = p;
%
load('out_scenario314');
out_314    = o;
params_314 = p;
%
load('out_scenario315');
out_315    = o;
params_315 = p;
%
load('out_scenario316');
out_316    = o;
params_316 = p;
%
%% some downstream calculations comparing alkalinity fluxes for time-dependent runs
J_Alk_in_temp  = (2.*out_311.J_w_sil) + (2.*out_311.J_w_carb) ...
                +(2*params_311.f_py_CaCO3_0.*out_311.J_w_py) + (2.*out_311.J_b_py);
J_Alk_out_temp = (2.*out_311.J_b_carb) + (2*params_311.f_py_H_0.*out_311.J_w_py); 
r_Alk_311      = J_Alk_in_temp./J_Alk_out_temp;
%
J_Alk_in_temp  = (2.*out_312.J_w_sil) + (2.*out_312.J_w_carb) ...
                +(2*params_312.f_py_CaCO3_0.*out_312.J_w_py) + (2.*out_312.J_b_py);
J_Alk_out_temp = (2.*out_312.J_b_carb) + (2*params_312.f_py_H_0.*out_312.J_w_py); 
r_Alk_312      = J_Alk_in_temp./J_Alk_out_temp;
%
J_Alk_in_temp  = (2.*out_313.J_w_sil) + (2.*out_313.J_w_carb) ...
                +(2*params_313.f_py_CaCO3_0.*out_313.J_w_py) + (2.*out_313.J_b_py);
J_Alk_out_temp = (2.*out_313.J_b_carb) + (2*params_313.f_py_H_0.*out_313.J_w_py); 
r_Alk_313      = J_Alk_in_temp./J_Alk_out_temp;
%
J_Alk_in_temp  = (2.*out_314.J_w_sil) + (2.*out_314.J_w_carb) ...
                +(2*params_314.f_py_CaCO3_0.*out_314.J_w_py) + (2.*out_314.J_b_py);
J_Alk_out_temp = (2.*out_314.J_b_carb) + (2*params_314.f_py_H_0.*out_314.J_w_py); 
r_Alk_314      = J_Alk_in_temp./J_Alk_out_temp;
%
J_Alk_in_temp  = (2.*out_315.J_w_sil) + (2.*out_315.J_w_carb) ...
                +(2*params_315.f_py_CaCO3_0.*out_315.J_w_py) + (2.*out_315.J_b_py);
J_Alk_out_temp = (2.*out_315.J_b_carb) + (2*params_315.f_py_H_0.*out_315.J_w_py); 
r_Alk_315      = J_Alk_in_temp./J_Alk_out_temp;
%
J_Alk_in_temp  = (2.*out_316.J_w_sil) + (2.*out_316.J_w_carb) ...
                +(2*params_316.f_py_CaCO3_0.*out_316.J_w_py) + (2.*out_316.J_b_py);
J_Alk_out_temp = (2.*out_316.J_b_carb) + (2*params_316.f_py_H_0.*out_316.J_w_py); 
r_Alk_316      = J_Alk_in_temp./J_Alk_out_temp;
%
% make figure
figure(4);
%
subplot(3,2,[1,3]);
plot(out_311.pO2,out_311.pH,'b-','linewidth',1.5);
hold on;
plot(out_312.pO2,out_312.pH,'b--','linewidth',1.5);
plot(out_313.pO2,out_313.pH,'b-.','linewidth',1.5);
xlim([0.01 0.1]);
ylim([7.5 7.8]);
xlabel('pO_2 [atm]');
ylabel('pH');
box on;
grid on;
%
subplot(3,2,[2,4]);
plot(out_314.pO2,out_314.pH,'r-','linewidth',1.5);
hold on;
plot(out_315.pO2,out_315.pH,'r--','linewidth',1.5);
plot(out_316.pO2,out_316.pH,'r-.','linewidth',1.5);
xlim([0.01 0.1]);
ylim([7.6 7.9]);
xlabel('pO_2 [atm]');
ylabel('pH');
box on;
grid on;
%
subplot(3,2,5);
plot(out_311.t./1.e6,r_Alk_311,'b-','linewidth',1.5);
hold on;
plot(out_312.t./1.e6,r_Alk_312,'b--','linewidth',1.5);
plot(out_313.t./1.e6,r_Alk_313,'b-.','linewidth',1.5);
xlim([1 6]);
ylim([0.85 1.15]);
xlabel('pO_2 [atm]');
ylabel('pH');
box on;
grid on;
%
subplot(3,2,6);
plot(out_314.t./1.e6,r_Alk_314,'r-','linewidth',1.5);
hold on;
plot(out_315.t./1.e6,r_Alk_315,'r--','linewidth',1.5);
plot(out_316.t./1.e6,r_Alk_316,'r-.','linewidth',1.5);
xlim([1 6]);
ylim([0.9 1.1]);
xlabel('pO_2 [atm]');
ylabel('pH');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Main_Text_Figure_3.pdf');
%
% clear output
clear;
%
% directory control
cd ~/Documents/MATLAB/carb.chem/output/
%
%% Main Text Figure 4
% load and parse data
load('out_scenario411');
out_411    = o;
params_411 = p;
%
load('out_scenario412');
out_412    = o;
params_412 = p;
%
load('out_scenario413');
out_413    = o;
params_413 = p;
%
load('out_scenario414');
out_414    = o;
params_414 = p;
%
load('out_scenario415');
out_415    = o;
params_415 = p;
%
load('out_scenario416');
out_416    = o;
params_416 = p;
%
load('out_scenario421');
out_421    = o;
params_421 = p;
%
load('out_scenario422');
out_422    = o;
params_422 = p;
%
load('out_scenario423');
out_423    = o;
params_423 = p;
%
load('out_scenario424');
out_424    = o;
params_424 = p;
%
load('out_scenario425');
out_425    = o;
params_425 = p;
%
load('out_scenario426');
out_426    = o;
params_426 = p;
%
% make figure
figure(5);
%
subplot(2,3,1);
plot(out_421.t./1.e6,out_421.pH-out_411.pH,'b-','linewidth',1.5);
hold on;
plot(out_422.t./1.e6,out_422.pH-out_412.pH,'b--','linewidth',1.5);
plot(out_423.t./1.e6,out_423.pH-out_413.pH,'b-.','linewidth',1.5);
xlim([0 20]);
ylim([-0.05 0.3]);
xlabel('time [Myr]');
ylabel('?_{pH}');
box on;
grid on;
%
subplot(2,3,2);
plot(out_421.t./1.e6,out_421.CO2g-out_411.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_422.t./1.e6,out_422.CO2g-out_412.CO2g,'b--','linewidth',1.5);
plot(out_423.t./1.e6,out_423.CO2g-out_413.CO2g,'b-.','linewidth',1.5);
xlim([0 20]);
ylim([-900 100]);
xlabel('time [Myr]');
ylabel('?_{CO2} [ppm]');
box on;
grid on;
%
subplot(2,3,3);
plot(out_421.t./1.e6,out_421.T_surf-out_411.T_surf,'b-','linewidth',1.5);
hold on;
plot(out_422.t./1.e6,out_422.T_surf-out_412.T_surf,'b--','linewidth',1.5);
plot(out_423.t./1.e6,out_423.T_surf-out_413.T_surf,'b-.','linewidth',1.5);
xlim([0 20]);
ylim([-1.75 0.25]);
xlabel('time [Myr]');
ylabel('?_{T} [ºC]');
box on;
grid on;
%
subplot(2,3,4);
plot(out_424.t./1.e6,out_424.pH-out_414.pH,'r-','linewidth',1.5);
hold on;
plot(out_425.t./1.e6,out_425.pH-out_415.pH,'r--','linewidth',1.5);
plot(out_426.t./1.e6,out_426.pH-out_416.pH,'r-.','linewidth',1.5);
xlim([0 20]);
ylim([-0.05 0.3]);
xlabel('time [Myr]');
ylabel('?_{pH}');
box on;
grid on;
%
subplot(2,3,5);
plot(out_424.t./1.e6,out_424.CO2g-out_414.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_425.t./1.e6,out_425.CO2g-out_415.CO2g,'r--','linewidth',1.5);
plot(out_426.t./1.e6,out_426.CO2g-out_416.CO2g,'r-.','linewidth',1.5);
xlim([0 20]);
ylim([-900 100]);
xlabel('time [Myr]');
ylabel('?_{CO2} [ppm]');
box on;
grid on;
%
subplot(2,3,6);
plot(out_424.t./1.e6,out_424.T_surf-out_414.T_surf,'r-','linewidth',1.5);
hold on;
plot(out_425.t./1.e6,out_425.T_surf-out_415.T_surf,'r--','linewidth',1.5);
plot(out_426.t./1.e6,out_426.T_surf-out_416.T_surf,'r-.','linewidth',1.5);
xlim([0 20]);
ylim([-1.75 0.25]);
xlabel('time [Myr]');
ylabel('?_{T} [ºC]');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Main_Text_Figure_4.pdf');
%
% clear output
clear;
%
% directory control
cd ~/Documents/MATLAB/carb.chem/output/
%
%% Supplementary Figure S3
% load and parse data
load('out_scenario211');
out_211    = o;
params_211 = p;
%
load('out_scenario212');
out_212    = o;
params_212 = p;
%
load('out_scenario213');
out_213    = o;
params_213 = p;
%
load('out_scenario214');
out_214    = o;
params_214 = p;
%
load('out_scenario215');
out_215    = o;
params_215 = p;
%
load('out_scenario216');
out_216    = o;
params_216 = p;
%
load('out_scenario221');
out_221    = o;
params_221 = p;
%
load('out_scenario222');
out_222    = o;
params_222 = p;
%
load('out_scenario223');
out_223    = o;
params_223 = p;
%
load('out_scenario224');
out_224    = o;
params_224 = p;
%
load('out_scenario225');
out_225    = o;
params_225 = p;
%
load('out_scenario226');
out_226    = o;
params_226 = p;
%
load('out_scenario231');
out_231    = o;
params_231 = p;
%
load('out_scenario232');
out_232    = o;
params_232 = p;
%
load('out_scenario233');
out_233    = o;
params_233 = p;
%
load('out_scenario234');
out_234    = o;
params_234 = p;
%
load('out_scenario235');
out_235    = o;
params_235 = p;
%
load('out_scenario236');
out_236    = o;
params_236 = p;
%
% make and save figure
figure(6);
%
subplot(2,3,1);
plot(out_211.t./1.e6,out_211.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_212.t./1.e6,out_212.CO2g,'b--','linewidth',1.5);
plot(out_213.t./1.e6,out_213.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(2,3,2);
plot(out_221.t./1.e6,out_221.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_222.t./1.e6,out_222.CO2g,'b--','linewidth',1.5);
plot(out_223.t./1.e6,out_223.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(2,3,3);
plot(out_231.t./1.e6,out_231.CO2g,'b-','linewidth',1.5);
hold on;
plot(out_232.t./1.e6,out_232.CO2g,'b--','linewidth',1.5);
plot(out_233.t./1.e6,out_233.CO2g,'b-.','linewidth',1.5);
xlim([2 10]);
ylim([500 4500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(2,3,4);
plot(out_214.t./1.e6,out_214.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_215.t./1.e6,out_215.CO2g,'r--','linewidth',1.5);
plot(out_216.t./1.e6,out_216.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(2,3,5);
plot(out_224.t./1.e6,out_224.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_225.t./1.e6,out_225.CO2g,'r--','linewidth',1.5);
plot(out_226.t./1.e6,out_226.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
subplot(2,3,6);
plot(out_234.t./1.e6,out_234.CO2g,'r-','linewidth',1.5);
hold on;
plot(out_235.t./1.e6,out_235.CO2g,'r--','linewidth',1.5);
plot(out_236.t./1.e6,out_236.CO2g,'r-.','linewidth',1.5);
xlim([2 10]);
ylim([0 2500]);
xlabel('time [Myr]');
ylabel('pCO_2 [ppm]');
box on;
grid on;
%
% print to figure
cd ~/Documents/MATLAB/carb.chem/figures/;
print(gcf, '-dpdf', 'Supplementary_Figure_S3.pdf');
%
% clear output
clear;
% return to home directory
cd ~/Documents/MATLAB/carb.chem/;
%
toc
%