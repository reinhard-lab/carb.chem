%
tic
%
% remove any existing output from spinup
cd ~/Documents/MATLAB/carb.chem/output/;
delete *.mat;
delete *.mat;
%
% directory control
cd ~/Documents/MATLAB/carb.chem/
%
% make directory for spinup files [if none exists]
if ~exist('spinup','dir')
    mkdir('spinup')
end
%
% make directory for output files [if none exists]
if ~exist('output','dir')
    mkdir('output')
end
%
% load and run configs for experimental runs
%
% set spinup parameters
spin = 0;
time = 30.e6;
%
% 2x pyrite weathering
for i=111:116
    carb_chem(i,spin,time);
end
%
% 2x pyrite + organic carbon weathering
for i=121:126
    carb_chem(i,spin,time);
end
%
% 0.5x all crustal phases
for i=131:136
    carb_chem(i,spin,time);
end
%
% 2x pyrite weathering, power function
for i=211:216
    carb_chem(i,spin,time);
end
%
% 2x pyrite + organic carbon weathering, power function
for i=221:226
    carb_chem(i,spin,time);
end
%
% 0.5x all crustal phases, power function
for i=231:236
    carb_chem(i,spin,time);
end
%
% transient perturbation scenarios
for i=311:316
    carb_chem(i,spin,time);
end
%
% volcanic outgassing, CO2 only
for i=411:416
    carb_chem(i,spin,time);
end
%
% volcanic outgassing, CO2 + red
for i=421:426
    carb_chem(i,spin,time);
end
%
clear
%
toc
%