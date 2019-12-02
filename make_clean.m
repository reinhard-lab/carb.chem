%
tic
%
% careful now...this will delete all spinups, output files, and figures
%
% navigate to home directory
cd ~/Documents/MATLAB/carb.chem/;
%
% remove spinups
if exist('~/Documents/MATLAB/carb.chem/spinup/','dir')
    rmdir('~/Documents/MATLAB/carb.chem/spinup/','s')
end
%
% remove output
if exist('~/Documents/MATLAB/carb.chem/output/','dir')
    rmdir('~/Documents/MATLAB/carb.chem/output/','s')
end
%
% remove figures
if exist('~/Documents/MATLAB/carb.chem/figures/','dir')
    rmdir('~/Documents/MATLAB/carb.chem/figures/','s')
end
%
cd ~/Documents/MATLAB/carb.chem/;
%
% remove spinup and output directories
toc
%