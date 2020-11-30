%specify input shapes
iF1 = 2;
iF2 = 55;

file1 = 'Data/SCAPE/mesh' + string(num2str(iF1,'%03d')) + '.off';
file2 = 'Data/SCAPE/mesh' + string(num2str(iF2,'%03d')) + '.off';

%initialize parameters (optional)
param = struct;
param = standardparams(param);

param.noPlot = false; %turn on/off for intermediate plots
param.GPUcorrespondences = true; %recommended option = true. Turn off, if no GPU is available

%method parameters with recommended settings
param.facFeat = 0.25;
param.kArrayLength = 50;
param.kMax = 500;
param.lambdaArapInit = 0.02;
param.lambdaFeat = 25;
param.normalDamping = 0.04;
param.numMCMC = 100;

%execute main script
[P,tau,C,X,Y] = smoothshells(file1,file2,param);