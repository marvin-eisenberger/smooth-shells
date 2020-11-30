%specify input shapes
file1 = 'path/to/first/shape';
file2 = 'path/to/second/shape';

%initialize parameters (optional)
param = struct;
param = standardparams(param);

param.noPlot = false; %turn on/off for intermediate plots
param.GPUcorrespondences = true; %recommended option = true. Turn off, if no GPU is available

%method parameters with recommended settings for interclass shapes
param.facFeat = 0.25;
param.kArrayLength = 50;
param.kMax = 500;
param.lambdaArapInit = 0.02;
param.lambdaFeat = 0.1;
param.normalDamping = 0.04;
param.numMCMC = 500;
param.areaPreservation = false;
param.noRigid = true; %set false, if rigid pre-alignment is needed

%execute main script
[P,tau,C,X,Y] = smoothshells(file1,file2,param);