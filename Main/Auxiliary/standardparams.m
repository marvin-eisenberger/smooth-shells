function param = standardparams(param)

if ~exist('param','var')
    param = struct;
end


if ~isfield(param,'areaPreservation')
    param.areaPreservation = true;
end

if ~isfield(param,'facFeat')
    param.facFeat = 0.25;
end

if ~isfield(param,'GPUcorrespondences')
    param.GPUcorrespondences = true;
end

if ~isfield(param,'kArrayLength')
    param.kArrayLength = 50;
end

if ~isfield(param,'kMax')
    param.kMax = 500;
end

if ~isfield(param,'lambdaArap')
    param.lambdaArap = 0;
end

if ~isfield(param,'lambdaArapInit')
    param.lambdaArapInit = 0.02;
end

if ~isfield(param,'lambdaFeat')
    param.lambdaFeat = 25;
end

if ~isfield(param,'lambdaLap')
    param.lambdaLap = 0;
end

if ~isfield(param,'matchingAlignment')
    param.matchingAlignment = false;
end

if ~isfield(param,'mode')
    param.mode = 2;
end

if ~isfield(param,'noRigid')
    param.noRigid = false;
end

if ~isfield(param,'normalDamping')
    param.normalDamping = 0.04;
end

if ~isfield(param,'numSub')
    param.numSub = 1;
end

if ~isfield(param,'problemSizeInit')
    param.problemSizeInit = 1000;
end

%output
if ~isfield(param,'intermediateOutput')
    param.intermediateOutput = true;
end

%plot
if ~isfield(param,'noPlot')
    param.noPlot = true;
end

%plot
if ~isfield(param,'noPlotInBetween')
    param.noPlotInBetween = true;
end
%plot
if ~isfield(param,'plotCorr')
    param.plotCorr = true;
end

%plot
if ~isfield(param,'showOff')
    param.showOff = true;
end

%plot
if ~isfield(param,'twoPlots')
    param.twoPlots = false;
end

end
