function [POut,tauOut,COut,X,Y] = smoothshells(input1,input2,param)
% shells - load a pair of input shapes and compute the features and eigenpairs
% input1,input2: either a mesh (as a struct) OR a file containing a mesh

if ~exist('param','var')
    param = struct;
end

param = standardparams(param);
kArray = [20,20+round(linspace(1,(param.kMax-20).^(1/4),param.kArrayLength).^4)];
feat = struct;

param.noPlot = false; %turn on/off for intermediate plots
% param.GPUcorrespondences = false; % turn on/off if the GPU should be used for speedup

%% load Shape

disp('Load shapes and compute features..')

[X,Y] = loadshapepair(input1,input2);

if ~param.noPlot
    surf_pair(X,Y);
    drawnow
end

%% rigid alignment

if ~param.noRigid
    disp('Computing rigid alignment...')
    [X] = alignrigid(X,Y,param);

    if ~param.noPlot
        surf_pair(X,Y);
        drawnow
    end
end

%% MCMC

disp('MCMC initialization...')

[feat.Cfix,tauInit] = initMCMC(X,Y,param);

%% full matching

disp('Full run...')

[tau,X,Y,featOut] = layeredmatching(X,Y,feat,kArray,tauInit,param);

vertCurrFull = X.vert + X.evecs(:,1:size(tau,1)) * tau;

 %% plot result

if ~param.noPlot

    %gt exists
    if size(vertCurrFull,1)==size(Y.vert,1)
        weightsPlotX = normv(vertCurrFull-Y.vert);
        weightsPlotY = weightsPlotX;
    else
        weightsPlotX = X.vert(:,2);
        weightsPlotY = Y.vert(:,2);
    end

    subplot(1,3,1)
    hold off
    trisurf(X.triv,vertCurrFull(:,1),vertCurrFull(:,2),vertCurrFull(:,3),weightsPlotX);
    axis equal
    colorbar
    title('Morphed shape $\hat{\mathcal{X}}$','interpreter','latex')

    subplot(1,3,2)
    hold off
    trisurf(X.triv,X.vert(:,1),X.vert(:,2),X.vert(:,3),weightsPlotX);
    axis equal
    colorbar
    title('Source shape $\mathcal{X}$','interpreter','latex')

    subplot(1,3,3)
    hold off
    trisurf(Y.triv,Y.vert(:,1),Y.vert(:,2),Y.vert(:,3),weightsPlotY);
    axis equal
    colorbar
    title('Reference shape $\mathcal{Y}$','interpreter','latex')

    if size(vertCurrFull,1)==size(Y.vert,1)
        disp('mean error: ' + string(mean(normv(vertCurrFull-Y.vert))));
    end

    drawnow
end

%% createOutput

tauOut = tau;
COut = featOut.C;
POut = struct;
POut.assignment = featOut.assignment;
POut.assignmentinv = featOut.assignmentinv;

end






















