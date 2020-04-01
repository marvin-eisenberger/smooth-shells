function [tau,X,Y,feat] = layeredmatching(X,Y,feat,kArray,aInit,param)
    % layeredmatching - main registration method for two input shapes X,Y
    %
    % X,Y: input shape
    % feat: dynamic exchange collection
    % kArray: coarse-to-fine shell sizes 
    % aInit: deformation initialization from initialization
    % param: parameter values

    if ~exist('param','var')
        param = standardparams();
    end

    param.noPlotInBetween = true;
    param.plotCorr = false;
    param.twoPlots = false;

    
    tau = zeros(kArray(end),3);
    if exist('aInit','var')
        tau(1:size(aInit,1),:) = aInit;
    end
    
    feat.wCurr = zeros(X.n,3);
    
    normalDamping = param.normalDamping;

    %% coarse-to-fine matching - all steps
    for iK = 1:length(kArray)
        
        if param.intermediateOutput
            disp('k = ' + string(kArray(iK)));
        end
        
        k = kArray(iK);
    
        param.lambdaArap = param.lambdaArapInit .* sqrt(max(k-20,0));
        param.normalDamping = normalDamping .* (k > 10);

        % if a functional map from the initialization is used, it should
        % only be used until the 20th step for maximum accuracy
        if k > 20 && isfield(feat,'Cfix')
            feat = rmfield(feat,'Cfix');
        end
        
        %call sub matching step
        [X,Y,tau,feat] = matchlayer(X,Y,feat,tau,k,param);
        
    end







    
    
    
    
        
        
        
        
