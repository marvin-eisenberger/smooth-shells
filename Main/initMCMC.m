function [Cfix,ainit,iMin,errorArrayGT,errorArray] = initMCMC(X,Y,param,kMin,kMax)

    % initMCMC - Markov chain Monte Carlo initialization strategy
    %
    % X,Y: input shapes
    % param: parameter values
    % kMin,kMax: range of shell sizes for surrogate runs
    
    if ~exist('kMin','var')
        kMin = 6;
    end
    
    if ~exist('kMax','var')
        kMax = 20;
    end
    
    kArray = kMin:kMax;
    
    param = standardparams(param);
    param.facFeat = 1.5;
    param.lambdaFeat = 0;
    param.plotCorr = false;
    param.twoPlots = false;
    param.mode = 2;
    param.intermediateOutput = false;
    
    
    numProp = param.numMCMC;
    
    feat = struct;
    feat.wCurr = zeros(X.n,3);
    
    
    [XSmooth,feat.basisWeights] = smoothshapesigmoid(X,kMax);
    feat.oldBasisWeights = feat.basisWeights;
    YSmooth = smoothshapesigmoid(Y,kMax);
    
    
    errorArray = zeros(numProp,1);
    errorArrayGT = zeros(numProp,1);
    CfixCollection = cell(numProp,1);
    aCollection = cell(numProp,1);
    
    problemSize = param.problemSizeInit;
    X.samples = fps_euclidean(X.vert, problemSize, randi(X.n));
    X.samples = sort(X.samples);
    Y.samples = fps_euclidean(Y.vert, problemSize, randi(Y.n));
    Y.samples = sort(Y.samples);
    X.vertSub = X.vert(X.samples,:);
    Y.vertSub = Y.vert(Y.samples,:);
    
    if kMax > 20
        X.neigh.mat = getNeighbors(X.vert');
        [X.neigh.row,X.neigh.col] = find(X.neigh.mat(X.samples,:));
    end
    
    %% MCMC
    for iSet = 1:numProp
        %% create proposal shell
        disp("Prop #" + string(iSet) + "...")
        
        
        modeCurr = randi(2);
        lambdaLapSet = [0,10];
        param.lambdaLap = lambdaLapSet(modeCurr);
        ainit = proposalshell(X,kMin);
        
        paramComp = param;
        paramComp.noPlot = true;
        [a,~,~,featOut] = layeredmatching(X,Y,feat,kArray,ainit,paramComp);
        
        vertCurrOrig = shiftedvertupsample(X,a);
        featOut.normalCurr = compute_normal(vertCurrOrig,X.triv',X.flipNormal)';
        
        %% determine whether to accept the current proposal
        
        paramEval = param;
        paramEval.matchingAlignment = true;
        paramEval.facFeat = 0.11;
        
        featOut = computeCorrCurr(X,Y,vertCurrOrig,paramEval,featOut);
        
        errorArray(iSet) = 1e7 .* mean([diag(X.A(X.samples,X.samples)) .* featOut.Dass;diag(Y.A(Y.samples,Y.samples)) .* featOut.Dassinv]);
        
        if size(vertCurrOrig,1) == size(YSmooth.vert,1)
            currErrorGT = mean(normv(vertCurrOrig-YSmooth.vert));
            errorArrayGT(iSet) = currErrorGT;
        end
        
        CfixCollection{iSet} = featOut.C;
        aCollection{iSet} = a;
        
        iMin = find(errorArray(1:iSet) == min(errorArray(1:iSet)),1);
        
        if ~param.noPlot
            
            if iMin == iSet
                plotskeletonlayered(X,Y,XSmooth,YSmooth,param,0,0,a,0,kMax);
            end
            
            title('MCMC - init: objective = ' + string(errorArray(iSet)) + ', min objective = ' + string(errorArray(iMin)) + ', iMin = ' + string(iMin) + ', iCurr = ' + string(iSet) + "/" + string(numProp));
            drawnow
        end
        
        disp("...objective value: " + string(errorArray(iSet)))
    end
    
    %% choose best configuration
    iMin = find(errorArray == min(errorArray),1);
    
    Cfix = CfixCollection{iMin};
    ainit = aCollection{iMin};
    

    
    
    function featOut = computeCorrCurr(X,Y,vertCurrOrig,paramEval,featOut)
        [assignment,~,featOut] = computecorrespondences(X,Y,vertCurrOrig(X.samples,:),vertCurrOrig,paramEval,featOut);
        featOut.assignment = assignment;
    end
end






