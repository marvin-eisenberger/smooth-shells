function [X] = alignrigid(X,Y,param)
    % alignrigid - compute the rigid alignment of two shapes
    %
    % X,Y: input shapes
    % param: parameter values
    
    
    kArrayTestSkeleton = 1:6;
    
    %% main axis alignment
    
    XSmooth2 = smoothshape(X,2);
    YSmooth2 = smoothshape(Y,2);
    
    
    
    XSmooth2.vert = (XSmooth2.vert-mean(XSmooth2.vert));
    YSmooth2.vert = (YSmooth2.vert-mean(YSmooth2.vert));
    XSmooth2.vert = XSmooth2.vert ./ mean(XSmooth2.vert.^2);
    YSmooth2.vert = YSmooth2.vert ./ mean(YSmooth2.vert.^2);
    
    
    dir1 = XSmooth2.xi(2,:)';
    dir2 = YSmooth2.xi(2,:)';
    
    dir1 = dir1 ./ norm(dir1);
    dir2 = dir2 ./ norm(dir2);
    
    if acos(dir1' * dir2) <= pi/2
        Rw = rotvectorpairs(dir1',dir2');
    else
        Rw = rotvectorpairs(dir1',-dir2');
    end
    
    midpointCurrX = mean(XSmooth2.vert,1);
    midpointCurrY = mean(YSmooth2.vert,1);
    
    xiTranslate = transform_SE3_se3_homCoords([eye(3),midpointCurrX';zeros(1,3),1]);
    xiMidpoints = transform_SE3_se3_homCoords([eye(3),midpointCurrY' - midpointCurrX';zeros(1,3),1]);
    
    xiRot = transform_SE3_se3_homCoords(Rw);
    xi = groupMult(groupMult(xiTranslate,xiMidpoints),groupMult(xiRot,-xiTranslate));

    X.vert = rigidTransform(X.vert',xi')';
    
    
    XSmooth.xi = (X.evecs(:,kArrayTestSkeleton)' * X.A * X.vert);
    XSmooth.vert = X.evecs(:,kArrayTestSkeleton) * XSmooth.xi;
    YSmooth.xi = (Y.evecs(:,kArrayTestSkeleton)' * Y.A * Y.vert);
    YSmooth.vert = Y.evecs(:,kArrayTestSkeleton) * YSmooth.xi;
    XSmooth.vert = (XSmooth.vert-mean(XSmooth.vert));
    YSmooth.vert = (YSmooth.vert-mean(YSmooth.vert));
    XSmooth.vert = XSmooth.vert ./ mean(XSmooth.vert.^2);
    YSmooth.vert = YSmooth.vert ./ mean(YSmooth.vert.^2);
    XSmooth.triv = X.triv;
    YSmooth.triv = Y.triv;
    
    
    %% ICP around main axis
    
    problemSize = 1000;
    
    samplesX = fps_euclidean(XSmooth.vert, problemSize, randi(X.n));
    samplesX = sort(samplesX);

    samplesY = fps_euclidean(YSmooth.vert, problemSize, randi(Y.n));
    samplesY = sort(samplesY);


    xiICP = computeAxisICP(X.vert(samplesX,:),Y.vert(samplesY,:),5,YSmooth.xi(2,:));
    X.vert = rigidTransform(X.vert',xiICP')';
    
    
    weightAX = diag(X.A) ./ mean(diag(X.A));
    weightAY = diag(Y.A) ./ mean(diag(Y.A));
    X.vert = X.vert - mean(X.vert .* weightAX,1) + mean(Y.vert .* weightAY,1);
    
    
    
    %% surrogate runs
    
    %params init
    kMax = 500;
    
    XSmooth = smoothshapesigmoid(X,kMax);
    YSmooth = smoothshapesigmoid(Y,kMax);
    
    [X.normal,~,X.flipNormal] = compute_normal(XSmooth.vert',XSmooth.triv');
    [Y.normal,~,Y.flipNormal] = compute_normal(YSmooth.vert',YSmooth.triv');
    X.normal = X.normal';
    Y.normal = Y.normal';
    
    
    
    XsamplesOLD = X.samples;
    YsamplesOLD = Y.samples;
    
    param = standardparams(param);
    param.facFeat = 1000;
    param.noPlot = true;
    ainit = zeros(1,3);
    
    
    problemSize = param.problemSizeInit;
    X.samples = fps_euclidean(X.vert, problemSize, randi(X.n));
    X.samples = sort(X.samples);
    Y.samples = fps_euclidean(Y.vert, problemSize, randi(Y.n));
    Y.samples = sort(Y.samples);
    X.vertSub = X.vert(X.samples,:);
    Y.vertSub = Y.vert(Y.samples,:);
    
    
    feat = struct;
    feat.wCurr = zeros(X.n,3);
    
    kTest = 20;
    
    [~,feat.basisWeights] = smoothshapesigmoid(X,kTest);
    feat.oldBasisWeights = feat.basisWeights;
    
    kArray = 2:kTest;
    [feat] = computefunctionalmap(X,Y,param,feat,kTest);
    
    
    
    
    %first surrogate runs for vertical swap
    [coeffX,~,~] = pca(X.vert);
    principalDirX = coeffX(:,3);
    
    X = surrogateRun(X,Y,feat,param,kArray,principalDirX,2);
    
    
    %init for second run
    kTest = 20;
    
    [~,feat.basisWeights] = smoothshapesigmoid(X,kTest);
    feat.oldBasisWeights = feat.basisWeights;
    
    kArray = 2:kTest;
    [feat] = computefunctionalmap(X,Y,param,feat,kTest);
    
    
    %second surrogate runs for main axis rotation
    [coeffX,~,~] = pca(X.vert);
    principalDirX = coeffX(:,1);
    
    X = surrogateRun(X,Y,feat,param,kArray,principalDirX,4);
    
    
    
    X.samples = XsamplesOLD;
    Y.samples = YsamplesOLD;
    
    
    kMax = 500;
    XSmooth = smoothshapesigmoid(X,kMax);
    YSmooth = smoothshapesigmoid(Y,kMax);

    [X.normal,~,X.flipNormal] = compute_normal(XSmooth.vert',XSmooth.triv');
    [Y.normal,~,Y.flipNormal] = compute_normal(YSmooth.vert',YSmooth.triv');
    X.normal = X.normal';
    Y.normal = Y.normal';
    
    
    function X = surrogateRun(X,Y,feat,param,kArray,principalDirX,numSurr)
    
        kMax = 500;
        
        param = standardparams(param);
        param.lambdaLap = 1;
        param.numSub = 5;
        param.mode = 1;
        param.facFeat = 0.11;
        param.normalDamping = 0.1;
        param.intermediateOutput = false;
        
        ainit = zeros(1,3);

        Xsaves = cell(numSurr,1);

        for iSurr = 1:numSurr

            X.vert = rigidTransform(X.vert',[0,0,0,principalDirX' .* 2 .* pi .* 1 ./ numSurr]')';
            
            weightAX = diag(X.A) ./ mean(diag(X.A));
            weightAY = diag(Y.A) ./ mean(diag(Y.A));
            X.vert = X.vert - mean(X.vert .* weightAX,1) + mean(Y.vert .* weightAY,1);

            Xsaves{iSurr} = X;


            XSmooth = smoothshapesigmoid(X,kMax);
            YSmooth = smoothshapesigmoid(Y,kMax);

            [X.normal,~,X.flipNormal] = compute_normal(XSmooth.vert',XSmooth.triv');
            [Y.normal,~,Y.flipNormal] = compute_normal(YSmooth.vert',YSmooth.triv');
            X.normal = X.normal';
            Y.normal = Y.normal';

            %% compute matching

            [a,~,~,featOut] = layeredmatching(X,Y,feat,kArray,ainit,param);


            vertCurrOrig = shiftedvertupsample(X,a);
            featOut.normalCurr = compute_normal(vertCurrOrig,X.triv',X.flipNormal)';


            %% compute error

            paramEval = param;
            paramEval.matchingAlignment = true;


            featOut = computeCorrCurr(X,Y,vertCurrOrig,paramEval,featOut);


            
            errorArray(iSurr) = mean([diag(X.A(X.samples,X.samples)) .* featOut.Dass;diag(Y.A(Y.samples,Y.samples)) .* featOut.Dassinv]);
            

            if size(vertCurrOrig,1) == size(Y.vert,1)
                currErrorGT = max(normv(vertCurrOrig-Y.vert));
                errorArrayGT(iSurr) = currErrorGT;
            end
        end

        iMin = find(errorArray == min(errorArray),1);

        X = Xsaves{iMin};
    end

    function featOut = computeCorrCurr(XSmooth,YSmooth,vertCurrFull,paramEval,featOut)
        [assignment,~,featOut] = computecorrespondences(XSmooth,YSmooth,vertCurrFull(XSmooth.samples,:),vertCurrFull,paramEval,featOut);
        featOut.assignment = assignment;
    end
end











