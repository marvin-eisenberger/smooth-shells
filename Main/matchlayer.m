
function [X,Y,tau,feat,avgError,assignment,assignmentinv] = matchlayer(X,Y,feat,tau,k,param)
    % matchlayer - substep of the main method for one shell pair
    %
    % X,Y: input shape
    % feat: dynamic exchange collection
    % a: current deformation coefficients
    % k: current shell size
    % param: parameter values


    [XSmooth,feat.basisWeights] = smoothshapesigmoid(X,k);
    [YSmooth,~] = smoothshapesigmoid(Y,k);

    avgError = 0;

    JleftCurr = speye(X.n*3);

    for iSub = 1:param.numSub

        %compute morphed shape from previous iteration
        vertCurrFull = XSmooth.vert + rotatevertices(JleftCurr,X.evecs(:,1:k) * tau(1:k,:));
        vertCurr = vertCurrFull(X.samples,:);

        feat.normalCurr = compute_normal(vertCurrFull,X.triv',X.flipNormal)';

        %compute correspondences
        feat.iSub = iSub;
        [assignment,assignmentinv,feat] = computecorrespondences(XSmooth,YSmooth,vertCurr,vertCurrFull,param,feat);

        feat.assignment = assignment;
        feat.assignmentinv = assignmentinv;

        %compute the linear system for the data energy term
        
        M = X.evecs(X.samples,1:k)' * X.evecs(X.samples,1:k) + ...
            X.evecs(assignmentinv,1:k)' * X.evecs(assignmentinv,1:k);
        z = -X.evecs(X.samples,1:k)' * (XSmooth.vertSub - YSmooth.vert(assignment,:)) ...
            -X.evecs(assignmentinv,1:k)' * (XSmooth.vert(assignmentinv,:) - YSmooth.vertSub);



        %compute the Gauss Newton system for the arap regularization term
        if param.lambdaArap > 0
            
            feat.RCurr = determinerotationsshifted(XSmooth,vertCurrFull);
            feat.wCurr = transform_SO3_so3(feat.RCurr);

            [Marap,zarap] = jacobianshiftarap(XSmooth.vert,feat.wCurr,X.evecs(:,1:k),X.neigh,(1:X.n)');

            lambdaArapCurr = param.lambdaArap * norm(M,inf) ./ norm(Marap,inf);

            M = M + lambdaArapCurr .* Marap;
            z = z + lambdaArapCurr .* zarap;

        end

        %solve the linear system
        tau(1:k,:) = M \ z;


        %compute the new functional map
        [feat,C,Cgt] = computefunctionalmap(X,Y,param,feat,k,assignment,assignmentinv);

        %plot intermediate results (if settings demand it)
        if ~param.noPlotInBetween
            avgError = mean(normv(Y.vert(X.samples,:)-Y.vert(assignment,:)));
            plotskeletonlayered(X,Y,XSmooth,YSmooth,param,assignment,assignmentinv,tau,avgError,k,C,Cgt);

        end
    end

    %output error, if gt is available
    if size(vertCurrFull,1)==size(YSmooth.vert,1) && param.intermediateOutput
        disp('max error: ' + string(max(normv(vertCurrFull-YSmooth.vert))) + ', mean error: ' + string(mean(normv(vertCurrFull-YSmooth.vert))))
    end

    

    %plot intermediate results (if settings demand it)
    if ~param.noPlot

        plotskeletonlayered(X,Y,XSmooth,YSmooth,param,assignment,assignmentinv,tau,avgError,k,C,Cgt);

    end
end