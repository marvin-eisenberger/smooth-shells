function [feat,C,Cgt] = computefunctionalmap(X,Y,param,feat,k,assignment,assignmentinv)
    % computefunctionalmap - compute the FM using all current information
    %
    % X,Y: input shapes
    % feat: dynamic exchange collection
    % k: current shell size
    % assignment/assignmentinv: current pointwise matchings

    %if correspondences are given, use them to get a more accurate FM
    if exist('assignment','var')

        samplesX = [X.samples;assignmentinv];
        samplesY = [assignment;Y.samples];
        
        numSamples = length(samplesX);
        
        
        
        kOverhang = length(feat.basisWeights)-k;
        
        Xfct = X.evecs(samplesX,1:k+kOverhang);
        Yfct = Y.evecs(samplesY,1:k+kOverhang);

        weightAX = diag(X.A);
        weightAX = sqrt(weightAX);
        
        weightAX = weightAX(samplesX);
        weightAX = spdiags(weightAX,0,numSamples,numSamples);
        
        
        Xfct = weightAX * Xfct;
        Yfct = weightAX * Yfct;
        
        M = Xfct' * Xfct;
        z = Xfct' * Yfct;
    else
        M = 0;
        z = 0;
        
        param.lambdaFeat = 1;
    end

    kCurr = length(feat.basisWeights);
    overhang = length(feat.basisWeights)-k;
    
    %incorporate descriptor information for the FM 
    if isfield(X,'basisfeatures')
        numEvecs = size(X.evecs,2);
        Mfeat = X.basisfeatures(1:numEvecs,:) * X.basisfeatures(1:numEvecs,:)';
        zfeat = X.basisfeatures(1:numEvecs,:) * Y.basisfeatures(1:numEvecs,:)';

        Mfeat = Mfeat(1:kCurr,1:kCurr);
        zfeat = zfeat(1:kCurr,1:kCurr);

        facNorm = norm(M) ./ norm(Mfeat);

        if facNorm == 0
            facNorm = 1;
        end

        M = M + facNorm .* param.lambdaFeat .* Mfeat;
        z = z + facNorm .* param.lambdaFeat .* zfeat;
    end

    %add a spectrum regularizer/preservation of the Laplacian
    if param.lambdaLap > 0
        L = sqrt(abs(X.evals(1:kCurr) - Y.evals(1:kCurr)'));
        L = param.lambdaLap ./ mean(L) .* mean(mean(abs(M))) .* L;
        L = spdiags(L(:),0,kCurr.^2,kCurr.^2);
        M = kron(speye(kCurr),M);
        z = z(:);

        M = M + L;

        C = M\z;
        C = full(reshape(C,kCurr,kCurr));
    else
        C = full(M\z);
    end

    %compute the GT FM (for reference/plotting of the output) if possible
    if size(X.vert,1) == size(Y.vert,1)
        Cgt = (X.evecs(:,1:kCurr)' * X.A * Y.evecs(:,1:kCurr)); 
        Cgt = Cgt(1:k+overhang,1:k+overhang);
    else
        Cgt = 0;
    end

    C = C(1:k+overhang,1:k+overhang);

    %project the FM on the space of orthogonal matrices
    if param.areaPreservation
        C = orthogonalizematrix(C,param.GPUcorrespondences);
    end

    %if some prior FM is given (e.g from the initialization), replace the
    %respective entries
    if isfield(feat,'Cfix')
        C(1:size(feat.Cfix,1),:) = 0;
        C(:,1:size(feat.Cfix,2)) = 0;

        C(1:size(feat.Cfix,1),1:size(feat.Cfix,2)) = feat.Cfix;

        if param.areaPreservation
            C(size(feat.Cfix,1)+1:end,size(feat.Cfix,2)+1:end) = orthogonalizematrix(C(size(feat.Cfix,2)+1:end,size(feat.Cfix,2)+1:end),param.GPUcorrespondences);
        end
    end


    feat.C = C;
    feat.k = k;
    feat.oldBasisWeights = feat.basisWeights;
        
end
        
        
        
        
        
        