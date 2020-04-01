function wCurr = determinerotationsshiftedshape(X,vertNew,neigh,samples,numIt,wCurr,normalNew)
    % determinerotationsshiftedshape - compute rigid transformation for
    % each vertex in X that describes the mapping to vertNew best
    %
    % X: input shape
    % vertNew: mapped shape coordinates
    % neigh: neighorhood/connectivity information
    % samples: subsamples of the shape for which the rotations are computed
    % numIt: number of refinement iterations
    % wCurr: old rotations from previous iteration
    % normalNew: outer normals of vertNew
    

    if ~exist('samples','var')
        samples = 1:X.n;
    end

    if ~isfield(neigh,'row')
        [neigh.row,neigh.col] = find(neigh.mat(samples,:));
    end

    if ~exist('numIt','var')
    	numIt = 50;
    end

    if ~exist('wCurr','var')
    	wCurr = zeros(size(X.vert(samples,:)));
    end



    normWCurr = normv(wCurr);

    Rw = rotationblock(wCurr,normWCurr,false);

    normalXrot = rotatevertices(Rw,X.normal(samples,:));

    [~,wNormals] = rotvectorpairs(normalXrot,normalNew,false);


    xiCurr = [zeros(size(wCurr)),wCurr];
    xiNormals = [zeros(size(wNormals)),wNormals];

    xiCurr = groupMult(xiNormals,xiCurr);

    wCurr = xiCurr(:,4:6);


    resHist = zeros(numIt,1);

    for it = 1:numIt

        [M,z] = jacobiannormalrotarap(X.vert,vertNew,normalNew,wCurr,neigh,samples);
        

        wDeltaFac = M\z;

        normalNewBlock = constructblockmat(normalNew,1:length(samples),1:length(samples),[length(samples) length(samples)])';
        wDelta = normalNewBlock * wDeltaFac;

        wDelta = reshape(wDelta,3,length(samples))';

        xiCurr = [zeros(size(wCurr)),wCurr];
        xiDelta = [zeros(size(wDelta)),wDelta];

        xiCurr = groupMult(xiDelta,xiCurr);

        wCurr = xiCurr(:,4:6);

        resHist(it) = mean(normv(wDelta));

        if it>1
            relError = abs(resHist(it)-resHist(it-1)) ./ mean(normv(wCurr));

            if relError < 0.01
                break;
            end
        end


    end

    wCurr = full(wCurr);
end