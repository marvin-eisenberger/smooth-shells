function xi = computeAxisICP(vert,vertRef,numIt,rotAxis,midPoint)
    % computeAxisICP - ICP around an axis for two pointclouds
    %
    % vert: vertices to transform
    % vertRef: reference vertices
    % numIt: number of iterations
    % rotAxis: the axis around which vert should be rotated
    % midPoint: center of mass of vert
    
    if ~exist('midPoint','var')
        midPoint = mean(vert,1);
    end

    if ~iscolumn(rotAxis)
        rotAxis = rotAxis';
    end

    rotAxis = rotAxis ./ norm(rotAxis);

    numSubsteps = 10;

    N = size(vert,1);

    xi = zeros(1,6);

    vert = vert - midPoint;
    vertRef = vertRef- midPoint;


    for it = 1:numIt*numSubsteps

        vertRot = rigidTransform(vert',xi')';

        searcher = KDTreeSearcher(vertRef);
        searcherinv = KDTreeSearcher(vertRot);

        assignment = knnsearch(searcher,vertRot);
        assignmentinv = knnsearch(searcherinv,vertRef);



        r = [vertRot - vertRef(assignment,:);vertRot(assignmentinv,:)-vertRef];

        T = se3_SE3(xi');

        r = reshape(T(1:3,1:3)' * reshape(r',3,size(r,1)),3*size(r,1),1);

        J = -hatOp(vert);
        Jinv = -hatOp(vert(assignmentinv,:));
        J = J * rotAxis;
        Jinv = Jinv * rotAxis;
        M = [J;Jinv]' * [J;Jinv];
        t = [J;Jinv]' * r;

        xi = groupMult(xi,[0,0,0,rotAxis' * (-M\t)']);

    end

    RTranslate = [eye(3),midPoint';zeros(1,3),1];
    RTranslateInv = [eye(3),-midPoint';zeros(1,3),1];

    R = se3_SE3_homCoords(xi);

    R = RTranslate * R * RTranslateInv;

    xi = transform_SE3_se3_homCoords(R);

end
