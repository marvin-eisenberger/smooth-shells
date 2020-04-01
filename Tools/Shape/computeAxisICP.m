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














% % % % function xi = computeDistICP(vert,vertRef,numIt)
% % % % 
% % % % vert = vert';
% % % % vertRef = vertRef';
% % % % 
% % % % assignCorresp = true;
% % % % 
% % % % numSubsteps = 10;
% % % % 
% % % % N = size(vert,1);
% % % % 
% % % % % constScale = 1 ./ mean(sqrt(sum(vert.^2,2)));
% % % % % 
% % % % % vert = vert .* constScale;
% % % % % vertRef = vertRef .* constScale;
% % % % 
% % % % J = jacobianRigid(vert);
% % % % 
% % % % xi = zeros(6,1);
% % % % 
% % % % quadRadius = 1;
% % % % sigma = 0.1;
% % % % 
% % % % W = zeros(size(vert,1));
% % % % 
% % % % for it = 1:numIt*numSubsteps
% % % % 
% % % % %     T = se3_SE3(xi);
% % % %     
% % % %     vertRot = rigidTransform(vert',xi)';
% % % %     
% % % %     rW = zeros(3*N,1);
% % % %     
% % % %     D = zeros(size(W));
% % % %     
% % % %     for d = 1:3
% % % %         Dcurr = vertRot(:,d)-vertRef(:,d)';
% % % %         D = D + Dcurr.^2;
% % % %     end
% % % %     
% % % %     clear('Dcurr')
% % % %     
% % % %     potRobust = 2;
% % % %     
% % % %     if mod(it,numSubsteps)==1
% % % % 
% % % %     %     W = (D'==min(D'));
% % % %         assignment = assignmentAlgs(norm(sqrt(D),inf)-sqrt(D) + 1e-10,'auction');
% % % %         W(:) = 0;
% % % %         for iAss = 1:length(assignment)
% % % %             W(iAss,assignment(iAss))=1;
% % % %         end
% % % % 
% % % %         
% % % % 
% % % % %         sigma = 0.01;
% % % % 
% % % % 
% % % %         huberLin = D > sigma.^2;
% % % %         invDWeights = 1 ./ (D ./ sigma.^2).^potRobust;
% % % %         invDWeights(~huberLin) = 1;
% % % %         
% % % % %         invDWeights = 1-huberLin;
% % % % 
% % % %         W = W .* invDWeights;
% % % %         
% % % %         quadRadius = sqrt(sum(D.*W,2) ./ (sum(W,2)+1e-3)) / 1.0;
% % % %         sigma = mean(quadRadius);
% % % % 
% % % %         clear('invDWeights')
% % % % 
% % % %         %determine M from W
% % % %         Wcols = sqrt(sum(W,2));
% % % %         Wcols = kron(Wcols,[1;1;1]);
% % % %         JW = J .* Wcols;
% % % %     %     M = JWevec' * JWevec + lambda*eye(6*K);
% % % %         M = JW' * JW;
% % % %     end
% % % %     
% % % %     
% % % %     for d = 1:3
% % % %         Dcurr = vertRot(:,d)-vertRef(:,d)';
% % % %         rW(d:3:3*N) = sum(W.*Dcurr,2);
% % % %     end
% % % %     
% % % %     T = se3_SE3(xi);
% % % %     
% % % %     rW = reshape(T(1:3,1:3)' * reshape(rW,3,N),3*N,1);
% % % %     
% % % %     t = J' * rW;
% % % %     
% % % %     xi = xi - M\t;
% % % %     
% % % % 
% % % % end
% % % % 
% % % % end
% % % % 
