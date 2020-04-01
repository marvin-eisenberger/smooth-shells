function xi = computeICP(vert,vertRef,W,numIt)
    % computeAxisICP - ICP for known correspondences
    %
    % vert: vertices to transform
    % vertRef: reference vertices
    % W: correspondence matrix
    % numIt: number of iterations


    vert = vert';
    vertRef = vertRef';

    numSubsteps = 10;

    N = size(vert,1);

    xi = zeros(6,1);



    %compute M from W
    J = jacobianrigid(vert);
    Wcols = sqrt(sum(W,2));
    Wcols = kron(Wcols,[1;1;1]);
    JW = J .* Wcols;
    M = JW' * JW;

    for it = 1:numIt*numSubsteps

        
        vertRot = rigidTransform(vert',xi)';
        
        rW = zeros(3*N,1);
        
        D = zeros(size(W));
        
        for d = 1:3
            Dcurr = vertRot(:,d)-vertRef(:,d)';
            D = D + Dcurr.^2;
        end
        
        clear('Dcurr')
        
        for d = 1:3
            Dcurr = vertRot(:,d)-vertRef(:,d)';
            rW(d:3:3*N) = sum(W.*Dcurr,2);
        end
        
        T = se3_SE3(xi);
        
        rW = reshape(T(1:3,1:3)' * reshape(rW,3,N),3*N,1);
        
        t = J' * rW;
        
        xi = groupMult(xi',(-M\t)')';

    end

end

