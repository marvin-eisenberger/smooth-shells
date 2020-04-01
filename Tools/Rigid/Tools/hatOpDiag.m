function vHat = hatOpDiag(v,homCoords)
    % hatOpDiag - hat operator for some vectors v, creates a block diagonal matrix
    %
    % v: vectors v
    % homCoords: output with homogeneous coordinates (4x4) or as a pure
    % rotation matrix (3x3)

    if nargin < 2 || ~homCoords
        vHat = kron(spdiags(v(:,1),0,size(v,1),size(v,1)),[0,0,0;0,0,-1;0,1,0]) + kron(spdiags(v(:,2),0,size(v,1),size(v,1)),[0,0,1;0,0,0;-1,0,0]) + kron(spdiags(v(:,3),0,size(v,1),size(v,1)),[0,-1,0;1,0,0;0,0,0]);
    else
        vHat = kron(spdiags(v(:,1),0,size(v,1),size(v,1)),extendMatrix4D([0,0,0;0,0,-1;0,1,0])) + kron(spdiags(v(:,2),0,size(v,1),size(v,1)),extendMatrix4D([0,0,1;0,0,0;-1,0,0])) + kron(spdiags(v(:,3),0,size(v,1),size(v,1)),extendMatrix4D([0,-1,0;1,0,0;0,0,0]));
    end

    function M = extendMatrix4D(M)
        M = [M,sparse(3,1);sparse(1,4)];
    end

end

