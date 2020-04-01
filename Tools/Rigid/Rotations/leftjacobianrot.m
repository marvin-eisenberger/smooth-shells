function Jleft = leftjacobianrot(w,Rw,normW,homCoords)
    % leftjacobianrot - left jacobian of a SO3 lie group transformation 
    %
    % w: so3 element
    % Rw: SO3 element
    % normW: norm of vectors w
    % homCoords: output with homogeneous coordinates (4x4) or as a pure
    % rotation matrix (3x3)


    if nargin < 4
       homCoords = true; 
    end

    if nargin < 3
       normW = normv(w); 
    end
    
    L = normW==0;
    normW(L) = 1e-10;
    
    L = kron(L,ones(3,1))==1;

    wHat = hatOpDiag(w ./ normW.^2,homCoords);
    
    if homCoords
        wCell = mat2cell([sparse(w ./ normW)';sparse(1,size(w,1))],4,ones(size(w,1),1));
    else
        wCell = mat2cell(sparse(w ./ normW)',3,ones(size(w,1),1));
    end
    wMat = blkdiag(wCell{:});
    
    Jleft = ((speye(size(wHat)) - Rw) * wHat + wMat * wMat');
    
    Jleft(L,L) = speye(sum(L));
end

