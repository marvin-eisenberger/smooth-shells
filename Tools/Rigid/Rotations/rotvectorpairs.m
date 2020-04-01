function [Rw,w] = rotvectorpairs(a,b,homCoords)
    % rotvectorpairs - compute direct rotation between the vectors a and b
    % 
    % a,b: input vectors
    % homCoords: output with homogeneous coordinates (4x4) or as a pure
    % rotation matrix (3x3)
    
    if nargin < 3
       homCoords = true;
    end

    a = (a ./ normv(a))';
    b = (b ./ normv(b))';

    w = (hatOpDiag(a')*b(:))';
    
    w = reshape(w,3,length(w) ./ 3)';
    
    normW = normv(w);
    
    L = normW<=1e-10;
    normW(L) = 1e-10;
    
    angle = acos(sum((a .* b),1))';
    
    w = w ./ normW .* angle;
    
    Rw = rotationblock(w,angle,homCoords);
end

