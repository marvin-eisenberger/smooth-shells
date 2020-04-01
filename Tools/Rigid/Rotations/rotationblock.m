function Rw = rotationblock(w,normW,homCoords)
    % rotationblock - exponential map for SO3 to block diagonal matrix
    % 
    % w: lie algebra element
    % normW: norm of vectors w
    % homCoords: output with homogeneous coordinates (4x4) or as a pure
    % rotation matrix (3x3)
    
    if nargin < 3
       homCoords = true;
    end

    if nargin < 2
       normW = normv(w); 
    end
    
    L = normW==0;
    normW(L) = 1e-10;

    wHat = hatOpDiag(w ./ normW,homCoords);    
    wHatSq = wHat^2;

    sinFac = scalartorotweight(sin(normW),size(w,1),homCoords);
    cosFac = scalartorotweight(1-cos(normW),size(w,1),homCoords);
    
    Rw = speye(size(wHat)) + wHat .* sinFac + wHatSq .* cosFac;
end

