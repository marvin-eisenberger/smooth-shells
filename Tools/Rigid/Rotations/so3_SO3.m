function Rw = so3_SO3(w,normW,homCoords)

    if ~exist('normW','var')
        normW = normv(w);
        normW(normW<=1e-10)=1e-10;
    end
    
    if ~exist('homCoords','var')
        homCoords = false;
    end
    
    wHat = hatOpDiag(w ./ normW .* sin(normW),homCoords);
    wHatSq = hatOpDiag((w ./ normW) .* sqrt(1-cos(normW)),homCoords)^2;
    Rw = speye(size(wHat)) + wHat + wHatSq;
end