function Sblock = scalartorotweight(scalarVec,sizeM,homCoords)
    % scalartorotweight - transform scalar vector to weights for block
    % diagonal roation matrices
    % 
    % scalarVec: rotation matrix weights
    % sizeM: size of the goal matrix
    % homCoords: output with homogeneous coordinates (4x4) or as a pure
    % rotation matrix (3x3)

    if nargin < 3 || homCoords
        Sblock = kron(spdiags(scalarVec,0,sizeM,sizeM),rotationweightmatrix());
    else
        Sblock = kron(spdiags(scalarVec,0,sizeM,sizeM),ones(3));
    end
    
end

