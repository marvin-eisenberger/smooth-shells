function [M,z] = jacobiannormalrotarap(vert,vertNew,normalNew,w,neigh,samples)
    % jacobiannormalrotarap - compute the arap jacobian for a fixed normal rotation
    %
    % vert: original vertices
    % vertNew: reference vertices
    % normalNew: morphed normals
    % w: SO3 lie algebra element of the normal rotation
    % neigh: connectivity information for vert
    % samples: samples of vert to compute the jacobian for

    neighRow = neigh.row;
    neighCol = neigh.col;
    
    normW = normv(w);
    L = normW<1e-10;
    normW(L) = 1e-10;
    
    nVertDiff = length(neighRow);
    nVert = size(vert,1);
    nSamples = length(samples);
    
    RwRow = rotationblock(w(neighRow,:),normW(neighRow),false);
     
    vertRow = vert(samples(neighRow),:);
    vertNewRow = vertNew(samples(neighRow),:);
    vertCol = vert(neighCol,:);
    vertNewCol = vertNew(neighCol,:);
    
    vertRow = vertRow';
    vertRow = vertRow(:);
    vertCol = vertCol';
    vertCol = vertCol(:);
    
    vertNewRow = vertNewRow';
    vertNewRow = vertNewRow(:);
    vertNewCol = vertNewCol';
    vertNewCol = vertNewCol(:);
    
    rotJac = RwRow*(vertCol-vertRow);
    rotJacHat = hatOpDiag(reshape(rotJac,3,nVertDiff)',false);
    

    normalNewBlock = constructblockmat(normalNew(neighRow,:),1:nVertDiff,1:nVertDiff,[nVertDiff nVertDiff])';
    rotJacHat = rotJacHat * normalNewBlock;
    
    
    z = blockstackmat(constructblockmat(rotJacHat' * (rotJac-(vertNewCol-vertNewRow)),neighRow,neighCol,[nSamples nVert]),nVert,1);
    
    M = constructblockmat(blockstackmat(constructblockmat(blockstackmat(rotJacHat' * rotJacHat,nVertDiff,1),neighRow,neighCol,[nSamples nVert]),nVert,1),1:nSamples,1:nSamples,[nSamples nSamples]);

end
