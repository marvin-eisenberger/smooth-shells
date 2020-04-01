function [M,z] = jacobianshiftarap(vert,w,V,neigh,samples)
    % jacobianshiftarap - compute the arap regularization linear system
    %
    % vert: original vertices
    % w: SO3 lie algebra element of the normal rotation
    % V: shifted vertices
    % neigh: connectivity information for vert
    % samples: samples of vert to compute the regularizer for

    
    nVert = size(vert,1);
    nVertDiff = length(neigh.row);
    nSamples = length(samples);
     
    vertRow = vert(samples(neigh.row),:);
    vertCol = vert(neigh.col,:);


    Vdiff = V(samples(neigh.row),:) - V(neigh.col,:);

    

    w = w(samples(neigh.row),:);

    normW = normv(w); 
    L = normW<1e-10;
    normW(L) = 1e-10;

    Rw = rotationblock(w,normW,false);
    
    vertDiff = vertRow-vertCol;
    vertDiff = rotatevertices(Rw-speye(size(Rw)),vertDiff);
    
    
    
    z = Vdiff'*vertDiff;
    
    M = Vdiff' * Vdiff;

    

end
