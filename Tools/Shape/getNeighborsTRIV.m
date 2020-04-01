function neigh = getNeighborsTRIV(TRIV,numPoints)
    
    TRIV = [TRIV(:,1),TRIV(:,2);TRIV(:,1),TRIV(:,3);TRIV(:,2),TRIV(:,3)];
    TRIV = [TRIV;TRIV(:,2),TRIV(:,1)];
    
    neigh = (sparse(TRIV(:,1)',TRIV(:,2)',ones(1,size(TRIV,1)),numPoints,numPoints)~=0);

end

