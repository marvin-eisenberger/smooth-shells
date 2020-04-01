function B = constructblockmat(blocks,neighRow,neighCol,sizeNeigh,colorDist)
    %writes all the blocks contained in block into a block matrix, 
    %where the block locations are specified by neighRow, neighCol.

    blockRow = size(blocks,1) / length(neighRow);
    blockCol = size(blocks,2);
    
    BRow = blockRow * sizeNeigh(1);
    BCol = blockCol * sizeNeigh(2);
    
    B = sparse(BRow,BCol);
    
    for iR = 1:blockRow
        selR = sparse(blockRow,1);
        selR(iR) = 1;
        selR = kron(speye(length(neighRow)),selR);
        
        for iC = 1:blockCol
            selC = sparse(blockCol,1);
            selC(iC) = 1;
            
            entriesCurr = selR' * blocks * selC;
            
            B = B + sparse((neighRow-1)*blockRow + iR,(neighCol-1)*blockCol + iC,entriesCurr,BRow,BCol);
        end
    end
    
    if nargin >= 5
        colorDist = kron(colorDist,ones(blockRow,blockCol));
        
        B = B .* colorDist;
    end
    
end

