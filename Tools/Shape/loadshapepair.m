function [X,Y]= loadshapepair(input1,input2)
    % loadshapepair - load a pair of input shapes and compute the features and eigenpairs
    % input1,input2: either a mesh (as a struct) OR a file containing a mesh

    %% load mesh
    
    X = loadshape(input1);
    Y = loadshape(input2);
   
    [X,Y] = basisfeatures(X,Y);

end

