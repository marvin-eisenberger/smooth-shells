function vert = rotatevertices(Rw,vert)
    % rotatevertices - apply rotation matrix to vertices
    % 
    % Rw: rotation matrix
    % vert: input vertices

    vert = reshape(Rw * reshape(vert',3*size(vert,1),1),3,size(vert,1))';

end

