function J = jacobianrigid(vert)
    % jacobianrigid - derivative of a SE3 lie group element applied to some vertex 
    %
    % vert: vertices to which the transformation is applied
    
    J = [repmat(eye(3),size(vert,1),1),-hatOp(vert)];

end

