function Rw = se3_SE3_homCoords(xi)
    % se3_SE3_homCoords - SE3 exponential map
    % 
    % xi: se3 lie algebra element
    
    v = [xi(:,1:3)';zeros(1,size(xi,1))];
    w = xi(:,4:6);
    
    
    %% Rotation component
    
    normW = normv(w);
    L = normW==0;
    normW(L) = 1e-10;
    
    Rw = rotationblock(w,normW);
    
    %% Translation component
    
    Jleft = leftjacobianrot(w,Rw,normW);
    
    Tw = Jleft * v(:);
    Tw(4:4:end) = [];
    Tw(kron(L,ones(3,1))==1) = v(1:3,L);
    
    Tw = reshape(Tw,3,length(Tw)/3);

    Rw = Rw + translationtoblocks(Tw,size(xi,1));
end