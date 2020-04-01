function xi = SE3_se3_homCoords(T)
    % SE3_se3_homCoords - SE3 log map
    % 
    % T: SE3 transformation matrix

    %invert rotation w
    Rw = T;
    Rw(4:4:end,:) = [];
    Rw(:,4:4:end) = [];

    Rwdiag = full(diag(Rw));
    Rwdiag = reshape(Rwdiag,3,length(Rwdiag)/3);
    
    theta = acos((sum(Rwdiag)-1)/2);
    
    lnRw = theta./(2.*sin(theta));
    lnRw(abs(theta)<1e-10) = 0.5;
    lnRw = spdiags(kron(lnRw,[1,1,1])',0,size(Rw,1),size(Rw,2))*(Rw-Rw');
    
    selVec = cell(3,1); 
    for d = 1:3
        ed = sparse(3,1);
        ed(d) = 1;
        selVec{d} = kron(speye(length(theta)),ed);
    end
    
    w=[-diag(selVec{2}' * lnRw * selVec{3}) diag(selVec{1}' * lnRw * selVec{3}) -diag(selVec{1}' * lnRw * selVec{2})];

    w = full(w);
    
    %invert translation u
    
    normW = normv(w);
    L = normW<1e-10;
    normW(L) = 1e-10;
    
    selVecTw = kron(ones(length(theta),1),[0;0;0;1]);
    
    Tw = T * selVecTw;
    Tw(4:4:end) = [];
    
    L = abs(kron((sum(Rwdiag)-3)',ones(3,1)))<1e-10;
    
    Jleft = leftjacobianrot(w,Rw,normW,false);
    
    u = Jleft \ Tw;

    L = L | isnan(u) | isinf(u);
    u(L) = Tw(L);
    u = reshape(u,3,length(u)/3);
    
    xi = [u;w']';
    
end

