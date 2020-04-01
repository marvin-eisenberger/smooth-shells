function w = SO3_so3(Rw)

    %invert rotation w
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
end

