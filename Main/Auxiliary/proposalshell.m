function a = proposalshell(X,k,aCurr)
    % proposalshell - sample a new shell from the space of shells
    %
    % X: input shape
    % k: shell smoothness parameter
    % aCurr: current iterate of a
    
    if ~exist('aCurr','var')
        aCurr = zeros(k,3);
    end
    
    sigma = normv(X.evecs(:,1:k)' * X.A * X.vert + aCurr);
    sigma(1) = 0;

    a = sigma .* randn(k,3) + aCurr;

end

