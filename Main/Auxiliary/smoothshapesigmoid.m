function [X,weights] = smoothshapesigmoid(X,k,truncthresh,smoothThresh)
    % smoothshapesigmoid - convert shape X to its smooth shell
    %
    % X: input shape
    % k: current shell size
    % truncthresh/smoothThresh: parameters resulting from the sigmoid 
    % truncation limit and sigma from the paper
    
    if ~exist('truncthresh','var')
        truncthresh = 1e-2;
    end
    
    if ~exist('smoothThresh','var')
        smoothThresh = 10;
    end
    
    %compute the sigmoid weights
    if k == 1
        weights = 1;
    elseif k >= smoothThresh
        t = -1./(smoothThresh-1) .* log(1/(1-truncthresh)-1);
        weights = 1 ./ (1+exp(t.*((-(k-1):-1+smoothThresh)')));
    else
        t = -1./(k-1) .* log(1/(1-truncthresh)-1);
        weights = 1 ./ (1+exp(t.*((-(k-1):k-1)')));
    end

    k = length(weights);

    %project the geometry on the eigenfunctions
    X.xi = (weights .* (X.evecs(:,1:k)' * X.A * X.vert));
    X.vert = X.evecs(:,1:k) * X.xi;

    if isfield(X,'normal')
        X.normal = compute_normal(X.vert',X.triv',X.flipNormal)';
    end
    
    if isfield(X,'vertSub')
        X.vertSub = X.vert(X.samples,:);
    end
    
        
end

