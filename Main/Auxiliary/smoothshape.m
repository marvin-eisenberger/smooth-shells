function X = smoothshape(X,k)
    % smoothshape - convert shape X to smoothed projected shape
    % like smoothshapesigmoid, but for the operator \mathcal{P} instead of
    % \mathcal{S} from the paper
    %
    % X: input shape
    % k: current shell size
    
    X.xi = (X.evecs(:,1:k)' * X.A * X.vert);
    X.vert = X.evecs(:,1:k) * X.xi;

    if isfield(X,'normal')
        X.normal = compute_normal(X.vert',X.triv',X.flipNormal)';
    end
    
    if isfield(X,'vertSub')
        X.vertSub = X.vert(X.samples,:);
    end
    
    
    if isfield(X,'vertSmall')
        X.vertSmall = X.vert(X.samplesSmall,:);
    end
    
    
        
end











