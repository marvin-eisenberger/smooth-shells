function R = determinerotationsshifted(X,vertNew)

if ~isfield(X,'samples')
    X.samples = (1:X.n)';
end

if ~isfield(X,'neigh')
    if X.n == size(samples,1)
        X.neigh.mat = getNeighborsTRIV(X.triv,X.n);
    else
        X.neigh.mat = getNeighbors(X.vert');
    end
end

if ~isfield(X.neigh,'row')
    [X.neigh.row,X.neigh.col] = find(X.neigh.mat(samples,:));
end

nVertDiff = length(X.neigh.row);

vertDiff = X.vert(X.neigh.col,:) - X.vert(X.neigh.row,:);
vertDiffNew = vertNew(X.neigh.col,:) - vertNew(X.neigh.row,:);

vertDiff = constructblockmat(vertDiff,1:nVertDiff,1:nVertDiff,[nVertDiff nVertDiff])';
vertDiffNew = constructblockmat(vertDiffNew,1:nVertDiff,1:nVertDiff,[nVertDiff nVertDiff])';

S = vertDiff * vertDiffNew';
S = blockstackmat(S,nVertDiff,3);
S = constructblockmat(S,X.neigh.row,X.neigh.col,[X.n X.n])';
S = blockstackmat(S,X.n,3);
% S = constructblockmat(S,(1:X.n)',(1:X.n)',[X.n X.n])';
S = full(S);
R = zeros(size(S));

Sigma = eye(3);

for iV = 1:X.n
    iCurr = iV*3-2:iV*3;
    [U,~,V] = svd(S(iCurr,:));
    Sigma(3,3) = det(V*U');
    R(iCurr,:) = V * Sigma * U';
end

R = constructblockmat(R,(1:X.n)',(1:X.n)',[X.n X.n])';

end