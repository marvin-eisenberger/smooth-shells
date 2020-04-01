function [ S ] = fps_euclidean(V, n, seed)
    % fps_euclidean - Samples K vertices from V by using farthest point sampling.


    if n > size(V, 1)
        n = size(V, 1);
    end

    S = zeros(n,1);
    S(1) = seed;
    d = pdist2(V,V(seed,:));

    for i=2:n
        [~,m] = max(d);
        S(i) = m(1);
        d = min(pdist2(V,V(S(i),:)),d);
    end

end

