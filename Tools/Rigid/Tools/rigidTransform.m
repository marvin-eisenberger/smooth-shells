function X1 = rigidTransform(X1,xi)
    % rigidTransform - rigid transformation of vertices X1 by a SE3 element
    %
    % X1: input coordinates
    % xi: se3 lie algebra element

    n = size(X1,2);

    X1 = [X1;ones(1,n)];

    T = se3_SE3(xi);

    X1 = T*X1;

    X1 = X1(1:3,:);

end

