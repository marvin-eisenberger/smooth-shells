function [X,ainit] = rotateainit(X,ainit)
    % rotateAInit - rigid transformation of X according to ainit
    % 
    % X: input shape
    % ainit: determined initial pose

    vertCurrOrig = shiftedvertupsample(X,ainit);
    kinit = size(ainit,1);
    XSmoothInit = smoothshapesigmoid(X,kinit);
    xiInit = XSmoothInit.xi;
    clear('XSmoothInit')
    R = determinerotationsshifted(X,vertCurrOrig);
    w = transform_SO3_so3(R);
    w = mean(w,1);
    R = so3_SO3(w);
    R = full(R);
    X.vert = X.vert*R';
    X.vertSub = X.vert(X.samples,:);
    ainit(1:kinit,:) = ainit + xiInit(1:kinit,:) * (eye(3)-R');
    
end
