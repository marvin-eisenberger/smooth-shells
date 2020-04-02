function [vertOrigFull,w,R] = shiftedvertupsample(X,a)
    %shiftedvertupsample - determine vertOrig from vertCurr
    %
    % X: input shape
    % a: shifting coefficients
    
    iLast = find(any(a~=0,2),1,'last');
    if isempty(iLast)
        iLast = 1;
    end
    
    a = a(1:iLast,:);
    
    k = size(a,1);
    
    XSmooth = smoothshapesigmoid(X,k);
    vertShifted = XSmooth.vert + X.evecs(:,1:size(a,1)) * a;
    
    R = determinerotationsshifted(XSmooth,vertShifted);

    w = real(transform_SO3_so3(R));

    R = so3_SO3(w);

    vertOrigFull = reshape(R * reshape((X.vert-XSmooth.vert)',3*X.n,1),3,X.n)' + XSmooth.vert + 1.*(vertShifted-XSmooth.vert);
end









