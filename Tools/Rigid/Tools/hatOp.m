function [vHat] = hatOp(v)
    % hatOp - hat operator for some vectors v

    vHat = kron(v(:,1),[0,0,0;0,0,-1;0,1,0]) + kron(v(:,2),[0,0,1;0,0,0;-1,0,0]) + kron(v(:,3),[0,-1,0;1,0,0;0,0,0]);

end

