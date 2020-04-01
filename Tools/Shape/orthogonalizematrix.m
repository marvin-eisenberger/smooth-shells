function C = orthogonalizematrix(C,useGPU)
    % orthogonalizematrix - map matrix on space of orthogonal matrices

    if useGPU
        C = gpuArray(C);
    end
    
    [orthleft,S,orthright] = svd(C);
    S = diag(S);
    S(:) = 1;
    S = diag(S);
    C = orthleft * S * orthright';
    
    if useGPU
        C = gather(C);
    end

end
