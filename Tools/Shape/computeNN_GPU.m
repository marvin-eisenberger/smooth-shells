function assignment = computeNN_GPU(pointsY,pointsSubX)
    % computeNN_GPU - Nearest neighbor computation
    %
    % pointsY: reference points
    % pointsSubX: points to map

    pointsY_GPU = gpuArray(pointsY);

    pointsSubX_GPU = gpuArray(pointsSubX);

    nSubX = size(pointsSubX,1);

    numSubsamples = round(1000 .* 50000 ./ size(pointsY_GPU,1));

    subDivX = [1:numSubsamples:nSubX,nSubX+1];

    assignment = gpuArray(zeros(nSubX,1));

    for iS = 1:length(subDivX)-1
        samplesSub = subDivX(iS):subDivX(iS+1)-1;


        Dcurr = sum(pointsSubX_GPU(samplesSub,:).^2,2)' + sum(pointsY_GPU.^2,2);
        Dcurr = Dcurr - 2 * pointsY_GPU * pointsSubX_GPU(samplesSub,:)';

        while true
            try
                [assignment(samplesSub),~] = find(Dcurr==min(Dcurr));
                break;
            catch
                Dcurr = Dcurr + randn(size(Dcurr)) .* mean(mean(Dcurr)) .* 1e-4;
            end
        end
    end
    assignment = gather(assignment);
        
end

    









