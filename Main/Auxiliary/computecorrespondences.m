function [assignment,assignmentinv,feat] = computecorrespondences(XSmooth,YSmooth,vertCurr,vertCurrFull,param,feat)
    % computecorrespondences - recover the pointwise correspondences
    % the difference between the two modes is that for param.mode = 1 the
    % assignment is computed from subvertices to full vertices and for
    % param.mode = 2 from subvertices to subvertices
    %
    % XSmooth,YSmooth: input shapes (not necessarily shells, but most of
    % the time)
    % vertCurr/vertCurrFull: current morphed shell position of
    % subsampled/full shape
    % param: parameter values
    % feat: dynamic exchange collection

    if ~isfield(feat,'iSub')
        feat.iSub = -1;
    end

    switch param.mode
        
        
    case{1}

        if isfield(feat,'C')


            overhang = length(feat.oldBasisWeights) - feat.k;

            kStart = 1;
            kCurr = feat.k;



            indexCheck = kStart:kCurr;



            basisSamplesX = XSmooth.evecs(:,1:feat.k+overhang) * feat.C(1:feat.k+overhang,indexCheck) ;
            basisSamplesY = YSmooth.evecs(:,indexCheck);

            basisSamplesSubX = basisSamplesX(XSmooth.samples,:);
            basisSamplesSubY = basisSamplesY(YSmooth.samples,:);


            facFeatX = param.facFeat .* norm(vertCurrFull,'F') ./ norm(basisSamplesX,'F') .* feat.oldBasisWeights(indexCheck)';
            facFeatY = param.facFeat .* norm(YSmooth.vert,'F') ./ norm(basisSamplesY,'F') .* feat.oldBasisWeights(indexCheck)';
        else

            basisSamplesX = [];
            basisSamplesY = [];
            basisSamplesSubX = [];
            basisSamplesSubY = [];

            facFeatX = 1;
            facFeatY = 1;
        end

        facNormalX = param.normalDamping .* norm(vertCurrFull,'F') ./ norm(feat.normalCurr,'F');
        facNormalY = param.normalDamping .* norm(YSmooth.vert,'F') ./ norm(YSmooth.normal,'F');



        pointsX = [vertCurrFull,facNormalX .* feat.normalCurr,facFeatX .* basisSamplesX];
        pointsY = [YSmooth.vert,facNormalY .* YSmooth.normal,facFeatY .* basisSamplesY];

        pointsSubX = [vertCurr,facNormalY .* feat.normalCurr(XSmooth.samples,:),facFeatY .* basisSamplesSubX];
        pointsSubY = [YSmooth.vertSub,facNormalX .* YSmooth.normal(YSmooth.samples,:),facFeatX .* basisSamplesSubY];


        if ~param.GPUcorrespondences
            searcherY = KDTreeSearcher(pointsY);
            assignment = knnsearch(searcherY,pointsSubX);


            searcherX = KDTreeSearcher(pointsX);
            assignmentinv = knnsearch(searcherX,pointsSubY);
        else
            assignment = computeNN_GPU(pointsY,pointsSubX);
            assignmentinv = computeNN_GPU(pointsX,pointsSubY);
        end


        if param.matchingAlignment

            facNormal = mean([facNormalX,facNormalY]);
            facFeat = mean([facFeatX;facFeatY],1);


            pointsWeightX = [facNormal .* feat.normalCurr,facFeat .* basisSamplesX];
            pointsWeightY = [facNormal .* YSmooth.normal,facFeat .* basisSamplesY];

            if isempty(basisSamplesX)
                feat.Dass = 0;
                feat.Dassinv = 0;
                feat.weightass = 1;
                feat.weightassinv = 1;
            else

                Dass = normv(pointsWeightX(XSmooth.samples,:)-pointsWeightY(assignment,:));
                Dassinv = normv(pointsWeightX(assignmentinv,:)-pointsWeightY(YSmooth.samples,:));

                feat.weightass = exp(-0.5 ./ mean(Dass.^2) .* Dass.^2);
                feat.weightassinv = exp(-0.5 ./ mean(Dassinv.^2) .* Dassinv.^2);

                pointsInitX = [vertCurrFull,facNormal .* feat.normalCurr];
                pointsInitY = [YSmooth.vert,facNormal .* YSmooth.normal];

                feat.Dass = normv(pointsInitX(XSmooth.samples,:)-pointsInitY(assignment,:));
                feat.Dassinv = normv(pointsInitX(assignmentinv,:)-pointsInitY(YSmooth.samples,:));
            end
        end


 %feat 
    case{2}

        if isfield(feat,'C')


            overhang = length(feat.oldBasisWeights) - feat.k;

            kStart = 1;
            kCurr = feat.k;



            indexCheck = kStart:kCurr;



            basisSamplesX = XSmooth.evecs(:,1:feat.k+overhang) * feat.C(1:feat.k+overhang,indexCheck) ;
            basisSamplesY = YSmooth.evecs(:,indexCheck);

            basisSamplesSubX = basisSamplesX(XSmooth.samples,:);
            basisSamplesSubY = basisSamplesY(YSmooth.samples,:);


            facFeatX = param.facFeat .* norm(vertCurrFull,'F') ./ norm(basisSamplesX,'F') .* feat.oldBasisWeights(indexCheck)';
            facFeatY = param.facFeat .* norm(YSmooth.vert,'F') ./ norm(basisSamplesY,'F') .* feat.oldBasisWeights(indexCheck)';
        else

            basisSamplesX = [];
            basisSamplesY = [];
            basisSamplesSubX = [];
            basisSamplesSubY = [];

            facFeatX = 1;
            facFeatY = 1;
        end

        facNormalX = param.normalDamping .* norm(vertCurrFull,'F') ./ norm(feat.normalCurr,'F');
        facNormalY = param.normalDamping .* norm(YSmooth.vert,'F') ./ norm(YSmooth.normal,'F');


        pointsSubX = [vertCurr,facNormalY .* feat.normalCurr(XSmooth.samples,:),facFeatY .* basisSamplesSubX];
        pointsSubY = [YSmooth.vertSub,facNormalX .* YSmooth.normal(YSmooth.samples,:),facFeatX .* basisSamplesSubY];


        if ~param.GPUcorrespondences
            searcherY = KDTreeSearcher(pointsSubY);
            assignment = knnsearch(searcherY,pointsSubX);


            searcherX = KDTreeSearcher(pointsSubX);
            assignmentinv = knnsearch(searcherX,pointsSubY);
        else
            assignment = computeNN_GPU(pointsSubY,pointsSubX);
            assignmentinv = computeNN_GPU(pointsSubX,pointsSubY);
        end

        assignmentinv = XSmooth.samples(assignmentinv);
        assignment = YSmooth.samples(assignment);


        if param.matchingAlignment

            facNormal = mean([facNormalX,facNormalY]);
            facFeat = mean([facFeatX;facFeatY],1);

            pointsWeightX = [facNormal .* feat.normalCurr,facFeat .* basisSamplesX];
            pointsWeightY = [facNormal .* YSmooth.normal,facFeat .* basisSamplesY];

            if isempty(basisSamplesX)
                feat.Dass = 0;
                feat.Dassinv = 0;
                feat.weightass = 1;
                feat.weightassinv = 1;
            else

                Dass = normv(pointsWeightX(XSmooth.samples,:)-pointsWeightY(assignment,:));
                Dassinv = normv(pointsWeightX(assignmentinv,:)-pointsWeightY(YSmooth.samples,:));

                feat.weightass = exp(-0.5 ./ mean(Dass.^2) .* Dass.^2);
                feat.weightassinv = exp(-0.5 ./ mean(Dassinv.^2) .* Dassinv.^2);

                pointsInitX = [vertCurrFull,facNormal .* feat.normalCurr];
                pointsInitY = [YSmooth.vert,facNormal .* YSmooth.normal];

                feat.Dass = normv(pointsInitX(XSmooth.samples,:)-pointsInitY(assignment,:));
                feat.Dassinv = normv(pointsInitX(assignmentinv,:)-pointsInitY(YSmooth.samples,:));
            end
        end
        
    end


end

    









