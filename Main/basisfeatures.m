function [X,Y] = basisfeatures(X,Y)
    % basisfeatures - compute the features for shapes X and Y, namely the
    % SHOT and HKS descriptors
    %
    % X,Y: input shapes

    
    rangeHKS = [0.002,0.05];


    %% SHOT
   
    Xnormalize = X;
    Ynormalize = Y;
    refarea = 1.93e+04;
    Xnormalize.vert = Xnormalize.vert ./ sqrt(sum(calc_tri_areas(Xnormalize))) .* sqrt(refarea);
    Ynormalize.vert = Ynormalize.vert ./ sqrt(sum(calc_tri_areas(Ynormalize))) .* sqrt(refarea);


    Xnormalize.area = sum(calc_tri_areas(Xnormalize));
    Ynormalize.area = sum(calc_tri_areas(Ynormalize));

    opts = struct;
    opts.shot_num_bins = 10; % number of bins for shot
    opts.shot_radius = 5; % percentage of the diameter used for shot

    X.SHOT = calc_shot(Xnormalize.vert', Xnormalize.triv', 1:Xnormalize.n, opts.shot_num_bins, opts.shot_radius*sqrt(Xnormalize.area)/100, 3)';
    Y.SHOT = calc_shot(Ynormalize.vert', Ynormalize.triv', 1:Ynormalize.n, opts.shot_num_bins, opts.shot_radius*sqrt(Ynormalize.area)/100, 3)';
 
    %% HKS

    tArray = exp(linspace(log(rangeHKS(1)),log(rangeHKS(2)),100));
    X.HKS = calc_HKS(X,200,tArray);
    Y.HKS = calc_HKS(Y,200,tArray);

    %%



    XbasisSHOT = X.evecs' * X.A * X.SHOT;
    YbasisSHOT = Y.evecs' * Y.A * Y.SHOT;

    normFacSHOT = sqrt(sum(sum(XbasisSHOT.^2)) + sum(sum(YbasisSHOT.^2)));

    XbasisSHOT = XbasisSHOT ./ normFacSHOT;
    YbasisSHOT = YbasisSHOT ./ normFacSHOT;

    XbasisHKS = X.evecs' * X.A * X.HKS;
    YbasisHKS = Y.evecs' * Y.A * Y.HKS;

    normFacHKS = sqrt(sum(sum(XbasisHKS.^2)) + sum(sum(YbasisHKS.^2)));

    XbasisHKS = XbasisHKS ./ normFacHKS;
    YbasisHKS = YbasisHKS ./ normFacHKS;
  

    
    X.basisfeatures = [XbasisSHOT,XbasisHKS];
    Y.basisfeatures = [YbasisSHOT,YbasisHKS];



    %normalize
    normFac = mean([sqrt(sum(X.basisfeatures.^2,2)),sqrt(sum(Y.basisfeatures.^2,2))],2);
    X.basisfeatures = X.basisfeatures ./ normFac;
    Y.basisfeatures = Y.basisfeatures ./ normFac;
end
