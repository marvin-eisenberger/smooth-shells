function plotskeletonlayered(X,Y,XSmooth,YSmooth,param,assignment,assignmentinv,tau,avgError,k,C,Cgt,normalCurr)
    % plotskeletonlayered - plot a current results from a subiteration
    %
    % X,Y,XSmooth,YSmooth: normal/smoothed input shapes
    % param: parameter values
    % assignment/assignmentinv: current matchings
    % a: morphing coefficients
    % avgError: errorEstimate
    % k: shell smoothing parameter
    % C/Cgt: current/GT FM
    % normalCurr: current outer normals of the morphed shell
  
    plotNormals = exist('normalCurr');

    functionalMap = exist('C','var');

    vertCurrFull = XSmooth.vert + X.evecs(:,1:k) * tau(1:k,:);
    vertCurrOrig = X.vert + X.evecs(:,1:k) * tau(1:k,:);

    if param.twoPlots
        figure(1)
    end
    
    if functionalMap
        subplot(2,2,1)
    else
        subplot(1,2,1)
    end
    hold off
    if size(vertCurrFull,1)==size(YSmooth.vert,1)
        trisurf(XSmooth.triv,vertCurrFull(:,1),vertCurrFull(:,2),vertCurrFull(:,3),normv(vertCurrFull-YSmooth.vert));
    else
        trisurf(XSmooth.triv,vertCurrFull(:,1),vertCurrFull(:,2),vertCurrFull(:,3));
    end
    
    if (plotNormals)
        hold on
        quiver3(vertCurrFull(X.samples,1),vertCurrFull(X.samples,2),vertCurrFull(X.samples,3),normalCurr(X.samples,1),normalCurr(X.samples,2),normalCurr(X.samples,3),2,'Color','black');
    end
    
    axis equal
    title('k = ' + string(k));
    shading interp

    hold on
    if size(vertCurrFull,1)==size(YSmooth.vert,1)
        trisurf(YSmooth.triv,YSmooth.vert(:,1),YSmooth.vert(:,2),YSmooth.vert(:,3),normv(vertCurrFull-YSmooth.vert));
    else
        trisurf(YSmooth.triv,YSmooth.vert(:,1),YSmooth.vert(:,2),YSmooth.vert(:,3));
    end
    axis equal
    title('k = ' + string(k));
    colorbar
    shading interp
    
    if param.plotCorr
        vertConnections = cell(3,1);
        vertConnectionsinv = cell(3,1);

        if length(X.samples) > 3000
            samplesPlot = randi(length(X.samples),1000,1);
            X.samples = X.samples(samplesPlot);
            Y.samples = Y.samples(samplesPlot);
            assignment = assignment(samplesPlot);
            assignmentinv = assignmentinv(samplesPlot);
        end
        
        for d = 1:3
            vertConnections{d} = [vertCurrFull(X.samples,d),YSmooth.vert(assignment,d)];
            vertConnectionsinv{d} = [vertCurrFull(assignmentinv,d),YSmooth.vert(Y.samples,d)];
        end

        plot3(vertConnections{1}',vertConnections{2}',vertConnections{3}','r')
        plot3(vertConnectionsinv{1}',vertConnectionsinv{2}',vertConnectionsinv{3}','g')
    end

    if functionalMap
        subplot(2,2,2)
    else
        subplot(1,2,2)
    end
    if param.showOff
        hold off
        if size(vertCurrOrig,1)==size(Y.vert,1)
            trisurf(X.triv,vertCurrOrig(:,1),vertCurrOrig(:,2),vertCurrOrig(:,3),normv(vertCurrOrig-Y.vert));
        else
            trisurf(X.triv,vertCurrOrig(:,1),vertCurrOrig(:,2),vertCurrOrig(:,3));
        end
        axis equal
        title('k = ' + string(k));
        colorbar
    else
        hold off
        trisurf(XSmooth.triv,XSmooth.vert(:,1),XSmooth.vert(:,2),XSmooth.vert(:,3),normv(vertCurrOrig-Y.vert));
        axis equal
        title('k = ' + string(k));
        shading interp

        hold on
        trisurf(YSmooth.triv,YSmooth.vert(:,1),YSmooth.vert(:,2),YSmooth.vert(:,3),normv(vertCurrOrig-Y.vert));
        axis equal
        title('k = ' + string(k));
        colorbar
        shading interp

        if param.plotCorr

            vertConnections = cell(3,1);
            vertConnectionsinv = cell(3,1);

            for d = 1:3
                vertConnections{d} = [XSmooth.vert(X.samples,d),YSmooth.vert(assignment,d)];
                vertConnectionsinv{d} = [XSmooth.vert(assignmentinv,d),YSmooth.vert(Y.samples,d)];
            end

            plot3(vertConnections{1}',vertConnections{2}',vertConnections{3}','r')
            plot3(vertConnectionsinv{1}',vertConnectionsinv{2}',vertConnectionsinv{3}','g')
        end
    end


    if param.twoPlots
    figure(2)

    subplot(1,3,1)
    hold off
    trisurf(XSmooth.triv,XSmooth.vert(:,1),XSmooth.vert(:,2),XSmooth.vert(:,3),normv(X.vert-XSmooth.vert));
    hold on
    plot3(X.vertSub(:,1),X.vertSub(:,2),X.vertSub(:,3),'xb');
    axis equal
    title('k = ' + string(k));


    subplot(1,3,2)
    hold off
    plot(assignment,'xb')
    hold on
    plot(X.samples,'xr')
    grid on
    title('avg error = ' + string(avgError));

    subplot(1,3,3)
    Dpoints = sort(normv(Y.vert(X.samples,:)-Y.vert(assignment,:)),'ascend');

    plot(Dpoints,linspace(0,1,length(Dpoints)));

    xlim([0 0.1])
    ylim([0 1])
    grid on;
    xlabel('% dist euclid','FontSize',20)
    ylabel('% matches','FontSize',20)
    
    end
    
    if functionalMap
        subplot(2,2,3);
        imagesc(C);
        title('Current FM $\mathbf{C}$','interpreter','latex')
        colorbar
        
        subplot(2,2,4)
        imagesc(Cgt);
        title('Ground-truth FM $\mathbf{C_\mathrm{gt}}$','interpreter','latex')
        colorbar
    end


    drawnow
end

