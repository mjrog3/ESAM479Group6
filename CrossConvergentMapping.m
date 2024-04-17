numVars = 3;
steps = 10000;
data = zeros(steps, numVars);
data(:,1) = 1:steps;
data(:,2) = steps:-1:1;
LMesh = [10 20 30 40 50];

testOut = crossConvergentMap(data, 2, 3, LMesh, 5)
%%

function [correlations] = crossConvergentMap(data, timeDelay, embedDimension, LMesh, nearestNeighbors)
    X = data(:, 1);
    Y = data(:, 2);

    if(max(LMesh) > length(X))
        error("L larger than length of data")
    end

    correlations = zeros(length(LMesh), 2);

    for i = 1:length(LMesh)
        L  = LMesh(i);

        % number of data points that don't have enough history to be predicted
        skip = timeDelay*(embedDimension-1); 

        xManifoldPoints = zeros(L-skip, embedDimension);
        yManifoldPoints = zeros(L-skip, embedDimension);
        for j = 1:embedDimension
            startPoint = skip+1 - timeDelay*(j-1);
            xManifoldPoints(:,j) = X(startPoint:startPoint+(L-skip)-1);
            yManifoldPoints(:,j) = Y(startPoint:startPoint+(L-skip)-1);
        end
        
        % nearestNeighbors+1 because we're going to manually remove the
        % point itself after the knnsearch. 
        [xIdx, xD] = knnsearch(xManifoldPoints, xManifoldPoints,'K',nearestNeighbors+1); 
        [yIdx, yD] = knnsearch(yManifoldPoints, yManifoldPoints,'K',nearestNeighbors+1); 
        weightCap = 1000;
        xWeights = min(weightCap, 1./xD(:,2:end)); %threshold & drop self
        xWeights = xWeights./sum(xWeights, 2); %normalize weights
        yWeights = min(weightCap, 1./yD(:,2:end));
        yWeights = yWeights./sum(yWeights, 2);
        xPred = zeros(L-skip, 1);
        yPred = zeros(L-skip, 1);
        
        for j = 1:L-skip % for every row, take the weighted average of the images
            yPred(j) = xWeights(j,:) * Y(xIdx(j,2:end)); % 1xn * nx1 = 1x1
            xPred(j) = yWeights(j,:) * X(yIdx(j,2:end));
        end
        
        xCorrelation = corr(xPred, X(skip+1:L)); % currently not saving p-values but could
        yCorrelation = corr(yPred, Y(skip+1:L));
        correlations(i,:) = [xCorrelation yCorrelation];
    end

    %anything else we want to output?
end