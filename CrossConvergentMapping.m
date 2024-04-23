%% Some things to keep in mind:

% When plotting, make sure that you manually rename the legend names
% to the appropriate "__ Predictions". For the asymmetric system, always
% "X" then "Y" but for Lorenz will need to be able to be "X" and "Z" etc.

% If using a non-zero noiseAmp factor then make sure that you generate
% fresh data every time you run the CCM block. Otherwise, every time you
% run it, it will add even more noise to the data.

% Make sure that your LMesh is fairly reasonable given the number of steps
% you took with your data. I think a good rule of thumb is to let LMesh be
% every 10% of your steps (so if steps is 1000, LMesh can be 100, 200, 300,
% ..., 1000). Too fine of a mesh and it will take too long to run but too
% coarse of a mesh and the figures won't be as insightful (or pretty)

% the CCM function now also outputs the time-delayed embeddings of the two
% variables so that we can plot those and include them in the report (see:
% the Manifold Visualization section of the file.
%% Generate Lorenz Data
numVars = 3;
steps = 4000;
tRange = [0 100];
out = generate_lorenz(10, 8/3, 28, [1 1 1], linspace(tRange(1), tRange(2), steps));
time = out(:,1);
x = out(:,2);
y = out(:,3);
z = out(:,4);
%% reselect Lorenz variables without having to re-run sim
var1 = y;
var2 = z;

data = [var1 var2];

%% Generate Coupled System Data
steps = 10000;
rx = 3.7;
ry = 3.7; 
betayx = 0.01;
betaxy = 0.0;
IC = [0.2 0.4];
data = coupled_system(IC, steps, rx, ry, betaxy, betayx);
%% CCM

noiseAmp = 0;
data = data + noiseAmp * randn(size(data));

% data(:,1) = 1:steps;
% data(:,2) = steps:-1:1;

LMesh = 100:100:steps;

[correlations, xM, yM] = crossConvergentMap(data, 30, 3, LMesh); 

%% plotting

figure()
% if the var1 predictions are converging with increasing time, that means
% that var2 can predict var1 which means that var1 causes var2
plot(LMesh, correlations(:,1), LMesh, correlations(:,2),'LineWidth',1.5)
legend("Y Predictions", "Z Predictions")
%ylim([-1 1])
ylabel("Predictive power corr$(\hat{Z}, Z)$",'Interpreter','latex')
% figure()
% plot(time, x, time, y, time, z);
% hold on
% %libraryTimes = time(LMesh);
% %xline(libraryTimes);
% legend("X", "Y", "Z")
% %xlim([20 30])
figure()
plot(data(:,1))
hold on
plot(data(:,2))
% testing if the step size is small enough to capture the dynamics
% figure()
% plot(time, x, '-o');
% hold on
% plot(time(1:10:end), x(1:10:end), '-o', 'MarkerFaceColor', [0.3 0.3 0.7])
% xlim([24 25])

%% Manifold Visualization
figure()
tiledlayout(1,2, "TileSpacing","tight")
nexttile()
plot3(xM(:,1),xM(:,2),xM(:,3));
xlabel("$X_t$",'Interpreter','latex');
ylabel("$X_{t-\tau}$", 'Interpreter','latex');
zlabel("$X_{t-2\tau}$", 'Interpreter','latex');
nexttile()
plot3(yM(:,1),yM(:,2),yM(:,3));
xlabel("$Y_t$",'Interpreter','latex');
ylabel("$Y_{t-\tau}$", 'Interpreter','latex');
zlabel("$Y_{t-2\tau}$", 'Interpreter','latex');

%%

function [correlations, var1Manifold, var2Manifold] = crossConvergentMap(data, timeDelay, embedDimension, LMesh)
    X = data(:, 1);
    Y = data(:, 2);

    % number of data points that don't have enough history to be predicted
    skip = timeDelay*(embedDimension-1); 
    
    if(min(LMesh) - skip < 1)
        error("Smallest L too short to make predictions with")
    end

    if(max(LMesh) > length(X))
        error("L larger than length of data")
    end

    correlations = zeros(length(LMesh), 2);
    %RMSE = zeros(length(LMesh), 2);

    for i = 1:length(LMesh)
        L  = LMesh(i);

        xManifoldPoints = zeros(L-skip, embedDimension);
        yManifoldPoints = zeros(L-skip, embedDimension);
        for j = 1:embedDimension
            startPoint = skip+1 - timeDelay*(j-1);
            xManifoldPoints(:,j) = X(startPoint:startPoint+(L-skip)-1);
            yManifoldPoints(:,j) = Y(startPoint:startPoint+(L-skip)-1);
        end
        
        % embedDimension+2 because we're going to manually remove the
        % point itself after the knnsearch and the paper says to do E + 1. 
        [xIdx, xD] = knnsearch(xManifoldPoints, xManifoldPoints,'K',embedDimension+2); 
        [yIdx, yD] = knnsearch(yManifoldPoints, yManifoldPoints,'K',embedDimension+2); 
        xWeights = exp(-xD(:,2:end)./max(xD(:,2),0.001));% no div by 0 in this house
        %xWeights = min(weightCap, 1./xD(:,2:end)); %threshold & drop self
        xWeights = xWeights./sum(xWeights, 2); %normalize weights
        %yWeights = min(weightCap, 1./yD(:,2:end));
        yWeights = exp(-yD(:,2:end)./max(yD(:,2),0.001));
        yWeights = yWeights./sum(yWeights, 2);
        xPred = zeros(L-skip, 1);
        yPred = zeros(L-skip, 1);
        
        for j = 1:L-skip % for every row, take the weighted average of the images
            yPred(j) = xWeights(j,:) * Y(xIdx(j,2:end)+skip); % 1xn * nx1 = 1x1
            xPred(j) = yWeights(j,:) * X(yIdx(j,2:end)+skip);
        end
        
        %xRMSE = rmse(xPred, X(skip+1:L));
        %yRMSE = rmse(yPred, Y(skip+1:L));
        %RMSE(i, :) = [xRMSE yRMSE];
        xCorrelation = corr(xPred, X(skip+1:L)); % currently not saving p-values but could
        yCorrelation = corr(yPred, Y(skip+1:L));
        correlations(i,:) = [xCorrelation yCorrelation];

        if i == length(LMesh)
            var1Manifold = xManifoldPoints;
            var2Manifold = yManifoldPoints;
        end
    end
    %anything else we want to output?

end