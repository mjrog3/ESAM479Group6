numVars = 3;
steps = 1000;
data = zeros(steps, numVars);

function [out] = crossConvergentMap(data, timeDelay, predictor, prediction)
    var1 = data(predictor, :);
    var2 = data(prediction, :);
end