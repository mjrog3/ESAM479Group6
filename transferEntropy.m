numRows = 50;
steps = 10000000;
epsilon = 0.1;
tic
data = tentMapData(numRows, steps, epsilon);
toc
tic
% frequency of each tuple in the data
possibleTuples = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
tupleFrequencies = zeros(length(possibleTuples), 1);
for i = 1:length(possibleTuples)
    tupleMatch = (data(1:end-1, 1:end-1) == possibleTuples(i, 1)).*...
                 (data(2:end  , 1:end-1) == possibleTuples(i, 2)).*...
                 (data(2:end  , 2:end  ) == possibleTuples(i, 3));
    tupleFrequencies(i) = sum(tupleMatch, "all");
end

toc