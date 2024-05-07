numRows = 100;
steps = 10000;
N = numRows*(steps-1);
%epsilon = 0.05;
numTrials = 10;
epsilonMesh = linspace(0.001, 0.05, 20);
tEntResults = zeros(numTrials, length(epsilonMesh));
noiseAmp = 0;
jModifier = +1; % for creating the alter data
tic
for epsilonIndex = 1:length(epsilonMesh)
    epsilon = epsilonMesh(epsilonIndex);
    for trialIndex = 1:numTrials
        data = tentMapData(numRows, steps, epsilon, noiseAmp);
        
        % frequency of each tuple in the data
        possibleTuples = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
        tupleFrequencies = zeros(length(possibleTuples), 1);
        alterPast = data([1+mod(1+jModifier, numRows):end-1 1:mod(1+jModifier, numRows)], 1:end-1);
        for i = 1:length(possibleTuples)
            tupleMatch = (alterPast == possibleTuples(i, 1)).*...
                         (data(2:end  , 1:end-1) == possibleTuples(i, 2)).*...
                         (data(2:end  , 2:end  ) == possibleTuples(i, 3));
            tupleFrequencies(i) = sum(tupleMatch, "all");
        end
        
        tupleProbablities = tupleFrequencies/N;
        tEnt = 0;
        for i = 1:length(possibleTuples)
            matchingAll = tupleProbablities(i);
            indices = find((possibleTuples(:,1) == possibleTuples(i,1)).*(possibleTuples(:,2) == possibleTuples(i,2)));
            matching12 = sum(tupleProbablities(indices));
            indices = find(possibleTuples(:,2) == possibleTuples(i,2));
            matching2 = sum(tupleProbablities(indices));
            indices = find((possibleTuples(:,2) == possibleTuples(i,2)).*(possibleTuples(:,3) == possibleTuples(i,3)));
            matching23 = sum(tupleProbablities(indices));
            tEntContribution = matchingAll*log(matchingAll/matching12*matching2/matching23);
            tEnt = tEnt + tEntContribution;
        end
        
        tEntResults(trialIndex,epsilonIndex) = tEnt;
    end
end
toc

%% Plotting


means = mean(tEntResults, 1);
stdErrors = std(tEntResults, 1);
% fit alpha for the quadratic
alpha = fminbnd(@(a) sum( (means-a^2*epsilonMesh.^2/log(2)).^2 ),0,0.8);

errorbar(epsilonMesh, means, stdErrors, 'LineWidth',1.5);
hold on
plot(epsilonMesh, alpha^2*epsilonMesh.^2/log(2), '--', 'LineWidth',1.5)
legend("Observed", "Theoretical")
xlabel('Coupling strength $\epsilon$', 'Interpreter','latex','FontSize',14)
ylabel('Transfer entropy [bits]', 'Interpreter', 'latex', 'FontSize',14)

%% Testing
a = 1:10;
b = 1:10;
testMat = a'*b;

testMat = [testMat(end,:) ; testMat];

modifier = -2;
slicendice = testMat([1+mod(1+modifier, 10):end-1 1:mod(1+modifier,10)], :);
