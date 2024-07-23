%pick the preset file to use 
qkdInput = FiniteQubitBB84FullyAdaptivePreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.


save("./fullyAdaptiveOutput.mat","results","qkdInput");
 
%plotting
plottingFullyAdaptive()

%% plot the result
%load("adaptiveOutput.mat");
%plotResults(results, qkdInput, 'dB-log')
% plotResults(results, qkdInput, 'linear')