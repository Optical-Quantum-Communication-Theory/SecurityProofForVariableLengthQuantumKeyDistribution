
for type = 0:1 %decides fixed vs variable-length
    %pick the preset file to use
    qkdInput = FiniteQubitBB84AdaptivePreset(type);

    %run the QKDSolver with this input
    results = MainIteration(qkdInput);
    
    if qkdInput.fixedParameters.use == 0 %(type = 0))
        filename = "./fixedOutput";
        save("../FullyAdaptiveQubitBB84/fixedOutputForComparison","results","qkdInput")
    elseif qkdInput.fixedParameters.use == 1 %(type = 1)
        filename = "./variableOutput";
    end

    save(filename,"results","qkdInput");

end

Plotting(); %call the plotting function.


%% plot the result
%load("adaptiveOutput.mat");
%plotResults(results, qkdInput, 'dB-log')
% plotResults(results, qkdInput, 'linear')