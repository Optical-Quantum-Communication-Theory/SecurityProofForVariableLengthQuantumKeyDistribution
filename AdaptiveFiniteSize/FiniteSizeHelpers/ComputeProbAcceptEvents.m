function [ProbAcceptEvents] = ComputeProbAcceptEvents(ProbDist, tlist, m, numSamples)
% Samples from ProbDist and computes accept probability, by computing the number of times
% |Fobs-ProbDist| < tlist(i)
% 
% 
% Input
% * ProbDist : The probability distriibution we sample from.
% * tlist : The list of values of t for which we wish to compute accept probability.
% * m : The number of signals used for testing.
% * numSamples : The number of times we run the sampling procedure to estimate the accept probability.
%   NumSamples number of times to obtain Fobs
% Output
% * ProbAcceptEvents : A list of size the same as tlist, storing the
% computed accept probabilities. 
arguments
    ProbDist(:,1) double;
    tlist(:,1) double;
    m double;
    numSamples double;
end




%counts the number of times the samples Fobs was in the acceptance test. 
count = zeros(numel(tlist),1);

for i=1:numSamples
    if mod(i,5000)==0
        fprintf("Generating samples : %i out of %i done\n",i,numSamples);
    end
    


    Fobs = GenerateSample(ProbDist,m); %generate a frequency by sampling m times from ProbDist
    L1distance = sum(abs(Fobs(:) - ProbDist(:))); %compute L1 distance between Fobs and ProbDist.

    for j=1:numel(tlist)
        if L1distance <= tlist(j)
            count(j)=count(j)+1;
        end
    end
end

ProbAcceptEvents = count/numSamples;
end







