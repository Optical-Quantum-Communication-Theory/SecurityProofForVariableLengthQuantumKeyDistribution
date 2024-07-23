function [Fobs] = GenerateSample(ProbDist,m)
% generates a sample by sampling from ProbDist m times.
% 
% * Input 
% * ProbDist : The probability distribution to sample from
% * m : Number of times to sample.
%
% Output
% * Fobs : The frequency of outcome in the sample.

numOutcomes = numel(ProbDist);
Fobs = zeros(numOutcomes,1);

Samples = randsample(numOutcomes, m, true, ProbDist);

for i = 1:numOutcomes
    Fobs(i) = sum(Samples==i)/m;
end

end

