function [keyRate, modParser, debugInfo] = FiniteQubitBB84AdaptiveKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% FiniteQubitBB84AdaptiveKeyRateFuncRenyi A finite size key rate function for a 
%  Adaptive Finite-size key BB84 protocol, which uses Renyi entropies. We take
% in a list of tExp values, and compute finitesize keyrate and accept prob
% for each value. We then compute expected keyrate. Accept probabilities
% and finite-size key rates for each value of tExp are stored in debugInfo.
%  
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * f: error correction effiency. A 1 means we are correcting at the
%   Shannon limit. A more practical value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments. These are organized in a 4x6 table. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * N: (finite size) the total number of signals Alice sent
% * eps: (finite size) a struct containing four epsilons used in finite key
%   rate analysis: PE (parameter estimation), bar (smoothing min entropy),
%   EC (error correction), PA (privacy amplification)
% * ptest: the fraction of signals sent that are used for testing.
%   Typically should be fairly low. This cuts into key rate, but also
%   entirely destroys key rate if it is too low.
% * tExpList : a finite size parameter list that decides the size of the acceptance
% test.
% * expectationsJoint : Joint expectations of Alice and bob's outcomes
% * ObservablesJoint : POVMs corresponding to Alice sending a given signal
%   and Bob obtaining a given outcome
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, PM46DescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


modParser = moduleParser(mfilename);

%% optionsParser = makeGlobalOptionsParser(mfilename);



modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJoint",@mustBeProbDist);
modParser.addRequiredParam("numRuns",@(x) mustBeGreaterThan(x, 0));
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

% modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsConditional"]);

%modParser.addRequiredParam("probSignalsA",@mustBeProbDist); %Add something for same number of announcementsA or something

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));
%modParser.addRequiredParam("probSignalsA",@mustBeProbDist); %Add something for same number of announcementsA or something



modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsConditional","announcementsA","announcementsB"])

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));
modParser.addRequiredParam("rhoA",@isDensityOperator)
%modParser.addOptionalParam("rhoA", nan, @isDensityOperator);
modParser.addRequiredParam("alphabet", @(x) mustBeInteger(x));
modParser.addRequiredParam("c",  @(x) mustBeGreaterThanOrEqual(x,0)); %alpha = 1 + c /sqrt( N);


%% finite key analysis parametersw
modParser.addRequiredParam("N", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("ptest", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("epsilon");
% modParser.addRequiredParam("t", @(x) mustBeInRange(x, 0, 1));
%modParser.addRequiredParam("tExpList", @(x) allCells(x, @(y) mustBeLessThanOrEqual(y, 0) ));
modParser.addRequiredParam("tExpList");
%modParser.addRequiredParam("tSiftExp", @(x) mustBeLessThanOrEqual(x,0)); %
%we set tExp = tSiftExp.



% modParser.addOptionalParam("blockDimsA", nan);
% modParser.addOptionalParam("blockDimsB", nan);
% modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
% modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
% modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

modParser.parse(params);

params = modParser.Results;

%% extra parameters
%modParser.addOptionalParam("acceptanceSetChoice", "none", @(x) mustBeMember(x, ["upper", "lower", "independent", "none"]));


%% parse parameters
modParser.parse(params);

params = modParser.Results;



% we need to do a key rate calculation for every value of tExp;
maxIndex = numel(params.tExpList);

keyRateList = zeros(maxIndex,1);
acceptProb = zeros(maxIndex,1);

%set up list of values of t.
for index=1:maxIndex
    tExp = params.tExpList{index};
    tlist(index) = 10^(tExp);
end

%compute accept probabiblity by sampling from params.expectationsJoint.
acceptProb = ComputeProbAcceptEvents(params.expectationsJoint(:),tlist,params.N*params.ptest,params.numRuns);

% This is the expected key rate keyrate over all accept events. This is
% unimportant, since we redo the calculation in a separate plotting
% function. 
keyRate = 0;


% for each value of t, we now do a key rate calculation. 
for index = 1: maxIndex
    fprintf(" \n Iteration : %d \n",index);
    debugKeyRateLeaves = debugInfo.addLeaves(strcat("FiniteKeyRateCalc",num2str(index)));

    params.tExp = params.tExpList{index}; % put the right value in tExp;
    [keyRateList(index), ~] = FiniteQubitBB84KeyRateFuncRenyi(params,options,mathSolverFunc,debugKeyRateLeaves);
   
end




% we compute the expected key rate 
for index = 1:(maxIndex - 1)
    keyRate = keyRate + acceptProb(index)*(keyRateList(index) - keyRateList(index+1));
end
keyRate = keyRate +acceptProb(maxIndex)*keyRateList(maxIndex);


debugInfo.storeInfo( "keyRateOnAcceptList", keyRateList);
debugInfo.storeInfo( "acceptProb",  acceptProb);

debugInfo.storeInfo("expectationsJoint", params.expectationsJoint); %actually computed in the channel model.
% needed for comparison when we compute fully adaptive 

end









function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end

function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("BasicKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end

function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end


function mapping = BB84StandardSquashingPostProccessingMap()
% squashes detector patterns from H,V,D,A (as bit string patterns 0000,
% 1000, ..., 1111) to qubit+vac values H,V,D,A,,vac.

    function mapping = quickMap(mapping,pattern,remaping)
        mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remaping;
    end

mapping = zeros(5,16);

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps.
mapping(5,:) = 1;

% vacume to vacume
mapping = quickMap(mapping,[0,0,0,0],[0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0],[1,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0],[0,1,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0],[0,0,1,0,0]); % D
mapping = quickMap(mapping,[0,0,0,1],[0,0,0,1,0]); % A

% double clicks
mapping = quickMap(mapping,[1,1,0,0],[0.5,0.5,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1],[0,0,0.5,0.5,0]); % X (DA)
end


function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="BasicBB84_WCPKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end


