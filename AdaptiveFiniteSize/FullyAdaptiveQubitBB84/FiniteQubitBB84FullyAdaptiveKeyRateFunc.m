function [keyRate, modParser, debugInfo] = FiniteQubitBB84FullyAdaptiveKeyRateFunc(params,options,mathSolverFunc,mathSolverOptions,debugInfo)
% FiniteQubitBB84FullyAdaptiveKeyRateFuncRenyi A finite size key rate function for a 
%  Adaptive Finite-size key BB84 protocol, which uses Renyi entropies.
%  This function takes in a channel model and computes expected keyrate for
%  variable-length protocol and also the probability of accept for the
%  fixed-length protocols.
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
% * tExp : a finite size parameter that decides the size of the acceptance
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
    mathSolverOptions (1,1) struct
    debugInfo (1,1) DebugInfo
end

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


modParser = moduleParser(mfilename);

%% optionsParser = makeGlobalOptionsParser(mfilename);



modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));



modParser.addRequiredParam("expectationsJoint",@mustBeProbDist);
modParser.addRequiredParam("FbarFixed",@mustBeProbDist); %the fixed-length acceptance set centre point


modParser.addRequiredParam("numRuns",@(x) mustBeGreaterThan(x, 0));
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);
modParser.addAdditionalConstraint(@isEqualSize,["FbarFixed","expectationsJoint"]);

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
modParser.addRequiredParam("c",  @(x) mustBeGreaterThanOrEqual(x,0)); %alpha = 1 + c /sqrt(nsift - Ntsift);


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
maxIndex = numel(params.tExpList);
numRuns = params.numRuns;


keyRateList = zeros(numRuns,1); %will store the key rate for each run.
AcceptCounts = zeros(maxIndex,1); %counts the number of times protocol is accepted for fixed-length


for index=1:numRuns
    %generate sample
    Fobs = GenerateSample(params.expectationsJoint(:),params.ptest*params.N);
    %check what fixed-length set is falls in, and update counts
    %acccordingly. 
    Fobs = reshape(Fobs, size(params.expectationsJoint));
   
    for indexT=1:maxIndex
        tExp = params.tExpList{indexT};
        t = 10^(tExp);
        if (sum(abs(Fobs(:)-params.FbarFixed(:) )) <t  )
            AcceptCounts(indexT) = AcceptCounts(indexT)+1;
        end     
    end

    params.Fobs = Fobs;
    % now call the key rate function and compute keyrate
    fprintf(" \n Computing Key rate for sample : %i out of %i \n",index,numRuns);
    debugKeyRateLeaves = debugInfo.addLeaves(strcat("FiniteKeyRateCalc",num2str(index)));


    [keyRateList(index),~] = FiniteQubitBB84KeyRateFuncRenyiForFullyAdaptive(params,options,mathSolverFunc,mathSolverOptions,debugKeyRateLeaves);
    
end



%for this channel model

keyRate = sum(keyRateList,'all') / numRuns; %this is the key rate expected from this channel model for variable-length case
acceptProb = AcceptCounts / numRuns; %accept probability of different channel 








debugInfo.storeInfo( "keyRateOnAcceptList", keyRateList);
debugInfo.storeInfo( "acceptProb",  acceptProb);
%%debugInfo.storeInfo("acceptProbUpper",acceptProbUpper);
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


