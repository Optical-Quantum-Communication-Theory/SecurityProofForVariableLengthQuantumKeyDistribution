function [keyRate, modParser, debugInfo] = FiniteQubitBB84KeyRateFuncRenyiForFullyAdaptive(params,options,mathSolverFunc,mathSolverOptions,debugInfo)
% FiniteQubitBB84KeyRateFuncRenyiForFullyAdaptive A finite size key rate function for a 
%  Finite-size key BB84 protocol, which uses Renyi entropies.  
% and compute key rate for each value. Basically t is set to zero, and the
% optimization set is "centred" around Fobs instead. 
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
% * 
% * eps: (finite size) a struct containing four epsilons used in finite key
%   rate analysis: PE (parameter estimation), bar (smoothing min entropy),
%   EC (error correction), PA (privacy amplification)
% * ptest: the fraction of signals sent that are used for testing.
%   Typically should be fairly low. This cuts into key rate, but also
%   entirely destroys key rate if it is too low.
% * tExp : a finite size parameter that decides the size of the acceptance
% test.
% * expectationsJoint : Joint expectations of Bob's outcome
%   Alice sent a given signal.
% * ObservablesJoint : POVMs corresponding to Alice sending a given signal
%   and Bob obtaining a given outcome,
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
modParser.addRequiredParam("Fobs",@mustBeProbDist);
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","Fobs"]);

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


%% finite key analysis parameters
modParser.addRequiredParam("N", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("ptest", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("epsilon");
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


debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();


%%% we start computations %%%



[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB, params.Fobs ,params.keyMap,params.f); 
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains",gains);



%% finite size calculations.
epsilon = params.epsilon;



m = params.N*params.ptest; % received testing signals
muBall = sqrt(2)*sqrt((log(1/epsilon.PE) + numel(params.Fobs)*log(m+1))/m);
params.muBall = muBall;



%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice


numObs = numel(params.observablesJoint);


mathSolverInput.vectorOneNormConstraints = VectorOneNormConstraint(params.observablesJoint(:),params.Fobs(:),muBall);


%% Translate for math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.rhoA = params.rhoA;


[relEnt,~] = mathSolverFunc(mathSolverInput,mathSolverOptions, debugMathSolver);


keyRate = finiteKeyRateRenyi(relEnt, deltaLeak, gains, params, options);
keyRate = max(keyRate,0);

if isfield(debugMathSolver.info,"relEntStep2Linearization")
    relEntStep2Linearization = debugMathSolver.info.relEntStep2Linearization; 
    
    keyRateStep2Linearization = finiteKeyRateRenyi( relEntStep2Linearization, deltaLeak, gains, params, options);
    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)
    
end




if options.verboseLevel>=1
    fprintf("Key rate: %e\n",keyRate);
end





debugInfo.storeInfo("keyRateOnAccept", keyRate);
debugInfo.storeInfo("renEnt",relEnt);
debugInfo.storeInfo("relEntStep2", debugMathSolver.info.relEntStep2Linearization);
debugInfo.storeInfo("Fobs",params.Fobs);

end








%%%%%%%%%%%  FINITE FADING CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [keyRate, debugInfo] = finiteKeyRateRenyi(relEnt, deltaLeak, gains, params, options, debugInfo)
%computes the finite size keyrate.

d = params.alphabet; %size of key register.
ptest = params.ptest;
pgen = 1 - ptest;

epsilon = params.epsilon;
N = params.N;

n=pgen*N;
alpha = 1 + params.c / sqrt(n);


RenyiTerm = ( sqrt(n) / N ) *params.c*(log2(2*d+1))^2;
privacyAmplification = (sqrt(n)/N)  *(alpha/params.c)*(2/alpha-log2(4*epsilon.PA));
ECLeakage = n/N*deltaLeak + ceil(log2(2/epsilon.EC))/N; %actually epsilon.EV this is!


keyRate = n/N*max(relEnt,0) - RenyiTerm - ECLeakage - privacyAmplification; 

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



