function qkdInput = FiniteQubitBB84FullyAdaptivePreset()
% FiniteQubitBB84FullyAdaptivePreset A preset for BB84 with qubits, where
% the key rate is computed in the fully adaptive case. 
qkdInput = QKDSolverInput();

%% Parameters



qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("f",1.16); %effiecienct of error-correction.

%Our channel has 20 possible values. 
qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(2*pi/180,10*pi/180,5)));
qkdInput.addScanParameter("depolarization",num2cell(linspace(0.02,0.05,4)));



%finite size parameters.

qkdInput.addFixedParameter("alphabet", 2); % encoding alphabet size; for qubits, this is 2

% store all epsilons in one struct
epsilon.PE = (1/4)*1e-12; %failure probability for parameter estimation 
epsilon.EC = (1/2)*1e-12; %failure probability for error-correction
epsilon.PA = (1/4)*1e-12; %failure probability for privacy amplification
qkdInput.addFixedParameter("epsilon", epsilon);
qkdInput.addFixedParameter("N",1e6);
qkdInput.addFixedParameter("ptest", 0.05); %ptest*N signals used for testing
qkdInput.addFixedParameter("c",sqrt(log2(1/epsilon.PA))/log2(2+1))% alpha = 1 + c/sqrt(pgen*N); c = sqrt(log2(1/epsPA)) / log(dz+1);



qkdInput.addFixedParameter("numRuns",50); % number of times we sample Fobs per channel behaviour.

%we take in the fixed-length protocol we are supposed to compare against.
%Make sure this is a "fair" comparison (i.e identical protocol)
datafixed = load('./fixedOutputForComparison.mat');
tExpList = datafixed.qkdInput.fixedParameters.tExpList;
FbarFixed = datafixed.results.debugInfo.keyRateModule.expectationsJoint;


%these are the fixed-length acceptance sets for which we must compute
%probability of accept. 
qkdInput.addFixedParameter("tExpList",tExpList);
qkdInput.addFixedParameter("FbarFixed",FbarFixed);



% description is the same as the lossy qubit description since we squash
% Bob's detector data down to a lossy qubit equivalent
descriptionModule = QKDDescriptionModule(@BasicBB84Alice2DDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model is very different from normal qubit
channelModule = QKDChannelModule(@BasicBB84Alice2DChannelFunc);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis
keyRateOptions = struct();
keyRateOptions.decoyTolerance = 1e-14;
keyRateOptions.decoySolver = "mosek";
keyRateOptions.decoyForceSep = true;
keyMod = QKDKeyRateModule(@FiniteQubitBB84FullyAdaptiveKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 25  ;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = false;
% mathSolverOptions.infiniteDecoyIntensities = true; % Good for testing
% with an infinite number of decoy intensities (Upper and lower decoy
% bounds converge.)
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek","cvxPrecision","high"));
