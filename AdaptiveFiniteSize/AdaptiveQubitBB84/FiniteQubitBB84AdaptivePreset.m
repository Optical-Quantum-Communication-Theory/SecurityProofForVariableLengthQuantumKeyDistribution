function qkdInput = FiniteQubitBB84AdaptivePreset(type)
% FiniteQubitBB84AdaptivePreset A preset for BB84 with qubits
% We also supply a fixed set of t values, and the key rate is computed at
% all the t values
qkdInput = QKDSolverInput();

%% Parameters



qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("f",1.16); %effieciency of error-correction.
qkdInput.addFixedParameter("misalignmentAngle", 2*pi/180);
qkdInput.addFixedParameter("depolarization",0.02);



%finite size parameters.

qkdInput.addFixedParameter("alphabet", 2); % encoding alphabet size; for qubits, this is 2

% store all epsilons in one struct
epsilon.PE = (1/4)*1e-10; %failure probability for parameter estimation 
epsilon.EC = (1/2)*1e-10; %failure probability for error-correction
epsilon.PA = (1/4)*1e-10; %failure probability for privacy amplification
qkdInput.addFixedParameter("epsilon", epsilon);
qkdInput.addFixedParameter("N",1e6);
qkdInput.addFixedParameter("ptest", 0.05); %ptest*N signals used for testing
qkdInput.addFixedParameter("c",sqrt(log2(1/epsilon.PA))/log2(2+1))% alpha = 1 + c/sqrt(pgen*N); c = sqrt(log2(1/epsPA)) / log(dz+1);


qkdInput.addFixedParameter("numRuns",100000); % number of times we sample Fobs to estimate expected key rates/
qkdInput.addFixedParameter("use",type) % 0 means we plot fixed-length keyrates
% 1 means we plot variable-length keyrates.



%qkdInput.addFixedParameter("tExp", -4);

qkdInput.addFixedParameter("tExpList", num2cell(linspace(-2.5,-0.5,30)));
%qkdInput.addFixedParameter("tSiftExp"); we set tSiftExp = tExp;

%qkdInput.addFixedParameter("physDimAB", 2*3); % physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection


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
keyMod = QKDKeyRateModule(@FiniteQubitBB84AdaptiveKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 20  ;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = false;
% mathSolverOptions.infiniteDecoyIntensities = true; % Good for testing
% with an infinite number of decoy intensities (Upper and lower decoy
% bounds converge.)
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek","cvxPrecision","high"));
