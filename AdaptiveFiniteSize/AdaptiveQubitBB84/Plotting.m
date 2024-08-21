fSize = 18;
lineWidth = 2;
markerSize = 6;

%% We start plotting the fixed-length keyrates

dataFix=load("./fixedOutput.mat");

tExpListFix = dataFix.qkdInput.fixedParameters.tExpList;
tList = 10.^(cell2mat(tExpListFix));


keyRateDataFix = dataFix.results.debugInfo.keyRateModule;
keyRateListAcceptFix = max(keyRateDataFix.keyRateOnAcceptList,0); % list of key rates when we accept;
%keyRateAdaptive = max(data.results.keyRate,0);

probAcceptFix = keyRateDataFix.acceptProb;


fixedLengthKeyRateExpected = probAcceptFix.*keyRateListAcceptFix;


if ~IsInNonIncreasingOrder(keyRateListAcceptFix,1e-3)
    warning('Key Rate List is not in decreasing order');
end


%% Now we start with variable-length keyrates
dataVar=load("./variableOutput.mat");

tExpListVar = dataVar.qkdInput.fixedParameters.tExpList;

if prod(cell2mat(tExpListVar) == cell2mat(tExpListFix))==0  %if a single value is unequal, this gives zero.
    warning('Value of t for variable and fixed-length are not identical!');
end



keyRateDataVar = dataVar.results.debugInfo.keyRateModule;
keyRateListAcceptVar = max(keyRateDataVar.keyRateOnAcceptList,0); % list of key rates when we accept;

probAcceptVar = keyRateDataVar.acceptProb;




if ~IsInNonIncreasingOrder(keyRateListAcceptVar,1e-2)
    warning('Key Rate List Var is not in decreasing order');
end



maxIndex = numel(keyRateListAcceptVar);
%i want to compute adaptive keyrate Lower bound
keyRateAdaptive = 0;
for index = 1:(maxIndex - 1)
    temp = (keyRateListAcceptVar(index) - keyRateListAcceptVar(index+1));
    keyRateAdaptive = keyRateAdaptive + probAcceptVar(index)*(temp);
end

keyRateAdaptive= keyRateAdaptive +probAcceptVar(maxIndex)*keyRateListAcceptVar(maxIndex);



% %i want to compute adaptive keyrate upper bound
% keyRateAdaptiveUpperBound = 0;
% for index = 1:(maxIndex - 1)
%     temp = (keyRateListAcceptVar(index) - keyRateListAcceptVar(index+1));
%     if temp>0
%         keyRateAdaptiveUpperBound = keyRateAdaptiveUpperBound + probAcceptUpperVar(index)*(temp);
%     else 
%          keyRateAdaptiveUpperBound = keyRateAdaptiveUpperBound + probAcceptLowerVar(index)*(temp);
%     end
% end
% 
% keyRateAdaptiveUpperBound= keyRateAdaptiveUpperBound +probAcceptUpperVar(index)*keyRateListAcceptVar(maxIndex);
% 


%axes('XScale', 'log', 'YScale', 'log')

%some of the early values in 

maxPoints = sum(keyRateListAcceptFix > 0,'all');

tList = tList(1:maxPoints);


fixedLengthKeyRateExpected = fixedLengthKeyRateExpected(1:maxPoints);
%fixedLengthKeyRateUpper = fixedLengthKeyRateUpper(1:maxPoints);
keyRateListAcceptFix = keyRateListAcceptFix(1:maxPoints);

semilogx(tList,fixedLengthKeyRateExpected,'-.<','Color',"#0072BD",...
    'DisplayName','$\bar{R}_{\mathrm{fixed},i}$');
hold on;
% loglog(tList,fixedLengthKeyRateUpper,'-.+','Color',"#D95319",...
%     'DisplayName','$\bar{R}^{\mathrm{u}}_{\mathrm{fixed},i}$');

semilogx(tList,keyRateListAcceptFix,'-.*','Color',"#EDB120",...
    'DisplayName','$R_{\mathrm{fixed},i}$');

yline(keyRateAdaptive,"--",'Color',"#7E2F8E","LineWidth",lineWidth,...
    'DisplayName','$\bar{R}_{\mathrm{variable}}$ ')

% yline(keyRateAdaptiveUpperBound,"-.",'Color',"#77AC30",'LineWidth',4,...
%     'DisplayName','$\bar{R}^{\mathrm{u}}_{\mathrm{variable}}$')


%Other things
xlabel('$t_i$','Interpreter','latex');
ylabel('Secure key bits / signal sent','Interpreter','latex');
%ylabel('Objective Function');
legend('Location','best','Interpreter','latex');
set(gca,"FontSize",fSize);
% fontsize(gca,fSize,"pixels")

set(findall(gcf,'Type','line'),'LineWidth',lineWidth);
set(findall(gcf,'Type','line'),'MarkerSize',markerSize);

xlim([tList(1), tList(numel(tList))])
savefig('./adaptivePlot.fig');
%saveas(gcf,'adaptivePlot.png');
%close();

hold off;


function output = IsInNonIncreasingOrder(list,tol)
output = true;

for i=1:(numel(list)-1)
    if list(i+1) > list(i) +tol
        output = false;
        return;
    end
end
return;


end