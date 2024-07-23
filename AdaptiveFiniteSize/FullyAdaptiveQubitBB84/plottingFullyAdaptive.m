fsize = 50;
lw = 4;
ms = 13;


% we will plot expected keyrate for fixed-length and fully-variable length
% protocol




dataVar = load("./fullyAdaptiveOutput.mat");

numChannels = numel(dataVar.results); %each keyrate is for a different channel

dataFixedLength=load("./fixedOutputForComparison.mat");

KeyRateListUponAcceptFixed = dataFixedLength.results.debugInfo.keyRateModule.keyRateOnAcceptList;

FixedLengthKeyRateExpected = zeros(numel(KeyRateListUponAcceptFixed),1);
AcceptProbFixed = zeros(numel(KeyRateListUponAcceptFixed),numChannels);
variableKeyRates = zeros(numChannels,1);


if (cell2mat(dataVar.qkdInput.fixedParameters.tExpList) ~= ...
        cell2mat(dataFixedLength.qkdInput.fixedParameters.tExpList))
    warning('Fixed length and variable length data is for different set of tExpList!');
end

tExpList = dataFixedLength.qkdInput.fixedParameters.tExpList;
tList = 10.^(cell2mat(tExpList));



for index=1:numChannels % sum over all channels.
    variableKeyRates(index) = dataVar.results(index).keyRate;
    AcceptProbFixed(:,index) = dataVar.results(index).debugInfo.keyRateModule.acceptProb;
    FixedLengthKeyRateExpected = FixedLengthKeyRateExpected + (1/numChannels) *AcceptProbFixed(:,index).*KeyRateListUponAcceptFixed;
end


maxPoints = sum(KeyRateListUponAcceptFixed > 0,'all');
tList = tList(1:maxPoints);

FixedLengthKeyRateExpected=  FixedLengthKeyRateExpected(1:maxPoints);
KeyRateListUponAcceptFixed  = KeyRateListUponAcceptFixed (1:maxPoints);
VariableLengthKeyRate = mean(variableKeyRates,'all');








semilogx(tList,FixedLengthKeyRateExpected,'-o','Color',"#0072BD",'DisplayName','$\bar{R}_{\mathrm{fixed},i}$');
hold on;
semilogx(tList,KeyRateListUponAcceptFixed,'-*','Color',"#EDB120",'DisplayName','$R_{\mathrm{fixed},i}$');

yline(VariableLengthKeyRate,'Color',"#7E2F8E",'LineWidth',lw,'DisplayName','$\bar{R}_{\mathrm{variable}}$ ')



%Other things
xlabel('$t_i$','Interpreter','latex');
ylabel("Secure Key Bits / Signal Sents","Interpreter","latex");
legend('Location','best','Interpreter','latex');
fontsize(gca,fsize,"pixels")
xlim([tList(1),tList(numel(tList))]);
set(findall(gcf,'Type','line'),'LineWidth',lw);
set(findall(gcf,'Type','line'),'MarkerSize',ms);


savefig('fullyAdaptivePlot.fig');
%close();

hold off;


function [output] = NumberInsideAcceptanceSet(FbarFixedLength, t, FbarList)
% computes the number of Fbars in FbarList that are all within an entrywise
% distance t of FbarFixedLength
count = 0;

for index = 1:numel(FbarList)
    currentFbar = FbarList{index};
    dist = max(abs(currentFbar - FbarFixedLength),[],'all');
    if dist < t
        count = count+1;
    end
end

output = count;
end



