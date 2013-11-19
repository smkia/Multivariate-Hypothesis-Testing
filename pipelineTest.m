function [] = pipelineTest(cfgFileAddress,resultFileAddress)
load(cfgFileAddress);
% % Simulation Config
% cfgSim = [];
% cfgSim.trialNum = 10;
% cfgSim.channelNum = 274;
% cfgSim.timeNum = 50;
% cfgSim.freqNum = 45;
% cfgSim.affectedChannel = [208:209];
% cfgSim.sigmaVar = [4,7];
% cfgSim.effectSize = 11;
% cfgSim.timeVar = [21,35];
% cfgSim.freqVar = [15,29];
% % Preparing Targets
% % Hierarchical Method config
% %cfgHrc.coefNum = 5;
% cfgHrc.criticalAlpha = 0.05;
% cfgHrc.iterations = 10000;
% %cfgHrc.featureExt = 'DCT';
% cfgHrc.featureExt = 'DCT';
% cfgSim.ampVar = [7,12];
% cfgHrc.MCPMethod = {'BF','BF','BF'};
data_tf = [];
targets = zeros(1,cfgSim.trialNum);
targets(1,1:round(cfgSim.trialNum/2)) = 1;
[data_tf,mask,SNR] = simulatingData(cfgSim,data_tf);
for j = 1 : 30
    % Hierarchy test Bonferroni
    cfgHrc.coefNum = j;
    [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
    [sensitivityHierarchyBF(j),specificityHierarchyBF(j)] = testEvaluation(hierarchyMask,mask);
    save(resultFileAddress,'SNR','sensitivityHierarchyBF','specificityHierarchyBF');
    disp(strcat(num2str(j)));
end