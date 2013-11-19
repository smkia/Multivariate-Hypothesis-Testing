%% Simulation 
% Simulation Config
cfgSim = [];
cfgSim.trialNum = 200;
cfgSim.channelNum = 274;
cfgSim.timeNum = 50;
cfgSim.freqNum = 45;
cfgSim.affectedChannel = [208:217];
cfgSim.sigmaVar = [4,7];
%cfgSim.ampVar = [7,12]; % 1
%cfgSim.ampVar = [70,120]; %2
%cfgSim.ampVar = [270,320]; %3
cfgSim.effectSize = 11;
cfgSim.timeVar = [21,35];
cfgSim.freqVar = [15,29];
% Preparing Targets
targets = zeros(1,cfgSim.trialNum);
targets(1,1:round(cfgSim.trialNum/2)) = 1;
% Hierarchical Method config
%cfgHrc.coefNum = 5;
cfgHrc.criticalAlpha = 0.05;
cfgHrc.iterations = 10000;
%cfgHrc.featureExt = 'DCT';

% Experiment 1: Different MCP correction methods for different SNRs
ampVar = [7,12;70,120;270,320];
cfgHrc.featureExt = 'DCT';
cfgHrc.coefNum = 5;
for j = 1 : size(ampVar,1)
    cfgSim.ampVar = ampVar(j,:);
    for i = 1 : 10
        % Data simulation
        [data_tf,mask,SNR(j,i)] = simulatingData(cfgSim,data_tf);
        % Cluster-based Test
        [clusterMask] = clusterBasedTest(data_tf,targets);
        [sensitivityCluster(j,i),specificityCluster(j,i)] = testEvaluation(clusterMask,mask);
        % Hierarchy test FDR-BH
        cfgHrc.MCPMethod = {'BH','BH','BH'};
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBH(j,i),specificityHierarchyBH(j,i)] = testEvaluation(hierarchyMask,mask);
        % Hierarchy test FDR-BR
        cfgHrc.MCPMethod = {'BR','BR','BR'};
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBR(j,i),specificityHierarchyBR(j,i)] = testEvaluation(hierarchyMask,mask);
        % Hierarchy test Bonferroni
        cfgHrc.MCPMethod = {'BF','BF','BF'};
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBF(j,i),specificityHierarchyBF(j,i)] = testEvaluation(hierarchyMask,mask);
        save('tempResult.mat','SNR','sensitivityCluster','specificityCluster','sensitivityHierarchyBH','specificityHierarchyBH', ...
            'sensitivityHierarchyBF','specificityHierarchyBF','sensitivityHierarchyBR','specificityHierarchyBR');
        disp(strcat(num2str(j),':',num2str(i)));
    end
end

% Experiment 2: Number of DCT coefficients
cfgHrc.featureExt = 'DCT';
cfgSim.ampVar = [7,12];
cfgHrc.MCPMethod = {'BF','BF','BF'};
for i = 1 : 10
    % Data simulation
    [data_tf,mask,SNR(i)] = simulatingData(cfgSim,data_tf);
    for j = 11 : 30
        % Hierarchy test Bonferroni
        cfgHrc.coefNum = j;
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBF(j,i),specificityHierarchyBF(j,i)] = testEvaluation(hierarchyMask,mask);
        save('tempResult.mat','SNR','sensitivityHierarchyBF','specificityHierarchyBF');
        disp(strcat(num2str(i),':',num2str(j)));
    end
end
save('Experiment2.mat','SNR','sensitivityHierarchyBF','specificityHierarchyBF');

% Experiment 3: KTST without DCT
ampVar = [7,12;14,25;35,60;70,120;270,320];
cfgHrc.MCPMethod = {'BF','BF','BF'};
for j = 1 : size(ampVar,1)
    cfgSim.ampVar = ampVar(j,:);
    for i = 1 : 10
        % Data simulation
        [data_tf,mask,SNR(j,i)] = simulatingData(cfgSim,data_tf);
%         % Just KTST
%         cfgHrc.featureExt = '';
%         [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
%         [sensitivityHierarchyKTST(j,i),specificityHierarchyKTST(j,i)] = testEvaluation(hierarchyMask,mask);
%         % KTST + DCT
%         cfgHrc.featureExt = 'DCT';
%         cfgHrc.coefNum = [cfgSim.freqNum,cfgSim.timeNum];
%         [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
%         [sensitivityHierarchyDCT(j,i),specificityHierarchyDCT(j,i)] = testEvaluation(hierarchyMask,mask);
        % KTST + DCT + Coef
        cfgHrc.featureExt = 'DCT';
        cfgHrc.coefNum = [5];
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyDCTCoef(j,i),specificityHierarchyDCTCoef(j,i)] = testEvaluation(hierarchyMask,mask);
        save('tempResult2.mat','SNR','sensitivityHierarchyKTST','specificityHierarchyKTST','sensitivityHierarchyDCT','specificityHierarchyDCT', ...
            'sensitivityHierarchyDCTCoef','specificityHierarchyDCTCoef');
        disp(strcat(num2str(i),':',num2str(j)));
    end
end

% Experiment 4: White noise checking
cfgSim.affectedChannel = [];
cfgHrc.featureExt = 'DCT';
coefNum = [5,15,30,45];
for j = 1 : length(coefNum)
    cfgHrc.coefNum = coefNum(j);
    for i = 1 : 10
        % Data simulation
        [data_tf,mask,SNR(j,i)] = simulatingData(cfgSim,data_tf);
        % Hierarchy test FDR-BH
        cfgHrc.MCPMethod = {'BH','BH','BH'};
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBH(j,i),specificityHierarchyBH(j,i)] = testEvaluation(hierarchyMask,mask);
        % Hierarchy test Bonferroni
        cfgHrc.MCPMethod = {'BF','BF','BF'};
        [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
        [sensitivityHierarchyBF(j,i),specificityHierarchyBF(j,i)] = testEvaluation(hierarchyMask,mask);
        save('tempResult3.mat','SNR','sensitivityHierarchyBH','specificityHierarchyBH', ...
            'sensitivityHierarchyBF','specificityHierarchyBF');
        disp(i);
    end
end
%% White noise checking
e = 0;
for k = 1 : 100
    data_tf.powspctrm = rand(600,1,45,60);
    [trialNum, channelNum, frequencyBinNum, timeBinNum] = size(data_tf.powspctrm);
    data_tf.powspctrm(isnan(data_tf.powspctrm)) = 0;
    coefNum = 5;
    criticalAlpha = 0.05;
    FDRMethod = {'BH','BH','BH'};
    ktstcfg = [];
    ktstcfg.iterations = 10000;
    channelACC = [];
    hChannels = [];
    pChannels = [];
    significantChannels = [];
    for z = 1 : channelNum
        features = zeros(coefNum*coefNum,trialNum);
        for i = 1 : trialNum
            D =[];
            D = dct2(squeeze(data_tf.powspctrm(i,z,:,:)));
            features(:,i) = reshape(D(1:coefNum,1:coefNum),coefNum*coefNum,1);
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        % KTST
        [pChannels(z)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        %disp(strcat(num2str(z),'/',num2str(channelNum)));
    end
    [hChannels] = FDR(pChannels,criticalAlpha,FDRMethod{1});
    significantChannels = find(hChannels);
    if ~isempty(significantChannels)
        e = e + 1;
    end
    disp(strcat(num2str(k),':',num2str(e)));
end
%% Hierarchi Permuted Labels
state = zeros(1,100);
for iter = 1 : 100
    [trialNum, channelNum, frequencyBinNum, timeBinNum] = size(data_tf.powspctrm);
    targets = zeros(trialNum,1);
    targets(randperm(600,300),1) = 1;
    data_tf.powspctrm(isnan(data_tf.powspctrm)) = 0;
    coefNum = 5;
    criticalAlpha = 0.05;
    FDRMethod = {'BH','BH','BH'};
    ktstcfg = [];
    ktstcfg.iterations = 10000;
    channelACC = [];
    hChannels = [];
    pChannels = [];
    for z = 1 : channelNum
        features = zeros(coefNum*coefNum,trialNum);
        for i = 1 : trialNum
            D =[];
            D = dct2(squeeze(data_tf.powspctrm(i,z,:,:)));
            features(:,i) = reshape(D(1:coefNum,1:coefNum),coefNum*coefNum,1);
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        % KTST
        [pChannels(z)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        %disp(strcat(num2str(z),'/',num2str(channelNum)));
    end
    [hChannels] = FDR(pChannels,criticalAlpha,FDRMethod{1});
    significantChannels{iter} = find(hChannels);
    if isempty(significantChannels{iter})
        state(1,iter) = 1;
        disp(strcat(num2str(iter),':',num2str(state(1,iter))));
        continue;
    end
    for i = 1 : length(significantChannels{iter})
        for j = 1 : frequencyBinNum
            features = zeros(coefNum,trialNum);
            for k = 1 : trialNum
                D =[];
                D = dct(squeeze(data_tf.powspctrm(k,significantChannels{iter}(i),j,:)));
                features(:,k) = reshape(D(1:coefNum),coefNum,1);
            end
            features(isnan(features)) = 0;
            features = mapstd(features);
            % KTST
            [pFreqs(i,j)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
            %disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j)));
        end
    end
    try
        [hFreqs{iter}] = FDR(pFreqs,criticalAlpha,FDRMethod{2});
        [temp, significantFreq{iter}] = find(hFreqs{iter});
        significantChannels2{iter} = significantChannels(temp);
        uSC{iter} = unique(significantChannels2{iter});
    catch
        hFreqs{iter} = [];
        significantFreq{iter} = [];
        significantChannels2{iter} = [];
        uSC{iter} = [];
    end
    if isempty(uSC{iter})
        state(1,iter) = 2;
        disp(strcat(num2str(iter),':',num2str(state(1,iter))));
        continue;
    end
    for i = 1 : length(uSC{iter})
        for j = 1 : timeBinNum
            features = zeros(coefNum,trialNum);
            for k = 1 : trialNum
                D =[];
                D = dct(squeeze(data_tf.powspctrm(k,uSC(i),:,j)));
                features(:,k) = reshape(D(1:coefNum),coefNum,1);
            end
            features(isnan(features)) = 0;
            features = mapstd(features);
            [pTime(i,j)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
            %disp(strcat(num2str(i),'/',num2str(length(uSC)),':',num2str(j)));
        end
    end
    [hTime{iter}] = FDR(pTime,criticalAlpha,FDRMethod{3});
    significantTime{iter} = find(hTime{iter});
    state(1,iter) = 3;
    disp(strcat(num2str(iter),':',num2str(state(1,iter))));
end



%% Evaluation
FP = sum(sum(sum(and(~m,stat.mask))))/sum(sum(sum(~m))); 
TP = sum(sum(sum(and(m,stat.mask))))/sum(sum(sum(m)));
FP = sum(sum(sum(and(~m,h))))/sum(sum(sum(~m))); 
TP = sum(sum(sum(and(m,h))))/sum(sum(sum(m)));
 