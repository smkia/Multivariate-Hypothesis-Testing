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
cfgSim.ampVar = [270,320]; %3
cfgSim.effectSize = 11;
cfgSim.timeVar = [21,35];
cfgSim.freqVar = [15,29];
% Preparing Targets
targets = zeros(1,cfgSim.trialNum);
targets(1,1:round(cfgSim.trialNum/2)) = 1;
% Hierarchical Method config
cfgHrc.coefNum = 5;
cfgHrc.criticalAlpha = 0.05;
cfgHrc.iterations = 10000;
% Loop
for i = 1 : 10
    % Data simulation
    [data_tf,mask,SNR(i)] = simulatingData(cfgSim,data_tf);
    % Cluster-based Test
    [clusterMask] = clusterBasedTest(data_tf,targets);
    [sensitivityCluster(i),specificityCluster(i)] = testEvaluation(clusterMask,mask);
    % Hierarchy test FDR-BH
    cfgHrc.MCPMethod = {'BH','BH','BH'};
    [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
    [sensitivityHierarchyBH(i),specificityHierarchyBH(i)] = testEvaluation(hierarchyMask,mask);
    % Hierarchy test FDR-BR
    cfgHrc.MCPMethod = {'BR','BR','BR'};
    [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
    [sensitivityHierarchyBR(i),specificityHierarchyBR(i)] = testEvaluation(hierarchyMask,mask);
    % Hierarchy test Bonferroni
    cfgHrc.MCPMethod = {'BF','BF','BF'};
    [hierarchyMask] = hierarchyTest(cfgHrc,data_tf,targets);
    [sensitivityHierarchyBF(i),specificityHierarchyBF(i)] = testEvaluation(hierarchyMask,mask);
    save('tempResult.mat','SNR','sensitivityCluster','specificityCluster','sensitivityHierarchyBH','specificityHierarchyBH', ...
        'sensitivityHierarchyBF','specificityHierarchyBF','sensitivityHierarchyBR','specificityHierarchyBR');
    disp(i);
end

%% Hierarchi-KTST

[trialNum, channelNum, frequencyBinNum, timeBinNum] = size(data_tf.powspctrm);
data_tf.powspctrm(isnan(data_tf.powspctrm)) = 0;
criticalAlpha = 0.05;
FDRMethod = {'BH','BH','BH'};
ktstcfg = [];
coefNum = 5;
ktstcfg.iterations = 10000;
hChannels = [];
pChannels = [];
for z = 1 : channelNum
    features = zeros(frequencyBinNum*timeBinNum,trialNum);
    for i = 1 : trialNum
        features(:,i) = reshape(data_tf.powspctrm(i,z,:,:),frequencyBinNum*timeBinNum,1);
    end
    features(isnan(features)) = 0;
    features = mapstd(features);
    % KTST
    [pChannels(z)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
    disp(strcat(num2str(z),'/',num2str(channelNum),':',num2str(pChannels(z))));
end
[hChannels] = FDR(pChannels,criticalAlpha,FDRMethod{1});
significantChannels = find(hChannels);

for i = 1 : length(significantChannels)
    data = squeeze(data_tf.powspctrm(:,significantChannels(i),:,:));
    data = padarray(data,[0,floor(coefNum/2),0],'replicate');
    for j = 1+floor(coefNum/2) : frequencyBinNum+floor(coefNum/2)
        features = zeros(coefNum*timeBinNum,trialNum);
        for k = 1 : trialNum
            features(:,k) = reshape(squeeze(data(k,j-2:j+2,:)),coefNum*timeBinNum,1);
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        % KTST
        [pFreqs(i,j-floor(coefNum/2))] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-floor(coefNum/2)),':',num2str(pFreqs(i,j-floor(coefNum/2)))));
    end
end
[hFreqs] = FDR(pFreqs,criticalAlpha,FDRMethod{2});
[temp, significantFreq] = find(hFreqs);
significantFreqChan = significantChannels(temp);

for i = 1 : length(significantChannels)
    data = squeeze(data_tf.powspctrm(:,significantChannels(i),:,:));
    data = padarray(data,[0,0,floor(coefNum/2)],'replicate');
    for j = 1+floor(coefNum/2) : timeBinNum+floor(coefNum/2) 
        features = zeros(coefNum*frequencyBinNum,trialNum);
        for k = 1 : trialNum
            features(:,k) = reshape(squeeze(data(k,:,j-2:j+2)),coefNum*frequencyBinNum,1);
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        [pTime(i,j-floor(coefNum/2))] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-floor(coefNum/2)),':',num2str(pTime(i,j-floor(coefNum/2)))));
    end
end
[hTime] = FDR(pTime,criticalAlpha,FDRMethod{3});
[temp, significantTime] = find(hTime);
significantTimeChan = significantChannels(temp);
h = zeros(channelNum,frequencyBinNum,timeBinNum);
for i = 1 : length(significantChannels)
    sc(i).freqs = significantFreq(significantFreqChan == significantChannels(i));
    sc(i).time = significantTime(significantTimeChan == significantChannels(i));
end
for i = 1 : length(significantChannels)
    h(significantChannels(i),sc(i).freqs,sc(i).time) = 1;
end
localizer = [];
l = 1;
for i = 1 : channelNum
    for j = 1 : frequencyBinNum
        for k = 1 : timeBinNum
            if h(i,j,k) == 1
                localizer(l,1) = i;
                localizer(l,2) = j;
                localizer(l,3) = k;
                l = l + 1;
            end
        end
    end
end

%% Plotting the results 2

temp = [];
temp{1} = localizer((localizer(:,2)>=3 & localizer(:,2)<=6),:);
temp{2} = localizer((localizer(:,2)>=7 & localizer(:,2)<=14),:);
temp{3} = localizer((localizer(:,2)>=15 & localizer(:,2)<=30),:);
temp{4} = localizer((localizer(:,2)>=31 & localizer(:,2)<=45),:);
TF1 = data_tf;
TF1.powspctrm(101:200,:,:,:) = [];
TF2 = data_tf;
TF2.powspctrm(1:100,:,:,:) = [];
for i = 1 : 4
    if ~isempty(temp{i})
        cfg = [];
        cfg.foilim = minmax(temp{i}(:,2)');
        data_TF_planar_cmb_left1 = ft_freqdescriptives(cfg, TF1);
        data_TF_planar_cmb_right1  = ft_freqdescriptives(cfg, TF2);
        raweffect = data_TF_planar_cmb_left1.powspctrm - data_TF_planar_cmb_right1.powspctrm;
        plotFormat.label = data_TF_planar_cmb_left1.label;
        plotFormat.freq = i;
        plotFormat.dimord = 'chan_freq_time';
        plotFormat.time = data_TF_planar_cmb_left1.time;
        plotFormat.cfg = data_TF_planar_cmb_left1.cfg;
        plotFormat.grad = data_TF_planar_cmb_left1.grad;
        plotFormat.powspctrm = raweffect;
        plotFormat.raweffect = raweffect;
        plotcfg = [];
        plotcfg.interpolation = 'v4';
        plotcfg.layout = 'CTF275.lay';
        plotcfg.comment = 'no';
        plotcfg.parameter = 'raweffect';
        timePoints = [];
        timePoints = unique(temp{i}(:,3));
        figNum = ceil(length(timePoints)/20);
        
        for k = 1 : figNum
            figPointer = figure;
            if k == figNum && figNum ~= 1 && mod(length(timePoints),20) ~= 0
                n = mod(length(timePoints),20);
            else
                n = 20;
            end
            for j = 1 : min(n,length(timePoints))
                plotcfg.highlight = 'on';
                chan = temp{i}(temp{i}(:,3)==timePoints(j+(k-1)*20),1);
                plotcfg.highlightchannel =  plotFormat.label(chan,1);
                plotcfg.highlightsymbol = '*';
                plotcfg.highlightcolor =  [0 0 0]; %(black))
                plotcfg.highlightsize = 6;
                plotcfg.highlightfontsize = 8;
                plotcfg.marker = 'off';
                subplot(4,5,j);
                ft_topoplotTFR(plotcfg,plotFormat);
                title(strcat('Time = ',num2str(plotFormat.time(timePoints(j+(k-1)*20))),' - Freq = ',num2str(cfg.foilim(1)),'-',num2str(cfg.foilim(2))));
            end
            %saveas(figPointer,strcat('H:\NILAB\Master Thesis\Figures\Hierarchy Cross-Validation final\sub3_ktst_',num2str(cfg.foilim(1)),'-',num2str(cfg.foilim(2)),'_',num2str(k)),'fig');
            %close(figPointer);
        end
    else
        continue;
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
 