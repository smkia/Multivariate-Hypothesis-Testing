function [h] = hierarchyTest(cfg,data_tf,targets)
[trialNum, channelNum, frequencyBinNum, timeBinNum] = size(data_tf.powspctrm);
data_tf.powspctrm(isnan(data_tf.powspctrm)) = 0;
if strcmp(cfg.featureExt,'DCT')
    if numel(cfg.coefNum) == 1
        coefNum(1,1) = cfg.coefNum;
        coefNum(1,2) = cfg.coefNum;
    else
        coefNum = cfg.coefNum;
    end
end
criticalAlpha = cfg.criticalAlpha;
FDRMethod = cfg.MCPMethod;
ktstcfg = [];
ktstcfg.iterations = cfg.iterations;
slidingWinSize = 2;
h = zeros(channelNum,frequencyBinNum,timeBinNum);
for z = 1 : channelNum
    if strcmp(cfg.featureExt,'DCT')
        features = zeros(coefNum(1,1)*coefNum(1,2),trialNum);
        for i = 1 : trialNum
            D =[];
            D = dct2(squeeze(data_tf.powspctrm(i,z,:,:)));
            features(:,i) = reshape(D(1:coefNum(1,1),1:coefNum(1,2)),coefNum(1,1)*coefNum(1,2),1);
        end
    else
        features = zeros(frequencyBinNum*timeBinNum,trialNum);
        for i = 1 : trialNum
            features(:,i) = reshape(data_tf.powspctrm(i,z,:,:),frequencyBinNum*timeBinNum,1);
        end
    end
    features(isnan(features)) = 0;
    features = mapstd(features);
    % KTST
    [pChannels(z)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
    disp(strcat(num2str(z),'/',num2str(channelNum),':',num2str(pChannels(z))));
end
if ~strcmp(FDRMethod{1},'BF')
    [hChannels] = FDR(pChannels,criticalAlpha,FDRMethod{1});
else
    correctedAlpha = criticalAlpha/channelNum;
    correctedAlpha = max(1/cfg.iterations,correctedAlpha);
    [hChannels] = double(pChannels <= correctedAlpha);
end
significantChannels = find(hChannels);
if isempty(significantChannels)
    return;
end
for i = 1 : length(significantChannels)
    data = squeeze(data_tf.powspctrm(:,significantChannels(i),:,:));
    data = padarray(data,[0,slidingWinSize,0],'replicate');
    for j = 1+slidingWinSize : frequencyBinNum+slidingWinSize
        if strcmp(cfg.featureExt,'DCT')
            features = zeros(coefNum(1,2),trialNum);
            for k = 1 : trialNum
                D =[];
                D = dct2(squeeze(data(k,j-slidingWinSize:j+slidingWinSize,:)));
                features(:,k) = reshape(D(1,1:coefNum(1,2)),coefNum(1,2),1);
            end
        else
            features = zeros(((2*slidingWinSize)+1)*timeBinNum,trialNum);
            for k = 1 : trialNum
                features(:,k) = reshape(squeeze(data(k,j-slidingWinSize:j+slidingWinSize,:)),((2*slidingWinSize)+1)*timeBinNum,1);
            end
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        % KTST
        [pFreqs(i,j-slidingWinSize)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-slidingWinSize),':',num2str(pFreqs(i,j-slidingWinSize))));
    end
end
if ~strcmp(FDRMethod{2},'BF')
    [hFreqs] = FDR(pFreqs,criticalAlpha,FDRMethod{2});
else
    correctedAlpha = criticalAlpha/(length(significantChannels)*frequencyBinNum);
    correctedAlpha = max(1/cfg.iterations,correctedAlpha);
    [hFreqs] = double(pFreqs <= correctedAlpha);
end
[temp, significantFreq] = find(hFreqs);
significantFreqChan = significantChannels(temp);
for i = 1 : length(significantChannels)
    data = squeeze(data_tf.powspctrm(:,significantChannels(i),:,:));
    data = padarray(data,[0,0,slidingWinSize],'replicate');
    for j = 1+slidingWinSize : timeBinNum+slidingWinSize
        if strcmp(cfg.featureExt,'DCT')
            features = zeros(coefNum(1,1),trialNum);
            for k = 1 : trialNum
                D =[];
                D = dct2(squeeze(data(k,:,j-slidingWinSize:j+slidingWinSize)));
                features(:,k) = reshape(D(1:coefNum(1,1),1),coefNum(1,1),1);
            end
        else
            features = zeros(((2*slidingWinSize)+1)*frequencyBinNum,trialNum);
            for k = 1 : trialNum
                features(:,k) = reshape(squeeze(data(k,:,j-slidingWinSize:j+slidingWinSize)),((2*slidingWinSize)+1)*frequencyBinNum,1);
            end
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        [pTime(i,j-slidingWinSize)] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-slidingWinSize),':',num2str(pTime(i,j-slidingWinSize))));
    end
end
if ~strcmp(FDRMethod{3},'BF')
    [hTime] = FDR(pTime,criticalAlpha,FDRMethod{3});
else
    correctedAlpha = criticalAlpha/(length(significantChannels)*timeBinNum);
    correctedAlpha = max(1/cfg.iterations,correctedAlpha);
    [hTime] = double(pTime <= correctedAlpha);
end
[temp, significantTime] = find(hTime);
significantTimeChan = significantChannels(temp);
for i = 1 : length(significantChannels)
    sc(i).freqs = significantFreq(significantFreqChan == significantChannels(i));
    sc(i).time = significantTime(significantTimeChan == significantChannels(i));
end
for i = 1 : length(significantChannels)
    h(significantChannels(i),sc(i).freqs,sc(i).time) = 1;
end