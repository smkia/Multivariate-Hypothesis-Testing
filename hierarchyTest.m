function [h] = hierarchyTest(cfg,data_tf,targets)
[trialNum, channelNum, frequencyBinNum, timeBinNum] = size(data_tf.powspctrm);
data_tf.powspctrm(isnan(data_tf.powspctrm)) = 0;
coefNum = cfg.coefNum;
criticalAlpha = cfg.criticalAlpha;
FDRMethod = cfg.MCPMethod;
ktstcfg = [];
ktstcfg.iterations = cfg.iterations;
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
for i = 1 : length(significantChannels)
    data = squeeze(data_tf.powspctrm(:,significantChannels(i),:,:));
    data = padarray(data,[0,floor(coefNum/2),0],'replicate');
    for j = 1+floor(coefNum/2) : frequencyBinNum+floor(coefNum/2)
        features = zeros(coefNum,trialNum);
        for k = 1 : trialNum
            D =[];
            D = dct2(squeeze(data(k,j-2:j+2,:)));
            features(:,k) = reshape(D(1,1:coefNum),coefNum,1);
%                       features(:,k) = squeeze(data_tf.powspctrm(k,significantChannels(i),j,:));
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        % KTST
        [pFreqs(i,j-floor(coefNum/2))] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-floor(coefNum/2)),':',num2str(pFreqs(i,j-floor(coefNum/2)))));
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
    data = padarray(data,[0,0,floor(coefNum/2)],'replicate');
    for j = 1+floor(coefNum/2) : timeBinNum+floor(coefNum/2) 
        features = zeros(coefNum,trialNum);
        for k = 1 : trialNum
            D =[];
            D = dct2(squeeze(data(k,:,j-2:j+2)));
            features(:,k) = reshape(D(1:coefNum,1),coefNum,1);
        end
        features(isnan(features)) = 0;
        features = mapstd(features);
        [pTime(i,j-floor(coefNum/2))] = KTST(features(:,targets == 0)',features(:,targets == 1)',ktstcfg);
        disp(strcat(num2str(i),'/',num2str(length(significantChannels)),':',num2str(j-floor(coefNum/2)),':',num2str(pTime(i,j-floor(coefNum/2)))));
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
h = zeros(channelNum,frequencyBinNum,timeBinNum);
for i = 1 : length(significantChannels)
    sc(i).freqs = significantFreq(significantFreqChan == significantChannels(i));
    sc(i).time = significantTime(significantTimeChan == significantChannels(i));
end
for i = 1 : length(significantChannels)
    h(significantChannels(i),sc(i).freqs,sc(i).time) = 1;
end
% localizer = [];
% l = 1;
% for i = 1 : channelNum
%     for j = 1 : frequencyBinNum
%         for k = 1 : timeBinNum
%             if h(i,j,k) == 1
%                 localizer(l,1) = i;
%                 localizer(l,2) = j;
%                 localizer(l,3) = k;
%                 l = l + 1;
%             end
%         end
%     end
% end