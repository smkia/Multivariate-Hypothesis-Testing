function [data_tf,mask,meanSNR] = simulatingData(cfg,data_tf)
data_tf.powspctrm = [];
data_tf.time = -.45:0.05:(cfg.timeNum-10)*0.05;
data_tf.powspctrm = zeros(cfg.trialNum,cfg.channelNum,cfg.freqNum,cfg.timeNum);
m = zeros(cfg.trialNum,cfg.channelNum,cfg.freqNum,cfg.timeNum);
for i = 1 : cfg.trialNum
    data_tf.powspctrm(i,:,:,:) = data_tf.powspctrm(i,:,:,:) + rand(1,cfg.channelNum,cfg.freqNum,cfg.timeNum);
end
SNR = zeros(round(cfg.trialNum/2),length(cfg.affectedChannel));
for i = 1 : round(cfg.trialNum/2)
    for j = 1: length(cfg.affectedChannel)
        g = GaussianFilter(cfg.effectSize,randi(cfg.sigmaVar,1))*randi(cfg.ampVar,1);
        tCenter = randi(cfg.timeVar,1);
        fCenter = randi(cfg.freqVar,1);
        v = floor(cfg.effectSize/2);
        noise  = data_tf.powspctrm(i,cfg.affectedChannel(j),:,:);
        signal = zeros(size(data_tf.powspctrm(1,1,:,:)));
        signal(1,1,fCenter - v : fCenter + v, tCenter - v : tCenter + v) = ...
            signal(1,1,fCenter - v : fCenter + v, tCenter - v : tCenter + v) + shiftdim(g,-2);
        data_tf.powspctrm(i,cfg.affectedChannel(j),:,:) = noise + signal;
        SNR(i,j) = snr(squeeze(noise),squeeze(signal));
        a = zeros(size(g));
        a(g ~= 0) = 1;
        m(i,cfg.affectedChannel(j),fCenter - v : fCenter + v, tCenter - v : tCenter + v) ...
            = m(i,cfg.affectedChannel(j),fCenter - v : fCenter + v, tCenter -v : tCenter + v) + shiftdim(a,-2);
    end
end
mask = zeros(cfg.channelNum,cfg.freqNum,cfg.timeNum);
for i = 1 : cfg.channelNum
    for j = 1 : 100
        mask(i,:,:) = shiftdim(or(squeeze(mask(i,:,:)),squeeze(m(j,i,:,:))),-1);
    end
end
meanSNR = mean2(SNR);