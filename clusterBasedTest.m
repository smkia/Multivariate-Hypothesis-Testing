function [mask] = clusterBasedTest(data_tf,targets)
%trialNum = size(data_tf.powspctrm,1);
cfg = [];
cfg.channel          = {'MEG'};
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, data_tf);
cfg.design(1,targets == 1) = 1;
cfg.design(1,targets == 0) = 2;
cfg.ivar = 1;
%cfg.avgoverchan = 'yes';
%cfg.avgovertime = 'yes';
%cfg.avgoverfreq = 'yes';
[stat] = ft_freqstatistics(cfg,data_tf);
mask = double(stat.mask);