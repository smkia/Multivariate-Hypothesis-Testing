function [] = plotresults(config,data_tf,targets,mask)
[m,n,o]=size(mask);
localizer = [];
l = 1;
for i = 1 : m
    for j = 1 : n
        for k = 1 : o
            if mask(i,j,k) == 1
                localizer(l,1) = i;
                localizer(l,2) = j;
                localizer(l,3) = k;
                l = l + 1;
            end
        end
    end
end
temp = [];
temp{1} = localizer(round(data_tf.freq((localizer(:,2))))>=3 & round(data_tf.freq(localizer(:,2)))<=6,:);
temp{2} = localizer(round(data_tf.freq((localizer(:,2))))>=7 & round(data_tf.freq(localizer(:,2)))<=14,:);
temp{3} = localizer(round((data_tf.freq(localizer(:,2))))>=15 & round(data_tf.freq(localizer(:,2)))<=30,:);
temp{4} = localizer(round((data_tf.freq(localizer(:,2))))>=31 & round(data_tf.freq(localizer(:,2)))<=45,:);
data_tf1 = data_tf;
data_tf1.powspctrm(targets==1,:,:,:) = [];
data_tf2 = data_tf;
data_tf2.powspctrm(targets==0,:,:,:) = [];
clear data_tf;
for i = 1 : 4
    if ~isempty(temp{i})
        cfg = [];
        cfg.foilim = minmax(round(data_tf1.freq(temp{i}(:,2))));
        data_tf1_avg = ft_freqdescriptives(cfg, data_tf1);
        data_tf2_avg = ft_freqdescriptives(cfg, data_tf2);
        raweffect = data_tf1_avg.powspctrm - data_tf2_avg.powspctrm;
        plotFormat.label = data_tf1.label;
        plotFormat.freq = i;
        plotFormat.dimord = 'chan_freq_time';
        plotFormat.time = data_tf1.time;
        plotFormat.cfg = data_tf1.cfg;
        plotFormat.grad = data_tf1.grad;
        plotFormat.powspctrm = raweffect;
        plotFormat.raweffect = raweffect;
        plotcfg = [];
        plotcfg.interpolation = 'v4';
        plotcfg.layout = 'neuromag306cmb.lay';
        plotcfg.comment = 'no';
        plotcfg.parameter = 'raweffect';
        timePoints = [];
        timePoints = data_tf1.time(unique(temp{i}(:,3)+min(config.timeIndex)-1));
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
                chan = temp{i}(data_tf1.time(temp{i}(:,3)+min(config.timeIndex)-1)==timePoints(j+(k-1)*20),1);
                plotcfg.highlightchannel =  plotFormat.label(chan,1);
                plotcfg.highlightsymbol = '*';
                plotcfg.highlightcolor =  [0 0 0]; %(black))
                plotcfg.highlightsize = 8;
                plotcfg.highlightfontsize = 8;
                plotcfg.marker = 'off';
                subplot(4,5,j);
                ft_topoplotTFR(plotcfg,plotFormat);
                title(strcat('Time = ',num2str(timePoints(j+(k-1)*20)),' - Freq = ',num2str(cfg.foilim(1)),'-',num2str(cfg.foilim(2))));
            end
            if ~isempty(config.path)
                saveas(figPointer,strcat(config.path,config.subName,num2str(cfg.foilim(1)),'-',num2str(cfg.foilim(2)),'_',num2str(k)),'fig');
                close(figPointer);
            end
        end
    else
        continue;
    end
end

