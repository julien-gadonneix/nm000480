% compute oscillatory peaks based on stimulus epoch of all correct trials
%
% example: peaks_irasa('S1')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = peaks_irasa(subid)

pth = ; % INPUT PATH TO DATA
datdir = [pth subid '/'];
savdir = [datdir 'irasa/'];
mkdir(savdir);

zthr = 1.96;

% load data
load([datdir subid], 'data');

srate = data.fsample;

cor = cell2mat(data.trialinfo(:,6));
trialinfo = data.trialinfo;

% select task-responsive elecs
load([datdir 'hfb/chans_stim']);
load([datdir 'hfb/chans_resp']);

cfg = [];
cfg.trials = find(cor == 1);
cfg.channel = unique([chans_stim chans_resp]);

data = ft_selectdata(cfg, data);
data.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo chans*

% select stim onset --> resp onset, no buffers
resp = str2double(data.trialinfo(:,18));
resp = resp./srate;

for r = 1:length(data.trial)
    t1 = nearest(data.time{r}, 0);
    t2 = nearest(data.time{r}, resp(r));
    data.trial{r} = data.trial{r}(:,t1:t2);
    data.time{r} = data.time{r}(1,t1:t2);
    clear t1 t2
end

data.trialinfo = 1:size(data.trialinfo,1);
data.trialinfo = data.trialinfo';

% partition the data into ten overlapping sub-segments
[~,idx] = min(cellfun('size', data.trial, 2)); % find shortest trial
w = data.time{idx}(end)-data.time{idx}(1); % window length
cfg = [];
cfg.length = w*.9;
cfg.overlap = 1-((w-cfg.length)/(10-1));

tmp = ft_redefinetrial(cfg, data);

% perform irasa and regular spectral analysis
cfg = [];
cfg.foilim = [3 30];
cfg.taper = 'hanning';
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';

cfg.method = 'irasa';
frac = ft_freqanalysis(cfg, tmp);

cfg.method = 'mtmfft';
orig = ft_freqanalysis(cfg, tmp);
clear tmp*

% average across the sub-segments
tmpf = {};
tmpo = {};
for rpt = unique(frac.trialinfo)'
    cfg = [];
    cfg.trials = find(frac.trialinfo == rpt);
    cfg.avgoverrpt = 'yes';
    
    tmpf{end+1} = ft_selectdata(cfg, frac);
    tmpo{end+1} = ft_selectdata(cfg, orig);
end

frac = ft_appendfreq([], tmpf{:});
orig = ft_appendfreq([], tmpo{:});
clear tmp*

% average across trials
cfg = [];
cfg.avgoverrpt = 'yes';

frac = ft_selectdata(cfg, frac);
orig = ft_selectdata(cfg, orig);

% subtract the fractal component from the power spectrum
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2-x1';
osci = ft_math(cfg, frac, orig);

% detect peaks per electrode
osci.pk.pk = cell(length(osci.label),1);
osci.pk.fr = osci.pk.pk;
osci.pk.wd = osci.pk.fr;
osci.theta = nan(length(osci.label),1);
osci.alpha = osci.theta;
osci.beta = osci.theta;

osci.powspctrm = zcont_5(osci.powspctrm); % z-score

for e = 1:length(osci.label)
    [tmp_1, tmp_fr, tmp_wd, tmp_pk] = findpeaks(osci.powspctrm(e,:), osci.freq);
    tmp1 = tmp_1 > 0; tmp2 = tmp_pk >= zthr; tmp = tmp1 + tmp2 == 2;
    tmp_pk = tmp_pk(tmp); tmp_fr = tmp_fr(tmp); tmp_wd = tmp_wd(tmp);
    [~,tmp] = sort(tmp_pk, 'descend');
    tmp_pk = tmp_pk(tmp); tmp_fr = tmp_fr(tmp); tmp_wd = tmp_wd(tmp);
    if length(tmp_pk) > 3
        tmp_pk = tmp_pk(1:3); tmp_fr = tmp_fr(1:3); tmp_wd = tmp_wd(1:3);
    end
    [~,tmp] = sort(tmp_fr, 'ascend');

    osci.pk.pk{e} = tmp_pk(tmp);
    osci.pk.fr{e} = tmp_fr(tmp);
    osci.pk.wd{e} = tmp_wd(tmp);
    clear tmp
        
    tmp = osci.pk.fr{e} > 4 & osci.pk.fr{e} < 8.5;
    if sum(tmp) == 1
        osci.theta(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.theta(e) = tmp(tmp2(1));
    end
    clear tmp tmp2
    
    tmp = osci.pk.fr{e} > 8 & osci.pk.fr{e} < 13.5;
    if sum(tmp) == 1
        osci.alpha(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.alpha(e) = tmp(tmp2(1));
    end
    clear tmp tmp2
    
    tmp = osci.pk.fr{e} > 13;
    if sum(tmp) == 1
        osci.beta(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.beta(e) = tmp(tmp2(1));
    end
    clear tmp*
end

save([savdir 'osci'], 'osci', 'orig', 'frac');

clear all

end
