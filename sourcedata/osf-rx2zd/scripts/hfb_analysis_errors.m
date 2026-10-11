% compute hfb (70-150 hz) on error trials. run stats on hfb activity:
% effect size. generate chance distributions and z-score
%
% example: hfb_analysis_errors('S1')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = hfb_analysis_errors(subid)

pth = ; % INPUT PATH TO DATA
datdir = [pth subid '/'];
savdir = [datdir 'hfb/'];
mkdir(savdir);

% load data
load([datdir subid], 'data');

srate = data.fsample;
times = str2double(data.trialinfo(:,16:17));

% 1000-ms data time windows (in sec, 500-ms on screen + 250-ms pre/post fixation)
s1 = [-0.25 0.75]; % 0 is onset
s2 = [times(:,1)/srate-0.25 times(:,1)/srate+0.75];
s3 = [times(:,2)/srate-0.25 times(:,2)/srate+0.75];

stim_max = max(s3(:,2));
stim = data;

% epoch stim-locked with buffers
for r = 1:length(data.trial)
    t1 = nearest(data.time{r}, -1);
    t2 = nearest(data.time{r}, stim_max+0.5);
    stim.trial{r} = data.trial{r}(:,t1:t2);
    stim.time{r} = -1:1/srate:stim_max+0.5;
    clear t1 t2
end
clear data

% set up filter
cfg = [];
cfg.pad = 'nextpow2';
cfg.method = 'mtmconvol';
cfg.taper = 'dpss';
cfg.foi = 75:10:145;
cfg.tapsmofrq = 10/2;
cfg.t_ftimwin = repmat(0.3, [1 length(cfg.foi)]); % length of sliding time window
cfg.output = 'pow'; % power
cfg.keeptrials = 'yes';

% compute hfb
cfg.toi = stim.time{1}(1):0.005:stim.time{1}(end);
data_stim = ft_freqanalysis(cfg, stim);
clear stim

t1 = nearest(data_stim.time, -0.25);
t2 = nearest(data_stim.time, -0.05);
base_mtx = data_stim.powspctrm(:,:,:,t1:t2);

% select error trials and balance trial counts for z-scoring by randomly 
% selecting correct trials to make up the difference
cor = cell2mat(data_stim.trialinfo(:,6));
trialinfo = data_stim.trialinfo;

ndiff = sum(cor == 1) - sum(cor == 0);
diff = randsample(find(cor == 1), ndiff);

cfg = [];
cfg.trials = [find(cor == 0); diff];

data_stim = ft_selectdata(cfg, data_stim);
data_stim.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo *diff

base_mtx = base_mtx(cfg.trials,:,:,:);

% z-score on baseline
data_stim = zbaseline_200(data_stim, base_mtx);
clear base_mtx

% select error trials
cor = cell2mat(data_stim.trialinfo(:,6));
trialinfo = data_stim.trialinfo;

cfg = [];
cfg.trials = find(cor == 0);

data_stim = ft_selectdata(cfg, data_stim);
data_stim.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo

% average over frequencies and remove buffers
cfg = [];
cfg.avgoverfreq = 'yes';

cfg.latency = [-0.25 max(s3(:,2))];
data_stim = ft_selectdata(cfg, data_stim);

data_stim.dimord = 'rpt_chan_time';
data_stim.trial = reshape(data_stim.powspctrm, [size(data_stim.powspctrm,1) ...
    size(data_stim.powspctrm,2) size(data_stim.powspctrm,4)]);
data_stim = rmfield(data_stim, {'powspctrm','freq'});

% select stim-locked time windows
dat_s1 = data_stim;
dat_s1.time = s1(1):0.005:s1(2);
dat_s1.trial = nan(size(data_stim.trial,1), size(data_stim.trial,2), length(dat_s1.time));
dat_s2 = dat_s1; dat_s3 = dat_s1;

for r = 1:size(data_stim.trial,1)
    tmp_on = nearest(data_stim.time, s1(1)); tmp_off = nearest(data_stim.time, s1(2));
    if tmp_off-tmp_on ~= length(dat_s1.time)-1
        tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
        tmp_off = tmp_off+tmp_dif;
    end
    dat_s1.trial(r,:,:) = data_stim.trial(r,:,tmp_on:tmp_off); clear tmp*
    tmp_on = nearest(data_stim.time, s2(r,1)); tmp_off = nearest(data_stim.time, s2(r,2));
    if tmp_off-tmp_on ~= length(dat_s1.time)-1
        tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
        tmp_off = tmp_off+tmp_dif;
    end
    dat_s2.trial(r,:,:) = data_stim.trial(r,:,tmp_on:tmp_off); clear tmp*
    tmp_on = nearest(data_stim.time, s3(r,1)); tmp_off = nearest(data_stim.time, s3(r,2));
    if tmp_off-tmp_on ~= length(dat_s1.time)-1
        tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
        tmp_off = tmp_off+tmp_dif;
    end
    dat_s3.trial(r,:,:) = data_stim.trial(r,:,tmp_on:tmp_off); clear tmp*
end

% concatenate all error targets and save by condition
if ~strcmp(subid,'CP34') % rename
    % no cl error trials
    [target, distractor, load1, load2] = split_cf_cl(dat_s1, dat_s2, dat_s3);
    data_stim = target;
    idx1 = cell2mat(load1.trialinfo(:,14)) == '1';
    idx2 = cell2mat(load2.trialinfo(:,14)) == '2';
    data_stim.trialinfo = cat(1, data_stim.trialinfo, load1.trialinfo(idx1,:), load2.trialinfo(idx2,:));
    data_stim.trial = cat(1, data_stim.trial, load1.trial(idx1,:,:), load2.trial(idx2,:,:));
    save([savdir 'data_stim_cl_inc'], 'load*');
else
    [target, distractor] = split_cf_cl(dat_s1, dat_s2, dat_s3);
    data_stim = target;
end
save([savdir 'data_stim_locked_inc'], 'data_stim');
save([savdir 'data_stim_cf_inc'], 'target', 'distractor');
clear data_stim

% stats
load([savdir 'chans_stim']);
load([savdir 'chans_resp']);

% select stimulus epoch and task-responsive elecs
cfg = [];
cfg.channel = unique([chans_stim, chans_resp]); clear chans*
cfg.latency = [0 0.75];
cfg.avgovertime = 'yes';
cfg.nanmean = 'yes';

% context first
target = ft_selectdata(cfg, target);
distractor = ft_selectdata(cfg, distractor);

hfb = [];
hfb.dimord = 'rpt_chan';
hfb.trialinfo = target.trialinfo;
hfb.label = target.label;
hfb.target = target.trial;
hfb.distractor = distractor.trial;

% mean difference
hfb.diff = squeeze(nanmean(hfb.target-hfb.distractor,1));
hfb.diff = hfb.diff';

save([savdir 'cf_hfb_inc'], 'hfb'); 
clear hfb target distractor

% context last
if ~strcmp(subid,'CP34') % rename
    load2 = ft_selectdata(cfg, load2);
    load1 = ft_selectdata(cfg, load1);
    
    hfb = [];
    hfb.dimord = 'rpt_chan';
    hfb.trialinfo = load2.trialinfo;
    hfb.label = load2.label;
    hfb.load2 = load2.trial;
    hfb.load1 = load1.trial;
    
    % mean difference
    hfb.diff = squeeze(nanmean(hfb.load2-hfb.load1,1));
    hfb.diff = hfb.diff';
    
    save([savdir 'cl_hfb_inc'], 'hfb');
end

clear all

end
