% compute hfb (70-150 hz) and select task-responsive electrodes
%
% example: hfb_analysis('S1')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = hfb_analysis(subid)

addpath([pwd '/toolboxes/fieldtrip/external/spm2']);

pth = ; % INPUT PATH TO DATA
datdir = [pth subid '/'];
savdir = [datdir 'hfb/'];
mkdir(savdir);

% set task-responsive elecs threshold
zthr = 1.96;
dur = 0.05;

% load data
load([datdir subid], 'data');

srate = data.fsample;

cor = cell2mat(data.trialinfo(:,6));
trialinfo = data.trialinfo;

cfg = [];
cfg.trials = find(cor == 1);

data = ft_selectdata(cfg, data);
data.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo 

rt = round(cell2mat(data.trialinfo(:,5)).*srate); % in srate
rt_min = min(rt);

times = str2double(data.trialinfo(:,16:17));

% 1000-ms data time windows (in sec, 500-ms on screen + 250-ms pre/post fixation)
s1 = [-0.25 0.75]; % 0 is onset
s2 = [times(:,1)/srate-0.25 times(:,1)/srate+0.75];
s3 = [times(:,2)/srate-0.25 times(:,2)/srate+0.75];

stim_max = max(s3(:,2));

stim = data;
resp = data;

% epoch stim- and response-locked + response with buffers (min rt for response-locked)
for r = 1:length(data.trial)
    t1 = nearest(data.time{r}, -1);
    t2 = nearest(data.time{r}, stim_max+0.5);
    stim.trial{r} = data.trial{r}(:,t1:t2);
    stim.time{r} = -1:1/srate:stim_max+0.5;
    clear t1 t2
    
    t2 = length(data.time{r}); % rt + 1-s buffer
    t1 = t2 - 2*srate - rt_min;
        
    resp.time{r} = nan(1,length(t1:t2));
    resp.time{r}(end) = 1;    
    if r == 1
        for t = length(t1:t2)-1:-1:1
            resp.time{r}(t) = resp.time{r}(t+1)-1/srate; % response at time 0
        end
    else
        resp.time{r} = resp.time{r-1};
    end
    resp.trial{r} = data.trial{r}(:,t1:t2);
    clear t1 t2
end

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

cfg.toi = resp.time{1}(1):0.005:resp.time{1}(end);
data_resp = ft_freqanalysis(cfg, resp);
clear resp

t1 = nearest(data_stim.time, -0.25);
t2 = nearest(data_stim.time, -0.05);
base_mtx = data_stim.powspctrm(:,:,:,t1:t2);

% z-score on baseline
data_stim = zbaseline_200(data_stim, base_mtx);
data_resp = zbaseline_200(data_resp, base_mtx);
clear base_mtx

% average over frequencies and remove buffers
cfg = [];
cfg.avgoverfreq = 'yes';

cfg.latency = [-0.25 max(s3(:,2))];
data_stim = ft_selectdata(cfg, data_stim);

cfg.latency = [data_resp.time(1)+1 0.05];
data_resp = ft_selectdata(cfg, data_resp);

data_stim.dimord = 'rpt_chan_time';
data_stim.trial = squeeze(data_stim.powspctrm);
data_stim = rmfield(data_stim, {'powspctrm','freq'});

data_resp.dimord = 'rpt_chan_time';
data_resp.trial = squeeze(data_resp.powspctrm);
data_resp = rmfield(data_resp, {'powspctrm','freq'});

save([savdir 'data_resp_locked'], 'data_resp');

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

% concatenate all correct targets and save by condition
[target, distractor, load1, load2] = split_cf_cl(dat_s1, dat_s2, dat_s3);
data_stim = target;
idx1 = cell2mat(load1.trialinfo(:,14)) == '1';
idx2 = cell2mat(load2.trialinfo(:,14)) == '2';
data_stim.trialinfo = cat(1, data_stim.trialinfo, load1.trialinfo(idx1,:), load2.trialinfo(idx2,:));
data_stim.trial = cat(1, data_stim.trial, load1.trial(idx1,:,:), load2.trial(idx2,:,:));

save([savdir 'data_stim_locked'], 'data_stim');

save([savdir 'data_stim_cf'], 'target', 'distractor');
clear target distractor

save([savdir 'data_stim_cl'], 'load*');
clear load*

% select task-responsive elecs
cfg = [];
cfg.avgoverrpt = 'yes';
cfg.nanmean = 'yes';

cfg.latency = [0 0.75];
data_stim = rmfield(data_stim, 'cumtapcnt');
data_stim = ft_selectdata(cfg, data_stim);

cfg.latency = [data_resp.time(1) 0];
data_resp = rmfield(data_resp, 'cumtapcnt');
data_resp = ft_selectdata(cfg, data_resp);

mask = data_stim.trial >= zthr;

c = 0;
for e = 1:size(mask,1)
    [tmp,tmpn] = spm_bwlabel(double(mask(e,:)),6);

    % min duration
    for b = 1:tmpn
        tmp_indx = find(tmp==b);
        tmp_dur = max(data_stim.time(tmp_indx)) - min(data_stim.time(tmp_indx));
        if tmp_dur >= dur
            c = c+1;
            chans_stim(c) = data_stim.label(e);
        end
        clear tmp_*
    end
    clear tmp*
end
if ~exist('chans_stim','var')
    chans_stim = [];
end
chans_stim = unique(chans_stim);

save([savdir 'chans_stim'], 'chans_stim');
clear mask

mask = data_resp.trial >= zthr;

c = 0;
for e = 1:size(mask,1)
    [tmp,tmpn] = spm_bwlabel(double(mask(e,:)),6);

    % min duration
    for b = 1:tmpn
        tmp_indx = find(tmp==b);
        tmp_dur = max(data_resp.time(tmp_indx)) - min(data_resp.time(tmp_indx));
        if tmp_dur >= dur
            c = c+1;
            chans_resp(c) = data_resp.label(e);
        end
        clear tmp_*
    end
    clear tmp*
end
if ~exist('chans_resp','var')
    chans_resp = [];
end
chans_resp = unique(chans_resp);

save([savdir 'chans_resp'], 'chans_resp');

clear all

end
