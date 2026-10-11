% run stats on hfb activity: effect size, correlation. generate chance 
% distributions and z-score
%
% example: hfb_stats('S1')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = hfb_stats(subid)

pth = ; % INPUT PATH TO DATA
datdir = [pth subid '/hfb/'];
savdir = datdir;

nperm = 1000;

% load data
load([datdir 'data_stim_cf']);
load([datdir 'data_stim_cl']);

load([datdir 'chans_stim']);
load([datdir 'chans_resp']);

% select stimulus epoch and task-responsive elecs
cfg = [];
cfg.channel = unique([chans_stim, chans_resp]); clear chans*
cfg.latency = [0 0.75];
cfg.avgovertime = 'yes';
cfg.nanmean = 'yes';

target = ft_selectdata(cfg, target);
distractor = ft_selectdata(cfg, distractor);
load2 = ft_selectdata(cfg, load2);
load1 = ft_selectdata(cfg, load1);
    
% stats
rt_cf = cell2mat(target.trialinfo(:,5));
rt_cl = cell2mat(load2.trialinfo(:,5));
jitter_cf = str2double(target.trialinfo(:,17))-str2double(target.trialinfo(:,16));
jitter_cl = str2double(load2.trialinfo(:,17))-str2double(load2.trialinfo(:,16));

% context first
hfb = [];
hfb.dimord = 'rpt_chan';
hfb.trialinfo = target.trialinfo;
hfb.label = target.label;
hfb.target = target.trial;
hfb.distractor = distractor.trial;

% mean difference
hfb.diff = squeeze(nanmean(hfb.target-hfb.distractor,1));
hfb.diff = hfb.diff';

% true correlation
hfb.r = nan(length(hfb.label),1);
for e1 = 1:size(hfb.r,1)
    hfb.r(e1) = partialcorr(hfb.target(:,e1)-hfb.distractor(:,e1), rt_cf, jitter_cf, 'Type', 'Spearman');
end

% chance distribution
hfb.r_perm = nan(length(hfb.label), nperm);
for z = 1:nperm
    x = randperm(length(rt_cf));
    for e1 = 1:size(hfb.r,1)
        hfb.r_perm(e1,z) = partialcorr(hfb.target(x,e1)-hfb.distractor(x,e1), rt_cf, jitter_cf, 'Type', 'Spearman');
    end
    clear x
end

% z-score
hfb.r_z = (hfb.r-nanmean(hfb.r_perm,2))./nanstd(hfb.r_perm,[],2);

save([savdir 'cf_hfb'], 'hfb'); 
clear hfb target distractor *cf

% context last
hfb = [];
hfb.dimord = 'rpt_chan';
hfb.trialinfo = load2.trialinfo;
hfb.label = load2.label;
hfb.load2 = load2.trial;
hfb.load1 = load1.trial;

% mean difference
hfb.diff = squeeze(nanmean(hfb.load2-hfb.load1,1));
hfb.diff = hfb.diff';

% true correlation
hfb.r = nan(length(hfb.label),1);
for e1 = 1:size(hfb.r,1)
    hfb.r(e1) = partialcorr(hfb.load2(:,e1)-hfb.load1(:,e1), rt_cl, jitter_cl, 'Type', 'Spearman');
end

% chance distribution
hfb.r_perm = nan(length(hfb.label), nperm);
for z = 1:nperm
    x = randperm(length(rt_cl));
    for e1 = 1:size(hfb.r,1)
        hfb.r_perm(e1,z) = partialcorr(hfb.load2(x,e1)-hfb.load1(x,e1), rt_cl, jitter_cl, 'Type', 'Spearman');
    end
    clear x
end

% z-score
hfb.r_z = (hfb.r-nanmean(hfb.r_perm,2))./nanstd(hfb.r_perm,[],2);

save([savdir 'cl_hfb'], 'hfb');

clear all

end
