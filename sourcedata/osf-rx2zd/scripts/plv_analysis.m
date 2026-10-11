% compute plv between task-responsive electrodes at oscillatory peaks. run
% stats (effect size, correlation), generate chance distributions, and 
% z-score. compute node strength statistics
%
% example: plv_analysis('S1')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = plv_analysis(subid)

pth = ; % INPUT PATH TO DATA
datdir = [pth subid '/'];
savdir = [datdir 'plv/'];
mkdir(savdir);

nperm = 1000;

% load data
load([datdir subid], 'data');

srate = data.fsample;

cor = cell2mat(data.trialinfo(:,6));
trialinfo = data.trialinfo;

% load peaks
load([datdir 'irasa/osci'], 'osci');

cfg = [];
cfg.trials = find(cor == 1);
cfg.channel = osci.label;

data = ft_selectdata(cfg, data);
data.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo 

times = str2double(data.trialinfo(:,16:17));

% 1000-ms data time windows (in sec, 500-ms on screen + 250-ms pre/post fixation)
s1 = [-0.25 0.75]; % 0 is onset
s2 = [times(:,1)/srate-0.25 times(:,1)/srate+0.75];
s3 = [times(:,2)/srate-0.25 times(:,2)/srate+0.75];

cfg = [];
cfg.latency = [min(s1)-1 max(s3(:,2))+1]; % 1-s buffers

data = ft_selectdata(cfg, data);

% set up filter
cfg = [];
cfg.pad = 'nextpow2';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.toi = data.time{1}(1):0.005:data.time{1}(end);

% filter at peaks
theta = [];
theta.dimord = 'rpt_chan_time';
theta.trialinfo = data.trialinfo;
alpha = theta;
beta = theta;

for e = 1:length(osci.label)
    if ~isnan(osci.theta(e))
        cfg.channel = osci.label(e);
        cfg.foi = osci.theta(e);
        cfg.tapsmofrq = ceil(cfg.foi/4);
        cfg.t_ftimwin = 3/cfg.foi;
        
        tmp = ft_freqanalysis(cfg, data);
        
        if ~isfield(theta, 'label')
            theta.label = tmp.label;
            theta.time = tmp.time;
            theta.trial = tmp.fourierspctrm;
            theta.freq = osci.theta(e);
        else
            theta.label = cat(1, theta.label, tmp.label);
            theta.trial = cat(2, theta.trial, tmp.fourierspctrm);
            theta.freq = cat(1, theta.freq, osci.theta(e));
        end
        clear tmp
    end
    
    if ~isnan(osci.alpha(e))
        cfg.channel = osci.label(e);
        cfg.foi = osci.alpha(e);
        cfg.tapsmofrq = ceil(cfg.foi/4); % 1/4 bandwidth
        cfg.t_ftimwin = 3/cfg.foi;
        
        tmp = ft_freqanalysis(cfg, data);
        
        if ~isfield(alpha, 'label')
            alpha.label = tmp.label;
            alpha.time = tmp.time;
            alpha.trial = tmp.fourierspctrm;
            alpha.freq = osci.alpha(e);
        else
            alpha.label = cat(1, alpha.label, tmp.label);
            alpha.trial = cat(2, alpha.trial, tmp.fourierspctrm);
            alpha.freq = cat(1, alpha.freq, osci.alpha(e));
        end
        clear tmp
    end
    
    if ~isnan(osci.beta(e))
        cfg.channel = osci.label(e);
        cfg.foi = osci.beta(e);
        cfg.tapsmofrq = ceil(cfg.foi/4); % 1/4 bandwidth
        cfg.t_ftimwin = 3/cfg.foi;
        
        tmp = ft_freqanalysis(cfg, data);
        
        if ~isfield(beta, 'label')
            beta.label = tmp.label;
            beta.time = tmp.time;
            beta.trial = tmp.fourierspctrm;
            beta.freq = osci.beta(e);
        else
            beta.label = cat(1, beta.label, tmp.label);
            beta.trial = cat(2, beta.trial, tmp.fourierspctrm);
            beta.freq = cat(1, beta.freq, osci.beta(e));
        end
        clear tmp
    end
end
clear data osci

if isfield(theta, 'trial')
    theta.trial = reshape(theta.trial, [size(theta.trial,1) size(theta.trial,2) size(theta.trial,4)]);

    % select stim-locked time windows
    dat_s1 = theta;
    dat_s1.time = s1(1):0.005:s1(2);
    dat_s1.trial = nan(size(theta.trial,1), size(theta.trial,2), length(dat_s1.time));
    dat_s1.cfg = 'dummy';
    dat_s2 = dat_s1; dat_s3 = dat_s1;

    for r = 1:size(theta.trial,1)
        tmp_on = nearest(theta.time, s1(1)); tmp_off = nearest(theta.time, s1(2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s1.trial(r,:,:) = theta.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(theta.time, s2(r,1)); tmp_off = nearest(theta.time, s2(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s2.trial(r,:,:) = theta.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(theta.time, s3(r,1)); tmp_off = nearest(theta.time, s3(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s3.trial(r,:,:) = theta.trial(r,:,tmp_on:tmp_off); clear tmp*
    end
    clear theta

    % concatenate by condition
    [target, distractor, load1, load2] = split_cf_cl(dat_s1, dat_s2, dat_s3);
    clear dat_*
    
    % compute single-trial plv and stats
    t1 = nearest(target.time, 0);
    t2 = length(target.time);
    rt_cf = cell2mat(target.trialinfo(:,5));
    rt_cl = cell2mat(load2.trialinfo(:,5));
    jitter_cf = str2double(target.trialinfo(:,17))-str2double(target.trialinfo(:,16));
    jitter_cl = str2double(load2.trialinfo(:,17))-str2double(load2.trialinfo(:,16));
    
    % context first
    if length(target.label) > 1
        tmp1 = nan(size(target.trial,1), size(target.trial,2), size(target.trial,2), size(target.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(target.label)
            for e2 = 1:length(target.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(target.trial(:,e1,t1:t2))-angle(target.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(distractor.trial(:,e1,t1:t2))-angle(distractor.trial(:,e2,t1:t2));
                end
            end
        end

        theta = [];
        theta.dimord = 'rpt_chan_chan';
        theta.trialinfo = target.trialinfo;
        theta.label = target.label;
        theta.target = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        theta.distractor = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        theta.diff = squeeze(nanmean(theta.target-theta.distractor,1));
        theta.diff = theta.diff';
                
        % true correlation
        theta.r = nan(length(theta.label),length(theta.label)); 
        for e1 = 1:size(theta.r,1)
            for e2 = 1:size(theta.r,2)
                if e1 ~= e2
                    theta.r(e1,e2) = partialcorr(theta.target(:,e1,e2)-theta.distractor(:,e1,e2), ...
                        rt_cf, jitter_cf, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        theta.r_perm = nan(length(theta.label),length(theta.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cf));
            for e1 = 1:size(theta.r,1)
                for e2 = 1:size(theta.r,2)
                    if e1 ~= e2
                        theta.r_perm(e1,e2,z) = partialcorr(theta.target(x,e1,e2)-theta.distractor(x,e1,e2), ...
                            rt_cf, jitter_cf, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        theta.r_z = (theta.r-nanmean(theta.r_perm,3))./nanstd(theta.r_perm,[],3);

        save([savdir 'cf_theta'], 'theta'); 
        clear target distractor
        
        % node strength stats
        hold = theta;
        
        theta = rmfield(theta, {'diff', 'r', 'r_perm', 'r_z'});
        theta.dimord = 'rpt_chan';
        theta.target = nansum(hold.target,3);
        theta.distractor = nansum(hold.distractor,3);
        theta.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        theta.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cf_theta_graph'], 'theta'); 
        clear theta hold
    end
    
    % context last
    if length(load2.label) > 1
        tmp1 = nan(size(load2.trial,1), size(load2.trial,2), size(load2.trial,2), size(load2.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(load2.label)
            for e2 = 1:length(load2.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(load2.trial(:,e1,t1:t2))-angle(load2.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(load1.trial(:,e1,t1:t2))-angle(load1.trial(:,e2,t1:t2));
                end
            end
        end

        theta = [];
        theta.dimord = 'rpt_chan_chan';
        theta.trialinfo = load2.trialinfo;
        theta.label = load2.label;
        theta.load2 = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        theta.load1 = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        theta.diff = squeeze(nanmean(theta.load2-theta.load1,1));
        theta.diff = theta.diff';
        
        % true correlation
        theta.r = nan(length(theta.label),length(theta.label)); 
        for e1 = 1:size(theta.r,1)
            for e2 = 1:size(theta.r,2)
                if e1 ~= e2
                    theta.r(e1,e2) = partialcorr(theta.load2(:,e1,e2)-theta.load1(:,e1,e2), ...
                        rt_cl, jitter_cl, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        theta.r_perm = nan(length(theta.label),length(theta.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cl));
            for e1 = 1:size(theta.r,1)
                for e2 = 1:size(theta.r,2)
                    if e1 ~= e2
                        theta.r_perm(e1,e2,z) = partialcorr(theta.load2(x,e1,e2)-theta.load1(x,e1,e2), ...
                            rt_cl, jitter_cl, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        theta.r_z = (theta.r-nanmean(theta.r_perm,3))./nanstd(theta.r_perm,[],3);

        save([savdir 'cl_theta'], 'theta'); 
        clear load2 load1
        
        % node strength stats
        hold = theta;
        
        theta = rmfield(theta, {'diff', 'r', 'r_perm', 'r_z'});
        theta.dimord = 'rpt_chan';
        theta.load2 = nansum(hold.load2,3);
        theta.load1 = nansum(hold.load1,3);
        theta.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        theta.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cl_theta_graph'], 'theta'); 
        clear theta hold
    end
end

if isfield(alpha, 'trial')
    alpha.trial = reshape(alpha.trial, [size(alpha.trial,1) size(alpha.trial,2) size(alpha.trial,4)]);

    % select stim-locked time windows
    dat_s1 = alpha;
    dat_s1.time = s1(1):0.005:s1(2);
    dat_s1.trial = nan(size(alpha.trial,1), size(alpha.trial,2), length(dat_s1.time));
    dat_s1.cfg = 'dummy';
    dat_s2 = dat_s1; dat_s3 = dat_s1;

    for r = 1:size(alpha.trial,1)
        tmp_on = nearest(alpha.time, s1(1)); tmp_off = nearest(alpha.time, s1(2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s1.trial(r,:,:) = alpha.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(alpha.time, s2(r,1)); tmp_off = nearest(alpha.time, s2(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s2.trial(r,:,:) = alpha.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(alpha.time, s3(r,1)); tmp_off = nearest(alpha.time, s3(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s3.trial(r,:,:) = alpha.trial(r,:,tmp_on:tmp_off); clear tmp*
    end
    clear alpha

    % concatenate by condition
    [target, distractor, load1, load2] = split_cf_cl(dat_s1, dat_s2, dat_s3);
    clear dat_*
    
    % compute single-trial plv and stats
    t1 = nearest(target.time, 0);
    t2 = length(target.time);
    rt_cf = cell2mat(target.trialinfo(:,5));
    rt_cl = cell2mat(load2.trialinfo(:,5));
    jitter_cf = str2double(target.trialinfo(:,17))-str2double(target.trialinfo(:,16));
    jitter_cl = str2double(load2.trialinfo(:,17))-str2double(load2.trialinfo(:,16));
    
    % context first
    if length(target.label) > 1
        tmp1 = nan(size(target.trial,1), size(target.trial,2), size(target.trial,2), size(target.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(target.label)
            for e2 = 1:length(target.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(target.trial(:,e1,t1:t2))-angle(target.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(distractor.trial(:,e1,t1:t2))-angle(distractor.trial(:,e2,t1:t2));
                end
            end
        end

        alpha = [];
        alpha.dimord = 'rpt_chan_chan';
        alpha.trialinfo = target.trialinfo;
        alpha.label = target.label;
        alpha.target = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        alpha.distractor = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        alpha.diff = squeeze(nanmean(alpha.target-alpha.distractor,1));
        alpha.diff = alpha.diff';
                
        % true correlation
        alpha.r = nan(length(alpha.label),length(alpha.label)); 
        for e1 = 1:size(alpha.r,1)
            for e2 = 1:size(alpha.r,2)
                if e1 ~= e2
                    alpha.r(e1,e2) = partialcorr(alpha.target(:,e1,e2)-alpha.distractor(:,e1,e2), ...
                        rt_cf, jitter_cf, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        alpha.r_perm = nan(length(alpha.label),length(alpha.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cf));
            for e1 = 1:size(alpha.r,1)
                for e2 = 1:size(alpha.r,2)
                    if e1 ~= e2
                        alpha.r_perm(e1,e2,z) = partialcorr(alpha.target(x,e1,e2)-alpha.distractor(x,e1,e2), ...
                            rt_cf, jitter_cf, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        alpha.r_z = (alpha.r-nanmean(alpha.r_perm,3))./nanstd(alpha.r_perm,[],3);

        save([savdir 'cf_alpha'], 'alpha'); 
        clear target distractor
        
        % node strength stats
        hold = alpha;
        
        alpha = rmfield(alpha, {'diff', 'r', 'r_perm', 'r_z'});
        alpha.dimord = 'rpt_chan';
        alpha.target = nansum(hold.target,3);
        alpha.distractor = nansum(hold.distractor,3);
        alpha.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        alpha.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cf_alpha_graph'], 'alpha'); 
        clear alpha hold
    end
    
    % context last
    if length(load2.label) > 1
        tmp1 = nan(size(load2.trial,1), size(load2.trial,2), size(load2.trial,2), size(load2.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(load2.label)
            for e2 = 1:length(load2.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(load2.trial(:,e1,t1:t2))-angle(load2.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(load1.trial(:,e1,t1:t2))-angle(load1.trial(:,e2,t1:t2));
                end
            end
        end

        alpha = [];
        alpha.dimord = 'rpt_chan_chan';
        alpha.trialinfo = load2.trialinfo;
        alpha.label = load2.label;
        alpha.load2 = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        alpha.load1 = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        alpha.diff = squeeze(nanmean(alpha.load2-alpha.load1,1));
        alpha.diff = alpha.diff';
        
        % true correlation
        alpha.r = nan(length(alpha.label),length(alpha.label)); 
        for e1 = 1:size(alpha.r,1)
            for e2 = 1:size(alpha.r,2)
                if e1 ~= e2
                    alpha.r(e1,e2) = partialcorr(alpha.load2(:,e1,e2)-alpha.load1(:,e1,e2), ...
                        rt_cl, jitter_cl, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        alpha.r_perm = nan(length(alpha.label),length(alpha.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cl));
            for e1 = 1:size(alpha.r,1)
                for e2 = 1:size(alpha.r,2)
                    if e1 ~= e2
                        alpha.r_perm(e1,e2,z) = partialcorr(alpha.load2(x,e1,e2)-alpha.load1(x,e1,e2), ...
                            rt_cl, jitter_cl, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        alpha.r_z = (alpha.r-nanmean(alpha.r_perm,3))./nanstd(alpha.r_perm,[],3);

        save([savdir 'cl_alpha'], 'alpha'); 
        clear load2 load1
        
        % node strength stats
        hold = alpha;
        
        alpha = rmfield(alpha, {'diff', 'r', 'r_perm', 'r_z'});
        alpha.dimord = 'rpt_chan';
        alpha.load2 = nansum(hold.load2,3);
        alpha.load1 = nansum(hold.load1,3);
        alpha.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        alpha.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cl_alpha_graph'], 'alpha'); 
        clear alpha hold
    end
end

if isfield(beta, 'trial')
    beta.trial = reshape(beta.trial, [size(beta.trial,1) size(beta.trial,2) size(beta.trial,4)]);

    % select stim-locked time windows
    dat_s1 = beta;
    dat_s1.time = s1(1):0.005:s1(2);
    dat_s1.trial = nan(size(beta.trial,1), size(beta.trial,2), length(dat_s1.time));
    dat_s1.cfg = 'dummy';
    dat_s2 = dat_s1; dat_s3 = dat_s1;

    for r = 1:size(beta.trial,1)
        tmp_on = nearest(beta.time, s1(1)); tmp_off = nearest(beta.time, s1(2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s1.trial(r,:,:) = beta.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(beta.time, s2(r,1)); tmp_off = nearest(beta.time, s2(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s2.trial(r,:,:) = beta.trial(r,:,tmp_on:tmp_off); clear tmp*
        tmp_on = nearest(beta.time, s3(r,1)); tmp_off = nearest(beta.time, s3(r,2));
        if tmp_off-tmp_on ~= length(dat_s1.time)-1
            tmp_dif = (length(dat_s1.time)-1) - (tmp_off-tmp_on);
            tmp_off = tmp_off+tmp_dif;
        end
        dat_s3.trial(r,:,:) = beta.trial(r,:,tmp_on:tmp_off); clear tmp*
    end
    clear beta

    % concatenate by condition
    [target, distractor, load1, load2] = split_cf_cl(dat_s1, dat_s2, dat_s3);
    clear dat_*
    
    % compute single-trial plv and stats
    t1 = nearest(target.time, 0);
    t2 = length(target.time);
    rt_cf = cell2mat(target.trialinfo(:,5));
    rt_cl = cell2mat(load2.trialinfo(:,5));
    jitter_cf = str2double(target.trialinfo(:,17))-str2double(target.trialinfo(:,16));
    jitter_cl = str2double(load2.trialinfo(:,17))-str2double(load2.trialinfo(:,16));
    
    % context first
    if length(target.label) > 1
        tmp1 = nan(size(target.trial,1), size(target.trial,2), size(target.trial,2), size(target.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(target.label)
            for e2 = 1:length(target.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(target.trial(:,e1,t1:t2))-angle(target.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(distractor.trial(:,e1,t1:t2))-angle(distractor.trial(:,e2,t1:t2));
                end
            end
        end

        beta = [];
        beta.dimord = 'rpt_chan_chan';
        beta.trialinfo = target.trialinfo;
        beta.label = target.label;
        beta.target = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        beta.distractor = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        beta.diff = squeeze(nanmean(beta.target-beta.distractor,1));
        beta.diff = beta.diff';
        
        % true correlation
        beta.r = nan(length(beta.label),length(beta.label)); 
        for e1 = 1:size(beta.r,1)
            for e2 = 1:size(beta.r,2)
                if e1 ~= e2
                    beta.r(e1,e2) = partialcorr(beta.target(:,e1,e2)-beta.distractor(:,e1,e2), ...
                        rt_cf, jitter_cf, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        beta.r_perm = nan(length(beta.label),length(beta.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cf));
            for e1 = 1:size(beta.r,1)
                for e2 = 1:size(beta.r,2)
                    if e1 ~= e2
                        beta.r_perm(e1,e2,z) = partialcorr(beta.target(x,e1,e2)-beta.distractor(x,e1,e2), ...
                            rt_cf, jitter_cf, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        beta.r_z = (beta.r-nanmean(beta.r_perm,3))./nanstd(beta.r_perm,[],3);

        save([savdir 'cf_beta'], 'beta'); 
        clear target distractor
        
        % node strength stats
        hold = beta;
        
        beta = rmfield(beta, {'diff', 'r', 'r_perm', 'r_z'});
        beta.dimord = 'rpt_chan';
        beta.target = nansum(hold.target,3);
        beta.distractor = nansum(hold.distractor,3);
        beta.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        beta.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cf_beta_graph'], 'beta'); 
        clear beta hold
    end
    
    % context last
    if length(load2.label) > 1
        tmp1 = nan(size(load2.trial,1), size(load2.trial,2), size(load2.trial,2), size(load2.trial(:,:,t1:t2),3)); 
        tmp2 = tmp1;

        for e1 = 1:length(load2.label)
            for e2 = 1:length(load2.label)
                if e1 ~= e2
                    tmp1(:,e1,e2,:) = angle(load2.trial(:,e1,t1:t2))-angle(load2.trial(:,e2,t1:t2));
                    tmp2(:,e1,e2,:) = angle(load1.trial(:,e1,t1:t2))-angle(load1.trial(:,e2,t1:t2));
                end
            end
        end

        beta = [];
        beta.dimord = 'rpt_chan_chan';
        beta.trialinfo = load2.trialinfo;
        beta.label = load2.label;
        beta.load2 = atanh(abs(mean(exp(1i.*(tmp1)),4))); % plv
        beta.load1 = atanh(abs(mean(exp(1i.*(tmp2)),4)));
        clear tmp*
        
        % mean difference
        beta.diff = squeeze(nanmean(beta.load2-beta.load1,1));
        beta.diff = beta.diff';
        
        % true correlation
        beta.r = nan(length(beta.label),length(beta.label)); 
        for e1 = 1:size(beta.r,1)
            for e2 = 1:size(beta.r,2)
                if e1 ~= e2
                    beta.r(e1,e2) = partialcorr(beta.load2(:,e1,e2)-beta.load1(:,e1,e2), ...
                        rt_cl, jitter_cl, 'Type', 'Spearman');
                end
            end
        end
        
        % chance distribution
        beta.r_perm = nan(length(beta.label),length(beta.label),nperm);
        for z = 1:nperm
            x = randperm(length(rt_cl));
            for e1 = 1:size(beta.r,1)
                for e2 = 1:size(beta.r,2)
                    if e1 ~= e2
                        beta.r_perm(e1,e2,z) = partialcorr(beta.load2(x,e1,e2)-beta.load1(x,e1,e2), ...
                            rt_cl, jitter_cl, 'Type', 'Spearman');
                    end
                end
            end
            clear x
        end

        % z-score
        beta.r_z = (beta.r-nanmean(beta.r_perm,3))./nanstd(beta.r_perm,[],3);

        save([savdir 'cl_beta'], 'beta'); 
        clear load2 load1
        
        % node strength stats
        hold = beta;
        
        beta = rmfield(beta, {'diff', 'r', 'r_perm', 'r_z'});
        beta.dimord = 'rpt_chan';
        beta.load2 = nansum(hold.load2,3);
        beta.load1 = nansum(hold.load1,3);
        beta.diff_sum = nansum(hold.diff,2);
        
        tmp_r = nansum(hold.r,2);
        tmp_perm = nansum(hold.r_perm,2);
        beta.r_sum = (tmp_r-nanmean(tmp_perm,3))./nanstd(tmp_perm,[],3);
        clear tmp*
                
        save([savdir 'cl_beta_graph'], 'beta'); 
        clear beta hold
    end
end

clear all

end
