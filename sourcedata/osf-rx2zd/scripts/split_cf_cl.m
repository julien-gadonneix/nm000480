% subfunction to split cf and cl trials
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [cf_target, cf_distractor, cl_load1, cl_load2] = split_cf_cl(dat_s1, dat_s2, dat_s3)

if strcmp(dat_s1.dimord, 'rpt_chan_time')
    param = 'trial';
    if isfield(dat_s1, 'cumtapcnt')
        dat_s1 = rmfield(dat_s1, 'cumtapcnt');
        dat_s2 = rmfield(dat_s2, 'cumtapcnt');
        dat_s3 = rmfield(dat_s3, 'cumtapcnt');
    end
else
    param = 'powspctrm';
end

context_dat = dat_s1.trialinfo(:,2);
cf_idx = cell2mat(context_dat) == 'CF';
cl_idx = cf_idx(:,2) == 0;
cf_idx = cf_idx(:,2);

% context first
if sum(cf_idx) > 0
    target_pos = cell2mat(dat_s1.trialinfo(:,14));
    t2_idx = target_pos == '2';
    t3_idx = target_pos == '3';

    t2_cf_idx = t2_idx + cf_idx == 2;
    t3_cf_idx = t3_idx + cf_idx == 2;

    distract_pos = cell2mat(dat_s1.trialinfo(:,15));
    d2_idx = distract_pos == '2';
    d3_idx = distract_pos == '3';

    d2_cf_idx = d2_idx + cf_idx == 2;
    d3_cf_idx = d3_idx + cf_idx == 2;

    clear *dat *pos

    cfg = []; cfg2 = [];

    cfg.trials = find(t2_cf_idx==1);
    if ~isempty(cfg.trials)
        one = ft_selectdata(cfg, dat_s2);
    else
        one = dat_s2;
        one.(param) = [];
    end
    
    cfg2.trials = find(t3_cf_idx==1);
    if ~isempty(cfg2.trials)
        two = ft_selectdata(cfg2, dat_s3);
    else
        two = dat_s3;
        two.(param) = [];
    end
    
    cf_target = dat_s1;
    cf_target.(param) = cat(1, one.(param), two.(param)); clear one two
    cf_target.trialinfo = dat_s1.trialinfo([cfg.trials; cfg2.trials],:);

    cfg.trials = find(d2_cf_idx==1);
    if ~isempty(cfg.trials)
        one = ft_selectdata(cfg, dat_s2);
    else
        one = dat_s2;
        one.(param) = [];
    end
    
    cfg2.trials = find(d3_cf_idx==1);
    if ~isempty(cfg2.trials)
        two = ft_selectdata(cfg2, dat_s3);
    else
        two = dat_s3;
        two.(param) = [];
    end
    
    cf_distractor = dat_s1;
    cf_distractor.(param) = cat(1, one.(param), two.(param)); clear one two
    cf_distractor.trialinfo = dat_s1.trialinfo([cfg.trials; cfg2.trials],:);

    clear t2* t3* d2* d3* cfg2
end

% context last
if sum(cl_idx) > 0
    cfg.trials = find(cl_idx==1);
    one = ft_selectdata(cfg, dat_s1);
    two = ft_selectdata(cfg, dat_s2);

    cl_load1 = dat_s1;
    cl_load1.(param) = one.(param); clear one
    cl_load1.trialinfo = dat_s1.trialinfo(cfg.trials,:);

    cl_load2 = dat_s2;
    cl_load2.(param) = two.(param); clear two
    cl_load2.trialinfo = dat_s2.trialinfo(cfg.trials,:);
end

clearvars -except cf_target cf_distractor cl_load1 cl_load2

end
