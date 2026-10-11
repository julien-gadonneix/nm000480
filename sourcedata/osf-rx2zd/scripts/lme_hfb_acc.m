% linear mixed effects model. iv: accuracy, network; dv: hfb effect size
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

pth = ; % INPUT PATH TO DATA
savdir = [pth '/hfb/'];
mkdir(savdir);

sbj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

c = 0;
for s = 1:length(sbj)
    datdir = [pth sbj{s} '/'];
    load([datdir 'yeo'], 'yeo_key');
    datdir = [pth sbj{s} '/hfb/'];
    
    cf_cor = load([datdir 'cf_hfb'], 'hfb');
    cf_inc = load([datdir 'cf_hfb_inc'], 'hfb');
    
    cl_cor = load([datdir 'cl_hfb'], 'hfb');
    if ~strcmp(sbj{s},'S3')
        cl_inc = load([datdir 'cl_hfb_inc'], 'hfb');
    end
    
    for e = 1:length(cf_cor.hfb.label)
        c = c+1;
        sub{c} = sbj{s};
        elec{c} = cf_cor.hfb.label{e};
        network{c} = yeo_key{e};
        cf_cor_m(c) = cf_cor.hfb.diff(e);
        cf_inc_m(c) = cf_inc.hfb.diff(e);
        cl_cor_m(c) = cl_cor.hfb.diff(e);
        if ~strcmp(sbj{s},'S3')
            cl_inc_m(c) = cl_inc.hfb.diff(e);
        else
            cl_inc_m(c) = nan;
        end
    end
    
    clear yeo_key *cor *inc datdir
end

acc = [zeros(length(sub),1); ones(length(sub),1)];
cf_m = [cf_inc_m'; cf_cor_m'];
cl_m = [cl_inc_m'; cl_cor_m'];
sub = [sub'; sub'];
elec = [elec'; elec'];
network = [network'; network'];

model_data = table(sub, elec, network, cf_m, cl_m, acc, ...
    'VariableNames', {'sub', 'elec', 'network', 'cf_m', 'cl_m', 'acc'});

cf_lme = fitlme(model_data, 'cf_m ~ acc * network + (1|sub) + (1|sub:elec)');
cf_lme_anova = anova(cf_lme);

cl_lme = fitlme(model_data, 'cl_m ~ acc * network + (1|sub) + (1|sub:elec)');
cl_lme_anova = anova(cl_lme);

save([savdir 'lme_acc'], 'model_data', '*lme*');

clear all
