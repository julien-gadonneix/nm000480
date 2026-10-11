% linear mixed effects model. iv: accuracy, network; dv: plv effect size.
% input freq = 'theta'/'alpha'/'beta'
%
% example: lme_plv_graph_acc('theta')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_graph_acc(freq)

pth = ; % INPUT PATH TO DATA
savdir = [pth '/plv/'];
mkdir(savdir);

sbj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

c = 0;
for s = 1:length(sbj)
    datdir = [pth sbj{s} '/'];
    load([datdir 'yeo'], 'yeo_key');
    load([datdir 'hfb/cf_hfb']);
    ch = hfb.label; clear hfb
    
    datdir = [pth sbj{s} '/revision/plv/'];
    
    if isfile([datdir 'cf_' freq '.mat'])
        cf_cor = load([datdir 'cf_' freq '_graph'], freq);
        cf_inc = load([datdir 'cf_' freq '_graph_inc'], freq);

        cl_cor = load([datdir 'cl_' freq '_graph'], freq);
        if ~strcmp(sbj{s},'S3')
            cl_inc = load([datdir 'cl_' freq '_graph_inc'], freq);
        end
        
        ch_idx = ismember(ch, cf_cor.(freq).label);
        yeo_key = yeo_key(ch_idx);      
        
        for e = 1:length(cf_cor.(freq).label)
            c = c+1;
            sub{c} = sbj{s};
            elec{c} = cf_cor.(freq).label{e};
            network{c} = yeo_key{e};
            cf_cor_m(c) = cf_cor.(freq).diff_sum(e);
            cf_inc_m(c) = cf_inc.(freq).diff_sum(e);
            cl_cor_m(c) = cl_cor.(freq).diff_sum(e);
            if ~strcmp(sbj{s},'S3')
                cl_inc_m(c) = cl_inc.(freq).diff_sum(e);
            else
                cl_inc_m(c) = nan;
            end
        end

        clear yeo_key *cor *inc datdir ch_idx
    end
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

save([savdir 'lme_' freq '_graph_acc'], 'model_data', '*lme*');

clear all

end
