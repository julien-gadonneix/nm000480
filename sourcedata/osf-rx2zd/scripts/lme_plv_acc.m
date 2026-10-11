% linear mixed effects model. iv: accuracy, network; dv: plv effect size.
% input freq = 'theta'/'alpha'/'beta'
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_acc(freq)

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
    
    datdir = [pth sbj{s} '/plv/'];
    
    if isfile([datdir 'cf_' freq '.mat'])
        cf_cor = load([datdir 'cf_' freq], freq);
        cf_inc = load([datdir 'cf_' freq '_inc'], freq);

        cl_cor = load([datdir 'cl_' freq], freq);
        if ~strcmp(sbj{s},'S3')
            cl_inc = load([datdir 'cl_' freq '_inc'], freq);
        end
        
        ch_idx = ismember(ch, cf_cor.(freq).label);
        yeo_key = yeo_key(ch_idx);      
        
        for e1 = 1:length(cf_cor.(freq).label)
            for e2 = 1:length(cf_cor.(freq).label)
                if e1 > e2
                    c = c+1;
                    sub{c} = sbj{s};
                    elec1{c} = cf_cor.(freq).label{e1};
                    elec2{c} = cf_cor.(freq).label{e2};
                    network{c} = strcat(yeo_key{e1}, '-', yeo_key{e2});
                    network1{c} = yeo_key{e1};
                    network2{c} = yeo_key{e2};
                    cf_cor_m(c) = cf_cor.(freq).diff(e1,e2);
                    cf_inc_m(c) = cf_inc.(freq).diff(e1,e2);
                    cl_cor_m(c) = cl_cor.(freq).diff(e1,e2);
                    if ~strcmp(sbj{s},'S3')
                        cl_inc_m(c) = cl_inc.(freq).diff(e1,e2);
                    else
                        cl_inc_m(c) = nan;
                    end
                end
            end
        end

        clear yeo_key *cor *inc datdir ch_idx
    end
end

acc = [zeros(length(sub),1); ones(length(sub),1)];
cf_m = [cf_inc_m'; cf_cor_m'];
cl_m = [cl_inc_m'; cl_cor_m'];
sub = [sub'; sub'];
elec1 = [elec1'; elec1'];
elec2 = [elec2'; elec2'];
network = [network'; network'];
network1 = [network1'; network1'];
network2 = [network2'; network2'];

model_data = table(sub, elec1, elec2, network, network1, network2, cf_m, cl_m, acc, ...
    'VariableNames', {'sub', 'elec1', 'elec2', 'network', 'network1', 'network2', 'cf_m', 'cl_m', 'acc'});

cf_lme = fitlme(model_data, 'cf_m ~ acc * network + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
cf_lme_anova = anova(cf_lme);

try
    cl_lme = fitlme(model_data, 'cl_m ~ acc * network + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
    cl_lme_anova = anova(cl_lme);
catch
    hold = model_data;
    model_data(strcmp(model_data.sub,'S3'),:) = [];
    
    cl_lme = fitlme(model_data, 'cl_m ~ acc * network + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
    cl_lme_anova = anova(cl_lme);

    model_data = hold;
end

save([savdir 'lme_' freq '_acc'], 'model_data', '*lme*');

clear all

end
