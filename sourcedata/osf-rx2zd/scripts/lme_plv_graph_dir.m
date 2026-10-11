% linear mixed effects model. iv: plv stat direction, network; dv: plv
% correlation. input freq = 'theta'/'alpha'/'beta'
%
% example: lme_plv_graph_dir('theta')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_graph_dir(freq)

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
    
    if isfile([datdir 'cf_' freq '_graph.mat'])
        cf = load([datdir 'cf_' freq '_graph'], freq);
        cl = load([datdir 'cl_' freq '_graph'], freq);
        
        ch_idx = ismember(ch, cf.(freq).label);
        yeo_key = yeo_key(ch_idx);

        for e = 1:length(cf.(freq).label)
            c = c+1;
            sub{c} = sbj{s};
            elec{c} = cf.(freq).label{e};
            network{c} = yeo_key{e};
            cf_z(c) = cf.(freq).r_sum(e);
            cl_z(c) = cl.(freq).r_sum(e);
        end
    end
    
    clear yeo_key cf cl datdir ch_idx
end

dir = [zeros(length(sub),1); ones(length(sub),1)];
cf_z = [zeros(length(sub),1); cf_z'];
cl_z = [zeros(length(sub),1); cl_z'];
sub = [sub'; sub'];
elec = [elec'; elec'];
network = [network'; network'];

model_data = table(sub, elec, network, cf_z, cl_z, dir, ...
    'VariableNames', {'sub', 'elec', 'network', 'cf_z', 'cl_z', 'dir'});

cf_lme = fitlme(model_data, 'cf_z ~ dir * network + (1|sub) + (1|sub:elec)');
cf_lme_anova = anova(cf_lme);

cl_lme = fitlme(model_data, 'cl_z ~ dir * network + (1|sub) + (1|sub:elec)');
cl_lme_anova = anova(cl_lme);

save([savdir 'lme_' freq '_graph_dir'], 'model_data', '*lme*');

clear all

end
