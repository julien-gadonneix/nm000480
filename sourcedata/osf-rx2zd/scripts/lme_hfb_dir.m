% linear mixed effects model. iv: hfb stat direction, network; dv: hfb
% correlation
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
    
    cf = load([datdir 'cf_hfb'], 'hfb');
    cl = load([datdir 'cl_hfb'], 'hfb');
    
    for e = 1:length(cf.hfb.label)
        c = c+1;
        sub{c} = sbj{s};
        elec{c} = cf.hfb.label{e};
        network{c} = yeo_key{e};
        cf_z(c) = cf.hfb.r_z(e);
        cl_z(c) = cl.hfb.r_z(e);
    end
    
    clear yeo_key cf cl datdir
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

save([savdir 'lme_dir'], 'model_data', '*lme*');

clear all
