% linear mixed effects model. iv: plv stat direction, network-network; 
% dv: plv correlation. input freq = 'theta'/'alpha'/'beta'
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_dir(freq)

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
        cf = load([datdir 'cf_' freq], freq);
        cl = load([datdir 'cl_' freq], freq);
        
        ch_idx = ismember(ch, cf.(freq).label);
        yeo_key = yeo_key(ch_idx);

        for e1 = 1:length(cf.(freq).label)
            for e2 = 1:length(cf.(freq).label)
                if e1 > e2
                    c = c+1;
                    sub{c} = sbj{s};
                    elec1{c} = cf.(freq).label{e1};
                    elec2{c} = cf.(freq).label{e2};
                    network{c} = strcat(yeo_key{e1}, '-', yeo_key{e2});
                    network1{c} = yeo_key{e1};
                    network2{c} = yeo_key{e2};
                    cf_z(c) = cf.(freq).r_z(e1,e2);
                    cl_z(c) = cl.(freq).r_z(e1,e2);
                end
            end
        end
    end
    
    clear yeo_key cf cl datdir ch_idx
end

dir = [zeros(length(sub),1); ones(length(sub),1)];
cf_z = [zeros(length(sub),1); cf_z'];
cl_z = [zeros(length(sub),1); cl_z'];
sub = [sub'; sub'];
elec1 = [elec1'; elec1'];
elec2 = [elec2'; elec2'];
network = [network'; network'];
network1 = [network1'; network1'];
network2 = [network2'; network2'];
    
model_data = table(sub, elec1, elec2, network, network1, network2, cf_z, cl_z, dir, ...
    'VariableNames', {'sub', 'elec1', 'elec2', 'network', 'network1', 'network2', 'cf_z', 'cl_z', 'dir'});

cf_lme = fitlme(model_data, 'cf_z ~ dir * network + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
cf_lme_anova = anova(cf_lme);

cl_lme = fitlme(model_data, 'cl_z ~ dir * network + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
cl_lme_anova = anova(cl_lme);

save([savdir 'lme_' freq '_dir'], 'model_data', '*lme*');

clear all

end
