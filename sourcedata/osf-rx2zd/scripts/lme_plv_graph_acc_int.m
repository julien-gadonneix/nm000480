% post hoc testing of accuracy x network interaction on plv node strength.
% tests accuracy per network. input freq = 'theta'/'alpha'/'beta'
%
% example: lme_plv_graph_acc_int('theta')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_graph_acc_int(freq)

pth = ; % INPUT PATH TO DATA
datdir = [pth '/plv/'];
savdir = datdir;

% load model data
load([datdir 'lme_' freq '_graph_acc'], 'model_data');

network{1} = strcmp(model_data.network, 'VIS');
network{2} = strcmp(model_data.network, 'SM');
network{3} = strcmp(model_data.network, 'DAN');
network{4} = strcmp(model_data.network, 'VAN');
network{5} = strcmp(model_data.network, 'LB');
network{6} = strcmp(model_data.network, 'FP');
network{7} = strcmp(model_data.network, 'DM');

cf_lme_acc = cell(length(network), 1);
cf_lme_acc_anova = cf_lme_acc;
cf_lme_acc_t = nan(length(network), 1);

cl_lme_acc = cf_lme_acc;
cl_lme_acc_anova = cf_lme_acc;
cl_lme_acc_t = cf_lme_acc_t;

lme_acc_network = cf_lme_acc;

for n1 = 1:length(network)
    tmp_model = model_data(network{n1},:);
    
    lme_acc_network{n1} = num2str(n1);
    
    try
        cf_lme_acc{n1} = fitlme(tmp_model, 'cf_m ~ acc + (1|sub) + (1|sub:elec)');
        cf_lme_acc_anova{n1} = anova(cf_lme_acc{n1});
        cf_lme_acc_t(n1) = cf_lme_acc{n1}.Coefficients.tStat(2);
    catch
        cf_lme_acc{n1} = nan;
        cf_lme_acc_anova{n1} = nan;
        cf_lme_acc_t(n1) = nan;
    end
    
    try
        cl_lme_acc{n1} = fitlme(tmp_model, 'cl_m ~ acc + (1|sub) + (1|sub:elec)');
        cl_lme_acc_anova{n1} = anova(cl_lme_acc{n1});
        cl_lme_acc_t(n1) = cl_lme_acc{n1}.Coefficients.tStat(2);
    catch
        cl_lme_acc{n1} = nan;
        cl_lme_acc_anova{n1} = nan;
        cl_lme_acc_t(n1) = nan;
    end
    
    clear tmp*
end

save([savdir 'lme_' freq '_graph_acc_int'], '*lme*');

clear all

end
