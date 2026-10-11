% post hoc testing of accuracy x network-network interaction on plv. tests
% accuracy per network-network. input freq = 'theta'/'alpha'/'beta'
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_acc_int(freq)

pth = ; % INPUT PATH TO DATA
datdir = [pth '/plv/'];
savdir = datdir;

% load model data
load([datdir 'lme_' freq '_acc'], 'model_data');

network1{1} = strcmp(model_data.network1, 'VIS');
network1{2} = strcmp(model_data.network1, 'SM');
network1{3} = strcmp(model_data.network1, 'DAN');
network1{4} = strcmp(model_data.network1, 'VAN');
network1{5} = strcmp(model_data.network1, 'LB');
network1{6} = strcmp(model_data.network1, 'FP');
network1{7} = strcmp(model_data.network1, 'DM');

network2{1} = strcmp(model_data.network2, 'VIS');
network2{2} = strcmp(model_data.network2, 'SM');
network2{3} = strcmp(model_data.network2, 'DAN');
network2{4} = strcmp(model_data.network2, 'VAN');
network2{5} = strcmp(model_data.network2, 'LB');
network2{6} = strcmp(model_data.network2, 'FP');
network2{7} = strcmp(model_data.network2, 'DM');

cf_lme_acc = cell(length(network1), length(network2));
cf_lme_acc_anova = cf_lme_acc;
cf_lme_acc_t = nan(length(network1), length(network2));

cl_lme_acc = cf_lme_acc;
cl_lme_acc_anova = cf_lme_acc;
cl_lme_acc_t = cf_lme_acc_t;

lme_acc_network = cf_lme_acc;

for n1 = 1:length(network1)
    for n2 = 1:length(network2)
        if n1 >= n2
            tmp_idx = network1{n1} + network2{n2} == 2;
            tmp_model = model_data(tmp_idx,:);
            
            lme_acc_network{n1,n2} = strcat(num2str(n1), '-', num2str(n2));
            
            try
                cf_lme_acc{n1,n2} = fitlme(tmp_model, 'cf_m ~ acc + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cf_lme_acc_anova{n1,n2} = anova(cf_lme_acc{n1,n2});
                cf_lme_acc_t(n1,n2) = cf_lme_acc{n1,n2}.Coefficients.tStat(2);
            catch
                cf_lme_acc{n1,n2} = nan;
                cf_lme_acc_anova{n1,n2} = nan;
                cf_lme_acc_t(n1,n2) = nan;
            end
            
            try
                cl_lme_acc{n1,n2} = fitlme(tmp_model, 'cl_m ~ acc + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cl_lme_acc_anova{n1,n2} = anova(cl_lme_acc{n1,n2});
                cl_lme_acc_t(n1,n2) = cl_lme_acc{n1,n2}.Coefficients.tStat(2);
            catch
                cl_lme_acc{n1,n2} = nan;
                cl_lme_acc_anova{n1,n2} = nan;
                cl_lme_acc_t(n1,n2) = nan;
            end

            clear tmp*
        end
    end
end

save([savdir 'lme_' freq '_acc_int'], '*lme*');

clear all

end
