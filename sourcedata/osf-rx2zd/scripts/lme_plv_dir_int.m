% post hoc testing of direction x network-network interaction on plv. tests
% direction per network-network. input freq = 'theta'/'alpha'/'beta'
%
% example: lme_plv_dir_int('theta')
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function [] = lme_plv_dir_int(freq)

pth = ; % INPUT PATH TO DATA
datdir = [pth '/plv/'];
savdir = datdir;

% load model data
load([datdir 'lme_' freq '_dir'], 'model_data');

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

cf_lme_dir = cell(length(network1), length(network2));
cf_lme_dir_anova = cf_lme_dir;
cf_lme_dir_t = nan(length(network1), length(network2));

cl_lme_dir = cf_lme_dir;
cl_lme_dir_anova = cf_lme_dir;
cl_lme_dir_t = cf_lme_dir_t;

lme_dir_network = cf_lme_dir;

for n1 = 1:length(network1)
    for n2 = 1:length(network2)
        if n1 >= n2
            tmp_idx = network1{n1} + network2{n2} == 2;
            tmp_model = model_data(tmp_idx,:);
            
            lme_dir_network{n1,n2} = strcat(num2str(n1), '-', num2str(n2));
            
            try
                cf_lme_dir{n1,n2} = fitlme(tmp_model, 'cf_z ~ dir + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cf_lme_dir_anova{n1,n2} = anova(cf_lme_dir{n1,n2});
                cf_lme_dir_t(n1,n2) = cf_lme_dir{n1,n2}.Coefficients.tStat(2);
            catch
                cf_lme_dir{n1,n2} = nan;
                cf_lme_dir_anova{n1,n2} = nan;
                cf_lme_dir_t(n1,n2) = nan;
            end
            
            try
                cl_lme_dir{n1,n2} = fitlme(tmp_model, 'cl_z ~ dir + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cl_lme_dir_anova{n1,n2} = anova(cl_lme_dir{n1,n2});
                cl_lme_dir_t(n1,n2) = cl_lme_dir{n1,n2}.Coefficients.tStat(2);
            catch
                cl_lme_dir{n1,n2} = nan;
                cl_lme_dir_anova{n1,n2} = nan;
                cl_lme_dir_t(n1,n2) = nan;
            end

            clear tmp*
        end
    end
end

save([savdir 'lme_' freq '_dir_int'], '*lme*');

clear all

end
