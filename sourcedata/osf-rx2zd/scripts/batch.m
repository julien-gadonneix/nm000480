%% subject-level analyses

sbj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

for s = 1:length(sbj)
    hfb_analysis(sbj{s});
    hfb_stats(sbj{s});
    hfb_analysis_errors(sbj{s});

    peaks_irasa(sbj{s});
    plv_analysis(sbj{s});
    plv_analysis_errors(sbj{s});
end

% subfunctions: zbaseline_200, zcont_5, split_cf_cl

%% group-level models

lme_hfb_acc;
lme_hfb_dir;

freqs = {'theta','alpha','beta'};

for f = 1:length(freqs)
    lme_plv_acc(freqs{f});
    lme_plv_dir(freqs{f});
    
    lme_plv_graph_acc(freqs{f});
    lme_plv_graph_dir(freqs{f}); 
end

% post hoc models

freqs = {'theta','alpha','beta'};

for f = 1:length(freqs)
    lme_plv_acc_int(freqs{f});    
    lme_plv_dir_int(freqs{f});
    
    lme_plv_graph_acc_int(freqs{f});
    lme_plv_graph_dir_int(freqs{f});
end
