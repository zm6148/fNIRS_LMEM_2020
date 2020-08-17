clear;

%% preprocess fnirs data
run('data_preprocess_exp_1.m');

%% process behavior data
run('reponse_analysis_calculate_exp1.m');

%% generate data for R
% HbO
data_analysis_LMER_exp1_build_data_for_R('HbO');
% HbR
data_analysis_LMER_exp1_build_data_for_R('HbR');

%% combine HbO, HbR to generate data for R
run('stich_together_HbO_HbR_exp1.m');
