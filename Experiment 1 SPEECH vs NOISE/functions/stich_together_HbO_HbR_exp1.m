clear;
close all;

% load csv as table
% HbO 1; HbR 2
HbO_M = csvread('\LMEM R data exp1\R_data_speech_noise_HbO.csv');
HbR_M = csvread('\LMEM R data exp1\R_data_speech_noise_HbR.csv');

% stich HbO HbR horizontally
% dublocate column 11 ~ 16
% make half of 11 ~ 16 zero
% add one column mark HbO or HbR (1,2)
combined = cat(1, HbO_M, HbR_M);
% dublicate column 11 ~ 16
combined_dup = cat(2, combined, combined(:, 10:16));
rows = length(combined_dup);
% fisrt half is HbO, second half is HbR
% make the second half of original 10 to 16 column zero
% make the first half of added 20 to 26 column zero
combined_dup(rows/2 + 1 : end, 10 : 16) =  0;
combined_dup(1 : rows/2, 20 : 26) =  0;
% add label column for HbO 1, HbR 2
label_column = [ones(rows/2, 1); ones(rows/2, 1) * 2];
% final data for R
data_for_R = cat(2, combined_dup, label_column);

% save as cvs
% go up one directory
csvwrite([pwd, '\LMEM R data exp1\R_data_speech_noise_HbO_HbR.csv'], data_for_R);