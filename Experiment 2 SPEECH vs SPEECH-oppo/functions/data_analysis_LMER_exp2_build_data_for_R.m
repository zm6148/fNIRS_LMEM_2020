function data_analysis_LMER_exp2_build_data_for_R(Hb_type)%% load data folder
% data analysis using GLM model
%clear;
%close all;

%% load processed fnirs data
path = [pwd, '\fNIRS data exp2\processed data'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

%% behavior data
load([pwd, '\behavioral data exp2\behavioral_data.mat']);

%% listener PTA data and left/right handedness
% left/right hand (left:1; right:2)
L_R_hand = [2
    2
    2
    1
    2
    2
    2
    2
    2
    2
    2
    2
    2
    2];
% PTA data
R_audio = [10
    8
    10
    7
    7
    12
    5
    7
    8
    13
    8
    13
    5
    3];

L_audio = [5
    8
    7
    10
    5
    7
    2
    3
    8
    12
    8
    15
    3
    3];

%% define fnirs parameters
% Difine block average window
window_b=[-5,40];

% which Hb type to process
switch Hb_type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

%% build data for R analysis using LMEM
data_for_R = [];
for index=1:length(files)
    
    % load breathholding data
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    
    
    % loop through ROI 
    for jj = 1 : 4 % 4 regions of interest
        
        % switch different channel number for differnet ROI
        % jj 1:E right cIFS; 2:A left cIFS; 3:F rigth STG; 4:B left STG
        % further divide into left/right, and location (STG, cIFS)
        % left/right: 1:left; 2:right
        % STG : 2; cIFS: 1
        
        % block number
        speech_d_prime = speech_d(index);
        speechoppo_d_prime = speechoppo_d(index);
        [blk_num, condition, behavior] = block_number_R(s_speech_task, s_speechoppo_task, speech_d_prime, speechoppo_d_prime);
        task_only_index =  find(blk_num ~= 0);
        blk_num_no_gap = blk_num(task_only_index);
        condition_no_gap = condition(task_only_index);
        behavior_no_gap = behavior(task_only_index);
        
        switch jj
            case 1 %E 
                RC = 14;
                channels = [11,12,13,15];
                hemsphere = ones(size(blk_num_no_gap,1),1) * 2; %right
                location = ones(size(blk_num_no_gap,1),1) * 1; % cIFS
            case 2 %A
                RC = 4;
                channels = [1,2,3,5];
                hemsphere = ones(size(blk_num_no_gap,1),1) * 1; %left
                location = ones(size(blk_num_no_gap,1),1) * 1; % cIFS
            case 3 %F 
                RC = 17;
                channels = [16,18,19,20];
                hemsphere = ones(size(blk_num_no_gap,1),1) * 2; %right
                location = ones(size(blk_num_no_gap,1),1) * 2; % STG
            case 4 %B 
                RC = 17;
                channels = [6,8,9,10];
                hemsphere = ones(size(blk_num_no_gap,1),1) * 1; %left
                location = ones(size(blk_num_no_gap,1),1) * 2; % STG
        end
        
        % block average at breahhodling part
        breath_blockAvg = hmrBlockAvg(dc(:, :, channels), s_b, t, window_b);

        % max breath used for normalize later
        max_b_HbO = max(abs(mean(breath_blockAvg(:,1,:),3)));
        max_b_HbR = max(abs(mean(breath_blockAvg(:,2,:),3)));
        max_b = (max_b_HbO + max_b_HbR)/2;
        
        % build new dc based on each block pieced togethor
        % piece together observed fnirs data
        fnirs_data_pieced = piece_blk(dc_task(:, number, channels), s_speech_task + s_speechoppo_task)/ max_b;
        
        % piece together RC channel fnirs data
        RC_fnirs_data_pieced = piece_blk(dc_task(:, number, RC), s_speech_task + s_speechoppo_task)/ max_b;
        
        % subject ID
        sid = ones(size(blk_num_no_gap,1),1) * index;
        
        % L R hand
        L_R_hand_S = ones(size(blk_num_no_gap,1),1) * L_R_hand(index);
        
        % R audiometric
        R_audio_S = ones(size(blk_num_no_gap,1),1) * R_audio(index);
        
        % L audiometric
        L_audio_S = ones(size(blk_num_no_gap,1),1) * L_audio(index);
        
        %roi_code E,A,F,B: 1,2,3,4
        roi_code = ones(size(blk_num_no_gap,1), 1) * jj;
        
        % task_type 1:SPEECH; 2: SPEECH-oppo
        task_type  = condition_no_gap;
        
        % first derivative
        % HRF function for SPEECH and SPEECH-oppo
        time_window  = 0 : 1/50 : 45;
        [target_hrf_speech, target_hrf_speechoppo, time] = target_HRF_cat_R(condition_no_gap, blk_num_no_gap, Hb_type, time_window);
        first_deriv_speech = [0; diff(target_hrf_speech)];
        first_deriv_speechoppo = [0; diff(target_hrf_speechoppo)];
        second_deriv_speech = [0; diff(first_deriv_speech)];
        second_deriv_speechoppo = [0; diff(first_deriv_speechoppo)];
        
        dummy = cat(2, sid, hemsphere, location, roi_code, task_type, ...
            L_R_hand_S, R_audio_S, L_audio_S, ...
            fnirs_data_pieced, RC_fnirs_data_pieced, ...
            target_hrf_speech, target_hrf_speechoppo, ...
            first_deriv_speech, first_deriv_speechoppo, ...
            second_deriv_speech, second_deriv_speechoppo, ...
            blk_num_no_gap, behavior_no_gap, time);
        
        data_for_R = cat(1, data_for_R, dummy);
    end
end

% save as cvs
csvwrite([pwd, '\LMEM R data exp2\R_data_speech_speechoppo_', Hb_type, '.csv'], data_for_R);
end
