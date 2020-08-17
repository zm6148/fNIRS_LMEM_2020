%% analysis on behavioral response
clear;

% load FAR lookup table
load('FAR_table.mat');

order = {'TNB_experiment_1', ...
    'TNE_experiment_1', ...
    'TOY_experiment_1', ...
    'TOZ_experiment_1', ...
    'TPA_experiment_1', ...
    'TPB_experiment_1', ...
    'TPC_experiment_1', ...
    'TPQ_experiment_1', ...
    'TPS_experiment_1', ...
    'TPT_experiment_1', ...
    'TPW_experiment_1', ...
    'TPZ_experiment_1', ...
    'TQA_experiment_1', ...
    'TQE_experiment_1'};

% 2 conditions
% 1: SPEECH
% 2: NOISE

out_put = [];
for ii = 1:length(order)
    name = order{ii};
    
    disp(name);
    load(name);
    
    
    all_block_hitrate_SNR_cal = [];
    for jj = 1 : length(all_response_timestamp)
        % loop through each target play time and try to find a match within
        % n s (1.2 s for now) in the response timestamp
        correct_window = 1.2;
        
        target_timestamp = all_target_timestamp{jj};
        response_timestamp = all_response_timestamp{jj};
        block_total_color = length(target_timestamp);
        block_total_response = length(response_timestamp);
        
        block_correct_total = 0;
        
        for k = 1 : size(target_timestamp, 1)
            target_time = target_timestamp(k,2);
            % loop through each response time to find a match
            for l = 1 : size(response_timestamp,1)
                response_time = response_timestamp(l,3);
                % if within the window
                if (target_time <= response_time) && (response_time <= (target_time + correct_window))
                    block_correct_total = block_correct_total + 1;
                    response_timestamp(l,:) = [];
                    break;
                end
            end
        end
        
        % calculate hitrate
        hitrate = block_correct_total/block_total_color;
        %calcualte d'
        %find hits and false alarm
        %simulated farate, depends on the number of keywords and button
        %pushes for that block and listener

        hits = block_correct_total;
        switch block_total_color
            case 3
                FAR_index= 1;
            case 4
                FAR_index = 2;
            case 5
                FAR_index = 3;
        end
        if block_total_response == 0
            farate = 0;
        else
            if block_total_response > 9
                block_total_response = 9;
            end
            farate = FAR(FAR_index, block_total_response);
        end
        Ntrials = block_total_color;
        
        % calcualte d'
        hitrate = hits/Ntrials;% observed percent correct for that block and listener
        % clipping
        hitrate(hitrate > .999) = .999;
        hitrate(hitrate < .001) = .001;
        farate(farate > .999) = .999;
        farate(farate < .001) = .001;
        
        zh = erfinv(1-2.*(1-hitrate)).*sqrt(2);%z-score of the observed hit rate
        zf = erfinv(1-2.*(1-farate)).*sqrt(2);% z-score of the simulated false alarm rate
        dprime = zh - zf;%d prime score for that listener and block
        beta = -(zh+zf)/2;% bias for that listener and block
        
        % total color and total hits
        %      correct                                d'      hits                 total color        total response
        temp =[block_correct_total/block_total_color, dprime, block_correct_total, block_total_color, block_total_response];
        all_block_hitrate_SNR_cal = [all_block_hitrate_SNR_cal; temp];
        
    end
    
    index_1 = find(condition_sequence_rand == 1);
    index_2 = find(condition_sequence_rand == 2);
    
    % select from all_block_hitrate_SNR
    % correct%  d'  hits    total color   total response
    condition_speech = mean(all_block_hitrate_SNR_cal(index_1,:),1);
    condition_noise = mean(all_block_hitrate_SNR_cal(index_2,:),1);
    
    % final matrix
    final_matrix = cat(1, condition_speech,  condition_noise);
    temp = [condition_speech(1:2), condition_noise(1:2)];
    out_put = [out_put; temp];
end

speech_d = out_put(:, 2);
noise_d = out_put(:, 4);
save([pwd, '\behavioral data exp1\behavioral_data.mat'], 'speech_d', 'noise_d');
