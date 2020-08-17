function [target_hrf_ITD,target_hrf_noITD, time_speech, time_SSN] = target_HRF_cat_R(condition_no_gap, blk_num_no_gap, Hb_type, time_window)

%function to build HRF for task 

data_length = length(condition_no_gap);

switch Hb_type
    case 'HbO'
        RT = 0.02;
        %p = [10 35 1 1 100 0 45];
        p = [6 16 1 1 100 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

        
    case 'HbR'
        RT = 0.02;
        %p = [43 45 1 1 100 0 45];
        p = [6 16 1 1 100 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

    case 'HbT'
        RT = 0.02;
        %p = [5 25 1 1 80 0 45];
        p = [6 16 1 1 80 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

end

% 2 tasks
s_ITD_task_square = zeros(length(blk_num_no_gap),1);
s_noITD_task_square = zeros(length(blk_num_no_gap),1);
time_speech = zeros(length(condition_no_gap),1);
unique_blk_num = unique(blk_num_no_gap);
start_indices = [];
for ii = 1 : length(unique_blk_num)
    blk_start = find(blk_num_no_gap == unique_blk_num(ii), 1, 'first');
    start_indices = [start_indices, blk_start];
end

for idx = 1:length(start_indices)
    element = start_indices(idx);
    %disp(element)
    time_speech(element:element+45*50) =  time_window;
    if condition_no_gap(element) == 1
        s_ITD_task_square(element:element+15*50) = 1;
    end
    if condition_no_gap(element) == 2
        s_noITD_task_square(element:element+15*50) = 1;
    end 
        
end

target_hrf_ITD = conv(h,s_ITD_task_square);
target_hrf_ITD=target_hrf_ITD(1:data_length)/max(abs(target_hrf_ITD));

target_hrf_noITD = conv(h,s_noITD_task_square);
target_hrf_noITD=target_hrf_noITD(1:data_length)/max(abs(target_hrf_noITD));


end

