function dc_pieced = piece_blk_new(dc_task, s_task)

dc_pieced = [];

blk_start_index =  find(s_task == 1);

% for each block substract the mean of first 5 seconds before task onset
for ii =  1 : length(blk_start_index)
    
    s = zeros(length(s_task),1);
    s(blk_start_index(ii)) = 1;
    t = [0 0.02];
    window_t = [-5, 45];
    blk = hmrBlockAvg(dc_task,s,t,window_t);
    blk = mean(blk(50 * 5 + 1 : end, :, :), 3);
   
    dc_pieced = [dc_pieced; blk];
    
end
end