function [blk_num, condition, behavior]= block_number_R(s_ITD_task, s_noITD_task, ITD_d_prime, noITD_d_prime)

%function to build block num
s_task = s_ITD_task + s_noITD_task;

tick = find(s_task==1);
tick_ITD = find(s_ITD_task == 1);
tick_noITD = find(s_noITD_task == 1);


blk_num = s_task;
condition = s_task;
behavior = s_task;
for idx = 1:numel(tick)
    
    element = tick(idx);
    %disp(element);
    blk_num(element:element+45*50)=idx;
    
    if ismember(element,tick_ITD)
        condition(element:element+45*50)= 1;
        behavior(element:element+45*50)= ITD_d_prime;
    end
    
    if ismember(element,tick_noITD)
        condition(element:element+45*50)= 2;
        behavior(element:element+45*50)= noITD_d_prime;
    end
        
end


end

