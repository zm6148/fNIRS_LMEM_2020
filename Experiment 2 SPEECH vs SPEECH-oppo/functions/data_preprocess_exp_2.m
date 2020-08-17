%% load data
% data preprocessing steps
clear all;
close all;

path = [pwd,'\fNIRS data exp2\unprocessed data'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

%% Difine window
window_b=[-5,40];
window_t=[-5,40];

%% load individual file and do processing

for index=1:length(files)
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    
    
    %based on aux(:,1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    stim_b_speech = aux(:,1);
    [pks,locs] = findpeaks(stim_b_speech);
    
    %find peak above certain threshold
    peak_index_aux1=[];
    for ii = 1:length(locs)
        if aux(locs(ii),1)>1
            peak_index_aux1=[peak_index_aux1,locs(ii)];
        end
    end
    %remove close by peaks
    index_to_remove=[];
    for ii = 2:length(peak_index_aux1)
        if peak_index_aux1(ii)-peak_index_aux1(ii-1)<1000
            index_to_remove=[index_to_remove,ii];
        end           
    end
    peak_index_aux1(index_to_remove)=[];
    
    %based on aux(:,2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stim_speechoppo = aux(:,2);
    [pks,locs] = findpeaks(stim_speechoppo);
    
    %find peak above certain threshold
    peak_index_aux2=[];
    for ii = 1:length(locs)
        if aux(locs(ii),2)>1
            peak_index_aux2=[peak_index_aux2,locs(ii)];
        end
    end
    %remove close by peaks
    index_to_remove=[];
    for ii = 2:length(peak_index_aux2)
        if peak_index_aux2(ii)-peak_index_aux2(ii-1)<1000
            index_to_remove=[index_to_remove,ii];
        end           
    end
    peak_index_aux2(index_to_remove)=[];
    
    %set s_b, s_speech, s_speechoppo at peaks to 1
    %first 10 is breath holding, last 10 is task
    s_b=zeros(length(s),1);
    s_speech=zeros(length(s),1);
    s_speechoppo=zeros(length(s),1);
       
    peak_b = peak_index_aux1(1:10);
    peak_speech = peak_index_aux1(end-6+1:end); % change based on subject 11 before 6 now
    peak_speechoppo = peak_index_aux2(end-6+1:end);% change based on subject
 
    s_b(peak_b)=1;
    s_speech(peak_speech)=1;
    s_speechoppo(peak_speechoppo)=1;
    
    % create new S and d data 
    % by cut from 10s before first s till last data point
    % 10s is 10*50=500 data points
    
    s_b=s_b(1:end);
    b_last_index = max(find(s_b==1));
    
    s_speech=s_speech(1:end);
    speech_first_index = min(find(s_speech==1));
    speech_last_index = max(find(s_speech==1));
    
    s_speechoppo=s_speechoppo(1:end);
    speechoppo_first_index = min(find(s_speechoppo==1));
    speechoppo_last_index = max(find(s_speechoppo==1));
    
    first_task_index = min([speech_first_index, speechoppo_first_index]);
    last_task_index = max([speech_last_index, speechoppo_last_index]);
    
    % data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert to dod and filter
    tIncMan=ones(length(s_b),1);
    dod = hmrIntensity2OD(d);
    dod_raw = dod;
    
    % bandpass filter
    [dod,svs,nSV] = enPCAFilter(dod,SD,tIncMan,0);
    % smaller low pass to get smoother curves
    dod = hmrBandpassFilt(dod, t, 0.01, 0.1);
    
    % 20 th ploynomial detrend
    dod = detrend( dod, 20 );
    
    % motion correction 
    % Molavi for motion correction
    iqr = 1;
    dod = hmrMotionCorrectWavelet(dod,SD,iqr);
    
    % convert to dc
    dc = hmrOD2Conc(dod,SD,[6 6]);
    
    % extract task portion--no breath holding data
    % by remove the breath holding part
    start_index = first_task_index-50*5;
    end_index = last_task_index+50*45;
    
    dc_task = dc(start_index:end_index,:,:);
    s_speech_task = s_speech(start_index:end_index);
    s_speechoppo_task = s_speechoppo(start_index:end_index);
    
    % go up one directory
    [upperPath, deepestFolder, ~] = fileparts(path);
        
    % save this subject data
    subject_data_name = [upperPath, '/processed data/', name, '_processed_data.mat'];
    save(subject_data_name,'dc_task','s_speech_task','s_speechoppo_task','dc','s_b','aux','t','first_task_index','last_task_index', 'dod_raw');
    
end

