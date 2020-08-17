clear
close all;

% load each csv file
% path = pwd;
% filepattern=fullfile(path,'/*.csv');
% files=dir(filepattern);

% file_name = {'L_STG_ITD_HbO_peak_d.csv', ...
%     'R_STG_ITD_HbO_peak_d.csv', ...
%     'L_STG_opposite_HbO_peak_d.csv', ...
%     'R_STG_opposite_HbO_peak_d.csv'};

file_name = {'R_STG_opposite_HbO_peak_d.csv'};



%% load individual file
% STG region
figure;
hold on;
for index=1:length(file_name)
    
%     base_name = files(index).name;
%     [folder, name, extension] = fileparts(base_name);
%     new_name = [path,'\',name,'.csv'];
%     disp(name);
    
    
    new_name = file_name{index};
    peak_d= csvread(new_name, 1, 1);
    
   
    %3 ROIs%
    %Z score%
    a = peak_d(1, :);
    %MI%
    b = peak_d(2, :);
    
    %Correlation Coefficient%
    r = corrcoef(a,b);
    rsq = r(1,2).^2;

    subplot(2,2,index);
    Mdl1=LinearModel.fit(b,a)
    plot(Mdl1,'marker','o','MarkerSize',8)
    xlim([0,3.5])
    ylim([0,1])
    axis square
    %     xlabel('d prime');
    %     ylabel('model peak');
    title(new_name);
    %title('')
    xlabel('')
    ylabel('')
    
    text(min(b)*1.25,max(a)*0.98,['R^2 = ',num2str(round(Mdl1.Rsquared.Adjusted*100)/100)])
    text(min(b)*1.25,max(a)*0.90,['p = ',num2str(Mdl1.Coefficients.pValue(2))]);
    figure_l = gca; legend(figure_l,'off');
    %saveas(gcf,[name, '.png']);
end


%% load individual file
% cIFS region
% figure;
% hold on;
% for index=1:length(file_name2)
%     
% %     base_name = files(index).name;
% %     [folder, name, extension] = fileparts(base_name);
% %     new_name = [path,'\',name,'.csv'];
% %     disp(name);
%     
%     
%     new_name = file_name2{index};
%     peak_d= csvread(new_name, 1, 1);
%     
%    
%     %3 ROIs%
%     %Z score%
%     a = peak_d(1, :);
%     %MI%
%     b = peak_d(2, :);
%     
%     %Correlation Coefficient%
%     r = corrcoef(a,b);
%     rsq = r(1,2).^2;
% 
%     subplot(2,2,index);
%     Mdl1=LinearModel.fit(b,a);
%     plot(Mdl1,'marker','o','MarkerSize',8)
%     xlim([0,3.5])
%     ylim([0,1])
%     axis square
%     %     xlabel('d prime');
%     %     ylabel('model peak');
%     title(new_name);
%     %title('')
%     xlabel('')
%     ylabel('')
%     
%     text(min(b)*1.25,max(a)*0.98,['R^2 = ',num2str(round(Mdl1.Rsquared.Adjusted*100)/100)])
%     text(min(b)*1.25,max(a)*0.90,['p = ',num2str(Mdl1.Coefficients.pValue(2))]);
%     figure_l = gca; legend(figure_l,'off');
%     %saveas(gcf,[name, '.png']);
% end

%%
% ITD_peak_d = csvread('R_STG_ITD_HbO_peak_d.csv', 1, 1);
% opposite_peak_d =  csvread('R_STG_opposite_HbO_peak_d.csv', 1, 1);
% 
% ITD_d = ITD_peak_d(2, :);
% ITD_peak= ITD_peak_d(1, :);
% opposite_peak = opposite_peak_d(1, :);
% 
% Mdl1=LinearModel.fit(ITD_peak, opposite_peak)
% plot(Mdl1,'marker','o','MarkerSize',8)
