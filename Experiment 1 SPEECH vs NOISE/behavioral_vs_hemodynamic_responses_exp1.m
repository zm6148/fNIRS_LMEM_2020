clear
close all;

%% load csv file for left and right STG under SPEECH and NOISE condition
file_name = {'L_STG_speech_HbO_peak_d.csv', ...
    'R_STG_speech_HbO_peak_d.csv', ...
    'L_STG_noise_HbO_peak_d.csv', ...
    'R_STG_noise_HbO_peak_d.csv'};
 
%% load individual file
for index=1:length(file_name)

    subject_name = file_name{index};
    peak_d= csvread(['behavioral and hemodynamic responses data exp1\', subject_name], 1, 1);
    hemodynamic_peak = peak_d(1, :);
    behavioral = peak_d(2, :);
    
    %Correlation Coefficient%
    r = corrcoef(hemodynamic_peak,behavioral);
    rsq = r(1,2).^2;
    
    % plot linear model fit
    subplot(2,2,index);
    Mdl1=LinearModel.fit(behavioral,hemodynamic_peak);
    plot(Mdl1,'marker','s','MarkerSize',8)
    xlim([0,3.5]);
    ylim([0,1]);
    axis square;
    title('')
    xlabel('d prime');
    ylabel('model peak');
    
    text(min(behavioral)*0.8,max(hemodynamic_peak)*0.78,['R^2 = ',num2str(round(Mdl1.Rsquared.Adjusted*100)/100)])
    text(min(behavioral)*0.8,max(hemodynamic_peak)*0.7,['p = ',num2str(Mdl1.Coefficients.pValue(2))]);
    figure = gca; 
    legend(figure,'off');
   
end
