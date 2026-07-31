% 
% Compute zIP from Santos and Duarte's data + artificial measurement noise
% Take data from participants with lowest Fx PSD:
% OLD subject 22 (in sagittal and frontal planes), 26 (in frontal plane)
% YOUNG subject 1 (in sagittal plane), 9 (in frontal plane)
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-06-13

clear;clc;close all;
ifig = 5;
%%
setup;

CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;

bool_plot_zIP = true;

first_data_file_nums = [1,115,1048,1129]; % Data file number of first trials of participants of interest (with lowest Fx PSD)
height_m = [1.575, 1.627, 1.40, 1.485]; % heights of particicpants
com_height_m = 0.56*height_m; % estimate of CoM height

plane = 'sgt';% Frontal plane: 'frt', sagittal plane: 'sgt';

%%
switch plane
    case 'sgt'
        % for sgt, keeps zIP/zCoM error < 0.1 most of the time (depends on noise)
        noise_Fx_N = 0.008*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.007-0.04
        noise_Fz_N = 0.1*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.014-0.3
        noise_My_Nm = 0.005*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.004-0.01
    case 'frt'
        % for frt, keeps zIP/zCoM error < 0.1 most of the time (depends on noise)
        noise_Fx_N = 0.005*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.007-0.04
        noise_Fz_N = 0.3*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.014-0.3
        noise_My_Nm = 0.003*sqrt(SAMPFREQ_HUMAN_DATA_HZ); % 0.004-0.01
end

%%
for iS = 1:4
for iTrial = 1:3 % for each trial of data
    % Load human data
    human_data = ...
        readtable(['Santos2016Data\BDS0' num2str(first_data_file_nums(iS) + iTrial - 1,'%04i') '.txt']); 

    if strcmp(plane,'sgt')
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_size = size(human_data.Fx_N_);
        data_noisy{iTrial}.FootForce = ...
            [human_data.Fx_N_ + noise_Fx_N*randn(data_size),...
             human_data.Fz_N_ + noise_Fz_N*randn(data_size)]';
        data_noisy{iTrial}.FootMoment = (-human_data.My_Nm_ + noise_My_Nm*randn(data_size))';
        data_noisy{iTrial}.COP = data_noisy{iTrial}.FootMoment./data_noisy{iTrial}.FootForce(2,:);
        data_noisy{iTrial}.COM = com_height_m(iS);
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fx_N_,...
             human_data.Fz_N_]';
        data_original{iTrial}.FootMoment = -human_data.My_Nm_';
        data_original{iTrial}.COP = human_data.COPx_cm_'*CM2M;
        data_original{iTrial}.COM = com_height_m(iS);
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;

    elseif strcmp(plane,'frt')
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_size = size(human_data.Fy_N_);
        data_noisy{iTrial}.FootForce = ...
            [human_data.Fy_N_ + noise_Fx_N*randn(data_size),...
             human_data.Fz_N_ + noise_Fz_N*randn(data_size)]';
        data_noisy{iTrial}.FootMoment = (human_data.Mx_Nm_ + noise_My_Nm*randn(data_size))';
        data_noisy{iTrial}.COP = data_noisy{iTrial}.FootMoment./data_noisy{iTrial}.FootForce(2,:);
        data_noisy{iTrial}.COM = com_height_m(iS);
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fy_N_,...
             human_data.Fz_N_]';
        data_original{iTrial}.FootMoment = human_data.Mx_Nm_';
        data_original{iTrial}.COP = human_data.COPy_cm_'*CM2M;
        data_original{iTrial}.COM = com_height_m(iS);
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;

    end

end % iTrial = 1:3

%% Compute and plot zIP from original and noisy data

params.method = 'cpsd_gain';
[f_Hz, zIP_original, GOF_original] = getZIPfromData(data_original,params);
[~, zIP_noisy, GOF_noisy] = getZIPfromData(data_noisy,params);

if bool_plot_zIP
    figure;hold on;
    plot(f_Hz,zIP_original);
    plot(f_Hz,zIP_noisy);
    xlim([0,6]);xticks(0:6);
    ylim([0,3]);
    xlabel('Frequency (Hz)');
    ylabel({'Intersection point height','normalized by CoM height'});
    yline(1,'k--');
    legend('Original','Noise-added');
end

%% Plot difference between clean and noisy zIP
figure(ifig);hold on;
plot(f_Hz,abs(zIP_noisy - zIP_original));
xlim([0,6]);xticks(0:6);
xlabel('Frequency (Hz)');
ylabel('Absolute error');

end % iS

yline(0.1,'k--');