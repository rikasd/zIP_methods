% 
% Example of zIP computed from Santos and Duarte's data with default
% parameters
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-06-15

clear;clc;close all;
%%
setup;

CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;

height_m = 1.575; % height of preliminary subject
com_height_m = 0.56*height_m; % estimate of CoM height

plane = 'sgt';% Frontal plane: 'frt', sagittal plane: 'sgt';

trial_duration_s = 60;

%%
dt_s = 1/SAMPFREQ_HUMAN_DATA_HZ;
time_s = dt_s:dt_s:trial_duration_s;
indices = 1:length(time_s);
for iTrial = 1:3 % for each trial of data
    % Load human data
    human_data = ...
        readtable(['Santos2016Data\BDS0000' num2str(iTrial) '.txt']);

    if strcmp(plane,'sgt')
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fx_N_(indices),...
             human_data.Fz_N_(indices)]';
        data_original{iTrial}.COP = human_data.COPx_cm_(indices)'*CM2M;
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;

    elseif strcmp(plane,'frt')
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fy_N_(indices),...
             human_data.Fz_N_(indices)]';
        data_original{iTrial}.COP = human_data.COPy_cm_(indices)'*CM2M;
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_HUMAN_DATA_HZ;

    end

end % iTrial = 1:3

%% Plot first trial of human data as example representative plots

figure;
subplot(3,1,1);
plot(human_data.Time_s_(indices)',data_original{1}.FootForce(1,:));
ylabel('F_x (N)');ylim([-5,1]);box off;
subplot(3,1,2);
plot(human_data.Time_s_(indices)',data_original{1}.FootForce(2,:));
ylabel('F_z (N)');ylim([535,541]);yticks(536:2:540);box off;
subplot(3,1,3);
plot(human_data.Time_s_(indices)',data_original{1}.COP);
ylabel('x_{CoP} (m)');ylim([-0.091,-0.069]);box off;
xlabel('Time (s)');

%% Compute and plot zIP with default filter parameters

[f_Hz, zIP] = getZIPfromData(data_original);

figure;hold on;
plot(f_Hz,zIP);
xlim([0,6]);xticks(0:6);
ylim([0,3]);yline(1,'--k');
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','(normalized by CoM height)'});