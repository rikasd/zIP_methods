% 
% Compute zIP from Santos and Duarte's data 
% and compare with double inverted pendulum (DIP) model
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-07-30

clear;clc;close all;
%%
setup;

G_MPS = 9.81; % acceleration due to gravity
CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;

height_m = 1.575; % height of preliminary subject
mass_kg = 54.2;   % mass of preliminary subject
com_height_m = 0.56*height_m; % estimate of CoM height
gender = 'M';     % Female: 'F', Male: 'M'
pose = 'pose_I';
agegroup = 1;     % Young: 1, Old: 2

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


%% Double-inverted-pendulum model parameters
lumped_params = getLumpedParams_DIP(mass_kg,height_m,gender,plane,pose,agegroup);
lumped_params.h_foot = 0.07;        % distance of ankle joint from floor
input_struct.lumped_params = lumped_params;
input_struct.controller_params.K_x =... % feedback gains
    [ 530  350  180   75;
        5  230   16   33];
input_struct.motorNoiseRatio = 0.01;     % motor noise ratio: ankle noise level over hip noise level

%% Human vs. model
[f_data_Hz, zIP_data] = getZIPfromData(data_original);
[f_model_Hz, zIP_model] = predictZIPfromModel(input_struct);

figure;hold on;colors = get(gca,'colororder');
plot(f_data_Hz,zIP_data,'--','Color',colors(1,:));
plot(f_model_Hz,zIP_model,'Color',colors(2,:));
xlim([0,6]);xticks(0:6);
ylim([0,3]);
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','normalized by CoM height'});
yline(1,'k--');
legend('Human','DIP model');
