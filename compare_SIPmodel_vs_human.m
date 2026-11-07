% 
% Compute zIP from Santos and Duarte's data 
% and compare with single inverted pendulum (SIP) model
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-06-20

clear;clc;close all;
%%
setup;

G_MPS = 9.81; % acceleration due to gravity
CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;

height_m = 1.575; % height of preliminary subject
mass_kg = 54.2;   % mass of preliminary subject
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


%% Single-inverted-pendulum model parameters
Ic = 5.7;         % moment of inertia about center of mass
m = mass_kg;      % mass
l = com_height_m; % center of mass distance from ankle
hf = 0.07;        % distance of ankle joint from floor

%% Human vs. model
[f_Hz, zIP_original] = getZIPfromData(data_original);

zIP_model = (Ic/m/l + l + hf + G_MPS/4/pi^2./f_Hz.^2)/l;

figure;hold on;colors = get(gca,'colororder');
plot(f_Hz,zIP_original,'--','Color',colors(1,:));
plot(f_Hz,zIP_model,'Color',colors(2,:));
xlim([0,6]);xticks(0:6);
ylim([0,3]);
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','normalized by CoM height'});
yline(1,'k--');
legend('Human','SIP model');
