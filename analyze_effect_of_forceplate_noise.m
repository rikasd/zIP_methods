% 
% Compute zIP from Santos and Duarte's data + force plate noise
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-06-13
% Updated: 2026-11-06

clear;clc;close all;
%%
setup;

G_MPS = 9.81; % acceleration due to gravity
CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;
SAMPFREQ_FP1_HZ = 100;
SAMPFREQ_FP2_HZ = 1200;

height_m = 1.575; % height of preliminary subject
mass_kg = 54.2;   % mass of preliminary subject
com_height_m = 0.56*height_m; % estimate of CoM height

plane = 'sgt';% Frontal plane: 'frt', sagittal plane: 'sgt';

%%
%% Test force plate 1
%%
for iTrial = 1:3 % for each trial of data
    % Load human data
    human_data = ...
        readtable(['Santos2016Data\BDS0000' num2str(iTrial) '.txt']);

    % Load measured device noise from force plate 1...
    forceplate_data = ...
        load(['Forceplate_noise\forceplate1_noise' num2str(iTrial) '.mat']);    

    if strcmp(plane,'sgt')

        noise_Fx_N = forceplate_data.Fx_N;
        noise_Fz_N = forceplate_data.Fz_N;
        noise_xCOP_m = -forceplate_data.My_Nm./(mass_kg*G_MPS + forceplate_data.Fz_N);
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_noisy{iTrial}.FootForce = ...
            [human_data.Fx_N_ + noise_Fx_N,...
             human_data.Fz_N_ + noise_Fz_N]';
        data_noisy{iTrial}.COP = (human_data.COPx_cm_*CM2M + noise_xCOP_m)';
        data_noisy{iTrial}.COM = com_height_m;
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_FP1_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fx_N_,...
             human_data.Fz_N_]';
        data_original{iTrial}.COP = human_data.COPx_cm_'*CM2M;
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_FP1_HZ;

    elseif strcmp(plane,'frt')

        noise_Fy_N = forceplate_data.Fy_N;
        noise_Fz_N = forceplate_data.Fz_N;
        noise_yCOP_m = forceplate_data.Mx_Nm./(mass_kg*G_MPS + forceplate_data.Fz_N);
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_noisy{iTrial}.FootForce = ...
            [human_data.Fy_N_ + noise_Fy_N,...
             human_data.Fz_N_ + noise_Fz_N]';
        data_noisy{iTrial}.COP = (human_data.COPy_cm_*CM2M + noise_yCOP_m)';
        data_noisy{iTrial}.COM = com_height_m;
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_FP1_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data.Fy_N_,...
             human_data.Fz_N_]';
        data_original{iTrial}.COP = human_data.COPy_cm_'*CM2M;
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_FP1_HZ;

    end

end % iTrial = 1:3

%% Plot first trial of human data and noise as example representative plots

figure;
subplot(3,1,1);
plot(human_data.Time_s_',data_original{1}.FootForce(1,:));
ylabel('F_x (N)');ylim([-10,10]);
subplot(3,1,2);
plot(human_data.Time_s_',data_original{1}.FootForce(2,:));
ylabel('F_z (N)');ylim([530,550]);
subplot(3,1,3);
plot(human_data.Time_s_',data_original{1}.COP);
ylabel('x_{CoP} (m)');ylim([-0.091,-0.069]);
xlabel('Time (s)');
sgtitle('Human data');

forceplate_data = ...
        load(['Forceplate_noise\forceplate1_noise' num2str(1) '.mat']); 
figure;
subplot(3,1,1);
plot(forceplate_data.time_s',forceplate_data.Fx_N);
ylabel('F_x (N)');ylim([-10,10]);
subplot(3,1,2);
plot(forceplate_data.time_s',forceplate_data.Fz_N);
ylabel('F_z (N)');ylim([-10,10]);
subplot(3,1,3);
noise_xCOP_m = -forceplate_data.My_Nm./(mass_kg*G_MPS + forceplate_data.Fz_N);
plot(forceplate_data.time_s',noise_xCOP_m);
ylabel('x_{CoP} (m)');ylim([-0.011,0.011]);
xlabel('Time (s)');
sgtitle('Force plate 1 measurement noise');

figure;
subplot(3,1,1);
plot(human_data.Time_s_',data_noisy{1}.FootForce(1,:));
ylabel('F_x (N)');
subplot(3,1,2);
plot(human_data.Time_s_',data_noisy{1}.FootForce(2,:));
ylabel('F_z (N)');
subplot(3,1,3);
plot(human_data.Time_s_',data_noisy{1}.COP);
ylabel('x_{CoP} (m)');
xlabel('Time (s)');
sgtitle('Human data + force plate 1 noise')

%% Compute and plot zIP from original and noisy data

[f_Hz, zIP_original, GOF_original] = getZIPfromData(data_original);
[~, zIP_noisy, GOF_noisy] = getZIPfromData(data_noisy);

fig = figure('Units', 'centimeters', 'Position', [0 0 26 10.5]);
colors = get(gca,'colororder');
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
hold on;
plot(f_Hz,zIP_original,'k--');
plot(f_Hz,zIP_noisy,'Color',colors(1,:));
xlim([0,6]);xticks(0:6);
ylim([0,3]);
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','normalized by CoM height'});
yline(1,'k--');
legend('Original','Noise-added');
text(-0.2, 1.12, 'A', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left','FontSize', 20)
title('Force plate 1');

%%
%% Test force plate 2
%%
%% Repeat above analysis
%%
for iTrial = 1:3 % for each trial of data
    % Load human data
    human_data = ...
        readtable(['Santos2016Data\BDS0000' num2str(iTrial) '.txt']);

    % Load measured device noise from force plate 1...
    forceplate_data = ...
        load(['Forceplate_noise\forceplate2_noise' num2str(iTrial) '.mat']);    

    if strcmp(plane,'sgt')

        noise_Fx_N = forceplate_data.Fx_N;
        noise_Fz_N = forceplate_data.Fz_N;
        noise_xCOP_m = -forceplate_data.My_Nm./(mass_kg*G_MPS + forceplate_data.Fz_N);
    
        % Interpolate human data to 1200 Hz for force plate 2 noise.
        human_data_interp.Fx_N = interp1(...
            human_data.Time_s_,...
            human_data.Fx_N_, ...
            forceplate_data.time_s(12:end));
        human_data_interp.Fz_N = interp1(...
            human_data.Time_s_,...
            human_data.Fz_N_, ...
            forceplate_data.time_s(12:end));
        human_data_interp.xCOP_m = interp1(...
            human_data.Time_s_,...
            human_data.COPx_cm_*CM2M, ...
            forceplate_data.time_s(12:end));
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_noisy{iTrial}.FootForce = ...
            [human_data_interp.Fx_N + noise_Fx_N(12:end),...
             human_data_interp.Fz_N + noise_Fz_N(12:end)]';
        data_noisy{iTrial}.COP = (human_data_interp.xCOP_m + noise_xCOP_m(12:end))';
        data_noisy{iTrial}.COM = com_height_m;
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_FP2_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data_interp.Fx_N,...
             human_data_interp.Fz_N]';
        data_original{iTrial}.COP = human_data_interp.xCOP_m';
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_FP2_HZ;

    elseif strcmp(plane,'frt')
        noise_Fy_N = forceplate_data.Fy_N;
        noise_Fz_N = forceplate_data.Fz_N;
        noise_yCOP_m = forceplate_data.Mx_Nm./(mass_kg*G_MPS + forceplate_data.Fz_N);
    
        % Interpolate human data to 1200 Hz for force plate 2 noise.
        human_data_interp.Fy_N = interp1(...
            human_data.Time_s_,...
            human_data.Fy_N_, ...
            forceplate_data.time_s(12:end));
        human_data_interp.Fz_N = interp1(...
            human_data.Time_s_,...
            human_data.Fz_N_, ...
            forceplate_data.time_s(12:end));
        human_data_interp.yCOP_m = interp1(...
            human_data.Time_s_,...
            human_data.COPy_cm_*CM2M, ...
            forceplate_data.time_s(12:end));
    
        % Add shear force noise, vertical force noise, and CoP noise to human data
        data_noisy{iTrial}.FootForce = ...
            [human_data_interp.Fy_N + noise_Fy_N(12:end),...
             human_data_interp.Fz_N + noise_Fz_N(12:end)]';
        data_noisy{iTrial}.COP = (human_data_interp.yCOP_m + noise_yCOP_m(12:end))';
        data_noisy{iTrial}.COM = com_height_m;
        data_noisy{iTrial}.sampFreq_Hz = SAMPFREQ_FP2_HZ;
    
        % original human data
        data_original{iTrial}.FootForce = ...
            [human_data_interp.Fy_N,...
             human_data_interp.Fz_N]';
        data_original{iTrial}.COP = human_data_interp.yCOP_m';
        data_original{iTrial}.COM = com_height_m;
        data_original{iTrial}.sampFreq_Hz = SAMPFREQ_FP2_HZ;

    end

end % iTrial = 1:3
%% Compute and plot zIP from original and noisy data

params.window_size = 2^13;
[f2_Hz, zIP_noisy] = getZIPfromData(data_noisy,params);
[f_Hz, zIP_original] = getZIPfromData(data_original,params);

nexttile
hold on;
plot(f_Hz,zIP_original,'k--');
plot(f_Hz,zIP_noisy,'Color',colors(2,:));
xlim([0,6]);xticks(0:6);
ylim([0,3]);
xlabel('Frequency (Hz)');
yline(1,'k--');
legend('Original','Noise-added');
text(-0.14, 1.12, 'B', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left','FontSize', 20)
title('Force plate 2');