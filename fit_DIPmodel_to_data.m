% 
% Compute zIP from Santos and Duarte's data 
% and compare with double inverted pendulum (DIP) model
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2026-10-20

clear;clc;close all;
%%
setup;

G_MPS = 9.81; % acceleration due to gravity
CM2M = 0.01;  % centimeters to meters

SAMPFREQ_HUMAN_DATA_HZ = 100;

height_m = 1.575; % height of preliminary subject [meters]
mass_kg = 54.2;   % mass of preliminary subject [kilograms]
ankle_height_m = 0.07; % estimated height of ankle joint from floor [meters]
com_height_m = 0.56*height_m; % estimated CoM height [meters]
sex = 'M';    % Female: 'F', Male: 'M'
pose = 'pose_I';
agegroup = 1; % Young: 1, Old: 2

plane = 'sgt';% Frontal plane: 'frt', sagittal plane: 'sgt';

trial_duration_s = 60;

max_freq_Hz = 6;

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

[f_data_Hz_temp, zIP_data_temp] = getZIPfromData(data_original);
index_max_freq = find(f_data_Hz_temp > max_freq_Hz, 1);
f_data_Hz = f_data_Hz_temp(1:index_max_freq);
zIP_data = zIP_data_temp(1:index_max_freq);

%% Double-inverted-pendulum model parameters
lumped_params = getLumpedParams_DIP(mass_kg,height_m,sex,plane,pose,agegroup);
lumped_params.h_foot = ankle_height_m;
input_struct.lumped_params = lumped_params;
input_struct.f = f_data_Hz;

%% Fit model
%  Search through LQR cost function parameters to find a model that best
%  reproduces the human zIP.
%  LQR cost function is given by: sum of x'*Q*x + u'*R*u
%  where x = states, u = torque commands, 
%  Q = state cost = eye(4) , R = input cost = alpha*[beta omega; omega 1/beta]
disp(['Searching through LQR cost function parameters ' ...
    'to find a model that best reproduces human zIP.']);
disp('This may take a few minutes.')
input_struct.controller_params.Q = eye(4);
beta_vec = 0.15:0.01:0.45; 
omega_vec = -0.9:0.1:0.9;     
sigma_r_vec = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];% motor noise ratio: ankle noise level over hip noise level

for b = 1:length(beta_vec)
for o = 1:length(omega_vec)
    input_struct.controller_params.R = ...
        1e6*[beta_vec(b) omega_vec(o); ...
                     omega_vec(o) 1/beta_vec(b)];
    for s = 1:length(sigma_r_vec)
        input_struct.motorNoiseRatio = sigma_r_vec(s);
        [~, zIP_model] = predictZIPfromModel(input_struct);
        rmse_results(b,o,s) = rmse(zIP_data,zIP_model);
    end
end
end

%% Plot best fit model
[v_min,loc_min] = min(abs(rmse_results(:)));
[b_min,o_min,s_min] = ind2sub(size(rmse_results),loc_min);
bestfit_params.beta = beta_vec(b_min);
bestfit_params.omega = omega_vec(o_min);
bestfit_params.sigma_r = sigma_r_vec(s_min);

input_struct.controller_params.R = ...
    1e6*[bestfit_params.beta bestfit_params.omega; ...
                         bestfit_params.omega 1/bestfit_params.beta];
input_struct.motorNoiseRatio = bestfit_params.sigma_r;

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

%% Analyze feedback gains

% Obtained feedback gains based on LQR cost function
[~,~,~,~,K_x] = getModel_linDIP(lumped_params,input_struct.controller_params);

% Extract stiffness matrix
K = K_x(1:2,1:2);

% Compute asymmetric and symmetric components of stiffness
K_A = (K-K')/2;
K_S = (K+K')/2;

% Compute level of asymmetry
K_ratio = sqrt(det(K_A))/sqrt(det(K_S));