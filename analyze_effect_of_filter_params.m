% 
% Compute zIP from Santos and Duarte's data 
% with variations in the zIP computation parameters to illustrate how these
% parameters affect the zIP
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

%% Effect of window size
params.window_taper_option = 2;
params.detrend_option = 3;

params.window_size = 2^8;
[f_Hz{1}, zIP{1}] = getZIPfromData(data_original,params);
params.window_size = 2^9;
[f_Hz{2}, zIP{2}] = getZIPfromData(data_original,params);
params.window_size = 2^10;
[f_Hz{3}, zIP{3}] = getZIPfromData(data_original,params);
params.window_size = 2^11;
[f_Hz{4}, zIP{4}] = getZIPfromData(data_original,params);

figure;hold on;colors = get(gca,'colororder');
plot(f_Hz{4},zIP{4}+1.5,'.','Color',colors(5,:),'MarkerSize',10);
plot(f_Hz{3},zIP{3}+1,'.','Color',colors(3,:),'MarkerSize',10);
plot(f_Hz{1},zIP{1}+0.5,'.','Color',colors(1,:),'MarkerSize',10);
plot(f_Hz{2},zIP{2},'.','Color',colors(2,:),'MarkerSize',10);
xlim([0,6]);xticks(0:6);
ylim([0,4]);yticks('');
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','(curves offset for clarity)'});
legend('2^{11} (20.5 s)','2^{10} (10.2 s)','2^9 (5.12 s)','2^8 (2.56 s)');

for i = 1:4
index_freq = find((f_Hz{i} > 1).*(f_Hz{i} < 6));
dzIP{i} = ((zIP{i}(index_freq)-zIP{i}(index_freq-1))/(f_Hz{i}(2)-f_Hz{i}(1)));
smoothness_win(i) = var(dzIP{i});
end

%% Effect of window taper
params.window_size = 2^9;
params.detrend_option = 3;

params.window_taper_option = 2;
[f_Hz{2}, zIP{2}] = getZIPfromData(data_original,params);
params.window_taper_option = 1;
[f_Hz{1}, zIP{1}] = getZIPfromData(data_original,params);
params.window_taper_option = 3;
[f_Hz{3}, zIP{3}] = getZIPfromData(data_original,params);

figure;hold on;colors = get(gca,'colororder');
plot(f_Hz{1},zIP{1}+1.5,'-','Color',colors(5,:));
plot(f_Hz{3},zIP{3}+1,'-.','Color',colors(3,:));
plot(f_Hz{2},zIP{2}+0.5,'--','Color',colors(1,:));
xlim([0,6]);xticks(0:6);
ylim([0,4]);yticks('');
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','(curves offset for clarity)'});
legend('Rectangular','Hamming','Hann');

for i = 1:3
index_freq = find((f_Hz{i} > 1).*(f_Hz{i} < 6));
dzIP{i} = ((zIP{i}(index_freq)-zIP{i}(index_freq-1))/(f_Hz{i}(2)-f_Hz{i}(1)));
smoothness_taper(i) = var(dzIP{i});
end

%% Effect of detrend
params.window_size = 2^9;
params.window_taper_option = 2;

params.detrend_option = 3;
[f_Hz{3}, zIP{3}] = getZIPfromData(data_original,params);
params.detrend_option = 1;
[f_Hz{1}, zIP{1}] = getZIPfromData(data_original,params);
params.detrend_option = 2;
[f_Hz{2}, zIP{2}] = getZIPfromData(data_original,params);

figure;hold on;colors = get(gca,'colororder');
plot(f_Hz{1},zIP{1}+1.5,'-','Color',colors(5,:));
plot(f_Hz{2},zIP{2}+1,'-.','Color',colors(3,:));
plot(f_Hz{3},zIP{3}+0.5,'--','Color',colors(1,:));
xlim([0,6]);xticks(0:6);
ylim([0,4]);yticks('');
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','(curves offset for clarity)'});
legend('no detrend','zero-mean','linear detrend');

for i = 1:2
index_freq = find((f_Hz{i} > 0).*(f_Hz{i} < 6));
rms_detrend(i) = rms(zIP{i}(index_freq)-zIP{3}(index_freq));
end