% 
% Compute zIP from Santos and Duarte's data 
% with variations in the experimental design (number of trials + duration)
% to illustrate how these parameters affect the zIP
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


%% Effect of number of trials

[f_Hz{1}, zIP{1}] = getZIPfromData(data_original);
[f_Hz{2}, zIP{2}] = getZIPfromData(data_original{1});
[f_Hz{3}, zIP{3}] = getZIPfromData({data_original{1},data_original{2}});

figure;hold on;colors = get(gca,'colororder');
plot(f_Hz{2},zIP{2},'-','Color',colors(3,:));
plot(f_Hz{3},zIP{3},'-.','Color',colors(2,:));
plot(f_Hz{1},zIP{1},'--','Color',colors(1,:));
xlim([0,6]);xticks(0:6);
ylim([0,3]);
xlabel('Frequency (Hz)');
ylabel({'Intersection point height','normalized by CoM height'});
yline(1,'k--');
legend('1 trial','2 trials','3 trials');

for i = 2:3
index_freq = find((f_Hz{i} > 0).*(f_Hz{i} < 6));
rms_numTrials(i) = rms(zIP{i}(index_freq)-zIP{1}(index_freq));
end

%% Effect of trial duration

figure;hold on;
colors = get(gca,'colororder');
line_type = {'-','-.','--'};
colors_reordered(1,:) = colors(3,:);
colors_reordered(2,:) = colors(2,:);
colors_reordered(3,:) = colors(1,:);
xlim([0,6]);xticks(0:6);
ylim([0,3]);

i = 1;

for trial_duration_s = [10 30 60]

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


[f_Hz{i}, zIP{i}] = getZIPfromData(data_original);


plot(f_Hz{i},zIP{i},line_type{i},'Color',colors_reordered(i,:));

i = i + 1;

end

xlabel('Frequency (Hz)');
ylabel({'Intersection point height','normalized by CoM height'});
yline(1,'k--');
legend('10 s x 3','30 s x 3','60 s x 3');


for i = 2:3
index_freq = find((f_Hz{i} > 0).*(f_Hz{i} < 6));
rms_trialDur(i) = rms(zIP{i}(index_freq)-zIP{1}(index_freq));
end