%
% Effect of force plate noise on zIP, across all participants
% (Santos and Duarte 2016 data; vision open, firm surface; 3 trials per
% participant combined into one zIP estimate).
%
% For each participant, measured noise from force plate 1 (100 Hz) and
% force plate 2 (1200 Hz) is added to the human data (trial k o-f the
% participant gets noise recording k). zIP is computed from the original and
% noise-added data, and the absolute error |zIP_noisy - zIP_original| is
% computed. The figure shows the mean +/- SD of this error across
% participants as a function of frequency, for each force plate.
%
% Code generated using Claude Sonnet 5.5 on 2026-10-21
% based on reference code provided by Rika Sugimoto Dimitrova

clear; clc; close all;

%% ------------------------------------------------------------------------
%  Settings (same values as analyze_effect_of_forceplate_noise.m)
%  ------------------------------------------------------------------------
setup;

bool_compute_errors = false; % If true, compute zIP from data and compute 
                             % errors (takes a few mintues)
                             % If false, load saved results

if bool_compute_errors

cfg.G_MPS  = 9.81;       % acceleration due to gravity
cfg.CM2M   = 0.01;       % centimeters to meters
HEIGHT2COM = 0.56;       % CoM height / body height

SAMPFREQ_FP1_HZ = 100;
SAMPFREQ_FP2_HZ = 1200;
FP2_START_IDX   = 12;    % first sample of the 1200 Hz noise that is used

plane = 'sgt';           % Frontal plane: 'frt', sagittal plane: 'sgt'

datapath  = 'Santos2016Data';
noisepath = 'Forceplate_noise';
FIG_DIR   = 'figures';

N_TRIALS_EXPECTED = 3;

MASS_COLUMN_CANDIDATES = {'Weight','Mass','BodyMass','Weight_kg_'};

% zIP parameters for force plate 2 (force plate 1 uses the defaults, as before)
params2.window_size = 2^13;

% Plot settings
F_XLIM_HZ = [0 6];
Y_MAX     = [];          % y-limit (upper) for the error plots; [] = automatic
LABEL_X   = -0.20;       % panel-letter x position (normalized axes units)
LABEL_Y   = 1.08;        % panel-letter y position

warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%% ------------------------------------------------------------------------
%  Select participants / trials: vision open, firm surface
%  ------------------------------------------------------------------------
info = readtable(fullfile(datapath, 'BDSinfo.txt'));

isSel    = strcmp(info.Vision, 'Open') & strcmp(info.Surface, 'Firm');
subjects = unique(info.Subject(isSel))';
nSubj    = numel(subjects);

fprintf('Found %d participants (eyes open, firm surface).\n', nSubj);

% Find the mass column (if any)
massCol = '';
for c = MASS_COLUMN_CANDIDATES
    if ismember(c{1}, info.Properties.VariableNames)
        massCol = c{1}; break;
    end
end
if isempty(massCol)
    warning(['No mass column found in BDSinfo.txt; using %.1f kg for all ' ...
             'participants. Edit MASS_COLUMN_CANDIDATES if needed.'], MASS_DEFAULT_KG);
end

%% ------------------------------------------------------------------------
%  Compute original and noisy zIP for every participant and force plate
%  ------------------------------------------------------------------------
err1 = [];  err2 = [];            % [participant x frequency] absolute errors
f1_Hz = []; f2_Hz = [];
subjUsed = false(1, nSubj);

for iS = 1:nSubj
    rows = find(isSel & info.Subject == subjects(iS));

    % Sort trials by name so trial k pairs with noise recording k
    [~, ord] = sort(info.Trial(rows));
    rows = rows(ord);

    if numel(rows) ~= N_TRIALS_EXPECTED
        warning('Subject %d has %d trials (expected %d). Skipped.', ...
                subjects(iS), numel(rows), N_TRIALS_EXPECTED);
        continue;
    end

    % Anthropometrics
    height_m = info.Height(rows(1)) * cfg.CM2M;
    if ~isempty(massCol), mass_kg = info.(massCol)(rows(1));
    else,                 mass_kg = MASS_DEFAULT_KG; end
    com_height_m = HEIGHT2COM * height_m;

    % Containers for the 3 trials
    d1_noisy = cell(1, N_TRIALS_EXPECTED); d1_orig = d1_noisy;
    d2_noisy = d1_noisy;                   d2_orig = d1_noisy;

    for iT = 1:N_TRIALS_EXPECTED
        human_data = readtable(fullfile(datapath, [info.Trial{rows(iT)} '.txt']));

        % ---- Force plate 1 (100 Hz, same sampling as human data)
        fp1 = load(fullfile(noisepath, sprintf('forceplate1_noise%d.mat', iT)));
        [d1_noisy{iT}, d1_orig{iT}] = buildTrialData(human_data, fp1, plane, ...
            mass_kg, com_height_m, SAMPFREQ_FP1_HZ, 1, cfg);

        % ---- Force plate 2 (1200 Hz; human data interpolated)
        fp2 = load(fullfile(noisepath, sprintf('forceplate2_noise%d.mat', iT)));
        [d2_noisy{iT}, d2_orig{iT}] = buildTrialData(human_data, fp2, plane, ...
            mass_kg, com_height_m, SAMPFREQ_FP2_HZ, FP2_START_IDX, cfg);
    end

    try
        % Force plate 1: default zIP parameters
        [f1, zIP1_orig] = getZIPfromData(d1_orig);
        [~,  zIP1_noisy] = getZIPfromData(d1_noisy);

        % Force plate 2: window size 2^13
        [f2o, zIP2_orig]  = getZIPfromData(d2_orig,  params2);
        [f2n, zIP2_noisy] = getZIPfromData(d2_noisy, params2);

        assert(numel(f2o) == numel(f2n), 'FP2 frequency vectors differ.');

        e1 = abs(zIP1_noisy(:)' - zIP1_orig(:)');
        e2 = abs(zIP2_noisy(:)' - zIP2_orig(:)');

        % Initialise on first success; then check consistent frequency axes
        if isempty(err1)
            f1_Hz = f1(:)';  f2_Hz = f2o(:)';
            err1 = nan(nSubj, numel(f1_Hz));
            err2 = nan(nSubj, numel(f2_Hz));
        end
        assert(numel(e1) == numel(f1_Hz) && numel(e2) == numel(f2_Hz), ...
               'Frequency axis differs from the first participant.');

        err1(iS,:) = e1;
        err2(iS,:) = e2;
        subjUsed(iS) = true;
    catch ME
        warning('Subject %d failed: %s', subjects(iS), ME.message);
    end

    fprintf('  Subject %2d done (%d/%d)\n', subjects(iS), iS, nSubj);
end

fprintf('%d of %d participants processed successfully.\n', nnz(subjUsed), nSubj);

%% ------------------------------------------------------------------------
%  Mean and SD across participants
%  ------------------------------------------------------------------------
mean1 = mean(err1, 1, 'omitnan');   sd1 = std(err1, 0, 1, 'omitnan');
mean2 = mean(err2, 1, 'omitnan');   sd2 = std(err2, 0, 1, 'omitnan');

%% ------------------------------------------------------------------------
%  Summary + save
%  ------------------------------------------------------------------------
inBand1 = f1_Hz >= F_XLIM_HZ(1) & f1_Hz <= F_XLIM_HZ(2);
inBand2 = f2_Hz >= F_XLIM_HZ(1) & f2_Hz <= F_XLIM_HZ(2);
fprintf('\nMean absolute error over %g-%g Hz (mean across participants):\n', F_XLIM_HZ);
fprintf('  Force plate 1: %.4f\n', mean(mean1(inBand1), 'omitnan'));
fprintf('  Force plate 2: %.4f\n', mean(mean2(inBand2), 'omitnan'));

save('./Results/zIP_forceplate_noise_allSubj.mat', 'err1','err2','f1_Hz','f2_Hz', ...
     'mean1','sd1','mean2','sd2','subjects','subjUsed','plane');

else 
load('./Results/zIP_forceplate_noise_allSubj.mat');
end % bool_compute_errors

%% ------------------------------------------------------------------------
%  Figure: mean +/- SD of absolute error vs. frequency, both force plates
%  ------------------------------------------------------------------------

figure;hold on;colors = get(gca,'colororder');
plotMeanSD(f1_Hz, mean1, sd1, colors(1,:));
plotMeanSD(f2_Hz, mean2, sd2, colors(2,:));
xlabel('Frequency (Hz)');
ylabel('Absolute error');
xlim([0,6]);
legend('Force plate 1','Force plate 2','Location','Northwest');


%% ========================================================================
%  Local functions
%  ========================================================================
function [d_noisy, d_orig] = buildTrialData(human, fp, plane, mass_kg, zCoM, ...
                                            fs, startIdx, cfg)
% Build "noisy" and "original" data structs for getZIPfromData, following
% analyze_effect_of_forceplate_noise.m.
%   startIdx = 1  : noise and human data share the same time base (100 Hz)
%   startIdx = 12 : human data are interpolated onto the noise time base
%                   (1200 Hz), using noise samples startIdx:end.

    switch plane
        case 'sgt'
            hF   = human.Fx_N_(:);
            hCOP = human.COPx_cm_(:) * cfg.CM2M;
            nF   = fp.Fx_N(:);
            nCOP = -fp.My_Nm(:) ./ (mass_kg*cfg.G_MPS + fp.Fz_N(:));
        case 'frt'
            hF   = human.Fy_N_(:);
            hCOP = human.COPy_cm_(:) * cfg.CM2M;
            nF   = fp.Fy_N(:);
            nCOP =  fp.Mx_Nm(:) ./ (mass_kg*cfg.G_MPS + fp.Fz_N(:));
        otherwise
            error('plane must be ''sgt'' or ''frt''');
    end
    hFz = human.Fz_N_(:);
    nFz = fp.Fz_N(:);

    hF = hF - mean(hF);

    if startIdx == 1
        % Same sampling rate: add noise directly (truncate to common length)
        n = min(numel(hF), numel(nF));
        hF = hF(1:n); hFz = hFz(1:n); hCOP = hCOP(1:n);
        nF = nF(1:n); nFz = nFz(1:n); nCOP = nCOP(1:n);
    else
        % Interpolate human data to the force plate noise time base
        t  = fp.time_s(:);  t = t(startIdx:end);
        hF   = interp1(human.Time_s_, hF,   t);
        hFz  = interp1(human.Time_s_, hFz,  t);
        hCOP = interp1(human.Time_s_, hCOP, t);
        nF   = nF(startIdx:end);
        nFz  = nFz(startIdx:end);
        nCOP = nCOP(startIdx:end);
    end

    % Noise-added data
    d_noisy.FootForce   = [hF + nF, hFz + nFz]';
    d_noisy.COP         = (hCOP + nCOP)';
    d_noisy.COM         = zCoM;
    d_noisy.sampFreq_Hz = fs;

    % Original human data (same samples)
    d_orig.FootForce    = [hF, hFz]';
    d_orig.COP          = hCOP';
    d_orig.COM          = zCoM;
    d_orig.sampFreq_Hz  = fs;
end

function plotMeanSD(f, m, s, col)
% Mean line with a shaded +/- 1 SD band (lower band clipped at 0, since the
% quantity is an absolute error).
    lo = max(m - s, 0);
    hi = m + s;
    ok = isfinite(m) & isfinite(lo) & isfinite(hi);
    fill([f(ok) fliplr(f(ok))], [hi(ok) fliplr(lo(ok))], col, ...
         'FaceAlpha',0.25, 'EdgeColor','none', 'HandleVisibility','off');
    plot(f, m, 'Color', col, 'LineWidth', 1.5);
    grid on; box on;
end