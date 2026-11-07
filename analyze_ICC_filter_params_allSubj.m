% analyze_ICC_filter_params_allSubj.m
%
% Sensitivity analysis of the zIP whole-curve ICC across preprocessing 
% parameters to address reviewer comments.
% Evaluates test-retest reliability using full 60-second trials.
%
% Parameters evaluated:
%   - Window size: 2^8, 2^9, 2^10, 2^11
%   - Taper: 1 (Rectangular), 2 (Hann), 3 (Hamming)
%   - Detrend: 1 (None), 2 (Zero mean), 3 (Linear)
%
% Based on: analyze_ICC_trial_length_allSubj.m
%
% Code generated using Gemini 3.1 Pro on 2026-10-26

clear; clc;close all;

%% ------------------------------------------------------------------------
%  User settings
%  ------------------------------------------------------------------------
setup;

bool_compute_ICC = false; % If true, compute zIP from data and compute ICC
                          % (takes a few mintues) 
                          % If false, load saved results

datapath = 'Santos2016Data';
plane    = 'sgt';                 % 'sgt' (sagittal) or 'frt' (frontal)

CM2M       = 0.01;
HEIGHT2COM = 0.56;
FS_HZ      = 100;
DEMEAN_FH  = true;

TRIAL_LEN_S = 60;                 % Fixed at 60s
N_TRIALS    = 3;                  % trials per participant (k)

F_BAND_HZ   = [0 6];              % bins used: f > lo & f < hi 

ALPHA       = 0.05;               % 95% intervals
N_BOOT      = 5000;
RNG_SEED    = 1;

CURVE_WEIGHTING = 'none';         % 'none' or 'invVar'

% Parameter grids
WIN_SIZES = [2^8, 2^9, 2^10, 2^11];
TAPERS    = [1, 2, 3];
DETRENDS  = [1, 2, 3];

nW = numel(WIN_SIZES);
nT = numel(TAPERS);
nD = numel(DETRENDS);

% Baseline/Recommended parameters for plotting
idx_base_w = 2; % 2 = 2^9
idx_base_t = 2; % 2 = Hann
idx_base_d = 3; % 3 = Linear

FIG_DIR = 'figures';

warning('off','MATLAB:table:ModifiedAndSavedVarnames');

if bool_compute_ICC
%% ------------------------------------------------------------------------
%  Select participants: vision open, firm surface
%  ------------------------------------------------------------------------
info = readtable(fullfile(datapath, 'BDSinfo.txt'));
isSel    = strcmp(info.Vision, 'Open') & strcmp(info.Surface, 'Firm');
subjects = unique(info.Subject(isSel))';
nSubj    = numel(subjects);
k        = N_TRIALS;
fprintf('Found %d participants (vision open, firm surface).\n', nSubj);

%% ------------------------------------------------------------------------
%  Compute single-trial zIP for all parameter combinations
%  ------------------------------------------------------------------------
fAll = cell(nSubj, k, nW, nT, nD);
zAll = cell(nSubj, k, nW, nT, nD);

nReq = round(TRIAL_LEN_S * FS_HZ);

for iS = 1:nSubj
    rows = find(isSel & info.Subject == subjects(iS));
    height_m = info.Height(rows(1)) * CM2M;

    if numel(rows) < k
        warning('Subject %d has only %d trials; excluded.', subjects(iS), numel(rows));
        continue;
    end
    rows = rows(1:k);

    tables = cell(1, k);
    for iT = 1:k
        tables{iT} = readtable(fullfile(datapath, [info.Trial{rows(iT)} '.txt']));
    end

    for iT = 1:k
        nSamp = min(nReq, height(tables{iT}));
        if nReq - nSamp > 0.5*FS_HZ
            warning('Subject %d trial %d shorter than %g s; skipped.', ...
                    subjects(iS), iT, TRIAL_LEN_S);
            continue;
        end
        d = formatTrial(tables{iT}, plane, height_m, HEIGHT2COM, ...
                        CM2M, FS_HZ, nSamp, DEMEAN_FH);
                    
        for w = 1:nW
            for t = 1:nT
                for dt = 1:nD
                    params.window_size         = WIN_SIZES(w);
                    params.window_taper_option = TAPERS(t);
                    params.detrend_option      = DETRENDS(dt);
                    try
                        [f, z] = getZIPfromData({d}, params);
                        fAll{iS,iT,w,t,dt} = f(:);
                        zAll{iS,iT,w,t,dt} = z(:);
                    catch ME
                        warning('Subject %d, trial %d, Win %d, Tap %d, Det %d failed: %s', ...
                                subjects(iS), iT, WIN_SIZES(w), TAPERS(t), DETRENDS(dt), ME.message);
                    end
                end
            end
        end
    end
    fprintf('  Subject %3d done (%d/%d)\n', subjects(iS), iS, nSubj);
end

%% ------------------------------------------------------------------------
%  Assemble array Z: participants x trials x bins x Win x Tap x Det
%  ------------------------------------------------------------------------
% Extract fRef correctly from the baseline parameters 
slice_fAll = fAll(:,:,idx_base_w, idx_base_t, idx_base_d);
[valid_row, valid_col] = find(~cellfun(@isempty, slice_fAll), 1);
fRef = fAll{valid_row, valid_col, idx_base_w, idx_base_t, idx_base_d};

inBand  = fRef > F_BAND_HZ(1) & fRef < F_BAND_HZ(2);
fBand   = fRef(inBand);
nF      = numel(fBand);

Z = nan(nSubj, k, nF, nW, nT, nD);
for iS = 1:nSubj
    for iT = 1:k
        for w = 1:nW
            for t = 1:nT
                for dt = 1:nD
                    z = zAll{iS,iT,w,t,dt};
                    if isempty(z), continue; end
                    f = fAll{iS,iT,w,t,dt};
                    
                    % Interpolate to align frequency bins
                    if numel(f) ~= numel(fRef) || any(abs(f - fRef) > 1e-9)
                        % FIXED: Use 'extrap' to prevent NaNs at the lowest frequency bins
                        % where 2^8 has lower resolution than the 2^9 reference.
                        z = interp1(f, z, fRef, 'linear', 'extrap');
                    end
                    Z(iS,iT,:,w,t,dt) = z(inBand);
                end
            end
        end
    end
end

% FIXED: Base participant completeness strictly on the baseline parameters 
Z_base = Z(:,:,:,idx_base_w, idx_base_t, idx_base_d);
isComplete = all(isfinite(Z_base), [2 3]); 

% Apply filter to full array
Z = Z(isComplete,:,:,:,:,:);
n = nnz(isComplete);
fprintf('\nUsing %d of %d participants (complete data), %d bins in %g-%g Hz.\n', ...
        n, nSubj, nF, F_BAND_HZ(1), F_BAND_HZ(2));

% Whole-curve weights based on baseline condition
switch CURVE_WEIGHTING
    case 'none'
        w_curve = ones(nF, 1);
    case 'invVar'
        % Use the same valid participants for inverse-variance weighting
        Z_base_complete = Z(:,:,:,idx_base_w, idx_base_t, idx_base_d);
        w_curve = 1 ./ var(reshape(Z_base_complete, [], nF), 0, 1)';
    otherwise
        error('Unknown CURVE_WEIGHTING');
end

%% ------------------------------------------------------------------------
%  ICC point estimates and bootstrap CIs for all combinations
%  ------------------------------------------------------------------------
fprintf('Bootstrapping ICCs for %d combinations (%d iterations)...\n', nW*nT*nD, N_BOOT);
rng(RNG_SEED);
bootIdx = randi(n, n, N_BOOT); 

iccCurve = nan(nW, nT, nD);
curveLo  = nan(nW, nT, nD);
curveHi  = nan(nW, nT, nD);

qLo = 100*ALPHA/2;  
qHi = 100*(1 - ALPHA/2);

for w = 1:nW
    for t = 1:nT
        for dt = 1:nD
            Y = Z(:,:,:,w,t,dt); % n x k x nF
            
            % Point Estimate
            S = iccStats(Y);
            iccCurve(w,t,dt) = sum(w_curve.*S.num) / sum(w_curve.*S.den);
            
            % Bootstrap
            Bc = nan(N_BOOT, 1);
            for b = 1:N_BOOT
                Sb = iccStats(Y(bootIdx(:,b),:,:));
                Bc(b) = sum(w_curve.*Sb.num) / sum(w_curve.*Sb.den);
            end
            
            pc = prctile(Bc, [qLo qHi]);
            curveLo(w,t,dt) = pc(1);
            curveHi(w,t,dt) = pc(2);
        end
    end
end
fprintf('Done.\n');

%% ------------------------------------------------------------------------
%  Save outputs
%  ------------------------------------------------------------------------

save('./Results/zIP_ICC_Sensitivity.mat', 'Z', 'fBand', 'WIN_SIZES', 'TAPERS', 'DETRENDS', ...
     'iccCurve', 'curveLo', 'curveHi', 'w_curve', 'F_BAND_HZ', 'N_BOOT');

else
load('./Results/zIP_ICC_Sensitivity.mat');
end

%% ------------------------------------------------------------------------
%  Figure Plotting
%  ------------------------------------------------------------------------
fig = figure('Units','centimeters','Position',[0 0 32 10]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
y_limits = [max(0, min(curveLo(:)) - 0.05), min(1, max(curveHi(:)) + 0.05)]; 
y_limits = [0.4,0.9];

% 1. ICC vs Window Size (Holding Taper and Detrend constant)
nexttile; hold on;
x_win = 1:nW;
y_win = squeeze(iccCurve(:, idx_base_t, idx_base_d));
l_win = squeeze(curveLo(:, idx_base_t, idx_base_d));
h_win = squeeze(curveHi(:, idx_base_t, idx_base_d));

errorbar(x_win, y_win, y_win - l_win, h_win - y_win, ...
    'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto', 'Color', '#0072BD');
xticks(x_win);
xticklabels({'2^8', '2^9', '2^{10}', '2^{11}'});
xlabel('Window Size', 'FontSize', 18);
ylabel('ICC_{curve}(2,1)', 'FontSize', 18);
text(-0.24, 1.08, 'A', 'Units','normalized', 'FontWeight','bold', ...
         'FontSize',20, 'VerticalAlignment','top');
ylim(y_limits); grid on; box on;

% 2. ICC vs Taper Option (Holding Window Size and Detrend constant)
nexttile; hold on;
x_tap = 1:nT;
y_tap = squeeze(iccCurve(idx_base_w, :, idx_base_d))';
l_tap = squeeze(curveLo(idx_base_w, :, idx_base_d))';
h_tap = squeeze(curveHi(idx_base_w, :, idx_base_d))';

errorbar(x_tap, y_tap, y_tap - l_tap, h_tap - y_tap, ...
    's-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto', 'Color', '#D95319');
xticks(x_tap);
xticklabels({'Rectangular', 'Hann', 'Hamming'});
xtickangle(0);
xlabel('Taper Option', 'FontSize', 18);
text(-0.24, 1.08, 'B', 'Units','normalized', 'FontWeight','bold', ...
         'FontSize',20, 'VerticalAlignment','top');
ylim(y_limits); grid on; box on;

% 3. ICC vs Detrend Option (Holding Window Size and Taper constant)
nexttile; hold on;
x_det = 1:nD;
y_det = squeeze(iccCurve(idx_base_w, idx_base_t, :));
l_det = squeeze(curveLo(idx_base_w, idx_base_t, :));
h_det = squeeze(curveHi(idx_base_w, idx_base_t, :));

errorbar(x_det, y_det, y_det - l_det, h_det - y_det, ...
    '^-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto', 'Color', '#EDB120');
xticks(x_det);
xticklabels({'None', 'Zero Mean', 'Linear'});
xtickangle(0);
xlabel('Detrend Option', 'FontSize', 18);
text(-0.24, 1.08, 'C', 'Units','normalized', 'FontWeight','bold', ...
         'FontSize',20, 'VerticalAlignment','top');
ylim(y_limits); grid on; box on;


%% ========================================================================
%  Local functions
%  ========================================================================
function S = iccStats(Y)
% Two-way ANOVA mean squares and ICCs, per bin.
% Y: n (participants) x k (trials) x nF (bins). Outputs are nF x 1.
    [n, k, ~] = size(Y);
    m   = mean(Y, [1 2]);                    % grand mean     (1x1xnF)
    mi  = mean(Y, 2);                        % participant means (nx1xnF)
    mj  = mean(Y, 1);                        % trial means    (1xkxnF)
    SST = sum((Y - m).^2,  [1 2]);
    SSR = k * sum((mi - m).^2, 1);           % between participants
    SSC = n * sum((mj - m).^2, 2);           % between trials (systematic)
    SSE = SST - SSR - SSC;                   % residual noise

    S.MSR = squeeze(SSR) / (n - 1);
    S.MSC = squeeze(SSC) / (k - 1);
    S.MSE = squeeze(SSE) / ((n - 1)*(k - 1));

    % ICC(A,1) = num/den
    S.num = S.MSR - S.MSE;
    S.den = S.MSR + (k - 1)*S.MSE + (k/n)*(S.MSC - S.MSE);
end

function d = formatTrial(T, plane, height_m, height2com, cm2m, fs, nSamp, demeanFh)
% Format the first nSamp samples of a trial table for getZIPfromData.
    idx = 1:nSamp;
    switch plane
        case 'sgt'
            cop_m = T.COPx_cm_(idx)' * cm2m;
            fh    = T.Fx_N_(idx)';
        case 'frt'
            cop_m = T.COPy_cm_(idx)' * cm2m;
            fh    = T.Fy_N_(idx)';
        otherwise
            error('plane must be ''sgt'' or ''frt''');
    end
    fz   = T.Fz_N_(idx)';
    zCoM = height2com * height_m;

    if demeanFh, fh = fh - mean(fh); end

    d.FootForce   = [fh; fz];
    d.COM         = zCoM;
    d.COP         = cop_m;
    d.sampFreq_Hz = fs;
end

% Label helper functions
function n = getWinName(idx)
    names = {'2^8', '2^9', '2^{10}', '2^{11}'};
    n = names{idx};
end
function n = getTaperName(idx)
    names = {'Rect.', 'Hann', 'Hamming'};
    n = names{idx};
end
function n = getDetrendName(idx)
    names = {'None', 'Zero Mean', 'Linear'};
    n = names{idx};
end