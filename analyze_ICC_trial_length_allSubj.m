% analyze_ICC_trial_length_allSubj_v2.m
%
% Test-retest reliability of single-trial zIP curves vs. trial length
% (Santos and Duarte 2016; vision open, firm surface).
%
% For each trial length L (first L s of each of N_TRIALS trials):
%   - one zIP curve per trial (trials are NOT pooled)
%   - per-bin ICC(A,1) [two-way random, absolute agreement, single trial]
%       * point estimate
%       * analytic 95% CI (McGraw & Wong 1996, Case 2A), as a cross-check
%       * bootstrap pointwise 95% CI
%       * bootstrap simultaneous 95% band (max-t, Fisher-z scale)
%   - whole-curve ICC(A,1) = sum_f w_f*num_f / sum_f w_f*den_f
%       * point estimate, bootstrap 95% CI
%   - Spearman-Brown: number of L-s trials needed to reach target ICC
%
% The bootstrap resamples PARTICIPANTS (each keeps all trials and the
% full curve), which preserves the correlation between frequency bins.
%
% Based on: analyze_effect_of_trial_num_dur_allSubj.m
%
% Code generated using Claude Opus 5.5 on 2026-10-25

clear; clc; close all;

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

% Fixed zIP parameters (same as reference script)
params.window_size         = 2^9;
params.window_taper_option = 2;   % Hann
params.detrend_option      = 3;   % linear detrend

TRIAL_LEN_S = [10 15 20 30 40 50 60];   % keep sorted ascending
N_TRIALS    = 3;                        % trials per participant (k)

F_BAND_HZ   = [0 6];      % bins used: f > lo & f < hi  (choose a priori!)

ALPHA       = 0.05;       % 95% intervals
N_BOOT      = 5000;
RNG_SEED    = 1;

ICC_TARGET  = [0.7 0.8 0.9];   % Spearman-Brown targets (good, excellent)

% Whole-curve weighting across bins:
%   'none'   : raw zIP units; bins with large variance dominate
%   'invVar' : each bin scaled by 1/total variance at the longest length
%              (~ average of per-bin ICCs; every bin counts equally)
CURVE_WEIGHTING = 'none';

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
nL       = numel(TRIAL_LEN_S);
k        = N_TRIALS;
fprintf('Found %d participants (vision open, firm surface).\n', nSubj);

%% ------------------------------------------------------------------------
%  Compute single-trial zIP for each participant, trial length and trial
%  ------------------------------------------------------------------------
fAll = cell(nSubj, nL, k);
zAll = cell(nSubj, nL, k);
ageGroup = cell(nSubj, 1);

for iS = 1:nSubj
    rows = find(isSel & info.Subject == subjects(iS));
    ageGroup{iS} = info.AgeGroup{rows(1)};
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

    for iL = 1:nL
        nReq = round(TRIAL_LEN_S(iL) * FS_HZ);
        for iT = 1:k
            nSamp = min(nReq, height(tables{iT}));
            if nReq - nSamp > 0.5*FS_HZ
                warning('Subject %d trial %d shorter than %g s; skipped.', ...
                        subjects(iS), iT, TRIAL_LEN_S(iL));
                continue;
            end
            d = formatTrial(tables{iT}, plane, height_m, HEIGHT2COM, ...
                            CM2M, FS_HZ, nSamp, DEMEAN_FH);
            try
                [f, z] = getZIPfromData({d}, params);   % ONE trial only
                fAll{iS,iL,iT} = f(:);
                zAll{iS,iL,iT} = z(:);
            catch ME
                warning('Subject %d, %g s, trial %d failed: %s', ...
                        subjects(iS), TRIAL_LEN_S(iL), iT, ME.message);
            end
        end
    end
    fprintf('  Subject %3d done (%d/%d)\n', subjects(iS), iS, nSubj);
end

%% ------------------------------------------------------------------------
%  Assemble array Z: participants x trials x bins x lengths
%  ------------------------------------------------------------------------
iFirst  = find(~cellfun(@isempty, fAll), 1);
fRef    = fAll{iFirst};
inBand  = fRef > F_BAND_HZ(1) & fRef < F_BAND_HZ(2);
fBand   = fRef(inBand);
nF      = numel(fBand);

Z = nan(nSubj, k, nF, nL);
for iS = 1:nSubj
    for iL = 1:nL
        for iT = 1:k
            z = zAll{iS,iL,iT};
            if isempty(z), continue; end
            f = fAll{iS,iL,iT};
            if numel(f) ~= numel(fRef) || any(abs(f - fRef) > 1e-9)
                z = interp1(f, z, fRef, 'linear', NaN);
            end
            Z(iS,iT,:,iL) = z(inBand);
        end
    end
end

% Complete cases only (same participants at every length)
isComplete = all(isfinite(Z), [2 3 4]);
Z          = Z(isComplete,:,:,:);
subjUsed   = subjects(isComplete);
n          = nnz(isComplete);
fprintf('\nUsing %d of %d participants (complete data), %d bins in %g-%g Hz.\n', ...
        n, nSubj, nF, F_BAND_HZ(1), F_BAND_HZ(2));

% Whole-curve weights (fixed across lengths and bootstrap replicates)
switch CURVE_WEIGHTING
    case 'none'
        w = ones(nF, 1);
    case 'invVar'
        [~, iLong] = max(TRIAL_LEN_S);
        w = 1 ./ var(reshape(Z(:,:,:,iLong), [], nF), 0, 1)';
    otherwise
        error('Unknown CURVE_WEIGHTING');
end

%% ------------------------------------------------------------------------
%  ICC point estimates, analytic CI, bootstrap CIs
%  ------------------------------------------------------------------------
rng(RNG_SEED);
bootIdx = randi(n, n, N_BOOT);       % same resamples for every length

iccBin   = nan(nF, nL);  iccBin11 = nan(nF, nL);
loA      = nan(nF, nL);  hiA      = nan(nF, nL);   % analytic
loBoot   = nan(nF, nL);  hiBoot   = nan(nF, nL);   % bootstrap pointwise
loSim    = nan(nF, nL);  hiSim    = nan(nF, nL);   % bootstrap simultaneous
qSim     = nan(1, nL);
iccCurve = nan(1, nL);   curveLo = nan(1, nL);   curveHi = nan(1, nL);
iccCurveBoot = nan(N_BOOT, nL);

qLo = 100*ALPHA/2;  qHi = 100*(1 - ALPHA/2);
fz  = @(r) atanh(min(max(r, -0.999), 0.999));      % Fisher-z (clipped)

for iL = 1:nL
    Y = Z(:,:,:,iL);                                % n x k x nF

    S = iccStats(Y);
    iccBin(:,iL)   = S.icc;
    iccBin11(:,iL) = S.icc11;
    [loA(:,iL), hiA(:,iL)] = iccA1_CI(S.MSR, S.MSC, S.MSE, n, k, ALPHA);
    iccCurve(iL) = sum(w.*S.num) / sum(w.*S.den);

    Bb = nan(N_BOOT, nF);  Bc = nan(N_BOOT, 1);
    for b = 1:N_BOOT
        Sb = iccStats(Y(bootIdx(:,b),:,:));
        Bb(b,:) = Sb.icc';
        Bc(b)   = sum(w.*Sb.num) / sum(w.*Sb.den);
    end

    % Pointwise percentile CI
    pr = prctile(Bb, [qLo qHi]);
    loBoot(:,iL) = pr(1,:)';  hiBoot(:,iL) = pr(2,:)';

    % Simultaneous band: max-|t| over bins, on Fisher-z scale
    G   = fz(Bb);  g0 = fz(S.icc');
    seG = std(G, 0, 1);
    T   = max(abs(G - g0) ./ seG, [], 2);
    qSim(iL)    = prctile(T, 100*(1 - ALPHA));
    loSim(:,iL) = tanh(g0 - qSim(iL)*seG)';
    hiSim(:,iL) = tanh(g0 + qSim(iL)*seG)';

    % Whole-curve percentile CI
    iccCurveBoot(:,iL) = Bc;
    pc = prctile(Bc, [qLo qHi]);
    curveLo(iL) = pc(1);  curveHi(iL) = pc(2);

    fprintf('  %2g s: whole-curve ICC = %.3f [%.3f, %.3f], q_sim = %.2f\n', ...
            TRIAL_LEN_S(iL), iccCurve(iL), curveLo(iL), curveHi(iL), qSim(iL));
end

% Sanity checks
fprintf('\nMax |ICC(A,1) - ICC(1,1)| over bins/lengths: %.3f\n', ...
        max(abs(iccBin(:) - iccBin11(:))));
fprintf('Max |analytic - bootstrap| lower limit:      %.3f\n', ...
        max(abs(loA(:) - loBoot(:))));
fprintf('Simultaneous critical value q (1.96 = bins fully correlated, ');
fprintf('%.2f = Bonferroni):\n  %s\n', norminv(1 - ALPHA/(2*nF)), num2str(qSim, '%.2f  '));

%% ------------------------------------------------------------------------
%  Spearman-Brown: number of L-s trials needed for target ICC
%  ------------------------------------------------------------------------
nTgt       = numel(ICC_TARGET);
minLoPoint = min(loBoot, [], 1)';      % worst-bin pointwise lower limit
minLoSim   = min(loSim,  [], 1)';      % worst-bin simultaneous lower limit
mCurve     = nan(nL, nTgt);
mAllBins   = nan(nL, nTgt);
for t = 1:nTgt
    mCurve(:,t)   = ceil(sbNeeded(curveLo(:),  ICC_TARGET(t)) - 1e-9);
    mAllBins(:,t) = ceil(sbNeeded(minLoPoint,  ICC_TARGET(t)) - 1e-9);
end

L = TRIAL_LEN_S(:);
summaryT = table(L, iccCurve(:), curveLo(:), curveHi(:), ...
                 min(iccBin,[],1)', minLoPoint, minLoSim, ...
    'VariableNames', {'TrialLen_s','ICCcurve','ICCcurve_lo','ICCcurve_hi', ...
                      'minBin_ICC','minBin_lo','minBin_loSimul'});
for t = 1:nTgt
    tag = sprintf('%02d', round(100*ICC_TARGET(t)));
    summaryT.(['nTrials_curve_'   tag]) = mCurve(:,t);
    summaryT.(['totTime_curve_'   tag]) = mCurve(:,t) .* L;
    summaryT.(['nTrials_allBins_' tag]) = mAllBins(:,t);
    summaryT.(['totTime_allBins_' tag]) = mAllBins(:,t) .* L;
end
fprintf('\n'); disp(summaryT);

% Shortest single-trial length whose lower limit stays >= target for all
% longer lengths (guards against a lucky crossing)
for t = 1:nTgt
    ok  = curveLo(:) >= ICC_TARGET(t);
    iOK = find(arrayfun(@(i) all(ok(i:end)), 1:nL), 1);
    if isempty(iOK)
        fprintf('Target %.2f: no single-trial length reaches it (whole-curve lower limit).\n', ...
                ICC_TARGET(t));
    else
        fprintf('Target %.2f: single trial of >= %g s suffices (whole-curve lower limit).\n', ...
                ICC_TARGET(t), TRIAL_LEN_S(iOK));
    end
end


%% ------------------------------------------------------------------------
%  Save
%  ------------------------------------------------------------------------
save('./Results/zIP_ICC_trialLength.mat', 'Z', 'fBand', 'subjUsed', 'TRIAL_LEN_S', ...
     'iccBin', 'iccBin11', 'loA', 'hiA', 'loBoot', 'hiBoot', 'loSim', 'hiSim', ...
     'qSim', 'iccCurve', 'curveLo', 'curveHi', 'iccCurveBoot', 'w', ...
     'mCurve', 'mAllBins', 'summaryT', 'params', 'F_BAND_HZ', 'ICC_TARGET', ...
     'CURVE_WEIGHTING', 'N_BOOT', 'RNG_SEED');

else

load('./Results/zIP_ICC_trialLength.mat');
L = TRIAL_LEN_S(:);
nTgt       = numel(ICC_TARGET);

end

% %% ------------------------------------------------------------------------
% %  Figure 1: per-bin ICC vs. frequency, one panel per trial length
% %  ------------------------------------------------------------------------
% nTiles = nL + 1;  nCols = 4;  nRows = ceil(nTiles / nCols);
% cmap   = parula(nL + 1);  cmap = cmap(1:nL,:);
% yLo    = max(-0.5, floor(10*min(loSim(:)))/10);
% lenLab = arrayfun(@(x) sprintf('%g s', x), TRIAL_LEN_S, 'UniformOutput', false);
% 
% fig1 = figure('Units','centimeters','Position',[0 0 32 8*nRows]);
% tl1  = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
% for iL = 1:nL
%     nexttile; hold on;
%     fill([fBand; flipud(fBand)], [loBoot(:,iL); flipud(hiBoot(:,iL))], ...
%          cmap(iL,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     plot(fBand, loSim(:,iL), '--', 'Color', cmap(iL,:), 'LineWidth', 1);
%     plot(fBand, hiSim(:,iL), '--', 'Color', cmap(iL,:), 'LineWidth', 1, ...
%          'HandleVisibility', 'off');
%     plot(fBand, iccBin(:,iL), '-', 'Color', cmap(iL,:), 'LineWidth', 1.5);
%     for t = 1:nTgt
%         yline(ICC_TARGET(t), ':k', 'HandleVisibility', 'off');
%     end
%     title(lenLab{iL}); ylim([yLo 1]); xlim([fBand(1) fBand(end)]);
%     grid on; box on;
%     if iL == 1
%         legend({'Pointwise 95% CI','Simultaneous 95% band','ICC(A,1)'}, ...
%                'Location','southeast');
%     end
% end
% nexttile; hold on;
% for iL = 1:nL
%     plot(fBand, iccBin(:,iL), '-', 'Color', cmap(iL,:), 'LineWidth', 1.5);
% end
% for t = 1:nTgt, yline(ICC_TARGET(t), ':k', 'HandleVisibility', 'off'); end
% legend(lenLab, 'Location', 'southeast', 'NumColumns', 2);
% title('All lengths (estimates)'); ylim([yLo 1]); xlim([fBand(1) fBand(end)]);
% grid on; box on;
% xlabel(tl1, 'Frequency (Hz)', 'FontSize', 18);
% ylabel(tl1, 'ICC(A,1) of single-trial z_{IP}', 'FontSize', 18);

%% ------------------------------------------------------------------------
%  Figure 2: whole-curve ICC and worst-bin lower limits vs. trial length
%  ------------------------------------------------------------------------
fig2 = figure('Units','centimeters','Position',[0 0 16 11]); hold on;
errorbar(L, iccCurve(:), iccCurve(:) - curveLo(:), curveHi(:) - iccCurve(:), ...
         'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
for t = 1:nTgt, yline(ICC_TARGET(t), ':k', 'HandleVisibility', 'off'); end
xlabel('Trial length (s)', 'FontSize', 18);
ylabel('ICC(A,1), single trial', 'FontSize', 18);
xticks(L); grid on; box on;

%% ------------------------------------------------------------------------
%  Figure 3: Spearman-Brown trials needed and total recording time
%  ------------------------------------------------------------------------
fig3 = figure('Units','centimeters','Position',[0 0 26 10]);
tl3  = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
co   = get(groot, 'defaultAxesColorOrder');
legSB = {};
for p = 1:2
    nexttile; hold on;
    for t = 1:nTgt
        scale = 1; if p == 2, scale = L; end
        plot(L, mCurve(:,t)   .* scale, 'o-',  'Color', co(t,:), 'LineWidth', 1.5);
        if p == 1
            legSB = [legSB, sprintf('ICC > %.2f', ICC_TARGET(t))]; %#ok<AGROW>
        end
    end
    xticks(L); grid on; box on;
    if p == 1
        ylabel('Trials needed (lower 95% limit)', 'FontSize', 18);
        legend(legSB, 'Location', 'northeast');
        title('Number of trials');
    else
        ylabel('Total recording time (s)', 'FontSize', 18);
        title('Trials \times trial length');
    end
end
xlabel(tl3, 'Trial length (s)', 'FontSize', 18);


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
    MSW   = squeeze(SSC + SSE) / (n*(k - 1));    % one-way "within"

    % ICC(A,1) = num/den; den/k = sigma_b^2 + sigma_c^2 + sigma_e^2
    S.num   = S.MSR - S.MSE;
    S.den   = S.MSR + (k - 1)*S.MSE + (k/n)*(S.MSC - S.MSE);
    S.icc   = S.num ./ S.den;
    S.icc11 = (S.MSR - MSW) ./ (S.MSR + (k - 1)*MSW);
end

function [lo, hi] = iccA1_CI(MSR, MSC, MSE, n, k, alpha)
% Analytic CI for ICC(A,1) (McGraw & Wong 1996, Case 2A; Satterthwaite df).
    icc = (MSR - MSE) ./ (MSR + (k - 1)*MSE + (k/n)*(MSC - MSE));
    Fj  = MSC ./ MSE;
    a   = n*(1 + (k - 1)*icc) - k*icc;
    vn  = (k - 1)*(n - 1) * (k*icc.*Fj + a).^2;
    vd  = (n - 1)*k^2 * icc.^2 .* Fj.^2 + a.^2;
    v   = vn ./ vd;
    Fu  = finv(1 - alpha/2, n - 1, v);
    Fl  = finv(1 - alpha/2, v, n - 1);
    lo  = n*(MSR - Fu.*MSE) ./ (Fu.*(k*MSC + (k*n - k - n)*MSE) + n*MSR);
    hi  = n*(Fl.*MSR - MSE) ./ (k*MSC + (k*n - k - n)*MSE + n*Fl.*MSR);
end

function m = sbNeeded(r, rTarget)
% Spearman-Brown: number of trials to average for reliability rTarget,
% given single-trial reliability r.
    m = rTarget .* (1 - r) ./ (r .* (1 - rTarget));
    m(r <= 0) = Inf;
    m = max(m, 1);
end

function d = formatTrial(T, plane, height_m, height2com, cm2m, fs, nSamp, ...
                         demeanFh)
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