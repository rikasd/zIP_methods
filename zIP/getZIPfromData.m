function [f_zIP, zIP_ratio, GOF] = getZIPfromData(data,params)
% GETZIPFROMDATA Get zIP over freq from data
% Obtain intersection-point-height (zIP) at different frequencies from
% standing data
%
% Inputs:
% data = a struct containing a single trial of data, or
%        a cell with M structs (one struct per trial, for M trials of data)
%        containing the following info:
%      - COM = scalar center-of-mass (COM) height, or
%              2xN center-of-mass position in the 2D plane of interest 
%              (sagittal, frontal, etc.)
%              with first row corresponding to the horizontal component
%              and the second row to the vertical component.
%      - COP = 1xN center-of-pressure (COP) position in the 2D plane of 
%              interest.
%      - FootForce = 2xN foot-force vector in the 2D plane of interest
%              with first row corresponding to the horizontal component
%              and the second row to the vertical component.
%      - sampFreq_Hz = sampling frequency of data, in Hz
% params = struct containing parameters for zIP analysis
%      - method = string that specifies the method for computing zIP from
%                 experimental data; choose from:
%                 'cpsd_gain' (default) = gain definition
%                 'cpsd_sdratio' = original definition 
%                 'cpsd_geom' = geometric definition
%                 'cpsd_pca'  = unit-sensitive (PCA-based) definition           
%                 'bpf_geom'= band-pass-filter (BPF) method from
%                               Boehm et al. 2019, used to
%                               compute zIP as cov(thetaF,xCoP)/var(thetaF)  
%                 'bpf_pca'   = unit-sensitive (PCA-based) definition, 
%                               computed using band-pass filtering
%      - window_size = CPSD window size (default: 2^9)
%      - detrend_option = 1: no detrend, 2: zeromean, 3: linear detrend (default)
%      - window_taper_option = 1: rectangular, 2: Hann (default), 3: Hamming
%      - f_int = BPF bandwidth (default: 0.2)
%      - f_i = min frequency to compute zIP for, using BPF method
%      - f_f = max frequency to compute zIP for, using BPF method
%
% Outputs:
% f_zIP     = frequencies at which zIP is evaluated
% zIP_ratio = zIP/zCOM; intersection-point height normalized by COM height
% GOF       = goodness of fit, 
%             for gain: coherence
%             for geom method: squared correlation coefficient
%             for PCA-based methods:variance accounted for by principal eigenvector
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2024-02-06
% Last updated: 2026-07-30

% CPSD detrend settings
iRAW = 1;
iZEROMEAN = 2;
iDETREND = 3;

% CPSD window taper settings
iRECTWIN = 1;
iHANNWIN = 2;
iHAMMWIN = 3;

iX = 1; % index of horizontal direction
iZ = 2; % index of vertical direction

% Default parameters
method = 'cpsd_gain';
window_size = 2^9;
detrend_option = iDETREND;
window_taper_option = iHANNWIN;
f_int = 0.2; % BPF frequency interval

if iscell(data)
    data_cell = data;
elseif isstruct(data)
    data_cell{1} = data;
else
    ME = MException('getZIPfromData:incorrectInputType', ...
        'data argument should be of type cell or struct');
    throw(ME)
end

if isscalar(data_cell{1}.COM)
    for i = 1:length(data_cell)
        data_cell{i}.COM = ...
            [zeros(1,length(data_cell{i}.COP));...
            data_cell{i}.COM*ones(1,length(data_cell{i}.COP))];
    end
end

if nargin > 1 && isfield(params,'method')
    method = params.method;
end

if strcmp(method,'bpf_pca') || strcmp(method,'bpf_geom')
    window_size = length(data_cell{1}.COP);   % Use full data length
elseif nargin > 1 && isfield(params,'window_size')
    window_size = params.window_size;
end

if nargin > 1 && isfield(params,'detrend_option')
    detrend_option = params.detrend_option;
end

if nargin > 1 && isfield(params,'window_taper_option')
    window_taper_option = params.window_taper_option;
end

if nargin > 1 && isfield(params,'f_int')
    f_int = params.f_int;
end

if nargin > 1 && isfield(params,'f_i')
    f_i = params.f_i;
else
    f_i = 0.5-f_int/2;
end

if nargin > 1 && isfield(params,'f_f')
    f_f = params.f_f;
else
    f_f = 8;
end

switch window_taper_option
    case iRECTWIN
        windowtaper = rectwin(window_size);
    case iHANNWIN
        windowtaper = hann(window_size);
    case iHAMMWIN
        windowtaper = hamming(window_size);
end

nfft = window_size*2;
noverlap = 0.5*window_size;

% BPF frequency bins
f_bpf = (f_i:f_int:f_f)'; % frequency bins
        
Fs_Hz = data_cell{1}.sampFreq_Hz;
dt = 1/Fs_Hz;

switch method
    case {'cpsd_gain','cpsd_sdratio','cpsd_geom','cpsd_pca'}

        for iTrial = 1:length(data_cell)
            COM_z = data_cell{iTrial}.COM(iZ,:);    % vertical CoM
            Fx = data_cell{iTrial}.FootForce(iX,:); % horizontal ground reaction force            
            COP_x = data_cell{iTrial}.COP;          % 1D CoP
            Fz = data_cell{iTrial}.FootForce(iZ,:); % vertical ground reaction force
            theta_F = -Fx./Fz;
            switch detrend_option
                case iRAW
                    y = [theta_F' COP_x']; 
                case iZEROMEAN
                    y = [detrend(theta_F,0)' detrend(COP_x,0)']; 
                case iDETREND
                    y = [detrend(theta_F,1)' detrend(COP_x,1)']; 
            end   
            
            [Gyy_temp, f_zIP] = ...
                cpsd_custom(y,y,windowtaper,noverlap,nfft,'mimo',Fs_Hz,...
                detrend_option);            

            f_zIP = f_zIP(6:end);
            Gyy11(:,iTrial) = Gyy_temp(6:end,1,1);
            Gyy12(:,iTrial) = Gyy_temp(6:end,1,2);
            Gyy21(:,iTrial) = Gyy_temp(6:end,2,1);
            Gyy22(:,iTrial) = Gyy_temp(6:end,2,2);
        end
            
        N = length(f_zIP);
        zIP_ratio = zeros(N,1);
        Gcpsd_phase = zeros(N,1);

        switch method
            case 'cpsd_gain'
                GOF = zeros(N,1);
                for i = 1:(N)
                    Gyy(1,1) = mean(Gyy11(i,:),2);
                    Gyy(1,2) = mean(Gyy12(i,:),2);
                    Gyy(2,1) = mean(Gyy21(i,:),2);
                    Gyy(2,2) = mean(Gyy22(i,:),2);
                    zIP_ratio(i) = abs(Gyy(1,2))/Gyy(1,1) / mean(COM_z);
                    GOF(i) = Gyy(1,2)*Gyy(2,1)/Gyy(1,1)/Gyy(2,2);
                end
            case 'cpsd_sdratio'
                GOF = zeros(N,1);
                for i = 1:(N)
                    Gyy(1,1) = mean(Gyy11(i,:),2);
                    Gyy(2,2) = mean(Gyy22(i,:),2);
                    zIP_ratio(i) = sqrt(Gyy(2,2)/Gyy(1,1)) / mean(COM_z);
                end
            case 'cpsd_geom'
                GOF = zeros(N,1);
                for i = 1:(N)
                    Gyy(1,1) = mean(Gyy11(i,:),2);
                    Gyy(1,2) = mean(Gyy12(i,:),2);
                    Gyy(2,1) = mean(Gyy21(i,:),2);
                    Gyy(2,2) = mean(Gyy22(i,:),2);
                    zIP_ratio(i) = real(Gyy(1,2))/real(Gyy(1,1)) / mean(COM_z);
                    GOF(i) = real(Gyy(1,2))^2/real(Gyy(1,1))/real(Gyy(2,2));
                end
            case 'cpsd_pca'
                GOF = zeros(N,1);
                for i = 1:(N)
                    Gyy(1,1) = mean(Gyy11(i,:),2);
                    Gyy(1,2) = mean(Gyy12(i,:),2);
                    Gyy(2,1) = mean(Gyy21(i,:),2);
                    Gyy(2,2) = mean(Gyy22(i,:),2);
                    [V,D] = eig((real(Gyy(:,:))));
                    [~,i_max] = max(diag(abs(D)));
                    zIP_ratio(i) = V(2,i_max)/V(1,i_max) / mean(COM_z);
                    GOF(i) = max(diag(abs(D)))/sum(diag(abs(D)));
                end
        end

    case {'bpf_pca', 'bpf_geom'}
        
        for iTrial = 1:length(data_cell)
            COM_z = data_cell{iTrial}.COM(iZ,:);
            COP_x = data_cell{iTrial}.COP(iX,:);
            Fx = data_cell{iTrial}.FootForce(iX,:);
            Fz = data_cell{iTrial}.FootForce(iZ,:);
            theta_F = -Fx./Fz;
    
            zIP = zeros(length(f_bpf),1);

            switch detrend_option
                case iRAW
                    y = [theta_F' COP_x']; 
                case iZEROMEAN
                    y = [detrend(theta_F,0)' detrend(COP_x,0)']; 
                case iDETREND
                    y = [detrend(theta_F,1)' detrend(COP_x,1)']; 
            end 

            cop_tapered = y(:,2).*windowtaper; 
            theta_tapered = y(:,1).*windowtaper;
            cop_padded = [ zeros(floor(length(cop_tapered)/2),1); 
                   cop_tapered;
                   zeros(floor(length(cop_tapered)/2),1)];
            theta_padded = [ zeros(floor(length(theta_tapered)/2),1); 
                   theta_tapered;
                   zeros(floor(length(theta_tapered)/2),1)];
            for f = 1:length(f_bpf)
                [B,A] = butter(2,[f_bpf(f) f_bpf(f)+f_int]/(1/dt/2)); % BPF
                COP_f(f,iTrial,:) = filtfilt(B,A,cop_padded);
                theta_f(f,iTrial,:) = filtfilt(B,A,theta_padded);
            end % f

            COM_z_allTrials(iTrial) = mean(COM_z);
        end % iTrials
        
        switch method
            case 'bpf_pca'
                for f = 1:length(f_bpf)
                    COP_f_allTrials = [];
                    theta_f_allTrials = [];
                    for iTrial = 1:length(data_cell)
                        COP_f_allTrials = [COP_f_allTrials; squeeze(COP_f(f,iTrial,:))];
                        theta_f_allTrials = [theta_f_allTrials; squeeze(theta_f(f,iTrial,:))];
                    end % iTrials
                    [coeff,~,~,~,explained,~] = pca([theta_f_allTrials COP_f_allTrials]);
                    zIP(f) = coeff(2,1)/coeff(1,1);
                    GOF(f) = explained(1)/100;
                end % f
            case 'bpf_geom'
                for f = 1:length(f_bpf)
                    COP_f_allTrials = [];
                    theta_f_allTrials = [];
                    for iTrial = 1:length(data_cell)
                        COP_f_allTrials = [COP_f_allTrials; squeeze(COP_f(f,iTrial,:))];
                        theta_f_allTrials = [theta_f_allTrials; squeeze(theta_f(f,iTrial,:))];
                    end % iTrials
                    COV_MAT = cov(theta_f_allTrials, COP_f_allTrials);
                    zIP(f) = COV_MAT(2,1)/COV_MAT(1,1);  
                    GOF(f) = COV_MAT(1,2)^2/COV_MAT(1,1)/COV_MAT(2,2);
                    [xcov_r,xcov_lag]=xcov(COP_f_allTrials,theta_f_allTrials);
                    lag_max = 1/(f_bpf(f) + f_int/2)/4/dt;
                    indices_to_search = round(length(xcov_r)/2+(-lag_max:lag_max));
                    lag_to_search = xcov_lag(indices_to_search);
                    [~,ind_max] = max(xcov_r(indices_to_search));
                    Gcpsd_phase(f) = (dt*lag_to_search(ind_max))*(f_bpf(f) + f_int/2)*360;

                end % f
        end

        zIP_ratio = zIP/mean(COM_z_allTrials);
        f_zIP = f_bpf + f_int/2;

end % switch method

end