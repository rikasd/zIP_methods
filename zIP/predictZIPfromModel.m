function [f, zIP_ratio, GOF] = predictZIPfromModel(input_struct)
% PREDICTZIPFROMMODEL Compute closed-form zIP estimate from model
%
% Compute the closed-form estimate of the intersection-point height (zIP)
% over frequency from a double-inverted-pendulum model of bipedal stance
% based on the cross-spectral-density (CPSD) matrix between 
% center of pressure (COP) and foot-force orientation (theta_F)
%
% Inputs:
% input_struct = struct of simulation parameters, including:
%   'lumped_params'     : struct of lumped parameters for DIP
%     'm1' : mass of link 1 (ankle to hip)
%     'm2' : mass of link 2 (above hip)
%     'c1' : distance from ankle joint to center of mass of link 1
%     'c2' : distance from ankle joint to center of mass of link 2
%     'j1' : moment of inertia of link 1 about its center of mass
%     'j2' : moment of inertia of link 2 about its center of mass
%     'L1' : length of link 1
%     'L_COM': height of net center of mass 
%     'h_foot' (optional): height of foot (distance of ankle joint off the
%            ground)
%   'controller_params' : struct containing controller parameters,
%       either
%       'K_x': 2x4 feedback gains
%       or
%       LQR cost function matrices
%       'Q'  : 4x4 state penalty matrix
%       'R'  : 2x2 control-input penalty matrix
%   'motorNoiseRatio' : ankle actuator noise / hip actuator noise,
%                        default 1.6
%   'f'  : frequencies at which to evaluate zIP
%   'ip_method : method for computing zIP; 'bpf','cpsd','cpsd_covvar',etc.
%
% Outputs:
% f         = frequencies at which zIP is evaluated
% zIP_ratio = zIP/zCOM; intersection-point height normalized by COM height
% H         = closed-loop transfer function of system
% GOF       = goodness of fit, 
%             for PCA-based methods: VAF - variance accounted for by principal eigenvector
%             for covvar method: SCC - squared correlation coefficient
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2024-02-06
% Last updated: 2026-07-30


%% Load params

if isfield(input_struct,'lumped_params')
    lumped_params = input_struct.lumped_params;
else
    lumped_params = getLumpedParams_DIP(70,170,'F','sgt','pose_I',1);
end

if isfield(input_struct,'controller_params')
    controller_params = input_struct.controller_params;
else
    alpha = 1e6; beta = 0.7; omega = -0.87;
    controller_params.Q = eye(4);
    controller_params.R = alpha*[beta omega; omega 1/beta];
end

if isfield(input_struct,'motorNoiseRatio')
    motorNoiseRatio = input_struct.motorNoiseRatio;
else
    motorNoiseRatio = 0.6;
end

if isfield(input_struct,'f')
    f = input_struct.f;
else
    f = 0.1:0.1:6;
end

if isfield(input_struct,'ip_method')
    ip_method = input_struct.ip_method;
else
    ip_method = 'cpsd_gain';
end

w = 2*pi*f;

%% Generate model dynamics
% linearize dynamics
[A_CL,B_CL,C_CL,D_CL] = getModel_linDIP(lumped_params,controller_params);

s = tf('s');
H = C_CL * inv(s*eye(size(A_CL)) - A_CL) * B_CL + D_CL;

Gww = [motorNoiseRatio^2 0; 0 1];

%% Get frequency response

[H_mag,H_phase,~] = bode(H,w);
[HT_mag,HT_phase,~] = bode(Gww*(H.'),w);
HT_frf = HT_mag.*exp(1j*pi/180*HT_phase);
HC_frf = H_mag.*exp(-1j*pi/180*H_phase);

%% Compute zIP
zIP_ratio = f*0;

switch ip_method
    case 'cpsd_gain'      
        SCC = f*0;
        for i = 1:length(w)
            Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i);
            zIP_ratio(i) = abs(Gcpsd_frf(1,2,i))/abs(Gcpsd_frf(1,1,i));
            SCC(i) = abs(Gcpsd_frf(1,2,i))^2/(Gcpsd_frf(1,1,i))/(Gcpsd_frf(2,2,i));
        end
        GOF = SCC;
    case 'cpsd_sdratio'
        for i = 1:length(w)
            Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i);
            zIP_ratio(i) = sqrt(real(Gcpsd_frf(2,2,i))/real(Gcpsd_frf(1,1,i))) / lumped_params.L_COM;   
        end
        GOF = [];
    case {'cpsd_geom','bpf_geom'}
        SCC = f*0;
        for i = 1:length(w)
            Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i);
            zIP_ratio(i) = real(Gcpsd_frf(1,2,i))/real(Gcpsd_frf(1,1,i)) / lumped_params.L_COM;
            SCC(i) = real(Gcpsd_frf(1,2,i))^2/real(Gcpsd_frf(1,1,i))/real(Gcpsd_frf(2,2,i));
        end
        GOF = SCC;
    case {'cpsd_pca','bpf_pca'} 
        VAF = f*0;
        for i = 1:length(w)
            Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i);
            [V,D] = eig(real(Gcpsd_frf(:,:,i)));
            [~,i_max] = max(diag(abs(D)));
            zIP_ratio(i) = V(2,i_max)/V(1,i_max) / lumped_params.L_COM;
            VAF(i) = max(diag(abs(D)))/sum(diag(abs(D)));
        end
        GOF = VAF;
end

end