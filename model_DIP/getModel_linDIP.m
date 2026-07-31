function [A_cl,B_cl,C_cl,D_cl,K_x] = getModel_linDIP(lumped_params,controller_params)
%% GETMODEL_LINDIP Get full-state-feedback closed-loop (cl) state-space model 
%               of linearized double-inverted pendulum (DIP) with 
%               linear quadratic regulator (LQR) controller 
%
% Inputs:
% lumped_params = struct containing lumped parameters for DIP
%   'm1' : mass of link 1 (ankle to hip)
%   'm2' : mass of link 2 (above hip)
%   'c1' : distance from ankle joint to center of mass of link 1
%   'c2' : distance from ankle joint to center of mass of link 2
%   'j1' : moment of inertia of link 1 about its center of mass
%   'j2' : moment of inertia of link 2 about its center of mass
%   'L1' : length of link 1
%   'h_foot' (optional): height of foot (distance of ankle joint off the
%          ground)
% controller_params = struct containing controller parameters,
%   either
%   'K_x': 2x4 feedback gains
%   or
%   LQR cost function matrices
%   'Q'  : 4x4 state penalty matrix
%   'R'  : 2x2 control-input penalty matrix
%
% Outputs:
% State-space model matrices A_cl, B_cl, C_cl, D_cl for closed-loop system
%   A = dynamics matrix
%   B = input matrix
%   C = output matrix
%   D = direct feedthrough matrix
%
% Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2024-02-10
% Reference: Shiozawa et al. 2021, Appendices 2-4
% Last updated: 2026-06-05

G_MPS = 9.81; % acceleration due to gravity

m1 = lumped_params.m1;
m2 = lumped_params.m2;
c1 = lumped_params.c1;
c2 = lumped_params.c2;
j1 = lumped_params.j1;
j2 = lumped_params.j2;
L1 = lumped_params.L1;

if isfield(lumped_params,'h_foot')
    h_foot = lumped_params.h_foot;
else
    h_foot = 0;
end

M = [j1 + m1*c1^2 + j2 + m2*c2^2 + m2*L1*(L1 + 2*c2), j2 + m2*c2^2 + m2*L1*c2;
    j2 + m2*c2^2 + m2*L1*c2, j2 + m2*c2^2];
J_G = -G_MPS*[m1*c1 + m2*(L1 + c2), m2*c2;
              m2*c2, m2*c2];

% open-loop system state-space matrices
A_ol = [zeros(2) eye(2); -M\J_G zeros(2)];
B_ol = [zeros(2); M\eye(2)];

% LQR controller gains
if isfield(controller_params,'K_x')
    K_x = controller_params.K_x;
else
    K_x = lqr(A_ol,B_ol,controller_params.Q,controller_params.R); 
end

% closed-loop system state-space matrices
A_cl = A_ol - B_ol*K_x;
B_cl = B_ol;

C1 = 1/G_MPS/(m1 + m2)*[0 0 (m1*c1 + m2*L1 + m2*c2) m2*c2]*A_cl;
C2 = 1/G_MPS/(m1+m2)*(-[1 0]*K_x + h_foot*[0 0 (m1*c1 + m2*L1 + m2*c2) m2*c2]*A_cl);
D1 = 1/G_MPS/(m1 + m2)*[0 0 (m1*c1 + m2*L1 + m2*c2) m2*c2]*B_cl;
D2 = 1/G_MPS/(m1+m2)*([1 0] + h_foot*[0 0 (m1*c1 + m2*L1 + m2*c2) m2*c2]*B_cl);

C_cl = [C1; C2];
D_cl = [D1; D2];

end