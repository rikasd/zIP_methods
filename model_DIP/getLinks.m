function [ClassUpperarm,ClassForearm,ClassHand,ClassThigh,ClassHead,ClassTrunk,ClassShank]...
    = getLinks(agegroup)
%GETLINKS Create links for the different body segments
%   For younger adults based on de Leva 1996, or 
%   for older adults based on Ho Hoang and Mombaur 2015.
% 
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% Reference: original code from Jongwoo Lee, PhD
% 2026-05-27
YOUNG = 1;
OLD = 2;

if agegroup == YOUNG
    HeightF_mm = 1735;
    HeightM_mm = 1741;
    % Anthropomorphic measurements from De Leva 1996, Table 4
    ClassUpperarm = ClassLink(275.1/HeightF_mm*100, 281.7/HeightM_mm*100, 2.55, 2.71, 57.54, 57.72, 27.8, 28.5, 26.0, 26.9);
    ClassForearm = ClassLink(264.3/HeightF_mm*100, 268.9/HeightM_mm*100, 1.38, 1.62, 45.59, 45.74, 26.1, 27.6, 25.7, 26.5);
    ClassHand = ClassLink(78.0/HeightF_mm*100, 86.2/HeightM_mm*100, 0.56, 0.61, 74.74, 79.00, 53.1, 62.8, 45.4, 51.3);
    
    ClassThigh = ClassLink(368.5/HeightF_mm*100, 422.2/HeightM_mm*100, 14.78, 14.16, 100-36.12, 100-40.95, 36.9, 32.9, 36.4, 32.9);
    
    ClassHead = ClassLink(243.7/HeightF_mm*100, 242.9/HeightM_mm*100, 6.68, 6.94, 100-48.41, 100-50.02, 27.1, 30.3, 29.5, 31.5);
    ClassTrunk = ClassLink(614.8/HeightF_mm*100, 603.3/HeightM_mm*100, 42.57, 43.46, 100-49.64, 100-51.38, 30.7, 32.8, 29.2, 30.6);
    
    ClassShank = ClassLink(438.6/HeightF_mm*100, 440.3/HeightM_mm*100, 4.81, 4.33, 100-43.52, 100-43.95, 26.7, 25.1, 26.3, 24.6);

elseif agegroup == OLD
    % Anthropomorphic measurements from Ho Hoang and Mombaur 2015, Table 4
    % Note: female arm length percentages recomputed from mm values because
    %       there seemed to be a discrepancy... 
    ClassUpperarm = ClassLink(16.33, 16.67, 3.24, 3.22, 6.94/16.33*100, 7.05/16.67*100, 4.54/16.33*100, 4.75/16.67*100, 4.25/16.33*100, 4.48/16.67*100);
    ClassForearm = ClassLink(15.69, 15.91, 1.12, 1.76, 8.54/15.69*100, 8.63/15.91*100, 4.10/15.69*100, 4.39/15.91*100, 4.03/15.69*100, 4.22/15.91*100);
    ClassHand = ClassLink(4.63, 5.10, 0.50, 0.65, 1.17/4.63*100, 1.07/5.10*100, 2.47/4.63*100, 3.20/5.10*100, 2.10/4.63*100, 2.62/5.10*100);
    
    ClassThigh = ClassLink(21.71, 24.47, 12.75, 9.48, (1 - 7.84/21.71)*100, (1 - 10.02/24.47)*100, 8.01/21.71*100, 8.05/24.47*100, 7.90/21.71*100, 8.05/24.47*100);
    ClassShank = ClassLink(25.47, 25.15, 4.35, 4.27, (1 - 11.25/25.47)*100, (1 - 11.22/25.15)*100, 6.90/25.47*100, 6.41/25.15*100, 6.80/25.47*100, 6.26/25.15*100);

    ClassHead = ClassLink(13.73, 13.82, 7.39, 6.94, 7.08/13.73*100, 6.91/13.82*100, 3.72/13.73*100, 4.19/13.82*100, 4.05/13.73*100, 4.35/13.82*100);
    
    L_UPT_rel_F = 12.84;
    L_MPT_rel_F = 11.56;
    L_pelvis_rel_F = 10.22;
    L_UPT_rel_M = 13.77;
    L_MPT_rel_M = 12.26;
    L_pelvis_rel_M = 8.29;
    m_UPT_rel_F = 15.58;
    m_MPT_rel_F = 16.55;
    m_pelvis_rel_F = 14.09;
    m_UPT_rel_M = 14.49;
    m_MPT_rel_M = 21.99;
    m_pelvis_rel_M = 15.04;
    L_com_UPT_rel_F = 6.36;
    L_com_MPT_rel_F = 6.35;
    L_com_pelvis_rel_F = 5.19;
    L_com_UPT_rel_M = 6.80;
    L_com_MPT_rel_M = 6.74;
    L_com_pelvis_rel_M = 3.22;
    r_UPT_rel_F_frt = 5.99;
    r_MPT_rel_F_frt = 5.01;
    r_pelvis_rel_F_frt = 4.43;
    r_UPT_rel_M_frt = 6.96;
    r_MPT_rel_M_frt = 5.91;
    r_pelvis_rel_M_frt = 5.10;
    r_UPT_rel_F_sgt = 4.03;
    r_MPT_rel_F_sgt = 4.09;
    r_pelvis_rel_F_sgt = 4.11;
    r_UPT_rel_M_sgt = 4.41;
    r_MPT_rel_M_sgt = 4.70;
    r_pelvis_rel_M_sgt = 4.57;
    L_trunk_rel_F = L_UPT_rel_F + L_MPT_rel_F + L_pelvis_rel_F;
    L_trunk_rel_M = L_UPT_rel_M + L_MPT_rel_M + L_pelvis_rel_M;
    m_trunk_rel_F = m_UPT_rel_F + m_MPT_rel_F + m_pelvis_rel_F;
    m_trunk_rel_M = m_UPT_rel_M + m_MPT_rel_M + m_pelvis_rel_M;
    L_com_trunk_rel_F = (m_UPT_rel_F*(L_com_UPT_rel_F + L_MPT_rel_F + L_pelvis_rel_F)...
        + m_MPT_rel_F*(L_com_MPT_rel_F + L_pelvis_rel_F)...
        + m_pelvis_rel_F*L_com_pelvis_rel_F)/m_trunk_rel_F;
    L_com_trunk_rel_M = (m_UPT_rel_M*(L_com_UPT_rel_M + L_MPT_rel_M + L_pelvis_rel_M)...
        + m_MPT_rel_M*(L_com_MPT_rel_M + L_pelvis_rel_M)...
        + m_pelvis_rel_M*L_com_pelvis_rel_M)/m_trunk_rel_M;
    j_trunk_rel_F_frt = m_UPT_rel_F*(r_UPT_rel_F_frt^2 + (L_com_trunk_rel_F - (L_com_UPT_rel_F + L_MPT_rel_F + L_pelvis_rel_F))^2)...
        + m_MPT_rel_F*(r_MPT_rel_F_frt^2 + (L_com_trunk_rel_F - (L_com_MPT_rel_F + L_pelvis_rel_F))^2) ...
        + m_pelvis_rel_F*(r_pelvis_rel_F_frt^2 + (L_com_trunk_rel_F - L_com_pelvis_rel_F)^2);
    j_trunk_rel_M_frt = m_UPT_rel_M*(r_UPT_rel_M_frt^2 + (L_com_trunk_rel_M - (L_com_UPT_rel_M + L_MPT_rel_M + L_pelvis_rel_M))^2)...
        + m_MPT_rel_M*(r_MPT_rel_M_frt^2 + (L_com_trunk_rel_M - (L_com_MPT_rel_M + L_pelvis_rel_M))^2) ...
        + m_pelvis_rel_M*(r_pelvis_rel_M_frt^2 + (L_com_trunk_rel_M - L_com_pelvis_rel_M)^2);
    j_trunk_rel_F_sgt = m_UPT_rel_F*(r_UPT_rel_F_sgt^2 + (L_com_trunk_rel_F - (L_com_UPT_rel_F + L_MPT_rel_F + L_pelvis_rel_F))^2)...
        + m_MPT_rel_F*(r_MPT_rel_F_sgt^2 + (L_com_trunk_rel_F - (L_com_MPT_rel_F + L_pelvis_rel_F))^2) ...
        + m_pelvis_rel_F*(r_pelvis_rel_F_sgt^2 + (L_com_trunk_rel_F - L_com_pelvis_rel_F)^2);
    j_trunk_rel_M_sgt = m_UPT_rel_M*(r_UPT_rel_M_sgt^2 + (L_com_trunk_rel_M - (L_com_UPT_rel_M + L_MPT_rel_M + L_pelvis_rel_M))^2)...
        + m_MPT_rel_M*(r_MPT_rel_M_sgt^2 + (L_com_trunk_rel_M - (L_com_MPT_rel_M + L_pelvis_rel_M))^2) ...
        + m_pelvis_rel_M*(r_pelvis_rel_M_sgt^2 + (L_com_trunk_rel_M - L_com_pelvis_rel_M)^2);
    ClassTrunk = ClassLink(L_trunk_rel_F, L_trunk_rel_M, ...
        m_trunk_rel_F, m_trunk_rel_M, ...
        L_com_trunk_rel_F/L_trunk_rel_F*100, L_com_trunk_rel_M/L_trunk_rel_M*100, ...
        sqrt(j_trunk_rel_F_frt/m_trunk_rel_F)/L_trunk_rel_F*100,...
        sqrt(j_trunk_rel_M_frt/m_trunk_rel_M)/L_trunk_rel_M*100,...
        sqrt(j_trunk_rel_F_sgt/m_trunk_rel_F)/L_trunk_rel_F*100,...
        sqrt(j_trunk_rel_M_sgt/m_trunk_rel_M)/L_trunk_rel_M*100);
    
    
end

end

