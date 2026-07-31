function lumped_params = getLumpedParams_DIP(totalMass_kg,totalHeight_m,gender,plane,pose,agegroup)
%% Compute lumped parameters for double inverted pendulum model
%
% Input: - totalMass_kg  = total mass of human, in kilograms
%        - totalHeight_m = total height of human, in meters
%        - gender        = 'F': female, 'M': male
%        - plane         = 'frt': frontal plane, 'sgt': sagittal plane
%        - pose          = 'pose_I': arms adducted, 'pose_T': arms abducted
%        - agegroup      = 1: young, 2: old
% Output:- lumped_params = struct containing the lumped parameters of the
%                          two links, including:
%                          m1: mass of link 1 (lower body)
%                          m2: mass of link 2 (upper body)
%                          c1: center of mass position of link 1 from joint
%                          c2: center of mass position of link 2 from joint
%                          j1: mass moment of inertia of link 1 about c1
%                          j2: mass moment of inertia of link 2 about c2
%                          L1: length of link 1
%                          L2: length of link 2
%
% References:
% - Used Dynamics_gen_JW/sub_Gen_Parameters.m and ClassLink.m from code by
%   Jongwoo Lee as reference.
% - All parameters are based on 
%   De Leva, Paolo. "Adjustments to Zatsiorsky-Seluyanov's segment inertia
%   parameters." Journal of biomechanics 29.9 (1996): 1223-1230.
% - Shoulder width for frontal plane model parameters were obtained from 
%   NASA man-systems integration standards,
%   https://msis.jsc.nasa.gov/sections/section03.htm
%
% Author: Rika Sugimoto Dimitrova (rikasd@mit.edu)
% Date:   2025-04-08
% Modified:
% 2026-05-28 - anthropomorphic data differentiate for different age groups
%              by adding agegroup as input, and getLinks function to
%              get the corresponding link parameters
%              (Note: l_sjc and w_sjc assumed to be the same for young and
%              old)

YOUNG = 1;
OLD = 2;

if nargin < 6
    agegroup = YOUNG;
end
[ClassUpperarm,ClassForearm,ClassHand,ClassThigh,ClassHead,ClassTrunk,ClassShank]...
    = getLinks(agegroup);

Upperarm = ClassUpperarm.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Forearm = ClassForearm.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Hand = ClassHand.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Thigh = ClassThigh.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Head = ClassHead.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Trunk = ClassTrunk.assignLink(totalMass_kg, totalHeight_m, gender, plane);
Shank = ClassShank.assignLink(totalMass_kg, totalHeight_m, gender, plane);

% shoulder longitudinal and transverse distances from the hip joint
l_sjc = L_SJC(totalHeight_m, gender);
w_sjc = w_SJC(totalHeight_m, gender, plane);

link1.m = 2*(Shank.m + Thigh.m);
link2.m = 2*(Upperarm.m + Forearm.m + Hand.m) + Trunk.m + Head.m;

link1.L = Shank.L + Thigh.L;
link2.L = Trunk.L + Head.L;

link1.c = (2*Shank.m*Shank.c + 2*Thigh.m*(Shank.L+Thigh.c))/link1.m;
link1.j = 2*(Shank.j + Shank.m*(Shank.c-link1.c)^2 + ...
    Thigh.j + Thigh.m*((Shank.L+Thigh.c)-link1.c)^2);

switch pose
    case 'pose_I'
        link2.c = (2*Upperarm.m*(l_sjc-Upperarm.c) + ...
            2*Forearm.m*(l_sjc-Upperarm.L-Forearm.c) + ...
            2*Hand.m*(l_sjc-Upperarm.L-Forearm.L-Hand.c) + ...
            Trunk.m*Trunk.c + Head.m*(Trunk.L+Head.c) )/link2.m;
        link2.j = 2*Upperarm.j + 2*Upperarm.m*((l_sjc-Upperarm.c-link2.c)^2+w_sjc^2) + ...
            2*Forearm.j + 2*Forearm.m*((l_sjc-Upperarm.L-Forearm.c-link2.c)^2+w_sjc^2) + ...
            2*Hand.j + 2*Hand.m*((l_sjc-Upperarm.L-Forearm.L-Hand.c-link2.c)^2+w_sjc^2) + ...
            Trunk.j + Trunk.m*(Trunk.c-link2.c)^2 + ...
            Head.j + Head.m*(Trunk.L+Head.c-link2.c)^2;
    case 'pose_T'
        link2.c = (2*(Upperarm.m + Forearm.m + Hand.m)*(l_sjc) + ...
            Trunk.m*Trunk.c + Head.m*(Trunk.L+Head.c) )/link2.m;
        switch plane
            case 'frt'
                link2.j = 2*Upperarm.j + 2*Upperarm.m*((l_sjc-link2.c)^2+(w_sjc+Upperarm.c)^2) + ...
                    2*Forearm.j + 2*Forearm.m*((l_sjc-link2.c)^2+(w_sjc+Upperarm.L+Forearm.c)^2) + ...
                    2*Hand.j + 2*Hand.m*((l_sjc-link2.c)^2+(w_sjc+Upperarm.L+Forearm.L+Hand.c)^2) + ...
                    Trunk.j + Trunk.m*(Trunk.c-link2.c)^2 + ...
                    Head.j + Head.m*(Trunk.L+Head.c-link2.c)^2;
            case 'sgt'
                link2.j = 2*Upperarm.j + 2*Upperarm.m*(l_sjc-link2.c)^2 + ...
                    2*Forearm.j + 2*Forearm.m*(l_sjc-link2.c)^2 + ...
                    2*Hand.j + 2*Hand.m*(l_sjc-link2.c)^2 + ...
                    Trunk.j + Trunk.m*(Trunk.c-link2.c)^2 + ...
                    Head.j + Head.m*(Trunk.L+Head.c-link2.c)^2;
        end
end

lumped_params.m1 = link1.m;
lumped_params.m2 = link2.m;
lumped_params.c1 = link1.c;
lumped_params.c2 = link2.c;
lumped_params.j1 = link1.j;
lumped_params.j2 = link2.j;
lumped_params.L1 = link1.L;
lumped_params.L2 = link2.L;

%%
function l_sjc = L_SJC(totalHeight, gender)
    switch gender
        case 'F'
            l_sjc = 497.9/1735 * totalHeight;
        case 'M'
            l_sjc = 515.5/1741 * totalHeight;
    end
end
function w_sjc = w_SJC(totalHeight, gender, plane)
% obtained from NASA standards, averaged 5th, 50th, 95th percentile data
    switch plane
        case 'sgt'
            w_sjc = 0;
        case 'frt'
            switch gender
                case 'F'
                    w_sjc = 0.23/2 * totalHeight;
                    % Japanese, 40 yo, 2000
                case 'M'
                    w_sjc = 0.22/2 * totalHeight;
                    % American, 40 yo, 2000
            end
    end
end

end