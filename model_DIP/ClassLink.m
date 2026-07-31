% Author: Jongwoo Lee, PhD
classdef ClassLink
    properties
        L_ratio_F
        L_ratio_M
        m_ratio_F
        m_ratio_M
        c_ratio_F
        c_ratio_M
        r_ratio_frt_F
        r_ratio_frt_M
        r_ratio_sgt_F
        r_ratio_sgt_M
        
        m % mass
        L % length
        c % com position from joint
        j % moment of inertia about com
        
        KE % kinetic energy
        PE % potential energy
        SE2_joint % SE(2) representation       
        SE2_com        
    end
    methods
        function obj = ClassLink(L_ratio_F, L_ratio_M, m_ratio_F, m_ratio_M, c_ratio_F, c_ratio_M, r_ratio_frt_F, r_ratio_frt_M, r_ratio_sgt_F, r_ratio_sgt_M)
            if nargin > 1
                obj.L_ratio_F = L_ratio_F/100;
                obj.L_ratio_M = L_ratio_M/100;
                obj.m_ratio_F = m_ratio_F/100;
                obj.m_ratio_M = m_ratio_M/100;
                obj.c_ratio_F = c_ratio_F/100;
                obj.c_ratio_M = c_ratio_M/100;
                obj.r_ratio_frt_F = r_ratio_frt_F/100;
                obj.r_ratio_frt_M = r_ratio_frt_M/100;
                obj.r_ratio_sgt_F = r_ratio_sgt_F/100;
                obj.r_ratio_sgt_M = r_ratio_sgt_M/100;
            end
        end
        function obj = assignLink(obj, totalMass, totalHeight, gender, plane)
            switch gender
                case 'F'
                    obj.m = obj.m_ratio_F*totalMass;
                    obj.L = obj.L_ratio_F*totalHeight;
                    obj.c = obj.c_ratio_F*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_F)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_F)^2;
                    end
                case 'M'
                    obj.m = obj.m_ratio_M*totalMass;
                    obj.L = obj.L_ratio_M*totalHeight;
                    obj.c = obj.c_ratio_M*obj.L;
                    switch plane
                        case 'sgt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_sgt_M)^2;
                        case 'frt'
                            obj.j = obj.m*(obj.L*obj.r_ratio_frt_M)^2;
                    end
            end
        end
        function assignSE2(SE2_joint, SE2_com)
            obj.SE2_joint = SE2_joint;
            obj.SE2_com = SE2_com;            
        end        
    end
end
