% Calculate homogenous joint transformation matrices for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_mdh [4x4x16]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh1m2TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 19:58:44
% EndTime: 2020-05-01 19:58:44
% DurationCPUTime: 0.25s
% Computational Cost: add. (206->43), mult. (336->44), div. (0->0), fcn. (552->24), ass. (0->43)
t156 = sin(pkin(20));
t159 = cos(pkin(20));
t165 = sin(pkin(18));
t171 = cos(pkin(18));
t141 = -t171 * t156 + t165 * t159;
t145 = t165 * t156 + t171 * t159;
t162 = sin(qJ(3));
t168 = cos(qJ(3));
t135 = t141 * t168 - t162 * t145;
t138 = t162 * t141 + t145 * t168;
t163 = sin(qJ(2));
t169 = cos(qJ(2));
t132 = t135 * t169 - t138 * t163;
t133 = t135 * t163 + t138 * t169;
t154 = pkin(22) + pkin(21);
t152 = sin(t154);
t153 = cos(t154);
t128 = t132 * t153 - t152 * t133;
t130 = t152 * t132 + t133 * t153;
t155 = sin(pkin(22));
t158 = cos(pkin(22));
t143 = t165 * t155 + t158 * t171;
t144 = t171 * t155 - t165 * t158;
t137 = t143 * t169 - t163 * t144;
t166 = sin(pkin(17));
t172 = cos(pkin(17));
t146 = t165 * t172 - t171 * t166;
t147 = t166 * t165 + t172 * t171;
t173 = t163 * t146 + t147 * t169;
t170 = cos(qJ(1));
t167 = cos(qJ(4));
t164 = sin(qJ(1));
t161 = sin(qJ(4));
t160 = cos(pkin(19));
t157 = sin(pkin(19));
t151 = pkin(18) - pkin(22) - qJ(3) - qJ(2);
t149 = cos(t151);
t148 = sin(t151);
t142 = -t162 * t157 + t168 * t160;
t140 = t168 * t157 + t162 * t160;
t139 = t146 * t169 - t147 * t163;
t136 = t163 * t143 + t144 * t169;
t1 = [t170, -t164, 0, 0; t164, t170, 0, 0; 0, 0, 1, pkin(13); 0, 0, 0, 1; -t163, -t169, 0, pkin(15); 0, 0, -1, 0; t169, -t163, 0, 0; 0, 0, 0, 1; t162, t168, 0, pkin(1); -t168, t162, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t130, t128, 0, pkin(5); -t128, -t130, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t167, -t161, 0, pkin(9); 0, 0, -1, -pkin(11); t161, t167, 0, 0; 0, 0, 0, 1; t139, -t173, 0, -pkin(14); 0, 0, -1, 0; t173, t139, 0, -pkin(16); 0, 0, 0, 1; t136, -t137, 0, pkin(1); t137, t136, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t140, t142, 0, 0; -t142, t140, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t140, -t142, 0, pkin(2); t142, t140, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t149, -t148, 0, cos(pkin(21)) * pkin(4); t148, -t149, 0, -sin(pkin(21)) * pkin(4); 0, 0, 1, 0; 0, 0, 0, 1; t136, t137, 0, t158 * pkin(3); -t137, t136, 0, t155 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t162, -t168, 0, t160 * pkin(6); t168, t162, 0, t157 * pkin(6); 0, 0, 1, 0; 0, 0, 0, 1; t130, t128, 0, t159 * pkin(10); -t128, t130, 0, t156 * pkin(10); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(7); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(12); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
