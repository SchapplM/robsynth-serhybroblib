% Calculate homogenous joint transformation matrices for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_mdh [4x4x12]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(12+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 21:14
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh3m2DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 19:26:00
% EndTime: 2020-06-30 19:26:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (194->35), mult. (318->38), div. (0->0), fcn. (514->22), ass. (0->39)
t148 = sin(pkin(16));
t149 = cos(pkin(16));
t156 = sin(pkin(15));
t162 = cos(pkin(15));
t136 = t148 * t162 + t149 * t156;
t137 = -t148 * t156 + t149 * t162;
t153 = sin(qJ(3));
t159 = cos(qJ(3));
t130 = t136 * t159 + t153 * t137;
t154 = sin(qJ(2));
t160 = cos(qJ(2));
t165 = t153 * t136 - t137 * t159;
t128 = t130 * t160 - t154 * t165;
t147 = pkin(17) + pkin(18);
t145 = sin(t147);
t146 = cos(t147);
t168 = t130 * t154 + t165 * t160;
t127 = t145 * t128 + t168 * t146;
t124 = t128 * t146 - t145 * t168;
t150 = sin(pkin(18));
t151 = cos(pkin(18));
t138 = t162 * t150 + t156 * t151;
t139 = -t156 * t150 + t162 * t151;
t133 = t160 * t138 + t154 * t139;
t157 = sin(pkin(14));
t163 = cos(pkin(14));
t140 = t163 * t156 - t157 * t162;
t141 = t157 * t156 + t163 * t162;
t164 = t160 * t140 + t154 * t141;
t161 = cos(qJ(1));
t158 = cos(qJ(4));
t155 = sin(qJ(1));
t152 = sin(qJ(4));
t144 = pkin(15) + qJ(3) + qJ(2) + pkin(18);
t143 = cos(t144);
t142 = sin(t144);
t135 = -t140 * t154 + t141 * t160;
t132 = t154 * t138 - t160 * t139;
t1 = [t161, -t155, 0, 0; t155, t161, 0, 0; 0, 0, 1, pkin(11); t160, -t154, 0, pkin(12); 0, 0, -1, 0; t154, t160, 0, 0; -t159, t153, 0, pkin(1); -t153, -t159, 0, 0; 0, 0, 1, 0; -t127, t124, 0, pkin(4); -t124, -t127, 0, 0; 0, 0, 1, 0; t158, -t152, 0, pkin(8); 0, 0, -1, -pkin(10); t152, t158, 0, 0; t135, -t164, 0, -pkin(6); 0, 0, -1, 0; t164, t135, 0, pkin(13); t132, -t133, 0, pkin(1); t133, t132, 0, 0; 0, 0, 1, 0; t143, -t142, 0, cos(pkin(17)) * pkin(3); t142, t143, 0, -sin(pkin(17)) * pkin(3); 0, 0, 1, 0; t132, t133, 0, t151 * pkin(2); -t133, t132, 0, t150 * pkin(2); 0, 0, 1, 0; t127, t124, 0, t149 * pkin(9); -t124, t127, 0, t148 * pkin(9); 0, 0, 1, 0; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0; -1, 0, 0, pkin(7); 0, -1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
