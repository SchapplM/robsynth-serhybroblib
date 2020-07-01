% Calculate homogenous joint transformation matrices for
% palh3m2DE1
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
% Datum: 2020-06-30 19:21
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh3m2DE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:56:53
% EndTime: 2020-06-30 18:56:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (194->35), mult. (318->38), div. (0->0), fcn. (514->22), ass. (0->39)
t133 = sin(pkin(16));
t134 = cos(pkin(16));
t141 = sin(pkin(15));
t147 = cos(pkin(15));
t121 = t133 * t147 + t134 * t141;
t122 = -t133 * t141 + t134 * t147;
t138 = sin(qJ(3));
t144 = cos(qJ(3));
t115 = t121 * t144 + t138 * t122;
t139 = sin(qJ(2));
t145 = cos(qJ(2));
t150 = t138 * t121 - t122 * t144;
t113 = t115 * t145 - t139 * t150;
t132 = pkin(17) + pkin(18);
t130 = sin(t132);
t131 = cos(t132);
t153 = t139 * t115 + t150 * t145;
t112 = t113 * t130 + t153 * t131;
t109 = t113 * t131 - t130 * t153;
t135 = sin(pkin(18));
t136 = cos(pkin(18));
t123 = t135 * t147 + t136 * t141;
t124 = -t135 * t141 + t136 * t147;
t118 = t145 * t123 + t139 * t124;
t142 = sin(pkin(14));
t148 = cos(pkin(14));
t125 = t148 * t141 - t142 * t147;
t126 = t142 * t141 + t148 * t147;
t149 = t145 * t125 + t139 * t126;
t146 = cos(qJ(1));
t143 = cos(qJ(4));
t140 = sin(qJ(1));
t137 = sin(qJ(4));
t129 = pkin(18) + pkin(15) + qJ(3) + qJ(2);
t128 = cos(t129);
t127 = sin(t129);
t120 = -t125 * t139 + t126 * t145;
t117 = t139 * t123 - t145 * t124;
t1 = [t146, -t140, 0, 0; t140, t146, 0, 0; 0, 0, 1, pkin(11); t145, -t139, 0, pkin(12); 0, 0, -1, 0; t139, t145, 0, 0; -t144, t138, 0, pkin(1); -t138, -t144, 0, 0; 0, 0, 1, 0; -t112, t109, 0, pkin(4); -t109, -t112, 0, 0; 0, 0, 1, 0; t143, -t137, 0, pkin(8); 0, 0, -1, -pkin(10); t137, t143, 0, 0; t120, -t149, 0, -pkin(6); 0, 0, -1, 0; t149, t120, 0, pkin(13); t117, -t118, 0, pkin(1); t118, t117, 0, 0; 0, 0, 1, 0; t128, -t127, 0, cos(pkin(17)) * pkin(3); t127, t128, 0, -sin(pkin(17)) * pkin(3); 0, 0, 1, 0; t117, t118, 0, t136 * pkin(2); -t118, t117, 0, t135 * pkin(2); 0, 0, 1, 0; t112, t109, 0, t134 * pkin(9); -t109, t112, 0, t133 * pkin(9); 0, 0, 1, 0; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0; -1, 0, 0, pkin(7); 0, -1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
