% Calculate homogenous joint transformation matrices for
% palh3m2TE
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

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh3m2TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:30
% EndTime: 2020-05-07 01:41:31
% DurationCPUTime: 0.21s
% Computational Cost: add. (194->35), mult. (318->38), div. (0->0), fcn. (514->22), ass. (0->39)
t129 = sin(pkin(16));
t130 = cos(pkin(16));
t137 = sin(pkin(15));
t143 = cos(pkin(15));
t117 = t129 * t143 + t130 * t137;
t118 = -t129 * t137 + t130 * t143;
t134 = sin(qJ(3));
t140 = cos(qJ(3));
t111 = t117 * t140 + t134 * t118;
t135 = sin(qJ(2));
t141 = cos(qJ(2));
t146 = t134 * t117 - t118 * t140;
t109 = t111 * t141 - t135 * t146;
t128 = pkin(17) + pkin(18);
t126 = sin(t128);
t127 = cos(t128);
t149 = t111 * t135 + t146 * t141;
t108 = t109 * t126 + t149 * t127;
t105 = t109 * t127 - t126 * t149;
t131 = sin(pkin(18));
t132 = cos(pkin(18));
t119 = t131 * t143 + t132 * t137;
t120 = -t131 * t137 + t132 * t143;
t114 = t141 * t119 + t135 * t120;
t138 = sin(pkin(14));
t144 = cos(pkin(14));
t121 = t144 * t137 - t143 * t138;
t122 = t138 * t137 + t144 * t143;
t145 = t141 * t121 + t135 * t122;
t142 = cos(qJ(1));
t139 = cos(qJ(4));
t136 = sin(qJ(1));
t133 = sin(qJ(4));
t125 = pkin(18) + pkin(15) + qJ(3) + qJ(2);
t124 = cos(t125);
t123 = sin(t125);
t116 = -t121 * t135 + t122 * t141;
t113 = t135 * t119 - t141 * t120;
t1 = [t142, -t136, 0, 0; t136, t142, 0, 0; 0, 0, 1, pkin(11); 0, 0, 0, 1; t141, -t135, 0, pkin(12); 0, 0, -1, 0; t135, t141, 0, 0; 0, 0, 0, 1; -t140, t134, 0, pkin(1); -t134, -t140, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t108, t105, 0, pkin(4); -t105, -t108, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t139, -t133, 0, pkin(8); 0, 0, -1, -pkin(10); t133, t139, 0, 0; 0, 0, 0, 1; t116, -t145, 0, -pkin(6); 0, 0, -1, 0; t145, t116, 0, pkin(13); 0, 0, 0, 1; t113, -t114, 0, pkin(1); t114, t113, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t124, -t123, 0, cos(pkin(17)) * pkin(3); t123, t124, 0, -sin(pkin(17)) * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t113, t114, 0, t132 * pkin(2); -t114, t113, 0, t131 * pkin(2); 0, 0, 1, 0; 0, 0, 0, 1; t108, t105, 0, t130 * pkin(9); -t105, t108, 0, t129 * pkin(9); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(7); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
