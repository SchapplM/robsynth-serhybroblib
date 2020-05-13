% Calculate homogenous joint transformation matrices for
% mg10hlTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% T_mdh [4x4x15]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:29
% Revision: f5729c120b58d1b0137ade1aa9321d1ea2b3cc0a (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = mg10hlTE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlTE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlTE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [17x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 12:25:05
% EndTime: 2020-04-11 12:25:05
% DurationCPUTime: 0.20s
% Computational Cost: add. (944->83), mult. (954->99), div. (144->7), fcn. (448->23), ass. (0->313)
unknown=NaN(60,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t5 = cos(pkin(15));
t7 = sin(pkin(15));
t9 = -t7 * t3 + t5 * t4;
t10 = t9 * pkin(2);
t13 = t5 * t3 + t7 * t4;
t14 = pkin(1) * t13;
t16 = 0.2e1 * pkin(2) * t14;
t17 = pkin(1) ^ 2;
t27 = sqrt(-(-t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4))) * (-t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4))));
t30 = -t13 * pkin(2) + pkin(1);
t31 = pkin(2) ^ 2;
t32 = pkin(5) ^ 2;
t33 = pkin(4) ^ 2;
t34 = -t16 + t17 + t31 - t32 + t33;
t36 = -t27 * t10 + t34 * t30;
t38 = 0.1e1 / pkin(4);
t40 = 0.1e1 / (-t16 + t17 + t31);
t41 = t40 * t38;
t45 = t34 * t10 + t27 * t30;
t48 = -t41 * t36 * t5 + t41 * t45 * t7;
t53 = t41 * t36 * t7 + t41 * t45 * t5;
t56 = cos(pkin(17));
t57 = 1 / pkin(7);
t58 = qJ(6) + pkin(8);
t60 = 0.1e1 / t58 * t57;
t61 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
t62 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
t64 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
t65 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
t68 = sqrt(-t65 * t64 * t62 * t61);
t69 = t68 * t60;
t70 = (pkin(6) ^ 2);
t71 = (pkin(7) ^ 2);
t72 = (pkin(8) ^ 2);
t74 = 2 * pkin(8) * qJ(6);
t75 = qJ(6) ^ 2;
t76 = t70 - t71 - t72 - t74 - t75;
t77 = t76 * t60;
t78 = atan2(t69, t77);
t79 = t78 + pkin(16);
t80 = cos(t79);
t82 = sin(pkin(17));
t83 = sin(t79);
t85 = -t80 * t56 - t83 * t82;
t88 = -t83 * t56 + t80 * t82;
t92 = t58 ^ 2;
t94 = 0.1e1 / t92 / pkin(6);
t96 = t57 * (t70 - t71 + t72 + t74 + t75);
t104 = -t57 * t65 * t64 * t62 * t61 * t94 + t76 * t96 * t94;
t105 = cos(pkin(16));
t112 = -t76 * t57 * t68 * t94 + t68 * t96 * t94;
t113 = sin(pkin(16));
t115 = t105 * t104 / 0.4e1 - t113 * t112 / 0.4e1;
t118 = t105 * t112 / 0.4e1 + t113 * t104 / 0.4e1;
t119 = cos(qJ(3));
t120 = sin(qJ(3));
t121 = sin(qJ(4));
t122 = cos(qJ(4));
t123 = cos(qJ(5));
t124 = sin(qJ(5));
t125 = t77 / 0.2e1;
t126 = t69 / 0.2e1;
t129 = t14 - pkin(2);
t130 = -t16 + t17 + t31 + t32 - t33;
t133 = 0.1e1 / pkin(5);
t136 = t40 * t133 * (-pkin(1) * t27 * t9 - t130 * t129) / 0.2e1;
t143 = t40 * t133 * (t130 * t9 * pkin(1) - t27 * t129) / 0.2e1;
t147 = t38 * t133 * (t16 - t17 - t31 + t32 + t33) / 0.2e1;
t150 = t27 * t38 * t133 / 0.2e1;
unknown(1,1) = t1;
unknown(1,2) = -t2;
unknown(1,3) = 0.0e0;
unknown(1,4) = 0.0e0;
unknown(2,1) = t2;
unknown(2,2) = t1;
unknown(2,3) = 0.0e0;
unknown(2,4) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.1e1;
unknown(3,4) = pkin(11);
unknown(4,1) = 0.0e0;
unknown(4,2) = 0.0e0;
unknown(4,3) = 0.0e0;
unknown(4,4) = 0.1e1;
unknown(5,1) = -t3;
unknown(5,2) = -t4;
unknown(5,3) = 0.0e0;
unknown(5,4) = pkin(13);
unknown(6,1) = 0.0e0;
unknown(6,2) = 0.0e0;
unknown(6,3) = -0.1e1;
unknown(6,4) = 0.0e0;
unknown(7,1) = t4;
unknown(7,2) = -t3;
unknown(7,3) = 0.0e0;
unknown(7,4) = 0.0e0;
unknown(8,1) = 0.0e0;
unknown(8,2) = 0.0e0;
unknown(8,3) = 0.0e0;
unknown(8,4) = 0.1e1;
unknown(9,1) = t48 / 0.2e1;
unknown(9,2) = t53 / 0.2e1;
unknown(9,3) = 0.0e0;
unknown(9,4) = pkin(1) * t5;
unknown(10,1) = -t53 / 0.2e1;
unknown(10,2) = t48 / 0.2e1;
unknown(10,3) = 0.0e0;
unknown(10,4) = pkin(1) * t7;
unknown(11,1) = 0.0e0;
unknown(11,2) = 0.0e0;
unknown(11,3) = 0.1e1;
unknown(11,4) = 0.0e0;
unknown(12,1) = 0.0e0;
unknown(12,2) = 0.0e0;
unknown(12,3) = 0.0e0;
unknown(12,4) = 0.1e1;
unknown(13,1) = t85;
unknown(13,2) = t88;
unknown(13,3) = 0.0e0;
unknown(13,4) = pkin(3) * t56;
unknown(14,1) = -t88;
unknown(14,2) = t85;
unknown(14,3) = 0.0e0;
unknown(14,4) = pkin(3) * t82;
unknown(15,1) = 0.0e0;
unknown(15,2) = 0.0e0;
unknown(15,3) = 0.1e1;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = t115;
unknown(17,2) = t118;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(7);
unknown(18,1) = -t118;
unknown(18,2) = t115;
unknown(18,3) = 0.0e0;
unknown(18,4) = 0.0e0;
unknown(19,1) = 0.0e0;
unknown(19,2) = 0.0e0;
unknown(19,3) = 0.1e1;
unknown(19,4) = 0.0e0;
unknown(20,1) = 0.0e0;
unknown(20,2) = 0.0e0;
unknown(20,3) = 0.0e0;
unknown(20,4) = 0.1e1;
unknown(21,1) = 0.0e0;
unknown(21,2) = 0.0e0;
unknown(21,3) = 0.1e1;
unknown(21,4) = pkin(12);
unknown(22,1) = t119;
unknown(22,2) = -t120;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = t120;
unknown(23,2) = t119;
unknown(23,3) = 0.0e0;
unknown(23,4) = 0.0e0;
unknown(24,1) = 0.0e0;
unknown(24,2) = 0.0e0;
unknown(24,3) = 0.0e0;
unknown(24,4) = 0.1e1;
unknown(25,1) = t121;
unknown(25,2) = t122;
unknown(25,3) = 0.0e0;
unknown(25,4) = 0.0e0;
unknown(26,1) = 0.0e0;
unknown(26,2) = 0.0e0;
unknown(26,3) = 0.1e1;
unknown(26,4) = 0.0e0;
unknown(27,1) = t122;
unknown(27,2) = -t121;
unknown(27,3) = 0.0e0;
unknown(27,4) = pkin(10);
unknown(28,1) = 0.0e0;
unknown(28,2) = 0.0e0;
unknown(28,3) = 0.0e0;
unknown(28,4) = 0.1e1;
unknown(29,1) = 0.0e0;
unknown(29,2) = 0.0e0;
unknown(29,3) = 0.1e1;
unknown(29,4) = pkin(9);
unknown(30,1) = t123;
unknown(30,2) = -t124;
unknown(30,3) = 0.0e0;
unknown(30,4) = 0.0e0;
unknown(31,1) = t124;
unknown(31,2) = t123;
unknown(31,3) = 0.0e0;
unknown(31,4) = 0.0e0;
unknown(32,1) = 0.0e0;
unknown(32,2) = 0.0e0;
unknown(32,3) = 0.0e0;
unknown(32,4) = 0.1e1;
unknown(33,1) = t125;
unknown(33,2) = -t126;
unknown(33,3) = 0.0e0;
unknown(33,4) = 0.0e0;
unknown(34,1) = t126;
unknown(34,2) = t125;
unknown(34,3) = 0.0e0;
unknown(34,4) = 0.0e0;
unknown(35,1) = 0.0e0;
unknown(35,2) = 0.0e0;
unknown(35,3) = 0.1e1;
unknown(35,4) = 0.0e0;
unknown(36,1) = 0.0e0;
unknown(36,2) = 0.0e0;
unknown(36,3) = 0.0e0;
unknown(36,4) = 0.1e1;
unknown(37,1) = t136;
unknown(37,2) = -t143;
unknown(37,3) = 0.0e0;
unknown(37,4) = -pkin(14);
unknown(38,1) = 0.0e0;
unknown(38,2) = 0.0e0;
unknown(38,3) = -0.1e1;
unknown(38,4) = 0.0e0;
unknown(39,1) = t143;
unknown(39,2) = t136;
unknown(39,3) = 0.0e0;
unknown(39,4) = 0.0e0;
unknown(40,1) = 0.0e0;
unknown(40,2) = 0.0e0;
unknown(40,3) = 0.0e0;
unknown(40,4) = 0.1e1;
unknown(41,1) = -t147;
unknown(41,2) = t150;
unknown(41,3) = 0.0e0;
unknown(41,4) = pkin(6);
unknown(42,1) = -t150;
unknown(42,2) = -t147;
unknown(42,3) = 0.0e0;
unknown(42,4) = 0.0e0;
unknown(43,1) = 0.0e0;
unknown(43,2) = 0.0e0;
unknown(43,3) = 0.1e1;
unknown(43,4) = 0.0e0;
unknown(44,1) = 0.0e0;
unknown(44,2) = 0.0e0;
unknown(44,3) = 0.0e0;
unknown(44,4) = 0.1e1;
unknown(45,1) = -t104 / 0.4e1;
unknown(45,2) = t112 / 0.4e1;
unknown(45,3) = 0.0e0;
unknown(45,4) = pkin(5);
unknown(46,1) = -t112 / 0.4e1;
unknown(46,2) = -t104 / 0.4e1;
unknown(46,3) = 0.0e0;
unknown(46,4) = 0.0e0;
unknown(47,1) = 0.0e0;
unknown(47,2) = 0.0e0;
unknown(47,3) = 0.1e1;
unknown(47,4) = 0.0e0;
unknown(48,1) = 0.0e0;
unknown(48,2) = 0.0e0;
unknown(48,3) = 0.0e0;
unknown(48,4) = 0.1e1;
unknown(49,1) = 0.0e0;
unknown(49,2) = 0.0e0;
unknown(49,3) = 0.1e1;
unknown(49,4) = qJ(6);
unknown(50,1) = -0.1e1;
unknown(50,2) = 0.0e0;
unknown(50,3) = 0.0e0;
unknown(50,4) = 0.0e0;
unknown(51,1) = 0.0e0;
unknown(51,2) = -0.1e1;
unknown(51,3) = 0.0e0;
unknown(51,4) = 0.0e0;
unknown(52,1) = 0.0e0;
unknown(52,2) = 0.0e0;
unknown(52,3) = 0.0e0;
unknown(52,4) = 0.1e1;
unknown(53,1) = 0.1e1;
unknown(53,2) = 0.0e0;
unknown(53,3) = 0.0e0;
unknown(53,4) = -pkin(4);
unknown(54,1) = 0.0e0;
unknown(54,2) = 0.1e1;
unknown(54,3) = 0.0e0;
unknown(54,4) = 0.0e0;
unknown(55,1) = 0.0e0;
unknown(55,2) = 0.0e0;
unknown(55,3) = 0.1e1;
unknown(55,4) = 0.0e0;
unknown(56,1) = 0.0e0;
unknown(56,2) = 0.0e0;
unknown(56,3) = 0.0e0;
unknown(56,4) = 0.1e1;
unknown(57,1) = 0.0e0;
unknown(57,2) = 0.0e0;
unknown(57,3) = 0.1e1;
unknown(57,4) = -pkin(8);
unknown(58,1) = -0.1e1;
unknown(58,2) = 0.0e0;
unknown(58,3) = 0.0e0;
unknown(58,4) = 0.0e0;
unknown(59,1) = 0.0e0;
unknown(59,2) = -0.1e1;
unknown(59,3) = 0.0e0;
unknown(59,4) = 0.0e0;
unknown(60,1) = 0.0e0;
unknown(60,2) = 0.0e0;
unknown(60,3) = 0.0e0;
unknown(60,4) = 0.1e1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,15);             % numerisch
else,                         T_mdh = sym('xx', [4,4,15]); end % symbolisch

for i = 1:15
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
