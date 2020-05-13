% Calculate homogenous joint transformation matrices for
% mg10hlDE2
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
% Datum: 2020-04-11 13:01
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = mg10hlDE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [17x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:00:02
% EndTime: 2020-04-11 13:00:03
% DurationCPUTime: 0.51s
% Computational Cost: add. (2507->95), mult. (2778->133), div. (368->17), fcn. (1240->35), ass. (0->346)
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
t26 = (-t16 + t17 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4))) * (-t16 + t17 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4)));
t27 = sqrt(-t26);
t30 = -t13 * pkin(2) + pkin(1);
t31 = pkin(2) ^ 2;
t32 = pkin(5) ^ 2;
t33 = pkin(4) ^ 2;
t34 = -t16 + t17 + t31 - t32 + t33;
t36 = -t27 * t10 + t34 * t30;
t38 = 0.1e1 / pkin(4);
t39 = -t16 + t17 + t31;
t40 = 0.1e1 / t39;
t44 = t34 * t10 + t27 * t30;
t45 = t44 ^ 2;
t46 = 0.1e1 / t33;
t48 = t39 ^ 2;
t49 = 0.1e1 / t48;
t51 = t36 ^ 2;
t55 = sqrt(t49 * t46 * t45 + t49 * t46 * t51);
t57 = 0.1e1 / t55 * t40 * t38;
t61 = -t57 * t36 * t5 + t57 * t44 * t7;
t66 = t57 * t36 * t7 + t57 * t44 * t5;
t69 = cos(pkin(17));
t70 = 1 / pkin(7);
t71 = qJ(6) + pkin(8);
t73 = 0.1e1 / t71 * t70;
t74 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
t75 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
t77 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
t78 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
t81 = sqrt(-t78 * t77 * t75 * t74);
t83 = (pkin(6) ^ 2);
t84 = (pkin(7) ^ 2);
t85 = (pkin(8) ^ 2);
t87 = 2 * pkin(8) * qJ(6);
t88 = qJ(6) ^ 2;
t89 = t83 - t84 - t85 - t87 - t88;
t91 = atan2(t81 * t73, t89 * t73);
t92 = t91 + pkin(16);
t93 = cos(t92);
t95 = sin(pkin(17));
t96 = sin(t92);
t98 = -t93 * t69 - t96 * t95;
t101 = -t96 * t69 + t93 * t95;
t105 = t71 ^ 2;
t106 = 1 / t105;
t107 = t106 / pkin(6);
t109 = t70 * (t83 - t84 + t85 + t87 + t88);
t113 = t77 * t75;
t117 = -t70 * t78 * t113 * t74 * t107 + t89 * t109 * t107;
t118 = cos(pkin(16));
t125 = -t89 * t70 * t81 * t107 + t81 * t109 * t107;
t126 = sin(pkin(16));
t128 = -t118 * t117 / 0.4e1 + t126 * t125 / 0.4e1;
t131 = t118 * t125 / 0.4e1 + t126 * t117 / 0.4e1;
t132 = t131 ^ 2;
t133 = t128 ^ 2;
t135 = sqrt(t132 + t133);
t136 = 0.1e1 / t135;
t137 = t136 * t128;
t138 = t136 * t131;
t139 = cos(qJ(3));
t140 = sin(qJ(3));
t141 = sin(qJ(4));
t142 = cos(qJ(4));
t143 = cos(qJ(5));
t144 = sin(qJ(5));
t146 = t106 / t84;
t150 = t89 ^ 2;
t153 = sqrt(-t78 * t113 * t74 * t146 + (t150 * t146));
t154 = 0.1e1 / t153;
t156 = t154 * t89 * t73;
t158 = t154 * t81 * t73;
t161 = t14 - pkin(2);
t162 = -t16 + t17 + t31 + t32 - t33;
t164 = -pkin(1) * t27 * t9 - t162 * t161;
t165 = 0.1e1 / pkin(5);
t170 = t162 * t9 * pkin(1) - t27 * t161;
t171 = t170 ^ 2;
t172 = 0.1e1 / t32;
t175 = t164 ^ 2;
t179 = sqrt(t49 * t172 * t171 + t49 * t172 * t175);
t181 = 0.1e1 / t179 * t40;
t182 = t181 * t165 * t164;
t184 = t181 * t165 * t170;
t185 = t16 - t17 - t31 + t32 + t33;
t189 = t185 ^ 2;
t193 = sqrt(t46 * t172 * t189 - t26 * t46 * t172);
t194 = 0.1e1 / t193;
t196 = t194 * t38 * t165 * t185;
t199 = t194 * t27 * t38 * t165;
t203 = sqrt(t125 ^ 2 / 0.16e2 + t117 ^ 2 / 0.16e2);
t204 = 0.1e1 / t203;
t205 = t204 * t117 / 0.4e1;
t206 = t204 * t125 / 0.4e1;
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
unknown(9,1) = t61;
unknown(9,2) = t66;
unknown(9,3) = 0.0e0;
unknown(9,4) = pkin(1) * t5;
unknown(10,1) = -t66;
unknown(10,2) = t61;
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
unknown(13,1) = t98;
unknown(13,2) = t101;
unknown(13,3) = 0.0e0;
unknown(13,4) = pkin(3) * t69;
unknown(14,1) = -t101;
unknown(14,2) = t98;
unknown(14,3) = 0.0e0;
unknown(14,4) = pkin(3) * t95;
unknown(15,1) = 0.0e0;
unknown(15,2) = 0.0e0;
unknown(15,3) = 0.1e1;
unknown(15,4) = 0.0e0;
unknown(16,1) = 0.0e0;
unknown(16,2) = 0.0e0;
unknown(16,3) = 0.0e0;
unknown(16,4) = 0.1e1;
unknown(17,1) = -t137;
unknown(17,2) = t138;
unknown(17,3) = 0.0e0;
unknown(17,4) = pkin(7);
unknown(18,1) = -t138;
unknown(18,2) = -t137;
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
unknown(22,1) = t139;
unknown(22,2) = -t140;
unknown(22,3) = 0.0e0;
unknown(22,4) = 0.0e0;
unknown(23,1) = t140;
unknown(23,2) = t139;
unknown(23,3) = 0.0e0;
unknown(23,4) = 0.0e0;
unknown(24,1) = 0.0e0;
unknown(24,2) = 0.0e0;
unknown(24,3) = 0.0e0;
unknown(24,4) = 0.1e1;
unknown(25,1) = t141;
unknown(25,2) = t142;
unknown(25,3) = 0.0e0;
unknown(25,4) = 0.0e0;
unknown(26,1) = 0.0e0;
unknown(26,2) = 0.0e0;
unknown(26,3) = 0.1e1;
unknown(26,4) = 0.0e0;
unknown(27,1) = t142;
unknown(27,2) = -t141;
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
unknown(30,1) = t143;
unknown(30,2) = -t144;
unknown(30,3) = 0.0e0;
unknown(30,4) = 0.0e0;
unknown(31,1) = t144;
unknown(31,2) = t143;
unknown(31,3) = 0.0e0;
unknown(31,4) = 0.0e0;
unknown(32,1) = 0.0e0;
unknown(32,2) = 0.0e0;
unknown(32,3) = 0.0e0;
unknown(32,4) = 0.1e1;
unknown(33,1) = t156;
unknown(33,2) = -t158;
unknown(33,3) = 0.0e0;
unknown(33,4) = 0.0e0;
unknown(34,1) = t158;
unknown(34,2) = t156;
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
unknown(37,1) = t182;
unknown(37,2) = -t184;
unknown(37,3) = 0.0e0;
unknown(37,4) = -pkin(14);
unknown(38,1) = 0.0e0;
unknown(38,2) = 0.0e0;
unknown(38,3) = -0.1e1;
unknown(38,4) = 0.0e0;
unknown(39,1) = t184;
unknown(39,2) = t182;
unknown(39,3) = 0.0e0;
unknown(39,4) = 0.0e0;
unknown(40,1) = 0.0e0;
unknown(40,2) = 0.0e0;
unknown(40,3) = 0.0e0;
unknown(40,4) = 0.1e1;
unknown(41,1) = -t196;
unknown(41,2) = t199;
unknown(41,3) = 0.0e0;
unknown(41,4) = pkin(6);
unknown(42,1) = -t199;
unknown(42,2) = -t196;
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
unknown(45,1) = -t205;
unknown(45,2) = t206;
unknown(45,3) = 0.0e0;
unknown(45,4) = pkin(5);
unknown(46,1) = -t206;
unknown(46,2) = -t205;
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
