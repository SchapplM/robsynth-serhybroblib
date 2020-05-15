% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_c_mdh [4x4x(12+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   9:  mdh base (link 0) -> mdh frame (9-1), link (9-1)
%   ...
%   12+1:  mdh base (link 0) -> mdh frame (12)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m1DE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:23:32
% EndTime: 2020-04-18 10:23:45
% DurationCPUTime: 13.11s
% Computational Cost: add. (327724->168), mult. (492493->246), div. (23544->6), fcn. (312110->42), ass. (0->147)
t152 = sin(qJ(3));
t153 = sin(qJ(2));
t155 = cos(qJ(3));
t156 = cos(qJ(2));
t130 = t156 * t152 + t153 * t155;
t147 = sin(pkin(17));
t148 = cos(pkin(17));
t154 = sin(pkin(16));
t157 = cos(pkin(16));
t129 = t153 * t154 - t156 * t157;
t127 = pkin(5) * t129;
t122 = (0.2e1 * t127 - pkin(1)) * pkin(1);
t164 = -pkin(6) + pkin(2);
t165 = -pkin(6) - pkin(2);
t116 = sqrt(-((pkin(5) - t164) * (pkin(5) + t164) - t122) * ((pkin(5) - t165) * (pkin(5) + t165) - t122));
t131 = t153 * t157 + t156 * t154;
t115 = t116 * t131;
t121 = pkin(5) ^ 2 - t122;
t80 = pkin(2) ^ 2;
t81 = pkin(6) ^ 2;
t149 = t81 - t80;
t118 = t121 - t149;
t120 = 0.1e1 / t121;
t124 = -t127 + pkin(1);
t169 = 0.1e1 / pkin(2);
t111 = t169 * t120 * (-pkin(5) * t115 + t118 * t124);
t109 = -t111 / 0.2e1;
t119 = t120 / 0.2e1;
t110 = t169 * (pkin(5) * t118 * t131 + t116 * t124) * t119;
t106 = t155 * t109 + t152 * t110;
t108 = t111 / 0.2e1;
t107 = t152 * t108 + t155 * t110;
t146 = pkin(18) + pkin(19);
t138 = sin(t146);
t139 = cos(t146);
t105 = -t106 * t139 - t107 * t138;
t103 = pkin(3) * t105;
t100 = -t103 + pkin(4);
t104 = t106 * t138 - t107 * t139;
t166 = 0.1e1 / pkin(10);
t162 = -pkin(8) + pkin(10);
t163 = -pkin(8) - pkin(10);
t99 = (-0.2e1 * t103 + pkin(4)) * pkin(4);
t92 = sqrt(-((pkin(3) - t162) * (pkin(3) + t162) + t99) * ((pkin(3) - t163) * (pkin(3) + t163) + t99));
t102 = pkin(4) * t105;
t98 = pkin(4) ^ 2 + (-0.2e1 * t102 + pkin(3)) * pkin(3);
t96 = 0.1e1 / t98;
t94 = t96 / 0.2e1;
t82 = pkin(8) ^ 2;
t83 = pkin(10) ^ 2;
t95 = -t82 + t83 + t98;
t89 = t166 * (pkin(3) * t104 * t95 + t100 * t92) * t94;
t91 = t92 * t104;
t90 = t166 * t96 * (-pkin(3) * t91 + t100 * t95);
t86 = atan2(t148 * t89 + t147 * t90 / 0.2e1, -t148 * t90 / 0.2e1 + t147 * t89);
t85 = cos(t86);
t84 = t130 * t85;
t69 = sin(pkin(18));
t161 = pkin(3) * t69;
t71 = cos(pkin(18));
t160 = pkin(3) * t71;
t78 = 0.1e1 / pkin(8);
t159 = t78 * t166;
t79 = 0.1e1 / pkin(6);
t158 = t79 * t169;
t117 = cos(pkin(15)) * t119;
t51 = t121 + t149;
t56 = pkin(1) * t129 - pkin(5);
t47 = -pkin(1) * t115 - t56 * t51;
t48 = pkin(1) * t131 * t51 - t116 * t56;
t74 = sin(pkin(15));
t44 = atan2((t48 * t117 - t47 * t120 * t74 / 0.2e1) * t79, (t119 * t48 * t74 + t117 * t47) * t79);
t41 = sin(t44);
t73 = sin(qJ(1));
t151 = t73 * t41;
t76 = cos(qJ(1));
t150 = t76 * t41;
t67 = pkin(12) + 0;
t145 = t73 * pkin(13) + 0;
t144 = t76 * pkin(13) + 0;
t143 = t73 * t156;
t142 = t76 * t156;
t42 = cos(t44);
t141 = pkin(6) * t42 - pkin(7);
t58 = pkin(1) * t143 + t145;
t59 = pkin(1) * t142 + t144;
t61 = t153 * pkin(1) + t67;
t140 = pkin(14) + t67;
t128 = t153 * t152 - t156 * t155;
t54 = t128 * t73;
t137 = t54 * pkin(4) + t58;
t55 = t128 * t76;
t136 = t55 * pkin(4) + t59;
t135 = -pkin(4) * t130 + t61;
t68 = sin(pkin(19));
t70 = cos(pkin(19));
t43 = atan2(t108 * t68 + t110 * t70, t109 * t70 + t110 * t68);
t39 = sin(t43);
t40 = cos(t43);
t33 = -t153 * t39 + t156 * t40;
t34 = t153 * t40 + t156 * t39;
t134 = t34 * t160 - t33 * t161 + t61;
t28 = t34 * t73;
t29 = t33 * t73;
t133 = t29 * t160 + t28 * t161 + t58;
t30 = t34 * t76;
t31 = t33 * t76;
t132 = t31 * t160 + t30 * t161 + t59;
t126 = t130 * t73;
t125 = t130 * t76;
t114 = atan2(t116 * t158 / 0.2e1, -(t80 + t81 - t121) * t158 / 0.2e1);
t113 = cos(t114);
t112 = sin(t114);
t97 = -t83 + t98;
t93 = t78 * t94;
t35 = t82 + t97;
t36 = t102 - pkin(3);
t88 = atan2((pkin(4) * t104 * t35 - t36 * t92) * t93, (-pkin(4) * t91 - t36 * t35) * t93);
t87 = sin(t88);
t75 = cos(qJ(4));
t72 = sin(qJ(4));
t46 = t112 * t68 - t113 * t70;
t45 = -t112 * t70 - t113 * t68;
t38 = t76 * t42;
t37 = t73 * t42;
t25 = atan2(t92 * t159 / 0.2e1, -(t82 - t97) * t159 / 0.2e1);
t24 = cos(t25);
t23 = sin(t25);
t22 = -t147 * t23 + t148 * t24;
t21 = t147 * t24 + t148 * t23;
t20 = cos(t88);
t18 = sin(t86);
t17 = -t71 * t20 - t69 * t87;
t16 = t20 * t69 - t71 * t87;
t12 = t128 * t18 - t84;
t11 = -t128 * t85 - t130 * t18;
t10 = t125 * t18 + t55 * t85;
t9 = t18 * t55 - t76 * t84;
t8 = t126 * t18 + t54 * t85;
t7 = t18 * t54 - t73 * t84;
t6 = t16 * t33 + t17 * t34;
t5 = -t16 * t34 + t17 * t33;
t4 = -t16 * t30 + t17 * t31;
t3 = -t16 * t31 - t17 * t30;
t2 = -t16 * t28 + t17 * t29;
t1 = -t16 * t29 - t17 * t28;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t76, -t73, 0, 0; t73, t76, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t142, -t76 * t153, t73, t144; t143, -t73 * t153, -t76, t145; t153, t156, 0, t67; 0, 0, 0, 1; t55, t125, t73, t59; t54, t126, -t76, t58; -t130, t128, 0, t61; 0, 0, 0, 1; t10, -t9, t73, t136; t8, -t7, -t76, t137; t12, -t11, 0, t135; 0, 0, 0, 1; t10 * t75 + t72 * t73, -t10 * t72 + t73 * t75, t9, pkin(9) * t10 + pkin(11) * t9 + t136; -t72 * t76 + t75 * t8, -t72 * t8 - t75 * t76, t7, pkin(9) * t8 + pkin(11) * t7 + t137; t12 * t75, -t12 * t72, t11, pkin(9) * t12 + pkin(11) * t11 + t135; 0, 0, 0, 1; t38, -t150, t73, -pkin(7) * t76 + 0; t37, -t151, -t76, -pkin(7) * t73 + 0; t41, t42, 0, t140; 0, 0, 0, 1; t31, -t30, t73, t59; t29, -t28, -t76, t58; t34, t33, 0, t61; 0, 0, 0, 1; t4, t3, t73, t132; t2, t1, -t76, t133; t6, t5, 0, t134; 0, 0, 0, 1; t38, -t150, t73, t141 * t76 + 0; t37, -t151, -t76, t141 * t73 + 0; t41, t42, 0, pkin(6) * t41 + t140; 0, 0, 0, 1; -t4, -t3, t73, pkin(8) * t4 + t132; -t2, -t1, -t76, pkin(8) * t2 + t133; -t6, -t5, 0, pkin(8) * t6 + t134; 0, 0, 0, 1; -t30 * t45 + t31 * t46, -t30 * t46 - t31 * t45, t73, (-t30 * t68 + t31 * t70) * pkin(2) + t59; -t28 * t45 + t29 * t46, -t28 * t46 - t29 * t45, -t76, (-t28 * t68 + t29 * t70) * pkin(2) + t58; t33 * t45 + t34 * t46, t33 * t46 - t34 * t45, 0, (t33 * t68 + t34 * t70) * pkin(2) + t61; 0, 0, 0, 1; t10 * t22 - t21 * t9, -t10 * t21 - t22 * t9, t73, (t148 * t10 - t147 * t9) * pkin(10) + t136; -t21 * t7 + t22 * t8, -t21 * t8 - t22 * t7, -t76, (-t147 * t7 + t148 * t8) * pkin(10) + t137; -t11 * t21 + t12 * t22, -t11 * t22 - t12 * t21, 0, (-t147 * t11 + t148 * t12) * pkin(10) + t135; 0, 0, 0, 1;];
T_ges = t13;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,12+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,12+1]); end % symbolisch
for i = 1:12+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
