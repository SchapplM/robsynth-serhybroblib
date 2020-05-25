% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:17:14
% EndTime: 2020-04-17 15:17:20
% DurationCPUTime: 6.35s
% Computational Cost: add. (164061->168), mult. (246709->253), div. (11772->6), fcn. (156206->24), ass. (0->142)
t166 = 0.1e1 / pkin(2);
t165 = 0.1e1 / pkin(8);
t164 = 0.1e1 / pkin(10);
t163 = t165 * t164;
t162 = -pkin(6) - pkin(2);
t161 = -pkin(6) + pkin(2);
t160 = -pkin(8) - pkin(10);
t159 = -pkin(8) + pkin(10);
t78 = sin(pkin(18));
t158 = pkin(3) * t78;
t80 = cos(pkin(18));
t157 = pkin(3) * t80;
t83 = sin(qJ(2));
t85 = sin(pkin(16));
t88 = cos(qJ(2));
t90 = cos(pkin(16));
t65 = t83 * t85 - t88 * t90;
t156 = pkin(5) * t65;
t155 = sin(pkin(15));
t93 = pkin(1) ^ 2;
t149 = pkin(5) ^ 2 + t93;
t147 = pkin(1) * t156;
t60 = -0.2e1 * t147;
t143 = t60 + t149;
t52 = 0.1e1 / t143;
t91 = 0.1e1 / pkin(6);
t154 = t52 * t91;
t150 = t60 + t93;
t49 = sqrt(-((pkin(5) - t161) * (pkin(5) + t161) + t150) * ((pkin(5) - t162) * (pkin(5) + t162) + t150));
t67 = t83 * t90 + t85 * t88;
t153 = t67 * t49;
t84 = sin(qJ(1));
t152 = t84 * t88;
t89 = cos(qJ(1));
t151 = t89 * t88;
t148 = t166 / 0.2e1;
t146 = pkin(18) + pkin(19);
t74 = pkin(12) + 0;
t145 = t84 * pkin(13) + 0;
t144 = t89 * pkin(13) + 0;
t95 = pkin(6) ^ 2;
t142 = -t95 + t149;
t141 = cos(pkin(15)) / 0.2e1;
t140 = pkin(1) - t156;
t94 = pkin(2) ^ 2;
t51 = -t94 + t95 + t143;
t59 = pkin(1) * t65 - pkin(5);
t129 = -pkin(1) * t153 - t51 * t59;
t130 = pkin(1) * t51 * t67 - t49 * t59;
t40 = (t129 * t141 + t130 * t155 / 0.2e1) * t154;
t139 = pkin(6) * t40 - pkin(7);
t62 = pkin(1) * t152 + t145;
t63 = pkin(1) * t151 + t144;
t68 = t83 * pkin(1) + t74;
t138 = pkin(14) + t74;
t137 = cos(t146);
t136 = sin(t146);
t82 = sin(qJ(3));
t87 = cos(qJ(3));
t64 = t82 * t83 - t87 * t88;
t56 = t64 * t84;
t135 = t56 * pkin(4) + t62;
t58 = t64 * t89;
t134 = t58 * pkin(4) + t63;
t131 = t82 * t88 + t83 * t87;
t133 = -pkin(4) * t131 + t68;
t132 = t60 + t94 + t142;
t124 = t166 * (-pkin(5) * t153 + t140 * t132);
t122 = -t124 / 0.2e1;
t123 = (pkin(5) * t67 * t132 + t140 * t49) * t148;
t77 = sin(pkin(19));
t79 = cos(pkin(19));
t38 = (t79 * t122 + t77 * t123) * t52;
t121 = t124 / 0.2e1;
t39 = (t77 * t121 + t79 * t123) * t52;
t33 = t38 * t88 - t39 * t83;
t32 = t38 * t83 + t39 * t88;
t128 = t32 * t157 - t33 * t158 + t68;
t27 = t33 * t84;
t28 = t32 * t84;
t127 = t27 * t157 + t28 * t158 + t62;
t29 = t33 * t89;
t30 = t32 * t89;
t126 = t29 * t157 + t30 * t158 + t63;
t125 = (-t142 + t94 + 0.2e1 * t147) * t148;
t120 = t82 * t121 + t87 * t123;
t119 = t87 * t122 + t82 * t123;
t118 = t52 * (-t137 * t119 - t136 * t120);
t117 = t52 * (t136 * t119 - t137 * t120);
t116 = pkin(4) * t118;
t115 = pkin(3) - t116;
t114 = -pkin(3) * t118 + pkin(4);
t113 = -0.2e1 * pkin(3) * t116 + pkin(4) ^ 2;
t112 = pkin(3) ^ 2 + t113;
t111 = 0.1e1 / t112;
t97 = pkin(10) ^ 2;
t110 = -t97 + t112;
t109 = t165 * t111;
t108 = t164 * t111;
t96 = pkin(8) ^ 2;
t107 = -t96 + t97 + t112;
t106 = t96 + t110;
t105 = -(t96 - t110) * t163 / 0.2e1;
t104 = sqrt(-((pkin(3) - t159) * (pkin(3) + t159) + t113) * ((pkin(3) - t160) * (pkin(3) + t160) + t113));
t103 = t104 * t163;
t102 = t104 * t117;
t101 = (-pkin(4) * t102 + t115 * t106) * t109;
t100 = (-pkin(3) * t102 + t114 * t107) * t108;
t99 = -(pkin(4) * t106 * t117 + t115 * t104) * t109 / 0.2e1;
t98 = (pkin(3) * t107 * t117 + t114 * t104) * t108 / 0.2e1;
t86 = cos(qJ(4));
t81 = sin(qJ(4));
t76 = cos(pkin(17));
t75 = sin(pkin(17));
t57 = t131 * t89;
t55 = t131 * t84;
t47 = (t77 * t49 * t148 + t79 * t125) * t91;
t46 = (t77 * t125 - t79 * t166 * t49 / 0.2e1) * t91;
t41 = (t130 * t141 - t129 * t155 / 0.2e1) * t154;
t37 = t89 * t41;
t36 = t89 * t40;
t35 = t84 * t41;
t34 = t84 * t40;
t22 = t76 * t105 - t75 * t103 / 0.2e1;
t21 = t75 * t105 + t76 * t103 / 0.2e1;
t16 = -t80 * t101 / 0.2e1 + t78 * t99;
t15 = t78 * t101 / 0.2e1 + t80 * t99;
t14 = t76 * t98 + t75 * t100 / 0.2e1;
t13 = -t76 * t100 / 0.2e1 + t75 * t98;
t12 = t13 * t64 + t131 * t14;
t11 = -t13 * t131 + t14 * t64;
t10 = t13 * t57 - t14 * t58;
t9 = t13 * t58 + t14 * t57;
t8 = t13 * t55 - t14 * t56;
t7 = t13 * t56 + t14 * t55;
t6 = -t15 * t32 + t16 * t33;
t5 = t15 * t33 + t16 * t32;
t4 = -t15 * t29 - t16 * t30;
t3 = -t15 * t30 + t16 * t29;
t2 = -t15 * t27 - t16 * t28;
t1 = -t15 * t28 + t16 * t27;
t17 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t89, -t84, 0, 0; t84, t89, 0, 0; 0, 0, 1, t74; 0, 0, 0, 1; t151, -t89 * t83, t84, t144; t152, -t84 * t83, -t89, t145; t83, t88, 0, t74; 0, 0, 0, 1; t58, t57, t84, t63; t56, t55, -t89, t62; -t131, t64, 0, t68; 0, 0, 0, 1; t9, t10, t84, t134; t7, t8, -t89, t135; t11, t12, 0, t133; 0, 0, 0, 1; t81 * t84 + t86 * t9, -t81 * t9 + t84 * t86, -t10, pkin(9) * t9 - pkin(11) * t10 + t134; t7 * t86 - t81 * t89, -t7 * t81 - t86 * t89, -t8, pkin(9) * t7 - pkin(11) * t8 + t135; t11 * t86, -t11 * t81, -t12, pkin(9) * t11 - pkin(11) * t12 + t133; 0, 0, 0, 1; t36, -t37, t84, -pkin(7) * t89 + 0; t34, -t35, -t89, -pkin(7) * t84 + 0; t41, t40, 0, t138; 0, 0, 0, 1; t29, -t30, t84, t63; t27, -t28, -t89, t62; t32, t33, 0, t68; 0, 0, 0, 1; t3, t4, t84, t126; t1, t2, -t89, t127; t5, t6, 0, t128; 0, 0, 0, 1; t36, -t37, t84, t139 * t89 + 0; t34, -t35, -t89, t139 * t84 + 0; t41, t40, 0, pkin(6) * t41 + t138; 0, 0, 0, 1; -t3, -t4, t84, pkin(8) * t3 + t126; -t1, -t2, -t89, pkin(8) * t1 + t127; -t5, -t6, 0, pkin(8) * t5 + t128; 0, 0, 0, 1; t29 * t47 - t30 * t46, -t29 * t46 - t30 * t47, t84, (t29 * t79 - t30 * t77) * pkin(2) + t63; t27 * t47 - t28 * t46, -t27 * t46 - t28 * t47, -t89, (t27 * t79 - t28 * t77) * pkin(2) + t62; t32 * t47 + t33 * t46, -t32 * t46 + t33 * t47, 0, (t32 * t79 + t33 * t77) * pkin(2) + t68; 0, 0, 0, 1; t10 * t21 + t22 * t9, t10 * t22 - t21 * t9, t84, (t10 * t75 + t76 * t9) * pkin(10) + t134; t21 * t8 + t22 * t7, -t21 * t7 + t22 * t8, -t89, (t7 * t76 + t75 * t8) * pkin(10) + t135; t11 * t22 + t12 * t21, -t11 * t21 + t12 * t22, 0, (t11 * t76 + t12 * t75) * pkin(10) + t133; 0, 0, 0, 1;];
T_ges = t17;
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
