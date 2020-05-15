% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh3m1DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:32:59
% EndTime: 2020-04-19 19:33:06
% DurationCPUTime: 5.81s
% Computational Cost: add. (113094->155), mult. (169525->186), div. (8100->6), fcn. (107348->48), ass. (0->137)
t102 = pkin(1) ^ 2;
t81 = sin(qJ(2));
t83 = sin(pkin(16));
t87 = cos(qJ(2));
t89 = cos(pkin(16));
t59 = t81 * t83 - t87 * t89;
t139 = pkin(5) * t59;
t114 = pkin(1) * t139;
t56 = -0.2e1 * t114;
t117 = t102 + t56;
t143 = -pkin(6) + pkin(2);
t144 = -pkin(6) - pkin(2);
t50 = sqrt(-((pkin(5) - t143) * (pkin(5) + t143) + t117) * ((pkin(5) - t144) * (pkin(5) + t144) + t117));
t60 = t81 * t89 + t87 * t83;
t134 = t50 * t60;
t100 = pkin(2) ^ 2;
t116 = pkin(5) ^ 2 + t102;
t109 = t56 + t116;
t95 = pkin(6) ^ 2;
t52 = t100 - t95 + t109;
t54 = pkin(1) - t139;
t47 = -pkin(5) * t134 + t54 * t52;
t150 = -t47 / 0.2e1;
t48 = pkin(5) * t60 * t52 + t54 * t50;
t149 = t48 / 0.2e1;
t148 = sin(pkin(17)) / 0.2e1;
t147 = sin(pkin(19)) / 0.2e1;
t146 = sin(qJ(3)) / 0.2e1;
t145 = cos(pkin(15)) / 0.2e1;
t142 = -pkin(8) - pkin(10);
t141 = -pkin(8) + pkin(10);
t101 = 0.1e1 / pkin(2);
t53 = 0.1e1 / t109;
t119 = t101 * t53;
t86 = cos(qJ(3));
t44 = (t48 * t146 + t86 * t150) * t119;
t45 = (t47 * t146 + t86 * t149) * t119;
t73 = pkin(18) + pkin(19);
t66 = sin(t73);
t67 = cos(t73);
t29 = -t67 * t44 - t66 * t45;
t140 = pkin(4) * t29;
t98 = pkin(4) ^ 2;
t120 = pkin(3) ^ 2 + t98;
t115 = pkin(3) * t140;
t25 = -0.2e1 * t115;
t111 = t25 + t120;
t22 = 0.1e1 / t111;
t92 = 0.1e1 / pkin(10);
t135 = t22 * t92;
t121 = t25 + t98;
t16 = sqrt(-((pkin(3) - t141) * (pkin(3) + t141) + t121) * ((pkin(3) - t142) * (pkin(3) + t142) + t121));
t28 = t66 * t44 - t67 * t45;
t136 = t16 * t28;
t93 = pkin(8) ^ 2;
t110 = -t93 + t120;
t91 = pkin(10) ^ 2;
t21 = t25 + t91 + t110;
t24 = -pkin(3) * t29 + pkin(4);
t14 = -pkin(3) * t136 + t24 * t21;
t15 = pkin(3) * t28 * t21 + t24 * t16;
t74 = qJ(2) + qJ(3);
t76 = cos(pkin(17));
t8 = atan2((t15 * t76 / 0.2e1 + t14 * t148) * t135, (-t14 * t76 / 0.2e1 + t15 * t148) * t135) + t74;
t6 = sin(t8);
t82 = sin(qJ(1));
t138 = t82 * t6;
t88 = cos(qJ(1));
t137 = t88 * t6;
t96 = 0.1e1 / pkin(6);
t133 = t53 * t96;
t94 = 0.1e1 / pkin(8);
t108 = t22 * t94 / 0.2e1;
t20 = -t91 + t93 + t111;
t23 = -pkin(3) + t140;
t78 = cos(pkin(19));
t39 = qJ(2) + atan2((t47 * t147 + t78 * t149) * t119, (t48 * t147 + t78 * t150) * t119);
t38 = pkin(18) - t39;
t13 = -atan2((pkin(4) * t28 * t20 - t23 * t16) * t108, (-pkin(4) * t136 - t23 * t20) * t108) + t38;
t11 = sin(t13);
t132 = t82 * t11;
t12 = cos(t13);
t131 = t82 * t12;
t107 = -t100 + t116;
t51 = t56 + t95 + t107;
t55 = pkin(1) * t59 - pkin(5);
t46 = -pkin(1) * t134 - t55 * t51;
t49 = pkin(1) * t60 * t51 - t55 * t50;
t84 = sin(pkin(15));
t43 = atan2((t49 * t145 - t46 * t84 / 0.2e1) * t133, (t46 * t145 + t49 * t84 / 0.2e1) * t133);
t40 = sin(t43);
t130 = t82 * t40;
t79 = sin(qJ(4));
t129 = t82 * t79;
t85 = cos(qJ(4));
t128 = t82 * t85;
t127 = t88 * t11;
t126 = t88 * t12;
t125 = t88 * t40;
t124 = t88 * t79;
t123 = t88 * t85;
t122 = t92 * t94;
t65 = t87 * pkin(1) + pkin(13);
t118 = t101 * t96;
t72 = pkin(12) + 0;
t69 = cos(t74);
t61 = -pkin(4) * t69 + t65;
t113 = t82 * t61 + 0;
t112 = t88 * t61 + 0;
t26 = pkin(3) * cos(t38) + t65;
t37 = pkin(19) + t39;
t64 = t81 * pkin(1) + t72;
t106 = pkin(14) + t72;
t5 = pkin(17) + t8;
t7 = cos(t8);
t105 = -pkin(9) * t7 - pkin(11) * t6;
t104 = -pkin(3) * sin(t38) + t64;
t68 = sin(t74);
t103 = -pkin(4) * t68 + t64;
t63 = t88 * t65 + 0;
t62 = t82 * t65 + 0;
t41 = cos(t43);
t36 = cos(t39);
t35 = sin(t39);
t34 = t88 * t41;
t33 = t82 * t41;
t32 = t41 * pkin(6) - pkin(7);
t27 = pkin(2) * cos(t37) + t65;
t19 = atan2(t50 * t118 / 0.2e1, -(-t107 + t95 + 0.2e1 * t114) * t118 / 0.2e1) + t37;
t18 = cos(t19);
t17 = sin(t19);
t10 = -pkin(8) * t12 + t26;
t4 = -pkin(10) * cos(t5) + t61;
t3 = atan2(t16 * t122 / 0.2e1, -(-t110 + t91 + 0.2e1 * t115) * t122 / 0.2e1) + t5;
t2 = cos(t3);
t1 = sin(t3);
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t88, -t82, 0, 0; t82, t88, 0, 0; 0, 0, 1, t72; 0, 0, 0, 1; t88 * t87, -t88 * t81, t82, t88 * pkin(13) + 0; t82 * t87, -t82 * t81, -t88, t82 * pkin(13) + 0; t81, t87, 0, t72; 0, 0, 0, 1; -t88 * t69, t88 * t68, t82, t63; -t82 * t69, t82 * t68, -t88, t62; -t68, -t69, 0, t64; 0, 0, 0, 1; -t88 * t7, t137, t82, t112; -t82 * t7, t138, -t88, t113; -t6, -t7, 0, t103; 0, 0, 0, 1; -t7 * t123 + t129, t7 * t124 + t128, -t137, t105 * t88 + t112; -t7 * t128 - t124, t7 * t129 - t123, -t138, t105 * t82 + t113; -t6 * t85, t6 * t79, t7, -t6 * pkin(9) + t7 * pkin(11) + t103; 0, 0, 0, 1; t34, -t125, t82, -t88 * pkin(7) + 0; t33, -t130, -t88, -t82 * pkin(7) + 0; t40, t41, 0, t106; 0, 0, 0, 1; t88 * t36, -t88 * t35, t82, t63; t82 * t36, -t82 * t35, -t88, t62; t35, t36, 0, t64; 0, 0, 0, 1; -t126, -t127, t82, t88 * t26 + 0; -t131, -t132, -t88, t82 * t26 + 0; t11, -t12, 0, t104; 0, 0, 0, 1; t34, -t125, t82, t88 * t32 + 0; t33, -t130, -t88, t82 * t32 + 0; t40, t41, 0, t40 * pkin(6) + t106; 0, 0, 0, 1; t126, t127, t82, t88 * t10 + 0; t131, t132, -t88, t82 * t10 + 0; -t11, t12, 0, pkin(8) * t11 + t104; 0, 0, 0, 1; -t88 * t18, t88 * t17, t82, t88 * t27 + 0; -t82 * t18, t82 * t17, -t88, t82 * t27 + 0; -t17, -t18, 0, pkin(2) * sin(t37) + t64; 0, 0, 0, 1; -t88 * t2, t88 * t1, t82, t88 * t4 + 0; -t82 * t2, t82 * t1, -t88, t82 * t4 + 0; -t1, -t2, 0, -pkin(10) * sin(t5) + t103; 0, 0, 0, 1;];
T_ges = t9;
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
