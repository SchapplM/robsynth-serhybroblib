% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_c_mdh [4x4x(16+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   11:  mdh base (link 0) -> mdh frame (11-1), link (11-1)
%   ...
%   16+1:  mdh base (link 0) -> mdh frame (16)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh1m1DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:00:28
% EndTime: 2020-04-14 20:00:34
% DurationCPUTime: 6.24s
% Computational Cost: add. (115280->207), mult. (172675->248), div. (8292->9), fcn. (109393->63), ass. (0->186)
t140 = pkin(1) ^ 2;
t114 = sin(qJ(2));
t116 = sin(pkin(19));
t120 = cos(qJ(2));
t122 = cos(pkin(19));
t86 = t114 * t122 - t120 * t116;
t192 = pkin(7) * t86;
t157 = pkin(1) * t192;
t80 = -0.2e1 * t157;
t169 = t140 + t80;
t197 = -pkin(8) + pkin(3);
t198 = -pkin(8) - pkin(3);
t68 = sqrt(-((pkin(7) - t197) * (pkin(7) + t197) + t169) * ((pkin(7) - t198) * (pkin(7) + t198) + t169));
t87 = t114 * t116 + t120 * t122;
t182 = t68 * t87;
t130 = pkin(8) ^ 2;
t136 = pkin(3) ^ 2;
t162 = pkin(7) ^ 2 + t140;
t154 = t80 + t162;
t72 = -t130 + t136 + t154;
t78 = pkin(1) - t192;
t64 = -pkin(7) * t182 + t78 * t72;
t202 = -t64 / 0.2e1;
t65 = pkin(7) * t87 * t72 + t78 * t68;
t201 = t65 / 0.2e1;
t200 = -pkin(2) - pkin(13);
t199 = -pkin(2) + pkin(13);
t196 = -pkin(9) - pkin(11);
t195 = -pkin(9) + pkin(11);
t113 = sin(qJ(3));
t137 = 0.1e1 / pkin(3);
t74 = 0.1e1 / t154;
t172 = t137 * t74;
t119 = cos(qJ(3));
t189 = -t119 / 0.2e1;
t61 = (t113 * t202 + t65 * t189) * t172;
t62 = (t113 * t201 + t64 * t189) * t172;
t103 = pkin(23) + pkin(22);
t97 = sin(t103);
t98 = cos(t103);
t29 = t98 * t61 + t97 * t62;
t194 = pkin(5) * t29;
t108 = sin(pkin(20));
t111 = cos(pkin(20));
t84 = -t119 * t108 - t113 * t111;
t193 = pkin(6) * t84;
t191 = sin(pkin(23)) / 0.2e1;
t190 = sin(pkin(21)) / 0.2e1;
t188 = cos(pkin(18)) / 0.2e1;
t187 = 0.1e1 / pkin(2) / 0.2e1;
t115 = sin(qJ(1));
t105 = qJ(2) + qJ(3);
t110 = cos(pkin(21));
t134 = pkin(5) ^ 2;
t159 = pkin(4) * t194;
t26 = -0.2e1 * t159;
t170 = t134 + t26;
t16 = sqrt(-((pkin(4) - t195) * (pkin(4) + t195) + t170) * ((pkin(4) - t196) * (pkin(4) + t196) + t170));
t30 = -t97 * t61 + t98 * t62;
t184 = t16 * t30;
t126 = pkin(11) ^ 2;
t128 = pkin(9) ^ 2;
t160 = pkin(4) ^ 2 + t134;
t147 = -t128 + t160;
t22 = t126 + t26 + t147;
t25 = -pkin(4) * t29 + pkin(5);
t14 = -pkin(4) * t184 + t25 * t22;
t15 = pkin(4) * t30 * t22 + t25 * t16;
t127 = 0.1e1 / pkin(11);
t152 = t26 + t160;
t23 = 0.1e1 / t152;
t174 = t127 * t23;
t8 = atan2((t14 * t190 + t15 * t110 / 0.2e1) * t174, (-t14 * t110 / 0.2e1 + t15 * t190) * t174) + t105;
t6 = sin(t8);
t186 = t115 * t6;
t121 = cos(qJ(1));
t185 = t121 * t6;
t133 = pkin(6) ^ 2;
t158 = pkin(1) * t193;
t77 = -0.2e1 * t158;
t171 = t133 + t77;
t67 = sqrt(-((pkin(1) - t199) * (pkin(1) + t199) + t171) * ((pkin(1) - t200) * (pkin(1) + t200) + t171));
t85 = t113 * t108 - t119 * t111;
t183 = t67 * t85;
t129 = 0.1e1 / pkin(9);
t151 = t129 * t23 / 0.2e1;
t21 = -t126 + t128 + t152;
t24 = -pkin(4) + t194;
t109 = cos(pkin(23));
t40 = qJ(2) + atan2((t109 * t201 + t64 * t191) * t172, (t109 * t202 + t65 * t191) * t172);
t34 = pkin(22) - t40;
t13 = -atan2((pkin(5) * t30 * t21 - t24 * t16) * t151, (-pkin(5) * t184 - t24 * t21) * t151) + t34;
t11 = sin(t13);
t181 = t115 * t11;
t12 = cos(t13);
t180 = t115 * t12;
t117 = sin(pkin(18));
t131 = 0.1e1 / pkin(8);
t173 = t131 * t74;
t149 = -t136 + t162;
t71 = t130 + t80 + t149;
t79 = pkin(1) * t86 - pkin(7);
t63 = -pkin(1) * t182 - t79 * t71;
t66 = pkin(1) * t87 * t71 - t79 * t68;
t44 = atan2((t66 * t188 + t63 * t117 / 0.2e1) * t173, (t63 * t188 - t117 * t66 / 0.2e1) * t173);
t41 = sin(t44);
t179 = t115 * t41;
t178 = t121 * t11;
t177 = t121 * t12;
t176 = t121 * t41;
t125 = 0.1e1 / pkin(13);
t161 = t133 + t140;
t153 = t77 + t161;
t73 = 0.1e1 / t153;
t175 = t125 * t73;
t150 = t73 * t187;
t124 = pkin(13) ^ 2;
t138 = pkin(2) ^ 2;
t69 = -t124 + t138 + t153;
t75 = -pkin(1) + t193;
t59 = qJ(2) + atan2((pkin(6) * t85 * t69 - t75 * t67) * t150, (-pkin(6) * t183 - t75 * t69) * t150);
t112 = sin(qJ(4));
t168 = t115 * t112;
t118 = cos(qJ(4));
t167 = t115 * t118;
t166 = t121 * t112;
t165 = t121 * t118;
t164 = t127 * t129;
t163 = t131 * t137;
t104 = pkin(14) + 0;
t101 = cos(t105);
t96 = -t114 * pkin(1) + pkin(16);
t88 = pkin(5) * t101 + t96;
t156 = t115 * t88 + 0;
t155 = t121 * t88 + 0;
t148 = -t138 + t161;
t57 = sin(t59);
t52 = -pkin(2) * t57 + pkin(16);
t146 = t125 * t187;
t33 = pkin(23) + t40;
t99 = pkin(20) + t105;
t58 = cos(t59);
t145 = pkin(2) * t58 + t104;
t144 = -pkin(17) + t104;
t91 = t120 * pkin(1) + t104;
t5 = pkin(21) + t8;
t7 = cos(t8);
t143 = pkin(10) * t7 + pkin(12) * t6;
t27 = pkin(4) * sin(t34) + t96;
t142 = pkin(4) * cos(t34) + t91;
t100 = sin(t105);
t141 = pkin(5) * t100 + t91;
t93 = t121 * pkin(16) + 0;
t92 = t115 * pkin(16) + 0;
t90 = t121 * t96 + 0;
t89 = t115 * t96 + 0;
t83 = pkin(6) * cos(t99) + t96;
t76 = -pkin(1) * t84 + pkin(6);
t70 = t124 + t77 + t148;
t56 = atan2((pkin(1) * t85 * t70 + t76 * t67) * t175 / 0.2e1, -(-pkin(1) * t183 + t76 * t70) * t175 / 0.2e1) + t99;
t54 = cos(t56);
t53 = sin(t56);
t51 = atan2(t67 * t146, (t124 - t148 + 0.2e1 * t158) * t146) + t59;
t50 = cos(t51);
t49 = sin(t51);
t48 = t121 * t50;
t47 = t121 * t49;
t46 = t115 * t50;
t45 = t115 * t49;
t42 = cos(t44);
t39 = t121 * t42;
t38 = t115 * t42;
t37 = t42 * pkin(8) - pkin(15);
t36 = cos(t40);
t35 = sin(t40);
t28 = -pkin(3) * sin(t33) + t96;
t20 = atan2(t68 * t163 / 0.2e1, -(t130 - t149 + 0.2e1 * t157) * t163 / 0.2e1) + t33;
t19 = cos(t20);
t18 = sin(t20);
t17 = pkin(13) * t49 + t52;
t10 = -pkin(9) * t11 + t27;
t4 = pkin(11) * cos(t5) + t88;
t3 = atan2(t16 * t164 / 0.2e1, -(t126 - t147 + 0.2e1 * t159) * t164 / 0.2e1) + t5;
t2 = cos(t3);
t1 = sin(t3);
t9 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t121, -t115, 0, 0; t115, t121, 0, 0; 0, 0, 1, t104; 0, 0, 0, 1; -t121 * t114, -t121 * t120, t115, t93; -t115 * t114, -t115 * t120, -t121, t92; t120, -t114, 0, t104; 0, 0, 0, 1; t121 * t101, -t121 * t100, t115, t90; t115 * t101, -t115 * t100, -t121, t89; t100, t101, 0, t91; 0, 0, 0, 1; t121 * t7, -t185, t115, t155; t115 * t7, -t186, -t121, t156; t6, t7, 0, t141; 0, 0, 0, 1; t7 * t165 + t168, -t7 * t166 + t167, t185, t143 * t121 + t155; t7 * t167 - t166, -t7 * t168 - t165, t186, t143 * t115 + t156; t6 * t118, -t6 * t112, -t7, t6 * pkin(10) - t7 * pkin(12) + t141; 0, 0, 0, 1; t39, -t176, t115, -t121 * pkin(15) + 0; t38, -t179, -t121, -t115 * pkin(15) + 0; t41, t42, 0, t144; 0, 0, 0, 1; -t121 * t35, -t121 * t36, t115, t90; -t115 * t35, -t115 * t36, -t121, t89; t36, -t35, 0, t91; 0, 0, 0, 1; -t121 * t57, -t121 * t58, t115, t93; -t115 * t57, -t115 * t58, -t121, t92; t58, -t57, 0, t104; 0, 0, 0, 1; t47, t48, t115, t121 * t52 + 0; t45, t46, -t121, t115 * t52 + 0; -t50, t49, 0, t145; 0, 0, 0, 1; -t178, t177, t115, t121 * t27 + 0; -t181, t180, -t121, t115 * t27 + 0; -t12, -t11, 0, t142; 0, 0, 0, 1; t121 * t18, t121 * t19, t115, t121 * t28 + 0; t115 * t18, t115 * t19, -t121, t115 * t28 + 0; -t19, t18, 0, pkin(3) * cos(t33) + t91; 0, 0, 0, 1; -t121 * t54, t121 * t53, t115, t121 * t83 + 0; -t115 * t54, t115 * t53, -t121, t115 * t83 + 0; -t53, -t54, 0, pkin(6) * sin(t99) + t91; 0, 0, 0, 1; t121 * t2, -t121 * t1, t115, t121 * t4 + 0; t115 * t2, -t115 * t1, -t121, t115 * t4 + 0; t1, t2, 0, pkin(11) * sin(t5) + t141; 0, 0, 0, 1; t39, -t176, t115, t121 * t37 + 0; t38, -t179, -t121, t115 * t37 + 0; t41, t42, 0, t41 * pkin(8) + t144; 0, 0, 0, 1; t47, t48, t115, t121 * t17 + 0; t45, t46, -t121, t115 * t17 + 0; -t50, t49, 0, -pkin(13) * t50 + t145; 0, 0, 0, 1; t178, -t177, t115, t121 * t10 + 0; t181, -t180, -t121, t115 * t10 + 0; t12, t11, 0, -pkin(9) * t12 + t142; 0, 0, 0, 1;];
T_ges = t9;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,16+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,16+1]); end % symbolisch
for i = 1:16+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
