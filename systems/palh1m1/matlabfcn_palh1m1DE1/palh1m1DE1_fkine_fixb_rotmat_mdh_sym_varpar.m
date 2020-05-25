% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh1m1DE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:45:45
% EndTime: 2020-04-13 14:45:58
% DurationCPUTime: 13.84s
% Computational Cost: add. (334286->229), mult. (502459->338), div. (24096->9), fcn. (318709->55), ass. (0->198)
t211 = sin(qJ(2));
t212 = sin(pkin(19));
t214 = cos(qJ(2));
t215 = cos(pkin(19));
t174 = t211 * t215 - t214 * t212;
t171 = pkin(7) * t174;
t167 = (0.2e1 * t171 - pkin(1)) * pkin(1);
t166 = pkin(7) ^ 2 - t167;
t220 = -pkin(8) + pkin(3);
t221 = -pkin(8) - pkin(3);
t157 = sqrt(-((pkin(7) - t220) * (pkin(7) + t220) - t167) * ((pkin(7) - t221) * (pkin(7) + t221) - t167));
t175 = t211 * t212 + t214 * t215;
t156 = t157 * t175;
t122 = pkin(3) ^ 2;
t123 = pkin(8) ^ 2;
t197 = t123 - t122;
t162 = t166 - t197;
t165 = 0.1e1 / t166;
t227 = 0.1e1 / pkin(3);
t164 = t227 * t165;
t169 = -t171 + pkin(1);
t151 = (-pkin(7) * t156 + t169 * t162) * t164;
t149 = -t151 / 0.2e1;
t152 = (pkin(7) * t175 * t162 + t169 * t157) * t164;
t210 = sin(qJ(3));
t213 = cos(qJ(3));
t147 = t210 * t149 - t213 * t152 / 0.2e1;
t150 = t152 / 0.2e1;
t148 = t213 * t149 + t210 * t150;
t196 = pkin(23) + pkin(22);
t183 = sin(t196);
t184 = cos(t196);
t145 = t184 * t147 + t183 * t148;
t143 = pkin(5) * t145;
t139 = pkin(5) ^ 2 + (-0.2e1 * t143 + pkin(4)) * pkin(4);
t224 = 0.1e1 / pkin(11);
t223 = -pkin(2) - pkin(13);
t222 = -pkin(2) + pkin(13);
t219 = -pkin(9) - pkin(11);
t218 = -pkin(9) + pkin(11);
t105 = sin(pkin(20));
t108 = cos(pkin(20));
t92 = -t213 * t105 - t210 * t108;
t217 = pkin(6) * t92;
t216 = 0.1e1 / pkin(2) / 0.2e1;
t104 = sin(pkin(22));
t209 = pkin(4) * t104;
t107 = cos(pkin(22));
t208 = pkin(4) * t107;
t117 = 0.1e1 / pkin(9);
t207 = t117 * t224;
t118 = 0.1e1 / pkin(8);
t206 = t118 * t227;
t119 = pkin(6) ^ 2;
t201 = -0.2e1 * pkin(1) * t217 + t119;
t76 = sqrt(-((pkin(1) - t222) * (pkin(1) + t222) + t201) * ((pkin(1) - t223) * (pkin(1) + t223) + t201));
t93 = t210 * t105 - t213 * t108;
t205 = t93 * t76;
t110 = sin(qJ(1));
t111 = sin(pkin(18));
t163 = t165 / 0.2e1;
t161 = cos(pkin(18)) * t163;
t78 = t166 + t197;
t88 = pkin(1) * t174 - pkin(7);
t69 = -pkin(1) * t156 - t88 * t78;
t70 = pkin(1) * t175 * t78 - t88 * t157;
t61 = atan2((t69 * t111 * t163 + t70 * t161) * t118, (t69 * t161 - t111 * t70 * t165 / 0.2e1) * t118);
t58 = sin(t61);
t204 = t110 * t58;
t113 = cos(qJ(1));
t203 = t113 * t58;
t116 = 0.1e1 / pkin(13);
t189 = pkin(1) ^ 2 + t201;
t79 = 0.1e1 / t189;
t202 = t116 * t79;
t200 = cos(pkin(21));
t199 = sin(pkin(21));
t115 = pkin(13) ^ 2;
t120 = pkin(2) ^ 2;
t198 = t120 - t115;
t102 = pkin(14) + 0;
t97 = t110 * pkin(16) + 0;
t98 = t113 * pkin(16) + 0;
t195 = -pkin(1) * t92 + pkin(6);
t194 = t79 * t216;
t59 = cos(t61);
t193 = pkin(8) * t59 - pkin(15);
t192 = t110 * t211;
t191 = t113 * t211;
t190 = t116 * t216;
t77 = t189 + t198;
t86 = -pkin(1) + t217;
t66 = atan2((pkin(6) * t77 * t93 - t76 * t86) * t194, (-pkin(6) * t205 - t77 * t86) * t194);
t64 = sin(t66);
t65 = cos(t66);
t51 = -t211 * t65 - t214 * t64;
t44 = t51 * t110;
t188 = t44 * pkin(2) + t97;
t179 = t211 * t64 - t214 * t65;
t187 = -pkin(2) * t179 + t102;
t186 = -pkin(17) + t102;
t96 = t214 * pkin(1) + t102;
t46 = t51 * t113;
t185 = t46 * pkin(2) + t98;
t173 = -t214 * t210 - t211 * t213;
t182 = -pkin(5) * t173 + t96;
t181 = t189 - t198;
t90 = -pkin(1) * t192 + t97;
t91 = -pkin(1) * t191 + t98;
t103 = sin(pkin(23));
t106 = cos(pkin(23));
t60 = atan2(t106 * t150 + t103 * t151 / 0.2e1, t103 * t150 + t106 * t149);
t56 = sin(t60);
t57 = cos(t60);
t39 = -t211 * t57 - t214 * t56;
t180 = t211 * t56 - t214 * t57;
t178 = -t180 * t208 - t39 * t209 + t96;
t95 = -t211 * t210 + t214 * t213;
t82 = t95 * t110;
t177 = t82 * pkin(5) + t90;
t84 = t95 * t113;
t176 = t84 * pkin(5) + t91;
t34 = t180 * t110;
t35 = t39 * t110;
t172 = t35 * t208 - t34 * t209 + t90;
t36 = t180 * t113;
t37 = t39 * t113;
t170 = t37 * t208 - t36 * t209 + t91;
t160 = atan2((pkin(1) * t93 * t181 + t195 * t76) * t202 / 0.2e1, -(-pkin(1) * t205 + t195 * t181) * t202 / 0.2e1);
t159 = cos(t160);
t158 = sin(t160);
t155 = atan2(t157 * t206 / 0.2e1, -(t122 + t123 - t166) * t206 / 0.2e1);
t154 = cos(t155);
t153 = sin(t155);
t146 = -t183 * t147 + t184 * t148;
t144 = pkin(4) * t145;
t141 = -t144 + pkin(5);
t140 = (-0.2e1 * t144 + pkin(5)) * pkin(5);
t125 = pkin(11) ^ 2;
t138 = -t125 + t139;
t137 = 0.1e1 / t139;
t124 = pkin(9) ^ 2;
t136 = -t124 + t125 + t139;
t135 = t137 / 0.2e1;
t134 = t117 * t135;
t133 = sqrt(-((pkin(4) - t218) * (pkin(4) + t218) + t140) * ((pkin(4) - t219) * (pkin(4) + t219) + t140));
t132 = t133 * t146;
t131 = t224 * t137 * (-pkin(4) * t132 + t141 * t136);
t130 = t224 * (pkin(4) * t146 * t136 + t141 * t133) * t135;
t48 = t124 + t138;
t53 = t143 - pkin(4);
t129 = atan2((pkin(5) * t146 * t48 - t53 * t133) * t134, (-pkin(5) * t132 - t53 * t48) * t134);
t128 = sin(t129);
t127 = atan2(t199 * t131 / 0.2e1 + t200 * t130, -t200 * t131 / 0.2e1 + t199 * t130);
t126 = cos(t127);
t112 = cos(qJ(4));
t109 = sin(qJ(4));
t85 = t173 * t113;
t83 = t173 * t110;
t75 = atan2(t76 * t190, (t115 - t119 + t120 + (-pkin(1) + 0.2e1 * t217) * pkin(1)) * t190);
t74 = cos(t75);
t73 = sin(t75);
t68 = t103 * t153 - t106 * t154;
t67 = -t103 * t154 - t106 * t153;
t55 = t113 * t59;
t54 = t110 * t59;
t50 = t105 * t158 - t108 * t159;
t49 = -t105 * t159 - t108 * t158;
t45 = t179 * t113;
t43 = t179 * t110;
t31 = t179 * t74 - t51 * t73;
t30 = -t179 * t73 - t51 * t74;
t29 = -t45 * t73 - t46 * t74;
t28 = -t45 * t74 + t46 * t73;
t27 = -t43 * t73 - t44 * t74;
t26 = -t43 * t74 + t44 * t73;
t25 = atan2(t133 * t207 / 0.2e1, -(t124 - t138) * t207 / 0.2e1);
t24 = cos(t25);
t23 = sin(t25);
t22 = -t199 * t23 + t200 * t24;
t21 = t199 * t24 + t200 * t23;
t20 = cos(t129);
t18 = sin(t127);
t17 = -t104 * t128 - t107 * t20;
t16 = t104 * t20 - t107 * t128;
t12 = -t126 * t173 + t95 * t18;
t11 = -t95 * t126 - t173 * t18;
t10 = t84 * t126 + t85 * t18;
t9 = -t85 * t126 + t18 * t84;
t8 = t82 * t126 + t83 * t18;
t7 = -t83 * t126 + t18 * t82;
t6 = t16 * t39 - t17 * t180;
t5 = t16 * t180 + t17 * t39;
t4 = t16 * t36 + t17 * t37;
t3 = -t16 * t37 + t17 * t36;
t2 = t16 * t34 + t17 * t35;
t1 = -t16 * t35 + t17 * t34;
t13 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t113, -t110, 0, 0; t110, t113, 0, 0; 0, 0, 1, t102; 0, 0, 0, 1; -t191, -t113 * t214, t110, t98; -t192, -t110 * t214, -t113, t97; t214, -t211, 0, t102; 0, 0, 0, 1; t84, t85, t110, t91; t82, t83, -t113, t90; -t173, t95, 0, t96; 0, 0, 0, 1; t10, -t9, t110, t176; t8, -t7, -t113, t177; t12, -t11, 0, t182; 0, 0, 0, 1; t10 * t112 + t109 * t110, -t10 * t109 + t110 * t112, t9, t10 * pkin(10) + t9 * pkin(12) + t176; -t109 * t113 + t112 * t8, -t109 * t8 - t112 * t113, t7, t8 * pkin(10) + t7 * pkin(12) + t177; t12 * t112, -t12 * t109, t11, pkin(10) * t12 + pkin(12) * t11 + t182; 0, 0, 0, 1; t55, -t203, t110, -pkin(15) * t113 + 0; t54, -t204, -t113, -pkin(15) * t110 + 0; t58, t59, 0, t186; 0, 0, 0, 1; t37, t36, t110, t91; t35, t34, -t113, t90; -t180, t39, 0, t96; 0, 0, 0, 1; t46, t45, t110, t98; t44, t43, -t113, t97; -t179, t51, 0, t102; 0, 0, 0, 1; t29, t28, t110, t185; t27, t26, -t113, t188; t31, t30, 0, t187; 0, 0, 0, 1; t4, t3, t110, t170; t2, t1, -t113, t172; t6, t5, 0, t178; 0, 0, 0, 1; t36 * t67 + t37 * t68, t36 * t68 - t37 * t67, t110, (t103 * t36 + t106 * t37) * pkin(3) + t91; t34 * t67 + t35 * t68, t34 * t68 - t35 * t67, -t113, (t103 * t34 + t106 * t35) * pkin(3) + t90; -t180 * t68 + t39 * t67, t180 * t67 + t39 * t68, 0, (t103 * t39 - t106 * t180) * pkin(3) + t96; 0, 0, 0, 1; t49 * t85 + t50 * t84, -t49 * t84 + t50 * t85, t110, (t105 * t85 + t108 * t84) * pkin(6) + t91; t49 * t83 + t50 * t82, -t49 * t82 + t50 * t83, -t113, (t105 * t83 + t108 * t82) * pkin(6) + t90; -t173 * t50 + t49 * t95, t173 * t49 + t50 * t95, 0, (t105 * t95 - t108 * t173) * pkin(6) + t96; 0, 0, 0, 1; t10 * t22 - t21 * t9, -t10 * t21 - t22 * t9, t110, (t10 * t200 - t199 * t9) * pkin(11) + t176; -t21 * t7 + t22 * t8, -t21 * t8 - t22 * t7, -t113, (-t199 * t7 + t200 * t8) * pkin(11) + t177; -t11 * t21 + t12 * t22, -t11 * t22 - t12 * t21, 0, (-t11 * t199 + t12 * t200) * pkin(11) + t182; 0, 0, 0, 1; t55, -t203, t110, t113 * t193 + 0; t54, -t204, -t113, t110 * t193 + 0; t58, t59, 0, pkin(8) * t58 + t186; 0, 0, 0, 1; t29, t28, t110, pkin(13) * t29 + t185; t27, t26, -t113, pkin(13) * t27 + t188; t31, t30, 0, pkin(13) * t31 + t187; 0, 0, 0, 1; -t4, -t3, t110, t4 * pkin(9) + t170; -t2, -t1, -t113, t2 * pkin(9) + t172; -t6, -t5, 0, pkin(9) * t6 + t178; 0, 0, 0, 1;];
T_ges = t13;
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
