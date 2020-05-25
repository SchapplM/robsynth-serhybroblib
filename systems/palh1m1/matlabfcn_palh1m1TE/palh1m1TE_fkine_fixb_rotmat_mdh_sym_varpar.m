% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:08:52
% EndTime: 2020-04-12 20:08:59
% DurationCPUTime: 7.32s
% Computational Cost: add. (167407->230), mult. (251842->355), div. (12048->9), fcn. (159511->28), ass. (0->196)
t119 = cos(qJ(2));
t128 = 0.1e1 / pkin(2);
t126 = pkin(6) ^ 2;
t129 = pkin(1) ^ 2;
t194 = t126 + t129;
t107 = sin(pkin(20));
t111 = cos(pkin(20));
t113 = sin(qJ(3));
t118 = cos(qJ(3));
t91 = -t107 * t118 - t111 * t113;
t213 = pkin(6) * t91;
t193 = pkin(1) * t213;
t85 = -0.2e1 * t193;
t188 = t85 + t194;
t76 = 0.1e1 / t188;
t202 = t128 * t76;
t114 = sin(qJ(2));
t210 = -t114 / 0.2e1;
t200 = t126 + t85;
t218 = -pkin(2) + pkin(13);
t219 = -pkin(2) - pkin(13);
t70 = sqrt(-((pkin(1) - t218) * (pkin(1) + t218) + t200) * ((pkin(1) - t219) * (pkin(1) + t219) + t200));
t92 = t107 * t113 - t111 * t118;
t206 = t92 * t70;
t122 = pkin(13) ^ 2;
t127 = pkin(2) ^ 2;
t74 = -t122 + t127 + t188;
t84 = -pkin(1) + t213;
t65 = -pkin(6) * t206 - t74 * t84;
t66 = pkin(6) * t74 * t92 - t70 * t84;
t227 = (t119 * t65 / 0.2e1 + t66 * t210) * t202;
t226 = 0.1e1 / pkin(3);
t225 = 0.1e1 / pkin(9);
t224 = 0.1e1 / pkin(11);
t223 = t225 * t224;
t222 = -t70 / 0.2e1;
t221 = t70 / 0.2e1;
t184 = -t127 + t194;
t220 = -t122 / 0.2e1 + t184 / 0.2e1 - t193;
t217 = -pkin(8) - pkin(3);
t216 = -pkin(8) + pkin(3);
t215 = -pkin(9) - pkin(11);
t214 = -pkin(9) + pkin(11);
t116 = sin(pkin(19));
t121 = cos(pkin(19));
t93 = t114 * t121 - t116 * t119;
t212 = pkin(7) * t93;
t211 = t107 / 0.2e1;
t209 = sin(pkin(18));
t105 = sin(pkin(22));
t208 = pkin(4) * t105;
t109 = cos(pkin(22));
t207 = pkin(4) * t109;
t192 = pkin(1) * t212;
t87 = -0.2e1 * t192;
t199 = t129 + t87;
t71 = sqrt(-((pkin(7) - t216) * (pkin(7) + t216) + t199) * ((pkin(7) - t217) * (pkin(7) + t217) + t199));
t96 = t114 * t116 + t119 * t121;
t205 = t96 * t71;
t123 = 0.1e1 / pkin(13);
t204 = t123 * t76;
t124 = 0.1e1 / pkin(8);
t195 = pkin(7) ^ 2 + t129;
t189 = t87 + t195;
t77 = 0.1e1 / t189;
t203 = t124 * t77;
t201 = t226 / 0.2e1;
t115 = sin(qJ(1));
t198 = t115 * t114;
t120 = cos(qJ(1));
t197 = t120 * t114;
t196 = t123 * t128;
t103 = pkin(14) + 0;
t191 = pkin(23) + pkin(22);
t98 = t115 * pkin(16) + 0;
t99 = t120 * pkin(16) + 0;
t190 = cos(pkin(18)) / 0.2e1;
t187 = -pkin(1) * t91 + pkin(6);
t186 = pkin(1) - t212;
t131 = pkin(8) ^ 2;
t185 = -t131 + t195;
t130 = pkin(3) ^ 2;
t75 = -t130 + t131 + t189;
t86 = pkin(1) * t93 - pkin(7);
t171 = -pkin(1) * t205 - t75 * t86;
t172 = pkin(1) * t75 * t96 - t71 * t86;
t57 = (t171 * t190 - t209 * t172 / 0.2e1) * t203;
t183 = pkin(8) * t57 - pkin(15);
t182 = pkin(2) * t227 + t103;
t181 = -pkin(17) + t103;
t97 = t119 * pkin(1) + t103;
t50 = (-t119 * t66 / 0.2e1 + t65 * t210) * t202;
t42 = t115 * t50;
t180 = t42 * pkin(2) + t98;
t44 = t120 * t50;
t179 = t44 * pkin(2) + t99;
t178 = cos(t191);
t177 = sin(t191);
t94 = t113 * t119 + t114 * t118;
t176 = t94 * pkin(5) + t97;
t175 = t122 + t85 + t184;
t174 = t130 + t87 + t185;
t104 = sin(pkin(23));
t108 = cos(pkin(23));
t159 = t226 * (-pkin(7) * t205 + t186 * t174);
t157 = -t159 / 0.2e1;
t160 = t226 * (pkin(7) * t96 * t174 + t186 * t71);
t158 = t160 / 0.2e1;
t55 = (t104 * t158 + t108 * t157) * t77;
t56 = (t108 * t158 + t104 * t159 / 0.2e1) * t77;
t173 = t114 * t56 - t119 * t55;
t39 = -t114 * t55 - t119 * t56;
t95 = -t113 * t114 + t118 * t119;
t89 = -pkin(1) * t198 + t98;
t90 = -pkin(1) * t197 + t99;
t80 = t95 * t115;
t169 = t80 * pkin(5) + t89;
t82 = t95 * t120;
t168 = t82 * pkin(5) + t90;
t167 = -t173 * t207 - t39 * t208 + t97;
t165 = (t130 - t185 + 0.2e1 * t192) * t201;
t33 = t39 * t115;
t34 = t173 * t115;
t164 = t33 * t207 - t34 * t208 + t89;
t35 = t39 * t120;
t36 = t173 * t120;
t163 = t35 * t207 - t36 * t208 + t90;
t162 = pkin(1) * t92 * t175 + t187 * t70;
t161 = -pkin(1) * t206 + t187 * t175;
t156 = t113 * t157 - t118 * t160 / 0.2e1;
t155 = t113 * t158 + t118 * t157;
t154 = t77 * (t177 * t155 + t178 * t156);
t153 = t77 * (t178 * t155 - t177 * t156);
t152 = pkin(5) * t154;
t151 = pkin(4) - t152;
t150 = -pkin(4) * t154 + pkin(5);
t149 = -0.2e1 * pkin(4) * t152 + pkin(5) ^ 2;
t148 = pkin(4) ^ 2 + t149;
t147 = 0.1e1 / t148;
t133 = pkin(11) ^ 2;
t146 = -t133 + t148;
t145 = t225 * t147;
t144 = t224 * t147;
t132 = pkin(9) ^ 2;
t143 = -t132 + t133 + t148;
t142 = t132 + t146;
t141 = -(t132 - t146) * t223 / 0.2e1;
t140 = sqrt(-((pkin(4) - t214) * (pkin(4) + t214) + t149) * ((pkin(4) - t215) * (pkin(4) + t215) + t149));
t139 = t140 * t223;
t138 = t140 * t153;
t137 = (-pkin(5) * t138 + t151 * t142) * t145;
t136 = (-pkin(4) * t138 + t150 * t143) * t144;
t135 = -(pkin(5) * t142 * t153 + t151 * t140) * t145 / 0.2e1;
t134 = (pkin(4) * t143 * t153 + t150 * t140) * t144 / 0.2e1;
t117 = cos(qJ(4));
t112 = sin(qJ(4));
t110 = cos(pkin(21));
t106 = sin(pkin(21));
t83 = t94 * t120;
t81 = t94 * t115;
t68 = (t104 * t71 * t201 + t108 * t165) * t124;
t67 = (t104 * t165 - t108 * t226 * t71 / 0.2e1) * t124;
t58 = (t172 * t190 + t171 * t209 / 0.2e1) * t203;
t54 = t120 * t58;
t53 = t120 * t57;
t52 = t115 * t58;
t51 = t115 * t57;
t48 = (t111 * t161 / 0.2e1 + t162 * t211) * t204;
t47 = (t161 * t211 - t111 * t162 / 0.2e1) * t204;
t45 = t120 * t227;
t43 = t115 * t227;
t29 = (t50 * t220 + t221 * t227) * t196;
t28 = (t220 * t227 + t50 * t222) * t196;
t27 = (-t45 * t220 + t44 * t221) * t196;
t26 = (t44 * t220 - t45 * t222) * t196;
t25 = (-t43 * t220 + t42 * t221) * t196;
t24 = (t42 * t220 - t43 * t222) * t196;
t22 = t110 * t141 - t106 * t139 / 0.2e1;
t21 = t106 * t141 + t110 * t139 / 0.2e1;
t16 = -t110 * t136 / 0.2e1 + t106 * t134;
t15 = t106 * t136 / 0.2e1 + t110 * t134;
t14 = -t109 * t137 / 0.2e1 + t105 * t135;
t13 = t105 * t137 / 0.2e1 + t109 * t135;
t12 = t15 * t95 + t16 * t94;
t11 = -t15 * t94 + t16 * t95;
t10 = -t15 * t83 + t16 * t82;
t9 = -t15 * t82 - t16 * t83;
t8 = -t15 * t81 + t16 * t80;
t7 = -t15 * t80 - t16 * t81;
t6 = t13 * t173 + t14 * t39;
t5 = t13 * t39 - t14 * t173;
t4 = -t13 * t35 + t14 * t36;
t3 = t13 * t36 + t14 * t35;
t2 = -t13 * t33 + t14 * t34;
t1 = t13 * t34 + t14 * t33;
t17 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t120, -t115, 0, 0; t115, t120, 0, 0; 0, 0, 1, t103; 0, 0, 0, 1; -t197, -t120 * t119, t115, t99; -t198, -t115 * t119, -t120, t98; t119, -t114, 0, t103; 0, 0, 0, 1; t82, -t83, t115, t90; t80, -t81, -t120, t89; t94, t95, 0, t97; 0, 0, 0, 1; t10, t9, t115, t168; t8, t7, -t120, t169; t12, t11, 0, t176; 0, 0, 0, 1; t10 * t117 + t112 * t115, -t10 * t112 + t115 * t117, -t9, pkin(10) * t10 - pkin(12) * t9 + t168; -t112 * t120 + t117 * t8, -t112 * t8 - t117 * t120, -t7, pkin(10) * t8 - pkin(12) * t7 + t169; t12 * t117, -t12 * t112, -t11, pkin(10) * t12 - pkin(12) * t11 + t176; 0, 0, 0, 1; t53, -t54, t115, -pkin(15) * t120 + 0; t51, -t52, -t120, -pkin(15) * t115 + 0; t58, t57, 0, t181; 0, 0, 0, 1; t35, t36, t115, t90; t33, t34, -t120, t89; -t173, t39, 0, t97; 0, 0, 0, 1; t44, -t45, t115, t99; t42, -t43, -t120, t98; t227, t50, 0, t103; 0, 0, 0, 1; t26, t27, t115, t179; t24, t25, -t120, t180; t28, t29, 0, t182; 0, 0, 0, 1; t3, t4, t115, t163; t1, t2, -t120, t164; t5, t6, 0, t167; 0, 0, 0, 1; t35 * t68 + t36 * t67, -t35 * t67 + t36 * t68, t115, (t104 * t36 + t108 * t35) * pkin(3) + t90; t33 * t68 + t34 * t67, -t33 * t67 + t34 * t68, -t120, (t104 * t34 + t108 * t33) * pkin(3) + t89; -t173 * t68 + t39 * t67, t173 * t67 + t39 * t68, 0, (t104 * t39 - t108 * t173) * pkin(3) + t97; 0, 0, 0, 1; -t47 * t83 + t48 * t82, -t47 * t82 - t48 * t83, t115, (-t107 * t83 + t111 * t82) * pkin(6) + t90; -t47 * t81 + t48 * t80, -t47 * t80 - t48 * t81, -t120, (-t107 * t81 + t111 * t80) * pkin(6) + t89; t47 * t95 + t48 * t94, -t47 * t94 + t48 * t95, 0, (t107 * t95 + t111 * t94) * pkin(6) + t97; 0, 0, 0, 1; t10 * t22 + t21 * t9, -t10 * t21 + t22 * t9, t115, (t10 * t110 + t106 * t9) * pkin(11) + t168; t21 * t7 + t22 * t8, -t21 * t8 + t22 * t7, -t120, (t106 * t7 + t110 * t8) * pkin(11) + t169; t11 * t21 + t12 * t22, t11 * t22 - t12 * t21, 0, (t106 * t11 + t110 * t12) * pkin(11) + t176; 0, 0, 0, 1; t53, -t54, t115, t120 * t183 + 0; t51, -t52, -t120, t115 * t183 + 0; t58, t57, 0, pkin(8) * t58 + t181; 0, 0, 0, 1; t26, t27, t115, pkin(13) * t26 + t179; t24, t25, -t120, pkin(13) * t24 + t180; t28, t29, 0, pkin(13) * t28 + t182; 0, 0, 0, 1; -t3, -t4, t115, pkin(9) * t3 + t163; -t1, -t2, -t120, pkin(9) * t1 + t164; -t5, -t6, 0, pkin(9) * t5 + t167; 0, 0, 0, 1;];
T_ges = t17;
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
