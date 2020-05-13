% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh3m1DE1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_inertiaJ_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t287 = 2 * pkin(3);
t143 = sin(qJ(3));
t144 = sin(qJ(2));
t148 = cos(qJ(3));
t149 = cos(qJ(2));
t123 = t143 * t144 - t148 * t149;
t125 = -t143 * t149 - t144 * t148;
t136 = sin(pkin(17));
t158 = pkin(4) ^ 2;
t157 = pkin(5) ^ 2;
t162 = pkin(1) ^ 2;
t145 = sin(pkin(16));
t249 = cos(pkin(16));
t124 = t144 * t145 - t149 * t249;
t248 = pkin(5) * t124;
t281 = -2 * pkin(1);
t215 = t248 * t281 + t162;
t118 = t157 + t215;
t116 = 0.1e1 / t118;
t161 = 0.1e1 / pkin(2);
t219 = t116 * t161;
t212 = pkin(2) ^ 2 - pkin(6) ^ 2;
t114 = t118 + t212;
t119 = pkin(1) - t248;
t268 = -pkin(6) - pkin(2);
t111 = (pkin(5) - t268) * (pkin(5) + t268) + t215;
t267 = -pkin(6) + pkin(2);
t112 = (pkin(5) - t267) * (pkin(5) + t267) + t215;
t163 = sqrt(-t112 * t111);
t126 = t144 * t249 + t149 * t145;
t247 = pkin(5) * t126;
t106 = t114 * t247 + t119 * t163;
t223 = t106 * t148;
t217 = t126 * t163;
t105 = -pkin(5) * t217 + t114 * t119;
t228 = t105 * t143;
t100 = (t223 / 0.2e1 + t228 / 0.2e1) * t219;
t133 = pkin(18) + pkin(19);
t131 = sin(t133);
t132 = cos(t133);
t224 = t106 * t143;
t227 = t105 * t148;
t99 = (-t227 / 0.2e1 + t224 / 0.2e1) * t219;
t77 = t100 * t131 + t132 * t99;
t256 = t77 * pkin(4);
t231 = t256 * t287 + t158;
t266 = (-pkin(8) - pkin(10));
t66 = ((pkin(3) - t266) * (pkin(3) + t266)) + t231;
t265 = (pkin(10) - pkin(8));
t67 = ((pkin(3) - t265) * (pkin(3) + t265)) + t231;
t164 = sqrt(-t67 * t66);
t177 = -t100 * t132 + t131 * t99;
t213 = pkin(8) ^ 2 - pkin(10) ^ 2;
t159 = pkin(3) ^ 2;
t72 = t159 + t231;
t69 = t72 - t213;
t284 = t177 * t69;
t74 = pkin(3) * t77 + pkin(4);
t52 = pkin(3) * t284 + t164 * t74;
t233 = t52 * t136;
t137 = cos(pkin(17));
t283 = t164 * t177;
t50 = -pkin(3) * t283 + t69 * t74;
t234 = t50 * t137;
t152 = 0.1e1 / pkin(10);
t70 = 0.1e1 / t72;
t238 = t152 * t70;
t37 = (-t234 / 0.2e1 + t233 / 0.2e1) * t238;
t232 = t52 * t137;
t254 = t136 / 0.2e1;
t38 = (t232 / 0.2e1 + t50 * t254) * t238;
t34 = atan2(t38, t37);
t31 = sin(t34);
t32 = cos(t34);
t24 = t123 * t31 + t125 * t32;
t286 = -0.2e1 * t24;
t180 = -pkin(3) * t137 * t238 / 0.2e1;
t71 = 0.1e1 / t72 ^ 2;
t258 = pkin(4) * t71;
t282 = (-0.2e1 * pkin(4) * t74 - t69) * t180 + (t233 - t234) * pkin(3) * t152 * t258;
t22 = -t123 * t32 + t125 * t31;
t280 = t22 ^ 2;
t139 = sin(pkin(18));
t141 = cos(pkin(18));
t154 = 0.1e1 / pkin(8);
t271 = t70 / 0.2e1;
t200 = t154 * t271;
t68 = t72 + t213;
t73 = -pkin(3) - t256;
t49 = -pkin(4) * t283 - t68 * t73;
t257 = pkin(4) * t177;
t51 = -t164 * t73 + t68 * t257;
t40 = atan2(t51 * t200, t49 * t200);
t264 = sin(t40);
t39 = cos(t40);
t29 = t139 * t39 - t141 * t264;
t30 = -t139 * t264 - t141 * t39;
t138 = sin(pkin(19));
t226 = t106 * t138;
t140 = cos(pkin(19));
t229 = t105 * t140;
t96 = (-t229 / 0.2e1 + t226 / 0.2e1) * t219;
t225 = t106 * t140;
t230 = t105 * t138;
t97 = (t225 / 0.2e1 + t230 / 0.2e1) * t219;
t87 = atan2(t97, t96);
t81 = sin(t87);
t82 = cos(t87);
t60 = t144 * t81 - t149 * t82;
t62 = t144 * t82 + t149 * t81;
t20 = -t29 * t60 + t30 * t62;
t279 = 0.2e1 * t20;
t113 = t118 - t212;
t120 = pkin(1) * t124 - pkin(5);
t104 = -pkin(1) * t217 - t113 * t120;
t107 = pkin(1) * t126 * t113 - t120 * t163;
t146 = sin(pkin(15));
t156 = 0.1e1 / pkin(6);
t220 = t116 * t156;
t150 = cos(pkin(15));
t250 = t150 / 0.2e1;
t101 = (t107 * t250 - t104 * t146 / 0.2e1) * t220;
t98 = (t104 * t250 + t107 * t146 / 0.2e1) * t220;
t88 = atan2(t101, t98);
t84 = sin(t88);
t278 = 0.2e1 * t84;
t277 = 0.1e1 / t37 ^ 2;
t128 = -pkin(1) * t149 - pkin(13);
t115 = -pkin(4) * t123 + t128;
t276 = 0.2e1 * t115;
t275 = -0.2e1 * t126 ^ 2;
t274 = 0.2e1 * t128;
t273 = 0.2e1 * t149;
t59 = 0.1e1 / t164;
t272 = -t59 / 0.2e1;
t211 = pkin(1) * t247;
t222 = 0.2e1 / t163 * (t111 + t112) * t211;
t198 = -t222 / 0.2e1;
t218 = t124 * t163;
t90 = (t218 + (t119 * t281 - t114 + t198) * t126) * pkin(5);
t270 = -t90 / 0.2e1;
t92 = t119 * t222 / 0.2e1 + t157 * pkin(1) * t275 + (-t114 * t124 - t217) * pkin(5);
t269 = t92 / 0.2e1;
t263 = pkin(1) * t81;
t262 = pkin(1) * t82;
t192 = 0.1e1 / t118 ^ 2 * t211;
t252 = t138 / 0.2e1;
t94 = 0.1e1 / t96 ^ 2;
t45 = 0.1e1 + (((t140 * t269 + t90 * t252) * t116 + (t225 + t230) * t192) / t96 - ((t140 * t270 + t92 * t252) * t116 + (t226 - t229) * t192) * t97 * t94) * t161 / (t94 * t97 ^ 2 + 0.1e1);
t261 = pkin(3) * t45;
t260 = pkin(4) * t31;
t259 = pkin(4) * t32;
t166 = pkin(4) * (t71 * t232 + (t50 * t71 - t70 * t74) * t136);
t196 = -0.2e1 * t159 * t257;
t203 = t59 * t74 / 0.2e1;
t174 = pkin(4) * (t66 + t67) * t287;
t53 = t177 * t174;
t167 = t177 * t196 + t53 * t203;
t202 = t177 * t272;
t170 = -t77 * t164 + t53 * t202;
t178 = t69 * t77 - t283;
t188 = t238 * t254;
t253 = t137 / 0.2e1;
t201 = t70 * t253;
t33 = 0.1e1 / (t38 ^ 2 * t277 + 0.1e1);
t206 = t33 * t38 * t277;
t207 = t152 * t33 / t37;
t10 = 0.1e1 - (t170 * t180 + (t178 * pkin(3) + t167) * t188 + t282 * t177) * t206 + (t167 * t201 + ((t178 * t253 + (t170 - t284) * t254) * t70 + t177 * t166) * pkin(3)) * t207;
t251 = t143 / 0.2e1;
t64 = ((t148 * t269 + t90 * t251) * t116 + (t223 + t228) * t192) * t161;
t65 = ((t148 * t270 + t92 * t251) * t116 + (t224 - t227) * t192) * t161;
t57 = -t131 * t64 - t132 * t65;
t42 = t57 * t174;
t168 = t57 * t196 + t42 * t203;
t56 = t131 * t65 - t132 * t64;
t171 = -t56 * t164 + t42 * t202;
t237 = t164 * t57;
t179 = t56 * t69 - t237;
t5 = 0.1e1 - (t171 * t180 + (t179 * pkin(3) + t168) * t188 + t282 * t57) * t206 + (t168 * t201 + ((t179 * t253 + (-t57 * t69 + t171) * t254) * t70 + t57 * t166) * pkin(3)) * t207;
t255 = t5 * t10;
t142 = sin(qJ(4));
t246 = t142 * t5;
t245 = t143 * pkin(1);
t147 = cos(qJ(4));
t244 = t147 * t5;
t243 = t148 * pkin(1);
t242 = t142 * t10;
t241 = t142 * t22;
t240 = t147 * t10;
t239 = t147 * t22;
t221 = t116 * t146;
t216 = t142 * t147;
t134 = t142 ^ 2;
t135 = t147 ^ 2;
t214 = t134 + t135;
t210 = -0.2e1 * t241;
t209 = 0.2e1 * t239;
t208 = 0.2e1 * t216;
t205 = t24 * t216;
t204 = t73 * t272;
t197 = t116 * t250;
t7 = pkin(11) * t10 + t260;
t195 = t214 * t7;
t193 = (-t134 + t135) * t24;
t48 = 0.1e1 / t49 ^ 2;
t191 = pkin(8) * t154 / (t48 * t51 ^ 2 + 0.1e1) * t72;
t190 = t5 * t205;
t189 = t10 * t205;
t130 = pkin(4) - t243;
t26 = t31 * t130 - t32 * t245;
t187 = 0.2e1 * t214;
t25 = t130 * t32 + t31 * t245;
t2 = -pkin(9) * t5 - t25;
t8 = -pkin(9) * t10 - t259;
t186 = t10 * t2 + t5 * t8;
t185 = 0.1e1 / t49 * t191;
t3 = pkin(11) * t5 + t26;
t184 = t2 * t24 - t22 * t3;
t183 = -t22 * t7 + t24 * t8;
t182 = t104 * t192;
t181 = t107 * t192;
t173 = pkin(3) * (t49 * t71 + t70 * t73);
t172 = pkin(4) * t48 * t51 * t191;
t169 = pkin(3) * (-t158 * t177 * t70 + t51 * t258);
t127 = t128 ^ 2;
t95 = 0.1e1 / t98 ^ 2;
t91 = t120 * t198 + t162 * pkin(5) * t275 + (-t113 * t124 - t217) * pkin(1);
t89 = (t218 + (0.2e1 * t120 * pkin(5) - t113 + t198) * t126) * pkin(1);
t85 = cos(t88);
t54 = (-t139 * t62 + t141 * t60) * pkin(3) + t128;
t46 = ((t91 * t197 + t150 * t181 - t89 * t221 / 0.2e1 - t146 * t182) / t98 - (t89 * t197 + t150 * t182 + t91 * t221 / 0.2e1 + t146 * t181) * t101 * t95) / (t101 ^ 2 * t95 + 0.1e1) * t156;
t44 = t141 * t261 + t262;
t43 = t139 * t261 + t263;
t21 = t24 ^ 2;
t19 = -t29 * t62 - t30 * t60;
t18 = t29 * t44 + t30 * t43;
t17 = -t29 * t43 + t30 * t44;
t16 = 0.2e1 * ((t53 * t204 + (t77 * t68 - t283) * pkin(4)) * t271 + t177 * t169) * t185 - 0.2e1 * ((-t177 * t68 + t170) * t271 + t177 * t173) * t172;
t14 = pkin(9) * t22 - pkin(11) * t24 + t115;
t12 = 0.2e1 * ((t42 * t204 + (t56 * t68 - t237) * pkin(4)) * t271 + t57 * t169) * t185 - 0.2e1 * ((-t57 * t68 + t171) * t271 + t57 * t173) * t172 + t45;
t9 = t10 ^ 2;
t4 = t5 ^ 2;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t144 ^ 2, t144 * t273, 0, t149 ^ 2, 0, 0, pkin(13) * t273, -0.2e1 * pkin(13) * t144, 0, pkin(13) ^ 2, t125 ^ 2, 0.2e1 * t123 * t125, 0, t123 ^ 2, 0, 0, -0.2e1 * t128 * t123, t125 * t274, 0, t127, t21, t22 * t286, 0, t280, 0, 0, t22 * t276, t24 * t276, 0, t115 ^ 2, t135 * t21, -0.2e1 * t21 * t216, t24 * t209, t134 * t21, t24 * t210, t280, t14 * t209, t14 * t210, t214 * t14 * t286, t214 * t14 ^ 2, t84 ^ 2, t85 * t278, 0, t85 ^ 2, 0, 0, -0.2e1 * pkin(7) * t85, pkin(7) * t278, 0, pkin(7) ^ 2, t62 ^ 2, -0.2e1 * t60 * t62, 0, t60 ^ 2, 0, 0, t60 * t274, t62 * t274, 0, t127, t20 ^ 2, t19 * t279, 0, t19 ^ 2, 0, 0, -0.2e1 * t54 * t19, t54 * t279, 0, t54 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, t149, 0, 0, 0, 0, 0, 0, 0, t125, 0, t123, 0, 0, 0, (-t123 * t143 + t125 * t148) * pkin(1), 0, 0, 0, t5 * t24, 0, -t22 * t5, 0, 0, 0, -t22 * t26 - t24 * t25, 0, t190, t5 * t193, t5 * t241, -t190, t5 * t239, 0, t184 * t142, t184 * t147, 0, 0, 0, 0, t46 * t84, 0, t46 * t85, 0, 0, 0, 0, 0, 0, 0, t45 * t62, 0, -t45 * t60, 0, 0, 0, (-t60 * t81 - t62 * t82) * pkin(1), 0, 0, 0, t12 * t20, 0, t12 * t19, 0, 0, 0, -t17 * t20 + t18 * t19, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -0.2e1 * t243, 0.2e1 * t245, 0, (t143 ^ 2 + t148 ^ 2) * t162, 0, 0, 0, 0, 0, t4, 0.2e1 * t25 * t5, -0.2e1 * t26 * t5, 0, t25 ^ 2 + t26 ^ 2, t134 * t4, t4 * t208, 0, t135 * t4, 0, 0, -0.2e1 * t2 * t244, 0.2e1 * t2 * t246, t5 * t3 * t187, t214 * t3 ^ 2 + t2 ^ 2, 0, 0, 0, 0, 0, t46 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 ^ 2, 0.2e1 * t45 * t262, -0.2e1 * t45 * t263, 0, (t81 ^ 2 + t82 ^ 2) * t162, 0, 0, 0, 0, 0, t12 ^ 2, 0.2e1 * t17 * t12, -0.2e1 * t18 * t12, 0, t17 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, t123, 0, 0, 0, 0, 0, 0, 0, t10 * t24, 0, -t22 * t10, 0, 0, 0, (-t22 * t31 - t24 * t32) * pkin(4), 0, t189, t10 * t193, t10 * t241, -t189, t10 * t239, 0, t183 * t142, t183 * t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t20, 0, t16 * t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t243, t245, 0, 0, 0, 0, 0, 0, 0, t255, t10 * t25 + t5 * t259, -t10 * t26 - t5 * t260, 0, (t25 * t32 + t26 * t31) * pkin(4), t134 * t255, t208 * t255, 0, t135 * t255, 0, 0, -t186 * t147, t186 * t142, t10 * t214 * t3 + t195 * t5, t195 * t3 + t2 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t16, t17 * t16, -t18 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0.2e1 * t10 * t259, -0.2e1 * t10 * t260, 0, (t31 ^ 2 + t32 ^ 2) * t158, t134 * t9, t9 * t208, 0, t135 * t9, 0, 0, -0.2e1 * t8 * t240, 0.2e1 * t8 * t242, t7 * t10 * t187, t214 * t7 ^ 2 + t8 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t24, 0, -t142 * t24, t22, t147 * t14, -t142 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, 0, t244, 0, -t142 * t3, -t147 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, 0, t240, 0, -t142 * t7, -t147 * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MM_reg = t1;
