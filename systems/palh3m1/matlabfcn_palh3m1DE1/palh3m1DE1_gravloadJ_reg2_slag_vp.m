% Calculate inertial parameters regressor of gravitation load for
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh3m1DE1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_gravloadJ_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t167 = sin(qJ(4));
t172 = cos(qJ(4));
t175 = cos(qJ(1));
t169 = sin(qJ(1));
t168 = sin(qJ(2));
t173 = cos(qJ(3));
t174 = cos(qJ(2));
t299 = sin(qJ(3));
t214 = t168 * t173 + t174 * t299;
t137 = t214 * t169;
t250 = t168 * t299;
t146 = -t174 * t173 + t250;
t138 = t146 * t169;
t162 = cos(pkin(17));
t177 = 0.1e1 / pkin(10);
t182 = pkin(3) ^ 2;
t181 = pkin(4) ^ 2;
t180 = pkin(5) ^ 2;
t185 = pkin(1) ^ 2;
t170 = sin(pkin(16));
t300 = cos(pkin(16));
t147 = t168 * t170 - t174 * t300;
t298 = pkin(5) * t147;
t316 = -2 * pkin(1);
t272 = t298 * t316 + t185;
t136 = t180 + t272;
t269 = pkin(2) ^ 2 - pkin(6) ^ 2;
t128 = t136 + t269;
t141 = pkin(1) - t298;
t312 = -pkin(6) - pkin(2);
t125 = (pkin(5) - t312) * (pkin(5) + t312) + t272;
t311 = -pkin(6) + pkin(2);
t126 = (pkin(5) - t311) * (pkin(5) + t311) + t272;
t186 = sqrt(-t126 * t125);
t149 = t168 * t300 + t174 * t170;
t297 = pkin(5) * t149;
t117 = t128 * t297 + t141 * t186;
t251 = t117 * t299;
t130 = 0.1e1 / t136;
t184 = 0.1e1 / pkin(2);
t279 = t130 * t184;
t277 = t149 * t186;
t116 = -pkin(5) * t277 + t128 * t141;
t285 = t116 * t173;
t110 = (-t285 / 0.2e1 + t251 / 0.2e1) * t279;
t252 = t116 * t299;
t284 = t117 * t173;
t111 = (t284 / 0.2e1 + t252 / 0.2e1) * t279;
t160 = pkin(18) + pkin(19);
t156 = sin(t160);
t157 = cos(t160);
t89 = t110 * t157 + t111 * t156;
t212 = pkin(4) * t89;
t286 = 0.2e1 * pkin(3) * t212 + t181;
t85 = t182 + t286;
t83 = 0.1e1 / t85;
t288 = t177 * t83;
t161 = sin(pkin(17));
t303 = t161 / 0.2e1;
t310 = -pkin(8) - pkin(10);
t80 = (pkin(3) - t310) * (pkin(3) + t310) + t286;
t309 = -pkin(8) + pkin(10);
t81 = (pkin(3) - t309) * (pkin(3) + t309) + t286;
t187 = sqrt(-t81 * t80);
t225 = t110 * t156 - t111 * t157;
t287 = t187 * t225;
t268 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t82 = t85 + t268;
t86 = pkin(3) * t89 + pkin(4);
t65 = -pkin(3) * t287 + t82 * t86;
t307 = pkin(3) * t225;
t66 = t187 * t86 + t82 * t307;
t59 = (-t65 * t162 / 0.2e1 + t66 * t303) * t288;
t60 = (t66 * t162 / 0.2e1 + t65 * t303) * t288;
t47 = atan2(t60, t59);
t45 = cos(t47);
t35 = t138 * t45;
t44 = sin(t47);
t26 = -t137 * t44 - t35;
t327 = t175 * t167 + t172 * t26;
t326 = t167 * t26 - t172 * t175;
t313 = pkin(3) * pkin(4);
t325 = 0.2e1 * t313;
t71 = 0.1e1 / t187;
t324 = -t71 / 0.2e1;
t267 = 0.1e1 / t85 ^ 2 * t313;
t244 = t162 * t267;
t266 = pkin(1) * t297;
t283 = 0.2e1 / t186 * (t125 + t126) * t266;
t249 = -t283 / 0.2e1;
t278 = t147 * t186;
t102 = (t278 + (t141 * t316 - t128 + t249) * t149) * pkin(5);
t243 = 0.1e1 / t136 ^ 2 * t266;
t264 = t299 / 0.2e1;
t314 = -0.2e1 * t149 ^ 2;
t104 = t141 * t283 / 0.2e1 + t180 * pkin(1) * t314 + (-t147 * t128 - t277) * pkin(5);
t304 = t104 / 0.2e1;
t78 = ((t102 * t264 + t173 * t304) * t130 + (t252 + t284) * t243) * t184;
t305 = -t102 / 0.2e1;
t79 = ((t104 * t264 + t173 * t305) * t130 + (t251 - t285) * t243) * t184;
t69 = -t156 * t78 - t157 * t79;
t239 = t69 * t244;
t245 = t161 * t267;
t241 = t69 * t245;
t293 = t162 * t83;
t255 = t293 / 0.2e1;
t256 = -t293 / 0.2e1;
t257 = t83 * t303;
t58 = 0.1e1 / t59 ^ 2;
t289 = t177 / (t58 * t60 ^ 2 + 0.1e1);
t296 = t58 * t60;
t260 = t225 * t324;
t224 = (t80 + t81) * t325;
t62 = t69 * t224;
t68 = t156 * t79 - t157 * t78;
t218 = -t68 * t187 + t62 * t260;
t315 = -0.2e1 * pkin(4);
t253 = t86 * t315 - t82;
t31 = (t253 * t69 + t218) * pkin(3);
t246 = t182 * t225 * t315;
t261 = t71 * t86 / 0.2e1;
t32 = t62 * t261 + t69 * t246 + (-t187 * t69 + t68 * t82) * pkin(3);
t57 = 0.1e1 / t59;
t15 = ((t66 * t239 + t65 * t241 + t32 * t255 + t31 * t257) * t57 - (-t65 * t239 + t66 * t241 + t31 * t256 + t32 * t257) * t296) * t289;
t323 = t15 + 0.1e1;
t238 = t225 * t244;
t240 = t225 * t245;
t67 = t225 * t224;
t217 = -t89 * t187 + t67 * t260;
t54 = (t225 * t253 + t217) * pkin(3);
t55 = t67 * t261 + t225 * t246 + (t89 * t82 - t287) * pkin(3);
t18 = ((t66 * t238 + t65 * t240 + t55 * t255 + t54 * t257) * t57 - (-t65 * t238 + t66 * t240 + t54 * t256 + t55 * t257) * t296) * t289;
t322 = t18 + 0.1e1;
t295 = t83 * t85;
t165 = cos(pkin(19));
t163 = sin(pkin(19));
t302 = t163 / 0.2e1;
t107 = (-t116 * t165 / 0.2e1 + t117 * t302) * t279;
t108 = (t117 * t165 / 0.2e1 + t116 * t302) * t279;
t99 = atan2(t108, t107);
t93 = sin(t99);
t290 = t168 * t93;
t94 = cos(t99);
t92 = t174 * t94;
t76 = t92 - t290;
t211 = -t212 - pkin(3);
t228 = t85 - t268;
t206 = -pkin(4) * t287 - t211 * t228;
t205 = 0.1e1 / t206 ^ 2;
t220 = pkin(4) * t228;
t207 = -t211 * t187 + t220 * t225;
t204 = 0.1e1 / (t207 ^ 2 * t205 + 0.1e1);
t198 = t204 * t205 * t207;
t199 = t204 / t206;
t208 = t211 * t325 - t220;
t209 = t211 * t324;
t223 = -pkin(4) * t187 - 0.2e1 * t181 * t307;
t321 = (t67 * t209 + t89 * t220 + t223 * t225) * t199 - (t217 * pkin(4) + t208 * t225) * t198;
t320 = (t62 * t209 + t68 * t220 + t223 * t69) * t199 - (t218 * pkin(4) + t208 * t69) * t198;
t164 = sin(pkin(18));
t166 = cos(pkin(18));
t265 = t83 / pkin(8) / 0.2e1;
t203 = atan2(t207 * t265, t206 * t265);
t202 = sin(t203);
t61 = cos(t203);
t319 = -t164 * t61 + t166 * t202;
t43 = -t164 * t202 - t166 * t61;
t318 = t295 * t43;
t317 = t319 * t295;
t105 = 0.1e1 / t107 ^ 2;
t232 = t117 * t243;
t233 = t116 * t243;
t248 = t130 * t302;
t282 = t130 * t165;
t63 = ((t102 * t248 + t163 * t233 + t165 * t232 + t282 * t304) / t107 - (t104 * t248 + t163 * t232 - t165 * t233 + t282 * t305) * t108 * t105) / (t105 * t108 ^ 2 + 0.1e1) * t184;
t308 = t63 + 0.1e1;
t41 = t214 * t45;
t306 = g(3) * (t146 * t44 - t41);
t176 = cos(pkin(15));
t301 = t176 / 0.2e1;
t159 = t174 * pkin(1);
t273 = t174 * t175;
t139 = t173 * t273 - t175 * t250;
t140 = t214 * t175;
t227 = t139 * t44 + t140 * t45;
t226 = t146 * t45 + t214 * t44;
t33 = t137 * t45;
t37 = t139 * t45;
t171 = sin(pkin(15));
t281 = t130 * t171;
t179 = 0.1e1 / pkin(6);
t280 = t130 * t179;
t276 = t168 * t169;
t144 = t146 * pkin(4);
t271 = t144 + t159;
t270 = pkin(1) * t273 + t175 * pkin(13);
t263 = t323 * t44;
t262 = t322 * t44;
t259 = t168 * t308;
t258 = t174 * t308;
t254 = -t139 * pkin(4) + t270;
t247 = t130 * t301;
t242 = t94 * t258;
t132 = t137 * pkin(4);
t237 = -pkin(1) * t276 + t132;
t134 = t140 * pkin(4);
t236 = -pkin(1) * t168 * t175 + t134;
t25 = -t138 * t44 + t33;
t235 = g(1) * t25 - g(2) * t227;
t127 = t136 - t269;
t142 = pkin(1) * t147 - pkin(5);
t115 = -pkin(1) * t277 - t127 * t142;
t234 = t115 * t243;
t118 = pkin(1) * t149 * t127 - t142 * t186;
t231 = t118 * t243;
t150 = g(1) * t175 + g(2) * t169;
t230 = g(1) * t169 - g(2) * t175;
t229 = (-pkin(13) - t159) * t169;
t77 = t168 * t94 + t174 * t93;
t1 = t137 * t263 + t15 * t35 + t35;
t3 = t140 * t263 - t15 * t37 - t37;
t5 = t146 * t263 - t15 * t41 - t41;
t222 = g(1) * t3 + g(2) * t1 + g(3) * t5;
t2 = -t138 * t263 + t15 * t33 + t33;
t4 = t323 * t227;
t6 = t323 * t226;
t221 = g(1) * t4 + g(2) * t2 + g(3) * t6;
t219 = t150 * t168;
t11 = t146 * t262 - t18 * t41 - t41;
t7 = t137 * t262 + t18 * t35 + t35;
t9 = t140 * t262 - t18 * t37 - t37;
t216 = g(1) * t9 + g(2) * t7 + g(3) * t11;
t215 = -t138 * pkin(4) + t229;
t10 = t322 * t227;
t12 = t322 * t226;
t8 = -t138 * t262 + t18 * t33 + t33;
t213 = g(1) * t10 + g(2) * t8 + g(3) * t12;
t120 = -g(1) * t140 - g(2) * t137 - g(3) * t146;
t52 = -t93 * t258 - t94 * t259;
t210 = -g(3) * t174 + t219;
t124 = -g(1) * t229 - g(2) * t270;
t129 = t210 * pkin(1);
t121 = -g(1) * t139 + g(2) * t138 - g(3) * t214;
t112 = (t118 * t301 - t115 * t171 / 0.2e1) * t280;
t109 = (t115 * t301 + t118 * t171 / 0.2e1) * t280;
t106 = 0.1e1 / t109 ^ 2;
t103 = t142 * t249 + t185 * pkin(5) * t314 + (-t147 * t127 - t277) * pkin(1);
t101 = (t278 + (0.2e1 * t142 * pkin(5) - t127 + t249) * t149) * pkin(1);
t100 = atan2(t112, t109);
t97 = cos(t100);
t96 = sin(t100);
t91 = t93 * t276;
t75 = t76 * t175;
t74 = t77 * t175;
t73 = t169 * t92 - t91;
t72 = t77 * t169;
t64 = ((t103 * t247 + t176 * t231 - t101 * t281 / 0.2e1 - t171 * t234) / t109 - (t101 * t247 + t176 * t234 + t103 * t281 / 0.2e1 + t171 * t231) * t112 * t106) / (t106 * t112 ^ 2 + 0.1e1) * t179;
t53 = -t93 * t259 + t63 * t92 + t92;
t51 = t52 * t175;
t50 = (t308 * t290 - t242) * t175;
t49 = t52 * t169;
t48 = t91 + (t63 * t290 - t242) * t169;
t29 = t140 * t44 - t37;
t24 = t167 * t169 + t172 * t29;
t23 = -t167 * t29 + t169 * t172;
t17 = t321 * t317;
t16 = t321 * t318;
t14 = t320 * t317;
t13 = t320 * t318;
t19 = [0, 0, 0, 0, 0, 0, t230, t150, 0, 0, 0, 0, 0, 0, 0, 0, t230 * t174, -t230 * t168, -t150, t230 * pkin(13), 0, 0, 0, 0, 0, 0, g(1) * t138 + g(2) * t139, g(1) * t137 - g(2) * t140, -t150, t124, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t29, t235, -t150, -g(1) * t215 - g(2) * t254, 0, 0, 0, 0, 0, 0, -g(1) * t327 - g(2) * t24, g(1) * t326 - g(2) * t23, -t235, -g(2) * (pkin(9) * t29 - pkin(11) * t227 + t254) - g(1) * (pkin(9) * t26 + pkin(11) * t25 + t215), 0, 0, 0, 0, 0, 0, t230 * t97, -t230 * t96, -t150, -t230 * pkin(7), 0, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t75, -g(1) * t72 + g(2) * t74, -t150, t124, 0, 0, 0, 0, 0, 0, -g(2) * (t319 * t74 + t43 * t75) - g(1) * (-t319 * t72 - t43 * t73), -g(2) * (t319 * t75 - t43 * t74) - g(1) * (-t319 * t73 + t43 * t72), -t150, (-g(2) * (t164 * t74 + t166 * t75) - g(1) * (-t164 * t72 - t166 * t73)) * pkin(3) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, g(3) * t168 + t150 * t174, 0, 0, 0, 0, 0, 0, 0, 0, t120, t121, 0, t129, 0, 0, 0, 0, 0, 0, -t221, t222, 0, -g(1) * t236 - g(2) * t237 - g(3) * t271, 0, 0, 0, 0, 0, 0, -t221 * t172, t221 * t167, -t222, -g(3) * (pkin(9) * t6 + pkin(11) * t5 + t271) - g(2) * (pkin(9) * t2 + pkin(11) * t1 + t237) - g(1) * (pkin(9) * t4 + pkin(11) * t3 + t236), 0, 0, 0, 0, 0, 0, (-g(3) * t97 + t150 * t96) * t64, (g(3) * t96 + t150 * t97) * t64, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t51 - g(2) * t49 - g(3) * t53, -g(1) * t50 - g(2) * t48 - g(3) * t52, 0, t129, 0, 0, 0, 0, 0, 0, -g(3) * (t13 * t76 + t14 * t77 - t319 * t52 + t43 * t53) - g(2) * (-t13 * t72 + t14 * t73 - t319 * t48 + t43 * t49) - g(1) * (-t13 * t74 + t14 * t75 - t319 * t50 + t43 * t51), -g(3) * (-t13 * t77 + t14 * t76 + t319 * t53 + t43 * t52) - g(2) * (-t13 * t73 - t14 * t72 + t319 * t49 + t43 * t48) - g(1) * (-t13 * t75 - t14 * t74 + t319 * t51 + t43 * t50), 0, -g(3) * t159 + pkin(1) * t219 + (-g(3) * (-t164 * t52 + t166 * t53) - g(2) * (-t164 * t48 + t166 * t49) - g(1) * (-t164 * t50 + t166 * t51)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t121, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t216, 0, t120 * pkin(4), 0, 0, 0, 0, 0, 0, -t213 * t172, t213 * t167, -t216, -g(3) * (pkin(9) * t12 + pkin(11) * t11 + t144) - g(2) * (pkin(9) * t8 + pkin(11) * t7 + t132) - g(1) * (pkin(9) * t10 + pkin(11) * t9 + t134), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * (t16 * t76 + t17 * t77) - g(2) * (-t16 * t72 + t17 * t73) - g(1) * (-t16 * t74 + t17 * t75), -g(3) * (-t16 * t77 + t17 * t76) - g(2) * (-t16 * t73 - t17 * t72) - g(1) * (-t16 * t75 - t17 * t74), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t326 + t167 * t306, g(1) * t24 - g(2) * t327 + t172 * t306, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t19;
