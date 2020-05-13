% Calculate inertial parameters regressor of gravitation load for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh1m1DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_gravloadJ_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t184 = pkin(5) ^ 2;
t155 = pkin(23) + pkin(22);
t151 = sin(t155);
t152 = cos(t155);
t182 = pkin(7) ^ 2;
t190 = pkin(1) ^ 2;
t165 = sin(qJ(2));
t170 = cos(qJ(2));
t172 = cos(pkin(19));
t300 = sin(pkin(19));
t142 = t165 * t172 - t170 * t300;
t293 = pkin(7) * t142;
t325 = -2 * pkin(1);
t256 = t293 * t325 + t190;
t133 = t182 + t256;
t253 = pkin(3) ^ 2 - pkin(8) ^ 2;
t125 = t133 + t253;
t136 = pkin(1) - t293;
t316 = -pkin(8) - pkin(3);
t117 = (pkin(7) - t316) * (pkin(7) + t316) + t256;
t315 = -pkin(8) + pkin(3);
t118 = (pkin(7) - t315) * (pkin(7) + t315) + t256;
t192 = sqrt(-t118 * t117);
t143 = t165 * t300 + t170 * t172;
t292 = pkin(7) * t143;
t106 = t125 * t292 + t136 * t192;
t169 = cos(qJ(3));
t258 = t169 * t106;
t269 = t143 * t192;
t105 = -pkin(7) * t269 + t125 * t136;
t164 = sin(qJ(3));
t263 = t164 * t105;
t131 = 0.1e1 / t133;
t187 = 0.1e1 / pkin(3);
t272 = t131 * t187;
t94 = (t263 / 0.2e1 + t258 / 0.2e1) * t272;
t259 = t169 * t105;
t262 = t164 * t106;
t95 = (-t259 / 0.2e1 + t262 / 0.2e1) * t272;
t63 = t151 * t95 - t152 * t94;
t310 = pkin(5) * t63;
t279 = -0.2e1 * pkin(4) * t310 + t184;
t314 = -pkin(9) - pkin(11);
t51 = (pkin(4) - t314) * (pkin(4) + t314) + t279;
t313 = -pkin(9) + pkin(11);
t52 = (pkin(4) - t313) * (pkin(4) + t313) + t279;
t193 = sqrt(-t52 * t51);
t229 = t151 * t94 + t152 * t95;
t327 = t193 * t229;
t166 = sin(qJ(1));
t171 = cos(qJ(1));
t216 = g(1) * t166 - g(2) * t171;
t326 = -g(1) * t171 - g(2) * t166;
t156 = qJ(2) + qJ(3);
t161 = cos(pkin(21));
t177 = 0.1e1 / pkin(11);
t185 = pkin(4) ^ 2;
t57 = t185 + t279;
t55 = 0.1e1 / t57;
t283 = t177 * t55;
t158 = sin(pkin(21));
t304 = t158 / 0.2e1;
t255 = pkin(9) ^ 2 - pkin(11) ^ 2;
t54 = t57 - t255;
t59 = -pkin(4) * t63 + pkin(5);
t35 = -pkin(4) * t327 + t54 * t59;
t37 = pkin(4) * t229 * t54 + t193 * t59;
t25 = (t35 * t304 + t37 * t161 / 0.2e1) * t283;
t26 = (-t35 * t161 / 0.2e1 + t37 * t304) * t283;
t16 = atan2(t25, t26) + t156;
t14 = sin(t16);
t15 = cos(t16);
t205 = -g(3) * t15 - t14 * t326;
t160 = cos(pkin(23));
t267 = t160 * t105;
t157 = sin(pkin(23));
t268 = t157 * t106;
t91 = (-t267 / 0.2e1 + t268 / 0.2e1) * t272;
t324 = 0.1e1 / t91 ^ 2;
t323 = -0.2e1 * t143 ^ 2;
t45 = 0.1e1 / t193;
t322 = -t45 / 0.2e1;
t321 = t55 / 0.2e1;
t251 = pkin(1) * t292;
t276 = 0.2e1 / t192 * (t117 + t118) * t251;
t238 = -t276 / 0.2e1;
t270 = t142 * t192;
t81 = (t270 + (t136 * t325 - t125 + t238) * t143) * pkin(7);
t320 = -t81 / 0.2e1;
t86 = t136 * t276 / 0.2e1 + t182 * pkin(1) * t323 + (-t142 * t125 - t269) * pkin(7);
t319 = t86 / 0.2e1;
t318 = (-pkin(2) - pkin(13));
t317 = (-pkin(2) + pkin(13));
t266 = t160 * t106;
t278 = t105 * t157;
t92 = (t266 / 0.2e1 + t278 / 0.2e1) * t272;
t70 = qJ(2) + atan2(t92, t91);
t67 = pkin(22) - t70;
t312 = pkin(4) * sin(t67);
t56 = 0.1e1 / t57 ^ 2;
t311 = pkin(5) * t56;
t309 = pkin(5) * t229;
t308 = g(3) * t14;
t159 = sin(pkin(20));
t162 = cos(pkin(20));
t139 = t159 * t169 + t162 * t164;
t295 = pkin(6) * t139;
t250 = pkin(1) * t295;
t135 = 0.2e1 * t250;
t183 = pkin(6) ^ 2;
t254 = t183 + t190;
t130 = t135 + t254;
t128 = 0.1e1 / t130;
t306 = t128 / 0.2e1;
t305 = t157 / 0.2e1;
t303 = -t169 / 0.2e1;
t173 = cos(pkin(18));
t302 = t173 / 0.2e1;
t189 = 0.1e1 / pkin(2);
t301 = t189 / 0.2e1;
t140 = t159 * t164 - t162 * t169;
t299 = pkin(1) * t140;
t298 = pkin(1) * t165;
t153 = sin(t156);
t297 = pkin(5) * t153;
t188 = pkin(2) ^ 2;
t241 = -pkin(13) ^ 2 + t254;
t123 = t135 + t188 + t241;
t134 = -pkin(1) - t295;
t257 = t135 + t183;
t115 = ((pkin(1) - t318) * (pkin(1) + t318)) + t257;
t116 = ((pkin(1) - t317) * (pkin(1) + t317)) + t257;
t275 = t116 * t115;
t191 = sqrt(-t275);
t294 = pkin(6) * t140;
t103 = t123 * t294 - t134 * t191;
t296 = pkin(6) * t103;
t287 = t170 * pkin(1);
t24 = 0.1e1 / t26 ^ 2;
t286 = t24 * t25;
t271 = t140 * t191;
t249 = pkin(6) * t271;
t102 = -t123 * t134 - t249;
t101 = 0.1e1 / t102 ^ 2;
t129 = 0.1e1 / t130 ^ 2;
t277 = 0.2e1 / t191 * (t115 + t116) * pkin(1) * t294;
t239 = -t277 / 0.2e1;
t43 = 0.2e1 * (((t134 * t239 + (t139 * t123 - t271) * pkin(6)) * t306 + (-t128 * t140 * t183 + t129 * t296) * t299) / t102 - ((-t139 * t191 + (t239 - t123) * t140) * t306 + (t102 * t129 + t128 * t134) * t299) * t101 * t296) * pkin(2) * t130 * t189 / (t101 * t103 ^ 2 + 0.1e1);
t285 = t161 * t55;
t284 = 0.1e1 / (t24 * t25 ^ 2 + 0.1e1) * t177;
t228 = 0.1e1 / t133 ^ 2 * t251;
t48 = ((t164 * t320 + t86 * t303) * t131 + (-t258 - t263) * t228) * t187;
t49 = ((t164 * t319 + t81 * t303) * t131 + (-t259 + t262) * t228) * t187;
t41 = t151 * t49 + t152 * t48;
t282 = t193 * t41;
t237 = t128 * t301;
t87 = qJ(2) + atan2(t103 * t237, t102 * t237);
t167 = sin(pkin(18));
t274 = t131 * t167;
t181 = 0.1e1 / pkin(8);
t273 = t131 * t181;
t163 = sin(qJ(4));
t265 = t163 * t166;
t264 = t163 * t171;
t168 = cos(qJ(4));
t261 = t166 * t168;
t260 = t168 * t171;
t252 = pkin(4) * t311;
t58 = -pkin(4) + t310;
t248 = t58 * t322;
t247 = t45 * t59 / 0.2e1;
t246 = t229 * t322;
t245 = t55 * t304;
t244 = -t285 / 0.2e1;
t243 = t285 / 0.2e1;
t179 = 0.1e1 / pkin(9);
t242 = t179 * t321;
t240 = -0.2e1 * t59 * pkin(5) - t54;
t236 = t131 * t302;
t235 = 0.1e1 / pkin(13) * t301;
t234 = -0.2e1 * t185 * t309;
t233 = pkin(16) - t298;
t232 = t158 * t252;
t231 = t161 * t252;
t154 = cos(t156);
t148 = pkin(5) * t154;
t230 = t148 - t298;
t53 = t57 + t255;
t34 = -pkin(5) * t327 - t53 * t58;
t33 = 0.1e1 / t34 ^ 2;
t36 = -t193 * t58 + t53 * t309;
t227 = pkin(9) * t179 / (t33 * t36 ^ 2 + 0.1e1) * t57;
t225 = t41 * t232;
t224 = t229 * t232;
t223 = t41 * t231;
t222 = t229 * t231;
t221 = 0.1e1 / t34 * t227;
t220 = pkin(10) * t15 + pkin(12) * t14;
t219 = -pkin(10) * t14 + pkin(12) * t15;
t124 = t133 - t253;
t137 = pkin(1) * t142 - pkin(7);
t104 = -pkin(1) * t269 - t124 * t137;
t218 = t104 * t228;
t107 = pkin(1) * t143 * t124 - t137 * t192;
t217 = t107 * t228;
t30 = 0.1e1 + (((t160 * t319 + t81 * t305) * t131 + (t266 + t278) * t228) / t91 - ((t160 * t320 + t86 * t305) * t131 + (-t267 + t268) * t228) * t92 * t324) * t187 / (t92 ^ 2 * t324 + 0.1e1);
t215 = 0.2e1 * pkin(4) * pkin(5) * (t51 + t52);
t213 = pkin(4) * (t34 * t56 + t55 * t58);
t212 = pkin(5) * t33 * t36 * t227;
t28 = t41 * t215;
t42 = -t151 * t48 + t152 * t49;
t211 = -t42 * t193 + t246 * t28;
t38 = t229 * t215;
t210 = t193 * t63 + t246 * t38;
t208 = pkin(4) * (-t184 * t229 * t55 + t36 * t311);
t206 = -t15 * t326 + t308;
t22 = -atan2(t36 * t242, t34 * t242) + t67;
t20 = sin(t22);
t21 = cos(t22);
t204 = -g(3) * t20 - t21 * t326;
t203 = g(3) * t21 - t20 * t326;
t122 = t188 - t241 - 0.2e1 * t250;
t79 = atan2(t191 * t235, t122 * t235) + t87;
t77 = sin(t79);
t78 = cos(t79);
t202 = -g(3) * t77 + t326 * t78;
t201 = -g(3) * t78 - t326 * t77;
t83 = sin(t87);
t84 = cos(t87);
t200 = g(3) * t83 - t326 * t84;
t199 = g(3) * t84 + t326 * t83;
t120 = -g(3) * t154 - t153 * t326;
t198 = g(3) * t165 - t170 * t326;
t197 = t200 * t43;
t196 = t163 * t205;
t195 = t168 * t205;
t146 = -t287 - t297;
t145 = t216 * pkin(16);
t144 = pkin(16) + t230;
t127 = t216 * t233;
t126 = t198 * pkin(1);
t121 = 0.1e1 / t122 ^ 2;
t119 = g(3) * t153 - t154 * t326;
t97 = (t107 * t302 + t104 * t167 / 0.2e1) * t273;
t96 = (t104 * t302 - t167 * t107 / 0.2e1) * t273;
t93 = 0.1e1 / t96 ^ 2;
t85 = t137 * t238 + t190 * pkin(7) * t323 + (-t142 * t124 - t269) * pkin(1);
t80 = (t270 + (0.2e1 * t137 * pkin(7) - t124 + t238) * t143) * pkin(1);
t76 = atan2(t97, t96);
t73 = cos(t76);
t72 = sin(t76);
t69 = cos(t70);
t68 = sin(t70);
t39 = (0.1e1 / t122 * t277 / 0.2e1 + t121 * t249 * t325) / (-t121 * t275 + 0.1e1) + t43;
t31 = ((t85 * t236 + t173 * t217 + t80 * t274 / 0.2e1 + t167 * t218) / t96 - (t80 * t236 + t173 * t218 - t85 * t274 / 0.2e1 - t167 * t217) * t97 * t93) / (t93 * t97 ^ 2 + 0.1e1) * t181;
t23 = 0.1e1 / t26;
t19 = t38 * t247 + t229 * t234 + (-t54 * t63 - t327) * pkin(4);
t18 = (t229 * t240 + t210) * pkin(4);
t13 = t15 * t260 + t265;
t12 = -t15 * t264 + t261;
t11 = -t15 * t261 + t264;
t10 = t15 * t265 + t260;
t9 = t28 * t247 + t41 * t234 + (t42 * t54 - t282) * pkin(4);
t8 = (t240 * t41 + t211) * pkin(4);
t7 = t216 * t14;
t6 = -0.2e1 * ((t38 * t248 + (-t53 * t63 - t327) * pkin(5)) * t321 + t229 * t208) * t221 + 0.2e1 * ((-t229 * t53 + t210) * t321 + t229 * t213) * t212;
t5 = -0.2e1 * ((t28 * t248 + (t42 * t53 - t282) * pkin(5)) * t321 + t41 * t208) * t221 + 0.2e1 * ((-t41 * t53 + t211) * t321 + t41 * t213) * t212 - t30;
t4 = 0.1e1 + ((t18 * t245 + t19 * t243 + t222 * t37 + t224 * t35) * t23 - (t18 * t244 + t19 * t245 - t222 * t35 + t224 * t37) * t286) * t284;
t3 = 0.1e1 + ((t223 * t37 + t225 * t35 + t243 * t9 + t245 * t8) * t23 - (-t223 * t35 + t225 * t37 + t244 * t8 + t245 * t9) * t286) * t284;
t2 = t206 * t4;
t1 = t206 * t3;
t17 = [0, 0, 0, 0, 0, 0, t216, -t326, 0, 0, 0, 0, 0, 0, 0, 0, -t216 * t165, -t216 * t170, t326, t145, 0, 0, 0, 0, 0, 0, t216 * t154, -t216 * t153, t326, t127, 0, 0, 0, 0, 0, 0, t216 * t15, -t7, t326, t216 * t144, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t7, t216 * (t144 + t220), 0, 0, 0, 0, 0, 0, t216 * t73, -t216 * t72, t326, -t216 * pkin(15), 0, 0, 0, 0, 0, 0, -t216 * t68, -t216 * t69, t326, t127, 0, 0, 0, 0, 0, 0, -t216 * t83, -t216 * t84, t326, t145, 0, 0, 0, 0, 0, 0, t216 * t77, t216 * t78, t326, t216 * (-pkin(2) * t83 + pkin(16)), 0, 0, 0, 0, 0, 0, -t216 * t20, t216 * t21, t326, t216 * (t233 + t312); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, g(3) * t170 + t165 * t326, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, 0, t126, 0, 0, 0, 0, 0, 0, t205 * t3, t1, 0, -g(3) * t230 + t326 * t146, 0, 0, 0, 0, 0, 0, t3 * t195, -t3 * t196, -t1, -g(3) * (t220 * t3 + t230) + t326 * (t219 * t3 + t146), 0, 0, 0, 0, 0, 0, (-g(3) * t73 - t326 * t72) * t31, (g(3) * t72 - t326 * t73) * t31, 0, 0, 0, 0, 0, 0, 0, 0, (g(3) * t68 - t326 * t69) * t30, (g(3) * t69 + t326 * t68) * t30, 0, t126, 0, 0, 0, 0, 0, 0, t200, t199, 0, 0, 0, 0, 0, 0, 0, 0, t202, t201, 0, t200 * pkin(2), 0, 0, 0, 0, 0, 0, t204 * t5, t203 * t5, 0, -g(3) * (t30 * t312 - t298) + t326 * (-pkin(4) * t30 * cos(t67) - t287); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, 0, 0, 0, 0, 0, 0, 0, 0, t205 * t4, t2, 0, t120 * pkin(5), 0, 0, 0, 0, 0, 0, t4 * t195, -t4 * t196, -t2, -g(3) * (t220 * t4 + t148) - t326 * (-t219 * t4 + t297), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t199 * t43, 0, 0, 0, 0, 0, 0, 0, 0, t202 * t39, t201 * t39, 0, pkin(2) * t197, 0, 0, 0, 0, 0, 0, t204 * t6, t203 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t10 + t163 * t308, g(1) * t13 - g(2) * t11 + t168 * t308, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t17;
