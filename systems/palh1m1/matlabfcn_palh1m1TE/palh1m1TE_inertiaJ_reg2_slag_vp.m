% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh1m1TE_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_inertiaJ_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t204 = pkin(2) ^ 2;
t199 = pkin(6) ^ 2;
t205 = pkin(1) ^ 2;
t274 = t199 + t205;
t370 = pkin(13) ^ 2;
t253 = t274 - t370;
t174 = sin(pkin(20));
t178 = cos(pkin(20));
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t150 = t174 * t184 + t178 * t180;
t318 = t150 * pkin(6);
t271 = pkin(1) * t318;
t134 = t204 - t253 - 0.2e1 * t271;
t371 = t134 ^ 2;
t369 = 2 * pkin(4);
t151 = t174 * t180 - t178 * t184;
t322 = pkin(1) * t151;
t368 = 0.2e1 * t322;
t198 = pkin(7) ^ 2;
t181 = sin(qJ(2));
t185 = cos(qJ(2));
t186 = cos(pkin(19));
t325 = sin(pkin(19));
t154 = t181 * t186 - t185 * t325;
t320 = pkin(7) * t154;
t362 = -0.2e1 * pkin(1);
t277 = t320 * t362 + t205;
t144 = t198 + t277;
t273 = pkin(3) ^ 2 - pkin(8) ^ 2;
t137 = t144 + t273;
t147 = pkin(1) - t320;
t342 = -pkin(8) - pkin(3);
t131 = (pkin(7) - t342) * (pkin(7) + t342) + t277;
t341 = -pkin(8) + pkin(3);
t132 = (pkin(7) - t341) * (pkin(7) + t341) + t277;
t207 = sqrt(-t132 * t131);
t157 = t181 * t325 + t185 * t186;
t319 = pkin(7) * t157;
t121 = t137 * t319 + t147 * t207;
t280 = t184 * t121;
t288 = t157 * t207;
t120 = -pkin(7) * t288 + t137 * t147;
t283 = t180 * t120;
t142 = 0.1e1 / t144;
t203 = 0.1e1 / pkin(3);
t292 = t142 * t203;
t105 = (t283 / 0.2e1 + t280 / 0.2e1) * t292;
t281 = t184 * t120;
t282 = t180 * t121;
t106 = (-t281 / 0.2e1 + t282 / 0.2e1) * t292;
t168 = pkin(23) + pkin(22);
t164 = sin(t168);
t165 = cos(t168);
t241 = t105 * t164 + t165 * t106;
t275 = pkin(9) ^ 2 - pkin(11) ^ 2;
t201 = pkin(4) ^ 2;
t200 = pkin(5) ^ 2;
t79 = t105 * t165 - t106 * t164;
t334 = t79 * pkin(5);
t300 = t334 * t369 + t200;
t75 = t201 + t300;
t72 = t75 - t275;
t366 = t241 * t72;
t340 = (-pkin(9) - pkin(11));
t69 = ((pkin(4) - t340) * (pkin(4) + t340)) + t300;
t339 = (-pkin(9) + pkin(11));
t70 = ((pkin(4) - t339) * (pkin(4) + t339)) + t300;
t208 = sqrt(-t70 * t69);
t365 = t208 * t241;
t193 = 0.1e1 / pkin(11);
t177 = cos(pkin(21));
t73 = 0.1e1 / t75;
t310 = t193 * t73;
t229 = -pkin(4) * t177 * t310 / 0.2e1;
t173 = sin(pkin(21));
t77 = pkin(4) * t79 + pkin(5);
t50 = pkin(4) * t366 + t208 * t77;
t304 = t50 * t173;
t48 = -pkin(4) * t365 + t72 * t77;
t305 = t48 * t177;
t74 = 0.1e1 / t75 ^ 2;
t335 = pkin(5) * t74;
t364 = (-0.2e1 * t77 * pkin(5) - t72) * t229 + (t304 - t305) * pkin(4) * t193 * t335;
t363 = 0.1e1 / pkin(2);
t153 = -t185 * t180 - t181 * t184;
t156 = -t180 * t181 + t184 * t185;
t303 = t50 * t177;
t330 = t173 / 0.2e1;
t33 = (t48 * t330 + t303 / 0.2e1) * t310;
t34 = (-t305 / 0.2e1 + t304 / 0.2e1) * t310;
t22 = t153 * t33 + t156 * t34;
t361 = t22 ^ 2;
t360 = -2 * pkin(16);
t172 = sin(pkin(22));
t195 = 0.1e1 / pkin(9);
t309 = t195 * t73;
t176 = cos(pkin(22));
t329 = -t176 / 0.2e1;
t71 = t75 + t275;
t76 = -pkin(4) - t334;
t47 = -pkin(5) * t365 - t71 * t76;
t333 = pkin(5) * t241;
t49 = -t208 * t76 + t71 * t333;
t31 = (t172 * t47 / 0.2e1 + t49 * t329) * t309;
t32 = (t47 * t329 - t172 * t49 / 0.2e1) * t309;
t175 = cos(pkin(23));
t286 = t175 * t120;
t171 = sin(pkin(23));
t287 = t171 * t121;
t101 = (-t286 / 0.2e1 + t287 / 0.2e1) * t292;
t285 = t175 * t121;
t299 = t120 * t171;
t102 = (t285 / 0.2e1 + t299 / 0.2e1) * t292;
t83 = t101 * t185 - t102 * t181;
t84 = -t101 * t181 - t102 * t185;
t19 = t31 * t84 + t32 * t83;
t359 = 0.2e1 * t19;
t24 = -t153 * t34 + t156 * t33;
t358 = 0.2e1 * t24;
t357 = 0.2e1 * t83;
t146 = 0.2e1 * t271;
t135 = t146 + t204 + t253;
t145 = -pkin(1) - t318;
t278 = t146 + t199;
t344 = -pkin(2) - pkin(13);
t129 = (pkin(1) - t344) * (pkin(1) + t344) + t278;
t343 = -pkin(2) + pkin(13);
t130 = (pkin(1) - t343) * (pkin(1) + t343) + t278;
t296 = t130 * t129;
t206 = sqrt(-t296);
t290 = t151 * t206;
t266 = pkin(6) * t290;
t117 = -t135 * t145 - t266;
t291 = t151 * t135;
t118 = pkin(6) * t291 - t145 * t206;
t141 = t146 + t274;
t138 = 0.1e1 / t141;
t301 = t363 / 0.2e1;
t259 = t185 * t301;
t260 = t138 * t301;
t95 = t181 * t117 * t260 + t118 * t138 * t259;
t91 = pkin(2) * t95 - pkin(16);
t356 = -0.2e1 * t91;
t355 = -0.2e1 * t95;
t354 = 0.1e1 / t371;
t211 = t34 ^ 2;
t353 = 0.1e1 / t211;
t136 = t144 - t273;
t148 = pkin(1) * t154 - pkin(7);
t119 = -pkin(1) * t288 - t136 * t148;
t122 = pkin(1) * t157 * t136 - t148 * t207;
t182 = sin(pkin(18));
t197 = 0.1e1 / pkin(8);
t293 = t142 * t197;
t187 = cos(pkin(18));
t326 = t187 / 0.2e1;
t108 = (t122 * t326 + t119 * t182 / 0.2e1) * t293;
t352 = 0.2e1 * t108;
t351 = -0.2e1 * t157 ^ 2;
t159 = t181 * pkin(1) - pkin(16);
t350 = -0.2e1 * t159;
t349 = -0.2e1 * t185;
t59 = 0.1e1 / t208;
t348 = -t59 / 0.2e1;
t347 = t73 / 0.2e1;
t272 = pkin(1) * t319;
t297 = 0.2e1 / t207 * (t131 + t132) * t272;
t250 = -t297 / 0.2e1;
t289 = t154 * t207;
t90 = (t289 + (t147 * t362 - t137 + t250) * t157) * pkin(7);
t346 = -t90 / 0.2e1;
t93 = t147 * t297 / 0.2e1 + t198 * pkin(1) * t351 + (-t137 * t154 - t288) * pkin(7);
t345 = t93 / 0.2e1;
t100 = t102 ^ 2;
t244 = 0.1e1 / t144 ^ 2 * t272;
t331 = t171 / 0.2e1;
t209 = t101 ^ 2;
t99 = 0.1e1 / t209;
t39 = 0.1e1 + (((t175 * t345 + t90 * t331) * t142 + (t285 + t299) * t244) / t101 - ((t175 * t346 + t93 * t331) * t142 + (-t286 + t287) * t244) * t102 * t99) * t203 / (t100 * t99 + 0.1e1);
t338 = pkin(4) * t39;
t337 = pkin(5) * t33;
t336 = pkin(5) * t34;
t214 = pkin(5) * (t74 * t303 + (t48 * t74 - t73 * t77) * t173);
t248 = -0.2e1 * t201 * t333;
t257 = t59 * t77 / 0.2e1;
t224 = pkin(5) * (t69 + t70) * t369;
t51 = t241 * t224;
t215 = t241 * t248 + t51 * t257;
t256 = t241 * t348;
t220 = -t79 * t208 + t51 * t256;
t227 = t79 * t72 - t365;
t238 = t310 * t330;
t328 = t177 / 0.2e1;
t254 = t73 * t328;
t29 = t33 ^ 2;
t27 = 0.1e1 / (t29 * t353 + 0.1e1);
t263 = t27 * t33 * t353;
t264 = t193 * t27 / t34;
t10 = 0.1e1 - (t220 * t229 + (t227 * pkin(4) + t215) * t238 + t364 * t241) * t263 + (t215 * t254 + (((t220 - t366) * t330 + t227 * t328) * t73 + t241 * t214) * pkin(4)) * t264;
t327 = -t184 / 0.2e1;
t63 = ((t180 * t346 + t93 * t327) * t142 + (-t280 - t283) * t244) * t203;
t64 = ((t180 * t345 + t90 * t327) * t142 + (-t281 + t282) * t244) * t203;
t54 = t164 * t64 + t165 * t63;
t36 = t54 * t224;
t216 = t54 * t248 + t36 * t257;
t55 = -t164 * t63 + t165 * t64;
t221 = -t55 * t208 + t36 * t256;
t308 = t208 * t54;
t228 = t55 * t72 - t308;
t5 = 0.1e1 - (t221 * t229 + (t228 * pkin(4) + t216) * t238 + t364 * t54) * t263 + (t216 * t254 + (((-t54 * t72 + t221) * t330 + t228 * t328) * t73 + t54 * t214) * pkin(4)) * t264;
t332 = t5 * t10;
t324 = pkin(1) * t101;
t323 = pkin(1) * t102;
t321 = pkin(6) * t118;
t179 = sin(qJ(4));
t317 = t179 * t5;
t183 = cos(qJ(4));
t316 = t183 * t5;
t315 = t184 * pkin(1);
t314 = t179 * t10;
t313 = t179 * t22;
t312 = t183 * t10;
t311 = t183 * t22;
t302 = -t363 / 0.2e1;
t298 = 0.1e1 / t206 * (t129 + t130) * pkin(6) * t368;
t190 = 0.1e1 / pkin(13);
t295 = t134 * t190;
t294 = t142 * t182;
t284 = t179 * t183;
t279 = t190 * t206;
t169 = t179 ^ 2;
t170 = t183 ^ 2;
t276 = t169 + t170;
t270 = 0.2e1 * t313;
t269 = -0.2e1 * t311;
t268 = 0.2e1 * t284;
t116 = 0.1e1 / t117 ^ 2;
t267 = pkin(2) / (t116 * t118 ^ 2 + 0.1e1) * t141;
t243 = 0.1e1 / t141 ^ 2 * t363 * t322;
t251 = -t298 / 0.2e1;
t212 = ((t145 * t251 + t199 * t151 ^ 2 * t362 + (t150 * t135 - t290) * pkin(6)) * t260 + t243 * t321) / t117 * t267;
t213 = ((t145 * t368 - t150 * t206 + t151 * t251 - t291) * t260 + t117 * t243) * t116 * t267 * t321;
t52 = 0.2e1 * t212 - 0.2e1 * t213 + (0.1e1 / t134 * t298 / 0.2e1 + t354 * t266 * t362) / (-t354 * t296 + 0.1e1);
t57 = 0.2e1 * t212 - 0.2e1 * t213;
t265 = t190 * t52 * t57;
t262 = t24 * t284;
t261 = t134 * t301;
t258 = t76 * t348;
t255 = t57 / 0.2e1 + t52 / 0.2e1;
t249 = t142 * t326;
t7 = pkin(12) * t10 + t337;
t247 = t276 * t7;
t245 = (-t169 + t170) * t24;
t46 = 0.1e1 / t47 ^ 2;
t242 = pkin(9) * t195 / (t46 * t49 ^ 2 + 0.1e1) * t75;
t240 = t5 * t262;
t239 = t10 * t262;
t166 = t180 * pkin(1);
t163 = t166 + pkin(5);
t26 = t33 * t163 - t34 * t315;
t237 = 0.2e1 * t276;
t25 = t163 * t34 + t33 * t315;
t2 = -pkin(10) * t5 - t25;
t8 = -pkin(10) * t10 - t336;
t236 = t10 * t2 + t5 * t8;
t140 = -pkin(5) * t156 + t159;
t235 = 0.1e1 / t47 * t242;
t3 = pkin(12) * t5 + t26;
t234 = t2 * t24 + t22 * t3;
t233 = t22 * t7 + t24 * t8;
t232 = t119 * t244;
t231 = t122 * t244;
t223 = pkin(4) * (t47 * t74 + t73 * t76);
t222 = pkin(5) * t46 * t49 * t242;
t219 = (-t296 / 0.4e1 + t371 / 0.4e1) / t370;
t218 = pkin(4) * (-t200 * t241 * t73 + t49 * t335);
t94 = (t181 * t118 * t302 + t117 * t259) * t138;
t66 = t94 * t190 * t261 + t95 * t279 * t302;
t68 = (t94 * t206 * t301 + t95 * t261) * t190;
t217 = (-t134 * t66 / 0.2e1 - t206 * t68 / 0.2e1) * t190;
t107 = (t119 * t326 - t182 * t122 / 0.2e1) * t293;
t210 = t107 ^ 2;
t188 = pkin(16) ^ 2;
t158 = t159 ^ 2;
t104 = t108 ^ 2;
t103 = 0.1e1 / t210;
t92 = t148 * t250 + t205 * pkin(7) * t351 + (-t136 * t154 - t288) * pkin(1);
t89 = (t289 + (0.2e1 * t148 * pkin(7) - t136 + t250) * t157) * pkin(1);
t60 = (-t172 * t83 - t176 * t84) * pkin(4) + t159;
t56 = t57 ^ 2;
t40 = ((t92 * t249 + t187 * t231 + t89 * t294 / 0.2e1 + t182 * t232) / t107 - (t89 * t249 + t187 * t232 - t92 * t294 / 0.2e1 - t182 * t231) * t108 * t103) / (t103 * t104 + 0.1e1) * t197;
t38 = t172 * t338 + t323;
t37 = t176 * t338 + t324;
t21 = t24 ^ 2;
t20 = -t31 * t83 + t32 * t84;
t18 = -t31 * t38 + t32 * t37;
t17 = t31 * t37 + t32 * t38;
t16 = -pkin(10) * t22 - pkin(12) * t24 + t140;
t14 = 0.2e1 * ((t51 * t258 + (t79 * t71 - t365) * pkin(5)) * t347 + t241 * t218) * t235 - 0.2e1 * ((-t241 * t71 + t220) * t347 + t241 * t223) * t222;
t12 = 0.2e1 * ((t36 * t258 + (t55 * t71 - t308) * pkin(5)) * t347 + t54 * t218) * t235 - 0.2e1 * ((-t54 * t71 + t221) * t347 + t54 * t223) * t222 + t39;
t9 = t10 ^ 2;
t4 = t5 ^ 2;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t185 ^ 2, t181 * t349, 0, t181 ^ 2, 0, 0, t181 * t360, pkin(16) * t349, 0, t188, t153 ^ 2, -0.2e1 * t153 * t156, 0, t156 ^ 2, 0, 0, t156 * t350, t153 * t350, 0, t158, t21, t22 * t358, 0, t361, 0, 0, -0.2e1 * t140 * t22, t140 * t358, 0, t140 ^ 2, t170 * t21, -0.2e1 * t21 * t284, t24 * t269, t169 * t21, t24 * t270, t361, t16 * t269, t16 * t270, -0.2e1 * t276 * t24 * t16, t276 * t16 ^ 2, t104, t107 * t352, 0, t210, 0, 0, -0.2e1 * pkin(15) * t107, pkin(15) * t352, 0, pkin(15) ^ 2, t83 ^ 2, t84 * t357, 0, t84 ^ 2, 0, 0, t84 * t350, t159 * t357, 0, t158, t94 ^ 2, t94 * t355, 0, t95 ^ 2, 0, 0, pkin(16) * t355, t94 * t360, 0, t188, t66 ^ 2, -0.2e1 * t68 * t66, 0, t68 ^ 2, 0, 0, t68 * t356, t66 * t356, 0, t91 ^ 2, t19 ^ 2, t20 * t359, 0, t20 ^ 2, 0, 0, -0.2e1 * t60 * t20, t60 * t359, 0, t60 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, 0, -t181, 0, 0, 0, 0, 0, 0, 0, -t153, 0, t156, 0, 0, 0, (t153 * t180 - t156 * t184) * pkin(1), 0, 0, 0, t5 * t24, 0, t5 * t22, 0, 0, 0, t22 * t26 - t24 * t25, 0, t240, t5 * t245, -t5 * t313, -t240, -t5 * t311, 0, t234 * t179, t234 * t183, 0, 0, 0, 0, t108 * t40, 0, t107 * t40, 0, 0, 0, 0, 0, 0, 0, t83 * t39, 0, t39 * t84, 0, 0, 0, (-t101 * t83 + t102 * t84) * pkin(1), 0, 0, 0, t94, 0, -t95, 0, 0, 0, 0, 0, 0, 0, -t66, 0, t68, 0, 0, 0, t217, 0, 0, 0, t12 * t19, 0, t12 * t20, 0, 0, 0, t17 * t20 - t18 * t19, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t166, 0.2e1 * t315, 0, (t180 ^ 2 + t184 ^ 2) * t205, 0, 0, 0, 0, 0, t4, 0.2e1 * t25 * t5, -0.2e1 * t26 * t5, 0, t25 ^ 2 + t26 ^ 2, t169 * t4, t4 * t268, 0, t170 * t4, 0, 0, -0.2e1 * t2 * t316, 0.2e1 * t2 * t317, t5 * t3 * t237, t276 * t3 ^ 2 + t2 ^ 2, 0, 0, 0, 0, 0, t40 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 ^ 2, 0.2e1 * t39 * t324, -0.2e1 * t39 * t323, 0, (t100 + t209) * t205, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t295, t279, 0, t219, 0, 0, 0, 0, 0, t12 ^ 2, 0.2e1 * t18 * t12, -0.2e1 * t17 * t12, 0, t17 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, 0, t156, 0, 0, 0, 0, 0, 0, 0, t10 * t24, 0, t10 * t22, 0, 0, 0, (t22 * t33 - t24 * t34) * pkin(5), 0, t239, t10 * t245, -t10 * t313, -t239, -t10 * t311, 0, t233 * t179, t233 * t183, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t94, 0, -t57 * t95, 0, 0, 0, 0, 0, 0, 0, -t52 * t66, 0, t52 * t68, 0, 0, 0, t57 * t217, 0, 0, 0, t14 * t19, 0, t14 * t20, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t166, t315, 0, 0, 0, 0, 0, 0, 0, t332, t10 * t25 + t5 * t336, -t10 * t26 - t5 * t337, 0, (t25 * t34 + t26 * t33) * pkin(5), t169 * t332, t268 * t332, 0, t170 * t332, 0, 0, -t236 * t183, t236 * t179, t10 * t276 * t3 + t247 * t5, t2 * t8 + t247 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t255 * t295, t255 * t279, 0, t57 * t219, 0, 0, 0, 0, 0, t12 * t14, t18 * t14, -t17 * t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0.2e1 * t10 * t336, -0.2e1 * t10 * t337, 0, (t211 + t29) * t200, t169 * t9, t9 * t268, 0, t170 * t9, 0, 0, -0.2e1 * t8 * t312, 0.2e1 * t8 * t314, t7 * t10 * t237, t276 * t7 ^ 2 + t8 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 ^ 2, -t134 * t265, t206 * t265, 0, t56 * t219, 0, 0, 0, 0, 0, t14 ^ 2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183 * t24, 0, -t179 * t24, -t22, t183 * t16, -t179 * t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, 0, t316, 0, -t179 * t3, -t183 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, 0, t312, 0, -t179 * t7, -t183 * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MM_reg = t1;