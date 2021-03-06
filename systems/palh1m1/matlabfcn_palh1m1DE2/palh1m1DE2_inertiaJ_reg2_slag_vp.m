% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh1m1DE2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_inertiaJ_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t371 = 2 * pkin(4);
t186 = sin(qJ(3));
t187 = sin(qJ(2));
t190 = cos(qJ(3));
t191 = cos(qJ(2));
t162 = t186 * t191 + t187 * t190;
t164 = -t186 * t187 + t190 * t191;
t183 = cos(pkin(21));
t205 = pkin(5) ^ 2;
t203 = pkin(7) ^ 2;
t211 = pkin(1) ^ 2;
t192 = cos(pkin(19));
t325 = sin(pkin(19));
t163 = t187 * t192 - t191 * t325;
t322 = pkin(7) * t163;
t365 = -2 * pkin(1);
t281 = t322 * t365 + t211;
t153 = t203 + t281;
t277 = pkin(3) ^ 2 - pkin(8) ^ 2;
t146 = t153 + t277;
t156 = pkin(1) - t322;
t348 = -pkin(8) - pkin(3);
t142 = (pkin(7) - t348) * (pkin(7) + t348) + t281;
t347 = -pkin(8) + pkin(3);
t143 = (pkin(7) - t347) * (pkin(7) + t347) + t281;
t213 = sqrt(-t143 * t142);
t165 = t187 * t325 + t191 * t192;
t321 = pkin(7) * t165;
t127 = t146 * t321 + t156 * t213;
t283 = t190 * t127;
t291 = t165 * t213;
t126 = -pkin(7) * t291 + t146 * t156;
t286 = t186 * t126;
t151 = 0.1e1 / t153;
t208 = 0.1e1 / pkin(3);
t294 = t151 * t208;
t113 = (t286 / 0.2e1 + t283 / 0.2e1) * t294;
t284 = t190 * t126;
t285 = t186 * t127;
t114 = (-t284 / 0.2e1 + t285 / 0.2e1) * t294;
t174 = pkin(23) + pkin(22);
t170 = sin(t174);
t171 = cos(t174);
t88 = t113 * t171 - t114 * t170;
t337 = t88 * pkin(5);
t301 = t337 * t371 + t205;
t346 = (-pkin(9) - pkin(11));
t74 = ((pkin(4) - t346) * (pkin(4) + t346)) + t301;
t345 = (-pkin(9) + pkin(11));
t75 = ((pkin(4) - t345) * (pkin(4) + t345)) + t301;
t214 = sqrt(-t75 * t74);
t246 = t113 * t170 + t171 * t114;
t279 = pkin(9) ^ 2 - pkin(11) ^ 2;
t206 = pkin(4) ^ 2;
t84 = t206 + t301;
t78 = t84 - t279;
t368 = t246 * t78;
t86 = pkin(4) * t88 + pkin(5);
t52 = pkin(4) * t368 + t214 * t86;
t302 = t52 * t183;
t198 = 0.1e1 / pkin(11);
t82 = 0.1e1 / t84;
t308 = t198 * t82;
t179 = sin(pkin(21));
t330 = t179 / 0.2e1;
t367 = t214 * t246;
t50 = -pkin(4) * t367 + t78 * t86;
t37 = (t50 * t330 + t302 / 0.2e1) * t308;
t303 = t52 * t179;
t304 = t50 * t183;
t38 = (-t304 / 0.2e1 + t303 / 0.2e1) * t308;
t34 = atan2(t37, t38);
t31 = sin(t34);
t32 = cos(t34);
t24 = t162 * t32 + t164 * t31;
t370 = -0.2e1 * t24;
t234 = -pkin(4) * t183 * t308 / 0.2e1;
t83 = 0.1e1 / t84 ^ 2;
t338 = pkin(5) * t83;
t366 = (-0.2e1 * t86 * pkin(5) - t78) * t234 + (t303 - t304) * pkin(4) * t198 * t338;
t22 = t162 * t31 - t32 * t164;
t364 = t22 ^ 2;
t363 = -2 * pkin(16);
t178 = sin(pkin(22));
t182 = cos(pkin(22));
t200 = 0.1e1 / pkin(9);
t351 = t82 / 0.2e1;
t262 = t200 * t351;
t77 = t84 + t279;
t85 = -pkin(4) - t337;
t49 = -pkin(5) * t367 - t77 * t85;
t336 = pkin(5) * t246;
t51 = -t214 * t85 + t77 * t336;
t40 = atan2(t51 * t262, t49 * t262);
t344 = sin(t40);
t39 = cos(t40);
t29 = t178 * t39 - t182 * t344;
t30 = -t178 * t344 - t182 * t39;
t181 = cos(pkin(23));
t289 = t181 * t126;
t177 = sin(pkin(23));
t290 = t177 * t127;
t110 = (-t289 / 0.2e1 + t290 / 0.2e1) * t294;
t288 = t181 * t127;
t300 = t126 * t177;
t111 = (t288 / 0.2e1 + t300 / 0.2e1) * t294;
t97 = atan2(t111, t110);
t92 = sin(t97);
t93 = cos(t97);
t66 = t187 * t93 + t191 * t92;
t68 = -t187 * t92 + t191 * t93;
t20 = -t29 * t66 + t30 * t68;
t362 = 0.2e1 * t20;
t180 = sin(pkin(20));
t184 = cos(pkin(20));
t159 = t180 * t190 + t184 * t186;
t318 = t159 * pkin(6);
t274 = pkin(1) * t318;
t155 = 0.2e1 * t274;
t209 = pkin(2) ^ 2;
t204 = pkin(6) ^ 2;
t278 = t204 + t211;
t261 = -pkin(13) ^ 2 + t278;
t144 = t155 + t209 + t261;
t154 = -pkin(1) - t318;
t160 = t180 * t186 - t184 * t190;
t282 = t155 + t204;
t350 = -pkin(2) - pkin(13);
t140 = (pkin(1) - t350) * (pkin(1) + t350) + t282;
t349 = -pkin(2) + pkin(13);
t141 = (pkin(1) - t349) * (pkin(1) + t349) + t282;
t297 = t141 * t140;
t212 = sqrt(-t297);
t293 = t160 * t212;
t270 = pkin(6) * t293;
t123 = -t144 * t154 - t270;
t317 = t160 * pkin(6);
t124 = t144 * t317 - t154 * t212;
t150 = t155 + t278;
t147 = 0.1e1 / t150;
t210 = 0.1e1 / pkin(2);
t326 = t210 / 0.2e1;
t257 = t147 * t326;
t107 = atan2(t124 * t257, t123 * t257);
t104 = sin(t107);
t105 = cos(t107);
t79 = t104 * t191 + t105 * t187;
t76 = pkin(2) * t79 - pkin(16);
t361 = -0.2e1 * t76;
t81 = -t104 * t187 + t105 * t191;
t360 = -0.2e1 * t81;
t145 = t153 - t277;
t157 = pkin(1) * t163 - pkin(7);
t125 = -pkin(1) * t291 - t145 * t157;
t128 = pkin(1) * t165 * t145 - t157 * t213;
t188 = sin(pkin(18));
t202 = 0.1e1 / pkin(8);
t295 = t151 * t202;
t193 = cos(pkin(18));
t327 = t193 / 0.2e1;
t115 = (t125 * t327 - t188 * t128 / 0.2e1) * t295;
t116 = (t128 * t327 + t125 * t188 / 0.2e1) * t295;
t99 = atan2(t116, t115);
t95 = sin(t99);
t359 = 0.2e1 * t95;
t358 = 0.1e1 / t38 ^ 2;
t223 = t209 - t261 - 0.2e1 * t274;
t357 = 0.1e1 / t223 ^ 2;
t167 = t187 * pkin(1) - pkin(16);
t149 = -pkin(5) * t164 + t167;
t356 = 0.2e1 * t149;
t355 = -0.2e1 * t165 ^ 2;
t354 = 0.2e1 * t167;
t353 = -0.2e1 * t191;
t65 = 0.1e1 / t214;
t352 = -t65 / 0.2e1;
t343 = pkin(1) * t92;
t342 = pkin(1) * t93;
t276 = pkin(1) * t321;
t298 = 0.2e1 / t213 * (t142 + t143) * t276;
t258 = -t298 / 0.2e1;
t292 = t163 * t213;
t101 = (t292 + (t156 * t365 - t146 + t258) * t165) * pkin(7);
t103 = t156 * t298 / 0.2e1 + t203 * pkin(1) * t355 + (-t163 * t146 - t291) * pkin(7);
t109 = 0.1e1 / t110 ^ 2;
t250 = 0.1e1 / t153 ^ 2 * t276;
t331 = t177 / 0.2e1;
t333 = t103 / 0.2e1;
t334 = -t101 / 0.2e1;
t45 = 0.1e1 + (((t101 * t331 + t181 * t333) * t151 + (t288 + t300) * t250) / t110 - ((t103 * t331 + t181 * t334) * t151 + (-t289 + t290) * t250) * t111 * t109) * t208 / (t109 * t111 ^ 2 + 0.1e1);
t341 = pkin(4) * t45;
t340 = pkin(5) * t31;
t339 = pkin(5) * t32;
t219 = pkin(5) * (t83 * t302 + (t50 * t83 - t82 * t86) * t179);
t254 = -0.2e1 * t206 * t336;
t265 = t65 * t86 / 0.2e1;
t229 = pkin(5) * (t74 + t75) * t371;
t53 = t246 * t229;
t220 = t246 * t254 + t53 * t265;
t264 = t246 * t352;
t224 = -t88 * t214 + t53 * t264;
t232 = t88 * t78 - t367;
t243 = t308 * t330;
t329 = t183 / 0.2e1;
t263 = t82 * t329;
t33 = 0.1e1 / (t37 ^ 2 * t358 + 0.1e1);
t268 = t33 * t37 * t358;
t269 = t198 * t33 / t38;
t10 = 0.1e1 - (t224 * t234 + (t232 * pkin(4) + t220) * t243 + t366 * t246) * t268 + (t220 * t263 + (((t224 - t368) * t330 + t232 * t329) * t82 + t246 * t219) * pkin(4)) * t269;
t328 = -t190 / 0.2e1;
t71 = ((t103 * t328 + t186 * t334) * t151 + (-t283 - t286) * t250) * t208;
t72 = ((t101 * t328 + t186 * t333) * t151 + (-t284 + t285) * t250) * t208;
t57 = t170 * t72 + t171 * t71;
t42 = t57 * t229;
t221 = t57 * t254 + t42 * t265;
t58 = -t170 * t71 + t171 * t72;
t225 = -t58 * t214 + t42 * t264;
t307 = t214 * t57;
t233 = t58 * t78 - t307;
t5 = 0.1e1 - (t225 * t234 + (t233 * pkin(4) + t221) * t243 + t366 * t57) * t268 + (t221 * t263 + (((-t57 * t78 + t225) * t330 + t233 * t329) * t82 + t57 * t219) * pkin(4)) * t269;
t335 = t5 * t10;
t332 = t147 / 0.2e1;
t324 = pkin(1) * t160;
t323 = pkin(6) * t124;
t255 = 0.1e1 / pkin(13) * t326;
t133 = atan2(t212 * t255, t223 * t255);
t131 = sin(t133);
t320 = t131 * pkin(2);
t132 = cos(t133);
t319 = t132 * pkin(2);
t185 = sin(qJ(4));
t316 = t185 * t5;
t189 = cos(qJ(4));
t315 = t189 * t5;
t314 = t190 * pkin(1);
t148 = 0.1e1 / t150 ^ 2;
t122 = 0.1e1 / t123 ^ 2;
t245 = pkin(2) / (t122 * t124 ^ 2 + 0.1e1) * t150 * t210;
t299 = 0.2e1 / t212 * (t140 + t141) * pkin(1) * t317;
t259 = -t299 / 0.2e1;
t217 = ((t154 * t259 + (t159 * t144 - t293) * pkin(6)) * t332 + (-t147 * t160 * t204 + t148 * t323) * t324) / t123 * t245;
t218 = ((-t159 * t212 + (t259 - t144) * t160) * t332 + (t123 * t148 + t147 * t154) * t324) * t122 * t245 * t323;
t55 = 0.2e1 * t217 - 0.2e1 * t218 + (0.1e1 / t223 * t299 / 0.2e1 + t357 * t270 * t365) / (-t357 * t297 + 0.1e1);
t60 = 0.2e1 * t217 - 0.2e1 * t218;
t313 = t55 + t60;
t312 = t185 * t10;
t311 = t185 * t22;
t310 = t189 * t10;
t309 = t189 * t22;
t296 = t151 * t188;
t287 = t185 * t189;
t175 = t185 ^ 2;
t176 = t189 ^ 2;
t280 = t175 + t176;
t275 = pkin(2) * t55 * t60;
t273 = -0.2e1 * t311;
t272 = 0.2e1 * t309;
t271 = 0.2e1 * t287;
t267 = t24 * t287;
t266 = t85 * t352;
t256 = t151 * t327;
t7 = pkin(12) * t10 + t340;
t253 = t280 * t7;
t251 = (-t175 + t176) * t24;
t48 = 0.1e1 / t49 ^ 2;
t249 = pkin(9) * t200 / (t48 * t51 ^ 2 + 0.1e1) * t84;
t248 = (t131 ^ 2 + t132 ^ 2) * t209;
t247 = t5 * t267;
t244 = t10 * t267;
t172 = t186 * pkin(1);
t169 = t172 + pkin(5);
t26 = t31 * t169 - t32 * t314;
t242 = 0.2e1 * t280;
t25 = t169 * t32 + t31 * t314;
t2 = -pkin(10) * t5 - t25;
t8 = -pkin(10) * t10 - t339;
t241 = t10 * t2 + t5 * t8;
t240 = 0.1e1 / t49 * t249;
t3 = pkin(12) * t5 + t26;
t239 = t2 * t24 - t22 * t3;
t238 = -t22 * t7 + t24 * t8;
t237 = t125 * t250;
t236 = t128 * t250;
t228 = pkin(4) * (t49 * t83 + t82 * t85);
t227 = pkin(5) * t48 * t51 * t249;
t61 = t131 * t81 + t132 * t79;
t62 = -t131 * t79 + t132 * t81;
t226 = (-t131 * t61 - t132 * t62) * pkin(2);
t222 = pkin(4) * (-t205 * t246 * t82 + t51 * t338);
t194 = pkin(16) ^ 2;
t166 = t167 ^ 2;
t112 = 0.1e1 / t115 ^ 2;
t102 = t157 * t258 + t211 * pkin(7) * t355 + (-t163 * t145 - t291) * pkin(1);
t100 = (t292 + (0.2e1 * t157 * pkin(7) - t145 + t258) * t165) * pkin(1);
t96 = cos(t99);
t59 = t60 ^ 2;
t54 = (-t178 * t68 + t182 * t66) * pkin(4) + t167;
t46 = ((t102 * t256 + t193 * t236 + t100 * t296 / 0.2e1 + t188 * t237) / t115 - (t100 * t256 + t193 * t237 - t102 * t296 / 0.2e1 - t188 * t236) * t116 * t112) / (t112 * t116 ^ 2 + 0.1e1) * t202;
t44 = t182 * t341 + t342;
t43 = t178 * t341 + t343;
t21 = t24 ^ 2;
t19 = -t29 * t68 - t30 * t66;
t18 = t29 * t44 + t30 * t43;
t17 = -t29 * t43 + t30 * t44;
t16 = 0.2e1 * ((t53 * t266 + (t88 * t77 - t367) * pkin(5)) * t351 + t246 * t222) * t240 - 0.2e1 * ((-t246 * t77 + t224) * t351 + t246 * t228) * t227;
t14 = pkin(10) * t22 - pkin(12) * t24 + t149;
t12 = 0.2e1 * ((t42 * t266 + (t58 * t77 - t307) * pkin(5)) * t351 + t57 * t222) * t240 - 0.2e1 * ((-t57 * t77 + t225) * t351 + t57 * t228) * t227 + t45;
t9 = t10 ^ 2;
t4 = t5 ^ 2;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t191 ^ 2, t187 * t353, 0, t187 ^ 2, 0, 0, t187 * t363, pkin(16) * t353, 0, t194, t162 ^ 2, 0.2e1 * t162 * t164, 0, t164 ^ 2, 0, 0, -0.2e1 * t167 * t164, t162 * t354, 0, t166, t21, t22 * t370, 0, t364, 0, 0, t22 * t356, t24 * t356, 0, t149 ^ 2, t176 * t21, -0.2e1 * t21 * t287, t24 * t272, t175 * t21, t24 * t273, t364, t14 * t272, t14 * t273, t280 * t14 * t370, t280 * t14 ^ 2, t95 ^ 2, t96 * t359, 0, t96 ^ 2, 0, 0, -0.2e1 * pkin(15) * t96, pkin(15) * t359, 0, pkin(15) ^ 2, t68 ^ 2, -0.2e1 * t66 * t68, 0, t66 ^ 2, 0, 0, t66 * t354, t68 * t354, 0, t166, t81 ^ 2, t79 * t360, 0, t79 ^ 2, 0, 0, t79 * t363, pkin(16) * t360, 0, t194, t62 ^ 2, -0.2e1 * t61 * t62, 0, t61 ^ 2, 0, 0, t61 * t361, t62 * t361, 0, t76 ^ 2, t20 ^ 2, t19 * t362, 0, t19 ^ 2, 0, 0, -0.2e1 * t54 * t19, t54 * t362, 0, t54 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, -t187, 0, 0, 0, 0, 0, 0, 0, t162, 0, t164, 0, 0, 0, (-t162 * t186 - t164 * t190) * pkin(1), 0, 0, 0, t24 * t5, 0, -t22 * t5, 0, 0, 0, -t22 * t26 - t24 * t25, 0, t247, t5 * t251, t5 * t311, -t247, t5 * t309, 0, t239 * t185, t239 * t189, 0, 0, 0, 0, t46 * t95, 0, t46 * t96, 0, 0, 0, 0, 0, 0, 0, t45 * t68, 0, -t45 * t66, 0, 0, 0, (-t66 * t92 - t68 * t93) * pkin(1), 0, 0, 0, t81, 0, -t79, 0, 0, 0, 0, 0, 0, 0, -t62, 0, t61, 0, 0, 0, t226, 0, 0, 0, t12 * t20, 0, t12 * t19, 0, 0, 0, -t17 * t20 + t18 * t19, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t172, 0.2e1 * t314, 0, (t186 ^ 2 + t190 ^ 2) * t211, 0, 0, 0, 0, 0, t4, 0.2e1 * t25 * t5, -0.2e1 * t26 * t5, 0, t25 ^ 2 + t26 ^ 2, t175 * t4, t4 * t271, 0, t176 * t4, 0, 0, -0.2e1 * t2 * t315, 0.2e1 * t2 * t316, t5 * t3 * t242, t280 * t3 ^ 2 + t2 ^ 2, 0, 0, 0, 0, 0, t46 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 ^ 2, 0.2e1 * t45 * t342, -0.2e1 * t45 * t343, 0, (t92 ^ 2 + t93 ^ 2) * t211, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -0.2e1 * t319, 0.2e1 * t320, 0, t248, 0, 0, 0, 0, 0, t12 ^ 2, 0.2e1 * t17 * t12, -0.2e1 * t18 * t12, 0, t17 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, t164, 0, 0, 0, 0, 0, 0, 0, t24 * t10, 0, -t22 * t10, 0, 0, 0, (-t22 * t31 - t24 * t32) * pkin(5), 0, t244, t10 * t251, t10 * t311, -t244, t10 * t309, 0, t238 * t185, t238 * t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t60, 0, -t60 * t79, 0, 0, 0, 0, 0, 0, 0, -t55 * t62, 0, t55 * t61, 0, 0, 0, t60 * t226, 0, 0, 0, t16 * t20, 0, t16 * t19, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t172, t314, 0, 0, 0, 0, 0, 0, 0, t335, t10 * t25 + t5 * t339, -t10 * t26 - t5 * t340, 0, (t25 * t32 + t26 * t31) * pkin(5), t175 * t335, t271 * t335, 0, t176 * t335, 0, 0, -t241 * t189, t241 * t185, t10 * t280 * t3 + t253 * t5, t2 * t8 + t253 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t313 * t319, t313 * t320, 0, t60 * t248, 0, 0, 0, 0, 0, t12 * t16, t17 * t16, -t18 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0.2e1 * t10 * t339, -0.2e1 * t10 * t340, 0, (t31 ^ 2 + t32 ^ 2) * t205, t175 * t9, t9 * t271, 0, t176 * t9, 0, 0, -0.2e1 * t8 * t310, 0.2e1 * t8 * t312, t7 * t10 * t242, t280 * t7 ^ 2 + t8 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 ^ 2, -0.2e1 * t132 * t275, 0.2e1 * t131 * t275, 0, t59 * t248, 0, 0, 0, 0, 0, t16 ^ 2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189 * t24, 0, -t185 * t24, t22, t189 * t14, -t185 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, 0, t315, 0, -t185 * t3, -t189 * t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, 0, t310, 0, -t185 * t7, -t189 * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MM_reg = t1;
