% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% tau_reg [2x(2*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnTE_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_reg2_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:19:54
% EndTime: 2020-04-12 19:20:11
% DurationCPUTime: 8.89s
% Computational Cost: add. (74527->350), mult. (105423->819), div. (3238->18), fcn. (28337->6), ass. (0->311)
t133 = 0.1e1 / pkin(4);
t139 = pkin(2) ^ 2;
t140 = pkin(1) ^ 2;
t127 = cos(qJ(2));
t334 = pkin(2) * t127;
t280 = -0.2e1 * pkin(1) * t334 + t140;
t119 = t139 + t280;
t116 = 0.1e1 / t119 ^ 2;
t125 = sin(qJ(2));
t275 = qJD(2) * t125;
t244 = pkin(2) * t275;
t219 = pkin(1) * t244;
t196 = t116 * t219;
t383 = pkin(3) ^ 2;
t384 = pkin(4) ^ 2;
t278 = t383 - t384;
t114 = t119 - t278;
t120 = pkin(1) - t334;
t355 = -pkin(3) - pkin(4);
t111 = (pkin(2) - t355) * (pkin(2) + t355) + t280;
t354 = pkin(4) - pkin(3);
t112 = (pkin(2) - t354) * (pkin(2) + t354) + t280;
t302 = t111 * t112;
t141 = sqrt(-t302);
t293 = t125 * t141;
t95 = -pkin(2) * t293 + t114 * t120;
t179 = t95 * t196;
t115 = 0.1e1 / t119;
t351 = -t115 / 0.2e1;
t212 = (-t111 - t112) * pkin(1) * pkin(2);
t101 = t125 * t212;
t100 = qJD(2) * t101;
t105 = 0.1e1 / t141;
t307 = t100 * t105;
t237 = t125 * t307;
t264 = 0.2e1 * t120 * pkin(1);
t288 = t127 * t141;
t42 = (-t237 + (-t288 + (t114 + t264) * t125) * qJD(2)) * pkin(2);
t27 = (t42 * t351 + t179) * t133;
t335 = pkin(2) * t125;
t109 = t114 * t335;
t96 = t120 * t141 + t109;
t178 = t96 * t196;
t350 = t115 / 0.2e1;
t123 = t125 ^ 2;
t296 = t123 * t139;
t247 = pkin(1) * t296;
t216 = qJD(2) * t247;
t234 = t141 * t275;
t274 = qJD(2) * t127;
t243 = pkin(2) * t274;
t283 = pkin(2) * t234 + t114 * t243;
t305 = t105 * t120;
t43 = t100 * t305 + 0.2e1 * t216 + t283;
t29 = (t43 * t350 - t178) * t133;
t90 = 0.1e1 / t95 ^ 2;
t322 = t90 * t96;
t92 = t96 ^ 2;
t57 = t90 * t92 + 0.1e1;
t53 = 0.1e1 / t57;
t89 = 0.1e1 / t95;
t163 = t53 * (t27 * t322 + t29 * t89);
t113 = t119 + t278;
t339 = pkin(1) * t125;
t110 = t113 * t339;
t121 = pkin(1) * t127 - pkin(2);
t97 = -t121 * t141 + t110;
t357 = -t97 / 0.2e1;
t370 = -0.2e1 * t121;
t266 = pkin(2) * t370;
t306 = t105 * t101;
t45 = t110 + (-t288 + (t266 - t306) * t125) * pkin(1);
t363 = t45 / 0.2e1;
t241 = t357 + t363;
t94 = -pkin(1) * t293 - t113 * t121;
t360 = -t94 / 0.2e1;
t295 = t123 * t140;
t304 = t105 * t121;
t49 = -t101 * t304 + 0.2e1 * pkin(2) * t295 + (t113 * t127 + t293) * pkin(1);
t242 = t360 - t49 / 0.2e1;
t152 = (t241 * t125 - t242 * t127) * t115;
t336 = pkin(2) * t116;
t260 = pkin(1) * t336;
t294 = t125 * t127;
t157 = (t123 * t94 + t97 * t294) * t260;
t388 = t157 - t152;
t299 = t116 * t125;
t248 = pkin(2) * t299;
t226 = pkin(1) * t248;
t199 = t97 * t226;
t315 = t116 * t94;
t252 = pkin(2) * t315;
t136 = 0.1e1 / pkin(3);
t289 = t127 * t136;
t297 = t123 * t136;
t387 = -pkin(1) * t252 * t297 - t199 * t289;
t86 = 0.1e1 / t94 ^ 2;
t93 = t97 ^ 2;
t58 = t86 * t93 + 0.1e1;
t55 = 0.1e1 / t58;
t386 = 0.2e1 * t55;
t382 = -0.4e1 * t95;
t258 = 0.2e1 * t116;
t381 = 0.6e1 * t125;
t332 = pkin(4) * t119;
t265 = 0.2e1 * t332;
t380 = 2 * qJDD(2);
t256 = t53 * t332;
t222 = t90 * t256;
t379 = t133 * t96 * t222;
t308 = t127 * t97;
t349 = -t125 / 0.2e1;
t170 = -t308 / 0.2e1 + t94 * t349;
t159 = t170 * t115;
t378 = t387 * qJD(1);
t273 = qJD(1) * qJD(2);
t229 = t125 * t273;
t211 = pkin(2) * t229;
t268 = t127 * qJDD(1);
t377 = pkin(2) * t268 - t211;
t130 = qJD(2) ^ 2;
t290 = t127 * t130;
t376 = qJDD(2) * t125 + t290;
t181 = -0.2e1 * t196;
t227 = t113 + t266;
t41 = (-t237 + (t227 * t125 - t288) * qJD(2)) * pkin(1);
t28 = (t115 * t41 + t94 * t181) * t136;
t235 = qJD(2) * t295;
t215 = pkin(2) * t235;
t282 = (t113 * t274 + t234) * pkin(1);
t44 = -t100 * t304 + 0.2e1 * t215 + t282;
t30 = (t115 * t44 + t97 * t181) * t136;
t324 = t86 * t97;
t85 = 0.1e1 / t94;
t375 = -t28 * t324 + t30 * t85;
t373 = -0.2e1 * t90;
t372 = 0.4e1 * t96;
t371 = 0.2e1 * t97;
t369 = 0.2e1 * qJD(2);
t367 = -t41 / 0.2e1;
t366 = -t42 / 0.2e1;
t365 = t43 / 0.2e1;
t364 = -t44 / 0.2e1;
t46 = t109 + (-t288 + (t264 - t306) * t125) * pkin(2);
t362 = t46 / 0.2e1;
t50 = t101 * t305 + 0.2e1 * t247 + (t114 * t127 + t293) * pkin(2);
t361 = -t50 / 0.2e1;
t359 = t94 / 0.2e1;
t358 = -t95 / 0.2e1;
t356 = t97 / 0.2e1;
t259 = pkin(1) * t335;
t190 = t259 * t386;
t84 = t94 ^ 2;
t87 = t85 / t84;
t262 = t87 * t371;
t319 = t97 * t49;
t323 = t87 * t93;
t117 = t115 * t116;
t298 = t117 * t139;
t185 = t235 * t298;
t175 = 0.8e1 * t185;
t303 = 0.4e1 * t105 / t302;
t210 = t100 * t101 * t303;
t182 = -t210 / 0.4e1;
t239 = t115 * t307;
t267 = -0.2e1 * t336;
t311 = t125 * t49;
t218 = -0.4e1 * t139 * t295;
t98 = (t127 * t212 + t218) * qJD(2);
t331 = ((t140 * t243 * t381 + t121 * t182 - t98 * t304) * t115 + t97 * t175 + ((t44 * t267 + t239) * t125 + ((t288 + (-t113 + t306) * t125) * t115 + (-t308 - t311) * pkin(2) * t258) * qJD(2)) * pkin(1)) * t136 * t85;
t56 = 0.1e1 / t58 ^ 2;
t150 = (t125 * t182 + 0.2e1 * (-t100 * t127 + t98 * t349) * t105) * t115;
t333 = pkin(3) * t119;
t223 = t55 * t86 * t333;
t203 = t97 * t223;
t291 = t127 * t115;
t309 = t127 * t94;
t312 = t125 * t45;
t338 = pkin(1) * t136;
t8 = (((0.4e1 * t215 + t282) * t115 + t94 * t175) * t136 + (t150 + (-0.2e1 * t41 * t299 + (t291 * t370 + (-t309 - t312) * t258) * qJD(2)) * pkin(2)) * t338) * t203;
t2 = -t8 + (t375 * t190 + (-0.2e1 * t375 * t56 * (t86 * t319 - t45 * t323) + (-t30 * t45 * t86 + t331 + (t45 * t262 - t49 * t86) * t28) * t55) * t119) * pkin(3);
t32 = (t45 * t115 - 0.2e1 * t94 * t226) * t136;
t34 = (t49 * t115 - 0.2e1 * t199) * t136;
t176 = t32 * t324 - t34 * t85;
t320 = t97 * t44;
t330 = (t86 * t320 - t41 * t323) * t56;
t4 = -t8 + (-t176 * qJD(2) * t190 + (0.2e1 * t176 * t330 + (-t34 * t41 * t86 + t331 + (t41 * t262 - t44 * t86) * t32) * t55) * t119) * pkin(3);
t353 = t2 - t4;
t197 = t133 * t89 * t256;
t250 = t139 * t339;
t198 = 0.6e1 * t127 * t250;
t340 = pkin(1) * t116;
t352 = -0.2e1 * ((0.4e1 * t216 + t283) * t351 + t185 * t382 + (-t150 / 0.2e1 + (t42 * t299 + (-t120 * t291 + (t125 * t46 + t127 * t95) * t116) * qJD(2)) * pkin(1)) * pkin(2)) * t379 - 0.2e1 * ((t120 * t210 / 0.4e1 + t98 * t305 + qJD(2) * t198) * t350 + t185 * t372 + ((t239 / 0.2e1 - t43 * t340) * t125 + ((t288 + (-t114 + t306) * t125) * t350 + (-t125 * t50 - t127 * t96) * t340) * qJD(2)) * pkin(2)) * t197;
t348 = t125 / 0.2e1;
t347 = -t127 / 0.2e1;
t346 = t127 / 0.2e1;
t126 = sin(qJ(1));
t345 = g(1) * t126;
t128 = cos(qJ(1));
t344 = g(1) * t128;
t343 = g(2) * t126;
t342 = g(2) * t128;
t341 = pkin(1) * t115;
t337 = pkin(2) * t115;
t51 = t136 * t159;
t327 = t51 * t94;
t310 = t125 * t97;
t171 = t309 / 0.2e1 - t310 / 0.2e1;
t160 = t115 * t171;
t52 = t136 * t160;
t326 = t52 * t97;
t54 = 0.1e1 / t57 ^ 2;
t325 = t54 * t89;
t88 = t95 ^ 2;
t91 = t89 / t88;
t321 = t91 * t92;
t300 = t115 * t136;
t231 = t300 / 0.2e1;
t205 = t127 * t231;
t232 = -t300 / 0.2e1;
t318 = (t205 * t94 + t232 * t310) * t126;
t206 = t125 * t231;
t317 = (t205 * t97 + t206 * t94) * t128;
t316 = pkin(1) * qJD(1);
t314 = t116 * t95;
t313 = t116 * t96;
t301 = t115 * t121;
t292 = t126 * t136;
t287 = t128 * t136;
t286 = t130 * t140;
t131 = qJD(1) ^ 2;
t134 = 0.1e1 / t384;
t285 = t131 * t134;
t284 = t141 * t130;
t124 = t127 ^ 2;
t279 = t123 - t124;
t277 = qJD(1) * t127;
t276 = qJD(2) * t116;
t272 = qJDD(1) * t124;
t271 = qJDD(2) * t113;
t270 = qJDD(2) * t114;
t263 = 0.2e1 * t53 * t90;
t261 = t91 * t372;
t257 = t85 * t333;
t21 = -t42 * t321 + t43 * t322;
t255 = t21 * t325;
t23 = -t46 * t321 + t50 * t322;
t230 = 0.4e1 * t54 * t322;
t1 = (-0.4e1 * t163 * t259 + ((0.4e1 * t23 * t325 + t46 * t263) * t29 + (t23 * t230 + (t46 * t261 + t50 * t373) * t53) * t27) * t119) * pkin(4) + t352;
t201 = t95 * t226;
t31 = (t46 * t351 + t201) * t133;
t200 = t96 * t226;
t33 = (t50 * t350 - t200) * t133;
t162 = t53 * (-t31 * t322 - t33 * t89);
t3 = (0.4e1 * t162 * t219 + ((t42 * t263 + 0.4e1 * t255) * t33 + (t21 * t230 + (t42 * t261 + t43 * t373) * t53) * t31) * t119) * pkin(4) + t352;
t251 = t3 / 0.2e1 - t1 / 0.2e1;
t249 = pkin(2) * t300;
t246 = qJD(1) * t335;
t245 = pkin(2) * t277;
t238 = t130 * t296;
t236 = t131 * t294;
t233 = t131 * t350;
t228 = qJDD(1) * t358;
t225 = t117 * t259;
t224 = t55 * t257;
t221 = qJDD(1) * t116 / 0.4e1;
t12 = t163 * t265;
t13 = t162 * t265;
t220 = qJD(2) * t13 + t12;
t217 = t130 * t250;
t214 = pkin(2) * t232;
t213 = pkin(2) * t231;
t209 = t94 * t232;
t204 = t127 * t229;
t195 = t117 * t219;
t194 = t225 / 0.2e1;
t118 = t343 + t344;
t193 = -t342 + t345;
t192 = t116 * t217;
t191 = t131 * t140 * t248;
t14 = -t32 * t203 + t34 * t224 + 0.1e1;
t189 = qJDD(2) * t14 * t249;
t188 = qJD(2) * t214;
t187 = qJD(2) * t213;
t186 = -t115 * t120 + t314;
t183 = t96 * t95 * t225;
t173 = t344 / 0.2e1 + t343 / 0.2e1;
t167 = t100 ^ 2 * t303 / 0.4e1 + t105 * (t130 * t218 + t376 * t212);
t166 = t136 * t14 * t192;
t165 = t140 * t211 * t258;
t164 = qJDD(2) * t141 + t307 * t369;
t158 = -qJDD(1) * pkin(1) - t345 / 0.2e1 + t342 / 0.2e1;
t156 = (-t123 * t97 + t94 * t294) * t260;
t155 = -g(3) * t127 + t118 * t125;
t154 = -t167 + t284;
t153 = -t114 * t130 + t164;
t151 = t171 * qJD(2) + t44 * t346 + t41 * t348;
t149 = t136 * (t156 + (-t242 * t125 - t241 * t127) * t115);
t148 = (t41 * t347 + t44 * t348) * t115 + (-t159 + t156) * qJD(2);
t48 = t97 * qJD(1) * t206 + t209 * t277;
t47 = qJD(1) * t51;
t20 = qJD(1) * t149;
t19 = t136 * qJD(1) * t152 + t378;
t18 = t148 * t136;
t17 = t387 * qJD(2) + t151 * t300;
t16 = (t148 * qJD(1) - qJDD(1) * t160) * t136;
t15 = (t151 * qJD(1) - t170 * qJDD(1)) * t300 + t378 * qJD(2);
t11 = -t28 * t203 + t30 * t224 + qJD(2);
t6 = qJDD(2) + ((-t167 * t301 + (0.8e1 * t97 * t117 * t238 + (t123 * t380 + t290 * t381) * t337) * t140) * t136 + (((t271 + t284) * t115 + t97 * t130 * t267) * t127 + ((-t113 * t130 + t164) * t115 + (-0.4e1 * t44 * qJD(2) - 0.2e1 * t97 * qJDD(2)) * t336) * t125) * t338) * t224 - 0.2e1 * t30 * t257 * t330 - ((0.8e1 * t94 * t298 + 0.4e1 * t337) * t286 * t297 + ((-t164 * t115 + (t227 * t115 - 0.2e1 * t252) * t130) * t289 + ((t154 + t271) * t115 + (-0.4e1 * t41 * t276 + (-t301 - t315) * t380) * pkin(2)) * t136 * t125) * pkin(1)) * t203 + (-t28 * t44 - t30 * t41) * t223 + (t87 * t55 * t41 + t86 * t330) * t28 * t333 * t371 + t375 * pkin(3) * t219 * t386;
t5 = -0.2e1 * ((t167 * t120 + t130 * t198) * t350 + (t117 * t286 * t372 + qJDD(2) * t341) * t296 + (((t270 + t284) * t127 + t153 * t125) * t350 + (-t96 * t290 + (-0.2e1 * qJD(2) * t43 - qJDD(2) * t96) * t125) * t340) * pkin(2)) * t197 - 0.2e1 * ((t117 * t140 * t382 - 0.2e1 * t341) * t238 + ((t186 * t130 * pkin(1) + t153 * t350) * t127 + ((t154 + t270) * t351 + (t186 * qJDD(2) + 0.2e1 * t42 * t276) * pkin(1)) * t125) * pkin(2)) * t379 + 0.2e1 * (-t27 * t43 + t29 * t42) * t222 - 0.4e1 * pkin(4) * t219 * t163 + 0.2e1 * (t29 * t255 + (t90 * t54 * t21 + t91 * t53 * t42) * t27 * t96) * t265;
t7 = [0, 0, 0, 0, 0, qJDD(1), t193, t118, 0, 0, qJDD(1) * t123 + 0.2e1 * t204, 0.2e1 * t125 * t268 - 0.2e1 * t279 * t273, t376, -0.2e1 * t204 + t272, qJDD(2) * t127 - t125 * t130, 0, t193 * t127, -t193 * t125, -t118, 0, -t15 * t51 - t17 * t47, t15 * t52 + t16 * t51 - t17 * t48 + t18 * t47, -t11 * t17 + t51 * t6, -t16 * t52 + t18 * t48, t11 * t18 - t52 * t6, 0, -g(1) * t318 + t16 * t334 + t18 * t245 - t48 * t244 - (-t342 + t377) * t52, -g(1) * t292 * t159 - g(2) * t317 + t15 * t334 + t17 * t245 + t244 * t47 - t377 * t51, ((-t326 - t327) * t192 + ((t326 / 0.2e1 + t327 / 0.2e1) * qJDD(2) + (t17 * t360 + t18 * t357 - t52 * t364 - t51 * t367) * qJD(2)) * t337) * t136 - t118, t139 * t272 + (pkin(2) * t193 - 0.2e1 * t139 * t229) * t127, (t92 * t221 + (-t92 * t195 + t313 * t365) * qJD(1)) * t134, (t228 * t313 + (t183 * t369 + (t43 * t358 + t96 * t366) * t116) * qJD(1)) * t134, (t12 * t178 + (-t12 * t365 + t96 * t5 / 0.2e1) * t115) * t133, (t88 * t221 + (t42 * t314 / 0.2e1 - t88 * t195) * qJD(1)) * t134, (-t12 * t179 + (-t12 * t366 + t5 * t358) * t115) * t133, 0, (t95 * t165 + (t158 * t95 - t316 * t42) * t115) * t133, (t96 * t165 + (t158 * t96 - t316 * t43) * t115) * t133, -t118, pkin(1) * t193 + qJDD(1) * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t279 * t131, t125 * qJDD(1), t236, t268, qJDD(2), t155, g(3) * t125 + t118 * t127, 0, 0, t47 * t19, t19 * t48 - t20 * t47, t11 * t19 - t14 * t15 - t353 * t47, -t48 * t20, -t11 * t20 + t14 * t16 - t353 * t48, -t353 * t11 + t14 * t6, pkin(2) * t6 * t209 + t41 * t14 * t188 + t189 * t360 + t48 * t246 - t20 * t245 - g(1) * (((t45 * t347 + t311 / 0.2e1) * t115 + t156) * t287 + t317) - t149 * t343 - g(3) * t388 * t136 + (t45 * t187 + t41 * t214) * t11 + (t2 * t187 + t4 * t188 + t166) * t94, t44 * t14 * t187 + t189 * t356 - t47 * t246 - t19 * t245 + g(1) * t388 * t287 - g(2) * (((t49 * t346 + t312 / 0.2e1) * t115 - t157) * t292 + t318) - g(3) * t149 + (t49 * t188 + t44 * t213) * t11 + (t4 * t187 + t2 * t188 + t6 * t213 - t166) * t97, (t48 * t364 + t16 * t357 - t47 * t367 + t15 * t360 + (t49 * t48 / 0.2e1 + t20 * t356 - t47 * t363 + t19 * t359) * qJD(2)) * t249, t139 * t236 + t155 * pkin(2) + ((-t93 / 0.2e1 - t84 / 0.2e1) * pkin(2) * t117 * t217 + ((-t319 / 0.4e1 - t94 * t45 / 0.4e1) * t130 + (t93 / 0.4e1 + t84 / 0.4e1) * qJDD(2) + (t320 / 0.2e1 + t41 * t359) * qJD(2)) * t139 * t116) / t383, (-t50 * t313 / 0.4e1 + t92 * t194) * t285, (-t183 + (t50 * t95 / 0.4e1 + t96 * t46 / 0.4e1) * t116) * t285, (t96 * qJDD(1) * t13 * t350 + (-t220 * t200 + (-t12 * t361 + t13 * t365 + t251 * t96) * t115) * qJD(1)) * t133, (-t46 * t314 / 0.4e1 + t88 * t194) * t285, (t115 * t13 * t228 + (t220 * t201 + (-t12 * t362 + t13 * t366 - t251 * t95) * t115) * qJD(1)) * t133, t13 * t5 - (-t1 + t3) * t12, (-t95 * t191 + (g(3) * t361 + t173 * t46) * t115 + (t46 * t233 + (g(3) * t96 - t118 * t95) * t248) * pkin(1)) * t133, (-t96 * t191 + (g(3) * t362 + t173 * t50) * t115 + (t50 * t233 + (-g(3) * t95 - t118 * t96) * t248) * pkin(1)) * t133, 0, 0;];
tau_reg = t7;