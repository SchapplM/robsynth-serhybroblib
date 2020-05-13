% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% picker2Dm1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% tau_reg [2x(12*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = picker2Dm1IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:55:03
% EndTime: 2020-05-11 05:55:08
% DurationCPUTime: 4.91s
% Computational Cost: add. (7132->480), mult. (8791->603), div. (257->10), fcn. (5393->70), ass. (0->320)
t257 = sin(qJ(2));
t348 = qJDD(1) * t257;
t265 = cos(qJ(2));
t359 = qJD(2) * t265;
t283 = qJD(1) * t359 + t348;
t421 = t283 * pkin(1);
t393 = pkin(1) * qJD(1);
t219 = t265 * t393;
t244 = qJD(1) + qJD(2);
t135 = pkin(3) * t244 + t219;
t263 = cos(qJ(4));
t420 = (qJD(4) * t135 + t421) * t263;
t242 = qJDD(1) + qJDD(2);
t347 = qJDD(1) * t265;
t215 = pkin(1) * t347;
t360 = qJD(1) * t257;
t341 = pkin(1) * t360;
t323 = qJD(2) * t341;
t289 = t215 - t323;
t419 = -pkin(2) * t242 + qJD(3) * t341 - t289;
t339 = qJD(6) * t393;
t416 = t265 * t339 + t421;
t368 = qJ(2) + qJ(4);
t229 = sin(t368);
t255 = sin(qJ(4));
t137 = -pkin(1) * t229 - pkin(3) * t255;
t246 = 0.1e1 / t255;
t269 = 0.1e1 / pkin(3);
t379 = t246 * t269;
t337 = t137 * t379;
t415 = t337 + 0.1e1;
t248 = qJ(1) + qJ(2);
t231 = sin(t248);
t233 = cos(t248);
t253 = sin(qJ(7));
t261 = cos(qJ(7));
t414 = (t231 * t253 + t233 * t261) * pkin(3);
t254 = sin(qJ(6));
t344 = qJD(2) + qJD(6);
t413 = t254 * t360 * t344;
t237 = qJ(1) + t368;
t204 = sin(t237);
t207 = cos(t237);
t412 = -g(1) * t204 + g(2) * t207;
t238 = qJ(3) + t248;
t205 = sin(t238);
t208 = cos(t238);
t411 = g(1) * t205 - g(2) * t208;
t304 = -g(1) * t231 + g(2) * t233;
t410 = pkin(1) * t344;
t367 = qJ(3) - qJ(6);
t333 = -qJ(9) - t367;
t94 = 0.1e1 / ((-cos(qJ(4) - t367) + cos(qJ(4) + t367)) * pkin(6) + (cos(qJ(4) + t333) - cos(qJ(4) - t333)) * pkin(2));
t409 = pkin(2) * t94;
t252 = sin(qJ(8));
t408 = pkin(1) * t252;
t266 = cos(qJ(1));
t407 = pkin(1) * t266;
t268 = 0.1e1 / pkin(4);
t406 = pkin(1) * t268;
t405 = pkin(2) * t233;
t221 = qJDD(4) + t242;
t404 = pkin(3) * t221;
t403 = pkin(3) * t233;
t258 = sin(qJ(1));
t400 = g(1) * t258;
t398 = g(2) * t266;
t239 = t258 * pkin(1);
t370 = t257 * t263;
t292 = -t255 * t265 - t370;
t117 = t292 * t393;
t374 = t255 * t257;
t291 = t263 * t265 - t374;
t119 = t291 * t393;
t210 = pkin(3) * t263 + pkin(4);
t249 = sin(qJ(10));
t250 = cos(qJ(10));
t349 = qJD(10) * t255;
t350 = qJD(10) * t250;
t378 = t249 * t255;
t397 = -t117 * t249 - t119 * t250 + t210 * t350 - (t249 * t349 + (-t250 * t263 + t378) * qJD(4)) * pkin(3);
t351 = qJD(10) * t249;
t377 = t250 * t255;
t396 = t117 * t250 - t119 * t249 + t210 * t351 + (t250 * t349 + (t249 * t263 + t377) * qJD(4)) * pkin(3);
t256 = sin(qJ(3));
t264 = cos(qJ(3));
t369 = t257 * t264;
t290 = t256 * t265 + t369;
t118 = t290 * t393;
t165 = t256 * t341;
t120 = -t264 * t219 + t165;
t211 = -pkin(2) * t264 + pkin(6);
t251 = sin(qJ(9));
t259 = cos(qJ(9));
t352 = qJD(9) * t259;
t353 = qJD(9) * t251;
t376 = t251 * t256;
t395 = -t118 * t251 - t120 * t259 + t211 * t352 - (-t256 * t353 + (t259 * t264 - t376) * qJD(3)) * pkin(2);
t372 = t256 * t259;
t394 = t118 * t259 - t120 * t251 + t211 * t353 + (-t256 * t352 + (-t251 * t264 - t372) * qJD(3)) * pkin(2);
t267 = 0.1e1 / pkin(5);
t247 = qJ(1) + qJ(8);
t245 = pkin(8) + qJ(5);
t311 = t245 + t367;
t285 = t311 - t247;
t312 = t245 - t367;
t286 = t312 - t247;
t392 = t267 / (pkin(6) * (cos(t286) - cos(t285)) + (-cos(qJ(9) - t286) + cos(qJ(9) + t285)) * pkin(2));
t296 = t204 * t253 + t207 * t261;
t77 = t296 * pkin(4) + t414;
t391 = t268 * t77;
t390 = t269 * t94;
t102 = t135 * t255 + t263 * t341;
t389 = t102 * t249;
t388 = t102 * t250;
t136 = pkin(2) * t244 + t219;
t103 = t136 * t256 + t264 * t341;
t387 = t103 * t251;
t386 = t103 * t259;
t297 = t204 * t233 - t207 * t231;
t104 = 0.1e1 / t297;
t385 = t104 * t296;
t384 = t137 * t246;
t383 = t137 * t269;
t138 = pkin(3) * t257 + pkin(4) * t229;
t382 = t138 * t269;
t330 = -qJ(1) + t245;
t381 = -0.1e1 / sin(qJ(8) - t330) * t252;
t243 = qJD(1) + qJD(8);
t380 = t243 * t252;
t107 = pkin(3) * t242 + t289;
t375 = t255 * t107;
t373 = t256 * t257;
t262 = cos(qJ(6));
t371 = t257 * t262;
t217 = -qJ(9) + t237;
t166 = t217 - t367;
t334 = qJ(9) + t248;
t216 = qJ(4) + t334;
t167 = t216 + t367;
t366 = sin(t167) - sin(t166);
t365 = cos(t167) - cos(t166);
t364 = -cos(0.2e1 * t247) + cos(0.2e1 * t245);
t363 = sin(t217) - sin(t216);
t362 = cos(t217) - cos(t216);
t201 = pkin(3) * t231;
t361 = t201 + t239;
t358 = qJD(3) * t256;
t357 = qJD(3) * t264;
t355 = qJD(4) * t255;
t354 = qJD(4) * t263;
t346 = qJDD(1) * pkin(1) ^ 2;
t260 = cos(qJ(8));
t345 = t260 * qJDD(1);
t343 = (-t365 * t253 + t366 * t261) * t409;
t241 = qJDD(1) + qJDD(8);
t342 = pkin(1) * t243 * t260;
t340 = 0.1e1 / t364 * t267 * (-t364 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1));
t338 = qJD(8) * t393;
t336 = t268 * t382;
t331 = t257 * t354;
t328 = t262 * t347;
t202 = pkin(2) * t231;
t327 = -pkin(6) * t205 + t202;
t226 = qJD(3) + t244;
t225 = qJD(4) + t244;
t313 = qJ(4) + t330;
t300 = -qJ(8) + t313;
t314 = -qJ(4) + t330;
t301 = -qJ(8) + t314;
t326 = pkin(1) / (cos(t301) - cos(t300)) * t267 * (pkin(3) * (cos(t314) - cos(t313)) + (cos(qJ(2) - t301) - cos(qJ(2) + t300)) * pkin(5));
t101 = -t136 * t264 + t165;
t240 = t265 * pkin(1);
t213 = t240 + pkin(2);
t125 = pkin(1) * t373 - t213 * t264;
t325 = qJD(1) * (-qJD(2) + t244);
t324 = qJD(2) * (-qJD(1) - t244);
t222 = qJDD(3) + t242;
t220 = qJDD(6) + t242;
t321 = qJD(4) * t341;
t318 = t246 * t336;
t317 = -t254 * t215 - t416 * t262;
t230 = sin(t247);
t232 = cos(t247);
t316 = g(1) * t232 + g(2) * t230 + qJDD(1) * t408 + t260 * t338;
t223 = pkin(1) * t398;
t315 = -g(1) * t239 + t223;
t212 = t240 + pkin(3);
t124 = -pkin(1) * t374 + t263 * t212;
t310 = t215 + t304;
t309 = g(1) * t230 - g(2) * t232 + t252 * t338;
t178 = pkin(4) * t204;
t110 = t178 + t361;
t306 = -pkin(4) * t207 - t403;
t111 = -t306 + t407;
t139 = -pkin(5) * t232 + t407;
t140 = -pkin(5) * t230 + t239;
t234 = qJ(9) + t245;
t235 = -qJ(9) + t245;
t308 = ((t363 * t110 + t362 * t111) * t390 + (-(sin(t235) - sin(t234)) * t140 - (cos(t235) - cos(t234)) * t139) * t392) * pkin(2) + t337;
t307 = -t405 - t407;
t305 = -g(1) * t233 - g(2) * t231;
t106 = t263 * t107;
t303 = -t135 * t355 + t106;
t302 = t416 * t254 + t262 * t323 + t339 * t371;
t272 = -t255 * t348 + (-t255 * t359 - t331) * qJD(1);
t53 = t272 * pkin(1) + t303;
t45 = t221 * pkin(4) + t53;
t149 = t255 * t321;
t52 = -t149 + t375 + t420;
t100 = t263 * t135 - t255 * t341;
t96 = pkin(4) * t225 + t100;
t24 = t102 * t350 + t249 * t52 - t250 * t45 + t96 * t351;
t55 = t136 * t358 + t256 * t421 + t419 * t264;
t46 = pkin(6) * t222 + t55;
t54 = -t136 * t357 + t419 * t256 - t264 * t421;
t97 = pkin(6) * t226 + t101;
t26 = -t103 * t352 + t251 * t54 - t259 * t46 + t97 * t353;
t191 = qJDD(9) + t222;
t113 = pkin(4) + t124;
t127 = pkin(1) * t370 + t212 * t255;
t70 = -t113 * t250 + t127 * t249;
t299 = t113 * t249 + t127 * t250;
t114 = pkin(6) + t125;
t128 = -pkin(1) * t369 - t213 * t256;
t74 = -t114 * t259 + t128 * t251;
t298 = t114 * t251 + t128 * t259;
t294 = t254 * t265 + t371;
t293 = t254 * t257 - t262 * t265;
t182 = qJDD(10) + t221;
t288 = (t362 * t253 - t363 * t261) * t409 - t385;
t287 = -g(1) * t207 - g(2) * t204 + t149;
t284 = t103 * t353 + t251 * t46 + t259 * t54 + t97 * t352;
t133 = t293 * pkin(1);
t23 = t102 * t351 - t249 * t45 - t250 * t52 - t96 * t350;
t281 = t287 - t375;
t218 = qJ(3) + t334;
t185 = sin(t218);
t188 = cos(t218);
t280 = -g(1) * t185 + g(2) * t188 + t26;
t209 = qJ(10) + t237;
t179 = sin(t209);
t180 = cos(t209);
t279 = g(1) * t179 - g(2) * t180 + t24;
t278 = g(1) * t180 + g(2) * t179 - t23;
t277 = t55 + t411;
t276 = -g(1) * t188 - g(2) * t185 + t284;
t275 = g(1) * t208 + g(2) * t205 - t54;
t236 = qJ(6) + t248;
t203 = sin(t236);
t206 = cos(t236);
t274 = -pkin(1) * t328 + g(1) * t203 - g(2) * t206 + t302;
t273 = -pkin(1) * t413 + g(1) * t206 + g(2) * t203 - t317;
t271 = (-pkin(3) * t225 - t135) * qJD(4) - t421;
t228 = cos(t245);
t227 = sin(t245);
t224 = qJD(1) + t344;
t196 = qJD(9) + t226;
t195 = qJD(10) + t225;
t190 = qJ(9) + t311;
t189 = -qJ(9) + t312;
t177 = pkin(6) * t208;
t134 = t294 * pkin(1);
t126 = pkin(2) * t372 - t211 * t251;
t123 = -pkin(2) * t376 - t211 * t259;
t122 = -pkin(3) * t377 - t210 * t249;
t121 = pkin(3) * t378 - t210 * t250;
t116 = t294 * t393;
t115 = qJD(1) * t133;
t93 = t294 * t410;
t92 = t293 * t410;
t89 = (t265 * t325 - t348) * pkin(1) + t305;
t88 = t257 * pkin(1) * t325 + t310;
t85 = t213 * t358 + (t290 * qJD(2) + t257 * t357) * pkin(1);
t84 = -t213 * t357 + (t257 * t358 + (-t264 * t265 + t373) * qJD(2)) * pkin(1);
t83 = -t212 * t355 + (t292 * qJD(2) - t331) * pkin(1);
t82 = t212 * t354 + (t291 * qJD(2) - t257 * t355) * pkin(1);
t66 = t297 * pkin(4) - t414;
t63 = -t101 * t259 - t387;
t62 = t101 * t251 - t386;
t61 = -t100 * t250 + t389;
t60 = t100 * t249 + t388;
t59 = t251 * t97 - t386;
t58 = -t259 * t97 - t387;
t57 = t249 * t96 + t388;
t56 = -t250 * t96 + t389;
t42 = t115 * t224 + t273;
t41 = -t116 * t224 + t274;
t38 = -t103 * t226 + t277;
t37 = t101 * t226 + t275;
t36 = t102 * t225 + t412 + t53;
t35 = t100 * t225 + t281 - t420;
t34 = t298 * qJD(9) + t251 * t84 - t259 * t85;
t33 = t74 * qJD(9) - t251 * t85 - t259 * t84;
t32 = t299 * qJD(10) + t249 * t82 - t250 * t83;
t31 = t70 * qJD(10) - t249 * t83 - t250 * t82;
t30 = t120 * t226 + (t222 * t256 + t226 * t357) * pkin(2) + t275;
t29 = -t118 * t226 + (-t222 * t264 + t226 * t358) * pkin(2) + t277;
t28 = t119 * t225 + (-t107 - t404) * t255 + t271 * t263 + t287;
t27 = -t117 * t225 + t106 + (-t321 + t404) * t263 + t271 * t255 + t412;
t19 = -t196 * t59 + t280;
t18 = t196 * t58 + t276;
t17 = -t101 * t118 + t103 * t120 + (-t256 * t54 - t264 * t55 + (t101 * t256 + t103 * t264) * qJD(3) + t304) * pkin(2);
t16 = -t100 * t117 - t102 * t119 + (t255 * t52 + t263 * t53 + (-t100 * t255 + t102 * t263) * qJD(4) + t304) * pkin(3);
t15 = -t195 * t57 + t279;
t14 = t195 * t56 + t278;
t13 = ((-t366 * t110 - t365 * t111) * t390 + ((sin(t190) - sin(t189)) * t140 + (cos(t190) - cos(t189)) * t139) * t392) * pkin(2);
t12 = t63 * t196 + (t191 * t251 + t196 * t352) * pkin(6) + t276;
t11 = -t62 * t196 + (-t191 * t259 + t196 * t353) * pkin(6) + t280;
t10 = t61 * t195 + (t182 * t249 + t195 * t350) * pkin(4) + t278;
t9 = -t60 * t195 + (-t182 * t250 + t195 * t351) * pkin(4) + t279;
t8 = t123 * t191 + t394 * t196 + t280;
t7 = -t126 * t191 + t395 * t196 + t276;
t6 = t121 * t182 + t396 * t195 + t279;
t5 = -t122 * t182 + t397 * t195 + t278;
t4 = -t58 * t62 + t59 * t63 + (t284 * t251 - t259 * t26 + (t251 * t58 + t259 * t59) * qJD(9) + t411) * pkin(6);
t3 = -t56 * t60 + t57 * t61 + (-t23 * t249 - t24 * t250 + (t249 * t56 + t250 * t57) * qJD(10) + t412) * pkin(4);
t2 = -t284 * t126 + t26 * t123 - g(2) * (t177 - t405) - g(1) * t327 + t395 * t59 + t394 * t58;
t1 = t23 * t122 + t24 * t121 - g(2) * t306 - g(1) * (t201 + t178) + t397 * t57 + t396 * t56;
t20 = [0, 0, 0, 0, 0, qJDD(1), t398 - t400, -g(1) * t266 - g(2) * t258, 0, 0, 0, 0, 0, 0, 0, t415 * t242, t88 * t337 + (t242 * t265 + t257 * t324) * pkin(1) + t310, t89 * t337 + ((-qJDD(1) - t242) * t257 + t265 * t324) * pkin(1) + t305, 0, (t257 ^ 2 + t265 ^ 2) * t346 + t315, 0, 0, 0, 0, 0, (t13 + t415) * t222, t125 * t222 + t13 * t38 + t226 * t85 + t29 * t337 + t277, -t128 * t222 + t13 * t37 - t226 * t84 + t30 * t337 + t275, 0, t54 * t128 - t103 * t84 + t55 * t125 + t101 * t85 - g(2) * t307 - g(1) * (t202 + t239) + t17 * t337, 0, 0, 0, 0, 0, ((t138 * t406 + t137) * t379 + 0.1e1) * t221, t27 * t337 + t124 * t221 + t83 * t225 + (t36 * t318 + t272) * pkin(1) + t303 + t412, t28 * t337 - t135 * t354 - t127 * t221 - t82 * t225 + (-t283 * t263 + t35 * t318) * pkin(1) + t281, 0, t52 * t127 + t102 * t82 + t53 * t124 + t100 * t83 - g(2) * (-t403 - t407) - g(1) * t361 + t16 * t337, 0, 0, 0, 0, 0, qJDD(5) * t381, (g(1) * t227 - g(2) * t228) * t381, (g(1) * t228 + g(2) * t227) * t381, 0, 0, 0, 0, 0, 0, 0, (t308 + 0.1e1) * t220, t133 * t220 + t93 * t224 + t308 * t41 + t274, t134 * t220 - t92 * t224 + t308 * t42 + t273, 0, -t317 * t134 - t116 * t92 + t302 * t133 + t115 * t93 + t223 + (-t133 * t328 - t134 * t413 - t400) * pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t340 + 0.1e1) * t241, t309 * t340 + (qJD(8) * t380 - t260 * t241 - t345 + (-qJD(1) * t380 - t345) * t340) * pkin(1) + t309, qJD(8) * t342 + t241 * t408 + (-qJD(1) * t342 + t316) * t340 + t316, 0, (t252 ^ 2 + t260 ^ 2) * t346 + t315, 0, 0, 0, 0, 0, (t13 + (-t326 + t384) * t269 + 0.1e1) * t191, t13 * t11 + t74 * t191 + t34 * t196 + (-t19 * t326 + t8 * t384) * t269 + t280, t13 * t12 + t298 * t191 - t33 * t196 + (-t18 * t326 + t7 * t384) * t269 + t276, 0, t284 * t298 - t59 * t33 + t26 * t74 + t58 * t34 - g(2) * (t177 + t307) - g(1) * (t239 + t327) + t13 * t4 + t2 * t337, 0, 0, 0, 0, 0, ((t383 + (-t257 + t382) * t406) * t246 + 0.1e1) * t182, t70 * t182 + t32 * t195 + (t6 * t383 + (-t15 * t257 + t9 * t382) * t406) * t246 + t279, t299 * t182 - t31 * t195 + (t5 * t383 + (t10 * t382 - t14 * t257) * t406) * t246 + t278, 0, -g(1) * t110 + g(2) * t111 - t23 * t299 + t24 * t70 - t57 * t31 + t56 * t32 + (pkin(1) * t3 * t336 + t1 * t383) * t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242 * t385, -t88 * t385, -t89 * t385, 0, 0, 0, 0, 0, 0, 0, (t343 - t385) * t222, -t29 * t385 + t38 * t343, -t30 * t385 + t37 * t343, 0, -t17 * t385, 0, 0, 0, 0, 0, (-t296 + t391) * t221 * t104, (-t27 * t296 + t36 * t391) * t104, (-t28 * t296 + t35 * t391) * t104, 0, -t16 * t385, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288 * t220, t288 * t41, t288 * t42, 0, 0, 0, 0, 0, 0, 0, qJDD(7), -g(1) * t261 - g(2) * t253, g(1) * t253 - g(2) * t261, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t343 - 0.2e1 * t385) * t191, t11 * t343 - (t19 + t8) * t385, t12 * t343 - (t18 + t7) * t385, 0, -t2 * t385 + t343 * t4, 0, 0, 0, 0, 0, (-t296 + (t66 + t77) * t268) * t182 * t104, (-t296 * t6 + (t15 * t66 + t77 * t9) * t268) * t104, (-t296 * t5 + (t10 * t77 + t14 * t66) * t268) * t104, 0, (-t1 * t296 + t3 * t391) * t104;];
tau_reg = t20;
