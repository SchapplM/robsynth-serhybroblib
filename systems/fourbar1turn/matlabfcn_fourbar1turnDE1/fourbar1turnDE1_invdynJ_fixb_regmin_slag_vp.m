% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbar1turnDE1
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
% tau_reg [2x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:24
% EndTime: 2020-06-27 16:35:59
% DurationCPUTime: 10.40s
% Computational Cost: add. (110564->307), mult. (158145->727), div. (5768->21), fcn. (42667->10), ass. (0->293)
t389 = pkin(1) * pkin(2);
t378 = pkin(4) ^ 2;
t120 = 0.1e1 / t378;
t124 = pkin(2) ^ 2;
t125 = pkin(1) ^ 2;
t113 = cos(qJ(2));
t334 = pkin(2) * t113;
t365 = -2 * pkin(1);
t276 = t334 * t365 + t125;
t105 = t124 + t276;
t102 = 0.1e1 / t105;
t103 = 0.1e1 / t105 ^ 2;
t104 = t102 * t103;
t111 = sin(qJ(2));
t335 = pkin(2) * t111;
t258 = pkin(1) * t335;
t208 = 0.4e1 * t258;
t184 = t104 * t208;
t377 = pkin(3) ^ 2;
t274 = t377 - t378;
t101 = t105 - t274;
t106 = pkin(1) - t334;
t355 = -pkin(3) - pkin(4);
t98 = (pkin(2) - t355) * (pkin(2) + t355) + t276;
t354 = pkin(4) - pkin(3);
t99 = (pkin(2) - t354) * (pkin(2) + t354) + t276;
t312 = t98 * t99;
t128 = sqrt(-t312);
t284 = t111 * t128;
t82 = -pkin(2) * t284 + t101 * t106;
t75 = t82 ^ 2;
t96 = t101 * t335;
t83 = t106 * t128 + t96;
t79 = t83 ^ 2;
t307 = t75 + t79;
t150 = t307 * t184;
t256 = 0.2e1 * t103;
t262 = 0.2e1 * t106 * pkin(1);
t282 = t113 * t128;
t217 = (-t98 - t99) * t389;
t88 = t111 * t217;
t92 = 0.1e1 / t128;
t313 = t92 * t88;
t53 = t96 + (-t282 + (t262 - t313) * t111) * pkin(2);
t109 = t111 ^ 2;
t288 = t109 * t124;
t244 = pkin(1) * t288;
t300 = t106 * t92;
t55 = t88 * t300 + 0.2e1 * t244 + (t101 * t113 + t284) * pkin(2);
t23 = ((t53 * t82 + t55 * t83) * t256 - t150) * t120;
t388 = t23 / 0.2e1;
t116 = qJD(2) ^ 2;
t280 = t116 * t124;
t287 = t109 * t125;
t387 = 0.2e1 * pkin(2) * t287;
t386 = 0.8e1 * t104;
t100 = t105 + t274;
t107 = pkin(1) * t113 - pkin(2);
t81 = -pkin(1) * t284 - t100 * t107;
t71 = t81 ^ 2;
t72 = 0.1e1 / t81;
t74 = t72 / t71;
t342 = pkin(1) * t100;
t97 = t111 * t342;
t84 = -t107 * t128 + t97;
t80 = t84 ^ 2;
t318 = t74 * t80;
t73 = 0.1e1 / t81 ^ 2;
t319 = t73 * t84;
t264 = -0.2e1 * pkin(2) * t107;
t229 = t100 + t264;
t87 = qJD(2) * t88;
t314 = t87 * t92;
t245 = t111 * t314;
t48 = (-t245 + (t229 * t111 - t282) * qJD(2)) * pkin(1);
t272 = qJD(2) * t111;
t233 = t128 * t272;
t271 = qJD(2) * t113;
t299 = t107 * t92;
t51 = pkin(1) * t233 + qJD(2) * t387 + t271 * t342 - t87 * t299;
t69 = t73 * t80 + 0.1e1;
t67 = 0.1e1 / t69 ^ 2;
t327 = (-t48 * t318 + t51 * t319) * t67;
t385 = -0.2e1 * t327;
t77 = 0.1e1 / t82 ^ 2;
t317 = t77 * t83;
t119 = 0.1e1 / pkin(4);
t240 = pkin(2) * t272;
t214 = t103 * t240;
t195 = pkin(1) * t214;
t348 = -t102 / 0.2e1;
t49 = (-t245 + (-t282 + (t101 + t262) * t111) * qJD(2)) * pkin(2);
t35 = (t82 * t195 + t49 * t348) * t119;
t347 = t102 / 0.2e1;
t382 = 0.2e1 * qJD(2);
t50 = t244 * t382 + t87 * t300 + (t101 * t271 + t233) * pkin(2);
t37 = (-t83 * t195 + t50 * t347) * t119;
t68 = t77 * t79 + 0.1e1;
t64 = 0.1e1 / t68;
t76 = 0.1e1 / t82;
t156 = t64 * (t35 * t317 + t37 * t76);
t286 = t111 * t113;
t381 = 0.6e1 * t116 * t286;
t333 = pkin(3) * t105;
t66 = 0.1e1 / t69;
t226 = t66 * t73 * t333;
t202 = t84 * t226;
t255 = t72 * t333;
t227 = t66 * t255;
t122 = 0.1e1 / pkin(3);
t179 = -0.2e1 * t195;
t36 = (t102 * t48 + t179 * t81) * t122;
t38 = (t102 * t51 + t179 * t84) * t122;
t17 = -t202 * t36 + t38 * t227 + qJD(2);
t291 = t102 * t122;
t247 = pkin(2) * t291;
t215 = qJD(2) * t247;
t123 = 0.1e1 / t377;
t308 = t71 + t80;
t63 = t308 * t123 * t103;
t58 = t63 ^ (-0.1e1 / 0.2e1);
t192 = t58 * t215;
t246 = t103 * t335;
t228 = pkin(1) * t246;
t206 = -0.2e1 * t228;
t52 = t97 + (-t282 + (t264 - t313) * t111) * pkin(1);
t40 = (t102 * t52 + t206 * t81) * t122;
t54 = -t88 * t299 + t387 + (t100 * t113 + t284) * pkin(1);
t42 = (t102 * t54 + t206 * t84) * t122;
t20 = -t202 * t40 + t42 * t227 + 0.1e1;
t239 = t58 * t291;
t221 = pkin(2) * t239;
t151 = t308 * t184;
t24 = ((t52 * t81 + t54 * t84) * t256 - t151) * t123;
t59 = t58 / t63;
t328 = t24 * t59;
t177 = t40 * t319 - t42 * t72;
t209 = 0.2e1 * t258;
t361 = 0.2e1 * t84;
t368 = -t36 * t319 + t38 * t72;
t375 = ((t177 * qJD(2) + t368) * t66 * t209 + (-0.2e1 * t368 * t67 * (-t52 * t318 + t54 * t319) + t177 * t385 + ((t52 * t36 - t48 * t40) * t74 * t361 + (-t54 * t36 - t38 * t52 + t51 * t40 + t42 * t48) * t73) * t66) * t105) * pkin(3);
t283 = t113 * t116;
t369 = qJDD(2) * t111 + t283;
t163 = (t87 ^ 2 / t312 + t369 * t217 - 0.4e1 * t287 * t280) * t92;
t278 = t116 * t128;
t145 = -t163 + t278;
t161 = qJDD(2) * t128 + t314 * t382;
t223 = pkin(1) * t240;
t236 = t109 * t280;
t336 = pkin(2) * t103;
t265 = -0.2e1 * t336;
t269 = qJDD(2) * t100;
t273 = qJD(2) * t103;
t279 = t116 * t125;
t285 = t111 * t122;
t293 = t102 * t107;
t337 = pkin(2) * t102;
t359 = -0.4e1 * qJD(2);
t362 = -0.2e1 * t84;
t376 = 0.2e1 * qJDD(2);
t6 = qJDD(2) + (-t163 * t293 + (t84 * t236 * t386 + (t109 * t376 + t381) * t337) * t125 + (((t269 + t278) * t102 + t84 * t116 * t265) * t113 + ((-t100 * t116 + t161) * t102 + (qJDD(2) * t362 + t51 * t359) * t336) * t111) * pkin(1)) * t122 * t227 + t38 * t255 * t385 - ((t124 * t81 * t386 + 0.4e1 * t337) * t122 * t109 * t279 + ((-t161 * t102 + (t229 * t102 + t81 * t265) * t116) * t122 * t113 + ((t145 + t269) * t102 + (-0.4e1 * t48 * t273 + (-t103 * t81 - t293) * t376) * pkin(2)) * t285) * pkin(1)) * t202 + (-t36 * t51 - t38 * t48) * t226 + (t74 * t66 * t48 + t73 * t327) * t36 * t333 * t361 + 0.2e1 * t368 * pkin(3) * t66 * t223;
t379 = -t375 * t192 + t17 * t215 * t328 / 0.2e1 + (qJDD(2) * t20 + t6) * t221;
t332 = pkin(4) * t105;
t263 = 0.2e1 * t332;
t306 = -t84 + t52;
t374 = t306 * t58;
t237 = -t328 / 0.2e1;
t164 = t84 * t237 + t54 * t58;
t373 = t164 * t113;
t188 = t58 * t256 * t389;
t372 = (t109 * t81 + t84 * t286) * t188;
t18 = t156 * t263;
t39 = (t82 * t228 + t53 * t348) * t119;
t41 = (-t83 * t228 + t55 * t347) * t119;
t155 = t64 * (-t39 * t317 - t41 * t76);
t19 = t155 * t263;
t62 = t307 * t120 * t103;
t56 = t62 ^ (-0.1e1 / 0.2e1);
t371 = (qJD(2) * t19 + t18) * t103 * t56 * t209;
t370 = t17 * t221 + t20 * t192;
t22 = ((t48 * t81 + t51 * t84) * t256 - qJD(2) * t151) * t123;
t330 = t22 * t59;
t366 = t247 * t330 * (qJD(2) * t20 + t17);
t60 = 0.1e1 / t62;
t364 = -0.2e1 * t77;
t363 = 0.4e1 * t83;
t21 = ((t49 * t82 + t50 * t83) * t256 - qJD(2) * t150) * t120;
t357 = t21 / 0.2e1;
t356 = t81 / 0.2e1;
t65 = 0.1e1 / t68 ^ 2;
t231 = 0.4e1 * t65 * t317;
t78 = t76 / t75;
t316 = t78 * t79;
t25 = -t49 * t316 + t50 * t317;
t320 = t65 * t76;
t253 = t25 * t320;
t259 = t78 * t363;
t261 = 0.2e1 * t64 * t77;
t27 = -t53 * t316 + t55 * t317;
t353 = ((-qJD(2) * t155 - t156) * t208 + ((t53 * t261 + 0.4e1 * t27 * t320) * t37 + (t27 * t231 + (t53 * t259 + t55 * t364) * t64) * t35 - (t49 * t261 + 0.4e1 * t253) * t41 - (t25 * t231 + (t49 * t259 + t50 * t364) * t64) * t39) * t105) * pkin(4);
t351 = g(3) * t82;
t350 = g(3) * t83;
t144 = -t101 * t116 + t161;
t183 = -t102 * t106 + t103 * t82;
t254 = t64 * t332;
t225 = t77 * t254;
t268 = qJDD(2) * t101;
t340 = pkin(1) * t103;
t341 = pkin(1) * t102;
t5 = 0.2e1 * (-t35 * t50 + t37 * t49) * t225 - 0.4e1 * pkin(4) * t223 * t156 + 0.2e1 * (t37 * t253 + (t77 * t65 * t25 + t78 * t64 * t49) * t35 * t83) * t263 + 0.2e1 * (-((pkin(1) * t124 * t381 + t106 * t163) * t347 + (t104 * t279 * t363 + qJDD(2) * t341) * t288 + (((t268 + t278) * t113 + t144 * t111) * t347 + (-t83 * t283 + (-0.2e1 * qJD(2) * t50 - qJDD(2) * t83) * t111) * t340) * pkin(2)) * t76 * t254 - ((-0.4e1 * t104 * t125 * t82 - 0.2e1 * t341) * t236 + ((t183 * t116 * pkin(1) + t144 * t347) * t113 + ((t145 + t268) * t348 + (qJDD(2) * t183 + 0.2e1 * t49 * t273) * pkin(1)) * t111) * pkin(2)) * t83 * t225) * t119;
t349 = t5 * t56;
t112 = sin(qJ(1));
t346 = g(1) * t112;
t114 = cos(qJ(1));
t345 = g(1) * t114;
t344 = g(2) * t112;
t343 = g(2) * t114;
t117 = (qJD(1) ^ 2);
t339 = pkin(1) * t117;
t57 = t56 * t60;
t329 = t23 * t57;
t324 = t57 * t82;
t323 = t57 * t83;
t322 = t60 * t83;
t61 = 0.1e1 / t62 ^ 2;
t321 = t61 * t79;
t315 = t82 * t83;
t212 = t113 * t239;
t238 = t58 * t285;
t213 = t102 * t238;
t310 = (t212 * t84 + t213 * t81) * t114;
t309 = t54 + t81;
t302 = pkin(1) * qJD(1);
t301 = pkin(2) * qJD(1);
t298 = t111 * t81;
t297 = t111 * t84;
t296 = t113 * t81;
t295 = t113 * t84;
t185 = t295 + t298;
t152 = t102 * t58 * t185;
t33 = t122 * t152;
t31 = qJD(1) * t33;
t281 = t114 * t122;
t277 = t117 * t120;
t275 = -t113 ^ 2 + t109;
t270 = qJDD(1) * t102;
t266 = t113 * qJDD(1);
t257 = -2 * t302;
t250 = t61 * t315;
t248 = -0.2e1 * t117 * t125;
t242 = t111 * t301;
t241 = t113 * t301;
t232 = qJDD(1) * t103 * t60;
t222 = t21 * t57 * t302;
t220 = t56 * t246;
t216 = t328 * t356;
t207 = t19 * t56 * t270;
t203 = t104 * t60 * t258;
t194 = t344 + t345;
t193 = -t343 + t346;
t191 = t56 * t214;
t190 = t84 * t213;
t47 = t81 * t212;
t186 = -t296 + t297;
t180 = t79 * t203;
t176 = t17 * t192;
t169 = t18 * t388 + t19 * t357;
t166 = t203 * t315;
t165 = -t295 / 0.2e1 - t298 / 0.2e1;
t162 = pkin(1) * t18 * t191;
t157 = 0.2e1 * t194;
t154 = t194 + t339;
t153 = qJDD(1) * t365 - t193;
t149 = t59 * (t296 / 0.2e1 - t297 / 0.2e1);
t148 = t20 * t238 * t280 * t340;
t147 = 0.4e1 * qJD(1) * t125 * t191;
t146 = -t339 / 0.2e1 - t345 / 0.2e1 - t344 / 0.2e1;
t141 = (-t109 * t84 + t81 * t286) * t188;
t139 = (t372 + (-t373 + (t216 - t374) * t111) * t102) * t122;
t138 = t122 * (t141 + (t24 * t149 + (t309 * t111 - t306 * t113) * t58) * t102);
t137 = qJD(2) * t372 + (-t165 * t330 + (qJD(2) * t186 - t111 * t48 - t113 * t51) * t58) * t102;
t136 = qJD(2) * t141 + (t22 * t149 + (qJD(2) * t185 + t111 * t51 - t113 * t48) * t58) * t102;
t46 = qJD(1) * t47;
t43 = t112 * t47;
t34 = -t47 + t190;
t32 = qJD(1) * t190 - t46;
t12 = qJD(1) * t138;
t11 = qJD(1) * t139 - t46;
t10 = t136 * t122;
t9 = t137 * t122;
t8 = (t186 * t58 * t270 + qJD(1) * t136) * t122;
t7 = (qJD(1) * t137 - qJDD(1) * t152) * t122;
t1 = [qJDD(1), t193, t194, 0.2e1 * qJD(1) * t111 * t271 + qJDD(1) * t109, -0.2e1 * t275 * qJD(2) * qJD(1) + 0.2e1 * t111 * t266, t369, qJDD(2) * t113 - t111 * t116, 0, t193 * t113, -t193 * t111, -t31 * t9 - t33 * t7, -t10 * t31 + t32 * t9 - t33 * t8 + t34 * t7, t17 * t9 - t33 * t6, t10 * t32 + t34 * t8, t10 * t17 + t34 * t6, 0, -g(1) * t43 + (-t186 * t343 + t297 * t346) * t239 + ((-qJD(1) * t34 - t32) * t272 + (qJD(1) * t10 + qJDD(1) * t34 + t8) * t113) * pkin(2), -g(2) * t310 + t33 * t346 + (-0.2e1 * t31 * t272 + (-qJD(1) * t9 + qJDD(1) * t33 - t7) * t113) * pkin(2), (t79 * t232 + (t180 * t359 + (-t21 * t321 + 0.2e1 * t50 * t322) * t103) * qJD(1)) * t120, (-0.2e1 * t232 * t315 + (0.8e1 * qJD(2) * t166 + (t21 * t250 + (-t49 * t83 - t50 * t82) * t60) * t256) * qJD(1)) * t120, (0.2e1 * t83 * t162 + (t83 * t349 - (t50 * t56 - t21 * t323 / 0.2e1) * t18) * t102) * t119, (-0.2e1 * t82 * t162 + (-t82 * t349 - (t324 * t357 - t49 * t56) * t18) * t102) * t119, 0, (t82 * t147 + (t82 * t222 + (t153 * t82 + t257 * t49) * t56) * t102) * t119, (t83 * t147 + (t83 * t222 + (t153 * t83 + t257 * t50) * t56) * t102) * t119; 0, 0, 0, -t117 * t286, t275 * t117, t111 * qJDD(1), t266, qJDD(2), -g(3) * t113 + t111 * t194, g(3) * t111 + t113 * t194, t31 * t11, -t11 * t32 + t12 * t31, -t11 * t17 + t20 * t7 + t31 * t375, -t32 * t12, -t12 * t17 + t20 * t8 - t32 * t375, -t17 * t375 + t20 * t6, t32 * t242 - t12 * t241 + t52 * t176 - g(1) * ((t141 + ((-t52 * t58 + t216) * t113 + t164 * t111) * t102) * t281 + t310) - t138 * t344 - g(3) * (-t47 + t139) - t370 * t48 + t356 * t366 + (0.2e1 * t148 - t379) * t81, t148 * t362 + t31 * t242 + t11 * t241 - t54 * t176 - g(1) * (-t372 + (t165 * t328 + (t306 * t111 + t309 * t113) * t58) * t102) * t281 - g(2) * (t43 + (-t372 + (t373 + (t81 * t237 + t374) * t111) * t102) * t122 * t112) - g(3) * t138 + t370 * t51 + (-t366 / 0.2e1 + t379) * t84, (0.2e1 * t180 + (t321 * t388 - t55 * t322) * t103) * t277, (-0.4e1 * t166 + (-t23 * t250 + (t53 * t83 + t55 * t82) * t60) * t103) * t277, (t83 * t207 + (-t83 * t371 + (-t169 * t323 + (t18 * t55 + t19 * t50 - t353 * t83) * t56) * t102) * qJD(1)) * t119, (-t82 * t207 + (t82 * t371 + (t169 * t324 + (-t53 * t18 - t49 * t19 + t353 * t82) * t56) * t102) * qJD(1)) * t119, t353 * t18 + t19 * t5, ((t82 * t248 + (-t157 * t82 + 0.2e1 * t350) * pkin(1)) * t220 + ((-g(3) * t55 + t154 * t53) * t56 + (t350 / 0.2e1 + t146 * t82) * t329) * t102) * t119, ((t83 * t248 + (-t157 * t83 - 0.2e1 * t351) * pkin(1)) * t220 + ((g(3) * t53 + t154 * t55) * t56 + (-t351 / 0.2e1 + t146 * t83) * t329) * t102) * t119;];
tau_reg = t1;
