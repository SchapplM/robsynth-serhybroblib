% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% tauJ_reg [10x(10*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = palh3m2OL_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_reg2_snew_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_reg2_snew_vp: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_reg2_snew_vp: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynJ_fixb_reg2_snew_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:40:41
% EndTime: 2020-05-07 04:41:48
% DurationCPUTime: 9.74s
% Computational Cost: add. (20705->555), mult. (46268->811), div. (0->0), fcn. (36912->18), ass. (0->364)
t342 = cos(qJ(3));
t379 = qJDD(2) + qJDD(3);
t334 = sin(qJ(3));
t343 = cos(qJ(2));
t387 = qJD(1) * t343;
t335 = sin(qJ(2));
t388 = qJD(1) * t335;
t276 = t334 * t388 - t342 * t387;
t389 = t335 * t342;
t278 = (-t343 * t334 - t389) * qJD(1);
t404 = t276 * t278;
t433 = t379 + t404;
t442 = t342 * t433;
t346 = qJD(1) ^ 2;
t336 = sin(qJ(1));
t424 = cos(qJ(1));
t353 = t424 * g(1) + t336 * g(2);
t288 = -t346 * pkin(12) - t353;
t261 = g(3) * t343 + t288 * t335;
t301 = t343 * t346 * t335;
t381 = qJDD(2) + t301;
t251 = t381 * pkin(1) - t261;
t263 = -g(3) * t335 + t288 * t343;
t318 = t343 ^ 2 * t346;
t345 = qJD(2) ^ 2;
t254 = (-t345 - t318) * pkin(1) + t263;
t198 = t334 * t251 + t342 * t254;
t271 = t276 ^ 2;
t324 = qJD(2) + qJD(3);
t321 = t324 ^ 2;
t426 = -t271 - t321;
t185 = t426 * pkin(4) - t198;
t333 = sin(qJ(4));
t341 = cos(qJ(4));
t196 = -t342 * t251 + t334 * t254;
t348 = pkin(4) * t433 + t196;
t132 = t333 * t185 - t341 * t348;
t133 = t341 * t185 + t333 * t348;
t65 = t132 * t341 - t133 * t333;
t441 = t342 * t65;
t329 = sin(pkin(15));
t337 = cos(qJ(8));
t417 = cos(pkin(15));
t423 = sin(qJ(8));
t284 = t329 * t337 - t417 * t423;
t330 = sin(qJ(7));
t338 = cos(qJ(7));
t274 = t330 * t388 - t338 * t387;
t390 = t335 * t338;
t277 = (t343 * t330 + t390) * qJD(1);
t285 = -t329 * t423 - t417 * t337;
t216 = -t274 * t285 - t277 * t284;
t217 = -t274 * t284 + t277 * t285;
t174 = t216 * t217;
t322 = qJDD(2) + qJDD(7);
t307 = qJDD(8) + t322;
t429 = t174 + t307;
t440 = t284 * t429;
t439 = t285 * t429;
t332 = sin(qJ(5));
t241 = t276 * t333 + t278 * t341;
t312 = t335 * qJDD(1);
t383 = qJD(1) * qJD(2);
t371 = t343 * t383;
t292 = t312 + t371;
t314 = t343 * qJDD(1);
t372 = t335 * t383;
t294 = t314 - t372;
t363 = t334 * t292 - t294 * t342;
t226 = -qJD(3) * t278 + t363;
t352 = -t276 * qJD(3) + t292 * t342 + t334 * t294;
t364 = -t341 * t226 - t333 * t352;
t157 = -qJD(4) * t241 - t364;
t156 = qJDD(5) - t157;
t310 = qJD(4) + t324;
t340 = cos(qJ(5));
t221 = t241 * t332 - t340 * t310;
t223 = t241 * t340 + t310 * t332;
t179 = t223 * t221;
t430 = t156 - t179;
t438 = t332 * t430;
t239 = -t341 * t276 + t278 * t333;
t192 = t241 * t239;
t308 = qJDD(4) + t379;
t428 = -t192 + t308;
t437 = t333 * t428;
t436 = t340 * t430;
t435 = t341 * t428;
t396 = t324 * t278;
t202 = t226 - t396;
t434 = t294 - t372;
t355 = t226 * t333 - t341 * t352;
t158 = -qJD(4) * t239 + t355;
t233 = t310 * t239;
t141 = t158 - t233;
t370 = t336 * g(1) - t424 * g(2);
t286 = qJDD(1) * pkin(12) + t370;
t250 = t434 * pkin(1) + t286;
t177 = t202 * pkin(4) + t250;
t54 = -t141 * pkin(10) + (t241 * t310 - t157) * pkin(8) - t177;
t190 = pkin(8) * t239 - pkin(10) * t241;
t425 = t310 ^ 2;
t86 = -t425 * pkin(8) + t308 * pkin(10) - t239 * t190 + t133;
t21 = t332 * t86 - t340 * t54;
t22 = t332 * t54 + t340 * t86;
t14 = t332 * t21 + t340 * t22;
t252 = t277 * t274;
t242 = -t252 + t322;
t432 = t242 * t338;
t362 = t292 * t330 - t338 * t294;
t225 = -qJD(7) * t277 - t362;
t227 = -t274 * qJD(7) + t338 * t292 + t330 * t294;
t146 = -t217 * qJD(8) + t225 * t285 - t227 * t284;
t323 = qJD(2) + qJD(7);
t309 = qJD(8) + t323;
t212 = t309 * t217;
t431 = t146 + t212;
t265 = t323 * t274;
t427 = -t265 + t227;
t236 = qJD(5) + t239;
t365 = t158 * t332 - t340 * t308;
t96 = (qJD(5) - t236) * t223 + t365;
t214 = t216 ^ 2;
t215 = t217 ^ 2;
t219 = t221 ^ 2;
t220 = t223 ^ 2;
t235 = t236 ^ 2;
t237 = t239 ^ 2;
t238 = t241 ^ 2;
t270 = t274 ^ 2;
t272 = t277 ^ 2;
t273 = t278 ^ 2;
t306 = t309 ^ 2;
t320 = t323 ^ 2;
t422 = pkin(3) * t329;
t421 = pkin(8) * t333;
t376 = pkin(3) * t417;
t195 = -t338 * t251 + t330 * t254;
t234 = (t417 * t274 - t277 * t329) * pkin(3);
t161 = -t277 * t234 + (t320 * t329 + t417 * t322) * pkin(3) - t195;
t197 = t330 * t251 + t338 * t254;
t162 = -t274 * t234 + (-t417 * t320 + t322 * t329) * pkin(3) + t197;
t88 = t285 * t161 - t284 * t162;
t89 = t284 * t161 + t285 * t162;
t39 = t284 * t89 + t285 * t88;
t40 = -t284 * t88 + t285 * t89;
t420 = t39 * t376 + t40 * t422;
t356 = t225 * t284 + t227 * t285;
t147 = qJD(8) * t216 + t356;
t211 = t309 * t216;
t116 = t147 - t211;
t50 = -t116 * t285 + t284 * t431;
t52 = t116 * t284 + t285 * t431;
t419 = t50 * t376 + t52 * t422;
t85 = -t308 * pkin(8) - t425 * pkin(10) + t190 * t241 + t132;
t418 = -pkin(8) * t85 + pkin(10) * t14;
t82 = t332 * t85;
t83 = t340 * t85;
t119 = t156 + t179;
t416 = t119 * t332;
t415 = t119 * t340;
t171 = -t174 + t307;
t414 = t171 * t284;
t413 = t171 * t285;
t412 = t177 * t333;
t411 = t177 * t341;
t410 = t177 * t342;
t188 = t192 + t308;
t409 = t188 * t333;
t408 = t188 * t341;
t407 = t236 * t332;
t406 = t236 * t340;
t405 = t250 * t330;
t403 = t276 * t324;
t148 = -t277 * t323 * t376 + (t417 * t225 + t427 * t329) * pkin(3) + t250;
t402 = t284 * t148;
t401 = t285 * t148;
t400 = t310 * t333;
t399 = t310 * t341;
t398 = t323 * t330;
t397 = t323 * t338;
t395 = t324 * t334;
t394 = t324 * t342;
t243 = t252 + t322;
t393 = t330 * t243;
t245 = -t404 + t379;
t392 = t334 * t245;
t391 = t335 * t334;
t386 = qJD(4) + t310;
t384 = qJD(5) + t236;
t382 = qJD(1) * qJD(6);
t9 = t14 * t333 - t341 * t85;
t380 = pkin(4) * t9 + t418;
t359 = -t158 * t340 - t308 * t332;
t101 = t384 * t221 + t359;
t168 = -t220 - t235;
t71 = -t168 * t332 - t415;
t378 = pkin(8) * t101 + pkin(10) * t71 + t82;
t159 = -t235 - t219;
t68 = t159 * t340 - t438;
t97 = -t384 * t223 - t365;
t377 = pkin(8) * t97 + pkin(10) * t68 - t83;
t375 = t333 * t179;
t374 = t341 * t179;
t373 = -pkin(8) * t341 - pkin(4);
t32 = t101 * t341 + t333 * t71;
t369 = pkin(4) * t32 + t378;
t30 = t333 * t68 + t341 * t97;
t368 = pkin(4) * t30 + t377;
t155 = t219 + t220;
t130 = -qJD(5) * t221 - t359;
t184 = t236 * t221;
t100 = t130 + t184;
t45 = t100 * t332 - t340 * t96;
t367 = pkin(8) * t155 + pkin(10) * t45 + t14;
t366 = t132 * t333 + t341 * t133;
t24 = t155 * t341 + t333 * t45;
t361 = pkin(4) * t24 + t367;
t186 = -t425 - t237;
t151 = t186 * t333 + t435;
t360 = pkin(4) * t151 - t132;
t331 = sin(qJ(6));
t311 = t331 * qJDD(1);
t339 = cos(qJ(6));
t291 = 0.2e1 * t339 * t382 + t311;
t313 = t339 * qJDD(1);
t293 = -0.2e1 * t331 * t382 + t313;
t13 = -t21 * t340 + t22 * t332;
t358 = t195 * t338 - t197 * t330;
t357 = t196 * t342 - t198 * t334;
t169 = -t306 - t214;
t110 = t169 * t284 + t439;
t111 = t169 * t285 - t440;
t351 = t110 * t376 + t111 * t422 + t88;
t194 = -t215 - t306;
t123 = t194 * t285 - t414;
t124 = -t194 * t284 - t413;
t350 = t123 * t376 + t124 * t422 - t89;
t349 = (-qJD(4) + t310) * t241 - t364;
t200 = (-qJD(7) + t323) * t277 - t362;
t208 = -t352 + t403;
t224 = -t238 - t425;
t163 = t224 * t341 - t409;
t347 = pkin(4) * t163 - t133;
t344 = qJD(6) ^ 2;
t317 = t339 ^ 2 * t346;
t316 = t335 ^ 2 * t346;
t315 = t331 ^ 2 * t346;
t299 = t339 * t346 * t331;
t295 = t314 - 0.2e1 * t372;
t290 = -t312 - 0.2e1 * t371;
t289 = t346 * pkin(6) - t353;
t287 = -qJDD(1) * pkin(6) + t370;
t262 = -g(3) * t331 + t289 * t339;
t260 = g(3) * t339 + t289 * t331;
t259 = -t273 + t321;
t258 = -t272 + t320;
t257 = t271 - t321;
t256 = t270 - t320;
t249 = t273 - t271;
t248 = t272 - t270;
t232 = -t238 + t425;
t231 = t237 - t425;
t230 = -t271 - t273;
t229 = -t270 - t272;
t210 = -t215 + t306;
t209 = t214 - t306;
t207 = t352 + t403;
t205 = t265 + t227;
t201 = t226 + t396;
t199 = (qJD(7) + t323) * t277 + t362;
t193 = (t343 * pkin(1) + pkin(12)) * t250;
t191 = t238 - t237;
t183 = -t220 + t235;
t182 = t219 - t235;
t181 = (-t239 * t341 + t241 * t333) * t310;
t180 = (-t239 * t333 - t241 * t341) * t310;
t178 = t220 - t219;
t175 = -t237 - t238;
t173 = t215 - t214;
t167 = t231 * t341 - t409;
t166 = -t232 * t333 + t435;
t165 = t231 * t333 + t408;
t164 = t232 * t341 + t437;
t153 = (t216 * t285 + t217 * t284) * t309;
t152 = (t216 * t284 - t217 * t285) * t309;
t149 = -t214 - t215;
t145 = (-t221 * t340 + t223 * t332) * t236;
t144 = (-t221 * t332 - t223 * t340) * t236;
t143 = -t386 * t239 + t355;
t142 = t158 + t233;
t138 = t386 * t241 + t364;
t137 = t158 * t341 - t241 * t400;
t136 = t158 * t333 + t241 * t399;
t135 = -t157 * t333 + t239 * t399;
t134 = t157 * t341 + t239 * t400;
t129 = -qJD(5) * t223 - t365;
t128 = t209 * t285 - t414;
t127 = -t210 * t284 + t439;
t126 = t209 * t284 + t413;
t125 = t210 * t285 + t440;
t117 = (qJD(8) + t309) * t216 + t356;
t115 = t147 + t211;
t113 = t146 - t212;
t109 = t147 * t285 - t284 * t212;
t108 = t147 * t284 + t285 * t212;
t107 = -t146 * t284 - t285 * t211;
t106 = t146 * t285 - t284 * t211;
t103 = -pkin(4) * t138 + t411;
t102 = -pkin(4) * t143 - t412;
t99 = t130 - t184;
t93 = t130 * t340 - t223 * t407;
t92 = t130 * t332 + t223 * t406;
t91 = -t129 * t332 + t221 * t406;
t90 = t129 * t340 + t221 * t407;
t81 = t145 * t341 + t156 * t333;
t80 = t145 * t333 - t156 * t341;
t79 = t182 * t340 - t416;
t78 = -t183 * t332 + t436;
t77 = t182 * t332 + t415;
t76 = t183 * t340 + t438;
t75 = -t138 * t341 - t141 * t333;
t74 = -t142 * t341 + t333 * t349;
t73 = -t138 * t333 + t141 * t341;
t72 = pkin(4) * t74;
t70 = t168 * t340 - t416;
t67 = t159 * t332 + t436;
t63 = t341 * t93 + t375;
t62 = t341 * t91 - t375;
t61 = t333 * t93 - t374;
t60 = t333 * t91 + t374;
t59 = pkin(4) * t65;
t58 = -t117 * t422 - t401;
t57 = t113 * t376 + t401;
t56 = t113 * t422 - t402;
t55 = -t117 * t376 - t402;
t51 = t113 * t285 - t115 * t284;
t49 = t113 * t284 + t115 * t285;
t46 = -pkin(4) * t175 + t366;
t44 = -t332 * t99 + t340 * t97;
t43 = -t100 * t340 - t332 * t96;
t42 = t332 * t97 + t340 * t99;
t36 = -t333 * t96 + t341 * t79;
t35 = t100 * t333 + t341 * t78;
t34 = t333 * t79 + t341 * t96;
t33 = -t100 * t341 + t333 * t78;
t28 = t178 * t333 + t341 * t44;
t27 = -t178 * t341 + t333 * t44;
t26 = -pkin(10) * t70 + t83;
t25 = -pkin(10) * t67 + t82;
t20 = -t149 * t422 - t39;
t19 = -t149 * t376 + t40;
t16 = -pkin(8) * t70 + t22;
t15 = -pkin(8) * t67 + t21;
t11 = -t16 * t333 + t26 * t341;
t10 = -t15 * t333 + t25 * t341;
t7 = -pkin(10) * t43 - t13;
t6 = -pkin(4) * t70 + t16 * t341 + t26 * t333;
t5 = -pkin(4) * t67 + t15 * t341 + t25 * t333;
t4 = t341 * t7 + t43 * t421;
t3 = (-pkin(10) * t341 + t421) * t13;
t2 = t333 * t7 + t373 * t43;
t1 = (-pkin(10) * t333 + t373) * t13;
t8 = [0, 0, 0, 0, 0, qJDD(1), t370, t353, 0, 0, (t292 + t371) * t335, -t290 * t343 + t295 * t335, t335 * t381 + t343 * (-t316 + t345), t434 * t343, t335 * (t318 - t345) + t343 * (qJDD(2) - t301), 0, pkin(12) * t295 + t286 * t343, pkin(12) * t290 - t286 * t335, t335 * t261 + t343 * t263 + pkin(12) * (t318 + t316), pkin(12) * t286, t335 * (t278 * t395 + t342 * t352) + t343 * (-t278 * t394 + t334 * t352), t335 * (-t202 * t342 + t208 * t334) + t343 * (-t202 * t334 - t208 * t342), t335 * (t259 * t334 - t442) + t343 * (-t259 * t342 - t334 * t433), t335 * (t226 * t334 + t276 * t394) + t343 * (-t226 * t342 + t276 * t395), t335 * (-t257 * t342 + t392) + t343 * (-t245 * t342 - t257 * t334), (t335 * (-t276 * t342 - t278 * t334) + t343 * (-t276 * t334 + t278 * t342)) * t324, t250 * t391 + t343 * (pkin(1) * t202 - t250 * t342) + pkin(12) * t202, t250 * t389 + t343 * (-pkin(1) * t208 + t250 * t334) - pkin(12) * t208, t335 * t357 + t343 * (-pkin(1) * t230 + t196 * t334 + t198 * t342) - pkin(12) * t230, t193, t335 * (t136 * t334 - t137 * t342) + t343 * (-t136 * t342 - t137 * t334), t335 * (t334 * t73 - t342 * t75) + t343 * (-t334 * t75 - t342 * t73), t335 * (t164 * t334 - t166 * t342) + t343 * (-t164 * t342 - t166 * t334), t335 * (t134 * t334 - t135 * t342) + t343 * (-t134 * t342 - t135 * t334), t335 * (t165 * t334 - t167 * t342) + t343 * (-t165 * t342 - t167 * t334), t335 * (t180 * t334 - t181 * t342) + t343 * (-t180 * t342 - t181 * t334), t335 * (t103 * t334 + t333 * t410) + t343 * (-pkin(1) * t138 - t103 * t342 + t334 * t412) - pkin(12) * t138, t335 * (t102 * t334 + t341 * t410) + t343 * (-pkin(1) * t143 - t102 * t342 + t334 * t411) - pkin(12) * t143, t335 * (t334 * t46 - t441) + t343 * (-pkin(1) * t175 - t334 * t65 - t342 * t46) - pkin(12) * t175, (pkin(4) * t391 + t343 * (-pkin(4) * t342 + pkin(1)) + pkin(12)) * t177, t335 * (t334 * t61 - t342 * t63) + t343 * (-t334 * t63 - t342 * t61), t335 * (t27 * t334 - t28 * t342) + t343 * (-t27 * t342 - t28 * t334), t335 * (t33 * t334 - t342 * t35) + t343 * (-t33 * t342 - t334 * t35), t335 * (t334 * t60 - t342 * t62) + t343 * (-t334 * t62 - t342 * t60), t335 * (t334 * t34 - t342 * t36) + t343 * (-t334 * t36 - t34 * t342), t335 * (t334 * t80 - t342 * t81) + t343 * (-t334 * t81 - t342 * t80), t335 * (-t10 * t342 + t334 * t5) + t343 * (-pkin(1) * t67 - t10 * t334 - t342 * t5) - pkin(12) * t67, t335 * (-t11 * t342 + t334 * t6) + t343 * (-pkin(1) * t70 - t11 * t334 - t342 * t6) - pkin(12) * t70, t335 * (t2 * t334 - t342 * t4) + t343 * (-pkin(1) * t43 - t2 * t342 - t334 * t4) - pkin(12) * t43, t335 * (t1 * t334 - t3 * t342) + t343 * (-pkin(1) * t13 - t1 * t342 - t3 * t334) - pkin(12) * t13, t291 * t331, t291 * t339 + t293 * t331, t331 * (qJDD(6) + t299) + t339 * (-t315 + t344), t293 * t339, t331 * (t317 - t344) + t339 * (qJDD(6) - t299), 0, -pkin(6) * t293 + t287 * t339, pkin(6) * t291 - t287 * t331, t331 * t260 + t339 * t262 - pkin(6) * (t315 + t317), -pkin(6) * t287, t335 * (t227 * t338 - t277 * t398) + t343 * (t227 * t330 + t277 * t397), t335 * (-t199 * t338 - t330 * t427) + t343 * (-t199 * t330 + t338 * t427), t335 * (-t258 * t330 + t432) + t343 * (t242 * t330 + t258 * t338), t335 * (-t225 * t330 + t274 * t397) + t343 * (t225 * t338 + t274 * t398), t335 * (t256 * t338 - t393) + t343 * (t243 * t338 + t256 * t330), (t335 * (-t274 * t338 + t277 * t330) + t343 * (-t274 * t330 - t277 * t338)) * t323, -t335 * t405 + t343 * (-pkin(1) * t199 + t250 * t338) - pkin(12) * t199, -t250 * t390 + t343 * (-pkin(1) * t427 - t405) - pkin(12) * t427, t335 * t358 + t343 * (-pkin(1) * t229 + t195 * t330 + t197 * t338) - pkin(12) * t229, t193, t335 * (-t108 * t330 + t109 * t338) + t343 * (t108 * t338 + t109 * t330), t335 * (-t330 * t49 + t338 * t51) + t343 * (t330 * t51 + t338 * t49), t335 * (-t125 * t330 + t127 * t338) + t343 * (t125 * t338 + t127 * t330), t335 * (-t106 * t330 + t107 * t338) + t343 * (t106 * t338 + t107 * t330), t335 * (-t126 * t330 + t128 * t338) + t343 * (t126 * t338 + t128 * t330), t335 * (-t152 * t330 + t153 * t338) + t343 * (t152 * t338 + t153 * t330), t335 * (-t330 * t57 + t338 * t56) + t343 * (pkin(1) * t113 + t330 * t56 + t338 * t57) + pkin(12) * t113, t335 * (-t330 * t55 + t338 * t58) + t343 * (-pkin(1) * t117 + t330 * t58 + t338 * t55) - pkin(12) * t117, t335 * (-t19 * t330 + t20 * t338) + t343 * (-pkin(1) * t149 + t19 * t338 + t20 * t330) - pkin(12) * t149, (t335 * (t329 * t338 - t417 * t330) * pkin(3) + t343 * (pkin(1) + (t329 * t330 + t417 * t338) * pkin(3)) + pkin(12)) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, -t318 + t316, t312, t301, t314, qJDD(2), -t261, -t263, 0, 0, -t404, t249, -t207, t404, t201, t379, pkin(1) * (-t334 * t426 - t442) + t196, pkin(1) * (t392 - t342 * (-t273 - t321)) + t198, pkin(1) * (-t334 * (-(qJD(3) - t324) * t278 + t363) - t342 * t207), -pkin(1) * t357, t192, t191, t142, -t192, t349, t308, pkin(1) * (-t334 * (t186 * t341 - t437) - t342 * t151) + t360, pkin(1) * (-t334 * (-t224 * t333 - t408) - t342 * t163) + t347, pkin(1) * (-t334 * (t142 * t333 + t341 * t349) - t342 * t74) + t72, pkin(1) * (-t334 * t366 + t441) - t59, t92, t42, t76, t90, t77, t144, pkin(1) * (-t334 * (-t333 * t97 + t341 * t68) - t342 * t30) + t368, pkin(1) * (-t334 * (-t101 * t333 + t341 * t71) - t342 * t32) + t369, pkin(1) * (-t334 * (-t155 * t333 + t341 * t45) - t342 * t24) + t361, pkin(1) * (-t334 * (t14 * t341 + t333 * t85) - t342 * t9) + t380, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t248, t205, -t252, t200, t322, pkin(1) * (t330 * (-t320 - t270) + t432) - t195, pkin(1) * (-t393 + t338 * (-t272 - t320)) - t197, pkin(1) * (t200 * t330 - t338 * t205), -pkin(1) * t358, -t174, t173, t116, t174, t431, t307, pkin(1) * (t110 * t338 + t111 * t330) + t351, pkin(1) * (t123 * t338 + t124 * t330) + t350, pkin(1) * (t330 * t52 + t338 * t50) + t419, pkin(1) * (t330 * t40 + t338 * t39) + t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t404, t249, -t207, t404, t201, t379, t196, t198, 0, 0, t192, t191, t142, -t192, t349, t308, t360, t347, t72, -t59, t92, t42, t76, t90, t77, t144, t368, t369, t361, t380, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t191, t142, -t192, t349, t308, -t132, -t133, 0, 0, t92, t42, t76, t90, t77, t144, t377, t378, t367, t418, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t178, t100, -t179, -t96, t156, -t21, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t315 - t317, t311, t299, t313, qJDD(6), -t260, -t262, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t248, t205, -t252, t200, t322, -t195, -t197, 0, 0, -t174, t173, t116, t174, t431, t307, t351, t350, t419, t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t173, t116, t174, t431, t307, t88, -t89, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauJ_reg = t8;
