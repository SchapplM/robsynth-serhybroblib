% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% palh3m1IC
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
% tau_reg [4x(10*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh3m1IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:31:47
% EndTime: 2020-04-20 17:32:02
% DurationCPUTime: 13.92s
% Computational Cost: add. (14092->695), mult. (28238->997), div. (234->6), fcn. (22733->36), ass. (0->376)
t304 = -qJ(7) + pkin(15);
t295 = -qJ(8) + t304;
t270 = sin(t295);
t194 = -t270 * pkin(7) + pkin(3) * sin(t304);
t273 = cos(t295);
t195 = -t273 * pkin(7) + pkin(3) * cos(t304);
t474 = qJ(3) + qJ(4);
t290 = pkin(14) + t474;
t265 = sin(t290);
t266 = cos(t290);
t151 = 0.1e1 / (t265 * t273 + t266 * t270) / pkin(9) / pkin(7);
t309 = qJ(2) + qJ(7);
t296 = pkin(16) + t309;
t271 = sin(t296);
t274 = cos(t296);
t313 = sin(qJ(6));
t320 = cos(qJ(6));
t159 = 0.1e1 / (-t271 * t320 + t274 * t313) / pkin(5) / pkin(2);
t317 = sin(qJ(2));
t515 = pkin(1) * t317;
t209 = -pkin(2) * t271 - t515;
t323 = cos(qJ(2));
t299 = t323 * pkin(1);
t210 = pkin(2) * t274 + t299;
t74 = (t209 * t320 + t210 * t313) * pkin(5) * t159;
t352 = t151 * t74;
t23 = (-t194 * t273 + t195 * t270) * pkin(7) * t352;
t558 = 0.1e1 + t23;
t316 = sin(qJ(3));
t204 = -pkin(4) * t316 - pkin(9) * t265;
t322 = cos(qJ(3));
t205 = pkin(4) * t322 + pkin(9) * t266;
t63 = (t204 * t273 - t205 * t270) * t151 * pkin(7);
t557 = t63 + 0.1e1;
t321 = cos(qJ(5));
t459 = qJD(5) * t321;
t466 = qJD(1) * t323;
t467 = qJD(1) * t317;
t190 = t316 * t467 - t322 * t466;
t203 = t316 * t323 + t317 * t322;
t192 = t203 * qJD(1);
t315 = sin(qJ(4));
t518 = cos(qJ(4));
t127 = t518 * t190 + t192 * t315;
t545 = t127 * t321;
t560 = t459 - t545;
t22 = (-t194 * t266 - t195 * t265) * pkin(9) * t352;
t519 = 0.1e1 - t74;
t559 = t22 + t519;
t303 = qJD(2) + qJD(3);
t502 = pkin(1) * qJD(2);
t443 = t322 * t502;
t212 = pkin(4) * t303 - t443;
t415 = t518 * t502;
t396 = t316 * t415;
t158 = -t315 * t212 + t396;
t428 = qJD(4) + t303;
t155 = pkin(10) * t428 - t158;
t314 = sin(qJ(5));
t276 = t299 + pkin(12);
t227 = t276 * qJD(1);
t156 = -pkin(4) * t190 - t227;
t374 = -t315 * t190 + t518 * t192;
t58 = -pkin(8) * t127 + pkin(10) * t374 + t156;
t48 = -t155 * t314 + t321 * t58;
t540 = qJD(5) - t127;
t556 = t540 * t48;
t49 = t155 * t321 + t314 * t58;
t555 = t540 * t49;
t553 = t127 * t374;
t319 = cos(qJ(7));
t434 = t319 * t466;
t312 = sin(qJ(7));
t436 = t312 * t467;
t188 = -t434 + t436;
t202 = t312 * t323 + t319 * t317;
t191 = t202 * qJD(1);
t311 = sin(pkin(15));
t503 = cos(pkin(15));
t516 = sin(qJ(8));
t517 = cos(qJ(8));
t200 = t311 * t516 + t503 * t517;
t349 = -t311 * t517 + t503 * t516;
t389 = -t188 * t349 + t191 * t200;
t390 = t188 * t200 + t191 * t349;
t552 = t389 * t390;
t397 = t321 * t428;
t108 = -t314 * t374 - t397;
t387 = t316 * t317 - t322 * t323;
t367 = t387 * qJD(3);
t332 = qJD(2) * t387 + t367;
t106 = qJD(1) * t332 - qJDD(1) * t203;
t333 = t303 * t203;
t107 = qJD(1) * t333 + qJDD(1) * t387;
t432 = qJD(4) * t518;
t463 = qJD(4) * t315;
t370 = t518 * t106 + t315 * t107 + t190 * t432 + t192 * t463;
t301 = qJDD(2) + qJDD(3);
t416 = qJDD(4) + t301;
t460 = qJD(5) * t314;
t29 = -qJD(5) * t397 - t314 * t416 - t321 * t370 - t374 * t460;
t110 = t314 * t428 - t321 * t374;
t462 = qJD(5) * t110;
t30 = t314 * t370 - t321 * t416 + t462;
t550 = -t560 * t108 - t321 * t29 - t30 * t314;
t26 = t29 * t314;
t549 = t560 * t110 - t26;
t423 = t315 * t106 - t518 * t107;
t45 = -qJD(4) * t374 + t423;
t43 = qJDD(5) + t45;
t548 = t110 * t374 + t43 * t314 + t560 * t540;
t298 = qJ(2) + t474;
t279 = cos(t298);
t259 = g(3) * t279;
t441 = qJD(3) * t502;
t262 = t316 * t441;
t451 = qJDD(2) * t322;
t168 = -pkin(1) * t451 + pkin(4) * t301 + t262;
t497 = pkin(1) * qJDD(2);
t542 = -t316 * t497 - t322 * t441;
t358 = -qJD(4) * t396 - t518 * t168 + t212 * t463 + t542 * t315;
t78 = -pkin(8) * t416 + t358;
t547 = t78 - t259;
t546 = t108 * t374;
t477 = t315 * t316;
t445 = pkin(1) * t477;
t244 = qJD(2) * t445;
t157 = t518 * t212 + t244;
t154 = -pkin(8) * t428 - t157;
t493 = t127 * t154;
t62 = (t204 * t266 + t205 * t265) * t151 * pkin(9);
t544 = t62 * t552;
t302 = qJD(7) + qJD(2);
t153 = t302 * t202;
t448 = t323 * qJDD(1);
t449 = t317 * qJDD(1);
t105 = qJD(1) * t153 + t312 * t449 - t319 * t448;
t543 = t302 * t191 - t105;
t492 = t127 * t314;
t541 = t460 - t492;
t318 = sin(qJ(1));
t324 = cos(qJ(1));
t230 = g(1) * t324 + g(2) * t318;
t539 = -t127 ^ 2 + t374 ^ 2;
t278 = sin(t298);
t478 = t314 * t324;
t479 = t314 * t318;
t346 = -t374 * t49 + t154 * t459 + t547 * t314 + (g(1) * t478 + g(2) * t479) * t278;
t408 = t154 * t460 + t321 * t259 + t374 * t48;
t258 = g(3) * t278;
t484 = t279 * t324;
t485 = t279 * t318;
t439 = g(1) * t484 + g(2) * t485 + t258;
t395 = t322 * t415;
t437 = t518 * t316;
t418 = pkin(1) * t437;
t82 = -qJD(3) * t395 + qJD(4) * t244 - qJDD(2) * t418 + t315 * t168 + t212 * t432;
t343 = -t156 * t127 - t439 - t82;
t537 = -t127 * t428 + t370;
t380 = t230 * t278;
t329 = t156 * t374 + t259 - t358 - t380;
t536 = -t303 * t374 - t423;
t73 = -pkin(8) * t374 - pkin(10) * t127;
t535 = t559 * t552;
t41 = t43 * t321;
t533 = t460 * t540 - t41;
t280 = -qJ(2) + t295;
t261 = cos(t280);
t532 = t230 * t261;
t529 = -g(1) * t318 + g(2) * t324;
t531 = t529 * t278;
t530 = t519 * t188 * t191;
t80 = -t192 * t303 + t107;
t310 = qJ(2) + qJ(3);
t292 = sin(t310);
t294 = cos(t310);
t348 = g(3) * t294 - t230 * t292;
t398 = -t314 * t48 + t321 * t49;
t528 = t316 ^ 2 + t322 ^ 2;
t186 = t200 * qJD(8);
t187 = t349 * qJD(8);
t514 = pkin(3) * t311;
t228 = pkin(1) * t312 + t514;
t442 = t503 * pkin(3);
t229 = t319 * pkin(1) + t442;
t372 = pkin(1) * (t200 * t312 + t319 * t349);
t526 = -qJD(7) * t372 - t186 * t228 - t187 * t229 + ((t186 * t311 + t503 * t187) * pkin(3) - qJD(2) * t372) * t74;
t371 = pkin(1) * (-t200 * t319 + t312 * t349);
t525 = -qJD(7) * t371 + t186 * t229 - t187 * t228 + ((-t503 * t186 + t187 * t311) * pkin(3) - qJD(2) * t371) * t74;
t289 = qJD(8) + t302;
t481 = t312 * t317;
t403 = t302 * t481;
t454 = qJD(1) * qJD(2);
t431 = t323 * t454;
t413 = -qJD(7) * t434 - t312 * t448 + (-t431 - t449) * t319;
t104 = qJD(1) * t403 + t413;
t31 = t200 * t104 + t105 * t349 + t186 * t188 + t187 * t191;
t524 = -t289 * t390 + t31;
t32 = -t104 * t349 + t200 * t105 + t186 * t191 - t187 * t188;
t523 = -t289 * t389 + t32;
t111 = -t227 + (t503 * t188 - t191 * t311) * pkin(3);
t260 = sin(t280);
t285 = t319 * t497;
t300 = qJDD(2) + qJDD(7);
t444 = t312 * t502;
t163 = -qJD(7) * t444 + t300 * t442 + t285;
t458 = qJD(7) * t319;
t164 = t300 * t514 + (qJD(2) * t458 + t312 * qJDD(2)) * pkin(1);
t196 = t302 * t514 + t444;
t197 = t302 * t442 + t319 * t502;
t50 = -t163 * t349 - t200 * t164 - t186 * t197 + t187 * t196;
t335 = g(3) * t260 - t111 * t390 - t50;
t51 = -t200 * t163 + t164 * t349 + t186 * t196 + t187 * t197;
t347 = g(3) * t261 + t111 * t389 + t230 * t260 + t51;
t522 = t389 ^ 2 - t390 ^ 2;
t513 = pkin(4) * t192;
t511 = pkin(8) * t278;
t506 = g(3) * t323;
t286 = pkin(1) * t467;
t180 = qJD(2) * t286 - qJDD(1) * t276;
t81 = -pkin(4) * t107 + t180;
t11 = pkin(8) * t45 - pkin(10) * t370 + t81;
t77 = pkin(10) * t416 + t82;
t7 = qJD(5) * t48 + t11 * t314 + t321 * t77;
t6 = t7 * t321;
t28 = t30 * t321;
t327 = qJD(1) ^ 2;
t71 = (-t209 * t274 - t210 * t271) * t159 * pkin(2);
t499 = t327 * t71;
t496 = t108 * t110;
t495 = t108 * t314;
t144 = -t518 * t203 + t315 * t387;
t490 = t144 * t314;
t489 = t144 * t321;
t487 = t190 * t192;
t483 = t302 * t188;
t480 = t313 * t320;
t476 = t318 * t321;
t475 = t321 * t324;
t297 = -qJ(2) + t304;
t275 = cos(t297);
t470 = pkin(3) * t275 + t299;
t305 = t313 ^ 2;
t307 = t320 ^ 2;
t469 = t305 - t307;
t306 = t317 ^ 2;
t308 = t323 ^ 2;
t468 = t306 - t308;
t465 = qJD(2) * t323;
t464 = qJD(3) * t303;
t461 = qJD(5) * t540;
t457 = qJDD(1) * pkin(6);
t456 = qJDD(1) * pkin(12);
t453 = qJD(1) * qJD(6);
t452 = qJDD(1) * t320;
t450 = qJDD(2) * pkin(1) ^ 2;
t287 = t317 * t502;
t440 = t317 * t327 * t323;
t430 = -pkin(4) * t294 + t299;
t215 = pkin(4) * t292 - t515;
t429 = -pkin(10) * t279 + t215;
t427 = t558 * t374;
t426 = t557 * t374;
t283 = -pkin(1) * t322 + pkin(4);
t182 = t315 * t283 - t418;
t176 = pkin(10) + t182;
t64 = -t513 + t73;
t59 = t286 + t64;
t424 = qJD(5) * t176 + t59;
t422 = t314 * t540;
t421 = 0.2e1 * pkin(6) * t453;
t420 = -0.2e1 * pkin(12) * t454;
t419 = t6 + t439;
t417 = t480 * t499;
t288 = qJDD(8) + t300;
t414 = pkin(4) * t432;
t412 = t317 * t431;
t411 = t453 * t480;
t373 = t315 * t322 + t437;
t178 = t373 * t502;
t409 = pkin(4) * t463 + t178;
t407 = pkin(8) * t279 + pkin(10) * t278;
t94 = t196 * t349 - t197 * t200;
t95 = t196 * t200 + t197 * t349;
t404 = t389 * t95 + t390 * t94;
t399 = t314 * t49 + t321 * t48;
t394 = -t176 * t43 - t493;
t281 = pkin(4) * t315 + pkin(10);
t393 = -t281 * t43 - t493;
t391 = t157 * t127 + t158 * t374;
t385 = t48 * t545 + t49 * t492 + t419;
t384 = t541 * t108 - t28;
t291 = sin(t309);
t293 = cos(t309);
t382 = g(3) * t291 - t227 * t188 + t230 * t293;
t381 = -qJD(5) * t58 - t258 - t77;
t357 = t518 * t387;
t60 = -qJD(4) * t357 - t203 * t463 - t315 * t333 - t518 * t332;
t378 = t144 * t459 - t314 * t60;
t377 = -t144 * t460 - t321 * t60;
t369 = -t413 + t483;
t368 = t227 * t203;
t181 = t518 * t283 + t445;
t365 = -pkin(6) * t327 + t230;
t364 = pkin(12) * t327 + t230;
t362 = t529 + 0.2e1 * t457;
t361 = -t529 + 0.2e1 * t456;
t359 = -t78 - t380;
t121 = (t188 * t311 + t503 * t191) * pkin(3);
t356 = -g(3) * t293 + t227 * t191 + t230 * t291 + t285;
t354 = t492 * t540 - t533 - t546;
t129 = -pkin(4) * t333 + t287;
t61 = t144 * qJD(4) + t315 * t332 - t518 * t333;
t21 = t61 * pkin(8) + t60 * pkin(10) + t129;
t143 = -t203 * t315 - t357;
t166 = -pkin(4) * t387 - t276;
t70 = t143 * pkin(8) - t144 * pkin(10) + t166;
t351 = qJD(5) * t144 * t154 + t21 * t540 + t43 * t70;
t350 = t144 * t78 - t154 * t60 - t461 * t70;
t345 = (t227 * t467 - t506) * pkin(1) + t230 * t515;
t344 = -qJD(7) * t302 + (-qJD(7) - t74 * (-qJD(7) + t302)) * qJD(2);
t342 = -t541 * t110 + t550;
t341 = -t192 * t227 + t262 + t348;
t10 = t321 * t11;
t8 = -qJD(5) * t49 - t314 * t77 + t10;
t339 = -qJD(5) * t399 - t8 * t314 + t6;
t336 = -g(3) * t292 + t190 * t227 - t230 * t294 - t542;
t326 = qJD(2) ^ 2;
t325 = qJD(6) ^ 2;
t282 = -t518 * pkin(4) - pkin(8);
t272 = sin(t297);
t218 = t324 * t511;
t217 = t318 * t511;
t211 = pkin(12) + t430;
t201 = -t319 * t323 + t481;
t179 = -t395 + t244;
t175 = -pkin(8) - t181;
t174 = t279 * t475 - t479;
t173 = t279 * t478 + t476;
t172 = t279 * t476 + t478;
t171 = t279 * t479 - t475;
t160 = t286 - t513;
t152 = -t319 * t465 - t323 * t458 + t403;
t135 = (-t200 * t311 - t349 * t503) * pkin(3);
t134 = (-t503 * t200 + t311 * t349) * pkin(3);
t133 = -t283 * t463 + (qJD(3) * t373 + t316 * t432) * pkin(1);
t132 = t283 * t432 + (t316 * t463 + (-t518 * t322 + t477) * qJD(3)) * pkin(1);
t122 = (t503 * t201 - t202 * t311) * pkin(3) - t276;
t115 = t286 + t121;
t114 = -t200 * t228 - t229 * t349;
t113 = -t200 * t229 + t228 * t349;
t112 = -t190 ^ 2 + t192 ^ 2;
t103 = -t200 * t202 + t201 * t349;
t102 = t200 * t201 + t202 * t349;
t88 = -t227 * t287 + (-t180 - t529) * t276;
t79 = -t190 * t303 + t106;
t76 = t287 + (t152 * t311 + t503 * t153) * pkin(3);
t57 = t157 * t321 + t314 * t73;
t56 = -t157 * t314 + t321 * t73;
t54 = t179 * t321 + t314 * t64;
t53 = -t179 * t314 + t321 * t64;
t52 = (t104 * t311 + t503 * t105) * pkin(3) + t180;
t47 = -t152 * t349 + t153 * t200 + t186 * t202 - t187 * t201;
t46 = t152 * t200 + t153 * t349 + t186 * t201 + t187 * t202;
t36 = -t158 * t428 + t329;
t35 = t157 * t428 + t343;
t17 = -t289 * t95 + t347;
t16 = t289 * t94 + t335 - t532;
t15 = t495 * t540 - t28;
t13 = -t422 * t540 + t41 - t546;
t5 = -t110 * t422 + t550;
t4 = -pkin(8) * t30 + t108 * t158 - t540 * t56 + (-pkin(10) * t43 - t493) * t314 + (-pkin(10) * t461 + t359) * t321 + t408;
t3 = pkin(8) * t29 + t533 * pkin(10) + t110 * t158 - t154 * t545 + t540 * t57 + t346;
t2 = -g(1) * t218 - g(2) * t217 + t154 * t158 - t48 * t56 - t49 * t57 - t547 * pkin(8) + (t230 * t279 + t258 + t339) * pkin(10);
t1 = t108 * t57 + t110 * t56 + (-t556 + (-t30 + t462) * pkin(10)) * t321 + (-t8 - t555 + (qJD(5) * t108 - t29) * pkin(10)) * t314 + t419;
t9 = [0, 0, 0, 0, 0, qJDD(1), -t529, t230, 0, 0, qJDD(1) * t306 + 0.2e1 * t412, 0.2e1 * t317 * t448 - 0.2e1 * t454 * t468, qJDD(2) * t317 + t323 * t326, qJDD(1) * t308 - 0.2e1 * t412, qJDD(2) * t323 - t317 * t326, 0, t317 * t420 + t323 * t361, -t317 * t361 + t323 * t420, -t230, (-t529 + t456) * pkin(12), -t106 * t203 - t192 * t332, t106 * t387 - t107 * t203 + t303 * (t190 * t387 - t192 * t203), -t301 * t203 + t303 * t332, t107 * t387 + t190 * t333, t301 * t387 + t303 * t333, 0, t276 * t107 - t180 * t387 + t529 * t294 + qJD(3) * t368 + (-t190 * t515 + t368) * qJD(2), -t276 * t106 - t180 * t203 - t529 * t292 - t227 * t367 + (-t192 * t515 - t227 * t387) * qJD(2), ((-t322 * t203 - t316 * t387) * qJDD(2) - t528 * qJD(2) * t465) * pkin(1) - t230, t88, t144 * t370 + t374 * t60, -t127 * t60 - t143 * t370 - t144 * t45 + t374 * t61, t144 * t416 - t428 * t60, -t127 * t61 + t143 * t45, -t143 * t416 - t428 * t61, 0, -t127 * t129 + t143 * t81 + t156 * t61 + t166 * t45 + t279 * t529, -t129 * t374 + t144 * t81 - t156 * t60 + t166 * t370 - t531, -t143 * t82 + t144 * t358 + t157 * t60 + t158 * t61 - t230, t129 * t156 + t166 * t81 - t211 * t529, t110 * t377 - t29 * t489, (t108 * t321 + t110 * t314) * t60 + (t26 - t28 + (-t110 * t321 + t495) * qJD(5)) * t144, t110 * t61 - t143 * t29 + t377 * t540 + t43 * t489, t108 * t378 + t30 * t490, -t108 * t61 - t143 * t30 - t378 * t540 - t43 * t490, t143 * t43 + t540 * t61, -g(1) * t172 + g(2) * t174 + t143 * t8 + t314 * t350 + t321 * t351 + t48 * t61, g(1) * t171 - g(2) * t173 - t143 * t7 - t314 * t351 + t321 * t350 - t49 * t61, t531 + (-t110 * t21 - t144 * t8 + t29 * t70 + t48 * t60 + (-t108 * t70 - t144 * t49) * qJD(5)) * t321 + (-t108 * t21 - t144 * t7 - t30 * t70 + t49 * t60 + (t110 * t70 + t144 * t48) * qJD(5)) * t314, t399 * t21 + (qJD(5) * t398 + t7 * t314 + t8 * t321) * t70 + t529 * (-t211 + t407), qJDD(1) * t305 + 0.2e1 * t411, 0.2e1 * t313 * t452 - 0.2e1 * t453 * t469, qJDD(6) * t313 + t320 * t325, qJDD(1) * t307 - 0.2e1 * t411, qJDD(6) * t320 - t313 * t325, 0, t313 * t421 - t320 * t362, t313 * t362 + t320 * t421, -t230, (t529 + t457) * pkin(6), -t104 * t202 - t152 * t191, t104 * t201 - t105 * t202 + t152 * t188 - t153 * t191, -t152 * t302 + t202 * t300, t105 * t201 + t153 * t188, -t153 * t302 - t201 * t300, 0, -t105 * t276 - t153 * t227 + t180 * t201 + t188 * t287 - t293 * t529, t104 * t276 + t152 * t227 + t180 * t202 + t191 * t287 + t291 * t529, ((-t201 * t312 - t202 * t319) * qJDD(2) + (t152 * t319 - t153 * t312 + (-t201 * t319 + t202 * t312) * qJD(7)) * qJD(2)) * pkin(1) - t230, t88, t103 * t31 - t389 * t46, t102 * t31 + t103 * t32 - t389 * t47 + t390 * t46, t103 * t288 + t289 * t46, t102 * t32 + t390 * t47, t102 * t288 + t289 * t47, 0, -t102 * t52 - t111 * t47 - t122 * t32 + t261 * t529 - t390 * t76, t103 * t52 + t111 * t46 + t122 * t31 + t260 * t529 - t389 * t76, t102 * t50 - t103 * t51 - t46 * t94 - t47 * t95 - t230, t111 * t76 + t122 * t52 - t529 * (pkin(12) + t470); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t440, t468 * t327, t449, t440, t448, qJDD(2), t317 * t364 - t506, g(3) * t317 + t323 * t364, 0, 0, t487, t112, t79, -t487, t80, t301, (t190 * t467 + t316 * t464 + (-qJDD(2) - t301) * t322) * pkin(1) + t341, (t192 * t467 + t316 * t301 + t322 * t464) * pkin(1) + t336, -t190 * t443 + ((-qJD(3) * t190 + t106) * t322 - t80 * t316) * pkin(1), t450 * t528 + t345, t558 * t553, t558 * t539, t558 * t537, -t427 * t127, t558 * t536, t558 * t416, t127 * t160 + t133 * t428 + t181 * t416 + t23 * t36 + t329, -t132 * t428 + t160 * t374 - t182 * t416 + t23 * t35 + t343, t127 * t132 + t133 * t374 - t181 * t370 - t182 * t45 + t391, -g(3) * t430 - t158 * t132 + t157 * t133 - t156 * t160 - t181 * t358 + t82 * t182 - t215 * t230, t558 * t549, t23 * t5 + t342, t558 * t548, t15 * t23 + t384, t13 * t23 + t354, t427 * t540, -t108 * t133 + t175 * t30 + t23 * t4 + (-t132 * t540 + t394) * t314 + (-t424 * t540 + t359) * t321 + t408, -t110 * t133 - t175 * t29 + t23 * t3 + t394 * t321 + (-t132 * t321 + t314 * t424) * t540 + t346, t1 * t23 + (-t108 * t132 + t110 * t59 - t176 * t30 + (t110 * t176 - t48) * qJD(5)) * t321 + (t108 * t59 + t110 * t132 - t176 * t29 - t8 + (t108 * t176 - t49) * qJD(5)) * t314 + t385, t78 * t175 - t154 * t133 - g(3) * (-t407 + t430) - g(2) * (t318 * t429 + t217) - g(1) * (t324 * t429 + t218) + t23 * t2 - t399 * t59 + t398 * t132 + t339 * t176, -t417, t469 * t499, t71 * t313 * qJDD(1), t417, t71 * t452, t71 * qJDD(6), t71 * (-g(3) * t320 + t313 * t365), t71 * (g(3) * t313 + t320 * t365), 0, 0, t530, t519 * (-t188 ^ 2 + t191 ^ 2), -t519 * t436 * t302 - t74 * t369 + t369, -t530, t543 * t519, t519 * t300, -t74 * t356 + (-t188 * t467 + t319 * t300 + t312 * t344) * pkin(1) + t356, -t74 * t382 + (-t191 * t467 + (qJDD(2) * t74 - qJDD(2) - t300) * t312 + t344 * t319) * pkin(1) + t382, ((t104 - t483) * t319 + t543 * t312) * pkin(1), (t312 ^ 2 + t319 ^ 2) * t450 + t345, t535, t559 * t522, t559 * t524, -t535, t559 * t523, t559 * t288, t113 * t288 + t115 * t390 - t74 * (t121 * t390 + t134 * t288 + t347) + t22 * t17 - t526 * t289 + t347, -t114 * t288 + t115 * t389 - t74 * (t121 * t389 - t135 * t288 + t335) + t22 * t16 + t525 * t289 + t335 - t519 * t532, -t113 * t31 + t114 * t32 - t74 * (-t134 * t31 + t135 * t32 + t404) - t526 * t389 - t525 * t390 + t404, t50 * t114 + t51 * t113 - t111 * t115 - g(3) * t470 + t525 * t95 - t526 * t94 - t230 * (pkin(3) * t272 - t515) - (-t111 * t121 + t134 * t51 + t135 * t50 + (-g(3) * t275 - t230 * t272) * pkin(3)) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, t112, t79, -t487, t80, t301, (-qJD(2) * t303 * t316 - t451) * pkin(1) + t341, -t303 * t443 + t336, 0, 0, t557 * t553, t557 * t539, t557 * t537, -t426 * t127, t557 * t536, t557 * t416, -t178 * t428 + t63 * t36 + (-t127 * t192 + t518 * t416 - t428 * t463) * pkin(4) + t329, t179 * t428 + t63 * t35 + (-t192 * t374 - t315 * t416 - t428 * t432) * pkin(4) + t343, -t179 * t127 - t178 * t374 + (-t518 * t370 - t315 * t45 + (t127 * t518 - t315 * t374) * qJD(4)) * pkin(4) + t391, -t157 * t178 + t158 * t179 + (-t518 * t358 + t156 * t192 + t315 * t82 + (-t157 * t315 - t518 * t158) * qJD(4) + t348) * pkin(4), t557 * t549, t5 * t63 + t342, t557 * t548, t15 * t63 + t384, t13 * t63 + t354, t426 * t540, -t53 * t540 + t282 * t30 + t63 * t4 + t409 * t108 + (-t414 * t540 + t393) * t314 + (-t281 * t461 + t359) * t321 + t408, -t282 * t29 + t63 * t3 + t393 * t321 + t409 * t110 + (t281 * t460 - t321 * t414 + t54) * t540 + t346, t63 * t1 + t54 * t108 + t53 * t110 + (-t108 * t414 - t281 * t30 + (t110 * t281 - t48) * qJD(5)) * t321 + (t110 * t414 - t281 * t29 - t8 + (t108 * t281 - t49) * qJD(5)) * t314 + t385, t78 * t282 - t49 * t54 - t48 * t53 + t154 * t178 + g(3) * t407 - g(2) * (-pkin(10) * t485 + t217) - g(1) * (-pkin(10) * t484 + t218) + t63 * t2 + ((t154 * t315 + t398 * t518) * qJD(4) + t348) * pkin(4) + t339 * t281, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t544, t62 * t522, t62 * t524, -t544, t62 * t523, t62 * t288, t62 * t17, t62 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t496, -t108 ^ 2 + t110 ^ 2, t108 * t540 - t29, -t496, t110 * t540 - t30, t43, -g(1) * t173 - g(2) * t171 - t110 * t154 - t155 * t459 + t314 * t381 + t10 + t555, -g(1) * t174 - g(2) * t172 + t108 * t154 + t556 + (qJD(5) * t155 - t11) * t314 + t381 * t321, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t9;