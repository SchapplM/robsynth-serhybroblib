% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% palh3m2IC
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
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = palh3m2IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:59:53
% EndTime: 2020-05-07 05:00:09
% DurationCPUTime: 14.23s
% Computational Cost: add. (12895->740), mult. (26324->1006), div. (204->8), fcn. (21348->49), ass. (0->388)
t317 = 0.1e1 / pkin(9);
t459 = (qJ(4) + pkin(14));
t437 = (qJ(3) + t459);
t414 = (pkin(15) + t437);
t489 = (qJ(2) - qJ(6));
t369 = (t414 - t489);
t342 = -2 * qJ(7) - pkin(16) + t369;
t370 = (t414 + t489);
t351 = pkin(16) + t370;
t544 = -cos((qJ(8) - t342)) + cos((qJ(8) - t351));
t150 = 0.1e1 / t544;
t320 = 0.1e1 / pkin(2);
t467 = qJ(7) + pkin(16);
t419 = t467 + t489;
t538 = pkin(3) * (pkin(1) * (cos((qJ(8) - t489)) - cos((qJ(8) + t489))) + (cos((qJ(8) - t419)) - cos((qJ(8) + t419))) * pkin(2)) * t150 * t320;
t404 = t317 * t538;
t574 = t404 + 0.1e1;
t468 = pkin(15) - qJ(7);
t441 = -qJ(8) + t468;
t485 = cos((2 * t437)) - cos((2 * t441));
t547 = (t485 * pkin(9) + (cos((2 * qJ(3) + t459)) - cos((-2 * qJ(7) + 2 * pkin(15) - 2 * qJ(8) - t459))) * pkin(4)) / t485;
t450 = t317 * t547;
t573 = t450 - 0.1e1;
t309 = cos(qJ(5));
t473 = qJD(5) * t309;
t304 = sin(qJ(3));
t310 = cos(qJ(3));
t311 = cos(qJ(2));
t481 = qJD(1) * t311;
t305 = sin(qJ(2));
t482 = qJD(1) * t305;
t186 = t304 * t482 - t310 * t481;
t199 = t304 * t311 + t305 * t310;
t188 = t199 * qJD(1);
t303 = sin(qJ(4));
t532 = cos(qJ(4));
t121 = t532 * t186 + t188 * t303;
t560 = t121 * t309;
t576 = t473 - t560;
t318 = 0.1e1 / pkin(7);
t354 = -qJ(7) + t370;
t355 = -qJ(7) + t369;
t453 = t150 * t318 * (((cos(t355) - cos(t354)) * pkin(1) + (cos(t342) - cos(t351)) * pkin(2)) * pkin(3) + ((-cos((qJ(8) - t355)) + cos((qJ(8) - t354))) * pkin(1) + t544 * pkin(2)) * pkin(7));
t256 = sin(t419);
t241 = 0.1e1 / t256;
t504 = (-pkin(2) * t256 - pkin(1) * sin(t489)) * t241;
t575 = (t453 + t504) * t320 + 0.1e1;
t292 = qJD(2) + qJD(3);
t517 = pkin(1) * qJD(2);
t454 = t310 * t517;
t204 = pkin(4) * t292 - t454;
t422 = t532 * t517;
t398 = t304 * t422;
t153 = -t303 * t204 + t398;
t434 = qJD(4) + t292;
t149 = pkin(10) * t434 - t153;
t302 = sin(qJ(5));
t520 = t311 * pkin(1);
t267 = pkin(12) + t520;
t220 = t267 * qJD(1);
t151 = -pkin(4) * t186 - t220;
t374 = -t303 * t186 + t532 * t188;
t57 = -pkin(8) * t121 + pkin(10) * t374 + t151;
t46 = -t149 * t302 + t309 * t57;
t555 = qJD(5) - t121;
t572 = t555 * t46;
t47 = t149 * t309 + t302 * t57;
t571 = t555 * t47;
t570 = t121 * t374;
t449 = t320 * t504;
t307 = cos(qJ(7));
t442 = t307 * t481;
t300 = sin(qJ(7));
t444 = t300 * t482;
t184 = -t442 + t444;
t198 = t300 * t311 + t305 * t307;
t187 = t198 * qJD(1);
t299 = sin(pkin(15));
t518 = cos(pkin(15));
t530 = sin(qJ(8));
t531 = cos(qJ(8));
t195 = t299 * t530 + t518 * t531;
t341 = -t299 * t531 + t518 * t530;
t390 = -t184 * t341 + t187 * t195;
t391 = t184 * t195 + t187 * t341;
t568 = t390 * t391;
t399 = t309 * t434;
t103 = -t302 * t374 - t399;
t389 = t304 * t305 - t310 * t311;
t363 = t389 * qJD(3);
t328 = qJD(2) * t389 + t363;
t101 = qJD(1) * t328 - qJDD(1) * t199;
t329 = t292 * t199;
t102 = qJD(1) * t329 + qJDD(1) * t389;
t440 = qJD(4) * t532;
t477 = qJD(4) * t303;
t368 = t532 * t101 + t303 * t102 + t186 * t440 + t188 * t477;
t290 = qJDD(2) + qJDD(3);
t423 = qJDD(4) + t290;
t474 = qJD(5) * t302;
t27 = -qJD(5) * t399 - t302 * t423 - t309 * t368 - t374 * t474;
t105 = t302 * t434 - t309 * t374;
t476 = qJD(5) * t105;
t28 = t302 * t368 - t309 * t423 + t476;
t566 = -t576 * t103 - t27 * t309 - t302 * t28;
t24 = t27 * t302;
t565 = t576 * t105 - t24;
t433 = t303 * t101 - t532 * t102;
t43 = -qJD(4) * t374 + t433;
t41 = qJDD(5) + t43;
t564 = t105 * t374 + t302 * t41 + t576 * t555;
t298 = qJ(2) + qJ(3);
t288 = qJ(4) + t298;
t270 = cos(t288);
t255 = g(3) * t270;
t451 = qJD(3) * t517;
t259 = t304 * t451;
t463 = qJDD(2) * t310;
t162 = -pkin(1) * t463 + pkin(4) * t290 + t259;
t513 = pkin(1) * qJDD(2);
t558 = -t304 * t513 - t310 * t451;
t350 = -qJD(4) * t398 - t532 * t162 + t204 * t477 + t303 * t558;
t73 = -pkin(8) * t423 + t350;
t562 = -t73 + t255;
t561 = t103 * t374;
t492 = t303 * t304;
t456 = pkin(1) * t492;
t239 = qJD(2) * t456;
t152 = t532 * t204 + t239;
t148 = -pkin(8) * t434 - t152;
t509 = t121 * t148;
t428 = -pkin(4) / sin((qJ(7) + qJ(8) - t414)) * sin(t459) * t318;
t559 = t428 * t568;
t508 = t121 * t302;
t556 = t474 - t508;
t306 = sin(qJ(1));
t312 = cos(qJ(1));
t223 = g(1) * t312 + g(2) * t306;
t554 = -t121 ^ 2 + t374 ^ 2;
t269 = sin(t288);
t493 = t302 * t312;
t494 = t302 * t306;
t338 = -t374 * t47 + t148 * t473 - t562 * t302 + (g(1) * t493 + g(2) * t494) * t269;
t412 = t148 * t474 + t309 * t255 + t374 * t46;
t254 = g(3) * t269;
t497 = t270 * t312;
t498 = t270 * t306;
t447 = g(1) * t497 + g(2) * t498 + t254;
t397 = t310 * t422;
t445 = t532 * t304;
t424 = pkin(1) * t445;
t77 = -qJD(3) * t397 + qJD(4) * t239 - qJDD(2) * t424 + t303 * t162 + t204 * t440;
t336 = -t151 * t121 - t447 - t77;
t552 = -t121 * t434 + t368;
t380 = t223 * t269;
t323 = t151 * t374 + t255 - t350 - t380;
t551 = -t374 * t292 - t433;
t69 = -pkin(8) * t374 - pkin(10) * t121;
t550 = t575 * t568;
t418 = 0.1e1 + t449;
t39 = t309 * t41;
t548 = t474 * t555 - t39;
t542 = -g(1) * t306 + g(2) * t312;
t546 = t542 * t269;
t545 = t418 * t184 * t187;
t75 = -t188 * t292 + t102;
t284 = sin(t298);
t286 = cos(t298);
t340 = g(3) * t286 - t223 * t284;
t479 = qJD(2) * t311;
t543 = -qJD(7) * t311 - t479;
t400 = -t302 * t46 + t309 * t47;
t541 = t304 ^ 2 + t310 ^ 2;
t291 = qJD(2) + qJD(7);
t281 = qJD(8) + t291;
t147 = t291 * t198;
t460 = t311 * qJDD(1);
t461 = t305 * qJDD(1);
t402 = t300 * t461 - t307 * t460;
t100 = qJD(1) * t147 + t402;
t181 = t195 * qJD(8);
t182 = t341 * qJD(8);
t496 = t300 * t305;
t403 = t291 * t496;
t465 = qJD(1) * qJD(2);
t438 = t311 * t465;
t420 = -qJD(7) * t442 - t300 * t460 + (-t438 - t461) * t307;
t99 = qJD(1) * t403 + t420;
t29 = t100 * t341 + t181 * t184 + t182 * t187 + t195 * t99;
t540 = -t281 * t391 + t29;
t30 = t195 * t100 + t181 * t187 - t182 * t184 - t341 * t99;
t539 = -t281 * t390 + t30;
t106 = -t220 + (t518 * t184 - t187 * t299) * pkin(3);
t271 = -qJ(2) + t441;
t257 = sin(t271);
t258 = cos(t271);
t276 = t307 * t513;
t289 = qJDD(2) + qJDD(7);
t452 = t518 * pkin(3);
t455 = t300 * t517;
t157 = -qJD(7) * t455 + t289 * t452 + t276;
t472 = qJD(7) * t307;
t364 = qJD(2) * t472 + t300 * qJDD(2);
t528 = pkin(3) * t299;
t158 = pkin(1) * t364 + t289 * t528;
t190 = t291 * t528 + t455;
t480 = qJD(2) * t307;
t191 = pkin(1) * t480 + t291 * t452;
t48 = -t157 * t341 - t195 * t158 - t181 * t191 + t182 * t190;
t322 = g(3) * t257 - t106 * t391 - t223 * t258 - t48;
t49 = -t195 * t157 + t158 * t341 + t181 * t190 + t182 * t191;
t339 = g(3) * t258 + t106 * t390 + t223 * t257 + t49;
t537 = t390 ^ 2 - t391 ^ 2;
t529 = pkin(1) * t305;
t527 = pkin(4) * t188;
t526 = pkin(8) * t269;
t521 = g(3) * t311;
t439 = t305 * t465;
t175 = pkin(1) * t439 - qJDD(1) * t267;
t76 = -pkin(4) * t102 + t175;
t11 = pkin(8) * t43 - pkin(10) * t368 + t76;
t72 = pkin(10) * t423 + t77;
t7 = qJD(5) * t46 + t11 * t302 + t309 * t72;
t6 = t7 * t309;
t26 = t28 * t309;
t512 = t103 * t302;
t511 = t105 * t103;
t139 = -t532 * t199 + t303 * t389;
t506 = t139 * t302;
t505 = t139 * t309;
t502 = t184 * t291;
t501 = t186 * t188;
t500 = t187 * t291;
t301 = sin(qJ(6));
t308 = cos(qJ(6));
t495 = t301 * t308;
t491 = t306 * t309;
t490 = t309 * t312;
t371 = pkin(1) * (-t195 * t307 + t300 * t341);
t488 = (-t518 * t181 + t182 * t299) * pkin(3) - qJD(2) * t371;
t372 = pkin(1) * (t195 * t300 + t307 * t341);
t487 = (t181 * t299 + t518 * t182) * pkin(3) - qJD(2) * t372;
t293 = t301 ^ 2;
t295 = t308 ^ 2;
t484 = t293 - t295;
t294 = t305 ^ 2;
t296 = t311 ^ 2;
t483 = t294 - t296;
t478 = qJD(3) * t292;
t475 = qJD(5) * t555;
t470 = qJDD(1) * pkin(6);
t469 = qJDD(1) * pkin(12);
t464 = qJD(1) * qJD(6);
t462 = qJDD(2) * pkin(1) ^ 2;
t278 = t305 * t517;
t315 = qJD(1) ^ 2;
t448 = t305 * t315 * t311;
t287 = -qJ(2) + t468;
t266 = cos(t287);
t436 = pkin(3) * t266 + t520;
t207 = pkin(4) * t284 - t529;
t435 = -pkin(10) * t270 + t207;
t274 = -pkin(1) * t310 + pkin(4);
t177 = t303 * t274 - t424;
t171 = pkin(10) + t177;
t277 = pkin(1) * t482;
t61 = -t527 + t69;
t58 = t277 + t61;
t431 = qJD(5) * t171 + t58;
t430 = t302 * t555;
t429 = pkin(1) * t241 * sin(t467) / pkin(5);
t427 = 0.2e1 * pkin(6) * t464;
t426 = -0.2e1 * pkin(12) * t465;
t425 = t6 + t447;
t279 = qJDD(8) + t289;
t421 = pkin(4) * t440;
t417 = t305 * t438;
t416 = t464 * t495;
t373 = t303 * t310 + t445;
t173 = t373 * t517;
t413 = pkin(4) * t477 + t173;
t410 = -pkin(4) * t286 + t520;
t409 = pkin(8) * t270 + pkin(10) * t269;
t89 = t190 * t341 - t191 * t195;
t90 = t190 * t195 + t191 * t341;
t406 = t390 * t90 + t391 * t89;
t405 = t315 * t429;
t401 = t302 * t47 + t309 * t46;
t396 = qJDD(1) * t429;
t395 = -t171 * t41 - t509;
t272 = pkin(4) * t303 + pkin(10);
t394 = -t272 * t41 - t509;
t392 = t152 * t121 + t153 * t374;
t387 = t317 * t423;
t386 = t46 * t560 + t47 * t508 + t425;
t385 = t103 * t556 - t26;
t384 = t573 * t374;
t297 = qJ(2) + qJ(7);
t283 = sin(t297);
t285 = cos(t297);
t382 = g(3) * t283 - t220 * t184 + t223 * t285;
t381 = -qJD(5) * t57 - t254 - t72;
t349 = t532 * t389;
t59 = -qJD(4) * t349 - t199 * t477 - t303 * t329 - t532 * t328;
t378 = t139 * t473 - t302 * t59;
t377 = -t139 * t474 - t309 * t59;
t367 = -t420 + t502;
t366 = t220 * t199;
t365 = t405 * t495;
t176 = t532 * t274 + t456;
t361 = -pkin(6) * t315 + t223;
t360 = pkin(12) * t315 + t223;
t358 = t574 * t374;
t357 = t542 + 0.2e1 * t470;
t356 = -t542 + 0.2e1 * t469;
t353 = -t402 + t500;
t352 = -t73 - t380;
t114 = (t184 * t299 + t518 * t187) * pkin(3);
t348 = -g(3) * t285 + t220 * t187 + t223 * t283 + t276;
t346 = t508 * t555 - t548 - t561;
t123 = -pkin(4) * t329 + t278;
t60 = t139 * qJD(4) + t303 * t328 - t532 * t329;
t21 = t60 * pkin(8) + t59 * pkin(10) + t123;
t138 = -t199 * t303 - t349;
t160 = -pkin(4) * t389 - t267;
t67 = t138 * pkin(8) - t139 * pkin(10) + t160;
t344 = qJD(5) * t139 * t148 + t21 * t555 + t41 * t67;
t343 = t139 * t73 - t148 * t59 - t475 * t67;
t337 = (t220 * t482 - t521) * pkin(1) + t223 * t529;
t335 = -t105 * t556 + t566;
t334 = -t188 * t220 + t259 + t340;
t10 = t309 * t11;
t8 = -qJD(5) * t47 - t302 * t72 + t10;
t333 = -qJD(5) * t401 - t8 * t302 + t6;
t331 = -g(3) * t284 + t186 * t220 - t223 * t286 - t558;
t314 = qJD(2) ^ 2;
t313 = qJD(6) ^ 2;
t273 = -t532 * pkin(4) - pkin(8);
t265 = sin(t287);
t222 = t307 * pkin(1) + t452;
t221 = pkin(1) * t300 + t528;
t211 = t312 * t526;
t210 = t306 * t526;
t203 = pkin(12) + t410;
t197 = -t307 * t311 + t496;
t174 = -t397 + t239;
t170 = -pkin(8) - t176;
t169 = t270 * t490 - t494;
t168 = t270 * t493 + t491;
t167 = t270 * t491 + t493;
t166 = t270 * t494 - t490;
t154 = t277 - t527;
t146 = t307 * t543 + t403;
t129 = (-t195 * t299 - t341 * t518) * pkin(3);
t128 = (-t518 * t195 + t299 * t341) * pkin(3);
t127 = -t274 * t477 + (qJD(3) * t373 + t304 * t440) * pkin(1);
t126 = t274 * t440 + (t304 * t477 + (-t532 * t310 + t492) * qJD(3)) * pkin(1);
t115 = (t518 * t197 - t198 * t299) * pkin(3) - t267;
t110 = t277 + t114;
t109 = -t195 * t221 - t222 * t341;
t108 = -t195 * t222 + t221 * t341;
t107 = -t186 ^ 2 + t188 ^ 2;
t98 = -t195 * t198 + t197 * t341;
t97 = t195 * t197 + t198 * t341;
t83 = -t220 * t278 + (-t175 - t542) * t267;
t74 = -t186 * t292 + t101;
t71 = t278 + (t146 * t299 + t518 * t147) * pkin(3);
t66 = qJD(7) * t372 + t181 * t221 + t182 * t222;
t65 = qJD(7) * t371 - t181 * t222 + t182 * t221;
t56 = t152 * t309 + t302 * t69;
t55 = -t152 * t302 + t309 * t69;
t52 = t174 * t309 + t302 * t61;
t51 = -t174 * t302 + t309 * t61;
t50 = (t518 * t100 + t299 * t99) * pkin(3) + t175;
t45 = -t146 * t341 + t147 * t195 + t181 * t198 - t182 * t197;
t44 = t146 * t195 + t147 * t341 + t181 * t197 + t182 * t198;
t34 = -t153 * t434 + t323;
t33 = t152 * t434 + t336;
t17 = -t281 * t90 + t339;
t16 = t281 * t89 + t322;
t15 = t512 * t555 - t26;
t12 = -t430 * t555 + t39 - t561;
t5 = -t105 * t430 + t566;
t4 = -pkin(8) * t28 + t103 * t153 - t555 * t55 + (-pkin(10) * t41 - t509) * t302 + (-pkin(10) * t475 + t352) * t309 + t412;
t3 = pkin(8) * t27 + pkin(10) * t548 + t105 * t153 - t148 * t560 + t555 * t56 + t338;
t2 = -g(1) * t211 - g(2) * t210 + t148 * t153 - t46 * t55 - t47 * t56 + t562 * pkin(8) + (t223 * t270 + t254 + t333) * pkin(10);
t1 = t103 * t56 + t105 * t55 + (-t572 + (-t28 + t476) * pkin(10)) * t309 + (-t8 - t571 + (qJD(5) * t103 - t27) * pkin(10)) * t302 + t425;
t9 = [0, 0, 0, 0, 0, qJDD(1), -t542, t223, 0, 0, qJDD(1) * t294 + 0.2e1 * t417, 0.2e1 * t305 * t460 - 0.2e1 * t465 * t483, qJDD(2) * t305 + t311 * t314, qJDD(1) * t296 - 0.2e1 * t417, qJDD(2) * t311 - t305 * t314, 0, t305 * t426 + t311 * t356, -t305 * t356 + t311 * t426, -t223, (-t542 + t469) * pkin(12), -t101 * t199 - t188 * t328, t101 * t389 - t199 * t102 + t292 * (t186 * t389 - t188 * t199), -t199 * t290 + t292 * t328, t102 * t389 + t186 * t329, t290 * t389 + t292 * t329, 0, t267 * t102 - t175 * t389 + t542 * t286 + qJD(3) * t366 + (-t186 * t529 + t366) * qJD(2), -t267 * t101 - t175 * t199 - t542 * t284 - t220 * t363 + (-t188 * t529 - t220 * t389) * qJD(2), ((-t310 * t199 - t304 * t389) * qJDD(2) - t541 * qJD(2) * t479) * pkin(1) - t223, t83, t139 * t368 + t374 * t59, -t121 * t59 - t138 * t368 - t139 * t43 + t374 * t60, t139 * t423 - t434 * t59, -t121 * t60 + t138 * t43, -t138 * t423 - t434 * t60, 0, -t121 * t123 + t138 * t76 + t151 * t60 + t160 * t43 + t270 * t542, -t123 * t374 + t139 * t76 - t151 * t59 + t160 * t368 - t546, -t138 * t77 + t139 * t350 + t152 * t59 + t153 * t60 - t223, t123 * t151 + t160 * t76 - t203 * t542, t105 * t377 - t27 * t505, (t103 * t309 + t105 * t302) * t59 + (t24 - t26 + (-t105 * t309 + t512) * qJD(5)) * t139, t105 * t60 - t138 * t27 + t377 * t555 + t41 * t505, t103 * t378 + t28 * t506, -t103 * t60 - t138 * t28 - t378 * t555 - t41 * t506, t138 * t41 + t555 * t60, -g(1) * t167 + g(2) * t169 + t138 * t8 + t302 * t343 + t309 * t344 + t46 * t60, g(1) * t166 - g(2) * t168 - t138 * t7 - t302 * t344 + t309 * t343 - t47 * t60, t546 + (-t105 * t21 - t139 * t8 + t27 * t67 + t46 * t59 + (-t103 * t67 - t139 * t47) * qJD(5)) * t309 + (-t103 * t21 - t139 * t7 - t28 * t67 + t47 * t59 + (t105 * t67 + t139 * t46) * qJD(5)) * t302, t401 * t21 + (qJD(5) * t400 + t7 * t302 + t8 * t309) * t67 + t542 * (-t203 + t409), qJDD(1) * t293 + 0.2e1 * t416, 0.2e1 * qJDD(1) * t495 - 0.2e1 * t464 * t484, qJDD(6) * t301 + t308 * t313, qJDD(1) * t295 - 0.2e1 * t416, qJDD(6) * t308 - t301 * t313, 0, t301 * t427 - t308 * t357, t301 * t357 + t308 * t427, -t223, (t542 + t470) * pkin(6), -t146 * t187 - t198 * t99, -t100 * t198 + t146 * t184 - t147 * t187 + t197 * t99, -t146 * t291 + t198 * t289, t100 * t197 + t147 * t184, -t147 * t291 - t197 * t289, 0, -t100 * t267 - t147 * t220 + t175 * t197 + t184 * t278 - t285 * t542, t146 * t220 + t175 * t198 + t187 * t278 + t267 * t99 + t283 * t542, ((-t197 * t300 - t198 * t307) * qJDD(2) + (t146 * t307 - t147 * t300 + (-t197 * t307 + t198 * t300) * qJD(7)) * qJD(2)) * pkin(1) - t223, t83, t29 * t98 - t390 * t44, t29 * t97 + t30 * t98 - t390 * t45 + t391 * t44, t279 * t98 + t281 * t44, t30 * t97 + t391 * t45, t279 * t97 + t281 * t45, 0, -t106 * t45 - t115 * t30 + t258 * t542 - t391 * t71 - t50 * t97, t106 * t44 + t115 * t29 + t257 * t542 - t390 * t71 + t50 * t98, -t44 * t89 - t45 * t90 + t48 * t97 - t49 * t98 - t223, t106 * t71 + t115 * t50 - t542 * (pkin(12) + t436); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, t483 * t315, t461, t448, t460, qJDD(2), t305 * t360 - t521, g(3) * t305 + t311 * t360, 0, 0, t501, t107, t74, -t501, t75, t290, (t186 * t482 + t304 * t478 + (-qJDD(2) - t290) * t310) * pkin(1) + t334, (t188 * t482 + t290 * t304 + t310 * t478) * pkin(1) + t331, -t186 * t454 + ((-qJD(3) * t186 + t101) * t310 - t75 * t304) * pkin(1), t462 * t541 + t337, t574 * t570, t574 * t554, t574 * t552, -t358 * t121, t574 * t551, t387 * t538 + t423, t121 * t154 + t127 * t434 + t176 * t423 + t34 * t404 + t323, -t126 * t434 + t154 * t374 - t177 * t423 + t33 * t404 + t336, t121 * t126 + t127 * t374 - t176 * t368 - t177 * t43 + t392, -g(3) * t410 - t153 * t126 + t152 * t127 - t151 * t154 - t176 * t350 + t77 * t177 - t207 * t223, t574 * t565, t404 * t5 + t335, t574 * t564, t15 * t404 + t385, t12 * t404 + t346, t358 * t555, t4 * t404 - t103 * t127 + t170 * t28 + (-t126 * t555 + t395) * t302 + (-t431 * t555 + t352) * t309 + t412, t3 * t404 - t105 * t127 - t170 * t27 + t395 * t309 + (-t126 * t309 + t302 * t431) * t555 + t338, t1 * t404 + (-t103 * t126 + t105 * t58 - t171 * t28 + (t105 * t171 - t46) * qJD(5)) * t309 + (t103 * t58 + t105 * t126 - t171 * t27 - t8 + (t103 * t171 - t47) * qJD(5)) * t302 + t386, t73 * t170 - t148 * t127 - g(1) * (t312 * t435 + t211) - g(2) * (t306 * t435 + t210) - g(3) * (-t409 + t410) + t2 * t404 - t401 * t58 + t400 * t126 + t333 * t171, -t365, t484 * t405, t301 * t396, t365, t308 * t396, qJDD(6) * t429, (-g(3) * t308 + t301 * t361) * t429, (g(3) * t301 + t308 * t361) * t429, 0, 0, t545, t418 * (-t184 ^ 2 + t187 ^ 2), -qJD(7) * t444 - t300 * t439 + (-t291 * t444 + t367) * t449 + t367, -t545, t353 * t449 + t353 + t418 * qJD(1) * ((-t472 - t480) * t305 + t543 * t300), t418 * t289, t348 * t449 + (-t184 * t482 + t307 * t289 + (-qJD(7) * t291 + (-qJD(7) + (-qJD(7) + t291) * t449) * qJD(2)) * t300) * pkin(1) + t348, t382 * t449 + (-t291 * t472 - t300 * t289 - t187 * t482 + (t291 * t480 - t364) * t449 - t364) * pkin(1) + t382, ((t99 - t502) * t307 + (-t100 + t500) * t300) * pkin(1), (t300 ^ 2 + t307 ^ 2) * t462 + t337, t550, t575 * t537, t575 * t540, -t550, t575 * t539, t575 * t279, t108 * t279 + t110 * t391 + t66 * t281 + ((t114 * t391 + t128 * t279 + t281 * t487 + t339) * t504 + t17 * t453) * t320 + t339, -t109 * t279 + t110 * t390 - t65 * t281 + ((t114 * t390 - t129 * t279 - t281 * t488 + t322) * t504 + t16 * t453) * t320 + t322, -t108 * t29 + t109 * t30 + t65 * t391 + t66 * t390 + (-t128 * t29 + t129 * t30 + t390 * t487 + t391 * t488 + t406) * t449 + t406, t48 * t109 - t90 * t65 + t49 * t108 + t89 * t66 - t106 * t110 - g(3) * t436 + (-t106 * t114 + t128 * t49 + t129 * t48 - t488 * t90 + t487 * t89 + (-g(3) * t266 - t223 * t265) * pkin(3)) * t449 - t223 * (pkin(3) * t265 - t529); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, t107, t74, -t501, t75, t290, (-qJD(2) * t292 * t304 - t463) * pkin(1) + t334, -t292 * t454 + t331, 0, 0, -t573 * t570, -t573 * t554, -t573 * t552, t384 * t121, -t573 * t551, -t387 * t547 + t423, -t173 * t434 - t34 * t450 + (-t121 * t188 + t423 * t532 - t434 * t477) * pkin(4) + t323, t174 * t434 - t33 * t450 + (-t188 * t374 - t303 * t423 - t434 * t440) * pkin(4) + t336, -t174 * t121 - t173 * t374 + (-t532 * t368 - t303 * t43 + (t121 * t532 - t303 * t374) * qJD(4)) * pkin(4) + t392, -t152 * t173 + t153 * t174 + (-t532 * t350 + t151 * t188 + t303 * t77 + (-t152 * t303 - t153 * t532) * qJD(4) + t340) * pkin(4), -t573 * t565, -t450 * t5 + t335, -t573 * t564, -t15 * t450 + t385, -t12 * t450 + t346, -t384 * t555, -t4 * t450 - t51 * t555 + t273 * t28 + t413 * t103 + (-t421 * t555 + t394) * t302 + (-t272 * t475 + t352) * t309 + t412, -t3 * t450 - t273 * t27 + t394 * t309 + t413 * t105 + (t272 * t474 - t309 * t421 + t52) * t555 + t338, -t1 * t450 + t52 * t103 + t51 * t105 + (-t103 * t421 - t272 * t28 + (t105 * t272 - t46) * qJD(5)) * t309 + (t105 * t421 - t27 * t272 - t8 + (t103 * t272 - t47) * qJD(5)) * t302 + t386, t73 * t273 - t47 * t52 - t46 * t51 + t148 * t173 - g(1) * (-pkin(10) * t497 + t211) - g(2) * (-pkin(10) * t498 + t210) + g(3) * t409 - t2 * t450 + ((t148 * t303 + t400 * t532) * qJD(4) + t340) * pkin(4) + t333 * t272, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t559, t537 * t428, t540 * t428, -t559, t539 * t428, t279 * t428, t17 * t428, t16 * t428, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t511, -t103 ^ 2 + t105 ^ 2, t103 * t555 - t27, -t511, t105 * t555 - t28, t41, -g(1) * t168 - g(2) * t166 - t105 * t148 - t149 * t473 + t302 * t381 + t10 + t571, -g(1) * t169 - g(2) * t167 + t103 * t148 + t572 + (qJD(5) * t149 - t11) * t302 + t381 * t309, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t9;
