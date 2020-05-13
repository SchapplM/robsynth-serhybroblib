% Calculate vector of inverse dynamics joint torques with ic for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnIC_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnIC_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:13
% EndTime: 2020-05-07 11:33:30
% DurationCPUTime: 15.08s
% Computational Cost: add. (6403->665), mult. (11311->956), div. (4->3), fcn. (9044->11), ass. (0->385)
t291 = cos(qJ(2));
t289 = sin(qJ(1));
t292 = cos(qJ(1));
t418 = t289 ^ 2 + t292 ^ 2;
t533 = t291 * t418;
t286 = qJ(2) + qJ(3);
t269 = sin(t286);
t270 = cos(t286);
t345 = Icges(4,5) * t269 + Icges(4,6) * t270;
t170 = t345 * t289;
t171 = t345 * t292;
t283 = qJD(2) + qJD(3);
t219 = t283 * t289;
t365 = rSges(4,1) * t269 + rSges(4,2) * t270;
t532 = t365 * t219;
t456 = t269 * t289;
t383 = -rSges(4,2) * t456 + rSges(4,3) * t292;
t454 = t270 * t289;
t140 = rSges(4,1) * t454 + t383;
t220 = t283 * t292;
t280 = t289 * rSges(4,3);
t455 = t269 * t292;
t422 = rSges(4,2) * t455 + t280;
t453 = t270 * t292;
t327 = rSges(4,1) * t453 - t422;
t493 = pkin(2) * qJD(2);
t55 = -t219 * t140 - t220 * t327 + t493 * t533;
t531 = t55 * qJD(2);
t472 = Icges(4,6) * t292;
t134 = Icges(4,4) * t454 - Icges(4,2) * t456 + t472;
t481 = Icges(4,4) * t270;
t352 = Icges(4,1) * t269 + t481;
t530 = -t352 * t289 - t134;
t349 = -Icges(4,2) * t269 + t481;
t135 = Icges(4,6) * t289 - t292 * t349;
t529 = -t352 * t292 + t135;
t177 = rSges(4,1) * t455 + rSges(4,2) * t453;
t482 = Icges(4,4) * t269;
t348 = Icges(4,2) * t270 + t482;
t353 = Icges(4,1) * t270 - t482;
t528 = -t348 + t353;
t527 = -t349 - t352;
t404 = qJD(1) * qJD(2);
t264 = t292 * t404;
t216 = qJDD(2) * t289 + t264;
t403 = qJD(1) * qJD(3);
t145 = qJDD(3) * t289 + t292 * t403 + t216;
t263 = t289 * t404;
t146 = t289 * t403 + t263 + (-qJDD(2) - qJDD(3)) * t292;
t217 = -qJDD(2) * t292 + t263;
t293 = qJD(2) ^ 2;
t288 = sin(qJ(2));
t378 = t418 * t288;
t497 = rSges(4,1) * t270;
t397 = qJD(1) * t497;
t409 = qJD(1) * t292;
t368 = rSges(4,3) * t409 + t177 * t283 + t289 * t397;
t262 = t269 * rSges(4,2);
t396 = qJD(1) * t262;
t79 = -t289 * t396 + t368;
t525 = t262 - t497;
t80 = t532 + (t292 * t525 + t280) * qJD(1);
t13 = -t145 * t140 + t219 * t80 + t146 * t327 + t220 * t79 + ((t216 * t289 - t217 * t292) * t291 - t293 * t378) * pkin(2);
t153 = t525 * t283;
t176 = rSges(4,1) * t456 + rSges(4,2) * t454;
t317 = qJD(1) * t327;
t401 = qJDD(1) * t289;
t294 = qJD(1) ^ 2;
t417 = -t293 - t294;
t33 = -qJD(1) * t80 + qJDD(1) * t140 - t146 * t365 - t153 * t220 + ((t217 + t263) * t288 + (t292 * t417 - t401) * t291) * pkin(2);
t400 = qJDD(1) * t292;
t32 = t145 * t365 - t219 * t153 - qJDD(1) * t327 + qJD(1) * t79 + ((-t216 - t264) * t288 + (t289 * t417 + t400) * t291) * pkin(2);
t487 = t32 * t289;
t500 = pkin(2) * t291;
t251 = t409 * t500;
t398 = t288 * t493;
t64 = t289 * t398 - t251 + t317 - t532;
t407 = qJD(2) * t292;
t395 = t288 * t407;
t410 = qJD(1) * t289;
t65 = -pkin(2) * (t291 * t410 + t395) + qJD(1) * t140 + t220 * t365;
t526 = t13 * (-t289 * t140 - t292 * t327) - (-t487 - t33 * t292 + (t289 * t65 + t292 * t64) * qJD(1)) * t365 + (t289 * t64 - t292 * t65) * t153 + t64 * (qJD(1) * t177 - t219 * t525) - g(1) * t177 - (-qJD(1) * t176 - t220 * t525) * t65 + (-t140 * t409 + t292 * t79 + (t317 + t80) * t289 - t219 * t176 - t177 * t220) * t55;
t290 = cos(qJ(4));
t277 = Icges(5,4) * t290;
t287 = sin(qJ(4));
t347 = -Icges(5,2) * t287 + t277;
t229 = Icges(5,1) * t287 + t277;
t278 = Icges(3,4) * t291;
t350 = -Icges(3,2) * t288 + t278;
t231 = Icges(3,1) * t288 + t278;
t224 = Icges(3,5) * t291 - Icges(3,6) * t288;
t223 = Icges(3,5) * t288 + Icges(3,6) * t291;
t321 = qJD(2) * t223;
t483 = Icges(3,4) * t288;
t232 = Icges(3,1) * t291 - t483;
t165 = Icges(3,5) * t289 + t232 * t292;
t161 = Icges(3,6) * t289 + t292 * t350;
t461 = t161 * t288;
t336 = -t165 * t291 + t461;
t470 = Icges(3,3) * t292;
t524 = -t292 * t321 + (-t224 * t289 + t336 + t470) * qJD(1);
t450 = t288 * t289;
t257 = Icges(3,4) * t450;
t447 = t289 * t291;
t479 = Icges(3,5) * t292;
t164 = Icges(3,1) * t447 - t257 - t479;
t473 = Icges(3,6) * t292;
t160 = Icges(3,4) * t447 - Icges(3,2) * t450 - t473;
t462 = t160 * t288;
t337 = -t164 * t291 + t462;
t157 = Icges(3,3) * t289 + t224 * t292;
t414 = qJD(1) * t157;
t523 = qJD(1) * t337 - t289 * t321 + t414;
t222 = Icges(5,5) * t290 - Icges(5,6) * t287;
t221 = Icges(5,5) * t287 + Icges(5,6) * t290;
t318 = qJD(4) * t221;
t480 = Icges(5,4) * t287;
t230 = Icges(5,1) * t290 - t480;
t163 = Icges(5,5) * t289 + t230 * t292;
t159 = Icges(5,6) * t289 + t292 * t347;
t463 = t159 * t287;
t338 = -t163 * t290 + t463;
t468 = Icges(5,3) * t292;
t522 = -t292 * t318 + (-t222 * t289 + t338 + t468) * qJD(1);
t452 = t287 * t289;
t256 = Icges(5,4) * t452;
t448 = t289 * t290;
t476 = Icges(5,5) * t292;
t162 = Icges(5,1) * t448 - t256 - t476;
t471 = Icges(5,6) * t292;
t158 = Icges(5,4) * t448 - Icges(5,2) * t452 - t471;
t464 = t158 * t287;
t339 = -t162 * t290 + t464;
t155 = Icges(5,3) * t289 + t222 * t292;
t415 = qJD(1) * t155;
t521 = qJD(1) * t339 - t289 * t318 + t415;
t346 = Icges(4,5) * t270 - Icges(4,6) * t269;
t242 = Icges(4,4) * t455;
t478 = Icges(4,5) * t289;
t137 = -Icges(4,1) * t453 + t242 + t478;
t467 = t137 * t270;
t469 = Icges(4,3) * t292;
t520 = -t283 * t171 + (t135 * t269 - t289 * t346 - t467 - t469) * qJD(1);
t241 = Icges(4,4) * t456;
t477 = Icges(4,5) * t292;
t136 = Icges(4,1) * t454 - t241 + t477;
t341 = t134 * t269 - t136 * t270;
t133 = Icges(4,3) * t289 - t292 * t346;
t416 = qJD(1) * t133;
t519 = qJD(1) * t341 + t170 * t283 + t416;
t156 = Icges(3,5) * t447 - Icges(3,6) * t450 - t470;
t58 = -t156 * t292 - t289 * t337;
t154 = Icges(5,5) * t448 - Icges(5,6) * t452 - t468;
t56 = -t154 * t292 - t289 * t339;
t333 = -t269 * t348 + t270 * t352;
t518 = qJD(1) * t333 + t346 * t283;
t225 = Icges(5,2) * t290 + t480;
t331 = t225 * t287 - t229 * t290;
t517 = t331 * qJD(1) + t222 * qJD(4);
t227 = Icges(3,2) * t291 + t483;
t329 = t227 * t288 - t231 * t291;
t516 = t329 * qJD(1) + t224 * qJD(2);
t515 = t289 * (-t227 * t292 + t165) - t292 * (-Icges(3,2) * t447 + t164 - t257);
t514 = t289 * (-t225 * t292 + t163) - t292 * (-Icges(5,2) * t448 + t162 - t256);
t513 = qJD(1) * t527 + t219 * (Icges(4,2) * t453 + t137 + t242) + t220 * (-Icges(4,2) * t454 + t136 - t241);
t512 = t145 / 0.2e1;
t511 = t146 / 0.2e1;
t402 = qJD(1) * qJD(4);
t214 = qJDD(4) * t289 + t292 * t402;
t510 = t214 / 0.2e1;
t215 = -qJDD(4) * t292 + t289 * t402;
t509 = t215 / 0.2e1;
t508 = t216 / 0.2e1;
t507 = t217 / 0.2e1;
t506 = -t219 / 0.2e1;
t505 = t219 / 0.2e1;
t504 = -t220 / 0.2e1;
t503 = t220 / 0.2e1;
t502 = t289 / 0.2e1;
t501 = -t292 / 0.2e1;
t499 = -qJD(1) / 0.2e1;
t498 = qJD(1) / 0.2e1;
t496 = rSges(5,1) * t290;
t495 = rSges(5,2) * t287;
t494 = rSges(5,2) * t290;
t233 = rSges(5,1) * t287 + t494;
t192 = t233 * t292;
t279 = t289 * rSges(5,3);
t446 = t290 * t292;
t419 = rSges(5,1) * t446 + t279;
t451 = t287 * t292;
t168 = -rSges(5,2) * t451 + t419;
t406 = qJD(4) * t289;
t94 = -t233 * t406 + (pkin(1) * t292 + t168) * qJD(1);
t492 = t192 * t94;
t258 = rSges(5,2) * t452;
t420 = t292 * rSges(5,3) + t258;
t166 = rSges(5,1) * t448 - t420;
t405 = qJD(4) * t292;
t93 = -t233 * t405 + (-pkin(1) * t289 - t166) * qJD(1);
t490 = t289 * t93;
t489 = t292 * rSges(3,3);
t488 = t292 * t93;
t486 = qJDD(1) / 0.2e1;
t460 = t221 * t289;
t459 = t221 * t292;
t458 = t223 * t289;
t457 = t223 * t292;
t449 = t288 * t292;
t445 = t291 * t292;
t132 = Icges(4,5) * t454 - Icges(4,6) * t456 + t469;
t444 = t292 * t132;
t67 = t348 * t456 - t352 * t454 - t171;
t443 = t67 * qJD(1);
t89 = -t289 * t331 - t459;
t442 = t89 * qJD(1);
t90 = -t289 * t329 - t457;
t441 = t90 * qJD(1);
t440 = -qJ(4) + qJ(2);
t439 = t289 * t132 + t134 * t455;
t438 = t289 * t133 + t135 * t455;
t437 = -t289 * t154 - t162 * t446;
t436 = t289 * t155 + t163 * t446;
t435 = -t289 * t156 - t164 * t445;
t434 = t289 * t157 + t165 * t445;
t426 = -t225 + t230;
t425 = t229 + t347;
t424 = -t227 + t232;
t423 = t231 + t350;
t421 = rSges(5,3) * t409 + qJD(1) * t258;
t237 = rSges(3,1) * t291 - rSges(3,2) * t288;
t169 = rSges(3,3) * t289 + t237 * t292;
t413 = qJD(1) * t169;
t412 = qJD(1) * t222;
t411 = qJD(1) * t224;
t408 = qJD(2) * t289;
t394 = -pkin(1) - t496;
t393 = t410 / 0.2e1;
t392 = t409 / 0.2e1;
t391 = -t408 / 0.2e1;
t390 = t408 / 0.2e1;
t389 = -t407 / 0.2e1;
t388 = t407 / 0.2e1;
t387 = -t406 / 0.2e1;
t386 = t406 / 0.2e1;
t385 = -t405 / 0.2e1;
t384 = t405 / 0.2e1;
t382 = -(-t292 * t353 + t478) * qJD(1) + t530 * t283;
t381 = -(t289 * t353 + t477) * qJD(1) + t529 * t283;
t325 = t283 * t348;
t380 = -qJD(1) * t135 + t136 * t283 - t289 * t325;
t379 = t137 * t283 + t292 * t325 + (t289 * t349 + t472) * qJD(1);
t113 = t135 * t456;
t377 = t133 * t292 - t113;
t124 = t163 * t448;
t376 = t155 * t292 - t124;
t125 = t165 * t447;
t375 = t157 * t292 - t125;
t374 = t132 + t467;
t373 = t527 * t283;
t372 = t528 * t283;
t371 = -t154 + t463;
t370 = -t156 + t461;
t367 = t497 - t500;
t238 = rSges(2,1) * t292 - rSges(2,2) * t289;
t235 = rSges(2,1) * t289 + rSges(2,2) * t292;
t234 = rSges(3,1) * t288 + rSges(3,2) * t291;
t236 = -t495 + t496;
t319 = qJD(4) * t225;
t100 = qJD(1) * t159 - t289 * t319;
t320 = qJD(4) * t229;
t104 = qJD(1) * t163 - t289 * t320;
t85 = t158 * t290 + t162 * t287;
t303 = qJD(1) * t154 - qJD(4) * t85 - t100 * t287 + t104 * t290;
t103 = -t292 * t320 + (-t230 * t289 + t476) * qJD(1);
t86 = t159 * t290 + t163 * t287;
t99 = -t292 * t319 + (-t289 * t347 + t471) * qJD(1);
t304 = -qJD(4) * t86 + t103 * t290 - t287 * t99 + t415;
t364 = -(t289 * t521 + t303 * t292) * t292 + (t289 * t522 + t304 * t292) * t289;
t322 = qJD(2) * t227;
t101 = -t292 * t322 + (-t289 * t350 + t473) * qJD(1);
t323 = qJD(2) * t231;
t105 = -t292 * t323 + (-t232 * t289 + t479) * qJD(1);
t88 = t161 * t291 + t165 * t288;
t301 = -qJD(2) * t88 - t101 * t288 + t105 * t291 + t414;
t102 = qJD(1) * t161 - t289 * t322;
t106 = qJD(1) * t165 - t289 * t323;
t87 = t160 * t291 + t164 * t288;
t302 = qJD(1) * t156 - qJD(2) * t87 - t102 * t288 + t106 * t291;
t363 = -(t289 * t523 + t302 * t292) * t292 + (t289 * t524 + t301 * t292) * t289;
t362 = -(t303 * t289 - t292 * t521) * t292 + (t304 * t289 - t292 * t522) * t289;
t361 = -(t302 * t289 - t292 * t523) * t292 + (t301 * t289 - t292 * t524) * t289;
t57 = -t159 * t452 - t376;
t360 = t289 * t57 - t292 * t56;
t59 = -t161 * t450 - t375;
t359 = t289 * t59 - t292 * t58;
t60 = -t158 * t451 - t437;
t61 = -t159 * t451 + t436;
t358 = t289 * t61 - t292 * t60;
t62 = -t160 * t449 - t435;
t63 = -t161 * t449 + t434;
t357 = t289 * t63 - t292 * t62;
t355 = -t289 * t94 - t488;
t107 = -t405 * t494 + (-t287 * t405 - t290 * t410) * rSges(5,1) + t421;
t190 = t233 * t289;
t109 = -qJD(4) * t190 + (t236 * t292 + t279) * qJD(1);
t344 = t107 * t292 + t109 * t289;
t326 = qJD(2) * t234;
t108 = -t292 * t326 + (-t237 * t289 + t489) * qJD(1);
t110 = -t289 * t326 + t413;
t343 = t108 * t292 + t110 * t289;
t111 = t234 * t408 - t413;
t167 = rSges(3,1) * t447 - rSges(3,2) * t450 - t489;
t112 = -qJD(1) * t167 - t234 * t407;
t342 = t111 * t289 - t112 * t292;
t69 = t134 * t270 + t136 * t269;
t335 = t166 * t289 + t168 * t292;
t334 = t167 * t289 + t169 * t292;
t332 = t225 * t290 + t229 * t287;
t330 = t227 * t291 + t231 * t288;
t324 = t341 * t289;
t315 = -qJD(1) * t346 - t170 * t220 + t171 * t219;
t314 = t158 * t292 - t159 * t289;
t313 = t160 * t292 - t161 * t289;
t312 = (-t287 * t425 + t290 * t426) * qJD(1);
t311 = (-t288 * t423 + t291 * t424) * qJD(1);
t307 = -qJD(1) * t132 - t269 * t380 + t270 * t382;
t10 = t307 * t289 - t292 * t519;
t306 = t269 * t379 + t270 * t381 + t416;
t11 = t306 * t289 + t292 * t520;
t51 = -t324 + t444;
t52 = -t137 * t454 - t377;
t24 = t219 * t52 - t220 * t51 - t443;
t53 = t136 * t453 - t439;
t54 = -t137 * t453 + t438;
t68 = t292 * t333 - t170;
t66 = t68 * qJD(1);
t25 = t219 * t54 - t220 * t53 + t66;
t308 = t528 * qJD(1) + t529 * t219 - t530 * t220;
t296 = t269 * t513 + t308 * t270;
t305 = -qJD(1) * t345 + t269 * t373 + t270 * t372;
t34 = -t289 * t518 + t305 * t292;
t35 = t305 * t289 + t292 * t518;
t36 = t269 * t382 + t270 * t380;
t37 = t269 * t381 - t270 * t379;
t70 = -t135 * t270 - t137 * t269;
t8 = t289 * t519 + t307 * t292;
t9 = -t289 * t520 + t306 * t292;
t310 = (qJD(1) * t34 + qJDD(1) * t68 + t145 * t54 + t146 * t53 + t219 * t9 - t220 * t8) * t502 + (t308 * t269 - t270 * t513) * t499 + t24 * t393 + t25 * t392 + (qJD(1) * t35 - qJDD(1) * t67 - t10 * t220 + t11 * t219 + t145 * t52 + t146 * t51) * t501 + (t289 * t54 - t292 * t53) * t512 + (t289 * t52 - t292 * t51) * t511 + (t289 * t9 - t292 * t8 + (t289 * t53 + t292 * t54) * qJD(1)) * t505 + (-t10 * t292 + t11 * t289 + (t289 * t51 + t292 * t52) * qJD(1)) * t504 + (t289 * t70 - t292 * t69) * t486 + (t289 * t37 - t292 * t36 + (t289 * t69 + t292 * t70) * qJD(1)) * t498 + (t289 * t315 + t292 * t296) * t506 + (t289 * t296 - t292 * t315) * t503;
t204 = t347 * qJD(4);
t206 = t230 * qJD(4);
t300 = qJD(1) * t221 - qJD(4) * t332 - t204 * t287 + t206 * t290;
t205 = t350 * qJD(2);
t207 = t232 * qJD(2);
t299 = qJD(1) * t223 - qJD(2) * t330 - t205 * t288 + t207 * t291;
t298 = -t287 * t514 + t314 * t290;
t297 = -t288 * t515 + t313 * t291;
t265 = sin(qJ(3) + t440);
t211 = t237 * qJD(2);
t210 = t236 * qJD(4);
t193 = t234 * t292;
t191 = t234 * t289;
t92 = -t292 * t329 + t458;
t91 = -t292 * t331 + t460;
t84 = t92 * qJD(1);
t83 = t91 * qJD(1);
t82 = t334 * qJD(2);
t81 = t335 * qJD(4);
t50 = -qJD(1) * t110 - qJDD(1) * t167 - t211 * t407 + t217 * t234;
t49 = qJD(1) * t108 + qJDD(1) * t169 - t211 * t408 - t216 * t234;
t48 = -t210 * t406 + qJD(1) * t107 + qJDD(1) * t168 - t214 * t233 + (-t289 * t294 + t400) * pkin(1);
t47 = -t210 * t405 - qJD(1) * t109 - qJDD(1) * t166 + t215 * t233 + (-t292 * t294 - t401) * pkin(1);
t45 = t299 * t289 - t292 * t516;
t44 = t300 * t289 - t292 * t517;
t43 = t289 * t516 + t299 * t292;
t42 = t289 * t517 + t300 * t292;
t41 = -qJD(2) * t336 + t101 * t291 + t105 * t288;
t40 = -t337 * qJD(2) + t102 * t291 + t106 * t288;
t39 = -qJD(4) * t338 + t103 * t287 + t290 * t99;
t38 = -t339 * qJD(4) + t100 * t290 + t104 * t287;
t31 = qJD(2) * t357 + t84;
t30 = qJD(4) * t358 + t83;
t29 = qJD(2) * t359 + t441;
t28 = qJD(4) * t360 + t442;
t1 = [(-t108 * t111 - t110 * t112 + (-g(2) + t49) * t169 + (g(1) - t50) * t167) * m(3) + (t70 + t68) * t512 + (t86 + t91) * t510 + (t85 + t89) * t509 + (t88 + t92) * t508 + (t87 + t90) * t507 + (t37 + t34) * t505 + (t443 + (t54 - t324 - t438) * t220 + (t289 * t374 - t113 + t53) * t219 + ((t133 + t341) * t219 + t374 * t220) * t292 + t24) * t506 + (t94 * t421 + (t233 * t490 - t492) * qJD(4) + (t48 - g(2)) * ((pkin(1) - t495) * t292 + t419) + (t47 - g(1)) * (t289 * t394 + t420)) * m(5) + (m(2) * (t235 ^ 2 + t238 ^ 2) + t330 + t332 + Icges(2,3) + t348 * t270 + t352 * t269) * qJDD(1) + (t66 + (t52 + (-t136 * t292 + t137 * t289) * t270 + t377 + t439) * t220 + (t134 * t456 - t444 + t51 + (-t136 * t289 - t137 * t292) * t270 + t438) * t219) * t503 + t69 * t511 + (t65 * (-t251 + (-t396 + t397) * t292) - t64 * (-pkin(2) * t395 + t368) + (t65 * (-t283 * t365 + t398) + (-t65 * rSges(4,3) - t64 * (-t262 - t500)) * qJD(1)) * t289 + (t32 - g(2)) * (-t292 * t367 + t422) + (t33 - g(1)) * (t289 * t367 + t383)) * m(4) + (-t442 + ((t292 * t371 - t436 + t61) * t292 + (t289 * t371 + t376 + t60) * t289) * qJD(4) + t28) * t387 + (-t441 + ((t292 * t370 - t434 + t63) * t292 + (t289 * t370 + t375 + t62) * t289) * qJD(2) + t29) * t391 + (t83 + ((t57 - t124 + (t155 + t464) * t292 + t437) * t292 + t436 * t289) * qJD(4)) * t384 + (t84 + ((t59 - t125 + (t157 + t462) * t292 + t435) * t292 + t434 * t289) * qJD(2)) * t388 + (t36 + t35 + t25) * t504 + (t39 + t42) * t386 + (t269 * t372 - t270 * t373 + m(5) * ((-pkin(1) - t236) * t488 + (-t93 * rSges(5,3) + t394 * t94) * t289) - qJD(2) * t329 + t205 * t291 + t207 * t288 - qJD(4) * t331 + t204 * t290 + t206 * t287) * qJD(1) + (t31 + t40 + t45) * t389 + (t38 + t44 + t30) * t385 + (t41 + t43) * t390 - t146 * t67 / 0.2e1 - m(2) * (-g(1) * t235 + g(2) * t238); ((t62 * t289 + t63 * t292) * qJD(1) + t363) * t390 + ((t58 * t289 + t59 * t292) * qJD(1) + t361) * t389 + ((-t407 * t458 - t411) * t292 + (t311 + (t292 * t457 + t297) * qJD(2)) * t289) * t388 + ((-t408 * t457 + t411) * t289 + (t311 + (t289 * t458 + t297) * qJD(2)) * t292) * t391 + (qJD(1) * t45 + qJD(2) * t361 + qJDD(1) * t90 + t216 * t59 + t217 * t58) * t501 + t310 + ((t288 * t424 + t291 * t423) * qJD(1) + (t313 * t288 + t291 * t515) * qJD(2)) * t499 + (t289 * t41 - t292 * t40 + (t87 * t289 + t292 * t88) * qJD(1)) * t498 + (qJD(1) * t43 + qJD(2) * t363 + qJDD(1) * t92 + t216 * t63 + t217 * t62) * t502 + (t289 * t88 - t292 * t87) * t486 + t357 * t508 + t359 * t507 + t29 * t393 + t31 * t392 + (-g(2) * (-pkin(2) * t450 + t176) - g(3) * (t262 - t367) + (t378 * t531 + t13 * t533 + g(1) * t449 + (-t487 + (qJD(1) * t64 - t33) * t292 - t418 * t531 - t409 * t64) * t288) * pkin(2) + t526) * m(4) + (sin(qJ(3)) * pkin(2) / pkin(4) * (t30 * t392 + (qJD(1) * t42 + qJD(4) * t364 + qJDD(1) * t91 + t214 * t61 + t215 * t60) * t502 + t358 * t510 + ((t60 * t289 + t61 * t292) * qJD(1) + t364) * t386 + t28 * t393 + (qJD(1) * t44 + qJD(4) * t362 + qJDD(1) * t89 + t214 * t57 + t215 * t56) * t501 + t360 * t509 + ((t56 * t289 + t57 * t292) * qJD(1) + t362) * t385 + (t289 * t86 - t292 * t85) * t486 + (t289 * t39 - t292 * t38 + (t85 * t289 + t292 * t86) * qJD(1)) * t498 + ((-t406 * t459 + t412) * t289 + (t312 + (t289 * t460 + t298) * qJD(4)) * t292) * t387 + ((-t405 * t460 - t412) * t292 + (t312 + (t292 * t459 + t298) * qJD(4)) * t289) * t384 + ((t287 * t426 + t290 * t425) * qJD(1) + (t314 * t287 + t290 * t514) * qJD(4)) * t499 + ((qJD(4) * t344 + t166 * t214 - t168 * t215) * t335 + t81 * ((t166 * t292 - t168 * t289) * qJD(1) + t344) + t355 * t210 + (-t48 * t289 - t47 * t292 + (-t292 * t94 + t490) * qJD(1)) * t233 - (t190 * t93 - t492) * qJD(1) - (t81 * (-t190 * t289 - t192 * t292) + t355 * t236) * qJD(4) + g(1) * t192 + g(2) * t190 - g(3) * t236) * m(5)) + (-pkin(3) * t265 + pkin(2) * sin(t440)) / pkin(3) * (t310 + (-g(2) * t176 - g(3) * t525 + t526) * m(4))) / t265 + ((qJD(2) * t343 + t167 * t216 - t169 * t217) * t334 + t82 * ((t167 * t292 - t169 * t289) * qJD(1) + t343) + t342 * t211 + (-t49 * t289 - t50 * t292 + (t111 * t292 + t112 * t289) * qJD(1)) * t234 - (t111 * t193 + t112 * t191) * qJD(1) - (t82 * (-t191 * t289 - t193 * t292) + t342 * t237) * qJD(2) + g(1) * t193 + g(2) * t191 - g(3) * t237) * m(3);];
tau = t1(:);
