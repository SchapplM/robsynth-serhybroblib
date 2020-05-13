% Calculate vector of inverse dynamics joint torques for
% fourbar1turnOL
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnOL_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:41:03
% DurationCPUTime: 11.16s
% Computational Cost: add. (6393->663), mult. (11306->954), div. (0->0), fcn. (9039->8), ass. (0->385)
t287 = sin(qJ(1));
t284 = qJ(2) + qJ(3);
t267 = sin(t284);
t268 = cos(t284);
t343 = Icges(4,5) * t267 + Icges(4,6) * t268;
t170 = t343 * t287;
t290 = cos(qJ(1));
t171 = t343 * t290;
t281 = qJD(2) + qJD(3);
t219 = t281 * t287;
t363 = rSges(4,1) * t267 + rSges(4,2) * t268;
t530 = t363 * t219;
t453 = t268 * t287;
t455 = t267 * t287;
t469 = Icges(4,6) * t290;
t134 = Icges(4,4) * t453 - Icges(4,2) * t455 + t469;
t478 = Icges(4,4) * t268;
t350 = Icges(4,1) * t267 + t478;
t529 = -t350 * t287 - t134;
t347 = -Icges(4,2) * t267 + t478;
t135 = Icges(4,6) * t287 - t290 * t347;
t528 = -t350 * t290 + t135;
t452 = t268 * t290;
t454 = t267 * t290;
t177 = rSges(4,1) * t454 + rSges(4,2) * t452;
t479 = Icges(4,4) * t267;
t346 = Icges(4,2) * t268 + t479;
t351 = Icges(4,1) * t268 - t479;
t527 = -t346 + t351;
t526 = -t347 - t350;
t381 = -rSges(4,2) * t455 + rSges(4,3) * t290;
t140 = rSges(4,1) * t453 + t381;
t220 = t281 * t290;
t278 = t287 * rSges(4,3);
t420 = rSges(4,2) * t454 + t278;
t325 = rSges(4,1) * t452 - t420;
t490 = pkin(2) * qJD(2);
t289 = cos(qJ(2));
t416 = t287 ^ 2 + t290 ^ 2;
t524 = t289 * t416;
t55 = -t219 * t140 - t220 * t325 + t490 * t524;
t525 = qJD(2) * t55;
t261 = t267 * rSges(4,2);
t494 = rSges(4,1) * t268;
t523 = t261 - t494;
t288 = cos(qJ(4));
t275 = Icges(5,4) * t288;
t285 = sin(qJ(4));
t345 = -Icges(5,2) * t285 + t275;
t229 = Icges(5,1) * t285 + t275;
t276 = Icges(3,4) * t289;
t286 = sin(qJ(2));
t348 = -Icges(3,2) * t286 + t276;
t231 = Icges(3,1) * t286 + t276;
t405 = qJD(2) * t290;
t393 = t286 * t405;
t408 = qJD(1) * t287;
t65 = -pkin(2) * (t289 * t408 + t393) + qJD(1) * t140 + t220 * t363;
t224 = Icges(3,5) * t289 - Icges(3,6) * t286;
t223 = Icges(3,5) * t286 + Icges(3,6) * t289;
t319 = qJD(2) * t223;
t480 = Icges(3,4) * t286;
t232 = Icges(3,1) * t289 - t480;
t165 = Icges(3,5) * t287 + t232 * t290;
t161 = Icges(3,6) * t287 + t290 * t348;
t460 = t161 * t286;
t334 = -t165 * t289 + t460;
t467 = Icges(3,3) * t290;
t522 = -t290 * t319 + (-t224 * t287 + t334 + t467) * qJD(1);
t449 = t286 * t287;
t257 = Icges(3,4) * t449;
t446 = t287 * t289;
t476 = Icges(3,5) * t290;
t164 = Icges(3,1) * t446 - t257 - t476;
t470 = Icges(3,6) * t290;
t160 = Icges(3,4) * t446 - Icges(3,2) * t449 - t470;
t461 = t160 * t286;
t335 = -t164 * t289 + t461;
t157 = Icges(3,3) * t287 + t224 * t290;
t412 = qJD(1) * t157;
t521 = qJD(1) * t335 - t287 * t319 + t412;
t222 = Icges(5,5) * t288 - Icges(5,6) * t285;
t221 = Icges(5,5) * t285 + Icges(5,6) * t288;
t316 = qJD(4) * t221;
t477 = Icges(5,4) * t285;
t230 = Icges(5,1) * t288 - t477;
t163 = Icges(5,5) * t287 + t230 * t290;
t159 = Icges(5,6) * t287 + t290 * t345;
t462 = t159 * t285;
t336 = -t163 * t288 + t462;
t465 = Icges(5,3) * t290;
t520 = -t290 * t316 + (-t222 * t287 + t336 + t465) * qJD(1);
t451 = t285 * t287;
t256 = Icges(5,4) * t451;
t447 = t287 * t288;
t473 = Icges(5,5) * t290;
t162 = Icges(5,1) * t447 - t256 - t473;
t468 = Icges(5,6) * t290;
t158 = Icges(5,4) * t447 - Icges(5,2) * t451 - t468;
t463 = t158 * t285;
t337 = -t162 * t288 + t463;
t155 = Icges(5,3) * t287 + t222 * t290;
t413 = qJD(1) * t155;
t519 = qJD(1) * t337 - t287 * t316 + t413;
t344 = Icges(4,5) * t268 - Icges(4,6) * t267;
t242 = Icges(4,4) * t454;
t475 = Icges(4,5) * t287;
t137 = -Icges(4,1) * t452 + t242 + t475;
t464 = t137 * t268;
t466 = Icges(4,3) * t290;
t518 = -t281 * t171 + (t135 * t267 - t287 * t344 - t464 - t466) * qJD(1);
t241 = Icges(4,4) * t455;
t474 = Icges(4,5) * t290;
t136 = Icges(4,1) * t453 - t241 + t474;
t339 = t134 * t267 - t136 * t268;
t133 = Icges(4,3) * t287 - t290 * t344;
t414 = qJD(1) * t133;
t517 = qJD(1) * t339 + t170 * t281 + t414;
t156 = Icges(3,5) * t446 - Icges(3,6) * t449 - t467;
t58 = -t290 * t156 - t287 * t335;
t154 = Icges(5,5) * t447 - Icges(5,6) * t451 - t465;
t56 = -t290 * t154 - t287 * t337;
t331 = -t267 * t346 + t268 * t350;
t516 = qJD(1) * t331 + t344 * t281;
t225 = Icges(5,2) * t288 + t477;
t329 = t225 * t285 - t229 * t288;
t515 = t329 * qJD(1) + t222 * qJD(4);
t227 = Icges(3,2) * t289 + t480;
t327 = t227 * t286 - t231 * t289;
t514 = t327 * qJD(1) + t224 * qJD(2);
t402 = qJD(1) * qJD(2);
t263 = t290 * t402;
t216 = qJDD(2) * t287 + t263;
t401 = qJD(1) * qJD(3);
t145 = qJDD(3) * t287 + t290 * t401 + t216;
t262 = t287 * t402;
t146 = t287 * t401 + t262 + (-qJDD(2) - qJDD(3)) * t290;
t217 = -qJDD(2) * t290 + t262;
t291 = qJD(2) ^ 2;
t376 = t286 * t416;
t395 = qJD(1) * t494;
t407 = qJD(1) * t290;
t366 = rSges(4,3) * t407 + t177 * t281 + t287 * t395;
t394 = qJD(1) * t261;
t79 = -t287 * t394 + t366;
t80 = t530 + (t290 * t523 + t278) * qJD(1);
t13 = -t145 * t140 + t219 * t80 + t146 * t325 + t220 * t79 + ((t216 * t287 - t217 * t290) * t289 - t291 * t376) * pkin(2);
t176 = rSges(4,1) * t455 + rSges(4,2) * t453;
t315 = qJD(1) * t325;
t153 = t523 * t281;
t399 = qJDD(1) * t287;
t292 = qJD(1) ^ 2;
t415 = -t291 - t292;
t33 = -qJD(1) * t80 + qJDD(1) * t140 - t146 * t363 - t153 * t220 + ((t217 + t262) * t286 + (t290 * t415 - t399) * t289) * pkin(2);
t398 = qJDD(1) * t290;
t32 = t145 * t363 - t219 * t153 - qJDD(1) * t325 + qJD(1) * t79 + ((-t216 - t263) * t286 + (t287 * t415 + t398) * t289) * pkin(2);
t484 = t32 * t287;
t497 = pkin(2) * t289;
t251 = t407 * t497;
t396 = t286 * t490;
t64 = t287 * t396 - t251 + t315 - t530;
t513 = -g(1) * t177 + t13 * (-t287 * t140 - t290 * t325) - (-t484 - t33 * t290 + (t287 * t65 + t290 * t64) * qJD(1)) * t363 + t64 * (qJD(1) * t177 - t219 * t523) + (-t140 * t407 + t290 * t79 + (t315 + t80) * t287 - t219 * t176 - t177 * t220) * t55;
t512 = t287 * (-t227 * t290 + t165) - t290 * (-Icges(3,2) * t446 + t164 - t257);
t511 = t287 * (-t225 * t290 + t163) - t290 * (-Icges(5,2) * t447 + t162 - t256);
t510 = qJD(1) * t526 + t219 * (Icges(4,2) * t452 + t137 + t242) + t220 * (-Icges(4,2) * t453 + t136 - t241);
t509 = t145 / 0.2e1;
t508 = t146 / 0.2e1;
t400 = qJD(1) * qJD(4);
t214 = qJDD(4) * t287 + t290 * t400;
t507 = t214 / 0.2e1;
t215 = -qJDD(4) * t290 + t287 * t400;
t506 = t215 / 0.2e1;
t505 = t216 / 0.2e1;
t504 = t217 / 0.2e1;
t503 = -t219 / 0.2e1;
t502 = t219 / 0.2e1;
t501 = -t220 / 0.2e1;
t500 = t220 / 0.2e1;
t499 = t287 / 0.2e1;
t498 = -t290 / 0.2e1;
t496 = -qJD(1) / 0.2e1;
t495 = qJD(1) / 0.2e1;
t493 = rSges(5,1) * t288;
t492 = rSges(5,2) * t285;
t491 = rSges(5,2) * t288;
t233 = rSges(5,1) * t285 + t491;
t192 = t233 * t290;
t277 = t287 * rSges(5,3);
t445 = t288 * t290;
t417 = rSges(5,1) * t445 + t277;
t450 = t285 * t290;
t168 = -rSges(5,2) * t450 + t417;
t404 = qJD(4) * t287;
t94 = -t233 * t404 + (pkin(1) * t290 + t168) * qJD(1);
t489 = t192 * t94;
t488 = t287 * t64;
t258 = rSges(5,2) * t451;
t418 = t290 * rSges(5,3) + t258;
t166 = rSges(5,1) * t447 - t418;
t403 = qJD(4) * t290;
t93 = -t233 * t403 + (-pkin(1) * t287 - t166) * qJD(1);
t487 = t287 * t93;
t486 = t290 * rSges(3,3);
t485 = t290 * t93;
t483 = qJDD(1) / 0.2e1;
t459 = t221 * t287;
t458 = t221 * t290;
t457 = t223 * t287;
t456 = t223 * t290;
t448 = t286 * t290;
t444 = t289 * t290;
t132 = Icges(4,5) * t453 - Icges(4,6) * t455 + t466;
t443 = t290 * t132;
t67 = t346 * t455 - t350 * t453 - t171;
t440 = t67 * qJD(1);
t89 = -t287 * t329 - t458;
t439 = t89 * qJD(1);
t90 = -t287 * t327 - t456;
t438 = t90 * qJD(1);
t437 = t287 * t132 + t134 * t454;
t436 = t287 * t133 + t135 * t454;
t435 = -t287 * t154 - t162 * t445;
t434 = t287 * t155 + t163 * t445;
t433 = -t287 * t156 - t164 * t444;
t432 = t287 * t157 + t165 * t444;
t424 = -t225 + t230;
t423 = t229 + t345;
t422 = -t227 + t232;
t421 = t231 + t348;
t419 = rSges(5,3) * t407 + qJD(1) * t258;
t237 = rSges(3,1) * t289 - rSges(3,2) * t286;
t169 = rSges(3,3) * t287 + t237 * t290;
t411 = qJD(1) * t169;
t410 = qJD(1) * t222;
t409 = qJD(1) * t224;
t406 = qJD(2) * t287;
t392 = -pkin(1) - t493;
t391 = t408 / 0.2e1;
t390 = t407 / 0.2e1;
t389 = -t406 / 0.2e1;
t388 = t406 / 0.2e1;
t387 = -t405 / 0.2e1;
t386 = t405 / 0.2e1;
t385 = -t404 / 0.2e1;
t384 = t404 / 0.2e1;
t383 = -t403 / 0.2e1;
t382 = t403 / 0.2e1;
t380 = -(-t290 * t351 + t475) * qJD(1) + t529 * t281;
t379 = -(t287 * t351 + t474) * qJD(1) + t528 * t281;
t323 = t281 * t346;
t378 = -qJD(1) * t135 + t136 * t281 - t287 * t323;
t377 = t137 * t281 + t290 * t323 + (t287 * t347 + t469) * qJD(1);
t113 = t135 * t455;
t375 = t290 * t133 - t113;
t124 = t163 * t447;
t374 = t290 * t155 - t124;
t125 = t165 * t446;
t373 = t290 * t157 - t125;
t372 = t132 + t464;
t371 = t526 * t281;
t370 = t527 * t281;
t369 = -t154 + t462;
t368 = -t156 + t460;
t367 = -qJD(1) * t176 - t220 * t523;
t365 = t494 - t497;
t238 = rSges(2,1) * t290 - rSges(2,2) * t287;
t235 = rSges(2,1) * t287 + rSges(2,2) * t290;
t234 = rSges(3,1) * t286 + rSges(3,2) * t289;
t236 = -t492 + t493;
t317 = qJD(4) * t225;
t100 = qJD(1) * t159 - t287 * t317;
t318 = qJD(4) * t229;
t104 = qJD(1) * t163 - t287 * t318;
t85 = t158 * t288 + t162 * t285;
t301 = qJD(1) * t154 - qJD(4) * t85 - t100 * t285 + t104 * t288;
t103 = -t290 * t318 + (-t230 * t287 + t473) * qJD(1);
t86 = t159 * t288 + t163 * t285;
t99 = -t290 * t317 + (-t287 * t345 + t468) * qJD(1);
t302 = -qJD(4) * t86 + t103 * t288 - t285 * t99 + t413;
t362 = -(t519 * t287 + t301 * t290) * t290 + (t520 * t287 + t302 * t290) * t287;
t320 = qJD(2) * t227;
t101 = -t290 * t320 + (-t287 * t348 + t470) * qJD(1);
t321 = qJD(2) * t231;
t105 = -t290 * t321 + (-t232 * t287 + t476) * qJD(1);
t88 = t161 * t289 + t165 * t286;
t299 = -qJD(2) * t88 - t101 * t286 + t105 * t289 + t412;
t102 = qJD(1) * t161 - t287 * t320;
t106 = qJD(1) * t165 - t287 * t321;
t87 = t160 * t289 + t164 * t286;
t300 = qJD(1) * t156 - qJD(2) * t87 - t102 * t286 + t106 * t289;
t361 = -(t521 * t287 + t300 * t290) * t290 + (t522 * t287 + t299 * t290) * t287;
t360 = -(t301 * t287 - t519 * t290) * t290 + (t302 * t287 - t520 * t290) * t287;
t359 = -(t300 * t287 - t521 * t290) * t290 + (t299 * t287 - t522 * t290) * t287;
t57 = -t159 * t451 - t374;
t358 = t287 * t57 - t290 * t56;
t59 = -t161 * t449 - t373;
t357 = t287 * t59 - t290 * t58;
t60 = -t158 * t450 - t435;
t61 = -t159 * t450 + t434;
t356 = t287 * t61 - t290 * t60;
t62 = -t160 * t448 - t433;
t63 = -t161 * t448 + t432;
t355 = t287 * t63 - t290 * t62;
t353 = -t287 * t94 - t485;
t107 = -t403 * t491 + (-t285 * t403 - t288 * t408) * rSges(5,1) + t419;
t190 = t233 * t287;
t109 = -qJD(4) * t190 + (t236 * t290 + t277) * qJD(1);
t342 = t107 * t290 + t109 * t287;
t324 = qJD(2) * t234;
t108 = -t290 * t324 + (-t237 * t287 + t486) * qJD(1);
t110 = -t287 * t324 + t411;
t341 = t108 * t290 + t110 * t287;
t111 = t234 * t406 - t411;
t167 = rSges(3,1) * t446 - rSges(3,2) * t449 - t486;
t112 = -qJD(1) * t167 - t234 * t405;
t340 = t111 * t287 - t112 * t290;
t69 = t134 * t268 + t136 * t267;
t333 = t166 * t287 + t168 * t290;
t332 = t167 * t287 + t169 * t290;
t330 = t225 * t288 + t229 * t285;
t328 = t227 * t289 + t231 * t286;
t322 = t339 * t287;
t313 = -qJD(1) * t344 - t170 * t220 + t171 * t219;
t312 = t158 * t290 - t159 * t287;
t311 = t160 * t290 - t161 * t287;
t310 = (-t285 * t423 + t288 * t424) * qJD(1);
t309 = (-t286 * t421 + t289 * t422) * qJD(1);
t305 = -qJD(1) * t132 - t267 * t378 + t268 * t380;
t10 = t305 * t287 - t517 * t290;
t304 = t267 * t377 + t268 * t379 + t414;
t11 = t304 * t287 + t518 * t290;
t51 = -t322 + t443;
t52 = -t137 * t453 - t375;
t24 = t219 * t52 - t220 * t51 - t440;
t53 = t136 * t452 - t437;
t54 = -t137 * t452 + t436;
t68 = t290 * t331 - t170;
t66 = t68 * qJD(1);
t25 = t219 * t54 - t220 * t53 + t66;
t306 = t527 * qJD(1) + t528 * t219 - t529 * t220;
t294 = t510 * t267 + t306 * t268;
t303 = -qJD(1) * t343 + t267 * t371 + t268 * t370;
t34 = -t516 * t287 + t303 * t290;
t35 = t303 * t287 + t516 * t290;
t36 = t267 * t380 + t268 * t378;
t37 = t267 * t379 - t268 * t377;
t70 = -t135 * t268 - t137 * t267;
t8 = t517 * t287 + t305 * t290;
t9 = -t518 * t287 + t304 * t290;
t308 = (qJD(1) * t34 + qJDD(1) * t68 + t145 * t54 + t146 * t53 + t219 * t9 - t220 * t8) * t499 + (t306 * t267 - t510 * t268) * t496 + t24 * t391 + t25 * t390 + (qJD(1) * t35 - qJDD(1) * t67 - t10 * t220 + t11 * t219 + t145 * t52 + t146 * t51) * t498 + (t287 * t54 - t290 * t53) * t509 + (t287 * t52 - t290 * t51) * t508 + (t287 * t9 - t290 * t8 + (t287 * t53 + t290 * t54) * qJD(1)) * t502 + (-t10 * t290 + t11 * t287 + (t287 * t51 + t290 * t52) * qJD(1)) * t501 + (t287 * t70 - t290 * t69) * t483 + (t287 * t37 - t290 * t36 + (t287 * t69 + t290 * t70) * qJD(1)) * t495 + (t287 * t313 + t290 * t294) * t503 + (t287 * t294 - t290 * t313) * t500;
t204 = t345 * qJD(4);
t206 = t230 * qJD(4);
t298 = qJD(1) * t221 - qJD(4) * t330 - t204 * t285 + t206 * t288;
t205 = t348 * qJD(2);
t207 = t232 * qJD(2);
t297 = qJD(1) * t223 - qJD(2) * t328 - t205 * t286 + t207 * t289;
t296 = -t511 * t285 + t312 * t288;
t295 = -t512 * t286 + t311 * t289;
t211 = t237 * qJD(2);
t210 = t236 * qJD(4);
t193 = t234 * t290;
t191 = t234 * t287;
t92 = -t290 * t327 + t457;
t91 = -t290 * t329 + t459;
t84 = t92 * qJD(1);
t83 = t91 * qJD(1);
t82 = t332 * qJD(2);
t81 = t333 * qJD(4);
t50 = -qJD(1) * t110 - qJDD(1) * t167 - t211 * t405 + t217 * t234;
t49 = qJD(1) * t108 + qJDD(1) * t169 - t211 * t406 - t216 * t234;
t48 = -t210 * t404 + qJD(1) * t107 + qJDD(1) * t168 - t214 * t233 + (-t287 * t292 + t398) * pkin(1);
t47 = -t210 * t403 - qJD(1) * t109 - qJDD(1) * t166 + t215 * t233 + (-t290 * t292 - t399) * pkin(1);
t45 = t297 * t287 - t514 * t290;
t44 = t298 * t287 - t515 * t290;
t43 = t514 * t287 + t297 * t290;
t42 = t515 * t287 + t298 * t290;
t41 = -qJD(2) * t334 + t101 * t289 + t105 * t286;
t40 = -t335 * qJD(2) + t102 * t289 + t106 * t286;
t39 = -qJD(4) * t336 + t103 * t285 + t288 * t99;
t38 = -t337 * qJD(4) + t100 * t288 + t104 * t285;
t31 = qJD(2) * t355 + t84;
t30 = qJD(4) * t356 + t83;
t29 = qJD(2) * t357 + t438;
t28 = qJD(4) * t358 + t439;
t1 = [(-t108 * t111 - t110 * t112 + (-g(2) + t49) * t169 + (g(1) - t50) * t167) * m(3) + (t70 + t68) * t509 + (t86 + t91) * t507 + (t85 + t89) * t506 + (t88 + t92) * t505 + (t87 + t90) * t504 + (t37 + t34) * t502 + (t36 + t35 + t25) * t501 + (t346 * t268 + t350 * t267 + t328 + t330 + m(2) * (t235 ^ 2 + t238 ^ 2) + Icges(2,3)) * qJDD(1) + (-t438 + ((t290 * t368 - t432 + t63) * t290 + (t287 * t368 + t373 + t62) * t287) * qJD(2) + t29) * t389 + (t41 + t43) * t388 + (t40 + t45 + t31) * t387 + (t28 - t439 + ((t290 * t369 - t434 + t61) * t290 + (t287 * t369 + t374 + t60) * t287) * qJD(4)) * t385 + (t440 + (t54 - t322 - t436) * t220 + (t372 * t287 - t113 + t53) * t219 + ((t133 + t339) * t219 + t372 * t220) * t290 + t24) * t503 + (t94 * t419 + (t233 * t487 - t489) * qJD(4) + (t48 - g(2)) * ((pkin(1) - t492) * t290 + t417) + (t47 - g(1)) * (t287 * t392 + t418)) * m(5) + (t39 + t42) * t384 + (t38 + t44 + t30) * t383 + (t65 * (-t251 + (-t394 + t395) * t290) - t64 * (-pkin(2) * t393 + t366) + (t65 * (-t281 * t363 + t396) + (-t65 * rSges(4,3) - t64 * (-t261 - t497)) * qJD(1)) * t287 + (t32 - g(2)) * (-t290 * t365 + t420) + (t33 - g(1)) * (t287 * t365 + t381)) * m(4) + (t83 + ((t57 - t124 + (t155 + t463) * t290 + t435) * t290 + t434 * t287) * qJD(4)) * t382 + (t84 + ((t59 - t125 + (t157 + t461) * t290 + t433) * t290 + t432 * t287) * qJD(2)) * t386 + (m(5) * ((-pkin(1) - t236) * t485 + (-t93 * rSges(5,3) + t392 * t94) * t287) - qJD(2) * t327 + t205 * t289 + t207 * t286 - qJD(4) * t329 + t204 * t288 + t206 * t285 + t267 * t370 - t268 * t371) * qJD(1) + (t66 + (t52 + (-t136 * t290 + t137 * t287) * t268 + t375 + t437) * t220 + (t134 * t455 - t443 + t51 + (-t136 * t287 - t137 * t290) * t268 + t436) * t219) * t500 + t69 * t508 - t146 * t67 / 0.2e1 - m(2) * (-g(1) * t235 + g(2) * t238); ((t62 * t287 + t63 * t290) * qJD(1) + t361) * t388 + ((-t405 * t457 - t409) * t290 + (t309 + (t290 * t456 + t295) * qJD(2)) * t287) * t386 + ((-t406 * t456 + t409) * t287 + (t309 + (t287 * t457 + t295) * qJD(2)) * t290) * t389 + ((t286 * t422 + t289 * t421) * qJD(1) + (t311 * t286 + t512 * t289) * qJD(2)) * t496 + (t287 * t41 - t290 * t40 + (t87 * t287 + t290 * t88) * qJD(1)) * t495 + t31 * t390 + t355 * t505 + t357 * t504 + (t287 * t88 - t290 * t87) * t483 + ((t58 * t287 + t59 * t290) * qJD(1) + t359) * t387 + (qJD(1) * t45 + qJD(2) * t359 + qJDD(1) * t90 + t216 * t59 + t217 * t58) * t498 + (qJD(1) * t43 + qJD(2) * t361 + qJDD(1) * t92 + t216 * t63 + t217 * t62) * t499 + t308 + t29 * t391 + (-g(2) * (-pkin(2) * t449 + t176) - g(3) * (t261 - t365) + t153 * t488 + (t376 * t525 + t13 * t524 + g(1) * t448 + (-t64 * t407 - t484 + (qJD(1) * t64 - t33) * t290 - t416 * t525) * t286) * pkin(2) + (-t153 * t290 - t367) * t65 + t513) * m(4) + (-(t111 * t193 + t112 * t191) * qJD(1) - (t82 * (-t191 * t287 - t193 * t290) + t340 * t237) * qJD(2) + (qJD(2) * t341 + t167 * t216 - t169 * t217) * t332 + t82 * ((t167 * t290 - t169 * t287) * qJD(1) + t341) + t340 * t211 + (-t49 * t287 - t50 * t290 + (t111 * t290 + t112 * t287) * qJD(1)) * t234 + g(1) * t193 + g(2) * t191 - g(3) * t237) * m(3); t308 + (-g(2) * t176 - g(3) * t523 + (-t290 * t65 + t488) * t153 - t65 * t367 + t513) * m(4); t30 * t390 + (qJD(1) * t42 + qJD(4) * t362 + qJDD(1) * t91 + t214 * t61 + t215 * t60) * t499 + t356 * t507 + ((t60 * t287 + t61 * t290) * qJD(1) + t362) * t384 + t28 * t391 + (qJD(1) * t44 + qJD(4) * t360 + qJDD(1) * t89 + t214 * t57 + t215 * t56) * t498 + t358 * t506 + ((t56 * t287 + t57 * t290) * qJD(1) + t360) * t383 + (t287 * t86 - t290 * t85) * t483 + (t287 * t39 - t290 * t38 + (t85 * t287 + t290 * t86) * qJD(1)) * t495 + ((-t404 * t458 + t410) * t287 + (t310 + (t287 * t459 + t296) * qJD(4)) * t290) * t385 + ((-t403 * t459 - t410) * t290 + (t310 + (t290 * t458 + t296) * qJD(4)) * t287) * t382 + ((t285 * t424 + t288 * t423) * qJD(1) + (t312 * t285 + t511 * t288) * qJD(4)) * t496 + ((qJD(4) * t342 + t166 * t214 - t168 * t215) * t333 + t81 * ((t166 * t290 - t168 * t287) * qJD(1) + t342) + t353 * t210 + (-t48 * t287 - t47 * t290 + (-t290 * t94 + t487) * qJD(1)) * t233 - (t190 * t93 - t489) * qJD(1) - (t81 * (-t190 * t287 - t192 * t290) + t353 * t236) * qJD(4) + g(1) * t192 + g(2) * t190 - g(3) * t236) * m(5); 0;];
tau = t1;
