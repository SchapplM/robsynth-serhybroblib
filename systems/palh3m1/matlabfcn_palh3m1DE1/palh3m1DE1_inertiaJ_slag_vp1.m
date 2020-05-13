% Calculate joint inertia matrix for
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1DE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_inertiaJ_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1DE1_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 11:17:17
% EndTime: 2020-04-18 12:07:06
% DurationCPUTime: 501.57s
% Computational Cost: add. (13022311->653), mult. (19603973->1095), div. (905272->22), fcn. (12417032->36), ass. (0->448)
t350 = sin(qJ(3));
t351 = sin(qJ(2));
t356 = cos(qJ(3));
t357 = cos(qJ(2));
t315 = t350 * t351 - t356 * t357;
t393 = t350 * t357 + t351 * t356;
t288 = -rSges(4,1) * t393 + rSges(4,2) * t315;
t508 = m(4) * t288;
t366 = pkin(5) ^ 2;
t371 = pkin(1) ^ 2;
t353 = sin(pkin(16));
t519 = cos(pkin(16));
t316 = t351 * t353 - t357 * t519;
t512 = pkin(5) * t316;
t547 = -2 * pkin(1);
t468 = t512 * t547 + t371;
t308 = t366 + t468;
t464 = pkin(2) ^ 2 - pkin(6) ^ 2;
t293 = t308 - t464;
t312 = pkin(1) * t316 - pkin(5);
t318 = t351 * t519 + t357 * t353;
t541 = -pkin(6) - pkin(2);
t289 = (pkin(5) - t541) * (pkin(5) + t541) + t468;
t542 = pkin(2) - pkin(6);
t290 = (pkin(5) - t542) * (pkin(5) + t542) + t468;
t373 = sqrt(-t290 * t289);
t474 = t318 * t373;
t263 = -pkin(1) * t474 - t293 * t312;
t266 = pkin(1) * t318 * t293 - t312 * t373;
t354 = sin(pkin(15));
t305 = 0.1e1 / t308;
t365 = 0.1e1 / pkin(6);
t477 = t305 * t365;
t359 = cos(pkin(15));
t521 = t359 / 0.2e1;
t252 = (t263 * t521 + t266 * t354 / 0.2e1) * t477;
t255 = (t266 * t521 - t263 * t354 / 0.2e1) * t477;
t241 = atan2(t255, t252);
t237 = sin(t241);
t563 = Icges(7,5) * t237;
t238 = cos(t241);
t562 = Icges(7,6) * t238;
t561 = t562 / 0.2e1 + t563 / 0.2e1;
t294 = t308 + t464;
t311 = pkin(1) - t512;
t511 = pkin(5) * t318;
t462 = pkin(1) * t511;
t480 = 0.2e1 / t373 * (t289 + t290) * t462;
t441 = -t480 / 0.2e1;
t475 = t316 * t373;
t245 = (t475 + (t311 * t547 - t294 + t441) * t318) * pkin(5);
t545 = -0.2e1 * t318 ^ 2;
t247 = t311 * t480 / 0.2e1 + t366 * pkin(1) * t545 + (-t294 * t316 - t474) * pkin(5);
t264 = -pkin(5) * t474 + t294 * t311;
t265 = t294 * t511 + t311 * t373;
t347 = cos(pkin(19));
t370 = 0.1e1 / pkin(2);
t476 = t305 * t370;
t345 = sin(pkin(19));
t525 = t345 / 0.2e1;
t250 = (-t264 * t347 / 0.2e1 + t265 * t525) * t476;
t248 = 0.1e1 / t250 ^ 2;
t251 = (t265 * t347 / 0.2e1 + t264 * t525) * t476;
t434 = 0.1e1 / t308 ^ 2 * t462;
t413 = t265 * t434;
t414 = t264 * t434;
t440 = t305 * t525;
t479 = t305 * t347;
t527 = t247 / 0.2e1;
t528 = -t245 / 0.2e1;
t560 = ((t245 * t440 + t345 * t414 + t347 * t413 + t479 * t527) / t250 - (t247 * t440 + t345 * t413 - t347 * t414 + t479 * t528) * t251 * t248) / (t248 * t251 ^ 2 + 0.1e1) * t370 + 0.1e1;
t368 = pkin(3) ^ 2;
t367 = pkin(4) ^ 2;
t482 = t265 * t350;
t483 = t264 * t356;
t253 = (-t483 / 0.2e1 + t482 / 0.2e1) * t476;
t481 = t265 * t356;
t484 = t264 * t350;
t254 = (t481 / 0.2e1 + t484 / 0.2e1) * t476;
t340 = pkin(18) + pkin(19);
t331 = sin(t340);
t332 = cos(t340);
t230 = t253 * t332 + t254 * t331;
t515 = t230 * pkin(4);
t558 = 2 * pkin(3);
t471 = t515 * t558 + t367;
t225 = t368 + t471;
t465 = pkin(8) ^ 2 - pkin(10) ^ 2;
t222 = t225 - t465;
t227 = pkin(3) * t230 + pkin(4);
t540 = (-pkin(8) - pkin(10));
t219 = ((pkin(3) - t540) * (pkin(3) + t540)) + t471;
t539 = (pkin(10) - pkin(8));
t220 = ((pkin(3) - t539) * (pkin(3) + t539)) + t471;
t372 = sqrt(-t220 * t219);
t396 = t253 * t331 - t254 * t332;
t548 = t396 * t372;
t174 = -pkin(3) * t548 + t222 * t227;
t176 = pkin(3) * t222 * t396 + t227 * t372;
t344 = cos(pkin(17));
t223 = 0.1e1 / t225;
t361 = 0.1e1 / pkin(10);
t487 = t223 * t361;
t343 = sin(pkin(17));
t526 = t343 / 0.2e1;
t158 = (-t174 * t344 / 0.2e1 + t176 * t526) * t487;
t159 = (t176 * t344 / 0.2e1 + t174 * t526) * t487;
t452 = atan2(t159, t158);
t422 = cos(t452);
t376 = t393 * t422;
t352 = sin(qJ(1));
t341 = t352 ^ 2;
t358 = cos(qJ(1));
t342 = t358 ^ 2;
t559 = t341 + t342;
t221 = t225 + t465;
t226 = -pkin(3) - t515;
t173 = -pkin(4) * t548 - t221 * t226;
t516 = pkin(4) * t396;
t175 = t221 * t516 - t226 * t372;
t363 = 0.1e1 / pkin(8);
t529 = t223 / 0.2e1;
t442 = t363 * t529;
t375 = atan2(t175 * t442, t173 * t442);
t160 = cos(t375);
t346 = sin(pkin(18));
t348 = cos(pkin(18));
t374 = sin(t375);
t147 = t160 * t346 - t348 * t374;
t148 = -t348 * t160 - t346 * t374;
t240 = atan2(t251, t250);
t234 = sin(t240);
t235 = cos(t240);
t210 = -t234 * t351 + t235 * t357;
t211 = t234 * t357 + t235 * t351;
t124 = -t147 * t211 + t148 * t210;
t125 = t147 * t210 + t148 * t211;
t96 = Icges(9,4) * t125 + Icges(9,2) * t124;
t557 = t96 / 0.2e1;
t97 = Icges(9,1) * t125 + Icges(9,4) * t124;
t556 = t97 / 0.2e1;
t149 = sin(t452);
t140 = -t149 * t393 - t315 * t422;
t141 = t315 * t149 - t376;
t112 = Icges(5,4) * t141 - Icges(5,2) * t140;
t555 = -t112 / 0.2e1;
t113 = Icges(5,1) * t141 - Icges(5,4) * t140;
t554 = t113 / 0.2e1;
t553 = t124 / 0.2e1;
t552 = t125 / 0.2e1;
t551 = -t140 / 0.2e1;
t550 = t141 / 0.2e1;
t523 = t352 / 0.2e1;
t522 = -t358 / 0.2e1;
t421 = 0.2e1 * t559;
t392 = t421 / 0.2e1;
t549 = t352 * t358;
t309 = t315 * t352;
t136 = t149 * t309 - t352 * t376;
t310 = t315 * t358;
t138 = t149 * t310 - t358 * t376;
t383 = t393 * t352;
t137 = t149 * t383 + t309 * t422;
t349 = sin(qJ(4));
t355 = cos(qJ(4));
t128 = -t137 * t349 - t355 * t358;
t129 = t137 * t355 - t349 * t358;
t69 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t136;
t71 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t136;
t73 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t136;
t23 = t128 * t71 + t129 * t73 + t136 * t69;
t382 = t393 * t358;
t139 = t149 * t382 + t310 * t422;
t130 = -t139 * t349 + t352 * t355;
t131 = t139 * t355 + t349 * t352;
t70 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t138;
t72 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t138;
t74 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t138;
t24 = t128 * t72 + t129 * t74 + t136 * t70;
t77 = Icges(6,3) * t140 + (Icges(6,5) * t355 - Icges(6,6) * t349) * t141;
t78 = Icges(6,6) * t140 + (Icges(6,4) * t355 - Icges(6,2) * t349) * t141;
t79 = Icges(6,5) * t140 + (Icges(6,1) * t355 - Icges(6,4) * t349) * t141;
t29 = t128 * t78 + t129 * t79 + t136 * t77;
t1 = t136 * t23 + t138 * t24 + t140 * t29;
t544 = t1 / 0.2e1;
t25 = t130 * t71 + t131 * t73 + t138 * t69;
t26 = t130 * t72 + t131 * t74 + t138 * t70;
t30 = t130 * t78 + t131 * t79 + t138 * t77;
t2 = t136 * t25 + t138 * t26 + t140 * t30;
t543 = t2 / 0.2e1;
t524 = t350 / 0.2e1;
t217 = ((t245 * t524 + t356 * t527) * t305 + (t481 + t484) * t434) * t370;
t218 = ((t247 * t524 + t356 * t528) * t305 + (t482 - t483) * t434) * t370;
t192 = -t217 * t331 - t218 * t332;
t390 = pkin(4) * (t219 + t220) * t558;
t164 = t192 * t390;
t191 = -t217 * t332 + t218 * t331;
t194 = 0.1e1 / t372;
t530 = -t194 / 0.2e1;
t446 = t396 * t530;
t381 = t164 * t446 - t191 * t372;
t435 = -0.2e1 * pkin(4) * t227 - t222;
t142 = (t192 * t435 + t381) * pkin(3);
t437 = -0.2e1 * t368 * t516;
t447 = t194 * t227 / 0.2e1;
t490 = t192 * t372;
t143 = t164 * t447 + t192 * t437 + (t191 * t222 - t490) * pkin(3);
t156 = 0.1e1 / t158;
t224 = 0.1e1 / t225 ^ 2;
t517 = pkin(4) * t224;
t461 = pkin(3) * t517;
t432 = t344 * t461;
t417 = t192 * t432;
t433 = t343 * t461;
t418 = t192 * t433;
t488 = t223 * t344;
t443 = t488 / 0.2e1;
t444 = -t488 / 0.2e1;
t445 = t223 * t526;
t157 = 0.1e1 / t158 ^ 2;
t491 = t157 * t159;
t492 = 0.1e1 / (t157 * t159 ^ 2 + 0.1e1) * t361;
t538 = ((t142 * t445 + t143 * t443 + t174 * t418 + t176 * t417) * t156 - (t142 * t444 + t143 * t445 - t174 * t417 + t176 * t418) * t491) * t492 + 0.1e1;
t177 = t396 * t390;
t380 = t177 * t446 - t230 * t372;
t153 = (t396 * t435 + t380) * pkin(3);
t154 = t177 * t447 + t396 * t437 + (t222 * t230 - t548) * pkin(3);
t415 = t396 * t432;
t416 = t396 * t433;
t537 = ((t153 * t445 + t154 * t443 + t174 * t416 + t176 * t415) * t156 - (t153 * t444 + t154 * t445 - t174 * t415 + t176 * t416) * t491) * t492 + 0.1e1;
t98 = rSges(9,1) * t125 + rSges(9,2) * t124;
t536 = m(9) * t98;
t535 = t136 / 0.2e1;
t534 = t138 / 0.2e1;
t533 = t140 / 0.2e1;
t518 = pkin(1) * t351;
t514 = pkin(4) * t393;
t513 = pkin(4) * t352;
t510 = pkin(9) * t137;
t323 = rSges(3,1) * t351 + rSges(3,2) * t357;
t509 = m(3) * t323;
t216 = rSges(7,1) * t237 + rSges(7,2) * t238;
t507 = m(7) * t216;
t307 = t310 * pkin(4);
t506 = t357 * pkin(1);
t505 = t141 * t355 * t79 + t140 * t77;
t504 = rSges(7,1) * t238;
t503 = rSges(7,2) * t237;
t502 = t349 * t78;
t501 = t358 * rSges(3,3);
t500 = t358 * rSges(7,3);
t407 = -rSges(6,1) * t129 - rSges(6,2) * t128;
t75 = rSges(6,3) * t136 - t407;
t499 = pkin(11) * t136 + t510 + t75;
t76 = t131 * rSges(6,1) + t130 * rSges(6,2) + t138 * rSges(6,3);
t498 = -t139 * pkin(9) - pkin(11) * t138 - t76;
t80 = rSges(6,3) * t140 + (rSges(6,1) * t355 - rSges(6,2) * t349) * t141;
t497 = pkin(9) * t141 + pkin(11) * t140 + t80;
t496 = Icges(3,4) * t351;
t495 = Icges(3,4) * t357;
t208 = t211 * t358;
t489 = t208 * t346;
t478 = t305 * t354;
t473 = t351 * t358;
t472 = t357 * t358;
t470 = t352 * rSges(7,3) + t358 * t504;
t408 = -t309 * rSges(4,1) + t358 * rSges(4,3);
t450 = t310 * rSges(4,1) + rSges(4,2) * t382 + t352 * rSges(4,3);
t267 = t352 * (rSges(4,2) * t383 - t408) + t358 * t450;
t469 = t358 * t307 + t309 * t513;
t467 = t559 * t506;
t339 = t358 * pkin(13);
t466 = pkin(1) * t472 + t339;
t165 = t560 * t352;
t460 = t352 * t518;
t459 = pkin(1) * t473;
t458 = t393 * t513;
t457 = t358 * t514;
t27 = t140 * t69 + (-t349 * t71 + t355 * t73) * t141;
t456 = t27 / 0.2e1 + t29 / 0.2e1;
t28 = t140 * t70 + (-t349 * t72 + t355 * t74) * t141;
t455 = t28 / 0.2e1 + t30 / 0.2e1;
t206 = t211 * t352;
t207 = t210 * t352;
t120 = -t147 * t207 - t148 * t206;
t121 = -t147 * t206 + t148 * t207;
t88 = Icges(9,4) * t121 + Icges(9,2) * t120 - Icges(9,6) * t358;
t90 = Icges(9,1) * t121 + Icges(9,4) * t120 - Icges(9,5) * t358;
t95 = Icges(9,5) * t125 + Icges(9,6) * t124;
t454 = t120 * t557 + t121 * t556 + t95 * t522 + t90 * t552 + t88 * t553;
t209 = t210 * t358;
t122 = -t147 * t209 - t148 * t208;
t123 = -t147 * t208 + t148 * t209;
t89 = Icges(9,4) * t123 + Icges(9,2) * t122 + Icges(9,6) * t352;
t91 = Icges(9,1) * t123 + Icges(9,4) * t122 + Icges(9,5) * t352;
t453 = t122 * t557 + t123 * t556 + t95 * t523 + t91 * t552 + t89 * t553;
t93 = t123 * rSges(9,1) + t122 * rSges(9,2) + t352 * rSges(9,3);
t108 = t139 * rSges(5,1) - t138 * rSges(5,2) + t352 * rSges(5,3);
t451 = t209 * rSges(8,1) - t208 * rSges(8,2) + t352 * rSges(8,3);
t449 = t307 + t466;
t448 = t226 * t530;
t439 = t305 * t521;
t438 = -pkin(13) - t506;
t436 = -t288 - t518;
t276 = Icges(4,5) * t309 + Icges(4,6) * t383 - Icges(4,3) * t358;
t277 = Icges(4,5) * t310 + Icges(4,6) * t382 + Icges(4,3) * t352;
t278 = Icges(4,4) * t309 + Icges(4,2) * t383 - Icges(4,6) * t358;
t279 = Icges(4,4) * t310 + Icges(4,2) * t382 + Icges(4,6) * t352;
t280 = Icges(4,1) * t309 + Icges(4,4) * t383 - Icges(4,5) * t358;
t281 = Icges(4,1) * t310 + Icges(4,4) * t382 + Icges(4,5) * t352;
t426 = -t358 * (-(-t358 * t276 + t309 * t280) * t358 + (t309 * t281 + t279 * t383 + (-t278 * t393 - t277) * t358) * t352) + t352 * ((t352 * t277 + t310 * t281) * t352 + (-t310 * t280 - t278 * t382 + (t279 * t393 - t276) * t352) * t358);
t172 = 0.1e1 / t173 ^ 2;
t425 = pkin(8) / (t172 * t175 ^ 2 + 0.1e1) * t225 * t363;
t285 = -Icges(4,5) * t393 + Icges(4,6) * t315;
t286 = -Icges(4,4) * t393 + Icges(4,2) * t315;
t287 = -Icges(4,1) * t393 + Icges(4,4) * t315;
t424 = (t279 * t315 - t281 * t393 + t352 * t285 + t286 * t382 + t310 * t287) * t523 + (t278 * t315 - t280 * t393 - t358 * t285 + t286 * t383 + t309 * t287) * t522;
t423 = t467 + t469;
t420 = t354 * t434;
t419 = t359 * t434;
t412 = t514 - t518;
t411 = t438 * t352;
t410 = 0.1e1 / t173 * t425;
t409 = rSges(3,1) * t357 - rSges(3,2) * t351;
t406 = -t503 + t504;
t405 = Icges(3,1) * t357 - t496;
t403 = -Icges(3,2) * t351 + t495;
t401 = Icges(3,5) * t357 - Icges(3,6) * t351;
t400 = Icges(7,5) * t238 - Icges(7,6) * t237;
t397 = -t206 * t346 - t207 * t348;
t391 = rSges(3,1) * t472 - rSges(3,2) * t473 + t352 * rSges(3,3);
t103 = Icges(5,4) * t137 - Icges(5,2) * t136 - Icges(5,6) * t358;
t105 = Icges(5,1) * t137 - Icges(5,4) * t136 - Icges(5,5) * t358;
t111 = Icges(5,5) * t141 - Icges(5,6) * t140;
t389 = t103 * t551 + t105 * t550 + t111 * t522 + t136 * t555 + t137 * t554 + t456;
t104 = Icges(5,4) * t139 - Icges(5,2) * t138 + Icges(5,6) * t352;
t106 = Icges(5,1) * t139 - Icges(5,4) * t138 + Icges(5,5) * t352;
t388 = t104 * t551 + t106 * t550 + t111 * t523 + t138 * t555 + t139 * t554 + t455;
t387 = t412 * t352;
t386 = t412 * t358;
t385 = pkin(3) * (t173 * t224 + t223 * t226);
t384 = pkin(4) * t172 * t175 * t425;
t107 = rSges(5,1) * t137 - rSges(5,2) * t136 - rSges(5,3) * t358;
t379 = rSges(8,1) * t207 - rSges(8,2) * t206 - rSges(8,3) * t358;
t92 = rSges(9,1) * t121 + rSges(9,2) * t120 - rSges(9,3) * t358;
t378 = -pkin(4) * t309 + t411;
t377 = pkin(3) * (-t223 * t367 * t396 + t175 * t517);
t325 = rSges(2,1) * t358 - rSges(2,2) * t352;
t324 = -rSges(2,1) * t352 - rSges(2,2) * t358;
t320 = Icges(3,5) * t351 + Icges(3,6) * t357;
t298 = Icges(3,3) * t352 + t358 * t401;
t297 = -Icges(3,3) * t358 + t352 * t401;
t292 = t339 + t391;
t291 = t501 + (-pkin(13) - t409) * t352;
t283 = t436 * t358;
t282 = t436 * t352;
t275 = t358 * t391 + (t352 * t409 - t501) * t352;
t270 = t450 + t466;
t269 = (-rSges(4,2) * t393 + t438) * t352 + t408;
t262 = t467 + t267;
t249 = 0.1e1 / t252 ^ 2;
t246 = t312 * t441 + t371 * pkin(5) * t545 + (-t293 * t316 - t474) * pkin(1);
t244 = (t475 + (0.2e1 * t312 * pkin(5) - t293 + t441) * t318) * pkin(1);
t213 = t562 + t563;
t201 = Icges(7,3) * t352 + t358 * t400;
t200 = -Icges(7,3) * t358 + t352 * t400;
t197 = (-pkin(7) - t503) * t358 + t470;
t196 = t500 + (pkin(7) - t406) * t352;
t195 = t209 * t348 * pkin(3);
t190 = rSges(8,1) * t211 + rSges(8,2) * t210;
t189 = Icges(8,1) * t211 + Icges(8,4) * t210;
t188 = Icges(8,4) * t211 + Icges(8,2) * t210;
t187 = Icges(8,5) * t211 + Icges(8,6) * t210;
t186 = (-t210 * t346 + t211 * t348) * pkin(3);
t185 = Icges(8,1) * t209 - Icges(8,4) * t208 + Icges(8,5) * t352;
t184 = Icges(8,1) * t207 - Icges(8,4) * t206 - Icges(8,5) * t358;
t183 = Icges(8,4) * t209 - Icges(8,2) * t208 + Icges(8,6) * t352;
t182 = Icges(8,4) * t207 - Icges(8,2) * t206 - Icges(8,6) * t358;
t181 = Icges(8,5) * t209 - Icges(8,6) * t208 + Icges(8,3) * t352;
t180 = Icges(8,5) * t207 - Icges(8,6) * t206 - Icges(8,3) * t358;
t179 = t451 + t466;
t178 = t411 - t379;
t170 = ((t246 * t439 + t266 * t419 - t244 * t478 / 0.2e1 - t263 * t420) / t252 - (t244 * t439 + t263 * t419 + t246 * t478 / 0.2e1 + t266 * t420) * t255 * t249) / (t249 * t255 ^ 2 + 0.1e1) * t365;
t166 = t560 * t358;
t162 = -t166 * t190 - t459;
t161 = -t165 * t190 - t460;
t152 = (t358 * (-t358 * t503 + t470) + (t352 * t406 - t500) * t352) * t170;
t151 = t165 * t379 + t166 * t451 + t467;
t117 = 0.2e1 * ((t177 * t448 + (t221 * t230 - t548) * pkin(4)) * t529 + t396 * t377) * t410 - 0.2e1 * ((-t221 * t396 + t380) * t529 + t396 * t385) * t384;
t114 = rSges(5,1) * t141 - rSges(5,2) * t140;
t102 = Icges(5,5) * t139 - Icges(5,6) * t138 + Icges(5,3) * t352;
t101 = Icges(5,5) * t137 - Icges(5,6) * t136 - Icges(5,3) * t358;
t100 = t449 + t108;
t99 = -t107 + t378;
t87 = Icges(9,5) * t123 + Icges(9,6) * t122 + Icges(9,3) * t352;
t86 = Icges(9,5) * t121 + Icges(9,6) * t120 - Icges(9,3) * t358;
t85 = pkin(3) * t489 + t195 + t466 + t93;
t84 = pkin(3) * t397 + t411 - t92;
t83 = 0.2e1 * ((t164 * t448 + (t191 * t221 - t490) * pkin(4)) * t529 + t192 * t377) * t410 - 0.2e1 * ((-t192 * t221 + t381) * t529 + t192 * t385) * t384;
t82 = (-t83 - t560) * t358;
t81 = t352 * t83 + t165;
t67 = t537 * t358;
t66 = t537 * t352;
t62 = t538 * t358;
t61 = t538 * t352;
t60 = -t166 * t186 + t82 * t98 - t459;
t59 = -t165 * t186 - t81 * t98 - t460;
t58 = t449 - t498;
t57 = -t510 + (-pkin(11) - rSges(6,3)) * t136 + t378 + t407;
t56 = -t114 * t67 + t457;
t55 = -t114 * t66 + t458;
t50 = -t114 * t62 + t386;
t49 = -t114 * t61 + t387;
t48 = (t352 * t92 + t358 * t93) * t117;
t45 = -t138 * t80 + t140 * t76;
t44 = t136 * t80 - t140 * t75;
t43 = t102 * t352 - t104 * t138 + t106 * t139;
t42 = t101 * t352 - t103 * t138 + t105 * t139;
t41 = -t102 * t358 - t104 * t136 + t106 * t137;
t40 = -t101 * t358 - t103 * t136 + t105 * t137;
t39 = -t136 * t76 + t138 * t75;
t36 = t122 * t89 + t123 * t91 + t352 * t87;
t35 = t122 * t88 + t123 * t90 + t352 * t86;
t34 = t120 * t89 + t121 * t91 - t358 * t87;
t33 = t120 * t88 + t121 * t90 - t358 * t86;
t32 = t166 * t195 + t81 * t92 - t82 * t93 + (-t165 * t397 + t166 * t489) * pkin(3) + t467;
t31 = t107 * t66 + t108 * t67 + t469;
t22 = -t497 * t67 + t457;
t21 = -t497 * t66 + t458;
t20 = (-t141 * t502 + t505) * t140;
t19 = t107 * t61 + t108 * t62 + t423;
t18 = -t497 * t62 + t386;
t17 = -t497 * t61 + t387;
t16 = (-t35 * t358 + t352 * t36) * t117;
t15 = (-t33 * t358 + t34 * t352) * t117;
t14 = t35 * t82 + t36 * t81;
t13 = t33 * t82 + t34 * t81;
t12 = -t42 * t67 + t43 * t66;
t11 = -t40 * t67 + t41 * t66;
t10 = -t498 * t67 + t499 * t66 + t469;
t9 = -t42 * t62 + t43 * t61;
t8 = -t40 * t62 + t41 * t61;
t7 = -t498 * t62 + t499 * t61 + t423;
t6 = -t25 * t67 + t26 * t66;
t5 = -t23 * t67 + t24 * t66;
t4 = -t25 * t62 + t26 * t61;
t3 = -t23 * t62 + t24 * t61;
t37 = [Icges(7,2) * t238 ^ 2 + (t57 ^ 2 + t58 ^ 2) * m(6) + t351 * (Icges(3,1) * t351 + t495) + t357 * (Icges(3,2) * t357 + t496) + (t100 ^ 2 + t99 ^ 2) * m(5) + (t269 ^ 2 + t270 ^ 2) * m(4) + t505 + t315 * t286 - t393 * t287 + t125 * t97 - t140 * t112 + (t84 ^ 2 + t85 ^ 2) * m(9) + (t178 ^ 2 + t179 ^ 2) * m(8) + t211 * t189 + t210 * t188 + t124 * t96 + (t113 - t502) * t141 + Icges(2,3) + m(7) * (t196 ^ 2 + t197 ^ 2) + m(3) * (t291 ^ 2 + t292 ^ 2) + m(2) * (t324 ^ 2 + t325 ^ 2) + (Icges(7,1) * t237 + 0.2e1 * Icges(7,4) * t238) * t237; ((-Icges(3,6) * t358 + t352 * t403) * t357 + (-Icges(3,5) * t358 + t352 * t405) * t351) * t522 + t342 * t320 / 0.2e1 + t453 * t81 + ((-t197 * t507 + t213 * t523 + t561 * t352) * t352 + (-t196 * t507 + (t213 / 0.2e1 + t561) * t358) * t358) * t170 + t454 * t82 + (t59 * t85 + t60 * t84) * m(9) + t424 + (-t292 * t509 + (Icges(3,6) * t352 + t358 * t403) * t357 / 0.2e1 + (Icges(3,5) * t352 + t358 * t405) * t351 / 0.2e1 + t320 * t523) * t352 - t291 * t358 * t509 + (t100 * t49 + t50 * t99) * m(5) + (t161 * t179 + t162 * t178) * m(8) + (t17 * t58 + t18 * t57) * m(6) + (t269 * t283 + t270 * t282) * m(4) + t388 * t61 - t389 * t62 + (t183 * t210 + t185 * t211 + t187 * t352 - t188 * t208 + t189 * t209) * t165 / 0.2e1 - (t182 * t210 + t184 * t211 - t187 * t358 - t188 * t206 + t189 * t207) * t166 / 0.2e1; (t32 ^ 2 + t59 ^ 2 + t60 ^ 2) * m(9) + t82 * t13 + (t262 ^ 2 + t282 ^ 2 + t283 ^ 2) * m(4) - t358 * (t342 * t297 - t298 * t549) + t165 * ((t181 * t352 - t183 * t208 + t185 * t209) * t165 - (t180 * t352 - t182 * t208 + t184 * t209) * t166) + (t17 ^ 2 + t18 ^ 2 + t7 ^ 2) * m(6) + (t19 ^ 2 + t49 ^ 2 + t50 ^ 2) * m(5) + t81 * t14 + (t151 ^ 2 + t161 ^ 2 + t162 ^ 2) * m(8) + m(7) * t152 ^ 2 - t166 * ((-t181 * t358 - t183 * t206 + t185 * t207) * t165 - (-t180 * t358 - t182 * t206 + t184 * t207) * t166) + t352 * (-t297 * t549 + t341 * t298) + m(3) * (t323 ^ 2 * t421 + 0.2e1 * t275 ^ 2) / 0.2e1 + t426 - (t3 + t8) * t62 + (t4 + t9) * t61 + (t352 * (-t200 * t549 + t341 * t201) - t358 * (t342 * t200 - t201 * t549) + m(7) * t216 ^ 2 * t392) * t170 ^ 2; (t21 * t58 + t22 * t57) * m(6) + (t100 * t55 + t56 * t99) * m(5) - t389 * t67 + t388 * t66 + (-t269 * t508 + (-t536 * t84 - t454) * t117) * t358 + (-t270 * t508 + (-t536 * t85 + t453) * t117) * t352 + t424; (t19 * t31 + t49 * t55 + t50 * t56) * m(5) + m(9) * t48 * t32 + t82 * t15 / 0.2e1 + m(4) * t267 * t262 + t81 * t16 / 0.2e1 + (t10 * t7 + t17 * t21 + t18 * t22) * m(6) - (t3 / 0.2e1 + t8 / 0.2e1) * t67 + (t4 / 0.2e1 + t9 / 0.2e1) * t66 - (t11 / 0.2e1 + t5 / 0.2e1) * t62 + (t6 / 0.2e1 + t12 / 0.2e1) * t61 + (-t282 * t352 - t283 * t358) * t508 + ((-t60 * t536 - t13 / 0.2e1) * t358 + (t14 / 0.2e1 - t59 * t536) * t352) * t117 + t426; (t31 ^ 2 + t55 ^ 2 + t56 ^ 2) * m(5) + (t10 ^ 2 + t21 ^ 2 + t22 ^ 2) * m(6) - (t5 + t11) * t67 + (t6 + t12) * t66 + (-t15 * t358 + t16 * t352) * t117 + t426 + (t117 ^ 2 * t392 * t98 ^ 2 + t48 ^ 2) * m(9) + (t288 ^ 2 * t392 + t267 ^ 2) * m(4); t20 + (t44 * t57 + t45 * t58) * m(6) + t455 * t138 + t456 * t136; t61 * t543 + t3 * t535 + t4 * t534 + (-t27 * t62 + t28 * t61) * t533 + (t17 * t45 + t18 * t44 + t39 * t7) * m(6) - t62 * t544; -t67 * t544 + t66 * t543 + (t10 * t39 + t21 * t45 + t22 * t44) * m(6) + t5 * t535 + t6 * t534 + (-t27 * t67 + t28 * t66) * t533; (t39 ^ 2 + t44 ^ 2 + t45 ^ 2) * m(6) + t140 * (t27 * t136 + t28 * t138 + t20) + t136 * t1 + t138 * t2;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t37(1), t37(2), t37(4), t37(7); t37(2), t37(3), t37(5), t37(8); t37(4), t37(5), t37(6), t37(9); t37(7), t37(8), t37(9), t37(10);];
Mq = res;
