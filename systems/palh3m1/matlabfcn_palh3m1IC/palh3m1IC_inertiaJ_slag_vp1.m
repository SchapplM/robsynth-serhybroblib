% Calculate joint inertia matrix for
% palh3m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_inertiaJ_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_inertiaJ_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1IC_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1IC_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:31:37
% EndTime: 2020-04-20 17:31:46
% DurationCPUTime: 7.31s
% Computational Cost: add. (24009->549), mult. (21436->806), div. (128->6), fcn. (21524->28), ass. (0->323)
t448 = sin(qJ(1));
t452 = cos(qJ(1));
t444 = qJ(2) + qJ(3);
t426 = sin(t444);
t428 = cos(t444);
t489 = -Icges(4,5) * t428 + Icges(4,6) * t426;
t327 = -Icges(4,3) * t452 + t489 * t448;
t328 = Icges(4,3) * t448 + t489 * t452;
t442 = t452 ^ 2;
t549 = Icges(4,4) * t428;
t495 = Icges(4,2) * t426 - t549;
t332 = Icges(4,6) * t448 + t495 * t452;
t550 = Icges(4,4) * t426;
t501 = -Icges(4,1) * t428 + t550;
t336 = Icges(4,5) * t448 + t501 * t452;
t476 = t332 * t426 - t336 * t428;
t331 = -Icges(4,6) * t452 + t495 * t448;
t335 = -Icges(4,5) * t452 + t501 * t448;
t477 = -t331 * t426 + t335 * t428;
t530 = qJ(3) + qJ(4);
t432 = qJ(2) + t530;
t422 = cos(t432);
t450 = cos(qJ(5));
t531 = t450 * t452;
t446 = sin(qJ(5));
t533 = t448 * t446;
t357 = t422 * t533 - t531;
t532 = t448 * t450;
t534 = t446 * t452;
t358 = -t422 * t532 - t534;
t421 = sin(t432);
t537 = t421 * t448;
t254 = Icges(6,5) * t358 + Icges(6,6) * t357 - Icges(6,3) * t537;
t256 = Icges(6,4) * t358 + Icges(6,2) * t357 - Icges(6,6) * t537;
t258 = Icges(6,1) * t358 + Icges(6,4) * t357 - Icges(6,5) * t537;
t195 = -t254 * t537 + t256 * t357 + t258 * t358;
t359 = t422 * t534 + t532;
t360 = -t422 * t531 + t533;
t536 = t421 * t452;
t255 = Icges(6,5) * t360 + Icges(6,6) * t359 - Icges(6,3) * t536;
t257 = Icges(6,4) * t360 + Icges(6,2) * t359 - Icges(6,6) * t536;
t259 = Icges(6,1) * t360 + Icges(6,4) * t359 - Icges(6,5) * t536;
t196 = -t255 * t537 + t257 * t357 + t259 * t358;
t178 = -t195 * t452 + t196 * t448;
t488 = -Icges(5,5) * t422 + Icges(5,6) * t421;
t309 = -Icges(5,3) * t452 + t488 * t448;
t310 = Icges(5,3) * t448 + t488 * t452;
t547 = Icges(5,4) * t422;
t494 = Icges(5,2) * t421 - t547;
t312 = Icges(5,6) * t448 + t494 * t452;
t548 = Icges(5,4) * t421;
t500 = -Icges(5,1) * t422 + t548;
t314 = Icges(5,5) * t448 + t500 * t452;
t480 = t312 * t421 - t314 * t422;
t311 = -Icges(5,6) * t452 + t494 * t448;
t313 = -Icges(5,5) * t452 + t500 * t448;
t481 = -t311 * t421 + t313 * t422;
t528 = -t178 - t442 * t309 - (t480 * t448 + (-t310 + t481) * t452) * t448;
t581 = -t442 * t327 - (t476 * t448 + (-t328 + t477) * t452) * t448 + t528;
t443 = qJ(2) + qJ(7);
t430 = pkin(16) + t443;
t416 = sin(t430);
t419 = cos(t430);
t445 = sin(qJ(6));
t449 = cos(qJ(6));
t337 = 0.1e1 / (-t416 * t449 + t419 * t445) / pkin(5) / pkin(2);
t447 = sin(qJ(2));
t570 = pkin(1) * t447;
t386 = -pkin(2) * t416 - t570;
t451 = cos(qJ(2));
t439 = t451 * pkin(1);
t387 = pkin(2) * t419 + t439;
t228 = (-t386 * t419 - t387 * t416) * t337 * pkin(2);
t580 = t228 ^ 2;
t441 = t448 ^ 2;
t520 = t441 + t442;
t440 = -qJ(7) + pkin(15);
t429 = -qJ(8) + t440;
t423 = -qJ(2) + t429;
t409 = sin(t423);
t410 = cos(t423);
t504 = -rSges(9,1) * t410 - rSges(9,2) * t409;
t459 = t448 * rSges(9,3) + t504 * t452;
t555 = rSges(9,3) * t452;
t242 = t448 * (t504 * t448 - t555) + t452 * t459;
t356 = rSges(9,1) * t409 - rSges(9,2) * t410;
t485 = -Icges(9,5) * t410 - Icges(9,6) * t409;
t292 = -Icges(9,3) * t452 + t485 * t448;
t293 = Icges(9,3) * t448 + t485 * t452;
t541 = Icges(9,4) * t410;
t491 = -Icges(9,2) * t409 - t541;
t295 = Icges(9,6) * t448 + t491 * t452;
t542 = Icges(9,4) * t409;
t497 = -Icges(9,1) * t410 - t542;
t297 = Icges(9,5) * t448 + t497 * t452;
t482 = -t295 * t409 - t297 * t410;
t294 = -Icges(9,6) * t452 + t491 * t448;
t296 = -Icges(9,5) * t452 + t497 * t448;
t483 = t294 * t409 + t296 * t410;
t186 = t448 * (t441 * t293 + (t483 * t452 + (-t292 + t482) * t448) * t452);
t187 = t442 * t292 + (t482 * t448 + (-t293 + t483) * t452) * t448;
t514 = -t452 * t187 + t186;
t168 = m(9) * (t520 * t356 ^ 2 + t242 ^ 2) + t514;
t415 = sin(t429);
t371 = -t415 * pkin(7) + pkin(3) * sin(t440);
t418 = cos(t429);
t372 = -t418 * pkin(7) + pkin(3) * cos(t440);
t424 = pkin(14) + t530;
t413 = sin(t424);
t414 = cos(t424);
t229 = (t386 * t449 + t387 * t445) * pkin(5) * t337;
t286 = 0.1e1 / (t413 * t418 + t414 * t415) / pkin(9) / pkin(7);
t455 = t286 * t229;
t190 = (-t371 * t414 - t372 * t413) * pkin(9) * t455;
t420 = t439 + pkin(12);
t431 = -qJ(2) + t440;
t381 = pkin(3) * cos(t431) + t420;
t363 = t452 * t381;
t404 = t452 * t420;
t214 = t441 * (t381 - t420) + t452 * (-t404 + t363) + t242;
t523 = t441 * (-pkin(12) + t420) + t452 * (-pkin(12) * t452 + t404);
t207 = t214 + t523;
t569 = pkin(3) * sin(t431);
t389 = t448 * t569;
t518 = -t356 - t570;
t272 = t518 * t448 + t389;
t391 = t452 * t569;
t273 = t518 * t452 + t391;
t278 = -t356 * t448 + t389;
t279 = -t356 * t452 + t391;
t578 = -t229 * (m(9) * (t242 * t214 + (-t278 * t448 - t279 * t452) * t356) + t514) + m(9) * (t242 * t207 + (-t272 * t448 - t273 * t452) * t356) + t514;
t157 = t190 * t168 + t578;
t507 = -rSges(6,1) * t358 - rSges(6,2) * t357;
t262 = -rSges(6,3) * t537 - t507;
t524 = t360 * rSges(6,1) + t359 * rSges(6,2);
t263 = -rSges(6,3) * t536 + t524;
t567 = pkin(8) * t422;
t208 = t448 * t262 + t452 * t263 + t520 * (-pkin(10) * t421 - t567);
t388 = -pkin(4) * t428 + t420;
t370 = t452 * t388;
t526 = t441 * (t388 - t420) + t452 * (-t404 + t370);
t192 = t208 + t526;
t185 = t192 + t523;
t564 = rSges(5,1) * t422;
t464 = rSges(5,2) * t536 + t448 * rSges(5,3) - t452 * t564;
t508 = rSges(5,2) * t421 - t564;
t557 = rSges(5,3) * t452;
t249 = t448 * (t508 * t448 - t557) + t452 * t464;
t217 = t526 + t249;
t209 = t217 + t523;
t384 = -pkin(9) * t413 - sin(qJ(3)) * pkin(4);
t385 = pkin(9) * t414 + cos(qJ(3)) * pkin(4);
t218 = (t384 * t414 + t385 * t413) * t286 * pkin(9);
t535 = t426 * t452;
t565 = rSges(4,1) * t428;
t465 = rSges(4,2) * t535 + t448 * rSges(4,3) - t452 * t565;
t509 = rSges(4,2) * t426 - t565;
t558 = rSges(4,3) * t452;
t253 = t448 * (t509 * t448 - t558) + t452 * t465;
t223 = t523 + t253;
t407 = t448 * pkin(4) * t426;
t302 = rSges(6,3) * t422 + (-rSges(6,1) * t450 + rSges(6,2) * t446) * t421;
t525 = pkin(8) * t421 - pkin(10) * t422 - t302;
t512 = t525 - t570;
t236 = t512 * t448 + t407;
t408 = pkin(4) * t535;
t237 = t512 * t452 + t408;
t260 = t525 * t448;
t247 = t407 + t260;
t261 = t525 * t452;
t248 = t408 + t261;
t368 = -rSges(5,1) * t421 - rSges(5,2) * t422;
t517 = -t368 - t570;
t276 = t517 * t448 + t407;
t277 = t517 * t452 + t408;
t303 = -t368 * t448 + t407;
t304 = -t368 * t452 + t408;
t383 = -rSges(4,1) * t426 - rSges(4,2) * t428;
t515 = -t383 - t570;
t316 = t515 * t448;
t318 = t515 * t452;
t197 = -t254 * t536 + t359 * t256 + t360 * t258;
t198 = -t255 * t536 + t359 * t257 + t360 * t259;
t179 = -t197 * t452 + t198 * t448;
t529 = (t441 * t310 + t179 + (t481 * t452 + (-t309 + t480) * t448) * t452) * t448;
t519 = t448 * (t441 * t328 + (t477 * t452 + (-t327 + t476) * t448) * t452) + t529;
t456 = t581 * t452 + t519;
t579 = m(6) * (t185 * t192 + t236 * t247 + t237 * t248) + m(5) * (t209 * t217 + t276 * t303 + t277 * t304) + m(4) * (t253 * t223 + (-t316 * t448 - t318 * t452) * t383) + t456 + t157 * t218;
t425 = sin(t443);
t427 = cos(t443);
t505 = rSges(8,1) * t427 - rSges(8,2) * t425;
t576 = t448 * t452;
t574 = -(pkin(10) + rSges(6,3)) * t421 - t567;
t573 = t448 / 0.2e1;
t572 = -t452 / 0.2e1;
t566 = rSges(3,1) * t451;
t563 = rSges(7,1) * t449;
t561 = rSges(3,2) * t447;
t560 = rSges(7,2) * t445;
t556 = rSges(8,3) * t452;
t554 = t452 * rSges(3,3);
t553 = t452 * rSges(7,3);
t552 = Icges(3,4) * t447;
t551 = Icges(3,4) * t451;
t546 = Icges(7,4) * t445;
t545 = Icges(7,4) * t449;
t544 = Icges(8,4) * t425;
t543 = Icges(8,4) * t427;
t205 = t254 * t422 + (t256 * t446 - t258 * t450) * t421;
t540 = t205 * t452;
t206 = t255 * t422 + (t257 * t446 - t259 * t450) * t421;
t539 = t206 * t448;
t291 = Icges(6,5) * t422 + (-Icges(6,1) * t450 + Icges(6,4) * t446) * t421;
t538 = t291 * t450;
t289 = Icges(6,3) * t422 + (-Icges(6,5) * t450 + Icges(6,6) * t446) * t421;
t290 = Icges(6,6) * t422 + (-Icges(6,4) * t450 + Icges(6,2) * t446) * t421;
t527 = t421 * t446 * t290 + t422 * t289;
t463 = t448 * rSges(8,3) + t505 * t452;
t252 = t448 * (t505 * t448 - t556) + t452 * t463;
t522 = t448 * rSges(7,3) + t452 * t563;
t521 = t448 * rSges(3,3) + t452 * t566;
t382 = rSges(8,1) * t425 + rSges(8,2) * t427;
t516 = -t382 - t570;
t353 = Icges(9,5) * t409 - Icges(9,6) * t410;
t354 = -Icges(9,2) * t410 + t542;
t355 = Icges(9,1) * t409 - t541;
t471 = -t354 * t409 - t355 * t410;
t513 = (-t295 * t410 + t297 * t409 + t448 * t353 + t471 * t452) * t573 + (-t294 * t410 + t296 * t409 - t353 * t452 + t471 * t448) * t572;
t510 = -t561 + t566;
t506 = -t560 + t563;
t212 = -t289 * t537 + t290 * t357 + t291 * t358;
t171 = t212 * t422 + (-t195 * t448 - t196 * t452) * t421;
t213 = -t289 * t536 + t359 * t290 + t360 * t291;
t172 = t213 * t422 + (-t197 * t448 - t198 * t452) * t421;
t503 = t171 * t572 + t172 * t573 - t178 * t537 / 0.2e1 - t179 * t536 / 0.2e1 + t422 * (t539 - t540) / 0.2e1;
t502 = Icges(3,1) * t451 - t552;
t499 = Icges(7,1) * t449 - t546;
t498 = Icges(8,1) * t427 - t544;
t496 = -Icges(3,2) * t447 + t551;
t493 = -Icges(7,2) * t445 + t545;
t492 = -Icges(8,2) * t425 + t543;
t490 = Icges(3,5) * t451 - Icges(3,6) * t447;
t487 = Icges(7,5) * t449 - Icges(7,6) * t445;
t486 = Icges(8,5) * t427 - Icges(8,6) * t425;
t329 = -Icges(8,6) * t452 + t492 * t448;
t333 = -Icges(8,5) * t452 + t498 * t448;
t479 = t329 * t425 - t333 * t427;
t330 = Icges(8,6) * t448 + t492 * t452;
t334 = Icges(8,5) * t448 + t498 * t452;
t478 = -t330 * t425 + t334 * t427;
t366 = -Icges(5,2) * t422 - t548;
t367 = -Icges(5,1) * t421 - t547;
t470 = t366 * t421 - t367 * t422;
t377 = Icges(8,2) * t427 + t544;
t379 = Icges(8,1) * t425 + t543;
t469 = -t377 * t425 + t379 * t427;
t378 = -Icges(4,2) * t428 - t550;
t380 = -Icges(4,1) * t426 - t549;
t468 = t378 * t426 - t380 * t428;
t397 = Icges(3,2) * t451 + t552;
t399 = Icges(3,1) * t447 + t551;
t466 = -t397 * t447 + t399 * t451;
t462 = t528 * t452 + t529;
t325 = -Icges(8,3) * t452 + t486 * t448;
t326 = Icges(8,3) * t448 + t486 * t452;
t199 = t448 * (t441 * t326 + (t479 * t452 + (-t325 + t478) * t448) * t452);
t201 = t442 * t325 + (t478 * t448 + (-t326 + t479) * t452) * t448;
t461 = t186 + t199 + (-t187 - t201) * t452;
t460 = -t452 * t201 + t199 + t514;
t365 = -Icges(5,5) * t421 - Icges(5,6) * t422;
t458 = -t540 / 0.2e1 + t539 / 0.2e1 + (-t312 * t422 - t314 * t421 + t448 * t365 + t470 * t452 + t213) * t573 + (-t311 * t422 - t313 * t421 - t365 * t452 + t470 * t448 + t212) * t572;
t375 = Icges(8,5) * t425 + Icges(8,6) * t427;
t457 = t513 + (t330 * t427 + t334 * t425 + t448 * t375 + t469 * t452) * t573 + (t329 * t427 + t333 * t425 - t375 * t452 + t469 * t448) * t572;
t158 = m(6) * (t185 * t208 + t236 * t260 + t237 * t261) + m(5) * (t249 * t209 + (-t276 * t448 - t277 * t452) * t368) + t462;
t159 = m(6) * (t192 * t208 + t247 * t260 + t248 * t261) + m(5) * (t249 * t217 + (-t303 * t448 - t304 * t452) * t368) + t462;
t376 = -Icges(4,5) * t426 - Icges(4,6) * t428;
t453 = t458 + (-t332 * t428 - t336 * t426 + t448 * t376 + t468 * t452) * t573 + (-t331 * t428 - t335 * t426 - t376 * t452 + t468 * t448) * t572;
t403 = rSges(2,1) * t452 - t448 * rSges(2,2);
t402 = -t448 * rSges(2,1) - rSges(2,2) * t452;
t401 = rSges(3,1) * t447 + rSges(3,2) * t451;
t400 = rSges(7,1) * t445 + rSges(7,2) * t449;
t395 = Icges(3,5) * t447 + Icges(3,6) * t451;
t343 = Icges(3,3) * t448 + t490 * t452;
t342 = -Icges(3,3) * t452 + t490 * t448;
t341 = Icges(7,3) * t448 + t487 * t452;
t340 = -Icges(7,3) * t452 + t487 * t448;
t324 = (-pkin(6) - t560) * t452 + t522;
t323 = (pkin(12) - t561) * t452 + t521;
t322 = t553 + (pkin(6) - t506) * t448;
t321 = t554 + (-pkin(12) - t510) * t448;
t317 = t516 * t452;
t315 = t516 * t448;
t283 = t404 + t465;
t282 = t404 + t463;
t281 = t558 + (-t420 - t509) * t448;
t280 = t556 + (-t420 - t505) * t448;
t271 = t370 + t464;
t270 = t557 + (-t388 - t508) * t448;
t268 = t363 + t459;
t267 = t555 + (-t381 - t504) * t448;
t265 = t452 * (-t452 * t561 + t521) + (t510 * t448 - t554) * t448;
t264 = t452 * (-t452 * t560 + t522) + (t506 * t448 - t553) * t448;
t225 = t574 * t452 + t370 + t524;
t224 = (-t388 - t574) * t448 + t507;
t222 = t523 + t252;
t221 = t422 * t263 + t302 * t536;
t220 = -t262 * t422 - t302 * t537;
t219 = (t384 * t418 - t385 * t415) * t286 * pkin(7);
t216 = (-t421 * t538 + t527) * t422;
t215 = (-t262 * t452 + t263 * t448) * t421;
t191 = (-t371 * t418 + t372 * t415) * pkin(7) * t455;
t173 = m(9) * (-t267 * t452 - t268 * t448) * t356 + t513;
t162 = t216 + m(6) * (t220 * t224 + t221 * t225) + ((-t206 / 0.2e1 - t213 / 0.2e1) * t452 + (-t205 / 0.2e1 - t212 / 0.2e1) * t448) * t421;
t161 = m(6) * (t224 * t261 + t225 * t260) + m(5) * (-t270 * t452 - t271 * t448) * t368 + t458;
t160 = m(6) * (t208 ^ 2 + t260 ^ 2 + t261 ^ 2) + m(5) * (t520 * t368 ^ 2 + t249 ^ 2) + t462;
t156 = m(6) * (t208 * t215 + t220 * t261 + t221 * t260) + t503;
t155 = t453 + m(6) * (t224 * t248 + t225 * t247) + m(5) * (t270 * t304 + t271 * t303) + t219 * t161 + t218 * t173 + m(4) * (-t281 * t452 - t283 * t448) * t383;
t154 = t219 * t160 + t159;
t153 = t191 * t160 + t158;
t152 = m(6) * (t192 * t215 + t220 * t248 + t221 * t247) + t219 * t156 + t503;
t151 = m(6) * (t185 * t215 + t220 * t237 + t221 * t236) + t191 * t156 + t503;
t150 = t453 + t457 + m(3) * (-t321 * t452 - t323 * t448) * t401 + t191 * t161 + t190 * t173 + t228 * (m(7) * (-t322 * t452 - t324 * t448) * t400 + (t442 / 0.2e1 + t441 / 0.2e1) * (Icges(7,5) * t445 + Icges(7,6) * t449)) - t229 * (m(9) * (t267 * t279 + t268 * t278) + m(8) * (-t280 * t452 - t282 * t448) * t382 + t457) + m(8) * (t280 * t317 + t282 * t315) + m(4) * (t281 * t318 + t283 * t316) + m(5) * (t270 * t277 + t271 * t276) + m(9) * (t267 * t273 + t268 * t272) + m(6) * (t224 * t237 + t225 * t236) + ((Icges(3,6) * t448 + t496 * t452) * t451 + (Icges(3,5) * t448 + t502 * t452) * t447 + t448 * t395 + t452 * t466 + t228 * (t449 * (Icges(7,6) * t448 + t493 * t452) + t445 * (Icges(7,5) * t448 + t499 * t452))) * t573 + ((-Icges(3,6) * t452 + t496 * t448) * t451 + (-Icges(3,5) * t452 + t502 * t448) * t447 + t228 * (t449 * (-Icges(7,6) * t452 + t493 * t448) + t445 * (-Icges(7,5) * t452 + t499 * t448)) - t395 * t452 + t448 * t466) * t572;
t1 = [t449 * (Icges(7,2) * t449 + t546) + t445 * (Icges(7,1) * t445 + t545) + t451 * t397 + t447 * t399 + t425 * t379 - t426 * t380 + t427 * t377 - t428 * t378 - t422 * t366 - t410 * t354 + t409 * t355 + (-t367 - t538) * t421 + m(4) * (t281 ^ 2 + t283 ^ 2) + m(3) * (t321 ^ 2 + t323 ^ 2) + m(2) * (t402 ^ 2 + t403 ^ 2) + m(9) * (t267 ^ 2 + t268 ^ 2) + m(8) * (t280 ^ 2 + t282 ^ 2) + m(7) * (t322 ^ 2 + t324 ^ 2) + m(6) * (t224 ^ 2 + t225 ^ 2) + m(5) * (t270 ^ 2 + t271 ^ 2) + t527 + Icges(2,3), t150, t155, t162; t150, t460 + t519 + t448 * t441 * t343 + t580 * (t448 * (-t340 * t576 + t441 * t341) + m(7) * (t520 * t400 ^ 2 + t264 ^ 2)) + m(8) * (t222 ^ 2 + t315 ^ 2 + t317 ^ 2) + m(4) * (t223 ^ 2 + t316 ^ 2 + t318 ^ 2) + m(9) * (t207 ^ 2 + t272 ^ 2 + t273 ^ 2) + m(5) * (t209 ^ 2 + t276 ^ 2 + t277 ^ 2) + m(6) * (t185 ^ 2 + t236 ^ 2 + t237 ^ 2) + m(3) * (t520 * t401 ^ 2 + t265 ^ 2) + (-t460 - 0.2e1 * m(9) * (t207 * t214 + t272 * t278 + t273 * t279) - 0.2e1 * m(8) * (t252 * t222 + (-t315 * t448 - t317 * t452) * t382) - t461 + t229 * (m(9) * (t214 ^ 2 + t278 ^ 2 + t279 ^ 2) + m(8) * (t520 * t382 ^ 2 + t252 ^ 2) + t461)) * t229 + (t158 + t153) * t191 + (t157 + t578) * t190 + (-t442 * t342 + t580 * (-t442 * t340 + t341 * t576) + (-t448 * t342 + t452 * t343) * t448 + t581) * t452, t153 * t219 + t191 * t159 + t579, t151; t155, t154 * t191 + t219 * t158 + t579, t218 ^ 2 * t168 + (t159 + t154) * t219 + m(6) * (t192 ^ 2 + t247 ^ 2 + t248 ^ 2) + m(5) * (t217 ^ 2 + t303 ^ 2 + t304 ^ 2) + m(4) * (t520 * t383 ^ 2 + t253 ^ 2) + t456, t152; t162, t151, t152, t422 * t216 + m(6) * (t215 ^ 2 + t220 ^ 2 + t221 ^ 2) + (t422 * (-t205 * t448 - t206 * t452) - t452 * t172 - t448 * t171) * t421;];
Mq = t1;
