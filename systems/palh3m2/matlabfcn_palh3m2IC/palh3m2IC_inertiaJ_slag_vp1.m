% Calculate joint inertia matrix for
% palh3m2IC
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
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_inertiaJ_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_inertiaJ_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2IC_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2IC_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:59:45
% EndTime: 2020-05-07 04:59:53
% DurationCPUTime: 6.67s
% Computational Cost: add. (23443->579), mult. (20395->822), div. (109->12), fcn. (20764->39), ass. (0->336)
t443 = sin(qJ(1));
t447 = cos(qJ(1));
t601 = t443 * t447;
t436 = t443 ^ 2;
t437 = t447 ^ 2;
t538 = t436 + t437;
t438 = qJ(2) + qJ(7);
t423 = sin(t438);
t425 = cos(t438);
t508 = rSges(8,1) * t425 - rSges(8,2) * t423;
t537 = pkin(15) - qJ(7);
t529 = (-qJ(8) + t537);
t420 = -qJ(2) + t529;
t412 = sin(t420);
t413 = cos(t420);
t507 = -rSges(9,1) * t413 - rSges(9,2) * t412;
t459 = t443 * rSges(9,3) + t507 * t447;
t578 = rSges(9,3) * t447;
t242 = t443 * (t507 * t443 - t578) + t447 * t459;
t446 = cos(qJ(2));
t417 = t446 * pkin(1) + pkin(12);
t427 = -qJ(2) + t537;
t381 = pkin(3) * cos(t427) + t417;
t363 = t447 * t381;
t403 = t447 * t417;
t217 = t436 * (t381 - t417) + t447 * (-t403 + t363) + t242;
t542 = t436 * (-pkin(12) + t417) + t447 * (-pkin(12) * t447 + t403);
t210 = t217 + t542;
t591 = pkin(3) * sin(t427);
t387 = t443 * t591;
t356 = rSges(9,1) * t412 - rSges(9,2) * t413;
t442 = sin(qJ(2));
t592 = pkin(1) * t442;
t527 = -t356 - t592;
t270 = t527 * t443 + t387;
t389 = t447 * t591;
t271 = t527 * t447 + t389;
t488 = -Icges(9,5) * t413 - Icges(9,6) * t412;
t291 = -Icges(9,3) * t447 + t488 * t443;
t292 = Icges(9,3) * t443 + t488 * t447;
t564 = Icges(9,4) * t413;
t494 = -Icges(9,2) * t412 - t564;
t294 = Icges(9,6) * t443 + t494 * t447;
t565 = Icges(9,4) * t412;
t500 = -Icges(9,1) * t413 - t565;
t296 = Icges(9,5) * t443 + t500 * t447;
t486 = -t294 * t412 - t296 * t413;
t293 = -Icges(9,6) * t447 + t494 * t443;
t295 = -Icges(9,5) * t447 + t500 * t443;
t487 = t293 * t412 + t295 * t413;
t191 = t443 * (t436 * t292 + (t487 * t447 + (-t291 + t486) * t443) * t447);
t192 = t437 * t291 + (t486 * t443 + (-t292 + t487) * t447) * t443;
t523 = -t447 * t192 + t191;
t168 = m(9) * (t210 * t242 + (-t270 * t443 - t271 * t447) * t356) + t523;
t278 = -t356 * t443 + t387;
t279 = -t356 * t447 + t389;
t169 = m(9) * (t217 * t242 + (-t278 * t443 - t279 * t447) * t356) + t523;
t452 = 0.1e1 / pkin(2);
t536 = qJ(7) + pkin(16);
t549 = qJ(2) - qJ(6);
t518 = t536 + t549;
t411 = sin(t518);
t410 = 0.1e1 / t411;
t557 = (-pkin(2) * t411 - pkin(1) * sin(t549)) * t410;
t531 = t452 * t557;
t602 = t169 * t531 + t168;
t535 = qJ(4) + pkin(14);
t528 = (qJ(3) + t535);
t516 = pkin(15) + t528;
t464 = t516 - t549;
t455 = -(2 * qJ(7)) - pkin(16) + t464;
t465 = t516 + t549;
t460 = pkin(16) + t465;
t599 = -cos(qJ(8) - t455) + cos(qJ(8) - t460);
t439 = qJ(2) + qJ(3);
t428 = qJ(4) + t439;
t418 = sin(t428);
t419 = cos(t428);
t590 = pkin(8) * t419;
t598 = -(pkin(10) + rSges(6,3)) * t418 - t590;
t595 = t443 / 0.2e1;
t594 = -t447 / 0.2e1;
t589 = rSges(3,1) * t446;
t426 = cos(t439);
t588 = rSges(4,1) * t426;
t587 = rSges(5,1) * t419;
t444 = cos(qJ(6));
t586 = rSges(7,1) * t444;
t584 = rSges(3,2) * t442;
t440 = sin(qJ(6));
t583 = rSges(7,2) * t440;
t581 = rSges(4,3) * t447;
t580 = rSges(5,3) * t447;
t579 = rSges(8,3) * t447;
t577 = t447 * rSges(3,3);
t576 = t447 * rSges(7,3);
t575 = Icges(3,4) * t442;
t574 = Icges(3,4) * t446;
t424 = sin(t439);
t573 = Icges(4,4) * t424;
t572 = Icges(4,4) * t426;
t571 = Icges(5,4) * t418;
t570 = Icges(5,4) * t419;
t569 = Icges(7,4) * t440;
t568 = Icges(7,4) * t444;
t567 = Icges(8,4) * t423;
t566 = Icges(8,4) * t425;
t173 = m(9) * (t538 * t356 ^ 2 + t242 ^ 2) + t523;
t563 = t173 / pkin(7) ^ 2;
t445 = cos(qJ(5));
t550 = t445 * t447;
t441 = sin(qJ(5));
t553 = t441 * t443;
t357 = t419 * t553 - t550;
t551 = t443 * t445;
t552 = t441 * t447;
t358 = -t419 * t551 - t552;
t556 = t418 * t443;
t254 = Icges(6,5) * t358 + Icges(6,6) * t357 - Icges(6,3) * t556;
t256 = Icges(6,4) * t358 + Icges(6,2) * t357 - Icges(6,6) * t556;
t258 = Icges(6,1) * t358 + Icges(6,4) * t357 - Icges(6,5) * t556;
t208 = t254 * t419 + (t256 * t441 - t258 * t445) * t418;
t562 = t208 * t447;
t359 = t419 * t552 + t551;
t360 = -t419 * t550 + t553;
t555 = t418 * t447;
t255 = Icges(6,5) * t360 + Icges(6,6) * t359 - Icges(6,3) * t555;
t257 = Icges(6,4) * t360 + Icges(6,2) * t359 - Icges(6,6) * t555;
t259 = Icges(6,1) * t360 + Icges(6,4) * t359 - Icges(6,5) * t555;
t209 = t255 * t419 + (t257 * t441 - t259 * t445) * t418;
t561 = t209 * t443;
t541 = cos((2 * t528)) - cos((2 * t529));
t560 = (t541 * pkin(9) + (cos(0.2e1 * qJ(3) + t535) - cos(-(2 * qJ(7)) + (2 * pkin(15)) - 0.2e1 * qJ(8) - t535)) * pkin(4)) / t541;
t290 = Icges(6,5) * t419 + (-Icges(6,1) * t445 + Icges(6,4) * t441) * t418;
t559 = t290 * t445;
t304 = 0.1e1 / t599;
t558 = t304 * t452;
t554 = t424 * t447;
t200 = -t254 * t555 + t256 * t359 + t258 * t360;
t201 = -t255 * t555 + t257 * t359 + t259 * t360;
t184 = -t200 * t447 + t201 * t443;
t491 = -Icges(5,5) * t419 + Icges(5,6) * t418;
t309 = -Icges(5,3) * t447 + t491 * t443;
t310 = Icges(5,3) * t443 + t491 * t447;
t497 = Icges(5,2) * t418 - t570;
t312 = Icges(5,6) * t443 + t497 * t447;
t503 = -Icges(5,1) * t419 + t571;
t314 = Icges(5,5) * t443 + t503 * t447;
t484 = t312 * t418 - t314 * t419;
t311 = -Icges(5,6) * t447 + t497 * t443;
t313 = -Icges(5,5) * t447 + t503 * t443;
t485 = -t311 * t418 + t313 * t419;
t548 = (t436 * t310 + t184 + (t485 * t447 + (-t309 + t484) * t443) * t447) * t443;
t198 = -t254 * t556 + t256 * t357 + t258 * t358;
t199 = -t255 * t556 + t257 * t357 + t259 * t358;
t183 = -t198 * t447 + t199 * t443;
t197 = t437 * t309 + (t484 * t443 + (-t310 + t485) * t447) * t443;
t547 = -t197 - t183;
t288 = Icges(6,3) * t419 + (-Icges(6,5) * t445 + Icges(6,6) * t441) * t418;
t289 = Icges(6,6) * t419 + (-Icges(6,4) * t445 + Icges(6,2) * t441) * t418;
t546 = t418 * t441 * t289 + t419 * t288;
t469 = rSges(5,2) * t555 + t443 * rSges(5,3) - t447 * t587;
t511 = rSges(5,2) * t418 - t587;
t249 = t443 * (t511 * t443 - t580) + t447 * t469;
t385 = -pkin(4) * t426 + t417;
t371 = t447 * t385;
t545 = t436 * (t385 - t417) + t447 * (-t403 + t371);
t301 = rSges(6,3) * t419 + (-rSges(6,1) * t445 + rSges(6,2) * t441) * t418;
t544 = pkin(8) * t418 - pkin(10) * t419 - t301;
t468 = t443 * rSges(8,3) + t508 * t447;
t252 = t443 * (t508 * t443 - t579) + t447 * t468;
t470 = rSges(4,2) * t554 + t443 * rSges(4,3) - t447 * t588;
t512 = rSges(4,2) * t424 - t588;
t253 = t443 * (t512 * t443 - t581) + t447 * t470;
t543 = t360 * rSges(6,1) + t359 * rSges(6,2);
t540 = t443 * rSges(7,3) + t447 * t586;
t539 = t443 * rSges(3,3) + t447 * t589;
t392 = -sin(qJ(7) + qJ(8) - t516);
t421 = sin(t535);
t534 = pkin(4) / t392 * t421;
t461 = -qJ(7) + t465;
t462 = -qJ(7) + t464;
t220 = ((cos(t462) - cos(t461)) * pkin(1) + (cos(t455) - cos(t460)) * pkin(2)) * pkin(3) + ((-cos(qJ(8) - t462) + cos(qJ(8) - t461)) * pkin(1) + t599 * pkin(2)) * pkin(7);
t533 = t220 * t558;
t449 = 0.1e1 / pkin(9);
t532 = t449 * t560;
t492 = -Icges(4,5) * t426 + Icges(4,6) * t424;
t327 = -Icges(4,3) * t447 + t492 * t443;
t328 = Icges(4,3) * t443 + t492 * t447;
t498 = Icges(4,2) * t424 - t572;
t332 = Icges(4,6) * t443 + t498 * t447;
t504 = -Icges(4,1) * t426 + t573;
t336 = Icges(4,5) * t443 + t504 * t447;
t480 = t332 * t424 - t336 * t426;
t331 = -Icges(4,6) * t447 + t498 * t443;
t335 = -Icges(4,5) * t447 + t504 * t443;
t481 = -t331 * t424 + t335 * t426;
t530 = t443 * (t436 * t328 + (t481 * t447 + (-t327 + t480) * t443) * t447) + t548;
t368 = -rSges(5,1) * t418 - rSges(5,2) * t419;
t526 = -t368 - t592;
t382 = rSges(8,1) * t423 + rSges(8,2) * t425;
t525 = -t382 - t592;
t383 = -rSges(4,1) * t424 - rSges(4,2) * t426;
t524 = -t383 - t592;
t260 = t544 * t443;
t261 = t544 * t447;
t522 = pkin(3) * (pkin(1) * (cos(qJ(8) - t549) - cos(qJ(8) + t549)) + (cos(qJ(8) - t518) - cos(qJ(8) + t518)) * pkin(2)) * t558;
t450 = 0.1e1 / pkin(7);
t521 = t450 * t534;
t520 = t450 * t533;
t352 = Icges(9,5) * t412 - Icges(9,6) * t413;
t353 = -Icges(9,2) * t413 + t565;
t354 = Icges(9,1) * t412 - t564;
t475 = -t353 * t412 - t354 * t413;
t519 = (-t294 * t413 + t296 * t412 + t352 * t443 + t475 * t447) * t595 + (-t293 * t413 + t295 * t412 - t352 * t447 + t475 * t443) * t594;
t221 = t545 + t249;
t510 = -rSges(6,1) * t358 - rSges(6,2) * t357;
t262 = -rSges(6,3) * t556 - t510;
t263 = -rSges(6,3) * t555 + t543;
t211 = t443 * t262 + t447 * t263 + t538 * (-pkin(10) * t418 - t590);
t517 = t544 - t592;
t514 = t449 * t522;
t513 = -t584 + t589;
t509 = -t583 + t586;
t215 = -t288 * t556 + t289 * t357 + t290 * t358;
t176 = t215 * t419 + (-t198 * t443 - t199 * t447) * t418;
t216 = -t288 * t555 + t289 * t359 + t290 * t360;
t177 = t216 * t419 + (-t200 * t443 - t201 * t447) * t418;
t506 = t176 * t594 + t177 * t595 - t183 * t556 / 0.2e1 - t184 * t555 / 0.2e1 + t419 * (t561 - t562) / 0.2e1;
t505 = Icges(3,1) * t446 - t575;
t502 = Icges(7,1) * t444 - t569;
t501 = Icges(8,1) * t425 - t567;
t499 = -Icges(3,2) * t442 + t574;
t496 = -Icges(7,2) * t440 + t568;
t495 = -Icges(8,2) * t423 + t566;
t493 = Icges(3,5) * t446 - Icges(3,6) * t442;
t490 = Icges(7,5) * t444 - Icges(7,6) * t440;
t489 = Icges(8,5) * t425 - Icges(8,6) * t423;
t329 = -Icges(8,6) * t447 + t495 * t443;
t333 = -Icges(8,5) * t447 + t501 * t443;
t483 = t329 * t423 - t333 * t425;
t330 = Icges(8,6) * t443 + t495 * t447;
t334 = Icges(8,5) * t443 + t501 * t447;
t482 = -t330 * t423 + t334 * t425;
t366 = -Icges(5,2) * t419 - t571;
t367 = -Icges(5,1) * t418 - t570;
t474 = t366 * t418 - t367 * t419;
t377 = Icges(8,2) * t425 + t567;
t379 = Icges(8,1) * t423 + t566;
t473 = -t377 * t423 + t379 * t425;
t378 = -Icges(4,2) * t426 - t573;
t380 = -Icges(4,1) * t424 - t572;
t472 = t378 * t424 - t380 * t426;
t396 = Icges(3,2) * t446 + t575;
t398 = Icges(3,1) * t442 + t574;
t471 = -t396 * t442 + t398 * t446;
t195 = t211 + t545;
t467 = t547 * t447 + t548;
t325 = -Icges(8,3) * t447 + t489 * t443;
t326 = Icges(8,3) * t443 + t489 * t447;
t202 = t443 * (t436 * t326 + (t483 * t447 + (-t325 + t482) * t443) * t447);
t204 = t437 * t325 + (t482 * t443 + (-t326 + t483) * t447) * t443;
t466 = t191 + t202 + (-t204 - t192) * t447;
t463 = -t447 * t204 + t202 + t523;
t365 = -Icges(5,5) * t418 - Icges(5,6) * t419;
t458 = -t562 / 0.2e1 + t561 / 0.2e1 + (-t312 * t419 - t314 * t418 + t365 * t443 + t474 * t447 + t216) * t595 + (-t311 * t419 - t313 * t418 - t365 * t447 + t474 * t443 + t215) * t594;
t375 = Icges(8,5) * t423 + Icges(8,6) * t425;
t457 = t519 + (t330 * t425 + t334 * t423 + t375 * t443 + t473 * t447) * t595 + (t329 * t425 + t333 * t423 - t375 * t447 + t473 * t443) * t594;
t205 = t437 * t327 + (t480 * t443 + (-t328 + t481) * t447) * t443;
t456 = (-t205 + t547) * t447 + t530;
t190 = t195 + t542;
t212 = t221 + t542;
t407 = t443 * pkin(4) * t424;
t236 = t517 * t443 + t407;
t408 = pkin(4) * t554;
t237 = t517 * t447 + t408;
t275 = t526 * t443 + t407;
t276 = t526 * t447 + t408;
t163 = m(6) * (t190 * t211 + t236 * t260 + t237 * t261) + m(5) * (t212 * t249 + (-t275 * t443 - t276 * t447) * t368) + t467;
t247 = t407 + t260;
t248 = t408 + t261;
t302 = -t368 * t443 + t407;
t303 = -t368 * t447 + t408;
t164 = m(6) * (t195 * t211 + t247 * t260 + t248 * t261) + m(5) * (t221 * t249 + (-t302 * t443 - t303 * t447) * t368) + t467;
t225 = t542 + t253;
t316 = t524 * t443;
t318 = t524 * t447;
t454 = m(6) * (t190 * t195 + t236 * t247 + t237 * t248) + m(5) * (t212 * t221 + t275 * t302 + t276 * t303) + m(4) * (t225 * t253 + (-t316 * t443 - t318 * t447) * t383) + t456;
t376 = -Icges(4,5) * t424 - Icges(4,6) * t426;
t453 = t458 + (-t332 * t426 - t336 * t424 + t376 * t443 + t472 * t447) * t595 + (-t331 * t426 - t335 * t424 - t376 * t447 + t472 * t443) * t594;
t422 = sin(t536);
t402 = rSges(2,1) * t447 - rSges(2,2) * t443;
t401 = -rSges(2,1) * t443 - rSges(2,2) * t447;
t400 = rSges(3,1) * t442 + rSges(3,2) * t446;
t399 = rSges(7,1) * t440 + rSges(7,2) * t444;
t394 = Icges(3,5) * t442 + Icges(3,6) * t446;
t342 = Icges(3,3) * t443 + t493 * t447;
t341 = -Icges(3,3) * t447 + t493 * t443;
t340 = Icges(7,3) * t443 + t490 * t447;
t339 = -Icges(7,3) * t447 + t490 * t443;
t324 = (-pkin(6) - t583) * t447 + t540;
t323 = (pkin(12) - t584) * t447 + t539;
t322 = t576 + (pkin(6) - t509) * t443;
t321 = t577 + (-pkin(12) - t513) * t443;
t317 = t525 * t447;
t315 = t525 * t443;
t283 = t403 + t470;
t282 = t403 + t468;
t281 = t581 + (-t417 - t512) * t443;
t280 = t579 + (-t417 - t508) * t443;
t269 = t371 + t469;
t268 = t580 + (-t385 - t511) * t443;
t267 = t363 + t459;
t266 = t578 + (-t381 - t507) * t443;
t265 = t447 * (-t447 * t584 + t539) + (t513 * t443 - t577) * t443;
t264 = t447 * (-t447 * t583 + t540) + (t509 * t443 - t576) * t443;
t227 = t447 * t598 + t371 + t543;
t226 = (-t385 - t598) * t443 + t510;
t224 = t542 + t252;
t223 = t263 * t419 + t301 * t555;
t222 = -t262 * t419 - t301 * t556;
t219 = (-t418 * t559 + t546) * t419;
t218 = (-t262 * t447 + t263 * t443) * t418;
t178 = m(9) * (-t266 * t447 - t267 * t443) * t356 + t519;
t167 = t219 + m(6) * (t222 * t226 + t223 * t227) + ((-t209 / 0.2e1 - t216 / 0.2e1) * t447 + (-t215 / 0.2e1 - t208 / 0.2e1) * t443) * t418;
t166 = m(6) * (t226 * t261 + t227 * t260) + m(5) * (-t268 * t447 - t269 * t443) * t368 + t458;
t165 = m(5) * (t538 * t368 ^ 2 + t249 ^ 2) + m(6) * (t211 ^ 2 + t260 ^ 2 + t261 ^ 2) + t467;
t162 = (t173 * t220 * t304 * t450 + t169 * t557) * t452 + t168;
t161 = m(6) * (t211 * t218 + t222 * t261 + t223 * t260) + t506;
t160 = t453 + t178 * t521 + m(6) * (t226 * t248 + t227 * t247) + m(5) * (t268 * t303 + t269 * t302) + m(4) * (-t281 * t447 - t283 * t443) * t383 - t166 * t532;
t159 = -t165 * t532 + t164;
t158 = t165 * t514 + t163;
t157 = m(6) * (t195 * t218 + t222 * t248 + t223 * t247) - t161 * t532 + t506;
t156 = m(6) * (t190 * t218 + t222 * t237 + t223 * t236) + t161 * t514 + t506;
t155 = m(6) * (t226 * t237 + t227 * t236) + m(9) * (t266 * t271 + t267 * t270) + m(5) * (t268 * t276 + t269 * t275) + m(8) * (t280 * t317 + t282 * t315) + m(4) * (t281 * t318 + t283 * t316) + t453 + m(3) * (-t321 * t447 - t323 * t443) * t400 + t166 * t514 + (m(9) * (t266 * t279 + t267 * t278) + m(8) * (-t280 * t447 - t282 * t443) * t382 + t457) * t531 + pkin(1) * t422 / pkin(5) * t410 * ((t444 * (Icges(7,6) * t443 + t496 * t447) + t440 * (Icges(7,5) * t443 + t502 * t447)) * t595 + (t444 * (-Icges(7,6) * t447 + t496 * t443) + t440 * (-Icges(7,5) * t447 + t502 * t443)) * t594 + m(7) * (-t322 * t447 - t324 * t443) * t399 + (t436 / 0.2e1 + t437 / 0.2e1) * (Icges(7,5) * t440 + Icges(7,6) * t444)) + t457 + t178 * t520 + (t394 * t443 + t447 * t471) * t595 + (-t394 * t447 + t443 * t471) * t594 + ((-Icges(3,6) * t447 + t499 * t443) * t446 + (-Icges(3,5) * t447 + t505 * t443) * t442) * t594 + ((Icges(3,6) * t443 + t499 * t447) * t446 + (Icges(3,5) * t443 + t505 * t447) * t442) * t595;
t1 = [t444 * (Icges(7,2) * t444 + t569) + t440 * (Icges(7,1) * t440 + t568) + m(9) * (t266 ^ 2 + t267 ^ 2) + m(8) * (t280 ^ 2 + t282 ^ 2) + m(7) * (t322 ^ 2 + t324 ^ 2) + m(6) * (t226 ^ 2 + t227 ^ 2) + m(5) * (t268 ^ 2 + t269 ^ 2) + m(4) * (t281 ^ 2 + t283 ^ 2) + m(3) * (t321 ^ 2 + t323 ^ 2) + m(2) * (t401 ^ 2 + t402 ^ 2) + (-t367 - t559) * t418 + t446 * t396 + t442 * t398 - t426 * t378 - t419 * t366 + t423 * t379 - t424 * t380 + t425 * t377 + t412 * t354 - t413 * t353 + Icges(2,3) + t546, t155, t160, t167; t155, pkin(1) ^ 2 * t422 ^ 2 / pkin(5) ^ 2 / t411 ^ 2 * (t443 * (-t339 * t601 + t436 * t340) - t447 * (t437 * t339 - t340 * t601) + m(7) * (t538 * t399 ^ 2 + t264 ^ 2)) + m(5) * (t212 ^ 2 + t275 ^ 2 + t276 ^ 2) + m(8) * (t224 ^ 2 + t315 ^ 2 + t317 ^ 2) + m(4) * (t225 ^ 2 + t316 ^ 2 + t318 ^ 2) + m(6) * (t190 ^ 2 + t236 ^ 2 + t237 ^ 2) + m(9) * (t210 ^ 2 + t270 ^ 2 + t271 ^ 2) + m(3) * (t538 * t400 ^ 2 + t265 ^ 2) + t463 - t447 * (t437 * t341 - t342 * t601) + t443 * (-t341 * t601 + t436 * t342) - t447 * t197 - t447 * t183 - t447 * t205 + t530 + (t463 + 0.2e1 * m(9) * (t210 * t217 + t270 * t278 + t271 * t279) + 0.2e1 * m(8) * (t224 * t252 + (-t315 * t443 - t317 * t447) * t382) + t466 + (m(8) * (t538 * t382 ^ 2 + t252 ^ 2) + m(9) * (t217 ^ 2 + t278 ^ 2 + t279 ^ 2) + t466) * t531) * t531 + (t162 + t602) * t520 + (t158 + t163) * t514, t162 * t521 + (-t158 * t560 + t164 * t522) * t449 + t454, t156; t160, (t159 * t522 - t163 * t560) * t449 + (t602 * t450 + t533 * t563) * t534 + t454, pkin(4) ^ 2 * t421 ^ 2 / t392 ^ 2 * t563 - (t159 + t164) * t532 + m(6) * (t195 ^ 2 + t247 ^ 2 + t248 ^ 2) + m(5) * (t221 ^ 2 + t302 ^ 2 + t303 ^ 2) + m(4) * (t538 * t383 ^ 2 + t253 ^ 2) + t456, t157; t167, t156, t157, t419 * t219 + m(6) * (t218 ^ 2 + t222 ^ 2 + t223 ^ 2) + (-t447 * t177 - t443 * t176 + t419 * (-t208 * t443 - t209 * t447)) * t418;];
Mq = t1;
