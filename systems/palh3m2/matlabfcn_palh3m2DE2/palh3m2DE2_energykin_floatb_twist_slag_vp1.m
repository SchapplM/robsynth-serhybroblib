% Calculate kinetic energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:04
% EndTime: 2020-05-07 02:14:15
% DurationCPUTime: 9.99s
% Computational Cost: add. (12054->596), mult. (17350->789), div. (0->0), fcn. (25861->87), ass. (0->351)
t617 = Icges(3,3) + Icges(7,3);
t454 = sin(pkin(16));
t455 = cos(pkin(16));
t462 = sin(pkin(15));
t468 = cos(pkin(15));
t348 = t454 * t468 + t455 * t462;
t349 = -t454 * t462 + t455 * t468;
t459 = sin(qJ(3));
t465 = cos(qJ(3));
t276 = t348 * t465 + t349 * t459;
t277 = -t348 * t459 + t349 * t465;
t460 = sin(qJ(2));
t466 = cos(qJ(2));
t252 = -t276 * t460 + t277 * t466;
t448 = pkin(17) + pkin(18);
t421 = sin(t448);
t422 = cos(t448);
t449 = qJ(2) + qJ(3);
t497 = t276 * t466 + t460 * t277;
t219 = atan2(-t252 * t421 - t497 * t422, t252 * t422 - t421 * t497) + t449;
t217 = sin(t219);
t218 = cos(t219);
t456 = sin(pkin(18));
t457 = cos(pkin(18));
t350 = t456 * t468 + t457 * t462;
t351 = -t456 * t462 + t457 * t468;
t255 = qJ(2) + atan2(t350 * t466 + t351 * t460, t350 * t460 - t351 * t466);
t253 = sin(t255);
t254 = cos(t255);
t616 = Icges(5,5) * t218 - Icges(8,5) * t254 - Icges(5,6) * t217 + Icges(8,6) * t253;
t463 = sin(pkin(14));
t469 = cos(pkin(14));
t355 = t462 * t469 - t463 * t468;
t359 = t462 * t463 + t468 * t469;
t278 = t355 * t466 + t359 * t460;
t531 = -t355 * t460 + t359 * t466;
t615 = Icges(3,5) * t466 + Icges(7,5) * t531 - Icges(3,6) * t460 - Icges(7,6) * t278;
t614 = Icges(2,2) + Icges(5,3) + Icges(8,3);
t613 = Icges(7,4) * t531;
t436 = qJ(2) + pkin(15) + pkin(18);
t405 = pkin(17) + qJ(3) + t436;
t400 = pkin(16) + t405;
t312 = atan2(-sin(t400), cos(t400)) + t449;
t305 = qJ(4) + t312;
t298 = qJ(1) + t305;
t286 = cos(t298);
t306 = -qJ(4) + t312;
t299 = qJ(1) - t306;
t287 = cos(t299);
t424 = V_base(6) + qJD(1);
t408 = qJD(4) + t424;
t450 = qJ(1) + qJ(4);
t432 = cos(t450);
t590 = t432 / 0.2e1;
t612 = (t590 + t286 / 0.4e1 + t287 / 0.4e1) * t408;
t300 = qJ(1) + t306;
t284 = sin(t300);
t301 = qJ(1) - t305;
t285 = sin(t301);
t409 = -qJD(4) + t424;
t451 = qJ(1) - qJ(4);
t428 = sin(t451);
t591 = t428 / 0.2e1;
t611 = (t591 - t284 / 0.4e1 - t285 / 0.4e1) * t409;
t282 = sin(t298);
t283 = sin(t299);
t610 = t283 + t282;
t609 = t285 + t284;
t608 = t287 + t286;
t288 = cos(t300);
t289 = cos(t301);
t607 = t288 + t289;
t461 = sin(qJ(1));
t467 = cos(qJ(1));
t507 = -Icges(7,2) * t278 + t613;
t248 = Icges(7,6) * t461 + t507 * t467;
t249 = -Icges(7,6) * t467 + t507 * t461;
t573 = Icges(7,4) * t278;
t513 = Icges(7,1) * t531 - t573;
t250 = Icges(7,5) * t461 + t513 * t467;
t251 = -Icges(7,5) * t467 + t513 * t461;
t257 = Icges(7,2) * t531 + t573;
t258 = Icges(7,1) * t278 + t613;
t578 = Icges(3,4) * t466;
t510 = -Icges(3,2) * t460 + t578;
t325 = -Icges(3,6) * t467 + t510 * t461;
t326 = Icges(3,6) * t461 + t510 * t467;
t579 = Icges(3,4) * t460;
t516 = Icges(3,1) * t466 - t579;
t327 = -Icges(3,5) * t467 + t516 * t461;
t328 = Icges(3,5) * t461 + t516 * t467;
t379 = Icges(3,2) * t466 + t579;
t382 = Icges(3,1) * t460 + t578;
t401 = -qJD(2) * t467 + V_base(5);
t402 = qJD(2) * t461 + V_base(4);
t606 = (-t257 * t278 + t258 * t531 - t379 * t460 + t382 * t466) * t424 + (-t248 * t278 + t250 * t531 - t326 * t460 + t328 * t466) * t402 + (-t249 * t278 + t251 * t531 - t325 * t460 + t327 * t466) * t401;
t605 = (Icges(3,5) * t460 + Icges(7,5) * t278 + Icges(3,6) * t466 + Icges(7,6) * t531) * t424 + (t617 * t461 + t615 * t467) * t402 + (t615 * t461 - t617 * t467) * t401;
t447 = qJD(2) + qJD(3);
t370 = -t447 * t467 + V_base(5);
t521 = rSges(3,1) * t466 - rSges(3,2) * t460;
t491 = pkin(12) + t521;
t604 = t461 * rSges(3,3) + t491 * t467;
t437 = Icges(2,4) * t467;
t580 = Icges(2,4) * t461;
t603 = (Icges(5,5) * t217 - Icges(8,5) * t253 + Icges(5,6) * t218 - Icges(8,6) * t254) * t424 + (t616 * t461 + t614 * t467 + t580) * V_base(5) + (-t614 * t461 + t616 * t467 + t437) * V_base(4);
t574 = Icges(5,4) * t218;
t508 = Icges(5,2) * t217 - t574;
t191 = -Icges(5,6) * t467 + t508 * t461;
t192 = Icges(5,6) * t461 + t508 * t467;
t575 = Icges(5,4) * t217;
t514 = -Icges(5,1) * t218 + t575;
t193 = -Icges(5,5) * t467 + t514 * t461;
t194 = Icges(5,5) * t461 + t514 * t467;
t196 = -Icges(5,2) * t218 - t575;
t197 = -Icges(5,1) * t217 - t574;
t571 = Icges(8,4) * t254;
t506 = -Icges(8,2) * t253 + t571;
t233 = -Icges(8,6) * t467 + t506 * t461;
t234 = Icges(8,6) * t461 + t506 * t467;
t572 = Icges(8,4) * t253;
t512 = Icges(8,1) * t254 - t572;
t235 = -Icges(8,5) * t467 + t512 * t461;
t236 = Icges(8,5) * t461 + t512 * t467;
t238 = Icges(8,2) * t254 + t572;
t239 = Icges(8,1) * t253 + t571;
t602 = (t196 * t217 - t197 * t218 - t238 * t253 + t239 * t254) * t424 + (Icges(2,1) * t461 + t191 * t217 - t193 * t218 - t233 * t253 + t235 * t254 + t437) * V_base(5) + (Icges(2,1) * t467 + t192 * t217 - t194 * t218 - t234 * t253 + t236 * t254 - t580) * V_base(4);
t307 = qJ(1) + t312;
t292 = sin(t307);
t601 = t292 / 0.2e1;
t308 = qJ(1) - t312;
t293 = sin(t308);
t600 = -t293 / 0.2e1;
t296 = cos(t307);
t599 = -t296 / 0.2e1;
t440 = t461 * rSges(6,2);
t589 = t440 / 0.2e1;
t588 = pkin(1) * t460;
t587 = pkin(1) * t466;
t431 = cos(t449);
t586 = pkin(4) * t431;
t585 = pkin(1) * qJD(2);
t584 = t461 * rSges(6,1);
t582 = t467 * rSges(6,1);
t581 = t467 * rSges(3,3);
t426 = sin(t449);
t577 = Icges(4,4) * t426;
t576 = Icges(4,4) * t431;
t356 = t459 * t468 + t462 * t465;
t357 = -t459 * t462 + t465 * t468;
t280 = t356 * t460 - t357 * t466;
t494 = t356 * t466 + t357 * t460;
t225 = pkin(17) - atan2(t280 * t421 - t494 * t422, t280 * t422 + t421 * t494) - t255;
t223 = sin(t225);
t570 = Icges(9,4) * t223;
t224 = cos(t225);
t569 = Icges(9,4) * t224;
t568 = t217 * t461;
t567 = t217 * t467;
t458 = sin(qJ(4));
t565 = t458 * t461;
t564 = t458 * t467;
t563 = t461 * t424;
t464 = cos(qJ(4));
t562 = t461 * t464;
t561 = t464 * t467;
t412 = pkin(12) + t587;
t399 = t467 * t412;
t560 = t467 * t424;
t559 = qJD(1) * t412;
t558 = qJD(4) * t217;
t343 = atan2(sin(t436), -cos(t436));
t336 = pkin(17) - t343;
t557 = t466 * t585 + V_base(3);
t556 = V_base(5) * pkin(11) + V_base(1);
t555 = V_base(4) / 0.2e1;
t554 = V_base(5) / 0.2e1;
t553 = pkin(1) * V_base(5);
t552 = V_base(4) * pkin(8);
t551 = t460 * t585;
t550 = rSges(6,1) * V_base(4);
t549 = rSges(6,1) * V_base(5);
t548 = rSges(9,1) * V_base(4);
t547 = rSges(6,2) * V_base(4);
t546 = rSges(6,2) * V_base(5);
t545 = rSges(9,2) * V_base(4);
t544 = t461 * V_base(4);
t543 = t467 * V_base(5);
t470 = -rSges(6,3) - pkin(10);
t542 = t470 * V_base(4);
t541 = sin(t305) / 0.2e1 + sin(t306) / 0.2e1;
t540 = t601 + t600;
t539 = t293 / 0.2e1 + t601;
t538 = cos(t305) / 0.2e1 - cos(t306) / 0.2e1;
t297 = cos(t308);
t537 = -t297 / 0.2e1 + t599;
t536 = t297 / 0.2e1 + t599;
t535 = -qJ(2) + t336;
t275 = -atan2(-sin(t405), -cos(t405)) + t336;
t534 = -pkin(11) - t588;
t533 = t544 * t587 + t557;
t532 = t460 * t553 + t556;
t530 = t467 * t551;
t274 = -qJ(2) + t275;
t411 = -qJD(2) + t424;
t410 = qJD(2) + t424;
t406 = pkin(12) * t544;
t525 = t406 + t533;
t303 = sin(t312);
t304 = cos(t312);
t522 = -pkin(8) * t303 - t470 * t304;
t520 = -rSges(4,1) * t431 + rSges(4,2) * t426;
t337 = qJ(2) + t343;
t519 = rSges(8,1) * cos(t337) - rSges(8,2) * sin(t337);
t518 = rSges(9,1) * cos(t274) + rSges(9,2) * sin(t274);
t517 = -V_base(4) * pkin(11) + V_base(2);
t515 = -Icges(4,1) * t431 + t577;
t511 = -Icges(9,1) * t224 - t570;
t509 = Icges(4,2) * t426 - t576;
t505 = -Icges(9,2) * t223 - t569;
t503 = -Icges(4,5) * t431 + Icges(4,6) * t426;
t499 = -Icges(9,5) * t224 - Icges(9,6) * t223;
t498 = -sin(t275) * (rSges(9,1) * t466 - rSges(9,2) * t460) + cos(t275) * (rSges(9,1) * t460 + rSges(9,2) * t466);
t496 = sin(t336) * t466 - cos(t336) * t460;
t339 = sin(t343);
t340 = cos(t343);
t393 = rSges(8,1) * t460 + rSges(8,2) * t466;
t395 = rSges(8,1) * t466 - rSges(8,2) * t460;
t495 = t339 * t395 + t340 * t393;
t391 = -t468 * rSges(7,1) + rSges(7,2) * t462;
t397 = rSges(7,1) * t462 + rSges(7,2) * t468;
t493 = t391 * t463 + t397 * t469;
t361 = t412 * V_base(6);
t341 = t361 * t467;
t492 = t399 * qJD(1) + t341 + V_base(2);
t490 = t412 * V_base(5);
t427 = sin(t450);
t489 = (t427 / 0.2e1 + t283 / 0.4e1 + t282 / 0.4e1) * t408;
t433 = cos(t451);
t488 = (t433 / 0.2e1 - t288 / 0.4e1 - t289 / 0.4e1) * t409;
t452 = qJ(1) + qJ(2);
t429 = sin(t452);
t453 = qJ(1) - qJ(2);
t430 = sin(t453);
t434 = cos(t452);
t435 = cos(t453);
t487 = t557 + (t429 + t430) * pkin(1) * t555 - (t434 + t435) * t553 / 0.2e1;
t486 = t543 - t544;
t371 = t447 * t461 + V_base(4);
t485 = (-Icges(9,3) * t467 + t499 * t461) * t370 + (Icges(9,3) * t461 + t499 * t467) * t371 + (Icges(9,5) * t223 - Icges(9,6) * t224) * t424;
t353 = qJD(3) * t461 + t402;
t483 = -(-Icges(4,3) * t467 + t503 * t461) * t370 - (Icges(4,3) * t461 + t503 * t467) * t353 - (-Icges(4,5) * t426 - Icges(4,6) * t431) * t424;
t481 = t532 + (-t361 - t559) * t461;
t480 = t424 * rSges(8,3) - t551;
t479 = t424 * rSges(9,3) - t551;
t201 = -Icges(9,6) * t467 + t505 * t461;
t202 = Icges(9,6) * t461 + t505 * t467;
t203 = -Icges(9,5) * t467 + t511 * t461;
t204 = Icges(9,5) * t461 + t511 * t467;
t206 = -Icges(9,2) * t224 + t570;
t207 = Icges(9,1) * t223 - t569;
t476 = (-t202 * t223 - t204 * t224) * t371 + (-t201 * t223 - t203 * t224) * t370 + (-t206 * t223 - t207 * t224) * t424;
t315 = -Icges(4,6) * t467 + t509 * t461;
t316 = Icges(4,6) * t461 + t509 * t467;
t317 = -Icges(4,5) * t467 + t515 * t461;
t318 = Icges(4,5) * t461 + t515 * t467;
t345 = -Icges(4,2) * t431 - t577;
t346 = -Icges(4,1) * t426 - t576;
t474 = (t316 * t426 - t318 * t431) * t353 + (t315 * t426 - t317 * t431) * t370 + (t345 * t426 - t346 * t431) * t424;
t446 = rSges(9,1) * V_base(5);
t445 = rSges(9,2) * V_base(5);
t444 = V_base(4) * pkin(12);
t441 = t467 * rSges(6,2);
t439 = qJ(1) - t449;
t438 = qJ(1) + t449;
t423 = 0.2e1 * V_base(5) * pkin(8);
t419 = -t584 / 0.2e1;
t417 = cos(t439);
t416 = cos(t438);
t415 = sin(t439);
t414 = sin(t438);
t413 = V_base(5) * t470;
t404 = -qJD(3) + t411;
t403 = qJD(3) + t410;
t396 = rSges(2,1) * t467 - rSges(2,2) * t461;
t390 = rSges(2,1) * t461 + rSges(2,2) * t467;
t389 = rSges(3,1) * t460 + rSges(3,2) * t466;
t378 = Icges(2,5) * t467 - Icges(2,6) * t461;
t377 = Icges(2,5) * t461 + Icges(2,6) * t467;
t369 = t467 * t559;
t364 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t363 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t362 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t358 = -t459 * t460 + t465 * t466;
t354 = t459 * t466 + t460 * t465;
t347 = -rSges(4,1) * t426 - rSges(4,2) * t431;
t331 = -qJ(1) + t535;
t330 = qJ(1) + t535;
t329 = cos(t535);
t320 = rSges(4,3) * t461 + t520 * t467;
t319 = -rSges(4,3) * t467 + t520 * t461;
t311 = -t391 * t469 + t397 * t463;
t310 = V_base(5) * rSges(2,3) - t390 * t424 + t556;
t309 = t396 * t424 + V_base(2) + (-pkin(11) - rSges(2,3)) * V_base(4);
t302 = t390 * V_base(4) - t396 * V_base(5) + V_base(3);
t271 = -qJ(1) + t274;
t270 = qJ(1) + t274;
t267 = -rSges(5,1) * t303 - rSges(5,2) * t304;
t266 = t401 * t389 + (-t491 * t461 + t581) * t424 + t556;
t265 = -t402 * t389 + t604 * t424 + t517;
t264 = t311 * t466 - t460 * t493;
t263 = t460 * t311 + t493 * t466;
t262 = t521 * qJD(2) + V_base(3) + t444 * t461 + (t521 * t461 - t581) * V_base(4) - t604 * V_base(5);
t261 = rSges(7,3) * t461 + t264 * t467;
t260 = -rSges(7,3) * t467 + t264 * t461;
t245 = t320 * t424 - t347 * t353 - t402 * t588 + t341 + t369 + t517;
t244 = -t319 * t424 + t347 * t370 + t481 - t530;
t243 = t319 * t353 - t320 * t370 - t467 * t490 + t525;
t242 = (-V_base(4) * rSges(8,3) - t490) * t467 + (-V_base(5) * rSges(8,3) + t444) * t461 + (t339 * t393 - t340 * t395) * t486 + t533;
t241 = t537 * rSges(5,1) + t540 * rSges(5,2) + t461 * rSges(5,3);
t240 = -t539 * rSges(5,1) + t536 * rSges(5,2) - t467 * rSges(5,3);
t230 = V_base(5) * pkin(13) + t263 * t401 + (pkin(6) * t461 - t260) * t424 + t556;
t229 = -t263 * t402 + V_base(2) + (-pkin(11) - pkin(13)) * V_base(4) + (-pkin(6) * t467 + t261) * t424;
t228 = t486 * pkin(6) + t260 * t402 - t261 * t401 + V_base(3);
t227 = t480 * t461 + t519 * t560 + (-t495 + t534) * V_base(4) + t492;
t226 = t480 * t467 + t495 * V_base(5) - t519 * t563 + t481;
t222 = -t461 * t551 + V_base(6) * t399 + t241 * t424 + t369 + V_base(2) + (t353 * t354 - t358 * t560) * pkin(4) + (-t267 + t534) * V_base(4);
t221 = -pkin(4) * t370 * t354 + t267 * V_base(5) - t530 + t532 + (-t240 + (pkin(4) * t358 - t412) * t461) * t424;
t220 = -(t544 + t447) * t586 + V_base(4) * t240 + (-t241 + (-t412 + t586) * t467) * V_base(5) + t525;
t216 = qJD(4) * t218 + t424;
t215 = -t467 * t558 + V_base(4);
t214 = -t461 * t558 + V_base(5);
t213 = -t218 * t561 + t565;
t212 = t218 * t564 + t562;
t211 = -t218 * t562 - t564;
t210 = t218 * t565 - t561;
t209 = t479 * t461 + t498 * t353 + (pkin(3) * t329 - t518) * t560 + (t496 * pkin(3) + t534) * V_base(4) + t492;
t208 = t479 * t467 + t518 * t563 - t498 * t370 + (-t329 * t563 - t496 * V_base(5)) * pkin(3) + t481;
t198 = (t446 - t545) * cos(t271) / 0.2e1 + (t445 + t548) * sin(t271) / 0.2e1 + (t446 + t545) * cos(t270) / 0.2e1 + (t445 - t548) * sin(t270) / 0.2e1 - t518 * t447 + (-rSges(9,3) * V_base(4) - V_base(5) * pkin(12)) * t467 + (-rSges(9,3) * V_base(5) + t444) * t461 + ((-cos(t331) / 0.2e1 - cos(t330) / 0.2e1) * V_base(5) + (-sin(t331) / 0.2e1 + sin(t330) / 0.2e1) * V_base(4)) * pkin(3) + t487;
t188 = V_base(2) + pkin(12) * t560 + (-t403 * t416 / 0.2e1 - t404 * t417 / 0.2e1) * pkin(4) + (t410 * t434 / 0.2e1 + t411 * t435 / 0.2e1) * pkin(1) + (t537 * pkin(8) + t540 * t470) * t424 + (t489 + t611) * rSges(6,2) + (t488 - t612) * rSges(6,1) + (pkin(4) * t426 + t541 * rSges(6,1) + t538 * rSges(6,2) - t522 + t534) * V_base(4);
t187 = t522 * V_base(5) - pkin(12) * t563 + (-t410 * t429 / 0.2e1 - t411 * t430 / 0.2e1) * pkin(1) + (-V_base(5) * t426 + t403 * t414 / 0.2e1 + t404 * t415 / 0.2e1) * pkin(4) + (t539 * pkin(8) - t536 * t470) * t424 + (-t538 * V_base(5) + t488 + t612) * rSges(6,2) + (-t541 * V_base(5) + t489 - t611) * rSges(6,1) + t532;
t186 = Icges(6,5) * t218 + (-Icges(6,1) * t464 + Icges(6,4) * t458) * t217;
t185 = Icges(6,6) * t218 + (-Icges(6,4) * t464 + Icges(6,2) * t458) * t217;
t184 = Icges(6,3) * t218 + (-Icges(6,5) * t464 + Icges(6,6) * t458) * t217;
t183 = Icges(6,1) * t213 + Icges(6,4) * t212 - Icges(6,5) * t567;
t182 = Icges(6,1) * t211 + Icges(6,4) * t210 - Icges(6,5) * t568;
t181 = Icges(6,4) * t213 + Icges(6,2) * t212 - Icges(6,6) * t567;
t180 = Icges(6,4) * t211 + Icges(6,2) * t210 - Icges(6,6) * t568;
t179 = Icges(6,5) * t213 + Icges(6,6) * t212 - Icges(6,3) * t567;
t178 = Icges(6,5) * t211 + Icges(6,6) * t210 - Icges(6,3) * t568;
t177 = t406 + (t423 + 0.2e1 * t542) * t297 / 0.4e1 + (t423 - 0.2e1 * t542) * t296 / 0.4e1 - pkin(12) * t543 + t487 + (-t413 + t552) * t600 + (-t413 - t552) * t601 - t447 * t586 + ((t441 + t584) * t433 + (t440 - t582) * t428 + (t441 - t584) * t432 + t427 * (t440 + t582) + t607 * (t419 - t441 / 0.2e1) + t608 * (t419 + t441 / 0.2e1) - t609 * (t589 - t582 / 0.2e1) + t610 * (t589 + t582 / 0.2e1) + ((-t296 + t297) * t467 + (-t292 + t293) * t461) * rSges(6,3)) * t558 / 0.2e1 + (t590 + t608 / 0.4e1) * (-t547 + t549) + (-t433 / 0.2e1 + t607 / 0.4e1) * (t547 + t549) + (t591 - t609 / 0.4e1) * (-t546 + t550) + (-t427 / 0.2e1 - t610 / 0.4e1) * (t546 + t550) + (-(t415 + t414) * V_base(4) / 0.2e1 + (t417 + t416) * t554) * pkin(4);
t1 = t353 * (-t483 * t461 + t474 * t467) / 0.2e1 + t371 * (t485 * t461 + t476 * t467) / 0.2e1 + m(8) * (t226 ^ 2 + t227 ^ 2 + t242 ^ 2) / 0.2e1 + m(4) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + m(5) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(7) * (t228 ^ 2 + t229 ^ 2 + t230 ^ 2) / 0.2e1 + t215 * ((-t179 * t567 + t212 * t181 + t213 * t183) * t215 + (-t178 * t567 + t180 * t212 + t182 * t213) * t214 + (-t184 * t567 + t185 * t212 + t186 * t213) * t216) / 0.2e1 + t214 * ((-t179 * t568 + t181 * t210 + t183 * t211) * t215 + (-t178 * t568 + t210 * t180 + t211 * t182) * t214 + (-t184 * t568 + t185 * t210 + t186 * t211) * t216) / 0.2e1 + t216 * ((t178 * t214 + t179 * t215 + t184 * t216) * t218 + ((t181 * t458 - t183 * t464) * t215 + (t180 * t458 - t182 * t464) * t214 + (t185 * t458 - t186 * t464) * t216) * t217) / 0.2e1 + m(6) * (t177 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(9) * (t198 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + V_base(6) * (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6)) / 0.2e1 + m(1) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(3) * (t262 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(2) * (t302 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + (t606 * t461 - t605 * t467) * t401 / 0.2e1 + (t605 * t461 + t606 * t467) * t402 / 0.2e1 + (Icges(1,1) * V_base(4) + Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + t378 * t424 - t603 * t461 + t602 * t467) * t555 + (Icges(1,4) * V_base(4) + Icges(1,2) * V_base(5) + Icges(1,6) * V_base(6) + t377 * t424 + t602 * t461 + t603 * t467) * t554 + ((-t316 * t431 - t318 * t426) * t353 + (-t202 * t224 + t204 * t223) * t371 + (t248 * t531 + t250 * t278 + t326 * t466 + t328 * t460) * t402 + (t249 * t531 + t251 * t278 + t325 * t466 + t327 * t460) * t401 + (-t191 * t218 - t193 * t217 + t233 * t254 + t235 * t253 + t377) * V_base(5) + (-t192 * t218 - t194 * t217 + t234 * t254 + t236 * t253 + t378) * V_base(4) + (-t196 * t218 - t197 * t217 - t206 * t224 + t207 * t223 + t238 * t254 + t239 * t253 + t257 * t531 + t258 * t278 - t345 * t431 - t346 * t426 + t379 * t466 + t382 * t460 + Icges(2,3)) * t424 + (-t201 * t224 + t203 * t223 - t315 * t431 - t317 * t426) * t370) * t424 / 0.2e1 + ((t483 - t485) * t467 + (t474 + t476) * t461) * t370 / 0.2e1;
T = t1;
