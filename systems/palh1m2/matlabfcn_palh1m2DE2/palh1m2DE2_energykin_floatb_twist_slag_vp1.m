% Calculate kinetic energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:54
% EndTime: 2020-05-02 20:58:10
% DurationCPUTime: 14.41s
% Computational Cost: add. (11823->717), mult. (17971->941), div. (0->0), fcn. (26578->92), ass. (0->395)
t631 = Icges(2,2) + Icges(8,3);
t630 = Icges(4,3) + Icges(9,3);
t430 = sin(pkin(19));
t433 = cos(pkin(19));
t435 = sin(qJ(3));
t441 = cos(qJ(3));
t292 = t430 * t441 + t433 * t435;
t294 = t430 * t435 - t433 * t441;
t245 = qJ(2) + atan2(t294, t292);
t243 = sin(t245);
t244 = cos(t245);
t423 = qJ(2) + qJ(3);
t399 = sin(t423);
t404 = cos(t423);
t629 = -Icges(4,5) * t404 + Icges(9,5) * t243 + Icges(4,6) * t399 + Icges(9,6) * t244;
t437 = sin(qJ(1));
t443 = cos(qJ(1));
t305 = -t437 * V_base(4) + V_base(5) * t443;
t472 = -qJD(2) + t305;
t628 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t220 = atan2(t294, -t292) + t245;
t217 = sin(t220);
t218 = cos(t220);
t438 = sin(pkin(18));
t439 = sin(pkin(17));
t444 = cos(pkin(18));
t445 = cos(pkin(17));
t301 = t438 * t445 - t439 * t444;
t304 = t438 * t439 + t444 * t445;
t436 = sin(qJ(2));
t442 = cos(qJ(2));
t241 = t301 * t436 + t304 * t442;
t517 = t301 * t442 - t304 * t436;
t627 = -Icges(3,5) * t436 + Icges(7,5) * t517 + Icges(10,5) * t217 - Icges(3,6) * t442 - Icges(7,6) * t241 + Icges(10,6) * t218;
t397 = V_base(6) + qJD(1);
t416 = V_base(6) * pkin(15);
t420 = pkin(22) + pkin(21);
t519 = pkin(18) - t420;
t378 = -pkin(20) + t519;
t362 = sin(t378);
t363 = cos(t378);
t607 = -rSges(5,1) * t363 + rSges(5,2) * t362;
t551 = t441 * t442;
t376 = pkin(5) * t551;
t387 = pkin(5) * t435 + pkin(1);
t613 = -t387 * t436 + t376;
t626 = (pkin(15) + t613) * qJD(1) + t607 * t397 + t613 * V_base(6) + t416;
t588 = pkin(1) * t436;
t306 = t588 * V_base(6) - t416;
t542 = V_base(5) * pkin(13) + V_base(1);
t587 = pkin(1) * t442;
t518 = V_base(5) * t587 + t542;
t385 = -pkin(15) + t588;
t546 = t385 * qJD(1);
t514 = t437 * t546 + t518;
t624 = t306 * t437 + t514;
t584 = pkin(5) * qJD(3);
t622 = qJD(2) * t613 + (-t435 * t436 + t551) * t584 + V_base(3);
t621 = Icges(7,4) * t517;
t366 = qJ(4) + t378;
t357 = qJ(1) + t366;
t322 = cos(t357);
t367 = -qJ(4) + t378;
t360 = -qJ(1) + t367;
t325 = cos(t360);
t381 = qJD(4) + t397;
t424 = qJ(1) + qJ(4);
t405 = cos(t424);
t589 = t405 / 0.2e1;
t620 = t381 * (t322 / 0.4e1 + t325 / 0.4e1 + t589);
t382 = -qJD(4) + t397;
t358 = qJ(1) + t367;
t319 = sin(t358);
t359 = -qJ(1) + t366;
t320 = sin(t359);
t425 = qJ(1) - qJ(4);
t401 = sin(t425);
t512 = t319 / 0.4e1 - t320 / 0.4e1 - t401 / 0.2e1;
t619 = t382 * t512;
t510 = t397 * t443;
t422 = pkin(18) - pkin(22);
t396 = -qJ(1) + t422;
t373 = sin(t396) / 0.2e1;
t395 = qJ(1) + t422;
t379 = sin(t395);
t280 = -t379 / 0.2e1 + t373;
t374 = -cos(t395) / 0.2e1;
t380 = cos(t396);
t281 = t380 / 0.2e1 + t374;
t573 = Icges(2,4) * t437;
t618 = Icges(8,5) * t280 + Icges(8,6) * t281 - t631 * t443 - t573;
t279 = t379 / 0.2e1 + t373;
t282 = -t380 / 0.2e1 + t374;
t409 = Icges(2,4) * t443;
t617 = Icges(8,5) * t282 + Icges(8,6) * t279 + t631 * t437 - t409;
t323 = cos(t358);
t324 = cos(t359);
t616 = t324 + t323;
t615 = t325 + t322;
t368 = qJ(1) + t378;
t351 = sin(t368);
t369 = -qJ(1) + t378;
t352 = sin(t369);
t614 = t351 + t352;
t565 = Icges(9,4) * t243;
t490 = -Icges(9,2) * t244 - t565;
t202 = -Icges(9,6) * t443 + t437 * t490;
t203 = Icges(9,6) * t437 + t443 * t490;
t564 = Icges(9,4) * t244;
t496 = -Icges(9,1) * t243 - t564;
t204 = -Icges(9,5) * t443 + t437 * t496;
t205 = Icges(9,5) * t437 + t443 * t496;
t213 = -Icges(9,2) * t243 + t564;
t214 = Icges(9,1) * t244 - t565;
t569 = Icges(4,4) * t404;
t493 = -Icges(4,2) * t399 + t569;
t253 = -Icges(4,6) * t443 + t437 * t493;
t254 = Icges(4,6) * t437 + t443 * t493;
t570 = Icges(4,4) * t399;
t499 = Icges(4,1) * t404 - t570;
t255 = -Icges(4,5) * t443 + t437 * t499;
t256 = Icges(4,5) * t437 + t443 * t499;
t289 = Icges(4,2) * t404 + t570;
t290 = Icges(4,1) * t399 + t569;
t604 = qJD(2) * t437 + V_base(4);
t299 = qJD(3) * t437 + t604;
t421 = qJD(2) + qJD(3);
t326 = -t421 * t443 + V_base(5);
t612 = (-t213 * t244 - t214 * t243 - t289 * t399 + t290 * t404) * t397 + (-t203 * t244 - t205 * t243 - t254 * t399 + t256 * t404) * t299 + (-t202 * t244 - t204 * t243 - t253 * t399 + t255 * t404) * t326;
t611 = (-Icges(4,5) * t399 - Icges(9,5) * t244 - Icges(4,6) * t404 + Icges(9,6) * t243) * t397 + (-t630 * t437 + t629 * t443) * t299 + (t629 * t437 + t630 * t443) * t326;
t308 = t387 * V_base(5);
t610 = (-pkin(15) - t376) * V_base(5) + t308 * t436;
t609 = rSges(10,1) * t217 + rSges(10,2) * t218;
t608 = rSges(10,1) * t218 - rSges(10,2) * t217;
t563 = Icges(10,4) * t217;
t489 = Icges(10,2) * t218 + t563;
t180 = -Icges(10,6) * t443 + t437 * t489;
t181 = Icges(10,6) * t437 + t443 * t489;
t562 = Icges(10,4) * t218;
t495 = Icges(10,1) * t217 + t562;
t182 = -Icges(10,5) * t443 + t437 * t495;
t183 = Icges(10,5) * t437 + t443 * t495;
t186 = Icges(10,2) * t217 - t562;
t187 = -Icges(10,1) * t218 + t563;
t491 = -Icges(7,2) * t241 + t621;
t196 = Icges(7,6) * t437 + t443 * t491;
t197 = -Icges(7,6) * t443 + t437 * t491;
t566 = Icges(7,4) * t241;
t497 = Icges(7,1) * t517 - t566;
t198 = Icges(7,5) * t437 + t443 * t497;
t199 = -Icges(7,5) * t443 + t437 * t497;
t207 = Icges(7,2) * t517 + t566;
t208 = Icges(7,1) * t241 + t621;
t572 = Icges(3,4) * t436;
t494 = -Icges(3,2) * t442 - t572;
t262 = -Icges(3,6) * t443 + t437 * t494;
t263 = Icges(3,6) * t437 + t443 * t494;
t571 = Icges(3,4) * t442;
t500 = -Icges(3,1) * t436 - t571;
t264 = -Icges(3,5) * t443 + t437 * t500;
t265 = Icges(3,5) * t437 + t443 * t500;
t332 = -Icges(3,2) * t436 + t571;
t335 = Icges(3,1) * t442 - t572;
t507 = qJD(2) * t443 - V_base(5);
t606 = (t181 * t218 + t183 * t217 - t196 * t241 + t198 * t517 - t263 * t442 - t265 * t436) * t604 + (-t180 * t218 - t182 * t217 + t197 * t241 - t199 * t517 + t262 * t442 + t264 * t436) * t507 + (t186 * t218 + t187 * t217 - t207 * t241 + t208 * t517 - t332 * t442 - t335 * t436) * t397;
t605 = (t628 * t437 + t627 * t443) * t604 + (-t627 * t437 + t628 * t443) * t507 + (Icges(3,5) * t442 + Icges(7,5) * t241 - Icges(10,5) * t218 - Icges(3,6) * t436 + Icges(7,6) * t517 + Icges(10,6) * t217) * t397;
t327 = t421 * t437 + V_base(4);
t343 = rSges(3,1) * t436 + rSges(3,2) * t442;
t474 = -pkin(15) + t343;
t603 = -t437 * rSges(3,3) + t443 * t474;
t318 = sin(t357);
t321 = sin(t360);
t400 = sin(t424);
t602 = t318 / 0.4e1 - t321 / 0.4e1 + t400 / 0.2e1;
t597 = -t352 / 0.2e1;
t355 = cos(t368);
t596 = -t355 / 0.2e1;
t356 = cos(t369);
t595 = -t356 / 0.2e1;
t383 = qJD(2) + t397;
t592 = -t383 / 0.2e1;
t384 = -qJD(2) + t397;
t591 = t384 / 0.2e1;
t586 = pkin(2) * t243;
t585 = pkin(2) * t244;
t577 = t437 * rSges(6,1);
t575 = t443 * rSges(6,1);
t574 = t443 * rSges(3,3);
t414 = t443 * pkin(15);
t390 = sin(t420);
t391 = cos(t420);
t429 = sin(pkin(20));
t432 = cos(pkin(20));
t293 = -t429 * t444 + t432 * t438;
t297 = t429 * t438 + t432 * t444;
t236 = t293 * t441 - t297 * t435;
t237 = t293 * t435 + t297 * t441;
t480 = t436 * t236 + t237 * t442;
t481 = t236 * t442 - t436 * t237;
t163 = atan2(t390 * t480 - t391 * t481, -t390 * t481 - t391 * t480) + t423;
t161 = sin(t163);
t568 = Icges(5,4) * t161;
t162 = cos(t163);
t567 = Icges(5,4) * t162;
t561 = t161 * t437;
t560 = t161 * t443;
t509 = t397 * t437;
t434 = sin(qJ(4));
t556 = t434 * t437;
t555 = t434 * t443;
t554 = t436 * t441;
t440 = cos(qJ(4));
t553 = t437 * t440;
t552 = t440 * t443;
t302 = -t435 * t444 + t438 * t441;
t303 = t435 * t438 + t441 * t444;
t239 = -t302 * t442 + t303 * t436;
t240 = t302 * t436 + t303 * t442;
t428 = sin(pkin(22));
t431 = cos(pkin(22));
t295 = t428 * t438 + t431 * t444;
t296 = t428 * t444 - t431 * t438;
t168 = -qJ(2) - atan2(t295 * t442 - t296 * t436, t295 * t436 + t296 * t442) + pkin(21) - atan2(t239 * t391 + t240 * t390, -t239 * t390 + t240 * t391);
t166 = sin(t168);
t549 = Icges(11,4) * t166;
t167 = cos(t168);
t548 = Icges(11,4) * t167;
t547 = qJD(4) * t161;
t545 = t436 * qJD(2);
t394 = -qJ(2) + t422;
t271 = pkin(21) - atan2(cos(t394), -sin(t394));
t541 = t397 * t586;
t540 = pkin(5) * t554;
t538 = V_base(5) / 0.2e1;
t537 = qJD(2) * t587;
t536 = rSges(6,1) * V_base(4);
t535 = rSges(6,1) * V_base(5);
t534 = rSges(6,2) * V_base(4);
t533 = rSges(6,2) * V_base(5);
t532 = V_base(4) * pkin(13);
t530 = V_base(5) * pkin(15);
t529 = t575 / 0.2e1;
t528 = rSges(11,1) * V_base(4);
t527 = rSges(11,2) * V_base(4);
t526 = rSges(11,2) * V_base(5);
t524 = t437 * V_base(6);
t523 = sin(t366) / 0.2e1 + sin(t367) / 0.2e1;
t522 = cos(t366) / 0.2e1 - cos(t367) / 0.2e1;
t521 = -pkin(13) - t587;
t520 = -qJ(2) + t271;
t370 = t519 - t423;
t235 = -atan2(-sin(t370), cos(t370)) + t271;
t228 = -qJ(2) + t235;
t513 = t397 * pkin(15);
t446 = pkin(11) + rSges(6,3);
t506 = -pkin(9) * t362 + t446 * t363;
t505 = rSges(4,1) * t404 - rSges(4,2) * t399;
t392 = sin(t422);
t393 = cos(t422);
t504 = rSges(8,1) * t393 - rSges(8,2) * t392;
t502 = V_base(2) - t532;
t501 = rSges(11,1) * sin(t228) - rSges(11,2) * cos(t228);
t498 = Icges(5,1) * t162 - t568;
t492 = -Icges(5,2) * t161 + t567;
t486 = Icges(5,5) * t162 - Icges(5,6) * t161;
t482 = sin(t235) * (rSges(11,1) * t436 + rSges(11,2) * t442) + cos(t235) * (rSges(11,1) * t442 - rSges(11,2) * t436);
t479 = sin(t271) * t436 + cos(t271) * t442;
t347 = -rSges(7,1) * t438 + rSges(7,2) * t444;
t348 = rSges(7,1) * t444 + rSges(7,2) * t438;
t478 = t347 * t439 - t348 * t445;
t477 = -Icges(11,1) * t166 + t548;
t476 = Icges(11,2) * t167 - t549;
t475 = -Icges(11,5) * t166 + Icges(11,6) * t167;
t473 = rSges(8,3) * t397 - t537;
t375 = t443 * t416;
t471 = t375 + t502;
t470 = t602 * t381;
t406 = cos(t425);
t469 = (-t323 / 0.4e1 - t324 / 0.4e1 + t406 / 0.2e1) * t382;
t468 = V_base(2) + (-t306 - t546) * t443;
t465 = (-Icges(11,3) * t443 + t437 * t475) * t326 + (Icges(11,3) * t437 + t443 * t475) * t327 + (-Icges(11,5) * t167 - Icges(11,6) * t166) * t397;
t459 = -rSges(5,1) * t362 - rSges(5,2) * t363 + t540;
t458 = rSges(11,3) * t397 - t537;
t309 = t387 * V_base(4);
t418 = V_base(4) * pkin(15);
t456 = -t309 * t436 + t376 * V_base(4) + t418;
t455 = (-Icges(5,3) * t443 + t437 * t486) * V_base(5) + (Icges(5,3) * t437 + t443 * t486) * V_base(4) + (Icges(5,5) * t161 + Icges(5,6) * t162) * t397;
t454 = rSges(5,3) * t397 - (t387 * t442 + t540) * qJD(2) - (t435 * t442 + t554) * t584;
t142 = -Icges(11,6) * t443 + t437 * t476;
t143 = Icges(11,6) * t437 + t443 * t476;
t144 = -Icges(11,5) * t443 + t437 * t477;
t145 = Icges(11,5) * t437 + t443 * t477;
t148 = -Icges(11,2) * t166 - t548;
t149 = -Icges(11,1) * t167 - t549;
t453 = (t143 * t167 - t145 * t166) * t327 + (t142 * t167 - t144 * t166) * t326 + (t148 * t167 - t149 * t166) * t397;
t133 = -Icges(5,6) * t443 + t437 * t492;
t134 = Icges(5,6) * t437 + t443 * t492;
t135 = -Icges(5,5) * t443 + t437 * t498;
t136 = Icges(5,5) * t437 + t443 * t498;
t138 = Icges(5,2) * t162 + t568;
t139 = Icges(5,1) * t161 + t567;
t447 = (-t134 * t161 + t136 * t162) * V_base(4) + (-t133 * t161 + t135 * t162) * V_base(5) + (-t138 * t161 + t139 * t162) * t397;
t427 = qJ(1) - qJ(2);
t426 = qJ(1) + qJ(2);
t415 = t443 * rSges(6,2);
t413 = t437 * rSges(6,2);
t412 = qJ(1) - t423;
t411 = qJ(1) + t423;
t410 = rSges(11,1) * V_base(5);
t408 = cos(t427);
t407 = cos(t426);
t403 = sin(t427);
t402 = sin(t426);
t388 = -t577 / 0.2e1;
t372 = -qJD(3) + t384;
t371 = qJD(3) + t383;
t346 = rSges(2,1) * t443 - rSges(2,2) * t437;
t345 = rSges(3,1) * t442 - rSges(3,2) * t436;
t344 = rSges(2,1) * t437 + rSges(2,2) * t443;
t337 = Icges(2,1) * t443 - t573;
t336 = Icges(2,1) * t437 + t409;
t331 = Icges(2,5) * t443 - Icges(2,6) * t437;
t330 = Icges(2,5) * t437 + Icges(2,6) * t443;
t312 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t311 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t310 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t291 = rSges(4,1) * t399 + rSges(4,2) * t404;
t287 = rSges(8,1) * t392 + rSges(8,2) * t393;
t285 = -Icges(8,1) * t392 - Icges(8,4) * t393;
t284 = -Icges(8,4) * t392 - Icges(8,2) * t393;
t283 = -Icges(8,5) * t392 - Icges(8,6) * t393;
t275 = -qJD(3) + t472;
t268 = -qJ(1) + t520;
t267 = qJ(1) + t520;
t266 = sin(t520);
t258 = rSges(4,3) * t437 + t443 * t505;
t257 = -rSges(4,3) * t443 + t437 * t505;
t250 = -t347 * t445 - t348 * t439;
t249 = V_base(5) * rSges(2,3) - t344 * t397 + t542;
t248 = t346 * t397 + V_base(2) + (-pkin(13) - rSges(2,3)) * V_base(4);
t247 = t344 * V_base(4) - t346 * V_base(5) + V_base(3);
t234 = Icges(8,1) * t282 + Icges(8,4) * t279 + Icges(8,5) * t437;
t233 = Icges(8,1) * t280 + Icges(8,4) * t281 - Icges(8,5) * t443;
t232 = Icges(8,4) * t282 + Icges(8,2) * t279 + Icges(8,6) * t437;
t231 = Icges(8,4) * t280 + Icges(8,2) * t281 - Icges(8,6) * t443;
t225 = -qJ(1) + t228;
t224 = qJ(1) + t228;
t221 = -t507 * t345 + (t437 * t474 + t574) * t397 + t542;
t219 = -t345 * t604 - t397 * t603 + t502;
t216 = t250 * t436 - t442 * t478;
t215 = t250 * t442 + t436 * t478;
t211 = -t343 * qJD(2) + t418 * t437 + V_base(3) + (-t343 * t437 - t574) * V_base(4) + t603 * V_base(5);
t210 = rSges(7,3) * t437 + t215 * t443;
t209 = -rSges(7,3) * t443 + t215 * t437;
t193 = V_base(3) - pkin(1) * t545 + (-V_base(4) * rSges(8,3) + t385 * V_base(5)) * t443 + (-V_base(5) * rSges(8,3) - t588 * V_base(4) + t418) * t437 + t504 * t305;
t192 = V_base(2) + t473 * t437 + (t287 + t521) * V_base(4) + (-t385 - t504) * t510;
t191 = -V_base(5) * t287 + t473 * t443 + (V_base(6) * t385 + t397 * t504) * t437 + t514;
t190 = t258 * t397 - t291 * t299 - t587 * t604 + t468 - t532;
t189 = -t257 * t397 + t291 * t326 - t443 * t537 + t624;
t188 = -t305 * pkin(15) + t257 * t299 - t258 * t326 + t472 * t588 + V_base(3);
t184 = V_base(3) + (-rSges(9,3) * V_base(4) - t530) * t443 + (-rSges(9,3) * V_base(5) + t418) * t437 + (rSges(9,1) * t243 + rSges(9,2) * t244) * t275;
t177 = (-t243 * t326 + t244 * t509) * rSges(9,2) + (t243 * t509 + t244 * t326) * rSges(9,1) + t542 + t397 * (rSges(9,3) * t443 - pkin(15) * t437);
t176 = (rSges(9,3) * t437 + t414) * qJD(1) + rSges(9,3) * t524 + (t243 * t327 - t244 * t510) * rSges(9,2) + (-t243 * t510 - t244 * t327) * rSges(9,1) + t471;
t175 = (-rSges(5,3) * V_base(4) + t610) * t443 + (-V_base(5) * rSges(5,3) + t456) * t437 - t607 * t305 + t622;
t174 = -V_base(5) * pkin(16) - t216 * t507 + (pkin(14) * t437 - t209) * t397 + t542;
t173 = -t216 * t604 + V_base(2) + (-pkin(13) + pkin(16)) * V_base(4) + (-pkin(14) * t443 + t210) * t397;
t172 = pkin(14) * t305 + t209 * t604 + t210 * t507 + V_base(3);
t171 = -t309 * t442 + V_base(2) + (-pkin(13) - t459) * V_base(4) + t626 * t443 + t454 * t437;
t170 = t308 * t442 - t626 * t437 + t454 * t443 + t459 * V_base(5) + t542;
t169 = t275 * t586 + V_base(3) + (-rSges(10,3) * V_base(4) - t530) * t443 + (-rSges(10,3) * V_base(5) + t418) * t437 - t609 * t472;
t165 = -t299 * t585 + (rSges(10,3) * t437 + t414) * qJD(1) + rSges(10,3) * t524 + (t397 * t609 - t541) * t443 + t471 + t608 * t604;
t164 = t326 * t585 + (-t513 + t541) * t437 + t542 + t608 * t507 + (t443 * rSges(10,3) - t437 * t609) * t397;
t160 = -qJD(4) * t162 + t397;
t159 = t443 * t547 + V_base(4);
t158 = t437 * t547 + V_base(5);
t157 = t162 * t552 + t556;
t156 = -t162 * t555 + t553;
t155 = t162 * t553 - t555;
t154 = -t162 * t556 - t552;
t153 = t458 * t437 + t482 * t299 + (pkin(4) * t266 - t501) * t510 + (-pkin(4) * t479 + t521) * V_base(4) + t468;
t152 = t458 * t443 + t501 * t509 - t482 * t326 + (-t266 * t509 + t479 * V_base(5)) * pkin(4) + t624;
t151 = qJD(1) * t414 + t375 + V_base(2) + (t371 * cos(t411) / 0.2e1 + t372 * cos(t412) / 0.2e1) * pkin(5) + (t402 * t592 + t403 * t591) * pkin(1) + ((-t351 / 0.2e1 + t597) * t446 + (t596 + t595) * pkin(9)) * t397 + (t470 - t619) * rSges(6,2) + (t469 - t620) * rSges(6,1) + (-pkin(5) * t399 + rSges(6,1) * t523 + rSges(6,2) * t522 - t506 + t521) * V_base(4);
t150 = t506 * V_base(5) - t437 * t513 + (t407 * t592 + t408 * t591) * pkin(1) + (-t371 * sin(t411) / 0.2e1 - t372 * sin(t412) / 0.2e1 + V_base(5) * t399) * pkin(5) + ((t596 + t356 / 0.2e1) * t446 + (t351 / 0.2e1 + t597) * pkin(9)) * t397 + (-t522 * V_base(5) + t469 + t620) * rSges(6,2) + (-t523 * V_base(5) + t470 + t619) * rSges(6,1) + t518;
t146 = (-t526 - t528) * cos(t225) / 0.2e1 + (t410 - t527) * sin(t225) / 0.2e1 + (-t526 + t528) * cos(t224) / 0.2e1 + (t410 + t527) * sin(t224) / 0.2e1 + V_base(3) - t501 * t421 + (-rSges(11,3) * V_base(4) - t530) * t443 + (-rSges(11,3) * V_base(5) + t418) * t437 + ((-sin(t268) / 0.2e1 - sin(t267) / 0.2e1) * V_base(5) + (cos(t268) / 0.2e1 - cos(t267) / 0.2e1) * V_base(4)) * pkin(4) + (-t545 + (-t403 / 0.2e1 + t402 / 0.2e1) * V_base(5) + (-t408 / 0.2e1 + t407 / 0.2e1) * V_base(4)) * pkin(1);
t130 = -Icges(6,5) * t162 + (Icges(6,1) * t440 - Icges(6,4) * t434) * t161;
t129 = -Icges(6,6) * t162 + (Icges(6,4) * t440 - Icges(6,2) * t434) * t161;
t128 = -Icges(6,3) * t162 + (Icges(6,5) * t440 - Icges(6,6) * t434) * t161;
t127 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t560;
t126 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t561;
t125 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t560;
t124 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t561;
t123 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t560;
t122 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t561;
t121 = -((t415 + t577) * t406 + (t413 - t575) * t401 + (t415 - t577) * t405 + t400 * (t413 + t575) + (-t320 + t319) * (t529 - t413 / 0.2e1) + (-t321 + t318) * (t529 + t413 / 0.2e1) + t616 * (-t415 / 0.2e1 + t388) + t615 * (t415 / 0.2e1 + t388) + ((-t355 + t356) * t443 - t614 * t437) * rSges(6,3)) * t547 / 0.2e1 + t610 * t443 + t456 * t437 + (t595 + t355 / 0.2e1) * V_base(4) * rSges(6,3) + t614 * rSges(6,3) * t538 + (pkin(9) * t363 + pkin(11) * t362) * t305 - t512 * (-t533 + t536) - t602 * (t533 + t536) + (t589 + t615 / 0.4e1) * (-t534 + t535) + (-t406 / 0.2e1 + t616 / 0.4e1) * (t534 + t535) + t622;
t1 = (-t611 * t437 + t612 * t443) * t299 / 0.2e1 + (Icges(1,6) * V_base(6) + t437 * t447 - t443 * t455 + (t280 * t285 + t281 * t284 - t283 * t443 + t330) * t397 + (t231 * t281 + t233 * t280 + t336 * t437 - t443 * t618 + Icges(1,2)) * V_base(5) + (t232 * t281 + t234 * t280 + t337 * t437 - t443 * t617 + Icges(1,4)) * V_base(4)) * t538 + t160 * ((-t122 * t158 - t123 * t159 - t128 * t160) * t162 + ((-t125 * t434 + t127 * t440) * t159 + (-t124 * t434 + t126 * t440) * t158 + (-t129 * t434 + t130 * t440) * t160) * t161) / 0.2e1 + V_base(6) * V_base(4) * Icges(1,5) + ((-t143 * t166 - t145 * t167) * t327 + (-t203 * t243 + t205 * t244 + t254 * t404 + t256 * t399) * t299 + (t133 * t162 + t135 * t161 - t231 * t393 - t233 * t392 + t330) * V_base(5) + (t134 * t162 + t136 * t161 - t232 * t393 - t234 * t392 + t331) * V_base(4) + (t181 * t217 - t183 * t218 + t196 * t517 + t198 * t241 - t263 * t436 + t265 * t442) * t604 - (t180 * t217 - t182 * t218 + t197 * t517 + t199 * t241 - t262 * t436 + t264 * t442) * t507 + (Icges(2,3) - t332 * t436 + t335 * t442 + t289 * t404 + t290 * t399 - t284 * t393 - t285 * t392 + t186 * t217 - t187 * t218 + t207 * t517 + t208 * t241 - t148 * t166 - t149 * t167 - t213 * t243 + t214 * t244 + t138 * t162 + t139 * t161) * t397 + (-t142 * t166 - t144 * t167 - t202 * t243 + t204 * t244 + t253 * t404 + t255 * t399) * t326) * t397 / 0.2e1 + t327 * (t437 * t465 + t453 * t443) / 0.2e1 + (t437 * t455 + t443 * t447 + (t279 * t284 + t282 * t285 + t283 * t437 + t331) * t397 + (t231 * t279 + t233 * t282 + t336 * t443 + t437 * t618 + Icges(1,4)) * V_base(5) + (t232 * t279 + t234 * t282 + t337 * t443 + t617 * t437 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + V_base(6) * (Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6)) / 0.2e1 + t159 * ((t123 * t560 + t156 * t125 + t157 * t127) * t159 + (t122 * t560 + t124 * t156 + t126 * t157) * t158 + (t128 * t560 + t129 * t156 + t130 * t157) * t160) / 0.2e1 + t158 * ((t123 * t561 + t125 * t154 + t127 * t155) * t159 + (t122 * t561 + t154 * t124 + t155 * t126) * t158 + (t128 * t561 + t129 * t154 + t130 * t155) * t160) / 0.2e1 + m(1) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 - (t606 * t437 - t605 * t443) * t507 / 0.2e1 + (t605 * t437 + t606 * t443) * t604 / 0.2e1 + m(6) * (t121 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(11) * (t146 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(10) * (t164 ^ 2 + t165 ^ 2 + t169 ^ 2) / 0.2e1 + m(7) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(5) * (t170 ^ 2 + t171 ^ 2 + t175 ^ 2) / 0.2e1 + m(9) * (t176 ^ 2 + t177 ^ 2 + t184 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(8) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(3) * (t211 ^ 2 + t219 ^ 2 + t221 ^ 2) / 0.2e1 + m(2) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + ((-t465 + t611) * t443 + (t453 + t612) * t437) * t326 / 0.2e1;
T = t1;
