% Calculate kinetic energy for
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:41
% EndTime: 2020-04-19 19:36:04
% DurationCPUTime: 78.72s
% Computational Cost: add. (1803541->543), mult. (2768798->898), div. (119448->22), fcn. (1735032->38), ass. (0->368)
t515 = sin(qJ(1));
t521 = cos(qJ(1));
t679 = -t515 * V_base(4) + V_base(5) * t521;
t678 = -2 * pkin(1);
t677 = -pkin(6) - pkin(2);
t676 = -pkin(6) + pkin(2);
t675 = -pkin(8) - pkin(10);
t674 = -pkin(8) + pkin(10);
t532 = pkin(3) ^ 2;
t531 = pkin(4) ^ 2;
t530 = pkin(5) ^ 2;
t535 = pkin(1) ^ 2;
t514 = sin(qJ(2));
t516 = sin(pkin(16));
t520 = cos(qJ(2));
t522 = cos(pkin(16));
t472 = t514 * t516 - t520 * t522;
t658 = pkin(5) * t472;
t618 = t658 * t678 + t535;
t457 = t530 + t618;
t615 = pkin(2) ^ 2 - pkin(6) ^ 2;
t445 = t457 + t615;
t459 = pkin(1) - t658;
t473 = t514 * t522 + t516 * t520;
t433 = (pkin(5) - t677) * (pkin(5) + t677) + t618;
t434 = (pkin(5) - t676) * (pkin(5) + t676) + t618;
t537 = sqrt(-t434 * t433);
t416 = pkin(5) * t445 * t473 + t459 * t537;
t513 = sin(qJ(3));
t633 = t416 * t513;
t625 = t473 * t537;
t415 = -pkin(5) * t625 + t445 * t459;
t519 = cos(qJ(3));
t634 = t415 * t519;
t559 = -t633 / 0.2e1 + t634 / 0.2e1;
t453 = 0.1e1 / t457;
t534 = 0.1e1 / pkin(2);
t627 = t453 * t534;
t410 = t559 * t627;
t632 = t416 * t519;
t635 = t415 * t513;
t558 = t632 / 0.2e1 + t635 / 0.2e1;
t411 = t558 * t627;
t506 = pkin(18) + pkin(19);
t496 = sin(t506);
t497 = cos(t506);
t389 = t410 * t497 - t411 * t496;
t660 = pkin(4) * t389;
t620 = -0.2e1 * pkin(3) * t660 + t531;
t383 = t532 + t620;
t381 = 0.1e1 / t383;
t673 = t381 / 0.2e1;
t462 = t473 * qJD(2);
t463 = t472 * qJD(2);
t611 = pkin(1) * pkin(5) * t462;
t631 = 0.2e1 * (t433 + t434) * t611 / t537;
t600 = -t631 / 0.2e1;
t556 = t463 * t537 + t473 * t600;
t402 = ((t459 * t678 - t445) * t462 + t556) * pkin(5);
t672 = -t402 / 0.2e1;
t609 = -0.2e1 * t462 * t473;
t626 = t462 * t537;
t403 = t459 * t631 / 0.2e1 + t530 * pkin(1) * t609 + (-t445 * t463 - t626) * pkin(5);
t671 = t403 / 0.2e1;
t508 = sin(pkin(17));
t669 = t508 / 0.2e1;
t510 = sin(pkin(19));
t668 = t510 / 0.2e1;
t667 = t513 / 0.2e1;
t523 = cos(pkin(15));
t666 = t523 / 0.2e1;
t665 = pkin(1) * t514;
t597 = 0.1e1 / t457 ^ 2 * t611;
t347 = ((t632 + t635) * t597 + (qJD(3) * t559 + t402 * t667 + t519 * t671) * t453) * t534;
t348 = ((t633 - t634) * t597 + (qJD(3) * t558 + t403 * t667 + t519 * t672) * t453) * t534;
t346 = -t347 * t496 - t348 * t497;
t664 = pkin(3) * t346;
t511 = cos(pkin(19));
t407 = (-t415 * t511 / 0.2e1 + t416 * t668) * t627;
t408 = (t416 * t511 / 0.2e1 + t415 * t668) * t627;
t394 = qJ(2) + atan2(t408, t407);
t393 = pkin(18) - t394;
t663 = pkin(3) * sin(t393);
t616 = pkin(8) ^ 2 - pkin(10) ^ 2;
t379 = t383 + t616;
t384 = -pkin(3) + t660;
t377 = (pkin(3) - t675) * (pkin(3) + t675) + t620;
t378 = (pkin(3) - t674) * (pkin(3) + t674) + t620;
t536 = sqrt(-t378 * t377);
t388 = -t410 * t496 - t411 * t497;
t661 = pkin(4) * t388;
t343 = t379 * t661 - t384 * t536;
t662 = pkin(4) * t343;
t507 = qJ(2) + qJ(3);
t501 = sin(t507);
t659 = pkin(4) * t501;
t657 = pkin(13) * t515;
t656 = t520 * pkin(1);
t655 = Icges(2,4) * t515;
t654 = Icges(3,4) * t514;
t653 = Icges(3,4) * t520;
t652 = Icges(4,4) * t501;
t502 = cos(t507);
t651 = Icges(4,4) * t502;
t380 = t383 - t616;
t385 = -pkin(3) * t389 + pkin(4);
t636 = t388 * t536;
t342 = -pkin(3) * t636 + t380 * t385;
t344 = pkin(3) * t380 * t388 + t385 * t536;
t509 = cos(pkin(17));
t525 = 0.1e1 / pkin(10);
t637 = t381 * t525;
t325 = (-t342 * t509 / 0.2e1 + t344 * t669) * t637;
t326 = (t344 * t509 / 0.2e1 + t342 * t669) * t637;
t317 = atan2(t326, t325) + t507;
t315 = sin(t317);
t650 = Icges(5,4) * t315;
t316 = cos(t317);
t649 = Icges(5,4) * t316;
t444 = t457 - t615;
t460 = pkin(1) * t472 - pkin(5);
t414 = -pkin(1) * t625 - t444 * t460;
t417 = pkin(1) * t444 * t473 - t460 * t537;
t517 = sin(pkin(15));
t529 = 0.1e1 / pkin(6);
t628 = t453 * t529;
t409 = (t414 * t666 + t417 * t517 / 0.2e1) * t628;
t412 = (t417 * t666 - t414 * t517 / 0.2e1) * t628;
t400 = atan2(t412, t409);
t396 = sin(t400);
t648 = Icges(7,4) * t396;
t397 = cos(t400);
t647 = Icges(7,4) * t397;
t391 = sin(t394);
t646 = Icges(8,4) * t391;
t392 = cos(t394);
t645 = Icges(8,4) * t392;
t341 = -pkin(4) * t636 - t379 * t384;
t527 = 0.1e1 / pkin(8);
t601 = t527 * t673;
t323 = -atan2(t343 * t601, t341 * t601) + t393;
t321 = sin(t323);
t644 = Icges(9,4) * t321;
t322 = cos(t323);
t643 = Icges(9,4) * t322;
t642 = t315 * t515;
t641 = t315 * t521;
t610 = pkin(4) * t664;
t640 = 0.2e1 * (t377 + t378) * t610 / t536;
t639 = t346 * t536;
t638 = t381 * t509;
t630 = t453 * t511;
t629 = t453 * t517;
t512 = sin(qJ(4));
t624 = t512 * t515;
t623 = t512 * t521;
t518 = cos(qJ(4));
t622 = t515 * t518;
t621 = t518 * t521;
t619 = pkin(3) * cos(t393);
t617 = pkin(4) * t502;
t614 = qJD(4) * t315;
t405 = 0.1e1 / t407 ^ 2;
t588 = t415 * t597;
t592 = t416 * t597;
t599 = t453 * t668;
t338 = ((t402 * t599 + t510 * t588 + t511 * t592 + t630 * t671) / t407 - (t403 * t599 + t510 * t592 - t511 * t588 + t630 * t672) * t408 * t405) / (t405 * t408 ^ 2 + 0.1e1) * t534;
t613 = -qJD(2) - t338;
t612 = -qJD(2) - qJD(3);
t608 = V_base(5) * pkin(12) + V_base(1);
t603 = -t640 / 0.2e1;
t602 = t381 * t669;
t598 = t453 * t666;
t493 = qJD(2) * t515 + V_base(4);
t498 = V_base(6) + qJD(1);
t382 = 0.1e1 / t383 ^ 2;
t596 = t382 * t610;
t464 = t656 * t515;
t595 = -t464 - t657;
t492 = -qJD(2) * t521 + V_base(5);
t594 = t492 * t665 + t608;
t334 = t338 * t515 + t493;
t471 = qJD(3) * t515 + t493;
t591 = t344 * t596;
t590 = t342 * t596;
t589 = t414 * t597;
t587 = t417 * t597;
t431 = t617 * t515;
t586 = t431 + t595;
t585 = -pkin(9) * t316 - pkin(11) * t315;
t584 = rSges(3,1) * t520 - rSges(3,2) * t514;
t583 = -rSges(4,1) * t502 + rSges(4,2) * t501;
t582 = -rSges(5,1) * t316 + rSges(5,2) * t315;
t581 = rSges(7,1) * t397 - rSges(7,2) * t396;
t580 = rSges(8,1) * t392 - rSges(8,2) * t391;
t579 = -rSges(9,1) * t322 - rSges(9,2) * t321;
t345 = -t347 * t497 + t348 * t496;
t557 = -t345 * t536 + t388 * t603;
t296 = ((-0.2e1 * pkin(4) * t385 - t380) * t346 + t557) * pkin(3);
t297 = t385 * t640 / 0.2e1 - 0.2e1 * t532 * t346 * t661 + (t345 * t380 - t639) * pkin(3);
t324 = 0.1e1 / t325 ^ 2;
t263 = ((t297 * t638 / 0.2e1 + t509 * t591 + t296 * t602 + t508 * t590) / t325 - (-t296 * t638 / 0.2e1 - t509 * t590 + t297 * t602 + t508 * t591) * t326 * t324) / (t324 * t326 ^ 2 + 0.1e1) * t525;
t261 = t263 * t515 + t471;
t578 = Icges(3,1) * t520 - t654;
t577 = -Icges(4,1) * t502 + t652;
t576 = -Icges(5,1) * t316 + t650;
t575 = Icges(7,1) * t397 - t648;
t574 = Icges(8,1) * t392 - t646;
t573 = -Icges(9,1) * t322 - t644;
t572 = -Icges(3,2) * t514 + t653;
t571 = Icges(4,2) * t501 - t651;
t570 = Icges(5,2) * t315 - t649;
t569 = -Icges(7,2) * t396 + t647;
t568 = -Icges(8,2) * t391 + t645;
t567 = -Icges(9,2) * t321 - t643;
t566 = Icges(3,5) * t520 - Icges(3,6) * t514;
t565 = -Icges(4,5) * t502 + Icges(4,6) * t501;
t564 = -Icges(5,5) * t316 + Icges(5,6) * t315;
t563 = Icges(7,5) * t397 - Icges(7,6) * t396;
t562 = Icges(8,5) * t392 - Icges(8,6) * t391;
t561 = -Icges(9,5) * t322 - Icges(9,6) * t321;
t560 = t498 * t521 * pkin(13) - V_base(4) * pkin(12) + V_base(2);
t470 = t521 * t612 + V_base(5);
t555 = -t470 * t659 + t594;
t554 = -t679 * pkin(13) + V_base(3);
t260 = V_base(5) + (-t263 + t612) * t521;
t553 = (-Icges(5,3) * t521 + t515 * t564) * t260 + (Icges(5,3) * t515 + t521 * t564) * t261 + (-Icges(5,5) * t315 - Icges(5,6) * t316) * t498;
t340 = 0.1e1 / t341 ^ 2;
t268 = 0.2e1 * (((t384 * t603 + (t345 * t379 - t639) * pkin(4)) * t673 + (-t381 * t388 * t531 + t382 * t662) * t664) / t341 - ((-t346 * t379 + t557) * t673 + (t341 * t382 + t381 * t384) * t664) * t340 * t662) * pkin(8) / (t340 * t343 ^ 2 + 0.1e1) * t383 * t527;
t266 = V_base(5) + (-t268 + t613) * t521;
t267 = t268 * t515 + t334;
t552 = (-Icges(9,3) * t521 + t515 * t561) * t266 + (Icges(9,3) * t515 + t521 * t561) * t267 + (Icges(9,5) * t321 - Icges(9,6) * t322) * t498;
t333 = t521 * t613 + V_base(5);
t551 = (-Icges(8,3) * t521 + t515 * t562) * t333 + (Icges(8,3) * t515 + t521 * t562) * t334 + (Icges(8,5) * t391 + Icges(8,6) * t392) * t498;
t401 = ((0.2e1 * pkin(5) * t460 - t444) * t462 + t556) * pkin(1);
t404 = t460 * t600 + t535 * pkin(5) * t609 + (-t444 * t463 - t626) * pkin(1);
t406 = 0.1e1 / t409 ^ 2;
t339 = ((t404 * t598 + t523 * t587 - t401 * t629 / 0.2e1 - t517 * t589) / t409 - (t401 * t598 + t523 * t589 + t404 * t629 / 0.2e1 + t517 * t587) * t412 * t406) / (t406 * t412 ^ 2 + 0.1e1) * t529;
t336 = -t339 * t521 + V_base(5);
t337 = t339 * t515 + V_base(4);
t550 = (-Icges(7,3) * t521 + t515 * t563) * t336 + (Icges(7,3) * t515 + t521 * t563) * t337 + (Icges(7,5) * t396 + Icges(7,6) * t397) * t498;
t549 = (-Icges(4,3) * t521 + t515 * t565) * t470 + (Icges(4,3) * t515 + t521 * t565) * t471 + (-Icges(4,5) * t501 - Icges(4,6) * t502) * t498;
t548 = (-Icges(3,3) * t521 + t515 * t566) * t492 + (Icges(3,3) * t515 + t521 * t566) * t493 + (Icges(3,5) * t514 + Icges(3,6) * t520) * t498;
t465 = t656 * t521;
t547 = t498 * t465 - t493 * t665 + t560;
t546 = t493 * t464 - t465 * t492 + t554;
t432 = t617 * t521;
t545 = -t498 * t432 + t471 * t659 + t547;
t544 = -t471 * t431 + t432 * t470 + t546;
t283 = -Icges(5,6) * t521 + t515 * t570;
t284 = Icges(5,6) * t515 + t521 * t570;
t285 = -Icges(5,5) * t521 + t515 * t576;
t286 = Icges(5,5) * t515 + t521 * t576;
t292 = -Icges(5,2) * t316 - t650;
t293 = -Icges(5,1) * t315 - t649;
t543 = (t284 * t315 - t286 * t316) * t261 + (t283 * t315 - t285 * t316) * t260 + (t292 * t315 - t293 * t316) * t498;
t300 = -Icges(9,6) * t521 + t515 * t567;
t301 = Icges(9,6) * t515 + t521 * t567;
t302 = -Icges(9,5) * t521 + t515 * t573;
t303 = Icges(9,5) * t515 + t521 * t573;
t307 = -Icges(9,2) * t322 + t644;
t308 = Icges(9,1) * t321 - t643;
t542 = (-t301 * t321 - t303 * t322) * t267 + (-t300 * t321 - t302 * t322) * t266 + (-t307 * t321 - t308 * t322) * t498;
t353 = -Icges(8,6) * t521 + t515 * t568;
t354 = Icges(8,6) * t515 + t521 * t568;
t355 = -Icges(8,5) * t521 + t515 * t574;
t356 = Icges(8,5) * t515 + t521 * t574;
t368 = Icges(8,2) * t392 + t646;
t369 = Icges(8,1) * t391 + t645;
t541 = (-t354 * t391 + t356 * t392) * t334 + (-t353 * t391 + t355 * t392) * t333 + (-t368 * t391 + t369 * t392) * t498;
t361 = -Icges(7,6) * t521 + t515 * t569;
t362 = Icges(7,6) * t515 + t521 * t569;
t363 = -Icges(7,5) * t521 + t515 * t575;
t364 = Icges(7,5) * t515 + t521 * t575;
t372 = Icges(7,2) * t397 + t648;
t373 = Icges(7,1) * t396 + t647;
t540 = (-t362 * t396 + t364 * t397) * t337 + (-t361 * t396 + t363 * t397) * t336 + (-t372 * t396 + t373 * t397) * t498;
t437 = -Icges(4,6) * t521 + t515 * t571;
t438 = Icges(4,6) * t515 + t521 * t571;
t439 = -Icges(4,5) * t521 + t515 * t577;
t440 = Icges(4,5) * t515 + t521 * t577;
t467 = -Icges(4,2) * t502 - t652;
t468 = -Icges(4,1) * t501 - t651;
t539 = (t438 * t501 - t440 * t502) * t471 + (t437 * t501 - t439 * t502) * t470 + (t467 * t501 - t468 * t502) * t498;
t449 = -Icges(3,6) * t521 + t515 * t572;
t450 = Icges(3,6) * t515 + t521 * t572;
t451 = -Icges(3,5) * t521 + t515 * t578;
t452 = Icges(3,5) * t515 + t521 * t578;
t482 = Icges(3,2) * t520 + t654;
t485 = Icges(3,1) * t514 + t653;
t538 = (-t450 * t514 + t452 * t520) * t493 + (-t449 * t514 + t451 * t520) * t492 + (-t482 * t514 + t485 * t520) * t498;
t503 = Icges(2,4) * t521;
t490 = rSges(2,1) * t521 - rSges(2,2) * t515;
t489 = rSges(2,1) * t515 + rSges(2,2) * t521;
t488 = rSges(3,1) * t514 + rSges(3,2) * t520;
t487 = Icges(2,1) * t521 - t655;
t486 = Icges(2,1) * t515 + t503;
t484 = -Icges(2,2) * t515 + t503;
t483 = Icges(2,2) * t521 + t655;
t478 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t477 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t476 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t469 = -rSges(4,1) * t501 - rSges(4,2) * t502;
t456 = rSges(3,3) * t515 + t521 * t584;
t455 = -rSges(3,3) * t521 + t515 * t584;
t443 = rSges(4,3) * t515 + t521 * t583;
t442 = -rSges(4,3) * t521 + t515 * t583;
t430 = V_base(5) * rSges(2,3) - t489 * t498 + t608;
t429 = t490 * t498 + V_base(2) + (-pkin(12) - rSges(2,3)) * V_base(4);
t428 = t489 * V_base(4) - t490 * V_base(5) + V_base(3);
t423 = t488 * t492 + (-t455 - t657) * t498 + t608;
t422 = t456 * t498 - t488 * t493 + t560;
t421 = t455 * t493 - t456 * t492 + t554;
t419 = t469 * t470 + (-t442 + t595) * t498 + t594;
t418 = t443 * t498 - t469 * t471 + t547;
t413 = t442 * t471 - t443 * t470 + t546;
t376 = t619 * t521;
t375 = t619 * t515;
t374 = rSges(7,1) * t396 + rSges(7,2) * t397;
t370 = rSges(8,1) * t391 + rSges(8,2) * t392;
t366 = rSges(7,3) * t515 + t521 * t581;
t365 = -rSges(7,3) * t521 + t515 * t581;
t358 = rSges(8,3) * t515 + t521 * t580;
t357 = -rSges(8,3) * t521 + t515 * t580;
t331 = V_base(5) * pkin(14) + t336 * t374 + (pkin(7) * t515 - t365) * t498 + t608;
t330 = -t337 * t374 + V_base(2) + (-pkin(12) - pkin(14)) * V_base(4) + (-pkin(7) * t521 + t366) * t498;
t329 = t333 * t370 + (-t357 + t595) * t498 + t594;
t328 = -t334 * t370 + t358 * t498 + t547;
t320 = t679 * pkin(7) - t336 * t366 + t337 * t365 + V_base(3);
t319 = -t333 * t358 + t334 * t357 + t546;
t314 = qJD(4) * t316 + t498;
t313 = -t316 * t621 + t624;
t312 = t316 * t623 + t622;
t311 = -t316 * t622 - t623;
t310 = t316 * t624 - t621;
t309 = rSges(9,1) * t321 - rSges(9,2) * t322;
t305 = rSges(9,3) * t515 + t521 * t579;
t304 = -rSges(9,3) * t521 + t515 * t579;
t295 = -pkin(9) * t315 + pkin(11) * t316;
t294 = -rSges(5,1) * t315 - rSges(5,2) * t316;
t290 = t585 * t521;
t289 = t585 * t515;
t288 = rSges(5,3) * t515 + t521 * t582;
t287 = -rSges(5,3) * t521 + t515 * t582;
t280 = rSges(6,3) * t316 + (-rSges(6,1) * t518 + rSges(6,2) * t512) * t315;
t279 = Icges(6,5) * t316 + (-Icges(6,1) * t518 + Icges(6,4) * t512) * t315;
t278 = Icges(6,6) * t316 + (-Icges(6,4) * t518 + Icges(6,2) * t512) * t315;
t277 = Icges(6,3) * t316 + (-Icges(6,5) * t518 + Icges(6,6) * t512) * t315;
t276 = rSges(6,1) * t313 + rSges(6,2) * t312 - rSges(6,3) * t641;
t275 = rSges(6,1) * t311 + rSges(6,2) * t310 - rSges(6,3) * t642;
t274 = Icges(6,1) * t313 + Icges(6,4) * t312 - Icges(6,5) * t641;
t273 = Icges(6,1) * t311 + Icges(6,4) * t310 - Icges(6,5) * t642;
t272 = Icges(6,4) * t313 + Icges(6,2) * t312 - Icges(6,6) * t641;
t271 = Icges(6,4) * t311 + Icges(6,2) * t310 - Icges(6,6) * t642;
t270 = Icges(6,5) * t313 + Icges(6,6) * t312 - Icges(6,3) * t641;
t269 = Icges(6,5) * t311 + Icges(6,6) * t310 - Icges(6,3) * t642;
t265 = -t333 * t663 + t266 * t309 + (-t304 - t375 + t595) * t498 + t594;
t264 = t334 * t663 - t267 * t309 + (t305 + t376) * t498 + t547;
t259 = -t521 * t614 + t261;
t258 = -t515 * t614 + t260;
t257 = -t266 * t305 + t267 * t304 - t333 * t376 + t334 * t375 + t546;
t256 = t260 * t294 + (-t287 + t586) * t498 + t555;
t255 = -t261 * t294 + t288 * t498 + t545;
t254 = -t260 * t288 + t261 * t287 + t544;
t253 = t258 * t280 + t260 * t295 - t275 * t314 + (-t289 + t586) * t498 + t555;
t252 = -t259 * t280 - t261 * t295 + t276 * t314 + t290 * t498 + t545;
t251 = -t258 * t276 + t259 * t275 - t260 * t290 + t261 * t289 + t544;
t1 = ((t450 * t520 + t452 * t514) * t493 + (t449 * t520 + t451 * t514) * t492 + (-t438 * t502 - t440 * t501) * t471 + (-t437 * t502 - t439 * t501) * t470 + (t362 * t397 + t364 * t396) * t337 + (t361 * t397 + t363 * t396) * t336 + (t354 * t392 + t356 * t391) * t334 + (t353 * t392 + t355 * t391) * t333 + (-t301 * t322 + t303 * t321) * t267 + (-t300 * t322 + t302 * t321) * t266 + (-t284 * t316 - t286 * t315) * t261 + (-t283 * t316 - t285 * t315) * t260 + (-t292 * t316 - t293 * t315 - t307 * t322 + t308 * t321 + t368 * t392 + t369 * t391 + t372 * t397 + t373 * t396 - t467 * t502 - t468 * t501 + t482 * t520 + t485 * t514 + Icges(2,3)) * t498) * t498 / 0.2e1 + ((-t483 * t515 + t486 * t521 + Icges(1,4)) * V_base(5) + (-t484 * t515 + t487 * t521 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t483 * t521 + t486 * t515 + Icges(1,2)) * V_base(5) + (t484 * t521 + t487 * t515 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t498 * V_base(5) * (Icges(2,5) * t515 + Icges(2,6) * t521) + t498 * V_base(4) * (Icges(2,5) * t521 - Icges(2,6) * t515) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(2) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + t493 * (t515 * t548 + t521 * t538) / 0.2e1 + t492 * (t515 * t538 - t521 * t548) / 0.2e1 + t471 * (t515 * t549 + t521 * t539) / 0.2e1 + t470 * (t515 * t539 - t521 * t549) / 0.2e1 + t337 * (t515 * t550 + t521 * t540) / 0.2e1 + t336 * (t515 * t540 - t521 * t550) / 0.2e1 + t334 * (t515 * t551 + t521 * t541) / 0.2e1 + t333 * (t515 * t541 - t521 * t551) / 0.2e1 + t267 * (t515 * t552 + t521 * t542) / 0.2e1 + t266 * (t515 * t542 - t521 * t552) / 0.2e1 + m(6) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + t261 * (t515 * t553 + t521 * t543) / 0.2e1 + t260 * (t515 * t543 - t521 * t553) / 0.2e1 + m(4) * (t413 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + m(3) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + t259 * ((-t270 * t641 + t312 * t272 + t313 * t274) * t259 + (-t269 * t641 + t271 * t312 + t273 * t313) * t258 + (-t277 * t641 + t278 * t312 + t279 * t313) * t314) / 0.2e1 + t258 * ((-t270 * t642 + t272 * t310 + t274 * t311) * t259 + (-t269 * t642 + t310 * t271 + t311 * t273) * t258 + (-t277 * t642 + t278 * t310 + t279 * t311) * t314) / 0.2e1 + m(9) * (t257 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t314 * ((t269 * t258 + t270 * t259 + t277 * t314) * t316 + ((t272 * t512 - t274 * t518) * t259 + (t271 * t512 - t273 * t518) * t258 + (t278 * t512 - t279 * t518) * t314) * t315) / 0.2e1 + m(7) * (t320 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(1) * (t476 ^ 2 + t477 ^ 2 + t478 ^ 2) / 0.2e1 + m(5) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + m(8) * (t319 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1;
T = t1;
