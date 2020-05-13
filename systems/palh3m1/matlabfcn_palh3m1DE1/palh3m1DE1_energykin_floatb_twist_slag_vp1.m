% Calculate kinetic energy for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1DE1_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:25:55
% EndTime: 2020-04-18 10:28:00
% DurationCPUTime: 119.96s
% Computational Cost: add. (2822013->606), mult. (4301224->986), div. (192288->22), fcn. (2706380->36), ass. (0->361)
t686 = -2 * pkin(1);
t685 = -pkin(6) - pkin(2);
t684 = -pkin(6) + pkin(2);
t683 = -pkin(8) - pkin(10);
t682 = -pkin(8) + pkin(10);
t576 = pkin(3) ^ 2;
t575 = pkin(4) ^ 2;
t574 = pkin(5) ^ 2;
t579 = pkin(1) ^ 2;
t558 = sin(qJ(2));
t560 = sin(pkin(16));
t564 = cos(qJ(2));
t566 = cos(pkin(16));
t518 = t558 * t560 - t564 * t566;
t664 = pkin(5) * t518;
t640 = t664 * t686 + t579;
t504 = t574 + t640;
t638 = pkin(2) ^ 2 - pkin(6) ^ 2;
t493 = t504 + t638;
t509 = pkin(1) - t664;
t520 = t558 * t566 + t560 * t564;
t489 = (pkin(5) - t685) * (pkin(5) + t685) + t640;
t490 = (pkin(5) - t684) * (pkin(5) + t684) + t640;
t581 = sqrt(-t490 * t489);
t464 = pkin(5) * t493 * t520 + t509 * t581;
t557 = sin(qJ(3));
t651 = t464 * t557;
t643 = t520 * t581;
t463 = -pkin(5) * t643 + t493 * t509;
t563 = cos(qJ(3));
t652 = t463 * t563;
t596 = -t651 / 0.2e1 + t652 / 0.2e1;
t500 = 0.1e1 / t504;
t578 = 0.1e1 / pkin(2);
t645 = t500 * t578;
t456 = t596 * t645;
t650 = t464 * t563;
t653 = t463 * t557;
t595 = t650 / 0.2e1 + t653 / 0.2e1;
t457 = t595 * t645;
t549 = pkin(18) + pkin(19);
t542 = sin(t549);
t543 = cos(t549);
t438 = t456 * t543 - t457 * t542;
t668 = pkin(4) * t438;
t641 = -0.2e1 * pkin(3) * t668 + t575;
t433 = t576 + t641;
t431 = 0.1e1 / t433;
t681 = t431 / 0.2e1;
t512 = t520 * qJD(2);
t513 = t518 * qJD(2);
t635 = pkin(1) * pkin(5) * t512;
t649 = 0.2e1 * (t489 + t490) * t635 / t581;
t622 = -t649 / 0.2e1;
t593 = t513 * t581 + t520 * t622;
t448 = ((t509 * t686 - t493) * t512 + t593) * pkin(5);
t680 = -t448 / 0.2e1;
t633 = -0.2e1 * t512 * t520;
t644 = t512 * t581;
t449 = t509 * t649 / 0.2e1 + t574 * pkin(1) * t633 + (-t493 * t513 - t644) * pkin(5);
t679 = t449 / 0.2e1;
t550 = sin(pkin(17));
t677 = t550 / 0.2e1;
t552 = sin(pkin(19));
t676 = t552 / 0.2e1;
t675 = t557 / 0.2e1;
t567 = cos(pkin(15));
t674 = t567 / 0.2e1;
t673 = pkin(1) * t558;
t672 = pkin(1) * t564;
t619 = 0.1e1 / t504 ^ 2 * t635;
t405 = ((t650 + t653) * t619 + (qJD(3) * t596 + t448 * t675 + t563 * t679) * t500) * t578;
t406 = ((t651 - t652) * t619 + (qJD(3) * t595 + t449 * t675 + t563 * t680) * t500) * t578;
t389 = -t405 * t542 - t406 * t543;
t671 = pkin(3) * t389;
t639 = pkin(8) ^ 2 - pkin(10) ^ 2;
t429 = t433 + t639;
t434 = -pkin(3) + t668;
t427 = (pkin(3) - t683) * (pkin(3) + t683) + t641;
t428 = (pkin(3) - t682) * (pkin(3) + t682) + t641;
t580 = sqrt(-t428 * t427);
t437 = -t456 * t542 - t457 * t543;
t669 = pkin(4) * t437;
t386 = t429 * t669 - t434 * t580;
t670 = pkin(4) * t386;
t517 = t557 * t558 - t563 * t564;
t559 = sin(qJ(1));
t506 = t517 * t559;
t667 = pkin(4) * t506;
t565 = cos(qJ(1));
t508 = t517 * t565;
t666 = pkin(4) * t508;
t598 = t557 * t564 + t558 * t563;
t665 = pkin(4) * t598;
t663 = Icges(2,4) * t559;
t662 = Icges(3,4) * t558;
t661 = Icges(3,4) * t564;
t492 = t504 - t638;
t510 = pkin(1) * t518 - pkin(5);
t462 = -pkin(1) * t643 - t492 * t510;
t465 = pkin(1) * t492 * t520 - t510 * t581;
t561 = sin(pkin(15));
t573 = 0.1e1 / pkin(6);
t646 = t500 * t573;
t455 = (t462 * t674 + t465 * t561 / 0.2e1) * t646;
t458 = (t465 * t674 - t462 * t561 / 0.2e1) * t646;
t446 = atan2(t458, t455);
t442 = sin(t446);
t660 = Icges(7,4) * t442;
t443 = cos(t446);
t659 = Icges(7,4) * t443;
t634 = pkin(4) * t671;
t658 = 0.2e1 * (t427 + t428) * t634 / t580;
t657 = t389 * t580;
t551 = cos(pkin(17));
t656 = t431 * t551;
t569 = 0.1e1 / pkin(10);
t655 = t431 * t569;
t654 = t437 * t580;
t554 = cos(pkin(19));
t648 = t500 * t554;
t647 = t500 * t561;
t544 = V_base(6) + qJD(1);
t642 = t544 * t565;
t453 = (-t463 * t554 / 0.2e1 + t464 * t676) * t645;
t451 = 0.1e1 / t453 ^ 2;
t454 = (t464 * t554 / 0.2e1 + t463 * t676) * t645;
t609 = t464 * t619;
t610 = t463 * t619;
t621 = t500 * t676;
t381 = ((t448 * t621 + t552 * t610 + t554 * t609 + t648 * t679) / t453 - (t449 * t621 + t552 * t609 - t554 * t610 + t648 * t680) * t454 * t451) / (t451 * t454 ^ 2 + 0.1e1) * t578;
t637 = -qJD(2) - t381;
t636 = -qJD(2) - qJD(3);
t627 = t559 * V_base(4);
t632 = pkin(13) * t627 + V_base(3);
t631 = V_base(5) * pkin(12) + V_base(1);
t430 = t433 - t639;
t435 = -pkin(3) * t438 + pkin(4);
t385 = -pkin(3) * t654 + t430 * t435;
t387 = pkin(3) * t430 * t437 + t435 * t580;
t369 = (-t385 * t551 / 0.2e1 + t387 * t677) * t655;
t370 = (t387 * t551 / 0.2e1 + t385 * t677) * t655;
t628 = atan2(t370, t369);
t626 = V_base(5) * t565;
t625 = -t658 / 0.2e1;
t624 = t431 * t677;
t571 = 0.1e1 / pkin(8);
t623 = t571 * t681;
t620 = t500 * t674;
t540 = qJD(2) * t559 + V_base(4);
t432 = 0.1e1 / t433 ^ 2;
t618 = t432 * t634;
t539 = -qJD(2) * t565 + V_base(5);
t617 = t539 * t673 + t631;
t615 = cos(t628);
t377 = t381 * t559 + t540;
t516 = qJD(3) * t559 + t540;
t614 = t465 * t619;
t613 = t385 * t618;
t612 = t387 * t618;
t611 = t462 * t619;
t515 = t565 * t636 + V_base(5);
t608 = -t515 * t665 + t617;
t607 = (-pkin(13) - t672) * t559;
t606 = rSges(3,1) * t564 - rSges(3,2) * t558;
t605 = rSges(7,1) * t443 - rSges(7,2) * t442;
t388 = -t405 * t543 + t406 * t542;
t594 = -t388 * t580 + t437 * t625;
t354 = ((-0.2e1 * pkin(4) * t435 - t430) * t389 + t594) * pkin(3);
t355 = t435 * t658 / 0.2e1 - 0.2e1 * t576 * t389 * t669 + (t388 * t430 - t657) * pkin(3);
t368 = 0.1e1 / t369 ^ 2;
t294 = ((t355 * t656 / 0.2e1 + t551 * t612 + t354 * t624 + t550 * t613) / t369 - (-t354 * t656 / 0.2e1 - t551 * t613 + t355 * t624 + t550 * t612) * t370 * t368) / (t368 * t370 ^ 2 + 0.1e1) * t569;
t292 = t294 * t559 + t516;
t604 = Icges(3,1) * t564 - t662;
t603 = Icges(7,1) * t443 - t660;
t602 = -Icges(3,2) * t558 + t661;
t601 = -Icges(7,2) * t442 + t659;
t600 = Icges(3,5) * t564 - Icges(3,6) * t558;
t599 = Icges(7,5) * t443 - Icges(7,6) * t442;
t445 = atan2(t454, t453);
t439 = sin(t445);
t440 = cos(t445);
t422 = t439 * t564 + t440 * t558;
t421 = -t439 * t558 + t440 * t564;
t597 = -V_base(4) * pkin(12) + pkin(13) * t642 + V_base(2);
t291 = V_base(5) + (-t294 + t636) * t565;
t592 = t607 - t667;
t447 = ((0.2e1 * pkin(5) * t510 - t492) * t512 + t593) * pkin(1);
t450 = t510 * t622 + t579 * pkin(5) * t633 + (-t492 * t513 - t644) * pkin(1);
t452 = 0.1e1 / t455 ^ 2;
t382 = ((t450 * t620 + t567 * t614 - t447 * t647 / 0.2e1 - t561 * t611) / t455 - (t447 * t620 + t567 * t611 + t450 * t647 / 0.2e1 + t561 * t614) * t458 * t452) / (t452 * t458 ^ 2 + 0.1e1) * t573;
t379 = -t382 * t565 + V_base(5);
t380 = t382 * t559 + V_base(4);
t591 = (-Icges(7,3) * t565 + t559 * t599) * t379 + (Icges(7,3) * t559 + t565 * t599) * t380 + (Icges(7,5) * t442 + Icges(7,6) * t443) * t544;
t590 = (-Icges(3,3) * t565 + t559 * t600) * t539 + (Icges(3,3) * t559 + t565 * t600) * t540 + (Icges(3,5) * t558 + Icges(3,6) * t564) * t544;
t384 = -pkin(4) * t654 - t429 * t434;
t589 = atan2(t386 * t623, t384 * t623);
t588 = -t540 * t673 + t642 * t672 + t597;
t587 = sin(t589);
t586 = t540 * t559 * t672 + (-pkin(13) * V_base(5) - t539 * t672) * t565 + t632;
t585 = t516 * t665 + t544 * t666 + t588;
t584 = -t515 * t666 + t516 * t667 + t586;
t411 = -Icges(7,6) * t565 + t559 * t601;
t412 = Icges(7,6) * t559 + t565 * t601;
t413 = -Icges(7,5) * t565 + t559 * t603;
t414 = Icges(7,5) * t559 + t565 * t603;
t424 = Icges(7,2) * t443 + t660;
t425 = Icges(7,1) * t442 + t659;
t583 = (-t412 * t442 + t414 * t443) * t380 + (-t411 * t442 + t413 * t443) * t379 + (-t424 * t442 + t425 * t443) * t544;
t496 = -Icges(3,6) * t565 + t559 * t602;
t497 = Icges(3,6) * t559 + t565 * t602;
t498 = -Icges(3,5) * t565 + t559 * t604;
t499 = Icges(3,5) * t559 + t565 * t604;
t529 = Icges(3,2) * t564 + t662;
t532 = Icges(3,1) * t558 + t661;
t582 = (-t497 * t558 + t499 * t564) * t540 + (-t496 * t558 + t498 * t564) * t539 + (-t529 * t558 + t532 * t564) * t544;
t562 = cos(qJ(4));
t556 = sin(qJ(4));
t555 = cos(pkin(18));
t553 = sin(pkin(18));
t547 = Icges(2,4) * t565;
t537 = rSges(2,1) * t565 - rSges(2,2) * t559;
t536 = rSges(2,1) * t559 + rSges(2,2) * t565;
t535 = rSges(3,1) * t558 + rSges(3,2) * t564;
t534 = Icges(2,1) * t565 - t663;
t533 = Icges(2,1) * t559 + t547;
t531 = -Icges(2,2) * t559 + t547;
t530 = Icges(2,2) * t565 + t663;
t525 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t524 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t523 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t507 = t598 * t565;
t505 = t598 * t559;
t503 = rSges(3,3) * t559 + t565 * t606;
t502 = -rSges(3,3) * t565 + t559 * t606;
t488 = V_base(5) * rSges(2,3) - t536 * t544 + t631;
t487 = t537 * t544 + V_base(2) + (-pkin(12) - rSges(2,3)) * V_base(4);
t485 = t536 * V_base(4) - t537 * V_base(5) + V_base(3);
t484 = -rSges(4,1) * t598 + rSges(4,2) * t517;
t483 = -Icges(4,1) * t598 + Icges(4,4) * t517;
t482 = -Icges(4,4) * t598 + Icges(4,2) * t517;
t481 = -Icges(4,5) * t598 + Icges(4,6) * t517;
t479 = rSges(4,1) * t508 + rSges(4,2) * t507 + rSges(4,3) * t559;
t478 = rSges(4,1) * t506 + rSges(4,2) * t505 - rSges(4,3) * t565;
t477 = Icges(4,1) * t508 + Icges(4,4) * t507 + Icges(4,5) * t559;
t476 = Icges(4,1) * t506 + Icges(4,4) * t505 - Icges(4,5) * t565;
t475 = Icges(4,4) * t508 + Icges(4,2) * t507 + Icges(4,6) * t559;
t474 = Icges(4,4) * t506 + Icges(4,2) * t505 - Icges(4,6) * t565;
t473 = Icges(4,5) * t508 + Icges(4,6) * t507 + Icges(4,3) * t559;
t472 = Icges(4,5) * t506 + Icges(4,6) * t505 - Icges(4,3) * t565;
t469 = t535 * t539 + (-pkin(13) * t559 - t502) * t544 + t631;
t468 = t503 * t544 - t535 * t540 + t597;
t467 = -pkin(13) * t626 + t502 * t540 - t503 * t539 + t632;
t461 = t484 * t515 + (-t478 + t607) * t544 + t617;
t460 = t479 * t544 - t484 * t516 + t588;
t459 = t478 * t516 - t479 * t515 + t586;
t426 = rSges(7,1) * t442 + rSges(7,2) * t443;
t420 = rSges(7,3) * t559 + t565 * t605;
t419 = -rSges(7,3) * t565 + t559 * t605;
t418 = t421 * t565;
t417 = t422 * t565;
t416 = t421 * t559;
t415 = t422 * t559;
t404 = rSges(8,1) * t422 + rSges(8,2) * t421;
t403 = Icges(8,1) * t422 + Icges(8,4) * t421;
t402 = Icges(8,4) * t422 + Icges(8,2) * t421;
t401 = Icges(8,5) * t422 + Icges(8,6) * t421;
t400 = (-t421 * t553 + t422 * t555) * pkin(3);
t399 = rSges(8,1) * t418 - rSges(8,2) * t417 + rSges(8,3) * t559;
t398 = rSges(8,1) * t416 - rSges(8,2) * t415 - rSges(8,3) * t565;
t397 = Icges(8,1) * t418 - Icges(8,4) * t417 + Icges(8,5) * t559;
t396 = Icges(8,1) * t416 - Icges(8,4) * t415 - Icges(8,5) * t565;
t395 = Icges(8,4) * t418 - Icges(8,2) * t417 + Icges(8,6) * t559;
t394 = Icges(8,4) * t416 - Icges(8,2) * t415 - Icges(8,6) * t565;
t393 = Icges(8,5) * t418 - Icges(8,6) * t417 + Icges(8,3) * t559;
t392 = Icges(8,5) * t416 - Icges(8,6) * t415 - Icges(8,3) * t565;
t391 = (t417 * t553 + t418 * t555) * pkin(3);
t390 = (t415 * t553 + t416 * t555) * pkin(3);
t383 = 0.1e1 / t384 ^ 2;
t376 = t565 * t637 + V_base(5);
t374 = V_base(5) * pkin(14) + t379 * t426 + (pkin(7) * t559 - t419) * t544 + t631;
t373 = -t380 * t426 + V_base(2) + (-pkin(12) - pkin(14)) * V_base(4) + (-pkin(7) * t565 + t420) * t544;
t371 = cos(t589);
t366 = t376 * t404 + (-t398 + t607) * t544 + t617;
t365 = -t377 * t404 + t399 * t544 + t588;
t364 = -t379 * t420 + t380 * t419 + V_base(3) + (t626 - t627) * pkin(7);
t363 = -t376 * t399 + t377 * t398 + t586;
t361 = sin(t628);
t360 = -t555 * t371 - t553 * t587;
t359 = t371 * t553 - t555 * t587;
t353 = t517 * t361 - t598 * t615;
t352 = -t361 * t598 - t517 * t615;
t351 = qJD(4) * t352 + t544;
t350 = t507 * t361 + t508 * t615;
t349 = t361 * t508 - t507 * t615;
t348 = t505 * t361 + t506 * t615;
t347 = t361 * t506 - t505 * t615;
t346 = t350 * t562 + t556 * t559;
t345 = -t350 * t556 + t559 * t562;
t344 = t348 * t562 - t556 * t565;
t343 = -t348 * t556 - t562 * t565;
t342 = t359 * t421 + t360 * t422;
t341 = -t359 * t422 + t360 * t421;
t340 = -t359 * t417 + t360 * t418;
t339 = -t359 * t418 - t360 * t417;
t338 = -t359 * t415 + t360 * t416;
t337 = -t359 * t416 - t360 * t415;
t336 = pkin(9) * t353 + pkin(11) * t352;
t335 = rSges(5,1) * t353 - rSges(5,2) * t352;
t334 = Icges(5,1) * t353 - Icges(5,4) * t352;
t333 = Icges(5,4) * t353 - Icges(5,2) * t352;
t332 = Icges(5,5) * t353 - Icges(5,6) * t352;
t331 = pkin(9) * t350 + pkin(11) * t349;
t330 = pkin(9) * t348 + pkin(11) * t347;
t329 = rSges(5,1) * t350 - rSges(5,2) * t349 + rSges(5,3) * t559;
t328 = rSges(5,1) * t348 - rSges(5,2) * t347 - rSges(5,3) * t565;
t327 = Icges(5,1) * t350 - Icges(5,4) * t349 + Icges(5,5) * t559;
t326 = Icges(5,1) * t348 - Icges(5,4) * t347 - Icges(5,5) * t565;
t325 = Icges(5,4) * t350 - Icges(5,2) * t349 + Icges(5,6) * t559;
t324 = Icges(5,4) * t348 - Icges(5,2) * t347 - Icges(5,6) * t565;
t323 = Icges(5,5) * t350 - Icges(5,6) * t349 + Icges(5,3) * t559;
t322 = Icges(5,5) * t348 - Icges(5,6) * t347 - Icges(5,3) * t565;
t321 = rSges(9,1) * t342 + rSges(9,2) * t341;
t320 = Icges(9,1) * t342 + Icges(9,4) * t341;
t319 = Icges(9,4) * t342 + Icges(9,2) * t341;
t318 = Icges(9,5) * t342 + Icges(9,6) * t341;
t317 = rSges(9,1) * t340 + rSges(9,2) * t339 + rSges(9,3) * t559;
t316 = rSges(9,1) * t338 + rSges(9,2) * t337 - rSges(9,3) * t565;
t315 = Icges(9,1) * t340 + Icges(9,4) * t339 + Icges(9,5) * t559;
t314 = Icges(9,1) * t338 + Icges(9,4) * t337 - Icges(9,5) * t565;
t313 = Icges(9,4) * t340 + Icges(9,2) * t339 + Icges(9,6) * t559;
t312 = Icges(9,4) * t338 + Icges(9,2) * t337 - Icges(9,6) * t565;
t311 = Icges(9,5) * t340 + Icges(9,6) * t339 + Icges(9,3) * t559;
t310 = Icges(9,5) * t338 + Icges(9,6) * t337 - Icges(9,3) * t565;
t309 = 0.2e1 * (((t434 * t625 + (t388 * t429 - t657) * pkin(4)) * t681 + (-t431 * t437 * t575 + t432 * t670) * t671) / t384 - ((-t389 * t429 + t594) * t681 + (t384 * t432 + t431 * t434) * t671) * t383 * t670) * pkin(8) / (t383 * t386 ^ 2 + 0.1e1) * t433 * t571;
t308 = t309 * t559 + t377;
t307 = V_base(5) + (-t309 + t637) * t565;
t306 = rSges(6,3) * t352 + (rSges(6,1) * t562 - rSges(6,2) * t556) * t353;
t305 = Icges(6,5) * t352 + (Icges(6,1) * t562 - Icges(6,4) * t556) * t353;
t304 = Icges(6,6) * t352 + (Icges(6,4) * t562 - Icges(6,2) * t556) * t353;
t303 = Icges(6,3) * t352 + (Icges(6,5) * t562 - Icges(6,6) * t556) * t353;
t302 = rSges(6,1) * t346 + rSges(6,2) * t345 + rSges(6,3) * t349;
t301 = rSges(6,1) * t344 + rSges(6,2) * t343 + rSges(6,3) * t347;
t300 = Icges(6,1) * t346 + Icges(6,4) * t345 + Icges(6,5) * t349;
t299 = Icges(6,1) * t344 + Icges(6,4) * t343 + Icges(6,5) * t347;
t298 = Icges(6,4) * t346 + Icges(6,2) * t345 + Icges(6,6) * t349;
t297 = Icges(6,4) * t344 + Icges(6,2) * t343 + Icges(6,6) * t347;
t296 = Icges(6,5) * t346 + Icges(6,6) * t345 + Icges(6,3) * t349;
t295 = Icges(6,5) * t344 + Icges(6,6) * t343 + Icges(6,3) * t347;
t290 = qJD(4) * t349 + t292;
t289 = qJD(4) * t347 + t291;
t288 = t307 * t321 + t376 * t400 + (-t316 - t390 + t607) * t544 + t617;
t287 = -t308 * t321 - t377 * t400 + (t317 + t391) * t544 + t588;
t286 = t291 * t335 + (-t328 + t592) * t544 + t608;
t285 = -t292 * t335 + t329 * t544 + t585;
t284 = -t307 * t317 + t308 * t316 - t376 * t391 + t377 * t390 + t586;
t283 = -t291 * t329 + t292 * t328 + t584;
t282 = t289 * t306 + t291 * t336 - t301 * t351 + (-t330 + t592) * t544 + t608;
t281 = -t290 * t306 - t292 * t336 + t302 * t351 + t331 * t544 + t585;
t280 = -t289 * t302 + t290 * t301 - t291 * t331 + t292 * t330 + t584;
t1 = ((t530 * t565 + t533 * t559 + Icges(1,2)) * V_base(5) + (t531 * t565 + t534 * t559 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t540 * (t559 * t590 + t565 * t582) / 0.2e1 + t539 * (t559 * t582 - t565 * t590) / 0.2e1 + t380 * (t559 * t591 + t565 * t583) / 0.2e1 + t379 * (t559 * t583 - t591 * t565) / 0.2e1 + t376 * ((-t393 * t565 - t395 * t415 + t397 * t416) * t377 + (-t392 * t565 - t394 * t415 + t396 * t416) * t376 + (-t401 * t565 - t402 * t415 + t403 * t416) * t544) / 0.2e1 + t377 * ((t393 * t559 - t395 * t417 + t397 * t418) * t377 + (t392 * t559 - t394 * t417 + t396 * t418) * t376 + (t401 * t559 - t402 * t417 + t403 * t418) * t544) / 0.2e1 + ((-t530 * t559 + t533 * t565 + Icges(1,4)) * V_base(5) + (-t531 * t559 + t534 * t565 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(2) * (t485 ^ 2 + t487 ^ 2 + t488 ^ 2) / 0.2e1 + m(3) * (t467 ^ 2 + t468 ^ 2 + t469 ^ 2) / 0.2e1 + m(9) * (t284 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + t544 * V_base(5) * (Icges(2,5) * t559 + Icges(2,6) * t565) + t544 * V_base(4) * (Icges(2,5) * t565 - Icges(2,6) * t559) + m(1) * (t523 ^ 2 + t524 ^ 2 + t525 ^ 2) / 0.2e1 + t351 * ((t289 * t295 + t290 * t296 + t303 * t351) * t352 + ((-t298 * t556 + t300 * t562) * t290 + (-t297 * t556 + t299 * t562) * t289 + (-t304 * t556 + t305 * t562) * t351) * t353) / 0.2e1 + m(4) * (t459 ^ 2 + t460 ^ 2 + t461 ^ 2) / 0.2e1 + t290 * ((t296 * t349 + t298 * t345 + t300 * t346) * t290 + (t295 * t349 + t297 * t345 + t299 * t346) * t289 + (t303 * t349 + t304 * t345 + t305 * t346) * t351) / 0.2e1 + t289 * ((t296 * t347 + t298 * t343 + t300 * t344) * t290 + (t295 * t347 + t297 * t343 + t299 * t344) * t289 + (t303 * t347 + t304 * t343 + t305 * t344) * t351) / 0.2e1 + m(6) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(5) * (t283 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + ((t497 * t564 + t499 * t558) * t540 + (t496 * t564 + t498 * t558) * t539 + (t475 * t517 - t477 * t598) * t516 + (t474 * t517 - t476 * t598) * t515 + (t412 * t443 + t414 * t442) * t380 + (t411 * t443 + t413 * t442) * t379 + (t395 * t421 + t397 * t422) * t377 + (t394 * t421 + t396 * t422) * t376 + (-t325 * t352 + t327 * t353) * t292 + (-t324 * t352 + t326 * t353) * t291 + (t313 * t341 + t315 * t342) * t308 + (t312 * t341 + t314 * t342) * t307 + (t319 * t341 + t320 * t342 - t333 * t352 + t334 * t353 + t402 * t421 + t403 * t422 + t424 * t443 + t425 * t442 + t482 * t517 - t483 * t598 + t529 * t564 + t532 * t558 + Icges(2,3)) * t544) * t544 / 0.2e1 + m(8) * (t363 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + t308 * ((t311 * t559 + t313 * t339 + t315 * t340) * t308 + (t310 * t559 + t312 * t339 + t314 * t340) * t307 + (t318 * t559 + t319 * t339 + t320 * t340) * t544) / 0.2e1 + t292 * ((t323 * t559 - t325 * t349 + t327 * t350) * t292 + (t322 * t559 - t324 * t349 + t326 * t350) * t291 + (t332 * t559 - t333 * t349 + t334 * t350) * t544) / 0.2e1 + t516 * ((t473 * t559 + t475 * t507 + t477 * t508) * t516 + (t472 * t559 + t474 * t507 + t476 * t508) * t515 + (t481 * t559 + t482 * t507 + t483 * t508) * t544) / 0.2e1 + t307 * ((-t311 * t565 + t313 * t337 + t315 * t338) * t308 + (-t310 * t565 + t312 * t337 + t314 * t338) * t307 + (-t318 * t565 + t319 * t337 + t320 * t338) * t544) / 0.2e1 + t291 * ((-t323 * t565 - t325 * t347 + t327 * t348) * t292 + (-t322 * t565 - t324 * t347 + t326 * t348) * t291 + (-t332 * t565 - t333 * t347 + t334 * t348) * t544) / 0.2e1 + t515 * ((-t473 * t565 + t475 * t505 + t477 * t506) * t516 + (-t472 * t565 + t474 * t505 + t476 * t506) * t515 + (-t481 * t565 + t482 * t505 + t483 * t506) * t544) / 0.2e1;
T = t1;
