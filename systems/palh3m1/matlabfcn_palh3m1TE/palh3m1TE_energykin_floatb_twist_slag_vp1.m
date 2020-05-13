% Calculate kinetic energy for
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:19:09
% EndTime: 2020-04-17 15:20:42
% DurationCPUTime: 84.60s
% Computational Cost: add. (1917843->607), mult. (2943256->983), div. (127608->22), fcn. (1845788->24), ass. (0->348)
t633 = sin(qJ(2));
t634 = sin(pkin(16));
t636 = cos(qJ(2));
t637 = cos(pkin(16));
t561 = t633 * t634 - t636 * t637;
t559 = pkin(5) * t561;
t553 = (-0.2e1 * t559 + pkin(1)) * pkin(1);
t645 = -pkin(6) - pkin(2);
t646 = pkin(5) + pkin(6);
t647 = pkin(5) - pkin(6);
t535 = sqrt(-((-pkin(2) + t646) * (pkin(2) + t647) + t553) * ((pkin(5) - t645) * (pkin(5) + t645) + t553));
t536 = pkin(2) ^ 2;
t562 = t633 * t637 + t634 * t636;
t479 = t562 * qJD(2);
t648 = pkin(1) * pkin(5);
t605 = t479 * t648;
t616 = 0.4e1 / t535 * (t646 * t647 - t536 + t553) * t605;
t653 = -t616 / 0.2e1;
t652 = pkin(5) ^ 2;
t651 = 0.1e1 / pkin(6);
t650 = -0.2e1 * pkin(1);
t473 = t553 + t652;
t649 = 0.1e1 / t473;
t644 = -pkin(8) - pkin(10);
t643 = -pkin(8) + pkin(10);
t532 = pkin(3) ^ 2;
t531 = pkin(4) ^ 2;
t611 = pkin(6) ^ 2 - t536;
t463 = t473 - t611;
t478 = -t559 + pkin(1);
t558 = t535 * t562;
t438 = -pkin(5) * t558 + t478 * t463;
t560 = pkin(5) * t562;
t439 = t463 * t560 + t478 * t535;
t525 = cos(qJ(3));
t533 = 0.1e1 / pkin(2);
t522 = sin(qJ(3));
t609 = t649 / 0.2e1;
t589 = t522 * t609;
t610 = -t649 / 0.2e1;
t428 = (t438 * t525 * t610 + t439 * t589) * t533;
t429 = (t439 * t525 * t609 + t438 * t589) * t533;
t514 = pkin(18) + pkin(19);
t507 = sin(t514);
t508 = cos(t514);
t410 = -t428 * t508 - t429 * t507;
t629 = pkin(4) * t410;
t613 = -0.2e1 * pkin(3) * t629 + t531;
t393 = t532 + t613;
t391 = 0.1e1 / t393;
t642 = t391 / 0.2e1;
t515 = sin(pkin(17));
t640 = t515 / 0.2e1;
t520 = cos(pkin(18));
t639 = -t520 / 0.2e1;
t638 = cos(pkin(15));
t635 = sin(pkin(15));
t480 = t561 * qJD(2);
t549 = t480 * t535 + t562 * t653;
t420 = ((t478 * t650 - t463) * t479 + t549) * pkin(5);
t470 = 0.1e1 / t473 ^ 2;
t588 = t470 * t605;
t581 = t438 * t588;
t554 = qJD(3) * t439 * t610 + t420 * t609 + t581;
t551 = (t560 * t650 - t535) * t479;
t421 = t478 * t616 / 0.2e1 + (-t480 * t463 + t551) * pkin(5);
t580 = t439 * t588;
t555 = t580 + (qJD(3) * t438 + t421) * t609;
t368 = (t522 * t554 + t525 * t555) * t533;
t369 = (t522 * t555 - t525 * t554) * t533;
t367 = -t368 * t507 - t369 * t508;
t632 = pkin(3) * t367;
t612 = pkin(8) ^ 2 - pkin(10) ^ 2;
t389 = t393 + t612;
t394 = -pkin(3) + t629;
t387 = (pkin(3) - t644) * (pkin(3) + t644) + t613;
t388 = (pkin(3) - t643) * (pkin(3) + t643) + t613;
t534 = sqrt(-t388 * t387);
t409 = t428 * t507 - t429 * t508;
t630 = pkin(4) * t409;
t364 = t389 * t630 - t394 * t534;
t631 = pkin(4) * t364;
t484 = t522 * t633 - t525 * t636;
t523 = sin(qJ(1));
t475 = t484 * t523;
t628 = pkin(4) * t475;
t526 = cos(qJ(1));
t477 = t484 * t526;
t627 = pkin(4) * t477;
t566 = t522 * t636 + t525 * t633;
t626 = pkin(4) * t566;
t625 = Icges(2,4) * t523;
t548 = t473 + t611;
t546 = pkin(1) * t548;
t556 = pkin(1) * t561 - pkin(5);
t543 = t651 * (-t535 * t556 + t546 * t562);
t538 = t543 * t609;
t542 = t651 * (-pkin(1) * t558 - t548 * t556);
t540 = t638 * t542;
t427 = t538 * t635 + t540 * t609;
t624 = Icges(7,4) * t427;
t539 = t635 * t542;
t430 = t538 * t638 + t539 * t610;
t623 = Icges(7,4) * t430;
t606 = pkin(4) * t632;
t622 = 0.2e1 * (t387 + t388) * t606 / t534;
t621 = t367 * t534;
t516 = cos(pkin(17));
t620 = t391 * t516;
t528 = 0.1e1 / pkin(10);
t619 = t391 * t528;
t530 = 0.1e1 / pkin(8);
t618 = t391 * t530;
t617 = t409 * t534;
t615 = t470 * t479;
t509 = V_base(6) + qJD(1);
t614 = t526 * t509;
t519 = cos(pkin(19));
t591 = t519 * t610;
t517 = sin(pkin(19));
t592 = t517 * t609;
t425 = (t438 * t591 + t439 * t592) * t533;
t423 = 0.1e1 / t425 ^ 2;
t590 = t519 * t609;
t426 = (t438 * t592 + t439 * t590) * t533;
t355 = ((t420 * t592 + t421 * t590 + t517 * t581 + t519 * t580) / t425 - (t420 * t591 + t421 * t592 + t517 * t580 - t519 * t581) * t426 * t423) / (t423 * t426 ^ 2 + 0.1e1) * t533;
t608 = -qJD(2) - t355;
t607 = -qJD(2) - qJD(3);
t596 = t523 * V_base(4);
t604 = pkin(13) * t596 + V_base(3);
t603 = V_base(5) * pkin(12) + V_base(1);
t600 = pkin(1) * t636;
t599 = t633 * pkin(1);
t598 = Icges(3,4) * t636;
t597 = Icges(3,4) * t633;
t595 = V_base(5) * t526;
t594 = -t622 / 0.2e1;
t593 = t391 * t640;
t505 = qJD(2) * t523 + V_base(4);
t392 = 0.1e1 / t393 ^ 2;
t587 = t392 * t606;
t504 = -qJD(2) * t526 + V_base(5);
t586 = t504 * t599 + t603;
t584 = t651 * t649 * t635;
t351 = t355 * t523 + t505;
t483 = qJD(3) * t523 + t505;
t390 = t393 - t612;
t395 = -pkin(3) * t410 + pkin(4);
t363 = -pkin(3) * t617 + t390 * t395;
t583 = t363 * t587;
t365 = pkin(3) * t390 * t409 + t395 * t534;
t582 = t365 * t587;
t482 = t526 * t607 + V_base(5);
t579 = -t482 * t626 + t586;
t578 = rSges(7,1) * t427 - rSges(7,2) * t430;
t577 = t651 * t638 * t609;
t366 = -t368 * t508 + t369 * t507;
t571 = -t366 * t534 + t409 * t594;
t317 = ((-0.2e1 * pkin(4) * t395 - t390) * t367 + t571) * pkin(3);
t318 = t395 * t622 / 0.2e1 - 0.2e1 * t532 * t367 * t630 + (t366 * t390 - t621) * pkin(3);
t340 = (-t363 * t516 / 0.2e1 + t365 * t640) * t619;
t339 = 0.1e1 / t340 ^ 2;
t341 = (t365 * t516 / 0.2e1 + t363 * t640) * t619;
t272 = ((t318 * t620 / 0.2e1 + t516 * t582 + t317 * t593 + t515 * t583) / t340 - (-t317 * t620 / 0.2e1 - t516 * t583 + t318 * t593 + t515 * t582) * t341 * t339) / (t339 * t341 ^ 2 + 0.1e1) * t528;
t270 = t272 * t523 + t483;
t576 = Icges(7,1) * t427 - t623;
t575 = -Icges(7,2) * t430 + t624;
t574 = Icges(7,5) * t427 - Icges(7,6) * t430;
t573 = (-t600 - pkin(13)) * t523;
t572 = -V_base(4) * pkin(12) + pkin(13) * t614 + V_base(2);
t570 = rSges(3,1) * t636 - rSges(3,2) * t633;
t569 = Icges(3,1) * t636 - t597;
t568 = -Icges(3,2) * t633 + t598;
t567 = Icges(3,5) * t636 - Icges(3,6) * t633;
t411 = t425 * t633 + t426 * t636;
t412 = t425 * t636 - t426 * t633;
t269 = V_base(5) + (-t272 + t607) * t526;
t419 = ((-0.3e1 * t652 + (0.4e1 * t559 - pkin(1)) * pkin(1) - t611) * t479 + t549) * pkin(1);
t422 = pkin(1) * t551 - t480 * t546 + t556 * t653;
t424 = 0.1e1 / t427 ^ 2;
t541 = t543 * t648;
t356 = ((t422 * t577 - t419 * t584 / 0.2e1 + (-t539 * t648 + t541 * t638) * t615) / t427 - (t419 * t577 + t422 * t584 / 0.2e1 + (t540 * t648 + t541 * t635) * t615) * t430 * t424) / (t424 * t430 ^ 2 + 0.1e1);
t353 = -t356 * t526 + V_base(5);
t354 = t356 * t523 + V_base(4);
t565 = (-Icges(7,3) * t526 + t523 * t574) * t353 + (Icges(7,3) * t523 + t526 * t574) * t354 + (Icges(7,5) * t430 + Icges(7,6) * t427) * t509;
t564 = (-Icges(3,3) * t526 + t523 * t567) * t504 + (Icges(3,3) * t523 + t526 * t567) * t505 + (Icges(3,5) * t633 + Icges(3,6) * t636) * t509;
t563 = t573 - t628;
t557 = -t505 * t599 + t600 * t614 + t572;
t552 = t505 * t523 * t600 + (-pkin(13) * V_base(5) - t504 * t600) * t526 + t604;
t550 = t483 * t626 + t509 * t627 + t557;
t547 = -t482 * t627 + t483 * t628 + t552;
t399 = -Icges(7,6) * t526 + t523 * t575;
t400 = Icges(7,6) * t523 + t526 * t575;
t401 = -Icges(7,5) * t526 + t523 * t576;
t402 = Icges(7,5) * t523 + t526 * t576;
t414 = Icges(7,2) * t427 + t623;
t415 = Icges(7,1) * t430 + t624;
t545 = (-t400 * t430 + t402 * t427) * t354 + (-t399 * t430 + t401 * t427) * t353 + (-t414 * t430 + t415 * t427) * t509;
t466 = -Icges(3,6) * t526 + t523 * t568;
t467 = Icges(3,6) * t523 + t526 * t568;
t468 = -Icges(3,5) * t526 + t523 * t569;
t469 = Icges(3,5) * t523 + t526 * t569;
t494 = Icges(3,2) * t636 + t597;
t497 = Icges(3,1) * t633 + t598;
t544 = (-t467 * t633 + t469 * t636) * t505 + (-t466 * t633 + t468 * t636) * t504 + (-t494 * t633 + t497 * t636) * t509;
t524 = cos(qJ(4));
t521 = sin(qJ(4));
t518 = sin(pkin(18));
t512 = Icges(2,4) * t526;
t502 = rSges(2,1) * t526 - rSges(2,2) * t523;
t501 = rSges(2,1) * t523 + rSges(2,2) * t526;
t500 = rSges(3,1) * t633 + rSges(3,2) * t636;
t499 = Icges(2,1) * t526 - t625;
t498 = Icges(2,1) * t523 + t512;
t496 = -Icges(2,2) * t523 + t512;
t495 = Icges(2,2) * t526 + t625;
t490 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t489 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t488 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t476 = t566 * t526;
t474 = t566 * t523;
t472 = t523 * rSges(3,3) + t526 * t570;
t471 = -t526 * rSges(3,3) + t523 * t570;
t461 = V_base(5) * rSges(2,3) - t501 * t509 + t603;
t460 = t502 * t509 + V_base(2) + (-pkin(12) - rSges(2,3)) * V_base(4);
t458 = t501 * V_base(4) - t502 * V_base(5) + V_base(3);
t457 = -rSges(4,1) * t566 + rSges(4,2) * t484;
t456 = -Icges(4,1) * t566 + Icges(4,4) * t484;
t455 = -Icges(4,4) * t566 + Icges(4,2) * t484;
t454 = -Icges(4,5) * t566 + Icges(4,6) * t484;
t452 = rSges(4,1) * t477 + rSges(4,2) * t476 + rSges(4,3) * t523;
t451 = rSges(4,1) * t475 + rSges(4,2) * t474 - rSges(4,3) * t526;
t450 = Icges(4,1) * t477 + Icges(4,4) * t476 + Icges(4,5) * t523;
t449 = Icges(4,1) * t475 + Icges(4,4) * t474 - Icges(4,5) * t526;
t448 = Icges(4,4) * t477 + Icges(4,2) * t476 + Icges(4,6) * t523;
t447 = Icges(4,4) * t475 + Icges(4,2) * t474 - Icges(4,6) * t526;
t446 = Icges(4,5) * t477 + Icges(4,6) * t476 + Icges(4,3) * t523;
t445 = Icges(4,5) * t475 + Icges(4,6) * t474 - Icges(4,3) * t526;
t443 = t500 * t504 + (-pkin(13) * t523 - t471) * t509 + t603;
t442 = t472 * t509 - t500 * t505 + t572;
t441 = -pkin(13) * t595 + t471 * t505 - t472 * t504 + t604;
t433 = t482 * t457 + (-t451 + t573) * t509 + t586;
t432 = t509 * t452 - t483 * t457 + t557;
t431 = t483 * t451 - t482 * t452 + t552;
t416 = rSges(7,1) * t430 + rSges(7,2) * t427;
t408 = t411 * t526;
t407 = t412 * t526;
t406 = t411 * t523;
t405 = t412 * t523;
t404 = rSges(7,3) * t523 + t526 * t578;
t403 = -rSges(7,3) * t526 + t523 * t578;
t386 = rSges(8,1) * t411 + rSges(8,2) * t412;
t385 = Icges(8,1) * t411 + Icges(8,4) * t412;
t384 = Icges(8,4) * t411 + Icges(8,2) * t412;
t383 = Icges(8,5) * t411 + Icges(8,6) * t412;
t382 = (t411 * t520 - t412 * t518) * pkin(3);
t381 = rSges(8,1) * t407 - rSges(8,2) * t408 + rSges(8,3) * t523;
t380 = rSges(8,1) * t405 - rSges(8,2) * t406 - rSges(8,3) * t526;
t379 = Icges(8,1) * t407 - Icges(8,4) * t408 + Icges(8,5) * t523;
t378 = Icges(8,1) * t405 - Icges(8,4) * t406 - Icges(8,5) * t526;
t377 = Icges(8,4) * t407 - Icges(8,2) * t408 + Icges(8,6) * t523;
t376 = Icges(8,4) * t405 - Icges(8,2) * t406 - Icges(8,6) * t526;
t375 = Icges(8,5) * t407 - Icges(8,6) * t408 + Icges(8,3) * t523;
t374 = Icges(8,5) * t405 - Icges(8,6) * t406 - Icges(8,3) * t526;
t373 = (t407 * t520 + t408 * t518) * pkin(3);
t372 = (t405 * t520 + t406 * t518) * pkin(3);
t362 = -pkin(4) * t617 - t389 * t394;
t361 = 0.1e1 / t362 ^ 2;
t350 = t526 * t608 + V_base(5);
t349 = V_base(5) * pkin(14) + t353 * t416 + (pkin(7) * t523 - t403) * t509 + t603;
t348 = -t354 * t416 + V_base(2) + (-pkin(12) - pkin(14)) * V_base(4) + (-pkin(7) * t526 + t404) * t509;
t346 = t350 * t386 + (-t380 + t573) * t509 + t586;
t345 = -t351 * t386 + t509 * t381 + t557;
t343 = (t362 * t639 - t518 * t364 / 0.2e1) * t618;
t342 = (t518 * t362 / 0.2e1 + t364 * t639) * t618;
t338 = -t353 * t404 + t354 * t403 + V_base(3) + (t595 - t596) * pkin(7);
t337 = -t350 * t381 + t351 * t380 + t552;
t335 = t340 * t484 + t341 * t566;
t334 = -t340 * t566 + t341 * t484;
t333 = -qJD(4) * t335 + t509;
t332 = t340 * t476 - t341 * t477;
t331 = t340 * t477 + t341 * t476;
t330 = t340 * t474 - t341 * t475;
t329 = t340 * t475 + t341 * t474;
t328 = t331 * t524 + t521 * t523;
t327 = -t331 * t521 + t523 * t524;
t326 = t329 * t524 - t521 * t526;
t325 = -t329 * t521 - t524 * t526;
t324 = -t342 * t411 + t343 * t412;
t323 = t342 * t412 + t343 * t411;
t322 = -t342 * t407 - t343 * t408;
t321 = -t342 * t408 + t343 * t407;
t320 = -t342 * t405 - t343 * t406;
t319 = -t342 * t406 + t343 * t405;
t316 = pkin(9) * t334 - pkin(11) * t335;
t315 = rSges(5,1) * t334 + rSges(5,2) * t335;
t314 = Icges(5,1) * t334 + Icges(5,4) * t335;
t313 = Icges(5,4) * t334 + Icges(5,2) * t335;
t312 = Icges(5,5) * t334 + Icges(5,6) * t335;
t311 = pkin(9) * t331 - pkin(11) * t332;
t310 = pkin(9) * t329 - pkin(11) * t330;
t309 = rSges(5,1) * t331 + rSges(5,2) * t332 + rSges(5,3) * t523;
t308 = rSges(5,1) * t329 + rSges(5,2) * t330 - rSges(5,3) * t526;
t307 = Icges(5,1) * t331 + Icges(5,4) * t332 + Icges(5,5) * t523;
t306 = Icges(5,1) * t329 + Icges(5,4) * t330 - Icges(5,5) * t526;
t305 = Icges(5,4) * t331 + Icges(5,2) * t332 + Icges(5,6) * t523;
t304 = Icges(5,4) * t329 + Icges(5,2) * t330 - Icges(5,6) * t526;
t303 = Icges(5,5) * t331 + Icges(5,6) * t332 + Icges(5,3) * t523;
t302 = Icges(5,5) * t329 + Icges(5,6) * t330 - Icges(5,3) * t526;
t301 = rSges(9,1) * t323 + rSges(9,2) * t324;
t300 = Icges(9,1) * t323 + Icges(9,4) * t324;
t299 = Icges(9,4) * t323 + Icges(9,2) * t324;
t298 = Icges(9,5) * t323 + Icges(9,6) * t324;
t297 = rSges(9,1) * t321 + rSges(9,2) * t322 + rSges(9,3) * t523;
t296 = rSges(9,1) * t319 + rSges(9,2) * t320 - rSges(9,3) * t526;
t295 = Icges(9,1) * t321 + Icges(9,4) * t322 + Icges(9,5) * t523;
t294 = Icges(9,1) * t319 + Icges(9,4) * t320 - Icges(9,5) * t526;
t293 = Icges(9,4) * t321 + Icges(9,2) * t322 + Icges(9,6) * t523;
t292 = Icges(9,4) * t319 + Icges(9,2) * t320 - Icges(9,6) * t526;
t291 = Icges(9,5) * t321 + Icges(9,6) * t322 + Icges(9,3) * t523;
t290 = Icges(9,5) * t319 + Icges(9,6) * t320 - Icges(9,3) * t526;
t289 = -rSges(6,3) * t335 + (rSges(6,1) * t524 - rSges(6,2) * t521) * t334;
t288 = -Icges(6,5) * t335 + (Icges(6,1) * t524 - Icges(6,4) * t521) * t334;
t287 = -Icges(6,6) * t335 + (Icges(6,4) * t524 - Icges(6,2) * t521) * t334;
t286 = -Icges(6,3) * t335 + (Icges(6,5) * t524 - Icges(6,6) * t521) * t334;
t285 = rSges(6,1) * t328 + rSges(6,2) * t327 - rSges(6,3) * t332;
t284 = rSges(6,1) * t326 + rSges(6,2) * t325 - rSges(6,3) * t330;
t283 = Icges(6,1) * t328 + Icges(6,4) * t327 - Icges(6,5) * t332;
t282 = Icges(6,1) * t326 + Icges(6,4) * t325 - Icges(6,5) * t330;
t281 = Icges(6,4) * t328 + Icges(6,2) * t327 - Icges(6,6) * t332;
t280 = Icges(6,4) * t326 + Icges(6,2) * t325 - Icges(6,6) * t330;
t279 = Icges(6,5) * t328 + Icges(6,6) * t327 - Icges(6,3) * t332;
t278 = Icges(6,5) * t326 + Icges(6,6) * t325 - Icges(6,3) * t330;
t277 = 0.2e1 * (((t394 * t594 + (t366 * t389 - t621) * pkin(4)) * t642 + (-t391 * t409 * t531 + t392 * t631) * t632) / t362 - ((-t367 * t389 + t571) * t642 + (t362 * t392 + t391 * t394) * t632) * t361 * t631) * pkin(8) / (t361 * t364 ^ 2 + 0.1e1) * t393 * t530;
t276 = t277 * t523 + t351;
t275 = V_base(5) + (-t277 + t608) * t526;
t274 = t275 * t301 + t350 * t382 + (-t296 - t372 + t573) * t509 + t586;
t273 = -t276 * t301 - t351 * t382 + (t297 + t373) * t509 + t557;
t268 = -qJD(4) * t332 + t270;
t267 = -qJD(4) * t330 + t269;
t266 = t269 * t315 + (-t308 + t563) * t509 + t579;
t265 = -t270 * t315 + t509 * t309 + t550;
t264 = -t275 * t297 + t276 * t296 - t350 * t373 + t351 * t372 + t552;
t263 = -t269 * t309 + t270 * t308 + t547;
t262 = t267 * t289 + t269 * t316 - t333 * t284 + (-t310 + t563) * t509 + t579;
t261 = -t268 * t289 - t270 * t316 + t333 * t285 + t509 * t311 + t550;
t260 = -t267 * t285 + t268 * t284 - t269 * t311 + t270 * t310 + t547;
t1 = ((-t495 * t523 + t498 * t526 + Icges(1,4)) * V_base(5) + (-t496 * t523 + t499 * t526 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t350 * ((-t375 * t526 - t377 * t406 + t379 * t405) * t351 + (-t374 * t526 - t376 * t406 + t378 * t405) * t350 + (-t383 * t526 - t384 * t406 + t385 * t405) * t509) / 0.2e1 + t351 * ((t375 * t523 - t377 * t408 + t379 * t407) * t351 + (t374 * t523 - t376 * t408 + t378 * t407) * t350 + (t383 * t523 - t384 * t408 + t385 * t407) * t509) / 0.2e1 + ((t467 * t636 + t469 * t633) * t505 + (t466 * t636 + t468 * t633) * t504 + (t448 * t484 - t450 * t566) * t483 + (t447 * t484 - t449 * t566) * t482 + (t400 * t427 + t402 * t430) * t354 + (t399 * t427 + t401 * t430) * t353 + (t377 * t412 + t379 * t411) * t351 + (t376 * t412 + t378 * t411) * t350 + (t305 * t335 + t307 * t334) * t270 + (t304 * t335 + t306 * t334) * t269 + (t293 * t324 + t295 * t323) * t276 + (t292 * t324 + t294 * t323) * t275 + (t299 * t324 + t300 * t323 + t313 * t335 + t314 * t334 + t384 * t412 + t385 * t411 + t414 * t427 + t415 * t430 + t455 * t484 - t456 * t566 + t494 * t636 + t497 * t633 + Icges(2,3)) * t509) * t509 / 0.2e1 + V_base(5) * t509 * (Icges(2,5) * t523 + Icges(2,6) * t526) + t509 * V_base(4) * (Icges(2,5) * t526 - Icges(2,6) * t523) + t268 * ((-t279 * t332 + t281 * t327 + t283 * t328) * t268 + (-t278 * t332 + t280 * t327 + t282 * t328) * t267 + (-t286 * t332 + t287 * t327 + t288 * t328) * t333) / 0.2e1 + t267 * ((-t279 * t330 + t281 * t325 + t283 * t326) * t268 + (-t278 * t330 + t280 * t325 + t282 * t326) * t267 + (-t286 * t330 + t287 * t325 + t288 * t326) * t333) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(8) * (t337 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(3) * (t441 ^ 2 + t442 ^ 2 + t443 ^ 2) / 0.2e1 + t333 * ((-t267 * t278 - t268 * t279 - t286 * t333) * t335 + ((-t281 * t521 + t283 * t524) * t268 + (-t280 * t521 + t282 * t524) * t267 + (-t287 * t521 + t288 * t524) * t333) * t334) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + m(4) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(2) * (t458 ^ 2 + t460 ^ 2 + t461 ^ 2) / 0.2e1 + ((t495 * t526 + t498 * t523 + Icges(1,2)) * V_base(5) + (t496 * t526 + t499 * t523 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t482 * ((-t446 * t526 + t448 * t474 + t450 * t475) * t483 + (-t445 * t526 + t447 * t474 + t449 * t475) * t482 + (-t454 * t526 + t455 * t474 + t456 * t475) * t509) / 0.2e1 + t275 * ((-t291 * t526 + t293 * t320 + t295 * t319) * t276 + (-t290 * t526 + t292 * t320 + t294 * t319) * t275 + (-t298 * t526 + t299 * t320 + t300 * t319) * t509) / 0.2e1 + t269 * ((-t303 * t526 + t305 * t330 + t307 * t329) * t270 + (-t302 * t526 + t304 * t330 + t306 * t329) * t269 + (-t312 * t526 + t313 * t330 + t314 * t329) * t509) / 0.2e1 + t505 * (t523 * t564 + t526 * t544) / 0.2e1 + t504 * (t523 * t544 - t526 * t564) / 0.2e1 + t354 * (t523 * t565 + t526 * t545) / 0.2e1 + t353 * (t523 * t545 - t526 * t565) / 0.2e1 + t483 * ((t446 * t523 + t448 * t476 + t450 * t477) * t483 + (t445 * t523 + t447 * t476 + t449 * t477) * t482 + (t454 * t523 + t455 * t476 + t456 * t477) * t509) / 0.2e1 + t276 * ((t291 * t523 + t293 * t322 + t295 * t321) * t276 + (t290 * t523 + t292 * t322 + t294 * t321) * t275 + (t298 * t523 + t299 * t322 + t300 * t321) * t509) / 0.2e1 + t270 * ((t303 * t523 + t305 * t332 + t307 * t331) * t270 + (t302 * t523 + t304 * t332 + t306 * t331) * t269 + (t312 * t523 + t313 * t332 + t314 * t331) * t509) / 0.2e1 + m(5) * (t263 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(1) * (t488 ^ 2 + t489 ^ 2 + t490 ^ 2) / 0.2e1 + m(9) * (t264 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(7) * (t338 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1;
T = t1;
