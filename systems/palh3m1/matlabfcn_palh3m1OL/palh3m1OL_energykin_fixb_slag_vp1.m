% Calculate kinetic energy for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1OL_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:03:52
% EndTime: 2020-04-20 17:03:57
% DurationCPUTime: 4.73s
% Computational Cost: add. (1817->371), mult. (1987->608), div. (0->0), fcn. (1782->18), ass. (0->223)
t531 = sin(qJ(1));
t618 = t531 ^ 2;
t535 = cos(qJ(1));
t617 = t535 ^ 2;
t526 = qJ(2) + qJ(7);
t523 = pkin(15) - t526;
t615 = pkin(3) * sin(t523);
t527 = qJ(2) + qJ(3);
t520 = sin(t527);
t614 = pkin(4) * t520;
t612 = pkin(12) * t531;
t534 = cos(qJ(2));
t611 = t534 * pkin(1);
t530 = sin(qJ(2));
t610 = Icges(3,4) * t530;
t609 = Icges(3,4) * t534;
t608 = Icges(4,4) * t520;
t522 = cos(t527);
t607 = Icges(4,4) * t522;
t524 = qJ(4) + t527;
t511 = sin(t524);
t606 = Icges(5,4) * t511;
t512 = cos(t524);
t605 = Icges(5,4) * t512;
t528 = sin(qJ(6));
t604 = Icges(7,4) * t528;
t532 = cos(qJ(6));
t603 = Icges(7,4) * t532;
t519 = sin(t526);
t602 = Icges(8,4) * t519;
t521 = cos(t526);
t601 = Icges(8,4) * t521;
t513 = -qJ(8) + t523;
t507 = sin(t513);
t600 = Icges(9,4) * t507;
t508 = cos(t513);
t599 = Icges(9,4) * t508;
t598 = t511 * t531;
t597 = t511 * t535;
t529 = sin(qJ(5));
t596 = t529 * t531;
t595 = t529 * t535;
t533 = cos(qJ(5));
t594 = t531 * t533;
t593 = t533 * t535;
t480 = t611 * t531;
t481 = t611 * t535;
t518 = qJD(2) * t531;
t589 = qJD(2) * t535;
t592 = t480 * t518 + t481 * t589;
t591 = pkin(3) * cos(t523);
t590 = pkin(4) * t522;
t492 = qJD(7) * t531 + t518;
t493 = qJD(3) * t531 + t518;
t588 = qJD(5) * t511;
t587 = qJD(6) * t531;
t586 = qJD(6) * t535;
t585 = -qJD(2) - qJD(3);
t584 = -qJD(2) - qJD(7);
t583 = pkin(1) * qJD(2) * t530;
t477 = qJD(4) * t531 + t493;
t582 = -t480 - t612;
t581 = t535 * t583;
t423 = t590 * t531;
t580 = t423 + t582;
t579 = -pkin(8) * t512 - pkin(10) * t511;
t578 = rSges(3,1) * t534 - rSges(3,2) * t530;
t577 = -rSges(4,1) * t522 + rSges(4,2) * t520;
t576 = -rSges(5,1) * t512 + rSges(5,2) * t511;
t575 = rSges(7,1) * t532 - rSges(7,2) * t528;
t574 = rSges(8,1) * t521 - rSges(8,2) * t519;
t573 = -rSges(9,1) * t508 - rSges(9,2) * t507;
t479 = (-qJD(4) + t585) * t535;
t572 = Icges(3,1) * t534 - t610;
t571 = -Icges(4,1) * t522 + t608;
t570 = -Icges(5,1) * t512 + t606;
t569 = Icges(7,1) * t532 - t604;
t568 = Icges(8,1) * t521 - t602;
t567 = -Icges(9,1) * t508 - t600;
t566 = -Icges(3,2) * t530 + t609;
t565 = Icges(4,2) * t520 - t607;
t564 = Icges(5,2) * t511 - t605;
t563 = -Icges(7,2) * t528 + t603;
t562 = -Icges(8,2) * t519 + t601;
t561 = -Icges(9,2) * t507 - t599;
t560 = Icges(3,5) * t534 - Icges(3,6) * t530;
t559 = -Icges(4,5) * t522 + Icges(4,6) * t520;
t558 = -Icges(5,5) * t512 + Icges(5,6) * t511;
t557 = Icges(7,5) * t532 - Icges(7,6) * t528;
t556 = Icges(8,5) * t521 - Icges(8,6) * t519;
t555 = -Icges(9,5) * t508 - Icges(9,6) * t507;
t447 = -Icges(7,6) * t535 + t531 * t563;
t451 = -Icges(7,5) * t535 + t531 * t569;
t554 = t447 * t528 - t451 * t532;
t448 = Icges(7,6) * t531 + t535 * t563;
t452 = Icges(7,5) * t531 + t535 * t569;
t553 = -t448 * t528 + t452 * t532;
t449 = -Icges(3,6) * t535 + t531 * t566;
t453 = -Icges(3,5) * t535 + t531 * t572;
t552 = t449 * t530 - t453 * t534;
t450 = Icges(3,6) * t531 + t535 * t566;
t454 = Icges(3,5) * t531 + t535 * t572;
t551 = -t450 * t530 + t454 * t534;
t499 = Icges(7,2) * t532 + t604;
t501 = Icges(7,1) * t528 + t603;
t550 = -t499 * t528 + t501 * t532;
t500 = Icges(3,2) * t534 + t610;
t502 = Icges(3,1) * t530 + t609;
t549 = -t500 * t530 + t502 * t534;
t424 = t590 * t535;
t495 = t585 * t535;
t548 = -t493 * t423 + t424 * t495 + t592;
t514 = qJD(1) * t535 * pkin(12);
t547 = qJD(1) * t481 - t531 * t583 + t514;
t546 = -t495 * t614 - t581;
t476 = qJD(8) * t531 + t492;
t478 = (-qJD(8) + t584) * t535;
t545 = (Icges(9,5) * t507 - Icges(9,6) * t508) * qJD(1) + (-Icges(9,3) * t535 + t531 * t555) * t478 + (Icges(9,3) * t531 + t535 * t555) * t476;
t544 = (-Icges(5,5) * t511 - Icges(5,6) * t512) * qJD(1) + (-Icges(5,3) * t535 + t531 * t558) * t479 + (Icges(5,3) * t531 + t535 * t558) * t477;
t494 = t584 * t535;
t543 = (Icges(8,5) * t519 + Icges(8,6) * t521) * qJD(1) + (-Icges(8,3) * t535 + t531 * t556) * t494 + (Icges(8,3) * t531 + t535 * t556) * t492;
t542 = (-Icges(4,5) * t520 - Icges(4,6) * t522) * qJD(1) + (-Icges(4,3) * t535 + t531 * t559) * t495 + (Icges(4,3) * t531 + t535 * t559) * t493;
t541 = -qJD(1) * t424 + t493 * t614 + t547;
t403 = -Icges(9,6) * t535 + t531 * t561;
t404 = Icges(9,6) * t531 + t535 * t561;
t405 = -Icges(9,5) * t535 + t531 * t567;
t406 = Icges(9,5) * t531 + t535 * t567;
t462 = -Icges(9,2) * t508 + t600;
t463 = Icges(9,1) * t507 - t599;
t540 = (-t404 * t507 - t406 * t508) * t476 + (-t403 * t507 - t405 * t508) * t478 + (-t462 * t507 - t463 * t508) * qJD(1);
t415 = -Icges(5,6) * t535 + t531 * t564;
t416 = Icges(5,6) * t531 + t535 * t564;
t417 = -Icges(5,5) * t535 + t531 * t570;
t418 = Icges(5,5) * t531 + t535 * t570;
t472 = -Icges(5,2) * t512 - t606;
t473 = -Icges(5,1) * t511 - t605;
t539 = (t416 * t511 - t418 * t512) * t477 + (t415 * t511 - t417 * t512) * t479 + (t472 * t511 - t473 * t512) * qJD(1);
t429 = -Icges(8,6) * t535 + t531 * t562;
t430 = Icges(8,6) * t531 + t535 * t562;
t433 = -Icges(8,5) * t535 + t531 * t568;
t434 = Icges(8,5) * t531 + t535 * t568;
t484 = Icges(8,2) * t521 + t602;
t486 = Icges(8,1) * t519 + t601;
t538 = (-t430 * t519 + t434 * t521) * t492 + (-t429 * t519 + t433 * t521) * t494 + (-t484 * t519 + t486 * t521) * qJD(1);
t431 = -Icges(4,6) * t535 + t531 * t565;
t432 = Icges(4,6) * t531 + t535 * t565;
t435 = -Icges(4,5) * t535 + t531 * t571;
t436 = Icges(4,5) * t531 + t535 * t571;
t485 = -Icges(4,2) * t522 - t608;
t487 = -Icges(4,1) * t520 - t607;
t537 = (t432 * t520 - t436 * t522) * t493 + (t431 * t520 - t435 * t522) * t495 + (t485 * t520 - t487 * t522) * qJD(1);
t506 = rSges(2,1) * t535 - rSges(2,2) * t531;
t505 = rSges(2,1) * t531 + rSges(2,2) * t535;
t504 = rSges(3,1) * t530 + rSges(3,2) * t534;
t503 = rSges(7,1) * t528 + rSges(7,2) * t532;
t498 = Icges(3,5) * t530 + Icges(3,6) * t534;
t497 = Icges(7,5) * t528 + Icges(7,6) * t532;
t496 = qJD(5) * t512 + qJD(1);
t490 = -rSges(4,1) * t520 - rSges(4,2) * t522;
t489 = rSges(8,1) * t519 + rSges(8,2) * t521;
t475 = -pkin(8) * t511 + pkin(10) * t512;
t474 = -rSges(5,1) * t511 - rSges(5,2) * t512;
t468 = -t512 * t593 + t596;
t467 = t512 * t595 + t594;
t466 = -t512 * t594 - t595;
t465 = t512 * t596 - t593;
t464 = rSges(9,1) * t507 - rSges(9,2) * t508;
t458 = rSges(3,3) * t531 + t535 * t578;
t457 = rSges(7,3) * t531 + t535 * t575;
t456 = -rSges(3,3) * t535 + t531 * t578;
t455 = -rSges(7,3) * t535 + t531 * t575;
t446 = Icges(3,3) * t531 + t535 * t560;
t445 = -Icges(3,3) * t535 + t531 * t560;
t444 = Icges(7,3) * t531 + t535 * t557;
t443 = -Icges(7,3) * t535 + t531 * t557;
t442 = t579 * t535;
t441 = t579 * t531;
t440 = rSges(4,3) * t531 + t535 * t577;
t439 = rSges(8,3) * t531 + t535 * t574;
t438 = -rSges(4,3) * t535 + t531 * t577;
t437 = -rSges(8,3) * t535 + t531 * t574;
t422 = rSges(5,3) * t531 + t535 * t576;
t421 = -rSges(5,3) * t535 + t531 * t576;
t420 = -t531 * t588 + t479;
t419 = -t535 * t588 + t477;
t412 = t591 * t535;
t411 = t591 * t531;
t409 = rSges(9,3) * t531 + t535 * t573;
t408 = -rSges(9,3) * t535 + t531 * t573;
t407 = rSges(6,3) * t512 + (-rSges(6,1) * t533 + rSges(6,2) * t529) * t511;
t400 = Icges(6,5) * t512 + (-Icges(6,1) * t533 + Icges(6,4) * t529) * t511;
t399 = Icges(6,6) * t512 + (-Icges(6,4) * t533 + Icges(6,2) * t529) * t511;
t398 = Icges(6,3) * t512 + (-Icges(6,5) * t533 + Icges(6,6) * t529) * t511;
t396 = -t503 * t587 + (-pkin(6) * t535 + t457) * qJD(1);
t395 = qJD(1) * t458 - t504 * t518 + t514;
t394 = -t503 * t586 + (pkin(6) * t531 - t455) * qJD(1);
t393 = -t504 * t589 + (-t456 - t612) * qJD(1);
t392 = (t456 * t531 + t458 * t535) * qJD(2);
t391 = (t455 * t531 + t457 * t535) * qJD(6);
t390 = rSges(6,1) * t468 + rSges(6,2) * t467 - rSges(6,3) * t597;
t389 = rSges(6,1) * t466 + rSges(6,2) * t465 - rSges(6,3) * t598;
t388 = Icges(6,1) * t468 + Icges(6,4) * t467 - Icges(6,5) * t597;
t387 = Icges(6,1) * t466 + Icges(6,4) * t465 - Icges(6,5) * t598;
t386 = Icges(6,4) * t468 + Icges(6,2) * t467 - Icges(6,6) * t597;
t385 = Icges(6,4) * t466 + Icges(6,2) * t465 - Icges(6,6) * t598;
t384 = Icges(6,5) * t468 + Icges(6,6) * t467 - Icges(6,3) * t597;
t383 = Icges(6,5) * t466 + Icges(6,6) * t465 - Icges(6,3) * t598;
t382 = qJD(1) * t440 - t490 * t493 + t547;
t381 = qJD(1) * t439 - t489 * t492 + t547;
t380 = -t581 + t490 * t495 + (-t438 + t582) * qJD(1);
t379 = -t581 + t489 * t494 + (-t437 + t582) * qJD(1);
t378 = t438 * t493 - t440 * t495 + t592;
t377 = t437 * t492 - t439 * t494 + t592;
t376 = qJD(1) * t422 - t474 * t477 + t541;
t375 = t474 * t479 + (-t421 + t580) * qJD(1) + t546;
t374 = t492 * t615 - t464 * t476 + (t409 + t412) * qJD(1) + t547;
t373 = -t581 - t494 * t615 + t464 * t478 + (-t408 - t411 + t582) * qJD(1);
t372 = t421 * t477 - t422 * t479 + t548;
t371 = t408 * t476 - t409 * t478 + t411 * t492 - t412 * t494 + t592;
t370 = qJD(1) * t442 + t390 * t496 - t419 * t407 - t475 * t477 + t541;
t369 = -t389 * t496 + t420 * t407 + t475 * t479 + (-t441 + t580) * qJD(1) + t546;
t368 = t419 * t389 - t420 * t390 + t441 * t477 - t442 * t479 + t548;
t1 = m(5) * (t372 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(8) * (t377 ^ 2 + t379 ^ 2 + t381 ^ 2) / 0.2e1 + m(4) * (t378 ^ 2 + t380 ^ 2 + t382 ^ 2) / 0.2e1 + m(6) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(9) * (t371 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + t496 * ((t383 * t420 + t384 * t419 + t398 * t496) * t512 + ((t386 * t529 - t388 * t533) * t419 + (t385 * t529 - t387 * t533) * t420 + (t399 * t529 - t400 * t533) * t496) * t511) / 0.2e1 + t493 * (t531 * t542 + t535 * t537) / 0.2e1 + t495 * (t531 * t537 - t535 * t542) / 0.2e1 + t492 * (t531 * t543 + t535 * t538) / 0.2e1 + t494 * (t531 * t538 - t535 * t543) / 0.2e1 + t477 * (t531 * t544 + t539 * t535) / 0.2e1 + t479 * (t539 * t531 - t535 * t544) / 0.2e1 + t476 * (t545 * t531 + t540 * t535) / 0.2e1 + t478 * (t540 * t531 - t545 * t535) / 0.2e1 + m(3) * (t392 ^ 2 + t393 ^ 2 + t395 ^ 2) / 0.2e1 + m(7) * (t391 ^ 2 + t394 ^ 2 + t396 ^ 2) / 0.2e1 - ((-t535 * t498 + t531 * t549) * qJD(1) + (t617 * t445 + (t551 * t531 + (-t446 + t552) * t535) * t531) * qJD(2)) * t589 / 0.2e1 - ((-t535 * t497 + t531 * t550) * qJD(1) + (t617 * t443 + (t553 * t531 + (-t444 + t554) * t535) * t531) * qJD(6)) * t586 / 0.2e1 + ((t531 * t498 + t535 * t549) * qJD(1) + (t618 * t446 + (t552 * t535 + (-t445 + t551) * t531) * t535) * qJD(2)) * t518 / 0.2e1 + ((t531 * t497 + t535 * t550) * qJD(1) + (t618 * t444 + (t554 * t535 + (-t443 + t553) * t531) * t535) * qJD(6)) * t587 / 0.2e1 + t419 * ((-t384 * t597 + t467 * t386 + t468 * t388) * t419 + (-t383 * t597 + t385 * t467 + t387 * t468) * t420 + (-t398 * t597 + t399 * t467 + t400 * t468) * t496) / 0.2e1 + t420 * ((-t384 * t598 + t386 * t465 + t388 * t466) * t419 + (-t383 * t598 + t465 * t385 + t466 * t387) * t420 + (-t398 * t598 + t399 * t465 + t400 * t466) * t496) / 0.2e1 + (Icges(2,3) + m(2) * (t505 ^ 2 + t506 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t450 * t534 + t454 * t530) * t531 - (t449 * t534 + t453 * t530) * t535) * qJD(2) + ((t448 * t532 + t452 * t528) * t531 - (t447 * t532 + t451 * t528) * t535) * qJD(6) + (-t432 * t522 - t436 * t520) * t493 + (-t431 * t522 - t435 * t520) * t495 + (t430 * t521 + t434 * t519) * t492 + (t429 * t521 + t433 * t519) * t494 + (-t416 * t512 - t418 * t511) * t477 + (-t415 * t512 - t417 * t511) * t479 + (-t404 * t508 + t406 * t507) * t476 + (-t403 * t508 + t405 * t507) * t478 + (-t508 * t462 + t507 * t463 - t512 * t472 - t511 * t473 + t521 * t484 - t522 * t485 + t519 * t486 - t520 * t487 + t499 * t532 + t534 * t500 + t528 * t501 + t530 * t502) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
