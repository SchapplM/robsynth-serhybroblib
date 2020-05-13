% Calculate kinetic energy for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2TE_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:57
% EndTime: 2020-05-07 01:42:02
% DurationCPUTime: 4.92s
% Computational Cost: add. (1766->297), mult. (3101->519), div. (0->0), fcn. (3659->22), ass. (0->180)
t611 = sin(qJ(1));
t687 = t611 ^ 2;
t617 = cos(qJ(1));
t686 = t617 ^ 2;
t700 = Icges(3,3) + Icges(7,3);
t612 = sin(pkin(15));
t613 = sin(pkin(14));
t618 = cos(pkin(15));
t619 = cos(pkin(14));
t569 = t612 * t619 - t613 * t618;
t570 = t612 * t613 + t618 * t619;
t610 = sin(qJ(2));
t616 = cos(qJ(2));
t528 = t569 * t616 + t570 * t610;
t529 = -t569 * t610 + t570 * t616;
t699 = Icges(3,5) * t616 + Icges(7,5) * t529 - Icges(3,6) * t610 - Icges(7,6) * t528;
t602 = qJ(3) + qJ(2);
t597 = sin(t602);
t598 = cos(t602);
t698 = Icges(9,5) * t597 + Icges(9,6) * t598;
t695 = Icges(7,4) * t528;
t600 = qJD(2) + qJD(3);
t572 = t600 * t611;
t573 = t600 * t617;
t694 = t611 * t617;
t693 = t611 * t699 - t617 * t700;
t692 = Icges(3,5) * t610 + Icges(7,5) * t528 + Icges(3,6) * t616 + Icges(7,6) * t529;
t691 = t611 * t700 + t617 * t699;
t496 = Icges(7,2) * t529 + t695;
t677 = Icges(7,4) * t529;
t497 = Icges(7,1) * t528 + t677;
t681 = Icges(3,4) * t610;
t577 = Icges(3,2) * t616 + t681;
t680 = Icges(3,4) * t616;
t578 = Icges(3,1) * t610 + t680;
t690 = -t496 * t528 + t497 * t529 - t577 * t610 + t578 * t616;
t649 = -Icges(7,2) * t528 + t677;
t491 = -Icges(7,6) * t617 + t611 * t649;
t653 = Icges(7,1) * t529 - t695;
t493 = -Icges(7,5) * t617 + t611 * t653;
t651 = -Icges(3,2) * t610 + t680;
t550 = -Icges(3,6) * t617 + t611 * t651;
t655 = Icges(3,1) * t616 - t681;
t552 = -Icges(3,5) * t617 + t611 * t655;
t689 = t491 * t528 - t493 * t529 + t550 * t610 - t552 * t616;
t492 = Icges(7,6) * t611 + t617 * t649;
t494 = Icges(7,5) * t611 + t617 * t653;
t551 = Icges(3,6) * t611 + t617 * t651;
t553 = Icges(3,5) * t611 + t617 * t655;
t688 = -t492 * t528 + t494 * t529 - t551 * t610 + t553 * t616;
t666 = t617 * qJD(1);
t683 = pkin(1) * qJD(2);
t682 = pkin(4) * qJD(3);
t615 = cos(qJ(3));
t609 = sin(qJ(3));
t672 = t609 * t610;
t567 = -t615 * t616 + t672;
t679 = Icges(4,4) * t567;
t671 = t609 * t616;
t568 = t610 * t615 + t671;
t678 = Icges(4,4) * t568;
t603 = sin(pkin(16));
t604 = cos(pkin(16));
t563 = t603 * t618 + t604 * t612;
t564 = -t603 * t612 + t604 * t618;
t601 = pkin(17) + pkin(18);
t595 = sin(t601);
t596 = cos(t601);
t524 = t563 * t596 + t564 * t595;
t674 = t524 * t611;
t673 = t524 * t617;
t593 = -pkin(4) * t615 + pkin(1);
t670 = pkin(4) * t672 + t593 * t616;
t669 = qJD(2) * t611;
t668 = qJD(2) * t617;
t667 = qJD(4) * t524;
t665 = t610 * t683;
t664 = t568 * t682;
t608 = sin(qJ(4));
t614 = cos(qJ(4));
t657 = rSges(6,1) * t614 - rSges(6,2) * t608;
t510 = rSges(6,3) * t564 + t563 * t657;
t511 = -rSges(6,3) * t563 + t564 * t657;
t663 = (t510 * t596 + t511 * t595) * t667;
t591 = pkin(1) * t616 + pkin(12);
t659 = t617 * t665;
t525 = qJD(2) * t670 + t567 * t682;
t558 = -pkin(4) * t671 + t593 * t610;
t658 = -t558 * t668 + t617 * t664;
t585 = rSges(3,1) * t616 - rSges(3,2) * t610;
t656 = -rSges(9,1) * t598 + rSges(9,2) * t597;
t654 = Icges(4,1) * t567 + t678;
t650 = Icges(4,2) * t568 + t679;
t646 = Icges(4,5) * t567 + Icges(4,6) * t568;
t644 = -Icges(9,5) * t598 + Icges(9,6) * t597;
t522 = t563 * t595 - t564 * t596;
t626 = t608 * t522;
t506 = -t611 * t626 - t614 * t617;
t624 = t522 * t614;
t507 = -t608 * t617 + t611 * t624;
t508 = t611 * t614 - t617 * t626;
t509 = t608 * t611 + t617 * t624;
t643 = (Icges(6,5) * t507 + Icges(6,6) * t506 + Icges(6,3) * t674) * t611 + (Icges(6,5) * t509 + Icges(6,6) * t508 + Icges(6,3) * t673) * t617;
t639 = t510 * t595 - t511 * t596;
t574 = rSges(5,1) * t603 - rSges(5,2) * t604;
t575 = rSges(5,1) * t604 + rSges(5,2) * t603;
t638 = -(t574 * t618 + t575 * t612) * t595 + (-t574 * t612 + t575 * t618) * t596;
t583 = -rSges(7,1) * t618 + rSges(7,2) * t612;
t587 = rSges(7,1) * t612 + rSges(7,2) * t618;
t631 = t583 * t613 + t587 * t619;
t605 = sin(pkin(18));
t607 = cos(pkin(18));
t565 = t605 * t618 + t607 * t612;
t630 = t591 * t666 - t611 * t665;
t629 = pkin(12) + t585;
t557 = pkin(12) + t670;
t628 = t557 * t666 - t558 * t669 + t611 * t664;
t627 = -(-rSges(9,1) * t597 - rSges(9,2) * t598) * t600 - t665;
t625 = t643 * t524;
t623 = (-Icges(4,5) * t568 + Icges(4,6) * t567) * qJD(1) - (-Icges(4,3) * t617 + t611 * t646) * t573 + (Icges(4,3) * t611 + t617 * t646) * t572;
t518 = -Icges(4,6) * t617 + t611 * t650;
t519 = Icges(4,6) * t611 + t617 * t650;
t520 = -Icges(4,5) * t617 + t611 * t654;
t521 = Icges(4,5) * t611 + t617 * t654;
t532 = Icges(4,2) * t567 - t678;
t533 = -Icges(4,1) * t568 + t679;
t622 = (t519 * t568 + t521 * t567) * t572 - (t518 * t568 + t520 * t567) * t573 + (t532 * t568 + t533 * t567) * qJD(1);
t620 = qJD(2) ^ 2;
t594 = t616 * t683;
t586 = rSges(2,1) * t617 - rSges(2,2) * t611;
t584 = rSges(6,1) * t608 + rSges(6,2) * t614;
t582 = rSges(2,1) * t611 + rSges(2,2) * t617;
t581 = rSges(3,1) * t610 + rSges(3,2) * t616;
t580 = pkin(8) * t604 - pkin(10) * t603;
t579 = pkin(8) * t603 + pkin(10) * t604;
t566 = -t605 * t612 + t607 * t618;
t546 = rSges(9,3) * t611 + t617 * t656;
t545 = -rSges(9,3) * t617 + t611 * t656;
t540 = Icges(9,3) * t611 + t617 * t644;
t539 = -Icges(9,3) * t617 + t611 * t644;
t538 = -t583 * t619 + t587 * t613;
t537 = (-rSges(4,1) * t616 + rSges(4,2) * t610) * t609 - t615 * (rSges(4,1) * t610 + rSges(4,2) * t616);
t536 = (-rSges(4,1) * t615 + rSges(4,2) * t609) * t616 + t610 * (rSges(4,1) * t609 + rSges(4,2) * t615);
t527 = rSges(4,3) * t611 + t536 * t617;
t526 = -rSges(4,3) * t617 + t536 * t611;
t515 = -t581 * t669 + (rSges(3,3) * t611 + t617 * t629) * qJD(1);
t514 = -t581 * t668 + (rSges(3,3) * t617 - t611 * t629) * qJD(1);
t513 = -qJD(4) * t522 + qJD(1);
t512 = (cos(pkin(17)) * t566 - t565 * sin(pkin(17))) * pkin(3) - t591;
t505 = qJD(1) * (((-rSges(8,1) * t618 - rSges(8,2) * t612) * t607 + t605 * (rSges(8,1) * t612 - rSges(8,2) * t618)) * t617 + t611 * rSges(8,3)) + t630;
t504 = -t659 + (t617 * rSges(8,3) + (-t591 + (rSges(8,1) * t607 + rSges(8,2) * t605) * t618 - t612 * (rSges(8,1) * t605 - rSges(8,2) * t607)) * t611) * qJD(1);
t503 = t538 * t616 - t610 * t631;
t502 = t538 * t610 + t616 * t631;
t501 = (-t579 * t612 + t580 * t618) * t596 - (t579 * t618 + t580 * t612) * t595;
t500 = t594 + (t545 * t611 + t546 * t617) * t600;
t499 = rSges(7,3) * t611 + t503 * t617;
t498 = -rSges(7,3) * t617 + t503 * t611;
t488 = qJD(1) * t527 - t537 * t572 + t630;
t487 = -t659 - t537 * t573 + (-t611 * t591 - t526) * qJD(1);
t486 = t526 * t572 + t527 * t573 + t594;
t485 = t627 * t611 + (-t512 * t617 + t546) * qJD(1);
t484 = t627 * t617 + (t512 * t611 - t545) * qJD(1);
t482 = -Icges(6,5) * t522 + (Icges(6,1) * t614 - Icges(6,4) * t608) * t524;
t481 = -Icges(6,6) * t522 + (Icges(6,4) * t614 - Icges(6,2) * t608) * t524;
t480 = -Icges(6,3) * t522 + (Icges(6,5) * t614 - Icges(6,6) * t608) * t524;
t479 = -t584 * t617 + t611 * t639;
t478 = t584 * t611 + t617 * t639;
t477 = -t502 * t669 + (-pkin(6) * t617 + t499) * qJD(1);
t476 = -t502 * t668 + (pkin(6) * t611 - t498) * qJD(1);
t475 = -qJD(1) * (-rSges(5,3) * t611 + t617 * t638) + t628;
t474 = (t617 * rSges(5,3) + (-t557 + t638) * t611) * qJD(1) + t658;
t473 = (t498 * t611 + t499 * t617) * qJD(2);
t472 = Icges(6,1) * t509 + Icges(6,4) * t508 + Icges(6,5) * t673;
t471 = Icges(6,1) * t507 + Icges(6,4) * t506 + Icges(6,5) * t674;
t470 = Icges(6,4) * t509 + Icges(6,2) * t508 + Icges(6,6) * t673;
t469 = Icges(6,4) * t507 + Icges(6,2) * t506 + Icges(6,6) * t674;
t466 = (-t478 * t611 + t479 * t617) * t667 + t525;
t465 = t513 * t478 + (-qJD(1) * t501 - t663) * t617 + t628;
t464 = -t513 * t479 + (t663 + (t501 - t557) * qJD(1)) * t611 + t658;
t1 = m(6) * (t464 ^ 2 + t465 ^ 2 + t466 ^ 2) / 0.2e1 + m(7) * (t473 ^ 2 + t476 ^ 2 + t477 ^ 2) / 0.2e1 + m(4) * (t486 ^ 2 + t487 ^ 2 + t488 ^ 2) / 0.2e1 + m(5) * (t474 ^ 2 + t475 ^ 2 + t525 ^ 2) / 0.2e1 + m(9) * (t484 ^ 2 + t485 ^ 2 + t500 ^ 2) / 0.2e1 + t513 * ((-t522 * t480 + (-t481 * t608 + t482 * t614) * t524) * t513 + (((-t470 * t608 + t472 * t614) * t617 + (-t469 * t608 + t471 * t614) * t611) * t524 - t643 * t522) * t667) / 0.2e1 + m(3) * (t585 ^ 2 * t620 + t514 ^ 2 + t515 ^ 2) / 0.2e1 + m(8) * (pkin(1) ^ 2 * t616 ^ 2 * t620 + t504 ^ 2 + t505 ^ 2) / 0.2e1 + ((t691 * t687 + (t689 * t617 + (t688 - t693) * t611) * t617) * qJD(2) + (t611 * t692 + t617 * t690) * qJD(1)) * t669 / 0.2e1 - ((t693 * t686 + (t688 * t611 + (t689 - t691) * t617) * t611) * qJD(2) + (t611 * t690 - t617 * t692) * qJD(1)) * t668 / 0.2e1 + (t611 * ((t480 * t674 + t481 * t506 + t482 * t507) * t513 + ((t470 * t506 + t472 * t507) * t617 + (t469 * t506 + t471 * t507 + t625) * t611) * t667) + t617 * ((t480 * t673 + t481 * t508 + t482 * t509) * t513 + ((t469 * t508 + t471 * t509) * t611 + (t470 * t508 + t472 * t509 + t625) * t617) * t667)) * t667 / 0.2e1 + (m(2) * (t582 ^ 2 + t586 ^ 2) + t566 ^ 2 * Icges(8,2) + (Icges(8,1) * t565 - 0.2e1 * Icges(8,4) * t566) * t565 + Icges(5,1) * t524 ^ 2 + (0.2e1 * Icges(5,4) * t524 + Icges(5,2) * t522) * t522 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t519 * t567 - t521 * t568) * t572 - (t518 * t567 - t520 * t568) * t573 + ((-t491 * t529 - t493 * t528 - t550 * t616 - t552 * t610) * t617 + (t492 * t529 + t494 * t528 + t551 * t616 + t553 * t610) * t611) * qJD(2) + (t598 ^ 2 * Icges(9,2) + t529 * t496 + t528 * t497 + t567 * t532 - t568 * t533 + t616 * t577 + t610 * t578 + (Icges(9,1) * t597 + 0.2e1 * Icges(9,4) * t598) * t597) * qJD(1) + (-t686 - t687) * t600 * t698) * qJD(1) / 0.2e1 + ((-t539 * t694 + t687 * t540) * t600 + t622 * t617 + (-qJD(1) * t698 + t623) * t611) * t572 / 0.2e1 - (t622 * t611 - t623 * t617 + t698 * t666 + (t686 * t539 - t540 * t694) * t600) * t573 / 0.2e1;
T = t1;
