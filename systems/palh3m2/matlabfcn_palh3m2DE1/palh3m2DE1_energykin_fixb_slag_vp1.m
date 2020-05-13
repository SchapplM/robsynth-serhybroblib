% Calculate kinetic energy for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE1_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:50
% EndTime: 2020-05-07 01:57:56
% DurationCPUTime: 5.53s
% Computational Cost: add. (1766->297), mult. (3101->516), div. (0->0), fcn. (3659->22), ass. (0->183)
t623 = cos(qJ(3));
t624 = cos(qJ(2));
t617 = sin(qJ(3));
t618 = sin(qJ(2));
t680 = t617 * t618;
t575 = -t623 * t624 + t680;
t679 = t617 * t624;
t576 = t618 * t623 + t679;
t608 = qJD(2) + qJD(3);
t619 = sin(qJ(1));
t580 = t608 * t619;
t625 = cos(qJ(1));
t581 = t608 * t625;
t653 = Icges(4,5) * t575 + Icges(4,6) * t576;
t610 = qJ(3) + qJ(2);
t605 = sin(t610);
t606 = cos(t610);
t709 = Icges(9,5) * t605 + Icges(9,6) * t606;
t712 = -(-Icges(4,3) * t625 + t619 * t653) * t581 + (Icges(4,3) * t619 + t625 * t653) * t580 + (-Icges(4,5) * t576 + Icges(4,6) * t575 - t709) * qJD(1);
t697 = t619 ^ 2;
t696 = t625 ^ 2;
t711 = Icges(3,3) + Icges(7,3);
t620 = sin(pkin(15));
t621 = sin(pkin(14));
t626 = cos(pkin(15));
t627 = cos(pkin(14));
t577 = t620 * t627 - t621 * t626;
t578 = t620 * t621 + t626 * t627;
t535 = t577 * t624 + t578 * t618;
t536 = -t577 * t618 + t578 * t624;
t710 = Icges(3,5) * t624 + Icges(7,5) * t536 - Icges(3,6) * t618 - Icges(7,6) * t535;
t706 = Icges(7,4) * t535;
t705 = t619 * t625;
t704 = Icges(3,5) * t618 + Icges(7,5) * t535 + Icges(3,6) * t624 + Icges(7,6) * t536;
t703 = t619 * t710 - t625 * t711;
t702 = t619 * t711 + t625 * t710;
t502 = Icges(7,2) * t536 + t706;
t687 = Icges(7,4) * t536;
t503 = Icges(7,1) * t535 + t687;
t691 = Icges(3,4) * t618;
t585 = Icges(3,2) * t624 + t691;
t690 = Icges(3,4) * t624;
t586 = Icges(3,1) * t618 + t690;
t701 = -t502 * t535 + t503 * t536 - t585 * t618 + t586 * t624;
t656 = -Icges(7,2) * t535 + t687;
t498 = Icges(7,6) * t619 + t625 * t656;
t660 = Icges(7,1) * t536 - t706;
t500 = Icges(7,5) * t619 + t625 * t660;
t658 = -Icges(3,2) * t618 + t690;
t558 = Icges(3,6) * t619 + t625 * t658;
t662 = Icges(3,1) * t624 - t691;
t560 = Icges(3,5) * t619 + t625 * t662;
t700 = -t498 * t535 + t500 * t536 - t558 * t618 + t560 * t624;
t497 = -Icges(7,6) * t625 + t619 * t656;
t499 = -Icges(7,5) * t625 + t619 * t660;
t557 = -Icges(3,6) * t625 + t619 * t658;
t559 = -Icges(3,5) * t625 + t619 * t662;
t699 = t497 * t535 - t499 * t536 + t557 * t618 - t559 * t624;
t693 = pkin(1) * qJD(2);
t692 = pkin(4) * qJD(3);
t689 = Icges(4,4) * t575;
t688 = Icges(4,4) * t576;
t611 = sin(pkin(16));
t612 = cos(pkin(16));
t571 = t611 * t626 + t612 * t620;
t572 = -t611 * t620 + t612 * t626;
t609 = pkin(17) + pkin(18);
t603 = sin(t609);
t604 = cos(t609);
t667 = t571 * t604 + t572 * t603;
t684 = t667 * t619;
t622 = cos(qJ(4));
t683 = t667 * t622;
t682 = t667 * t625;
t616 = sin(qJ(4));
t681 = t667 * t616;
t601 = -pkin(4) * t623 + pkin(1);
t678 = pkin(4) * t680 + t601 * t624;
t677 = qJD(1) * t625;
t676 = qJD(2) * t619;
t675 = qJD(2) * t625;
t674 = qJD(4) * t667;
t673 = t618 * t693;
t672 = t576 * t692;
t664 = rSges(6,1) * t622 - rSges(6,2) * t616;
t516 = rSges(6,3) * t572 + t571 * t664;
t517 = -rSges(6,3) * t571 + t572 * t664;
t671 = (t516 * t604 + t517 * t603) * t674;
t599 = pkin(1) * t624 + pkin(12);
t666 = t625 * t673;
t532 = t678 * qJD(2) + t575 * t692;
t566 = -pkin(4) * t679 + t601 * t618;
t665 = -t566 * t675 + t625 * t672;
t593 = rSges(3,1) * t624 - rSges(3,2) * t618;
t663 = -rSges(9,1) * t606 + rSges(9,2) * t605;
t661 = Icges(4,1) * t575 + t688;
t657 = Icges(4,2) * t576 + t689;
t651 = -Icges(9,5) * t606 + Icges(9,6) * t605;
t647 = t516 * t603 - t517 * t604;
t582 = rSges(5,1) * t611 - rSges(5,2) * t612;
t583 = rSges(5,1) * t612 + rSges(5,2) * t611;
t646 = -(t582 * t626 + t583 * t620) * t603 + (-t582 * t620 + t583 * t626) * t604;
t528 = t571 * t603 - t572 * t604;
t591 = -t626 * rSges(7,1) + rSges(7,2) * t620;
t595 = rSges(7,1) * t620 + rSges(7,2) * t626;
t639 = t591 * t621 + t595 * t627;
t613 = sin(pkin(18));
t615 = cos(pkin(18));
t573 = t613 * t626 + t615 * t620;
t638 = t599 * t677 - t619 * t673;
t637 = pkin(12) + t593;
t565 = pkin(12) + t678;
t636 = t565 * t677 - t566 * t676 + t619 * t672;
t635 = -(-rSges(9,1) * t605 - rSges(9,2) * t606) * t608 - t673;
t634 = t616 * t528;
t512 = -t619 * t634 - t622 * t625;
t632 = t528 * t622;
t513 = -t616 * t625 + t619 * t632;
t473 = Icges(6,5) * t513 + Icges(6,6) * t512 + Icges(6,3) * t684;
t514 = t619 * t622 - t625 * t634;
t515 = t616 * t619 + t625 * t632;
t474 = Icges(6,5) * t515 + Icges(6,6) * t514 + Icges(6,3) * t682;
t633 = (t473 * t619 + t474 * t625) * t667;
t524 = -Icges(4,6) * t625 + t619 * t657;
t525 = Icges(4,6) * t619 + t625 * t657;
t526 = -Icges(4,5) * t625 + t619 * t661;
t527 = Icges(4,5) * t619 + t625 * t661;
t539 = Icges(4,2) * t575 - t688;
t540 = -Icges(4,1) * t576 + t689;
t630 = (t525 * t576 + t527 * t575) * t580 - (t524 * t576 + t526 * t575) * t581 + (t539 * t576 + t540 * t575) * qJD(1);
t628 = qJD(2) ^ 2;
t602 = t624 * t693;
t594 = rSges(2,1) * t625 - rSges(2,2) * t619;
t592 = rSges(6,1) * t616 + rSges(6,2) * t622;
t590 = rSges(2,1) * t619 + rSges(2,2) * t625;
t589 = rSges(3,1) * t618 + rSges(3,2) * t624;
t588 = pkin(8) * t612 - pkin(10) * t611;
t587 = pkin(8) * t611 + pkin(10) * t612;
t574 = -t613 * t620 + t615 * t626;
t553 = rSges(9,3) * t619 + t625 * t663;
t552 = -rSges(9,3) * t625 + t619 * t663;
t547 = Icges(9,3) * t619 + t625 * t651;
t546 = -Icges(9,3) * t625 + t619 * t651;
t545 = -t591 * t627 + t595 * t621;
t544 = (-rSges(4,1) * t624 + rSges(4,2) * t618) * t617 - t623 * (rSges(4,1) * t618 + rSges(4,2) * t624);
t543 = (-rSges(4,1) * t623 + rSges(4,2) * t617) * t624 + t618 * (rSges(4,1) * t617 + rSges(4,2) * t623);
t534 = rSges(4,3) * t619 + t543 * t625;
t533 = -rSges(4,3) * t625 + t543 * t619;
t521 = -t589 * t676 + (t619 * rSges(3,3) + t625 * t637) * qJD(1);
t520 = -t589 * t675 + (t625 * rSges(3,3) - t619 * t637) * qJD(1);
t519 = -qJD(4) * t528 + qJD(1);
t518 = (cos(pkin(17)) * t574 - t573 * sin(pkin(17))) * pkin(3) - t599;
t511 = qJD(1) * (((-rSges(8,1) * t626 - rSges(8,2) * t620) * t615 + t613 * (rSges(8,1) * t620 - rSges(8,2) * t626)) * t625 + t619 * rSges(8,3)) + t638;
t510 = -t666 + (t625 * rSges(8,3) + (-t599 + (rSges(8,1) * t615 + rSges(8,2) * t613) * t626 - t620 * (rSges(8,1) * t613 - rSges(8,2) * t615)) * t619) * qJD(1);
t509 = t545 * t624 - t618 * t639;
t508 = t618 * t545 + t624 * t639;
t507 = (-t587 * t620 + t588 * t626) * t604 - (t587 * t626 + t588 * t620) * t603;
t506 = t602 + (t552 * t619 + t553 * t625) * t608;
t505 = rSges(7,3) * t619 + t509 * t625;
t504 = -rSges(7,3) * t625 + t509 * t619;
t494 = qJD(1) * t534 - t544 * t580 + t638;
t493 = -t666 - t544 * t581 + (-t619 * t599 - t533) * qJD(1);
t492 = t533 * t580 + t534 * t581 + t602;
t491 = t635 * t619 + (-t518 * t625 + t553) * qJD(1);
t490 = t635 * t625 + (t518 * t619 - t552) * qJD(1);
t488 = Icges(6,1) * t683 - Icges(6,4) * t681 - Icges(6,5) * t528;
t487 = Icges(6,4) * t683 - Icges(6,2) * t681 - Icges(6,6) * t528;
t486 = Icges(6,5) * t683 - Icges(6,6) * t681 - Icges(6,3) * t528;
t485 = t592 * t619 + t625 * t647;
t484 = -t592 * t625 + t619 * t647;
t483 = -t508 * t676 + (-pkin(6) * t625 + t505) * qJD(1);
t482 = -t508 * t675 + (pkin(6) * t619 - t504) * qJD(1);
t481 = -(-rSges(5,3) * t619 + t625 * t646) * qJD(1) + t636;
t480 = (t625 * rSges(5,3) + (-t565 + t646) * t619) * qJD(1) + t665;
t479 = (t504 * t619 + t505 * t625) * qJD(2);
t478 = Icges(6,1) * t515 + Icges(6,4) * t514 + Icges(6,5) * t682;
t477 = Icges(6,1) * t513 + Icges(6,4) * t512 + Icges(6,5) * t684;
t476 = Icges(6,4) * t515 + Icges(6,2) * t514 + Icges(6,6) * t682;
t475 = Icges(6,4) * t513 + Icges(6,2) * t512 + Icges(6,6) * t684;
t472 = (t484 * t625 - t485 * t619) * t674 + t532;
t471 = t519 * t485 + (-qJD(1) * t507 - t671) * t625 + t636;
t470 = -t519 * t484 + (t671 + (t507 - t565) * qJD(1)) * t619 + t665;
t1 = t519 * ((-t486 * t528 - t487 * t681 + t488 * t683) * t519 + ((-t474 * t528 - t476 * t681) * t625 + (-t473 * t528 - t475 * t681) * t619 + (t477 * t619 + t478 * t625) * t683) * t674) / 0.2e1 + m(8) * (pkin(1) ^ 2 * t624 ^ 2 * t628 + t510 ^ 2 + t511 ^ 2) / 0.2e1 + m(3) * (t593 ^ 2 * t628 + t520 ^ 2 + t521 ^ 2) / 0.2e1 + m(7) * (t479 ^ 2 + t482 ^ 2 + t483 ^ 2) / 0.2e1 + m(6) * (t470 ^ 2 + t471 ^ 2 + t472 ^ 2) / 0.2e1 + m(5) * (t480 ^ 2 + t481 ^ 2 + t532 ^ 2) / 0.2e1 + m(9) * (t490 ^ 2 + t491 ^ 2 + t506 ^ 2) / 0.2e1 + m(4) * (t492 ^ 2 + t493 ^ 2 + t494 ^ 2) / 0.2e1 + ((t702 * t697 + (t699 * t625 + (t700 - t703) * t619) * t625) * qJD(2) + (t704 * t619 + t701 * t625) * qJD(1)) * t676 / 0.2e1 - ((t703 * t696 + (t700 * t619 + (t699 - t702) * t625) * t619) * qJD(2) + (t701 * t619 - t704 * t625) * qJD(1)) * t675 / 0.2e1 + (t625 * ((t486 * t682 + t514 * t487 + t515 * t488) * t519 + ((t514 * t475 + t515 * t477) * t619 + (t514 * t476 + t515 * t478 + t633) * t625) * t674) + t619 * ((t486 * t684 + t512 * t487 + t513 * t488) * t519 + ((t512 * t476 + t513 * t478) * t625 + (t512 * t475 + t513 * t477 + t633) * t619) * t674)) * t674 / 0.2e1 + (m(2) * (t590 ^ 2 + t594 ^ 2) + t574 ^ 2 * Icges(8,2) + (Icges(8,1) * t573 - 0.2e1 * Icges(8,4) * t574) * t573 + Icges(5,1) * t667 ^ 2 + (0.2e1 * Icges(5,4) * t667 + Icges(5,2) * t528) * t528 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t525 * t575 - t527 * t576) * t580 - (t524 * t575 - t526 * t576) * t581 + ((-t497 * t536 - t499 * t535 - t557 * t624 - t559 * t618) * t625 + (t498 * t536 + t500 * t535 + t558 * t624 + t560 * t618) * t619) * qJD(2) + (t606 ^ 2 * Icges(9,2) + t536 * t502 + t535 * t503 + t575 * t539 - t576 * t540 + t624 * t585 + t618 * t586 + (Icges(9,1) * t605 + 0.2e1 * Icges(9,4) * t606) * t605) * qJD(1) + (-t696 - t697) * t608 * t709) * qJD(1) / 0.2e1 + (t625 * t630 + (-t546 * t705 + t697 * t547) * t608 + t712 * t619) * t580 / 0.2e1 - (t619 * t630 + (t696 * t546 - t547 * t705) * t608 - t712 * t625) * t581 / 0.2e1;
T = t1;
