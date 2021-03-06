% Calculate potential energy for
% picker2Dm1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 19:54
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE1_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE1_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1DE1_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 09:11:15
% EndTime: 2020-05-10 09:11:35
% DurationCPUTime: 12.25s
% Computational Cost: add. (173261->504), mult. (465802->630), div. (7180->12), fcn. (123736->38), ass. (0->286)
t534 = 2 * pkin(2);
t501 = sin(pkin(8));
t502 = cos(pkin(8));
t504 = sin(qJ(1));
t507 = cos(qJ(1));
t439 = -t501 * t504 - t502 * t507;
t705 = pkin(5) * t439;
t715 = -2 * pkin(1);
t436 = t705 * t715;
t538 = pkin(5) ^ 2;
t430 = t436 + 0.2e1 * t538;
t434 = -pkin(1) + t705;
t535 = 2 * pkin(1);
t680 = t436 + t538;
t424 = sqrt(-(-(t535 + pkin(5)) * pkin(5) + t680) * (pkin(5) * (t535 - pkin(5)) + t680));
t438 = t501 * t507 - t502 * t504;
t700 = t424 * t438;
t422 = -pkin(5) * t700 - t430 * t434;
t423 = pkin(5) * t430 * t438 - t424 * t434;
t549 = pkin(1) ^ 2;
t547 = t549 ^ 2;
t552 = pkin(7) ^ 2;
t473 = t507 ^ 2;
t556 = pkin(3) ^ 2;
t707 = sin(qJ(2));
t640 = pkin(3) * t707;
t616 = pkin(7) * t640;
t609 = 0.2e1 * t616;
t598 = t552 + t609;
t567 = t707 ^ 2;
t685 = t556 * t567;
t648 = 0.2e1 * t685;
t583 = -t556 + t598 + t648;
t581 = t473 * t583;
t541 = pkin(4) ^ 2;
t589 = t549 + t598;
t587 = -t541 + t589;
t605 = t640 + pkin(7);
t591 = t605 * t507;
t506 = cos(qJ(2));
t695 = t506 * t504;
t644 = pkin(3) * t695;
t615 = pkin(1) * t644;
t711 = 0.2e1 * t556;
t661 = t552 + t711;
t623 = t541 - t661;
t662 = -t552 + t541;
t625 = t549 - t662;
t646 = pkin(1) * t695;
t647 = t707 * pkin(7);
t665 = t549 - t552;
t716 = 0.1e1 / pkin(3);
t702 = t716 / 0.2e1;
t723 = 0.1e1 / t702;
t579 = sqrt(-0.4e1 * pkin(1) * (0.2e1 * (-t646 + t647) * pkin(3) + t625) * t591 + 0.4e1 * t587 * t615 + 0.4e1 * t665 * t685 - 0.4e1 * t625 * t616 - t547 - (t552 - (t723 + pkin(4)) * pkin(4)) * (t552 + (t723 - pkin(4)) * pkin(4)) + (-0.4e1 * t581 + 0.2e1 * t623) * t549);
t582 = t591 - t644;
t580 = 0.1e1 / (t535 * t582 + t556 + t589);
t713 = 3 * t549;
t602 = t713 - t623;
t694 = t506 * t507;
t643 = pkin(3) * t694;
t710 = 0.4e1 * t556;
t575 = t580 * ((t605 * t504 + t643) * t579 + t581 * t715 - ((-0.4e1 * t646 + 0.2e1 * t647) * pkin(3) + t602) * t591 + (t609 + t602) * t644 + (t648 - t609 - t710 - t625) * pkin(1));
t573 = t575 / 0.4e1;
t462 = pkin(3) * t506;
t584 = t711 + t587;
t463 = pkin(1) * t507;
t654 = 0.2e1 * t463;
t578 = t580 * ((pkin(1) + t582) * t579 + (t583 * t654 + t584 * t605) * t504 + (t507 * t584 + (0.4e1 * t473 - 0.2e1) * pkin(1) * t605) * t462);
t577 = t578 / 0.4e1;
t576 = t716 * t577;
t429 = 0.1e1 / (t549 + t680);
t699 = t429 / pkin(5);
t407 = (t422 * t573 * t716 + t423 * t576) * t699;
t574 = t716 * t575;
t408 = (-t423 * t574 / 0.4e1 + t422 * t576) * t699;
t505 = sin(pkin(9));
t508 = cos(pkin(9));
t399 = t407 * t508 + t408 * t505;
t704 = pkin(6) * t399;
t398 = t704 * t534;
t536 = pkin(6) ^ 2;
t684 = t398 + t536;
t385 = sqrt(-(-(t534 + pkin(6)) * pkin(6) + t684) * (pkin(6) * (t534 - pkin(6)) + t684));
t396 = t398 + 0.2e1 * t536;
t397 = -pkin(2) - t704;
t708 = 0.1e1 / pkin(6) / 0.2e1;
t620 = 0.1e1 / ((pkin(2) ^ 2) + t684) * t708;
t400 = t407 * t505 - t408 * t508;
t703 = pkin(6) * t400;
t373 = atan2((-t385 * t397 + t396 * t703) * t620, (-t385 * t703 - t396 * t397) * t620);
t371 = sin(t373);
t372 = cos(t373);
t572 = atan2(t578 * t702, t574 / 0.2e1);
t570 = sin(t572);
t571 = cos(t572);
t405 = -t504 * t571 - t507 * t570;
t406 = t504 * t570 - t507 * t571;
t369 = t371 * t406 + t372 * t405;
t383 = atan2(0.1e1 / pkin(2) * t385 * t708, -t399);
t381 = sin(t383);
t382 = cos(t383);
t595 = t371 * t405 - t372 * t406;
t727 = t369 * t381 + t382 * t595;
t726 = t369 * t382 - t381 * t595;
t542 = 0.1e1 / pkin(4);
t691 = t542 / t556;
t558 = t556 ^ 2;
t540 = t541 ^ 2;
t551 = t552 ^ 2;
t666 = t547 + t551;
t533 = 0.2e1 * t552;
t669 = t533 - t541;
t688 = t552 * t541;
t588 = t669 * t549 + t540 / 0.6e1 + t666 - t688;
t437 = -t558 / 0.6e1 + t588;
t495 = -t556 / 0.3e1;
t459 = t495 + t552;
t606 = -0.2e1 * t615;
t440 = t459 * t606;
t444 = t463 + t605;
t465 = t552 - 0.3e1 * t556;
t472 = t507 * t473;
t553 = pkin(1) * t549;
t686 = t553 * t472;
t651 = pkin(7) * t686;
t618 = 0.8e1 * t651;
t447 = t465 * t618;
t491 = 0.4e1 / 0.3e1 * t556;
t466 = t549 + t552;
t485 = -t541 / 0.3e1;
t629 = t485 + t466;
t448 = t491 + t629;
t486 = -t541 / 0.2e1;
t659 = t556 + t552;
t622 = t549 + t659;
t450 = t486 + t622;
t511 = 10 * t549;
t532 = 0.3e1 * t552;
t717 = t532 - t541 - t556;
t451 = t717 * t511;
t452 = -t541 + t622;
t455 = pkin(7) * t654;
t531 = 0.4e1 * t552;
t457 = (t531 + t541) * t549;
t460 = -t549 / 0.3e1 + t552;
t461 = t466 ^ 2;
t464 = -0.30e2 * t541 + 0.60e2 * t552;
t468 = -(3 * t549) + t552;
t474 = 0.10e2 / 0.3e1 * t549;
t475 = -0.20e2 / 0.3e1 * t549;
t483 = -t541 / 0.6e1;
t484 = -t541 / 0.4e1;
t488 = -0.3e1 / 0.2e1 * t541;
t492 = 0.2e1 / 0.3e1 * t556;
t497 = 0.4e1 / 0.3e1 * t549;
t499 = t549 / 0.2e1;
t509 = 15 * t547;
t516 = 0.18e2 * t552;
t517 = -0.2e1 * t541;
t518 = -0.5e1 * t541;
t519 = -0.6e1 * t541;
t520 = 0.5e1 * t558;
t523 = 7 * t547;
t524 = 5 * t547;
t525 = 6 * t549;
t526 = 2 * t549;
t528 = 0.3e1 * t551;
t529 = 0.8e1 * t552;
t530 = 0.6e1 * t552;
t557 = pkin(3) * t556;
t543 = t557 ^ 2;
t562 = pkin(7) * t552;
t585 = 0.5e1 / 0.6e1 * t558 + t588;
t645 = t549 * t462;
t586 = t473 * (-t504 * t553 + t645);
t668 = t540 - t558;
t590 = 0.6e1 * t551 + t668 - 0.6e1 * t688;
t596 = t552 - t615;
t678 = t540 / 0.2e1 - t558 / 0.2e1;
t603 = -0.3e1 * t688 + t528 + t678;
t608 = -0.6e1 * t615;
t673 = (15 * t549) + t532;
t681 = t466 * ((t488 + t533) * t549 - 0.3e1 / 0.2e1 * t688 + t666 + t678) + t543;
t597 = ((t474 + t669) * t556 + t585) * t608 + (t509 + (t516 - 0.9e1 * t541) * t549 + t603) * t556 + (t488 + t673) * t558 + t681;
t607 = -0.4e1 * t615;
t599 = t450 * t607;
t621 = t713 + t659;
t706 = pkin(1) * t504;
t600 = -(0.3e1 * t556 + t466) * t706 + t621 * t462;
t487 = -0.2e1 / 0.3e1 * t541;
t496 = -0.2e1 / 0.3e1 * t556;
t671 = t517 - 0.2e1 * t556;
t626 = t530 + t671;
t627 = t487 + t466;
t672 = t511 + t533;
t677 = t496 + t552;
t601 = -(t520 + (t511 + t626) * t556 + (t496 + t627) * t466) * t706 + (t558 + (t487 + t496 + t672) * t556 + t524 + t626 * t549 + t552 * (t487 + t677)) * t462;
t628 = t487 + t492 + t533;
t679 = (t492 + t627) * t466 + t558;
t604 = t448 * t607 + (t525 + t628) * t556 + t679;
t454 = t463 + pkin(7);
t610 = t454 * t640;
t611 = t686 * t462;
t690 = t547 * t473 ^ 2;
t612 = t690 * t462;
t635 = t454 * t707 * t567 * t557;
t614 = -0.8e1 * t635;
t641 = 0.16e2 * t686;
t617 = pkin(7) * t641;
t660 = -t556 + t552;
t624 = -t541 + t660;
t493 = t556 / 0.3e1;
t630 = t483 + t493 + t552;
t631 = t541 / 0.3e1 + t493 + t533;
t632 = 0.2e1 / 0.3e1 * t541 + t492 + t531;
t633 = 0.4e1 / 0.3e1 * t541 + t491 - 0.2e1 * t552;
t693 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t634 = t504 * t693;
t652 = 0.6e1 * t463;
t636 = pkin(7) * t652;
t653 = 0.4e1 * t463;
t637 = pkin(7) * t653;
t639 = -t706 / 0.2e1;
t689 = t549 * t473;
t642 = 0.12e2 * t689;
t649 = 0.4e1 * t689;
t650 = 0.8e1 * t690;
t655 = 0.2e1 * t706;
t656 = pkin(7) * t463;
t658 = 0.4e1 * pkin(7);
t663 = t551 + t558;
t664 = t551 - t547;
t670 = t518 - 0.5e1 * t556;
t674 = 0.4e1 / 0.7e1 * t552 - t541 / 0.7e1;
t675 = t499 + t552;
t676 = t549 / 0.3e1 + t552;
t687 = t552 * t549;
t692 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t696 = t567 ^ 2 * t558;
t697 = t461 * (t549 + t624);
t712 = -0.6e1 * t556;
t714 = 4 * t547;
t718 = t484 + t556 / 0.2e1;
t701 = ((-0.24e2 * (0.4e1 / 0.3e1 * t689 + t455 + t460) * t696 * t706 + (0.4e1 * t586 + (t462 + t655) * t455 + t468 * t462 + (t486 + t621) * t655) * t614 - 0.12e2 * (-0.8e1 / 0.3e1 * t612 + ((t497 + t630) * t462 - (0.7e1 / 0.6e1 * t556 + t483 + t675) * t706) * t649 + (-t556 * t665 - 0.5e1 / 0.3e1 * t547 + t631 * t549 + t552 * (t485 + t459)) * t462 + (-t558 + (t475 + t632) * t556 - (3 * t547) + t633 * t549 + t551) * t639 + (-t504 * t547 * t472 + ((t549 + t630) * t462 + (-t549 + t661) * t639) * t463) * t658) * t685 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t556 + t499 + t483) * t506 * t723 + pkin(1) * t634) * t689 + (-0.8e1 * t611 + ((t485 + t492 + t713 + t552) * t462 - (0.8e1 / 0.3e1 * t556 + t629) * t706) * t653) * pkin(7) + t601) * t610 + 0.24e2 * t459 * t612 + ((t552 + 0.5e1 / 0.2e1 * t556 + 0.3e1 / 0.2e1 * t549 + t486) * t462 + t465 * t706 / 0.2e1) * t617 - 0.6e1 * ((-0.3e1 * t558 + (t475 + t633) * t556 + t632 * t549 + t664) * t462 - 0.2e1 * (-0.5e1 / 0.3e1 * t558 + (-t549 + t631) * t556 + t552 * (t495 + t629)) * t706) * t689 - 0.6e1 * t601 * t656 - (t543 + ((21 * t549) + t717) * t558 + (t552 * t671 + t451 + t528 + (35 * t547)) * t556 + (t523 + (t529 + t670) * t549 + t552 * t624) * t466) * t462 + (0.7e1 * t543 + ((35 * t549) + 0.15e2 * t552 + t670) * t558 + ((21 * t547) + t451 + 0.9e1 * t551 + (t519 + t712) * t552) * t556 + t697) * t706) * t579 + (0.16e2 * (t650 + t617 + (-(8 * t547) + 0.12e2 * t687) * t473 + (0.4e1 * pkin(1) * t562 - 0.12e2 * pkin(7) * t553) * t507 - 0.6e1 * t687 + t666) * t696 + 0.32e2 * (t618 + (-0.4e1 * t553 * t644 + t714 + (t710 + t517 + t529) * t549) * t473 + (-t549 + t596 + t718) * t637 + t606 * t692 + t468 * t450) * t635 + 0.24e2 * (t677 * t650 + 0.14e2 * (-0.32e2 / 0.21e2 * (t552 + t556 / 0.4e1 + t549 / 0.4e1 - t541 / 0.8e1) * t615 + 0.5e1 / 0.42e2 * t558 + (0.16e2 / 0.21e2 * t549 + t674) * t556 + t547 / 0.7e1 + t674 * t549 + t551 - 0.3e1 / 0.7e1 * t688 + t540 / 0.42e2) * t689 + t460 * t599 - t665 * t558 + (t457 - 0.10e2 / 0.3e1 * t547 + 0.2e1 * t551 - t688) * t556 + t437 * t692 + ((-0.2e1 / 0.3e1 * t615 + t484 + t675) * t641 + (-0.8e1 / 0.3e1 * (t676 + t718) * t615 + 0.5e1 / 0.18e2 * t558 + (0.4e1 / 0.3e1 * t552 + t497 + t485) * t556 + t551 + 0.2e1 / 0.3e1 * t687 - 0.2e1 / 0.3e1 * t688 - t547 / 0.3e1 + t540 / 0.18e2) * t652) * pkin(7)) * t685 + 0.8e1 * (t447 + (t450 * t693 + t440) * t642 + (t599 + (t525 + t669) * t556 + t585) * t636 + t597) * t610 + 0.16e2 * (t552 * t712 + t663) * t690 + 0.32e2 * (t450 * t465 + t606 * t693) * t651 + 0.24e2 * (t459 * t599 - t543 + (-t474 + t662) * t558 + (t457 + t558 / 0.6e1 - t540 / 0.6e1 + t664) * t556 + t437 * t552) * t689 + 0.8e1 * t597 * t656 - 0.8e1 * ((t488 + t532 + (7 * t549)) * t558 + (t523 + (t518 + 0.10e2 * t552) * t549 + t603) * t556 + t681) * t615 + t558 ^ 2 + (t517 + t531 + (28 * t549)) * t543 + (t464 * t549 + (70 * t547) + t590) * t558 + (t464 * t547 + t519 * t551 + t590 * t525 + t668 * t533 + (28 * t553 ^ 2) + 0.4e1 * t562 ^ 2) * t556 + t452 * t697) * t444) / ((t614 * t706 - 0.4e1 * (0.2e1 * t586 - t665 * t462 + (0.2e1 * pkin(7) * t643 + t504 * (t526 + t556)) * pkin(1)) * t685 - 0.4e1 * (-0.2e1 * t473 * t645 + (t462 - t706) * t455 + t600) * t610 + 0.8e1 * pkin(7) * t611 + ((pkin(3) * t714 + 0.8e1 * t549 * t557) * t506 + 0.4e1 * t553 * t634) * t473 - 0.4e1 * t600 * t656 - (t556 * t672 + t524 + t663 + 0.6e1 * t687) * t462 + (t520 + (t511 + t530) * t556 + t461) * t706) * t579 + (0.8e1 * (t455 + t649 + t468) * t635 + 0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t615 + 0.4e1 / 0.9e1 * t556 - t541 / 0.9e1 + t676) * t689 + t460 * t606 + t448 * t692 + (t686 + (t483 + t492 + t596) * t463) * t658) * t685 + 0.6e1 * (t660 * t649 + (t448 + t606) * t637 + t604) * t610 + t447 + (t448 * t693 + t440) * t642 + t604 * t636 + ((t474 + t628) * t556 + t679) * t608 + t543 + (-t541 + t556 + t673) * t558 + (t509 + (t516 + t519 + 0.6e1 * t556) * t549 + t528 + (t517 + t711) * t552) * t556 + t461 * t452) * t444);
t389 = (t573 * t701 - t579 * t578 / 0.4e1) * t691;
t390 = (t573 * t579 + t577 * t701) * t691;
t441 = -t504 * t707 - t694;
t442 = -t507 * t707 + t695;
t376 = atan2(t389 * t441 + t390 * t442, -t389 * t442 + t390 * t441);
t374 = sin(t376);
t375 = cos(t376);
t638 = t542 * t702;
t393 = atan2(t579 * t638, t638 * t701);
t391 = sin(t393);
t392 = cos(t393);
t378 = t391 * t406 + t392 * t405;
t613 = -t391 * t405 + t406 * t392;
t725 = t374 * t378 - t375 * t613;
t724 = t374 * t613 + t375 * t378;
t709 = pkin(5) * m(6);
t698 = t429 / pkin(1);
t683 = t406 * pkin(3) - t463;
t682 = t406 * pkin(2) - t463;
t657 = m(3) + m(7) + m(9);
t619 = t699 / 0.2e1;
t388 = atan2(t400, t399);
t386 = sin(t388);
t387 = cos(t388);
t594 = t386 * t406 + t387 * t405;
t593 = t386 * t405 - t387 * t406;
t435 = -pkin(1) * t439 + pkin(5);
t431 = t436 + t526;
t421 = atan2((pkin(1) * t431 * t438 + t424 * t435) * t698 / 0.2e1, -(-pkin(1) * t700 + t431 * t435) * t698 / 0.2e1);
t420 = atan2(t423 * t619, t422 * t619);
t419 = cos(t421);
t418 = cos(t420);
t417 = sin(t421);
t416 = sin(t420);
t413 = -t416 * t504 + t418 * t507;
t412 = t416 * t507 + t418 * t504;
t411 = -t417 * t501 + t419 * t502;
t410 = t417 * t502 + t419 * t501;
t402 = t405 * pkin(2);
t401 = t405 * pkin(3);
t1 = (-mrSges(7,3) - mrSges(9,3) - mrSges(6,3) - mrSges(11,3) - mrSges(10,3) - mrSges(8,3) - mrSges(5,3) - mrSges(4,3) - mrSges(3,3) - mrSges(2,3) - mrSges(1,3)) * g(3) + (-mrSges(3,1) * t405 - mrSges(3,2) * t406 - mrSges(6,1) * t410 - mrSges(6,2) * t411 - mrSges(9,1) * t412 - mrSges(9,2) * t413 + mrSges(2,2) * t507 - t595 * mrSges(4,2) - t613 * mrSges(5,2) - m(11) * (pkin(4) * t378 + t401) - m(5) * t401 - m(4) * t402 + t369 * mrSges(4,1) - m(10) * (-pkin(6) * t369 + t402) + t724 * mrSges(11,1) - t725 * mrSges(11,2) - t726 * mrSges(10,1) + t727 * mrSges(10,2) + (mrSges(2,1) + (m(11) + m(4) + m(5) + m(10) + t657) * pkin(1)) * t504 - t707 * mrSges(8,2) - t501 * t709 + t506 * mrSges(8,1) - t593 * mrSges(7,2) + t594 * mrSges(7,1) - t378 * mrSges(5,1) - mrSges(1,2)) * g(2) + (-m(4) * t682 - m(5) * t683 - mrSges(6,1) * t411 + mrSges(6,2) * t410 - mrSges(9,1) * t413 + mrSges(9,2) * t412 - mrSges(3,1) * t406 + mrSges(3,2) * t405 - mrSges(2,2) * t504 - m(10) * (pkin(6) * t595 + t682) - t595 * mrSges(4,1) - m(11) * (pkin(4) * t613 + t683) - t613 * mrSges(5,1) - t724 * mrSges(11,2) - t725 * mrSges(11,1) + t726 * mrSges(10,2) + t727 * mrSges(10,1) + (pkin(1) * t657 + mrSges(2,1)) * t507 - t707 * mrSges(8,1) + t378 * mrSges(5,2) - t502 * t709 - t506 * mrSges(8,2) - t593 * mrSges(7,1) - t594 * mrSges(7,2) - m(8) * pkin(7) - t369 * mrSges(4,2) - mrSges(1,1)) * g(1);
U = t1;
