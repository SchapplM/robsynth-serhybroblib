% Calculate potential energy for
% picker2Dm1TE
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
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1TE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_energypot_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1TE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_energypot_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1TE_energypot_fixb_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:11:41
% EndTime: 2020-05-10 00:11:55
% DurationCPUTime: 10.28s
% Computational Cost: add. (86733->498), mult. (233078->680), div. (3590->12), fcn. (61716->14), ass. (0->272)
t517 = pkin(5) ^ 2;
t483 = sin(pkin(8));
t484 = cos(pkin(8));
t486 = sin(qJ(1));
t488 = cos(qJ(1));
t418 = -t483 * t486 - t484 * t488;
t689 = pkin(5) * t418;
t634 = pkin(1) * t689;
t409 = t517 - t634;
t413 = -pkin(1) + t689;
t514 = 0.2e1 * pkin(1);
t659 = -0.2e1 * t634 + t517;
t403 = sqrt(-(-(t514 + pkin(5)) * pkin(5) + t659) * (pkin(5) * (t514 - pkin(5)) + t659));
t417 = t483 * t488 - t484 * t486;
t683 = t403 * t417;
t399 = -pkin(5) * t683 - 0.2e1 * t409 * t413;
t702 = 0.2e1 * t417;
t401 = pkin(5) * t409 * t702 - t403 * t413;
t703 = 0.1e1 / pkin(3);
t686 = t703 / 0.2e1;
t512 = 0.1e1 / t686;
t530 = pkin(1) ^ 2;
t528 = t530 ^ 2;
t533 = pkin(7) ^ 2;
t456 = t488 ^ 2;
t525 = pkin(3) ^ 2;
t691 = sin(qJ(2));
t628 = t691 * pkin(7);
t601 = 0.2e1 * t628;
t591 = pkin(3) * t601;
t577 = t533 + t591;
t547 = t691 ^ 2;
t669 = t525 * t547;
t631 = 0.2e1 * t669;
t563 = -t525 + t577 + t631;
t560 = t456 * t563;
t520 = pkin(4) ^ 2;
t570 = t530 + t577;
t566 = -t520 + t570;
t621 = pkin(3) * t691;
t587 = t621 + pkin(7);
t572 = t587 * t488;
t596 = -0.4e1 * t621;
t487 = cos(qJ(2));
t674 = t487 * t486;
t625 = pkin(3) * t674;
t597 = pkin(1) * t625;
t647 = t520 - t533;
t606 = -0.2e1 * t525 + t647;
t449 = t530 + t533;
t607 = -t520 + t449;
t627 = pkin(1) * t674;
t644 = t530 - t533;
t698 = 0.2e1 * t530;
t398 = sqrt(-0.4e1 * t530 * t560 - 0.4e1 * pkin(1) * (0.2e1 * (-t627 + t628) * pkin(3) + t607) * t572 + 0.4e1 * t566 * t597 + 0.4e1 * t644 * t669 + pkin(7) * t607 * t596 - t528 + t606 * t698 - (t533 - (t512 + pkin(4)) * pkin(4)) * (t533 + (t512 - pkin(4)) * pkin(4)));
t561 = t572 - t625;
t559 = 0.1e1 / (t561 * t514 + t525 + t570);
t697 = 0.3e1 * t530;
t583 = t697 - t606;
t673 = t487 * t488;
t624 = pkin(3) * t673;
t700 = 0.4e1 * t525;
t554 = t559 * ((t486 * t587 + t624) * t398 - ((t601 - 0.4e1 * t627) * pkin(3) + t583) * t572 + (t591 + t583) * t625 + (-0.2e1 * t560 + t631 - t591 - t700 - t607) * pkin(1));
t553 = t703 * t554;
t444 = pkin(3) * t487;
t701 = 0.2e1 * t525;
t564 = t701 + t566;
t445 = pkin(1) * t488;
t637 = 0.2e1 * t445;
t558 = t559 * ((pkin(1) + t561) * t398 + (t563 * t637 + t564 * t587) * t486 + (t488 * t564 + (0.4e1 * t456 - 0.2e1) * pkin(1) * t587) * t444);
t556 = t558 / 0.4e1;
t555 = t703 * t556;
t408 = 0.1e1 / (t530 + t659);
t682 = t408 / pkin(5);
t384 = (-t401 * t553 / 0.4e1 + t399 * t555) * t682;
t489 = cos(pkin(9));
t552 = t554 / 0.4e1;
t550 = (t703 * t399 * t552 + t401 * t555) * t682;
t692 = sin(pkin(9));
t381 = t384 * t692 + t489 * t550;
t513 = 2 * pkin(2);
t688 = pkin(6) * t381;
t380 = t688 * t513;
t515 = pkin(6) ^ 2;
t663 = t380 + t515;
t371 = sqrt(-(-(t513 + pkin(6)) * pkin(6) + t663) * (pkin(6) * (t513 - pkin(6)) + t663));
t378 = t380 + 0.2e1 * t515;
t379 = -pkin(2) - t688;
t382 = t384 * t489 - t550 * t692;
t687 = pkin(6) * t382;
t363 = t371 * t687 - t378 * t379;
t364 = -t371 * t379 - t378 * t687;
t551 = -t553 / 0.2e1;
t557 = t703 * t558;
t694 = t486 / 0.2e1;
t391 = t488 * t551 + t557 * t694;
t516 = 0.1e1 / pkin(6);
t684 = 0.1e1 / ((pkin(2) ^ 2) + t663) * t516;
t390 = t486 * t551 - t488 * t557 / 0.2e1;
t696 = t390 / 0.2e1;
t705 = (t391 * t364 / 0.2e1 + t363 * t696) * t684;
t712 = t381 * t705;
t521 = 0.1e1 / pkin(4);
t670 = t521 / pkin(3) ^ 2;
t538 = t525 ^ 2;
t519 = t520 ^ 2;
t532 = t533 ^ 2;
t645 = t528 + t532;
t510 = 0.2e1 * t533;
t650 = t510 - t520;
t666 = t533 * t520;
t567 = t650 * t530 + t519 / 0.6e1 + t645 - t666;
t416 = -t538 / 0.6e1 + t567;
t477 = -t525 / 0.3e1;
t440 = t477 + t533;
t588 = -0.2e1 * t597;
t419 = t440 * t588;
t424 = t445 + t587;
t448 = -0.3e1 * t525 + t533;
t455 = t488 * t456;
t534 = pkin(1) * t530;
t664 = t534 * t455;
t633 = pkin(7) * t664;
t603 = 0.8e1 * t633;
t427 = t448 * t603;
t447 = -t520 - t525;
t509 = 0.3e1 * t533;
t436 = t509 + t447;
t679 = t436 * t530;
t428 = 0.10e2 * t679;
t473 = 0.4e1 / 0.3e1 * t525;
t467 = -t520 / 0.3e1;
t611 = t467 + t449;
t429 = t473 + t611;
t468 = -t520 / 0.2e1;
t605 = t525 + t449;
t431 = t468 + t605;
t432 = -t520 + t605;
t434 = t445 + pkin(7);
t435 = pkin(7) * t637;
t508 = 0.4e1 * t533;
t438 = (t508 + t520) * t530;
t441 = -t530 / 0.3e1 + t533;
t442 = 0.10e2 / 0.3e1 * t530;
t443 = t449 ^ 2;
t446 = -0.30e2 * t520 + 0.60e2 * t533;
t451 = -0.3e1 * t530 + t533;
t465 = -t520 / 0.6e1;
t466 = -t520 / 0.4e1;
t474 = 0.2e1 / 0.3e1 * t525;
t479 = 0.4e1 / 0.3e1 * t530;
t481 = t530 / 0.2e1;
t490 = 0.15e2 * t528;
t491 = 0.15e2 * t530;
t492 = 0.10e2 * t530;
t497 = -0.2e1 * t520;
t498 = -0.5e1 * t520;
t499 = 0.5e1 * t538;
t500 = 0.7e1 * t528;
t501 = 0.5e1 * t528;
t502 = 0.7e1 * t530;
t503 = 0.6e1 * t530;
t506 = 0.3e1 * t532;
t507 = 0.8e1 * t533;
t537 = pkin(3) * t525;
t522 = t537 ^ 2;
t542 = pkin(7) * t533;
t565 = 0.5e1 / 0.6e1 * t538 + t567;
t575 = t533 - t597;
t657 = t519 / 0.2e1 - t538 / 0.2e1;
t585 = -0.3e1 * t666 + t506 + t657;
t590 = -0.6e1 * t597;
t470 = -0.3e1 / 0.2e1 * t520;
t656 = t470 + t509;
t660 = t449 * ((t470 + t510) * t530 - 0.3e1 / 0.2e1 * t666 + t645 + t657) + t522;
t576 = ((t442 + t650) * t525 + t565) * t590 + (t490 + (-0.9e1 * t520 + 0.18e2 * t533) * t530 + t585) * t525 + (t491 + t656) * t538 + t660;
t589 = -0.4e1 * t597;
t578 = t431 * t589;
t641 = t533 + t697;
t604 = t525 + t641;
t690 = pkin(1) * t486;
t579 = -(0.3e1 * t525 + t449) * t690 + t604 * t444;
t469 = -0.2e1 / 0.3e1 * t520;
t478 = -0.2e1 / 0.3e1 * t525;
t609 = t469 + t449;
t651 = t492 + t510;
t655 = t478 + t533;
t580 = -(t499 + (0.5e1 * t530 + t436) * t701 + (t478 + t609) * t449) * t690 + (t538 + (t469 + t478 + t651) * t525 + t501 + 0.2e1 * t679 + t533 * (t469 + t655)) * t444;
t649 = t519 - t538;
t584 = -0.6e1 * t666 + 0.6e1 * t532 + t649;
t610 = t469 + t474 + t510;
t658 = (t474 + t609) * t449 + t538;
t586 = t429 * t589 + (t503 + t610) * t525 + t658;
t593 = t664 * t444;
t668 = t528 * t456 ^ 2;
t594 = t668 * t444;
t622 = 0.16e2 * t664;
t600 = pkin(7) * t622;
t602 = 0.20e2 / 0.3e1 * t530;
t648 = -t520 + t525;
t608 = t509 + t648;
t475 = t525 / 0.3e1;
t612 = t465 + t475 + t533;
t613 = t520 / 0.3e1 + t475 + t510;
t614 = 0.2e1 / 0.3e1 * t520 + t474 + t508;
t615 = 0.4e1 / 0.3e1 * t520 + t473 - 0.2e1 * t533;
t672 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t616 = t486 * t672;
t635 = 0.6e1 * t445;
t618 = pkin(7) * t635;
t636 = 0.4e1 * t445;
t619 = pkin(7) * t636;
t620 = -t690 / 0.2e1;
t667 = t530 * t456;
t623 = 0.12e2 * t667;
t626 = t530 * t444;
t629 = 0.4e1 * t667;
t630 = 0.8e1 * t668;
t675 = t691 * t547 * t537;
t632 = -0.8e1 * t675;
t638 = 0.2e1 * t690;
t639 = pkin(7) * t445;
t640 = 0.4e1 * pkin(7);
t642 = t532 + t538;
t643 = t532 - t528;
t646 = -t525 + t533;
t652 = 0.4e1 / 0.7e1 * t533 - t520 / 0.7e1;
t653 = t481 + t533;
t654 = t530 / 0.3e1 + t533;
t665 = t533 * t530;
t671 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t676 = t547 ^ 2 * t538;
t677 = t447 * t533;
t678 = t443 * (-t525 + t607);
t680 = (-t486 * t534 + t626) * t456;
t699 = 0.4e1 * t528;
t704 = t466 + t525 / 0.2e1;
t685 = ((-0.24e2 * (0.4e1 / 0.3e1 * t667 + t435 + t441) * t676 * t690 - 0.12e2 * (-0.8e1 / 0.3e1 * t594 + ((t479 + t612) * t444 - (0.7e1 / 0.6e1 * t525 + t465 + t653) * t690) * t629 + (-t525 * t644 - 0.5e1 / 0.3e1 * t528 + t613 * t530 + t533 * (t467 + t440)) * t444 + (-t538 + (-t602 + t614) * t525 - 0.3e1 * t528 + t615 * t530 + t532) * t620 + (-t486 * t528 * t455 + ((t530 + t612) * t444 + (t701 - t644) * t620) * t445) * t640) * t669 + 0.24e2 * t440 * t594 + ((t533 + 0.5e1 / 0.2e1 * t525 + 0.3e1 / 0.2e1 * t530 + t468) * t444 + t448 * t690 / 0.2e1) * t600 - 0.6e1 * ((-0.3e1 * t538 + (-t602 + t615) * t525 + t614 * t530 + t643) * t444 - 0.2e1 * (-0.5e1 / 0.3e1 * t538 + (-t530 + t613) * t525 + t533 * (t477 + t611)) * t690) * t667 - 0.6e1 * t580 * t639 - (t522 + (0.21e2 * t530 + t436) * t538 + (t428 + t506 + 0.35e2 * t528 + 0.2e1 * t677) * t525 + (t500 + (t498 + t507 - 0.5e1 * t525) * t530 + t533 * (-t520 + t646)) * t449) * t444 + (0.7e1 * t522 + (t502 + t436) * t499 + (t428 + 0.21e2 * t528 + 0.9e1 * t532 + 0.6e1 * t677) * t525 + t678) * t690) * t398 + (0.16e2 * (t630 + t600 + (-0.8e1 * t528 + 0.12e2 * t665) * t456 + (0.4e1 * pkin(1) * t542 - 0.12e2 * pkin(7) * t534) * t488 - 0.6e1 * t665 + t645) * t676 + 0.24e2 * (t655 * t630 + 0.14e2 * (-0.32e2 / 0.21e2 * (t533 + t525 / 0.4e1 + t530 / 0.4e1 - t520 / 0.8e1) * t597 + 0.5e1 / 0.42e2 * t538 + (0.16e2 / 0.21e2 * t530 + t652) * t525 + t528 / 0.7e1 + t652 * t530 + t532 - 0.3e1 / 0.7e1 * t666 + t519 / 0.42e2) * t667 + t441 * t578 - t644 * t538 + (t438 - 0.10e2 / 0.3e1 * t528 + 0.2e1 * t532 - t666) * t525 + t416 * t671 + ((-0.2e1 / 0.3e1 * t597 + t466 + t653) * t622 + (-0.8e1 / 0.3e1 * (t654 + t704) * t597 + 0.5e1 / 0.18e2 * t538 + (0.4e1 / 0.3e1 * t533 + t479 + t467) * t525 + t532 + 0.2e1 / 0.3e1 * t665 - 0.2e1 / 0.3e1 * t666 - t528 / 0.3e1 + t519 / 0.18e2) * t635) * pkin(7)) * t669 + 0.16e2 * (-0.6e1 * t533 * t525 + t642) * t668 + 0.32e2 * (t431 * t448 + t588 * t672) * t633 + 0.24e2 * (t440 * t578 - t522 + (-t442 + t647) * t538 + (t438 + t538 / 0.6e1 - t519 / 0.6e1 + t643) * t525 + t416 * t533) * t667 + 0.8e1 * t576 * t639 - 0.8e1 * ((t502 + t656) * t538 + (t500 + (t498 + 0.10e2 * t533) * t530 + t585) * t525 + t660) * t597 + t538 ^ 2 + (t497 + t508 + 0.28e2 * t530) * t522 + (t446 * t530 + 0.70e2 * t528 + t584) * t538 + (t446 * t528 + t503 * t584 + t510 * t649 - 0.6e1 * t532 * t520 + 0.28e2 * t534 ^ 2 + 0.4e1 * t542 ^ 2) * t525 + t432 * t678) * t424 + (((0.4e1 * t680 + (t444 + t638) * t435 + t451 * t444 + (t468 + t604) * t638) * t632 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t525 + t481 + t465) * t487 * t512 + pkin(1) * t616) * t667 + (-0.8e1 * t593 + ((t467 + t474 + t641) * t444 - (0.8e1 / 0.3e1 * t525 + t611) * t690) * t636) * pkin(7) + t580) * t621) * t398 + (0.32e2 * (t603 + (-0.4e1 * t534 * t625 + t699 + (t700 + t497 + t507) * t530) * t456 + (-t530 + t575 + t704) * t619 + t588 * t671 + t451 * t431) * t675 + 0.8e1 * (t427 + (t431 * t672 + t419) * t623 + (t578 + (t503 + t650) * t525 + t565) * t618 + t576) * t621) * t424) * t434) / ((-0.4e1 * (-t644 * t444 + 0.2e1 * t680 + (0.2e1 * pkin(7) * t624 + t486 * (t525 + t698)) * pkin(1)) * t669 + 0.8e1 * pkin(7) * t593 + ((pkin(3) * t699 + 0.8e1 * t530 * t537) * t487 + 0.4e1 * t534 * t616) * t456 - 0.4e1 * t579 * t639 - (t525 * t651 + t501 + t642 + 0.6e1 * t665) * t444 + (t499 + (t492 + 0.6e1 * t533) * t525 + t443) * t690) * t398 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t597 + 0.4e1 / 0.9e1 * t525 - t520 / 0.9e1 + t654) * t667 + t441 * t588 + t429 * t671 + (t664 + (t465 + t474 + t575) * t445) * t640) * t669 + t427 + (t429 * t672 + t419) * t623 + t586 * t618 + ((t442 + t610) * t525 + t658) * t590 + t522 + (t491 + t608) * t538 + (t503 * t608 + t510 * t648 + t490 + t506) * t525 + t443 * t432) * t424 + ((t632 * t690 + (-0.2e1 * t456 * t626 + (t444 - t690) * t435 + t579) * t596) * t398 + (0.8e1 * (t435 + t629 + t451) * t675 + 0.6e1 * (t646 * t629 + (t429 + t588) * t619 + t586) * t621) * t424) * t434);
t372 = (t552 * t685 - t398 * t558 / 0.4e1) * t670;
t373 = (t398 * t552 + t556 * t685) * t670;
t420 = -t486 * t691 - t673;
t421 = -t488 * t691 + t674;
t365 = t372 * t420 + t373 * t421;
t366 = t372 * t421 - t373 * t420;
t592 = t685 * t686;
t571 = (-t390 * t703 * t398 / 0.2e1 + t391 * t592) * t521;
t706 = (t391 * t398 * t686 + t390 * t592) * t521;
t711 = t365 * t706 + t366 * t571;
t710 = -t365 * t571 + t366 * t706;
t569 = (-t391 * t363 / 0.2e1 + t364 * t696) * t684;
t707 = t381 * t569;
t695 = -t483 / 0.2e1;
t693 = t488 / 0.2e1;
t681 = t408 / pkin(1);
t662 = t391 * pkin(3) - t445;
t661 = t391 * pkin(2) - t445;
t617 = t371 * t516 / pkin(2);
t599 = t390 * pkin(3) - t690;
t598 = t390 * pkin(2) - t690;
t595 = -t391 * t381 - t382 * t390;
t582 = -t617 / 0.2e1;
t581 = t617 / 0.2e1;
t574 = t381 * t390 - t382 * t391;
t414 = -pkin(1) * t418 + pkin(5);
t410 = t530 - t634;
t402 = pkin(1) * t410 * t702 + t403 * t414;
t400 = -pkin(1) * t683 + 0.2e1 * t410 * t414;
t396 = (t399 * t693 - t486 * t401 / 0.2e1) * t682;
t395 = (t399 * t694 + t401 * t693) * t682;
t394 = (-t484 * t400 / 0.2e1 + t402 * t695) * t681;
t393 = (t400 * t695 + t484 * t402 / 0.2e1) * t681;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (-rSges(2,1) * t488 + rSges(2,2) * t486) + g(2) * (-rSges(2,1) * t486 - rSges(2,2) * t488) + g(3) * rSges(2,3)) - m(3) * (g(1) * (rSges(3,1) * t391 - rSges(3,2) * t390 - t445) + g(2) * (rSges(3,1) * t390 + rSges(3,2) * t391 - t690) + g(3) * rSges(3,3)) - m(4) * (g(1) * (rSges(4,1) * t569 + rSges(4,2) * t705 + t661) + g(2) * (-rSges(4,1) * t705 + rSges(4,2) * t569 + t598) + g(3) * rSges(4,3)) - m(5) * (g(1) * (rSges(5,1) * t571 - rSges(5,2) * t706 + t662) + g(2) * (rSges(5,1) * t706 + rSges(5,2) * t571 + t599) + g(3) * rSges(5,3)) - m(6) * (g(1) * (pkin(5) * t484 + rSges(6,1) * t394 - rSges(6,2) * t393) + g(2) * (pkin(5) * t483 + rSges(6,1) * t393 + rSges(6,2) * t394) + g(3) * rSges(6,3)) - m(7) * (g(1) * (rSges(7,1) * t595 + rSges(7,2) * t574 - t445) + g(2) * (-rSges(7,1) * t574 + rSges(7,2) * t595 - t690) + g(3) * rSges(7,3)) - m(8) * (g(1) * (rSges(8,1) * t691 + t487 * rSges(8,2) + pkin(7)) + g(2) * (-t487 * rSges(8,1) + rSges(8,2) * t691) + g(3) * rSges(8,3)) - m(9) * (g(1) * (rSges(9,1) * t396 - rSges(9,2) * t395 - t445) + g(2) * (rSges(9,1) * t395 + rSges(9,2) * t396 - t690) + g(3) * rSges(9,3)) - m(10) * (g(1) * (t569 * pkin(6) + (t582 * t705 + t707) * rSges(10,1) + (t569 * t581 + t712) * rSges(10,2) + t661) + g(2) * (-t705 * pkin(6) + (t569 * t582 - t712) * rSges(10,1) + (-t581 * t705 + t707) * rSges(10,2) + t598) + g(3) * rSges(10,3)) - m(11) * (g(1) * (t571 * pkin(4) + t711 * rSges(11,1) - t710 * rSges(11,2) + t662) + g(2) * (t706 * pkin(4) + t710 * rSges(11,1) + t711 * rSges(11,2) + t599) + g(3) * rSges(11,3));
U = t1;