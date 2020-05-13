% Calculate joint inertia matrix with Newton Euler for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_inertiaJ_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_inertiaJ_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_inertiaJ_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_inertiaJ_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2IC_inertiaJ_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 05:00:17
% EndTime: 2020-05-07 05:00:21
% DurationCPUTime: 3.44s
% Computational Cost: add. (5218->441), mult. (7399->534), div. (109->12), fcn. (8004->37), ass. (0->230)
t601 = 0.1e1 / pkin(9);
t652 = (qJ(4) + pkin(14));
t643 = (qJ(3) + t652);
t665 = (qJ(7) + qJ(8));
t657 = cos((2 * t643)) - cos((2 * pkin(15) - 2 * t665));
t673 = (t657 * pkin(9) + (cos((2 * qJ(3) + t652)) - cos((-2 * qJ(7) + 2 * pkin(15) - 2 * qJ(8) - t652))) * pkin(4)) / t657;
t646 = t601 * t673;
t588 = sin(qJ(5));
t595 = cos(qJ(5));
t656 = Ifges(6,5) * t588 + Ifges(6,6) * t595;
t625 = t588 * mrSges(6,1) + t595 * mrSges(6,2);
t700 = t625 * pkin(10);
t701 = t656 - t700;
t703 = t646 * t701;
t633 = pkin(15) + t643;
t666 = (qJ(2) - qJ(6));
t617 = t633 - t666;
t610 = -2 * qJ(7) - pkin(16) + t617;
t618 = t633 + t666;
t612 = pkin(16) + t618;
t696 = -cos((qJ(8) - t610)) + cos((qJ(8) - t612));
t506 = 0.1e1 / t696;
t604 = 0.1e1 / pkin(2);
t672 = t506 * t604;
t653 = qJ(7) + pkin(16);
t634 = t653 + t666;
t688 = pkin(3) * (pkin(1) * (cos((qJ(8) - t666)) - cos((qJ(8) + t666))) + (cos((qJ(8) - t634)) - cos((qJ(8) + t634))) * pkin(2));
t638 = t672 * t688;
t628 = t601 * t638;
t702 = t701 * t628;
t597 = cos(qJ(3));
t682 = t597 * pkin(1);
t571 = pkin(4) - t682;
t589 = sin(qJ(4));
t590 = sin(qJ(3));
t596 = cos(qJ(4));
t667 = t590 * t596;
t531 = -pkin(1) * t667 + t589 * t571;
t527 = pkin(10) + t531;
t699 = t625 * t527;
t685 = t589 * pkin(4);
t569 = pkin(10) + t685;
t698 = t625 * t569;
t591 = sin(qJ(2));
t598 = cos(qJ(2));
t550 = t590 * t591 - t597 * t598;
t552 = -t590 * t598 - t597 * t591;
t494 = -t596 * t550 + t589 * t552;
t495 = t589 * t550 + t596 * t552;
t568 = -t598 * pkin(1) - pkin(12);
t520 = -t550 * pkin(4) + t568;
t430 = t494 * pkin(8) - t495 * pkin(10) + t520;
t626 = t595 * mrSges(6,1) - t588 * mrSges(6,2);
t697 = t626 * t430;
t586 = sin(qJ(7));
t593 = cos(qJ(7));
t549 = -t586 * t591 + t593 * t598;
t551 = t586 * t598 + t593 * t591;
t627 = m(6) * t430 - mrSges(6,3) * t495;
t403 = -t494 * mrSges(6,2) + t627 * t588;
t404 = t494 * mrSges(6,1) + t627 * t595;
t390 = t588 * t403 + t595 * t404;
t616 = m(5) * t520 + t494 * mrSges(5,1) + t495 * mrSges(5,2) + t390;
t585 = sin(pkin(15));
t592 = cos(qJ(8));
t681 = cos(pkin(15));
t690 = sin(qJ(8));
t545 = t585 * t592 - t681 * t690;
t546 = -t585 * t690 - t681 * t592;
t452 = -t545 * t551 + t546 * t549;
t453 = t545 * t549 + t546 * t551;
t477 = (-t681 * t549 - t551 * t585) * pkin(3) + t568;
t630 = -m(9) * t477 + t452 * mrSges(9,1) - t453 * mrSges(9,2);
t695 = t550 * mrSges(4,1) - t552 * mrSges(4,2) - t616 + t549 * mrSges(8,1) - t551 * mrSges(8,2) + t630 + (-m(4) - m(8)) * t568;
t694 = 0.2e1 * t598;
t691 = pkin(10) * m(6);
t689 = pkin(1) * t590;
t687 = t585 * pkin(3);
t686 = t586 * pkin(1);
t684 = t593 * pkin(1);
t683 = t596 * pkin(4);
t680 = mrSges(5,2) * t589;
t679 = mrSges(9,3) * t545;
t678 = mrSges(9,3) * t546;
t582 = Ifges(6,4) * t588;
t583 = Ifges(6,4) * t595;
t677 = Ifges(9,3) / pkin(7) ^ 2;
t675 = t588 * mrSges(6,3);
t577 = t595 * mrSges(6,3);
t613 = -qJ(7) + t618;
t614 = -qJ(7) + t617;
t420 = ((cos(t614) - cos(t613)) * pkin(1) + (cos(t610) - cos(t612)) * pkin(2)) * pkin(3) + ((-cos((qJ(8) - t614)) + cos((qJ(8) - t613))) * pkin(1) + t696 * pkin(2)) * pkin(7);
t602 = 0.1e1 / pkin(7);
t674 = t420 * t602;
t671 = t527 * t595;
t566 = sin(t634);
t565 = 0.1e1 / t566;
t670 = (-pkin(2) * t566 - pkin(1) * sin(t666)) * t565;
t557 = sin((t633 - t665));
t575 = sin(t652);
t669 = 0.1e1 / t557 * t575;
t668 = t569 * t595;
t642 = t595 * t403 - t588 * t404;
t387 = -t494 * mrSges(5,3) + t642;
t425 = (-mrSges(5,3) - t625) * t495;
t664 = t589 * t387 + t596 * t425;
t510 = (-m(6) * t527 - mrSges(6,3)) * t588;
t511 = m(6) * t671 + t577;
t641 = -t588 * t510 + t595 * t511;
t433 = m(5) * t531 - mrSges(5,2) + t641;
t530 = t596 * t571 + t589 * t689;
t526 = -pkin(8) - t530;
t623 = -m(6) * t526 + t626;
t470 = m(5) * t530 + mrSges(5,1) + t623;
t663 = t589 * t433 + t596 * t470;
t560 = t686 + t687;
t648 = t681 * pkin(3);
t561 = t648 + t684;
t472 = -t545 * t560 + t546 * t561;
t466 = m(9) * t472 + mrSges(9,1);
t473 = t545 * t561 + t546 * t560;
t467 = m(9) * t473 - mrSges(9,2);
t662 = t546 * t466 + t545 * t467;
t661 = -t545 * t466 + t546 * t467;
t422 = Ifges(9,5) * t453 + Ifges(9,6) * t452;
t486 = (-t545 * t585 + t681 * t546) * pkin(3);
t482 = m(9) * t486 + mrSges(9,1);
t487 = (t681 * t545 + t546 * t585) * pkin(3);
t483 = m(9) * t487 - mrSges(9,2);
t660 = t546 * t482 + t545 * t483;
t659 = -t545 * t482 + t546 * t483;
t534 = (-m(6) * t569 - mrSges(6,3)) * t588;
t535 = m(6) * t668 + t577;
t640 = -t588 * t534 + t595 * t535;
t474 = m(5) * t685 - mrSges(5,2) + t640;
t570 = -pkin(8) - t683;
t622 = -m(6) * t570 + t626;
t514 = m(5) * t683 + mrSges(5,1) + t622;
t658 = t589 * t474 + t596 * t514;
t655 = Ifges(6,2) * t595 + t582;
t654 = Ifges(6,1) * t588 + t583;
t651 = t601 * t688;
t650 = Ifges(5,6) - t656;
t499 = t546 * mrSges(9,1) - t545 * mrSges(9,2);
t501 = -t545 * mrSges(9,1) - t546 * mrSges(9,2);
t649 = t499 * t648 + t501 * t687 + Ifges(9,3);
t647 = t506 * t674;
t645 = t604 * t670;
t644 = t602 * t669;
t555 = (-mrSges(6,3) - t691) * t588;
t556 = t595 * t691 + t577;
t639 = -t588 * t555 + t595 * t556;
t637 = pkin(4) * t644;
t635 = -t625 * t685 + t701;
t429 = mrSges(9,1) * t472 - mrSges(9,2) * t473 + Ifges(9,3);
t435 = mrSges(9,1) * t486 - mrSges(9,2) * t487 + Ifges(9,3);
t632 = t452 * t679 - t453 * t678;
t631 = t452 * t678 + t453 * t679;
t629 = m(6) * pkin(8) + t626;
t528 = mrSges(6,1) * pkin(8) + pkin(10) * t577 + t655;
t529 = -mrSges(6,2) * pkin(8) + pkin(10) * t675 + t654;
t431 = pkin(8) * t629 + pkin(10) * t639 + t595 * t528 + t588 * t529 + Ifges(5,3);
t624 = t595 * Ifges(6,5) - t588 * Ifges(6,6);
t505 = -mrSges(5,2) + t639;
t544 = mrSges(5,1) + t629;
t468 = t589 * t505 + t596 * t544;
t621 = pkin(4) * t468 + t431;
t508 = -mrSges(6,1) * t570 + mrSges(6,3) * t668 + t655;
t509 = mrSges(6,2) * t570 + t569 * t675 + t654;
t620 = pkin(8) * t622 + pkin(10) * t640 + mrSges(5,1) * t683 + t595 * t508 + t588 * t509 + Ifges(5,3);
t396 = t430 * t675 + Ifges(6,6) * t494 + (-Ifges(6,2) * t588 + t583) * t495;
t397 = -t430 * t577 + Ifges(6,5) * t494 + (Ifges(6,1) * t595 - t582) * t495;
t376 = pkin(10) * t642 - Ifges(5,6) * t494 + t595 * t396 + t588 * t397 + (-pkin(8) * t625 + Ifges(5,5)) * t495;
t619 = Ifges(8,5) * t551 + Ifges(8,6) * t549 + t631 * t687 + t632 * t648 + t422;
t615 = t660 * t648 + t659 * t687 + Ifges(8,3) + t435;
t611 = pkin(4) * t658 + Ifges(4,3) + t620;
t460 = -mrSges(6,1) * t526 + mrSges(6,3) * t671 + t655;
t461 = mrSges(6,2) * t526 + t527 * t675 + t654;
t392 = pkin(8) * t623 + pkin(10) * t641 + mrSges(5,1) * t530 - mrSges(5,2) * t531 + t595 * t460 + t588 * t461 + Ifges(5,3);
t609 = pkin(4) * t664 + Ifges(4,5) * t552 + Ifges(4,6) * t550 + t376;
t606 = mrSges(8,1) * t684 - mrSges(8,2) * t686 + t662 * t648 + t661 * t687 + Ifges(8,3) + t429;
t605 = pkin(4) * t663 - mrSges(4,1) * t682 + mrSges(4,2) * t689 + Ifges(4,3) + t392;
t594 = cos(qJ(6));
t587 = sin(qJ(6));
t576 = sin(t653);
t547 = -pkin(8) * t626 - Ifges(6,3);
t518 = -pkin(10) * t626 + t624;
t513 = pkin(1) * t576 / pkin(5) * t565 * (t587 * Ifges(7,5) + t594 * Ifges(7,6));
t507 = t595 * t555 + t588 * t556;
t500 = t546 * Ifges(9,5) - t545 * Ifges(9,6);
t498 = t545 * Ifges(9,5) + t546 * Ifges(9,6);
t484 = t595 * t534 + t588 * t535;
t481 = -mrSges(9,3) * t486 + Ifges(9,5);
t480 = mrSges(9,3) * t487 + Ifges(9,6);
t476 = t596 * t518 - t589 * t547;
t465 = -mrSges(9,3) * t472 + Ifges(9,5);
t464 = mrSges(9,3) * t473 + Ifges(9,6);
t463 = -pkin(8) * t507 + t650 + t700;
t458 = -pkin(4) * t626 + t589 * t518 + t596 * t547;
t446 = t595 * t510 + t588 * t511;
t438 = -pkin(10) * t507 - t588 * t528 + t595 * t529 + Ifges(5,5);
t437 = -pkin(8) * t484 + mrSges(5,3) * t685 + t650 + t698;
t428 = -t545 * t480 + t546 * t481 + Ifges(8,5);
t427 = t546 * t480 + t545 * t481 + Ifges(8,6);
t426 = -pkin(10) * t484 - mrSges(5,3) * t683 - t588 * t508 + t595 * t509 + Ifges(5,5);
t419 = t596 * t438 - t589 * t463;
t418 = -pkin(4) * t680 + t620;
t415 = -mrSges(8,3) * t684 - t545 * t464 + t546 * t465 + Ifges(8,5);
t414 = mrSges(8,3) * t686 + t546 * t464 + t545 * t465 + Ifges(8,6);
t413 = -pkin(8) * t446 + mrSges(5,3) * t531 + t650 + t699;
t409 = t591 * (-t586 * t498 + t593 * t500) + t598 * (t593 * t498 + t586 * t500);
t407 = mrSges(9,2) * t477 + Ifges(9,1) * t453 + Ifges(9,4) * t452;
t406 = -mrSges(9,1) * t477 + Ifges(9,4) * t453 + Ifges(9,2) * t452;
t405 = -pkin(4) * t507 + t589 * t438 + t596 * t463;
t399 = t596 * t426 - t589 * t437 + Ifges(4,5);
t398 = -pkin(10) * t446 - mrSges(5,3) * t530 - t588 * t460 + t595 * t461 + Ifges(5,5);
t393 = -pkin(4) * t484 + t589 * t426 + t596 * t437 + Ifges(4,6);
t391 = -t431 * t646 + t621;
t388 = mrSges(4,3) * t682 + t596 * t398 - t589 * t413 + Ifges(4,5);
t385 = pkin(1) * (t593 * t499 + t586 * t501) + (Ifges(9,3) * t647 + t649 * t670) * t604 + t649;
t384 = pkin(1) * (-t590 * (t596 * t505 - t589 * t544) - t597 * t468) + t431 * t628 + t621;
t383 = -pkin(4) * t446 - mrSges(4,3) * t689 + t589 * t398 + t596 * t413 + Ifges(4,6);
t381 = mrSges(8,2) * t568 + Ifges(8,1) * t551 + Ifges(8,4) * t549 - t545 * t406 + t546 * t407 + t630 * t687;
t380 = -mrSges(8,1) * t568 + Ifges(8,4) * t551 + Ifges(8,2) * t549 + t546 * t406 + t545 * t407 + t630 * t648;
t379 = t591 * (t590 * t405 - t597 * t419) + t598 * (-pkin(1) * t507 - t597 * t405 - t590 * t419) - pkin(12) * t507;
t378 = -pkin(8) * t390 - mrSges(5,1) * t520 + (-Ifges(5,2) - Ifges(6,3)) * t494 - t697 + (Ifges(5,4) - t624) * t495;
t377 = -pkin(10) * t390 + mrSges(5,2) * t520 + Ifges(5,1) * t495 - Ifges(5,4) * t494 - t588 * t396 + t595 * t397;
t375 = mrSges(4,2) * t568 + Ifges(4,1) * t552 + Ifges(4,4) * t550 + t596 * t377 - t589 * t378;
t374 = -pkin(4) * t616 - mrSges(4,1) * t568 + Ifges(4,4) * t552 + Ifges(4,2) * t550 + t589 * t377 + t596 * t378;
t1 = [Ifges(2,3) + t598 * (t695 * pkin(1) + Ifges(3,2) * t598 - t597 * t374 - t590 * t375 + t593 * t380 + t586 * t381) + Ifges(7,2) * t594 ^ 2 + (Ifges(3,1) * t591 + Ifges(3,4) * t694 + t590 * t374 - t597 * t375 - t586 * t380 + t593 * t381) * t591 + (Ifges(7,1) * t587 + 0.2e1 * Ifges(7,4) * t594) * t587 + (m(3) * pkin(12) + mrSges(3,1) * t694 - 0.2e1 * t591 * mrSges(3,2) + t695) * pkin(12) + (m(7) * pkin(6) - 0.2e1 * t594 * mrSges(7,1) + 0.2e1 * t587 * mrSges(7,2)) * pkin(6), t591 * (t590 * t383 - t597 * t388 - t586 * t414 + t593 * t415 + Ifges(3,5)) + t598 * (-pkin(1) * t446 - t597 * t383 - t590 * t388 + t593 * t414 + t586 * t415 + Ifges(3,6)) - pkin(12) * t446 + t513 + ((t591 * (-t586 * t427 + t593 * t428) + t598 * (t593 * t427 + t586 * t428)) * t670 + (t379 * t651 + t409 * t674) * t506) * t604, t591 * (t590 * t393 - t597 * t399) + t598 * (-pkin(1) * t484 - t597 * t393 - t590 * t399) - pkin(12) * t484 - t379 * t646 + t409 * t637, t591 * (t590 * t458 - t597 * t476) + t598 * (-pkin(1) * t626 - t597 * t458 - t590 * t476) - pkin(12) * t626; t513 + t609 + (t619 * t670 + (t376 * t651 + t422 * t674) * t506) * t604 + (t586 * (t549 * mrSges(8,3) + t631) + t593 * (-t551 * mrSges(8,3) + t632) - t590 * (t550 * mrSges(4,3) + t596 * t387 - t589 * t425) - t597 * (-t552 * mrSges(4,3) + t664)) * pkin(1) + Ifges(3,5) * t591 + Ifges(3,6) * t598 + t619, t606 + t605 + Ifges(3,3) + ((t435 * t647 + t615 * t670) * t604 + t615 + t606) * t645 + (t385 + t429) * t604 * t647 + (t384 + t392) * t628 + (t586 * (m(8) * t686 - mrSges(8,2) + t661) + t593 * (m(8) * t684 + mrSges(8,1) + t662) + (t586 * (-mrSges(8,2) + t659) + t593 * (mrSges(8,1) + t660)) * t645 - t590 * (-m(4) * t689 + t596 * t433 - t589 * t470 - mrSges(4,2)) - t597 * (-m(4) * t682 + mrSges(4,1) + t663) + t576 ^ 2 / pkin(5) ^ 2 / t566 ^ 2 * Ifges(7,3) * pkin(1)) * pkin(1), pkin(1) * (-t590 * (t596 * t474 - t589 * t514 - mrSges(4,2)) - t597 * (mrSges(4,1) + t658)) + (-t384 * t673 + t418 * t638) * t601 + (t385 * t644 - t680) * pkin(4) + t611, t702 - pkin(1) * (-t589 * t597 - t667) * t625 + t635; -t376 * t646 + t422 * t637 + t609, t605 + (t391 * t638 - t392 * t673) * t601 + (t420 * t672 * t677 + (t435 * t645 + t429) * t602) * pkin(4) * t669, -(t391 + t418) * t646 + (-t680 + t575 ^ 2 / t557 ^ 2 * pkin(4) * t677) * pkin(4) + t611, t635 - t703; Ifges(6,3) * t494 + t624 * t495 + t697, t656 - t699 + t702, t656 - t698 - t703, Ifges(6,3);];
Mq = t1;
