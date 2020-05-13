% Calculate vector of centrifugal and Coriolis load on the joints for
% palh3m2DE2
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh3m2DE2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:16
% EndTime: 2020-05-07 04:22:22
% DurationCPUTime: 5.52s
% Computational Cost: add. (9404->340), mult. (18210->544), div. (0->0), fcn. (22369->20), ass. (0->194)
t562 = sin(qJ(3));
t669 = pkin(1) * qJD(2);
t612 = t562 * t669;
t554 = pkin(17) + pkin(18);
t549 = sin(t554);
t550 = cos(t554);
t563 = sin(qJ(2));
t566 = cos(qJ(2));
t558 = cos(pkin(16));
t567 = cos(pkin(15));
t647 = sin(pkin(16));
t675 = sin(pkin(15));
t527 = t558 * t675 + t647 * t567;
t528 = t558 * t567 - t647 * t675;
t565 = cos(qJ(3));
t586 = t527 * t562 - t528 * t565;
t604 = -t527 * t565 - t562 * t528;
t606 = t563 * t586 + t604 * t566;
t710 = t586 * t566;
t695 = -t563 * t604 + t710;
t714 = t549 * t695 + t550 * t606;
t440 = t606 * t549 - t695 * t550;
t553 = qJD(2) + qJD(3);
t672 = pkin(4) * t553;
t582 = t565 * t669 - t672;
t719 = t440 * t582;
t720 = -t714 * t612 + t719;
t611 = t440 * t669;
t436 = t562 * t611;
t715 = t714 * t582;
t423 = -t436 - t715;
t480 = t604 * qJD(3);
t717 = t695 * qJD(2) - t563 * t480;
t479 = t586 * qJD(3);
t441 = -t479 * t566 - t717;
t442 = t606 * qJD(2) + t479 * t563 + t480 * t566;
t431 = -t441 * t549 + t442 * t550;
t547 = -pkin(1) * t565 + pkin(4);
t637 = t442 * t549;
t638 = t714 * t565;
t577 = (-qJD(3) * t638 + (t441 * t550 + t637) * t562) * pkin(1);
t674 = pkin(1) * t562;
t613 = t440 * t674;
t623 = -qJD(3) * t613 - t431 * t547 + t423 + t577;
t713 = mrSges(9,1) + mrSges(4,1);
t712 = mrSges(4,2) + mrSges(9,2);
t568 = cos(pkin(14));
t676 = sin(pkin(14));
t531 = -t567 * t676 + t675 * t568;
t534 = t567 * t568 + t675 * t676;
t496 = t531 * t566 + t534 * t563;
t711 = t496 * Ifges(7,4);
t560 = cos(pkin(18));
t670 = sin(pkin(18));
t530 = t560 * t567 - t675 * t670;
t557 = qJ(2) + qJ(3);
t551 = sin(t557);
t552 = cos(t557);
t576 = -t675 * t560 - t567 * t670;
t708 = mrSges(8,1) * t530 + mrSges(9,1) * t552 - mrSges(8,2) * t576 - mrSges(9,2) * t551;
t585 = t551 * t562 + t552 * t565;
t471 = t527 * t550 + t528 * t549;
t466 = t471 * qJD(1);
t472 = t527 * t549 - t528 * t550;
t467 = t472 * qJD(1);
t533 = t562 * t563 - t565 * t566;
t523 = t533 * qJD(1);
t545 = -t566 * pkin(1) - pkin(12);
t503 = -pkin(4) * t523 + t545 * qJD(1);
t445 = -pkin(8) * t467 - pkin(10) * t466 + t503;
t561 = sin(qJ(4));
t564 = cos(qJ(4));
t417 = -t423 * t561 + t445 * t564;
t418 = t423 * t564 + t445 * t561;
t706 = (-t417 * t564 - t418 * t561) * qJD(4);
t603 = -t531 * t563 + t534 * t566;
t687 = -t603 / 0.2e1;
t705 = Ifges(7,4) * t603;
t464 = (cos(pkin(17)) * t530 + t576 * sin(pkin(17))) * pkin(3) + t545;
t704 = t464 * (-mrSges(9,1) * t551 - mrSges(9,2) * t552);
t584 = t562 * t566 + t563 * t565;
t524 = t584 * qJD(1);
t699 = (-mrSges(4,1) * t523 - mrSges(4,2) * t524 + qJD(1) * t708) * t563;
t616 = t551 * qJD(1);
t667 = mrSges(4,3) * t524;
t698 = (mrSges(9,3) * t616 + t553 * t713 + t667) * t562;
t697 = (-m(8) - m(4)) * t545;
t430 = (qJD(3) * t710 + t717) * t550 - t637;
t574 = -t431 * t562 + (-t440 * t565 + t562 * t714) * qJD(3);
t412 = pkin(1) * t574 + t430 * t547;
t696 = t412 - t720;
t609 = qJD(3) * t669;
t601 = t562 * t609;
t409 = qJD(2) * t577 + t431 * t582 - t440 * t601;
t426 = -t440 * t547 - t674 * t714;
t694 = t409 * t426 + t623 * t720;
t591 = -t417 * t561 + t418 * t564;
t410 = t430 * t672 + (-t430 * t565 + t574) * t669;
t500 = t553 * t584;
t484 = t500 * qJD(1);
t610 = qJD(1) * t669;
t463 = -pkin(4) * t484 + t563 * t610;
t407 = qJD(4) * t417 + t410 * t564 + t463 * t561;
t408 = -qJD(4) * t418 - t410 * t561 + t463 * t564;
t693 = t407 * t564 - t408 * t561;
t692 = -0.2e1 * pkin(12);
t691 = qJD(1) ^ 2;
t485 = t496 * qJD(2);
t689 = -0.3e1 / 0.2e1 * t485;
t486 = t603 * qJD(2);
t688 = 0.3e1 / 0.2e1 * t486;
t686 = t496 / 0.2e1;
t499 = t553 * t533;
t685 = t499 / 0.2e1;
t684 = t500 / 0.2e1;
t682 = t523 / 0.2e1;
t681 = -t524 / 0.2e1;
t679 = t561 / 0.2e1;
t678 = -t563 / 0.2e1;
t677 = t566 / 0.2e1;
t673 = pkin(4) * t524;
t671 = m(9) * t464;
t668 = mrSges(4,3) * t523;
t666 = mrSges(6,3) * t466;
t665 = Ifges(3,4) * t563;
t664 = Ifges(3,4) * t566;
t663 = Ifges(6,4) * t561;
t662 = Ifges(6,4) * t564;
t661 = Ifges(9,4) * t551;
t660 = Ifges(9,4) * t552;
t465 = qJD(4) - t467;
t659 = Ifges(6,5) * t465;
t658 = Ifges(7,5) * t486;
t657 = Ifges(9,5) * t552;
t656 = Ifges(9,5) * t553;
t655 = Ifges(6,6) * t465;
t654 = Ifges(7,6) * t485;
t653 = Ifges(9,6) * t551;
t652 = Ifges(9,6) * t553;
t651 = t496 * Ifges(7,1);
t650 = t603 * Ifges(7,2);
t649 = t524 * Ifges(4,4);
t648 = t566 * Ifges(3,2);
t416 = t715 - t436 + (t440 * t562 - t638) * t669;
t644 = t416 * t720;
t624 = t564 * t466;
t620 = qJD(1) * t552;
t622 = mrSges(9,3) * t620 + t553 * t712 - t668;
t583 = -pkin(4) * t533 + t545;
t446 = -pkin(8) * t472 - pkin(10) * t471 + t583;
t618 = qJD(4) * t446;
t617 = qJD(4) * t466;
t615 = t563 * qJD(1);
t607 = t472 * t617;
t599 = t408 * mrSges(6,1) - t407 * mrSges(6,2);
t598 = mrSges(6,1) * t564 - mrSges(6,2) * t561;
t592 = -t409 * t440 - t431 * t720;
t449 = -mrSges(6,2) * t465 - t561 * t666;
t450 = mrSges(6,1) * t465 - mrSges(6,3) * t624;
t590 = t564 * t449 - t561 * t450;
t493 = t530 * t563 - t566 * t576;
t481 = t493 * qJD(2);
t492 = t530 * t566 + t563 * t576;
t482 = t492 * qJD(2);
t587 = t481 * t576 - t482 * t530;
t581 = Ifges(3,5) * t677 + Ifges(3,6) * t678;
t548 = qJD(1) * t664;
t580 = (Ifges(3,6) * qJD(2) + (t648 + t665) * qJD(1)) * t678 + (Ifges(3,1) * t615 + Ifges(3,5) * qJD(2) + t548) * t677;
t578 = t585 * t553;
t575 = -qJD(1) * t704 - t545 * (-mrSges(4,1) * t524 + mrSges(4,2) * t523) - t552 * t656;
t573 = qJD(4) * (-t449 * t561 - t450 * t564 + (-t561 ^ 2 - t564 ^ 2) * t666);
t469 = t523 * Ifges(4,2) + t553 * Ifges(4,6) - t649;
t513 = Ifges(4,4) * t523;
t470 = -t524 * Ifges(4,1) + t553 * Ifges(4,5) + t513;
t483 = t499 * qJD(1);
t511 = t652 + (-Ifges(9,2) * t552 - t661) * qJD(1);
t512 = t656 + (-Ifges(9,1) * t551 - t660) * qJD(1);
t572 = t469 * t681 + t524 * (Ifges(4,1) * t523 + t649) / 0.2e1 + Ifges(4,6) * t484 + Ifges(4,5) * t483 + t512 * t620 / 0.2e1 + t612 * t667 - (Ifges(4,2) * t524 + t470 + t513) * t523 / 0.2e1 - (Ifges(4,5) * t523 + Ifges(4,6) * t524 + (t653 - t657) * qJD(1)) * t553 / 0.2e1 + (-t511 / 0.2e1 + t652) * t616 + (t551 * (-Ifges(9,1) * t552 + t661) + t552 * (Ifges(9,2) * t551 - t660)) * t691 / 0.2e1 + t585 * mrSges(9,3) * t610 + t713 * t601 + t712 * t565 * t609;
t506 = pkin(1) * t615 - t673;
t473 = -pkin(4) * t500 + t563 * t669;
t454 = Ifges(7,5) * qJD(2) + qJD(1) * (t651 + t705);
t453 = Ifges(7,6) * qJD(2) + qJD(1) * (t650 + t711);
t451 = -mrSges(5,1) * t467 + mrSges(5,2) * t466;
t448 = (mrSges(6,1) * t561 + mrSges(6,2) * t564) * t466;
t447 = t598 * t617;
t444 = t659 + (t564 * Ifges(6,1) - t663) * t466;
t443 = t655 + (-t561 * Ifges(6,2) + t662) * t466;
t427 = t547 * t714 - t613;
t420 = t506 * t561 + t564 * t720;
t419 = t506 * t564 - t561 * t720;
t415 = -t565 * t611 + t719;
t414 = t415 * t564 - t561 * t673;
t413 = -t415 * t561 - t564 * t673;
t1 = [-t485 * t453 / 0.2e1 + t486 * t454 / 0.2e1 + t469 * t684 + t473 * t451 + t470 * t685 + t545 * (-mrSges(4,1) * t484 + mrSges(4,2) * t483) + (t463 * mrSges(5,2) + t409 * mrSges(5,3)) * t471 + (t533 * t484 + t500 * t682) * Ifges(4,2) + (-t483 * t584 + t499 * t681) * Ifges(4,1) + (t463 * t583 + t503 * t473) * m(5) + (Ifges(4,5) * t685 + Ifges(4,6) * t684 + t551 * t511 / 0.2e1 - t552 * t512 / 0.2e1 + (-t657 / 0.2e1 + t653 / 0.2e1) * t553) * t553 + (-t463 * mrSges(5,1) + t410 * mrSges(5,3) - t599) * t472 + (t533 * t483 - t484 * t584 + t499 * t682 + t500 * t681) * Ifges(4,4) + (t658 / 0.2e1 - t654 / 0.2e1 + t581 * qJD(2) + (t587 * mrSges(8,3) + t699 + (t585 * qJD(3) - t578) * mrSges(9,3) + (t499 * t565 - t500 * t562 + (-t533 * t565 + t562 * t584) * qJD(3)) * mrSges(4,3)) * pkin(1) + t580) * qJD(2) + (t449 * t618 + t473 * t450 + m(6) * (t408 * t446 + t417 * t473 + t418 * t618) + Ifges(6,6) * t607 + (t409 * mrSges(6,2) - t408 * mrSges(6,3) + (t720 * mrSges(6,1) - t418 * mrSges(6,3) - t655 / 0.2e1 - t443 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,4) * t624) * qJD(4)) * t471) * t564 + (Ifges(6,5) * t607 + t473 * t449 - t450 * t618 + m(6) * (t407 * t446 - t417 * t618 + t418 * t473) + (t409 * mrSges(6,1) - t407 * mrSges(6,3) + (-t720 * mrSges(6,2) - t659 / 0.2e1 + t417 * mrSges(6,3) - t444 / 0.2e1 + (0.3e1 / 0.2e1 * t663 + (-0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t564) * t466) * qJD(4)) * t471) * t561 + (t545 * (-mrSges(4,1) * t500 + mrSges(4,2) * t499) + t651 * t688 + t650 * t689 + 0.2e1 * pkin(6) * (mrSges(7,1) * t485 + mrSges(7,2) * t486) + (t496 * t689 + t603 * t688) * Ifges(7,4) + ((mrSges(3,2) * t692 + 0.3e1 / 0.2e1 * t664) * t566 + (mrSges(3,1) * t692 - 0.3e1 / 0.2e1 * t665 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t566 + (-mrSges(4,1) * t533 - mrSges(4,2) * t584 + 0.2e1 * t671 - 0.2e1 * t697 + t708) * pkin(1)) * t563) * qJD(2) + (0.2e1 * t704 + (0.3e1 / 0.2e1 * Ifges(9,1) - 0.3e1 / 0.2e1 * Ifges(9,2)) * t551 * t552 + (0.3e1 / 0.2e1 * t552 ^ 2 - 0.3e1 / 0.2e1 * t551 ^ 2) * Ifges(9,4)) * t553) * qJD(1); -t506 * t451 + (t623 * t466 + t696 * t467) * mrSges(5,3) + t426 * t447 - t420 * t449 - t419 * t450 + t590 * t412 + (-t566 * t548 / 0.2e1 - t654 + t658 + t454 * t687 + t453 * t686 + (t563 * t648 / 0.2e1 + pkin(12) * (mrSges(3,1) * t563 + mrSges(3,2) * t566) + (Ifges(3,1) * t566 - t665) * t678 - pkin(6) * (mrSges(7,1) * t496 + mrSges(7,2) * t603) - t496 * (Ifges(7,1) * t603 - t711) / 0.2e1 + (-Ifges(7,2) * t496 + t705) * t687) * qJD(1) + (Ifges(7,5) * t687 + Ifges(7,6) * t686 + t581) * qJD(2) + t575 - t580) * qJD(1) + t623 * t448 + t427 * t573 + t572 + (t591 * t412 + (t706 + t693) * t427 - t417 * t419 - t418 * t420 + t694) * m(6) + (t410 * t427 + t696 * t423 - t503 * t506 + t694) * m(5) + ((-t671 + t697) * t691 * t563 + (-t562 * t484 + (-qJD(2) * t523 + t483) * t565) * mrSges(4,3) + (t622 * t565 + t698) * qJD(3) + (-mrSges(9,3) * t578 - t699 + ((t492 * t530 - t493 * t576) * qJD(2) + t587) * mrSges(8,3)) * qJD(1) + 0.2e1 * m(8) * (-t481 * t492 + t482 * t493) * t669) * pkin(1); -m(6) * (t413 * t417 + t414 * t418 + t644) - m(5) * (t415 * t423 + t644) + (-t440 * t447 + t524 * t451 + (-t466 * mrSges(5,3) - t448) * t431 + m(6) * (t706 * t714 + t592) + (m(6) * t693 + t573) * t714 + (m(6) * t591 + t467 * mrSges(5,3) + t590) * t430 + (t410 * t714 + t423 * t430 + t503 * t524 + t592) * m(5)) * pkin(4) - t416 * t448 - t414 * t449 - t413 * t450 + (-t415 * t467 - t416 * t466) * mrSges(5,3) + t572 + (-t698 + (-t622 - t668) * t565) * t669 + t575 * qJD(1); -t417 * t449 + t418 * t450 + (-t720 * t598 + t564 * t443 / 0.2e1 + t444 * t679 + ((-Ifges(6,2) * t564 - t663) * t679 - t564 * (-Ifges(6,1) * t561 - t662) / 0.2e1) * t466 + t591 * mrSges(6,3) + (-t465 / 0.2e1 + qJD(4)) * (-Ifges(6,5) * t561 - Ifges(6,6) * t564)) * t466 + t599;];
tauc = t1(:);
