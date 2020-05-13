% Calculate joint inertia matrix with Newton Euler for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_inertiaJ_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_inertiaJ_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_inertiaJ_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_inertiaJ_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1IC_inertiaJ_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:32:10
% EndTime: 2020-04-20 17:32:13
% DurationCPUTime: 3.32s
% Computational Cost: add. (5784->413), mult. (8440->519), div. (128->6), fcn. (8764->24), ass. (0->215)
t578 = -qJ(7) + pkin(15);
t568 = -qJ(8) + t578;
t556 = sin(t568);
t524 = -t556 * pkin(7) + pkin(3) * sin(t578);
t558 = cos(t568);
t525 = -t558 * pkin(7) + pkin(3) * cos(t578);
t569 = pkin(16) + qJ(7) + qJ(2);
t557 = sin(t569);
t559 = cos(t569);
t581 = sin(qJ(6));
t588 = cos(qJ(6));
t503 = 0.1e1 / (-t557 * t588 + t559 * t581) / pkin(5) / pkin(2);
t585 = sin(qJ(2));
t543 = -pkin(1) * t585 - pkin(2) * t557;
t592 = cos(qJ(2));
t644 = t592 * pkin(1);
t544 = pkin(2) * t559 + t644;
t423 = (t543 * t588 + t544 * t581) * pkin(5) * t503;
t567 = qJ(3) + qJ(4) + pkin(14);
t553 = sin(t567);
t554 = cos(t567);
t488 = 0.1e1 / (t553 * t558 + t554 * t556) / pkin(9) / pkin(7);
t599 = t488 * t423;
t379 = (-t524 * t558 + t525 * t556) * pkin(7) * t599;
t582 = sin(qJ(5));
t589 = cos(qJ(5));
t627 = Ifges(6,5) * t582 + Ifges(6,6) * t589;
t609 = mrSges(6,1) * t582 + mrSges(6,2) * t589;
t660 = t609 * pkin(10);
t661 = t627 - t660;
t663 = t379 * t661;
t584 = sin(qJ(3));
t541 = -pkin(4) * t584 - pkin(9) * t553;
t591 = cos(qJ(3));
t542 = pkin(4) * t591 + pkin(9) * t554;
t415 = (t541 * t558 - t542 * t556) * t488 * pkin(7);
t662 = t415 * t661;
t650 = pkin(1) * t591;
t563 = pkin(4) - t650;
t583 = sin(qJ(4));
t590 = cos(qJ(4));
t635 = t584 * t590;
t520 = -pkin(1) * t635 + t583 * t563;
t516 = pkin(10) + t520;
t659 = t609 * t516;
t648 = pkin(4) * t583;
t561 = pkin(10) + t648;
t658 = t609 * t561;
t538 = t584 * t585 - t591 * t592;
t540 = -t584 * t592 - t585 * t591;
t484 = -t590 * t538 + t540 * t583;
t485 = t538 * t583 + t540 * t590;
t560 = -pkin(12) - t644;
t510 = -pkin(4) * t538 + t560;
t418 = pkin(8) * t484 - pkin(10) * t485 + t510;
t610 = t589 * mrSges(6,1) - mrSges(6,2) * t582;
t657 = t610 * t418;
t580 = sin(qJ(7));
t587 = cos(qJ(7));
t537 = -t580 * t585 + t587 * t592;
t539 = t580 * t592 + t585 * t587;
t611 = m(6) * t418 - mrSges(6,3) * t485;
t389 = -t484 * mrSges(6,2) + t611 * t582;
t390 = t484 * mrSges(6,1) + t611 * t589;
t374 = t582 * t389 + t589 * t390;
t602 = m(5) * t510 + t484 * mrSges(5,1) + t485 * mrSges(5,2) + t374;
t579 = sin(pkin(15));
t586 = cos(qJ(8));
t643 = cos(pkin(15));
t653 = sin(qJ(8));
t534 = t579 * t586 - t643 * t653;
t535 = -t579 * t653 - t643 * t586;
t442 = -t534 * t539 + t535 * t537;
t443 = t534 * t537 + t535 * t539;
t469 = (-t643 * t537 - t539 * t579) * pkin(3) + t560;
t613 = -m(9) * t469 + t442 * mrSges(9,1) - t443 * mrSges(9,2);
t656 = t538 * mrSges(4,1) - t540 * mrSges(4,2) - t602 + t537 * mrSges(8,1) - t539 * mrSges(8,2) + t613 + (-m(4) - m(8)) * t560;
t655 = 0.2e1 * t592;
t654 = pkin(10) * m(6);
t652 = pkin(1) * t580;
t651 = pkin(1) * t584;
t649 = pkin(3) * t579;
t647 = pkin(4) * t590;
t645 = t587 * pkin(1);
t641 = mrSges(6,3) * t582;
t640 = mrSges(9,3) * t534;
t639 = mrSges(9,3) * t535;
t575 = Ifges(6,4) * t582;
t576 = Ifges(6,4) * t589;
t378 = (-t524 * t554 - t525 * t553) * pkin(9) * t599;
t638 = t378 * Ifges(9,3);
t570 = t589 * mrSges(6,3);
t637 = t516 * t589;
t636 = t561 * t589;
t621 = t589 * t389 - t390 * t582;
t371 = -mrSges(5,3) * t484 + t621;
t410 = (-mrSges(5,3) - t609) * t485;
t634 = t583 * t371 + t590 * t410;
t500 = (-m(6) * t516 - mrSges(6,3)) * t582;
t501 = m(6) * t637 + t570;
t619 = -t500 * t582 + t589 * t501;
t422 = m(5) * t520 - mrSges(5,2) + t619;
t519 = t563 * t590 + t583 * t651;
t515 = -pkin(8) - t519;
t606 = -m(6) * t515 + t610;
t460 = m(5) * t519 + mrSges(5,1) + t606;
t633 = t583 * t422 + t590 * t460;
t550 = t649 + t652;
t622 = t643 * pkin(3);
t551 = t622 + t645;
t462 = -t534 * t550 + t535 * t551;
t456 = m(9) * t462 + mrSges(9,1);
t463 = t534 * t551 + t535 * t550;
t457 = m(9) * t463 - mrSges(9,2);
t632 = t535 * t456 + t534 * t457;
t631 = -t534 * t456 + t535 * t457;
t407 = Ifges(9,5) * t443 + Ifges(9,6) * t442;
t477 = (-t534 * t579 + t643 * t535) * pkin(3);
t473 = m(9) * t477 + mrSges(9,1);
t478 = (t643 * t534 + t535 * t579) * pkin(3);
t474 = m(9) * t478 - mrSges(9,2);
t630 = t535 * t473 + t534 * t474;
t629 = -t534 * t473 + t535 * t474;
t522 = (-m(6) * t561 - mrSges(6,3)) * t582;
t523 = m(6) * t636 + t570;
t618 = -t522 * t582 + t589 * t523;
t466 = m(5) * t648 - mrSges(5,2) + t618;
t562 = -pkin(8) - t647;
t605 = -m(6) * t562 + t610;
t504 = m(5) * t647 + mrSges(5,1) + t605;
t628 = t583 * t466 + t590 * t504;
t626 = Ifges(6,2) * t589 + t575;
t625 = Ifges(6,1) * t582 + t576;
t624 = Ifges(5,6) - t627;
t490 = mrSges(9,1) * t535 - mrSges(9,2) * t534;
t492 = -mrSges(9,1) * t534 - mrSges(9,2) * t535;
t623 = t490 * t622 + t492 * t649 + Ifges(9,3);
t416 = mrSges(9,1) * t462 - mrSges(9,2) * t463 + Ifges(9,3);
t425 = mrSges(9,1) * t477 - mrSges(9,2) * t478 + Ifges(9,3);
t620 = -t423 * t425 + t416;
t546 = (-mrSges(6,3) - t654) * t582;
t547 = t589 * t654 + t570;
t617 = -t546 * t582 + t589 * t547;
t616 = -t609 * t648 + t661;
t615 = t442 * t640 - t443 * t639;
t614 = t442 * t639 + t443 * t640;
t612 = m(6) * pkin(8) + t610;
t517 = mrSges(6,1) * pkin(8) + pkin(10) * t570 + t626;
t518 = -mrSges(6,2) * pkin(8) + pkin(10) * t641 + t625;
t420 = pkin(8) * t612 + pkin(10) * t617 + t589 * t517 + t582 * t518 + Ifges(5,3);
t608 = Ifges(6,5) * t589 - Ifges(6,6) * t582;
t496 = -mrSges(5,2) + t617;
t533 = mrSges(5,1) + t612;
t458 = t496 * t583 + t533 * t590;
t604 = pkin(4) * t458 + t420;
t382 = t418 * t641 + Ifges(6,6) * t484 + (-Ifges(6,2) * t582 + t576) * t485;
t383 = -t418 * t570 + Ifges(6,5) * t484 + (Ifges(6,1) * t589 - t575) * t485;
t360 = pkin(10) * t621 - Ifges(5,6) * t484 + t589 * t382 + t582 * t383 + (-pkin(8) * t609 + Ifges(5,5)) * t485;
t603 = Ifges(8,5) * t539 + Ifges(8,6) * t537 + t614 * t649 + t615 * t622 + t407;
t600 = mrSges(8,1) * t645 + t632 * t622 + t631 * t649 + Ifges(8,3) + t416;
t450 = -mrSges(6,1) * t515 + mrSges(6,3) * t637 + t626;
t451 = mrSges(6,2) * t515 + t516 * t641 + t625;
t376 = pkin(8) * t606 + pkin(10) * t619 + mrSges(5,1) * t519 - mrSges(5,2) * t520 + t589 * t450 + t582 * t451 + Ifges(5,3);
t598 = pkin(4) * t634 + Ifges(4,5) * t540 + Ifges(4,6) * t538 + t360;
t498 = -mrSges(6,1) * t562 + mrSges(6,3) * t636 + t626;
t499 = mrSges(6,2) * t562 + t561 * t641 + t625;
t404 = pkin(8) * t605 + pkin(10) * t618 + mrSges(5,1) * t647 - mrSges(5,2) * t648 + t589 * t498 + t582 * t499 + Ifges(5,3);
t596 = pkin(4) * t628 + Ifges(4,3) + t404;
t595 = pkin(4) * t633 + mrSges(4,2) * t651 + Ifges(4,3) + t376;
t536 = -pkin(8) * t610 - Ifges(6,3);
t508 = -pkin(10) * t610 + t608;
t497 = t546 * t589 + t547 * t582;
t491 = Ifges(9,5) * t535 - Ifges(9,6) * t534;
t489 = Ifges(9,5) * t534 + Ifges(9,6) * t535;
t475 = t522 * t589 + t523 * t582;
t472 = -mrSges(9,3) * t477 + Ifges(9,5);
t471 = mrSges(9,3) * t478 + Ifges(9,6);
t468 = t508 * t590 - t536 * t583;
t455 = -mrSges(9,3) * t462 + Ifges(9,5);
t454 = mrSges(9,3) * t463 + Ifges(9,6);
t453 = -pkin(8) * t497 + t624 + t660;
t448 = -pkin(4) * t610 + t508 * t583 + t536 * t590;
t436 = t500 * t589 + t501 * t582;
t428 = -pkin(10) * t497 - t517 * t582 + t518 * t589 + Ifges(5,5);
t427 = -pkin(8) * t475 + mrSges(5,3) * t648 + t624 + t658;
t419 = (-t543 * t559 - t544 * t557) * t503 * pkin(2);
t417 = t419 * (Ifges(7,5) * t581 + Ifges(7,6) * t588);
t414 = (t541 * t554 + t542 * t553) * t488 * pkin(9);
t413 = -t471 * t534 + t472 * t535 + Ifges(8,5);
t412 = t471 * t535 + t472 * t534 + Ifges(8,6);
t411 = -pkin(10) * t475 - mrSges(5,3) * t647 - t498 * t582 + t499 * t589 + Ifges(5,5);
t405 = t428 * t590 - t453 * t583;
t401 = -mrSges(8,3) * t645 - t454 * t534 + t455 * t535 + Ifges(8,5);
t400 = mrSges(8,3) * t652 + t454 * t535 + t455 * t534 + Ifges(8,6);
t399 = -pkin(8) * t436 + mrSges(5,3) * t520 + t624 + t659;
t395 = t585 * (-t489 * t580 + t491 * t587) + t592 * (t489 * t587 + t491 * t580);
t393 = mrSges(9,2) * t469 + Ifges(9,1) * t443 + Ifges(9,4) * t442;
t392 = -mrSges(9,1) * t469 + Ifges(9,4) * t443 + Ifges(9,2) * t442;
t391 = -pkin(4) * t497 + t428 * t583 + t453 * t590;
t385 = t411 * t590 - t427 * t583 + Ifges(4,5);
t384 = -pkin(10) * t436 - mrSges(5,3) * t519 - t450 * t582 + t451 * t589 + Ifges(5,5);
t377 = -pkin(4) * t475 + t411 * t583 + t427 * t590 + Ifges(4,6);
t375 = t415 * t420 + t604;
t372 = mrSges(4,3) * t650 + t384 * t590 - t399 * t583 + Ifges(4,5);
t369 = -pkin(4) * t436 - mrSges(4,3) * t651 + t384 * t583 + t399 * t590 + Ifges(4,6);
t367 = pkin(1) * (t490 * t587 + t492 * t580) - t423 * t623 + t638 + t623;
t366 = mrSges(8,2) * t560 + Ifges(8,1) * t539 + Ifges(8,4) * t537 - t392 * t534 + t393 * t535 + t613 * t649;
t365 = -mrSges(8,1) * t560 + Ifges(8,4) * t539 + Ifges(8,2) * t537 + t535 * t392 + t534 * t393 + t613 * t622;
t364 = pkin(1) * (-t584 * (t496 * t590 - t533 * t583) - t591 * t458) + t379 * t420 + t604;
t363 = t585 * (t391 * t584 - t405 * t591) + t592 * (-pkin(1) * t497 - t391 * t591 - t405 * t584) - pkin(12) * t497;
t362 = -pkin(8) * t374 - mrSges(5,1) * t510 + (-Ifges(5,2) - Ifges(6,3)) * t484 - t657 + (Ifges(5,4) - t608) * t485;
t361 = -pkin(10) * t374 + mrSges(5,2) * t510 + Ifges(5,1) * t485 - Ifges(5,4) * t484 - t382 * t582 + t383 * t589;
t359 = mrSges(4,2) * t560 + Ifges(4,1) * t540 + Ifges(4,4) * t538 + t361 * t590 - t362 * t583;
t358 = -pkin(4) * t602 - mrSges(4,1) * t560 + Ifges(4,4) * t540 + Ifges(4,2) * t538 + t583 * t361 + t590 * t362;
t1 = [Ifges(2,3) + t592 * (t656 * pkin(1) + Ifges(3,2) * t592 - t591 * t358 - t584 * t359 + t587 * t365 + t580 * t366) + Ifges(7,2) * t588 ^ 2 + (Ifges(3,1) * t585 + Ifges(3,4) * t655 + t358 * t584 - t359 * t591 - t365 * t580 + t366 * t587) * t585 + (Ifges(7,1) * t581 + 0.2e1 * Ifges(7,4) * t588) * t581 + (m(3) * pkin(12) + mrSges(3,1) * t655 - 0.2e1 * t585 * mrSges(3,2) + t656) * pkin(12) + (m(7) * pkin(6) - 0.2e1 * mrSges(7,1) * t588 + 0.2e1 * mrSges(7,2) * t581) * pkin(6), t585 * (t584 * t369 - t591 * t372 - t580 * t400 + t587 * t401 + Ifges(3,5)) + t592 * (-pkin(1) * t436 - t591 * t369 - t584 * t372 + t587 * t400 + t580 * t401 + Ifges(3,6)) - pkin(12) * t436 + t363 * t379 + t417 - (t585 * (-t412 * t580 + t413 * t587) + t592 * (t412 * t587 + t413 * t580)) * t423 + t395 * t378, t585 * (t377 * t584 - t385 * t591) + t592 * (-pkin(1) * t475 - t377 * t591 - t385 * t584) - pkin(12) * t475 + t363 * t415 + t395 * t414, t585 * (t448 * t584 - t468 * t591) + t592 * (-pkin(1) * t610 - t448 * t591 - t468 * t584) - pkin(12) * t610; t598 + t378 * t407 + t379 * t360 + t417 - t423 * t603 + (t580 * (mrSges(8,3) * t537 + t614) + t587 * (-mrSges(8,3) * t539 + t615) - t584 * (mrSges(4,3) * t538 + t371 * t590 - t410 * t583) - t591 * (-mrSges(4,3) * t540 + t634)) * pkin(1) + Ifges(3,5) * t585 + Ifges(3,6) * t592 + t603, t595 + (t367 + t620) * t378 + t419 ^ 2 * Ifges(7,3) + (t364 + t376) * t379 - t423 * t600 + Ifges(3,3) + (-t584 * (t422 * t590 - t460 * t583 - mrSges(4,2)) + (-0.2e1 * mrSges(4,1) - t633) * t591 + (t584 ^ 2 + t591 ^ 2) * pkin(1) * m(4) + (mrSges(8,1) - (mrSges(8,1) + t630) * t423 + m(8) * t645 + t632) * t587 + (-t629 * t423 + m(8) * t652 + (0.2e1 * t423 - 0.2e1) * mrSges(8,2) + t631) * t580) * pkin(1) + t600 - (0.1e1 - t423) * t423 * (t630 * t622 + t629 * t649 + Ifges(8,3) + t425), pkin(1) * (-t584 * (t466 * t590 - t504 * t583 - mrSges(4,2)) - t591 * (mrSges(4,1) + t628)) + t379 * t404 + t364 * t415 + t367 * t414 + t596, t663 - pkin(1) * (-t583 * t591 - t635) * t609 + t616; t360 * t415 + t407 * t414 + t598, t595 - mrSges(4,1) * t650 + t375 * t379 + t376 * t415 + (t620 + t638) * t414, Ifges(9,3) * t414 ^ 2 + (t375 + t404) * t415 + t596, t616 + t662; Ifges(6,3) * t484 + t608 * t485 + t657, t627 - t659 + t663, t627 - t658 + t662, Ifges(6,3);];
Mq = t1;
