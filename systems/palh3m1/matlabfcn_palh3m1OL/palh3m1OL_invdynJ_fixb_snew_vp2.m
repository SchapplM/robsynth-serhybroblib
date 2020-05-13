% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tauJ [10x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = palh3m1OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1OL_invdynJ_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:05:20
% EndTime: 2020-04-20 17:05:35
% DurationCPUTime: 3.72s
% Computational Cost: add. (20116->412), mult. (44327->543), div. (0->0), fcn. (34142->18), ass. (0->162)
t634 = sin(qJ(7));
t639 = sin(qJ(2));
t642 = cos(qJ(7));
t647 = cos(qJ(2));
t603 = (t647 * t634 + t639 * t642) * qJD(1);
t665 = qJD(1) * qJD(2);
t612 = t639 * qJDD(1) + t647 * t665;
t663 = t639 * t665;
t614 = t647 * qJDD(1) - t663;
t560 = -t603 * qJD(7) - t634 * t612 + t642 * t614;
t638 = sin(qJ(3));
t646 = cos(qJ(3));
t604 = (-t647 * t638 - t639 * t646) * qJD(1);
t561 = -t604 * qJD(3) + t638 * t612 - t646 * t614;
t601 = (-t639 * t634 + t647 * t642) * qJD(1);
t562 = t601 * qJD(7) + t642 * t612 + t634 * t614;
t602 = (t639 * t638 - t647 * t646) * qJD(1);
t563 = t602 * qJD(3) - t646 * t612 - t638 * t614;
t640 = sin(qJ(1));
t648 = cos(qJ(1));
t662 = t640 * g(1) - t648 * g(2);
t607 = -qJDD(1) * pkin(12) - t662;
t580 = t607 + (-t614 + t663) * pkin(1);
t630 = qJD(2) + qJD(7);
t585 = -t630 * mrSges(8,2) + t601 * mrSges(8,3);
t631 = qJD(2) + qJD(3);
t586 = -t631 * mrSges(4,2) + t602 * mrSges(4,3);
t587 = t630 * mrSges(8,1) - t603 * mrSges(8,3);
t588 = t631 * mrSges(4,1) - t604 * mrSges(4,3);
t637 = sin(qJ(4));
t645 = cos(qJ(4));
t575 = t637 * t602 + t645 * t604;
t517 = -t575 * qJD(4) + t645 * t561 - t637 * t563;
t574 = t645 * t602 - t637 * t604;
t518 = t574 * qJD(4) + t637 * t561 + t645 * t563;
t531 = t580 + (t604 * t631 - t561) * pkin(4);
t625 = qJD(4) + t631;
t489 = (-t574 * t625 - t518) * pkin(10) + (t575 * t625 - t517) * pkin(8) + t531;
t649 = qJD(1) ^ 2;
t660 = -t648 * g(1) - t640 * g(2);
t609 = -t649 * pkin(12) + t660;
t590 = -t647 * g(3) - t639 * t609;
t583 = (t639 * t647 * t649 + qJDD(2)) * pkin(1) + t590;
t592 = -t639 * g(3) + t647 * t609;
t672 = t647 ^ 2;
t584 = (-qJD(2) ^ 2 - t672 * t649) * pkin(1) + t592;
t545 = -t646 * t583 + t638 * t584;
t629 = qJDD(2) + qJDD(3);
t528 = (t602 * t604 + t629) * pkin(4) + t545;
t547 = -t638 * t583 - t646 * t584;
t535 = (-t602 ^ 2 - t631 ^ 2) * pkin(4) + t547;
t500 = t637 * t528 + t645 * t535;
t541 = -t574 * pkin(8) - t575 * pkin(10);
t621 = t625 ^ 2;
t623 = qJDD(4) + t629;
t491 = -t621 * pkin(8) + t623 * pkin(10) + t574 * t541 + t500;
t636 = sin(qJ(5));
t644 = cos(qJ(5));
t480 = t644 * t489 - t636 * t491;
t554 = -t636 * t575 + t644 * t625;
t497 = t554 * qJD(5) + t644 * t518 + t636 * t623;
t515 = qJDD(5) - t517;
t555 = t644 * t575 + t636 * t625;
t529 = -t554 * mrSges(6,1) + t555 * mrSges(6,2);
t573 = qJD(5) - t574;
t532 = -t573 * mrSges(6,2) + t554 * mrSges(6,3);
t475 = m(6) * t480 + t515 * mrSges(6,1) - t497 * mrSges(6,3) - t555 * t529 + t573 * t532;
t481 = t636 * t489 + t644 * t491;
t496 = -t555 * qJD(5) - t636 * t518 + t644 * t623;
t533 = t573 * mrSges(6,1) - t555 * mrSges(6,3);
t476 = m(6) * t481 - t515 * mrSges(6,2) + t496 * mrSges(6,3) + t554 * t529 - t573 * t533;
t463 = t644 * t475 + t636 * t476;
t564 = -t625 * mrSges(5,2) + t574 * mrSges(5,3);
t565 = t625 * mrSges(5,1) - t575 * mrSges(5,3);
t657 = -m(5) * t531 + t517 * mrSges(5,1) - t518 * mrSges(5,2) + t574 * t564 - t575 * t565 - t463;
t632 = sin(pkin(15));
t633 = cos(pkin(15));
t641 = cos(qJ(8));
t671 = sin(qJ(8));
t605 = t632 * t641 - t633 * t671;
t606 = -t632 * t671 - t633 * t641;
t551 = t605 * t601 + t606 * t603;
t506 = -t551 * qJD(8) + t606 * t560 - t605 * t562;
t550 = t606 * t601 - t605 * t603;
t507 = t550 * qJD(8) + t605 * t560 + t606 * t562;
t508 = (-t560 * t633 - t562 * t632 + (-t601 * t632 + t603 * t633) * t630) * pkin(3) + t580;
t624 = qJD(8) + t630;
t548 = -t624 * mrSges(9,2) + t550 * mrSges(9,3);
t549 = t624 * mrSges(9,1) - t551 * mrSges(9,3);
t659 = -m(9) * t508 + t506 * mrSges(9,1) - t507 * mrSges(9,2) + t550 * t548 - t551 * t549;
t673 = t561 * mrSges(4,1) + t560 * mrSges(8,1) - t563 * mrSges(4,2) - t562 * mrSges(8,2) + t601 * t585 + t602 * t586 - t603 * t587 - t604 * t588 + t657 + t659 + (-m(4) - m(8)) * t580;
t670 = pkin(3) * t632;
t669 = pkin(3) * t633;
t540 = -t574 * mrSges(5,1) + t575 * mrSges(5,2);
t661 = -t636 * t475 + t644 * t476;
t461 = m(5) * t500 - t623 * mrSges(5,2) + t517 * mrSges(5,3) + t574 * t540 - t625 * t565 + t661;
t499 = t645 * t528 - t637 * t535;
t490 = -t623 * pkin(8) - t621 * pkin(10) + t575 * t541 - t499;
t658 = -m(6) * t490 + t496 * mrSges(6,1) - t497 * mrSges(6,2) + t554 * t532 - t555 * t533;
t469 = m(5) * t499 + t623 * mrSges(5,1) - t518 * mrSges(5,3) - t575 * t540 + t625 * t564 + t658;
t668 = t637 * t461 + t645 * t469;
t544 = t642 * t583 - t634 * t584;
t572 = (-t601 * t633 - t603 * t632) * pkin(3);
t627 = t630 ^ 2;
t628 = qJDD(2) + qJDD(7);
t519 = -t603 * t572 + (t627 * t632 + t628 * t633) * pkin(3) + t544;
t546 = t634 * t583 + t642 * t584;
t520 = t601 * t572 + (-t627 * t633 + t628 * t632) * pkin(3) + t546;
t493 = t606 * t519 - t605 * t520;
t526 = -t550 * mrSges(9,1) + t551 * mrSges(9,2);
t622 = qJDD(8) + t628;
t487 = m(9) * t493 + t622 * mrSges(9,1) - t507 * mrSges(9,3) - t551 * t526 + t624 * t548;
t494 = t605 * t519 + t606 * t520;
t488 = m(9) * t494 - t622 * mrSges(9,2) + t506 * mrSges(9,3) + t550 * t526 - t624 * t549;
t667 = t606 * t487 + t605 * t488;
t666 = -t605 * t487 + t606 * t488;
t664 = qJD(1) * qJD(6);
t522 = Ifges(9,4) * t551 + Ifges(9,2) * t550 + Ifges(9,6) * t624;
t523 = Ifges(9,1) * t551 + Ifges(9,4) * t550 + Ifges(9,5) * t624;
t656 = mrSges(9,1) * t493 - mrSges(9,2) * t494 + Ifges(9,5) * t507 + Ifges(9,6) * t506 + Ifges(9,3) * t622 + t551 * t522 - t550 * t523;
t509 = Ifges(6,5) * t555 + Ifges(6,6) * t554 + Ifges(6,3) * t573;
t511 = Ifges(6,1) * t555 + Ifges(6,4) * t554 + Ifges(6,5) * t573;
t466 = -mrSges(6,1) * t490 + mrSges(6,3) * t481 + Ifges(6,4) * t497 + Ifges(6,2) * t496 + Ifges(6,6) * t515 - t555 * t509 + t573 * t511;
t510 = Ifges(6,4) * t555 + Ifges(6,2) * t554 + Ifges(6,6) * t573;
t467 = mrSges(6,2) * t490 - mrSges(6,3) * t480 + Ifges(6,1) * t497 + Ifges(6,4) * t496 + Ifges(6,5) * t515 + t554 * t509 - t573 * t510;
t537 = Ifges(5,4) * t575 + Ifges(5,2) * t574 + Ifges(5,6) * t625;
t538 = Ifges(5,1) * t575 + Ifges(5,4) * t574 + Ifges(5,5) * t625;
t655 = pkin(8) * t658 + pkin(10) * t661 + mrSges(5,1) * t499 - mrSges(5,2) * t500 + Ifges(5,5) * t518 + Ifges(5,6) * t517 + Ifges(5,3) * t623 + t644 * t466 + t636 * t467 + t575 * t537 - t574 * t538;
t654 = mrSges(6,1) * t480 - mrSges(6,2) * t481 + Ifges(6,5) * t497 + Ifges(6,6) * t496 + Ifges(6,3) * t515 + t555 * t510 - t554 * t511;
t568 = Ifges(8,4) * t603 + Ifges(8,2) * t601 + Ifges(8,6) * t630;
t570 = Ifges(8,1) * t603 + Ifges(8,4) * t601 + Ifges(8,5) * t630;
t652 = mrSges(8,1) * t544 - mrSges(8,2) * t546 + Ifges(8,5) * t562 + Ifges(8,6) * t560 + Ifges(8,3) * t628 + t603 * t568 - t601 * t570 + t666 * t670 + t667 * t669 + t656;
t569 = Ifges(4,4) * t604 + Ifges(4,2) * t602 + Ifges(4,6) * t631;
t571 = Ifges(4,1) * t604 + Ifges(4,4) * t602 + Ifges(4,5) * t631;
t650 = pkin(4) * t668 + mrSges(4,1) * t545 - mrSges(4,2) * t547 + Ifges(4,5) * t563 + Ifges(4,6) * t561 + Ifges(4,3) * t629 + t604 * t569 - t602 * t571 + t655;
t643 = cos(qJ(6));
t635 = sin(qJ(6));
t613 = t643 * qJDD(1) - t635 * t664;
t611 = t635 * qJDD(1) + t643 * t664;
t610 = t649 * pkin(6) + t660;
t608 = qJDD(1) * pkin(6) - t662;
t600 = Ifges(3,5) * qJD(2) + (t639 * Ifges(3,1) + t647 * Ifges(3,4)) * qJD(1);
t599 = Ifges(7,5) * qJD(6) + (t635 * Ifges(7,1) + t643 * Ifges(7,4)) * qJD(1);
t598 = Ifges(3,6) * qJD(2) + (t639 * Ifges(3,4) + t647 * Ifges(3,2)) * qJD(1);
t597 = Ifges(7,6) * qJD(6) + (t635 * Ifges(7,4) + t643 * Ifges(7,2)) * qJD(1);
t591 = -t635 * g(3) + t643 * t610;
t589 = -t643 * g(3) - t635 * t610;
t579 = -t602 * mrSges(4,1) + t604 * mrSges(4,2);
t578 = -t601 * mrSges(8,1) + t603 * mrSges(8,2);
t567 = Ifges(4,5) * t604 + Ifges(4,6) * t602 + Ifges(4,3) * t631;
t566 = Ifges(8,5) * t603 + Ifges(8,6) * t601 + Ifges(8,3) * t630;
t536 = Ifges(5,5) * t575 + Ifges(5,6) * t574 + Ifges(5,3) * t625;
t521 = Ifges(9,5) * t551 + Ifges(9,6) * t550 + Ifges(9,3) * t624;
t478 = mrSges(9,2) * t508 - mrSges(9,3) * t493 + Ifges(9,1) * t507 + Ifges(9,4) * t506 + Ifges(9,5) * t622 + t550 * t521 - t624 * t522;
t477 = -mrSges(9,1) * t508 + mrSges(9,3) * t494 + Ifges(9,4) * t507 + Ifges(9,2) * t506 + Ifges(9,6) * t622 - t551 * t521 + t624 * t523;
t459 = mrSges(8,2) * t580 - mrSges(8,3) * t544 + Ifges(8,1) * t562 + Ifges(8,4) * t560 + Ifges(8,5) * t628 - t605 * t477 + t606 * t478 + t601 * t566 - t630 * t568 + t659 * t670;
t458 = -mrSges(8,1) * t580 + mrSges(8,3) * t546 + Ifges(8,4) * t562 + Ifges(8,2) * t560 + Ifges(8,6) * t628 + t606 * t477 + t605 * t478 - t603 * t566 + t630 * t570 + t659 * t669;
t456 = -pkin(8) * t463 - mrSges(5,1) * t531 + mrSges(5,3) * t500 + Ifges(5,4) * t518 + Ifges(5,2) * t517 + Ifges(5,6) * t623 - t575 * t536 + t625 * t538 - t654;
t455 = -pkin(10) * t463 + mrSges(5,2) * t531 - mrSges(5,3) * t499 + Ifges(5,1) * t518 + Ifges(5,4) * t517 + Ifges(5,5) * t623 - t636 * t466 + t644 * t467 + t574 * t536 - t625 * t537;
t454 = mrSges(4,2) * t580 - mrSges(4,3) * t545 + Ifges(4,1) * t563 + Ifges(4,4) * t561 + Ifges(4,5) * t629 + t645 * t455 - t637 * t456 + t602 * t567 - t631 * t569;
t453 = pkin(4) * t657 - mrSges(4,1) * t580 + mrSges(4,3) * t547 + Ifges(4,4) * t563 + Ifges(4,2) * t561 + Ifges(4,6) * t629 + t637 * t455 + t645 * t456 - t604 * t567 + t631 * t571;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t662 - mrSges(2,2) * t660 + t639 * (mrSges(3,2) * t607 - mrSges(3,3) * t590 + Ifges(3,1) * t612 + Ifges(3,4) * t614 + Ifges(3,5) * qJDD(2) - qJD(2) * t598 + t638 * t453 - t646 * t454 - t634 * t458 + t642 * t459) + t647 * (t673 * pkin(1) - mrSges(3,1) * t607 + mrSges(3,3) * t592 + Ifges(3,4) * t612 + Ifges(3,2) * t614 + Ifges(3,6) * qJDD(2) + qJD(2) * t600 - t646 * t453 - t638 * t454 + t642 * t458 + t634 * t459) + pkin(12) * (-m(3) * t607 + t614 * mrSges(3,1) - t612 * mrSges(3,2) + t673) + t635 * (mrSges(7,2) * t608 - mrSges(7,3) * t589 + Ifges(7,1) * t611 + Ifges(7,4) * t613 + Ifges(7,5) * qJDD(6) - qJD(6) * t597) + t643 * (-mrSges(7,1) * t608 + mrSges(7,3) * t591 + Ifges(7,4) * t611 + Ifges(7,2) * t613 + Ifges(7,6) * qJDD(6) + qJD(6) * t599) - pkin(6) * (-m(7) * t608 + t613 * mrSges(7,1) - t611 * mrSges(7,2)) + ((-pkin(6) * (t635 ^ 2 + t643 ^ 2) * mrSges(7,3) + pkin(12) * (t639 ^ 2 + t672) * mrSges(3,3)) * qJD(1) - pkin(6) * (-mrSges(7,1) * t635 - mrSges(7,2) * t643) * qJD(6) + pkin(12) * (-mrSges(3,1) * t639 - mrSges(3,2) * t647) * qJD(2)) * qJD(1); (-t638 * (m(4) * t547 - t629 * mrSges(4,2) + t561 * mrSges(4,3) + t645 * t461 - t637 * t469 + t602 * t579 - t631 * t588) - t646 * (m(4) * t545 + t629 * mrSges(4,1) - t563 * mrSges(4,3) - t604 * t579 + t631 * t586 + t668) + t634 * (m(8) * t546 - t628 * mrSges(8,2) + t560 * mrSges(8,3) + t601 * t578 - t630 * t587 + t666) + t642 * (m(8) * t544 + t628 * mrSges(8,1) - t562 * mrSges(8,3) - t603 * t578 + t630 * t585 + t667)) * pkin(1) + t652 + t650 + Ifges(3,3) * qJDD(2) + (t639 * t598 - t647 * t600) * qJD(1) + Ifges(3,5) * t612 + Ifges(3,6) * t614 + mrSges(3,1) * t590 - mrSges(3,2) * t592; t650; t655; t654; mrSges(7,1) * t589 - mrSges(7,2) * t591 + Ifges(7,5) * t611 + Ifges(7,6) * t613 + Ifges(7,3) * qJDD(6) + (t635 * t597 - t643 * t599) * qJD(1); t652; t656; 0; 0;];
tauJ = t1;
