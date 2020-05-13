% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+12)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = picker2Dm1OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_invdynJB_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:09
% EndTime: 2020-05-11 05:45:11
% DurationCPUTime: 1.74s
% Computational Cost: add. (14002->282), mult. (17680->341), div. (0->0), fcn. (9718->22), ass. (0->141)
t687 = sin(qJ(1));
t696 = cos(qJ(1));
t638 = -t687 * g(1) + t696 * g(2);
t633 = qJDD(1) * pkin(1) + t638;
t640 = t696 * g(1) + t687 * g(2);
t699 = qJD(1) ^ 2;
t634 = -t699 * pkin(1) + t640;
t686 = sin(qJ(2));
t695 = cos(qJ(2));
t617 = t695 * t633 - t686 * t634;
t619 = t686 * t633 + t695 * t634;
t682 = sin(qJ(6));
t691 = cos(qJ(6));
t592 = -t691 * t617 + t682 * t619;
t671 = qJD(1) + qJD(2);
t660 = qJD(6) + t671;
t654 = t660 ^ 2;
t669 = qJDD(1) + qJDD(2);
t657 = qJDD(6) + t669;
t583 = m(7) * t592 + t657 * mrSges(7,1) - t654 * mrSges(7,2);
t593 = -t682 * t617 - t691 * t619;
t584 = m(7) * t593 - t654 * mrSges(7,1) - t657 * mrSges(7,2);
t667 = t671 ^ 2;
t611 = t669 * pkin(2) + t617;
t613 = -t667 * pkin(2) + t619;
t685 = sin(qJ(3));
t694 = cos(qJ(3));
t588 = -t694 * t611 + t685 * t613;
t662 = qJD(3) + t671;
t656 = t662 ^ 2;
t659 = qJDD(3) + t669;
t580 = t659 * pkin(6) + t588;
t590 = -t685 * t611 - t694 * t613;
t582 = -t656 * pkin(6) + t590;
t679 = sin(qJ(9));
t688 = cos(qJ(9));
t574 = -t688 * t580 + t679 * t582;
t648 = qJDD(9) + t659;
t651 = qJD(9) + t662;
t649 = t651 ^ 2;
t568 = m(10) * t574 + t648 * mrSges(10,1) - t649 * mrSges(10,2);
t575 = -t679 * t580 - t688 * t582;
t569 = m(10) * t575 - t649 * mrSges(10,1) - t648 * mrSges(10,2);
t705 = -t688 * t568 - t679 * t569;
t558 = m(4) * t588 + t659 * mrSges(4,1) - t656 * mrSges(4,2) + t705;
t559 = m(4) * t590 - t656 * mrSges(4,1) - t659 * mrSges(4,2) + t679 * t568 - t688 * t569;
t707 = -t694 * t558 - t685 * t559;
t610 = t669 * pkin(3) + t617;
t612 = -t667 * pkin(3) + t619;
t684 = sin(qJ(4));
t693 = cos(qJ(4));
t587 = t693 * t610 - t684 * t612;
t661 = qJD(4) + t671;
t655 = t661 ^ 2;
t658 = qJDD(4) + t669;
t579 = t658 * pkin(4) + t587;
t589 = t684 * t610 + t693 * t612;
t581 = -t655 * pkin(4) + t589;
t675 = sin(qJ(10));
t677 = cos(qJ(10));
t572 = -t677 * t579 + t675 * t581;
t643 = qJDD(10) + t658;
t650 = qJD(10) + t661;
t647 = t650 ^ 2;
t566 = m(11) * t572 + t643 * mrSges(11,1) - t647 * mrSges(11,2);
t573 = -t675 * t579 - t677 * t581;
t567 = m(11) * t573 - t647 * mrSges(11,1) - t643 * mrSges(11,2);
t706 = -t677 * t566 - t675 * t567;
t556 = m(5) * t587 + t658 * mrSges(5,1) - t655 * mrSges(5,2) + t706;
t557 = m(5) * t589 - t655 * mrSges(5,1) - t658 * mrSges(5,2) + t675 * t566 - t677 * t567;
t718 = t693 * t556 + t684 * t557;
t546 = m(3) * t617 + t669 * mrSges(3,1) - t667 * mrSges(3,2) - t691 * t583 - t682 * t584 + t707 + t718;
t547 = m(3) * t619 - t667 * mrSges(3,1) - t669 * mrSges(3,2) - t684 * t556 + t693 * t557 + t685 * t558 - t694 * t559 + t682 * t583 - t691 * t584;
t680 = sin(qJ(8));
t689 = cos(qJ(8));
t616 = -t689 * t633 + t680 * t634;
t670 = qJD(1) + qJD(8);
t666 = t670 ^ 2;
t668 = qJDD(1) + qJDD(8);
t600 = m(9) * t616 + t668 * mrSges(9,1) - t666 * mrSges(9,2);
t618 = -t680 * t633 - t689 * t634;
t601 = m(9) * t618 - t666 * mrSges(9,1) - t668 * mrSges(9,2);
t723 = t695 * t546 + t686 * t547 - t689 * t600 - t680 * t601;
t722 = pkin(5) * m(6);
t721 = -m(4) - m(10);
t720 = -m(11) - m(5);
t676 = sin(pkin(8));
t678 = cos(pkin(8));
t683 = sin(qJ(5));
t692 = cos(qJ(5));
t631 = t676 * t692 + t678 * t683;
t632 = -t676 * t683 + t678 * t692;
t621 = t631 * g(1) - t632 * g(2);
t698 = qJD(5) ^ 2;
t608 = m(6) * t621 + qJDD(5) * mrSges(6,1) - t698 * mrSges(6,2);
t622 = -t632 * g(1) - t631 * g(2);
t609 = m(6) * t622 - t698 * mrSges(6,1) - qJDD(5) * mrSges(6,2);
t717 = t632 * t608 + t631 * t609;
t716 = -t631 * t608 + t632 * t609;
t681 = sin(qJ(7));
t690 = cos(qJ(7));
t637 = -t681 * g(1) + t690 * g(2);
t697 = qJD(7) ^ 2;
t626 = m(8) * t637 - t697 * mrSges(8,1) - qJDD(7) * mrSges(8,2);
t639 = -t690 * g(1) - t681 * g(2);
t627 = m(8) * t639 + qJDD(7) * mrSges(8,1) - t697 * mrSges(8,2);
t715 = -t690 * t626 + t681 * t627;
t714 = mrSges(6,1) * t621 - mrSges(6,2) * t622 + Ifges(6,3) * qJDD(5);
t713 = mrSges(7,1) * t592 - mrSges(7,2) * t593 + Ifges(7,3) * t657;
t712 = mrSges(8,1) * t639 - mrSges(8,2) * t637 + Ifges(8,3) * qJDD(7);
t711 = mrSges(9,1) * t616 - mrSges(9,2) * t618 + Ifges(9,3) * t668;
t710 = mrSges(10,1) * t574 - mrSges(10,2) * t575 + Ifges(10,3) * t648;
t709 = mrSges(11,1) * t572 - mrSges(11,2) * t573 + Ifges(11,3) * t643;
t708 = m(3) + m(7) + m(9) - t720 - t721;
t703 = pkin(6) * t705 + mrSges(4,1) * t588 - mrSges(4,2) * t590 + Ifges(4,3) * t659 + t710;
t702 = pkin(4) * t706 + mrSges(5,1) * t587 - mrSges(5,2) * t589 + Ifges(5,3) * t658 + t709;
t701 = pkin(2) * t707 + pkin(3) * t718 + mrSges(3,1) * t617 - mrSges(3,2) * t619 + Ifges(3,3) * t669 + t702 + t703 + t713;
t700 = t723 * pkin(1) + mrSges(2,1) * t638 - mrSges(2,2) * t640 + Ifges(2,3) * qJDD(1) + t701 + t711;
t625 = -mrSges(8,2) * g(3) - mrSges(8,3) * t639 + Ifges(8,5) * qJDD(7) - t697 * Ifges(8,6);
t624 = mrSges(8,1) * g(3) + mrSges(8,3) * t637 + t697 * Ifges(8,5) + Ifges(8,6) * qJDD(7);
t603 = -mrSges(6,2) * g(3) - mrSges(6,3) * t621 + Ifges(6,5) * qJDD(5) - t698 * Ifges(6,6);
t602 = mrSges(6,1) * g(3) + mrSges(6,3) * t622 + t698 * Ifges(6,5) + Ifges(6,6) * qJDD(5);
t599 = -mrSges(9,2) * g(3) - mrSges(9,3) * t616 + Ifges(9,5) * t668 - t666 * Ifges(9,6);
t598 = mrSges(9,1) * g(3) + mrSges(9,3) * t618 + t666 * Ifges(9,5) + Ifges(9,6) * t668;
t578 = -mrSges(7,2) * g(3) - mrSges(7,3) * t592 + Ifges(7,5) * t657 - t654 * Ifges(7,6);
t577 = mrSges(7,1) * g(3) + mrSges(7,3) * t593 + t654 * Ifges(7,5) + Ifges(7,6) * t657;
t565 = -mrSges(10,2) * g(3) - mrSges(10,3) * t574 + Ifges(10,5) * t648 - t649 * Ifges(10,6);
t564 = mrSges(10,1) * g(3) + mrSges(10,3) * t575 + t649 * Ifges(10,5) + Ifges(10,6) * t648;
t563 = -mrSges(11,2) * g(3) - mrSges(11,3) * t572 + Ifges(11,5) * t643 - t647 * Ifges(11,6);
t562 = mrSges(11,1) * g(3) + mrSges(11,3) * t573 + t647 * Ifges(11,5) + Ifges(11,6) * t643;
t553 = -mrSges(4,2) * g(3) - mrSges(4,3) * t588 + Ifges(4,5) * t659 - t656 * Ifges(4,6) + t679 * t564 - t688 * t565;
t552 = -mrSges(5,2) * g(3) - mrSges(5,3) * t587 + Ifges(5,5) * t658 - t655 * Ifges(5,6) + t675 * t562 - t677 * t563;
t551 = mrSges(4,3) * t590 + t656 * Ifges(4,5) + Ifges(4,6) * t659 - t688 * t564 - t679 * t565 + (pkin(6) * m(10) + mrSges(4,1)) * g(3);
t550 = mrSges(5,3) * t589 + t655 * Ifges(5,5) + Ifges(5,6) * t658 - t677 * t562 - t675 * t563 + (pkin(4) * m(11) + mrSges(5,1)) * g(3);
t543 = -mrSges(3,2) * g(3) - mrSges(3,3) * t617 + Ifges(3,5) * t669 - t667 * Ifges(3,6) - t684 * t550 + t685 * t551 + t693 * t552 - t694 * t553 + t682 * t577 - t691 * t578;
t542 = mrSges(3,3) * t619 + t667 * Ifges(3,5) + Ifges(3,6) * t669 + t693 * t550 - t694 * t551 + t684 * t552 - t685 * t553 - t691 * t577 - t682 * t578 + (-pkin(2) * t721 - pkin(3) * t720 + mrSges(3,1)) * g(3);
t540 = m(2) * t640 - t699 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t686 * t546 + t695 * t547 + t680 * t600 - t689 * t601;
t539 = m(2) * t638 + qJDD(1) * mrSges(2,1) - t699 * mrSges(2,2) + t723;
t538 = -mrSges(2,2) * g(3) - mrSges(2,3) * t638 + Ifges(2,5) * qJDD(1) - t699 * Ifges(2,6) - t686 * t542 + t695 * t543 + t680 * t598 - t689 * t599;
t537 = mrSges(2,3) * t640 + t699 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t695 * t542 + t686 * t543 - t689 * t598 - t680 * t599 + (t708 * pkin(1) + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t687 * t539 - t696 * t540 + t681 * t626 + t690 * t627 + t716; -m(1) * g(2) - t696 * t539 - t687 * t540 + t715 + t717; (-m(1) - m(2) - m(6) - m(8) - t708) * g(3); mrSges(1,3) * g(2) + t687 * t537 - t696 * t538 - t631 * t602 + t632 * t603 + t690 * t624 + t681 * t625 + (-t676 * t722 - mrSges(1,2)) * g(3); -mrSges(1,3) * g(1) - t696 * t537 - t687 * t538 + t632 * t602 + t631 * t603 + t681 * t624 - t690 * t625 + (m(8) * pkin(7) + t678 * t722 + mrSges(1,1)) * g(3); -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t700 + (-t676 * t716 + t678 * t717) * pkin(5) + pkin(7) * t715 + t712 + t714; t700; t701; t703; t702; t714; t713; t712; t711; t710; t709; 0; 0;];
tauJB = t1;
