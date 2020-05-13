% Calculate joint inertia matrix with Newton Euler for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_inertiaJ_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_inertiaJ_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_inertiaJ_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_inertiaJ_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1IC_inertiaJ_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:03:11
% EndTime: 2020-04-15 20:03:16
% DurationCPUTime: 4.14s
% Computational Cost: add. (6351->488), mult. (9241->616), div. (156->8), fcn. (9526->32), ass. (0->255)
t708 = -qJ(7) + pkin(19);
t694 = -qJ(10) + t708;
t677 = sin(t694);
t636 = t677 * pkin(8) - pkin(4) * sin(t708);
t678 = cos(t694);
t637 = t678 * pkin(8) - pkin(4) * cos(t708);
t696 = pkin(20) + qJ(7) + qJ(2);
t680 = sin(t696);
t682 = cos(t696);
t716 = sin(qJ(6));
t724 = cos(qJ(6));
t611 = 0.1e1 / (-t680 * t716 - t682 * t724) / pkin(7) / pkin(3);
t720 = sin(qJ(2));
t707 = t720 * pkin(1);
t660 = pkin(3) * t680 + t707;
t728 = cos(qJ(2));
t661 = t728 * pkin(1) + pkin(3) * t682;
t521 = (t660 * t716 + t661 * t724) * pkin(7) * t611;
t695 = qJ(3) + qJ(4) + pkin(18);
t679 = sin(t695);
t681 = cos(t695);
t595 = 0.1e1 / (-t677 * t679 + t678 * t681) / pkin(10) / pkin(8);
t735 = t595 * t521;
t469 = (t636 * t678 - t637 * t677) * pkin(8) * t735;
t717 = sin(qJ(5));
t725 = cos(qJ(5));
t767 = Ifges(6,5) * t717 + Ifges(6,6) * t725;
t748 = t717 * mrSges(6,1) + t725 * mrSges(6,2);
t808 = t748 * pkin(11);
t810 = t767 - t808;
t812 = t469 * t810;
t719 = sin(qJ(3));
t658 = t719 * pkin(5) + pkin(10) * t679;
t727 = cos(qJ(3));
t659 = t727 * pkin(5) + pkin(10) * t681;
t510 = (t658 * t677 - t659 * t678) * t595 * pkin(8);
t811 = t510 * t810;
t710 = qJ(8) + qJ(9);
t697 = sin(t710);
t714 = sin(qJ(8));
t662 = t714 * pkin(2) - t697 * pkin(12);
t698 = cos(t710);
t722 = cos(qJ(8));
t663 = -t722 * pkin(2) + t698 * pkin(12);
t709 = qJ(3) + pkin(17);
t692 = sin(t709);
t693 = cos(t709);
t797 = pkin(6) / (t662 * t698 + t663 * t697) / pkin(12);
t509 = (-t662 * t692 + t663 * t693) * t797;
t514 = (-t692 * t697 - t693 * t698) * pkin(12) * t797;
t713 = sin(qJ(9));
t721 = cos(qJ(9));
t798 = pkin(2) * t721;
t799 = pkin(2) * t713;
t809 = pkin(2) * (-t721 * mrSges(10,1) + t713 * mrSges(10,2)) + Ifges(10,3);
t738 = Ifges(9,3) + pkin(2) * (-t713 * (-m(10) * t799 - mrSges(10,2)) - t721 * (-m(10) * t798 + mrSges(10,1))) + t809;
t776 = t809 * t509 + t514 * t738;
t792 = t719 * pkin(1);
t685 = pkin(5) + t792;
t718 = sin(qJ(4));
t726 = cos(qJ(4));
t778 = t726 * t727;
t631 = -pkin(1) * t778 + t718 * t685;
t627 = pkin(11) + t631;
t807 = t748 * t627;
t793 = t718 * pkin(5);
t684 = pkin(11) + t793;
t806 = t748 * t684;
t654 = t719 * t728 + t727 * t720;
t657 = -t719 * t720 + t727 * t728;
t588 = t718 * t654 - t726 * t657;
t590 = t726 * t654 + t718 * t657;
t683 = -pkin(15) + t707;
t621 = -t657 * pkin(5) + t683;
t515 = t588 * pkin(9) - t590 * pkin(11) + t621;
t749 = t725 * mrSges(6,1) - t717 * mrSges(6,2);
t805 = t749 * t515;
t715 = sin(qJ(7));
t723 = cos(qJ(7));
t653 = -t715 * t728 - t723 * t720;
t656 = -t715 * t720 + t723 * t728;
t750 = m(6) * t515 - mrSges(6,3) * t590;
t482 = -t588 * mrSges(6,2) + t750 * t717;
t483 = t588 * mrSges(6,1) + t750 * t725;
t465 = t717 * t482 + t725 * t483;
t739 = m(5) * t621 + t588 * mrSges(5,1) + t590 * mrSges(5,2) + t465;
t711 = sin(pkin(19));
t712 = cos(qJ(10));
t784 = cos(pkin(19));
t789 = sin(qJ(10));
t644 = t711 * t712 - t784 * t789;
t645 = -t711 * t789 - t784 * t712;
t544 = -t644 * t656 + t645 * t653;
t545 = t644 * t653 + t645 * t656;
t573 = (-t784 * t653 - t656 * t711) * pkin(4) + t683;
t751 = -m(11) * t573 + t544 * mrSges(11,1) - t545 * mrSges(11,2);
t804 = -t657 * mrSges(4,1) - t653 * mrSges(8,1) + t654 * mrSges(4,2) + t656 * mrSges(8,2) + (m(4) + m(8)) * t683 + t739 - t751;
t803 = -0.2e1 * t728;
t802 = pkin(11) * m(6);
t800 = pkin(1) * t727;
t795 = t711 * pkin(4);
t794 = t715 * pkin(1);
t791 = t723 * pkin(1);
t790 = t726 * pkin(5);
t788 = mrSges(6,3) * t717;
t704 = Ifges(6,4) * t717;
t705 = Ifges(6,4) * t725;
t699 = t725 * mrSges(6,3);
t783 = mrSges(11,3) * t644;
t782 = mrSges(11,3) * t645;
t470 = (-t636 * t679 + t637 * t681) * pkin(10) * t735;
t781 = t470 * Ifges(11,3);
t780 = t627 * t725;
t779 = t684 * t725;
t761 = t725 * t482 - t717 * t483;
t462 = -t588 * mrSges(5,3) + t761;
t507 = (-mrSges(5,3) - t748) * t590;
t777 = t718 * t462 + t726 * t507;
t605 = (-m(6) * t627 - mrSges(6,3)) * t717;
t606 = m(6) * t780 + t699;
t759 = -t717 * t605 + t725 * t606;
t522 = m(5) * t631 - mrSges(5,2) + t759;
t630 = t726 * t685 + t718 * t800;
t626 = -pkin(9) - t630;
t743 = -m(6) * t626 + t749;
t562 = m(5) * t630 + mrSges(5,1) + t743;
t775 = t718 * t522 + t726 * t562;
t668 = t794 + t795;
t762 = t784 * pkin(4);
t671 = t762 + t791;
t560 = -t644 * t668 + t645 * t671;
t550 = m(11) * t560 + mrSges(11,1);
t561 = t644 * t671 + t645 * t668;
t551 = m(11) * t561 - mrSges(11,2);
t774 = t645 * t550 + t644 * t551;
t773 = -t644 * t550 + t645 * t551;
t503 = Ifges(11,5) * t545 + Ifges(11,6) * t544;
t577 = (-t644 * t711 + t784 * t645) * pkin(4);
t570 = m(11) * t577 + mrSges(11,1);
t578 = (t784 * t644 + t645 * t711) * pkin(4);
t571 = m(11) * t578 - mrSges(11,2);
t772 = t645 * t570 + t644 * t571;
t771 = -t644 * t570 + t645 * t571;
t634 = (-m(6) * t684 - mrSges(6,3)) * t717;
t635 = m(6) * t779 + t699;
t757 = -t717 * t634 + t725 * t635;
t566 = m(5) * t793 - mrSges(5,2) + t757;
t686 = -pkin(9) - t790;
t742 = -m(6) * t686 + t749;
t612 = m(5) * t790 + mrSges(5,1) + t742;
t770 = t718 * t566 + t726 * t612;
t652 = -t714 * t728 - t722 * t720;
t655 = -t714 * t720 + t722 * t728;
t587 = -t721 * t652 + t713 * t655;
t589 = -t713 * t652 - t721 * t655;
t769 = Ifges(10,5) * t589 + Ifges(10,6) * t587;
t674 = -mrSges(10,3) * t799 + Ifges(10,6);
t675 = mrSges(10,3) * t798 + Ifges(10,5);
t614 = -t721 * t674 - t713 * t675 + Ifges(9,6);
t615 = t713 * t674 - t721 * t675 + Ifges(9,5);
t768 = t722 * t614 + t714 * t615;
t766 = Ifges(6,2) * t725 + t704;
t765 = Ifges(6,1) * t717 + t705;
t764 = Ifges(5,6) - t767;
t592 = t645 * mrSges(11,1) - t644 * mrSges(11,2);
t594 = -t644 * mrSges(11,1) - t645 * mrSges(11,2);
t763 = t592 * t762 + t594 * t795 + Ifges(11,3);
t511 = mrSges(11,1) * t560 - mrSges(11,2) * t561 + Ifges(11,3);
t523 = mrSges(11,1) * t577 - mrSges(11,2) * t578 + Ifges(11,3);
t760 = t523 * t521 + t511;
t758 = -t714 * t614 + t722 * t615;
t665 = (-mrSges(6,3) - t802) * t717;
t666 = t725 * t802 + t699;
t756 = -t717 * t665 + t725 * t666;
t755 = -t748 * t793 + t810;
t754 = m(6) * pkin(9) + t749;
t753 = t544 * t783 - t545 * t782;
t752 = t544 * t782 + t545 * t783;
t628 = mrSges(6,1) * pkin(9) + pkin(11) * t699 + t766;
t629 = -mrSges(6,2) * pkin(9) + pkin(11) * t788 + t765;
t519 = pkin(9) * t754 + pkin(11) * t756 + t725 * t628 + t717 * t629 + Ifges(5,3);
t747 = t725 * Ifges(6,5) - t717 * Ifges(6,6);
t746 = Ifges(9,5) * t655 + Ifges(9,6) * t652 + pkin(2) * (-t587 * t713 + t589 * t721) * mrSges(10,3) + t769;
t601 = -mrSges(5,2) + t756;
t649 = mrSges(5,1) + t754;
t559 = t718 * t601 + t726 * t649;
t741 = pkin(5) * t559 + t519;
t475 = t515 * t788 + Ifges(6,6) * t588 + (-Ifges(6,2) * t717 + t705) * t590;
t476 = -t515 * t699 + Ifges(6,5) * t588 + (Ifges(6,1) * t725 - t704) * t590;
t451 = pkin(11) * t761 - Ifges(5,6) * t588 + t725 * t475 + t717 * t476 + (-pkin(9) * t748 + Ifges(5,5)) * t590;
t740 = Ifges(8,5) * t656 + Ifges(8,6) * t653 + t752 * t795 + t753 * t762 + t503;
t632 = -t652 * pkin(2) - pkin(15);
t737 = m(10) * t632 - t587 * mrSges(10,1) + t589 * mrSges(10,2);
t554 = -mrSges(6,1) * t626 + mrSges(6,3) * t780 + t766;
t555 = mrSges(6,2) * t626 + t627 * t788 + t765;
t468 = pkin(9) * t743 + pkin(11) * t759 + mrSges(5,1) * t630 - mrSges(5,2) * t631 + t725 * t554 + t717 * t555 + Ifges(5,3);
t734 = mrSges(8,1) * t791 + t774 * t762 + t773 * t795 + Ifges(8,3) + t511;
t733 = pkin(5) * t777 + Ifges(4,5) * t654 + Ifges(4,6) * t657 + t451;
t603 = -mrSges(6,1) * t686 + mrSges(6,3) * t779 + t766;
t604 = mrSges(6,2) * t686 + t684 * t788 + t765;
t498 = pkin(9) * t742 + pkin(11) * t757 + mrSges(5,1) * t790 - mrSges(5,2) * t793 + t725 * t603 + t717 * t604 + Ifges(5,3);
t732 = pkin(5) * t770 + Ifges(4,3) + t498;
t731 = pkin(5) * t775 + mrSges(4,1) * t792 + mrSges(4,2) * t800 + Ifges(4,3) + t468;
t670 = -t721 * Ifges(10,5) + t713 * Ifges(10,6);
t667 = -t713 * Ifges(10,5) - t721 * Ifges(10,6);
t651 = -pkin(9) * t749 - Ifges(6,3);
t619 = -pkin(11) * t749 + t747;
t602 = t725 * t665 + t717 * t666;
t593 = t645 * Ifges(11,5) - t644 * Ifges(11,6);
t591 = t644 * Ifges(11,5) + t645 * Ifges(11,6);
t576 = t725 * t634 + t717 * t635;
t572 = t726 * t619 - t718 * t651;
t569 = -mrSges(11,3) * t577 + Ifges(11,5);
t568 = mrSges(11,3) * t578 + Ifges(11,6);
t557 = -pkin(9) * t602 + t764 + t808;
t552 = -pkin(5) * t749 + t718 * t619 + t726 * t651;
t549 = -mrSges(11,3) * t560 + Ifges(11,5);
t548 = mrSges(11,3) * t561 + Ifges(11,6);
t536 = t725 * t605 + t717 * t606;
t528 = -pkin(11) * t602 - t717 * t628 + t725 * t629 + Ifges(5,5);
t527 = -pkin(9) * t576 + mrSges(5,3) * t793 + t764 + t806;
t518 = mrSges(10,2) * t632 + Ifges(10,1) * t589 + Ifges(10,4) * t587;
t517 = -mrSges(10,1) * t632 + Ifges(10,4) * t589 + Ifges(10,2) * t587;
t516 = (t660 * t682 - t661 * t680) * t611 * pkin(3);
t513 = t516 * (t716 * Ifges(7,5) + t724 * Ifges(7,6));
t512 = (-t658 * t681 + t659 * t679) * t595 * pkin(10);
t508 = -pkin(11) * t576 - mrSges(5,3) * t790 - t717 * t603 + t725 * t604 + Ifges(5,5);
t505 = -t644 * t568 + t645 * t569 + Ifges(8,5);
t504 = t645 * t568 + t644 * t569 + Ifges(8,6);
t499 = t726 * t528 - t718 * t557;
t496 = -pkin(9) * t536 + mrSges(5,3) * t631 + t764 + t807;
t492 = -mrSges(8,3) * t791 - t644 * t548 + t645 * t549 + Ifges(8,5);
t491 = mrSges(8,3) * t794 + t645 * t548 + t644 * t549 + Ifges(8,6);
t488 = t728 * (-t715 * t591 + t723 * t593) - t720 * (t723 * t591 + t715 * t593);
t486 = mrSges(11,2) * t573 + Ifges(11,1) * t545 + Ifges(11,4) * t544;
t485 = -mrSges(11,1) * t573 + Ifges(11,4) * t545 + Ifges(11,2) * t544;
t484 = -pkin(5) * t602 + t718 * t528 + t726 * t557;
t478 = t726 * t508 - t718 * t527 + Ifges(4,5);
t477 = -pkin(11) * t536 - mrSges(5,3) * t630 - t717 * t554 + t725 * t555 + Ifges(5,5);
t472 = -mrSges(9,2) * pkin(15) + Ifges(9,1) * t655 + Ifges(9,4) * t652 + t713 * t517 - t721 * t518;
t471 = -pkin(5) * t576 + t718 * t508 + t726 * t527 + Ifges(4,6);
t467 = -pkin(2) * t737 + mrSges(9,1) * pkin(15) + Ifges(9,4) * t655 + Ifges(9,2) * t652 - t721 * t517 - t713 * t518;
t466 = t510 * t519 + t741;
t463 = -mrSges(4,3) * t792 + t726 * t477 - t718 * t496 + Ifges(4,5);
t460 = -pkin(5) * t536 - mrSges(4,3) * t800 + t718 * t477 + t726 * t496 + Ifges(4,6);
t458 = mrSges(8,2) * t683 + Ifges(8,1) * t656 + Ifges(8,4) * t653 - t644 * t485 + t645 * t486 + t751 * t795;
t457 = -mrSges(8,1) * t683 + Ifges(8,4) * t656 + Ifges(8,2) * t653 + t645 * t485 + t644 * t486 + t751 * t762;
t456 = pkin(1) * (t723 * t592 + t715 * t594) + t521 * t763 - t781 + t763;
t455 = pkin(1) * (-t727 * (t726 * t601 - t718 * t649) + t719 * t559) - t469 * t519 + t741;
t454 = t728 * (t727 * t484 + t719 * t499) - t720 * (-pkin(1) * t602 + t719 * t484 - t727 * t499) - pkin(15) * t602;
t453 = -pkin(9) * t465 - mrSges(5,1) * t621 + (-Ifges(5,2) - Ifges(6,3)) * t588 - t805 + (Ifges(5,4) - t747) * t590;
t452 = -pkin(11) * t465 + mrSges(5,2) * t621 + Ifges(5,1) * t590 - Ifges(5,4) * t588 - t717 * t475 + t725 * t476;
t450 = mrSges(4,2) * t683 + Ifges(4,1) * t654 + Ifges(4,4) * t657 + t726 * t452 - t718 * t453;
t449 = -pkin(5) * t739 - mrSges(4,1) * t683 + Ifges(4,4) * t654 + Ifges(4,2) * t657 + t718 * t452 + t726 * t453;
t1 = [Ifges(2,3) + t728 * (Ifges(3,1) * t728 + t727 * t449 + t719 * t450 - t715 * t457 + t723 * t458 - t714 * t467 + t722 * t472) + Ifges(7,2) * t724 ^ 2 + (pkin(1) * t804 + Ifges(3,4) * t803 + Ifges(3,2) * t720 - t719 * t449 + t727 * t450 - t723 * t457 - t715 * t458 - t722 * t467 - t714 * t472) * t720 + (Ifges(7,1) * t716 + 0.2e1 * Ifges(7,4) * t724) * t716 + (m(7) * pkin(14) - 0.2e1 * t724 * mrSges(7,1) + 0.2e1 * t716 * mrSges(7,2)) * pkin(14) + ((m(9) + m(3)) * pkin(15) + mrSges(3,2) * t803 - 0.2e1 * t720 * mrSges(3,1) - t655 * mrSges(9,2) + t652 * mrSges(9,1) - t737 - t804) * pkin(15), t728 * (t727 * t460 + t719 * t463 - t715 * t491 + t723 * t492 + Ifges(3,5) + t758) - t720 * (-pkin(1) * t536 + t719 * t460 - t727 * t463 + t723 * t491 + t715 * t492 + Ifges(3,6) + t768) - pkin(15) * t536 - t454 * t469 + t513 + (t728 * (-t715 * t504 + t723 * t505) - t720 * (t723 * t504 + t715 * t505)) * t521 - t488 * t470, t728 * (t727 * t471 + t719 * t478) - t720 * (-pkin(1) * t576 + t719 * t471 - t727 * t478) - pkin(15) * t576 + t454 * t510 + (-t720 * t768 + t728 * t758) * t514 + (t728 * (-t714 * t667 + t722 * t670) - t720 * (t722 * t667 + t714 * t670)) * t509 + t488 * t512, t728 * (t727 * t552 + t719 * t572) - t720 * (-pkin(1) * t749 + t719 * t552 - t727 * t572) - pkin(15) * t749; (-t727 * (t657 * mrSges(4,3) + t726 * t462 - t718 * t507) + t719 * (-t654 * mrSges(4,3) + t777) + t715 * (t653 * mrSges(8,3) + t752) + t723 * (-t656 * mrSges(8,3) + t753)) * pkin(1) + t733 + t740 + Ifges(3,5) * t728 - Ifges(3,6) * t720 + t521 * t740 - t470 * t503 - t469 * t451 + t513 + t746, (-t727 * (t726 * t522 - t718 * t562 - mrSges(4,2)) + t719 * (mrSges(4,1) + t775) + (t719 ^ 2 + t727 ^ 2) * pkin(1) * m(4) + ((mrSges(8,1) + t772) * t521 + mrSges(8,1) + m(8) * t791 + t774) * t723 + (t771 * t521 + m(8) * t794 + (-0.2e1 * t521 - 0.2e1) * mrSges(8,2) + t773) * t715) * pkin(1) + t734 + t731 - (t455 + t468) * t469 + t521 * t734 - (t456 + t760) * t470 + t516 ^ 2 * Ifges(7,3) + Ifges(3,3) + t738 + (0.1e1 + t521) * t521 * (t772 * t762 + t771 * t795 + Ifges(8,3) + t523), t732 + pkin(1) * (-t727 * (t726 * t566 - t718 * t612 - mrSges(4,2)) + t719 * (mrSges(4,1) + t770)) + t455 * t510 + t456 * t512 - t469 * t498 + t776, -t812 - pkin(1) * (t718 * t719 - t778) * t748 + t755; t510 * t451 + t512 * t503 + t509 * t769 + t514 * t746 + t733, (t760 - t781) * t512 + t731 + t510 * t468 - t466 * t469 + t776, t776 * t514 + (t509 * Ifges(10,3) + t514 * t809) * t509 + t512 ^ 2 * Ifges(11,3) + (t498 + t466) * t510 + t732, t755 + t811; Ifges(6,3) * t588 + t747 * t590 + t805, t767 - t807 - t812, t767 - t806 + t811, Ifges(6,3);];
Mq = t1;
