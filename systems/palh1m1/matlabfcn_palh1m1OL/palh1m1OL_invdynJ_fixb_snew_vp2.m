% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tauJ [13x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = palh1m1OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynJ_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:29:51
% EndTime: 2020-04-15 19:30:05
% DurationCPUTime: 4.52s
% Computational Cost: add. (21962->522), mult. (48336->681), div. (0->0), fcn. (37204->22), ass. (0->201)
t796 = sin(qJ(7));
t801 = sin(qJ(2));
t805 = cos(qJ(7));
t810 = cos(qJ(2));
t755 = (-t796 * t801 + t805 * t810) * qJD(1);
t831 = qJD(1) * qJD(2);
t829 = t810 * t831;
t764 = -qJDD(1) * t801 - t829;
t766 = qJDD(1) * t810 - t801 * t831;
t696 = -qJD(7) * t755 + t764 * t805 - t766 * t796;
t800 = sin(qJ(3));
t809 = cos(qJ(3));
t756 = (-t800 * t801 + t809 * t810) * qJD(1);
t697 = qJD(3) * t756 - t764 * t809 + t766 * t800;
t752 = (-t796 * t810 - t801 * t805) * qJD(1);
t699 = qJD(7) * t752 + t764 * t796 + t766 * t805;
t753 = (t800 * t810 + t801 * t809) * qJD(1);
t700 = -qJD(3) * t753 + t764 * t800 + t766 * t809;
t802 = sin(qJ(1));
t811 = cos(qJ(1));
t828 = g(1) * t802 - t811 * g(2);
t759 = -qJDD(1) * pkin(15) - t828;
t727 = t759 + (-t764 + t829) * pkin(1);
t789 = qJD(2) + qJD(7);
t736 = -mrSges(8,2) * t789 + mrSges(8,3) * t752;
t790 = qJD(2) + qJD(3);
t737 = mrSges(4,1) * t790 - mrSges(4,3) * t753;
t738 = mrSges(8,1) * t789 - mrSges(8,3) * t755;
t739 = -mrSges(4,2) * t790 + mrSges(4,3) * t756;
t799 = sin(qJ(4));
t808 = cos(qJ(4));
t721 = t753 * t808 + t756 * t799;
t639 = -qJD(4) * t721 - t697 * t799 + t700 * t808;
t719 = -t753 * t799 + t756 * t808;
t641 = qJD(4) * t719 + t697 * t808 + t700 * t799;
t655 = t727 + (t753 * t790 - t700) * pkin(5);
t781 = qJD(4) + t790;
t605 = (-t719 * t781 - t641) * pkin(11) + (t721 * t781 - t639) * pkin(9) + t655;
t812 = qJD(1) ^ 2;
t826 = -g(1) * t811 - g(2) * t802;
t761 = -pkin(15) * t812 + t826;
t743 = t801 * g(3) - t761 * t810;
t731 = (-t801 * t810 * t812 + qJDD(2)) * pkin(1) + t743;
t741 = -g(3) * t810 - t761 * t801;
t838 = t801 ^ 2;
t732 = (-qJD(2) ^ 2 - t812 * t838) * pkin(1) + t741;
t680 = t800 * t731 + t809 * t732;
t787 = qJDD(2) + qJDD(3);
t652 = (t753 * t756 + t787) * pkin(5) + t680;
t678 = -t731 * t809 + t800 * t732;
t659 = (-t756 ^ 2 - t790 ^ 2) * pkin(5) + t678;
t616 = t799 * t652 + t808 * t659;
t671 = -pkin(9) * t719 - pkin(11) * t721;
t776 = t781 ^ 2;
t778 = qJDD(4) + t787;
t607 = -pkin(9) * t776 + pkin(11) * t778 + t671 * t719 + t616;
t798 = sin(qJ(5));
t807 = cos(qJ(5));
t595 = t605 * t807 - t607 * t798;
t686 = -t721 * t798 + t781 * t807;
t613 = qJD(5) * t686 + t641 * t807 + t778 * t798;
t636 = qJDD(5) - t639;
t687 = t721 * t807 + t781 * t798;
t653 = -mrSges(6,1) * t686 + mrSges(6,2) * t687;
t717 = qJD(5) - t719;
t656 = -mrSges(6,2) * t717 + mrSges(6,3) * t686;
t588 = m(6) * t595 + mrSges(6,1) * t636 - t613 * mrSges(6,3) - t653 * t687 + t656 * t717;
t596 = t605 * t798 + t607 * t807;
t612 = -qJD(5) * t687 - t641 * t798 + t778 * t807;
t657 = mrSges(6,1) * t717 - mrSges(6,3) * t687;
t589 = m(6) * t596 - mrSges(6,2) * t636 + t612 * mrSges(6,3) + t653 * t686 - t657 * t717;
t574 = t807 * t588 + t798 * t589;
t704 = -mrSges(5,2) * t781 + mrSges(5,3) * t719;
t706 = mrSges(5,1) * t781 - mrSges(5,3) * t721;
t823 = -m(5) * t655 + t639 * mrSges(5,1) - t641 * mrSges(5,2) + t719 * t704 - t721 * t706 - t574;
t791 = sin(pkin(19));
t792 = cos(pkin(19));
t793 = cos(qJ(10));
t835 = sin(qJ(10));
t757 = t791 * t793 - t792 * t835;
t758 = -t791 * t835 - t792 * t793;
t682 = t752 * t757 + t755 * t758;
t621 = -qJD(10) * t682 + t696 * t758 - t699 * t757;
t681 = t752 * t758 - t755 * t757;
t622 = qJD(10) * t681 + t696 * t757 + t699 * t758;
t624 = (-t696 * t792 - t699 * t791 + (-t752 * t791 + t755 * t792) * t789) * pkin(4) + t727;
t779 = qJD(10) + t789;
t675 = -mrSges(11,2) * t779 + mrSges(11,3) * t681;
t676 = mrSges(11,1) * t779 - mrSges(11,3) * t682;
t825 = -m(11) * t624 + t621 * mrSges(11,1) - t622 * mrSges(11,2) + t681 * t675 - t682 * t676;
t839 = t700 * mrSges(4,1) + t696 * mrSges(8,1) - t697 * mrSges(4,2) - t699 * mrSges(8,2) + t752 * t736 - t753 * t737 - t755 * t738 + t756 * t739 + t823 + t825 + (-m(4) - m(8)) * t727;
t837 = pkin(4) * t791;
t836 = pkin(4) * t792;
t670 = -mrSges(5,1) * t719 + mrSges(5,2) * t721;
t827 = -t588 * t798 + t807 * t589;
t572 = m(5) * t616 - mrSges(5,2) * t778 + mrSges(5,3) * t639 + t670 * t719 - t706 * t781 + t827;
t615 = t652 * t808 - t659 * t799;
t606 = -pkin(9) * t778 - pkin(11) * t776 + t671 * t721 - t615;
t824 = -m(6) * t606 + t612 * mrSges(6,1) - t613 * mrSges(6,2) + t686 * t656 - t657 * t687;
t582 = m(5) * t615 + mrSges(5,1) * t778 - mrSges(5,3) * t641 - t670 * t721 + t704 * t781 + t824;
t834 = t799 * t572 + t808 * t582;
t677 = t805 * t731 - t732 * t796;
t707 = (-t752 * t792 - t755 * t791) * pkin(4);
t784 = t789 ^ 2;
t786 = qJDD(2) + qJDD(7);
t642 = -t707 * t755 + (t784 * t791 + t786 * t792) * pkin(4) + t677;
t679 = t796 * t731 + t805 * t732;
t643 = t707 * t752 + (-t784 * t792 + t786 * t791) * pkin(4) + t679;
t609 = t642 * t758 - t643 * t757;
t650 = -mrSges(11,1) * t681 + mrSges(11,2) * t682;
t772 = qJDD(10) + t786;
t601 = m(11) * t609 + mrSges(11,1) * t772 - t622 * mrSges(11,3) - t650 * t682 + t675 * t779;
t610 = t642 * t757 + t643 * t758;
t602 = m(11) * t610 - mrSges(11,2) * t772 + t621 * mrSges(11,3) + t650 * t681 - t676 * t779;
t833 = t758 * t601 + t757 * t602;
t832 = -t757 * t601 + t758 * t602;
t795 = sin(qJ(8));
t804 = cos(qJ(8));
t702 = t804 * t741 + t795 * t743;
t830 = qJD(1) * qJD(6);
t788 = qJD(2) + qJD(8);
t785 = qJDD(2) + qJDD(8);
t701 = -t741 * t795 + t804 * t743;
t751 = (-t795 * t810 - t801 * t804) * qJD(1);
t754 = (-t795 * t801 + t804 * t810) * qJD(1);
t666 = (t751 * t754 + t785) * pkin(2) + t701;
t672 = (-t751 ^ 2 - t788 ^ 2) * pkin(2) + t702;
t794 = sin(qJ(9));
t803 = cos(qJ(9));
t626 = -t666 * t803 + t672 * t794;
t627 = -t666 * t794 - t672 * t803;
t695 = -qJD(8) * t754 + t764 * t804 - t766 * t795;
t698 = qJD(8) * t751 + t764 * t795 + t766 * t804;
t720 = -t751 * t794 - t754 * t803;
t638 = -qJD(9) * t720 - t695 * t803 + t698 * t794;
t718 = -t751 * t803 + t754 * t794;
t640 = qJD(9) * t718 - t695 * t794 - t698 * t803;
t780 = qJD(9) + t788;
t662 = Ifges(10,4) * t720 + Ifges(10,2) * t718 + Ifges(10,6) * t780;
t664 = Ifges(10,1) * t720 + Ifges(10,4) * t718 + Ifges(10,5) * t780;
t777 = qJDD(9) + t785;
t822 = mrSges(10,1) * t626 - mrSges(10,2) * t627 + Ifges(10,5) * t640 + Ifges(10,6) * t638 + Ifges(10,3) * t777 + t720 * t662 - t718 * t664;
t645 = Ifges(11,4) * t682 + Ifges(11,2) * t681 + Ifges(11,6) * t779;
t646 = Ifges(11,1) * t682 + Ifges(11,4) * t681 + Ifges(11,5) * t779;
t821 = mrSges(11,1) * t609 - mrSges(11,2) * t610 + Ifges(11,5) * t622 + Ifges(11,6) * t621 + Ifges(11,3) * t772 + t682 * t645 - t681 * t646;
t668 = (t754 * t788 - t695) * pkin(2) + t759;
t703 = -mrSges(10,2) * t780 + mrSges(10,3) * t718;
t705 = mrSges(10,1) * t780 - mrSges(10,3) * t720;
t820 = m(10) * t668 - t638 * mrSges(10,1) + t640 * mrSges(10,2) - t718 * t703 + t720 * t705;
t628 = Ifges(6,5) * t687 + Ifges(6,6) * t686 + Ifges(6,3) * t717;
t630 = Ifges(6,1) * t687 + Ifges(6,4) * t686 + Ifges(6,5) * t717;
t578 = -mrSges(6,1) * t606 + mrSges(6,3) * t596 + Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t636 - t628 * t687 + t630 * t717;
t629 = Ifges(6,4) * t687 + Ifges(6,2) * t686 + Ifges(6,6) * t717;
t579 = mrSges(6,2) * t606 - mrSges(6,3) * t595 + Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t636 + t628 * t686 - t629 * t717;
t663 = Ifges(5,4) * t721 + Ifges(5,2) * t719 + Ifges(5,6) * t781;
t665 = Ifges(5,1) * t721 + Ifges(5,4) * t719 + Ifges(5,5) * t781;
t819 = pkin(9) * t824 + pkin(11) * t827 + mrSges(5,1) * t615 - mrSges(5,2) * t616 + Ifges(5,5) * t641 + Ifges(5,6) * t639 + Ifges(5,3) * t778 + t807 * t578 + t798 * t579 + t721 * t663 - t719 * t665;
t818 = mrSges(6,1) * t595 - mrSges(6,2) * t596 + Ifges(6,5) * t613 + Ifges(6,6) * t612 + Ifges(6,3) * t636 + t629 * t687 - t630 * t686;
t669 = -mrSges(10,1) * t718 + mrSges(10,2) * t720;
t711 = Ifges(9,4) * t754 + Ifges(9,2) * t751 + Ifges(9,6) * t788;
t714 = Ifges(9,1) * t754 + Ifges(9,4) * t751 + Ifges(9,5) * t788;
t816 = -mrSges(9,2) * t702 - t751 * t714 + pkin(2) * (-t794 * (m(10) * t627 - mrSges(10,2) * t777 + mrSges(10,3) * t638 + t669 * t718 - t705 * t780) - t803 * (m(10) * t626 + mrSges(10,1) * t777 - mrSges(10,3) * t640 - t669 * t720 + t703 * t780)) + t754 * t711 + Ifges(9,6) * t695 + Ifges(9,5) * t698 + mrSges(9,1) * t701 + Ifges(9,3) * t785 + t822;
t712 = Ifges(8,4) * t755 + Ifges(8,2) * t752 + Ifges(8,6) * t789;
t715 = Ifges(8,1) * t755 + Ifges(8,4) * t752 + Ifges(8,5) * t789;
t815 = mrSges(8,1) * t677 - mrSges(8,2) * t679 + Ifges(8,5) * t699 + Ifges(8,6) * t696 + Ifges(8,3) * t786 + t755 * t712 - t752 * t715 + t832 * t837 + t833 * t836 + t821;
t713 = Ifges(4,4) * t753 + Ifges(4,2) * t756 + Ifges(4,6) * t790;
t716 = Ifges(4,1) * t753 + Ifges(4,4) * t756 + Ifges(4,5) * t790;
t813 = pkin(5) * t834 + mrSges(4,1) * t680 - mrSges(4,2) * t678 + Ifges(4,5) * t697 + Ifges(4,6) * t700 + Ifges(4,3) * t787 + t753 * t713 - t756 * t716 + t819;
t806 = cos(qJ(6));
t797 = sin(qJ(6));
t765 = qJDD(1) * t806 - t797 * t830;
t763 = qJDD(1) * t797 + t806 * t830;
t762 = pkin(14) * t812 + t826;
t760 = qJDD(1) * pkin(14) - t828;
t750 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t810 - Ifges(3,4) * t801) * qJD(1);
t749 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t797 + Ifges(7,4) * t806) * qJD(1);
t748 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t810 - Ifges(3,2) * t801) * qJD(1);
t747 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t797 + Ifges(7,2) * t806) * qJD(1);
t742 = -g(3) * t797 + t762 * t806;
t740 = -g(3) * t806 - t762 * t797;
t726 = -mrSges(8,1) * t752 + mrSges(8,2) * t755;
t725 = -mrSges(4,1) * t756 + mrSges(4,2) * t753;
t710 = Ifges(4,5) * t753 + Ifges(4,6) * t756 + Ifges(4,3) * t790;
t709 = Ifges(8,5) * t755 + Ifges(8,6) * t752 + Ifges(8,3) * t789;
t708 = Ifges(9,5) * t754 + Ifges(9,6) * t751 + Ifges(9,3) * t788;
t661 = Ifges(5,5) * t721 + Ifges(5,6) * t719 + Ifges(5,3) * t781;
t660 = Ifges(10,5) * t720 + Ifges(10,6) * t718 + Ifges(10,3) * t780;
t644 = Ifges(11,5) * t682 + Ifges(11,6) * t681 + Ifges(11,3) * t779;
t604 = mrSges(10,2) * t668 - mrSges(10,3) * t626 + Ifges(10,1) * t640 + Ifges(10,4) * t638 + Ifges(10,5) * t777 + t660 * t718 - t662 * t780;
t603 = -mrSges(10,1) * t668 + mrSges(10,3) * t627 + Ifges(10,4) * t640 + Ifges(10,2) * t638 + Ifges(10,6) * t777 - t660 * t720 + t664 * t780;
t591 = mrSges(11,2) * t624 - mrSges(11,3) * t609 + Ifges(11,1) * t622 + Ifges(11,4) * t621 + Ifges(11,5) * t772 + t644 * t681 - t645 * t779;
t590 = -mrSges(11,1) * t624 + mrSges(11,3) * t610 + Ifges(11,4) * t622 + Ifges(11,2) * t621 + Ifges(11,6) * t772 - t644 * t682 + t646 * t779;
t580 = mrSges(9,2) * t759 - mrSges(9,3) * t701 + Ifges(9,1) * t698 + Ifges(9,4) * t695 + Ifges(9,5) * t785 + t603 * t794 - t604 * t803 + t708 * t751 - t711 * t788;
t575 = -pkin(2) * t820 - mrSges(9,1) * t759 + mrSges(9,3) * t702 + Ifges(9,4) * t698 + Ifges(9,2) * t695 + Ifges(9,6) * t785 - t803 * t603 - t794 * t604 - t754 * t708 + t788 * t714;
t570 = mrSges(8,2) * t727 - mrSges(8,3) * t677 + Ifges(8,1) * t699 + Ifges(8,4) * t696 + Ifges(8,5) * t786 - t590 * t757 + t591 * t758 + t709 * t752 - t712 * t789 + t825 * t837;
t569 = -mrSges(8,1) * t727 + mrSges(8,3) * t679 + Ifges(8,4) * t699 + Ifges(8,2) * t696 + Ifges(8,6) * t786 + t590 * t758 + t591 * t757 - t709 * t755 + t715 * t789 + t825 * t836;
t567 = -pkin(9) * t574 - mrSges(5,1) * t655 + mrSges(5,3) * t616 + Ifges(5,4) * t641 + Ifges(5,2) * t639 + Ifges(5,6) * t778 - t661 * t721 + t665 * t781 - t818;
t566 = -pkin(11) * t574 + mrSges(5,2) * t655 - mrSges(5,3) * t615 + Ifges(5,1) * t641 + Ifges(5,4) * t639 + Ifges(5,5) * t778 - t578 * t798 + t579 * t807 + t661 * t719 - t663 * t781;
t565 = mrSges(4,2) * t727 - mrSges(4,3) * t680 + Ifges(4,1) * t697 + Ifges(4,4) * t700 + Ifges(4,5) * t787 + t566 * t808 - t567 * t799 + t710 * t756 - t713 * t790;
t564 = pkin(5) * t823 - mrSges(4,1) * t727 + mrSges(4,3) * t678 + Ifges(4,4) * t697 + Ifges(4,2) * t700 + Ifges(4,6) * t787 + t799 * t566 + t808 * t567 - t753 * t710 + t790 * t716;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t828 - mrSges(2,2) * t826 + t810 * (mrSges(3,2) * t759 - mrSges(3,3) * t743 + Ifges(3,1) * t766 + Ifges(3,4) * t764 + Ifges(3,5) * qJDD(2) - qJD(2) * t748 + t564 * t809 + t565 * t800 - t569 * t796 + t570 * t805 - t575 * t795 + t580 * t804) - t801 * (t839 * pkin(1) - mrSges(3,1) * t759 + mrSges(3,3) * t741 + Ifges(3,4) * t766 + Ifges(3,2) * t764 + Ifges(3,6) * qJDD(2) + qJD(2) * t750 + t800 * t564 - t809 * t565 + t805 * t569 + t796 * t570 + t804 * t575 + t795 * t580) + pkin(15) * (-t754 * (mrSges(9,1) * t788 - mrSges(9,3) * t754) + t751 * (-mrSges(9,2) * t788 + mrSges(9,3) * t751) + t764 * mrSges(3,1) - t766 * mrSges(3,2) + t695 * mrSges(9,1) - t698 * mrSges(9,2) - t820 + (-m(3) - m(9)) * t759 + t839) + t797 * (mrSges(7,2) * t760 - mrSges(7,3) * t740 + Ifges(7,1) * t763 + Ifges(7,4) * t765 + Ifges(7,5) * qJDD(6) - qJD(6) * t747) + t806 * (-mrSges(7,1) * t760 + mrSges(7,3) * t742 + Ifges(7,4) * t763 + Ifges(7,2) * t765 + Ifges(7,6) * qJDD(6) + qJD(6) * t749) - pkin(14) * (-m(7) * t760 + t765 * mrSges(7,1) - t763 * mrSges(7,2)) + ((-pkin(14) * (t797 ^ 2 + t806 ^ 2) * mrSges(7,3) + pkin(15) * (t810 ^ 2 + t838) * mrSges(3,3)) * qJD(1) - pkin(14) * (-mrSges(7,1) * t797 - mrSges(7,2) * t806) * qJD(6) + pkin(15) * (-mrSges(3,1) * t810 + mrSges(3,2) * t801) * qJD(2)) * qJD(1); (t810 * t748 + t801 * t750) * qJD(1) + t815 + t816 + t813 + Ifges(3,6) * t764 + Ifges(3,5) * t766 - mrSges(3,2) * t741 + mrSges(3,1) * t743 + (-t809 * (m(4) * t678 - mrSges(4,2) * t787 + mrSges(4,3) * t700 + t572 * t808 - t582 * t799 + t725 * t756 - t737 * t790) + t800 * (m(4) * t680 + mrSges(4,1) * t787 - mrSges(4,3) * t697 - t725 * t753 + t739 * t790 + t834) + t796 * (m(8) * t679 - mrSges(8,2) * t786 + mrSges(8,3) * t696 + t726 * t752 - t738 * t789 + t832) + t805 * (m(8) * t677 + mrSges(8,1) * t786 - mrSges(8,3) * t699 - t726 * t755 + t736 * t789 + t833)) * pkin(1) + Ifges(3,3) * qJDD(2); t813; t819; t818; mrSges(7,1) * t740 - mrSges(7,2) * t742 + Ifges(7,5) * t763 + Ifges(7,6) * t765 + Ifges(7,3) * qJDD(6) + (t797 * t747 - t806 * t749) * qJD(1); t815; t816; t822; t821; 0; 0; 0;];
tauJ = t1;
