% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2TE
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2TE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:47:40
% EndTime: 2020-05-07 01:47:41
% DurationCPUTime: 1.31s
% Computational Cost: add. (2434->236), mult. (2480->301), div. (0->0), fcn. (1310->80), ass. (0->152)
t780 = 2 * qJ(2);
t865 = -sin(t780) / 0.2e1;
t760 = qJ(3) + qJ(2);
t736 = 2 * t760;
t864 = -sin(t736) / 0.2e1;
t759 = pkin(17) + pkin(18);
t814 = (pkin(15) + t759);
t807 = (pkin(16) + t814);
t720 = qJ(2) + t807;
t712 = qJ(4) + t720;
t697 = qJ(3) + t712;
t721 = -qJ(2) + t807;
t713 = -qJ(4) + t721;
t698 = -qJ(3) + t713;
t863 = -sin(t697) / 0.4e1 + sin(t698) / 0.4e1;
t770 = sin(qJ(2));
t773 = cos(qJ(3));
t815 = -pkin(4) * t773 + pkin(1);
t769 = sin(qJ(3));
t774 = cos(qJ(2));
t832 = t769 * t774;
t653 = pkin(4) * t832 - t815 * t770;
t833 = t769 * t770;
t726 = pkin(4) * t833;
t654 = t815 * t774 + t726;
t771 = sin(pkin(15));
t775 = cos(pkin(15));
t632 = t653 * t775 - t654 * t771;
t764 = sin(pkin(16));
t765 = cos(pkin(16));
t813 = t653 * t771 + t654 * t775;
t862 = t764 * t632 + t765 * t813;
t861 = t632 * t765 - t764 * t813;
t830 = t773 * t774;
t685 = t830 - t833;
t686 = t773 * t770 + t832;
t834 = t686 * t775;
t638 = t771 * t685 + t834;
t831 = t771 * t686;
t812 = t685 * t775 - t831;
t860 = -t764 * t638 + t765 * t812;
t859 = t638 * t765 + t764 * t812;
t854 = sin(t712) / 0.4e1 - sin(t713) / 0.4e1;
t853 = 2 * qJ(4);
t776 = m(5) + m(6);
t851 = pkin(4) * mrSges(6,1);
t850 = pkin(4) * mrSges(6,2);
t849 = pkin(8) * mrSges(6,1);
t848 = pkin(8) * mrSges(6,2);
t847 = cos(t712);
t846 = cos(t713);
t841 = pkin(8) * t771;
t840 = pkin(8) * t775;
t839 = t776 * pkin(4);
t768 = sin(qJ(4));
t837 = t768 * mrSges(6,1);
t772 = cos(qJ(4));
t836 = t772 * mrSges(6,2);
t708 = t772 * mrSges(6,1) - t768 * mrSges(6,2);
t835 = t653 * t708;
t829 = qJ(4) - qJ(2);
t828 = qJD(4) * t708;
t755 = -pkin(15) + pkin(14) - qJ(2);
t724 = 0.2e1 * t755;
t725 = mrSges(9,1) + mrSges(4,1) + t839;
t732 = sin(t755);
t735 = cos(t755);
t737 = -pkin(10) * m(6) + mrSges(5,2) - mrSges(6,3);
t761 = t780 + qJ(3);
t742 = sin(t761);
t746 = cos(t761);
t820 = pkin(18) + pkin(15);
t749 = qJ(2) + t820;
t750 = -qJ(2) + t820;
t758 = pkin(8) * m(6) + mrSges(5,1);
t767 = mrSges(4,2) + mrSges(9,2);
t693 = cos(t697);
t694 = cos(t698);
t714 = qJ(3) + t720;
t715 = -qJ(3) + t721;
t808 = qJ(2) + t814;
t722 = qJ(3) + t808;
t809 = -qJ(2) + t814;
t723 = -qJ(3) + t809;
t800 = -qJ(4) + t720;
t796 = qJ(3) + t800;
t799 = qJ(4) + t721;
t797 = -qJ(3) + t799;
t785 = -sin(t796) / 0.4e1 + sin(t797) / 0.4e1;
t786 = cos(t797) / 0.4e1 + cos(t796) / 0.4e1;
t782 = (Ifges(9,4) + Ifges(4,4)) * cos(t736) + ((-cos(t722) / 0.2e1 - cos(t723) / 0.2e1) * mrSges(9,2) + (sin(t723) / 0.2e1 - sin(t722) / 0.2e1) * mrSges(9,1)) * pkin(3) + (t839 * t864 + (-sin(t714) / 0.2e1 + sin(t715) / 0.2e1) * t758 + (cos(t714) / 0.2e1 - cos(t715) / 0.2e1) * t737 + t786 * mrSges(6,2) + t785 * mrSges(6,1)) * pkin(4) + (-Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2)) * t864 + t863 * t851 - (t693 + t694) * t850 / 0.4e1;
t787 = -sin(t799) / 0.4e1 + sin(t800) / 0.4e1;
t788 = -cos(t799) / 0.4e1 - cos(t800) / 0.4e1;
t741 = sin(t760);
t745 = cos(t760);
t804 = t725 * t741 + t767 * t745;
t619 = (t725 * t742 + t767 * t746 + (sin(t720) / 0.2e1 - sin(t721) / 0.2e1) * t758 + (pkin(1) * t865 - pkin(12) * t770) * (m(4) + m(8) + m(9) + t776) + (-cos(t720) / 0.2e1 + cos(t721) / 0.2e1) * t737 + (-cos(t749) / 0.2e1 + cos(t750) / 0.2e1) * mrSges(8,2) + t788 * mrSges(6,2) + (sin(t749) / 0.2e1 - sin(t750) / 0.2e1) * mrSges(8,1) + t787 * mrSges(6,1) + (-sin(t809) / 0.2e1 + sin(t808) / 0.2e1) * m(9) * pkin(3)) * pkin(1) + t782 + (-Ifges(3,1) + Ifges(3,2)) * t865 + (-mrSges(7,1) * t732 + mrSges(7,2) * t735) * pkin(6) + (Ifges(7,2) - Ifges(7,1)) * sin(t724) / 0.2e1 + (-mrSges(3,1) * t770 - mrSges(3,2) * t774 + t804) * pkin(12) + Ifges(3,4) * cos(t780) + Ifges(7,4) * cos(t724) + t854 * pkin(1) * mrSges(6,1) + (t847 + t846) * pkin(1) * mrSges(6,2) / 0.4e1;
t827 = t619 * qJD(1);
t751 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t752 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t802 = qJ(4) + t807;
t794 = 2 * t802;
t803 = -qJ(4) + t807;
t795 = 2 * t803;
t810 = 2 * t807;
t805 = qJ(4) + t810;
t806 = -qJ(4) + t810;
t620 = -(t751 + t849) * sin(t806) / 0.4e1 - (-t752 - t848) * cos(t805) / 0.4e1 + (t837 / 0.2e1 + t836 / 0.2e1) * pkin(8) + (-t751 + t849) * sin(t805) / 0.4e1 + (-t752 + t848) * cos(t806) / 0.4e1 + (cos(t853) / 0.2e1 - cos(t795) / 0.4e1 - cos(t794) / 0.4e1) * Ifges(6,4) + ((-cos(t802) / 0.2e1 - cos(t803) / 0.2e1) * mrSges(6,2) + (sin(t803) / 0.2e1 - sin(t802) / 0.2e1) * mrSges(6,1)) * pkin(12) + ((t693 / 0.4e1 + t694 / 0.4e1 + t786) * mrSges(6,2) + (t785 - t863) * mrSges(6,1)) * pkin(4) + ((-t847 / 0.4e1 - t846 / 0.4e1 + t788) * mrSges(6,2) + (t787 - t854) * mrSges(6,1)) * pkin(1) + (-sin(t795) / 0.8e1 + sin(t794) / 0.8e1 - sin(t853) / 0.4e1) * (Ifges(6,2) - Ifges(6,1));
t826 = t620 * qJD(1);
t623 = t804 * pkin(12) + ((t773 / 0.2e1 + t746 / 0.2e1) * t767 + (t769 / 0.2e1 + t742 / 0.2e1) * t725) * pkin(1) + t782;
t825 = t623 * qJD(1);
t762 = qJ(2) + qJ(4);
t743 = sin(t762);
t744 = -sin(t829);
t747 = cos(t762);
t748 = cos(t829);
t753 = qJ(4) + t760;
t730 = sin(t753);
t754 = qJ(3) - t829;
t731 = sin(t754);
t733 = cos(t753);
t734 = cos(t754);
t792 = (-t734 / 0.4e1 + t733 / 0.4e1) * mrSges(6,2) + (t731 / 0.4e1 + t730 / 0.4e1) * mrSges(6,1);
t789 = t792 * pkin(4);
t784 = t789 + ((t748 / 0.4e1 - t747 / 0.4e1) * mrSges(6,2) + (-t744 / 0.4e1 - t743 / 0.4e1) * mrSges(6,1)) * pkin(1);
t627 = t835 / 0.2e1 + t784;
t824 = t627 * qJD(1);
t629 = t741 * (Ifges(4,6) + Ifges(9,6)) + (pkin(4) * mrSges(5,3) - Ifges(4,5) - Ifges(9,5)) * t745 + (t730 / 0.2e1 - t731 / 0.2e1) * t851 + (t733 + t734) * t850 / 0.2e1;
t823 = t629 * qJD(3);
t798 = t686 * t708;
t630 = (t798 / 0.2e1 + t792) * pkin(4);
t822 = t630 * qJD(1);
t707 = -t836 - t837;
t739 = sin(t759);
t740 = cos(t759);
t819 = pkin(4) * (-t859 * t739 + t860 * t740) * t707;
t811 = -t819 / 0.2e1;
t801 = mrSges(6,1) * t841 + t752 * t775;
t790 = -t775 * t830 + t831;
t791 = t771 * t830 + t834;
t783 = (((-t764 * t771 + t775 * t765) * t740 + (-t775 * t764 - t771 * t765) * t739) * t726 + ((t764 * t791 + t790 * t765) * t740 + (-t790 * t764 + t791 * t765) * t739) * pkin(4)) * t707;
t621 = t811 - t783 / 0.2e1;
t648 = (t725 * t769 + t767 * t773) * pkin(1);
t793 = -t648 * qJD(2) + t621 * qJD(4);
t680 = mrSges(6,2) * t841 + t751 * t775;
t679 = mrSges(6,1) * t840 - t771 * t752;
t678 = mrSges(6,2) * t840 - t771 * t751;
t644 = t648 * qJD(3);
t631 = -pkin(4) * t798 / 0.2e1 + t789;
t628 = -t835 / 0.2e1 + t784;
t622 = t783 / 0.2e1 + t811;
t1 = [t619 * qJD(2) + t623 * qJD(3) - t620 * qJD(4), t628 * qJD(4) + t823 + t827 + (Ifges(3,5) * t774 + Ifges(7,5) * t735 - Ifges(3,6) * t770 + Ifges(7,6) * t732 + t629 + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t774 + (-t748 / 0.2e1 - t747 / 0.2e1) * mrSges(6,2) + (t744 / 0.2e1 - t743 / 0.2e1) * mrSges(6,1)) * pkin(1)) * qJD(2), t629 * qJD(2) + t631 * qJD(4) + t823 + t825, -t826 + t628 * qJD(2) + t631 * qJD(3) + ((-(t679 * t765 - t764 * t801) * t768 + (-t678 * t765 + t680 * t764) * t772) * t740 + (-(-t764 * t679 - t801 * t765) * t768 + (t764 * t678 + t680 * t765) * t772) * t739 - (pkin(12) + t654) * t707) * qJD(4); t627 * qJD(4) - t827, t644, t644 - t793, t824 - t621 * qJD(3) - (t862 * t739 - t861 * t740) * t828; t630 * qJD(4) - t825, t793, 0, t822 + t621 * qJD(2) + pkin(4) * (t860 * t739 + t859 * t740) * t828; -t627 * qJD(2) - t630 * qJD(3) + t826, -t824 + (t861 * t739 + t862 * t740) * t707 * qJD(2) + t622 * qJD(3), t622 * qJD(2) - qJD(3) * t819 - t822, 0;];
Cq = t1;
