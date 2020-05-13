% Calculate vector of inverse dynamics base forces with Newton-Euler for
% palh3m2OL
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = palh3m2OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_invdynB_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:35:17
% EndTime: 2020-05-07 04:35:52
% DurationCPUTime: 7.83s
% Computational Cost: add. (70483->490), mult. (155057->635), div. (0->0), fcn. (119130->18), ass. (0->192)
t914 = sin(qJ(8));
t875 = sin(qJ(7));
t880 = sin(qJ(2));
t883 = cos(qJ(7));
t888 = cos(qJ(2));
t834 = (-t875 * t880 + t883 * t888) * qJD(1);
t836 = (t875 * t888 + t880 * t883) * qJD(1);
t873 = sin(pkin(15));
t874 = cos(pkin(15));
t882 = cos(qJ(8));
t842 = t873 * t882 - t874 * t914;
t843 = -t873 * t914 - t874 * t882;
t786 = t842 * t834 + t843 * t836;
t902 = qJD(1) * qJD(2);
t851 = t880 * qJDD(1) + t888 * t902;
t900 = t880 * t902;
t853 = t888 * qJDD(1) - t900;
t791 = -t836 * qJD(7) - t875 * t851 + t883 * t853;
t793 = t834 * qJD(7) + t883 * t851 + t875 * t853;
t735 = -t786 * qJD(8) + t843 * t791 - t842 * t793;
t785 = t843 * t834 - t842 * t836;
t736 = t785 * qJD(8) + t842 * t791 + t843 * t793;
t881 = sin(qJ(1));
t889 = cos(qJ(1));
t858 = t881 * g(1) - t889 * g(2);
t844 = -qJDD(1) * pkin(12) - t858;
t813 = t844 + (-t853 + t900) * pkin(1);
t871 = qJD(2) + qJD(7);
t737 = (-t791 * t874 - t793 * t873 + (-t834 * t873 + t836 * t874) * t871) * pkin(3) + t813;
t865 = qJD(8) + t871;
t783 = -t865 * mrSges(9,2) + t785 * mrSges(9,3);
t784 = t865 * mrSges(9,1) - t786 * mrSges(9,3);
t715 = m(9) * t737 - t735 * mrSges(9,1) + t736 * mrSges(9,2) - t785 * t783 + t786 * t784;
t913 = pkin(3) * t715;
t859 = -t889 * g(1) - t881 * g(2);
t890 = qJD(1) ^ 2;
t846 = -t890 * pkin(12) + t859;
t823 = -t888 * g(3) - t880 * t846;
t849 = (-mrSges(3,1) * t888 + mrSges(3,2) * t880) * qJD(1);
t903 = qJD(1) * t888;
t857 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t903;
t816 = (t880 * t888 * t890 + qJDD(2)) * pkin(1) + t823;
t825 = -t880 * g(3) + t888 * t846;
t817 = (-t888 ^ 2 * t890 - qJD(2) ^ 2) * pkin(1) + t825;
t879 = sin(qJ(3));
t887 = cos(qJ(3));
t780 = -t887 * t816 + t879 * t817;
t835 = (t879 * t880 - t887 * t888) * qJD(1);
t794 = t835 * qJD(3) - t887 * t851 - t879 * t853;
t837 = (-t879 * t888 - t880 * t887) * qJD(1);
t812 = -t835 * mrSges(4,1) + t837 * mrSges(4,2);
t872 = qJD(2) + qJD(3);
t819 = -t872 * mrSges(4,2) + t835 * mrSges(4,3);
t870 = qJDD(2) + qJDD(3);
t756 = (t835 * t837 + t870) * pkin(4) + t780;
t782 = -t879 * t816 - t887 * t817;
t765 = (-t835 ^ 2 - t872 ^ 2) * pkin(4) + t782;
t878 = sin(qJ(4));
t886 = cos(qJ(4));
t731 = t878 * t756 + t886 * t765;
t792 = -t837 * qJD(3) + t879 * t851 - t887 * t853;
t808 = t878 * t835 + t886 * t837;
t745 = -t808 * qJD(4) + t886 * t792 - t878 * t794;
t807 = t886 * t835 - t878 * t837;
t771 = -t807 * mrSges(5,1) + t808 * mrSges(5,2);
t866 = qJD(4) + t872;
t796 = t866 * mrSges(5,1) - t808 * mrSges(5,3);
t864 = qJDD(4) + t870;
t746 = t807 * qJD(4) + t878 * t792 + t886 * t794;
t761 = t813 + (t837 * t872 - t792) * pkin(4);
t722 = (-t807 * t866 - t746) * pkin(10) + (t808 * t866 - t745) * pkin(8) + t761;
t772 = -t807 * pkin(8) - t808 * pkin(10);
t862 = t866 ^ 2;
t724 = -t862 * pkin(8) + t864 * pkin(10) + t807 * t772 + t731;
t877 = sin(qJ(5));
t885 = cos(qJ(5));
t713 = t885 * t722 - t877 * t724;
t787 = -t877 * t808 + t885 * t866;
t729 = t787 * qJD(5) + t885 * t746 + t877 * t864;
t743 = qJDD(5) - t745;
t788 = t885 * t808 + t877 * t866;
t759 = -t787 * mrSges(6,1) + t788 * mrSges(6,2);
t804 = qJD(5) - t807;
t762 = -t804 * mrSges(6,2) + t787 * mrSges(6,3);
t709 = m(6) * t713 + t743 * mrSges(6,1) - t729 * mrSges(6,3) - t788 * t759 + t804 * t762;
t714 = t877 * t722 + t885 * t724;
t728 = -t788 * qJD(5) - t877 * t746 + t885 * t864;
t763 = t804 * mrSges(6,1) - t788 * mrSges(6,3);
t710 = m(6) * t714 - t743 * mrSges(6,2) + t728 * mrSges(6,3) + t787 * t759 - t804 * t763;
t898 = -t877 * t709 + t885 * t710;
t696 = m(5) * t731 - t864 * mrSges(5,2) + t745 * mrSges(5,3) + t807 * t771 - t866 * t796 + t898;
t730 = t886 * t756 - t878 * t765;
t795 = -t866 * mrSges(5,2) + t807 * mrSges(5,3);
t723 = -t864 * pkin(8) - t862 * pkin(10) + t808 * t772 - t730;
t893 = -m(6) * t723 + t728 * mrSges(6,1) - t729 * mrSges(6,2) + t787 * t762 - t788 * t763;
t705 = m(5) * t730 + t864 * mrSges(5,1) - t746 * mrSges(5,3) - t808 * t771 + t866 * t795 + t893;
t910 = t878 * t696 + t886 * t705;
t689 = m(4) * t780 + t870 * mrSges(4,1) - t794 * mrSges(4,3) - t837 * t812 + t872 * t819 + t910;
t821 = t872 * mrSges(4,1) - t837 * mrSges(4,3);
t690 = m(4) * t782 - t870 * mrSges(4,2) + t792 * mrSges(4,3) + t886 * t696 - t878 * t705 + t835 * t812 - t872 * t821;
t779 = t883 * t816 - t875 * t817;
t811 = -t834 * mrSges(8,1) + t836 * mrSges(8,2);
t818 = -t871 * mrSges(8,2) + t834 * mrSges(8,3);
t869 = qJDD(2) + qJDD(7);
t803 = (-t834 * t874 - t836 * t873) * pkin(3);
t868 = t871 ^ 2;
t747 = -t836 * t803 + (t868 * t873 + t869 * t874) * pkin(3) + t779;
t781 = t875 * t816 + t883 * t817;
t748 = t834 * t803 + (-t868 * t874 + t869 * t873) * pkin(3) + t781;
t725 = t843 * t747 - t842 * t748;
t754 = -t785 * mrSges(9,1) + t786 * mrSges(9,2);
t863 = qJDD(8) + t869;
t720 = m(9) * t725 + t863 * mrSges(9,1) - t736 * mrSges(9,3) - t786 * t754 + t865 * t783;
t726 = t842 * t747 + t843 * t748;
t721 = m(9) * t726 - t863 * mrSges(9,2) + t735 * mrSges(9,3) + t785 * t754 - t865 * t784;
t909 = t843 * t720 + t842 * t721;
t702 = m(8) * t779 + t869 * mrSges(8,1) - t793 * mrSges(8,3) - t836 * t811 + t871 * t818 + t909;
t820 = t871 * mrSges(8,1) - t836 * mrSges(8,3);
t908 = -t842 * t720 + t843 * t721;
t703 = m(8) * t781 - t869 * mrSges(8,2) + t791 * mrSges(8,3) + t834 * t811 - t871 * t820 + t908;
t894 = t887 * t689 + t879 * t690 - t883 * t702 - t875 * t703;
t905 = qJD(1) * t880;
t685 = m(3) * t823 + qJDD(2) * mrSges(3,1) - t851 * mrSges(3,3) + qJD(2) * t857 - t849 * t905 - t894;
t855 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t905;
t686 = m(3) * t825 - qJDD(2) * mrSges(3,2) + t853 * mrSges(3,3) - qJD(2) * t855 + t879 * t689 - t887 * t690 - t875 * t702 + t883 * t703 + t849 * t903;
t847 = t890 * pkin(6) + t859;
t876 = sin(qJ(6));
t884 = cos(qJ(6));
t822 = -t884 * g(3) - t876 * t847;
t848 = (-mrSges(7,1) * t884 + mrSges(7,2) * t876) * qJD(1);
t901 = qJD(1) * qJD(6);
t850 = t876 * qJDD(1) + t884 * t901;
t904 = qJD(1) * t884;
t856 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t904;
t906 = qJD(1) * t876;
t777 = m(7) * t822 + qJDD(6) * mrSges(7,1) - t850 * mrSges(7,3) + qJD(6) * t856 - t848 * t906;
t824 = -t876 * g(3) + t884 * t847;
t852 = t884 * qJDD(1) - t876 * t901;
t854 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t906;
t778 = m(7) * t824 - qJDD(6) * mrSges(7,2) + t852 * mrSges(7,3) - qJD(6) * t854 + t848 * t904;
t897 = -t876 * t777 + t884 * t778;
t680 = m(2) * t859 - t890 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t880 * t685 + t888 * t686 + t897;
t697 = t885 * t709 + t877 * t710;
t895 = m(5) * t761 - t745 * mrSges(5,1) + t746 * mrSges(5,2) - t807 * t795 + t808 * t796 + t697;
t892 = t835 * t819 + t791 * mrSges(8,1) + t792 * mrSges(4,1) + t834 * t818 + (-m(4) - m(8)) * t813 - t793 * mrSges(8,2) - t794 * mrSges(4,2) - t836 * t820 - t837 * t821 - t895 - t715;
t891 = -m(3) * t844 + t853 * mrSges(3,1) - t851 * mrSges(3,2) + t857 * t903 + t892;
t845 = qJDD(1) * pkin(6) - t858;
t896 = -m(7) * t845 + t852 * mrSges(7,1) - t850 * mrSges(7,2) + t856 * t904;
t692 = t891 + (-t876 * t854 - t880 * t855) * qJD(1) - t890 * mrSges(2,2) + qJDD(1) * mrSges(2,1) + m(2) * t858 + t896;
t912 = t881 * t680 + t889 * t692;
t911 = t888 * t685 + t880 * t686;
t907 = t884 * t777 + t876 * t778;
t899 = t889 * t680 - t881 * t692;
t833 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t880 + Ifges(3,4) * t888) * qJD(1);
t832 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t876 + Ifges(7,4) * t884) * qJD(1);
t831 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t880 + Ifges(3,2) * t888) * qJD(1);
t830 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t876 + Ifges(7,2) * t884) * qJD(1);
t829 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t880 + Ifges(3,6) * t888) * qJD(1);
t828 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t876 + Ifges(7,6) * t884) * qJD(1);
t802 = Ifges(4,1) * t837 + Ifges(4,4) * t835 + Ifges(4,5) * t872;
t801 = Ifges(8,1) * t836 + Ifges(8,4) * t834 + Ifges(8,5) * t871;
t800 = Ifges(4,4) * t837 + Ifges(4,2) * t835 + Ifges(4,6) * t872;
t799 = Ifges(8,4) * t836 + Ifges(8,2) * t834 + Ifges(8,6) * t871;
t798 = Ifges(4,5) * t837 + Ifges(4,6) * t835 + Ifges(4,3) * t872;
t797 = Ifges(8,5) * t836 + Ifges(8,6) * t834 + Ifges(8,3) * t871;
t776 = -t854 * t906 + t896;
t768 = Ifges(5,1) * t808 + Ifges(5,4) * t807 + Ifges(5,5) * t866;
t767 = Ifges(5,4) * t808 + Ifges(5,2) * t807 + Ifges(5,6) * t866;
t766 = Ifges(5,5) * t808 + Ifges(5,6) * t807 + Ifges(5,3) * t866;
t758 = mrSges(7,2) * t845 - mrSges(7,3) * t822 + Ifges(7,1) * t850 + Ifges(7,4) * t852 + Ifges(7,5) * qJDD(6) - qJD(6) * t830 + t828 * t904;
t757 = -mrSges(7,1) * t845 + mrSges(7,3) * t824 + Ifges(7,4) * t850 + Ifges(7,2) * t852 + Ifges(7,6) * qJDD(6) + qJD(6) * t832 - t828 * t906;
t751 = Ifges(9,1) * t786 + Ifges(9,4) * t785 + Ifges(9,5) * t865;
t750 = Ifges(9,4) * t786 + Ifges(9,2) * t785 + Ifges(9,6) * t865;
t749 = Ifges(9,5) * t786 + Ifges(9,6) * t785 + Ifges(9,3) * t865;
t740 = Ifges(6,1) * t788 + Ifges(6,4) * t787 + Ifges(6,5) * t804;
t739 = Ifges(6,4) * t788 + Ifges(6,2) * t787 + Ifges(6,6) * t804;
t738 = Ifges(6,5) * t788 + Ifges(6,6) * t787 + Ifges(6,3) * t804;
t712 = mrSges(9,2) * t737 - mrSges(9,3) * t725 + Ifges(9,1) * t736 + Ifges(9,4) * t735 + Ifges(9,5) * t863 + t785 * t749 - t865 * t750;
t711 = -mrSges(9,1) * t737 + mrSges(9,3) * t726 + Ifges(9,4) * t736 + Ifges(9,2) * t735 + Ifges(9,6) * t863 - t786 * t749 + t865 * t751;
t699 = mrSges(6,2) * t723 - mrSges(6,3) * t713 + Ifges(6,1) * t729 + Ifges(6,4) * t728 + Ifges(6,5) * t743 + t787 * t738 - t804 * t739;
t698 = -mrSges(6,1) * t723 + mrSges(6,3) * t714 + Ifges(6,4) * t729 + Ifges(6,2) * t728 + Ifges(6,6) * t743 - t788 * t738 + t804 * t740;
t694 = mrSges(8,2) * t813 - mrSges(8,3) * t779 + Ifges(8,1) * t793 + Ifges(8,4) * t791 + Ifges(8,5) * t869 - t842 * t711 + t843 * t712 + t834 * t797 - t871 * t799 - t873 * t913;
t693 = -mrSges(8,1) * t813 + mrSges(8,3) * t781 + Ifges(8,4) * t793 + Ifges(8,2) * t791 + Ifges(8,6) * t869 + t843 * t711 + t842 * t712 - t836 * t797 + t871 * t801 - t874 * t913;
t688 = -pkin(8) * t697 - mrSges(5,1) * t761 - mrSges(6,1) * t713 + mrSges(6,2) * t714 + mrSges(5,3) * t731 + Ifges(5,4) * t746 - Ifges(6,5) * t729 + Ifges(5,2) * t745 + Ifges(5,6) * t864 - Ifges(6,6) * t728 - Ifges(6,3) * t743 - t788 * t739 + t787 * t740 - t808 * t766 + t866 * t768;
t687 = -pkin(10) * t697 + mrSges(5,2) * t761 - mrSges(5,3) * t730 + Ifges(5,1) * t746 + Ifges(5,4) * t745 + Ifges(5,5) * t864 - t877 * t698 + t885 * t699 + t807 * t766 - t866 * t767;
t682 = mrSges(4,2) * t813 - mrSges(4,3) * t780 + Ifges(4,1) * t794 + Ifges(4,4) * t792 + Ifges(4,5) * t870 + t886 * t687 - t878 * t688 + t835 * t798 - t872 * t800;
t681 = -pkin(4) * t895 - mrSges(4,1) * t813 + mrSges(4,3) * t782 + Ifges(4,4) * t794 + Ifges(4,2) * t792 + Ifges(4,6) * t870 + t878 * t687 + t886 * t688 - t837 * t798 + t872 * t802;
t677 = mrSges(3,2) * t844 - mrSges(3,3) * t823 + Ifges(3,1) * t851 + Ifges(3,4) * t853 + Ifges(3,5) * qJDD(2) - qJD(2) * t831 + t879 * t681 - t887 * t682 - t875 * t693 + t883 * t694 + t829 * t903;
t676 = t892 * pkin(1) - mrSges(3,1) * t844 + mrSges(3,3) * t825 + Ifges(3,4) * t851 + Ifges(3,2) * t853 + Ifges(3,6) * qJDD(2) + qJD(2) * t833 - t887 * t681 - t879 * t682 + t883 * t693 + t875 * t694 - t829 * t905;
t675 = -pkin(12) * t911 - pkin(8) * t893 + t894 * pkin(1) + t890 * Ifges(2,5) + pkin(6) * t907 + (-t873 * t908 - t874 * t909) * pkin(3) - pkin(4) * t910 + pkin(13) * t897 - pkin(10) * t898 + (-t876 * t830 - t880 * t831 + t884 * t832 + t888 * t833) * qJD(1) - t885 * t698 - t877 * t699 - Ifges(7,5) * t850 - Ifges(3,5) * t851 - Ifges(7,6) * t852 - Ifges(3,6) * t853 - mrSges(7,1) * t822 - mrSges(3,1) * t823 + mrSges(7,2) * t824 + mrSges(3,2) * t825 + t807 * t768 - t808 * t767 - Ifges(8,5) * t793 - Ifges(4,5) * t794 - Ifges(8,6) * t791 - Ifges(4,6) * t792 + mrSges(8,2) * t781 + mrSges(4,2) * t782 + t785 * t751 - t786 * t750 - mrSges(8,1) * t779 - mrSges(4,1) * t780 + Ifges(2,6) * qJDD(1) - Ifges(5,6) * t745 - Ifges(5,5) * t746 - Ifges(9,6) * t735 - Ifges(9,5) * t736 + mrSges(5,2) * t731 - mrSges(5,1) * t730 - mrSges(9,1) * t725 + mrSges(9,2) * t726 - t836 * t799 - t837 * t800 + t835 * t802 + t834 * t801 - Ifges(3,3) * qJDD(2) - Ifges(7,3) * qJDD(6) + mrSges(2,1) * g(3) - Ifges(8,3) * t869 - Ifges(4,3) * t870 + mrSges(2,3) * t859 - Ifges(9,3) * t863 - Ifges(5,3) * t864;
t674 = -mrSges(2,2) * g(3) - mrSges(2,3) * t858 - pkin(13) * t776 + Ifges(2,5) * qJDD(1) - t890 * Ifges(2,6) - t880 * t676 + t888 * t677 - t876 * t757 + t884 * t758;
t1 = [-m(1) * g(1) + t899; -m(1) * g(2) + t912; (-m(1) - m(2)) * g(3) + t907 + t911; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(11) * t912 + t889 * t674 - t881 * t675; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(11) * t899 + t881 * t674 + t889 * t675; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t858 - mrSges(2,2) * t859 + t880 * t677 + t888 * t676 + pkin(12) * (-t855 * t905 + t891) + t876 * t758 + t884 * t757 - pkin(6) * t776;];
tauB = t1;
