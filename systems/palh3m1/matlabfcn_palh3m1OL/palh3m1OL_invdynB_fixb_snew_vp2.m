% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = palh3m1OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1OL_invdynB_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:05:35
% EndTime: 2020-04-20 17:06:10
% DurationCPUTime: 8.07s
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
t786 = t834 * t842 + t836 * t843;
t902 = qJD(1) * qJD(2);
t851 = qJDD(1) * t880 + t888 * t902;
t900 = t880 * t902;
t853 = qJDD(1) * t888 - t900;
t791 = -qJD(7) * t836 - t851 * t875 + t853 * t883;
t793 = qJD(7) * t834 + t851 * t883 + t853 * t875;
t735 = -qJD(8) * t786 + t791 * t843 - t793 * t842;
t785 = t834 * t843 - t836 * t842;
t736 = qJD(8) * t785 + t791 * t842 + t793 * t843;
t881 = sin(qJ(1));
t889 = cos(qJ(1));
t858 = g(1) * t881 - t889 * g(2);
t844 = -qJDD(1) * pkin(12) - t858;
t813 = t844 + (-t853 + t900) * pkin(1);
t871 = qJD(2) + qJD(7);
t737 = (-t791 * t874 - t793 * t873 + (-t834 * t873 + t836 * t874) * t871) * pkin(3) + t813;
t865 = qJD(8) + t871;
t783 = -mrSges(9,2) * t865 + mrSges(9,3) * t785;
t784 = mrSges(9,1) * t865 - mrSges(9,3) * t786;
t715 = m(9) * t737 - t735 * mrSges(9,1) + t736 * mrSges(9,2) - t785 * t783 + t786 * t784;
t913 = pkin(3) * t715;
t859 = -g(1) * t889 - g(2) * t881;
t890 = qJD(1) ^ 2;
t846 = -pkin(12) * t890 + t859;
t823 = -g(3) * t888 - t846 * t880;
t849 = (-mrSges(3,1) * t888 + mrSges(3,2) * t880) * qJD(1);
t903 = qJD(1) * t888;
t857 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t903;
t816 = (t880 * t888 * t890 + qJDD(2)) * pkin(1) + t823;
t825 = -g(3) * t880 + t888 * t846;
t817 = (-t888 ^ 2 * t890 - qJD(2) ^ 2) * pkin(1) + t825;
t879 = sin(qJ(3));
t887 = cos(qJ(3));
t780 = -t816 * t887 + t879 * t817;
t835 = (t879 * t880 - t887 * t888) * qJD(1);
t794 = qJD(3) * t835 - t851 * t887 - t853 * t879;
t837 = (-t879 * t888 - t880 * t887) * qJD(1);
t812 = -mrSges(4,1) * t835 + mrSges(4,2) * t837;
t872 = qJD(2) + qJD(3);
t819 = -mrSges(4,2) * t872 + mrSges(4,3) * t835;
t870 = qJDD(2) + qJDD(3);
t756 = (t835 * t837 + t870) * pkin(4) + t780;
t782 = -t816 * t879 - t817 * t887;
t765 = (-t835 ^ 2 - t872 ^ 2) * pkin(4) + t782;
t878 = sin(qJ(4));
t886 = cos(qJ(4));
t731 = t878 * t756 + t886 * t765;
t792 = -qJD(3) * t837 + t851 * t879 - t853 * t887;
t808 = t835 * t878 + t837 * t886;
t745 = -qJD(4) * t808 + t792 * t886 - t794 * t878;
t807 = t835 * t886 - t837 * t878;
t771 = -mrSges(5,1) * t807 + mrSges(5,2) * t808;
t866 = qJD(4) + t872;
t796 = mrSges(5,1) * t866 - mrSges(5,3) * t808;
t864 = qJDD(4) + t870;
t746 = qJD(4) * t807 + t792 * t878 + t794 * t886;
t761 = t813 + (t837 * t872 - t792) * pkin(4);
t722 = (-t807 * t866 - t746) * pkin(10) + (t808 * t866 - t745) * pkin(8) + t761;
t772 = -pkin(8) * t807 - pkin(10) * t808;
t862 = t866 ^ 2;
t724 = -pkin(8) * t862 + pkin(10) * t864 + t772 * t807 + t731;
t877 = sin(qJ(5));
t885 = cos(qJ(5));
t713 = t722 * t885 - t724 * t877;
t787 = -t808 * t877 + t866 * t885;
t729 = qJD(5) * t787 + t746 * t885 + t864 * t877;
t743 = qJDD(5) - t745;
t788 = t808 * t885 + t866 * t877;
t759 = -mrSges(6,1) * t787 + mrSges(6,2) * t788;
t804 = qJD(5) - t807;
t762 = -mrSges(6,2) * t804 + mrSges(6,3) * t787;
t709 = m(6) * t713 + mrSges(6,1) * t743 - mrSges(6,3) * t729 - t759 * t788 + t762 * t804;
t714 = t722 * t877 + t724 * t885;
t728 = -qJD(5) * t788 - t746 * t877 + t864 * t885;
t763 = mrSges(6,1) * t804 - mrSges(6,3) * t788;
t710 = m(6) * t714 - mrSges(6,2) * t743 + mrSges(6,3) * t728 + t759 * t787 - t763 * t804;
t898 = -t709 * t877 + t885 * t710;
t696 = m(5) * t731 - mrSges(5,2) * t864 + mrSges(5,3) * t745 + t771 * t807 - t796 * t866 + t898;
t730 = t756 * t886 - t765 * t878;
t795 = -mrSges(5,2) * t866 + mrSges(5,3) * t807;
t723 = -pkin(8) * t864 - pkin(10) * t862 + t772 * t808 - t730;
t893 = -m(6) * t723 + t728 * mrSges(6,1) - mrSges(6,2) * t729 + t787 * t762 - t763 * t788;
t705 = m(5) * t730 + mrSges(5,1) * t864 - mrSges(5,3) * t746 - t771 * t808 + t795 * t866 + t893;
t910 = t878 * t696 + t886 * t705;
t689 = m(4) * t780 + mrSges(4,1) * t870 - mrSges(4,3) * t794 - t812 * t837 + t819 * t872 + t910;
t821 = mrSges(4,1) * t872 - mrSges(4,3) * t837;
t690 = m(4) * t782 - mrSges(4,2) * t870 + mrSges(4,3) * t792 + t696 * t886 - t705 * t878 + t812 * t835 - t821 * t872;
t779 = t883 * t816 - t817 * t875;
t811 = -mrSges(8,1) * t834 + mrSges(8,2) * t836;
t818 = -mrSges(8,2) * t871 + mrSges(8,3) * t834;
t869 = qJDD(2) + qJDD(7);
t803 = (-t834 * t874 - t836 * t873) * pkin(3);
t868 = t871 ^ 2;
t747 = -t803 * t836 + (t868 * t873 + t869 * t874) * pkin(3) + t779;
t781 = t875 * t816 + t883 * t817;
t748 = t803 * t834 + (-t868 * t874 + t869 * t873) * pkin(3) + t781;
t725 = t747 * t843 - t748 * t842;
t754 = -mrSges(9,1) * t785 + mrSges(9,2) * t786;
t863 = qJDD(8) + t869;
t720 = m(9) * t725 + mrSges(9,1) * t863 - mrSges(9,3) * t736 - t754 * t786 + t783 * t865;
t726 = t747 * t842 + t748 * t843;
t721 = m(9) * t726 - mrSges(9,2) * t863 + mrSges(9,3) * t735 + t754 * t785 - t784 * t865;
t909 = t843 * t720 + t842 * t721;
t702 = m(8) * t779 + mrSges(8,1) * t869 - mrSges(8,3) * t793 - t811 * t836 + t818 * t871 + t909;
t820 = mrSges(8,1) * t871 - mrSges(8,3) * t836;
t908 = -t842 * t720 + t843 * t721;
t703 = m(8) * t781 - mrSges(8,2) * t869 + mrSges(8,3) * t791 + t811 * t834 - t820 * t871 + t908;
t894 = t689 * t887 + t690 * t879 - t883 * t702 - t875 * t703;
t905 = qJD(1) * t880;
t685 = m(3) * t823 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t851 + qJD(2) * t857 - t849 * t905 - t894;
t855 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t905;
t686 = m(3) * t825 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t853 - qJD(2) * t855 + t689 * t879 - t690 * t887 - t702 * t875 + t703 * t883 + t849 * t903;
t847 = pkin(6) * t890 + t859;
t876 = sin(qJ(6));
t884 = cos(qJ(6));
t822 = -g(3) * t884 - t847 * t876;
t848 = (-mrSges(7,1) * t884 + mrSges(7,2) * t876) * qJD(1);
t901 = qJD(1) * qJD(6);
t850 = qJDD(1) * t876 + t884 * t901;
t904 = qJD(1) * t884;
t856 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t904;
t906 = qJD(1) * t876;
t777 = m(7) * t822 + qJDD(6) * mrSges(7,1) - mrSges(7,3) * t850 + qJD(6) * t856 - t848 * t906;
t824 = -g(3) * t876 + t847 * t884;
t852 = qJDD(1) * t884 - t876 * t901;
t854 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t906;
t778 = m(7) * t824 - qJDD(6) * mrSges(7,2) + mrSges(7,3) * t852 - qJD(6) * t854 + t848 * t904;
t897 = -t777 * t876 + t884 * t778;
t680 = m(2) * t859 - mrSges(2,1) * t890 - qJDD(1) * mrSges(2,2) - t685 * t880 + t686 * t888 + t897;
t697 = t885 * t709 + t877 * t710;
t895 = m(5) * t761 - t745 * mrSges(5,1) + t746 * mrSges(5,2) - t807 * t795 + t808 * t796 + t697;
t892 = t792 * mrSges(4,1) + t791 * mrSges(8,1) + (-m(4) - m(8)) * t813 + t835 * t819 - t895 - t836 * t820 - t837 * t821 - t793 * mrSges(8,2) - t794 * mrSges(4,2) + t834 * t818 - t715;
t891 = -m(3) * t844 + t853 * mrSges(3,1) - t851 * mrSges(3,2) + t857 * t903 + t892;
t845 = qJDD(1) * pkin(6) - t858;
t896 = -m(7) * t845 + t852 * mrSges(7,1) - t850 * mrSges(7,2) + t856 * t904;
t692 = (-t854 * t876 - t855 * t880) * qJD(1) + m(2) * t858 + qJDD(1) * mrSges(2,1) - t890 * mrSges(2,2) + t891 + t896;
t912 = t881 * t680 + t889 * t692;
t911 = t888 * t685 + t880 * t686;
t907 = t884 * t777 + t876 * t778;
t899 = t889 * t680 - t692 * t881;
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
t712 = mrSges(9,2) * t737 - mrSges(9,3) * t725 + Ifges(9,1) * t736 + Ifges(9,4) * t735 + Ifges(9,5) * t863 + t749 * t785 - t750 * t865;
t711 = -mrSges(9,1) * t737 + mrSges(9,3) * t726 + Ifges(9,4) * t736 + Ifges(9,2) * t735 + Ifges(9,6) * t863 - t749 * t786 + t751 * t865;
t699 = mrSges(6,2) * t723 - mrSges(6,3) * t713 + Ifges(6,1) * t729 + Ifges(6,4) * t728 + Ifges(6,5) * t743 + t738 * t787 - t739 * t804;
t698 = -mrSges(6,1) * t723 + mrSges(6,3) * t714 + Ifges(6,4) * t729 + Ifges(6,2) * t728 + Ifges(6,6) * t743 - t738 * t788 + t740 * t804;
t694 = mrSges(8,2) * t813 - mrSges(8,3) * t779 + Ifges(8,1) * t793 + Ifges(8,4) * t791 + Ifges(8,5) * t869 - t711 * t842 + t712 * t843 + t797 * t834 - t799 * t871 - t873 * t913;
t693 = -mrSges(8,1) * t813 + mrSges(8,3) * t781 + Ifges(8,4) * t793 + Ifges(8,2) * t791 + Ifges(8,6) * t869 + t711 * t843 + t712 * t842 - t797 * t836 + t801 * t871 - t874 * t913;
t688 = -pkin(8) * t697 - mrSges(5,1) * t761 - mrSges(6,1) * t713 + mrSges(6,2) * t714 + mrSges(5,3) * t731 + Ifges(5,4) * t746 - Ifges(6,5) * t729 + Ifges(5,2) * t745 + Ifges(5,6) * t864 - Ifges(6,6) * t728 - Ifges(6,3) * t743 - t739 * t788 + t740 * t787 - t766 * t808 + t768 * t866;
t687 = -pkin(10) * t697 + mrSges(5,2) * t761 - mrSges(5,3) * t730 + Ifges(5,1) * t746 + Ifges(5,4) * t745 + Ifges(5,5) * t864 - t698 * t877 + t699 * t885 + t766 * t807 - t767 * t866;
t682 = mrSges(4,2) * t813 - mrSges(4,3) * t780 + Ifges(4,1) * t794 + Ifges(4,4) * t792 + Ifges(4,5) * t870 + t687 * t886 - t688 * t878 + t798 * t835 - t800 * t872;
t681 = -pkin(4) * t895 - mrSges(4,1) * t813 + mrSges(4,3) * t782 + Ifges(4,4) * t794 + Ifges(4,2) * t792 + Ifges(4,6) * t870 + t878 * t687 + t886 * t688 - t837 * t798 + t872 * t802;
t677 = mrSges(3,2) * t844 - mrSges(3,3) * t823 + Ifges(3,1) * t851 + Ifges(3,4) * t853 + Ifges(3,5) * qJDD(2) - qJD(2) * t831 + t681 * t879 - t682 * t887 - t693 * t875 + t694 * t883 + t829 * t903;
t676 = t892 * pkin(1) - mrSges(3,1) * t844 + mrSges(3,3) * t825 + Ifges(3,4) * t851 + Ifges(3,2) * t853 + Ifges(3,6) * qJDD(2) + qJD(2) * t833 - t887 * t681 - t879 * t682 + t883 * t693 + t875 * t694 - t829 * t905;
t675 = (-t830 * t876 - t831 * t880 + t832 * t884 + t833 * t888) * qJD(1) + pkin(13) * t897 - pkin(10) * t898 - pkin(4) * t910 - pkin(12) * t911 + pkin(6) * t907 + (-t873 * t908 - t874 * t909) * pkin(3) - pkin(8) * t893 + t894 * pkin(1) - Ifges(3,5) * t851 - Ifges(7,6) * t852 - Ifges(3,6) * t853 + mrSges(2,3) * t859 - Ifges(7,5) * t850 + t835 * t802 - t836 * t799 - t837 * t800 + t834 * t801 - mrSges(7,1) * t822 - mrSges(3,1) * t823 + mrSges(7,2) * t824 + mrSges(3,2) * t825 - t808 * t767 + t807 * t768 - Ifges(8,6) * t791 - Ifges(4,6) * t792 + Ifges(2,6) * qJDD(1) - Ifges(8,5) * t793 - Ifges(4,5) * t794 + mrSges(8,2) * t781 + mrSges(4,2) * t782 + t785 * t751 - t786 * t750 - mrSges(8,1) * t779 - mrSges(4,1) * t780 - Ifges(5,5) * t746 - Ifges(5,6) * t745 - Ifges(9,5) * t736 - Ifges(9,6) * t735 - mrSges(5,1) * t730 + mrSges(5,2) * t731 - mrSges(9,1) * t725 + mrSges(9,2) * t726 - Ifges(3,3) * qJDD(2) - Ifges(7,3) * qJDD(6) + mrSges(2,1) * g(3) - Ifges(9,3) * t863 - Ifges(5,3) * t864 + t890 * Ifges(2,5) - t885 * t698 - t877 * t699 - Ifges(8,3) * t869 - Ifges(4,3) * t870;
t674 = -mrSges(2,2) * g(3) - mrSges(2,3) * t858 - pkin(13) * t776 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t890 - t676 * t880 + t677 * t888 - t757 * t876 + t758 * t884;
t1 = [-m(1) * g(1) + t899; -m(1) * g(2) + t912; (-m(1) - m(2)) * g(3) + t907 + t911; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(11) * t912 + t889 * t674 - t881 * t675; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(11) * t899 + t881 * t674 + t889 * t675; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t858 - mrSges(2,2) * t859 + t880 * t677 + t888 * t676 + pkin(12) * (-t855 * t905 + t891) + t876 * t758 + t884 * t757 - pkin(6) * t776;];
tauB = t1;
