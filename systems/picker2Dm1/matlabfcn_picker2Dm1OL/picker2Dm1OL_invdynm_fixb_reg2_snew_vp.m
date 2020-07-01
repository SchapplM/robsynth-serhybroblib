% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
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
%
% Output:
% m_new_reg [(3*11)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = picker2Dm1OL_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_invdynm_fixb_reg2_snew_vp: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_invdynm_fixb_reg2_snew_vp: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1OL_invdynm_fixb_reg2_snew_vp: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:52
% EndTime: 2020-05-11 05:46:14
% DurationCPUTime: 12.90s
% Computational Cost: add. (13357->354), mult. (17803->437), div. (0->0), fcn. (11392->22), ass. (0->267)
t856 = sin(qJ(1));
t865 = cos(qJ(1));
t795 = g(1) * t856 - g(2) * t865;
t784 = qJDD(1) * pkin(1) - t795;
t797 = g(1) * t865 + g(2) * t856;
t871 = qJD(1) ^ 2;
t785 = -pkin(1) * t871 + t797;
t855 = sin(qJ(2));
t864 = cos(qJ(2));
t717 = -t784 * t864 + t785 * t855;
t841 = qJDD(1) + qJDD(2);
t703 = pkin(3) * t841 - t717;
t719 = t784 * t855 + t785 * t864;
t843 = qJD(1) + qJD(2);
t839 = t843 ^ 2;
t705 = -pkin(3) * t839 + t719;
t853 = sin(qJ(4));
t862 = cos(qJ(4));
t658 = -t703 * t862 + t705 * t853;
t823 = qJDD(4) + t841;
t653 = pkin(4) * t823 - t658;
t660 = t853 * t703 + t862 * t705;
t827 = qJD(4) + t843;
t820 = t827 ^ 2;
t655 = -pkin(4) * t820 + t660;
t844 = sin(qJ(10));
t846 = cos(qJ(10));
t623 = -t653 * t846 + t655 * t844;
t624 = t653 * t844 + t655 * t846;
t606 = t623 * t846 - t624 * t844;
t896 = t623 * t844 + t624 * t846;
t593 = t606 * t853 + t862 * t896;
t595 = t606 * t862 - t853 * t896;
t589 = -t593 * t855 + t595 * t864;
t704 = pkin(2) * t841 - t717;
t706 = -pkin(2) * t839 + t719;
t854 = sin(qJ(3));
t863 = cos(qJ(3));
t659 = -t704 * t863 + t706 * t854;
t824 = qJDD(3) + t841;
t654 = pkin(6) * t824 + t659;
t661 = t854 * t704 + t863 * t706;
t828 = qJD(3) + t843;
t821 = t828 ^ 2;
t656 = -pkin(6) * t821 - t661;
t848 = sin(qJ(9));
t857 = cos(qJ(9));
t625 = -t654 * t857 + t656 * t848;
t626 = t654 * t848 + t656 * t857;
t610 = t625 * t857 - t626 * t848;
t895 = t625 * t848 + t626 * t857;
t597 = -t610 * t854 - t863 * t895;
t598 = t610 * t863 - t854 * t895;
t985 = t597 * t855 + t598 * t864;
t806 = qJDD(10) + t823;
t817 = qJD(10) + t827;
t808 = t817 ^ 2;
t734 = t806 * t844 + t808 * t846;
t736 = t806 * t846 - t808 * t844;
t675 = t734 * t862 + t736 * t853;
t908 = t734 * t853 - t736 * t862;
t972 = t675 * t855 + t864 * t908;
t980 = t675 * t864 - t855 * t908;
t984 = -t856 * t980 - t865 * t972;
t983 = t856 * t972 - t865 * t980;
t809 = qJDD(9) + t824;
t818 = qJD(9) + t828;
t811 = t818 ^ 2;
t741 = t809 * t848 + t811 * t857;
t743 = t809 * t857 - t811 * t848;
t957 = t741 * t854 - t743 * t863;
t958 = -t741 * t863 - t743 * t854;
t650 = t855 * t957 + t864 * t958;
t971 = t855 * t958 - t864 * t957;
t982 = -t650 * t856 - t865 * t971;
t981 = t650 * t865 - t856 * t971;
t629 = t658 * t853 + t660 * t862;
t633 = t658 * t862 - t660 * t853;
t618 = -t629 * t855 + t633 * t864;
t631 = t659 * t854 + t661 * t863;
t634 = t659 * t863 - t661 * t854;
t619 = -t631 * t855 + t634 * t864;
t826 = qJD(6) + t843;
t819 = t826 ^ 2;
t822 = qJDD(6) + t841;
t851 = sin(qJ(6));
t860 = cos(qJ(6));
t752 = t819 * t851 - t822 * t860;
t902 = t819 * t860 + t822 * t851;
t907 = t752 * t864 + t855 * t902;
t959 = -t752 * t855 + t864 * t902;
t970 = -t856 * t959 - t865 * t907;
t969 = t856 * t907 - t865 * t959;
t753 = t820 * t862 + t823 * t853;
t881 = t820 * t853 - t823 * t862;
t690 = t753 * t864 - t855 * t881;
t906 = -t753 * t855 - t864 * t881;
t968 = t690 * t856 - t865 * t906;
t967 = t690 * t865 + t856 * t906;
t757 = t821 * t863 + t824 * t854;
t758 = t821 * t854 - t824 * t863;
t905 = t757 * t855 + t758 * t864;
t960 = t864 * t757 - t758 * t855;
t966 = -t856 * t960 - t865 * t905;
t965 = t856 * t905 - t865 * t960;
t928 = t862 * g(3);
t933 = t853 * g(3);
t764 = t844 * t933 - t846 * t928;
t765 = (t844 * t862 + t846 * t853) * g(3);
t707 = t764 * t864 + t765 * t855;
t904 = -t764 * t855 + t765 * t864;
t964 = t707 * t856 - t865 * t904;
t963 = t707 * t865 + t856 * t904;
t833 = t854 * g(3);
t927 = t863 * g(3);
t770 = -t833 * t848 + t857 * t927;
t776 = (-t848 * t863 - t854 * t857) * g(3);
t711 = t770 * t864 + t776 * t855;
t903 = -t770 * t855 + t776 * t864;
t962 = t711 * t856 - t865 * t903;
t961 = t711 * t865 + t856 * t903;
t938 = pkin(5) * g(3);
t845 = sin(pkin(8));
t847 = cos(pkin(8));
t852 = sin(qJ(5));
t861 = cos(qJ(5));
t782 = t845 * t861 + t847 * t852;
t937 = t782 * g(3);
t783 = -t845 * t852 + t847 * t861;
t936 = t783 * g(3);
t935 = t846 * g(3);
t850 = sin(qJ(7));
t934 = t850 * g(3);
t932 = t855 * g(3);
t834 = t856 * g(3);
t931 = t857 * g(3);
t858 = cos(qJ(8));
t930 = t858 * g(3);
t929 = t860 * g(3);
t926 = t864 * g(3);
t925 = t865 * g(3);
t920 = t855 * t856;
t919 = t855 * t865;
t602 = pkin(4) * t606;
t918 = -pkin(3) * t595 - t602;
t603 = pkin(6) * t610;
t917 = pkin(2) * t598 - t603;
t916 = pkin(4) * t933;
t915 = pkin(6) * t833;
t914 = pkin(1) * t925;
t913 = t854 * t932;
t912 = g(3) * t920;
t911 = g(3) * t919;
t910 = t853 * t932;
t859 = cos(qJ(7));
t794 = g(1) * t850 - g(2) * t859;
t724 = -g(1) * t782 + g(2) * t783;
t725 = g(1) * t783 + g(2) * t782;
t909 = t724 * t782 - t725 * t783;
t665 = t717 * t860 + t719 * t851;
t849 = sin(qJ(8));
t716 = -t784 * t858 + t785 * t849;
t842 = qJD(1) + qJD(8);
t838 = t842 ^ 2;
t840 = qJDD(1) + qJDD(8);
t901 = t838 * t858 + t840 * t849;
t796 = g(1) * t859 + g(2) * t850;
t900 = -pkin(4) * t736 + t623;
t899 = -pkin(6) * t743 + t625;
t898 = -pkin(3) * t881 - t658;
t897 = pkin(2) * t758 + t659;
t668 = -t851 * t717 + t860 * t719;
t635 = t665 * t851 + t668 * t860;
t636 = t665 * t860 - t668 * t851;
t622 = -t635 * t855 + t636 * t864;
t718 = t849 * t784 + t858 * t785;
t667 = t716 * t858 - t718 * t849;
t669 = t717 * t864 - t719 * t855;
t894 = t724 * t783 + t725 * t782;
t767 = t838 * t849 - t840 * t858;
t892 = t767 * t856 - t865 * t901;
t891 = -t767 * t865 - t856 * t901;
t768 = t839 * t864 + t841 * t855;
t769 = -t839 * t855 + t841 * t864;
t890 = t768 * t865 + t769 * t856;
t889 = t768 * t856 - t769 * t865;
t771 = t851 * t932 - t860 * t926;
t777 = (t851 * t864 + t855 * t860) * g(3);
t888 = t771 * t865 + t777 * t856;
t887 = t771 * t856 - t777 * t865;
t772 = t862 * t926 - t910;
t778 = (-t853 * t864 - t855 * t862) * g(3);
t886 = t772 * t865 + t778 * t856;
t885 = t772 * t856 - t778 * t865;
t773 = -t863 * t926 + t913;
t779 = (t854 * t864 + t855 * t863) * g(3);
t884 = t773 * t865 + t779 * t856;
t883 = t773 * t856 - t779 * t865;
t882 = t794 * t859 - t796 * t850;
t870 = qJD(5) ^ 2;
t880 = qJDD(5) * t783 - t782 * t870;
t720 = qJDD(5) * t782 + t783 * t870;
t879 = pkin(3) * t908 + t900;
t878 = -pkin(2) * t957 + t899;
t877 = pkin(4) * t734 + t624;
t876 = pkin(6) * t741 + t626;
t875 = -pkin(3) * t753 - t660;
t874 = pkin(2) * t757 + t661;
t873 = pkin(3) * t675 + t877;
t872 = pkin(2) * t958 + t876;
t869 = qJD(7) ^ 2;
t868 = pkin(1) * g(3);
t867 = pkin(2) * g(3);
t866 = pkin(3) * g(3);
t835 = t859 * g(3);
t832 = t851 * g(3);
t831 = t849 * g(3);
t830 = t848 * g(3);
t829 = t844 * g(3);
t825 = pkin(1) * t834;
t803 = pkin(2) * t926 + t868;
t802 = pkin(3) * t926 + t868;
t801 = -pkin(6) * t927 + t867;
t800 = pkin(4) * t928 + t866;
t790 = qJDD(7) * t850 + t859 * t869;
t789 = -qJDD(1) * t865 + t856 * t871;
t788 = -qJDD(1) * t856 - t865 * t871;
t787 = qJDD(7) * t859 - t850 * t869;
t781 = (t856 * t864 + t919) * g(3);
t780 = (-t849 * t865 - t856 * t858) * g(3);
t775 = (-t864 * t865 + t920) * g(3);
t774 = (-t849 * t856 + t858 * t865) * g(3);
t745 = -t801 * t855 + t864 * t915;
t744 = -t800 * t855 - t864 * t916;
t730 = pkin(6) * t913 + t801 * t864 + t868;
t729 = -pkin(4) * t910 + t800 * t864 + t868;
t687 = pkin(1) * t769 - t717;
t686 = -pkin(1) * t768 - t719;
t685 = pkin(1) * t767 + t716;
t684 = pkin(1) * t901 + t718;
t666 = t717 * t855 + t719 * t864;
t664 = t716 * t849 + t718 * t858;
t663 = pkin(1) * t669;
t662 = pkin(1) * t667;
t642 = pkin(1) * t907 + t665;
t641 = pkin(1) * t959 + t668;
t640 = pkin(1) * t905 + t897;
t639 = pkin(1) * t960 + t874;
t638 = pkin(1) * t906 + t898;
t637 = -pkin(1) * t690 + t875;
t628 = pkin(2) * t634;
t627 = pkin(3) * t633;
t621 = t635 * t864 + t636 * t855;
t620 = pkin(1) * t622;
t617 = t631 * t864 + t634 * t855;
t616 = t629 * t864 + t633 * t855;
t615 = pkin(1) * t971 + t878;
t614 = pkin(1) * t650 + t872;
t613 = pkin(1) * t972 + t879;
t612 = pkin(1) * t980 + t873;
t601 = -pkin(1) * t619 - t628;
t600 = -pkin(1) * t618 - t627;
t590 = t597 * t864 - t598 * t855;
t588 = t593 * t864 + t595 * t855;
t587 = pkin(1) * t985 + t917;
t586 = -pkin(1) * t589 + t918;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t789, 0, -t788, 0, t834, t925, -t795 * t865 + t797 * t856, 0, 0, 0, t889, 0, t890, 0, t781, -t775, t666 * t856 - t669 * t865, t825, 0, 0, t966, 0, t965, 0, t883, t884, t617 * t856 - t619 * t865, pkin(2) * t911 + t803 * t856, 0, 0, t968, 0, t967, 0, t885, t886, t616 * t856 - t618 * t865, pkin(3) * t911 + t802 * t856, 0, 0, t880, 0, -t720, 0, -t937, -t936, t894, -t845 * t938, 0, 0, t970, 0, t969, 0, t887, t888, t621 * t856 - t622 * t865, t825, 0, 0, t790, 0, t787, 0, t835, -t934, -t882, 0, 0, 0, t891, 0, t892, 0, t780, -t774, t664 * t856 - t667 * t865, t825, 0, 0, t982, 0, -t981, 0, t962, t961, t590 * t856 + t865 * t985, t730 * t856 - t745 * t865, 0, 0, t984, 0, t983, 0, t964, t963, t588 * t856 - t589 * t865, t729 * t856 - t744 * t865; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t788, 0, t789, 0, -t925, t834, -t795 * t856 - t797 * t865, 0, 0, 0, -t890, 0, t889, 0, t775, t781, -t666 * t865 - t669 * t856, -t914, 0, 0, -t965, 0, t966, 0, -t884, t883, -t617 * t865 - t619 * t856, pkin(2) * t912 - t803 * t865, 0, 0, -t967, 0, t968, 0, -t886, t885, -t616 * t865 - t618 * t856, pkin(3) * t912 - t802 * t865, 0, 0, t720, 0, t880, 0, t936, -t937, t909, t847 * t938, 0, 0, -t969, 0, t970, 0, -t888, t887, -t621 * t865 - t622 * t856, -t914, 0, 0, -t787, 0, t790, 0, t934, t835, -t794 * t850 - t796 * t859, pkin(7) * g(3), 0, 0, -t892, 0, t891, 0, t774, t780, -t664 * t865 - t667 * t856, -t914, 0, 0, t981, 0, t982, 0, -t961, t962, -t590 * t865 + t856 * t985, -t730 * t865 - t745 * t856, 0, 0, -t983, 0, t984, 0, -t963, t964, -t588 * t865 - t589 * t856, -t729 * t865 - t744 * t856; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t795, -t797, 0, 0, 0, 0, 0, 0, 0, t841, t687, t686, 0, -t663, 0, 0, 0, 0, 0, t824, t640, t639, 0, t601, 0, 0, 0, 0, 0, t823, t638, t637, 0, t600, 0, 0, 0, 0, 0, qJDD(5), (t720 * t845 + t847 * t880) * pkin(5) - t724, (-t720 * t847 + t845 * t880) * pkin(5) + t725, 0, (-t845 * t909 - t847 * t894) * pkin(5), 0, 0, 0, 0, 0, t822, t642, t641, 0, -t620, 0, 0, 0, 0, 0, qJDD(7), pkin(7) * t790 - t796, pkin(7) * t787 + t794, 0, pkin(7) * t882, 0, 0, 0, 0, 0, t840, t685, t684, 0, -t662, 0, 0, 0, 0, 0, t809, t615, t614, 0, t587, 0, 0, 0, 0, 0, t806, t613, t612, 0, t586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t871, 0, 0, -g(3), t795, 0, 0, 0, t769, 0, -t768, 0, -t932, -t926, t669, 0, 0, 0, t905, 0, t960, 0, t779, -t773, t619, -pkin(2) * t932, 0, 0, t906, 0, -t690, 0, t778, -t772, t618, -pkin(3) * t932, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t907, 0, t959, 0, t777, -t771, t622, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, 0, t901, 0, t831, t930, t667, 0, 0, 0, t971, 0, t650, 0, t903, -t711, -t985, t745, 0, 0, t972, 0, t980, 0, t904, -t707, t589, t744; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t871, 0, qJDD(1), 0, g(3), 0, t797, 0, 0, 0, t768, 0, t769, 0, t926, -t932, t666, t868, 0, 0, -t960, 0, t905, 0, t773, t779, t617, t803, 0, 0, t690, 0, t906, 0, t772, t778, t616, t802, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t959, 0, t907, 0, t771, t777, t621, t868, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t901, 0, t767, 0, -t930, t831, t664, t868, 0, 0, -t650, 0, t971, 0, t711, t903, t590, t730, 0, 0, -t980, 0, t972, 0, t707, t904, t588, t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t795, -t797, 0, 0, 0, 0, 0, 0, 0, t841, t687, t686, 0, -t663, 0, 0, 0, 0, 0, t824, t640, t639, 0, t601, 0, 0, 0, 0, 0, t823, t638, t637, 0, t600, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, t642, t641, 0, -t620, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t840, t685, t684, 0, -t662, 0, 0, 0, 0, 0, t809, t615, t614, 0, t587, 0, 0, 0, 0, 0, t806, t613, t612, 0, t586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t841, 0, -t839, 0, 0, -g(3), t717, 0, 0, 0, t758, 0, t757, 0, t833, t927, t634, 0, 0, 0, -t881, 0, -t753, 0, -t933, -t928, t633, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t752, 0, t902, 0, t832, t929, t636, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t957, 0, t958, 0, t776, -t770, -t598, t915, 0, 0, t908, 0, t675, 0, t765, -t764, t595, -t916; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t839, 0, t841, 0, g(3), 0, t719, 0, 0, 0, -t757, 0, t758, 0, -t927, t833, t631, t867, 0, 0, t753, 0, -t881, 0, t928, -t933, t629, t866, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t902, 0, t752, 0, -t929, t832, t635, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t958, 0, -t957, 0, t770, t776, t597, t801, 0, 0, -t675, 0, t908, 0, t764, t765, t593, t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t841, -t717, -t719, 0, 0, 0, 0, 0, 0, 0, t824, t897, t874, 0, -t628, 0, 0, 0, 0, 0, t823, t898, t875, 0, -t627, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, t665, t668, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t809, t878, t872, 0, t917, 0, 0, 0, 0, 0, t806, t879, t873, 0, t918; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t824, 0, -t821, 0, 0, -g(3), -t659, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t743, 0, t741, 0, t830, t931, t610, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t821, 0, t824, 0, g(3), 0, -t661, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t741, 0, -t743, 0, -t931, t830, t895, pkin(6) * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t824, t659, t661, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t809, t899, t876, 0, -t603, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, 0, -t820, 0, 0, -g(3), t658, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t736, 0, t734, 0, t829, t935, t606, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t820, 0, t823, 0, g(3), 0, t660, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t734, 0, -t736, 0, -t935, t829, t896, pkin(4) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, -t658, -t660, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806, t900, t877, 0, -t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(5), 0, -t870, 0, 0, -g(3), t724, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t870, 0, qJDD(5), 0, g(3), 0, -t725, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(5), -t724, t725, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, 0, -t819, 0, 0, -g(3), -t665, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t819, 0, t822, 0, g(3), 0, -t668, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, t665, t668, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(7), 0, -t869, 0, 0, -g(3), t796, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t869, 0, qJDD(7), 0, g(3), 0, -t794, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(7), -t796, t794, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t840, 0, -t838, 0, 0, -g(3), -t716, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t838, 0, t840, 0, g(3), 0, -t718, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t840, t716, t718, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t809, 0, -t811, 0, 0, -g(3), -t625, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t811, 0, t809, 0, g(3), 0, -t626, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t809, t625, t626, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806, 0, -t808, 0, 0, -g(3), -t623, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t808, 0, t806, 0, g(3), 0, -t624, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806, t623, t624, 0, 0;];
m_new_reg = t1;