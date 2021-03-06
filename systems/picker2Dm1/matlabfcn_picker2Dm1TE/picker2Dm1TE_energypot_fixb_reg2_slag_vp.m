% Calculate inertial parameters regressor of potential energy for
% picker2Dm1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = picker2Dm1TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 03:28:53
% EndTime: 2020-05-10 03:29:03
% DurationCPUTime: 9.79s
% Computational Cost: add. (86703->465), mult. (233051->647), div. (3590->12), fcn. (61716->14), ass. (0->274)
t789 = sin(qJ(1));
t791 = cos(qJ(1));
t1007 = 0.1e1 / pkin(3);
t828 = (pkin(3) ^ 2);
t1004 = 4 * t828;
t833 = (pkin(1) ^ 2);
t1002 = 2 * t833;
t990 = t1007 / 0.2e1;
t815 = 0.1e1 / t990;
t831 = t833 ^ 2;
t836 = pkin(7) ^ 2;
t759 = t791 ^ 2;
t995 = sin(qJ(2));
t932 = t995 * pkin(7);
t905 = 0.2e1 * t932;
t895 = pkin(3) * t905;
t880 = t895 + t836;
t850 = t995 ^ 2;
t973 = t828 * t850;
t935 = 0.2e1 * t973;
t866 = -t828 + t880 + t935;
t863 = t759 * t866;
t823 = pkin(4) ^ 2;
t873 = t833 + t880;
t869 = -t823 + t873;
t925 = pkin(3) * t995;
t891 = t925 + pkin(7);
t875 = t891 * t791;
t900 = -0.4e1 * t925;
t790 = cos(qJ(2));
t978 = t790 * t789;
t930 = pkin(3) * t978;
t901 = pkin(1) * t930;
t951 = t823 - t836;
t910 = -(2 * t828) + t951;
t752 = t833 + t836;
t911 = -t823 + t752;
t931 = pkin(1) * t978;
t948 = t833 - t836;
t700 = sqrt(-0.4e1 * t833 * t863 - 0.4e1 * pkin(1) * (0.2e1 * (-t931 + t932) * pkin(3) + t911) * t875 + 0.4e1 * t869 * t901 + 0.4e1 * t948 * t973 + pkin(7) * t911 * t900 - t831 + t910 * t1002 - (t836 - (t815 + pkin(4)) * pkin(4)) * (t836 + (t815 - pkin(4)) * pkin(4)));
t817 = 0.2e1 * pkin(1);
t864 = t875 - t930;
t862 = 0.1e1 / (t864 * t817 + t828 + t873);
t1001 = 3 * t833;
t887 = t1001 - t910;
t977 = t790 * t791;
t929 = pkin(3) * t977;
t857 = t862 * ((t789 * t891 + t929) * t700 - ((t905 - 0.4e1 * t931) * pkin(3) + t887) * t875 + (t895 + t887) * t930 + (-0.2e1 * t863 + t935 - t895 - t1004 - t911) * pkin(1));
t856 = t1007 * t857;
t854 = -t856 / 0.2e1;
t747 = pkin(3) * t790;
t1005 = 2 * t828;
t867 = t1005 + t869;
t748 = pkin(1) * t791;
t941 = 0.2e1 * t748;
t861 = t862 * ((pkin(1) + t864) * t700 + (t866 * t941 + t891 * t867) * t789 + (t867 * t791 + (0.4e1 * t759 - 0.2e1) * pkin(1) * t891) * t747);
t860 = t1007 * t861;
t692 = t789 * t854 - t791 * t860 / 0.2e1;
t1000 = t692 / 0.2e1;
t816 = 2 * pkin(2);
t820 = pkin(5) ^ 2;
t786 = sin(pkin(8));
t787 = cos(pkin(8));
t720 = -t786 * t789 - t787 * t791;
t993 = pkin(5) * t720;
t938 = pkin(1) * t993;
t711 = t820 - t938;
t715 = -pkin(1) + t993;
t963 = -0.2e1 * t938 + t820;
t705 = sqrt(-(-(t817 + pkin(5)) * pkin(5) + t963) * (pkin(5) * (t817 - pkin(5)) + t963));
t719 = t786 * t791 - t787 * t789;
t987 = t705 * t719;
t701 = -pkin(5) * t987 - 0.2e1 * t711 * t715;
t1006 = 0.2e1 * t719;
t703 = pkin(5) * t711 * t1006 - t705 * t715;
t859 = t861 / 0.4e1;
t858 = t1007 * t859;
t710 = 0.1e1 / (t833 + t963);
t986 = t710 / pkin(5);
t686 = (-t703 * t856 / 0.4e1 + t701 * t858) * t986;
t792 = cos(pkin(9));
t855 = t857 / 0.4e1;
t853 = (t1007 * t701 * t855 + t703 * t858) * t986;
t996 = sin(pkin(9));
t683 = t996 * t686 + t792 * t853;
t992 = pkin(6) * t683;
t682 = t992 * t816;
t818 = pkin(6) ^ 2;
t967 = t682 + t818;
t673 = sqrt(-(-(t816 + pkin(6)) * pkin(6) + t967) * (pkin(6) * (t816 - pkin(6)) + t967));
t680 = t682 + 0.2e1 * t818;
t681 = -pkin(2) - t992;
t684 = t686 * t792 - t996 * t853;
t991 = pkin(6) * t684;
t665 = t673 * t991 - t680 * t681;
t666 = -t673 * t681 - t680 * t991;
t998 = t789 / 0.2e1;
t693 = t791 * t854 + t860 * t998;
t819 = 0.1e1 / pkin(6);
t988 = 0.1e1 / ((pkin(2) ^ 2) + t967) * t819;
t1009 = (t693 * t666 / 0.2e1 + t665 * t1000) * t988;
t1016 = t1009 * t683;
t824 = 0.1e1 / pkin(4);
t1003 = 4 * t831;
t769 = -t823 / 0.4e1;
t1008 = t769 + t828 / 0.2e1;
t841 = t828 ^ 2;
t822 = t823 ^ 2;
t835 = t836 ^ 2;
t949 = t831 + t835;
t813 = 0.2e1 * t836;
t954 = t813 - t823;
t970 = t836 * t823;
t870 = t954 * t833 + t822 / 0.6e1 + t949 - t970;
t718 = -t841 / 0.6e1 + t870;
t780 = -t828 / 0.3e1;
t743 = t780 + t836;
t892 = -0.2e1 * t901;
t721 = t743 * t892;
t727 = t748 + t891;
t751 = -(3 * t828) + t836;
t758 = t791 * t759;
t837 = pkin(1) * t833;
t968 = t837 * t758;
t937 = pkin(7) * t968;
t907 = 0.8e1 * t937;
t730 = t751 * t907;
t750 = -t823 - t828;
t812 = 0.3e1 * t836;
t739 = t812 + t750;
t983 = t739 * t833;
t731 = 0.10e2 * t983;
t776 = 0.4e1 / 0.3e1 * t828;
t770 = -t823 / 0.3e1;
t915 = t770 + t752;
t732 = t776 + t915;
t771 = -t823 / 0.2e1;
t909 = t828 + t752;
t734 = t771 + t909;
t735 = -t823 + t909;
t737 = t748 + pkin(7);
t738 = pkin(7) * t941;
t811 = 0.4e1 * t836;
t741 = (t811 + t823) * t833;
t744 = -t833 / 0.3e1 + t836;
t745 = 0.10e2 / 0.3e1 * t833;
t746 = t752 ^ 2;
t749 = -0.30e2 * t823 + 0.60e2 * t836;
t754 = -(3 * t833) + t836;
t768 = -t823 / 0.6e1;
t777 = 0.2e1 / 0.3e1 * t828;
t782 = 0.4e1 / 0.3e1 * t833;
t784 = t833 / 0.2e1;
t793 = 15 * t831;
t794 = 15 * t833;
t795 = 10 * t833;
t800 = -0.2e1 * t823;
t801 = -0.5e1 * t823;
t802 = 5 * t841;
t803 = 7 * t831;
t804 = 5 * t831;
t805 = 7 * t833;
t806 = 6 * t833;
t809 = 0.3e1 * t835;
t810 = 0.8e1 * t836;
t840 = pkin(3) * t828;
t825 = t840 ^ 2;
t845 = pkin(7) * t836;
t868 = 0.5e1 / 0.6e1 * t841 + t870;
t878 = t836 - t901;
t961 = t822 / 0.2e1 - t841 / 0.2e1;
t889 = -0.3e1 * t970 + t809 + t961;
t894 = -0.6e1 * t901;
t773 = -0.3e1 / 0.2e1 * t823;
t960 = t773 + t812;
t964 = t752 * ((t773 + t813) * t833 - 0.3e1 / 0.2e1 * t970 + t949 + t961) + t825;
t879 = ((t745 + t954) * t828 + t868) * t894 + (t793 + (-0.9e1 * t823 + 0.18e2 * t836) * t833 + t889) * t828 + (t794 + t960) * t841 + t964;
t893 = -0.4e1 * t901;
t881 = t734 * t893;
t945 = t836 + t1001;
t908 = t828 + t945;
t994 = pkin(1) * t789;
t883 = -((3 * t828) + t752) * t994 + t908 * t747;
t772 = -0.2e1 / 0.3e1 * t823;
t781 = -0.2e1 / 0.3e1 * t828;
t913 = t772 + t752;
t955 = t795 + t813;
t959 = t781 + t836;
t884 = -(t802 + ((5 * t833) + t739) * t1005 + (t781 + t913) * t752) * t994 + (t841 + (t772 + t781 + t955) * t828 + t804 + 0.2e1 * t983 + t836 * (t772 + t959)) * t747;
t953 = t822 - t841;
t888 = -0.6e1 * t970 + 0.6e1 * t835 + t953;
t914 = t772 + t777 + t813;
t962 = (t777 + t913) * t752 + t841;
t890 = t732 * t893 + (t806 + t914) * t828 + t962;
t897 = t968 * t747;
t972 = t831 * t759 ^ 2;
t898 = t972 * t747;
t926 = 0.16e2 * t968;
t904 = pkin(7) * t926;
t906 = 0.20e2 / 0.3e1 * t833;
t952 = -t823 + t828;
t912 = t812 + t952;
t778 = t828 / 0.3e1;
t916 = t768 + t778 + t836;
t917 = t823 / 0.3e1 + t778 + t813;
t918 = 0.2e1 / 0.3e1 * t823 + t777 + t811;
t919 = 0.4e1 / 0.3e1 * t823 + t776 - 0.2e1 * t836;
t976 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t920 = t789 * t976;
t939 = 0.6e1 * t748;
t922 = pkin(7) * t939;
t940 = 0.4e1 * t748;
t923 = pkin(7) * t940;
t924 = -t994 / 0.2e1;
t971 = t833 * t759;
t927 = 0.12e2 * t971;
t928 = t833 * t747;
t933 = 0.4e1 * t971;
t934 = 0.8e1 * t972;
t979 = t995 * t850 * t840;
t936 = -0.8e1 * t979;
t942 = 0.2e1 * t994;
t943 = pkin(7) * t748;
t944 = 0.4e1 * pkin(7);
t946 = t835 + t841;
t947 = t835 - t831;
t950 = -t828 + t836;
t956 = 0.4e1 / 0.7e1 * t836 - t823 / 0.7e1;
t957 = t784 + t836;
t958 = t833 / 0.3e1 + t836;
t969 = t836 * t833;
t975 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t980 = t850 ^ 2 * t841;
t981 = t750 * t836;
t982 = t746 * (-t828 + t911);
t984 = (-t789 * t837 + t928) * t759;
t989 = ((-0.24e2 * (0.4e1 / 0.3e1 * t971 + t738 + t744) * t980 * t994 - 0.12e2 * (-0.8e1 / 0.3e1 * t898 + ((t782 + t916) * t747 - (0.7e1 / 0.6e1 * t828 + t768 + t957) * t994) * t933 + (-t828 * t948 - 0.5e1 / 0.3e1 * t831 + t917 * t833 + t836 * (t770 + t743)) * t747 + (-t841 + (-t906 + t918) * t828 - (3 * t831) + t919 * t833 + t835) * t924 + (-t789 * t831 * t758 + ((t833 + t916) * t747 + (t1005 - t948) * t924) * t748) * t944) * t973 + 0.24e2 * t743 * t898 + ((t836 + 0.5e1 / 0.2e1 * t828 + 0.3e1 / 0.2e1 * t833 + t771) * t747 + t751 * t994 / 0.2e1) * t904 - 0.6e1 * ((-(3 * t841) + (-t906 + t919) * t828 + t918 * t833 + t947) * t747 - 0.2e1 * (-0.5e1 / 0.3e1 * t841 + (-t833 + t917) * t828 + t836 * (t780 + t915)) * t994) * t971 - 0.6e1 * t884 * t943 - (t825 + ((21 * t833) + t739) * t841 + (t731 + t809 + (35 * t831) + 0.2e1 * t981) * t828 + (t803 + (t801 + t810 - (5 * t828)) * t833 + t836 * (-t823 + t950)) * t752) * t747 + (0.7e1 * t825 + (t805 + t739) * t802 + (t731 + (21 * t831) + 0.9e1 * t835 + 0.6e1 * t981) * t828 + t982) * t994) * t700 + (0.16e2 * (t934 + t904 + (-(8 * t831) + 0.12e2 * t969) * t759 + (0.4e1 * pkin(1) * t845 - 0.12e2 * pkin(7) * t837) * t791 - 0.6e1 * t969 + t949) * t980 + 0.24e2 * (t959 * t934 + 0.14e2 * (-0.32e2 / 0.21e2 * (t836 + t828 / 0.4e1 + t833 / 0.4e1 - t823 / 0.8e1) * t901 + 0.5e1 / 0.42e2 * t841 + (0.16e2 / 0.21e2 * t833 + t956) * t828 + t831 / 0.7e1 + t956 * t833 + t835 - 0.3e1 / 0.7e1 * t970 + t822 / 0.42e2) * t971 + t744 * t881 - t948 * t841 + (t741 - 0.10e2 / 0.3e1 * t831 + 0.2e1 * t835 - t970) * t828 + t718 * t975 + ((-0.2e1 / 0.3e1 * t901 + t769 + t957) * t926 + (-0.8e1 / 0.3e1 * (t958 + t1008) * t901 + 0.5e1 / 0.18e2 * t841 + (0.4e1 / 0.3e1 * t836 + t782 + t770) * t828 + t835 + 0.2e1 / 0.3e1 * t969 - 0.2e1 / 0.3e1 * t970 - t831 / 0.3e1 + t822 / 0.18e2) * t939) * pkin(7)) * t973 + 0.16e2 * (-0.6e1 * t836 * t828 + t946) * t972 + 0.32e2 * (t734 * t751 + t892 * t976) * t937 + 0.24e2 * (t743 * t881 - t825 + (-t745 + t951) * t841 + (t741 + t841 / 0.6e1 - t822 / 0.6e1 + t947) * t828 + t718 * t836) * t971 + 0.8e1 * t879 * t943 - 0.8e1 * ((t805 + t960) * t841 + (t803 + (t801 + 0.10e2 * t836) * t833 + t889) * t828 + t964) * t901 + (t841 ^ 2) + (t800 + t811 + (28 * t833)) * t825 + (t749 * t833 + (70 * t831) + t888) * t841 + (t749 * t831 + t888 * t806 + t953 * t813 - 0.6e1 * t835 * t823 + 0.28e2 * t837 ^ 2 + 0.4e1 * t845 ^ 2) * t828 + t735 * t982) * t727 + (((0.4e1 * t984 + (t747 + t942) * t738 + t754 * t747 + (t771 + t908) * t942) * t936 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t828 + t784 + t768) * t790 * t815 + pkin(1) * t920) * t971 + (-0.8e1 * t897 + ((t770 + t777 + t945) * t747 - (0.8e1 / 0.3e1 * t828 + t915) * t994) * t940) * pkin(7) + t884) * t925) * t700 + (0.32e2 * (t907 + (-0.4e1 * t837 * t930 + t1003 + (t1004 + t800 + t810) * t833) * t759 + (-t833 + t878 + t1008) * t923 + t892 * t975 + t754 * t734) * t979 + 0.8e1 * (t730 + (t734 * t976 + t721) * t927 + (t881 + (t806 + t954) * t828 + t868) * t922 + t879) * t925) * t727) * t737) / ((-0.4e1 * (-t948 * t747 + 0.2e1 * t984 + (0.2e1 * pkin(7) * t929 + t789 * (t828 + t1002)) * pkin(1)) * t973 + 0.8e1 * pkin(7) * t897 + ((pkin(3) * t1003 + 0.8e1 * t833 * t840) * t790 + 0.4e1 * t837 * t920) * t759 - 0.4e1 * t883 * t943 - (t955 * t828 + t804 + t946 + 0.6e1 * t969) * t747 + (t802 + (t795 + 0.6e1 * t836) * t828 + t746) * t994) * t700 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t901 + 0.4e1 / 0.9e1 * t828 - t823 / 0.9e1 + t958) * t971 + t744 * t892 + t732 * t975 + (t968 + (t768 + t777 + t878) * t748) * t944) * t973 + t730 + (t732 * t976 + t721) * t927 + t890 * t922 + ((t745 + t914) * t828 + t962) * t894 + t825 + (t794 + t912) * t841 + (t912 * t806 + t952 * t813 + t793 + t809) * t828 + t746 * t735) * t727 + ((t936 * t994 + (-0.2e1 * t759 * t928 + (t747 - t994) * t738 + t883) * t900) * t700 + (0.8e1 * (t738 + t933 + t754) * t979 + 0.6e1 * (t950 * t933 + (t732 + t892) * t923 + t890) * t925) * t727) * t737);
t896 = t989 * t990;
t1010 = (t693 * t700 * t990 + t692 * t896) * t824;
t974 = t824 / pkin(3) ^ 2;
t674 = (t855 * t989 - t700 * t861 / 0.4e1) * t974;
t675 = (t700 * t855 + t859 * t989) * t974;
t722 = -t995 * t789 - t977;
t723 = -t995 * t791 + t978;
t667 = t674 * t722 + t675 * t723;
t668 = t674 * t723 - t675 * t722;
t874 = (-t692 * t1007 * t700 / 0.2e1 + t693 * t896) * t824;
t1015 = t1010 * t667 + t668 * t874;
t1014 = t1010 * t668 - t667 * t874;
t872 = (t666 * t1000 - t693 * t665 / 0.2e1) * t988;
t1011 = t683 * t872;
t999 = -t786 / 0.2e1;
t997 = t791 / 0.2e1;
t985 = t710 / pkin(1);
t966 = t693 * pkin(3) - t748;
t965 = t693 * pkin(2) - t748;
t921 = t673 * t819 / pkin(2);
t903 = t692 * pkin(3) - t994;
t902 = t692 * pkin(2) - t994;
t899 = t693 * t683 + t684 * t692;
t886 = -t921 / 0.2e1;
t885 = t921 / 0.2e1;
t882 = g(1) * t791 + g(2) * t789;
t877 = t683 * t692 - t684 * t693;
t725 = t882 * pkin(1);
t716 = -pkin(1) * t720 + pkin(5);
t712 = t833 - t938;
t704 = pkin(1) * t712 * t1006 + t705 * t716;
t702 = -pkin(1) * t987 + 0.2e1 * t712 * t716;
t698 = (t701 * t997 - t789 * t703 / 0.2e1) * t986;
t697 = (t701 * t998 + t703 * t997) * t986;
t696 = (-t787 * t702 / 0.2e1 + t704 * t999) * t985;
t695 = (t702 * t999 + t787 * t704 / 0.2e1) * t985;
t1 = [0, 0, 0, 0, 0, 0, t882, -g(1) * t789 + g(2) * t791, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(1) * t693 - g(2) * t692, g(1) * t692 - g(2) * t693, -g(3), t725, 0, 0, 0, 0, 0, 0, -g(1) * t872 + g(2) * t1009, -g(1) * t1009 - g(2) * t872, -g(3), -g(1) * t965 - g(2) * t902, 0, 0, 0, 0, 0, 0, -g(1) * t874 - g(2) * t1010, g(1) * t1010 - g(2) * t874, -g(3), -g(1) * t966 - g(2) * t903, 0, 0, 0, 0, 0, 0, -g(1) * t696 - g(2) * t695, g(1) * t695 - g(2) * t696, -g(3), (-g(1) * t787 - g(2) * t786) * pkin(5), 0, 0, 0, 0, 0, 0, g(1) * t899 + g(2) * t877, -g(1) * t877 + g(2) * t899, -g(3), t725, 0, 0, 0, 0, 0, 0, -g(1) * t995 + g(2) * t790, -g(1) * t790 - g(2) * t995, -g(3), -g(1) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t698 - g(2) * t697, g(1) * t697 - g(2) * t698, -g(3), t725, 0, 0, 0, 0, 0, 0, -g(2) * (t872 * t886 - t1016) - g(1) * (t1009 * t886 + t1011), -g(2) * (-t1009 * t885 + t1011) - g(1) * (t872 * t885 + t1016), -g(3), -g(2) * (-pkin(6) * t1009 + t902) - g(1) * (pkin(6) * t872 + t965), 0, 0, 0, 0, 0, 0, -g(1) * t1015 - g(2) * t1014, g(1) * t1014 - g(2) * t1015, -g(3), -g(2) * (pkin(4) * t1010 + t903) - g(1) * (pkin(4) * t874 + t966);];
U_reg = t1;
