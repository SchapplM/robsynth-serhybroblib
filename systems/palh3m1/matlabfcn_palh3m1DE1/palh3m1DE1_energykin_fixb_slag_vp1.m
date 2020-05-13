% Calculate kinetic energy for
% palh3m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1DE1_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:28:08
% EndTime: 2020-04-18 10:30:05
% DurationCPUTime: 117.32s
% Computational Cost: add. (2821709->547), mult. (4301101->943), div. (192288->22), fcn. (2706320->36), ass. (0->354)
t1075 = -2 * pkin(1);
t946 = sin(qJ(1));
t1074 = t946 ^ 2;
t952 = cos(qJ(1));
t1073 = t952 ^ 2;
t1072 = -pkin(6) - pkin(2);
t1071 = -pkin(6) + pkin(2);
t1070 = -pkin(8) - pkin(10);
t1069 = -pkin(8) + pkin(10);
t945 = sin(qJ(2));
t947 = sin(pkin(16));
t951 = cos(qJ(2));
t953 = cos(pkin(16));
t917 = t945 * t947 - t951 * t953;
t1053 = pkin(5) * t917;
t967 = pkin(1) ^ 2;
t1027 = t1053 * t1075 + t967;
t962 = pkin(5) ^ 2;
t906 = t962 + t1027;
t902 = 0.1e1 / t906;
t966 = 0.1e1 / pkin(2);
t1031 = t902 * t966;
t1024 = pkin(2) ^ 2 - pkin(6) ^ 2;
t894 = t906 + t1024;
t911 = pkin(1) - t1053;
t919 = t945 * t953 + t947 * t951;
t890 = (pkin(5) - t1072) * (pkin(5) + t1072) + t1027;
t891 = (pkin(5) - t1071) * (pkin(5) + t1071) + t1027;
t969 = sqrt(-t891 * t890);
t869 = pkin(5) * t894 * t919 + t911 * t969;
t944 = sin(qJ(3));
t1037 = t869 * t944;
t1029 = t919 * t969;
t868 = -pkin(5) * t1029 + t894 * t911;
t950 = cos(qJ(3));
t1038 = t868 * t950;
t979 = -t1037 / 0.2e1 + t1038 / 0.2e1;
t861 = t979 * t1031;
t1036 = t869 * t950;
t1039 = t868 * t944;
t978 = t1036 / 0.2e1 + t1039 / 0.2e1;
t862 = t978 * t1031;
t936 = pkin(18) + pkin(19);
t932 = sin(t936);
t933 = cos(t936);
t843 = t861 * t933 - t862 * t932;
t1056 = pkin(4) * t843;
t963 = pkin(4) ^ 2;
t1028 = -0.2e1 * pkin(3) * t1056 + t963;
t964 = pkin(3) ^ 2;
t838 = t964 + t1028;
t836 = 0.1e1 / t838;
t1068 = t836 / 0.2e1;
t914 = t919 * qJD(2);
t1019 = pkin(1) * pkin(5) * t914;
t1035 = 0.2e1 * (t890 + t891) * t1019 / t969;
t1010 = -t1035 / 0.2e1;
t915 = t917 * qJD(2);
t974 = t1010 * t919 + t915 * t969;
t853 = ((t1075 * t911 - t894) * t914 + t974) * pkin(5);
t1067 = -t853 / 0.2e1;
t1017 = -0.2e1 * t914 * t919;
t1030 = t914 * t969;
t854 = t911 * t1035 / 0.2e1 + t962 * pkin(1) * t1017 + (-t894 * t915 - t1030) * pkin(5);
t1066 = t854 / 0.2e1;
t937 = sin(pkin(17));
t1065 = t937 / 0.2e1;
t939 = sin(pkin(19));
t1064 = t939 / 0.2e1;
t1063 = t944 / 0.2e1;
t954 = cos(pkin(15));
t1062 = t954 / 0.2e1;
t1060 = pkin(1) * t951;
t1007 = 0.1e1 / t906 ^ 2 * t1019;
t810 = ((t1036 + t1039) * t1007 + (qJD(3) * t979 + t1063 * t853 + t1066 * t950) * t902) * t966;
t811 = ((t1037 - t1038) * t1007 + (qJD(3) * t978 + t1063 * t854 + t1067 * t950) * t902) * t966;
t794 = -t810 * t932 - t811 * t933;
t1059 = pkin(3) * t794;
t842 = -t861 * t932 - t862 * t933;
t1057 = pkin(4) * t842;
t1025 = pkin(8) ^ 2 - pkin(10) ^ 2;
t834 = t838 + t1025;
t839 = -pkin(3) + t1056;
t832 = (pkin(3) - t1070) * (pkin(3) + t1070) + t1028;
t833 = (pkin(3) - t1069) * (pkin(3) + t1069) + t1028;
t968 = sqrt(-t833 * t832);
t791 = t1057 * t834 - t839 * t968;
t1058 = pkin(4) * t791;
t935 = qJD(2) * t946;
t920 = qJD(3) * t946 + t935;
t1055 = pkin(4) * t920;
t1020 = -qJD(2) - qJD(3);
t921 = t1020 * t952;
t1054 = pkin(4) * t921;
t1051 = pkin(1) * qJD(2);
t1050 = Icges(3,4) * t945;
t1049 = Icges(3,4) * t951;
t961 = 0.1e1 / pkin(6);
t1032 = t902 * t961;
t893 = t906 - t1024;
t912 = pkin(1) * t917 - pkin(5);
t867 = -pkin(1) * t1029 - t893 * t912;
t870 = pkin(1) * t893 * t919 - t912 * t969;
t948 = sin(pkin(15));
t860 = (t867 * t1062 + t870 * t948 / 0.2e1) * t1032;
t863 = (t870 * t1062 - t867 * t948 / 0.2e1) * t1032;
t851 = atan2(t863, t860);
t847 = sin(t851);
t1048 = Icges(7,4) * t847;
t848 = cos(t851);
t1047 = Icges(7,4) * t848;
t1018 = pkin(4) * t1059;
t1046 = 0.2e1 * (t832 + t833) * t1018 / t968;
t1000 = t948 * t1007;
t1008 = t902 * t1062;
t1033 = t902 * t948;
t852 = ((0.2e1 * pkin(5) * t912 - t893) * t914 + t974) * pkin(1);
t855 = t912 * t1010 + t967 * pkin(5) * t1017 + (-t893 * t915 - t1030) * pkin(1);
t857 = 0.1e1 / t860 ^ 2;
t999 = t954 * t1007;
t787 = ((t855 * t1008 + t870 * t999 - t852 * t1033 / 0.2e1 - t867 * t1000) / t860 - (t852 * t1008 + t867 * t999 + t855 * t1033 / 0.2e1 + t870 * t1000) * t863 * t857) / (t857 * t863 ^ 2 + 0.1e1) * t961;
t1045 = t787 * t946;
t1044 = t787 * t952;
t1043 = t794 * t968;
t938 = cos(pkin(17));
t1042 = t836 * t938;
t957 = 0.1e1 / pkin(10);
t1041 = t836 * t957;
t1040 = t842 * t968;
t941 = cos(pkin(19));
t1034 = t902 * t941;
t1026 = (t1073 + t1074) * t951 * t1051;
t1001 = t941 * t1007;
t1002 = t939 * t1007;
t1009 = t902 * t1064;
t858 = (-t868 * t941 / 0.2e1 + t869 * t1064) * t1031;
t856 = 0.1e1 / t858 ^ 2;
t859 = (t869 * t941 / 0.2e1 + t868 * t1064) * t1031;
t786 = ((t1001 * t869 + t1002 * t868 + t1009 * t853 + t1034 * t1066) / t858 - (-t1001 * t868 + t1002 * t869 + t1009 * t854 + t1034 * t1067) * t859 * t856) / (t856 * t859 ^ 2 + 0.1e1) * t966;
t783 = t786 * t946 + t935;
t1023 = qJD(1) * t952;
t1022 = qJD(2) * t952;
t1021 = -qJD(2) - t786;
t1016 = t945 * t1051;
t835 = t838 - t1025;
t840 = -pkin(3) * t843 + pkin(4);
t790 = -pkin(3) * t1040 + t835 * t840;
t792 = pkin(3) * t835 * t842 + t840 * t968;
t776 = (-t790 * t938 / 0.2e1 + t792 * t1065) * t1041;
t777 = (t792 * t938 / 0.2e1 + t790 * t1065) * t1041;
t1014 = atan2(t777, t776);
t1012 = t836 * t1065;
t1013 = -t1046 / 0.2e1;
t793 = -t810 * t933 + t811 * t932;
t975 = t1013 * t842 - t793 * t968;
t761 = ((-0.2e1 * pkin(4) * t840 - t835) * t794 + t975) * pkin(3);
t762 = t840 * t1046 / 0.2e1 - 0.2e1 * t964 * t794 * t1057 + (t793 * t835 - t1043) * pkin(3);
t775 = 0.1e1 / t776 ^ 2;
t837 = 0.1e1 / t838 ^ 2;
t1006 = t837 * t1018;
t997 = t938 * t1006;
t998 = t937 * t1006;
t701 = ((t762 * t1042 / 0.2e1 + t792 * t997 + t761 * t1012 + t790 * t998) / t776 - (-t761 * t1042 / 0.2e1 - t790 * t997 + t762 * t1012 + t792 * t998) * t777 * t775) / (t775 * t777 ^ 2 + 0.1e1) * t957;
t698 = t701 * t946 + t920;
t959 = 0.1e1 / pkin(8);
t1011 = t959 * t1068;
t1004 = t952 * t1016;
t1003 = cos(t1014);
t996 = (-pkin(13) - t1060) * t946;
t699 = (-t701 + t1020) * t952;
t995 = rSges(3,1) * t951 - rSges(3,2) * t945;
t994 = rSges(7,1) * t848 - rSges(7,2) * t847;
t993 = Icges(3,1) * t951 - t1050;
t992 = Icges(7,1) * t848 - t1048;
t991 = -Icges(3,2) * t945 + t1049;
t990 = -Icges(7,2) * t847 + t1047;
t989 = Icges(3,5) * t951 - Icges(3,6) * t945;
t988 = Icges(7,5) * t848 - Icges(7,6) * t847;
t816 = -Icges(7,6) * t952 + t946 * t990;
t818 = -Icges(7,5) * t952 + t946 * t992;
t987 = t816 * t847 - t818 * t848;
t817 = Icges(7,6) * t946 + t952 * t990;
t819 = Icges(7,5) * t946 + t952 * t992;
t986 = -t817 * t847 + t819 * t848;
t829 = Icges(7,2) * t848 + t1048;
t830 = Icges(7,1) * t847 + t1047;
t985 = -t829 * t847 + t830 * t848;
t850 = atan2(t859, t858);
t844 = sin(t850);
t845 = cos(t850);
t827 = t844 * t951 + t845 * t945;
t826 = -t844 * t945 + t845 * t951;
t898 = -Icges(3,6) * t952 + t946 * t991;
t900 = -Icges(3,5) * t952 + t946 * t993;
t984 = t898 * t945 - t900 * t951;
t899 = Icges(3,6) * t946 + t952 * t991;
t901 = Icges(3,5) * t946 + t952 * t993;
t983 = -t899 * t945 + t901 * t951;
t923 = Icges(3,2) * t951 + t1050;
t924 = Icges(3,1) * t945 + t1049;
t982 = -t923 * t945 + t924 * t951;
t981 = t944 * t951 + t945 * t950;
t916 = t944 * t945 - t950 * t951;
t980 = -t1054 * t981 - t1004;
t908 = t916 * t946;
t910 = t916 * t952;
t977 = -t1054 * t910 + t908 * t1055 + t1026;
t931 = pkin(13) * t1023;
t976 = -t1016 * t946 + t1023 * t1060 + t931;
t973 = -pkin(4) * t908 + t996;
t789 = -pkin(4) * t1040 - t834 * t839;
t972 = atan2(t791 * t1011, t789 * t1011);
t971 = qJD(1) * t910 * pkin(4) + t1055 * t981 + t976;
t970 = sin(t972);
t949 = cos(qJ(4));
t943 = sin(qJ(4));
t942 = cos(pkin(18));
t940 = sin(pkin(18));
t927 = rSges(2,1) * t952 - rSges(2,2) * t946;
t926 = rSges(2,1) * t946 + rSges(2,2) * t952;
t925 = rSges(3,1) * t945 + rSges(3,2) * t951;
t922 = Icges(3,5) * t945 + Icges(3,6) * t951;
t909 = t981 * t952;
t907 = t981 * t946;
t905 = rSges(3,3) * t946 + t952 * t995;
t904 = -rSges(3,3) * t952 + t946 * t995;
t897 = Icges(3,3) * t946 + t952 * t989;
t896 = -Icges(3,3) * t952 + t946 * t989;
t889 = -rSges(4,1) * t981 + rSges(4,2) * t916;
t887 = -Icges(4,1) * t981 + Icges(4,4) * t916;
t886 = -Icges(4,4) * t981 + Icges(4,2) * t916;
t885 = -Icges(4,5) * t981 + Icges(4,6) * t916;
t884 = qJD(1) * t905 - t925 * t935 + t931;
t883 = -t925 * t1022 + (-pkin(13) * t946 - t904) * qJD(1);
t882 = rSges(4,1) * t910 + rSges(4,2) * t909 + rSges(4,3) * t946;
t881 = rSges(4,1) * t908 + rSges(4,2) * t907 - rSges(4,3) * t952;
t880 = Icges(4,1) * t910 + Icges(4,4) * t909 + Icges(4,5) * t946;
t879 = Icges(4,1) * t908 + Icges(4,4) * t907 - Icges(4,5) * t952;
t878 = Icges(4,4) * t910 + Icges(4,2) * t909 + Icges(4,6) * t946;
t877 = Icges(4,4) * t908 + Icges(4,2) * t907 - Icges(4,6) * t952;
t876 = Icges(4,5) * t910 + Icges(4,6) * t909 + Icges(4,3) * t946;
t875 = Icges(4,5) * t908 + Icges(4,6) * t907 - Icges(4,3) * t952;
t874 = (t904 * t946 + t905 * t952) * qJD(2);
t866 = qJD(1) * t882 - t889 * t920 + t976;
t865 = -t1004 + t889 * t921 + (-t881 + t996) * qJD(1);
t864 = t881 * t920 - t882 * t921 + t1026;
t831 = rSges(7,1) * t847 + rSges(7,2) * t848;
t828 = Icges(7,5) * t847 + Icges(7,6) * t848;
t825 = rSges(7,3) * t946 + t952 * t994;
t824 = -rSges(7,3) * t952 + t946 * t994;
t823 = t826 * t952;
t822 = t827 * t952;
t821 = t826 * t946;
t820 = t827 * t946;
t815 = Icges(7,3) * t946 + t952 * t988;
t814 = -Icges(7,3) * t952 + t946 * t988;
t809 = rSges(8,1) * t827 + rSges(8,2) * t826;
t808 = Icges(8,1) * t827 + Icges(8,4) * t826;
t807 = Icges(8,4) * t827 + Icges(8,2) * t826;
t806 = Icges(8,5) * t827 + Icges(8,6) * t826;
t805 = (-t826 * t940 + t827 * t942) * pkin(3);
t804 = rSges(8,1) * t823 - rSges(8,2) * t822 + rSges(8,3) * t946;
t803 = rSges(8,1) * t821 - rSges(8,2) * t820 - rSges(8,3) * t952;
t802 = Icges(8,1) * t823 - Icges(8,4) * t822 + Icges(8,5) * t946;
t801 = Icges(8,1) * t821 - Icges(8,4) * t820 - Icges(8,5) * t952;
t800 = Icges(8,4) * t823 - Icges(8,2) * t822 + Icges(8,6) * t946;
t799 = Icges(8,4) * t821 - Icges(8,2) * t820 - Icges(8,6) * t952;
t798 = Icges(8,5) * t823 - Icges(8,6) * t822 + Icges(8,3) * t946;
t797 = Icges(8,5) * t821 - Icges(8,6) * t820 - Icges(8,3) * t952;
t796 = (t822 * t940 + t823 * t942) * pkin(3);
t795 = (t820 * t940 + t821 * t942) * pkin(3);
t788 = 0.1e1 / t789 ^ 2;
t784 = t1021 * t952;
t781 = -t831 * t1045 + (-pkin(7) * t952 + t825) * qJD(1);
t780 = -t831 * t1044 + (pkin(7) * t946 - t824) * qJD(1);
t778 = cos(t972);
t773 = qJD(1) * t804 - t783 * t809 + t976;
t772 = -t1004 + t784 * t809 + (-t803 + t996) * qJD(1);
t771 = (t824 * t946 + t825 * t952) * t787;
t770 = t783 * t803 - t784 * t804 + t1026;
t768 = sin(t1014);
t767 = -t942 * t778 - t940 * t970;
t766 = t778 * t940 - t942 * t970;
t760 = -t1003 * t981 + t916 * t768;
t759 = -t1003 * t916 - t768 * t981;
t758 = qJD(4) * t759 + qJD(1);
t757 = t1003 * t910 + t909 * t768;
t756 = -t1003 * t909 + t768 * t910;
t755 = t1003 * t908 + t907 * t768;
t754 = -t1003 * t907 + t768 * t908;
t753 = t757 * t949 + t943 * t946;
t752 = -t757 * t943 + t946 * t949;
t751 = t755 * t949 - t943 * t952;
t750 = -t755 * t943 - t949 * t952;
t749 = t766 * t826 + t767 * t827;
t748 = -t766 * t827 + t767 * t826;
t747 = -t766 * t822 + t767 * t823;
t746 = -t766 * t823 - t767 * t822;
t745 = -t766 * t820 + t767 * t821;
t744 = -t766 * t821 - t767 * t820;
t743 = pkin(9) * t760 + pkin(11) * t759;
t742 = rSges(5,1) * t760 - rSges(5,2) * t759;
t741 = Icges(5,1) * t760 - Icges(5,4) * t759;
t740 = Icges(5,4) * t760 - Icges(5,2) * t759;
t739 = Icges(5,5) * t760 - Icges(5,6) * t759;
t738 = pkin(9) * t757 + pkin(11) * t756;
t737 = pkin(9) * t755 + pkin(11) * t754;
t736 = rSges(5,1) * t757 - rSges(5,2) * t756 + rSges(5,3) * t946;
t735 = rSges(5,1) * t755 - rSges(5,2) * t754 - rSges(5,3) * t952;
t734 = Icges(5,1) * t757 - Icges(5,4) * t756 + Icges(5,5) * t946;
t733 = Icges(5,1) * t755 - Icges(5,4) * t754 - Icges(5,5) * t952;
t732 = Icges(5,4) * t757 - Icges(5,2) * t756 + Icges(5,6) * t946;
t731 = Icges(5,4) * t755 - Icges(5,2) * t754 - Icges(5,6) * t952;
t730 = Icges(5,5) * t757 - Icges(5,6) * t756 + Icges(5,3) * t946;
t729 = Icges(5,5) * t755 - Icges(5,6) * t754 - Icges(5,3) * t952;
t728 = rSges(9,1) * t749 + rSges(9,2) * t748;
t727 = Icges(9,1) * t749 + Icges(9,4) * t748;
t726 = Icges(9,4) * t749 + Icges(9,2) * t748;
t725 = Icges(9,5) * t749 + Icges(9,6) * t748;
t724 = rSges(9,1) * t747 + rSges(9,2) * t746 + rSges(9,3) * t946;
t723 = rSges(9,1) * t745 + rSges(9,2) * t744 - rSges(9,3) * t952;
t722 = Icges(9,1) * t747 + Icges(9,4) * t746 + Icges(9,5) * t946;
t721 = Icges(9,1) * t745 + Icges(9,4) * t744 - Icges(9,5) * t952;
t720 = Icges(9,4) * t747 + Icges(9,2) * t746 + Icges(9,6) * t946;
t719 = Icges(9,4) * t745 + Icges(9,2) * t744 - Icges(9,6) * t952;
t718 = Icges(9,5) * t747 + Icges(9,6) * t746 + Icges(9,3) * t946;
t717 = Icges(9,5) * t745 + Icges(9,6) * t744 - Icges(9,3) * t952;
t716 = 0.2e1 * (((t839 * t1013 + (t793 * t834 - t1043) * pkin(4)) * t1068 + (-t836 * t842 * t963 + t1058 * t837) * t1059) / t789 - ((-t794 * t834 + t975) * t1068 + (t789 * t837 + t836 * t839) * t1059) * t788 * t1058) * pkin(8) / (t788 * t791 ^ 2 + 0.1e1) * t838 * t959;
t715 = (-t716 + t1021) * t952;
t714 = t716 * t946 + t783;
t713 = rSges(6,3) * t759 + (rSges(6,1) * t949 - rSges(6,2) * t943) * t760;
t712 = Icges(6,5) * t759 + (Icges(6,1) * t949 - Icges(6,4) * t943) * t760;
t711 = Icges(6,6) * t759 + (Icges(6,4) * t949 - Icges(6,2) * t943) * t760;
t710 = Icges(6,3) * t759 + (Icges(6,5) * t949 - Icges(6,6) * t943) * t760;
t709 = rSges(6,1) * t753 + rSges(6,2) * t752 + rSges(6,3) * t756;
t708 = rSges(6,1) * t751 + rSges(6,2) * t750 + rSges(6,3) * t754;
t707 = Icges(6,1) * t753 + Icges(6,4) * t752 + Icges(6,5) * t756;
t706 = Icges(6,1) * t751 + Icges(6,4) * t750 + Icges(6,5) * t754;
t705 = Icges(6,4) * t753 + Icges(6,2) * t752 + Icges(6,6) * t756;
t704 = Icges(6,4) * t751 + Icges(6,2) * t750 + Icges(6,6) * t754;
t703 = Icges(6,5) * t753 + Icges(6,6) * t752 + Icges(6,3) * t756;
t702 = Icges(6,5) * t751 + Icges(6,6) * t750 + Icges(6,3) * t754;
t697 = qJD(4) * t754 + t699;
t696 = qJD(4) * t756 + t698;
t695 = -t714 * t728 - t783 * t805 + (t724 + t796) * qJD(1) + t976;
t694 = -t1004 + t715 * t728 + t784 * t805 + (-t723 - t795 + t996) * qJD(1);
t693 = qJD(1) * t736 - t698 * t742 + t971;
t692 = t699 * t742 + (-t735 + t973) * qJD(1) + t980;
t691 = t714 * t723 - t715 * t724 + t783 * t795 - t784 * t796 + t1026;
t690 = t698 * t735 - t699 * t736 + t977;
t689 = qJD(1) * t738 - t696 * t713 - t698 * t743 + t709 * t758 + t971;
t688 = t697 * t713 + t699 * t743 - t708 * t758 + (-t737 + t973) * qJD(1) + t980;
t687 = t696 * t708 - t697 * t709 + t698 * t737 - t699 * t738 + t977;
t1 = -((-t952 * t922 + t946 * t982) * qJD(1) + (t1073 * t896 + (t983 * t946 + (-t897 + t984) * t952) * t946) * qJD(2)) * t1022 / 0.2e1 - ((-t952 * t828 + t946 * t985) * qJD(1) + (t1073 * t814 + (t986 * t946 + (-t815 + t987) * t952) * t946) * t787) * t1044 / 0.2e1 + ((t946 * t922 + t952 * t982) * qJD(1) + (t1074 * t897 + (t984 * t952 + (-t896 + t983) * t946) * t952) * qJD(2)) * t935 / 0.2e1 + ((t946 * t828 + t952 * t985) * qJD(1) + (t1074 * t815 + (t987 * t952 + (-t814 + t986) * t946) * t952) * t787) * t1045 / 0.2e1 + m(5) * (t690 ^ 2 + t692 ^ 2 + t693 ^ 2) / 0.2e1 + m(9) * (t691 ^ 2 + t694 ^ 2 + t695 ^ 2) / 0.2e1 + t783 * ((t946 * t798 - t822 * t800 + t823 * t802) * t783 + (t797 * t946 - t799 * t822 + t801 * t823) * t784 + (t806 * t946 - t807 * t822 + t808 * t823) * qJD(1)) / 0.2e1 + m(7) * (t771 ^ 2 + t780 ^ 2 + t781 ^ 2) / 0.2e1 + m(8) * (t770 ^ 2 + t772 ^ 2 + t773 ^ 2) / 0.2e1 + t696 * ((t703 * t756 + t705 * t752 + t707 * t753) * t696 + (t702 * t756 + t704 * t752 + t706 * t753) * t697 + (t710 * t756 + t711 * t752 + t712 * t753) * t758) / 0.2e1 + t697 * ((t703 * t754 + t705 * t750 + t707 * t751) * t696 + (t702 * t754 + t704 * t750 + t706 * t751) * t697 + (t710 * t754 + t711 * t750 + t712 * t751) * t758) / 0.2e1 + t758 * ((t703 * t696 + t702 * t697 + t710 * t758) * t759 + ((-t705 * t943 + t707 * t949) * t696 + (-t704 * t943 + t706 * t949) * t697 + (-t711 * t943 + t712 * t949) * t758) * t760) / 0.2e1 + t714 * ((t946 * t718 + t746 * t720 + t747 * t722) * t714 + (t717 * t946 + t719 * t746 + t721 * t747) * t715 + (t725 * t946 + t726 * t746 + t727 * t747) * qJD(1)) / 0.2e1 + t920 * ((t946 * t876 + t909 * t878 + t910 * t880) * t920 + (t875 * t946 + t877 * t909 + t879 * t910) * t921 + (t885 * t946 + t886 * t909 + t887 * t910) * qJD(1)) / 0.2e1 + t698 * ((t946 * t730 - t756 * t732 + t757 * t734) * t698 + (t729 * t946 - t731 * t756 + t733 * t757) * t699 + (t739 * t946 - t740 * t756 + t741 * t757) * qJD(1)) / 0.2e1 + m(6) * (t687 ^ 2 + t688 ^ 2 + t689 ^ 2) / 0.2e1 + m(4) * (t864 ^ 2 + t865 ^ 2 + t866 ^ 2) / 0.2e1 + m(3) * (t874 ^ 2 + t883 ^ 2 + t884 ^ 2) / 0.2e1 + t715 * ((-t718 * t952 + t720 * t744 + t722 * t745) * t714 + (-t952 * t717 + t744 * t719 + t745 * t721) * t715 + (-t725 * t952 + t726 * t744 + t727 * t745) * qJD(1)) / 0.2e1 + t784 * ((-t798 * t952 - t800 * t820 + t802 * t821) * t783 + (-t952 * t797 - t820 * t799 + t821 * t801) * t784 + (-t806 * t952 - t807 * t820 + t808 * t821) * qJD(1)) / 0.2e1 + t921 * ((-t876 * t952 + t878 * t907 + t880 * t908) * t920 + (-t952 * t875 + t907 * t877 + t908 * t879) * t921 + (-t885 * t952 + t886 * t907 + t887 * t908) * qJD(1)) / 0.2e1 + t699 * ((-t730 * t952 - t732 * t754 + t734 * t755) * t698 + (-t952 * t729 - t754 * t731 + t755 * t733) * t699 + (-t739 * t952 - t740 * t754 + t741 * t755) * qJD(1)) / 0.2e1 + (Icges(2,3) + m(2) * (t926 ^ 2 + t927 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t899 * t951 + t901 * t945) * t946 - (t898 * t951 + t900 * t945) * t952) * qJD(2) + ((t817 * t848 + t819 * t847) * t946 - (t816 * t848 + t818 * t847) * t952) * t787 + (t878 * t916 - t880 * t981) * t920 + (t877 * t916 - t879 * t981) * t921 + (t800 * t826 + t802 * t827) * t783 + (t799 * t826 + t801 * t827) * t784 + (-t732 * t759 + t734 * t760) * t698 + (-t731 * t759 + t733 * t760) * t699 + (t720 * t748 + t722 * t749) * t714 + (t719 * t748 + t721 * t749) * t715 + (t726 * t748 + t727 * t749 - t740 * t759 + t741 * t760 + t826 * t807 + t827 * t808 + t829 * t848 + t830 * t847 + t916 * t886 - t887 * t981 + t951 * t923 + t945 * t924) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
