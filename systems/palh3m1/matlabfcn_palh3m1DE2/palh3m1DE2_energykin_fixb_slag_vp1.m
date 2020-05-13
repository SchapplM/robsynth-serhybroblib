% Calculate kinetic energy for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1DE2_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:36:09
% EndTime: 2020-04-19 19:37:25
% DurationCPUTime: 76.60s
% Computational Cost: add. (1803237->486), mult. (2768675->856), div. (119448->22), fcn. (1734982->38), ass. (0->362)
t1088 = -2 * pkin(1);
t924 = sin(qJ(1));
t1087 = t924 ^ 2;
t930 = cos(qJ(1));
t1086 = t930 ^ 2;
t1085 = -pkin(6) - pkin(2);
t1084 = -pkin(6) + pkin(2);
t1083 = -pkin(8) - pkin(10);
t1082 = -pkin(8) + pkin(10);
t923 = sin(qJ(2));
t925 = sin(pkin(16));
t929 = cos(qJ(2));
t931 = cos(pkin(16));
t895 = t923 * t925 - t929 * t931;
t1067 = pkin(5) * t895;
t945 = pkin(1) ^ 2;
t1024 = t1067 * t1088 + t945;
t940 = pkin(5) ^ 2;
t881 = t940 + t1024;
t875 = 0.1e1 / t881;
t944 = 0.1e1 / pkin(2);
t1034 = t875 * t944;
t1021 = pkin(2) ^ 2 - pkin(6) ^ 2;
t868 = t881 + t1021;
t884 = pkin(1) - t1067;
t896 = t923 * t931 + t925 * t929;
t857 = (pkin(5) - t1085) * (pkin(5) + t1085) + t1024;
t858 = (pkin(5) - t1084) * (pkin(5) + t1084) + t1024;
t947 = sqrt(-t858 * t857);
t843 = pkin(5) * t868 * t896 + t884 * t947;
t922 = sin(qJ(3));
t1040 = t843 * t922;
t1032 = t896 * t947;
t842 = -pkin(5) * t1032 + t868 * t884;
t928 = cos(qJ(3));
t1041 = t842 * t928;
t962 = -t1040 / 0.2e1 + t1041 / 0.2e1;
t837 = t962 * t1034;
t1039 = t843 * t928;
t1042 = t842 * t922;
t961 = t1039 / 0.2e1 + t1042 / 0.2e1;
t838 = t961 * t1034;
t915 = pkin(18) + pkin(19);
t908 = sin(t915);
t909 = cos(t915);
t816 = t837 * t909 - t838 * t908;
t1069 = pkin(4) * t816;
t941 = pkin(4) ^ 2;
t1027 = -0.2e1 * pkin(3) * t1069 + t941;
t942 = pkin(3) ^ 2;
t810 = t942 + t1027;
t808 = 0.1e1 / t810;
t1081 = t808 / 0.2e1;
t887 = t896 * qJD(2);
t1016 = pkin(1) * pkin(5) * t887;
t1038 = 0.2e1 * (t857 + t858) * t1016 / t947;
t1009 = -t1038 / 0.2e1;
t888 = t895 * qJD(2);
t958 = t896 * t1009 + t888 * t947;
t829 = ((t884 * t1088 - t868) * t887 + t958) * pkin(5);
t1080 = -t829 / 0.2e1;
t1014 = -0.2e1 * t887 * t896;
t1033 = t887 * t947;
t830 = t884 * t1038 / 0.2e1 + t940 * pkin(1) * t1014 + (-t868 * t888 - t1033) * pkin(5);
t1079 = t830 / 0.2e1;
t917 = sin(pkin(17));
t1078 = t917 / 0.2e1;
t919 = sin(pkin(19));
t1077 = t919 / 0.2e1;
t1076 = t922 / 0.2e1;
t932 = cos(pkin(15));
t1075 = t932 / 0.2e1;
t1006 = 0.1e1 / t881 ^ 2 * t1016;
t774 = ((t1039 + t1042) * t1006 + (t962 * qJD(3) + t829 * t1076 + t928 * t1079) * t875) * t944;
t775 = ((t1040 - t1041) * t1006 + (t961 * qJD(3) + t830 * t1076 + t928 * t1080) * t875) * t944;
t773 = -t774 * t908 - t775 * t909;
t1073 = pkin(3) * t773;
t920 = cos(pkin(19));
t834 = (-t842 * t920 / 0.2e1 + t843 * t1077) * t1034;
t835 = (t843 * t920 / 0.2e1 + t842 * t1077) * t1034;
t821 = qJ(2) + atan2(t835, t834);
t820 = pkin(18) - t821;
t1072 = pkin(3) * sin(t820);
t815 = -t837 * t908 - t838 * t909;
t1070 = pkin(4) * t815;
t1022 = pkin(8) ^ 2 - pkin(10) ^ 2;
t806 = t810 + t1022;
t811 = -pkin(3) + t1069;
t804 = (pkin(3) - t1083) * (pkin(3) + t1083) + t1027;
t805 = (pkin(3) - t1082) * (pkin(3) + t1082) + t1027;
t946 = sqrt(-t805 * t804);
t770 = t806 * t1070 - t811 * t946;
t1071 = pkin(4) * t770;
t916 = qJ(2) + qJ(3);
t912 = sin(t916);
t1068 = pkin(4) * t912;
t1065 = pkin(13) * t924;
t1064 = t929 * pkin(1);
t1063 = Icges(3,4) * t923;
t1062 = Icges(3,4) * t929;
t1061 = Icges(4,4) * t912;
t913 = cos(t916);
t1060 = Icges(4,4) * t913;
t935 = 0.1e1 / pkin(10);
t1044 = t808 * t935;
t1043 = t815 * t946;
t807 = t810 - t1022;
t812 = -pkin(3) * t816 + pkin(4);
t769 = -pkin(3) * t1043 + t807 * t812;
t771 = pkin(3) * t807 * t815 + t812 * t946;
t918 = cos(pkin(17));
t754 = (-t769 * t918 / 0.2e1 + t771 * t1078) * t1044;
t755 = (t771 * t918 / 0.2e1 + t769 * t1078) * t1044;
t746 = atan2(t755, t754) + t916;
t744 = sin(t746);
t1059 = Icges(5,4) * t744;
t745 = cos(t746);
t1058 = Icges(5,4) * t745;
t939 = 0.1e1 / pkin(6);
t1035 = t875 * t939;
t867 = t881 - t1021;
t885 = pkin(1) * t895 - pkin(5);
t841 = -pkin(1) * t1032 - t867 * t885;
t844 = pkin(1) * t867 * t896 - t885 * t947;
t926 = sin(pkin(15));
t836 = (t841 * t1075 + t844 * t926 / 0.2e1) * t1035;
t839 = (t844 * t1075 - t841 * t926 / 0.2e1) * t1035;
t827 = atan2(t839, t836);
t823 = sin(t827);
t1057 = Icges(7,4) * t823;
t824 = cos(t827);
t1056 = Icges(7,4) * t824;
t818 = sin(t821);
t1055 = Icges(8,4) * t818;
t819 = cos(t821);
t1054 = Icges(8,4) * t819;
t937 = 0.1e1 / pkin(8);
t1010 = t937 * t1081;
t768 = -pkin(4) * t1043 - t806 * t811;
t752 = -atan2(t770 * t1010, t768 * t1010) + t820;
t750 = sin(t752);
t1053 = Icges(9,4) * t750;
t751 = cos(t752);
t1052 = Icges(9,4) * t751;
t1051 = t744 * t924;
t1050 = t744 * t930;
t1015 = pkin(4) * t1073;
t1049 = 0.2e1 * (t804 + t805) * t1015 / t946;
t1007 = t875 * t1075;
t1036 = t875 * t926;
t828 = ((0.2e1 * pkin(5) * t885 - t867) * t887 + t958) * pkin(1);
t831 = t885 * t1009 + t945 * pkin(5) * t1014 + (-t867 * t888 - t1033) * pkin(1);
t833 = 0.1e1 / t836 ^ 2;
t998 = t932 * t1006;
t999 = t926 * t1006;
t766 = ((t831 * t1007 + t844 * t998 - t828 * t1036 / 0.2e1 - t841 * t999) / t836 - (t828 * t1007 + t841 * t998 + t831 * t1036 / 0.2e1 + t844 * t999) * t839 * t833) / (t833 * t839 ^ 2 + 0.1e1) * t939;
t1048 = t766 * t924;
t1047 = t766 * t930;
t1046 = t773 * t946;
t1045 = t808 * t918;
t1037 = t875 * t920;
t921 = sin(qJ(4));
t1031 = t921 * t924;
t1030 = t921 * t930;
t927 = cos(qJ(4));
t1029 = t924 * t927;
t1028 = t927 * t930;
t1026 = pkin(3) * cos(t820);
t1020 = qJD(2) * t930;
t889 = t1064 * t924;
t890 = t1064 * t930;
t911 = qJD(2) * t924;
t1025 = t890 * t1020 + t889 * t911;
t1023 = pkin(4) * t913;
t1000 = t920 * t1006;
t1001 = t919 * t1006;
t1008 = t875 * t1077;
t832 = 0.1e1 / t834 ^ 2;
t765 = ((t843 * t1000 + t842 * t1001 + t829 * t1008 + t1037 * t1079) / t834 - (-t842 * t1000 + t843 * t1001 + t830 * t1008 + t1037 * t1080) * t835 * t832) / (t832 * t835 ^ 2 + 0.1e1) * t944;
t762 = t765 * t924 + t911;
t898 = qJD(3) * t924 + t911;
t1019 = qJD(4) * t744;
t1018 = -qJD(2) - t765;
t1017 = -qJD(2) - qJD(3);
t1013 = pkin(1) * qJD(2) * t923;
t1011 = t808 * t1078;
t1012 = -t1049 / 0.2e1;
t772 = -t774 * t909 + t775 * t908;
t959 = t815 * t1012 - t772 * t946;
t725 = ((-0.2e1 * pkin(4) * t812 - t807) * t773 + t959) * pkin(3);
t726 = t812 * t1049 / 0.2e1 - 0.2e1 * t942 * t773 * t1070 + (t772 * t807 - t1046) * pkin(3);
t753 = 0.1e1 / t754 ^ 2;
t809 = 0.1e1 / t810 ^ 2;
t1005 = t809 * t1015;
t996 = t918 * t1005;
t997 = t917 * t1005;
t692 = ((t726 * t1045 / 0.2e1 + t771 * t996 + t725 * t1011 + t769 * t997) / t754 - (-t725 * t1045 / 0.2e1 - t769 * t996 + t726 * t1011 + t771 * t997) * t755 * t753) / (t753 * t755 ^ 2 + 0.1e1) * t935;
t689 = t692 * t924 + t898;
t1004 = -t889 - t1065;
t1002 = t930 * t1013;
t855 = t1023 * t924;
t995 = t1004 + t855;
t994 = -pkin(9) * t745 - pkin(11) * t744;
t690 = (-t692 + t1017) * t930;
t993 = rSges(3,1) * t929 - rSges(3,2) * t923;
t992 = -rSges(4,1) * t913 + rSges(4,2) * t912;
t991 = -rSges(5,1) * t745 + rSges(5,2) * t744;
t990 = rSges(7,1) * t824 - rSges(7,2) * t823;
t989 = rSges(8,1) * t819 - rSges(8,2) * t818;
t988 = -rSges(9,1) * t751 - rSges(9,2) * t750;
t987 = Icges(3,1) * t929 - t1063;
t986 = -Icges(4,1) * t913 + t1061;
t985 = -Icges(5,1) * t745 + t1059;
t984 = Icges(7,1) * t824 - t1057;
t983 = Icges(8,1) * t819 - t1055;
t982 = -Icges(9,1) * t751 - t1053;
t981 = -Icges(3,2) * t923 + t1062;
t980 = Icges(4,2) * t912 - t1060;
t979 = Icges(5,2) * t744 - t1058;
t978 = -Icges(7,2) * t823 + t1056;
t977 = -Icges(8,2) * t818 + t1054;
t976 = -Icges(9,2) * t750 - t1052;
t975 = Icges(3,5) * t929 - Icges(3,6) * t923;
t974 = -Icges(4,5) * t913 + Icges(4,6) * t912;
t973 = -Icges(5,5) * t745 + Icges(5,6) * t744;
t972 = Icges(7,5) * t824 - Icges(7,6) * t823;
t971 = Icges(8,5) * t819 - Icges(8,6) * t818;
t970 = -Icges(9,5) * t751 - Icges(9,6) * t750;
t788 = -Icges(7,6) * t930 + t978 * t924;
t790 = -Icges(7,5) * t930 + t984 * t924;
t969 = t788 * t823 - t790 * t824;
t789 = Icges(7,6) * t924 + t978 * t930;
t791 = Icges(7,5) * t924 + t984 * t930;
t968 = -t789 * t823 + t791 * t824;
t799 = Icges(7,2) * t824 + t1057;
t800 = Icges(7,1) * t823 + t1056;
t967 = -t799 * t823 + t800 * t824;
t871 = -Icges(3,6) * t930 + t981 * t924;
t873 = -Icges(3,5) * t930 + t987 * t924;
t966 = t871 * t923 - t873 * t929;
t872 = Icges(3,6) * t924 + t981 * t930;
t874 = Icges(3,5) * t924 + t987 * t930;
t965 = -t872 * t923 + t874 * t929;
t901 = Icges(3,2) * t929 + t1063;
t902 = Icges(3,1) * t923 + t1062;
t964 = -t901 * t923 + t902 * t929;
t856 = t1023 * t930;
t899 = t1017 * t930;
t963 = -t898 * t855 + t856 * t899 + t1025;
t907 = qJD(1) * t930 * pkin(13);
t960 = qJD(1) * t890 - t924 * t1013 + t907;
t957 = -t899 * t1068 - t1002;
t956 = qJD(1) * (-Icges(5,5) * t744 - Icges(5,6) * t745) + t689 * (Icges(5,3) * t924 + t973 * t930) + t690 * (-Icges(5,3) * t930 + t973 * t924);
t767 = 0.1e1 / t768 ^ 2;
t697 = 0.2e1 * (((t811 * t1012 + (t772 * t806 - t1046) * pkin(4)) * t1081 + (-t808 * t815 * t941 + t809 * t1071) * t1073) / t768 - ((-t773 * t806 + t959) * t1081 + (t768 * t809 + t808 * t811) * t1073) * t767 * t1071) * pkin(8) / (t767 * t770 ^ 2 + 0.1e1) * t810 * t937;
t695 = t697 * t924 + t762;
t696 = (-t697 + t1018) * t930;
t955 = qJD(1) * (Icges(9,5) * t750 - Icges(9,6) * t751) + t695 * (Icges(9,3) * t924 + t970 * t930) + t696 * (-Icges(9,3) * t930 + t970 * t924);
t763 = t1018 * t930;
t954 = qJD(1) * (Icges(8,5) * t818 + Icges(8,6) * t819) + t762 * (Icges(8,3) * t924 + t971 * t930) + t763 * (-Icges(8,3) * t930 + t971 * t924);
t953 = qJD(1) * (-Icges(4,5) * t912 - Icges(4,6) * t913) + (-Icges(4,3) * t930 + t974 * t924) * t899 + (Icges(4,3) * t924 + t974 * t930) * t898;
t952 = -qJD(1) * t856 + t898 * t1068 + t960;
t712 = -Icges(5,6) * t930 + t979 * t924;
t713 = Icges(5,6) * t924 + t979 * t930;
t714 = -Icges(5,5) * t930 + t985 * t924;
t715 = Icges(5,5) * t924 + t985 * t930;
t721 = -Icges(5,2) * t745 - t1059;
t722 = -Icges(5,1) * t744 - t1058;
t951 = (t713 * t744 - t715 * t745) * t689 + (t712 * t744 - t714 * t745) * t690 + (t721 * t744 - t722 * t745) * qJD(1);
t729 = -Icges(9,6) * t930 + t976 * t924;
t730 = Icges(9,6) * t924 + t976 * t930;
t731 = -Icges(9,5) * t930 + t982 * t924;
t732 = Icges(9,5) * t924 + t982 * t930;
t736 = -Icges(9,2) * t751 + t1053;
t737 = Icges(9,1) * t750 - t1052;
t950 = (-t730 * t750 - t732 * t751) * t695 + (-t729 * t750 - t731 * t751) * t696 + (-t736 * t750 - t737 * t751) * qJD(1);
t780 = -Icges(8,6) * t930 + t977 * t924;
t781 = Icges(8,6) * t924 + t977 * t930;
t782 = -Icges(8,5) * t930 + t983 * t924;
t783 = Icges(8,5) * t924 + t983 * t930;
t795 = Icges(8,2) * t819 + t1055;
t796 = Icges(8,1) * t818 + t1054;
t949 = (-t781 * t818 + t783 * t819) * t762 + (-t780 * t818 + t782 * t819) * t763 + (-t795 * t818 + t796 * t819) * qJD(1);
t861 = -Icges(4,6) * t930 + t980 * t924;
t862 = Icges(4,6) * t924 + t980 * t930;
t863 = -Icges(4,5) * t930 + t986 * t924;
t864 = Icges(4,5) * t924 + t986 * t930;
t892 = -Icges(4,2) * t913 - t1061;
t893 = -Icges(4,1) * t912 - t1060;
t948 = (t862 * t912 - t864 * t913) * t898 + (t861 * t912 - t863 * t913) * t899 + (t892 * t912 - t893 * t913) * qJD(1);
t905 = rSges(2,1) * t930 - rSges(2,2) * t924;
t904 = rSges(2,1) * t924 + rSges(2,2) * t930;
t903 = rSges(3,1) * t923 + rSges(3,2) * t929;
t900 = Icges(3,5) * t923 + Icges(3,6) * t929;
t894 = -rSges(4,1) * t912 - rSges(4,2) * t913;
t878 = rSges(3,3) * t924 + t993 * t930;
t877 = -rSges(3,3) * t930 + t993 * t924;
t870 = Icges(3,3) * t924 + t975 * t930;
t869 = -Icges(3,3) * t930 + t975 * t924;
t866 = rSges(4,3) * t924 + t992 * t930;
t865 = -rSges(4,3) * t930 + t992 * t924;
t852 = qJD(1) * t878 - t903 * t911 + t907;
t851 = -t903 * t1020 + (-t877 - t1065) * qJD(1);
t850 = (t877 * t924 + t878 * t930) * qJD(2);
t847 = qJD(1) * t866 - t894 * t898 + t960;
t846 = -t1002 + t894 * t899 + (t1004 - t865) * qJD(1);
t840 = t865 * t898 - t866 * t899 + t1025;
t803 = t1026 * t930;
t802 = t1026 * t924;
t801 = rSges(7,1) * t823 + rSges(7,2) * t824;
t798 = Icges(7,5) * t823 + Icges(7,6) * t824;
t797 = rSges(8,1) * t818 + rSges(8,2) * t819;
t793 = rSges(7,3) * t924 + t990 * t930;
t792 = -rSges(7,3) * t930 + t990 * t924;
t787 = Icges(7,3) * t924 + t972 * t930;
t786 = -Icges(7,3) * t930 + t972 * t924;
t785 = rSges(8,3) * t924 + t989 * t930;
t784 = -rSges(8,3) * t930 + t989 * t924;
t760 = -t801 * t1048 + (-pkin(7) * t930 + t793) * qJD(1);
t759 = -t801 * t1047 + (pkin(7) * t924 - t792) * qJD(1);
t758 = qJD(1) * t785 - t762 * t797 + t960;
t757 = -t1002 + t763 * t797 + (t1004 - t784) * qJD(1);
t749 = (t792 * t924 + t793 * t930) * t766;
t748 = t762 * t784 - t763 * t785 + t1025;
t743 = qJD(4) * t745 + qJD(1);
t742 = -t745 * t1028 + t1031;
t741 = t745 * t1030 + t1029;
t740 = -t745 * t1029 - t1030;
t739 = t745 * t1031 - t1028;
t738 = rSges(9,1) * t750 - rSges(9,2) * t751;
t734 = rSges(9,3) * t924 + t988 * t930;
t733 = -rSges(9,3) * t930 + t988 * t924;
t724 = -pkin(9) * t744 + pkin(11) * t745;
t723 = -rSges(5,1) * t744 - rSges(5,2) * t745;
t719 = t994 * t930;
t718 = t994 * t924;
t717 = rSges(5,3) * t924 + t991 * t930;
t716 = -rSges(5,3) * t930 + t991 * t924;
t709 = rSges(6,3) * t745 + (-rSges(6,1) * t927 + rSges(6,2) * t921) * t744;
t708 = Icges(6,5) * t745 + (-Icges(6,1) * t927 + Icges(6,4) * t921) * t744;
t707 = Icges(6,6) * t745 + (-Icges(6,4) * t927 + Icges(6,2) * t921) * t744;
t706 = Icges(6,3) * t745 + (-Icges(6,5) * t927 + Icges(6,6) * t921) * t744;
t705 = rSges(6,1) * t742 + rSges(6,2) * t741 - rSges(6,3) * t1050;
t704 = rSges(6,1) * t740 + rSges(6,2) * t739 - rSges(6,3) * t1051;
t703 = Icges(6,1) * t742 + Icges(6,4) * t741 - Icges(6,5) * t1050;
t702 = Icges(6,1) * t740 + Icges(6,4) * t739 - Icges(6,5) * t1051;
t701 = Icges(6,4) * t742 + Icges(6,2) * t741 - Icges(6,6) * t1050;
t700 = Icges(6,4) * t740 + Icges(6,2) * t739 - Icges(6,6) * t1051;
t699 = Icges(6,5) * t742 + Icges(6,6) * t741 - Icges(6,3) * t1050;
t698 = Icges(6,5) * t740 + Icges(6,6) * t739 - Icges(6,3) * t1051;
t694 = t762 * t1072 - t695 * t738 + (t734 + t803) * qJD(1) + t960;
t693 = -t1002 - t763 * t1072 + t696 * t738 + (t1004 - t733 - t802) * qJD(1);
t688 = -t924 * t1019 + t690;
t687 = -t930 * t1019 + t689;
t686 = t695 * t733 - t696 * t734 + t762 * t802 - t763 * t803 + t1025;
t685 = qJD(1) * t717 - t689 * t723 + t952;
t684 = t690 * t723 + (-t716 + t995) * qJD(1) + t957;
t683 = t689 * t716 - t690 * t717 + t963;
t682 = qJD(1) * t719 - t687 * t709 - t689 * t724 + t705 * t743 + t952;
t681 = t688 * t709 + t690 * t724 - t704 * t743 + (-t718 + t995) * qJD(1) + t957;
t680 = t687 * t704 - t688 * t705 + t689 * t718 - t690 * t719 + t963;
t1 = t688 * ((-t699 * t1051 + t701 * t739 + t703 * t740) * t687 + (-t698 * t1051 + t739 * t700 + t740 * t702) * t688 + (-t706 * t1051 + t707 * t739 + t708 * t740) * t743) / 0.2e1 + t687 * ((-t699 * t1050 + t741 * t701 + t742 * t703) * t687 + (-t698 * t1050 + t700 * t741 + t702 * t742) * t688 + (-t706 * t1050 + t707 * t741 + t708 * t742) * t743) / 0.2e1 + t898 * (t953 * t924 + t948 * t930) / 0.2e1 + t899 * (t948 * t924 - t953 * t930) / 0.2e1 + t762 * (t954 * t924 + t949 * t930) / 0.2e1 + m(8) * (t748 ^ 2 + t757 ^ 2 + t758 ^ 2) / 0.2e1 + m(4) * (t840 ^ 2 + t846 ^ 2 + t847 ^ 2) / 0.2e1 + m(3) * (t850 ^ 2 + t851 ^ 2 + t852 ^ 2) / 0.2e1 + m(7) * (t749 ^ 2 + t759 ^ 2 + t760 ^ 2) / 0.2e1 + t743 * ((t687 * t699 + t688 * t698 + t706 * t743) * t745 + ((t701 * t921 - t703 * t927) * t687 + (t700 * t921 - t702 * t927) * t688 + (t707 * t921 - t708 * t927) * t743) * t744) / 0.2e1 + m(6) * (t680 ^ 2 + t681 ^ 2 + t682 ^ 2) / 0.2e1 + t763 * (t949 * t924 - t954 * t930) / 0.2e1 + t695 * (t955 * t924 + t950 * t930) / 0.2e1 + t696 * (t950 * t924 - t955 * t930) / 0.2e1 + t689 * (t956 * t924 + t951 * t930) / 0.2e1 + t690 * (t951 * t924 - t956 * t930) / 0.2e1 + m(9) * (t686 ^ 2 + t693 ^ 2 + t694 ^ 2) / 0.2e1 + m(5) * (t683 ^ 2 + t684 ^ 2 + t685 ^ 2) / 0.2e1 - ((-t930 * t900 + t964 * t924) * qJD(1) + (t1086 * t869 + (t965 * t924 + (-t870 + t966) * t930) * t924) * qJD(2)) * t1020 / 0.2e1 - ((-t930 * t798 + t967 * t924) * qJD(1) + (t1086 * t786 + (t968 * t924 + (-t787 + t969) * t930) * t924) * t766) * t1047 / 0.2e1 + ((t924 * t900 + t964 * t930) * qJD(1) + (t1087 * t870 + (t966 * t930 + (-t869 + t965) * t924) * t930) * qJD(2)) * t911 / 0.2e1 + ((t924 * t798 + t967 * t930) * qJD(1) + (t1087 * t787 + (t969 * t930 + (-t786 + t968) * t924) * t930) * t766) * t1048 / 0.2e1 + (m(2) * (t904 ^ 2 + t905 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t872 * t929 + t874 * t923) * t924 - (t871 * t929 + t873 * t923) * t930) * qJD(2) + ((t789 * t824 + t791 * t823) * t924 - (t788 * t824 + t790 * t823) * t930) * t766 + (-t862 * t913 - t864 * t912) * t898 + (-t861 * t913 - t863 * t912) * t899 + (t781 * t819 + t783 * t818) * t762 + (t780 * t819 + t782 * t818) * t763 + (-t730 * t751 + t732 * t750) * t695 + (-t729 * t751 + t731 * t750) * t696 + (-t713 * t745 - t715 * t744) * t689 + (-t712 * t745 - t714 * t744) * t690 + (-t721 * t745 - t722 * t744 - t736 * t751 + t737 * t750 + t795 * t819 + t796 * t818 + t799 * t824 + t823 * t800 - t913 * t892 - t912 * t893 + t929 * t901 + t923 * t902) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
