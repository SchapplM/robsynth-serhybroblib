% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = palh1m1OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynB_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:30:05
% EndTime: 2020-04-15 19:30:35
% DurationCPUTime: 9.54s
% Computational Cost: add. (79038->616), mult. (173489->789), div. (0->0), fcn. (133204->22), ass. (0->237)
t1086 = sin(pkin(19));
t1087 = cos(pkin(19));
t1088 = cos(qJ(10));
t1132 = sin(qJ(10));
t1048 = t1086 * t1088 - t1087 * t1132;
t1049 = -t1086 * t1132 - t1087 * t1088;
t1091 = sin(qJ(7));
t1096 = sin(qJ(2));
t1100 = cos(qJ(7));
t1105 = cos(qJ(2));
t1042 = (-t1091 * t1105 - t1096 * t1100) * qJD(1);
t1045 = (-t1091 * t1096 + t1100 * t1105) * qJD(1);
t973 = t1048 * t1042 + t1049 * t1045;
t1123 = qJD(1) * t1105;
t1121 = qJD(2) * t1123;
t1059 = -t1096 * qJDD(1) - t1121;
t1125 = qJD(1) * t1096;
t1061 = -qJD(2) * t1125 + t1105 * qJDD(1);
t980 = -t1045 * qJD(7) + t1100 * t1059 - t1091 * t1061;
t983 = t1042 * qJD(7) + t1091 * t1059 + t1100 * t1061;
t908 = -t973 * qJD(10) - t1048 * t983 + t1049 * t980;
t972 = t1049 * t1042 - t1048 * t1045;
t909 = t972 * qJD(10) + t1048 * t980 + t1049 * t983;
t1097 = sin(qJ(1));
t1106 = cos(qJ(1));
t1066 = t1097 * g(1) - t1106 * g(2);
t1052 = -qJDD(1) * pkin(15) - t1066;
t1015 = t1052 + (-t1059 + t1121) * pkin(1);
t1084 = qJD(2) + qJD(7);
t911 = (-t1086 * t983 - t1087 * t980 + (-t1042 * t1086 + t1045 * t1087) * t1084) * pkin(4) + t1015;
t1074 = qJD(10) + t1084;
t964 = -t1074 * mrSges(11,2) + t972 * mrSges(11,3);
t965 = t1074 * mrSges(11,1) - t973 * mrSges(11,3);
t883 = m(11) * t911 - t908 * mrSges(11,1) + t909 * mrSges(11,2) - t972 * t964 + t973 * t965;
t1133 = pkin(4) * t883;
t1067 = -t1106 * g(1) - t1097 * g(2);
t1107 = qJD(1) ^ 2;
t1092 = sin(qJ(6));
t1101 = cos(qJ(6));
t1055 = t1107 * pkin(14) + t1067;
t1030 = -t1101 * g(3) - t1092 * t1055;
t1056 = (-mrSges(7,1) * t1101 + mrSges(7,2) * t1092) * qJD(1);
t1124 = qJD(1) * t1101;
t1058 = qJD(6) * t1124 + t1092 * qJDD(1);
t1064 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t1124;
t1126 = qJD(1) * t1092;
t966 = m(7) * t1030 + qJDD(6) * mrSges(7,1) - t1058 * mrSges(7,3) + qJD(6) * t1064 - t1056 * t1126;
t1032 = -t1092 * g(3) + t1101 * t1055;
t1060 = -qJD(6) * t1126 + t1101 * qJDD(1);
t1062 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t1126;
t967 = m(7) * t1032 - qJDD(6) * mrSges(7,2) + t1060 * mrSges(7,3) - qJD(6) * t1062 + t1056 * t1124;
t1120 = -t1092 * t966 + t1101 * t967;
t1054 = -t1107 * pkin(15) + t1067;
t1031 = -t1105 * g(3) - t1096 * t1054;
t1057 = (mrSges(3,1) * t1096 + mrSges(3,2) * t1105) * qJD(1);
t1065 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1123;
t1090 = sin(qJ(8));
t1095 = sin(qJ(3));
t1099 = cos(qJ(8));
t1104 = cos(qJ(3));
t1043 = (t1095 * t1105 + t1096 * t1104) * qJD(1);
t1046 = (-t1095 * t1096 + t1104 * t1105) * qJD(1);
t1012 = -t1046 * mrSges(4,1) + t1043 * mrSges(4,2);
t1085 = qJD(2) + qJD(3);
t1026 = t1085 * mrSges(4,1) - t1043 * mrSges(4,3);
t1082 = qJDD(2) + qJDD(3);
t1094 = sin(qJ(4));
t1103 = cos(qJ(4));
t1006 = -t1094 * t1043 + t1103 * t1046;
t1073 = qJDD(4) + t1082;
t1076 = qJD(4) + t1085;
t1093 = sin(qJ(5));
t1102 = cos(qJ(5));
t1001 = qJD(5) - t1006;
t1008 = t1103 * t1043 + t1094 * t1046;
t981 = t1046 * qJD(3) - t1104 * t1059 + t1095 * t1061;
t984 = -t1043 * qJD(3) + t1095 * t1059 + t1104 * t1061;
t923 = -t1008 * qJD(4) - t1094 * t981 + t1103 * t984;
t925 = t1006 * qJD(4) + t1094 * t984 + t1103 * t981;
t940 = t1015 + (t1043 * t1085 - t984) * pkin(5);
t896 = (-t1006 * t1076 - t925) * pkin(11) + (t1008 * t1076 - t923) * pkin(9) + t940;
t1071 = t1076 ^ 2;
t1033 = t1096 * g(3) - t1105 * t1054;
t1019 = (-t1096 * t1105 * t1107 + qJDD(2)) * pkin(1) + t1033;
t1020 = (-t1096 ^ 2 * t1107 - qJD(2) ^ 2) * pkin(1) + t1031;
t971 = t1095 * t1019 + t1104 * t1020;
t935 = (t1043 * t1046 + t1082) * pkin(5) + t971;
t969 = -t1104 * t1019 + t1095 * t1020;
t944 = (-t1046 ^ 2 - t1085 ^ 2) * pkin(5) + t969;
t905 = t1094 * t935 + t1103 * t944;
t958 = -t1006 * pkin(9) - t1008 * pkin(11);
t898 = -t1071 * pkin(9) + t1073 * pkin(11) + t1006 * t958 + t905;
t884 = -t1093 * t898 + t1102 * t896;
t974 = -t1093 * t1008 + t1102 * t1076;
t903 = t974 * qJD(5) + t1093 * t1073 + t1102 * t925;
t920 = qJDD(5) - t923;
t975 = t1102 * t1008 + t1093 * t1076;
t938 = -t974 * mrSges(6,1) + t975 * mrSges(6,2);
t941 = -t1001 * mrSges(6,2) + t974 * mrSges(6,3);
t877 = m(6) * t884 + t920 * mrSges(6,1) - t903 * mrSges(6,3) + t1001 * t941 - t975 * t938;
t885 = t1093 * t896 + t1102 * t898;
t902 = -t975 * qJD(5) + t1102 * t1073 - t1093 * t925;
t942 = t1001 * mrSges(6,1) - t975 * mrSges(6,3);
t878 = m(6) * t885 - t920 * mrSges(6,2) + t902 * mrSges(6,3) - t1001 * t942 + t974 * t938;
t1119 = -t1093 * t877 + t1102 * t878;
t957 = -t1006 * mrSges(5,1) + t1008 * mrSges(5,2);
t990 = t1076 * mrSges(5,1) - t1008 * mrSges(5,3);
t862 = m(5) * t905 - t1073 * mrSges(5,2) + t923 * mrSges(5,3) + t1006 * t957 - t1076 * t990 + t1119;
t904 = -t1094 * t944 + t1103 * t935;
t897 = -t1073 * pkin(9) - t1071 * pkin(11) + t1008 * t958 - t904;
t1111 = -m(6) * t897 + t902 * mrSges(6,1) - t903 * mrSges(6,2) + t974 * t941 - t975 * t942;
t988 = -t1076 * mrSges(5,2) + t1006 * mrSges(5,3);
t873 = m(5) * t904 + t1073 * mrSges(5,1) - t925 * mrSges(5,3) - t1008 * t957 + t1076 * t988 + t1111;
t857 = m(4) * t969 - t1082 * mrSges(4,2) + t984 * mrSges(4,3) + t1046 * t1012 - t1085 * t1026 - t1094 * t873 + t1103 * t862;
t1029 = -t1085 * mrSges(4,2) + t1046 * mrSges(4,3);
t1130 = t1094 * t862 + t1103 * t873;
t858 = m(4) * t971 + t1082 * mrSges(4,1) - t981 * mrSges(4,3) - t1043 * t1012 + t1085 * t1029 + t1130;
t1014 = -t1042 * mrSges(8,1) + t1045 * mrSges(8,2);
t1025 = -t1084 * mrSges(8,2) + t1042 * mrSges(8,3);
t1081 = qJDD(2) + qJDD(7);
t1070 = qJDD(10) + t1081;
t1079 = t1084 ^ 2;
t968 = t1100 * t1019 - t1091 * t1020;
t991 = (-t1042 * t1087 - t1045 * t1086) * pkin(4);
t926 = -t1045 * t991 + (t1079 * t1086 + t1081 * t1087) * pkin(4) + t968;
t970 = t1091 * t1019 + t1100 * t1020;
t927 = t1042 * t991 + (-t1079 * t1087 + t1081 * t1086) * pkin(4) + t970;
t899 = -t1048 * t927 + t1049 * t926;
t933 = -t972 * mrSges(11,1) + t973 * mrSges(11,2);
t890 = m(11) * t899 + t1070 * mrSges(11,1) - t909 * mrSges(11,3) + t1074 * t964 - t973 * t933;
t900 = t1048 * t926 + t1049 * t927;
t891 = m(11) * t900 - t1070 * mrSges(11,2) + t908 * mrSges(11,3) - t1074 * t965 + t972 * t933;
t1129 = t1048 * t891 + t1049 * t890;
t870 = m(8) * t968 + t1081 * mrSges(8,1) - t983 * mrSges(8,3) - t1045 * t1014 + t1084 * t1025 + t1129;
t1028 = t1084 * mrSges(8,1) - t1045 * mrSges(8,3);
t1128 = -t1048 * t890 + t1049 * t891;
t871 = m(8) * t970 - t1081 * mrSges(8,2) + t980 * mrSges(8,3) + t1042 * t1014 - t1084 * t1028 + t1128;
t1041 = (-t1090 * t1105 - t1096 * t1099) * qJD(1);
t1044 = (-t1090 * t1096 + t1099 * t1105) * qJD(1);
t1013 = -t1041 * mrSges(9,1) + t1044 * mrSges(9,2);
t1083 = qJD(2) + qJD(8);
t1024 = -t1083 * mrSges(9,2) + t1041 * mrSges(9,3);
t1080 = qJDD(2) + qJDD(8);
t1089 = sin(qJ(9));
t1098 = cos(qJ(9));
t1007 = -t1089 * t1041 - t1098 * t1044;
t1072 = qJDD(9) + t1080;
t1075 = qJD(9) + t1083;
t985 = -t1090 * t1031 + t1099 * t1033;
t951 = (t1041 * t1044 + t1080) * pkin(2) + t985;
t986 = t1099 * t1031 + t1090 * t1033;
t963 = (-t1041 ^ 2 - t1083 ^ 2) * pkin(2) + t986;
t912 = t1089 * t963 - t1098 * t951;
t1005 = -t1098 * t1041 + t1089 * t1044;
t979 = -t1044 * qJD(8) + t1099 * t1059 - t1090 * t1061;
t982 = t1041 * qJD(8) + t1090 * t1059 + t1099 * t1061;
t924 = t1005 * qJD(9) - t1089 * t979 - t1098 * t982;
t956 = -t1005 * mrSges(10,1) + t1007 * mrSges(10,2);
t987 = -t1075 * mrSges(10,2) + t1005 * mrSges(10,3);
t894 = m(10) * t912 + t1072 * mrSges(10,1) - t924 * mrSges(10,3) - t1007 * t956 + t1075 * t987;
t913 = -t1089 * t951 - t1098 * t963;
t922 = -t1007 * qJD(9) + t1089 * t982 - t1098 * t979;
t989 = t1075 * mrSges(10,1) - t1007 * mrSges(10,3);
t895 = m(10) * t913 - t1072 * mrSges(10,2) + t922 * mrSges(10,3) + t1005 * t956 - t1075 * t989;
t1116 = -t1089 * t895 - t1098 * t894;
t881 = m(9) * t985 + t1080 * mrSges(9,1) - t982 * mrSges(9,3) - t1044 * t1013 + t1083 * t1024 + t1116;
t1027 = t1083 * mrSges(9,1) - t1044 * mrSges(9,3);
t882 = m(9) * t986 - t1080 * mrSges(9,2) + t979 * mrSges(9,3) + t1041 * t1013 - t1083 * t1027 + t1089 * t894 - t1098 * t895;
t849 = m(3) * t1031 - qJDD(2) * mrSges(3,2) + t1059 * mrSges(3,3) - qJD(2) * t1065 - t1057 * t1125 - t1090 * t881 - t1091 * t870 + t1095 * t857 + t1099 * t882 + t1100 * t871 + t1104 * t858;
t1063 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t1125;
t1114 = -t1091 * t871 - t1095 * t858 - t1100 * t870 + t1104 * t857;
t850 = m(3) * t1033 + qJDD(2) * mrSges(3,1) - t1061 * mrSges(3,3) + qJD(2) * t1063 - t1057 * t1123 + t1090 * t882 + t1099 * t881 - t1114;
t846 = m(2) * t1067 - t1107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t1096 * t849 - t1105 * t850 + t1120;
t863 = t1093 * t878 + t1102 * t877;
t1113 = m(5) * t940 - t923 * mrSges(5,1) + t925 * mrSges(5,2) - t1006 * t988 + t1008 * t990 + t863;
t1109 = -t1113 + (-m(4) - m(8)) * t1015 - t1045 * t1028 - t1043 * t1026 - t981 * mrSges(4,2) - t983 * mrSges(8,2) + t984 * mrSges(4,1) + t1046 * t1029 + t1042 * t1025 + t980 * mrSges(8,1) - t883;
t955 = (t1044 * t1083 - t979) * pkin(2) + t1052;
t1110 = m(10) * t955 - t922 * mrSges(10,1) + t924 * mrSges(10,2) - t1005 * t987 + t1007 * t989;
t1108 = t1059 * mrSges(3,1) + (-m(3) - m(9)) * t1052 - t1061 * mrSges(3,2) - t1044 * t1027 - t982 * mrSges(9,2) + t1109 + t979 * mrSges(9,1) + t1041 * t1024 - t1110;
t1053 = qJDD(1) * pkin(14) - t1066;
t1112 = -m(7) * t1053 + t1060 * mrSges(7,1) - t1058 * mrSges(7,2) + t1064 * t1124;
t1115 = -t1096 * t1063 - t1105 * t1065;
t1122 = t1092 * t1062;
t855 = (t1115 - t1122) * qJD(1) - t1107 * mrSges(2,2) + m(2) * t1066 + qJDD(1) * mrSges(2,1) + t1108 + t1112;
t1131 = t1097 * t846 + t1106 * t855;
t1127 = t1092 * t967 + t1101 * t966;
t1118 = -t1096 * t850 + t1105 * t849;
t1117 = -t1097 * t855 + t1106 * t846;
t1040 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t1105 - Ifges(3,4) * t1096) * qJD(1);
t1039 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t1092 + Ifges(7,4) * t1101) * qJD(1);
t1038 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t1105 - Ifges(3,2) * t1096) * qJD(1);
t1037 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t1092 + Ifges(7,2) * t1101) * qJD(1);
t1036 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t1105 - Ifges(3,6) * t1096) * qJD(1);
t1035 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t1092 + Ifges(7,6) * t1101) * qJD(1);
t1000 = Ifges(4,1) * t1043 + Ifges(4,4) * t1046 + Ifges(4,5) * t1085;
t999 = Ifges(8,1) * t1045 + Ifges(8,4) * t1042 + Ifges(8,5) * t1084;
t998 = Ifges(9,1) * t1044 + Ifges(9,4) * t1041 + Ifges(9,5) * t1083;
t997 = Ifges(4,4) * t1043 + Ifges(4,2) * t1046 + Ifges(4,6) * t1085;
t996 = Ifges(8,4) * t1045 + Ifges(8,2) * t1042 + Ifges(8,6) * t1084;
t995 = Ifges(9,4) * t1044 + Ifges(9,2) * t1041 + Ifges(9,6) * t1083;
t994 = Ifges(4,5) * t1043 + Ifges(4,6) * t1046 + Ifges(4,3) * t1085;
t993 = Ifges(8,5) * t1045 + Ifges(8,6) * t1042 + Ifges(8,3) * t1084;
t992 = Ifges(9,5) * t1044 + Ifges(9,6) * t1041 + Ifges(9,3) * t1083;
t962 = -qJD(1) * t1122 + t1112;
t950 = Ifges(5,1) * t1008 + Ifges(5,4) * t1006 + Ifges(5,5) * t1076;
t949 = Ifges(10,1) * t1007 + Ifges(10,4) * t1005 + Ifges(10,5) * t1075;
t948 = Ifges(5,4) * t1008 + Ifges(5,2) * t1006 + Ifges(5,6) * t1076;
t947 = Ifges(10,4) * t1007 + Ifges(10,2) * t1005 + Ifges(10,6) * t1075;
t946 = Ifges(5,5) * t1008 + Ifges(5,6) * t1006 + Ifges(5,3) * t1076;
t945 = Ifges(10,5) * t1007 + Ifges(10,6) * t1005 + Ifges(10,3) * t1075;
t937 = mrSges(7,2) * t1053 - mrSges(7,3) * t1030 + Ifges(7,1) * t1058 + Ifges(7,4) * t1060 + Ifges(7,5) * qJDD(6) - qJD(6) * t1037 + t1035 * t1124;
t936 = -mrSges(7,1) * t1053 + mrSges(7,3) * t1032 + Ifges(7,4) * t1058 + Ifges(7,2) * t1060 + Ifges(7,6) * qJDD(6) + qJD(6) * t1039 - t1035 * t1126;
t930 = Ifges(11,1) * t973 + Ifges(11,4) * t972 + Ifges(11,5) * t1074;
t929 = Ifges(11,4) * t973 + Ifges(11,2) * t972 + Ifges(11,6) * t1074;
t928 = Ifges(11,5) * t973 + Ifges(11,6) * t972 + Ifges(11,3) * t1074;
t916 = Ifges(6,1) * t975 + Ifges(6,4) * t974 + Ifges(6,5) * t1001;
t915 = Ifges(6,4) * t975 + Ifges(6,2) * t974 + Ifges(6,6) * t1001;
t914 = Ifges(6,5) * t975 + Ifges(6,6) * t974 + Ifges(6,3) * t1001;
t893 = mrSges(10,2) * t955 - mrSges(10,3) * t912 + Ifges(10,1) * t924 + Ifges(10,4) * t922 + Ifges(10,5) * t1072 + t1005 * t945 - t1075 * t947;
t892 = -mrSges(10,1) * t955 + mrSges(10,3) * t913 + Ifges(10,4) * t924 + Ifges(10,2) * t922 + Ifges(10,6) * t1072 - t1007 * t945 + t1075 * t949;
t880 = mrSges(11,2) * t911 - mrSges(11,3) * t899 + Ifges(11,1) * t909 + Ifges(11,4) * t908 + Ifges(11,5) * t1070 - t1074 * t929 + t972 * t928;
t879 = -mrSges(11,1) * t911 + mrSges(11,3) * t900 + Ifges(11,4) * t909 + Ifges(11,2) * t908 + Ifges(11,6) * t1070 + t1074 * t930 - t973 * t928;
t867 = mrSges(9,2) * t1052 - mrSges(9,3) * t985 + Ifges(9,1) * t982 + Ifges(9,4) * t979 + Ifges(9,5) * t1080 + t1041 * t992 - t1083 * t995 + t1089 * t892 - t1098 * t893;
t866 = mrSges(6,2) * t897 - mrSges(6,3) * t884 + Ifges(6,1) * t903 + Ifges(6,4) * t902 + Ifges(6,5) * t920 - t1001 * t915 + t974 * t914;
t865 = -mrSges(6,1) * t897 + mrSges(6,3) * t885 + Ifges(6,4) * t903 + Ifges(6,2) * t902 + Ifges(6,6) * t920 + t1001 * t916 - t975 * t914;
t864 = -pkin(2) * t1110 - mrSges(9,1) * t1052 + mrSges(9,3) * t986 + Ifges(9,4) * t982 + Ifges(9,2) * t979 + Ifges(9,6) * t1080 - t1044 * t992 + t1083 * t998 - t1089 * t893 - t1098 * t892;
t860 = mrSges(8,2) * t1015 - mrSges(8,3) * t968 + Ifges(8,1) * t983 + Ifges(8,4) * t980 + Ifges(8,5) * t1081 + t1042 * t993 - t1048 * t879 + t1049 * t880 - t1084 * t996 - t1086 * t1133;
t859 = -mrSges(8,1) * t1015 + mrSges(8,3) * t970 + Ifges(8,4) * t983 + Ifges(8,2) * t980 + Ifges(8,6) * t1081 - t1045 * t993 + t1048 * t880 + t1049 * t879 + t1084 * t999 - t1087 * t1133;
t853 = -pkin(9) * t863 - mrSges(5,1) * t940 - mrSges(6,1) * t884 + mrSges(6,2) * t885 + mrSges(5,3) * t905 + Ifges(5,4) * t925 - Ifges(6,5) * t903 + Ifges(5,2) * t923 + Ifges(5,6) * t1073 - Ifges(6,6) * t902 - Ifges(6,3) * t920 - t1008 * t946 + t1076 * t950 - t975 * t915 + t974 * t916;
t852 = -pkin(11) * t863 + mrSges(5,2) * t940 - mrSges(5,3) * t904 + Ifges(5,1) * t925 + Ifges(5,4) * t923 + Ifges(5,5) * t1073 + t1006 * t946 - t1076 * t948 - t1093 * t865 + t1102 * t866;
t851 = mrSges(4,2) * t1015 - mrSges(4,3) * t971 + Ifges(4,1) * t981 + Ifges(4,4) * t984 + Ifges(4,5) * t1082 + t1046 * t994 - t1085 * t997 - t1094 * t853 + t1103 * t852;
t847 = -pkin(5) * t1113 - mrSges(4,1) * t1015 + mrSges(4,3) * t969 + Ifges(4,4) * t981 + Ifges(4,2) * t984 + Ifges(4,6) * t1082 + t1085 * t1000 - t1043 * t994 + t1094 * t852 + t1103 * t853;
t843 = mrSges(3,2) * t1052 - mrSges(3,3) * t1033 + Ifges(3,1) * t1061 + Ifges(3,4) * t1059 + Ifges(3,5) * qJDD(2) - qJD(2) * t1038 - t1036 * t1125 - t1090 * t864 - t1091 * t859 + t1095 * t851 + t1099 * t867 + t1100 * t860 + t1104 * t847;
t842 = t1109 * pkin(1) - mrSges(3,1) * t1052 + mrSges(3,3) * t1031 + Ifges(3,4) * t1061 + Ifges(3,2) * t1059 + Ifges(3,6) * qJDD(2) + qJD(2) * t1040 - t1036 * t1123 + t1090 * t867 + t1091 * t860 + t1095 * t847 + t1099 * t864 + t1100 * t859 - t1104 * t851;
t841 = (-t1092 * t1037 - t1105 * t1038 + t1101 * t1039 - t1096 * t1040) * qJD(1) + mrSges(2,1) * g(3) + pkin(14) * t1127 + (-t1086 * t1128 - t1087 * t1129) * pkin(4) - pkin(5) * t1130 + t1107 * Ifges(2,5) - t1102 * t865 - t1093 * t866 - Ifges(4,3) * t1082 - Ifges(5,3) * t1073 - Ifges(9,3) * t1080 - Ifges(8,3) * t1081 + mrSges(2,3) * t1067 - Ifges(11,3) * t1070 - Ifges(10,3) * t1072 - Ifges(7,5) * t1058 - Ifges(3,6) * t1059 - Ifges(7,6) * t1060 - Ifges(3,5) * t1061 - t1045 * t996 + t1046 * t1000 + t1041 * t998 + t1042 * t999 - t1043 * t997 - t1044 * t995 - mrSges(7,1) * t1030 + mrSges(3,2) * t1031 + mrSges(7,2) * t1032 - mrSges(3,1) * t1033 - pkin(15) * t1118 - pkin(11) * t1119 - pkin(16) * t1120 + t1006 * t950 - t1007 * t947 - t1008 * t948 - pkin(2) * t1116 + t1005 * t949 - Ifges(4,6) * t984 - mrSges(9,1) * t985 + mrSges(9,2) * t986 + t1114 * pkin(1) - Ifges(9,6) * t979 - Ifges(8,6) * t980 - Ifges(4,5) * t981 - Ifges(9,5) * t982 - Ifges(8,5) * t983 - t973 * t929 - mrSges(8,1) * t968 + mrSges(4,2) * t969 + mrSges(8,2) * t970 - mrSges(4,1) * t971 + t972 * t930 - pkin(9) * t1111 - Ifges(5,5) * t925 - Ifges(10,6) * t922 - Ifges(5,6) * t923 - Ifges(10,5) * t924 - mrSges(10,1) * t912 + mrSges(10,2) * t913 - Ifges(11,6) * t908 - Ifges(11,5) * t909 + Ifges(2,6) * qJDD(1) + mrSges(11,2) * t900 - mrSges(5,1) * t904 + mrSges(5,2) * t905 - mrSges(11,1) * t899 - Ifges(3,3) * qJDD(2) - Ifges(7,3) * qJDD(6);
t840 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1066 + pkin(16) * t962 + Ifges(2,5) * qJDD(1) - t1107 * Ifges(2,6) - t1092 * t936 - t1096 * t843 + t1101 * t937 - t1105 * t842;
t1 = [-m(1) * g(1) + t1117; -m(1) * g(2) + t1131; (-m(1) - m(2)) * g(3) + t1118 + t1127; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(13) * t1131 - t1097 * t841 + t1106 * t840; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(13) * t1117 + t1097 * t840 + t1106 * t841; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t1066 - mrSges(2,2) * t1067 + t1105 * t843 - t1096 * t842 + pkin(15) * (t1115 * qJD(1) + t1108) + t1092 * t937 + t1101 * t936 - pkin(14) * t962;];
tauB = t1;
