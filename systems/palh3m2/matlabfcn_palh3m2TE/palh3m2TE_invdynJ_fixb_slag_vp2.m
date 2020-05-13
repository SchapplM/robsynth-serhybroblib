% Calculate vector of inverse dynamics joint torques for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2TE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2TE_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:47:42
% EndTime: 2020-05-07 01:47:54
% DurationCPUTime: 11.70s
% Computational Cost: add. (10699->499), mult. (20079->730), div. (0->0), fcn. (24260->22), ass. (0->235)
t1074 = sin(qJ(3));
t1147 = qJD(2) * t1074;
t1132 = pkin(1) * t1147;
t1065 = qJD(3) + qJD(2);
t1080 = cos(qJ(3));
t1187 = pkin(1) * t1080;
t1131 = qJD(2) * t1187;
t1103 = -pkin(4) * t1065 + t1131;
t1066 = pkin(17) + pkin(18);
t1059 = sin(t1066);
t1060 = cos(t1066);
t1075 = sin(qJ(2));
t1081 = cos(qJ(2));
t1068 = sin(pkin(16));
t1069 = cos(pkin(16));
t1077 = sin(pkin(15));
t1083 = cos(pkin(15));
t1028 = t1068 * t1083 + t1069 * t1077;
t1029 = -t1068 * t1077 + t1069 * t1083;
t1106 = t1028 * t1074 - t1029 * t1080;
t1124 = -t1028 * t1080 - t1074 * t1029;
t1210 = t1075 * t1106 + t1081 * t1124;
t1225 = t1106 * t1081;
t1215 = -t1124 * t1075 + t1225;
t902 = t1210 * t1059 - t1060 * t1215;
t1236 = t1103 * t902;
t901 = t1059 * t1215 + t1060 * t1210;
t1237 = -t901 * t1132 + t1236;
t1231 = t1103 * t901;
t898 = t902 * t1132;
t882 = -t898 - t1231;
t1054 = pkin(4) - t1187;
t1188 = pkin(1) * t1074;
t1133 = t902 * t1188;
t1146 = qJD(3) * t1080;
t962 = t1106 * qJD(3);
t963 = t1124 * qJD(3);
t904 = qJD(2) * t1210 + t1075 * t962 + t1081 * t963;
t1158 = t904 * t1059;
t1234 = qJD(2) * t1215 - t963 * t1075;
t903 = -t1081 * t962 - t1234;
t888 = t1060 * t903 + t1158;
t890 = -t1059 * t903 + t1060 * t904;
t1165 = -t890 * t1054 - qJD(3) * t1133 + (t1074 * t888 - t1146 * t901) * pkin(1) + t882;
t1067 = qJ(3) + qJ(2);
t1061 = sin(t1067);
t1062 = cos(t1067);
t1233 = (t1061 * t1074 + t1062 * t1080) * mrSges(9,3);
t1229 = mrSges(4,1) + mrSges(9,1);
t1072 = mrSges(4,2) + mrSges(9,2);
t1189 = pkin(1) * qJD(2);
t1078 = sin(pkin(14));
t1084 = cos(pkin(14));
t1034 = t1077 * t1084 - t1083 * t1078;
t1035 = t1077 * t1078 + t1083 * t1084;
t978 = t1034 * t1081 + t1035 * t1075;
t1228 = Ifges(7,4) * t978;
t889 = (qJD(3) * t1225 + t1234) * t1060 - t1158;
t871 = t1054 * t889 + (-t1074 * t890 + (t1074 * t901 - t1080 * t902) * qJD(3)) * pkin(1);
t1227 = t871 - t1237;
t1134 = qJDD(2) * t1080;
t1026 = (qJD(3) * t1147 - t1134) * pkin(1);
t1027 = (-qJD(2) * t1146 - qJDD(2) * t1074) * pkin(1);
t1064 = qJDD(2) + qJDD(3);
t868 = t888 * t1132 + t890 * t1103 - t902 * (pkin(4) * t1064 + t1026) + t901 * t1027;
t885 = -t902 * t1054 - t1188 * t901;
t1224 = t1165 * t1237 + t868 * t885;
t1052 = pkin(1) * t1081 + pkin(12);
t1130 = (m(8) + m(4)) * t1052;
t1071 = cos(pkin(18));
t1190 = sin(pkin(18));
t1030 = t1071 * t1077 + t1190 * t1083;
t1031 = t1071 * t1083 - t1190 * t1077;
t948 = t1052 + (sin(pkin(17)) * t1030 - cos(pkin(17)) * t1031) * pkin(3);
t1199 = m(9) * t948;
t1223 = t1130 + t1199;
t1085 = m(5) + m(6);
t1123 = -t1034 * t1075 + t1035 * t1081;
t1220 = -t1123 / 0.2e1;
t1219 = pkin(4) * t1085;
t1218 = Ifges(7,4) * t1123;
t1217 = t948 * (-mrSges(9,1) * t1061 - mrSges(9,2) * t1062);
t959 = t1028 * t1059 - t1029 * t1060;
t954 = t959 * qJD(1);
t958 = t1028 * t1060 + t1029 * t1059;
t953 = t958 * qJD(1);
t919 = -mrSges(5,1) * t954 + mrSges(5,2) * t953;
t1032 = t1074 * t1075 - t1080 * t1081;
t1019 = t1032 * qJD(1);
t987 = -pkin(4) * t1019 - t1052 * qJD(1);
t1216 = m(5) * t987 + t919;
t952 = t959 * qJDD(1);
t1073 = sin(qJ(4));
t1079 = cos(qJ(4));
t910 = -pkin(8) * t954 - pkin(10) * t953 + t987;
t876 = -t1073 * t882 + t1079 * t910;
t877 = t1073 * t910 + t1079 * t882;
t1212 = -t1073 * t876 + t1079 * t877;
t869 = t902 * t1027 + (t1064 * t901 + t1065 * t889) * pkin(4) + (-t901 * t1134 + (-t889 * t1080 + (qJD(3) * t901 - t890) * t1074) * qJD(2)) * pkin(1);
t1138 = t1075 * qJD(1);
t1128 = qJD(2) * t1138;
t1046 = pkin(1) * t1128;
t1135 = -t1052 * qJDD(1) + t1046;
t1038 = qJDD(1) * t1081 - t1128;
t1137 = t1081 * qJD(1);
t1039 = qJD(2) * t1137 + qJDD(1) * t1075;
t1105 = t1074 * t1081 + t1075 * t1080;
t1098 = t1105 * qJD(3);
t944 = qJD(1) * t1098 - t1038 * t1080 + t1039 * t1074;
t933 = -pkin(4) * t944 + t1135;
t951 = t958 * qJDD(1);
t905 = -pkin(8) * t952 - pkin(10) * t951 + t933;
t866 = qJD(4) * t876 + t1073 * t905 + t1079 * t869;
t867 = -qJD(4) * t877 - t1073 * t869 + t1079 * t905;
t1211 = -t1073 * t867 + t1079 * t866;
t1145 = qJD(4) * t1073;
t915 = t1079 * t951 - t953 * t1145;
t1144 = qJD(4) * t1079;
t916 = -t1073 * t951 - t953 * t1144;
t947 = qJDD(4) - t952;
t1209 = t867 * mrSges(6,1) - t866 * mrSges(6,2) + Ifges(6,5) * t915 + Ifges(6,6) * t916 + Ifges(6,3) * t947;
t1121 = mrSges(6,1) * t1079 - mrSges(6,2) * t1073;
t949 = qJD(4) - t954;
t1208 = -mrSges(6,3) * t1212 + t1237 * t1121 + t949 * (-Ifges(6,5) * t1073 - Ifges(6,6) * t1079) + (Ifges(6,4) * t1073 ^ 2 - (Ifges(6,4) * t1079 + (Ifges(6,1) - Ifges(6,2)) * t1073) * t1079) * t953;
t1207 = qJD(1) ^ 2;
t1205 = pkin(6) * m(7);
t967 = t978 * qJD(2);
t1204 = -t967 / 0.2e1;
t968 = t1123 * qJD(2);
t1203 = t968 / 0.2e1;
t1097 = t1032 * qJD(3);
t981 = qJD(2) * t1032 + t1097;
t1202 = t981 / 0.2e1;
t982 = qJD(2) * t1105 + t1098;
t1201 = t982 / 0.2e1;
t1196 = t1019 / 0.2e1;
t1020 = t1105 * qJD(1);
t1195 = -t1020 / 0.2e1;
t1194 = t1075 / 0.2e1;
t1192 = t1081 / 0.2e1;
t1186 = pkin(4) * t1020;
t1185 = mrSges(6,3) * t953;
t1184 = Ifges(7,1) * t978;
t1183 = Ifges(7,2) * t1123;
t1182 = mrSges(4,3) * t1019;
t1181 = mrSges(4,3) * t1020;
t1180 = Ifges(9,1) * t1061;
t1179 = Ifges(3,4) * t1075;
t1176 = Ifges(9,4) * t1061;
t1175 = Ifges(9,4) * t1062;
t1174 = Ifges(3,5) * t1081;
t1173 = Ifges(9,5) * t1062;
t1172 = Ifges(9,2) * t1062;
t1171 = Ifges(3,6) * t1075;
t1170 = Ifges(9,6) * t1061;
t1169 = t1020 * Ifges(4,4);
t1168 = t1081 * Ifges(3,2);
t875 = t1231 - t898 + (t902 * t1074 - t1080 * t901) * t1189;
t1167 = t875 * t1237;
t1141 = t1061 * qJD(1);
t1152 = mrSges(9,3) * t1141 + t1229 * t1065 + t1181;
t1149 = qJD(1) * t1062;
t1151 = mrSges(9,3) * t1149 + t1072 * t1065 - t1182;
t1150 = m(4) + t1085;
t1148 = qJD(1) * t1065;
t1116 = mrSges(9,1) * t1062 - mrSges(9,2) * t1061;
t1117 = mrSges(8,1) * t1031 + mrSges(8,2) * t1030;
t1129 = mrSges(4,1) * t1019 + mrSges(4,2) * t1020 + (-t1116 - t1117) * qJD(1);
t1047 = t1229 + t1219;
t1076 = sin(qJ(1));
t1082 = cos(qJ(1));
t1043 = g(1) * t1082 + g(2) * t1076;
t1120 = mrSges(6,1) * t1073 + mrSges(6,2) * t1079;
t1119 = mrSges(7,1) * t1083 - mrSges(7,2) * t1077;
t1118 = mrSges(7,1) * t1077 + mrSges(7,2) * t1083;
t1112 = qJD(1) * t1217 + t1052 * (-mrSges(4,1) * t1020 + mrSges(4,2) * t1019);
t1110 = t1073 * t877 + t1079 * t876;
t917 = -mrSges(6,2) * t949 - t1073 * t1185;
t918 = mrSges(6,1) * t949 - t1079 * t1185;
t1108 = -t1073 * t918 + t1079 * t917;
t1107 = -t1073 * t917 - t1079 * t918;
t975 = t1030 * t1081 + t1031 * t1075;
t974 = -t1075 * t1030 + t1031 * t1081;
t1104 = -pkin(4) * t1032 - t1052;
t1102 = t1061 * (-Ifges(9,1) * t1062 + t1176);
t1101 = t1062 * (Ifges(9,2) * t1061 - t1175);
t1055 = Ifges(3,4) * t1137;
t1100 = -t1075 * (Ifges(3,6) * qJD(2) + (t1168 + t1179) * qJD(1)) / 0.2e1 + (Ifges(3,1) * t1138 + Ifges(3,5) * qJD(2) + t1055) * t1192;
t1096 = t1116 - t1199;
t1095 = -mrSges(8,3) - mrSges(9,3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - t1120;
t906 = mrSges(6,1) * t947 - mrSges(6,3) * t915;
t907 = -mrSges(6,2) * t947 + mrSges(6,3) * t916;
t1094 = qJD(4) * t1107 - t1073 * t906 + t1079 * t907;
t1093 = -pkin(12) * (mrSges(3,1) * t1075 + mrSges(3,2) * t1081) + (Ifges(3,1) * t1081 - t1179) * t1194;
t1058 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t1063 = pkin(8) * m(6) + mrSges(5,1);
t1003 = t1058 * t1069 + t1063 * t1068;
t1004 = -t1058 * t1068 + t1063 * t1069;
t1041 = mrSges(7,1) * t1078 - mrSges(7,2) * t1084;
t1042 = mrSges(7,1) * t1084 + mrSges(7,2) * t1078;
t1048 = mrSges(4,1) + t1219;
t1090 = m(8) * t1052 + (t1003 * t1083 + t1004 * t1077 + t1121 * t1028) * t1059 - (-t1003 * t1077 + t1004 * t1083 + t1121 * t1029) * t1060 + t1075 * (mrSges(4,2) * t1080 + t1041 * t1083 - t1042 * t1077 + t1048 * t1074 - mrSges(3,2)) + t1081 * (t1150 * pkin(1) + mrSges(4,2) * t1074 + t1041 * t1077 + t1042 * t1083 - t1048 * t1080 + mrSges(3,1)) - t1205 + mrSges(2,1) - (-m(3) - t1150) * pkin(12) - t1096 + (-mrSges(8,1) * t1083 - mrSges(8,2) * t1077) * t1071 + t1190 * (mrSges(8,1) * t1077 - mrSges(8,2) * t1083);
t1001 = -qJDD(1) * t1061 - t1062 * t1148;
t1002 = -qJDD(1) * t1062 + t1065 * t1141;
t1008 = Ifges(4,4) * t1019;
t943 = qJD(1) * t1097 - t1038 * t1074 - t1039 * t1080;
t956 = t1019 * Ifges(4,2) + t1065 * Ifges(4,6) - t1169;
t957 = -t1020 * Ifges(4,1) + t1065 * Ifges(4,5) + t1008;
t996 = Ifges(9,6) * t1065 + (-t1172 - t1176) * qJD(1);
t997 = Ifges(9,5) * t1065 + (-t1175 - t1180) * qJD(1);
t1089 = -t1072 * t1027 + t956 * t1195 + Ifges(9,5) * t1001 + t1020 * (Ifges(4,1) * t1019 + t1169) / 0.2e1 + Ifges(4,6) * t944 + Ifges(4,5) * t943 - t1065 * (Ifges(4,5) * t1019 + Ifges(4,6) * t1020) / 0.2e1 + t997 * t1149 / 0.2e1 - t996 * t1141 / 0.2e1 + t1132 * t1181 - (t1170 - t1173) * t1148 / 0.2e1 + Ifges(9,6) * t1002 - (Ifges(4,2) * t1020 + t1008 + t957) * t1019 / 0.2e1 + (t1102 + t1101) * t1207 / 0.2e1 + qJD(1) * t1189 * t1233 + (Ifges(9,3) + Ifges(4,3)) * t1064 + t1229 * t1026;
t1006 = g(3) * t1083 + t1043 * t1077;
t1005 = -g(3) * t1077 + t1043 * t1083;
t991 = g(3) * t1047 - t1043 * t1072;
t990 = pkin(1) * t1138 - t1186;
t985 = g(3) * t1072 + t1043 * t1047;
t966 = t974 * qJD(2);
t965 = t975 * qJD(2);
t946 = t1047 * t1074 + t1072 * t1080 + t1119 * t1078 - t1118 * t1084 - mrSges(3,2);
t939 = -t1047 * t1080 + t1072 * t1074 + t1119 * t1084 + (m(8) + m(9) + t1150) * pkin(1) + mrSges(3,1) + t1118 * t1078;
t932 = qJD(1) * t968 + qJDD(1) * t978;
t931 = -qJD(1) * t967 + qJDD(1) * t1123;
t926 = Ifges(7,5) * qJD(2) + (t1184 + t1218) * qJD(1);
t925 = Ifges(7,6) * qJD(2) + (t1183 + t1228) * qJD(1);
t924 = (qJD(2) * t965 - qJDD(2) * t974) * pkin(1);
t923 = (qJD(2) * t966 + qJDD(2) * t975) * pkin(1);
t914 = t1120 * t953;
t891 = -mrSges(6,1) * t916 + mrSges(6,2) * t915;
t886 = t1054 * t901 - t1133;
t879 = t1073 * t990 + t1079 * t1237;
t878 = -t1073 * t1237 + t1079 * t990;
t874 = -t902 * t1131 + t1236;
t873 = -t1073 * t1186 + t1079 * t874;
t872 = -t1073 * t874 - t1079 * t1186;
t1 = [t1038 * t1168 + t931 * t1183 + t932 * t1184 + (t1038 * t1075 + t1039 * t1081) * Ifges(3,4) + (-t1062 * t1001 - t1061 * t1002) * Ifges(9,4) + t1039 * t1075 * Ifges(3,1) + (m(6) * t1110 - t1107 + t1216) * (-pkin(4) * t982 + t1075 * t1189) + (Ifges(4,5) * t1202 + Ifges(4,6) * t1201 - t1062 * t997 / 0.2e1 + t1061 * t996 / 0.2e1 + (-t1173 / 0.2e1 + t1170 / 0.2e1) * t1065 + (-t1102 / 0.2e1 - t1101 / 0.2e1 - t1217) * qJD(1)) * t1065 + t956 * t1201 + t957 * t1202 + t926 * t1203 + t925 * t1204 + (Ifges(3,5) * t1075 + Ifges(7,5) * t978 + Ifges(3,6) * t1081 + Ifges(7,6) * t1123) * qJDD(2) + (t978 * t931 + t932 * t1123 + (t1123 * t1203 + t978 * t1204) * qJD(1)) * Ifges(7,4) + (Ifges(2,3) + (-t1052 * mrSges(8,1) + Ifges(8,2) * t1031) * t1031 + (-t1052 * mrSges(8,2) + Ifges(8,1) * t1030 - 0.2e1 * Ifges(8,4) * t1031) * t1030 + (m(3) * pkin(12) + mrSges(3,1) * t1081 - mrSges(3,2) * t1075) * pkin(12) + (-mrSges(7,1) * t1123 + mrSges(7,2) * t978 + t1205) * pkin(6)) * qJDD(1) + (-t1052 * (-mrSges(4,1) * t982 + mrSges(4,2) * t981) + pkin(6) * (mrSges(7,1) * t967 + mrSges(7,2) * t968) + t1184 * t1203 + t1183 * t1204) * qJD(1) + (t1076 * t1090 + t1082 * t1095) * g(1) + (t1076 * t1095 - t1082 * t1090) * g(2) + (t1026 * t1105 + t1027 * t1032) * mrSges(4,3) + (-Ifges(4,5) * t1105 - Ifges(9,5) * t1061 + Ifges(4,6) * t1032 - Ifges(9,6) * t1062) * t1064 + (-t1105 * t943 + t981 * t1195) * Ifges(4,1) + (t1032 * t943 - t1105 * t944 + t982 * t1195 + t981 * t1196) * Ifges(4,4) + (-mrSges(4,1) * t1032 - mrSges(4,2) * t1105 + t1117 - t1130) * t1135 + (t1026 * t1061 - t1027 * t1062) * mrSges(9,3) + (m(5) * t933 - mrSges(5,1) * t952 + mrSges(5,2) * t951) * t1104 + t1096 * (-qJDD(1) * t948 + t1046) + pkin(6) * (-mrSges(7,1) * t931 + mrSges(7,2) * t932) + (-t924 * t1030 - t923 * t1031) * mrSges(8,3) - t1052 * (-mrSges(4,1) * t944 + mrSges(4,2) * t943) + (t933 * mrSges(5,2) + t868 * mrSges(5,3) + Ifges(5,1) * t951 + Ifges(5,4) * t952 + (t868 * mrSges(6,2) - t867 * mrSges(6,3) + Ifges(6,1) * t915 + Ifges(6,4) * t916 + Ifges(6,5) * t947) * t1079 + (t868 * mrSges(6,1) - t866 * mrSges(6,3) - Ifges(6,4) * t915 - Ifges(6,2) * t916 - Ifges(6,6) * t947) * t1073 + t1208 * qJD(4)) * t958 + (-t933 * mrSges(5,1) + t869 * mrSges(5,3) + Ifges(5,4) * t951 + Ifges(5,2) * t952 - t1209) * t959 - pkin(12) * (-mrSges(3,1) * t1038 + mrSges(3,2) * t1039) - t948 * (-mrSges(9,1) * t1002 + mrSges(9,2) * t1001) + (t1032 * t944 + t982 * t1196) * Ifges(4,2) - t1001 * t1180 - t1002 * t1172 + (t917 * t1144 + m(6) * (t1073 * t866 + t1079 * t867 + t877 * t1144 - t876 * t1145) - t918 * t1145 + t1079 * t906 + t1073 * t907) * (-pkin(8) * t959 - pkin(10) * t958 + t1104) + (Ifges(7,5) * t1203 + Ifges(7,6) * t1204 + (t1174 / 0.2e1 - t1171 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t1081 - Ifges(3,2) * t1075) * t1192 + t1093) * qJD(1) + (-t1065 * t1233 + (-t1074 * t982 + t1080 * t981) * mrSges(4,3) + (-t1223 * qJD(1) - t1129) * t1075) * pkin(1) + t1100) * qJD(2); t885 * t891 - t879 * t917 - t878 * t918 - (g(3) * t939 + t1043 * t946) * t1081 + (-g(3) * t946 + t1043 * t939) * t1075 + t1089 + (Ifges(7,3) + Ifges(3,3)) * qJDD(2) + Ifges(7,6) * t931 + Ifges(7,5) * t932 + t1094 * t886 + (t1165 * t953 + t1227 * t954 + t885 * t951 + t886 * t952) * mrSges(5,3) + (t926 * t1220 - t1081 * t1055 / 0.2e1 + t978 * t925 / 0.2e1 + (t1168 * t1194 - pkin(6) * (mrSges(7,1) * t978 + mrSges(7,2) * t1123) - t978 * (Ifges(7,1) * t1123 - t1228) / 0.2e1 + (-Ifges(7,2) * t978 + t1218) * t1220 - t1093) * qJD(1) - t1100 + t1112 - (Ifges(7,5) * t1123 - Ifges(7,6) * t978 - t1171 + t1174) * qJD(2) / 0.2e1) * qJD(1) + Ifges(3,6) * t1038 + Ifges(3,5) * t1039 - t990 * t919 + t1165 * t914 + t1108 * t871 + (-t876 * t878 - t877 * t879 + t1212 * t871 + (-qJD(4) * t1110 + t1211) * t886 + t1224) * m(6) + (t1227 * t882 + t869 * t886 - t987 * t990 + t1224) * m(5) + ((t1030 * t974 - t1031 * t975) * qJDD(1) * mrSges(8,3) + t1223 * t1207 * t1075 + (-mrSges(4,3) * t944 - mrSges(9,3) * t1002 + t1152 * qJD(3) + t1072 * t1064) * t1074 + (mrSges(9,3) * t1001 - t1229 * t1064 + t1151 * qJD(3) + (-qJD(2) * t1019 + t943) * mrSges(4,3)) * t1080 + 0.2e1 * (m(9) / 0.2e1 + m(4) / 0.2e1) * (-t1026 * t1080 - t1027 * t1074) + (t1129 * t1075 + (-t1030 * t965 - t1031 * t966 + (t1030 * t975 + t1031 * t974) * qJD(2)) * mrSges(8,3)) * qJD(1) + (t923 * t975 - t924 * t974 + (-t965 * t974 + t966 * t975) * t1189) * m(8)) * pkin(1); (-t1152 * t1074 + (-t1151 - t1182) * t1080) * t1189 - m(6) * (t872 * t876 + t873 * t877 + t1167) - t873 * t917 - t872 * t918 - t875 * t914 + t1112 * qJD(1) - m(5) * (t874 * t882 + t1167) - (t1074 * t991 + t1080 * t985) * t1075 + (-t1074 * t985 + t1080 * t991) * t1081 + t1089 + (-t890 * t914 - t902 * t891 + (-t890 * t953 - t902 * t951) * mrSges(5,3) + t1085 * (-t1237 * t890 - t868 * t902) + (m(5) * t882 + m(6) * t1212 + t954 * mrSges(5,3) + t1108) * t889 + t1216 * t1020 + (t1094 + m(5) * t869 + m(6) * (-t1144 * t876 - t1145 * t877 + t1211) + t952 * mrSges(5,3)) * t901) * pkin(4) + (-t874 * t954 - t875 * t953) * mrSges(5,3); t877 * t918 - t876 * t917 - t1121 * (g(1) * t1076 - g(2) * t1082) + (-(t1005 * t1069 - t1006 * t1068) * t1060 + (t1005 * t1068 + t1006 * t1069) * t1059) * t1120 - t1208 * t953 + t1209;];
tau = t1;
