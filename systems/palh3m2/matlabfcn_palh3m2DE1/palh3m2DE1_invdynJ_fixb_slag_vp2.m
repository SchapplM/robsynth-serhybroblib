% Calculate vector of inverse dynamics joint torques for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2DE1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE1_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:03:22
% EndTime: 2020-05-07 02:03:34
% DurationCPUTime: 11.37s
% Computational Cost: add. (10699->498), mult. (20079->728), div. (0->0), fcn. (24260->22), ass. (0->235)
t1102 = sin(qJ(3));
t1180 = qJD(2) * t1102;
t1160 = pkin(1) * t1180;
t1093 = qJD(3) + qJD(2);
t1108 = cos(qJ(3));
t1219 = pkin(1) * t1108;
t1159 = qJD(2) * t1219;
t1131 = -pkin(4) * t1093 + t1159;
t1094 = pkin(17) + pkin(18);
t1087 = sin(t1094);
t1088 = cos(t1094);
t1103 = sin(qJ(2));
t1109 = cos(qJ(2));
t1096 = sin(pkin(16));
t1097 = cos(pkin(16));
t1105 = sin(pkin(15));
t1111 = cos(pkin(15));
t1056 = t1096 * t1111 + t1097 * t1105;
t1057 = -t1096 * t1105 + t1097 * t1111;
t1134 = t1056 * t1102 - t1057 * t1108;
t1152 = -t1056 * t1108 - t1102 * t1057;
t1241 = t1103 * t1134 + t1109 * t1152;
t1257 = t1134 * t1109;
t1246 = -t1152 * t1103 + t1257;
t930 = t1087 * t1241 - t1246 * t1088;
t1267 = t930 * t1131;
t929 = t1087 * t1246 + t1241 * t1088;
t1268 = -t929 * t1160 + t1267;
t1263 = t929 * t1131;
t926 = t930 * t1160;
t910 = -t926 - t1263;
t1082 = pkin(4) - t1219;
t1220 = pkin(1) * t1102;
t1161 = t930 * t1220;
t1178 = qJD(3) * t1108;
t990 = t1134 * qJD(3);
t991 = t1152 * qJD(3);
t932 = qJD(2) * t1241 + t1103 * t990 + t1109 * t991;
t1194 = t1087 * t932;
t1265 = t1246 * qJD(2) - t991 * t1103;
t931 = -t1109 * t990 - t1265;
t916 = t1088 * t931 + t1194;
t918 = -t1087 * t931 + t1088 * t932;
t1195 = -t918 * t1082 - qJD(3) * t1161 + (t1102 * t916 - t1178 * t929) * pkin(1) + t910;
t1261 = mrSges(4,1) + mrSges(9,1);
t1100 = mrSges(4,2) + mrSges(9,2);
t1221 = pkin(1) * qJD(2);
t1106 = sin(pkin(14));
t1112 = cos(pkin(14));
t1062 = t1105 * t1112 - t1106 * t1111;
t1063 = t1105 * t1106 + t1111 * t1112;
t1006 = t1062 * t1109 + t1063 * t1103;
t1260 = Ifges(7,4) * t1006;
t917 = (qJD(3) * t1257 + t1265) * t1088 - t1194;
t899 = t1082 * t917 + (-t1102 * t918 + (t1102 * t929 - t1108 * t930) * qJD(3)) * pkin(1);
t1259 = t899 - t1268;
t1162 = qJDD(2) * t1108;
t1054 = (qJD(3) * t1180 - t1162) * pkin(1);
t1055 = (-qJD(2) * t1178 - qJDD(2) * t1102) * pkin(1);
t1092 = qJDD(2) + qJDD(3);
t896 = t916 * t1160 + t918 * t1131 - t930 * (pkin(4) * t1092 + t1054) + t929 * t1055;
t913 = -t930 * t1082 - t1220 * t929;
t1256 = t1195 * t1268 + t896 * t913;
t1080 = pkin(1) * t1109 + pkin(12);
t1158 = (m(4) + m(8)) * t1080;
t1099 = cos(pkin(18));
t1222 = sin(pkin(18));
t1058 = t1099 * t1105 + t1222 * t1111;
t1059 = t1099 * t1111 - t1222 * t1105;
t976 = t1080 + (sin(pkin(17)) * t1058 - cos(pkin(17)) * t1059) * pkin(3);
t1233 = m(9) * t976;
t1255 = t1158 + t1233;
t1095 = qJ(3) + qJ(2);
t1089 = sin(t1095);
t1090 = cos(t1095);
t1254 = t1089 * t1102 + t1090 * t1108;
t1200 = Ifges(9,2) * t1090;
t1205 = Ifges(9,4) * t1090;
t1206 = Ifges(9,4) * t1089;
t1211 = Ifges(9,1) * t1089;
t1253 = t1090 * (Ifges(9,5) * t1093 + (-t1205 - t1211) * qJD(1)) / 0.2e1 - t1089 * (Ifges(9,6) * t1093 + (-t1200 - t1206) * qJD(1)) / 0.2e1;
t1113 = m(5) + m(6);
t1151 = -t1062 * t1103 + t1063 * t1109;
t1251 = -t1151 / 0.2e1;
t1250 = pkin(4) * t1113;
t1249 = Ifges(7,4) * t1151;
t1248 = t976 * (-mrSges(9,1) * t1089 - mrSges(9,2) * t1090);
t987 = t1056 * t1087 - t1057 * t1088;
t982 = t987 * qJD(1);
t1060 = t1102 * t1103 - t1108 * t1109;
t1047 = t1060 * qJD(1);
t1015 = -pkin(4) * t1047 - t1080 * qJD(1);
t986 = t1056 * t1088 + t1057 * t1087;
t981 = t986 * qJD(1);
t947 = -mrSges(5,1) * t982 + mrSges(5,2) * t981;
t1247 = m(5) * t1015 + t947;
t980 = t987 * qJDD(1);
t1101 = sin(qJ(4));
t1107 = cos(qJ(4));
t938 = -pkin(8) * t982 - pkin(10) * t981 + t1015;
t904 = -t1101 * t910 + t1107 * t938;
t905 = t1101 * t938 + t1107 * t910;
t1243 = -t1101 * t904 + t1107 * t905;
t897 = t930 * t1055 + (t1092 * t929 + t1093 * t917) * pkin(4) + (-t929 * t1162 + (-t917 * t1108 + (qJD(3) * t929 - t918) * t1102) * qJD(2)) * pkin(1);
t1179 = qJD(2) * t1103;
t1156 = qJD(1) * t1179;
t1074 = pkin(1) * t1156;
t1163 = -t1080 * qJDD(1) + t1074;
t1066 = qJDD(1) * t1109 - t1156;
t1181 = qJD(1) * t1109;
t1067 = qJD(2) * t1181 + qJDD(1) * t1103;
t1133 = t1102 * t1109 + t1103 * t1108;
t1126 = t1133 * qJD(3);
t972 = qJD(1) * t1126 - t1066 * t1108 + t1067 * t1102;
t961 = -pkin(4) * t972 + t1163;
t979 = t986 * qJDD(1);
t933 = -pkin(8) * t980 - pkin(10) * t979 + t961;
t894 = qJD(4) * t904 + t1101 * t933 + t1107 * t897;
t895 = -qJD(4) * t905 - t1101 * t897 + t1107 * t933;
t1242 = -t1101 * t895 + t1107 * t894;
t1177 = qJD(4) * t1101;
t943 = t1107 * t979 - t981 * t1177;
t1176 = qJD(4) * t1107;
t944 = -t1101 * t979 - t981 * t1176;
t975 = qJDD(4) - t980;
t1240 = t895 * mrSges(6,1) - t894 * mrSges(6,2) + Ifges(6,5) * t943 + Ifges(6,6) * t944 + Ifges(6,3) * t975;
t1149 = mrSges(6,1) * t1107 - mrSges(6,2) * t1101;
t977 = qJD(4) - t982;
t1239 = -t1243 * mrSges(6,3) + t1268 * t1149 - (Ifges(6,5) * t1101 + Ifges(6,6) * t1107) * t977 + (-Ifges(6,4) * t1107 ^ 2 + (Ifges(6,4) * t1101 + (-Ifges(6,1) + Ifges(6,2)) * t1107) * t1101) * t981;
t1238 = qJD(1) ^ 2;
t1236 = pkin(6) * m(7);
t995 = t1006 * qJD(2);
t1235 = -t995 / 0.2e1;
t996 = t1151 * qJD(2);
t1234 = t996 / 0.2e1;
t1125 = t1060 * qJD(3);
t1009 = t1060 * qJD(2) + t1125;
t1231 = t1009 / 0.2e1;
t1010 = t1133 * qJD(2) + t1126;
t1230 = t1010 / 0.2e1;
t1228 = t1047 / 0.2e1;
t1048 = t1133 * qJD(1);
t1227 = -t1048 / 0.2e1;
t1225 = t1103 / 0.2e1;
t1224 = t1109 / 0.2e1;
t1218 = pkin(4) * t1048;
t1216 = mrSges(6,3) * t981;
t1215 = mrSges(4,3) * t1047;
t1214 = mrSges(4,3) * t1048;
t1213 = mrSges(9,3) * qJD(1);
t1212 = Ifges(7,1) * t1006;
t1210 = Ifges(3,4) * t1103;
t1209 = Ifges(4,4) * t1048;
t1204 = Ifges(3,5) * t1109;
t1203 = Ifges(9,5) * t1090;
t1202 = Ifges(3,2) * t1109;
t1201 = Ifges(7,2) * t1151;
t1199 = Ifges(3,6) * t1103;
t1198 = Ifges(9,6) * t1089;
t903 = t1263 - t926 + (t930 * t1102 - t1108 * t929) * t1221;
t1197 = t903 * t1268;
t1184 = m(4) + t1113;
t1183 = qJD(1) * t1093;
t1182 = qJD(1) * t1103;
t1171 = t1089 * qJD(1);
t1165 = t1090 * t1213 + t1093 * t1100 - t1215;
t1164 = mrSges(9,3) * t1171 + t1093 * t1261 + t1214;
t1144 = mrSges(9,1) * t1090 - mrSges(9,2) * t1089;
t1145 = mrSges(8,1) * t1059 + mrSges(8,2) * t1058;
t1157 = mrSges(4,1) * t1047 + mrSges(4,2) * t1048 + (-t1144 - t1145) * qJD(1);
t1075 = t1261 + t1250;
t1104 = sin(qJ(1));
t1110 = cos(qJ(1));
t1071 = g(1) * t1110 + g(2) * t1104;
t1148 = mrSges(6,1) * t1101 + mrSges(6,2) * t1107;
t1147 = mrSges(7,1) * t1111 - mrSges(7,2) * t1105;
t1146 = mrSges(7,1) * t1105 + mrSges(7,2) * t1111;
t1140 = qJD(1) * t1248 + t1080 * (-mrSges(4,1) * t1048 + mrSges(4,2) * t1047);
t1138 = t1101 * t905 + t1107 * t904;
t945 = -mrSges(6,2) * t977 - t1101 * t1216;
t946 = mrSges(6,1) * t977 - t1107 * t1216;
t1136 = -t1101 * t946 + t1107 * t945;
t1135 = -t1101 * t945 - t1107 * t946;
t1003 = t1058 * t1109 + t1059 * t1103;
t1002 = -t1103 * t1058 + t1059 * t1109;
t1132 = -pkin(4) * t1060 - t1080;
t1130 = t1089 * (-Ifges(9,1) * t1090 + t1206);
t1129 = t1090 * (Ifges(9,2) * t1089 - t1205);
t1083 = Ifges(3,4) * t1181;
t1128 = -t1103 * (Ifges(3,6) * qJD(2) + (t1202 + t1210) * qJD(1)) / 0.2e1 + (Ifges(3,1) * t1182 + Ifges(3,5) * qJD(2) + t1083) * t1224;
t1124 = t1144 - t1233;
t1123 = -mrSges(8,3) - mrSges(9,3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - t1148;
t934 = mrSges(6,1) * t975 - mrSges(6,3) * t943;
t935 = -mrSges(6,2) * t975 + mrSges(6,3) * t944;
t1122 = t1135 * qJD(4) - t1101 * t934 + t1107 * t935;
t1121 = -pkin(12) * (mrSges(3,1) * t1103 + mrSges(3,2) * t1109) + (Ifges(3,1) * t1109 - t1210) * t1225;
t1086 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t1091 = pkin(8) * m(6) + mrSges(5,1);
t1031 = t1086 * t1097 + t1091 * t1096;
t1032 = -t1086 * t1096 + t1091 * t1097;
t1069 = mrSges(7,1) * t1106 - mrSges(7,2) * t1112;
t1070 = mrSges(7,1) * t1112 + mrSges(7,2) * t1106;
t1076 = mrSges(4,1) + t1250;
t1118 = m(8) * t1080 + (t1031 * t1111 + t1032 * t1105 + t1149 * t1056) * t1087 - (-t1031 * t1105 + t1032 * t1111 + t1149 * t1057) * t1088 + t1103 * (mrSges(4,2) * t1108 + t1069 * t1111 - t1070 * t1105 + t1076 * t1102 - mrSges(3,2)) + t1109 * (t1184 * pkin(1) + mrSges(4,2) * t1102 + t1069 * t1105 + t1070 * t1111 - t1076 * t1108 + mrSges(3,1)) + (-mrSges(8,1) * t1111 - mrSges(8,2) * t1105) * t1099 + t1222 * (mrSges(8,1) * t1105 - mrSges(8,2) * t1111) - t1236 + mrSges(2,1) - (-m(3) - t1184) * pkin(12) - t1124;
t1029 = -qJDD(1) * t1089 - t1090 * t1183;
t1030 = -qJDD(1) * t1090 + t1093 * t1171;
t1036 = Ifges(4,4) * t1047;
t971 = qJD(1) * t1125 - t1066 * t1102 - t1067 * t1108;
t984 = t1047 * Ifges(4,2) + Ifges(4,6) * t1093 - t1209;
t985 = -t1048 * Ifges(4,1) + t1093 * Ifges(4,5) + t1036;
t1117 = -t1100 * t1055 + t984 * t1227 + t1160 * t1214 - (t1198 - t1203) * t1183 / 0.2e1 + Ifges(9,6) * t1030 + Ifges(9,5) * t1029 + t1048 * (Ifges(4,1) * t1047 + t1209) / 0.2e1 + Ifges(4,6) * t972 + Ifges(4,5) * t971 - t1093 * (Ifges(4,5) * t1047 + Ifges(4,6) * t1048) / 0.2e1 - (Ifges(4,2) * t1048 + t1036 + t985) * t1047 / 0.2e1 + (t1130 + t1129) * t1238 / 0.2e1 + t1254 * t1213 * t1221 + (Ifges(9,3) + Ifges(4,3)) * t1092 + t1261 * t1054 + t1253 * qJD(1);
t1034 = g(3) * t1111 + t1071 * t1105;
t1033 = -g(3) * t1105 + t1071 * t1111;
t1019 = g(3) * t1075 - t1071 * t1100;
t1018 = pkin(1) * t1182 - t1218;
t1013 = g(3) * t1100 + t1071 * t1075;
t994 = t1002 * qJD(2);
t993 = t1003 * qJD(2);
t974 = t1075 * t1102 + t1100 * t1108 + t1147 * t1106 - t1146 * t1112 - mrSges(3,2);
t967 = -t1075 * t1108 + t1100 * t1102 + t1147 * t1112 + (m(8) + m(9) + t1184) * pkin(1) + mrSges(3,1) + t1146 * t1106;
t960 = qJD(1) * t996 + qJDD(1) * t1006;
t959 = -qJD(1) * t995 + qJDD(1) * t1151;
t954 = Ifges(7,5) * qJD(2) + (t1212 + t1249) * qJD(1);
t953 = Ifges(7,6) * qJD(2) + (t1201 + t1260) * qJD(1);
t952 = (qJD(2) * t994 + qJDD(2) * t1003) * pkin(1);
t951 = (qJD(2) * t993 - qJDD(2) * t1002) * pkin(1);
t942 = t1148 * t981;
t919 = -mrSges(6,1) * t944 + mrSges(6,2) * t943;
t914 = t1082 * t929 - t1161;
t907 = t1018 * t1101 + t1107 * t1268;
t906 = t1018 * t1107 - t1101 * t1268;
t902 = -t930 * t1159 + t1267;
t901 = -t1101 * t1218 + t1107 * t902;
t900 = -t1101 * t902 - t1107 * t1218;
t1 = [(t1054 * t1089 - t1055 * t1090) * mrSges(9,3) + t984 * t1230 + t985 * t1231 + t954 * t1234 + t953 * t1235 + t960 * t1212 + (Ifges(7,5) * t1234 + Ifges(7,6) * t1235 + (t1204 / 0.2e1 - t1199 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t1109 - Ifges(3,2) * t1103) * t1224 + t1121) * qJD(1) + (-t1254 * t1093 * mrSges(9,3) + (t1009 * t1108 - t1010 * t1102) * mrSges(4,3) + (-qJD(1) * t1255 - t1157) * t1103) * pkin(1) + t1128) * qJD(2) + t959 * t1201 + t1066 * t1202 + (Ifges(4,5) * t1231 + Ifges(4,6) * t1230 + (-t1203 / 0.2e1 + t1198 / 0.2e1) * t1093 + (-t1248 - t1129 / 0.2e1 - t1130 / 0.2e1) * qJD(1) - t1253) * t1093 + t1124 * (-qJDD(1) * t976 + t1074) - t976 * (-mrSges(9,1) * t1030 + mrSges(9,2) * t1029) + (t1066 * t1103 + t1067 * t1109) * Ifges(3,4) + t1067 * Ifges(3,1) * t1103 + (m(6) * t1138 - t1135 + t1247) * (pkin(1) * t1179 - pkin(4) * t1010) + (t1118 * t1104 + t1123 * t1110) * g(1) + (t1123 * t1104 - t1118 * t1110) * g(2) + (-Ifges(4,5) * t1133 - Ifges(9,5) * t1089 + Ifges(4,6) * t1060 - Ifges(9,6) * t1090) * t1092 + (t1054 * t1133 + t1055 * t1060) * mrSges(4,3) + (t1009 * t1228 + t1010 * t1227 + t971 * t1060 - t1133 * t972) * Ifges(4,4) + (t1009 * t1227 - t1133 * t971) * Ifges(4,1) + (-mrSges(4,1) * t1060 - mrSges(4,2) * t1133 + t1145 - t1158) * t1163 + (m(5) * t961 - mrSges(5,1) * t980 + mrSges(5,2) * t979) * t1132 + (-t1058 * t951 - t1059 * t952) * mrSges(8,3) + (Ifges(3,5) * t1103 + Ifges(7,5) * t1006 + Ifges(3,6) * t1109 + Ifges(7,6) * t1151) * qJDD(2) + (t1006 * t959 + t960 * t1151 + (t1006 * t1235 + t1151 * t1234) * qJD(1)) * Ifges(7,4) + (Ifges(2,3) + (-t1080 * mrSges(8,1) + Ifges(8,2) * t1059) * t1059 + (-t1080 * mrSges(8,2) + Ifges(8,1) * t1058 - 0.2e1 * Ifges(8,4) * t1059) * t1058 + (m(3) * pkin(12) + mrSges(3,1) * t1109 - mrSges(3,2) * t1103) * pkin(12) + (-mrSges(7,1) * t1151 + mrSges(7,2) * t1006 + t1236) * pkin(6)) * qJDD(1) + (pkin(6) * (mrSges(7,1) * t995 + mrSges(7,2) * t996) + t1212 * t1234 + t1201 * t1235 - t1080 * (-mrSges(4,1) * t1010 + mrSges(4,2) * t1009)) * qJD(1) - t1080 * (-mrSges(4,1) * t972 + mrSges(4,2) * t971) - pkin(12) * (-mrSges(3,1) * t1066 + mrSges(3,2) * t1067) + (t961 * mrSges(5,2) + t896 * mrSges(5,3) + Ifges(5,1) * t979 + Ifges(5,4) * t980 + (t896 * mrSges(6,2) - t895 * mrSges(6,3) + Ifges(6,1) * t943 + Ifges(6,4) * t944 + Ifges(6,5) * t975) * t1107 + (t896 * mrSges(6,1) - t894 * mrSges(6,3) - Ifges(6,4) * t943 - Ifges(6,2) * t944 - Ifges(6,6) * t975) * t1101 + t1239 * qJD(4)) * t986 + pkin(6) * (-mrSges(7,1) * t959 + mrSges(7,2) * t960) + (-t961 * mrSges(5,1) + t897 * mrSges(5,3) + Ifges(5,4) * t979 + Ifges(5,2) * t980 - t1240) * t987 + (t945 * t1176 + t1101 * t935 - t946 * t1177 + t1107 * t934 + m(6) * (t1101 * t894 + t1107 * t895 + t905 * t1176 - t904 * t1177)) * (-pkin(8) * t987 - pkin(10) * t986 + t1132) + (t1010 * t1228 + t972 * t1060) * Ifges(4,2) - t1029 * t1211 + (-t1029 * t1090 - t1030 * t1089) * Ifges(9,4) - t1030 * t1200; t913 * t919 - t1018 * t947 + t1136 * t899 + t1117 + Ifges(3,6) * t1066 + Ifges(3,5) * t1067 + t1195 * t942 + t1122 * t914 + (t1195 * t981 + t1259 * t982 + t913 * t979 + t914 * t980) * mrSges(5,3) + (t954 * t1251 + t1006 * t953 / 0.2e1 - t1109 * t1083 / 0.2e1 + (-pkin(6) * (mrSges(7,1) * t1006 + mrSges(7,2) * t1151) - t1006 * (Ifges(7,1) * t1151 - t1260) / 0.2e1 + (-Ifges(7,2) * t1006 + t1249) * t1251 + t1202 * t1225 - t1121) * qJD(1) - t1128 + t1140 - (Ifges(7,5) * t1151 - Ifges(7,6) * t1006 - t1199 + t1204) * qJD(2) / 0.2e1) * qJD(1) - t1109 * (g(3) * t967 + t1071 * t974) + (Ifges(7,3) + Ifges(3,3)) * qJDD(2) + Ifges(7,6) * t959 + Ifges(7,5) * t960 + (-g(3) * t974 + t1071 * t967) * t1103 - t907 * t945 - t906 * t946 + (t1243 * t899 + (-t1138 * qJD(4) + t1242) * t914 - t904 * t906 - t905 * t907 + t1256) * m(6) + (-t1015 * t1018 + t1259 * t910 + t897 * t914 + t1256) * m(5) + (0.2e1 * (m(9) / 0.2e1 + m(4) / 0.2e1) * (-t1054 * t1108 - t1055 * t1102) + (t1157 * t1103 + (-t1058 * t993 - t1059 * t994 + (t1002 * t1059 + t1003 * t1058) * qJD(2)) * mrSges(8,3)) * qJD(1) + (t1002 * t1058 - t1003 * t1059) * qJDD(1) * mrSges(8,3) + t1255 * t1238 * t1103 + (-mrSges(4,3) * t972 - mrSges(9,3) * t1030 + t1164 * qJD(3) + t1100 * t1092) * t1102 + (mrSges(9,3) * t1029 - t1261 * t1092 + t1165 * qJD(3) + (-qJD(2) * t1047 + t971) * mrSges(4,3)) * t1108 + (-t1002 * t951 + t1003 * t952 + (-t1002 * t993 + t1003 * t994) * t1221) * m(8)) * pkin(1); -t1103 * (t1013 * t1108 + t1019 * t1102) + (-t1013 * t1102 + t1019 * t1108) * t1109 + (-t1164 * t1102 + (-t1165 - t1215) * t1108) * t1221 + t1117 - m(6) * (t900 * t904 + t901 * t905 + t1197) + (-t902 * t982 - t903 * t981) * mrSges(5,3) - m(5) * (t902 * t910 + t1197) + (-t918 * t942 - t930 * t919 + (-t918 * t981 - t930 * t979) * mrSges(5,3) + t1113 * (-t1268 * t918 - t896 * t930) + (m(5) * t910 + m(6) * t1243 + t982 * mrSges(5,3) + t1136) * t917 + t1247 * t1048 + (t1122 + m(5) * t897 + m(6) * (-t1176 * t904 - t1177 * t905 + t1242) + t980 * mrSges(5,3)) * t929) * pkin(4) + t1140 * qJD(1) - t903 * t942 - t901 * t945 - t900 * t946; t905 * t946 - t904 * t945 - t1149 * (g(1) * t1104 - g(2) * t1110) + (-(t1033 * t1097 - t1034 * t1096) * t1088 + (t1033 * t1096 + t1034 * t1097) * t1087) * t1148 - t1239 * t981 + t1240;];
tau = t1;
