% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
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
% tauJB [(6+10)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = palh3m2OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_invdynJB_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:35:52
% EndTime: 2020-05-07 04:36:23
% DurationCPUTime: 7.88s
% Computational Cost: add. (90599->496), mult. (199384->635), div. (0->0), fcn. (153272->18), ass. (0->203)
t1163 = sin(qJ(1));
t1171 = cos(qJ(1));
t1136 = -t1171 * g(1) - t1163 * g(2);
t1172 = qJD(1) ^ 2;
t1122 = -t1172 * pkin(12) + t1136;
t1162 = sin(qJ(2));
t1170 = cos(qJ(2));
t1099 = -t1170 * g(3) - t1162 * t1122;
t1101 = -t1162 * g(3) + t1170 * t1122;
t1195 = qJD(1) * t1170;
t1127 = qJD(2) * t1195 + t1162 * qJDD(1);
t1197 = qJD(1) * t1162;
t1191 = qJD(2) * t1197;
t1129 = t1170 * qJDD(1) - t1191;
t1092 = (t1162 * t1170 * t1172 + qJDD(2)) * pkin(1) + t1099;
t1093 = (-t1170 ^ 2 * t1172 - qJD(2) ^ 2) * pkin(1) + t1101;
t1161 = sin(qJ(3));
t1169 = cos(qJ(3));
t1050 = -t1169 * t1092 + t1161 * t1093;
t1052 = -t1161 * t1092 - t1169 * t1093;
t1113 = (-t1161 * t1170 - t1162 * t1169) * qJD(1);
t1068 = -t1113 * qJD(3) + t1161 * t1127 - t1169 * t1129;
t1111 = (t1161 * t1162 - t1169 * t1170) * qJD(1);
t1070 = t1111 * qJD(3) - t1169 * t1127 - t1161 * t1129;
t1153 = qJD(2) + qJD(3);
t1076 = Ifges(4,4) * t1113 + Ifges(4,2) * t1111 + Ifges(4,6) * t1153;
t1078 = Ifges(4,1) * t1113 + Ifges(4,4) * t1111 + Ifges(4,5) * t1153;
t1151 = qJDD(2) + qJDD(3);
t1160 = sin(qJ(4));
t1168 = cos(qJ(4));
t1084 = t1160 * t1111 + t1168 * t1113;
t1010 = -t1084 * qJD(4) + t1168 * t1068 - t1160 * t1070;
t1083 = t1168 * t1111 - t1160 * t1113;
t1011 = t1083 * qJD(4) + t1160 * t1068 + t1168 * t1070;
t1147 = qJD(4) + t1153;
t1035 = Ifges(5,4) * t1084 + Ifges(5,2) * t1083 + Ifges(5,6) * t1147;
t1036 = Ifges(5,1) * t1084 + Ifges(5,4) * t1083 + Ifges(5,5) * t1147;
t1145 = qJDD(4) + t1151;
t1159 = sin(qJ(5));
t1167 = cos(qJ(5));
t1059 = -t1159 * t1084 + t1167 * t1147;
t1080 = qJD(5) - t1083;
t1030 = -t1080 * mrSges(6,2) + t1059 * mrSges(6,3);
t1060 = t1167 * t1084 + t1159 * t1147;
t1031 = t1080 * mrSges(6,1) - t1060 * mrSges(6,3);
t1040 = -t1083 * pkin(8) - t1084 * pkin(10);
t1143 = t1147 ^ 2;
t1024 = (t1111 * t1113 + t1151) * pkin(4) + t1050;
t1033 = (-t1111 ^ 2 - t1153 ^ 2) * pkin(4) + t1052;
t990 = t1168 * t1024 - t1160 * t1033;
t981 = -t1145 * pkin(8) - t1143 * pkin(10) + t1084 * t1040 - t990;
t987 = -t1060 * qJD(5) - t1159 * t1011 + t1167 * t1145;
t988 = t1059 * qJD(5) + t1167 * t1011 + t1159 * t1145;
t1183 = -m(6) * t981 + t987 * mrSges(6,1) - t988 * mrSges(6,2) + t1059 * t1030 - t1060 * t1031;
t1008 = qJDD(5) - t1010;
t1027 = -t1059 * mrSges(6,1) + t1060 * mrSges(6,2);
t1135 = t1163 * g(1) - t1171 * g(2);
t1120 = -qJDD(1) * pkin(12) - t1135;
t1089 = t1120 + (-t1129 + t1191) * pkin(1);
t1029 = t1089 + (t1113 * t1153 - t1068) * pkin(4);
t980 = (-t1083 * t1147 - t1011) * pkin(10) + (t1084 * t1147 - t1010) * pkin(8) + t1029;
t991 = t1160 * t1024 + t1168 * t1033;
t982 = -t1143 * pkin(8) + t1145 * pkin(10) + t1083 * t1040 + t991;
t971 = -t1159 * t982 + t1167 * t980;
t965 = m(6) * t971 + t1008 * mrSges(6,1) - t988 * mrSges(6,3) - t1060 * t1027 + t1080 * t1030;
t972 = t1159 * t980 + t1167 * t982;
t966 = m(6) * t972 - t1008 * mrSges(6,2) + t987 * mrSges(6,3) + t1059 * t1027 - t1080 * t1031;
t1190 = -t1159 * t965 + t1167 * t966;
t1001 = Ifges(6,5) * t1060 + Ifges(6,6) * t1059 + Ifges(6,3) * t1080;
t1003 = Ifges(6,1) * t1060 + Ifges(6,4) * t1059 + Ifges(6,5) * t1080;
t950 = -mrSges(6,1) * t981 + mrSges(6,3) * t972 + Ifges(6,4) * t988 + Ifges(6,2) * t987 + Ifges(6,6) * t1008 - t1060 * t1001 + t1080 * t1003;
t1002 = Ifges(6,4) * t1060 + Ifges(6,2) * t1059 + Ifges(6,6) * t1080;
t951 = mrSges(6,2) * t981 - mrSges(6,3) * t971 + Ifges(6,1) * t988 + Ifges(6,4) * t987 + Ifges(6,5) * t1008 + t1059 * t1001 - t1080 * t1002;
t1179 = -pkin(8) * t1183 - pkin(10) * t1190 - mrSges(5,1) * t990 + mrSges(5,2) * t991 - Ifges(5,5) * t1011 - Ifges(5,6) * t1010 - Ifges(5,3) * t1145 - t1084 * t1035 + t1083 * t1036 - t1159 * t951 - t1167 * t950;
t1039 = -t1083 * mrSges(5,1) + t1084 * mrSges(5,2);
t1072 = t1147 * mrSges(5,1) - t1084 * mrSges(5,3);
t944 = m(5) * t991 - t1145 * mrSges(5,2) + t1010 * mrSges(5,3) + t1083 * t1039 - t1147 * t1072 + t1190;
t1071 = -t1147 * mrSges(5,2) + t1083 * mrSges(5,3);
t957 = m(5) * t990 + t1145 * mrSges(5,1) - t1011 * mrSges(5,3) - t1084 * t1039 + t1147 * t1071 + t1183;
t1201 = t1160 * t944 + t1168 * t957;
t1176 = -pkin(4) * t1201 - mrSges(4,1) * t1050 + mrSges(4,2) * t1052 - Ifges(4,5) * t1070 - Ifges(4,6) * t1068 - Ifges(4,3) * t1151 - t1113 * t1076 + t1111 * t1078 + t1179;
t1157 = sin(qJ(7));
t1165 = cos(qJ(7));
t1049 = t1165 * t1092 - t1157 * t1093;
t1051 = t1157 * t1092 + t1165 * t1093;
t1112 = (t1157 * t1170 + t1162 * t1165) * qJD(1);
t1067 = -t1112 * qJD(7) - t1157 * t1127 + t1165 * t1129;
t1110 = (-t1157 * t1162 + t1165 * t1170) * qJD(1);
t1069 = t1110 * qJD(7) + t1165 * t1127 + t1157 * t1129;
t1152 = qJD(2) + qJD(7);
t1075 = Ifges(8,4) * t1112 + Ifges(8,2) * t1110 + Ifges(8,6) * t1152;
t1077 = Ifges(8,1) * t1112 + Ifges(8,4) * t1110 + Ifges(8,5) * t1152;
t1150 = qJDD(2) + qJDD(7);
t1155 = sin(pkin(15));
t1156 = cos(pkin(15));
t1164 = cos(qJ(8));
t1207 = sin(qJ(8));
t1118 = t1155 * t1164 - t1156 * t1207;
t1119 = -t1155 * t1207 - t1156 * t1164;
t1055 = t1119 * t1110 - t1118 * t1112;
t1056 = t1118 * t1110 + t1119 * t1112;
t1146 = qJD(8) + t1152;
t1015 = Ifges(9,4) * t1056 + Ifges(9,2) * t1055 + Ifges(9,6) * t1146;
t1016 = Ifges(9,1) * t1056 + Ifges(9,4) * t1055 + Ifges(9,5) * t1146;
t1144 = qJDD(8) + t1150;
t1079 = (-t1110 * t1156 - t1112 * t1155) * pkin(3);
t1149 = t1152 ^ 2;
t1012 = -t1112 * t1079 + (t1149 * t1155 + t1150 * t1156) * pkin(3) + t1049;
t1013 = t1110 * t1079 + (-t1149 * t1156 + t1150 * t1155) * pkin(3) + t1051;
t984 = t1119 * t1012 - t1118 * t1013;
t985 = t1118 * t1012 + t1119 * t1013;
t998 = -t1056 * qJD(8) + t1119 * t1067 - t1118 * t1069;
t999 = t1055 * qJD(8) + t1118 * t1067 + t1119 * t1069;
t1182 = -mrSges(9,1) * t984 + mrSges(9,2) * t985 - Ifges(9,5) * t999 - Ifges(9,6) * t998 - Ifges(9,3) * t1144 - t1056 * t1015 + t1055 * t1016;
t1020 = -t1055 * mrSges(9,1) + t1056 * mrSges(9,2);
t1053 = -t1146 * mrSges(9,2) + t1055 * mrSges(9,3);
t978 = m(9) * t984 + t1144 * mrSges(9,1) - t999 * mrSges(9,3) - t1056 * t1020 + t1146 * t1053;
t1054 = t1146 * mrSges(9,1) - t1056 * mrSges(9,3);
t979 = m(9) * t985 - t1144 * mrSges(9,2) + t998 * mrSges(9,3) + t1055 * t1020 - t1146 * t1054;
t1199 = -t1118 * t978 + t1119 * t979;
t1200 = t1118 * t979 + t1119 * t978;
t1205 = pkin(3) * t1156;
t1206 = pkin(3) * t1155;
t1177 = -mrSges(8,1) * t1049 + mrSges(8,2) * t1051 - Ifges(8,5) * t1069 - Ifges(8,6) * t1067 - Ifges(8,3) * t1150 - t1112 * t1075 + t1110 * t1077 - t1199 * t1206 - t1200 * t1205 + t1182;
t1088 = -t1111 * mrSges(4,1) + t1113 * mrSges(4,2);
t1095 = -t1153 * mrSges(4,2) + t1111 * mrSges(4,3);
t933 = m(4) * t1050 + t1151 * mrSges(4,1) - t1070 * mrSges(4,3) - t1113 * t1088 + t1153 * t1095 + t1201;
t1097 = t1153 * mrSges(4,1) - t1113 * mrSges(4,3);
t934 = m(4) * t1052 - t1151 * mrSges(4,2) + t1068 * mrSges(4,3) + t1111 * t1088 - t1153 * t1097 - t1160 * t957 + t1168 * t944;
t1087 = -t1110 * mrSges(8,1) + t1112 * mrSges(8,2);
t1094 = -t1152 * mrSges(8,2) + t1110 * mrSges(8,3);
t954 = m(8) * t1049 + t1150 * mrSges(8,1) - t1069 * mrSges(8,3) - t1112 * t1087 + t1152 * t1094 + t1200;
t1096 = t1152 * mrSges(8,1) - t1112 * mrSges(8,3);
t955 = m(8) * t1051 - t1150 * mrSges(8,2) + t1067 * mrSges(8,3) + t1110 * t1087 - t1152 * t1096 + t1199;
t1202 = t1157 * t955 - t1161 * t934 + t1165 * t954 - t1169 * t933;
t1208 = -t1202 * pkin(1) - mrSges(3,1) * t1099 + mrSges(3,2) * t1101 - Ifges(3,5) * t1127 - Ifges(3,6) * t1129 - Ifges(3,3) * qJDD(2) + t1176 + t1177;
t1123 = t1172 * pkin(6) + t1136;
t1158 = sin(qJ(6));
t1166 = cos(qJ(6));
t1098 = -t1166 * g(3) - t1158 * t1123;
t1124 = (-mrSges(7,1) * t1166 + mrSges(7,2) * t1158) * qJD(1);
t1196 = qJD(1) * t1166;
t1126 = qJD(6) * t1196 + t1158 * qJDD(1);
t1133 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t1196;
t1198 = qJD(1) * t1158;
t1047 = m(7) * t1098 + qJDD(6) * mrSges(7,1) - t1126 * mrSges(7,3) + qJD(6) * t1133 - t1124 * t1198;
t1100 = -t1158 * g(3) + t1166 * t1123;
t1128 = -qJD(6) * t1198 + t1166 * qJDD(1);
t1131 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t1198;
t1048 = m(7) * t1100 - qJDD(6) * mrSges(7,2) + t1128 * mrSges(7,3) - qJD(6) * t1131 + t1124 * t1196;
t1188 = -t1158 * t1047 + t1166 * t1048;
t1125 = (-mrSges(3,1) * t1170 + mrSges(3,2) * t1162) * qJD(1);
t1134 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1195;
t928 = m(3) * t1099 + qJDD(2) * mrSges(3,1) - t1127 * mrSges(3,3) + qJD(2) * t1134 - t1125 * t1197 + t1202;
t1132 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1197;
t929 = m(3) * t1101 - qJDD(2) * mrSges(3,2) + t1129 * mrSges(3,3) - qJD(2) * t1132 + t1125 * t1195 - t1157 * t954 + t1161 * t933 + t1165 * t955 - t1169 * t934;
t923 = m(2) * t1136 - t1172 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t1162 * t928 + t1170 * t929 + t1188;
t946 = t1159 * t966 + t1167 * t965;
t1184 = m(5) * t1029 - t1010 * mrSges(5,1) + t1011 * mrSges(5,2) - t1083 * t1071 + t1084 * t1072 + t946;
t1000 = (-t1067 * t1156 - t1069 * t1155 + (-t1110 * t1155 + t1112 * t1156) * t1152) * pkin(3) + t1089;
t973 = m(9) * t1000 - t998 * mrSges(9,1) + t999 * mrSges(9,2) - t1055 * t1053 + t1056 * t1054;
t1175 = t1110 * t1094 + t1111 * t1095 - t1184 + t1068 * mrSges(4,1) + (-m(4) - m(8)) * t1089 - t1112 * t1096 - t1113 * t1097 - t1069 * mrSges(8,2) - t1070 * mrSges(4,2) + t1067 * mrSges(8,1) - t973;
t1174 = -m(3) * t1120 + t1129 * mrSges(3,1) - t1127 * mrSges(3,2) + t1134 * t1195 + t1175;
t1121 = qJDD(1) * pkin(6) - t1135;
t1185 = -m(7) * t1121 + t1128 * mrSges(7,1) - t1126 * mrSges(7,2) + t1133 * t1196;
t1193 = t1162 * t1132;
t1194 = t1158 * t1131;
t938 = t1174 + qJDD(1) * mrSges(2,1) + (-t1193 - t1194) * qJD(1) - t1172 * mrSges(2,2) + m(2) * t1135 + t1185;
t1204 = t1163 * t923 + t1171 * t938;
t1203 = t1162 * t929 + t1170 * t928;
t1192 = t1166 * t1047 + t1158 * t1048;
t1189 = -t1163 * t938 + t1171 * t923;
t1106 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t1158 + Ifges(7,2) * t1166) * qJD(1);
t1108 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t1158 + Ifges(7,4) * t1166) * qJD(1);
t1187 = t1158 * t1106 - t1166 * t1108;
t1107 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t1162 + Ifges(3,2) * t1170) * qJD(1);
t1109 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t1162 + Ifges(3,4) * t1170) * qJD(1);
t1186 = t1162 * t1107 - t1170 * t1109;
t1104 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t1158 + Ifges(7,6) * t1166) * qJD(1);
t1025 = -mrSges(7,1) * t1121 + mrSges(7,3) * t1100 + Ifges(7,4) * t1126 + Ifges(7,2) * t1128 + Ifges(7,6) * qJDD(6) + qJD(6) * t1108 - t1104 * t1198;
t1026 = mrSges(7,2) * t1121 - mrSges(7,3) * t1098 + Ifges(7,1) * t1126 + Ifges(7,4) * t1128 + Ifges(7,5) * qJDD(6) - qJD(6) * t1106 + t1104 * t1196;
t1044 = -qJD(1) * t1194 + t1185;
t1105 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t1162 + Ifges(3,6) * t1170) * qJD(1);
t1074 = Ifges(4,5) * t1113 + Ifges(4,6) * t1111 + Ifges(4,3) * t1153;
t1034 = Ifges(5,5) * t1084 + Ifges(5,6) * t1083 + Ifges(5,3) * t1147;
t931 = -pkin(10) * t946 + mrSges(5,2) * t1029 - mrSges(5,3) * t990 + Ifges(5,1) * t1011 + Ifges(5,4) * t1010 + Ifges(5,5) * t1145 + t1083 * t1034 - t1147 * t1035 - t1159 * t950 + t1167 * t951;
t1178 = mrSges(6,1) * t971 - mrSges(6,2) * t972 + Ifges(6,5) * t988 + Ifges(6,6) * t987 + Ifges(6,3) * t1008 + t1060 * t1002 - t1059 * t1003;
t932 = -pkin(8) * t946 - mrSges(5,1) * t1029 + mrSges(5,3) * t991 + Ifges(5,4) * t1011 + Ifges(5,2) * t1010 + Ifges(5,6) * t1145 - t1084 * t1034 + t1147 * t1036 - t1178;
t924 = -pkin(4) * t1184 - mrSges(4,1) * t1089 + mrSges(4,3) * t1052 + Ifges(4,4) * t1070 + Ifges(4,2) * t1068 + Ifges(4,6) * t1151 - t1113 * t1074 + t1153 * t1078 + t1160 * t931 + t1168 * t932;
t925 = mrSges(4,2) * t1089 - mrSges(4,3) * t1050 + Ifges(4,1) * t1070 + Ifges(4,4) * t1068 + Ifges(4,5) * t1151 + t1111 * t1074 - t1153 * t1076 - t1160 * t932 + t1168 * t931;
t1073 = Ifges(8,5) * t1112 + Ifges(8,6) * t1110 + Ifges(8,3) * t1152;
t1014 = Ifges(9,5) * t1056 + Ifges(9,6) * t1055 + Ifges(9,3) * t1146;
t967 = -mrSges(9,1) * t1000 + mrSges(9,3) * t985 + Ifges(9,4) * t999 + Ifges(9,2) * t998 + Ifges(9,6) * t1144 - t1056 * t1014 + t1146 * t1016;
t968 = mrSges(9,2) * t1000 - mrSges(9,3) * t984 + Ifges(9,1) * t999 + Ifges(9,4) * t998 + Ifges(9,5) * t1144 + t1055 * t1014 - t1146 * t1015;
t940 = -mrSges(8,1) * t1089 + mrSges(8,3) * t1051 + Ifges(8,4) * t1069 + Ifges(8,2) * t1067 + Ifges(8,6) * t1150 - t1112 * t1073 + t1152 * t1077 + t1118 * t968 + t1119 * t967 - t973 * t1205;
t941 = mrSges(8,2) * t1089 - mrSges(8,3) * t1049 + Ifges(8,1) * t1069 + Ifges(8,4) * t1067 + Ifges(8,5) * t1150 + t1110 * t1073 - t1152 * t1075 - t1118 * t967 + t1119 * t968 - t973 * t1206;
t918 = t1175 * pkin(1) - mrSges(3,1) * t1120 + mrSges(3,3) * t1101 + Ifges(3,4) * t1127 + Ifges(3,2) * t1129 + Ifges(3,6) * qJDD(2) + qJD(2) * t1109 - t1105 * t1197 + t1157 * t941 - t1161 * t925 + t1165 * t940 - t1169 * t924;
t920 = mrSges(3,2) * t1120 - mrSges(3,3) * t1099 + Ifges(3,1) * t1127 + Ifges(3,4) * t1129 + Ifges(3,5) * qJDD(2) - qJD(2) * t1107 + t1105 * t1195 - t1157 * t940 + t1161 * t924 + t1165 * t941 - t1169 * t925;
t1181 = -pkin(6) * t1044 - mrSges(2,2) * t1136 + t1162 * t920 + t1170 * t918 + t1158 * t1026 + t1166 * t1025 + mrSges(2,1) * t1135 + Ifges(2,3) * qJDD(1) + pkin(12) * (-qJD(1) * t1193 + t1174);
t1180 = mrSges(7,1) * t1098 - mrSges(7,2) * t1100 + Ifges(7,5) * t1126 + Ifges(7,6) * t1128 + Ifges(7,3) * qJDD(6);
t916 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-t1186 - t1187) * qJD(1) + t1172 * Ifges(2,5) + pkin(13) * t1188 + mrSges(2,3) * t1136 + pkin(6) * t1192 - pkin(12) * t1203 - t1180 + t1208;
t915 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1135 - pkin(13) * t1044 + Ifges(2,5) * qJDD(1) - t1172 * Ifges(2,6) - t1158 * t1025 + t1166 * t1026 - t1162 * t918 + t1170 * t920;
t1 = [-m(1) * g(1) + t1189; -m(1) * g(2) + t1204; (-m(1) - m(2)) * g(3) + t1192 + t1203; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(11) * t1204 - t1163 * t916 + t1171 * t915; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(11) * t1189 + t1163 * t915 + t1171 * t916; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1181; t1181; t1186 * qJD(1) - t1208; -t1176; -t1179; t1178; t1187 * qJD(1) + t1180; -t1177; -t1182; 0; 0;];
tauJB = t1;
