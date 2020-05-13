% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
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
% 
% Output:
% tauB_reg [6x(9*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = palh3m2OL_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_reg2_snew_vp: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_reg2_snew_vp: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynB_fixb_reg2_snew_vp: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynB_fixb_reg2_snew_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:41:48
% EndTime: 2020-05-07 04:43:43
% DurationCPUTime: 38.33s
% Computational Cost: add. (72829->953), mult. (162254->1536), div. (0->0), fcn. (129842->18), ass. (0->755)
t1144 = sin(qJ(3));
t1190 = qJDD(2) + qJDD(3);
t1152 = cos(qJ(3));
t1153 = cos(qJ(2));
t1203 = t1152 * t1153;
t1145 = sin(qJ(2));
t1215 = t1144 * t1145;
t1056 = (-t1203 + t1215) * qJD(1);
t1210 = t1145 * t1152;
t1214 = t1144 * t1153;
t1058 = (-t1210 - t1214) * qJD(1);
t1232 = t1056 * t1058;
t1268 = t1190 + t1232;
t1285 = t1144 * t1268;
t1284 = t1152 * t1268;
t1143 = sin(qJ(4));
t1151 = cos(qJ(4));
t1146 = sin(qJ(1));
t1154 = cos(qJ(1));
t1095 = t1154 * g(1) + t1146 * g(2);
t1157 = qJD(1) ^ 2;
t1072 = -t1157 * pkin(12) - t1095;
t1039 = -t1145 * g(3) + t1153 * t1072;
t1138 = t1153 ^ 2;
t1127 = t1138 * t1157;
t1156 = qJD(2) ^ 2;
t1106 = -t1127 - t1156;
t1004 = pkin(1) * t1106 + t1039;
t1037 = t1153 * g(3) + t1145 * t1072;
t1110 = t1153 * t1157 * t1145;
t1092 = qJDD(2) + t1110;
t997 = pkin(1) * t1092 - t1037;
t892 = t1144 * t1004 - t1152 * t997;
t1158 = pkin(4) * t1268 + t892;
t894 = -t1152 * t1004 - t1144 * t997;
t1049 = t1056 ^ 2;
t1133 = qJD(2) + qJD(3);
t1130 = t1133 ^ 2;
t982 = -t1049 - t1130;
t866 = pkin(4) * t982 + t894;
t764 = t1143 * t866 - t1151 * t1158;
t765 = t1143 * t1158 + t1151 * t866;
t653 = t1143 * t765 - t1151 * t764;
t1283 = t1144 * t653;
t1140 = sin(qJ(7));
t1148 = cos(qJ(7));
t891 = t1140 * t1004 - t1148 * t997;
t893 = t1148 * t1004 + t1140 * t997;
t809 = t1140 * t893 - t1148 * t891;
t1282 = t1145 * t809;
t813 = t1144 * t894 + t1152 * t892;
t1281 = t1145 * t813;
t1280 = t1152 * t653;
t1279 = t1153 * t809;
t1278 = t1153 * t813;
t1225 = t1133 * t1058;
t1122 = t1145 * qJDD(1);
t1234 = qJD(1) * qJD(2);
t1175 = t1153 * t1234;
t1077 = t1122 + t1175;
t1124 = t1153 * qJDD(1);
t1179 = t1145 * t1234;
t1080 = t1124 - t1179;
t943 = -t1058 * qJD(3) + t1144 * t1077 - t1152 * t1080;
t903 = t943 - t1225;
t1139 = sin(pkin(15));
t1147 = cos(qJ(8));
t1259 = cos(pkin(15));
t1260 = sin(qJ(8));
t1065 = t1139 * t1147 - t1259 * t1260;
t1131 = qJDD(2) + qJDD(7);
t1170 = qJDD(8) + t1131;
t1206 = t1148 * t1153;
t1220 = t1140 * t1145;
t1054 = (-t1206 + t1220) * qJD(1);
t1211 = t1145 * t1148;
t1219 = t1140 * t1153;
t1057 = (t1211 + t1219) * qJD(1);
t1066 = -t1139 * t1260 - t1259 * t1147;
t932 = -t1066 * t1054 - t1065 * t1057;
t933 = -t1065 * t1054 + t1066 * t1057;
t846 = t932 * t933;
t1265 = t846 + t1170;
t1277 = t1065 * t1265;
t1276 = t1066 * t1265;
t998 = t1057 * t1054;
t1263 = -t998 + t1131;
t1275 = t1140 * t1263;
t1142 = sin(qJ(5));
t1161 = -t1056 * qJD(3) + t1152 * t1077 + t1144 * t1080;
t1171 = -t1143 * t1161 - t1151 * t943;
t980 = t1143 * t1056 + t1151 * t1058;
t807 = -t980 * qJD(4) - t1171;
t1160 = qJDD(5) - t807;
t1121 = qJD(4) + t1133;
t1150 = cos(qJ(5));
t938 = -t1150 * t1121 + t1142 * t980;
t940 = t1142 * t1121 + t1150 * t980;
t856 = t940 * t938;
t1266 = t1160 - t856;
t1274 = t1142 * t1266;
t1118 = qJDD(4) + t1190;
t978 = -t1151 * t1056 + t1143 * t1058;
t877 = t980 * t978;
t1264 = -t877 + t1118;
t1273 = t1143 * t1264;
t1272 = t1148 * t1263;
t1271 = t1150 * t1266;
t1270 = t1151 * t1264;
t1132 = qJD(2) + qJD(7);
t1042 = t1132 * t1054;
t1159 = t1054 * qJD(7) - t1148 * t1077 - t1140 * t1080;
t1262 = -t1042 - t1159;
t1269 = t1080 - t1179;
t1168 = t1140 * t1077 - t1148 * t1080;
t942 = -t1057 * qJD(7) - t1168;
t783 = t932 * qJD(8) + t1065 * t942 - t1066 * t1159;
t1120 = qJD(8) + t1132;
t917 = t1120 * t932;
t1267 = t783 + t917;
t808 = -t978 * qJD(4) + t1143 * t943 - t1151 * t1161;
t961 = t1121 * t978;
t777 = t808 - t961;
t762 = -t938 * qJD(5) + t1142 * t1118 + t1150 * t808;
t973 = qJD(5) + t978;
t865 = t973 * t938;
t708 = -t865 + t762;
t1162 = t933 * qJD(8) - t1065 * t1159 - t1066 * t942;
t918 = t1120 * t933;
t739 = t1162 - t918;
t1043 = t1133 * t1056;
t908 = t1161 - t1043;
t1169 = -t1150 * t1118 + t1142 * t808;
t705 = (qJD(5) - t973) * t940 + t1169;
t775 = (qJD(4) - t1121) * t980 + t1171;
t900 = (qJD(7) - t1132) * t1057 + t1168;
t930 = t932 ^ 2;
t931 = t933 ^ 2;
t936 = t938 ^ 2;
t937 = t940 ^ 2;
t972 = t973 ^ 2;
t976 = t978 ^ 2;
t977 = t980 ^ 2;
t1048 = t1054 ^ 2;
t1050 = t1057 ^ 2;
t1051 = t1058 ^ 2;
t1117 = t1120 ^ 2;
t1261 = t1121 ^ 2;
t1129 = t1132 ^ 2;
t1258 = pkin(3) * t1139;
t1257 = pkin(8) * t1143;
t949 = t1153 * t1037 - t1145 * t1039;
t1256 = pkin(12) * t949;
t1094 = t1146 * g(1) - t1154 * g(2);
t1069 = qJDD(1) * pkin(12) + t1094;
t995 = pkin(1) * t1269 + t1069;
t854 = pkin(4) * t903 + t995;
t634 = -t777 * pkin(10) + (t1121 * t980 - t807) * pkin(8) - t854;
t871 = t978 * pkin(8) - t980 * pkin(10);
t688 = -t1261 * pkin(8) + t1118 * pkin(10) - t978 * t871 + t765;
t580 = t1142 * t634 + t1150 * t688;
t971 = (t1259 * t1054 - t1057 * t1139) * pkin(3);
t828 = -t1057 * t971 + (t1129 * t1139 + t1259 * t1131) * pkin(3) - t891;
t829 = -t1054 * t971 + (-t1259 * t1129 + t1131 * t1139) * pkin(3) + t893;
t700 = t1065 * t828 + t1066 * t829;
t1188 = pkin(3) * t1259;
t784 = -t1057 * t1132 * t1188 + (t1139 * t1262 + t1259 * t942) * pkin(3) + t995;
t1255 = t1065 * t784;
t843 = -t846 + t1170;
t1254 = t1065 * t843;
t1253 = t1066 * t784;
t1252 = t1066 * t843;
t985 = t998 + t1131;
t1251 = t1140 * t985;
t687 = -t1118 * pkin(8) - t1261 * pkin(10) + t980 * t871 + t764;
t1250 = t1142 * t687;
t744 = t1160 + t856;
t1249 = t1142 * t744;
t1248 = t1142 * t973;
t1247 = t1143 * t854;
t869 = t877 + t1118;
t1246 = t1143 * t869;
t988 = -t1232 + t1190;
t1245 = t1144 * t988;
t1244 = t1146 * t995;
t1243 = t1148 * t985;
t1242 = t1150 * t687;
t1241 = t1150 * t744;
t1240 = t1150 * t973;
t1239 = t1151 * t854;
t1238 = t1151 * t869;
t1237 = t1152 * t854;
t1236 = t1152 * t988;
t1233 = qJD(1) * qJD(6);
t1231 = t1065 * t1120;
t1230 = t1066 * t1120;
t1229 = t1121 * t1143;
t1228 = t1121 * t1151;
t1227 = t1132 * t1140;
t1226 = t1132 * t1148;
t1224 = t1133 * t1144;
t1223 = t1133 * t1152;
t1141 = sin(qJ(6));
t1135 = t1141 ^ 2;
t1222 = t1135 * t1157;
t1136 = t1145 ^ 2;
t1221 = t1136 * t1157;
t1070 = -qJDD(1) * pkin(6) + t1094;
t1218 = t1141 * t1070;
t1149 = cos(qJ(6));
t1108 = t1149 * t1157 * t1141;
t1090 = qJDD(6) + t1108;
t1217 = t1141 * t1090;
t1091 = qJDD(6) - t1108;
t1216 = t1141 * t1091;
t1213 = t1145 * t1092;
t1093 = qJDD(2) - t1110;
t1212 = t1145 * t1093;
t1209 = t1146 * t1069;
t1208 = t1146 * t1070;
t1137 = t1149 ^ 2;
t1195 = t1135 + t1137;
t1082 = t1195 * qJDD(1);
t1207 = t1146 * t1082;
t1205 = t1149 * t1070;
t1204 = t1149 * t1091;
t1202 = t1153 * t1069;
t1201 = t1153 * t1093;
t1200 = t1154 * t1069;
t1199 = t1154 * t1070;
t1198 = t1154 * t1082;
t1194 = t1136 + t1138;
t1193 = t1141 * qJDD(1);
t1192 = t1146 * qJDD(1);
t1123 = t1149 * qJDD(1);
t1191 = t1154 * qJDD(1);
t1189 = pkin(1) * t1145 * t995;
t1187 = t1143 * t856;
t1186 = t1146 * t846;
t1185 = t1146 * t877;
t1184 = t1151 * t856;
t1183 = t1154 * t846;
t1182 = t1154 * t877;
t1181 = -pkin(8) * t1151 - pkin(4);
t1180 = t1141 * t1233;
t1178 = t1146 * t998;
t1177 = t1146 * t1232;
t1176 = t1149 * t1233;
t1174 = t1154 * t998;
t1173 = t1154 * t1232;
t1073 = t1157 * pkin(6) - t1095;
t1038 = -t1141 * g(3) + t1149 * t1073;
t699 = -t1065 * t829 + t1066 * t828;
t608 = -t1065 * t699 + t1066 * t700;
t812 = t1140 * t891 + t1148 * t893;
t579 = t1142 * t688 - t1150 * t634;
t655 = t1143 * t764 + t1151 * t765;
t1036 = t1149 * g(3) + t1141 * t1073;
t948 = t1141 * t1036 + t1149 * t1038;
t950 = t1145 * t1037 + t1153 * t1039;
t1010 = -t1146 * t1094 - t1154 * t1095;
t1167 = t1146 * t1108;
t1166 = t1146 * t1110;
t1165 = t1154 * t1108;
t1164 = t1154 * t1110;
t1085 = -t1146 * t1157 + t1191;
t1163 = -pkin(11) * t1085 - t1146 * g(3);
t607 = t1065 * t700 + t1066 * t699;
t520 = t1142 * t580 - t1150 * t579;
t521 = t1142 * t579 + t1150 * t580;
t814 = t1144 * t892 - t1152 * t894;
t946 = -t1149 * t1036 + t1141 * t1038;
t1009 = t1154 * t1094 - t1146 * t1095;
t902 = t943 + t1225;
t1155 = qJD(6) ^ 2;
t1126 = t1137 * t1157;
t1105 = t1127 - t1156;
t1104 = -t1126 - t1155;
t1103 = t1126 - t1155;
t1102 = -t1156 - t1221;
t1101 = t1156 - t1221;
t1100 = -t1155 - t1222;
t1099 = t1155 - t1222;
t1089 = t1127 - t1221;
t1088 = t1127 + t1221;
t1087 = t1126 - t1222;
t1086 = t1126 + t1222;
t1084 = t1154 * t1157 + t1192;
t1083 = t1194 * qJDD(1);
t1081 = t1124 - 0.2e1 * t1179;
t1079 = t1123 - 0.2e1 * t1180;
t1078 = t1123 - t1180;
t1076 = t1176 + t1193;
t1075 = 0.2e1 * t1176 + t1193;
t1074 = -t1122 - 0.2e1 * t1175;
t1068 = t1153 * t1092;
t1067 = t1149 * t1090;
t1064 = t1194 * t1234;
t1063 = t1195 * t1233;
t1047 = -pkin(11) * t1084 + t1154 * g(3);
t1035 = -t1051 + t1130;
t1034 = -t1050 + t1129;
t1033 = t1049 - t1130;
t1032 = t1048 - t1129;
t1031 = t1153 * t1077 - t1136 * t1234;
t1030 = t1149 * t1076 - t1135 * t1233;
t1029 = -t1145 * t1080 - t1138 * t1234;
t1028 = -t1141 * t1078 - t1137 * t1233;
t1024 = -t1051 - t1130;
t1023 = -t1050 - t1129;
t1022 = -t1145 * t1102 - t1201;
t1021 = -t1145 * t1101 + t1068;
t1020 = t1153 * t1106 - t1213;
t1019 = t1153 * t1105 - t1212;
t1018 = -t1141 * t1100 - t1204;
t1017 = -t1141 * t1099 + t1067;
t1016 = t1149 * t1104 - t1217;
t1015 = t1149 * t1103 - t1216;
t1014 = t1153 * t1102 - t1212;
t1013 = t1145 * t1106 + t1068;
t1012 = t1149 * t1100 - t1216;
t1011 = t1141 * t1104 + t1067;
t1008 = t1154 * t1083 - t1146 * t1088;
t1007 = -t1146 * t1086 + t1198;
t1006 = t1146 * t1083 + t1154 * t1088;
t1005 = t1154 * t1086 + t1207;
t1003 = pkin(13) * t1075 - t1205;
t1002 = -pkin(13) * t1079 - t1218;
t1001 = t1145 * t1074 + t1153 * t1081;
t1000 = -t1141 * t1075 + t1149 * t1079;
t994 = -t1051 + t1049;
t993 = -t1050 + t1048;
t983 = t1154 * t995;
t981 = -t1129 - t1048;
t975 = t1154 * t1189;
t974 = t1146 * t1189;
t970 = t1154 * t1022 - t1146 * t1074;
t969 = t1154 * t1020 - t1146 * t1081;
t968 = t1154 * t1018 + t1146 * t1075;
t967 = t1154 * t1016 - t1146 * t1079;
t966 = t1146 * t1022 + t1154 * t1074;
t965 = t1146 * t1020 + t1154 * t1081;
t964 = t1146 * t1018 - t1154 * t1075;
t963 = t1146 * t1016 + t1154 * t1079;
t960 = -t977 + t1261;
t959 = t976 - t1261;
t958 = (-t1056 * t1152 - t1058 * t1144) * t1133;
t957 = (-t1054 * t1148 + t1057 * t1140) * t1132;
t956 = (-t1056 * t1144 + t1058 * t1152) * t1133;
t955 = (-t1054 * t1140 - t1057 * t1148) * t1132;
t954 = -pkin(12) * t1014 + t1039;
t953 = -pkin(12) * t1013 + t1037;
t952 = -t1049 - t1051;
t951 = -t1048 - t1050;
t941 = -t977 - t1261;
t929 = -pkin(13) * t1086 - t946;
t927 = -t1152 * t1033 + t1245;
t926 = t1144 * t1035 - t1284;
t925 = t1148 * t1032 - t1251;
t924 = -t1140 * t1034 + t1272;
t923 = -t1144 * t1033 - t1236;
t922 = -t1152 * t1035 - t1285;
t921 = t1140 * t1032 + t1243;
t920 = t1148 * t1034 + t1275;
t916 = -t931 + t1117;
t915 = t930 - t1117;
t914 = t1144 * t1024 + t1236;
t913 = -t1140 * t1023 - t1243;
t912 = -t1152 * t1024 + t1245;
t911 = t1148 * t1023 - t1251;
t909 = t1043 + t1161;
t907 = -t1042 + t1159;
t899 = (qJD(7) + t1132) * t1057 + t1168;
t898 = t1154 * t950 - t1209;
t897 = t1154 * t948 - t1208;
t896 = t1146 * t950 + t1200;
t895 = t1146 * t948 + t1199;
t890 = -t931 - t1117;
t888 = (pkin(1) * t1153 + pkin(12)) * t995;
t887 = pkin(6) * t1012 + pkin(13) * t1018 + t1038;
t886 = pkin(6) * t1011 + pkin(13) * t1016 + t1036;
t885 = t1058 * t1224 + t1152 * t1161;
t884 = -t1057 * t1227 - t1148 * t1159;
t883 = -t1058 * t1223 + t1144 * t1161;
t882 = t1057 * t1226 - t1140 * t1159;
t881 = t1056 * t1223 + t1144 * t943;
t880 = t1054 * t1226 - t1140 * t942;
t879 = t1056 * t1224 - t1152 * t943;
t878 = t1054 * t1227 + t1148 * t942;
t876 = -t977 + t976;
t875 = -t1152 * t982 + t1285;
t874 = t1148 * t981 - t1275;
t873 = -t1144 * t982 - t1284;
t872 = t1140 * t981 + t1272;
t867 = -t1261 - t976;
t864 = -t937 + t972;
t863 = t936 - t972;
t862 = (t1143 * t980 - t1151 * t978) * t1121;
t861 = (-t1143 * t978 - t1151 * t980) * t1121;
t859 = -t1145 * t956 + t1153 * t958;
t858 = -t1145 * t955 + t1153 * t957;
t857 = pkin(6) * t946 + pkin(13) * t948;
t855 = -t937 + t936;
t853 = pkin(1) * t903 - t1152 * t995;
t852 = -pkin(1) * t899 + t1148 * t995;
t851 = pkin(1) * t908 + t1144 * t995;
t850 = -pkin(1) * t1262 - t1140 * t995;
t847 = -t976 - t977;
t845 = -t931 + t930;
t841 = -t1117 - t930;
t840 = -t937 - t972;
t839 = -t1145 * t923 + t1153 * t927;
t838 = -t1145 * t922 + t1153 * t926;
t837 = -t1145 * t921 + t1153 * t925;
t836 = -t1145 * t920 + t1153 * t924;
t835 = t1151 * t959 - t1246;
t834 = -t1143 * t960 + t1270;
t833 = t1143 * t959 + t1238;
t832 = t1151 * t960 + t1273;
t831 = -t1143 * t941 - t1238;
t830 = t1151 * t941 - t1246;
t827 = -t1145 * t912 + t1153 * t914;
t826 = -t1145 * t911 + t1153 * t913;
t825 = t1145 * t914 + t1153 * t912;
t824 = t1145 * t913 + t1153 * t911;
t823 = -t972 - t936;
t822 = -t1144 * t908 - t1152 * t903;
t821 = t1144 * t909 - t1152 * t902;
t820 = -t1140 * t907 - t1148 * t900;
t819 = -t1140 * t1262 - t1148 * t899;
t818 = -t1144 * t903 + t1152 * t908;
t817 = -t1144 * t902 - t1152 * t909;
t816 = -t1140 * t900 + t1148 * t907;
t815 = -t1140 * t899 + t1148 * t1262;
t804 = -t1145 * t851 + t1203 * t995;
t803 = -t1145 * t850 - t1206 * t995;
t802 = -t1145 * t853 + t1214 * t995;
t801 = -t1145 * t852 - t1219 * t995;
t800 = t936 + t937;
t797 = (t1065 * t933 + t1066 * t932) * t1120;
t796 = (t1065 * t932 - t1066 * t933) * t1120;
t795 = -t1145 * t883 + t1153 * t885;
t794 = -t1145 * t882 + t1153 * t884;
t793 = -t1145 * t879 + t1153 * t881;
t792 = -t1145 * t878 + t1153 * t880;
t791 = -t1145 * t873 + t1153 * t875;
t790 = -t1145 * t872 + t1153 * t874;
t789 = t1145 * t875 + t1153 * t873;
t788 = t1145 * t874 + t1153 * t872;
t787 = t1151 * t867 - t1273;
t786 = t1143 * t867 + t1270;
t785 = -t930 - t931;
t781 = (t1142 * t940 - t1150 * t938) * t973;
t780 = (t1142 * t938 + t1150 * t940) * t973;
t779 = -t808 - t961;
t774 = (qJD(4) + t1121) * t980 + t1171;
t773 = -pkin(1) * t952 + t814;
t772 = -pkin(1) * t951 + t812;
t771 = t1144 * t861 - t1152 * t862;
t770 = -t1144 * t862 - t1152 * t861;
t769 = t1151 * t808 - t1229 * t980;
t768 = t1143 * t808 + t1228 * t980;
t767 = -t1143 * t807 + t1228 * t978;
t766 = t1151 * t807 + t1229 * t978;
t761 = -t940 * qJD(5) - t1169;
t760 = (-pkin(4) * t1152 + pkin(1)) * t854;
t759 = -t1146 * t908 + t1154 * t827;
t758 = t1146 * t1262 + t1154 * t826;
t757 = t1146 * t827 + t1154 * t908;
t756 = t1146 * t826 - t1154 * t1262;
t755 = t1066 * t915 - t1254;
t754 = -t1065 * t916 + t1276;
t753 = t1065 * t915 + t1252;
t752 = t1066 * t916 + t1277;
t751 = -t1065 * t890 - t1252;
t750 = t1066 * t890 - t1254;
t749 = -t1146 * t903 + t1154 * t791;
t748 = t1146 * t899 + t1154 * t790;
t747 = t1146 * t791 + t1154 * t903;
t746 = t1146 * t790 - t1154 * t899;
t742 = -t783 + t917;
t737 = -t1162 - t918;
t736 = t1066 * t841 - t1277;
t735 = t1065 * t841 + t1276;
t734 = t1066 * t783 - t1231 * t933;
t733 = t1065 * t783 + t1230 * t933;
t732 = t1144 * t833 - t1152 * t835;
t731 = t1144 * t832 - t1152 * t834;
t730 = -t1144 * t835 - t1152 * t833;
t729 = -t1144 * t834 - t1152 * t832;
t728 = t1065 * t1162 - t1230 * t932;
t727 = -t1066 * t1162 - t1231 * t932;
t726 = -pkin(1) * t912 - pkin(12) * t825 + t894;
t725 = -pkin(1) * t911 - pkin(12) * t824 + t893;
t724 = t1144 * t830 - t1152 * t831;
t723 = -t1144 * t831 - t1152 * t830;
t722 = -t1145 * t818 + t1153 * t822;
t721 = -t1145 * t817 + t1153 * t821;
t720 = -t1145 * t816 + t1153 * t820;
t719 = -t1145 * t815 + t1153 * t819;
t718 = t1145 * t821 + t1153 * t817;
t717 = t1145 * t820 + t1153 * t816;
t716 = t1153 * t814 + t1281;
t715 = t1153 * t812 - t1282;
t714 = t1145 * t814 - t1278;
t713 = t1145 * t812 + t1279;
t712 = -pkin(4) * t774 + t1239;
t711 = -pkin(4) * t777 - t1247;
t709 = -t865 - t762;
t706 = (-qJD(5) - t973) * t940 - t1169;
t704 = t1150 * t762 - t1248 * t940;
t703 = -t1142 * t762 - t1240 * t940;
t702 = -t1142 * t761 + t1240 * t938;
t701 = -t1150 * t761 - t1248 * t938;
t698 = -pkin(1) * t873 - pkin(12) * t789 - t892;
t697 = -pkin(1) * t872 - pkin(12) * t788 + t891;
t696 = -t1140 * t796 + t1148 * t797;
t695 = t1140 * t797 + t1148 * t796;
t693 = pkin(4) * t1214 * t854 - t1145 * t760;
t692 = t1154 * t716 - t1244;
t691 = t1154 * t715 - t1244;
t690 = t1146 * t716 + t983;
t689 = t1146 * t715 + t983;
t685 = t1143 * t1160 + t1151 * t781;
t684 = t1143 * t781 - t1151 * t1160;
t683 = -t1145 * t773 + t1278;
t682 = -t1145 * t772 - t1279;
t681 = t1150 * t863 - t1249;
t680 = -t1142 * t864 + t1271;
t679 = -t1142 * t863 - t1241;
t678 = -t1150 * t864 - t1274;
t677 = t1146 * t952 + t1154 * t721;
t676 = t1146 * t951 + t1154 * t720;
t675 = t1146 * t721 - t1154 * t952;
t674 = t1146 * t720 - t1154 * t951;
t673 = t1144 * t786 - t1152 * t787;
t672 = -t1144 * t787 - t1152 * t786;
t671 = -t1143 * t779 - t1151 * t775;
t670 = -t1143 * t777 - t1151 * t774;
t669 = -t1143 * t775 + t1151 * t779;
t668 = -t1143 * t774 + t1151 * t777;
t667 = (t1139 * t1148 - t1259 * t1140) * t784 * pkin(3);
t666 = -t1142 * t840 - t1241;
t665 = t1150 * t840 - t1249;
t664 = -t1145 * t770 + t1153 * t771;
t663 = t1144 * t768 - t1152 * t769;
t662 = t1144 * t766 - t1152 * t767;
t661 = -t1144 * t769 - t1152 * t768;
t660 = -t1144 * t767 - t1152 * t766;
t659 = t1150 * t823 - t1274;
t658 = t1142 * t823 + t1271;
t657 = t1144 * t711 + t1151 * t1237;
t656 = t1143 * t1237 + t1144 * t712;
t652 = t1151 * t704 + t1187;
t651 = t1151 * t702 - t1187;
t650 = t1143 * t704 - t1184;
t649 = t1143 * t702 + t1184;
t648 = -pkin(1) * t817 - pkin(12) * t718;
t647 = -pkin(1) * t816 - pkin(12) * t717;
t646 = -t1140 * t753 + t1148 * t755;
t645 = -t1140 * t752 + t1148 * t754;
t644 = t1140 * t755 + t1148 * t753;
t643 = t1140 * t754 + t1148 * t752;
t642 = pkin(1) * t813 - pkin(12) * t714;
t641 = -pkin(1) * t809 - pkin(12) * t713;
t640 = -t1140 * t750 + t1148 * t751;
t639 = t1140 * t751 + t1148 * t750;
t638 = -t1258 * t1267 - t1253;
t637 = t1188 * t737 + t1253;
t636 = t737 * t1258 - t1255;
t635 = -t1188 * t1267 - t1255;
t631 = -t1065 * t742 - t1066 * t739;
t630 = -t1065 * t1267 + t1066 * t737;
t629 = -t1065 * t739 + t1066 * t742;
t628 = t1065 * t737 + t1066 * t1267;
t627 = -pkin(4) * t847 + t655;
t626 = -t1140 * t735 + t1148 * t736;
t625 = t1140 * t736 + t1148 * t735;
t624 = -t1145 * t730 + t1153 * t732;
t623 = -t1145 * t729 + t1153 * t731;
t622 = -t1140 * t733 + t1148 * t734;
t621 = -t1140 * t727 + t1148 * t728;
t620 = t1140 * t734 + t1148 * t733;
t619 = t1140 * t728 + t1148 * t727;
t618 = (pkin(1) + (t1139 * t1140 + t1259 * t1148) * pkin(3)) * t784;
t617 = -t1145 * t723 + t1153 * t724;
t616 = t1145 * t724 + t1153 * t723;
t615 = -pkin(1) * t777 + t1144 * t1239 - t1152 * t711;
t614 = -pkin(1) * t774 + t1144 * t1247 - t1152 * t712;
t613 = -t1142 * t709 - t1150 * t705;
t612 = -t1142 * t708 + t1150 * t706;
t611 = -t1142 * t705 + t1150 * t709;
t610 = -t1142 * t706 - t1150 * t708;
t609 = -t1145 * t695 + t1153 * t696;
t606 = -t1143 * t705 + t1151 * t681;
t605 = -t1143 * t709 + t1151 * t680;
t604 = t1143 * t681 + t1151 * t705;
t603 = t1143 * t680 + t1151 * t709;
t602 = t1144 * t684 - t1152 * t685;
t601 = -t1144 * t685 - t1152 * t684;
t600 = t1143 * t708 + t1151 * t666;
t599 = t1143 * t666 - t1151 * t708;
t598 = -t1143 * t706 + t1151 * t659;
t597 = t1143 * t659 + t1151 * t706;
t596 = -t1145 * t672 + t1153 * t673;
t595 = t1145 * t673 + t1153 * t672;
t594 = -t1143 * t855 + t1151 * t612;
t593 = t1143 * t612 + t1151 * t855;
t592 = -pkin(10) * t665 + t1242;
t591 = -pkin(10) * t658 + t1250;
t590 = -t1143 * t800 + t1151 * t613;
t589 = t1143 * t613 + t1151 * t800;
t588 = t1144 * t669 - t1152 * t671;
t587 = t1144 * t668 - t1152 * t670;
t586 = -t1144 * t671 - t1152 * t669;
t585 = -t1144 * t670 - t1152 * t668;
t584 = t1146 * t777 + t1154 * t617;
t583 = t1146 * t617 - t1154 * t777;
t582 = -t1145 * t661 + t1153 * t663;
t581 = -t1145 * t660 + t1153 * t662;
t578 = -t785 * t1258 - t607;
t577 = -t1188 * t785 + t608;
t576 = -t1152 * t655 + t1283;
t575 = -t1144 * t655 - t1280;
t574 = t1144 * t650 - t1152 * t652;
t573 = t1144 * t649 - t1152 * t651;
t572 = -t1144 * t652 - t1152 * t650;
t571 = -t1144 * t651 - t1152 * t649;
t570 = -t1145 * t644 + t1153 * t646;
t569 = -t1145 * t643 + t1153 * t645;
t568 = t1146 * t774 + t1154 * t596;
t567 = t1146 * t596 - t1154 * t774;
t566 = -t1145 * t639 + t1153 * t640;
t565 = t1145 * t640 + t1153 * t639;
t564 = -t1140 * t635 + t1148 * t638;
t563 = -t1140 * t637 + t1148 * t636;
t562 = t1144 * t627 + t1280;
t561 = -t1145 * t618 + t1153 * t667;
t560 = -t1145 * t615 + t1153 * t657;
t559 = -t1145 * t614 + t1153 * t656;
t558 = -t1140 * t629 + t1148 * t631;
t557 = -t1140 * t628 + t1148 * t630;
t556 = t1140 * t631 + t1148 * t629;
t555 = t1140 * t630 + t1148 * t628;
t554 = -pkin(1) * t847 - t1152 * t627 + t1283;
t553 = -t1145 * t625 + t1153 * t626;
t552 = t1145 * t626 + t1153 * t625;
t551 = -t1145 * t620 + t1153 * t622;
t550 = -t1145 * t619 + t1153 * t621;
t549 = -pkin(1) * t1267 + t1140 * t638 + t1148 * t635;
t548 = pkin(1) * t737 + t1140 * t636 + t1148 * t637;
t547 = t1146 * t1267 + t1154 * t566;
t546 = t1146 * t566 - t1154 * t1267;
t545 = -pkin(8) * t665 + t580;
t544 = -pkin(1) * t723 - pkin(4) * t830 - pkin(12) * t616 + t765;
t543 = -pkin(8) * t658 + t579;
t542 = -t1140 * t607 + t1148 * t608;
t541 = t1140 * t608 + t1148 * t607;
t540 = t1144 * t604 - t1152 * t606;
t539 = t1144 * t603 - t1152 * t605;
t538 = -t1144 * t606 - t1152 * t604;
t537 = -t1144 * t605 - t1152 * t603;
t536 = -t1145 * t601 + t1153 * t602;
t535 = -t1146 * t737 + t1154 * t553;
t534 = t1146 * t553 + t1154 * t737;
t533 = t1144 * t599 - t1152 * t600;
t532 = -t1144 * t600 - t1152 * t599;
t531 = t1144 * t597 - t1152 * t598;
t530 = -t1144 * t598 - t1152 * t597;
t529 = t1144 * t593 - t1152 * t594;
t528 = -t1144 * t594 - t1152 * t593;
t527 = -pkin(1) * t672 - pkin(4) * t786 - pkin(12) * t595 + t764;
t526 = t1144 * t589 - t1152 * t590;
t525 = -t1144 * t590 - t1152 * t589;
t524 = -t1145 * t586 + t1153 * t588;
t523 = -t1145 * t585 + t1153 * t587;
t522 = t1145 * t588 + t1153 * t586;
t519 = -t1140 * t577 + t1148 * t578;
t518 = t1146 * t847 + t1154 * t524;
t517 = t1146 * t524 - t1154 * t847;
t516 = -t1145 * t575 + t1153 * t576;
t515 = t1145 * t576 + t1153 * t575;
t514 = -t1145 * t572 + t1153 * t574;
t513 = -t1145 * t571 + t1153 * t573;
t512 = -t1146 * t854 + t1154 * t516;
t511 = t1146 * t516 + t1154 * t854;
t510 = -pkin(1) * t785 + t1140 * t578 + t1148 * t577;
t509 = -t1143 * t545 + t1151 * t592;
t508 = -t1143 * t543 + t1151 * t591;
t507 = t1143 * t687 + t1151 * t521;
t506 = t1143 * t521 - t1151 * t687;
t505 = -t1145 * t554 + t1153 * t562;
t504 = -t1145 * t549 + t1153 * t564;
t503 = -t1145 * t548 + t1153 * t563;
t502 = -pkin(1) * t639 - pkin(12) * t565 + (-t1139 * t751 - t1259 * t750) * pkin(3) + t700;
t501 = -t1145 * t556 + t1153 * t558;
t500 = -t1145 * t555 + t1153 * t557;
t499 = t1145 * t558 + t1153 * t556;
t498 = -pkin(10) * t611 - t520;
t497 = t1146 * t785 + t1154 * t501;
t496 = t1146 * t501 - t1154 * t785;
t495 = -pkin(4) * t665 + t1143 * t592 + t1151 * t545;
t494 = -pkin(4) * t658 + t1143 * t591 + t1151 * t543;
t493 = -pkin(1) * t625 - pkin(12) * t552 + (-t1139 * t736 - t1259 * t735) * pkin(3) - t699;
t492 = -t1145 * t541 + t1153 * t542;
t491 = t1145 * t542 + t1153 * t541;
t490 = -t1145 * t538 + t1153 * t540;
t489 = -t1145 * t537 + t1153 * t539;
t488 = -pkin(1) * t586 - pkin(4) * t669 - pkin(12) * t522;
t487 = -t1145 * t532 + t1153 * t533;
t486 = t1145 * t533 + t1153 * t532;
t485 = -t1146 * t784 + t1154 * t492;
t484 = t1146 * t492 + t1154 * t784;
t483 = -t1145 * t530 + t1153 * t531;
t482 = t1145 * t531 + t1153 * t530;
t481 = t1151 * t498 + t611 * t1257;
t480 = -t1145 * t528 + t1153 * t529;
t479 = -pkin(1) * t575 - pkin(4) * t653 - pkin(12) * t515;
t478 = -t1145 * t525 + t1153 * t526;
t477 = t1145 * t526 + t1153 * t525;
t476 = t1146 * t665 + t1154 * t487;
t475 = t1146 * t487 - t1154 * t665;
t474 = t1146 * t658 + t1154 * t483;
t473 = t1146 * t483 - t1154 * t658;
t472 = (-pkin(10) * t1151 + t1257) * t520;
t471 = t1143 * t498 + t1181 * t611;
t470 = -t1145 * t510 + t1153 * t519;
t469 = t1146 * t611 + t1154 * t478;
t468 = t1146 * t478 - t1154 * t611;
t467 = t1144 * t506 - t1152 * t507;
t466 = -t1144 * t507 - t1152 * t506;
t465 = -pkin(1) * t556 - pkin(12) * t499 + (-t1139 * t631 - t1259 * t629) * pkin(3);
t464 = t1144 * t495 - t1152 * t509;
t463 = t1144 * t494 - t1152 * t508;
t462 = -pkin(1) * t665 - t1144 * t509 - t1152 * t495;
t461 = -pkin(1) * t658 - t1144 * t508 - t1152 * t494;
t460 = (-pkin(10) * t1143 + t1181) * t520;
t459 = -pkin(1) * t541 - pkin(12) * t491 + (-t1139 * t608 - t1259 * t607) * pkin(3);
t458 = -pkin(1) * t532 - pkin(4) * t599 + pkin(8) * t708 - pkin(10) * t666 - pkin(12) * t486 - t1250;
t457 = -pkin(1) * t530 - pkin(4) * t597 - pkin(8) * t706 - pkin(10) * t659 - pkin(12) * t482 + t1242;
t456 = t1144 * t471 - t1152 * t481;
t455 = -pkin(1) * t611 - t1144 * t481 - t1152 * t471;
t454 = -pkin(1) * t525 - pkin(4) * t589 - pkin(8) * t800 - pkin(10) * t613 - pkin(12) * t477 - t521;
t453 = -t1145 * t466 + t1153 * t467;
t452 = t1145 * t467 + t1153 * t466;
t451 = t1144 * t460 - t1152 * t472;
t450 = -t1145 * t462 + t1153 * t464;
t449 = -t1145 * t461 + t1153 * t463;
t448 = t1146 * t520 + t1154 * t453;
t447 = t1146 * t453 - t1154 * t520;
t446 = -pkin(1) * t520 - t1144 * t472 - t1152 * t460;
t445 = -t1145 * t455 + t1153 * t456;
t444 = -pkin(1) * t466 - pkin(4) * t506 + pkin(8) * t687 - pkin(10) * t521 - pkin(12) * t452;
t443 = -t1145 * t446 + t1153 * t451;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1084, -t1085, 0, t1010, 0, 0, 0, 0, 0, 0, t969, t970, t1008, t898, 0, 0, 0, 0, 0, 0, t749, t759, t677, t692, 0, 0, 0, 0, 0, 0, t568, t584, t518, t512, 0, 0, 0, 0, 0, 0, t474, t476, t469, t448, 0, 0, 0, 0, 0, 0, t967, t968, t1007, t897, 0, 0, 0, 0, 0, 0, t748, t758, t676, t691, 0, 0, 0, 0, 0, 0, t535, t547, t497, t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1085, -t1084, 0, t1009, 0, 0, 0, 0, 0, 0, t965, t966, t1006, t896, 0, 0, 0, 0, 0, 0, t747, t757, t675, t690, 0, 0, 0, 0, 0, 0, t567, t583, t517, t511, 0, 0, 0, 0, 0, 0, t473, t475, t468, t447, 0, 0, 0, 0, 0, 0, t963, t964, t1005, t895, 0, 0, 0, 0, 0, 0, t746, t756, t674, t689, 0, 0, 0, 0, 0, 0, t534, t546, t496, t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1013, t1014, 0, -t949, 0, 0, 0, 0, 0, 0, t789, t825, t718, t714, 0, 0, 0, 0, 0, 0, t595, t616, t522, t515, 0, 0, 0, 0, 0, 0, t482, t486, t477, t452, 0, 0, 0, 0, 0, 0, t1011, t1012, 0, t946, 0, 0, 0, 0, 0, 0, t788, t824, t717, t713, 0, 0, 0, 0, 0, 0, t552, t565, t499, t491; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1085, 0, -t1084, 0, t1163, -t1047, -t1009, -pkin(11) * t1009, t1154 * t1031 - t1166, t1154 * t1001 - t1146 * t1089, t1154 * t1021 + t1145 * t1192, t1154 * t1029 + t1166, t1154 * t1019 + t1124 * t1146, t1146 * qJDD(2) + t1154 * t1064, -pkin(11) * t965 - t1145 * t1200 - t1146 * t953, -pkin(11) * t966 - t1146 * t954 - t1153 * t1200, -pkin(11) * t1006 + t1154 * t949, -pkin(11) * t896 - t1146 * t1256, t1154 * t795 - t1177, -t1146 * t994 + t1154 * t722, -t1146 * t909 + t1154 * t838, t1154 * t793 + t1177, t1146 * t902 + t1154 * t839, t1146 * t1190 + t1154 * t859, -pkin(11) * t747 - t1146 * t698 + t1154 * t802, -pkin(11) * t757 - t1146 * t726 + t1154 * t804, -pkin(11) * t675 - t1146 * t648 + t1154 * t683, -pkin(11) * t690 - t1146 * t642 - t975, t1154 * t582 + t1185, -t1146 * t876 + t1154 * t523, -t1146 * t779 + t1154 * t623, t1154 * t581 - t1185, -t1146 * t775 + t1154 * t624, t1146 * t1118 + t1154 * t664, -pkin(11) * t567 - t1146 * t527 + t1154 * t559, -pkin(11) * t583 - t1146 * t544 + t1154 * t560, -pkin(11) * t517 - t1146 * t488 + t1154 * t505, -pkin(11) * t511 - t1146 * t479 + t1154 * t693, -t1146 * t703 + t1154 * t514, -t1146 * t610 + t1154 * t480, -t1146 * t678 + t1154 * t489, -t1146 * t701 + t1154 * t513, -t1146 * t679 + t1154 * t490, -t1146 * t780 + t1154 * t536, -pkin(11) * t473 - t1146 * t457 + t1154 * t449, -pkin(11) * t475 - t1146 * t458 + t1154 * t450, -pkin(11) * t468 - t1146 * t454 + t1154 * t445, -pkin(11) * t447 - t1146 * t444 + t1154 * t443, t1154 * t1030 - t1167, t1154 * t1000 - t1146 * t1087, t1154 * t1017 + t1141 * t1192, t1154 * t1028 + t1167, t1154 * t1015 + t1123 * t1146, t1146 * qJDD(6) + t1154 * t1063, -pkin(11) * t963 + t1154 * t1002 - t1146 * t886, -pkin(11) * t964 + t1154 * t1003 - t1146 * t887, -pkin(11) * t1005 - pkin(13) * t1207 + t1154 * t929, -pkin(11) * t895 - pkin(13) * t1199 - t1146 * t857, t1154 * t794 + t1178, -t1146 * t993 + t1154 * t719, -t1146 * t907 + t1154 * t836, t1154 * t792 - t1178, -t1146 * t900 + t1154 * t837, t1146 * t1131 + t1154 * t858, -pkin(11) * t746 - t1146 * t697 + t1154 * t801, -pkin(11) * t756 - t1146 * t725 + t1154 * t803, -pkin(11) * t674 - t1146 * t647 + t1154 * t682, -pkin(11) * t689 - t1146 * t641 - t975, t1154 * t551 - t1186, -t1146 * t845 + t1154 * t500, -t1146 * t742 + t1154 * t569, t1154 * t550 + t1186, -t1146 * t739 + t1154 * t570, t1146 * t1170 + t1154 * t609, -pkin(11) * t534 - t1146 * t493 + t1154 * t503, -pkin(11) * t546 - t1146 * t502 + t1154 * t504, -pkin(11) * t496 - t1146 * t465 + t1154 * t470, -pkin(11) * t484 - t1146 * t459 + t1154 * t561; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1084, 0, t1085, 0, t1047, t1163, t1010, pkin(11) * t1010, t1146 * t1031 + t1164, t1146 * t1001 + t1154 * t1089, t1146 * t1021 - t1145 * t1191, t1146 * t1029 - t1164, t1146 * t1019 - t1153 * t1191, -t1154 * qJDD(2) + t1146 * t1064, pkin(11) * t969 - t1145 * t1209 + t1154 * t953, pkin(11) * t970 - t1146 * t1202 + t1154 * t954, pkin(11) * t1008 + t1146 * t949, pkin(11) * t898 + t1154 * t1256, t1146 * t795 + t1173, t1146 * t722 + t1154 * t994, t1146 * t838 + t1154 * t909, t1146 * t793 - t1173, t1146 * t839 - t1154 * t902, t1146 * t859 - t1154 * t1190, pkin(11) * t749 + t1146 * t802 + t1154 * t698, pkin(11) * t759 + t1146 * t804 + t1154 * t726, pkin(11) * t677 + t1146 * t683 + t1154 * t648, pkin(11) * t692 + t1154 * t642 - t974, t1146 * t582 - t1182, t1146 * t523 + t1154 * t876, t1146 * t623 + t1154 * t779, t1146 * t581 + t1182, t1146 * t624 + t1154 * t775, -t1154 * t1118 + t1146 * t664, pkin(11) * t568 + t1146 * t559 + t1154 * t527, pkin(11) * t584 + t1146 * t560 + t1154 * t544, pkin(11) * t518 + t1146 * t505 + t1154 * t488, pkin(11) * t512 + t1146 * t693 + t1154 * t479, t1146 * t514 + t1154 * t703, t1146 * t480 + t1154 * t610, t1146 * t489 + t1154 * t678, t1146 * t513 + t1154 * t701, t1146 * t490 + t1154 * t679, t1146 * t536 + t1154 * t780, pkin(11) * t474 + t1146 * t449 + t1154 * t457, pkin(11) * t476 + t1146 * t450 + t1154 * t458, pkin(11) * t469 + t1146 * t445 + t1154 * t454, pkin(11) * t448 + t1146 * t443 + t1154 * t444, t1146 * t1030 + t1165, t1146 * t1000 + t1154 * t1087, t1146 * t1017 - t1141 * t1191, t1146 * t1028 - t1165, t1146 * t1015 - t1149 * t1191, -t1154 * qJDD(6) + t1146 * t1063, pkin(11) * t967 + t1146 * t1002 + t1154 * t886, pkin(11) * t968 + t1146 * t1003 + t1154 * t887, pkin(11) * t1007 + pkin(13) * t1198 + t1146 * t929, pkin(11) * t897 - pkin(13) * t1208 + t1154 * t857, t1146 * t794 - t1174, t1146 * t719 + t1154 * t993, t1146 * t836 + t1154 * t907, t1146 * t792 + t1174, t1146 * t837 + t1154 * t900, -t1154 * t1131 + t1146 * t858, pkin(11) * t748 + t1146 * t801 + t1154 * t697, pkin(11) * t758 + t1146 * t803 + t1154 * t725, pkin(11) * t676 + t1146 * t682 + t1154 * t647, pkin(11) * t691 + t1154 * t641 - t974, t1146 * t551 + t1183, t1146 * t500 + t1154 * t845, t1146 * t569 + t1154 * t742, t1146 * t550 - t1183, t1146 * t570 + t1154 * t739, t1146 * t609 - t1154 * t1170, pkin(11) * t535 + t1146 * t503 + t1154 * t493, pkin(11) * t547 + t1146 * t504 + t1154 * t502, pkin(11) * t497 + t1146 * t470 + t1154 * t465, pkin(11) * t485 + t1146 * t561 + t1154 * t459; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1094, t1095, 0, 0, (t1077 + t1175) * t1145, -t1153 * t1074 + t1145 * t1081, t1153 * t1101 + t1213, t1269 * t1153, t1145 * t1105 + t1201, 0, pkin(12) * t1081 + t1202, pkin(12) * t1074 - t1145 * t1069, pkin(12) * t1088 + t950, pkin(12) * t1069, t1145 * t885 + t1153 * t883, t1145 * t822 + t1153 * t818, t1145 * t926 + t1153 * t922, t1145 * t881 + t1153 * t879, t1145 * t927 + t1153 * t923, t1145 * t958 + t1153 * t956, pkin(12) * t903 + t1153 * t853 + t1215 * t995, pkin(12) * t908 + t1153 * t851 + t1210 * t995, -pkin(12) * t952 + t1153 * t773 + t1281, t888, t1145 * t663 + t1153 * t661, t1145 * t587 + t1153 * t585, t1145 * t731 + t1153 * t729, t1145 * t662 + t1153 * t660, t1145 * t732 + t1153 * t730, t1145 * t771 + t1153 * t770, -pkin(12) * t774 + t1145 * t656 + t1153 * t614, -pkin(12) * t777 + t1145 * t657 + t1153 * t615, -pkin(12) * t847 + t1145 * t562 + t1153 * t554, t1153 * t760 + (pkin(4) * t1215 + pkin(12)) * t854, t1145 * t574 + t1153 * t572, t1145 * t529 + t1153 * t528, t1145 * t539 + t1153 * t537, t1145 * t573 + t1153 * t571, t1145 * t540 + t1153 * t538, t1145 * t602 + t1153 * t601, -pkin(12) * t658 + t1145 * t463 + t1153 * t461, -pkin(12) * t665 + t1145 * t464 + t1153 * t462, -pkin(12) * t611 + t1145 * t456 + t1153 * t455, -pkin(12) * t520 + t1145 * t451 + t1153 * t446, (t1076 + t1176) * t1141, t1149 * t1075 + t1141 * t1079, t1149 * t1099 + t1217, (t1078 - t1180) * t1149, t1141 * t1103 + t1204, 0, -pkin(6) * t1079 + t1205, pkin(6) * t1075 - t1218, -pkin(6) * t1086 + t948, -pkin(6) * t1070, t1145 * t884 + t1153 * t882, t1145 * t819 + t1153 * t815, t1145 * t924 + t1153 * t920, t1145 * t880 + t1153 * t878, t1145 * t925 + t1153 * t921, t1145 * t957 + t1153 * t955, -pkin(12) * t899 + t1153 * t852 - t1220 * t995, -pkin(12) * t1262 + t1153 * t850 - t1211 * t995, -pkin(12) * t951 + t1153 * t772 - t1282, t888, t1145 * t622 + t1153 * t620, t1145 * t557 + t1153 * t555, t1145 * t645 + t1153 * t643, t1145 * t621 + t1153 * t619, t1145 * t646 + t1153 * t644, t1145 * t696 + t1153 * t695, pkin(12) * t737 + t1145 * t563 + t1153 * t548, -pkin(12) * t1267 + t1145 * t564 + t1153 * t549, -pkin(12) * t785 + t1145 * t519 + t1153 * t510, pkin(12) * t784 + t1145 * t667 + t1153 * t618;];
tauB_reg = t1;
