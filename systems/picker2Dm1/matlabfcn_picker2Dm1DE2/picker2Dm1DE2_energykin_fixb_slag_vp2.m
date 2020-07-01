% Calculate kinetic energy for
% picker2Dm1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm1DE2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1DE2_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:23:17
% EndTime: 2020-05-10 20:24:16
% DurationCPUTime: 48.98s
% Computational Cost: add. (785024->732), mult. (2183357->1257), div. (27852->42), fcn. (593310->35), ass. (0->554)
t1082 = (pkin(3) ^ 2);
t1374 = -2 * t1082;
t1091 = (pkin(7) ^ 2);
t1067 = -2 * t1091;
t1034 = sin(qJ(2));
t1035 = sin(qJ(1));
t1277 = t1034 * t1035;
t1346 = pkin(1) * pkin(3);
t1170 = t1277 * t1346;
t1146 = qJD(2) * t1170;
t1037 = cos(qJ(2));
t1038 = cos(qJ(1));
t1273 = t1037 * t1038;
t1188 = qJD(1) * t1273;
t1147 = t1188 * t1346;
t1373 = -0.6e1 * t1147 + 0.6e1 * t1146;
t1372 = -t1147 + t1146;
t1077 = (pkin(4) ^ 2);
t1088 = (pkin(1) ^ 2);
t1086 = t1088 ^ 2;
t1275 = t1035 * t1037;
t1217 = pkin(3) * t1275;
t1357 = -4 * pkin(7);
t1243 = t1091 - t1077;
t976 = t1088 + t1243;
t1221 = pkin(3) * t976 * t1357;
t1290 = 4 * pkin(1);
t1002 = t1038 ^ 2;
t1247 = -t1082 + t1091;
t998 = t1034 ^ 2;
t1304 = t1082 * t998;
t1352 = 0.2e1 * t1034;
t1262 = pkin(7) * t1352;
t971 = pkin(3) * t1262;
t926 = t1247 + t971 + 0.2e1 * t1304;
t1296 = t926 * t1002;
t987 = pkin(3) * t1034;
t969 = t987 + pkin(7);
t1306 = t1038 * t969;
t1369 = 0.2e1 * pkin(3);
t935 = t971 + t976;
t1168 = pkin(1) * t1217;
t953 = -0.2e1 * t1168;
t901 = t953 + t935;
t1366 = 4 * t1088;
t1003 = t1082 * t1366;
t1266 = t1082 * t1091;
t974 = t1003 - 4 * t1266;
t832 = t974 * t998 + t1034 * t1221 - t1086 - (t1091 - (t1369 + pkin(4)) * pkin(4)) * (t1091 + (t1369 - pkin(4)) * pkin(4)) + (t1067 + (2 * t1077) - (4 * t1082) - 0.4e1 * t1296) * t1088 + (t935 * t1217 - t901 * t1306) * t1290;
t1092 = sqrt(t832);
t1040 = 0.15e2 * t1086;
t1047 = 18 * t1091;
t1048 = -2 * t1077;
t1050 = -6 * t1077;
t1052 = 2 * t1082;
t1090 = t1091 ^ 2;
t1061 = 3 * t1090;
t1099 = pkin(3) * t1082;
t1079 = t1099 ^ 2;
t1100 = t1082 ^ 2;
t1265 = t1088 * t1002;
t1213 = 0.12e2 * t1265;
t1225 = 0.12e2 * t1304;
t1065 = 3 * t1091;
t1255 = 15 * t1088 + t1065;
t1289 = 6 * pkin(1);
t997 = t1034 * t998;
t1298 = t1099 * t997;
t1335 = pkin(7) * t1038;
t1174 = -0.2e1 * t1217;
t1001 = t1038 * t1002;
t1095 = pkin(1) * t1088;
t1279 = t1001 * t1095;
t1228 = pkin(7) * t1279;
t1261 = 0.4e1 * t1335;
t1270 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t1258 = t1088 / 0.3e1 + t1091;
t906 = -0.4e1 / 0.9e1 * t1168 + 0.4e1 / 0.9e1 * t1082 - t1077 / 0.9e1 + t1258;
t1013 = -t1077 / 0.6e1;
t1022 = 0.2e1 / 0.3e1 * t1082;
t1135 = t1091 - t1168;
t920 = t1013 + t1022 + t1135;
t1021 = 0.4e1 / 0.3e1 * t1082;
t1015 = -t1077 / 0.3e1;
t993 = t1088 + t1091;
t1181 = t1015 + t993;
t955 = t1021 + t1181;
t985 = -t1088 / 0.3e1 + t1091;
t841 = 0.4e1 * t1228 + 0.6e1 * t906 * t1265 + t955 * t1270 + (t985 * t1174 + t920 * t1261) * pkin(1);
t989 = pkin(1) * t1038;
t1240 = 0.4e1 * t989;
t1057 = 6 * t1088;
t1157 = -0.4e1 * t1168;
t1017 = -0.2e1 / 0.3e1 * t1077;
t1066 = 2 * t1091;
t1180 = t1017 + t1022 + t1066;
t1179 = t1017 + t993;
t1292 = t1100 + (t1022 + t1179) * t993;
t856 = t955 * t1157 + t1292 + (t1057 + t1180) * t1082;
t914 = t953 + t955;
t973 = t1247 * t1366;
t843 = pkin(7) * t914 * t1240 + t1002 * t973 + t856;
t1271 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t1025 = -t1082 / 0.3e1;
t983 = t1025 + t1091;
t928 = t983 * t953;
t875 = t955 * t1271 + t928;
t1004 = 0.10e2 / 0.3e1 * t1088;
t876 = (t1004 + t1180) * t1082 + t1292;
t1222 = 0.4e1 * t1265;
t1241 = 0.2e1 * t989;
t972 = pkin(7) * t1241;
t995 = -3 * t1088 + t1091;
t927 = t972 + t1222 + t995;
t1186 = 0.8e1 * t1228;
t992 = -3 * t1082 + t1091;
t943 = t992 * t1186;
t1177 = t1082 + t993;
t966 = -t1077 + t1177;
t970 = t989 + pkin(7);
t986 = t993 ^ 2;
t813 = t841 * t1225 + t943 + t875 * t1213 + t1079 + ((-t1077 + t1082 + t1255) * t1100) + (t1040 + ((t1047 + t1050 + 6 * t1082) * t1088) + t1061 + ((t1048 + t1052) * t1091)) * t1082 + (t986 * t966) + (0.8e1 * t927 * t1298 + 0.6e1 * t843 * t987) * t970 + (-t876 * t1217 + t856 * t1335) * t1289;
t1209 = t970 * t1298;
t1167 = -0.8e1 * t1209;
t1233 = -0.4e1 * t1304;
t1309 = t1034 * t970;
t988 = pkin(3) * t1037;
t1215 = t1088 * t988;
t978 = 3 * t1082 + t993;
t1307 = t1035 * t978;
t1058 = 3 * t1088;
t1251 = t1058 + t1091;
t977 = t1082 + t1251;
t942 = t977 * t988;
t1331 = pkin(1) * t1035;
t946 = t988 - t1331;
t855 = -0.2e1 * t1002 * t1215 + t942 + (0.2e1 * t946 * t1335 - t1307) * pkin(1);
t1239 = 0.2e1 * t988;
t1212 = pkin(7) * t1239;
t1124 = -t1035 * t1095 + t1215;
t936 = 0.2e1 * t1124;
t1059 = 2 * t1088;
t991 = t1059 + t1082;
t994 = -t1088 + t1091;
t858 = t994 * t988 + t1002 * t936 + (t1035 * t991 + t1038 * t1212) * pkin(1);
t895 = -pkin(1) * t1307 + t942;
t1189 = t1035 * t1271;
t1347 = 4 * t1095;
t1348 = 0.4e1 * t1086;
t1365 = 8 * t1088;
t968 = pkin(3) * t1348 + t1099 * t1365;
t898 = t1037 * t968 + t1189 * t1347;
t1056 = 0.5e1 * t1086;
t1244 = t1090 + t1100;
t1042 = 10 * t1088;
t1254 = t1042 + t1066;
t1263 = t1091 * t1088;
t923 = (t1254 * t1082) + t1056 + t1244 + (6 * t1263);
t1051 = 5 * t1100;
t1063 = 6 * t1091;
t932 = t1051 + (t1042 + t1063) * t1082 + t986;
t823 = t858 * t1233 + t1002 * t898 + (-0.4e1 * t855 * t1309 + (-t923 + t1186) * t1037) * pkin(3) + (-0.4e1 * t895 * t1335 + (t932 + t1167) * t1035) * pkin(1);
t938 = t987 + t970;
t806 = t1092 * t823 + t813 * t938;
t1371 = t806 ^ 2;
t1311 = 2 * pkin(7);
t1370 = 2 * pkin(1);
t1242 = pkin(7) * t989;
t1211 = 0.6e1 * t1242;
t1018 = -0.3e1 / 0.2e1 * t1077;
t1076 = t1077 ^ 2;
t1259 = t1076 / 0.2e1 - t1100 / 0.2e1;
t1264 = t1091 * t1077;
t1149 = -(3 * t1264) + t1061 + t1259;
t1246 = t1086 + t1090;
t1293 = t1079 + t993 * ((t1018 + t1066) * t1088 - 0.3e1 / 0.2e1 * t1264 + t1246 + t1259);
t1250 = t1066 - t1077;
t1122 = (t1250 * t1088) + t1076 / 0.6e1 + t1246 - t1264;
t1118 = 0.5e1 / 0.6e1 * t1100 + t1122;
t884 = (t1004 + t1250) * t1082 + t1118;
t840 = -0.6e1 * t884 * t1168 + t1293 + (t1040 + ((t1047 - 9 * t1077) * t1088) + t1149) * t1082 + (t1018 + t1255) * t1100;
t1175 = -0.4e1 * t1217;
t1016 = -t1077 / 0.2e1;
t960 = t1016 + t1177;
t1150 = t960 * t1175;
t1139 = pkin(1) * t1150;
t859 = t1139 + ((t1057 + t1250) * t1082) + t1118;
t877 = t960 * t1271 + t928;
t827 = t859 * t1211 + t877 * t1213 + t840 + t943;
t1368 = 0.8e1 * t827;
t1023 = t1082 / 0.3e1;
t1183 = t1077 / 0.3e1 + t1023 + t1066;
t1367 = 0.10e2 / 0.3e1 * t1100 + (-t1088 + t1183) * t1374 + (t1025 + t1181) * t1067;
t1364 = t992 / 0.2e1;
t1070 = 2 * pkin(2);
t1072 = pkin(6) ^ 2;
t1036 = sin(pkin(9));
t1039 = cos(pkin(9));
t934 = t1052 + t935;
t860 = t926 * t1241 + t934 * t969;
t864 = t1038 * t934 + (0.4e1 * t1002 - 0.2e1) * t969 * pkin(1);
t952 = -pkin(1) + t1217;
t907 = -t952 + t1306;
t822 = t1035 * t860 + t1092 * t907 + t864 * t988;
t1074 = pkin(5) ^ 2;
t1031 = sin(pkin(8));
t1032 = cos(pkin(8));
t925 = -t1031 * t1035 - t1032 * t1038;
t1338 = pkin(5) * t925;
t912 = -0.2e1 * pkin(1) * t1338;
t1294 = t1074 + t912;
t889 = -(t1370 + pkin(5)) * pkin(5) + t1294;
t890 = pkin(5) * (t1370 - pkin(5)) + t1294;
t1094 = sqrt(-t889 * t890);
t903 = t912 + 0.2e1 * t1074;
t910 = -pkin(1) + t1338;
t924 = t1031 * t1038 - t1032 * t1035;
t837 = pkin(5) * t903 * t924 - t1094 * t910;
t1321 = t822 * t837;
t965 = t1052 + t1058 + t1243;
t896 = t965 + t971 + t1157;
t1216 = pkin(3) * t1273;
t909 = t1035 * t969 + t1216;
t820 = -t896 * t1306 + t909 * t1092 + (t952 * t1262 + t965 * t1275) * pkin(3) + (-0.2e1 * t1296 + (0.2e1 * t998 - 0.4e1) * t1082 - t976) * pkin(1);
t1299 = t1094 * t924;
t835 = -pkin(5) * t1299 - t903 * t910;
t1324 = t820 * t835;
t1127 = t1324 / 0.4e1 + t1321 / 0.4e1;
t1083 = 0.1e1 / pkin(3);
t1231 = pkin(1) * t1306;
t873 = t1177 + t953 + t971 + 0.2e1 * t1231;
t871 = 0.1e1 / t873;
t1303 = t1083 * t871;
t897 = t1088 + t1294;
t893 = 0.1e1 / t897;
t1305 = 0.1e1 / pkin(5) * t893;
t1158 = t1303 * t1305;
t799 = t1127 * t1158;
t1320 = t835 * t822;
t1323 = t820 * t837;
t1126 = t1323 / 0.4e1 - t1320 / 0.4e1;
t800 = t1126 * t1158;
t795 = -t1036 * t800 + t1039 * t799;
t1337 = pkin(6) * t795;
t792 = t1337 * t1070;
t1295 = t1072 + t792;
t785 = -(t1070 + pkin(6)) * pkin(6) + t1295;
t786 = pkin(6) * (t1070 - pkin(6)) + t1295;
t1325 = t785 * t786;
t1093 = sqrt(-t1325);
t1363 = t1093 / 0.2e1;
t1362 = qJD(1) - qJD(2);
t1014 = -t1077 / 0.4e1;
t1361 = t1014 + t1082 / 0.2e1;
t1360 = t1065 - t1077 - t1082;
t1148 = qJD(1) * t1170;
t1283 = qJD(2) * t1037;
t1199 = t970 * t1283;
t1163 = pkin(3) * t1199;
t1359 = t1148 - t1163;
t1192 = qJD(2) * t1277;
t1119 = t1188 - t1192;
t1358 = -8 * pkin(7);
t1356 = -4 * pkin(1);
t1355 = 0.1e1 / t795;
t1354 = 0.1e1 / t835;
t1353 = 0.1e1 / t795 ^ 2;
t1351 = -0.8e1 * t1038;
t1350 = -0.2e1 * t1038;
t1349 = 0.2e1 * t1072;
t804 = 0.1e1 / t806;
t1345 = t804 / 0.4e1;
t1344 = t820 / 0.4e1;
t1343 = t822 / 0.4e1;
t1342 = t835 / 0.4e1;
t1341 = t871 / 0.2e1;
t1005 = -0.20e2 / 0.3e1 * t1088;
t1064 = 4 * t1091;
t1184 = 0.2e1 / 0.3e1 * t1077 + t1022 + t1064;
t1185 = 0.4e1 / 0.3e1 * t1077 + t1021 + t1067;
t1340 = t1100 / 0.2e1 - (t1005 + t1184) * t1082 / 0.2e1 + 0.3e1 / 0.2e1 * t1086 - t1185 * t1088 / 0.2e1 - t1090 / 0.2e1;
t894 = 0.1e1 / t897 ^ 2;
t916 = t924 * qJD(1);
t1317 = t894 * t916;
t1169 = pkin(1) * t871 * t1317;
t1312 = 0.2e1 * t1372;
t963 = qJD(2) * t1212;
t1210 = t963 + t1312;
t1287 = qJD(1) * t1035;
t1318 = ((qJD(2) * t1216 - t969 * t1287) * t1370 + t1210) / t873 ^ 2;
t1285 = qJD(1) * t1038;
t1193 = t1035 * t1285;
t1156 = t926 * t1193;
t1282 = qJD(2) * t1082;
t1152 = t1034 * t1037 * t1282;
t1143 = 0.4e1 * t1152;
t918 = t963 + t1143;
t1310 = t1002 * t918;
t1322 = ((0.8e1 * t1156 - 0.4e1 * t1310) * t1088 + (t974 * t1352 + t1221) * t1283 + (-0.4e1 * t1210 * t1306 + (0.8e1 * pkin(7) * t1037 ^ 2 * t1282 + 0.4e1 * qJD(1) * t901 * t969) * t1035 + 0.4e1 * (t935 * t1188 + (-t901 * t1273 - t935 * t1277) * qJD(2)) * pkin(3)) * pkin(1)) / t1092;
t1204 = t1322 / 0.2e1;
t1272 = t1038 * t1092;
t1274 = t1035 * t1092;
t1286 = qJD(1) * t1037;
t801 = t907 * t1204 + (t1038 * t860 - t969 * t1274) * qJD(1) + (t1038 * t918 - t926 * t1287) * t1035 * t1370 + ((-t1272 + (-t934 - 0.8e1 * t1231) * t1035) * t1286 + ((-t864 + t1274) * t1034 + (t1272 + (t1311 * t969 + t934) * t1035 + (t1002 * t1370 - pkin(1) + t1335) * t1037 * t1369) * t1037) * qJD(2)) * pkin(3);
t1115 = 0.4e1 * t1119;
t1113 = t1115 * t989;
t1276 = t1034 * t1038;
t1134 = -t1275 + t1276;
t802 = t909 * t1204 + (-t963 * t1038 + (t1035 * t896 + t1272) * qJD(1)) * t969 + (t1143 + 0.4e1 * t1156 - 0.2e1 * t1310) * pkin(1) + ((t1038 * t965 - t1274) * t1286 + (-t1134 * t1092 - t896 * t1273 - t965 * t1277) * qJD(2) + (t1119 * t987 + t952 * t1283) * t1311 + t969 * t1113) * pkin(3);
t1201 = t916 * t1370;
t1319 = (t889 + t890) * pkin(5) * t1201 / t1094;
t1203 = t1319 / 0.2e1;
t917 = t925 * qJD(1);
t1123 = -t917 * t1094 + t924 * t1203;
t814 = (-(t1370 * t910 - t903) * t916 + t1123) * pkin(5);
t1238 = 0.2e1 * t916 * t924;
t1300 = t1094 * t916;
t815 = t910 * t1203 + t1074 * pkin(1) * t1238 + (t903 * t917 + t1300) * pkin(5);
t761 = (-(t1324 / 0.2e1 + t1321 / 0.2e1) * t1169 + (-t1127 * t1318 + (t802 * t1342 + t814 * t1344 + t801 * t837 / 0.4e1 + t815 * t1343) * t871) * t1305) * t1083;
t762 = (-(-t1323 / 0.2e1 + t1320 / 0.2e1) * t1169 + (t1126 * t1318 + (-t802 * t837 / 0.4e1 - t820 * t815 / 0.4e1 + t814 * t1343 + t801 * t1342) * t871) * t1305) * t1083;
t744 = t1036 * t762 + t1039 * t761;
t1339 = pkin(2) * t744;
t797 = t1036 * t799 + t1039 * t800;
t1336 = pkin(6) * t797;
t1334 = -t1092 / 0.4e1;
t1333 = t1092 / 0.4e1;
t1332 = pkin(1) * qJD(1);
t1330 = pkin(1) * t1002;
t1329 = (-t785 - t786) * pkin(6) * t1339 / t1363;
t1328 = t744 * t797;
t1138 = t1099 * t998 * t1199;
t1128 = -0.24e2 * t1138;
t1278 = t1002 * t1095;
t1229 = pkin(7) * t1278;
t1144 = t1229 * t1287;
t1129 = -0.24e2 * t1144;
t1131 = t985 * t1146;
t1132 = t985 * t1147;
t1151 = t1088 * t1193;
t1140 = -0.24e2 * t1151;
t1284 = qJD(2) * t1034;
t1218 = pkin(3) * t1284;
t1141 = t1218 * t1279;
t1200 = qJD(1) * t1298;
t1164 = pkin(1) * t1200;
t1145 = t1035 * t1164;
t1153 = t1271 * t1285;
t1220 = pkin(1) * t1287;
t1176 = pkin(7) * t1220;
t1159 = -0.6e1 * t1176;
t1161 = -0.2e1 * t1193;
t1219 = pkin(1) * t1285;
t1162 = t978 * t1219;
t1171 = t1265 * t1346;
t1191 = qJD(2) * t1276;
t1230 = pkin(3) * t1309;
t1313 = 0.4e1 * t1372 * t955;
t1116 = t1119 * pkin(3);
t1112 = pkin(1) * t983 * t1116 * t1265;
t1314 = t992 * t1129 - 0.24e2 * t1112;
t881 = t1088 * t1116 * t1261;
t933 = pkin(3) * t1283 - t1220;
t964 = -0.2e1 * t1176;
t999 = t1035 ^ 2;
t1327 = ((t970 * t1164 * t1351 + t999 * t1200 * t1365 + t1128 * t1331 + (0.2e1 * (-t1088 * t1218 - t1095 * t1285) * t1002 + t936 * t1161 - t994 * t1218 + (t991 * t1285 + (-qJD(1) * t1275 - t1191) * pkin(7) * t1369) * pkin(1)) * t1233 - 0.8e1 * t858 * t1152 - 0.4e1 * (-t1162 + (0.4e1 * t1037 * t1151 + (-t977 + 0.2e1 * t1265) * t1284) * pkin(3) + ((-t1218 - t1219) * t989 - t946 * t1220) * t1311) * t1230 + t1141 * t1358 + t1129 * t988 + (t1153 * t1347 - t968 * t1284) * t1002 + t898 * t1161 - 0.4e1 * (-t977 * t1218 - t1162) * t1242 + 0.4e1 * t895 * t1176 + t923 * t1218 + t932 * t1219 + 0.4e1 * t1359 * t855) * t1092 + t823 * t1204 + (0.8e1 * (t964 - 0.8e1 * t1151) * t1209 + (-0.12e2 * t1144 + 0.6e1 * (-0.4e1 / 0.9e1 * t1188 + 0.4e1 / 0.9e1 * t1192) * t1171 - 0.12e2 * t906 * t1151 - t881 - 0.4e1 * t920 * t1176 - 0.2e1 * t1132 + 0.2e1 * t1131) * t1225 + 0.24e2 * t841 * t1152 + 0.6e1 * (t973 * t1161 + (t1312 * t1038 - t914 * t1287) * pkin(7) * t1290 + t1313) * t1230 + t875 * t1140 + t1313 * t1211 + t856 * t1159 + t1314 + (-0.8e1 * t1145 + 0.24e2 * t1138) * t927 + t1373 * t876 - 0.6e1 * t1359 * t843) * t938 + t813 * t933) / t1371;
t1000 = t1002 ^ 2;
t1104 = pkin(7) * t1091;
t1249 = t1076 - t1100;
t1125 = 6 * t1090 + t1249 - 6 * t1264;
t1172 = 0.32e2 * t1228;
t1316 = t986 * (-t1082 + t976);
t1173 = 0.16e2 * t1228;
t1281 = t1000 * t1086;
t1227 = 0.8e1 * t1281;
t1260 = 0.6e1 * t1335;
t1256 = 0.4e1 / 0.7e1 * t1091 - t1077 / 0.7e1;
t949 = t1091 + t1082 / 0.4e1 + t1088 / 0.4e1 - t1077 / 0.8e1;
t852 = -0.32e2 / 0.21e2 * t949 * t1168 + 0.5e1 / 0.42e2 * t1100 + (0.16e2 / 0.21e2 * t1088 + t1256) * t1082 + t1086 / 0.7e1 + t1256 * t1088 + t1090 - 0.3e1 / 0.7e1 * t1264 + t1076 / 0.42e2;
t1027 = 0.4e1 / 0.3e1 * t1088;
t951 = t1258 + t1361;
t854 = -0.8e1 / 0.3e1 * t951 * t1168 + 0.5e1 / 0.18e2 * t1100 + (t1027 + t1015) * t1082 + t1090 - t1086 / 0.3e1 + t1076 / 0.18e2 + (t1021 + 0.2e1 / 0.3e1 * t1088 + t1017) * t1091;
t915 = -t1100 / 0.6e1 + t1122;
t1029 = t1088 / 0.2e1;
t1257 = t1029 + t1091;
t919 = -0.2e1 / 0.3e1 * t1168 + t1014 + t1257;
t980 = (t1064 + t1077) * t1088;
t1026 = -0.2e1 / 0.3e1 * t1082;
t984 = t1026 + t1091;
t826 = t984 * t1227 + t919 * t1173 + 0.14e2 * t852 * t1265 + (t994 * t1100) + (-0.10e2 / 0.3e1 * t1086 + (2 * t1090) - t1264 + t980) * t1082 + t915 * t1270 + (t985 * t1150 + t854 * t1260) * pkin(1);
t1062 = 8 * t1091;
t902 = t1095 * t1175 + t1003 + t1348 + ((t1048 + t1062) * t1088);
t913 = -t1088 + t1135 + t1361;
t842 = t1186 + t1002 * t902 + t960 * t995 + (t1174 * t1270 + t913 * t1261) * pkin(1);
t1245 = t1090 - t1086;
t844 = t983 * t1139 - t1079 + (-t1004 - t1243) * t1100 + (t980 + t1100 / 0.6e1 - t1076 / 0.6e1 + t1245) * t1082 + t915 * t1091;
t1049 = -5 * t1077;
t1055 = 0.7e1 * t1086;
t848 = (t1018 + t1065 + (7 * t1088)) * t1100 + (t1055 + ((t1049 + 10 * t1091) * t1088) + t1149) * t1082 + t1293;
t975 = -12 * pkin(7) * t1095 + t1104 * t1290;
t982 = -0.8e1 * t1086 + (12 * t1263);
t866 = t1002 * t982 + t1038 * t975 + t1173 + t1227 + t1246 - (6 * t1263);
t879 = t953 * t1271 + t960 * t992;
t944 = 0.16e2 * (t1244 - 6 * t1266) * t1086;
t990 = -30 * t1077 + 60 * t1091;
t996 = t998 ^ 2;
t798 = t944 * t1000 + t879 * t1172 + 0.24e2 * t844 * t1265 + (t1048 + t1064 + 28 * t1088) * t1079 + (t966 * t1316) + (0.32e2 * t842 * t1298 + t1368 * t987) * t970 + ((t1090 * t1050) + (t1057 * t1125) + (t1066 * t1249) + t990 * t1086 + (28 * t1095 ^ 2) + (4 * t1104 ^ 2) + 0.24e2 * t826 * t998) * t1082 + 0.8e1 * (-t848 * t1217 + t840 * t1335) * pkin(1) + ((t1088 * t990) + 0.16e2 * t866 * t996 + 0.70e2 * t1086 + t1100 + t1125) * t1100;
t1297 = t1100 * t996;
t921 = 0.4e1 / 0.3e1 * t1265 + t972 + t985;
t1165 = -0.24e2 * t921 * t1297;
t1224 = -0.6e1 * t1265;
t1226 = -0.12e2 * t1304;
t1280 = t1001 * t1086;
t1190 = t1035 * t1280;
t947 = 0.7e1 / 0.6e1 * t1082 + t1013 + t1257;
t1182 = t1013 + t1023 + t1091;
t950 = t1027 + t1182;
t874 = -t947 * t1331 + t950 * t988;
t880 = (t1082 * t994) - 0.5e1 / 0.3e1 * t1086 + t1183 * t1088 + t1091 * (t1015 + t983);
t957 = t1088 + t1182;
t979 = t1052 + t994;
t886 = t957 * t988 - t979 * t1331 / 0.2e1;
t828 = t1190 * t1357 + t874 * t1222 + (-0.8e1 / 0.3e1 * t1281 + t880) * t988 + (t1035 * t1340 + t886 * t1261) * pkin(1);
t1214 = t1095 * t988;
t1223 = -0.4e1 * t1265;
t1253 = t1048 + t1374;
t1178 = t1063 + t1253;
t885 = t1051 + ((t1042 + t1178) * t1082) + (t1026 + t1179) * t993;
t1308 = t1035 * t885;
t870 = t1100 + (t1017 + t1026 + t1254) * t1082 + t1056 + (t1178 * t1088) + t1091 * (t1017 + t984);
t863 = t870 * t988;
t954 = 0.8e1 / 0.3e1 * t1082 + t1181;
t956 = t1015 + t1022 + t1251;
t878 = -t954 * t1331 + t956 * t988;
t961 = 0.5e1 / 0.6e1 * t1082 + t1029 + t1013;
t892 = pkin(1) * t1189 + t961 * t1239;
t829 = t1001 * t1214 * t1358 + t892 * t1223 + t863 + (t878 * t1261 - t1308) * pkin(1);
t845 = -pkin(1) * t1308 + t863;
t883 = -(3 * t1100) + (t1005 + t1185) * t1082 + t1184 * t1088 + t1245;
t846 = t1331 * t1367 + t883 * t988;
t1252 = t1049 - 5 * t1082;
t962 = t1360 * t1042;
t847 = t1079 + ((21 * t1088 + t1360) * t1100) + ((t1253 * t1091) + t1061 + 0.35e2 * t1086 + t962) * t1082 + (t1055 + ((t1062 + t1252) * t1088) + (t1091 * (-t1082 + t1243))) * t993;
t937 = 0.4e1 * t1124;
t945 = t988 + 0.2e1 * t1331;
t959 = t1016 + t977;
t849 = t995 * t988 + t1002 * t937 + (t1035 * t959 + t945 * t1335) * t1370;
t857 = 0.7e1 * t1079 + ((35 * t1088 + 15 * t1091 + t1252) * t1100) + (0.21e2 * t1086 + t962 + (9 * t1090) + ((t1050 - 6 * t1082) * t1091)) * t1082 + t1316;
t948 = t1091 + 0.5e1 / 0.2e1 * t1082 + 0.3e1 / 0.2e1 * t1088 + t1016;
t888 = t1331 * t1364 + t948 * t988;
t807 = t888 * t1173 + t849 * t1167 + t846 * t1224 + t828 * t1226 + (-0.6e1 * t829 * t1309 + (0.24e2 * t983 * t1281 - t847) * t1037) * pkin(3) + (-0.6e1 * t845 * t1335 + (t857 + t1165) * t1035) * pkin(1);
t784 = t1092 * t807 + t798 * t938;
t1326 = 0.1e1 / t784 ^ 2 * t1371;
t1315 = t1373 * t884;
t1302 = 0.1e1 / pkin(1) * t893;
t1301 = t1093 * t744;
t1073 = 0.1e1 / pkin(6);
t1269 = t1073 / pkin(2);
t1078 = 0.1e1 / pkin(4);
t1268 = t1078 * t1083;
t1267 = t1078 / pkin(3) ^ 2;
t1098 = pkin(2) ^ 2;
t789 = t1098 + t1295;
t1236 = 0.1e1 / t789 ^ 2 * t1339;
t1235 = pkin(5) * t1317;
t1194 = t1303 / 0.2e1;
t811 = atan2(t822 * t1194, t820 * t1194);
t809 = sin(t811);
t1232 = t809 * t1332;
t1208 = -t1329 / 0.2e1;
t1207 = -t1327 / 0.4e1;
t1206 = t804 * t1344;
t1205 = t804 * t1343;
t1202 = -t1318 / 0.2e1;
t819 = 0.1e1 / t820 ^ 2;
t765 = qJD(1) + 0.2e1 * ((t822 * t1202 + t801 * t1341) / t820 - (t820 * t1202 + t802 * t1341) * t822 * t819) * pkin(3) * t1083 / (t819 * t822 ^ 2 + 0.1e1) * t873;
t1198 = t1035 * t1297;
t1197 = t871 * t1267;
t1196 = t1073 / t789 / 0.2e1;
t1195 = t1305 / 0.2e1;
t1187 = t1268 / 0.2e1;
t790 = t792 + t1349;
t791 = -pkin(2) - t1337;
t747 = -t1093 * t1336 - t790 * t791;
t748 = -t1093 * t791 + t790 * t1336;
t735 = atan2(t748 * t1196, t747 * t1196);
t733 = sin(t735);
t734 = cos(t735);
t810 = cos(t811);
t808 = t810 * t1332;
t764 = pkin(2) * t765 + t808;
t724 = t733 * t1232 - t734 * t764;
t1160 = pkin(3) * t1191;
t1155 = t1100 * t997 * t1283;
t1154 = qJD(1) * t1190;
t745 = t1036 * t761 - t1039 * t762;
t746 = 0.1e1 / t747 ^ 2;
t718 = t765 + 0.2e1 * (((t791 * t1208 + pkin(2) * t1328 * t1349 + (t745 * t790 + t1301) * pkin(6)) * t1196 - t748 * t1236) / t747 - (-t747 * t1236 + (t797 * t1208 - t1093 * t745 + (-0.2e1 * pkin(2) * t791 + t790) * t744) * pkin(6) * t1196) * t748 * t746) * pkin(6) / (t746 * t748 ^ 2 + 0.1e1) * t789;
t1114 = t1119 * t1370;
t1117 = -t870 * t1218 - t885 * t1219;
t1130 = -0.48e2 * t1144;
t1142 = t1218 * t1281;
t755 = (t1165 * t1219 - 0.24e2 * pkin(1) * (-0.8e1 / 0.3e1 * t1151 + t964) * t1198 - 0.96e2 * t921 * t1155 * t1331 + ((-t995 + t1223) * t1218 + 0.2e1 * (pkin(1) * t959 - t1035 * t937 - 0.2e1 * t1278) * t1285 + (-0.2e1 * t1160 + (-0.2e1 * t1035 * t945 + 0.4e1 * t1330) * qJD(1)) * pkin(1) * pkin(7)) * t1167 + (0.8e1 / 0.3e1 * t1142 + (-t950 * t1218 - t947 * t1219) * t1222 - t880 * t1218 + t1219 * t1340 + (0.32e2 / 0.3e1 * t1280 * t988 + t1088 * t874 * t1351) * t1287 + (t957 * t1160 * t1356 + ((0.12e2 * t1002 * t999 - 0.4e1 * t1000) * t1086 + (-0.4e1 * t1035 * t886 - 0.2e1 * t979 * t1330) * pkin(1)) * qJD(1)) * pkin(7)) * t1226 - 0.24e2 * t828 * t1152 - 0.6e1 * ((-0.4e1 * (pkin(1) * t1153 - 0.2e1 * t961 * t1218) * t1002 + 0.8e1 * t892 * t1193) * t1088 + (0.8e1 * t1141 + (-t956 * t1218 - t954 * t1219) * t1240 + (0.24e2 * t1002 * t1214 + t878 * t1356) * t1287) * pkin(7) + t1117) * t1230 + (-t948 * t1218 + t1219 * t1364) * t1173 + t888 * t1130 + (-t883 * t1218 + t1219 * t1367) * t1224 + 0.12e2 * t846 * t1151 - 0.6e1 * t1117 * t1242 + 0.6e1 * t845 * t1176 + t847 * t1218 + t857 * t1219 + (-0.96e2 * t1154 * t988 - 0.24e2 * t1142) * t983 + (0.8e1 * t1145 + t1128) * t849 + 0.6e1 * t1359 * t829) * t1092 + t807 * t1204 + (0.16e2 * (t982 * t1350 - 0.48e2 * t1229 - 0.32e2 * t1280 - t975) * qJD(1) * t1198 + 0.64e2 * t866 * t1155 + 0.32e2 * (-t881 + (t902 * t1350 + (t913 * t1356 - 0.24e2 * t1278) * pkin(7)) * t1287 + (-t1114 * t1270 - t1115 * t1278) * pkin(3)) * t1209 + 0.24e2 * (-0.32e2 * t984 * t1154 + (-0.2e1 / 0.3e1 * t1188 + 0.2e1 / 0.3e1 * t1192) * t1173 * t1346 + t919 * t1130 - 0.28e2 * t852 * t1151 + (-0.8e1 / 0.3e1 * t1188 + 0.8e1 / 0.3e1 * t1192) * t951 * pkin(3) * t1088 * t1260 + t854 * t1159 + 0.4e1 * (-t1132 + t1131) * t960 - 0.64e2 / 0.3e1 * t1119 * t949 * t1171) * t1304 + 0.48e2 * t826 * t1152 - 0.8e1 * t827 * t1148 + 0.8e1 * (t877 * t1140 + (-t960 * pkin(3) * t1113 - t859 * t1287) * pkin(7) * t1289 + t1314 + t1315) * t1230 + t1163 * t1368 - 0.4e1 * t944 * t1001 * t1287 - pkin(3) * t1114 * t1172 * t1271 - 0.96e2 * t879 * t1144 - 0.96e2 * t960 * t1112 - 0.48e2 * t844 * t1151 + 0.8e1 * t1315 * t1242 - 0.8e1 * t840 * t1176 + 0.8e1 * t1372 * t848 + (-0.32e2 * t1145 + 0.96e2 * t1138) * t842) * t938 + t798 * t933;
t728 = t765 + (0.1e1 / t784 * t806 * t1204 - 0.2e1 * (t755 * t804 / 0.2e1 - t784 * t1327 / 0.2e1) * pkin(4) * pkin(3) * t1092 * t1268 * t1326) / (t832 * t1326 + 0.1e1);
t763 = pkin(3) * t765 + t808;
t776 = atan2(t1092 * t1187, t784 * t804 * t1187);
t774 = sin(t776);
t775 = cos(t776);
t742 = -t774 * t1232 + t775 * t763;
t1133 = t1273 + t1277;
t1121 = t784 * t1205 + t820 * t1333;
t1120 = t784 * t1206 + t822 * t1334;
t911 = -pkin(1) * t925 + pkin(5);
t904 = t912 + t1059;
t869 = t1362 * t1133;
t868 = t1362 * t1134;
t838 = pkin(1) * t904 * t924 + t1094 * t911;
t836 = -pkin(1) * t1299 + t904 * t911;
t834 = 0.1e1 / t836 ^ 2;
t833 = 0.1e1 / t835 ^ 2;
t818 = atan2(t837 * t1195, t835 * t1195);
t817 = cos(t818);
t816 = sin(t818);
t782 = (-((-t911 * t1319 / 0.2e1 + t1088 * pkin(5) * t1238 + (t904 * t917 + t1300) * pkin(1)) * t1302 / 0.2e1 - t838 * t1235) / t836 - (t836 * t1235 - (-(-0.2e1 * pkin(5) * t911 - t904) * t916 + t1123) * pkin(1) * t1302 / 0.2e1) * t838 * t834) / (t834 * t838 ^ 2 + 0.1e1) * t897 * t1370;
t781 = qJD(1) + (t815 * t1354 * t1305 + (-t814 * t833 * t1305 - (-t835 * t833 + t1354) * t894 * t1201) * t837) * t897 / (t833 * t837 ^ 2 + 0.1e1) * pkin(5);
t773 = t1121 * t1197;
t772 = t1120 * t1197;
t771 = atan2(t797, t795);
t769 = cos(t771);
t768 = sin(t771);
t759 = atan2(t1269 * t1363, -t795);
t758 = cos(t759);
t757 = sin(t759);
t754 = -t1133 * t773 + t1134 * t772;
t753 = -t1133 * t772 - t1134 * t773;
t752 = 0.1e1 / t754 ^ 2;
t750 = (-t768 * t810 - t769 * t809) * t1332;
t749 = (t768 * t809 - t769 * t810) * t1332;
t743 = t775 * t1232 + t763 * t774;
t741 = atan2(t753, t754);
t739 = cos(t741);
t738 = sin(t741);
t731 = (-t1121 * t1318 + (t802 * t1333 + t820 * t1322 / 0.8e1 + t755 * t1205 + (t822 * t1207 + t801 * t1345) * t784) * t871) * t1267;
t730 = (-t1120 * t1318 + (t755 * t1206 + t801 * t1334 - t822 * t1322 / 0.8e1 + (t820 * t1207 + t802 * t1345) * t784) * t871) * t1267;
t727 = (-t1328 * t1353 + t745 * t1355) / (t1353 * t797 ^ 2 + 0.1e1) + t765;
t726 = pkin(4) * t728 + t742;
t725 = -t734 * t1232 - t733 * t764;
t723 = -t726 * t738 - t739 * t743;
t722 = -t726 * t739 + t738 * t743;
t719 = ((-t1133 * t730 - t1134 * t731 - t772 * t868 + t773 * t869) / t754 - (-t1133 * t731 + t1134 * t730 - t772 * t869 - t773 * t868) * t753 * t752) / (t752 * t753 ^ 2 + 0.1e1) + t728;
t717 = pkin(6) * t718 + t724;
t716 = (-t1355 * t1329 / 0.4e1 + t1353 * t1301 / 0.2e1) / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t1098 * t1353 * t1325 / 0.4e1) * t1269 + t718;
t715 = -t717 * t757 - t725 * t758;
t714 = -t717 * t758 + t725 * t757;
t1 = qJD(2) ^ 2 * Ifges(8,3) / 0.2e1 + m(11) * (t722 ^ 2 + t723 ^ 2) / 0.2e1 + m(4) * (t724 ^ 2 + t725 ^ 2) / 0.2e1 + m(10) * (t714 ^ 2 + t715 ^ 2) / 0.2e1 + m(5) * (t742 ^ 2 + t743 ^ 2) / 0.2e1 + t765 ^ 2 * Ifges(3,3) / 0.2e1 + m(7) * (t749 ^ 2 + t750 ^ 2) / 0.2e1 + t781 ^ 2 * Ifges(9,3) / 0.2e1 + t782 ^ 2 * Ifges(6,3) / 0.2e1 + (t742 * mrSges(5,1) - t743 * mrSges(5,2) + Ifges(5,3) * t728 / 0.2e1) * t728 + (t749 * mrSges(7,1) - t750 * mrSges(7,2) + Ifges(7,3) * t727 / 0.2e1) * t727 + (t722 * mrSges(11,1) - t723 * mrSges(11,2) + Ifges(11,3) * t719 / 0.2e1) * t719 + (t724 * mrSges(4,1) - t725 * mrSges(4,2) + Ifges(4,3) * t718 / 0.2e1) * t718 + (t714 * mrSges(10,1) - t715 * mrSges(10,2) + Ifges(10,3) * t716 / 0.2e1) * t716 + ((-mrSges(9,1) * t817 + mrSges(9,2) * t816) * t781 + (mrSges(3,1) * t810 - mrSges(3,2) * t809) * t765) * t1332 + (Ifges(2,3) / 0.2e1 + (m(9) * (t816 ^ 2 + t817 ^ 2) / 0.2e1 + m(3) * (t809 ^ 2 + t810 ^ 2) / 0.2e1) * t1088) * qJD(1) ^ 2;
T = t1;