% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh1m1DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh1m1DE1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE1_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_jacobiaD_rot_sym_varpar: pkin has to be [23x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:46:16
	% EndTime: 2020-04-14 18:46:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:50:36
	% EndTime: 2020-04-14 19:22:40
	% DurationCPUTime: 1575.69s
	% Computational Cost: add. (43503447->376), mult. (66897732->720), div. (2794967->26), fcn. (41969201->24), ass. (0->342)
	t1215 = pkin(7) ^ 2;
	t1219 = pkin(1) ^ 2;
	t1208 = sin(qJ(2));
	t1213 = cos(pkin(19));
	t1408 = sin(pkin(19));
	t1409 = cos(qJ(2));
	t1251 = -t1208 * t1213 + t1409 * t1408;
	t1391 = pkin(7) * t1251;
	t1416 = -2 * pkin(1);
	t1329 = -t1391 * t1416 + t1219;
	t1172 = t1215 + t1329;
	t1169 = 0.1e1 / t1172;
	t1184 = t1208 * t1408 + t1409 * t1213;
	t1180 = t1184 * qJD(2);
	t1218 = 0.1e1 / pkin(3);
	t1168 = pkin(3) ^ 2 - pkin(8) ^ 2 + t1172;
	t1177 = pkin(1) + t1391;
	t1413 = -pkin(8) - pkin(3);
	t1166 = (pkin(7) - t1413) * (pkin(7) + t1413) + t1329;
	t1412 = -pkin(8) + pkin(3);
	t1167 = (pkin(7) - t1412) * (pkin(7) + t1412) + t1329;
	t1354 = t1167 * t1166;
	t1221 = sqrt(-t1354);
	t1350 = t1184 * t1221;
	t1153 = -pkin(7) * t1350 + t1168 * t1177;
	t1211 = cos(qJ(3));
	t1333 = t1211 * t1153;
	t1390 = pkin(7) * t1184;
	t1154 = t1168 * t1390 + t1177 * t1221;
	t1207 = sin(qJ(3));
	t1342 = t1207 * t1154;
	t1272 = -t1333 + t1342;
	t1170 = 0.1e1 / t1172 ^ 2;
	t1414 = pkin(1) * pkin(7);
	t1327 = t1170 * t1414;
	t1243 = t1272 * t1327;
	t1332 = t1211 * t1154;
	t1343 = t1207 * t1153;
	t1264 = t1343 / 0.2e1 + t1332 / 0.2e1;
	t1270 = 0.2e1 * (t1166 + t1167) * t1414;
	t1155 = t1180 * t1270;
	t1162 = 0.1e1 / t1221;
	t1402 = -t1162 / 0.2e1;
	t1147 = t1155 * t1390 * t1402;
	t1328 = t1177 * t1416;
	t1302 = -t1168 + t1328;
	t1179 = t1251 * qJD(2);
	t1352 = t1179 * t1221;
	t1140 = t1147 + (t1302 * t1180 - t1352) * pkin(7);
	t1337 = t1211 * t1140;
	t1429 = t1180 * t1215;
	t1318 = t1184 * t1429;
	t1288 = pkin(1) * t1318;
	t1401 = t1162 / 0.2e1;
	t1307 = t1177 * t1401;
	t1351 = t1180 * t1221;
	t1141 = t1155 * t1307 - 0.2e1 * t1288 + (t1179 * t1168 - t1351) * pkin(7);
	t1346 = t1207 * t1141;
	t1116 = (t1180 * t1243 + (-t1337 / 0.2e1 + t1346 / 0.2e1 + t1264 * qJD(3)) * t1169) * t1218;
	t1262 = t1333 / 0.2e1 - t1342 / 0.2e1;
	t1271 = t1332 + t1343;
	t1336 = t1211 * t1141;
	t1347 = t1207 * t1140;
	t1117 = (t1271 * t1180 * t1327 + (t1347 / 0.2e1 + t1336 / 0.2e1 + t1262 * qJD(3)) * t1169) * t1218;
	t1203 = pkin(23) + pkin(22);
	t1201 = sin(t1203);
	t1202 = cos(t1203);
	t1111 = t1116 * t1201 - t1117 * t1202;
	t1216 = pkin(5) ^ 2;
	t1353 = t1169 * t1218;
	t1145 = t1264 * t1353;
	t1146 = t1262 * t1353;
	t1134 = t1145 * t1202 + t1146 * t1201;
	t1393 = pkin(4) * t1134;
	t1415 = -2 * pkin(5);
	t1330 = -t1393 * t1415 + t1216;
	t1411 = -pkin(9) - pkin(11);
	t1125 = (pkin(4) - t1411) * (pkin(4) + t1411) + t1330;
	t1410 = -pkin(9) + pkin(11);
	t1126 = (pkin(4) - t1410) * (pkin(4) + t1410) + t1330;
	t1384 = 2 * pkin(5);
	t1269 = pkin(4) * (t1125 + t1126) * t1384;
	t1102 = t1111 * t1269;
	t1217 = pkin(4) ^ 2;
	t1131 = t1217 + t1330;
	t1127 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t1131;
	t1132 = pkin(5) + t1393;
	t1301 = t1132 * t1415 - t1127;
	t1357 = t1126 * t1125;
	t1220 = sqrt(-t1357);
	t1119 = 0.1e1 / t1220;
	t1289 = t1145 * t1201 - t1202 * t1146;
	t1403 = -t1289 / 0.2e1;
	t1311 = t1119 * t1403;
	t1290 = t1202 * t1116 + t1117 * t1201;
	t1428 = t1290 * t1220;
	t1069 = (t1102 * t1311 + t1301 * t1111 - t1428) * pkin(4);
	t1303 = t1289 * t1217 * t1415;
	t1404 = t1119 / 0.2e1;
	t1312 = t1132 * t1404;
	t1433 = t1111 * t1220;
	t1070 = t1102 * t1312 + t1111 * t1303 + (t1127 * t1290 - t1433) * pkin(4);
	t1128 = 0.1e1 / t1131;
	t1214 = 0.1e1 / pkin(11);
	t1129 = 0.1e1 / t1131 ^ 2;
	t1326 = pkin(4) * pkin(5) * t1129;
	t1392 = pkin(4) * t1289;
	t1107 = t1127 * t1392 + t1132 * t1220;
	t1205 = cos(pkin(21));
	t1359 = t1107 * t1205;
	t1106 = t1127 * t1132 - t1220 * t1392;
	t1204 = sin(pkin(21));
	t1362 = t1106 * t1204;
	t1241 = (t1359 + t1362) * t1326;
	t1398 = t1205 / 0.2e1;
	t1400 = t1204 / 0.2e1;
	t1041 = ((t1069 * t1400 + t1070 * t1398) * t1128 + t1111 * t1241) * t1214;
	t1360 = t1107 * t1204;
	t1361 = t1106 * t1205;
	t1242 = (t1360 - t1361) * t1326;
	t1399 = -t1205 / 0.2e1;
	t1042 = ((t1069 * t1399 + t1070 * t1400) * t1128 + t1111 * t1242) * t1214;
	t1355 = t1128 * t1214;
	t1100 = (t1362 / 0.2e1 + t1359 / 0.2e1) * t1355;
	t1096 = t1100 ^ 2;
	t1101 = (-t1361 / 0.2e1 + t1360 / 0.2e1) * t1355;
	t1098 = 0.1e1 / t1101 ^ 2;
	t1089 = t1096 * t1098 + 0.1e1;
	t1087 = 0.1e1 / t1089;
	t1097 = 0.1e1 / t1101;
	t1364 = t1098 * t1100;
	t1025 = (t1041 * t1097 - t1042 * t1364) * t1087;
	t1090 = atan2(t1100, t1101);
	t1085 = sin(t1090);
	t1086 = cos(t1090);
	t1209 = sin(qJ(1));
	t1341 = t1208 * t1207;
	t1317 = t1209 * t1341;
	t1323 = t1409 * t1211;
	t1173 = -t1209 * t1323 + t1317;
	t1212 = cos(qJ(1));
	t1324 = t1409 * t1207;
	t1286 = t1212 * t1324;
	t1340 = t1208 * t1212;
	t1427 = qJD(2) + qJD(3);
	t1157 = t1173 * qJD(1) + t1427 * (-t1211 * t1340 - t1286);
	t1266 = t1208 * t1211 + t1324;
	t1174 = t1266 * t1209;
	t1183 = t1323 - t1341;
	t1158 = -t1427 * t1212 * t1183 + qJD(1) * t1174;
	t1176 = t1266 * t1212;
	t1082 = t1176 * t1086;
	t1418 = t1207 * t1340 - t1212 * t1323;
	t1421 = t1418 * t1085 - t1082;
	t1426 = t1025 * t1421 + t1158 * t1085 + t1086 * t1157;
	t1257 = t1184 * t1271;
	t1156 = t1184 * t1270;
	t1394 = pkin(1) * t1215;
	t1143 = t1156 * t1307 - 0.2e1 * t1184 ^ 2 * t1394 + (t1168 * t1251 - t1350) * pkin(7);
	t1334 = t1211 * t1143;
	t1281 = t1156 * t1402 - t1168;
	t1142 = (-t1251 * t1221 + (t1281 + t1328) * t1184) * pkin(7);
	t1345 = t1207 * t1142;
	t1265 = t1345 / 0.2e1 + t1334 / 0.2e1;
	t1123 = (-t1265 * t1169 - t1257 * t1327) * t1218;
	t1335 = t1211 * t1142;
	t1344 = t1207 * t1143;
	t1263 = -t1335 / 0.2e1 + t1344 / 0.2e1;
	t1124 = (t1263 * t1169 + t1184 * t1243) * t1218;
	t1113 = t1123 * t1202 + t1124 * t1201;
	t1103 = t1113 * t1269;
	t1114 = -t1123 * t1201 + t1124 * t1202;
	t1071 = (t1103 * t1311 + t1301 * t1113 - t1114 * t1220) * pkin(4);
	t1072 = t1103 * t1312 + t1113 * t1303 + (-t1113 * t1220 + t1114 * t1127) * pkin(4);
	t1044 = ((t1071 * t1400 + t1072 * t1398) * t1128 + t1113 * t1241) * t1214;
	t1045 = ((t1071 * t1399 + t1072 * t1400) * t1128 + t1113 * t1242) * t1214;
	t1026 = (t1044 * t1097 - t1045 * t1364) * t1087;
	t1396 = t1026 + 0.1e1;
	t1112 = t1289 * t1269;
	t1094 = (t1112 * t1311 - t1134 * t1220 + t1289 * t1301) * pkin(4);
	t1095 = t1112 * t1312 + t1289 * t1303 + (t1127 * t1134 - t1220 * t1289) * pkin(4);
	t1067 = ((t1094 * t1400 + t1095 * t1398) * t1128 + t1289 * t1241) * t1214;
	t1068 = ((t1094 * t1399 + t1095 * t1400) * t1128 + t1289 * t1242) * t1214;
	t1028 = (t1067 * t1097 - t1068 * t1364) * t1087;
	t1395 = t1028 + 0.1e1;
	t1206 = sin(qJ(4));
	t1210 = cos(qJ(4));
	t1339 = t1209 * t1210;
	t1367 = t1086 * t1418;
	t1422 = -t1176 * t1085 - t1367;
	t1052 = t1206 * t1422 - t1339;
	t1046 = t1052 ^ 2;
	t1349 = t1206 * t1209;
	t1053 = t1210 * t1422 + t1349;
	t1048 = 0.1e1 / t1053 ^ 2;
	t1033 = t1046 * t1048 + 0.1e1;
	t1031 = 0.1e1 / t1033;
	t1047 = 0.1e1 / t1053;
	t1373 = t1048 * t1052;
	t1258 = -t1206 * t1047 + t1210 * t1373;
	t1348 = t1206 * t1212;
	t1383 = qJD(4) * t1052;
	t999 = qJD(1) * t1348 + t1210 * t1426 - t1383;
	t1385 = t1047 * t1048 * t999;
	t1298 = 0.2e1 * t1052 * t1385;
	t1338 = t1210 * t1212;
	t998 = -qJD(1) * t1338 + t1053 * qJD(4) + t1206 * t1426;
	t1387 = 0.2e1 / t1033 ^ 2 * (-t1046 * t1385 + t998 * t1373);
	t1434 = (t1258 * t1387 + ((qJD(4) * t1047 + t1298) * t1210 + (-t1210 * t998 + (-t999 + t1383) * t1206) * t1048) * t1031) * t1421;
	t1165 = t1427 * t1266;
	t1273 = -t1085 * t1266 + t1086 * t1183;
	t1306 = t1273 * t1025 - t1165 * t1085;
	t1056 = t1085 * t1173 - t1174 * t1086;
	t1039 = atan2(t1056, -t1273);
	t1034 = sin(t1039);
	t1035 = cos(t1039);
	t1019 = t1034 * t1056 - t1035 * t1273;
	t1017 = 0.1e1 / t1019 ^ 2;
	t1381 = t1017 * t1421;
	t1016 = 0.1e1 / t1019;
	t1252 = t1427 * t1409;
	t1316 = qJD(3) * t1341;
	t1164 = -qJD(2) * t1341 + t1252 * t1211 - t1316;
	t1084 = t1183 * t1085;
	t1420 = t1086 * t1266 + t1084;
	t1004 = t1025 * t1420 + t1085 * t1164 + t1086 * t1165;
	t1159 = t1418 * qJD(1) + t1427 * t1174;
	t1160 = qJD(1) * t1286 - t1209 * t1316 - qJD(2) * t1317 + (qJD(1) * t1340 + t1252 * t1209) * t1211;
	t1423 = t1085 * t1174 + t1173 * t1086;
	t1417 = t1025 * t1423 + t1085 * t1159 - t1160 * t1086;
	t1054 = t1056 ^ 2;
	t1063 = 0.1e1 / t1273 ^ 2;
	t1038 = t1054 * t1063 + 0.1e1;
	t1036 = 0.1e1 / t1038;
	t1062 = 0.1e1 / t1273;
	t1371 = t1056 * t1063;
	t1261 = -t1004 * t1371 - t1062 * t1417;
	t992 = t1261 * t1036;
	t987 = (t1056 * t992 + t1004) * t1035 + (t1273 * t992 + t1417) * t1034;
	t1388 = t1016 * t1017 * t987;
	t1299 = 0.2e1 * t1421 * t1388;
	t1431 = t1004 * t1063;
	t1430 = t1042 * t1098;
	t1424 = t1157 * t1085 - t1158 * t1086;
	t1000 = t1025 * t1422 + t1424;
	t1407 = -t1102 / 0.2e1;
	t1406 = -t1103 / 0.2e1;
	t1405 = -t1112 / 0.2e1;
	t1397 = -t1211 / 0.2e1;
	t1055 = t1421 ^ 2;
	t1015 = t1017 * t1055 + 0.1e1;
	t1389 = 0.2e1 / t1015 ^ 2 * (-t1000 * t1381 - t1055 * t1388);
	t1382 = t1062 * t1431;
	t1386 = 0.2e1 / t1038 ^ 2 * (t1054 * t1382 + t1371 * t1417);
	t1380 = t1025 * t1173;
	t1379 = t1025 * t1174;
	t1378 = t1025 * t1176;
	t1374 = t1097 * t1430;
	t1377 = 0.2e1 * (t1041 * t1364 - t1096 * t1374) / t1089 ^ 2;
	t1376 = t1034 * t1421;
	t1375 = t1035 * t1421;
	t1372 = t1056 * t1062;
	t1363 = t1102 * t1119 / t1357;
	t1358 = t1111 * t1217;
	t1356 = t1128 * t1205;
	t1322 = t1062 * t1386;
	t1321 = t1217 * t1384;
	t1320 = t1216 * t1358;
	t1319 = t1162 / t1354 * t1156 * t1155;
	t1313 = t1363 / 0.4e1;
	t1310 = t1128 * t1400;
	t1309 = -t1356 / 0.2e1;
	t1308 = t1356 / 0.2e1;
	t1305 = t1025 * t1367 - t1424;
	t1304 = 0.4e1 * pkin(5) * t1358;
	t1300 = t1087 * t1326;
	t1005 = t1396 * t1423;
	t1008 = t1396 * t1420;
	t1260 = -t1005 * t1062 - t1008 * t1371;
	t994 = t1260 * t1036;
	t1297 = t1273 * t994 + t1005;
	t1296 = -t1056 * t994 - t1008;
	t1009 = t1395 * t1423;
	t1012 = t1395 * t1420;
	t1259 = -t1009 * t1062 - t1012 * t1371;
	t996 = t1259 * t1036;
	t1295 = t1273 * t996 + t1009;
	t1294 = -t1056 * t996 - t1012;
	t1293 = -0.2e1 * t1056 * t1382;
	t1292 = -0.8e1 * t1320;
	t1284 = t1128 * t1129 * t1320;
	t1283 = t1219 * t1318;
	t1282 = -t1289 * t1363 / 0.4e1;
	t1279 = t1205 * t1284;
	t1051 = t1210 * t1423 + t1348;
	t1050 = t1206 * t1423 - t1338;
	t1268 = 0.4e1 * t1204 * t1284;
	t1256 = -0.4e1 * t1106 * t1279;
	t1255 = 0.4e1 * t1107 * t1279;
	t1254 = t1113 * t1268;
	t1253 = t1289 * t1268;
	t1250 = -t1087 * t1430 - t1097 * t1377;
	t1249 = t1031 * t1258;
	t1248 = -t1034 + (t1035 * t1372 + t1034) * t1036;
	t1152 = t1179 * t1270 - 0.8e1 * t1283;
	t1121 = t1147 + (t1319 / 0.4e1 + t1152 * t1401) * t1177 + (-0.4e1 * t1179 * t1184 - 0.2e1 * t1180 * t1251) * t1394 + (t1281 * t1180 - t1352) * pkin(7);
	t1122 = 0.4e1 * t1288 + (t1351 - t1184 * t1319 / 0.4e1 + t1302 * t1179 + (-t1251 * t1155 / 0.2e1 - t1179 * t1156 / 0.2e1 - t1184 * t1152 / 0.2e1) * t1162) * pkin(7);
	t1171 = t1169 * t1170;
	t1104 = (0.4e1 * t1272 * t1171 * t1283 + (t1122 * t1397 + t1207 * t1121 / 0.2e1 + t1265 * qJD(3)) * t1169 + ((-t1335 + t1344) * t1180 + t1272 * t1179 + (t1271 * qJD(3) - t1337 + t1346) * t1184) * t1327) * t1218;
	t1105 = (-0.4e1 * t1219 * t1171 * t1257 * t1429 + (-t1207 * t1122 / 0.2e1 + t1121 * t1397 + t1263 * qJD(3)) * t1169 + ((-t1334 - t1345) * t1180 - t1271 * t1179 + (t1272 * qJD(3) - t1336 - t1347) * t1184) * t1327) * t1218;
	t1092 = t1104 * t1201 + t1105 * t1202;
	t1247 = t1069 * t1113 + t1071 * t1111 + t1092 * t1106;
	t1246 = t1069 * t1289 + t1094 * t1111 + t1106 * t1290;
	t1245 = t1070 * t1113 + t1072 * t1111 + t1092 * t1107;
	t1244 = t1070 * t1289 + t1095 * t1111 + t1107 * t1290;
	t1239 = t1086 * t1164 + t1306;
	t1238 = t1364 * t1377 + (-t1041 * t1098 + 0.2e1 * t1100 * t1374) * t1087;
	t1237 = t1016 * t1426;
	t1093 = t1269 * t1290 + t1289 * t1292;
	t1091 = t1104 * t1202 - t1105 * t1201;
	t1076 = t1159 * t1086;
	t1066 = t1092 * t1269 + t1113 * t1292;
	t1043 = t1289 * t1304 + (t1433 + t1112 * t1282 + t1301 * t1290 + (t1093 * t1403 + t1134 * t1407 + t1290 * t1405) * t1119) * pkin(4);
	t1040 = (t1093 * t1404 + t1112 * t1313) * t1132 + (-t1111 * t1134 - 0.2e1 * t1289 * t1290) * t1321 + (-t1428 - t1111 * t1127 + (t1111 * t1405 + t1289 * t1407) * t1119) * pkin(4);
	t1030 = t1113 * t1304 + (-t1091 * t1220 + t1103 * t1282 + t1301 * t1092 + (t1066 * t1403 + t1114 * t1407 + t1290 * t1406) * t1119) * pkin(4);
	t1029 = (t1066 * t1404 + t1103 * t1313) * t1132 + (-t1092 * t1289 - t1111 * t1114 - t1113 * t1290) * t1321 + (t1091 * t1127 - t1092 * t1220 + (t1111 * t1406 + t1113 * t1407) * t1119) * pkin(4);
	t1013 = 0.1e1 / t1015;
	t1010 = t1395 * t1422;
	t1006 = t1396 * t1422;
	t1002 = -t1056 * t1025 + t1085 * t1160 + t1076;
	t997 = t1248 * t1421;
	t991 = t1250 * t1067 + t1238 * t1068 + (((t1040 * t1308 + t1043 * t1310 + t1106 * t1253 + t1255 * t1289) * t1097 - (t1040 * t1310 + t1043 * t1309 + t1107 * t1253 + t1256 * t1289) * t1364) * t1087 + ((t1244 * t1097 + t1246 * t1364) * t1205 + (t1246 * t1097 - t1244 * t1364) * t1204) * t1300) * t1214;
	t990 = t1250 * t1044 + t1238 * t1045 + (((t1029 * t1308 + t1030 * t1310 + t1106 * t1254 + t1113 * t1255) * t1097 - (t1029 * t1310 + t1030 * t1309 + t1107 * t1254 + t1113 * t1256) * t1364) * t1087 + ((t1245 * t1097 + t1247 * t1364) * t1205 + (t1247 * t1097 - t1245 * t1364) * t1204) * t1300) * t1214;
	t989 = t1295 * t1034 - t1294 * t1035;
	t988 = t1297 * t1034 - t1296 * t1035;
	t986 = t991 * t1084 + (t1266 * t991 + t1164) * t1086 + t1239 * t1028 + t1306;
	t985 = -t1076 + (-t1028 * t1159 - t1173 * t991 - t1395 * t1379) * t1086 + (-t1028 * t1160 - t1174 * t991 + t1395 * t1380 - t1160) * t1085;
	t982 = t990 * t1084 + (t1266 * t990 + t1164) * t1086 + t1239 * t1026 + t1306;
	t981 = -t1076 + (-t1026 * t1159 - t1173 * t990 - t1396 * t1379) * t1086 + (-t1026 * t1160 - t1174 * t990 + t1396 * t1380 - t1160) * t1085;
	t979 = -t1259 * t1386 + (t1012 * t1293 + t1062 * t985 + (-t1004 * t1009 - t1012 * t1417 - t1056 * t986) * t1063) * t1036;
	t978 = -t1260 * t1386 + (t1008 * t1293 + t1062 * t981 + (-t1004 * t1005 - t1008 * t1417 - t1056 * t982) * t1063) * t1036;
	t1 = [t1421 * t1322 + (t1000 * t1062 - t1421 * t1431) * t1036, t978, t979, 0; (-t1056 * t1016 + t997 * t1381) * t1389 + (t997 * t1299 + t1417 * t1016 + (t997 * t1000 - t1056 * t987 - ((-t1036 * t992 * t1372 - t1386) * t1376 + (-t1056 * t1322 - t992 + (-t1261 + t992) * t1036) * t1375 - t1248 * t1000) * t1421) * t1017) * t1013, (-t1006 * t1016 - t988 * t1381) * t1389 + ((t1422 * t990 + t1426) * t1016 - t988 * t1299 + t1026 * t1237 + (-t988 * t1000 - t1006 * t987 + (t1056 * t978 + t1297 * t992 + t1417 * t994 + t982) * t1375 + (-t1004 * t994 + t1273 * t978 + t1296 * t992 - t981) * t1376) * t1017) * t1013, (-t1010 * t1016 - t989 * t1381) * t1389 + ((t1422 * t991 + t1426) * t1016 - t989 * t1299 + t1028 * t1237 + (-t989 * t1000 - t1010 * t987 + (t1056 * t979 + t1295 * t992 + t1417 * t996 + t986) * t1375 + (-t1004 * t996 + t1273 * t979 + t1294 * t992 - t985) * t1376) * t1017) * t1013, 0; (-t1047 * t1050 + t1051 * t1373) * t1387 + ((qJD(1) * t1339 + t1051 * qJD(4) + t1002 * t1206) * t1047 + t1051 * t1298 + (-t1050 * t999 - (-qJD(1) * t1349 - t1050 * qJD(4) + t1002 * t1210) * t1052 - t1051 * t998) * t1048) * t1031, -(-t990 * t1082 + (t1418 * t990 + t1378) * t1085 - t1000 * t1026 + t1305) * t1249 + t1396 * t1434, -(-t991 * t1082 + (t1418 * t991 + t1378) * t1085 - t1000 * t1028 + t1305) * t1249 + t1395 * t1434, -t1387 + (0.2e1 * t998 * t1048 * t1031 + (-0.2e1 * t1031 * t1385 - t1048 * t1387) * t1052) * t1052;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:06
	% EndTime: 2020-04-14 18:43:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:08
	% EndTime: 2020-04-14 18:43:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:05
	% EndTime: 2020-04-14 18:43:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:12
	% EndTime: 2020-04-14 18:43:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:45:52
	% EndTime: 2020-04-14 18:45:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end