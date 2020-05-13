% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m1DE1
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
%   Wie in palh3m1DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m1DE1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE1_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_jacobiaD_rot_sym_varpar: pkin has to be [19x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:21:05
	% EndTime: 2020-04-19 18:21:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:25:27
	% EndTime: 2020-04-19 18:56:46
	% DurationCPUTime: 1588.54s
	% Computational Cost: add. (43503447->369), mult. (66897732->715), div. (2794967->26), fcn. (41969201->24), ass. (0->338)
	t1225 = pkin(5) ^ 2;
	t1229 = pkin(1) ^ 2;
	t1422 = sin(qJ(2));
	t1424 = sin(pkin(16));
	t1425 = cos(qJ(2));
	t1427 = cos(pkin(16));
	t1191 = t1422 * t1424 - t1425 * t1427;
	t1408 = pkin(5) * t1191;
	t1358 = -0.2e1 * pkin(1) * t1408 + t1229;
	t1181 = t1225 + t1358;
	t1178 = 0.1e1 / t1181;
	t1193 = t1422 * t1427 + t1425 * t1424;
	t1188 = t1193 * qJD(2);
	t1228 = 0.1e1 / pkin(2);
	t1177 = pkin(2) ^ 2 - pkin(6) ^ 2 + t1181;
	t1186 = pkin(1) - t1408;
	t1431 = -pkin(6) - pkin(2);
	t1175 = (pkin(5) - t1431) * (pkin(5) + t1431) + t1358;
	t1430 = -pkin(6) + pkin(2);
	t1176 = (pkin(5) - t1430) * (pkin(5) + t1430) + t1358;
	t1366 = t1176 * t1175;
	t1231 = sqrt(-t1366);
	t1362 = t1193 * t1231;
	t1162 = -pkin(5) * t1362 + t1177 * t1186;
	t1421 = sin(qJ(3));
	t1348 = t1162 * t1421;
	t1407 = pkin(5) * t1193;
	t1163 = t1177 * t1407 + t1186 * t1231;
	t1223 = cos(qJ(3));
	t1368 = t1163 * t1223;
	t1262 = t1368 / 0.2e1 + t1348 / 0.2e1;
	t1347 = t1163 * t1421;
	t1369 = t1162 * t1223;
	t1276 = t1347 - t1369;
	t1432 = pkin(1) * pkin(5);
	t1286 = 0.2e1 * (t1175 + t1176) * t1432;
	t1164 = t1188 * t1286;
	t1189 = t1191 * qJD(2);
	t1445 = t1193 * t1225;
	t1334 = t1188 * t1445;
	t1313 = pkin(1) * t1334;
	t1171 = 0.1e1 / t1231;
	t1415 = t1171 / 0.2e1;
	t1327 = t1186 * t1415;
	t1364 = t1188 * t1231;
	t1150 = t1164 * t1327 - 0.2e1 * t1313 + (-t1189 * t1177 - t1364) * pkin(5);
	t1351 = t1150 * t1421;
	t1179 = 0.1e1 / t1181 ^ 2;
	t1357 = t1179 * t1432;
	t1446 = t1164 * t1171;
	t1155 = -t1407 * t1446 / 0.2e1;
	t1406 = t1186 * pkin(1);
	t1324 = t1177 + 0.2e1 * t1406;
	t1363 = t1189 * t1231;
	t1149 = t1155 + (-t1324 * t1188 + t1363) * pkin(5);
	t1373 = t1149 * t1223;
	t1125 = ((-t1351 / 0.2e1 + t1373 / 0.2e1 - t1262 * qJD(3)) * t1178 - t1276 * t1188 * t1357) * t1228;
	t1277 = t1348 + t1368;
	t1259 = t1188 * t1277;
	t1261 = -t1347 / 0.2e1 + t1369 / 0.2e1;
	t1352 = t1149 * t1421;
	t1372 = t1150 * t1223;
	t1126 = (t1259 * t1357 + (t1372 / 0.2e1 + t1352 / 0.2e1 + t1261 * qJD(3)) * t1178) * t1228;
	t1218 = pkin(18) + pkin(19);
	t1216 = sin(t1218);
	t1217 = cos(t1218);
	t1121 = t1125 * t1217 - t1126 * t1216;
	t1226 = pkin(4) ^ 2;
	t1365 = t1178 * t1228;
	t1153 = t1261 * t1365;
	t1154 = t1262 * t1365;
	t1144 = -t1153 * t1217 + t1154 * t1216;
	t1409 = pkin(3) * t1144;
	t1433 = -2 * pkin(4);
	t1360 = -t1409 * t1433 + t1226;
	t1429 = -pkin(8) - pkin(10);
	t1134 = (pkin(3) - t1429) * (pkin(3) + t1429) + t1360;
	t1428 = -pkin(8) + pkin(10);
	t1135 = (pkin(3) - t1428) * (pkin(3) + t1428) + t1360;
	t1403 = 2 * pkin(4);
	t1285 = pkin(3) * (t1134 + t1135) * t1403;
	t1112 = t1121 * t1285;
	t1227 = pkin(3) ^ 2;
	t1140 = t1227 + t1360;
	t1136 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t1140;
	t1141 = pkin(4) + t1409;
	t1323 = t1141 * t1433 - t1136;
	t1376 = t1135 * t1134;
	t1230 = sqrt(-t1376);
	t1128 = 0.1e1 / t1230;
	t1287 = -t1153 * t1216 - t1154 * t1217;
	t1416 = -t1287 / 0.2e1;
	t1331 = t1128 * t1416;
	t1288 = -t1125 * t1216 - t1126 * t1217;
	t1444 = t1288 * t1230;
	t1079 = (t1112 * t1331 + t1323 * t1121 - t1444) * pkin(3);
	t1325 = t1287 * t1227 * t1433;
	t1417 = t1128 / 0.2e1;
	t1332 = t1141 * t1417;
	t1453 = t1121 * t1230;
	t1080 = t1112 * t1332 + t1121 * t1325 + (t1136 * t1288 - t1453) * pkin(3);
	t1137 = 0.1e1 / t1140;
	t1224 = 0.1e1 / pkin(10);
	t1138 = 0.1e1 / t1140 ^ 2;
	t1356 = pkin(3) * pkin(4) * t1138;
	t1410 = pkin(3) * t1287;
	t1117 = t1136 * t1410 + t1141 * t1230;
	t1219 = sin(pkin(17));
	t1379 = t1117 * t1219;
	t1116 = t1136 * t1141 - t1230 * t1410;
	t1220 = cos(pkin(17));
	t1380 = t1116 * t1220;
	t1252 = (t1379 - t1380) * t1356;
	t1413 = -t1220 / 0.2e1;
	t1414 = t1219 / 0.2e1;
	t1051 = ((t1079 * t1413 + t1080 * t1414) * t1137 + t1121 * t1252) * t1224;
	t1378 = t1117 * t1220;
	t1381 = t1116 * t1219;
	t1251 = (t1378 + t1381) * t1356;
	t1412 = t1220 / 0.2e1;
	t1052 = ((t1079 * t1414 + t1080 * t1412) * t1137 + t1121 * t1251) * t1224;
	t1374 = t1137 * t1224;
	t1110 = (-t1380 / 0.2e1 + t1379 / 0.2e1) * t1374;
	t1107 = 0.1e1 / t1110 ^ 2;
	t1111 = (t1378 / 0.2e1 + t1381 / 0.2e1) * t1374;
	t1109 = t1111 ^ 2;
	t1099 = t1107 * t1109 + 0.1e1;
	t1097 = 0.1e1 / t1099;
	t1106 = 0.1e1 / t1110;
	t1383 = t1107 * t1111;
	t1035 = (-t1051 * t1383 + t1052 * t1106) * t1097;
	t1100 = atan2(t1111, t1110);
	t1095 = sin(t1100);
	t1096 = cos(t1100);
	t1305 = t1422 * t1421;
	t1423 = sin(qJ(1));
	t1278 = t1423 * t1305;
	t1309 = t1425 * t1423;
	t1183 = t1223 * t1309 - t1278;
	t1308 = t1425 * t1421;
	t1426 = cos(qJ(1));
	t1281 = t1426 * t1308;
	t1311 = t1426 * t1422;
	t1185 = t1223 * t1311 + t1281;
	t1442 = -qJD(3) - qJD(2);
	t1166 = -t1183 * qJD(1) + t1442 * t1185;
	t1279 = t1423 * t1308;
	t1280 = t1426 * t1305;
	t1307 = t1423 * t1422;
	t1312 = t1426 * t1425;
	t1167 = -qJD(1) * t1279 + (-qJD(1) * t1307 - t1442 * t1312) * t1223 + t1442 * t1280;
	t1184 = t1223 * t1312 - t1280;
	t1436 = t1184 * t1095 + t1185 * t1096;
	t1441 = t1035 * t1436 + t1167 * t1095 - t1166 * t1096;
	t1190 = -t1425 * t1223 + t1305;
	t1173 = t1190 * t1442;
	t1192 = -t1422 * t1223 - t1308;
	t1174 = t1442 * t1192;
	t1451 = -t1095 * t1192 + t1096 * t1190;
	t1459 = t1451 * t1035 + t1174 * t1095 - t1096 * t1173;
	t1165 = t1193 * t1286;
	t1301 = t1165 * t1415 + t1177;
	t1151 = (t1191 * t1231 + (-t1301 - 0.2e1 * t1406) * t1193) * pkin(5);
	t1350 = t1151 * t1421;
	t1411 = pkin(1) * t1225;
	t1152 = t1165 * t1327 - 0.2e1 * t1193 ^ 2 * t1411 + (-t1177 * t1191 - t1362) * pkin(5);
	t1370 = t1152 * t1223;
	t1264 = t1370 / 0.2e1 + t1350 / 0.2e1;
	t1321 = t1193 * t1357;
	t1132 = (t1264 * t1178 + t1277 * t1321) * t1228;
	t1349 = t1152 * t1421;
	t1371 = t1151 * t1223;
	t1263 = -t1349 / 0.2e1 + t1371 / 0.2e1;
	t1133 = (-t1263 * t1178 + t1276 * t1321) * t1228;
	t1124 = -t1132 * t1216 - t1133 * t1217;
	t1113 = t1124 * t1285;
	t1123 = -t1132 * t1217 + t1133 * t1216;
	t1081 = (t1113 * t1331 - t1123 * t1230 + t1323 * t1124) * pkin(3);
	t1082 = t1113 * t1332 + t1124 * t1325 + (t1123 * t1136 - t1124 * t1230) * pkin(3);
	t1054 = ((t1081 * t1413 + t1082 * t1414) * t1137 + t1124 * t1252) * t1224;
	t1055 = ((t1081 * t1414 + t1082 * t1412) * t1137 + t1124 * t1251) * t1224;
	t1036 = (-t1054 * t1383 + t1055 * t1106) * t1097;
	t1458 = t1036 + 0.1e1;
	t1122 = t1287 * t1285;
	t1104 = (t1122 * t1331 - t1144 * t1230 + t1287 * t1323) * pkin(3);
	t1105 = t1122 * t1332 + t1287 * t1325 + (t1136 * t1144 - t1230 * t1287) * pkin(3);
	t1077 = ((t1104 * t1413 + t1105 * t1414) * t1137 + t1287 * t1252) * t1224;
	t1078 = ((t1104 * t1414 + t1105 * t1412) * t1137 + t1287 * t1251) * t1224;
	t1039 = (-t1077 * t1383 + t1078 * t1106) * t1097;
	t1457 = t1039 + 0.1e1;
	t1222 = cos(qJ(4));
	t1221 = sin(qJ(4));
	t1344 = t1423 * t1221;
	t1435 = t1185 * t1095 - t1096 * t1184;
	t1063 = t1222 * t1435 + t1344;
	t1345 = t1426 * t1222;
	t1008 = -qJD(1) * t1345 + t1063 * qJD(4) + t1221 * t1441;
	t1346 = t1426 * t1221;
	t1343 = t1423 * t1222;
	t1274 = -t1221 * t1435 + t1343;
	t1452 = t1274 * qJD(4);
	t1009 = qJD(1) * t1346 + t1222 * t1441 + t1452;
	t1056 = t1274 ^ 2;
	t1058 = 0.1e1 / t1063 ^ 2;
	t1043 = t1056 * t1058 + 0.1e1;
	t1041 = 0.1e1 / t1043;
	t1057 = 0.1e1 / t1063;
	t1392 = t1058 * t1274;
	t1270 = -t1057 * t1221 - t1222 * t1392;
	t1400 = t1009 * t1057 * t1058;
	t1449 = -0.2e1 * t1400;
	t1318 = t1274 * t1449;
	t1402 = 0.2e1 * (-t1008 * t1392 - t1056 * t1400) / t1043 ^ 2;
	t1455 = (t1270 * t1402 + ((qJD(4) * t1057 + t1318) * t1222 + (-t1008 * t1222 + (-t1009 - t1452) * t1221) * t1058) * t1041) * t1436;
	t1265 = -t1035 * t1435 + t1166 * t1095 + t1167 * t1096;
	t1450 = -0.4e1 * t1178 * t1179;
	t1434 = t1190 * t1095 + t1096 * t1192;
	t1014 = t1035 * t1434 - t1095 * t1173 - t1096 * t1174;
	t1073 = 0.1e1 / t1451 ^ 2;
	t1448 = t1014 * t1073;
	t1182 = t1223 * t1307 + t1279;
	t1066 = t1095 * t1183 + t1182 * t1096;
	t1049 = atan2(t1066, -t1451);
	t1044 = sin(t1049);
	t1045 = cos(t1049);
	t1029 = t1044 * t1066 - t1045 * t1451;
	t1027 = 0.1e1 / t1029 ^ 2;
	t1398 = t1027 * t1436;
	t1447 = t1051 * t1107;
	t1437 = t1095 * t1182 - t1183 * t1096;
	t1168 = qJD(1) * t1281 + (qJD(1) * t1311 - t1442 * t1309) * t1223 + t1442 * t1278;
	t1169 = t1184 * qJD(1) + t1442 * t1182;
	t1013 = t1035 * t1437 - t1169 * t1095 - t1168 * t1096;
	t1026 = 0.1e1 / t1029;
	t1072 = 0.1e1 / t1451;
	t1420 = -t1112 / 0.2e1;
	t1419 = -t1113 / 0.2e1;
	t1418 = -t1122 / 0.2e1;
	t1065 = t1436 ^ 2;
	t1025 = t1027 * t1065 + 0.1e1;
	t1064 = t1066 ^ 2;
	t1048 = t1064 * t1073 + 0.1e1;
	t1046 = 0.1e1 / t1048;
	t1390 = t1066 * t1073;
	t1273 = t1013 * t1072 - t1014 * t1390;
	t1002 = t1273 * t1046;
	t1294 = t1044 * t1451 + t1045 * t1066;
	t997 = t1294 * t1002 - t1013 * t1044 + t1014 * t1045;
	t1404 = t1026 * t1027 * t997;
	t1405 = 0.2e1 / t1025 ^ 2 * (-t1065 * t1404 + t1265 * t1398);
	t1399 = t1072 * t1448;
	t1401 = 0.2e1 * (-t1013 * t1390 + t1064 * t1399) / t1048 ^ 2;
	t1393 = t1106 * t1447;
	t1397 = 0.2e1 * (t1052 * t1383 - t1109 * t1393) / t1099 ^ 2;
	t1396 = t1041 * t1058;
	t1395 = t1044 * t1436;
	t1394 = t1045 * t1436;
	t1391 = t1066 * t1072;
	t1382 = t1112 * t1128 / t1376;
	t1377 = t1121 * t1227;
	t1375 = t1137 * t1220;
	t1355 = 0.2e1 * t1404;
	t1354 = t1421 / 0.2e1;
	t1342 = t1035 * t1458;
	t1341 = t1035 * t1457;
	t1340 = t1227 * t1403;
	t1339 = t1058 * t1402;
	t1338 = t1072 * t1401;
	t1337 = t1008 * t1396;
	t1336 = t1226 * t1377;
	t1335 = 0.1e1 / t1366 * t1165 * t1446;
	t1333 = t1382 / 0.4e1;
	t1330 = t1137 * t1414;
	t1329 = -t1375 / 0.2e1;
	t1328 = t1375 / 0.2e1;
	t1326 = 0.4e1 * pkin(4) * t1377;
	t1320 = t1097 * t1356;
	t1319 = t1436 * t1355;
	t1317 = -0.2e1 * t1066 * t1399;
	t1316 = -0.8e1 * t1336;
	t1304 = t1137 * t1138 * t1336;
	t1303 = t1229 * t1334;
	t1302 = -t1287 * t1382 / 0.4e1;
	t1295 = t1220 * t1304;
	t1284 = 0.4e1 * t1219 * t1304;
	t1061 = -t1222 * t1437 + t1346;
	t1275 = t1221 * t1437 + t1345;
	t1015 = t1458 * t1437;
	t1018 = t1458 * t1434;
	t1272 = t1015 * t1072 - t1018 * t1390;
	t1019 = t1457 * t1437;
	t1022 = t1457 * t1434;
	t1271 = t1019 * t1072 - t1022 * t1390;
	t1269 = -0.4e1 * t1116 * t1295;
	t1268 = 0.4e1 * t1117 * t1295;
	t1267 = t1124 * t1284;
	t1266 = t1287 * t1284;
	t1260 = -t1097 * t1447 - t1106 * t1397;
	t1258 = t1041 * t1270;
	t1257 = -t1044 + (t1045 * t1391 + t1044) * t1046;
	t1161 = -t1189 * t1286 - 0.8e1 * t1303;
	t1130 = t1155 + (t1335 / 0.4e1 + t1161 * t1415) * t1186 + (0.2e1 * t1188 * t1191 + 0.4e1 * t1189 * t1193) * t1411 + (-t1301 * t1188 + t1363) * pkin(5);
	t1131 = 0.4e1 * t1313 + (t1364 - t1193 * t1335 / 0.4e1 + t1324 * t1189 + (t1191 * t1164 / 0.2e1 + t1189 * t1165 / 0.2e1 - t1193 * t1161 / 0.2e1) * t1171) * pkin(5);
	t1114 = (-t1276 * t1303 * t1450 + (-t1131 * t1223 / 0.2e1 + t1130 * t1354 + t1264 * qJD(3)) * t1178 + (-t1276 * t1189 - (-t1349 + t1371) * t1188 + (t1277 * qJD(3) + t1351 - t1373) * t1193) * t1357) * t1228;
	t1115 = (-t1229 * t1259 * t1445 * t1450 + (t1130 * t1223 / 0.2e1 + t1131 * t1354 + t1263 * qJD(3)) * t1178 + (-t1277 * t1189 - (-t1350 - t1370) * t1188 + (-t1276 * qJD(3) + t1352 + t1372) * t1193) * t1357) * t1228;
	t1102 = -t1114 * t1217 - t1115 * t1216;
	t1256 = t1079 * t1124 + t1081 * t1121 + t1102 * t1116;
	t1255 = t1079 * t1287 + t1104 * t1121 + t1116 * t1288;
	t1254 = t1080 * t1124 + t1082 * t1121 + t1102 * t1117;
	t1253 = t1080 * t1287 + t1105 * t1121 + t1117 * t1288;
	t1248 = t1383 * t1397 + (-t1052 * t1107 + 0.2e1 * t1111 * t1393) * t1097;
	t1247 = t1026 * t1441;
	t1103 = t1285 * t1288 + t1287 * t1316;
	t1101 = t1114 * t1216 - t1115 * t1217;
	t1087 = t1169 * t1096;
	t1076 = t1102 * t1285 + t1124 * t1316;
	t1053 = t1287 * t1326 + (t1453 + t1122 * t1302 + t1323 * t1288 + (t1103 * t1416 + t1144 * t1420 + t1288 * t1418) * t1128) * pkin(3);
	t1050 = (t1103 * t1417 + t1122 * t1333) * t1141 + (-t1121 * t1144 - 0.2e1 * t1287 * t1288) * t1340 + (-t1121 * t1136 - t1444 + (t1121 * t1418 + t1287 * t1420) * t1128) * pkin(3);
	t1040 = t1124 * t1326 + (-t1101 * t1230 + t1113 * t1302 + t1323 * t1102 + (t1076 * t1416 + t1123 * t1420 + t1288 * t1419) * t1128) * pkin(3);
	t1038 = (t1076 * t1417 + t1113 * t1333) * t1141 + (-t1102 * t1287 - t1121 * t1123 - t1124 * t1288) * t1340 + (t1101 * t1136 - t1102 * t1230 + (t1121 * t1419 + t1124 * t1420) * t1128) * pkin(3);
	t1023 = 0.1e1 / t1025;
	t1020 = t1457 * t1435;
	t1016 = t1458 * t1435;
	t1012 = -t1066 * t1035 - t1095 * t1168 + t1087;
	t1007 = t1257 * t1436;
	t1006 = t1271 * t1046;
	t1004 = t1272 * t1046;
	t1001 = t1260 * t1078 + t1248 * t1077 + (((t1050 * t1328 + t1053 * t1330 + t1116 * t1266 + t1268 * t1287) * t1106 - (t1050 * t1330 + t1053 * t1329 + t1117 * t1266 + t1269 * t1287) * t1383) * t1097 + ((t1253 * t1106 + t1255 * t1383) * t1220 + (t1255 * t1106 - t1253 * t1383) * t1219) * t1320) * t1224;
	t1000 = t1260 * t1055 + t1248 * t1054 + (((t1038 * t1328 + t1040 * t1330 + t1116 * t1267 + t1124 * t1268) * t1106 - (t1038 * t1330 + t1040 * t1329 + t1117 * t1267 + t1124 * t1269) * t1383) * t1097 + ((t1254 * t1106 + t1256 * t1383) * t1220 + (t1256 * t1106 - t1254 * t1383) * t1219) * t1320) * t1224;
	t999 = t1294 * t1006 - t1019 * t1044 + t1022 * t1045;
	t998 = t1294 * t1004 - t1015 * t1044 + t1018 * t1045;
	t996 = t1001 * t1434 + t1039 * t1459 + t1459;
	t995 = -t1087 + (-t1001 * t1183 - t1039 * t1169 + t1182 * t1341) * t1096 + (t1001 * t1182 + t1039 * t1168 + t1183 * t1341 + t1168) * t1095;
	t992 = t1000 * t1434 + t1036 * t1459 + t1459;
	t991 = -t1087 + (-t1000 * t1183 - t1036 * t1169 + t1182 * t1342) * t1096 + (t1000 * t1182 + t1036 * t1168 + t1183 * t1342 + t1168) * t1095;
	t989 = -t1271 * t1401 + (t1022 * t1317 + t1072 * t995 + (t1013 * t1022 + t1014 * t1019 - t1066 * t996) * t1073) * t1046;
	t988 = -t1272 * t1401 + (t1018 * t1317 + t1072 * t991 + (t1013 * t1018 + t1014 * t1015 - t1066 * t992) * t1073) * t1046;
	t1 = [t1436 * t1338 + (-t1072 * t1265 - t1436 * t1448) * t1046, t988, t989, 0; (t1007 * t1398 - t1066 * t1026) * t1405 + (-t1013 * t1026 + (-t1007 * t1265 - t1066 * t997) * t1027 - (((-t1002 * t1046 * t1391 - t1401) * t1044 + (-t1066 * t1338 - t1002 + (t1002 - t1273) * t1046) * t1045) * t1398 - t1355 * t1007 + t1257 * t1027 * t1265) * t1436) * t1023, (-t1016 * t1026 - t998 * t1398) * t1405 + ((t1435 * t1000 + t1441) * t1026 - t998 * t1319 + t1036 * t1247 + (t998 * t1265 - t1016 * t997 + (-t1004 * t1013 + t1066 * t988 + t992 + (t1004 * t1451 - t1015) * t1002) * t1394 + (-t1004 * t1014 + t1451 * t988 - t991 + (-t1004 * t1066 - t1018) * t1002) * t1395) * t1027) * t1023, (-t1020 * t1026 - t999 * t1398) * t1405 + ((t1435 * t1001 + t1441) * t1026 - t999 * t1319 + t1039 * t1247 + (t999 * t1265 - t1020 * t997 + (-t1006 * t1013 + t1066 * t989 + t996 + (t1006 * t1451 - t1019) * t1002) * t1394 + (-t1006 * t1014 + t1451 * t989 - t995 + (-t1006 * t1066 - t1022) * t1002) * t1395) * t1027) * t1023, 0; (-t1274 * t1339 - t1337) * t1061 - (-t1009 * t1396 - t1057 * t1402) * t1275 + ((qJD(1) * t1343 + qJD(4) * t1061 + t1012 * t1221) * t1057 + (-qJD(1) * t1344 + qJD(4) * t1275 + t1012 * t1222) * t1392 + t1061 * t1318) * t1041, -(t1000 * t1436 + t1036 * t1265 + t1265) * t1258 + t1458 * t1455, -(t1001 * t1436 + t1039 * t1265 + t1265) * t1258 + t1457 * t1455, -t1402 - (0.2e1 * t1337 - (t1041 * t1449 - t1339) * t1274) * t1274;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:20
	% EndTime: 2020-04-19 18:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:22
	% EndTime: 2020-04-19 18:18:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:20:45
	% EndTime: 2020-04-19 18:20:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end