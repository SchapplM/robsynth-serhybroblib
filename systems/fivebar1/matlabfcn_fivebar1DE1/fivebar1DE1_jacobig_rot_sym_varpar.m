% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fivebar1DE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% Jg_rot [3x2]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = fivebar1DE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_jacobig_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fivebar1DE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_jacobig_rot_sym_varpar: pkin has to be [5x1] (double)');
Jg_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 1, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:53:04
	% EndTime: 2020-04-27 03:53:14
	% DurationCPUTime: 10.17s
	% Computational Cost: add. (65317->621), mult. (196552->964), div. (464->9), fcn. (23660->6), ass. (0->422)
	t1300 = (pkin(2) ^ 2);
	t1296 = (pkin(3) ^ 2);
	t1302 = pkin(1) ^ 2;
	t1443 = t1296 - t1302;
	t1198 = t1443 * t1300;
	t1533 = 2 * pkin(1);
	t1532 = 2 * pkin(3);
	t1531 = 4 * pkin(3);
	t1259 = sin(qJ(2));
	t1260 = sin(qJ(1));
	t1476 = t1259 * t1260;
	t1187 = pkin(2) * t1476;
	t1350 = pkin(3) * t1187;
	t1182 = -0.2e1 * t1350;
	t1262 = cos(qJ(1));
	t1502 = pkin(2) * t1262;
	t1426 = pkin(1) * t1502;
	t1197 = -0.2e1 * t1426;
	t1441 = (t1300 + t1302);
	t1361 = t1296 + t1441;
	t1192 = pkin(1) - t1502;
	t1261 = cos(qJ(2));
	t1176 = t1192 * t1261;
	t1397 = pkin(3) * t1176;
	t1112 = t1182 + t1197 + t1361 + 0.2e1 * t1397;
	t1530 = -0.2e1 / t1112 ^ 2;
	t1215 = pkin(3) * t1261;
	t1164 = t1192 + t1215;
	t1520 = 0.8e1 * t1164;
	t1529 = 2 * t1296;
	t1528 = 8 * t1300;
	t1229 = t1262 ^ 2;
	t1464 = t1300 * t1229;
	t1425 = pkin(1) * t1215;
	t1195 = 0.2e1 * t1425;
	t1209 = -t1296 / 0.3e1 + t1302;
	t1226 = t1261 ^ 2;
	t1465 = t1296 * t1226;
	t1143 = 0.4e1 / 0.3e1 * t1465 + t1195 + t1209;
	t1289 = pkin(5) ^ 2;
	t1292 = pkin(4) ^ 2;
	t1217 = -t1289 - t1292;
	t1285 = 0.3e1 * t1302;
	t1201 = t1285 + t1217;
	t1486 = t1201 * t1296;
	t1175 = 0.10e2 * t1486;
	t1274 = 7 * t1296;
	t1308 = t1300 ^ 2;
	t1278 = 5 * t1308;
	t1294 = t1296 ^ 2;
	t1307 = pkin(2) * t1300;
	t1297 = t1307 ^ 2;
	t1301 = t1302 ^ 2;
	t1227 = t1229 ^ 2;
	t1480 = t1227 * t1308;
	t1483 = t1217 * t1302;
	t1219 = t1296 + t1302;
	t1212 = t1219 ^ 2;
	t1360 = -t1289 + t1219;
	t1484 = t1212 * (-t1292 + t1360);
	t1527 = 0.7e1 * t1297 + (t1274 + t1201) * t1278 + (t1175 + (21 * t1294) + 0.9e1 * t1301 + 0.6e1 * t1483) * t1300 + t1484 - 0.24e2 * t1143 * t1480;
	t1221 = -(3 * t1296) + t1302;
	t1412 = 0.4e1 * t1465;
	t1526 = t1221 + t1412;
	t1240 = -t1289 / 0.3e1;
	t1525 = t1240 - t1292 / 0.3e1;
	t1288 = t1289 ^ 2;
	t1291 = t1292 ^ 2;
	t1524 = -t1288 / 0.6e1 + t1291 / 0.6e1;
	t1523 = -0.4e1 * pkin(2);
	t1213 = 0.10e2 / 0.3e1 * t1296;
	t1444 = t1294 + t1301;
	t1286 = 0.2e1 * t1302;
	t1448 = t1286 - t1289;
	t1463 = t1302 * t1289;
	t1140 = t1448 * t1296 + t1444 - t1463 - t1524;
	t1322 = t1140 + t1308;
	t1114 = (t1213 + t1448) * t1300 + t1322;
	t1522 = -0.6e1 * t1114;
	t1410 = 0.2e1 * t1464;
	t1442 = -t1300 + t1302;
	t1148 = t1197 + t1410 + t1442;
	t1521 = -0.2e1 * t1148;
	t1225 = t1261 * t1226;
	t1519 = -0.8e1 * t1225;
	t1518 = 0.4e1 * t1226;
	t1517 = -0.8e1 * t1261;
	t1516 = 0.8e1 * t1261;
	t1515 = 4 * t1294;
	t1275 = 6 * t1296;
	t1514 = 2 * t1300;
	t1513 = pkin(1) * pkin(3);
	t1512 = -pkin(4) - pkin(5);
	t1511 = -pkin(4) + pkin(5);
	t1254 = t1296 / 0.3e1;
	t1130 = -0.4e1 / 0.9e1 * t1350 + t1302 + t1300 / 0.3e1 + t1254 + t1292 / 0.9e1 - t1289 / 0.9e1;
	t1238 = -t1289 / 0.6e1;
	t1257 = t1300 / 0.2e1;
	t1451 = t1257 + t1302;
	t1325 = -t1350 + t1451;
	t1138 = t1292 / 0.6e1 + t1238 + t1325;
	t1248 = t1292 / 0.3e1;
	t1163 = t1240 + t1248 + t1361;
	t1469 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
	t1311 = pkin(3) * t1296;
	t1481 = t1225 * t1311;
	t1494 = 4 * pkin(1);
	t1085 = t1209 * t1182 + 0.6e1 * t1130 * t1465 + t1163 * t1469 + (t1138 * t1215 + t1481) * t1494;
	t1338 = -0.4e1 * t1350;
	t1242 = -0.2e1 / 0.3e1 * t1289;
	t1247 = 0.2e1 / 0.3e1 * t1292;
	t1364 = t1242 + t1247 + t1286;
	t1363 = t1242 + t1219;
	t1457 = t1219 * (t1247 + t1363) + t1308;
	t1097 = t1163 * t1338 + (t1275 + t1364) * t1300 + t1457;
	t1136 = t1182 + t1163;
	t1199 = t1442 * t1296;
	t1386 = 0.4e1 * t1425;
	t1086 = t1136 * t1386 + t1199 * t1518 + t1097;
	t1210 = t1302 - t1300 / 0.3e1;
	t1149 = t1210 * t1182;
	t1468 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
	t1111 = t1163 * t1468 + t1149;
	t1113 = (t1213 + t1364) * t1300 + t1457;
	t1147 = t1195 + t1526;
	t1280 = -3 * t1300;
	t1222 = t1280 + t1302;
	t1407 = pkin(1) * t1481;
	t1358 = 0.8e1 * t1407;
	t1167 = t1222 * t1358;
	t1185 = t1292 + t1360;
	t1267 = 15 * t1294;
	t1268 = 15 * t1296;
	t1282 = 0.3e1 * t1301;
	t1193 = pkin(1) + t1215;
	t1228 = t1262 * t1229;
	t1479 = t1228 * t1307;
	t1381 = t1193 * t1479;
	t1344 = -0.8e1 * t1381;
	t1446 = -t1289 + t1292;
	t1362 = t1285 + t1446;
	t1385 = 0.6e1 * t1425;
	t1389 = 0.12e2 * t1465;
	t1395 = pkin(3) * t1476;
	t1487 = t1193 * t1262;
	t1071 = t1147 * t1344 + t1167 + t1111 * t1389 + t1097 * t1385 + t1297 + (t1268 + t1362) * t1308 + t1212 * t1185 + (0.12e2 * t1085 * t1229 + t1362 * t1275 + t1446 * t1286 + t1267 + t1282) * t1300 + 0.6e1 * (-t1086 * t1487 - t1113 * t1395) * pkin(2);
	t1276 = 3 * t1296;
	t1204 = t1276 + t1441;
	t1216 = pkin(2) * t1260;
	t1168 = t1204 * t1216;
	t1499 = pkin(3) * t1259;
	t1173 = t1216 - t1499;
	t1401 = t1296 * t1216;
	t1340 = t1226 * t1401;
	t1279 = 3 * t1300;
	t1202 = t1279 + t1219;
	t1485 = t1202 * t1259;
	t1504 = pkin(1) * t1261;
	t1098 = -0.2e1 * t1340 + t1168 + (0.2e1 * t1173 * t1504 - t1485) * pkin(3);
	t1218 = t1529 + t1300;
	t1474 = t1259 * t1311;
	t1165 = t1401 - t1474;
	t1488 = t1165 * t1226;
	t1099 = t1218 * t1499 + 0.2e1 * t1488 + (-t1443 + t1195) * t1216;
	t1129 = -pkin(3) * t1485 + t1168;
	t1200 = pkin(2) * t1515 + 0.8e1 * t1296 * t1307;
	t1377 = t1311 * t1468;
	t1131 = t1200 * t1260 + 0.4e1 * t1259 * t1377;
	t1273 = 5 * t1294;
	t1439 = t1301 + t1308;
	t1269 = 10 * t1296;
	t1449 = t1269 + t1286;
	t1462 = t1302 * t1296;
	t1146 = t1449 * t1300 + t1273 + t1439 + 0.6e1 * t1462;
	t1154 = t1278 + (t1269 + 0.6e1 * t1302) * t1300 + t1212;
	t1343 = 0.8e1 * t1381;
	t1411 = -0.4e1 * t1464;
	t1075 = t1099 * t1411 + t1131 * t1226 + (-0.4e1 * t1129 * t1504 + (t1154 + t1343) * t1259) * pkin(3) + (0.4e1 * t1098 * t1487 + (-t1146 + t1358) * t1260) * pkin(2);
	t1359 = -t1289 + t1441;
	t1329 = t1296 + t1359;
	t1179 = -t1292 + t1329;
	t1145 = t1197 + t1179;
	t1126 = t1182 + t1145;
	t1445 = t1289 - t1302;
	t1328 = t1226 * t1521 + t1445;
	t1303 = sqrt(0.4e1 * t1198 * t1229 + 0.4e1 * t1179 * t1426 - t1294 - (t1302 + (pkin(2) - t1511) * (pkin(2) + t1511)) * (t1302 + (pkin(2) - t1512) * (pkin(2) + t1512)) + (t1280 + t1292 + t1328) * t1529 + (-t1126 * t1176 + t1145 * t1187) * t1531);
	t1065 = t1071 * t1164 + t1075 * t1303;
	t1510 = 0.1e1 / t1065 / 0.4e1;
	t1509 = -0.1e1 / t1065 ^ 2 / 0.4e1;
	t1508 = 0.4e1 / 0.3e1 * t1300;
	t1507 = -t1303 / 0.4e1;
	t1506 = t1303 / 0.4e1;
	t1505 = pkin(1) * (t1187 - pkin(3));
	t1503 = pkin(1) * t1294;
	t1501 = pkin(3) * t1192;
	t1500 = pkin(3) * t1226;
	t1455 = t1288 / 0.2e1 - t1291 / 0.2e1;
	t1333 = -0.3e1 * t1463 + t1282 + t1455;
	t1243 = -0.3e1 / 0.2e1 * t1289;
	t1452 = t1243 + t1285;
	t1458 = t1219 * ((t1243 + t1286) * t1296 - 0.3e1 / 0.2e1 * t1463 + t1444 + t1455) + t1297;
	t1084 = t1350 * t1522 + (t1267 + (-0.9e1 * t1289 + 0.18e2 * t1302) * t1296 + t1333) * t1300 + (t1268 + t1452) * t1308 + t1458;
	t1498 = t1084 * pkin(3);
	t1271 = -0.5e1 * t1289;
	t1272 = 7 * t1294;
	t1094 = (t1274 + t1452) * t1308 + (t1272 + (t1271 + 0.10e2 * t1302) * t1296 + t1333) * t1300 + t1458;
	t1497 = t1094 * pkin(3);
	t1496 = t1209 * pkin(3);
	t1495 = -2 * t1513;
	t1171 = t1302 + t1300 / 0.4e1 + t1296 / 0.4e1 - t1289 / 0.8e1;
	t1450 = 0.4e1 / 0.7e1 * t1302 - t1289 / 0.7e1;
	t1091 = -0.32e2 / 0.21e2 * t1171 * t1350 + t1308 / 0.7e1 + (0.16e2 / 0.21e2 * t1296 + t1450) * t1300 + t1294 / 0.7e1 + t1450 * t1296 + t1301 - 0.3e1 / 0.7e1 * t1463 + t1288 / 0.42e2 - t1291 / 0.42e2;
	t1239 = -t1289 / 0.4e1;
	t1172 = t1239 + t1254 + t1451;
	t1253 = 0.4e1 / 0.3e1 * t1296;
	t1095 = -0.8e1 / 0.3e1 * t1172 * t1350 + t1308 / 0.3e1 + (t1253 + t1240) * t1300 + t1301 - t1294 / 0.3e1 + (t1508 + 0.2e1 / 0.3e1 * t1296 + t1242) * t1302 + t1288 / 0.18e2 - t1291 / 0.18e2;
	t1255 = t1296 / 0.2e1;
	t1142 = -0.2e1 / 0.3e1 * t1350 + t1302 + t1255 + t1239;
	t1284 = 0.4e1 * t1302;
	t1205 = (t1284 + t1289) * t1296;
	t1211 = t1302 - 0.2e1 / 0.3e1 * t1300;
	t1241 = -t1289 / 0.2e1;
	t1184 = t1241 + t1361;
	t1326 = t1184 * t1338;
	t1391 = 0.16e2 * t1481;
	t1224 = t1226 ^ 2;
	t1466 = t1294 * t1224;
	t1416 = 0.8e1 * t1466;
	t1074 = t1211 * t1416 + 0.14e2 * t1091 * t1465 + t1209 * t1326 - t1443 * t1308 + (t1205 - 0.10e2 / 0.3e1 * t1294 + 0.2e1 * t1301 - t1463) * t1300 + t1140 * t1469 + (0.6e1 * t1095 * t1215 + t1142 * t1391) * pkin(1);
	t1100 = t1326 + (t1275 + t1448) * t1300 + t1322;
	t1117 = t1184 * t1468 + t1149;
	t1076 = t1100 * t1385 + t1117 * t1389 + t1084 + t1167;
	t1270 = -0.2e1 * t1289;
	t1283 = 0.8e1 * t1302;
	t1354 = t1474 * t1523;
	t1135 = t1260 * t1354 + t1515 + ((4 * t1300) + t1270 + t1283) * t1296;
	t1141 = t1239 - t1296 + t1325;
	t1421 = 0.8e1 * t1481;
	t1429 = 0.4e1 * t1215;
	t1087 = t1182 * t1469 + t1135 * t1226 + t1184 * t1221 + (t1141 * t1429 + t1421) * pkin(1);
	t1440 = t1301 - t1294;
	t1089 = t1210 * t1326 - t1297 + (-t1213 + t1445) * t1308 + (t1205 + t1440 + t1524) * t1300 + t1302 * t1140;
	t1304 = pkin(1) * t1302;
	t1191 = -(12 * pkin(1) * t1311) + t1304 * t1531;
	t1208 = -(8 * t1294) + 0.12e2 * t1462;
	t1347 = pkin(1) * t1391;
	t1105 = t1191 * t1261 + t1208 * t1226 + t1347 + t1416 + t1444 - 0.6e1 * t1462;
	t1120 = t1182 * t1468 + t1184 * t1222;
	t1177 = (-0.6e1 * t1300 * t1302 + t1439) * t1294;
	t1214 = -0.30e2 * t1289 + 0.60e2 * t1302;
	t1447 = t1288 - t1291;
	t1332 = -0.6e1 * t1463 + 0.6e1 * t1301 + t1447;
	t1062 = -0.32e2 * t1087 * t1381 + 0.16e2 * t1177 * t1224 + 0.24e2 * t1089 * t1465 + (t1270 + t1284 + (28 * t1296)) * t1297 + t1185 * t1484 + (0.24e2 * t1074 * t1229 + t1214 * t1294 + t1332 * t1275 + t1447 * t1286 - 0.6e1 * t1301 * t1289 + 0.4e1 * t1304 ^ 2 + (28 * t1311 ^ 2)) * t1300 + 0.8e1 * (-t1076 * t1487 - t1094 * t1395) * pkin(2) + (0.32e2 * t1120 * t1481 + t1498 * t1516) * pkin(1) + (0.16e2 * t1105 * t1227 + t1214 * t1296 + (70 * t1294) + t1308 + t1332) * t1308;
	t1454 = t1238 - t1292 / 0.6e1;
	t1367 = t1302 + t1454;
	t1156 = t1508 + t1255 + t1367;
	t1331 = t1257 + t1367;
	t1157 = t1253 + t1331;
	t1106 = -t1156 * t1499 + t1157 * t1216;
	t1369 = 0.2e1 / 0.3e1 * t1289 + t1247 + t1284;
	t1370 = 0.4e1 / 0.3e1 * t1289 + 0.4e1 / 0.3e1 * t1292 - 0.2e1 * t1302;
	t1383 = 0.20e2 / 0.3e1 * t1296;
	t1115 = -t1308 + (t1369 - t1383) * t1300 - (3 * t1294) + t1370 * t1296 + t1301;
	t1159 = t1296 + t1331;
	t1203 = t1514 - t1443;
	t1393 = -t1499 / 0.2e1;
	t1118 = t1159 * t1216 + t1203 * t1393;
	t1366 = t1302 + t1525;
	t1368 = t1289 / 0.3e1 + t1248 + t1286;
	t1334 = -0.8e1 / 0.3e1 * t1466 - t1198 - 0.5e1 / 0.3e1 * t1294 + t1368 * t1296 + t1302 * t1366;
	t1482 = t1225 * t1294;
	t1077 = t1106 * t1412 + t1115 * t1393 + t1334 * t1216 + (t1118 * t1215 - t1259 * t1482) * t1494;
	t1252 = -0.2e1 / 0.3e1 * t1292;
	t1453 = t1242 + t1252;
	t1107 = t1308 + (t1449 + t1453) * t1300 + t1273 + 0.2e1 * t1486 + t1302 * (t1302 + t1453);
	t1119 = t1278 + ((5 * t1296) + t1201) * t1514 + t1219 * (t1252 + t1363);
	t1088 = t1107 * t1216 - t1119 * t1499;
	t1330 = t1296 + t1366;
	t1161 = t1279 + t1330;
	t1162 = t1204 + t1525;
	t1108 = -t1161 * t1499 + t1162 * t1216;
	t1170 = t1255 + t1300 + t1454;
	t1437 = 0.2e1 * t1216;
	t1123 = t1170 * t1437 + t1468 * t1499;
	t1400 = t1311 * t1216;
	t1414 = -0.4e1 * t1465;
	t1078 = t1123 * t1414 + (t1108 * t1429 + t1400 * t1519) * pkin(1) + t1088;
	t1116 = -(3 * t1308) + (t1370 - t1383) * t1300 + t1369 * t1296 + t1440;
	t1124 = -0.5e1 / 0.3e1 * t1308 + (-t1296 + t1368) * t1300 + t1302 * t1330;
	t1432 = -0.2e1 * t1499;
	t1090 = t1116 * t1216 + t1124 * t1432;
	t1365 = t1241 - t1292 / 0.2e1 + t1302;
	t1158 = 0.3e1 / 0.2e1 * t1300 + t1276 + t1365;
	t1174 = t1216 + 0.2e1 * t1499;
	t1093 = t1221 * t1216 + 0.4e1 * t1488 + (t1158 * t1259 + t1174 * t1504) * t1532;
	t1160 = t1279 + 0.3e1 / 0.2e1 * t1296 + t1365;
	t1122 = t1160 * t1216 + t1222 * t1499 / 0.2e1;
	t1327 = 0.24e2 * t1210 * t1466 - t1297 - ((21 * t1296) + t1201) * t1308 - (t1175 + t1282 + (35 * t1294) + 0.2e1 * t1483) * t1300 - (t1272 + (t1271 + t1283 - 0.5e1 * t1292) * t1296 + t1302 * (-t1292 - t1445)) * t1219;
	t1388 = -0.12e2 * t1464;
	t1415 = -0.6e1 * t1465;
	t1066 = t1122 * t1347 + t1093 * t1343 + t1077 * t1388 + t1090 * t1415 + (-0.6e1 * t1088 * t1504 + t1527 * t1259) * pkin(3) + (0.6e1 * t1078 * t1487 + t1327 * t1260) * pkin(2);
	t1059 = t1062 * t1164 + t1066 * t1303;
	t1144 = t1197 + t1292 + t1329;
	t1398 = t1148 * t1215;
	t1101 = t1144 * t1192 + 0.2e1 * t1398;
	t1104 = t1144 * t1261 + (t1518 - 0.2e1) * t1501;
	t1456 = -t1187 + t1176;
	t1137 = pkin(3) + t1456;
	t1073 = t1101 * t1259 + t1104 * t1216 + t1137 * t1303;
	t1178 = t1276 + t1292 + t1359;
	t1127 = t1178 + t1197 + t1338;
	t1473 = t1261 * t1260;
	t1139 = pkin(2) * t1473 + t1192 * t1259;
	t1490 = t1139 * t1303;
	t1072 = -t1127 * t1176 + t1490 + (t1178 * t1476 - 0.2e1 * t1262 * t1505) * pkin(2) + (-t1279 - t1292 - t1296 + t1328 + t1410) * pkin(3);
	t1375 = t1072 * t1510;
	t1323 = t1059 * t1375 + t1073 * t1507;
	t1109 = 0.1e1 / t1112;
	t1467 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
	t1382 = t1109 * t1467;
	t1057 = t1323 * t1382;
	t1472 = t1261 * t1262;
	t1152 = -t1472 - t1476;
	t1056 = t1057 * t1152;
	t1374 = t1073 * t1510;
	t1324 = t1059 * t1374 + t1072 * t1506;
	t1058 = t1324 * t1382;
	t1475 = t1259 * t1262;
	t1153 = -t1473 + t1475;
	t1052 = -t1058 * t1153 - t1056;
	t1050 = 0.1e1 / t1052 ^ 2;
	t1055 = t1057 * t1153;
	t1051 = -t1058 * t1152 + t1055;
	t1493 = t1050 * t1051;
	t1081 = 0.1e1 / t1303;
	t1207 = pkin(1) * t1216;
	t1196 = 0.2e1 * t1207;
	t1470 = t1262 * t1300;
	t1417 = -0.4e1 * t1470;
	t1342 = t1260 * t1417;
	t1151 = t1196 + t1342;
	t1402 = pkin(2) * t1475;
	t1349 = pkin(3) * t1402;
	t1352 = -0.8e1 * t1397;
	t1431 = -0.4e1 * t1215;
	t1492 = (t1151 * t1414 + (t1207 - t1349) * t1352 + 0.4e1 * t1145 * t1349 + (pkin(2) * t1126 * t1431 - 0.8e1 * t1198 * t1262 + (t1179 * t1523 + t1395 * t1528) * pkin(1)) * t1260) * t1081;
	t1413 = 0.2e1 * t1465;
	t1471 = t1261 * t1296;
	t1491 = t1081 * ((t1126 * t1501 + 0.2e1 * t1148 * t1471) * t1259 + (t1145 * t1215 + t1192 * t1413) * t1216);
	t1461 = t1311 * t1226;
	t1406 = pkin(1) * t1461;
	t1345 = -0.12e2 * t1406;
	t1489 = t1164 * t1259 * t1345;
	t1478 = t1228 * t1308;
	t1477 = t1229 * t1307;
	t1194 = pkin(1) * t1432;
	t1223 = t1259 ^ 2;
	t1335 = t1225 * t1377;
	t1384 = 0.32e2 / 0.3e1 * t1294;
	t1336 = t1225 * t1384;
	t1337 = 0.64e2 / 0.3e1 * t1171 * t1311;
	t1341 = t1468 * t1503;
	t1346 = -0.48e2 * t1406;
	t1348 = -0.16e2 * pkin(1) * t1172 * t1296;
	t1353 = -0.4e1 * t1184 * t1496;
	t1404 = pkin(1) * t1465;
	t1355 = 0.4e1 * t1404;
	t1356 = -0.2e1 * t1404;
	t1357 = -0.4e1 * t1404;
	t1376 = t1259 * t1471;
	t1378 = t1210 * t1481;
	t1390 = -0.64e2 * t1479;
	t1392 = -0.32e2 * t1482;
	t1394 = pkin(3) * t1469;
	t1396 = pkin(3) * t1479;
	t1403 = t1164 * t1216;
	t1408 = pkin(1) * t1466;
	t1409 = 0.3e1 * t1464;
	t1419 = 0.8e1 * t1479;
	t1420 = -0.8e1 * t1479;
	t1422 = -0.2e1 * t1481;
	t1423 = -0.4e1 * t1481;
	t1424 = 0.2e1 * t1491;
	t1430 = 0.2e1 * t1215;
	t1434 = 0.6e1 * t1502;
	t1436 = -0.6e1 * t1502;
	t1053 = ((0.12e2 * t1223 * t1226 * t1503 + t1156 * t1423 + t1203 * t1356 - 0.4e1 * t1408) * t1388 + 0.8e1 * t1222 * t1408 + 0.12e2 * t1124 * t1481 + 0.6e1 * t1119 * t1404 + (-t1115 * t1388 / 0.2e1 + t1527) * t1215) * t1303 + t1066 * t1424 + (((t1158 * t1430 + t1355 + t1423) * t1419 + (-t1119 * t1215 + t1161 * t1357 - 0.4e1 * t1335) * t1434) * t1303 + t1390 * t1489) * t1193 + (0.24e2 * (-pkin(1) * t1224 * t1384 - t1225 * t1337 + t1226 * t1348 + t1261 * t1353) * t1464 - 0.64e2 * t1224 * t1341 - 0.96e2 * t1184 * t1378 - 0.48e2 * t1114 * t1404 + t1497 * t1517 + ((-t1261 * t1394 + t1356 + t1422) * t1390 - 0.48e2 * (-t1114 * t1215 + t1184 * t1357 - 0.4e1 * t1378) * t1502) * t1193) * t1403 + (-pkin(3) * t1062 + (0.2e1 * (-0.2e1 * t1208 * t1261 - t1191 + t1346 + t1392) * t1480 + 0.4e1 * t1087 * t1396 + (-t1135 * t1261 + t1141 * t1495) * t1344 + (-0.28e2 * t1091 * t1471 - 0.6e1 * t1095 * t1513 + t1142 * t1346 + t1211 * t1392) * t1409 + pkin(3) * t1076 * t1502 + t1193 * (-t1100 * t1513 - 0.4e1 * t1117 * t1471 - 0.4e1 * t1222 * t1406) * t1436 + t1177 * t1519 + t1120 * t1345 - 0.6e1 * t1089 * t1471 - pkin(1) * t1498) * t1520 + ((-0.8e1 * t1106 * t1471 + t1336 * t1216) * t1388 - 0.96e2 * t1210 * t1482 * t1216 + t1122 * t1346 + 0.12e2 * t1090 * t1471 + (t1078 * t1436 + t1093 * t1420 + (-0.24e2 * t1194 + 0.64e2 * t1376) * t1480 + (0.48e2 * t1118 * t1464 + 0.6e1 * t1088) * pkin(1)) * pkin(3) + ((t1165 * t1517 + t1174 * t1495) * t1419 + (0.24e2 * pkin(1) * t1226 * t1400 - 0.4e1 * t1108 * t1513 + 0.8e1 * t1123 * t1471) * t1434) * t1193) * t1303) * t1259;
	t1387 = 0.12e2 * t1464;
	t1418 = t1260 * t1528;
	t1427 = -0.2e1 * t1496;
	t1428 = t1098 * t1523;
	t1433 = -0.4e1 * pkin(3) * t1163;
	t1435 = 0.4e1 * t1502;
	t1060 = (t1296 * t1223 * t1420 + (t1218 * t1215 + t1422) * t1411 + 0.4e1 * t1335 + t1202 * t1355 + t1154 * t1215) * t1303 + t1075 * t1424 + t1387 * t1489 + ((-0.8e1 / 0.3e1 * t1481 + t1357 + t1261 * t1427) * t1387 - 0.24e2 * t1378 - 0.24e2 * t1163 * t1404 - 0.6e1 * t1113 * t1215) * t1403 + (0.24e2 * (-t1164 * t1222 - t1303 * t1216) * t1406 + ((0.16e2 * t1165 * t1464 - 0.2e1 * t1131) * t1303 + t1164 * (-0.144e3 * t1130 * t1464 - 0.24e2 * t1111) * t1296) * t1261 + (t1262 * t1303 * t1428 - t1071 + t1164 * (t1086 * t1434 + t1147 * t1419) + ((pkin(2) * t1229 * t1418 + 0.4e1 * t1129) * t1303 + t1164 * (-0.48e2 * t1138 * t1464 - 0.6e1 * t1097)) * pkin(1)) * pkin(3)) * t1259 + ((t1396 * t1516 + ((-pkin(3) * t1202 + 0.4e1 * t1296 * t1187) * t1261 + (-t1173 * t1499 - t1465) * t1533) * t1435) * t1303 + t1164 * ((t1194 - 0.8e1 * t1376) * t1420 + ((-0.8e1 * t1199 * t1259 + t1433 * t1216) * t1261 + (-0.4e1 * t1136 * t1499 - 0.8e1 * t1340) * pkin(1)) * t1436)) * t1193;
	t1069 = -t1490 + t1137 * t1424 + pkin(3) * t1223 * t1521 + t1101 * t1261 + (-t1144 + t1352) * t1187;
	t1438 = -2 * pkin(1) * t1300;
	t1070 = t1456 * t1303 + t1139 * t1424 + (t1192 * t1127 + 0.4e1 * t1398) * t1259 + (t1438 * t1472 + (t1178 * t1261 + 0.4e1 * t1192 * t1500) * pkin(2)) * t1260;
	t1132 = t1139 * pkin(3);
	t1320 = t1323 * t1530;
	t1373 = t1072 * t1509;
	t1460 = (-t1132 * t1320 + (t1053 * t1375 + t1069 * t1507 - t1073 * t1491 / 0.2e1 + (t1060 * t1373 + t1070 * t1510) * t1059) * t1109) * t1467 - t1058;
	t1399 = pkin(2) * t1472;
	t1351 = pkin(1) * t1399;
	t1155 = -0.4e1 * t1296 * t1259 * t1351;
	t1169 = t1351 * t1532;
	t1371 = t1492 / 0.2e1;
	t1379 = t1210 * t1461;
	t1380 = t1193 * t1477;
	t1405 = pkin(1) * t1471;
	t1054 = t1066 * t1371 + (-0.32e2 * t1155 * t1164 + 0.8e1 * t1169 * t1303) * t1381 + (0.24e2 * (0.4e1 * t1143 * t1478 * t1499 + t1077 * t1470 - t1093 * t1380) * t1303 + (-0.6e1 * t1074 * t1470 + 0.12e2 * t1087 * t1380 - 0.8e1 * t1105 * t1478) * t1520) * t1260 + ((t1062 + (t1076 * t1520 - 0.6e1 * t1078 * t1303) * t1193) * t1260 + (((t1157 * t1412 + t1159 * t1386 + t1334) * t1388 + t1160 * t1347 + t1116 * t1415 - 0.6e1 * t1107 * t1425 + (t1526 * t1419 + (t1162 * t1386 - 0.8e1 * t1170 * t1465 + t1107 - 0.8e1 * t1407) * t1434) * t1193 + t1327) * t1303 + ((-pkin(1) * t1336 - t1226 * t1337 + t1261 * t1348 + t1353) * t1409 + t1341 * t1519 - 0.12e2 * t1184 * t1379 + t1405 * t1522 - t1497 + (-0.4e1 * (-(2 * t1394) - 0.4e1 * t1461) * t1479 - (pkin(3) * t1522 - 0.24e2 * t1184 * t1405 - 0.24e2 * t1379) * t1502) * t1193) * t1259 * t1520) * t1262) * pkin(2);
	t1061 = (t1169 * t1411 + (t1099 * t1418 + t1200 * t1226 + ((-t1443 + t1413) * t1411 - t1146 + (t1204 * t1431 + t1421) * pkin(1)) * pkin(2)) * t1262 + ((t1169 + (t1204 - 0.2e1 * t1465) * t1502) * t1435 + (-0.24e2 * t1477 * t1499 + t1428) * t1260) * t1193) * t1303 + t1075 * t1371 + t1071 * t1216 + 0.6e1 * t1164 * (0.4e1 * t1147 * t1260 * t1380 + (t1155 + (-0.8e1 / 0.3e1 * t1461 + t1427) * t1402) * t1410 + t1085 * t1342 + t1210 * t1262 * t1226 * t1354 + t1163 * t1155 - t1113 * t1349 + (-(-0.8e1 * t1405 + t1433) * t1259 * t1464 + t1086 * t1216) * t1193);
	t1067 = t1137 * t1371 + t1151 * t1259 * t1430 + ((-t1259 * t1303 + t1104) * t1262 + (t1261 * t1303 + (t1192 * t1533 + t1144) * t1259 + (-pkin(3) + 0.2e1 * t1500 + t1504) * t1437) * t1260) * pkin(2);
	t1068 = (t1187 + t1399) * t1303 + t1139 * t1371 - 0.2e1 * t1151 * t1500 - (t1196 - 0.4e1 * t1349) * t1176 + t1229 * t1259 * t1438 + t1178 * t1402 + (pkin(3) * t1417 + (-t1127 * t1261 + 0.2e1 * t1505) * pkin(2)) * t1260;
	t1128 = -t1153 * pkin(3) * pkin(2) + t1207;
	t1459 = (t1128 * t1320 + (t1054 * t1375 + t1067 * t1507 - t1073 * t1492 / 0.8e1 + (t1061 * t1373 + t1068 * t1510) * t1059) * t1109) * t1467 + t1058;
	t1372 = t1073 * t1509;
	t1321 = t1324 * t1530;
	t1049 = 0.1e1 / t1052;
	t1048 = 0.1e1 / (t1050 * t1051 ^ 2 + 0.1e1);
	t1047 = (t1128 * t1321 + (t1068 * t1506 + t1072 * t1492 / 0.8e1 + t1054 * t1374 + (t1061 * t1372 + t1067 * t1510) * t1059) * t1109) * t1467;
	t1045 = (-t1132 * t1321 + (t1070 * t1506 + t1072 * t1491 / 0.2e1 + t1053 * t1374 + (t1060 * t1372 + t1069 * t1510) * t1059) * t1109) * t1467;
	t1 = [0, 0; 0, 0; 0.1e1 + ((t1459 * t1153 + (-t1047 + t1057) * t1152) * t1049 - (-t1047 * t1153 - t1459 * t1152 + t1055) * t1493) * t1048, ((-t1045 * t1152 + t1153 * t1460 - t1056) * t1049 - ((-t1045 - t1057) * t1153 - t1460 * t1152) * t1493) * t1048;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 1;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:52:02
	% EndTime: 2020-04-27 03:52:03
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (1788->105), mult. (3903->181), div. (28->6), fcn. (958->6), ass. (0->81)
	t314 = -4 * pkin(2);
	t269 = pkin(2) ^ 2;
	t263 = cos(qJ(1));
	t305 = pkin(2) * t263;
	t288 = pkin(1) * t305;
	t255 = -0.2e1 * t288;
	t270 = pkin(1) ^ 2;
	t291 = t255 + t270;
	t259 = t263 ^ 2;
	t311 = 0.2e1 * t259;
	t245 = t269 * t311 - t269 + t291;
	t313 = -0.2e1 * t245;
	t253 = pkin(1) - t305;
	t312 = 0.2e1 * t253;
	t310 = (-pkin(4) - pkin(5));
	t309 = (-pkin(4) + pkin(5));
	t260 = sin(qJ(2));
	t261 = sin(qJ(1));
	t296 = t260 * t261;
	t252 = pkin(2) * t296;
	t279 = pkin(3) * t252;
	t251 = -0.2e1 * t279;
	t268 = pkin(3) ^ 2;
	t262 = cos(qJ(2));
	t247 = t253 * t262;
	t283 = pkin(3) * t247;
	t236 = t251 + t268 + t269 + 0.2e1 * t283 + t291;
	t234 = 0.1e1 / t236;
	t308 = t234 / 0.2e1;
	t307 = pkin(1) * (t252 - pkin(3));
	t306 = pkin(2) * t261;
	t304 = pkin(3) * t253;
	t258 = t262 ^ 2;
	t303 = pkin(3) * t258;
	t302 = pkin(3) * t262;
	t235 = 0.1e1 / t236 ^ 2;
	t257 = pkin(1) * t306;
	t294 = t261 * t262;
	t295 = t260 * t263;
	t301 = t235 * (t257 + (t294 - t295) * pkin(3) * pkin(2));
	t242 = pkin(2) * t294 + t260 * t253;
	t300 = t235 * t242 * pkin(3);
	t266 = pkin(4) ^ 2;
	t290 = -pkin(5) ^ 2 + t270;
	t282 = t269 + t290;
	t276 = t268 + t282;
	t249 = -t266 + t276;
	t244 = t255 + t249;
	t237 = t251 + t244;
	t256 = (t268 - t270) * t269;
	t275 = t258 * t313 - t290;
	t271 = sqrt(0.4e1 * t256 * t259 + 0.4e1 * t249 * t288 - (t270 + ((pkin(2) - t309) * (pkin(2) + t309))) * (t270 + ((pkin(2) - t310) * (pkin(2) + t310))) + 0.4e1 * (-t237 * t247 + t244 * t252) * pkin(3) + (0.2e1 * t266 - (6 * t269) + 0.2e1 * t275 - t268) * t268);
	t299 = t242 * t271;
	t298 = t245 * t262;
	t297 = t258 * t268;
	t293 = t262 * t263;
	t292 = -t252 + t247;
	t289 = -0.2e1 * pkin(1) * t269;
	t229 = 0.1e1 / t271;
	t287 = 0.2e1 * t229 * ((t237 * t304 + 0.2e1 * t268 * t298) * t260 + (t244 * t302 + t297 * t312) * t306);
	t286 = -0.4e1 * t263 * t269;
	t285 = pkin(2) * t295;
	t284 = pkin(3) * t298;
	t254 = 0.2e1 * t257;
	t246 = t261 * t286 + t254;
	t278 = pkin(3) * t285;
	t280 = -0.8e1 * t283;
	t281 = (-0.4e1 * t246 * t297 + (t257 - t278) * t280 + 0.4e1 * t244 * t278 + (t237 * t302 * t314 - 0.8e1 * t256 * t263 + (0.8e1 * pkin(3) * t269 * t296 + t249 * t314) * pkin(1)) * t261) * t229 / 0.2e1;
	t248 = 0.3e1 * t268 + t266 + t282;
	t238 = t248 + t255 - 0.4e1 * t279;
	t225 = -t238 * t247 + t299 + (t248 * t296 - 0.2e1 * t263 * t307) * pkin(2) + (-t266 - t268 + (t311 - 0.3e1) * t269 + t275) * pkin(3);
	t224 = 0.1e1 / t225 ^ 2;
	t243 = t255 + t266 + t276;
	t232 = t243 * t253 + 0.2e1 * t284;
	t233 = t243 * t262 + (0.4e1 * t258 - 0.2e1) * t304;
	t241 = pkin(3) + t292;
	t226 = t232 * t260 + t233 * t306 + t241 * t271;
	t277 = 0.1e1 / (t224 * t226 ^ 2 + 0.1e1) * t236;
	t274 = 0.1e1 / t225 * t277;
	t273 = t224 * t226 * t277;
	t1 = [0, 0; 0, 0; 0.2e1 * ((0.2e1 * t246 * t260 * t302 + t241 * t281) * t308 - t226 * t301 + ((-t260 * t271 + t233) * t263 * t308 + (t262 * t271 / 0.2e1 + (pkin(1) * t312 + t243) * t260 / 0.2e1 + (pkin(1) * t262 - pkin(3) + 0.2e1 * t303) * t306) * t234 * t261) * pkin(2)) * t274 - 0.2e1 * (((pkin(2) * t293 + t252) * t271 + t242 * t281 - 0.2e1 * t246 * t303 - (t254 - 0.4e1 * t278) * t247 + t259 * t260 * t289 + t248 * t285 + (pkin(3) * t286 + (-t238 * t262 + 0.2e1 * t307) * pkin(2)) * t261) * t308 - t225 * t301) * t273, 0.1e1 + 0.2e1 * ((-t299 + t241 * t287 + t232 * t262 + ((-t243 + t280) * t306 + pkin(3) * t260 * t313) * t260) * t308 + t226 * t300) * t274 - 0.2e1 * ((t292 * t271 + t242 * t287 + (t253 * t238 + 0.4e1 * t284) * t260 + (t289 * t293 + (t248 * t262 + 0.4e1 * t253 * t303) * pkin(2)) * t261) * t308 + t225 * t300) * t273;];
	Jg_rot = t1;
end