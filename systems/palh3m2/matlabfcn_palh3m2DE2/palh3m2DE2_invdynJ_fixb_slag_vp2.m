% Calculate vector of inverse dynamics joint torques for
% palh3m2DE2
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
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2DE2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:24
% EndTime: 2020-05-07 04:22:39
% DurationCPUTime: 14.71s
% Computational Cost: add. (11464->575), mult. (20187->801), div. (0->0), fcn. (24274->79), ass. (0->302)
t1262 = pkin(17) + pkin(18);
t1244 = sin(t1262);
t1245 = cos(t1262);
t1273 = sin(qJ(2));
t1276 = cos(qJ(2));
t1268 = cos(pkin(16));
t1277 = cos(pkin(15));
t1398 = sin(pkin(16));
t1426 = sin(pkin(15));
t1182 = t1268 * t1426 + t1277 * t1398;
t1183 = t1268 * t1277 - t1398 * t1426;
t1272 = sin(qJ(3));
t1275 = cos(qJ(3));
t1305 = t1182 * t1272 - t1183 * t1275;
t1320 = -t1182 * t1275 - t1272 * t1183;
t1322 = t1273 * t1305 + t1320 * t1276;
t1476 = t1305 * t1276;
t1457 = -t1273 * t1320 + t1476;
t1029 = t1322 * t1244 - t1457 * t1245;
t1261 = qJD(2) + qJD(3);
t1371 = qJD(2) * t1275;
t1333 = pkin(1) * t1371;
t1302 = -pkin(4) * t1261 + t1333;
t1492 = t1029 * t1302;
t1406 = pkin(1) * t1272;
t1331 = t1029 * t1406;
t1025 = qJD(2) * t1331;
t1485 = t1244 * t1457 + t1245 * t1322;
t1486 = t1485 * t1302;
t1009 = -t1025 - t1486;
t1102 = t1305 * qJD(3);
t1103 = t1320 * qJD(3);
t1489 = t1457 * qJD(2) - t1273 * t1103;
t1030 = -t1102 * t1276 - t1489;
t1031 = t1322 * qJD(2) + t1102 * t1273 + t1103 * t1276;
t1363 = t1031 * t1244;
t1015 = t1030 * t1245 + t1363;
t1017 = -t1030 * t1244 + t1031 * t1245;
t1236 = -pkin(1) * t1275 + pkin(4);
t1366 = t1485 * t1275;
t1374 = -t1017 * t1236 - qJD(3) * t1331 + (-qJD(3) * t1366 + t1015 * t1272) * pkin(1) + t1009;
t1270 = cos(pkin(18));
t1408 = sin(pkin(18));
t1187 = t1270 * t1277 - t1408 * t1426;
t1228 = -t1276 * pkin(1) - pkin(12);
t1294 = -t1270 * t1426 - t1277 * t1408;
t1074 = (cos(pkin(17)) * t1187 + t1294 * sin(pkin(17))) * pkin(3) + t1228;
t1446 = m(4) + m(8);
t1330 = t1446 * t1228;
t1491 = m(9) * t1074 + t1330;
t1376 = pkin(15) + qJ(2);
t1255 = pkin(18) + t1376;
t1206 = pkin(17) + qJ(3) + t1255;
t1375 = qJ(2) + atan2(sin(t1255), -cos(t1255));
t1328 = pkin(17) - t1375;
t1093 = -atan2(-sin(t1206), -cos(t1206)) + t1328;
t1089 = qJ(1) + t1093;
t1263 = qJ(2) + qJ(3);
t1258 = qJ(1) - t1263;
t1279 = mrSges(9,2) * g(2);
t1428 = mrSges(9,2) * g(1);
t1436 = mrSges(4,2) * g(1);
t1438 = mrSges(9,1) * g(2);
t1439 = mrSges(9,1) * g(1);
t1282 = m(5) + m(6);
t1453 = pkin(4) * t1282 + mrSges(4,1);
t1479 = g(2) * t1453;
t1488 = -(t1479 + t1436) * cos(t1258) / 0.2e1 - (t1428 + t1438) * cos(t1089) / 0.2e1 + (-t1279 + t1439) * sin(t1089) / 0.2e1;
t1090 = -qJ(1) + t1093;
t1257 = qJ(1) + t1263;
t1482 = (t1279 + t1439) * sin(t1090) / 0.2e1 + (-t1428 + t1438) * cos(t1090) / 0.2e1 + (t1479 - t1436) * cos(t1257) / 0.2e1;
t1481 = mrSges(4,1) + mrSges(9,1);
t1410 = mrSges(4,2) + mrSges(9,2);
t1480 = g(1) * t1453;
t1336 = qJDD(2) * t1275;
t1372 = qJD(2) * t1272;
t1298 = qJD(3) * t1372 - t1336;
t1180 = t1298 * pkin(1);
t1278 = cos(pkin(14));
t1427 = sin(pkin(14));
t1188 = -t1277 * t1427 + t1278 * t1426;
t1191 = t1277 * t1278 + t1426 * t1427;
t1117 = t1188 * t1276 + t1191 * t1273;
t1478 = Ifges(7,4) * t1117;
t1334 = pkin(1) * t1372;
t1007 = -t1485 * t1334 + t1492;
t1016 = (qJD(3) * t1476 + t1489) * t1245 - t1363;
t998 = t1016 * t1236 + (-t1017 * t1272 + (-t1029 * t1275 + t1272 * t1485) * qJD(3)) * pkin(1);
t1477 = -t1007 + t998;
t1332 = t1485 * t1406;
t1008 = -qJD(2) * t1332 + t1492;
t1012 = -t1029 * t1236 - t1332;
t1181 = (-qJD(3) * t1371 - qJDD(2) * t1272) * pkin(1);
t1260 = qJDD(2) + qJDD(3);
t995 = t1015 * t1334 + t1017 * t1302 - t1029 * (pkin(4) * t1260 + t1180) + t1485 * t1181;
t1475 = t1374 * t1008 + t1012 * t995;
t1246 = sin(t1263);
t1251 = cos(t1263);
t1474 = t1246 * t1272 + t1251 * t1275;
t1384 = Ifges(9,2) * t1251;
t1387 = Ifges(9,4) * t1251;
t1388 = Ifges(9,4) * t1246;
t1392 = Ifges(9,1) * t1246;
t1472 = -t1246 * (Ifges(9,6) * t1261 + (-t1384 - t1388) * qJD(1)) / 0.2e1 + t1251 * (Ifges(9,5) * t1261 + (-t1387 - t1392) * qJD(1)) / 0.2e1;
t1471 = -mrSges(4,2) * t1246 + t1453 * t1251;
t1199 = pkin(16) + t1206;
t1329 = atan2(-sin(t1199), cos(t1199)) + t1263;
t1133 = -qJ(4) + t1329;
t1130 = qJ(1) + t1133;
t1132 = qJ(4) + t1329;
t1131 = qJ(1) - t1132;
t1265 = qJ(1) - qJ(4);
t1433 = mrSges(6,2) * g(2);
t1434 = mrSges(6,2) * g(1);
t1443 = mrSges(6,1) * g(2);
t1444 = mrSges(6,1) * g(1);
t1470 = -(sin(t1265) / 0.2e1 - sin(t1130) / 0.4e1 - sin(t1131) / 0.4e1) * (-t1433 + t1444) + (cos(t1265) / 0.2e1 - cos(t1130) / 0.4e1 - cos(t1131) / 0.4e1) * (t1434 + t1443);
t1469 = -m(9) / 0.2e1;
t1319 = -t1188 * t1273 + t1191 * t1276;
t1468 = -t1319 / 0.2e1;
t1467 = Ifges(7,4) * t1319;
t1091 = t1182 * t1245 + t1183 * t1244;
t1080 = t1091 * qJD(1);
t1092 = t1182 * t1244 - t1183 * t1245;
t1081 = t1092 * qJD(1);
t1190 = t1272 * t1273 - t1275 * t1276;
t1172 = t1190 * qJD(1);
t1136 = -pkin(4) * t1172 + t1228 * qJD(1);
t1037 = -pkin(8) * t1081 - pkin(10) * t1080 + t1136;
t1271 = sin(qJ(4));
t1274 = cos(qJ(4));
t1003 = -t1009 * t1271 + t1037 * t1274;
t1466 = qJD(4) * t1003;
t1004 = t1009 * t1274 + t1037 * t1271;
t1465 = qJD(4) * t1004;
t1464 = t1074 * (-mrSges(9,1) * t1246 - mrSges(9,2) * t1251);
t1046 = -mrSges(5,1) * t1081 + mrSges(5,2) * t1080;
t1459 = m(5) * t1136 + t1046;
t1456 = -t1003 * t1271 + t1004 * t1274;
t1079 = t1092 * qJDD(1);
t1078 = t1091 * qJDD(1);
t1359 = t1080 * t1271;
t1042 = -qJD(4) * t1359 + t1078 * t1274;
t1358 = t1080 * t1274;
t1043 = -qJD(4) * t1358 - t1078 * t1271;
t1073 = qJDD(4) - t1079;
t1345 = t1273 * qJD(1);
t1327 = qJD(2) * t1345;
t1196 = qJDD(1) * t1276 - t1327;
t1341 = t1276 * qJD(1);
t1197 = qJD(2) * t1341 + qJDD(1) * t1273;
t1304 = t1272 * t1276 + t1273 * t1275;
t1296 = t1304 * qJD(3);
t1071 = qJD(1) * t1296 - t1196 * t1275 + t1197 * t1272;
t1211 = pkin(1) * t1327;
t1337 = t1228 * qJDD(1) + t1211;
t1058 = -pkin(4) * t1071 + t1337;
t1032 = -pkin(8) * t1079 - pkin(10) * t1078 + t1058;
t996 = t1029 * t1181 + (t1016 * t1261 + t1260 * t1485) * pkin(4) + (-t1485 * t1336 + (-t1016 * t1275 + (qJD(3) * t1485 - t1017) * t1272) * qJD(2)) * pkin(1);
t993 = t1032 * t1271 + t1274 * t996 + t1466;
t994 = t1032 * t1274 - t1271 * t996 - t1465;
t1451 = t994 * mrSges(6,1) - t993 * mrSges(6,2) + Ifges(6,5) * t1042 + Ifges(6,6) * t1043 + Ifges(6,3) * t1073;
t1075 = qJD(4) - t1081;
t1450 = -mrSges(6,3) * t1456 + t1008 * (mrSges(6,1) * t1274 - mrSges(6,2) * t1271) + t1075 * (-Ifges(6,5) * t1271 - Ifges(6,6) * t1274) + (Ifges(6,4) * t1271 ^ 2 - (Ifges(6,4) * t1274 + (Ifges(6,1) - Ifges(6,2)) * t1271) * t1274) * t1080;
t1449 = 0.2e1 * g(1);
t1284 = qJD(1) ^ 2;
t1448 = 0.2e1 * pkin(8) * m(6) + 0.2e1 * mrSges(5,1);
t1447 = m(8) / 0.2e1;
t1445 = m(7) * pkin(6);
t1442 = mrSges(7,1) * g(1);
t1441 = mrSges(7,1) * g(2);
t1440 = mrSges(8,1) * g(2);
t1437 = mrSges(3,2) * g(1);
t1435 = mrSges(4,2) * g(2);
t1432 = mrSges(7,2) * g(1);
t1431 = mrSges(7,2) * g(2);
t1430 = mrSges(8,2) * g(1);
t1429 = mrSges(8,2) * g(2);
t1106 = t1117 * qJD(2);
t1424 = -t1106 / 0.2e1;
t1107 = t1319 * qJD(2);
t1423 = t1107 / 0.2e1;
t1295 = t1190 * qJD(3);
t1120 = qJD(2) * t1190 + t1295;
t1422 = t1120 / 0.2e1;
t1121 = qJD(2) * t1304 + t1296;
t1421 = t1121 / 0.2e1;
t1419 = t1172 / 0.2e1;
t1173 = t1304 * qJD(1);
t1418 = -t1173 / 0.2e1;
t1413 = t1273 / 0.2e1;
t1411 = t1276 / 0.2e1;
t1135 = qJ(1) - t1329;
t1409 = cos(t1135);
t1407 = pkin(1) * qJD(2);
t1405 = pkin(4) * t1173;
t1243 = m(9) + t1282 + t1446;
t1200 = pkin(1) * t1243 + mrSges(3,1);
t1402 = g(1) * t1200;
t1401 = g(2) * t1200;
t1242 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t1399 = g(2) * t1242;
t1397 = mrSges(3,2) * t1273;
t1395 = mrSges(4,3) * t1172;
t1394 = mrSges(4,3) * t1173;
t1393 = mrSges(9,3) * qJD(1);
t1391 = Ifges(3,4) * t1273;
t1386 = Ifges(3,5) * t1276;
t1385 = Ifges(9,5) * t1251;
t1383 = Ifges(3,6) * t1273;
t1382 = Ifges(9,6) * t1246;
t1381 = t1117 * Ifges(7,1);
t1380 = t1319 * Ifges(7,2);
t1379 = t1173 * Ifges(4,4);
t1378 = t1298 * mrSges(9,1);
t1377 = t1276 * Ifges(3,2);
t1373 = qJD(1) * t1261;
t1002 = t1486 - t1025 + (t1029 * t1272 - t1366) * t1407;
t1370 = t1002 * t1008;
t1355 = t1181 * t1272;
t1350 = t1246 * qJD(1);
t1339 = t1251 * t1393 + t1410 * t1261 - t1395;
t1338 = mrSges(9,3) * t1350 + t1481 * t1261 + t1394;
t1134 = qJ(1) + t1329;
t1335 = cos(t1134) / 0.4e1;
t1256 = pkin(14) - t1376;
t1316 = mrSges(9,1) * t1251 - mrSges(9,2) * t1246;
t1317 = mrSges(8,1) * t1187 - mrSges(8,2) * t1294;
t1324 = mrSges(4,1) * t1172 + mrSges(4,2) * t1173 + (-t1316 - t1317) * qJD(1);
t1313 = -t1008 * t1017 - t1029 * t995;
t1312 = t1003 * t1274 + t1004 * t1271;
t1044 = -mrSges(6,2) * t1075 - mrSges(6,3) * t1359;
t1045 = mrSges(6,1) * t1075 - mrSges(6,3) * t1358;
t1310 = t1044 * t1274 - t1045 * t1271;
t1309 = -t1271 * t1044 - t1274 * t1045;
t1306 = -qJD(1) * t1464 - t1228 * (-mrSges(4,1) * t1173 + mrSges(4,2) * t1172);
t1114 = t1187 * t1273 - t1276 * t1294;
t1113 = t1187 * t1276 + t1273 * t1294;
t1303 = -pkin(4) * t1190 + t1228;
t1301 = t1246 * (-Ifges(9,1) * t1251 + t1388);
t1300 = t1251 * (Ifges(9,2) * t1246 - t1387);
t1237 = Ifges(3,4) * t1341;
t1299 = -t1273 * (Ifges(3,6) * qJD(2) + (t1377 + t1391) * qJD(1)) / 0.2e1 + (Ifges(3,1) * t1345 + Ifges(3,5) * qJD(2) + t1237) * t1411;
t1128 = qJ(1) + t1132;
t1129 = qJ(1) - t1133;
t1224 = t1433 + t1444;
t1227 = -t1434 + t1443;
t1264 = qJ(1) + qJ(4);
t1293 = t1227 * cos(t1264) / 0.2e1 - sin(t1264) * t1224 / 0.2e1 - (sin(t1128) + sin(t1129)) * t1224 / 0.4e1 + (cos(t1128) + cos(t1129)) * t1227 / 0.4e1;
t1238 = qJ(1) + t1256;
t1267 = qJ(1) - qJ(2);
t1292 = (t1401 + t1437) * cos(t1267) / 0.2e1 - (-t1431 + t1442) * sin(t1238) / 0.2e1 + (t1432 + t1441) * cos(t1238) / 0.2e1;
t1239 = -qJ(1) + t1256;
t1266 = qJ(1) + qJ(2);
t1291 = (-t1401 + t1437) * cos(t1266) / 0.2e1 + (-t1431 - t1442) * sin(t1239) / 0.2e1 + (t1432 - t1441) * cos(t1239) / 0.2e1;
t1033 = mrSges(6,1) * t1073 - mrSges(6,3) * t1042;
t1034 = -mrSges(6,2) * t1073 + mrSges(6,3) * t1043;
t1290 = qJD(4) * t1309 - t1271 * t1033 + t1274 * t1034;
t1289 = -pkin(12) * (mrSges(3,1) * t1273 + mrSges(3,2) * t1276) + (Ifges(3,1) * t1276 - t1391) * t1413;
t1070 = qJD(1) * t1295 - t1196 * t1272 - t1197 * t1275;
t1087 = t1172 * Ifges(4,2) + t1261 * Ifges(4,6) - t1379;
t1161 = Ifges(4,4) * t1172;
t1088 = -t1173 * Ifges(4,1) + t1261 * Ifges(4,5) + t1161;
t1157 = -qJDD(1) * t1246 - t1251 * t1373;
t1158 = -qJDD(1) * t1251 + t1261 * t1350;
t1230 = sin(t1257);
t1231 = sin(t1258);
t1286 = Ifges(9,5) * t1157 + Ifges(9,6) * t1158 + t1472 * qJD(1) + (Ifges(4,3) + Ifges(9,3)) * t1260 - t1410 * t1181 + (-t1435 - t1480) * t1230 / 0.2e1 + (-t1435 + t1480) * t1231 / 0.2e1 + Ifges(4,5) * t1070 + Ifges(4,6) * t1071 + t1474 * t1393 * t1407 + (mrSges(9,1) * cos(t1093) + mrSges(9,2) * sin(t1093)) * g(3) + t1334 * t1394 - t1261 * (Ifges(4,5) * t1172 + Ifges(4,6) * t1173) / 0.2e1 + t1087 * t1418 + t1482 + t1180 * mrSges(4,1) + t1488 + t1173 * (Ifges(4,1) * t1172 + t1379) / 0.2e1 - (t1382 - t1385) * t1373 / 0.2e1 + (t1301 + t1300) * t1284 / 0.2e1 - (Ifges(4,2) * t1173 + t1088 + t1161) * t1172 / 0.2e1;
t1281 = mrSges(8,1) * g(1);
t1280 = mrSges(3,2) * g(2);
t1250 = sin(t1267);
t1249 = sin(t1266);
t1205 = g(1) * t1448;
t1203 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t1198 = t1453 * t1449;
t1194 = 0.2e1 * t1402;
t1177 = t1445 - mrSges(2,1) + (-m(3) - t1243) * pkin(12);
t1160 = qJ(1) - t1375;
t1159 = qJ(1) + t1375;
t1152 = -qJ(1) + t1328;
t1151 = qJ(1) + t1328;
t1139 = pkin(1) * t1345 - t1405;
t1105 = t1113 * qJD(2);
t1104 = t1114 * qJD(2);
t1094 = -pkin(4) * t1121 + t1273 * t1407;
t1061 = qJDD(1) * t1074 + t1211;
t1057 = qJD(1) * t1107 + qJDD(1) * t1117;
t1056 = -qJD(1) * t1106 + qJDD(1) * t1319;
t1051 = Ifges(7,5) * qJD(2) + (t1381 + t1467) * qJD(1);
t1050 = Ifges(7,6) * qJD(2) + (t1380 + t1478) * qJD(1);
t1049 = (qJD(2) * t1105 + qJDD(2) * t1114) * pkin(1);
t1048 = (qJD(2) * t1104 - qJDD(2) * t1113) * pkin(1);
t1041 = (mrSges(6,1) * t1271 + mrSges(6,2) * t1274) * t1080;
t1038 = -pkin(8) * t1092 - pkin(10) * t1091 + t1303;
t1018 = -mrSges(6,1) * t1043 + mrSges(6,2) * t1042;
t1013 = t1236 * t1485 - t1331;
t1006 = t1007 * t1274 + t1139 * t1271;
t1005 = -t1007 * t1271 + t1139 * t1274;
t1001 = -t1029 * t1333 + t1492;
t1000 = t1001 * t1274 - t1271 * t1405;
t999 = -t1001 * t1271 - t1274 * t1405;
t1 = [pkin(6) * (-mrSges(7,1) * t1056 + mrSges(7,2) * t1057) + t1074 * (-mrSges(9,1) * t1158 + mrSges(9,2) * t1157) + (g(1) * t1203 + g(2) * t1177) * cos(qJ(1)) + (Ifges(7,5) * t1423 + Ifges(7,6) * t1424 + (t1386 / 0.2e1 - t1383 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t1276 - Ifges(3,2) * t1273) * t1411 + t1289) * qJD(1) + (-t1474 * t1261 * mrSges(9,3) + (t1120 * t1275 - t1121 * t1272) * mrSges(4,3) + (t1491 * qJD(1) - t1324) * t1273) * pkin(1) + t1299) * qJD(2) + (Ifges(4,5) * t1422 + Ifges(4,6) * t1421 + (-t1385 / 0.2e1 + t1382 / 0.2e1) * t1261 + (t1464 - t1300 / 0.2e1 - t1301 / 0.2e1) * qJD(1) - t1472) * t1261 - t1157 * t1392 + (t1071 * t1190 + t1121 * t1419) * Ifges(4,2) - (g(1) * t1177 - g(2) * t1203) * sin(qJ(1)) + t1196 * t1377 + t1056 * t1380 + t1057 * t1381 + (-t1157 * t1251 - t1158 * t1246) * Ifges(9,4) + (-t1309 + t1459) * t1094 - t1470 + t1087 * t1421 + t1088 * t1422 + t1051 * t1423 + t1050 * t1424 + t1228 * (-mrSges(4,1) * t1071 + mrSges(4,2) * t1070) + (t1194 + 0.2e1 * t1280) * t1249 / 0.4e1 + (t1194 - 0.2e1 * t1280) * t1250 / 0.4e1 + (t1281 - t1429) * sin(t1160) / 0.2e1 + (t1281 + t1429) * sin(t1159) / 0.2e1 + (-t1198 + 0.2e1 * t1435) * t1231 / 0.4e1 + (-t1198 - 0.2e1 * t1435) * t1230 / 0.4e1 + (-t1430 - t1440) * cos(t1160) / 0.2e1 + (t1430 - t1440) * cos(t1159) / 0.2e1 + (t1058 * mrSges(5,2) + t995 * mrSges(5,3) + Ifges(5,1) * t1078 + Ifges(5,4) * t1079 + (t995 * mrSges(6,2) - t994 * mrSges(6,3) + Ifges(6,1) * t1042 + Ifges(6,4) * t1043 + Ifges(6,5) * t1073) * t1274 + (t995 * mrSges(6,1) - t993 * mrSges(6,3) - Ifges(6,4) * t1042 - Ifges(6,2) * t1043 - Ifges(6,6) * t1073) * t1271 + t1450 * qJD(4)) * t1091 + (-t1070 * t1304 + t1120 * t1418) * Ifges(4,1) + (t1070 * t1190 - t1071 * t1304 + t1120 * t1419 + t1121 * t1418) * Ifges(4,4) + (t1180 * t1304 + t1181 * t1190) * mrSges(4,3) + (-mrSges(4,1) * t1190 - mrSges(4,2) * t1304 + t1317 + t1330) * t1337 - pkin(12) * (-mrSges(3,1) * t1196 + mrSges(3,2) * t1197) + (t1196 * t1273 + t1197 * t1276) * Ifges(3,4) + (t1061 * t1074 + ((-cos(t1151) / 0.2e1 - cos(t1152) / 0.2e1) * g(2) + (sin(t1151) / 0.2e1 - sin(t1152) / 0.2e1) * g(1)) * pkin(3)) * m(9) + t1482 + (t1409 / 0.4e1 + t1335) * g(2) * t1448 + (-t1409 / 0.4e1 + t1335) * t1242 * t1449 - t1488 + t1061 * t1316 + (t1310 * qJD(4) + t1033 * t1274 + t1034 * t1271) * t1038 + (Ifges(3,5) * t1273 + Ifges(7,5) * t1117 + Ifges(3,6) * t1276 + Ifges(7,6) * t1319) * qJDD(2) + (t1117 * t1056 + t1057 * t1319 + (t1117 * t1424 + t1319 * t1423) * qJD(1)) * Ifges(7,4) + (pkin(6) * (mrSges(7,1) * t1106 + mrSges(7,2) * t1107) + t1380 * t1424 + t1381 * t1423 + t1228 * (-mrSges(4,1) * t1121 + mrSges(4,2) * t1120)) * qJD(1) + ((qJD(4) * t1456 + t1271 * t993 + t1274 * t994) * t1038 + t1312 * t1094) * m(6) - t1158 * t1384 + (m(5) * t1058 - mrSges(5,1) * t1079 + mrSges(5,2) * t1078) * t1303 + (t1246 * t1180 - t1181 * t1251) * mrSges(9,3) + (-t1058 * mrSges(5,1) + t996 * mrSges(5,3) + Ifges(5,4) * t1078 + Ifges(5,2) * t1079 - t1451) * t1092 + (t1048 * t1294 - t1049 * t1187) * mrSges(8,3) + (Ifges(2,3) + (t1228 * mrSges(8,1) + Ifges(8,2) * t1187) * t1187 - (t1228 * mrSges(8,2) - Ifges(8,1) * t1294 - 0.2e1 * Ifges(8,4) * t1187) * t1294 + (m(3) * pkin(12) + mrSges(3,1) * t1276 - t1397) * pkin(12) + (-mrSges(7,1) * t1319 + mrSges(7,2) * t1117 + t1445) * pkin(6)) * qJDD(1) + (-Ifges(4,5) * t1304 - Ifges(9,5) * t1246 + Ifges(4,6) * t1190 - Ifges(9,6) * t1251) * t1260 + t1197 * t1273 * Ifges(3,1) + t1293 + t1291 - t1292 + (-t1205 + 0.2e1 * t1399) * sin(t1134) / 0.4e1 + (-t1205 - 0.2e1 * t1399) * sin(t1135) / 0.4e1; Ifges(7,6) * t1056 + Ifges(7,5) * t1057 + t1012 * t1018 + (-mrSges(3,1) * g(1) + t1280) * t1250 / 0.2e1 + t1286 - t1139 * t1046 + (t1051 * t1468 + t1117 * t1050 / 0.2e1 - t1276 * t1237 / 0.2e1 + (-pkin(6) * (mrSges(7,1) * t1117 + mrSges(7,2) * t1319) - t1117 * (Ifges(7,1) * t1319 - t1478) / 0.2e1 + (-Ifges(7,2) * t1117 + t1467) * t1468 + t1377 * t1413 - t1289) * qJD(1) - t1299 + t1306 - (Ifges(7,5) * t1319 - Ifges(7,6) * t1117 - t1383 + t1386) * qJD(2) / 0.2e1) * qJD(1) + (Ifges(7,3) + Ifges(3,3)) * qJDD(2) + t1290 * t1013 + (-t1271 * t998 - t1005) * t1045 + Ifges(3,6) * t1196 + Ifges(3,5) * t1197 + t1374 * t1041 + (t1012 * t1078 + t1013 * t1079 + t1374 * t1080 + t1477 * t1081) * mrSges(5,3) + (t1274 * t998 - t1006) * t1044 + (-mrSges(7,1) * cos(t1256) + t1397 - mrSges(7,2) * sin(t1256) - t1200 * t1276 + t1471) * g(3) + t1291 + t1292 + (t1280 + t1402) * t1249 / 0.2e1 + (t1456 * t998 + (-qJD(4) * t1312 - t1271 * t994 + t1274 * t993) * t1013 - t1003 * t1005 - t1004 * t1006 + t1475) * m(6) + (t1477 * t1009 + t1013 * t996 - t1136 * t1139 + t1475) * m(5) + (-t1243 * g(1) * t1250 / 0.2e1 + t1378 + (-t1113 * t1294 - t1114 * t1187) * qJDD(1) * mrSges(8,3) + (-mrSges(4,3) * t1071 - mrSges(9,3) * t1158 + qJD(3) * t1338 + t1260 * t1410) * t1272 + (mrSges(9,3) * t1157 - t1481 * t1260 + t1339 * qJD(3) + (-qJD(2) * t1172 + t1070) * mrSges(4,3)) * t1275 + (t1324 * t1273 + (t1104 * t1294 - t1105 * t1187 + (t1113 * t1187 - t1114 * t1294) * qJD(2)) * mrSges(8,3)) * qJD(1) + m(4) * (-t1180 * t1275 - t1355) + 0.2e1 * (-t1048 * t1113 + t1049 * t1114) * t1447 + 0.2e1 * t1355 * t1469 + 0.2e1 * (t1298 * t1275 * t1469 + (-t1104 * t1113 + t1105 * t1114) * qJD(2) * t1447) * pkin(1) - t1491 * t1273 * t1284) * pkin(1); (-t1001 * t1081 - t1002 * t1080) * mrSges(5,3) - m(6) * (t1000 * t1004 + t1003 * t999 + t1370) + t1286 + (-t1017 * t1041 - t1029 * t1018 + t1290 * t1485 + m(5) * (t1485 * t996 + t1313) + m(6) * (t1313 + ((-t1465 - t994) * t1271 + (-t1466 + t993) * t1274) * t1485) + (-t1017 * t1080 - t1029 * t1078 + t1079 * t1485) * mrSges(5,3) + t1459 * t1173 + (m(5) * t1009 + m(6) * t1456 + t1081 * mrSges(5,3) + t1310) * t1016) * pkin(4) - m(5) * (t1001 * t1009 + t1370) - t1002 * t1041 - t1000 * t1044 - t999 * t1045 + t1471 * g(3) + t1306 * qJD(1) + (t1378 + (-t1338 * t1272 + (-t1339 - t1395) * t1275) * qJD(2)) * pkin(1); -t1003 * t1044 + t1004 * t1045 + ((-sin(t1133) / 0.2e1 - sin(t1132) / 0.2e1) * mrSges(6,2) + (cos(t1132) / 0.2e1 - cos(t1133) / 0.2e1) * mrSges(6,1)) * g(3) - t1450 * t1080 + t1293 + t1451 + t1470;];
tau = t1;
