% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
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
% tauJB [(6+13)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = palh1m1OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynJB_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:30:35
% EndTime: 2020-04-15 19:31:06
% DurationCPUTime: 9.57s
% Computational Cost: add. (101000->624), mult. (221825->791), div. (0->0), fcn. (170408->22), ass. (0->250)
t1443 = sin(qJ(1));
t1452 = cos(qJ(1));
t1406 = -t1452 * g(1) - t1443 * g(2);
t1453 = qJD(1) ^ 2;
t1392 = -t1453 * pkin(15) + t1406;
t1442 = sin(qJ(2));
t1451 = cos(qJ(2));
t1369 = -t1451 * g(3) - t1442 * t1392;
t1371 = t1442 * g(3) - t1451 * t1392;
t1486 = qJD(1) * t1451;
t1478 = qJD(2) * t1486;
t1397 = -t1442 * qJDD(1) - t1478;
t1488 = qJD(1) * t1442;
t1399 = -qJD(2) * t1488 + t1451 * qJDD(1);
t1357 = (-t1442 * t1451 * t1453 + qJDD(2)) * pkin(1) + t1371;
t1358 = (-t1442 ^ 2 * t1453 - qJD(2) ^ 2) * pkin(1) + t1369;
t1441 = sin(qJ(3));
t1450 = cos(qJ(3));
t1297 = -t1450 * t1357 + t1441 * t1358;
t1299 = t1441 * t1357 + t1450 * t1358;
t1384 = (-t1441 * t1442 + t1450 * t1451) * qJD(1);
t1319 = t1384 * qJD(3) - t1450 * t1397 + t1441 * t1399;
t1381 = (t1441 * t1451 + t1442 * t1450) * qJD(1);
t1322 = -t1381 * qJD(3) + t1441 * t1397 + t1450 * t1399;
t1430 = qJD(2) + qJD(3);
t1335 = Ifges(4,4) * t1381 + Ifges(4,2) * t1384 + Ifges(4,6) * t1430;
t1338 = Ifges(4,1) * t1381 + Ifges(4,4) * t1384 + Ifges(4,5) * t1430;
t1427 = qJDD(2) + qJDD(3);
t1440 = sin(qJ(4));
t1449 = cos(qJ(4));
t1346 = t1449 * t1381 + t1440 * t1384;
t1245 = -t1346 * qJD(4) - t1440 * t1319 + t1449 * t1322;
t1344 = -t1440 * t1381 + t1449 * t1384;
t1247 = t1344 * qJD(4) + t1449 * t1319 + t1440 * t1322;
t1405 = t1443 * g(1) - t1452 * g(2);
t1390 = -qJDD(1) * pkin(15) - t1405;
t1353 = t1390 + (-t1397 + t1478) * pkin(1);
t1266 = t1353 + (t1381 * t1430 - t1322) * pkin(5);
t1421 = qJD(4) + t1430;
t1208 = (-t1344 * t1421 - t1247) * pkin(11) + (t1346 * t1421 - t1245) * pkin(9) + t1266;
t1261 = (t1381 * t1384 + t1427) * pkin(5) + t1299;
t1270 = (-t1384 ^ 2 - t1430 ^ 2) * pkin(5) + t1297;
t1219 = t1440 * t1261 + t1449 * t1270;
t1284 = -t1344 * pkin(9) - t1346 * pkin(11);
t1416 = t1421 ^ 2;
t1418 = qJDD(4) + t1427;
t1210 = -t1416 * pkin(9) + t1418 * pkin(11) + t1344 * t1284 + t1219;
t1439 = sin(qJ(5));
t1448 = cos(qJ(5));
t1197 = t1439 * t1208 + t1448 * t1210;
t1218 = t1449 * t1261 - t1440 * t1270;
t1209 = -t1418 * pkin(9) - t1416 * pkin(11) + t1346 * t1284 - t1218;
t1306 = t1448 * t1346 + t1439 * t1421;
t1215 = -t1306 * qJD(5) - t1439 * t1247 + t1448 * t1418;
t1305 = -t1439 * t1346 + t1448 * t1421;
t1216 = t1305 * qJD(5) + t1448 * t1247 + t1439 * t1418;
t1339 = qJD(5) - t1344;
t1232 = Ifges(6,5) * t1306 + Ifges(6,6) * t1305 + Ifges(6,3) * t1339;
t1234 = Ifges(6,1) * t1306 + Ifges(6,4) * t1305 + Ifges(6,5) * t1339;
t1242 = qJDD(5) - t1245;
t1169 = -mrSges(6,1) * t1209 + mrSges(6,3) * t1197 + Ifges(6,4) * t1216 + Ifges(6,2) * t1215 + Ifges(6,6) * t1242 - t1306 * t1232 + t1339 * t1234;
t1196 = t1448 * t1208 - t1439 * t1210;
t1233 = Ifges(6,4) * t1306 + Ifges(6,2) * t1305 + Ifges(6,6) * t1339;
t1170 = mrSges(6,2) * t1209 - mrSges(6,3) * t1196 + Ifges(6,1) * t1216 + Ifges(6,4) * t1215 + Ifges(6,5) * t1242 + t1305 * t1232 - t1339 * t1233;
t1274 = Ifges(5,4) * t1346 + Ifges(5,2) * t1344 + Ifges(5,6) * t1421;
t1276 = Ifges(5,1) * t1346 + Ifges(5,4) * t1344 + Ifges(5,5) * t1421;
t1267 = -t1339 * mrSges(6,2) + t1305 * mrSges(6,3);
t1268 = t1339 * mrSges(6,1) - t1306 * mrSges(6,3);
t1467 = -m(6) * t1209 + t1215 * mrSges(6,1) - t1216 * mrSges(6,2) + t1305 * t1267 - t1306 * t1268;
t1264 = -t1305 * mrSges(6,1) + t1306 * mrSges(6,2);
t1185 = m(6) * t1196 + t1242 * mrSges(6,1) - t1216 * mrSges(6,3) - t1306 * t1264 + t1339 * t1267;
t1186 = m(6) * t1197 - t1242 * mrSges(6,2) + t1215 * mrSges(6,3) + t1305 * t1264 - t1339 * t1268;
t1475 = -t1439 * t1185 + t1448 * t1186;
t1461 = -pkin(9) * t1467 - pkin(11) * t1475 - mrSges(5,1) * t1218 + mrSges(5,2) * t1219 - Ifges(5,5) * t1247 - Ifges(5,6) * t1245 - Ifges(5,3) * t1418 - t1448 * t1169 - t1439 * t1170 - t1346 * t1274 + t1344 * t1276;
t1283 = -t1344 * mrSges(5,1) + t1346 * mrSges(5,2);
t1328 = t1421 * mrSges(5,1) - t1346 * mrSges(5,3);
t1162 = m(5) * t1219 - t1418 * mrSges(5,2) + t1245 * mrSges(5,3) + t1344 * t1283 - t1421 * t1328 + t1475;
t1326 = -t1421 * mrSges(5,2) + t1344 * mrSges(5,3);
t1177 = m(5) * t1218 + t1418 * mrSges(5,1) - t1247 * mrSges(5,3) - t1346 * t1283 + t1421 * t1326 + t1467;
t1482 = t1440 * t1162 + t1449 * t1177;
t1457 = -pkin(5) * t1482 - mrSges(4,1) * t1299 + mrSges(4,2) * t1297 - Ifges(4,5) * t1319 - Ifges(4,6) * t1322 - Ifges(4,3) * t1427 - t1381 * t1335 + t1384 * t1338 + t1461;
t1437 = sin(qJ(7));
t1446 = cos(qJ(7));
t1296 = t1446 * t1357 - t1437 * t1358;
t1298 = t1437 * t1357 + t1446 * t1358;
t1383 = (-t1437 * t1442 + t1446 * t1451) * qJD(1);
t1318 = -t1383 * qJD(7) + t1446 * t1397 - t1437 * t1399;
t1380 = (-t1437 * t1451 - t1442 * t1446) * qJD(1);
t1321 = t1380 * qJD(7) + t1437 * t1397 + t1446 * t1399;
t1429 = qJD(2) + qJD(7);
t1334 = Ifges(8,4) * t1383 + Ifges(8,2) * t1380 + Ifges(8,6) * t1429;
t1337 = Ifges(8,1) * t1383 + Ifges(8,4) * t1380 + Ifges(8,5) * t1429;
t1426 = qJDD(2) + qJDD(7);
t1432 = sin(pkin(19));
t1433 = cos(pkin(19));
t1329 = (-t1380 * t1433 - t1383 * t1432) * pkin(4);
t1424 = t1429 ^ 2;
t1248 = -t1383 * t1329 + (t1424 * t1432 + t1426 * t1433) * pkin(4) + t1296;
t1249 = t1380 * t1329 + (-t1424 * t1433 + t1426 * t1432) * pkin(4) + t1298;
t1434 = cos(qJ(10));
t1492 = sin(qJ(10));
t1386 = t1432 * t1434 - t1433 * t1492;
t1387 = -t1432 * t1492 - t1433 * t1434;
t1212 = t1387 * t1248 - t1386 * t1249;
t1213 = t1386 * t1248 + t1387 * t1249;
t1301 = t1386 * t1380 + t1387 * t1383;
t1225 = -t1301 * qJD(10) + t1387 * t1318 - t1386 * t1321;
t1300 = t1387 * t1380 - t1386 * t1383;
t1226 = t1300 * qJD(10) + t1386 * t1318 + t1387 * t1321;
t1419 = qJD(10) + t1429;
t1251 = Ifges(11,4) * t1301 + Ifges(11,2) * t1300 + Ifges(11,6) * t1419;
t1252 = Ifges(11,1) * t1301 + Ifges(11,4) * t1300 + Ifges(11,5) * t1419;
t1412 = qJDD(10) + t1426;
t1464 = -mrSges(11,1) * t1212 + mrSges(11,2) * t1213 - Ifges(11,5) * t1226 - Ifges(11,6) * t1225 - Ifges(11,3) * t1412 - t1301 * t1251 + t1300 * t1252;
t1257 = -t1300 * mrSges(11,1) + t1301 * mrSges(11,2);
t1292 = -t1419 * mrSges(11,2) + t1300 * mrSges(11,3);
t1202 = m(11) * t1212 + t1412 * mrSges(11,1) - t1226 * mrSges(11,3) - t1301 * t1257 + t1419 * t1292;
t1293 = t1419 * mrSges(11,1) - t1301 * mrSges(11,3);
t1203 = m(11) * t1213 - t1412 * mrSges(11,2) + t1225 * mrSges(11,3) + t1300 * t1257 - t1419 * t1293;
t1480 = -t1386 * t1202 + t1387 * t1203;
t1481 = t1387 * t1202 + t1386 * t1203;
t1490 = pkin(4) * t1433;
t1491 = pkin(4) * t1432;
t1458 = -mrSges(8,1) * t1296 + mrSges(8,2) * t1298 - Ifges(8,5) * t1321 - Ifges(8,6) * t1318 - Ifges(8,3) * t1426 - t1383 * t1334 + t1380 * t1337 - t1480 * t1491 - t1481 * t1490 + t1464;
t1436 = sin(qJ(8));
t1445 = cos(qJ(8));
t1382 = (-t1436 * t1442 + t1445 * t1451) * qJD(1);
t1317 = -t1382 * qJD(8) + t1445 * t1397 - t1436 * t1399;
t1379 = (-t1436 * t1451 - t1442 * t1445) * qJD(1);
t1320 = t1379 * qJD(8) + t1436 * t1397 + t1445 * t1399;
t1323 = -t1436 * t1369 + t1445 * t1371;
t1324 = t1445 * t1369 + t1436 * t1371;
t1428 = qJD(2) + qJD(8);
t1333 = Ifges(9,4) * t1382 + Ifges(9,2) * t1379 + Ifges(9,6) * t1428;
t1336 = Ifges(9,1) * t1382 + Ifges(9,4) * t1379 + Ifges(9,5) * t1428;
t1425 = qJDD(2) + qJDD(8);
t1277 = (t1379 * t1382 + t1425) * pkin(2) + t1323;
t1289 = (-t1379 ^ 2 - t1428 ^ 2) * pkin(2) + t1324;
t1435 = sin(qJ(9));
t1444 = cos(qJ(9));
t1230 = -t1444 * t1277 + t1435 * t1289;
t1231 = -t1435 * t1277 - t1444 * t1289;
t1345 = -t1435 * t1379 - t1444 * t1382;
t1244 = -t1345 * qJD(9) - t1444 * t1317 + t1435 * t1320;
t1343 = -t1444 * t1379 + t1435 * t1382;
t1246 = t1343 * qJD(9) - t1435 * t1317 - t1444 * t1320;
t1420 = qJD(9) + t1428;
t1273 = Ifges(10,4) * t1345 + Ifges(10,2) * t1343 + Ifges(10,6) * t1420;
t1275 = Ifges(10,1) * t1345 + Ifges(10,4) * t1343 + Ifges(10,5) * t1420;
t1417 = qJDD(9) + t1425;
t1465 = -mrSges(10,1) * t1230 + mrSges(10,2) * t1231 - Ifges(10,5) * t1246 - Ifges(10,6) * t1244 - Ifges(10,3) * t1417 - t1345 * t1273 + t1343 * t1275;
t1282 = -t1343 * mrSges(10,1) + t1345 * mrSges(10,2);
t1325 = -t1420 * mrSges(10,2) + t1343 * mrSges(10,3);
t1206 = m(10) * t1230 + t1417 * mrSges(10,1) - t1246 * mrSges(10,3) - t1345 * t1282 + t1420 * t1325;
t1327 = t1420 * mrSges(10,1) - t1345 * mrSges(10,3);
t1207 = m(10) * t1231 - t1417 * mrSges(10,2) + t1244 * mrSges(10,3) + t1343 * t1282 - t1420 * t1327;
t1473 = -t1444 * t1206 - t1435 * t1207;
t1459 = -pkin(2) * t1473 - mrSges(9,1) * t1323 + mrSges(9,2) * t1324 - Ifges(9,5) * t1320 - Ifges(9,6) * t1317 - Ifges(9,3) * t1425 - t1382 * t1333 + t1379 * t1336 + t1465;
t1350 = -t1384 * mrSges(4,1) + t1381 * mrSges(4,2);
t1364 = t1430 * mrSges(4,1) - t1381 * mrSges(4,3);
t1154 = m(4) * t1297 - t1427 * mrSges(4,2) + t1322 * mrSges(4,3) + t1449 * t1162 - t1440 * t1177 + t1384 * t1350 - t1430 * t1364;
t1367 = -t1430 * mrSges(4,2) + t1384 * mrSges(4,3);
t1155 = m(4) * t1299 + t1427 * mrSges(4,1) - t1319 * mrSges(4,3) - t1381 * t1350 + t1430 * t1367 + t1482;
t1352 = -t1380 * mrSges(8,1) + t1383 * mrSges(8,2);
t1363 = -t1429 * mrSges(8,2) + t1380 * mrSges(8,3);
t1174 = m(8) * t1296 + t1426 * mrSges(8,1) - t1321 * mrSges(8,3) - t1383 * t1352 + t1429 * t1363 + t1481;
t1366 = t1429 * mrSges(8,1) - t1383 * mrSges(8,3);
t1175 = m(8) * t1298 - t1426 * mrSges(8,2) + t1318 * mrSges(8,3) + t1380 * t1352 - t1429 * t1366 + t1480;
t1483 = -t1450 * t1154 + t1441 * t1155 + t1446 * t1174 + t1437 * t1175;
t1493 = -t1483 * pkin(1) - mrSges(3,1) * t1371 + mrSges(3,2) * t1369 - Ifges(3,5) * t1399 - Ifges(3,6) * t1397 - Ifges(3,3) * qJDD(2) + t1457 + t1458 + t1459;
t1438 = sin(qJ(6));
t1489 = qJD(1) * t1438;
t1447 = cos(qJ(6));
t1487 = qJD(1) * t1447;
t1401 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t1489;
t1485 = t1438 * t1401;
t1351 = -t1379 * mrSges(9,1) + t1382 * mrSges(9,2);
t1362 = -t1428 * mrSges(9,2) + t1379 * mrSges(9,3);
t1189 = m(9) * t1323 + t1425 * mrSges(9,1) - t1320 * mrSges(9,3) - t1382 * t1351 + t1428 * t1362 + t1473;
t1365 = t1428 * mrSges(9,1) - t1382 * mrSges(9,3);
t1190 = m(9) * t1324 - t1425 * mrSges(9,2) + t1317 * mrSges(9,3) + t1435 * t1206 - t1444 * t1207 + t1379 * t1351 - t1428 * t1365;
t1395 = (mrSges(3,1) * t1442 + mrSges(3,2) * t1451) * qJD(1);
t1404 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1486;
t1144 = m(3) * t1369 - qJDD(2) * mrSges(3,2) + t1397 * mrSges(3,3) - qJD(2) * t1404 + t1441 * t1154 + t1450 * t1155 - t1437 * t1174 + t1446 * t1175 - t1436 * t1189 + t1445 * t1190 - t1395 * t1488;
t1402 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t1488;
t1145 = m(3) * t1371 + qJDD(2) * mrSges(3,1) - t1399 * mrSges(3,3) + qJD(2) * t1402 + t1445 * t1189 + t1436 * t1190 - t1395 * t1486 + t1483;
t1393 = t1453 * pkin(14) + t1406;
t1368 = -t1447 * g(3) - t1438 * t1393;
t1394 = (-mrSges(7,1) * t1447 + mrSges(7,2) * t1438) * qJD(1);
t1396 = qJD(6) * t1487 + t1438 * qJDD(1);
t1403 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t1487;
t1294 = m(7) * t1368 + qJDD(6) * mrSges(7,1) - t1396 * mrSges(7,3) + qJD(6) * t1403 - t1394 * t1489;
t1370 = -t1438 * g(3) + t1447 * t1393;
t1398 = -qJD(6) * t1489 + t1447 * qJDD(1);
t1295 = m(7) * t1370 - qJDD(6) * mrSges(7,2) + t1398 * mrSges(7,3) - qJD(6) * t1401 + t1394 * t1487;
t1474 = -t1438 * t1294 + t1447 * t1295;
t1141 = m(2) * t1406 - t1453 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t1442 * t1144 - t1451 * t1145 + t1474;
t1228 = (-t1318 * t1433 - t1321 * t1432 + (-t1380 * t1432 + t1383 * t1433) * t1429) * pkin(4) + t1353;
t1195 = m(11) * t1228 - t1225 * mrSges(11,1) + t1226 * mrSges(11,2) - t1300 * t1292 + t1301 * t1293;
t1164 = t1448 * t1185 + t1439 * t1186;
t1468 = m(5) * t1266 - t1245 * mrSges(5,1) + t1247 * mrSges(5,2) - t1344 * t1326 + t1346 * t1328 + t1164;
t1456 = t1384 * t1367 + t1380 * t1363 + t1322 * mrSges(4,1) + t1318 * mrSges(8,1) - t1468 + (-m(4) - m(8)) * t1353 - t1381 * t1364 - t1383 * t1366 - t1321 * mrSges(8,2) - t1319 * mrSges(4,2) - t1195;
t1281 = (t1382 * t1428 - t1317) * pkin(2) + t1390;
t1466 = m(10) * t1281 - t1244 * mrSges(10,1) + t1246 * mrSges(10,2) - t1343 * t1325 + t1345 * t1327;
t1455 = t1456 + t1379 * t1362 + t1397 * mrSges(3,1) + t1317 * mrSges(9,1) + (-m(3) - m(9)) * t1390 - t1399 * mrSges(3,2) - t1382 * t1365 - t1320 * mrSges(9,2) - t1466;
t1391 = qJDD(1) * pkin(14) - t1405;
t1469 = -m(7) * t1391 + t1398 * mrSges(7,1) - t1396 * mrSges(7,2) + t1403 * t1487;
t1470 = -t1442 * t1402 - t1451 * t1404;
t1151 = t1455 + qJDD(1) * mrSges(2,1) + (t1470 - t1485) * qJD(1) - t1453 * mrSges(2,2) + m(2) * t1405 + t1469;
t1484 = t1443 * t1141 + t1452 * t1151;
t1479 = t1447 * t1294 + t1438 * t1295;
t1477 = t1452 * t1141 - t1443 * t1151;
t1476 = t1451 * t1144 - t1442 * t1145;
t1375 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t1438 + Ifges(7,2) * t1447) * qJD(1);
t1377 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t1438 + Ifges(7,4) * t1447) * qJD(1);
t1472 = t1438 * t1375 - t1447 * t1377;
t1376 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t1451 - Ifges(3,2) * t1442) * qJD(1);
t1378 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t1451 - Ifges(3,4) * t1442) * qJD(1);
t1471 = t1451 * t1376 + t1442 * t1378;
t1463 = mrSges(7,1) * t1368 - mrSges(7,2) * t1370 + Ifges(7,5) * t1396 + Ifges(7,6) * t1398 + Ifges(7,3) * qJDD(6);
t1272 = Ifges(5,5) * t1346 + Ifges(5,6) * t1344 + Ifges(5,3) * t1421;
t1148 = -pkin(11) * t1164 + mrSges(5,2) * t1266 - mrSges(5,3) * t1218 + Ifges(5,1) * t1247 + Ifges(5,4) * t1245 + Ifges(5,5) * t1418 - t1439 * t1169 + t1448 * t1170 + t1344 * t1272 - t1421 * t1274;
t1460 = mrSges(6,1) * t1196 - mrSges(6,2) * t1197 + Ifges(6,5) * t1216 + Ifges(6,6) * t1215 + Ifges(6,3) * t1242 + t1306 * t1233 - t1305 * t1234;
t1149 = -pkin(9) * t1164 - mrSges(5,1) * t1266 + mrSges(5,3) * t1219 + Ifges(5,4) * t1247 + Ifges(5,2) * t1245 + Ifges(5,6) * t1418 - t1346 * t1272 + t1421 * t1276 - t1460;
t1332 = Ifges(4,5) * t1381 + Ifges(4,6) * t1384 + Ifges(4,3) * t1430;
t1142 = -pkin(5) * t1468 - mrSges(4,1) * t1353 + mrSges(4,3) * t1297 + Ifges(4,4) * t1319 + Ifges(4,2) * t1322 + Ifges(4,6) * t1427 + t1440 * t1148 + t1449 * t1149 - t1381 * t1332 + t1430 * t1338;
t1146 = mrSges(4,2) * t1353 - mrSges(4,3) * t1299 + Ifges(4,1) * t1319 + Ifges(4,4) * t1322 + Ifges(4,5) * t1427 + t1449 * t1148 - t1440 * t1149 + t1384 * t1332 - t1430 * t1335;
t1250 = Ifges(11,5) * t1301 + Ifges(11,6) * t1300 + Ifges(11,3) * t1419;
t1187 = -mrSges(11,1) * t1228 + mrSges(11,3) * t1213 + Ifges(11,4) * t1226 + Ifges(11,2) * t1225 + Ifges(11,6) * t1412 - t1301 * t1250 + t1419 * t1252;
t1188 = mrSges(11,2) * t1228 - mrSges(11,3) * t1212 + Ifges(11,1) * t1226 + Ifges(11,4) * t1225 + Ifges(11,5) * t1412 + t1300 * t1250 - t1419 * t1251;
t1331 = Ifges(8,5) * t1383 + Ifges(8,6) * t1380 + Ifges(8,3) * t1429;
t1158 = -mrSges(8,1) * t1353 + mrSges(8,3) * t1298 + Ifges(8,4) * t1321 + Ifges(8,2) * t1318 + Ifges(8,6) * t1426 + t1387 * t1187 + t1386 * t1188 - t1195 * t1490 - t1383 * t1331 + t1429 * t1337;
t1159 = mrSges(8,2) * t1353 - mrSges(8,3) * t1296 + Ifges(8,1) * t1321 + Ifges(8,4) * t1318 + Ifges(8,5) * t1426 - t1386 * t1187 + t1387 * t1188 - t1195 * t1491 + t1380 * t1331 - t1429 * t1334;
t1271 = Ifges(10,5) * t1345 + Ifges(10,6) * t1343 + Ifges(10,3) * t1420;
t1204 = -mrSges(10,1) * t1281 + mrSges(10,3) * t1231 + Ifges(10,4) * t1246 + Ifges(10,2) * t1244 + Ifges(10,6) * t1417 - t1345 * t1271 + t1420 * t1275;
t1205 = mrSges(10,2) * t1281 - mrSges(10,3) * t1230 + Ifges(10,1) * t1246 + Ifges(10,4) * t1244 + Ifges(10,5) * t1417 + t1343 * t1271 - t1420 * t1273;
t1330 = Ifges(9,5) * t1382 + Ifges(9,6) * t1379 + Ifges(9,3) * t1428;
t1166 = -pkin(2) * t1466 - mrSges(9,1) * t1390 + mrSges(9,3) * t1324 + Ifges(9,4) * t1320 + Ifges(9,2) * t1317 + Ifges(9,6) * t1425 - t1444 * t1204 - t1435 * t1205 - t1382 * t1330 + t1428 * t1336;
t1171 = mrSges(9,2) * t1390 - mrSges(9,3) * t1323 + Ifges(9,1) * t1320 + Ifges(9,4) * t1317 + Ifges(9,5) * t1425 + t1435 * t1204 - t1444 * t1205 + t1379 * t1330 - t1428 * t1333;
t1374 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t1451 - Ifges(3,6) * t1442) * qJD(1);
t1136 = t1456 * pkin(1) - mrSges(3,1) * t1390 + mrSges(3,3) * t1369 + Ifges(3,4) * t1399 + Ifges(3,2) * t1397 + Ifges(3,6) * qJDD(2) + qJD(2) * t1378 + t1441 * t1142 - t1450 * t1146 + t1446 * t1158 + t1437 * t1159 + t1445 * t1166 + t1436 * t1171 - t1374 * t1486;
t1138 = mrSges(3,2) * t1390 - mrSges(3,3) * t1371 + Ifges(3,1) * t1399 + Ifges(3,4) * t1397 + Ifges(3,5) * qJDD(2) - qJD(2) * t1376 + t1450 * t1142 + t1441 * t1146 - t1437 * t1158 + t1446 * t1159 - t1436 * t1166 + t1445 * t1171 - t1374 * t1488;
t1373 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t1438 + Ifges(7,6) * t1447) * qJD(1);
t1262 = -mrSges(7,1) * t1391 + mrSges(7,3) * t1370 + Ifges(7,4) * t1396 + Ifges(7,2) * t1398 + Ifges(7,6) * qJDD(6) + qJD(6) * t1377 - t1373 * t1489;
t1263 = mrSges(7,2) * t1391 - mrSges(7,3) * t1368 + Ifges(7,1) * t1396 + Ifges(7,4) * t1398 + Ifges(7,5) * qJDD(6) - qJD(6) * t1375 + t1373 * t1487;
t1288 = -qJD(1) * t1485 + t1469;
t1462 = -mrSges(2,2) * t1406 - pkin(14) * t1288 - t1442 * t1136 + t1451 * t1138 + pkin(15) * (t1470 * qJD(1) + t1455) + t1438 * t1263 + t1447 * t1262 + mrSges(2,1) * t1405 + Ifges(2,3) * qJDD(1);
t1135 = mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-t1471 - t1472) * qJD(1) + t1453 * Ifges(2,5) - pkin(16) * t1474 - pkin(15) * t1476 + mrSges(2,3) * t1406 + pkin(14) * t1479 - t1463 + t1493;
t1134 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1405 + pkin(16) * t1288 + Ifges(2,5) * qJDD(1) - t1453 * Ifges(2,6) - t1451 * t1136 - t1442 * t1138 - t1438 * t1262 + t1447 * t1263;
t1 = [-m(1) * g(1) + t1477; -m(1) * g(2) + t1484; (-m(1) - m(2)) * g(3) + t1476 + t1479; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(13) * t1484 + t1452 * t1134 - t1443 * t1135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(13) * t1477 + t1443 * t1134 + t1452 * t1135; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1462; t1462; t1471 * qJD(1) - t1493; -t1457; -t1461; t1460; t1472 * qJD(1) + t1463; -t1458; -t1459; -t1465; -t1464; 0; 0; 0;];
tauJB = t1;
