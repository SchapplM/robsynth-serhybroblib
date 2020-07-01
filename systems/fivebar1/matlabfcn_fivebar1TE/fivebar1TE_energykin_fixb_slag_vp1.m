% Calculate kinetic energy for
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1TE_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:13:55
% EndTime: 2020-04-27 09:14:11
% DurationCPUTime: 16.54s
% Computational Cost: add. (124300->616), mult. (375788->982), div. (866->13), fcn. (47089->6), ass. (0->426)
t1374 = (pkin(2) ^ 2);
t1612 = -6 * t1374;
t1330 = sin(qJ(2));
t1333 = cos(qJ(1));
t1540 = t1330 * t1333;
t1586 = pkin(2) * pkin(3);
t1452 = t1540 * t1586;
t1331 = sin(qJ(1));
t1332 = cos(qJ(2));
t1537 = t1332 * t1331;
t1478 = qJD(2) * t1537;
t1610 = -qJD(1) * t1452 - t1478 * t1586;
t1370 = pkin(3) ^ 2;
t1284 = 0.10e2 / 0.3e1 * t1370;
t1368 = t1370 ^ 2;
t1376 = (pkin(1) ^ 2);
t1375 = t1376 ^ 2;
t1514 = t1368 + t1375;
t1357 = 2 * t1376;
t1362 = (pkin(5) ^ 2);
t1517 = t1357 - t1362;
t1529 = t1376 * t1362;
t1361 = t1362 ^ 2;
t1365 = pkin(4) ^ 2;
t1364 = t1365 ^ 2;
t1598 = -t1361 / 0.6e1 + t1364 / 0.6e1;
t1210 = t1517 * t1370 + t1514 - t1529 - t1598;
t1382 = t1374 ^ 2;
t1399 = t1210 + t1382;
t1609 = (t1284 + t1517) * t1612 - 0.6e1 * t1399;
t1611 = -2 * t1376;
t1513 = t1370 - t1376;
t1269 = t1513 * t1374;
t1286 = pkin(3) * t1332;
t1266 = pkin(1) + t1286;
t1567 = qJD(1) * t1331;
t1482 = t1266 * t1567;
t1599 = pkin(2) * t1482 + qJD(2) * t1452;
t1571 = 0.2e1 * pkin(1);
t1568 = 0.4e1 * pkin(3);
t1311 = -t1362 / 0.3e1;
t1601 = t1311 - t1365 / 0.3e1;
t1463 = t1376 + t1601;
t1429 = t1370 + t1463;
t1319 = t1365 / 0.3e1;
t1465 = t1362 / 0.3e1 + t1319 + t1357;
t1608 = 0.10e2 / 0.3e1 * t1382 - 0.2e1 * (-t1370 + t1465) * t1374 + t1429 * t1611;
t1566 = qJD(1) * t1333;
t1272 = pkin(2) * t1566;
t1607 = 0.2e1 * t1272;
t1287 = pkin(2) * t1331;
t1606 = 0.2e1 * t1287;
t1605 = 0.2e1 * t1370;
t1473 = qJD(1) * t1540;
t1402 = t1473 + t1478;
t1604 = -0.4e1 * t1402;
t1541 = t1330 * t1331;
t1490 = pkin(3) * t1541;
t1453 = pkin(2) * t1490;
t1256 = -0.2e1 * t1453;
t1575 = pkin(2) * t1333;
t1506 = pkin(1) * t1575;
t1268 = -0.2e1 * t1506;
t1511 = t1374 + t1376;
t1458 = t1370 + t1511;
t1265 = pkin(1) - t1575;
t1553 = t1265 * t1332;
t1491 = pkin(3) * t1553;
t1183 = t1256 + t1268 + t1458 + 0.2e1 * t1491;
t1224 = -t1537 + t1540;
t1498 = pkin(2) * t1567;
t1263 = pkin(1) * t1498;
t1554 = t1265 * t1330;
t1560 = (t1263 + (-qJD(2) * t1554 + (-qJD(1) * t1224 - t1478) * pkin(2)) * pkin(3)) / t1183 ^ 2;
t1603 = -0.2e1 * t1560;
t1574 = pkin(2) * t1370;
t1297 = t1332 ^ 2;
t1573 = pkin(3) * t1297;
t1602 = qJD(1) - qJD(2);
t1280 = -t1370 / 0.3e1 + t1376;
t1600 = t1610 * t1280;
t1597 = 0.8e1 * pkin(1);
t1296 = t1332 * t1297;
t1596 = -0.8e1 * t1296;
t1595 = -0.2e1 * t1297;
t1594 = 0.4e1 * t1297;
t1593 = -0.4e1 * t1330;
t1592 = 0.4e1 * t1368;
t1591 = -0.8e1 * t1370;
t1346 = 0.6e1 * t1370;
t1590 = 2 * t1374;
t1349 = 5 * t1382;
t1385 = pkin(3) * t1370;
t1589 = 0.4e1 * t1385;
t1281 = t1376 - t1374 / 0.3e1;
t1588 = 0.24e2 * t1281;
t1587 = pkin(1) * pkin(2);
t1585 = -pkin(4) - pkin(5);
t1584 = -pkin(4) + pkin(5);
t1325 = t1370 / 0.3e1;
t1200 = -0.4e1 / 0.9e1 * t1453 + t1376 + t1374 / 0.3e1 + t1325 + t1365 / 0.9e1 - t1362 / 0.9e1;
t1309 = -t1362 / 0.6e1;
t1328 = t1374 / 0.2e1;
t1520 = t1328 + t1376;
t1405 = -t1453 + t1520;
t1208 = t1365 / 0.6e1 + t1309 + t1405;
t1236 = t1311 + t1319 + t1458;
t1530 = t1370 * t1297;
t1534 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t1546 = t1296 * t1385;
t1570 = 0.4e1 * pkin(1);
t1151 = t1280 * t1256 + 0.6e1 * t1200 * t1530 + t1236 * t1534 + (t1208 * t1286 + t1546) * t1570;
t1441 = -0.4e1 * t1453;
t1313 = -0.2e1 / 0.3e1 * t1362;
t1318 = 0.2e1 / 0.3e1 * t1365;
t1461 = t1313 + t1318 + t1357;
t1290 = t1370 + t1376;
t1460 = t1313 + t1290;
t1525 = (t1318 + t1460) * t1290 + t1382;
t1164 = t1236 * t1441 + (t1346 + t1461) * t1374 + t1525;
t1206 = t1256 + t1236;
t1512 = -t1374 + t1376;
t1270 = t1512 * t1370;
t1505 = pkin(1) * t1286;
t1152 = 0.4e1 * t1206 * t1505 + t1270 * t1594 + t1164;
t1221 = t1281 * t1256;
t1533 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t1182 = t1236 * t1533 + t1221;
t1184 = (t1284 + t1461) * t1374 + t1525;
t1267 = 0.2e1 * t1505;
t1292 = -0.3e1 * t1370 + t1376;
t1499 = 0.4e1 * t1530;
t1219 = t1267 + t1499 + t1292;
t1351 = -3 * t1374;
t1293 = t1351 + t1376;
t1504 = 0.8e1 * t1546;
t1455 = pkin(1) * t1504;
t1240 = t1293 * t1455;
t1508 = t1376 - t1362;
t1457 = t1370 + t1508;
t1261 = t1365 + t1457;
t1283 = t1290 ^ 2;
t1300 = t1333 ^ 2;
t1338 = 0.15e2 * t1368;
t1339 = 0.15e2 * t1370;
t1353 = 0.3e1 * t1375;
t1381 = pkin(2) * t1374;
t1371 = t1381 ^ 2;
t1299 = t1333 * t1300;
t1543 = t1299 * t1381;
t1476 = t1266 * t1543;
t1449 = -0.8e1 * t1476;
t1356 = 3 * t1376;
t1515 = -t1362 + t1365;
t1459 = t1356 + t1515;
t1484 = 0.6e1 * t1505;
t1485 = 0.12e2 * t1530;
t1552 = t1266 * t1333;
t1137 = t1219 * t1449 + t1240 + t1182 * t1485 + t1164 * t1484 + t1371 + (t1339 + t1459) * t1382 + t1283 * t1261 + (0.12e2 * t1151 * t1300 + t1346 * t1459 + t1357 * t1515 + t1338 + t1353) * t1374 + 0.6e1 * (-t1152 * t1552 - t1184 * t1490) * pkin(2);
t1347 = 0.3e1 * t1370;
t1276 = t1347 + t1511;
t1241 = t1276 * t1287;
t1572 = t1330 * pkin(3);
t1249 = t1287 - t1572;
t1495 = t1331 * t1574;
t1350 = 3 * t1374;
t1274 = t1350 + t1290;
t1550 = t1274 * t1330;
t1576 = pkin(1) * t1332;
t1165 = t1495 * t1595 + t1241 + (0.2e1 * t1249 * t1576 - t1550) * pkin(3);
t1289 = t1605 + t1374;
t1538 = t1330 * t1385;
t1238 = t1495 - t1538;
t1555 = t1238 * t1297;
t1166 = t1289 * t1572 + 0.2e1 * t1555 + (-t1513 + t1267) * t1287;
t1199 = -pkin(3) * t1550 + t1241;
t1271 = pkin(2) * t1592 + 0.8e1 * t1370 * t1381;
t1472 = t1330 * t1533;
t1202 = t1271 * t1331 + t1472 * t1589;
t1344 = 0.5e1 * t1368;
t1509 = t1375 + t1382;
t1340 = 0.10e2 * t1370;
t1518 = t1340 + t1357;
t1528 = t1376 * t1370;
t1218 = t1518 * t1374 + t1344 + t1509 + 0.6e1 * t1528;
t1226 = t1349 + (t1340 + (6 * t1376)) * t1374 + t1283;
t1448 = 0.8e1 * t1476;
t1542 = t1300 * t1374;
t1503 = -0.4e1 * t1542;
t1143 = t1166 * t1503 + t1202 * t1297 + (-0.4e1 * t1199 * t1576 + (t1226 + t1448) * t1330) * pkin(3) + (0.4e1 * t1165 * t1552 + (-t1218 + t1455) * t1331) * pkin(2);
t1237 = t1265 + t1286;
t1456 = t1374 + t1508;
t1428 = t1370 + t1456;
t1254 = -t1365 + t1428;
t1217 = t1268 + t1254;
t1197 = t1256 + t1217;
t1502 = 0.2e1 * t1542;
t1220 = t1268 + t1502 + t1512;
t1422 = t1220 * t1595 - t1508;
t1496 = pkin(2) * t1541;
t1377 = sqrt(0.4e1 * t1269 * t1300 + 0.4e1 * t1254 * t1506 - t1368 - (t1376 + (pkin(2) - t1585) * (pkin(2) + t1585)) * (t1376 + (pkin(2) - t1584) * (pkin(2) + t1584)) + (t1351 + t1365 + t1422) * t1605 + (-t1197 * t1553 + t1217 * t1496) * t1568);
t1132 = t1137 * t1237 + t1143 * t1377;
t1583 = 0.1e1 / t1132 / 0.4e1;
t1253 = t1347 + t1365 + t1456;
t1198 = t1253 + t1268 + t1441;
t1209 = pkin(2) * t1537 + t1554;
t1255 = -pkin(3) + t1496;
t1139 = -t1198 * t1553 + t1209 * t1377 + (-0.2e1 * pkin(1) * t1255 * t1333 + t1253 * t1541) * pkin(2) + (-t1350 - t1365 - t1370 + t1422 + t1502) * pkin(3);
t1582 = t1139 / 0.2e1;
t1180 = 0.1e1 / t1183;
t1581 = t1180 / 0.2e1;
t1580 = t1332 / 0.2e1;
t1579 = 0.4e1 / 0.3e1 * t1374;
t1578 = -t1377 / 0.4e1;
t1577 = t1377 / 0.4e1;
t1569 = 0.2e1 * pkin(3);
t1565 = qJD(1) * t1374;
t1564 = qJD(2) * t1330;
t1563 = qJD(2) * t1331;
t1562 = qJD(2) * t1332;
t1400 = t1332 * t1402;
t1398 = pkin(3) * t1400;
t1397 = t1265 * t1398;
t1479 = t1330 * t1562;
t1436 = t1220 * t1479;
t1260 = 0.2e1 * t1263;
t1434 = t1333 * t1331 * t1565;
t1423 = -0.4e1 * t1434;
t1212 = t1260 + t1423;
t1557 = t1212 * t1297;
t1561 = ((0.8e1 * t1436 - 0.4e1 * t1557) * t1370 + (-0.4e1 * t1254 * t1587 - 0.8e1 * t1269 * t1333) * t1567 + (t1330 * t1331 ^ 2 * t1565 * t1597 + (0.4e1 * t1197 * t1564 - 0.8e1 * t1263 * t1332) * t1265 + 0.4e1 * (-qJD(1) * t1197 * t1537 + t1217 * t1402 + 0.2e1 * t1397) * pkin(2)) * pkin(3)) / t1377;
t1366 = 0.1e1 / pkin(4);
t1559 = t1180 * t1366;
t1401 = pkin(3) * t1402;
t1201 = pkin(2) * t1401;
t1558 = t1201 * t1332;
t1556 = (t1272 * t1370 - t1385 * t1562) * t1297;
t1288 = -t1362 - t1365;
t1273 = t1356 + t1288;
t1551 = t1273 * t1370;
t1549 = t1283 * (-t1365 + t1457);
t1548 = t1288 * t1376;
t1547 = t1296 * t1368;
t1545 = t1297 * t1385;
t1298 = t1300 ^ 2;
t1544 = t1298 * t1382;
t1539 = t1330 * t1377;
t1536 = t1332 * t1333;
t1535 = t1332 * t1377;
t1532 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t1295 = t1297 ^ 2;
t1531 = t1368 * t1295;
t1395 = t1201 * t1530;
t1394 = t1281 * t1395;
t1433 = t1297 * qJD(2) * t1538;
t1416 = -0.24e2 * t1433;
t1407 = pkin(1) * t1416;
t1527 = t1293 * t1407 - 0.24e2 * t1394;
t1314 = -0.3e1 / 0.2e1 * t1362;
t1524 = t1361 / 0.2e1 - t1364 / 0.2e1;
t1526 = t1290 * ((t1314 + t1357) * t1370 - 0.3e1 / 0.2e1 * t1529 + t1514 + t1524) + t1371;
t1523 = t1309 - t1365 / 0.6e1;
t1323 = -0.2e1 / 0.3e1 * t1365;
t1522 = t1313 + t1323;
t1521 = t1314 + t1356;
t1519 = 0.4e1 / 0.7e1 * t1376 - t1362 / 0.7e1;
t1516 = t1361 - t1364;
t1510 = t1375 - t1368;
t1507 = 0.4e1 * t1286;
t1501 = 0.8e1 * t1531;
t1500 = -0.6e1 * t1530;
t1497 = pkin(2) * t1552;
t1494 = t1385 * t1287;
t1493 = pkin(3) * t1564;
t1492 = pkin(3) * t1562;
t1489 = pkin(3) * t1293 / 0.2e1;
t1488 = 0.16e2 * t1546;
t1487 = -0.12e2 * t1542;
t1486 = -t1572 / 0.2e1;
t1483 = 0.20e2 / 0.3e1 * t1370;
t1481 = qJD(1) * t1536;
t1480 = qJD(2) * t1543;
t1477 = t1180 * t1532;
t1475 = t1330 * t1547;
t1474 = t1330 * t1544;
t1227 = -t1493 + t1498;
t1454 = pkin(1) * t1493;
t1259 = -0.2e1 * t1454;
t1294 = t1330 ^ 2;
t1396 = t1370 * t1400 * t1587;
t1414 = t1300 * t1381 * t1482;
t1406 = -0.24e2 * t1414;
t1411 = t1276 * t1272 - t1274 * t1492;
t1435 = t1370 * t1479;
t1417 = -0.24e2 * t1435;
t1418 = t1272 * t1546;
t1445 = pkin(3) * t1480;
t1420 = t1330 * t1445;
t1421 = pkin(1) * t1433;
t1424 = -0.6e1 * t1435;
t1437 = t1533 * t1562;
t1439 = t1238 * t1479;
t1442 = -0.6e1 * t1454;
t1450 = -0.2e1 * t1479;
t1468 = t1561 / 0.2e1;
t1471 = -((0.8e1 * t1332 * t1266 * t1445 + t1294 * t1480 * t1591 + t1406 * t1572 + (t1289 * t1492 - 0.4e1 * t1439 + 0.2e1 * t1556 + (-t1513 * t1566 + (-qJD(2) * t1541 + t1481) * pkin(1) * t1569) * pkin(2)) * t1503 + 0.8e1 * t1166 * t1434 + 0.4e1 * ((0.4e1 * t1330 * t1478 + t1566 * t1595) * t1574 + ((t1272 - t1492) * t1286 - t1249 * t1493) * t1571 + t1411) * t1497 + t1418 * t1597 + t1407 * t1287 + (t1271 * t1566 + t1437 * t1589) * t1297 + t1202 * t1450 - 0.4e1 * t1411 * t1505 + 0.4e1 * t1199 * t1454 - t1218 * t1272 + t1226 * t1492 - 0.4e1 * t1599 * t1165) * t1377 + t1143 * t1468 + ((t1259 - 0.8e1 * t1435) * t1449 + 0.24e2 * (-0.6e1 * t1421 + 0.3e1 * (-0.4e1 / 0.9e1 * t1478 - 0.4e1 / 0.9e1 * t1473) * t1530 * t1586 + t1200 * t1424 - t1201 * t1267 + t1208 * t1259 + t1600) * t1542 - 0.24e2 * t1151 * t1434 - 0.6e1 * (-0.8e1 * t1270 * t1479 + (t1236 * pkin(2) * t1604 + (-0.4e1 * t1206 * t1564 - 0.8e1 * t1558) * pkin(1)) * pkin(3)) * t1497 + t1182 * t1417 - 0.24e2 * t1236 * t1396 + t1164 * t1442 + t1527 + (0.8e1 * t1420 + 0.24e2 * t1414) * t1219) * t1237 + t1137 * t1227 + 0.6e1 * (t1599 * t1152 + t1610 * t1184) * t1237) / t1132 ^ 2 / 0.4e1;
t1470 = t1139 * t1583;
t1216 = t1268 + t1365 + t1428;
t1168 = t1216 * t1265 + 0.2e1 * t1220 * t1286;
t1172 = t1216 * t1332 + (t1594 - 0.2e1) * t1265 * pkin(3);
t1207 = -t1255 + t1553;
t1141 = t1168 * t1330 + t1172 * t1287 + t1207 * t1377;
t1469 = t1141 * t1583;
t1467 = 0.4e1 / 0.3e1 * t1362 + 0.4e1 / 0.3e1 * t1365 + t1611;
t1355 = 4 * t1376;
t1466 = 0.2e1 / 0.3e1 * t1362 + t1318 + t1355;
t1464 = t1376 + t1523;
t1312 = -t1362 / 0.2e1;
t1462 = t1312 - t1365 / 0.2e1 + t1376;
t1451 = pkin(1) * t1488;
t1247 = t1376 + t1374 / 0.4e1 + t1370 / 0.4e1 - t1362 / 0.8e1;
t1157 = -0.32e2 / 0.21e2 * t1247 * t1453 + t1382 / 0.7e1 + (0.16e2 / 0.21e2 * t1370 + t1519) * t1374 + t1368 / 0.7e1 + t1519 * t1370 + t1375 - 0.3e1 / 0.7e1 * t1529 + t1361 / 0.42e2 - t1364 / 0.42e2;
t1310 = -t1362 / 0.4e1;
t1248 = t1310 + t1325 + t1520;
t1324 = 0.4e1 / 0.3e1 * t1370;
t1161 = -0.8e1 / 0.3e1 * t1248 * t1453 + t1382 / 0.3e1 + (t1324 + t1311) * t1374 + t1375 - t1368 / 0.3e1 + (t1579 + 0.2e1 / 0.3e1 * t1370 + t1313) * t1376 + t1361 / 0.18e2 - t1364 / 0.18e2;
t1326 = t1370 / 0.2e1;
t1213 = -0.2e1 / 0.3e1 * t1453 + t1376 + t1326 + t1310;
t1277 = (t1355 + t1362) * t1370;
t1282 = t1376 - 0.2e1 / 0.3e1 * t1374;
t1258 = t1312 + t1458;
t1415 = t1258 * t1441;
t1142 = t1282 * t1501 + 0.14e2 * t1157 * t1530 + t1280 * t1415 - t1513 * t1382 + (t1277 - 0.10e2 / 0.3e1 * t1368 + 0.2e1 * t1375 - t1529) * t1374 + t1210 * t1534 + (0.6e1 * t1161 * t1286 + t1213 * t1488) * pkin(1);
t1432 = -(3 * t1529) + t1353 + t1524;
t1150 = t1453 * t1609 + (t1338 + (-9 * t1362 + 18 * t1376) * t1370 + t1432) * t1374 + (t1339 + t1521) * t1382 + t1526;
t1167 = t1415 + (t1346 + t1517) * t1374 + t1399;
t1188 = t1258 * t1533 + t1221;
t1144 = t1167 * t1484 + t1188 * t1485 + t1150 + t1240;
t1341 = -2 * t1362;
t1354 = 8 * t1376;
t1205 = t1494 * t1593 + t1592 + (4 * t1374 + t1341 + t1354) * t1370;
t1211 = t1310 - t1370 + t1405;
t1153 = t1256 * t1534 + t1205 * t1297 + t1258 * t1292 + (t1211 * t1507 + t1504) * pkin(1);
t1155 = t1281 * t1415 - t1371 + (-t1284 - t1508) * t1382 + (t1277 + t1510 + t1598) * t1374 + t1210 * t1376;
t1342 = -5 * t1362;
t1343 = 0.7e1 * t1368;
t1345 = 0.7e1 * t1370;
t1160 = (t1345 + t1521) * t1382 + (t1343 + (t1342 + 10 * t1376) * t1370 + t1432) * t1374 + t1526;
t1378 = pkin(1) * t1376;
t1264 = -0.12e2 * pkin(1) * t1385 + t1378 * t1568;
t1279 = -0.8e1 * t1368 + 0.12e2 * t1528;
t1174 = t1264 * t1332 + t1279 * t1297 + t1451 + t1501 + t1514 - 0.6e1 * t1528;
t1191 = t1256 * t1533 + t1258 * t1293;
t1252 = ((t1376 * t1612) + t1509) * t1368;
t1285 = -30 * t1362 + 60 * t1376;
t1431 = -(6 * t1529) + 0.6e1 * t1375 + t1516;
t1127 = -0.32e2 * t1153 * t1476 + 0.16e2 * t1252 * t1295 + 0.24e2 * t1155 * t1530 + (t1341 + t1355 + 0.28e2 * t1370) * t1371 + t1261 * t1549 + (0.24e2 * t1142 * t1300 + t1285 * t1368 + t1346 * t1431 + t1357 * t1516 - 0.6e1 * t1375 * t1362 + 0.4e1 * t1378 ^ 2 + 0.28e2 * t1385 ^ 2) * t1374 + 0.8e1 * (-t1144 * t1552 - t1160 * t1490) * pkin(2) + (0.8e1 * t1150 * t1286 + 0.32e2 * t1191 * t1546) * pkin(1) + (0.16e2 * t1174 * t1298 + t1285 * t1370 + 0.70e2 * t1368 + t1382 + t1431) * t1382;
t1228 = t1579 + t1326 + t1464;
t1430 = t1328 + t1464;
t1229 = t1324 + t1430;
t1175 = -t1228 * t1572 + t1229 * t1287;
t1186 = -t1382 + (t1466 - t1483) * t1374 - 0.3e1 * t1368 + t1467 * t1370 + t1375;
t1232 = t1370 + t1430;
t1275 = t1590 - t1513;
t1189 = t1232 * t1287 + t1275 * t1486;
t1192 = -t1269 - 0.5e1 / 0.3e1 * t1368 + t1465 * t1370 + t1376 * t1463;
t1145 = t1175 * t1499 + t1186 * t1486 + (-0.8e1 / 0.3e1 * t1531 + t1192) * t1287 + (t1189 * t1286 - t1475) * t1570;
t1176 = t1382 + (t1518 + t1522) * t1374 + t1344 + 0.2e1 * t1551 + t1376 * (t1376 + t1522);
t1190 = t1349 + (0.5e1 * t1370 + t1273) * t1590 + (t1323 + t1460) * t1290;
t1154 = t1176 * t1287 - t1190 * t1572;
t1234 = t1350 + t1429;
t1235 = t1276 + t1601;
t1179 = -t1234 * t1572 + t1235 * t1287;
t1246 = t1326 + t1374 + t1523;
t1194 = pkin(3) * t1472 + t1246 * t1606;
t1146 = -0.4e1 * t1194 * t1530 + (t1179 * t1507 + t1494 * t1596) * pkin(1) + t1154;
t1187 = -(3 * t1382) + (t1467 - t1483) * t1374 + t1466 * t1370 + t1510;
t1156 = t1187 * t1287 + t1572 * t1608;
t1251 = 0.10e2 * t1551;
t1158 = t1371 + (0.21e2 * t1370 + t1273) * t1382 + (t1251 + t1353 + 0.35e2 * t1368 + 0.2e1 * t1548) * t1374 + (t1343 + (t1342 + t1354 - 0.5e1 * t1365) * t1370 + t1376 * (-t1365 + t1508)) * t1290;
t1231 = 0.3e1 / 0.2e1 * t1374 + t1347 + t1462;
t1250 = t1287 + 0.2e1 * t1572;
t1159 = t1292 * t1287 + 0.4e1 * t1555 + (t1231 * t1330 + t1250 * t1576) * t1569;
t1162 = 0.7e1 * t1371 + (t1345 + t1273) * t1349 + (t1251 + 0.21e2 * t1368 + 0.9e1 * t1375 + 0.6e1 * t1548) * t1374 + t1549;
t1233 = t1350 + 0.3e1 / 0.2e1 * t1370 + t1462;
t1193 = t1233 * t1287 + t1330 * t1489;
t1214 = 0.4e1 / 0.3e1 * t1530 + t1267 + t1280;
t1443 = -0.24e2 * t1214 * t1544;
t1133 = t1193 * t1451 + t1159 * t1448 + t1145 * t1487 + t1156 * t1500 + (-0.6e1 * t1154 * t1576 + (t1162 + t1443) * t1330) * pkin(3) + (0.6e1 * t1146 * t1552 + (t1531 * t1588 - t1158) * t1331) * pkin(2);
t1126 = t1127 * t1237 + t1133 * t1377;
t1403 = t1126 * t1470 + t1141 * t1578;
t1123 = t1403 * t1477;
t1404 = t1126 * t1469 + t1139 * t1577;
t1124 = t1404 * t1477;
t1225 = t1536 + t1541;
t1116 = t1123 * t1224 + t1124 * t1225;
t1117 = t1123 * t1225 - t1124 * t1224;
t1447 = t1116 * t1331 - t1333 * t1117;
t1440 = t1299 * t1382 * t1567;
t1438 = qJD(2) * t1475;
t1419 = t1272 * t1531;
t1413 = t1116 * t1333 + t1117 * t1331;
t1412 = t1176 * t1272 - t1190 * t1492;
t1408 = -0.48e2 * t1421;
t1245 = rSges(2,1) * t1333 - rSges(2,2) * t1331;
t1244 = rSges(4,1) * t1332 - rSges(4,2) * t1330;
t1243 = rSges(2,1) * t1331 + rSges(2,2) * t1333;
t1242 = rSges(4,1) * t1330 + rSges(4,2) * t1332;
t1178 = t1602 * t1225;
t1177 = t1602 * t1224;
t1138 = 0.1e1 / t1139 ^ 2;
t1135 = (t1139 * t1580 - t1330 * t1141 / 0.2e1) * t1559;
t1134 = (t1141 * t1580 + t1330 * t1582) * t1559;
t1129 = t1209 * t1468 + (-t1260 * t1332 + (t1198 * t1330 + t1535) * qJD(2)) * t1265 + (t1423 + 0.4e1 * t1436 - 0.2e1 * t1557) * pkin(3) + ((t1253 * t1332 - t1539) * t1563 + (-t1198 * t1537 + t1225 * t1377 + t1253 * t1540) * qJD(1) + 0.4e1 * t1397 + (t1255 * t1567 - t1402 * t1575) * t1571) * pkin(2);
t1128 = t1207 * t1468 + (t1168 * t1332 - t1265 * t1539) * qJD(2) + (t1212 * t1332 - t1220 * t1564) * t1330 * t1569 + ((-t1535 + (-t1216 - 0.8e1 * t1491) * t1330) * t1563 + ((t1172 - t1539) * t1333 + (t1535 + (t1265 * t1571 + t1216) * t1330 + (-pkin(3) + t1576 + 0.2e1 * t1573) * t1606) * t1331) * qJD(1)) * pkin(2);
t1121 = qJD(2) + ((t1128 * t1581 - t1141 * t1560) / t1582 - 0.2e1 * (t1129 * t1581 - t1139 * t1560) * t1141 * t1138) * pkin(4) / (t1138 * t1141 ^ 2 + 0.1e1) * t1183 * t1366;
t1120 = t1492 + t1121 * (rSges(5,1) * t1135 - rSges(5,2) * t1134);
t1119 = -t1493 - t1121 * (rSges(5,1) * t1134 + rSges(5,2) * t1135);
t1118 = (t1443 * t1492 - 0.24e2 * pkin(3) * (-0.8e1 / 0.3e1 * t1435 + t1259) * t1474 + 0.96e2 * t1214 * t1440 * t1572 + (t1292 * t1272 + 0.2e1 * t1231 * t1492 + (t1332 * t1607 + (-0.2e1 * t1250 * t1330 + 0.4e1 * t1573) * qJD(2)) * pkin(3) * pkin(1) - 0.8e1 * t1439 + 0.4e1 * t1556) * t1448 + (-0.8e1 / 0.3e1 * t1419 + (-t1228 * t1492 + t1229 * t1272) * t1499 + t1192 * t1272 - t1186 * t1492 / 0.2e1 + (0.32e2 / 0.3e1 * t1547 * t1287 + t1175 * t1332 * t1591) * t1564 + (0.4e1 * t1232 * t1481 * t1586 + ((0.12e2 * t1294 * t1297 - 0.4e1 * t1295) * t1368 + (t1189 * t1593 - 0.2e1 * t1275 * t1573) * pkin(3)) * qJD(2)) * pkin(1)) * t1487 + 0.24e2 * t1145 * t1434 + 0.6e1 * ((-0.4e1 * (pkin(3) * t1437 + t1246 * t1607) * t1297 + 0.8e1 * t1194 * t1479) * t1370 + (-0.8e1 * t1418 + (-t1234 * t1492 + t1235 * t1272) * t1507 + (-0.4e1 * pkin(3) * t1179 + 0.24e2 * t1297 * t1494) * t1564) * pkin(1) + t1412) * t1497 + t1419 * t1588 - 0.96e2 * t1281 * t1438 * t1287 + (t1233 * t1272 + t1489 * t1562) * t1451 + t1193 * t1408 + (t1187 * t1272 + t1492 * t1608) * t1500 + 0.12e2 * t1156 * t1435 - 0.6e1 * t1412 * t1505 + 0.6e1 * t1154 * t1454 - t1158 * t1272 + t1162 * t1492 + (-0.8e1 * t1420 + t1406) * t1159 - 0.6e1 * t1599 * t1146) * t1377 + t1133 * t1468 + 0.8e1 * (0.2e1 * (-0.48e2 * pkin(1) * t1545 - 0.2e1 * t1279 * t1332 - t1264 - 0.32e2 * t1547) * qJD(2) * t1474 - 0.8e1 * t1174 * t1440 - 0.4e1 * (t1205 * t1450 + (t1416 + (-t1211 * t1564 - t1558) * t1568) * pkin(1) + (-0.2e1 * t1401 * t1534 + t1545 * t1604) * pkin(2)) * t1476 + 0.3e1 * (-0.32e2 * t1282 * t1438 + (-0.2e1 / 0.3e1 * t1478 - 0.2e1 / 0.3e1 * t1473) * t1451 * t1586 + t1213 * t1408 - 0.64e2 / 0.3e1 * t1247 * t1395 - 0.28e2 * t1157 * t1435 + 0.6e1 * (-0.8e1 / 0.3e1 * t1478 - 0.8e1 / 0.3e1 * t1473) * t1248 * t1576 * t1574 + t1161 * t1442 + 0.4e1 * t1600 * t1258) * t1542 - 0.6e1 * t1142 * t1434 - (t1188 * t1417 + (-0.6e1 * pkin(1) * t1167 * t1564 + (-0.24e2 * pkin(1) * t1258 * t1398 + t1402 * t1609) * pkin(2)) * pkin(3) + t1527) * t1497 + t1252 * t1564 * t1596 - t1201 * t1455 * t1533 - 0.12e2 * t1191 * t1421 - 0.12e2 * t1258 * t1394 + t1155 * t1424 + t1396 * t1609 - t1150 * t1454 + t1610 * t1160 + (0.4e1 * t1420 + 0.12e2 * t1414) * t1153 + t1599 * t1144) * t1237 + t1127 * t1227;
t1115 = 0.1e1 / t1117 ^ 2;
t1112 = (t1404 * t1603 + (t1129 * t1577 + t1139 * t1561 / 0.8e1 + t1118 * t1469 + (t1128 * t1583 + t1141 * t1471) * t1126) * t1180) * t1532;
t1111 = (t1403 * t1603 + (t1118 * t1470 + t1128 * t1578 - t1141 * t1561 / 0.8e1 + (t1129 * t1583 + t1139 * t1471) * t1126) * t1180) * t1532;
t1110 = qJD(1) + ((t1111 * t1224 + t1112 * t1225 - t1123 * t1178 + t1124 * t1177) / t1117 - (t1111 * t1225 - t1112 * t1224 + t1123 * t1177 + t1124 * t1178) * t1116 * t1115) / (t1115 * t1116 ^ 2 + 0.1e1);
t1109 = t1272 + t1110 * (rSges(3,1) * t1447 + rSges(3,2) * t1413);
t1108 = -t1498 - t1110 * (-rSges(3,1) * t1413 + rSges(3,2) * t1447);
t1 = m(3) * (t1108 ^ 2 + t1109 ^ 2) / 0.2e1 + t1110 ^ 2 * Icges(3,3) / 0.2e1 + m(5) * (t1119 ^ 2 + t1120 ^ 2) / 0.2e1 + t1121 ^ 2 * Icges(5,3) / 0.2e1 + (m(2) * (t1243 ^ 2 + t1245 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (m(4) * (t1242 ^ 2 + t1244 ^ 2) / 0.2e1 + Icges(4,3) / 0.2e1) * qJD(2) ^ 2;
T = t1;