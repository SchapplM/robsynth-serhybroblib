% Calculate joint inertia matrix with Newton Euler for
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_inertiaJ_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_inertiaJ_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_inertiaJ_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2IC_inertiaJ_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2IC_inertiaJ_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:52:56
% EndTime: 2020-05-03 06:53:02
% DurationCPUTime: 6.35s
% Computational Cost: add. (2533->339), mult. (3772->443), div. (0->0), fcn. (2248->10), ass. (0->181)
t1624 = pkin(3) * m(7);
t1497 = sin(qJ(6));
t1502 = cos(qJ(6));
t1615 = t1502 * mrSges(7,1);
t1640 = -mrSges(7,2) * t1497 + t1615;
t1428 = mrSges(6,1) + t1640 + t1624;
t1489 = t1502 ^ 2;
t1507 = m(6) + m(7);
t1487 = m(5) + t1507;
t1648 = m(4) + t1487;
t1665 = t1648 * pkin(4) ^ 2;
t1500 = sin(qJ(3));
t1501 = sin(qJ(2));
t1590 = t1500 * t1501;
t1444 = pkin(1) * t1590 - pkin(2);
t1664 = 0.2e1 * t1444;
t1498 = sin(qJ(5));
t1503 = cos(qJ(5));
t1616 = Ifges(7,4) * t1497;
t1439 = (pkin(3) * mrSges(7,1) + t1616) * t1502;
t1495 = Ifges(7,1) - Ifges(7,2);
t1494 = pkin(3) ^ 2 * m(7);
t1571 = Ifges(6,2) + Ifges(7,3) + t1494;
t1538 = -t1495 * t1489 - Ifges(6,1) - Ifges(7,2) + 0.2e1 * t1439 + t1571;
t1609 = pkin(3) * t1497;
t1563 = mrSges(7,2) * t1609;
t1639 = -t1538 + 0.2e1 * t1563;
t1540 = 0.2e1 * t1639;
t1496 = mrSges(7,3) + mrSges(6,2);
t1622 = (pkin(5) * t1496);
t1663 = (t1540 * t1498 - (2 * t1622)) * t1503;
t1447 = mrSges(7,1) * t1497 + mrSges(7,2) * t1502;
t1662 = -mrSges(6,3) + t1447;
t1596 = t1496 * t1498;
t1437 = -pkin(5) * t1507 - mrSges(5,1) + t1596;
t1623 = pkin(4) * t1500;
t1661 = t1437 * t1623;
t1505 = cos(qJ(3));
t1463 = pkin(4) * t1505 + pkin(2);
t1504 = cos(qJ(4));
t1499 = sin(qJ(4));
t1593 = t1499 * t1500;
t1660 = -pkin(4) * t1593 + (-pkin(2) + t1463) * t1504;
t1414 = t1428 * t1503;
t1659 = t1437 - t1414;
t1565 = pkin(5) * t1596;
t1652 = -Ifges(5,1) + Ifges(5,2) - 0.2e1 * t1565;
t1658 = t1540 + 0.2e1 * t1652;
t1628 = -0.2e1 * t1503;
t1490 = t1503 ^ 2;
t1656 = t1639 * t1490;
t1441 = Ifges(7,5) * t1497 + Ifges(7,6) * t1502 - Ifges(6,6);
t1493 = pkin(3) * mrSges(7,2);
t1636 = -mrSges(7,1) * t1609 + 0.2e1 * Ifges(7,4) * t1489 + (t1495 * t1497 - t1493) * t1502 - Ifges(7,4) - Ifges(6,5);
t1362 = t1498 * t1441 - t1503 * t1636;
t1529 = pkin(5) ^ 2;
t1421 = pkin(5) * t1428;
t1476 = Ifges(7,6) * t1497;
t1477 = Ifges(7,5) * t1502;
t1446 = t1477 - t1476;
t1617 = mrSges(7,3) * pkin(3);
t1427 = -Ifges(6,4) - t1446 + t1617;
t1601 = t1427 * t1498;
t1605 = (t1421 - 0.2e1 * t1601) * t1503;
t1654 = t1605 - t1565 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + (m(6) / 0.2e1 + m(7) / 0.2e1) * t1529 + t1563;
t1637 = t1498 * t1639 - t1622;
t1592 = t1499 * t1503;
t1435 = t1498 * t1504 + t1592;
t1594 = t1498 * t1499;
t1546 = -t1503 * t1504 + t1594;
t1651 = -t1435 * t1505 + t1546 * t1500;
t1626 = Ifges(7,2) / 0.2e1;
t1650 = t1439 + (t1626 - Ifges(7,1) / 0.2e1) * t1489 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1 + t1494 / 0.2e1;
t1649 = 0.2e1 * t1499;
t1471 = t1487 * pkin(2) ^ 2;
t1647 = -Ifges(4,2) - t1471;
t1635 = -Ifges(7,2) / 0.2e1 + t1650;
t1373 = -t1563 + t1635;
t1606 = t1373 * t1490;
t1645 = t1499 * (0.2e1 * t1606 - t1635 + t1654);
t1644 = (-t1659 - t1414 + t1596) * t1504;
t1413 = t1498 * t1428;
t1408 = t1413 + mrSges(5,2);
t1466 = t1496 * t1503;
t1552 = -t1408 - t1466;
t1433 = pkin(3) * t1640 + Ifges(7,3);
t1383 = t1433 * t1498 - t1446 * t1503;
t1634 = -0.2e1 * pkin(4);
t1633 = -0.2e1 * t1408;
t1455 = pkin(2) * t1487 + mrSges(4,1);
t1632 = 0.2e1 * t1455;
t1631 = -0.2e1 * t1496;
t1629 = -0.4e1 * t1501;
t1627 = 0.2e1 * t1503;
t1625 = -(2 * Ifges(4,1)) + 0.2e1 * Ifges(4,2);
t1621 = pkin(5) * t1498;
t1614 = mrSges(4,2) * t1634;
t1613 = pkin(1) * t1501;
t1403 = pkin(2) * t1408;
t1612 = pkin(2) * t1496;
t1611 = pkin(2) * t1499;
t1610 = pkin(2) * t1500;
t1459 = -0.2e1 * t1563;
t1379 = t1459 + t1538;
t1491 = t1504 ^ 2;
t1567 = pkin(5) * t1413;
t1382 = -Ifges(5,4) - t1427 + t1567;
t1602 = t1427 * t1490;
t1535 = t1382 + 0.2e1 * t1602;
t1608 = ((t1379 * t1498 + t1622) * t1503 + t1535) * t1491;
t1360 = t1373 * t1498 + t1622 / 0.2e1;
t1607 = t1360 * t1500;
t1377 = Ifges(4,4) + t1382;
t1604 = t1377 * t1500;
t1603 = t1408 * t1499;
t1597 = t1455 * t1500;
t1595 = t1496 * t1499;
t1591 = t1500 * pkin(1);
t1589 = t1500 * t1504;
t1480 = t1507 * t1529;
t1461 = pkin(5) * t1503 + pkin(3);
t1588 = pkin(3) - t1461;
t1356 = pkin(5) * t1662 + Ifges(5,5) + t1362;
t1363 = -t1441 * t1503 - t1498 * t1636;
t1361 = Ifges(5,6) + t1363;
t1348 = t1356 * t1504 - t1361 * t1499;
t1349 = t1356 * t1499 + t1361 * t1504;
t1506 = cos(qJ(2));
t1587 = (t1348 * t1505 - t1349 * t1500) * t1501 + (t1348 * t1500 + t1349 * t1505) * t1506;
t1420 = pkin(2) * t1428;
t1586 = -0.4e1 * t1499 * t1637 - 0.2e1 * t1420;
t1583 = 0.8e1 * t1601 - 0.4e1 * t1421;
t1375 = -0.4e1 * t1601 + 0.2e1 * t1421;
t1412 = pkin(4) * t1589 + t1463 * t1499;
t1582 = -pkin(2) * t1592 + t1412 * t1503 + t1498 * t1660;
t1581 = t1466 + t1413;
t1579 = mrSges(7,3) * t1621;
t1578 = -0.2e1 * t1613;
t1577 = -0.4e1 * t1610;
t1576 = 0.4e1 * t1610;
t1575 = -0.4e1 * t1602;
t1442 = -t1610 + t1613;
t1574 = -0.2e1 * t1442 * t1499;
t1570 = t1428 * t1623;
t1568 = t1496 * t1623;
t1566 = t1447 * t1621;
t1564 = t1408 * t1623;
t1562 = pkin(2) * t1595;
t1561 = t1500 * t1608;
t1560 = pkin(2) * t1603;
t1558 = t1500 * t1602;
t1424 = pkin(2) * t1437;
t1557 = 0.2e1 * t1424 + (0.4e1 * t1382 + 0.8e1 * t1602) * t1499;
t1549 = -0.8e1 * t1382 - 0.16e2 * t1602;
t1366 = pkin(5) * t1640 + t1433 * t1503 + t1446 * t1498;
t1545 = 0.2e1 * t1552;
t1544 = -0.2e1 * t1659;
t1541 = 0.4e1 * t1639;
t1539 = 0.4e1 * t1645;
t1537 = t1480 + t1639 + t1652;
t1533 = t1541 * t1490 - 0.2e1 * t1480 - t1658;
t1532 = t1375 * t1503 + t1537;
t1368 = t1540 * t1490;
t1531 = t1537 + (t1583 * t1503 + t1533) * t1491 - t1368 + Ifges(4,1) + t1647;
t1492 = t1505 ^ 2;
t1399 = pkin(4) * t1631 + t1428 * t1577;
t1398 = t1420 - t1568;
t1397 = -t1570 + t1612;
t1396 = t1570 + t1612;
t1395 = -0.2e1 * t1564;
t1393 = mrSges(5,2) + t1581;
t1381 = -t1435 * t1500 - t1505 * t1546;
t1372 = pkin(4) * t1633 - t1437 * t1577;
t1371 = -t1403 + t1661;
t1370 = t1403 + t1661;
t1354 = t1381 * t1506 + t1501 * t1651;
t1353 = t1381 * t1501 - t1506 * t1651;
t1352 = t1362 * t1499 + t1363 * t1504;
t1351 = t1362 * t1504 - t1363 * t1499;
t1339 = Ifges(3,6) * t1506 + (Ifges(3,5) + (-mrSges(4,3) - mrSges(5,3) + t1662) * pkin(4)) * t1501;
t1 = [(t1659 * t1664 + (t1575 - 0.2e1 * t1382 + t1663) * t1499) * t1504 + (0.4e1 * t1561 + (t1442 * t1545 + t1500 * t1539) * t1504 - 0.4e1 * t1558 + (t1428 * t1574 - 0.4e1 * t1607) * t1503 - t1437 * t1574 + mrSges(4,2) * t1578 - 0.2e1 * t1604) * t1505 + ((0.2e1 * t1562 + t1375) * t1503 + t1531 + 0.2e1 * t1560 + (t1586 * t1503 + t1557) * t1504) * t1492 + t1459 + (t1444 * t1595 + t1601) * t1627 + t1603 * t1664 + t1656 + t1571 + (mrSges(3,2) + t1597) * t1578 + 0.2e1 * pkin(3) * t1615 + (-t1368 + t1532) * t1491 + Ifges(3,1) + Ifges(5,1) + Ifges(2,3) - t1647 + (m(3) + t1648) * pkin(1) ^ 2 + (((0.4e1 * t1602 + 0.2e1 * t1567 - 0.2e1 * t1617 + 0.2e1 * t1477 - 0.2e1 * t1476 - 0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,4) - t1663) * t1491 + (pkin(2) * t1466 + t1403 + (t1532 + 0.4e1 * t1606) * t1499) * t1504 + (t1428 * t1611 + t1637) * t1503 - t1437 * t1611 - Ifges(4,4) - t1535) * t1492 * t1629 + ((pkin(1) * t1544 + (t1399 * t1503 + t1372 + (0.16e2 * t1360 * t1503 - t1549) * t1593) * t1501) * t1504 + ((-t1625 + t1658) * t1500 + t1614 + 0.2e1 * (t1480 - t1471) * t1500) * t1501 + pkin(1) * t1632 + (((t1428 * t1634 + t1496 * t1576) * t1501 + pkin(1) * t1631) * t1503 + (t1408 * t1576 - t1437 * t1634) * t1501 + pkin(1) * t1633) * t1499 + (-0.8e1 * (t1379 * t1490 + t1626 - t1650 + t1654) * t1491 + 0.8e1 * t1606 + 0.4e1 * t1605) * t1590) * t1505 + 0.4e1 * t1501 * t1608 + (t1545 * t1591 + (t1397 * t1627 + 0.2e1 * t1370 + t1539) * t1501) * t1504 + t1501 * t1575 + (((t1420 + t1568) * t1501 - t1428 * t1591) * t1649 + t1360 * t1629) * t1503 + ((-t1424 + t1564) * t1501 + t1437 * t1591) * t1649 + 0.2e1 * (-pkin(4) * t1597 + Ifges(3,4) - t1377) * t1501 + 0.2e1 * (t1648 * pkin(4) - mrSges(4,2) * t1500 + mrSges(3,1)) * pkin(1) + ((-0.8e1 * t1561 + (pkin(4) * t1544 + (0.4e1 * pkin(2) * t1552 - 0.8e1 * t1645) * t1500) * t1504 + 0.8e1 * t1558 + (t1399 * t1499 + 0.8e1 * t1607) * t1503 + t1372 * t1499 + 0.4e1 * t1604 + pkin(4) * t1632) * t1505 + ((-0.8e1 * t1656 + (0.8e1 * t1421 - 0.16e2 * t1601) * t1503 - 0.8e1 * t1565 + 0.4e1 * t1480 - 0.4e1 * Ifges(5,1) + 0.4e1 * Ifges(5,2) + t1541) * t1491 + (0.4e1 * t1420 * t1503 - 0.4e1 * t1424 + (0.8e1 * t1503 * t1637 + t1549) * t1499) * t1504 + t1533 + 0.2e1 * t1471 - 0.4e1 * t1560 + (-0.4e1 * t1562 + t1583) * t1503 + t1625) * t1492 + (t1397 * t1649 + t1375) * t1503 + t1531 + t1500 * t1614 + ((-0.2e1 * t1568 + t1586) * t1503 + t1395 + t1557) * t1504 + t1370 * t1649 - Ifges(3,1) + Ifges(3,2) + t1665) * t1506) * t1506, t1339, -(t1351 * t1500 + t1352 * t1505) * t1506 - t1501 * (t1351 * t1505 - t1352 * t1500) + t1587, ((pkin(2) * t1640 + t1366 * t1504) * t1505 - t1383 * t1589 + pkin(4) * t1640 + (-t1366 * t1500 - t1383 * t1505) * t1499) * t1506 + (-pkin(2) * t1590 + pkin(1)) * t1640 + t1383 * t1499 * t1590 + (-(t1366 * t1499 + t1383 * t1504) * t1505 - t1366 * t1589) * t1501; t1339, Ifges(3,3) + t1665 + t1562 * t1628 + (-0.2e1 * t1661 + 0.2e1 * t1403 + 0.2e1 * t1371 + (-t1570 - 0.2e1 * t1612 + t1396) * t1628) * t1499 + (0.2e1 * t1564 + t1398 * t1627 + (0.2e1 * t1420 - t1568) * t1628 + t1395 + 0.2e1 * (-t1659 + t1437) * pkin(2)) * t1504, (t1398 * t1498 + t1371) * t1499 + (t1396 * t1498 - t1424) * t1504 + (-t1408 * t1589 + (t1644 + (-t1393 + t1581) * t1499) * t1505) * pkin(4) + (-t1644 + (t1408 - t1413) * t1499) * pkin(2), t1651 * t1447 * pkin(4); t1441 * t1354 + (-Ifges(6,5) + t1497 * (Ifges(7,1) * t1502 - t1616) + t1502 * (Ifges(7,4) * t1502 - Ifges(7,2) * t1497) - pkin(3) * t1447) * t1353 + t1587, t1661 * t1499 + (-(t1393 * t1499 + t1504 * t1659) * t1505 + (-t1428 * t1592 + t1504 * t1552) * t1500) * pkin(4) + t1428 * (-pkin(2) * t1594 + t1412 * t1498 - t1503 * t1660) + (mrSges(6,2) + (t1497 ^ 2 + t1489) * mrSges(7,3)) * t1582, Ifges(5,3) + t1480 + t1588 * t1624 + ((-mrSges(6,1) + t1428) * t1503 + (mrSges(6,2) - t1496) * t1498) * pkin(5) + (mrSges(7,2) * t1461 + t1497 * t1579 - t1493) * t1497 + (t1588 * mrSges(7,1) + t1502 * t1579) * t1502, -t1566; Ifges(7,3) * t1354 + t1446 * t1353 - t1640 * (-pkin(1) - t1506 * pkin(4) - (t1505 * t1506 - t1590) * pkin(2) - ((t1504 * t1505 - t1593) * t1506 - t1501 * (t1499 * t1505 + t1589)) * pkin(5) - t1354 * pkin(3)), -t1447 * t1582, -t1566, Ifges(7,3);];
Mq = t1;
