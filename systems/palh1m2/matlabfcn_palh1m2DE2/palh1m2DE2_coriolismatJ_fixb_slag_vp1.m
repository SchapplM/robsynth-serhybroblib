% Calculate matrix of centrifugal and coriolis load on the joints for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2DE2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE2_coriolismatJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:46
% EndTime: 2020-05-02 21:03:04
% DurationCPUTime: 123.82s
% Computational Cost: add. (746187->1015), mult. (1208471->1425), div. (24850->0), fcn. (1826307->74), ass. (0->709)
t1431 = sin(pkin(20));
t1493 = sin(pkin(18));
t1495 = cos(pkin(18));
t965 = cos(pkin(20));
t825 = -t1495 * t1431 + t1493 * t965;
t831 = t1493 * t1431 + t1495 * t965;
t968 = sin(qJ(3));
t972 = cos(qJ(3));
t1583 = t968 * t825 + t831 * t972;
t672 = t825 * t972 - t831 * t968;
t969 = sin(qJ(2));
t973 = cos(qJ(2));
t1571 = -t969 * t1583 + t672 * t973;
t1641 = t1583 * t973 + t969 * t672;
t958 = pkin(22) + pkin(21);
t924 = sin(t958);
t925 = cos(t958);
t961 = qJ(2) + qJ(3);
t359 = atan2(-t1571 * t925 + t1641 * t924, -t1571 * t924 - t1641 * t925) + t961;
t357 = sin(t359);
t967 = sin(qJ(4));
t971 = cos(qJ(4));
t271 = (-Icges(6,5) * t967 - Icges(6,6) * t971) * t357;
t358 = cos(t359);
t1387 = t358 * t271;
t1531 = -m(6) / 0.4e1;
t1416 = Icges(6,4) * t971;
t245 = -Icges(6,6) * t358 + (-Icges(6,2) * t967 + t1416) * t357;
t1417 = Icges(6,4) * t967;
t246 = -Icges(6,5) * t358 + (Icges(6,1) * t971 - t1417) * t357;
t272 = (-Icges(6,2) * t971 - t1417) * t357;
t273 = (-Icges(6,1) * t967 - t1416) * t357;
t1173 = pkin(18) - t958;
t906 = -pkin(20) + t1173;
t1078 = qJ(1) + t906;
t1046 = sin(t1078);
t1035 = t1046 / 0.2e1;
t1079 = -qJ(1) + t906;
t1049 = cos(t1079);
t1036 = t1049 / 0.2e1;
t1047 = sin(t1079);
t1048 = cos(t1078);
t1497 = pkin(11) + rSges(6,3);
t1076 = qJ(4) + t906;
t869 = qJ(1) + t1076;
t843 = sin(t869);
t1077 = -qJ(4) + t906;
t870 = -qJ(1) + t1077;
t844 = sin(t870);
t962 = qJ(1) + qJ(4);
t929 = sin(t962);
t1602 = t929 / 0.2e1 - t844 / 0.4e1 + t843 / 0.4e1;
t845 = cos(t869);
t846 = cos(t870);
t932 = cos(t962);
t1603 = t932 / 0.2e1 + t846 / 0.4e1 + t845 / 0.4e1;
t1056 = qJ(1) + t1077;
t1037 = sin(t1056);
t1057 = -qJ(1) + t1076;
t1038 = sin(t1057);
t1282 = qJ(1) - qJ(4);
t1192 = sin(t1282);
t984 = t1037 / 0.4e1 - t1192 / 0.2e1 - t1038 / 0.4e1;
t1039 = cos(t1056);
t1040 = cos(t1057);
t1193 = cos(t1282);
t985 = t1039 / 0.4e1 - t1193 / 0.2e1 + t1040 / 0.4e1;
t494 = (t984 - t1602) * rSges(6,2) + (t985 + t1603) * rSges(6,1);
t1283 = qJ(1) + qJ(2);
t1155 = sin(t1283) / 0.2e1;
t1209 = qJ(1) + t961;
t1086 = -cos(t1209) / 0.2e1;
t1210 = qJ(1) - t961;
t1156 = cos(t1210);
t1536 = -pkin(5) / 0.2e1;
t771 = pkin(5) * t1086 + t1156 * t1536;
t963 = qJ(1) - qJ(2);
t930 = sin(t963);
t700 = (t1155 - t930 / 0.2e1) * pkin(1) + t771;
t974 = cos(qJ(1));
t951 = t974 * pkin(15);
t421 = -t951 + (t1048 / 0.2e1 + t1036) * pkin(9) + t494 + t700 + (t1035 + t1047 / 0.2e1) * t1497;
t970 = sin(qJ(1));
t1437 = t970 * pkin(15);
t1073 = t1602 * rSges(6,1) + t1603 * rSges(6,2);
t1625 = rSges(6,1) * t984;
t988 = -t1040 - t1039;
t497 = t1073 + t1625 + (t1193 / 0.2e1 + t988 / 0.4e1) * rSges(6,2);
t1492 = pkin(1) * cos(t963);
t1484 = pkin(5) * sin(t1210);
t879 = sin(t1209) * t1536;
t769 = -t1484 / 0.2e1 + t879;
t897 = -pkin(1) * cos(t1283) / 0.2e1;
t698 = t1492 / 0.2e1 + t897 + t769;
t422 = (t1035 - t1047 / 0.2e1) * pkin(9) + t497 - t1437 + t698 + (-t1048 / 0.2e1 + t1036) * t1497;
t1582 = t846 + t845;
t495 = (t984 + t1602) * rSges(6,2) + (-t1582 / 0.4e1 - t932 / 0.2e1 + t985) * rSges(6,1);
t496 = rSges(6,2) * t985 + t1073 - t1625;
t1769 = 0.4e1 * (-t421 * t496 - t422 * t495) * t1531 + t1387 / 0.2e1 - ((t273 / 0.2e1 - t245 / 0.2e1) * t971 - (t246 / 0.2e1 + t272 / 0.2e1) * t967) * t357;
t1771 = t1769 * qJD(1);
t1225 = pkin(18) - pkin(22);
t927 = -qJ(2) + t1225;
t1224 = pkin(21) - atan2(cos(t927), -sin(t927));
t894 = t1173 - t961;
t662 = -atan2(-sin(t894), cos(t894)) + t1224;
t653 = sin(t662);
t654 = cos(t662);
t847 = -t969 * rSges(11,1) - rSges(11,2) * t973;
t848 = rSges(11,1) * t973 - rSges(11,2) * t969;
t1690 = -t653 * t848 - t654 * t847;
t1763 = t1690 * t970;
t655 = -qJ(2) + t662;
t650 = sin(t655);
t651 = cos(t655);
t540 = -rSges(11,1) * t650 + rSges(11,2) * t651;
t1569 = pkin(4) * sin(qJ(2) - t1224) - t540;
t1491 = pkin(1) * t969;
t913 = -pkin(15) + t1491;
t506 = rSges(11,3) * t970 + (-t913 - t1569) * t974;
t1770 = t1763 * t506;
t1141 = rSges(11,1) * t651 + rSges(11,2) * t650;
t1093 = -t653 * t847 + t654 * t848;
t523 = t1093 * t970;
t1768 = t1141 * t540 - t1763 * t523;
t1490 = pkin(1) * t973;
t519 = (t1093 - t1490) * t970;
t538 = t540 - t1491;
t1767 = t1141 * t538 - t1763 * t519;
t1762 = -t357 / 0.2e1;
t1586 = t357 / 0.2e1;
t1304 = t969 * t972;
t1306 = t968 * t973;
t833 = t1304 + t1306;
t1486 = pkin(5) * t833;
t1759 = Icges(11,3) + Icges(4,3) + Icges(9,3);
t1044 = t357 * t1486;
t1758 = Icges(7,3) + Icges(3,3);
t832 = -t968 * t1493 - t972 * t1495;
t834 = t972 * t1493 - t968 * t1495;
t1165 = -t832 * t973 + t834 * t969;
t1166 = t832 * t969 + t973 * t834;
t1433 = cos(pkin(22));
t964 = sin(pkin(22));
t829 = t1433 * t1495 + t1493 * t964;
t830 = -t1493 * t1433 + t1495 * t964;
t416 = -qJ(2) - atan2(t829 * t973 - t830 * t969, t829 * t969 + t830 * t973) + pkin(21) - atan2(t1165 * t924 - t1166 * t925, t1165 * t925 + t1166 * t924);
t414 = sin(t416);
t415 = cos(t416);
t1082 = Icges(11,5) * t414 - Icges(11,6) * t415;
t1757 = t1082 * t974;
t1671 = t843 - t844;
t1438 = t970 * rSges(6,1);
t1212 = t1438 / 0.2e1;
t956 = t974 * rSges(6,2);
t839 = t1212 + t956 / 0.2e1;
t840 = t1212 - t956 / 0.2e1;
t1435 = t974 * rSges(6,1);
t1211 = t1435 / 0.2e1;
t949 = t970 * rSges(6,2);
t841 = t1211 + t949 / 0.2e1;
t842 = t1211 - t949 / 0.2e1;
t858 = -t949 + t1435;
t859 = -t949 - t1435;
t862 = t956 + t1438;
t863 = t956 - t1438;
t987 = t1038 - t1037;
t411 = t862 * t1192 + t858 * t1193 + t1582 * t841 + t1671 * t840 + t987 * t839 + t988 * t842 - t859 * t932 - t863 * t929;
t1756 = t411 * t1762;
t1239 = Icges(11,4) * t415;
t1241 = Icges(11,1) * t414;
t1085 = -t1239 + t1241;
t1237 = Icges(11,5) * t970;
t325 = -t1085 * t974 + t1237;
t1395 = t325 * t414;
t1236 = Icges(11,2) * t415;
t1240 = Icges(11,4) * t414;
t1084 = -t1236 + t1240;
t1235 = Icges(11,6) * t970;
t323 = -t1084 * t974 + t1235;
t1397 = t323 * t415;
t328 = -Icges(11,1) * t415 - t1240;
t1658 = t328 * t974;
t327 = -Icges(11,2) * t414 - t1239;
t1659 = t327 * t974;
t1678 = t1658 * t415 + t1659 * t414 - t1395 + t1397;
t1677 = (-t1084 + t328) * t415 + (t1085 + t327) * t414;
t1704 = -t1082 * t970 + t1677 * t974;
t1432 = sin(pkin(19));
t966 = cos(pkin(19));
t823 = t972 * t1432 + t968 * t966;
t827 = -t968 * t1432 + t966 * t972;
t701 = qJ(2) + atan2(-t827, t823);
t691 = sin(t701);
t692 = cos(t701);
t1120 = Icges(9,5) * t691 + Icges(9,6) * t692;
t1617 = t1120 * t970;
t1413 = Icges(9,4) * t691;
t1129 = Icges(9,2) * t692 + t1413;
t1412 = Icges(9,4) * t692;
t1135 = Icges(9,1) * t691 + t1412;
t601 = -Icges(9,2) * t691 + t1412;
t603 = Icges(9,1) * t692 - t1413;
t1702 = (t1129 - t603) * t692 + (t601 + t1135) * t691;
t1732 = t1702 * t974 - t1617;
t1755 = t1732 + t1678 + t1704;
t1238 = Icges(11,4) * t970;
t324 = Icges(11,5) * t974 - t415 * t1238 + t970 * t1241;
t1396 = t324 * t414;
t322 = Icges(11,6) * t974 - t970 * t1236 + t414 * t1238;
t1106 = t322 * t415 - t1396;
t1656 = t328 * t970;
t1661 = t327 * t970;
t1679 = t1656 * t415 + t1661 * t414 - t1106;
t1705 = t1677 * t970 + t1757;
t1616 = t1120 * t974;
t1733 = t1702 * t970 + t1616;
t1754 = t1733 + t1679 + t1705;
t1753 = -t1678 / 0.4e1 - t1704 / 0.4e1 - t1732 / 0.4e1;
t1752 = t1678 / 0.2e1 + t1704 / 0.2e1 + t1732 / 0.2e1;
t1751 = -t1679 / 0.2e1 - t1705 / 0.2e1 - t1733 / 0.2e1;
t1750 = t1679 / 0.4e1 + t1705 / 0.4e1 + t1733 / 0.4e1;
t1749 = t415 / 0.2e1;
t1691 = Icges(11,5) * t415 + Icges(11,6) * t414;
t1657 = t970 * t1691;
t1660 = t974 * t1691;
t1105 = t322 * t414 + t324 * t415;
t1662 = t1105 * t974 - (-t323 * t414 - t325 * t415) * t970;
t1681 = t1658 * t414 - t1659 * t415;
t1682 = -t1656 * t414 + t1661 * t415;
t959 = t970 ^ 2;
t1700 = t959 * t1660 + (t1682 * t974 + (t1681 - t1657) * t970 + t1662) * t974;
t1748 = -t1700 / 0.2e1;
t960 = t974 ^ 2;
t1701 = t960 * t1657 + (t1681 * t970 + (t1682 - t1660) * t974 + t1662) * t970;
t1747 = t1701 / 0.2e1;
t1519 = -t414 / 0.2e1;
t931 = cos(t961);
t907 = Icges(4,4) * t931;
t928 = sin(t961);
t803 = -Icges(4,2) * t928 + t907;
t804 = Icges(4,1) * t928 + t907;
t1672 = t803 + t804;
t1421 = Icges(4,4) * t928;
t802 = Icges(4,2) * t931 + t1421;
t805 = Icges(4,1) * t931 - t1421;
t1636 = t327 * t1749 + t328 * t1519 + (-t802 / 0.2e1 + t805 / 0.2e1) * t928 + t1672 * t931 / 0.2e1;
t1729 = -t603 / 0.2e1;
t1730 = -t601 / 0.2e1;
t1744 = -t691 * t1729 - t692 * t1730 - t1636;
t1494 = sin(pkin(17));
t1496 = cos(pkin(17));
t1025 = t1493 * t1496 - t1495 * t1494;
t1026 = t1494 * t1493 + t1496 * t1495;
t1553 = t1025 * t973 - t1026 * t969;
t1554 = -t969 * t1025 - t1026 * t973;
t620 = atan2(-t827, -t823) + t701;
t618 = sin(t620);
t619 = cos(t620);
t1743 = -Icges(3,5) * t973 + Icges(7,5) * t1554 + Icges(10,5) * t619 + Icges(3,6) * t969 - Icges(7,6) * t1553 - Icges(10,6) * t618;
t1290 = t972 * t973;
t901 = pkin(5) * t1290;
t919 = pkin(5) * t968 + pkin(1);
t1581 = -t919 * t969 + t901;
t1742 = rSges(5,1) * cos(t906) - pkin(15) - t1581 - rSges(5,2) * sin(t906);
t1313 = t931 * t970;
t1317 = t928 * t970;
t1735 = -Icges(4,5) * t1313 + Icges(4,6) * t1317 - t415 * t1235 + t414 * t1237 + t1759 * t974 + t1617;
t801 = Icges(4,5) * t931 - Icges(4,6) * t928;
t1329 = t801 * t974;
t1734 = t1759 * t970 + t1329 - t1616 - t1757;
t1726 = t1129 / 0.2e1;
t1725 = -t1135 / 0.2e1;
t1276 = t246 + t272;
t1277 = -t245 + t273;
t1719 = (-t1387 + (-t1276 * t967 + t1277 * t971) * t357) * t358;
t1342 = t692 * t970;
t1344 = t691 * t970;
t587 = rSges(9,1) * t1344 + rSges(9,2) * t1342;
t543 = rSges(9,3) * t974 - t1437 + t587;
t1145 = rSges(9,1) * t691 + rSges(9,2) * t692;
t588 = t1145 * t974;
t1718 = t543 * t588;
t544 = -t970 * rSges(9,3) + t588 - t951;
t1146 = rSges(9,1) * t692 - rSges(9,2) * t691;
t586 = t1146 * t974;
t1717 = t544 * t586;
t1121 = Icges(9,5) * t692 - Icges(9,6) * t691;
t578 = t974 * t1121;
t1716 = t959 * t578;
t577 = t1121 * t970;
t1715 = t960 * t577;
t727 = Icges(4,5) * t970 + t805 * t974;
t1714 = -t727 * t1313 - t970 * t1397;
t1125 = Icges(3,5) * t969 + Icges(3,6) * t973;
t1409 = Icges(7,5) * t1553;
t1713 = (Icges(7,6) * t1554 - t1125 + t1409) * t974 + t1758 * t970;
t1348 = t1553 * t970;
t1350 = t1554 * t970;
t1673 = -Icges(7,5) * t1348 - Icges(7,6) * t1350 + t1125 * t970 + t1758 * t974;
t561 = Icges(9,6) * t970 - t1129 * t974;
t563 = Icges(9,5) * t970 - t1135 * t974;
t1094 = -t561 * t692 - t563 * t691;
t1710 = t1094 - t1395;
t1488 = pkin(2) * t692;
t1709 = rSges(10,2) * t618 + t1488;
t560 = Icges(9,6) * t974 + t1129 * t970;
t562 = Icges(9,5) * t974 + t1135 * t970;
t1095 = t560 * t691 - t562 * t692;
t1503 = t970 / 0.2e1;
t1504 = -t970 / 0.2e1;
t724 = Icges(4,4) * t1313 - Icges(4,2) * t1317 - Icges(4,6) * t974;
t873 = Icges(4,4) * t1317;
t726 = Icges(4,1) * t1313 - Icges(4,5) * t974 - t873;
t1708 = (t1503 + t1504) * (t724 * t931 + t726 * t928 + t1095 + t1105) + t1744;
t1706 = t1145 * t1146 - t586 * t588;
t632 = -rSges(5,3) * t970 + t1742 * t974;
t1703 = t692 * t1725 + t691 * t1726;
t1013 = pkin(2) * t691;
t1627 = t1553 / 0.2e1;
t1415 = Icges(7,4) * t1554;
t1597 = Icges(7,2) * t1553 - t1415;
t1650 = Icges(7,1) * t1553 + t1415;
t552 = Icges(7,5) * t970 + t1650 * t974;
t1699 = -t1597 * t974 + t552;
t1526 = m(10) / 0.4e1;
t1534 = m(4) / 0.4e1;
t1242 = t959 + t960;
t1162 = 0.4e1 * t1242;
t1307 = t968 * t969;
t835 = t1290 - t1307;
t1213 = pkin(5) ^ 2 * t835 * t833;
t770 = (t1156 / 0.2e1 + t1086) * pkin(5);
t1332 = t769 * t770;
t768 = t1484 / 0.2e1 + t879;
t1333 = t768 * t771;
t1529 = m(6) / 0.4e1;
t1532 = m(5) / 0.4e1;
t706 = -0.4e1 * t1213;
t1605 = (-0.4e1 * t1332 + t706 + 0.4e1 * t1333) * t1529 + (t1162 * t1213 + t706) * t1532;
t1220 = -0.1e1 + t1242;
t1151 = 0.4e1 * t1220;
t433 = pkin(2) ^ 2 * t692 * t691 * t1151;
t1081 = rSges(4,1) * t1313 - rSges(4,2) * t1317 - t974 * rSges(4,3);
t1312 = t931 * t974;
t1243 = -rSges(4,1) * t1312 - t970 * rSges(4,3);
t1316 = t928 * t974;
t622 = t970 * t1081 + t974 * (-rSges(4,2) * t1316 - t1243);
t822 = rSges(4,1) * t928 + rSges(4,2) * t931;
t761 = t822 * t970;
t762 = t822 * t974;
t644 = -t970 * t761 - t974 * t762;
t1447 = rSges(4,2) * t928;
t826 = rSges(4,1) * t931 - t1447;
t1567 = t1242 * t826 * t822 + t622 * t644;
t468 = 0.4e1 * t1567;
t1689 = -t433 * t1526 + t468 * t1534 + t1605;
t1381 = t415 * t974;
t1687 = -t726 * t1312 + t322 * t1381 + t1735 * t970;
t1644 = t727 * t1312 + t323 * t1381 + t1734 * t970;
t1686 = t1734 * t974 + t1714;
t1683 = t974 * rSges(5,3) + t1742 * t970;
t1680 = -t1084 * t1519 + t1085 * t1749;
t1676 = -0.2e1 * t496;
t1675 = -0.2e1 * t768;
t1453 = m(5) * qJD(1);
t1347 = t1553 * t974;
t553 = Icges(7,1) * t1348 + Icges(7,4) * t1350 - Icges(7,5) * t974;
t1669 = -t553 * t1347 + t1673 * t970;
t1423 = Icges(3,4) * t969;
t1133 = Icges(3,2) * t973 + t1423;
t749 = Icges(3,6) * t974 + t1133 * t970;
t1422 = Icges(3,4) * t973;
t1140 = Icges(3,1) * t969 + t1422;
t751 = Icges(3,5) * t974 + t1140 * t970;
t1088 = t749 * t973 + t751 * t969;
t551 = Icges(7,4) * t1348 + Icges(7,2) * t1350 - Icges(7,6) * t974;
t1370 = t551 * t1554;
t1655 = t1088 + t1370;
t516 = t552 * t1348;
t1653 = t1713 * t974 - t516;
t1096 = t560 * t692 + t562 * t691;
t1652 = t724 * t928 - t1096;
t1651 = -t1683 * t970 - t632 * t974;
t1414 = Icges(7,4) * t1553;
t1649 = -Icges(7,1) * t1554 + t1414;
t725 = Icges(4,6) * t970 + t803 * t974;
t1648 = -t725 * t1317 + t1710 * t970 - t1686;
t1383 = t414 * t974;
t1647 = -t1096 * t974 + t724 * t1316 - t324 * t1383 + t1687;
t1646 = t1094 * t974 - t725 * t1316 - t325 * t1383 + t1644;
t1359 = t619 * t974;
t1362 = t618 * t974;
t508 = Icges(10,5) * t618 + Icges(10,6) * t619;
t489 = Icges(10,3) * t970 + t508 * t974;
t612 = Icges(10,4) * t618;
t1572 = Icges(10,2) * t619 + t612;
t491 = Icges(10,6) * t970 + t1572 * t974;
t607 = Icges(10,4) * t1359;
t493 = Icges(10,1) * t1362 + Icges(10,5) * t970 + t607;
t1645 = t552 * t1347 + t491 * t1359 + t493 * t1362 + (t489 + t1713) * t970;
t1643 = t725 * t928 - t1710 + t1735;
t1512 = t1554 / 0.2e1;
t1515 = -t1554 / 0.2e1;
t1528 = m(7) / 0.4e1;
t857 = rSges(3,1) * t969 + rSges(3,2) * t973;
t1570 = t974 * rSges(3,3) + (-pkin(15) + t857) * t970;
t1045 = t822 + t1490;
t1579 = t1045 * t974;
t1580 = t1045 * t970;
t1587 = 0.4e1 * t1534;
t1588 = 0.4e1 * t1529;
t1590 = 0.4e1 * t1526;
t1524 = m(11) / 0.4e1;
t1591 = 0.4e1 * t1524;
t865 = -rSges(7,1) * t1493 + t1495 * rSges(7,2);
t866 = t1495 * rSges(7,1) + rSges(7,2) * t1493;
t710 = -t1494 * t866 - t865 * t1496;
t711 = -t1494 * t865 + t866 * t1496;
t1091 = t710 * t969 + t711 * t973;
t1592 = -0.4e1 * t1091;
t1598 = Icges(7,2) * t1554 + t1414;
t382 = t1141 - t1490;
t345 = t382 * t970;
t346 = t382 * t974;
t1143 = rSges(10,1) * t618 + rSges(10,2) * t619;
t666 = pkin(2) * t1344;
t453 = rSges(10,3) * t974 + t666 + (-pkin(15) - t1143) * t970;
t1487 = pkin(2) * t974;
t458 = -t1143 * t974 + t691 * t1487;
t454 = rSges(10,3) * t970 - t458 + t951;
t1360 = t619 * t970;
t455 = -rSges(10,1) * t1360 + t1709 * t970;
t472 = rSges(10,1) * t619 - t1709;
t456 = t472 * t974;
t868 = t913 * t970;
t505 = rSges(11,3) * t974 + t1569 * t970 + t868;
t1092 = t710 * t973 - t711 * t969;
t952 = t974 * rSges(7,3);
t583 = t952 + (pkin(14) - t1092) * t970;
t1259 = t970 * rSges(7,3) + t1092 * t974;
t584 = pkin(14) * t974 - t1259;
t585 = rSges(9,1) * t1342 - rSges(9,2) * t1344;
t1568 = -sin(t1225) * rSges(8,2) + rSges(8,1) * cos(t1225);
t677 = -t970 * rSges(8,3) + (t913 + t1568) * t974;
t678 = t974 * rSges(8,3) + t1568 * t970 + t868;
t688 = (t913 + t1447) * t974 + t1243;
t689 = -t1081 + t868;
t697 = t897 - t1492 / 0.2e1 + t768;
t699 = t770 + (t1155 + t930 / 0.2e1) * pkin(1);
t721 = rSges(3,3) * t970 - t857 * t974 + t951;
t856 = -t973 * rSges(3,1) + rSges(3,2) * t969;
t791 = t856 * t974;
t851 = -Icges(3,2) * t969 + t1422;
t1637 = -(t851 / 0.2e1 + t1140 / 0.2e1) * t973 + (-t345 * t505 + t346 * t506) * t1591 + (t1579 * t688 + t1580 * t689) * t1587 + (-t421 * t697 + t422 * t699) * t1588 + (-t583 * t970 - t584 * t974) * t1528 * t1592 + m(8) * (t677 * t974 + t678 * t970) * t1490 + (t543 * t585 + t1717) * m(9) + (t453 * t455 + t454 * t456) * t1590 + t1650 * t1515 + t1597 * t1512 + (t1649 + t1598) * t1627 + m(3) * (-t1570 * t856 * t970 + t721 * t791) + t1680;
t1391 = t357 * t970;
t1284 = t974 * t967;
t1295 = t970 * t971;
t353 = t358 * t1295 - t1284;
t350 = Icges(6,4) * t353;
t1291 = t971 * t974;
t1308 = t967 * t970;
t352 = t358 * t1308 + t1291;
t232 = Icges(6,2) * t352 - Icges(6,6) * t1391 - t350;
t349 = Icges(6,4) * t352;
t234 = Icges(6,1) * t353 + Icges(6,5) * t1391 - t349;
t354 = -t358 * t1284 + t1295;
t355 = t358 * t1291 + t1308;
t1634 = -t232 * t354 + t355 * t234;
t1633 = t699 * t1676 - 0.2e1 * t495 * t697;
t1632 = t495 * t1675 + t770 * t1676;
t1631 = t699 * t1675 + 0.2e1 * t697 * t770;
t1630 = -m(4) / 0.4e1;
t1500 = -t974 / 0.2e1;
t1498 = t974 / 0.2e1;
t580 = t601 * t974;
t582 = t603 * t974;
t1101 = t580 * t692 + t582 * t691;
t579 = t601 * t970;
t581 = t603 * t970;
t1102 = -t579 * t692 - t581 * t691;
t1549 = t1095 * t974 + (t561 * t691 - t563 * t692) * t970;
t148 = -t1715 + (t1101 * t970 + (t1102 + t578) * t974 + t1549) * t970;
t1624 = t148 + t1701;
t149 = -t1716 + (t1102 * t974 + (t1101 + t577) * t970 + t1549) * t974;
t1623 = t149 + t1700;
t1261 = -t563 + t580;
t1262 = t562 + t579;
t1263 = t561 + t582;
t1264 = t560 - t581;
t1543 = (t1261 * t970 - t1262 * t974) * t692 + (t1263 * t970 + t1264 * t974) * t691;
t226 = -t1715 + (t974 * t578 + t1543) * t970;
t1622 = t226 + t1701;
t227 = -t1716 + (t970 * t577 + t1543) * t974;
t1621 = t227 + t1700;
t1390 = t357 * t974;
t228 = Icges(6,5) * t353 - Icges(6,6) * t352 + Icges(6,3) * t1391;
t1612 = t228 * t1390;
t750 = Icges(3,6) * t970 - t1133 * t974;
t752 = Icges(3,5) * t970 - t1140 * t974;
t1087 = -t750 * t973 - t752 * t969;
t550 = Icges(7,6) * t970 + t1598 * t974;
t1609 = -t1554 * t550 - t1087 + t1673;
t1608 = t1683 * t974 - t632 * t970;
t1349 = t1554 * t974;
t1607 = t1087 * t974 + t550 * t1349 + t1645;
t1379 = (-Icges(10,3) * t974 + t508 * t970) * t974;
t1606 = t1379 + t1645;
t1596 = t1087 * t970 + t1088 * t974 + t551 * t1349 + t550 * t1350 - t1653 - t1669;
t1594 = 0.2e1 * m(9);
t1593 = 0.4e1 * m(9);
t490 = -Icges(10,6) * t974 + t1572 * t970;
t1411 = Icges(10,4) * t619;
t511 = Icges(10,1) * t618 + t1411;
t492 = -Icges(10,5) * t974 + t511 * t970;
t1578 = (t490 * t619 + t492 * t618) * t970;
t121 = t1612 + t1634;
t1573 = t121 - t1612;
t512 = -Icges(10,1) * t619 + t612;
t1560 = t1743 * t970;
t1559 = t1743 * t974;
t1557 = t1633 * t1529;
t1556 = t1632 * t1529;
t1164 = t692 * t455 * t1487;
t432 = -0.2e1 * t1164;
t1555 = (t432 + 0.2e1 * t1164) * t1526;
t795 = -pkin(5) * t1307 + t901;
t1223 = 0.2e1 * t795;
t1326 = t826 * t974;
t1221 = -0.2e1 * t1326;
t1327 = t826 * t970;
t1222 = 0.2e1 * t1327;
t1257 = t689 * t1221 + t688 * t1222;
t1538 = 0.2e1 * t771;
t1272 = t422 * t1538 - 0.2e1 * t769 * t421;
t400 = t1145 * t970;
t402 = t1146 * t970;
t1462 = (-t400 * t544 + t1718 + (-t402 + t585) * t586) * t1594;
t265 = t1013 * t970;
t266 = t1013 * t974;
t267 = t1488 * t970;
t268 = t1488 * t974;
t347 = -0.2e1 * t1770;
t420 = t974 * t1690;
t348 = -0.2e1 * t420 * t505;
t428 = t970 * t1141;
t429 = t974 * t1141;
t1489 = pkin(1) * t974;
t1323 = t848 * t974;
t1324 = t847 * t974;
t524 = t654 * t1323 - t653 * t1324;
t520 = -t973 * t1489 + t524;
t1546 = -m(11) * (-0.2e1 * t428 * t520 + 0.2e1 * t429 * t519 + t347 + t348) / 0.4e1 - m(10) * (t265 * t454 + t266 * t453 + t267 * t456 + t268 * t455) / 0.2e1 - t1462 / 0.4e1 + (t1631 + t1272) * t1531 + t1608 * t1223 * t1532 + (-0.2e1 * t1579 * t761 + 0.2e1 * t1580 * t762 + t1257) * t1630;
t1485 = pkin(5) * t835;
t1260 = -0.2e1 * t1608 * t1485;
t1460 = (-t544 * t587 + t1718) * t1594;
t342 = t653 * t1323 + t654 * t1324;
t1545 = -m(11) * (t342 * t505 - t345 * t524 + t346 * t523 - t1770) / 0.2e1 - m(10) * (t432 + 0.2e1 * (-t456 * t1342 + (t453 * t974 + t454 * t970) * t691) * pkin(2)) / 0.4e1 - t1460 / 0.4e1 + (-t1631 + t1272) * t1531 - m(5) * t1260 / 0.4e1 + (0.2e1 * (t1579 * t970 - t1580 * t974) * t822 + t1257) * t1630;
t1248 = -t802 * t974 + t727;
t1249 = -Icges(4,2) * t1313 + t726 - t873;
t1250 = -t804 * t974 - t725;
t1251 = t804 * t970 + t724;
t1544 = (-t1248 * t970 + t1249 * t974) * t928 + (t1250 * t970 + t1251 * t974) * t931;
t501 = -Icges(10,2) * t1362 + t607;
t503 = t512 * t974;
t570 = t1649 * t974;
t787 = t851 * t974;
t853 = Icges(3,1) * t973 - t1423;
t789 = t853 * t974;
t1542 = t1699 * t1350 + (-t550 - t570) * t1348 + (t493 + t501) * t1360 + ((-t752 + t787) * t973 + (t750 + t789) * t969 + (-t491 - t503) * t618) * t970;
t1128 = -Icges(10,2) * t618 + t1411;
t500 = t1128 * t970;
t502 = t512 * t970;
t569 = t1597 * t970;
t571 = t1649 * t970;
t786 = t851 * t970;
t788 = t853 * t970;
t1541 = -(t751 + t786) * t973 + (t749 - t788) * t969 - (t553 - t569) * t1554 + (t551 + t571) * t1553 - (t492 + t500) * t619 + (t490 + t502) * t618;
t24 = (t1644 * t970 + ((-t1396 + t1652 + t1734) * t974 + t1648 + t1687 + t1714) * t974) * t1500 + (t1646 * t970 + t1647 * t974) * t1498 + ((t1643 * t970 - t1647 + t1648 + t1686) * t970 + ((t1643 - t1735) * t974 + (-t726 * t931 + t1106 + t1652) * t970 - t1644 + t1646) * t974) * t1503;
t1539 = 0.2e1 * t495;
t1537 = 0.4e1 * t1581;
t1535 = m(4) / 0.2e1;
t1533 = m(5) / 0.2e1;
t1530 = m(6) / 0.2e1;
t1525 = m(8) * pkin(1);
t1502 = t970 / 0.4e1;
t1501 = t973 / 0.4e1;
t1499 = -t974 / 0.4e1;
t621 = t622 - t1491;
t1358 = t621 * t644;
t1481 = 0.4e1 * m(4) * (t1358 + (t1579 * t974 + t1580 * t970) * t826);
t1479 = m(4) * t822;
t1475 = m(5) * t1581;
t1468 = m(6) * t1633;
t1466 = m(6) * t1632;
t1463 = m(7) * t1092;
t1461 = (-t585 * t400 + t1706) * t1593;
t1459 = (-t585 * t587 + t1706) * t1593;
t1452 = m(6) * qJD(2);
t1451 = m(6) * qJD(4);
t1450 = m(10) * qJD(2);
t1430 = m(11) * qJD(2);
t1429 = m(11) * qJD(3);
t1418 = Icges(6,4) * t355;
t1369 = t553 * t1553;
t1340 = t697 * t771;
t1339 = t699 * t769;
t899 = pkin(5) * t1304;
t794 = -pkin(5) * t1306 - t899;
t1070 = t357 * t794;
t385 = -t858 * t1192 + t862 * t1193 - t929 * t859 + t863 * t932 - t987 * t842 + t1671 * t841 - t1582 * t840 + t988 * t839 + ((t1049 - t1048) * t974 + (-t1047 - t1046) * t970) * rSges(6,3);
t82 = (t1044 * t1529 - t1070 * t1531) * t385;
t1328 = t82 * qJD(4);
t1281 = -Icges(6,1) * t352 + t232 - t350;
t233 = Icges(6,2) * t354 + Icges(6,6) * t1390 + t1418;
t1280 = Icges(6,1) * t354 - t1418 - t233;
t1279 = -Icges(6,2) * t353 + t234 - t349;
t351 = Icges(6,4) * t354;
t236 = Icges(6,1) * t355 + Icges(6,5) * t1390 + t351;
t1278 = -Icges(6,2) * t355 + t236 + t351;
t1124 = Icges(4,5) * t928 + Icges(4,6) * t931;
t755 = t1124 * t970;
t756 = t974 * t1124;
t1273 = (-t959 * t756 + (t970 * t755 + t1544) * t974) * t1503 + (-t960 * t755 + (t974 * t756 + t1544) * t970) * t1500;
t1271 = t496 * t1538 + t769 * t1539;
t1233 = qJD(2) * t974;
t1232 = qJD(3) * t970;
t1231 = qJD(3) * t974;
t1228 = -0.2e1 * t1486;
t1226 = 0.2e1 * t1485;
t122 = (Icges(6,5) * t355 + Icges(6,6) * t354 + Icges(6,3) * t1390) * t1390 + t354 * t233 + t355 * t236;
t1216 = ((t228 * t1391 + t122) * t974 + t1573 * t970) * t1586 + (t121 * t970 + t122 * t974) * t1762;
t1215 = (t1573 - t1634) * t974 * t1586;
t773 = t919 * t973 + t899;
t1191 = t773 * t1226;
t1186 = t833 * t1242;
t1161 = pkin(2) * m(10);
t1152 = -t493 * t618 * t970 - t491 * t1360 + t489 * t974;
t1150 = -t1125 / 0.2e1 + Icges(7,6) * t1512 + t1409 / 0.2e1 + t508 / 0.2e1;
t1142 = t773 * t1162;
t1100 = -t455 * t970 + t456 * t974;
t1059 = t1273 + t1461 / 0.4e1 + t1623 * t1503 + t1624 * t1500;
t1058 = t1273 + t1459 / 0.4e1 + t1621 * t1503 + t1622 * t1500;
t1033 = (-t802 + t805) * t931 - t1672 * t928;
t308 = -t1379 + t1578;
t1030 = t1152 * t1504 + (t308 + t1673 * t974 + (t1369 + t1655) * t970) * t1500 + (t1609 * t974 - t1606 + t1607) * t1498 + (t490 * t1359 + t492 * t1362 + t1609 * t970 + t1152 + t1596 + t1653) * t1503;
t1029 = (t308 - t1578 + t1606) * t1504 + t1607 * t1503 + (-t516 + (-t1655 + t1713) * t974 + t1596 + t1669) * t1500;
t335 = t1262 * t691 + t1264 * t692;
t336 = t1261 * t691 - t1263 * t692;
t1000 = -t1545 + (t336 + t1755) * t1502 + (t335 + t1754) * t1499;
t991 = t24 + t1555;
t990 = -t24 + (t1033 * t974 + t1248 * t931 + t1250 * t928 + t801 * t970) * t1503 + (t1033 * t970 + t1249 * t931 - t1251 * t928 - t1329) * t1500;
t240 = t579 * t691 - t581 * t692 + t1096;
t241 = t580 * t691 - t582 * t692 + t1094;
t986 = -t1546 + (t241 + t1755) * t1502 + (t240 + t1754) * t1499;
t981 = (-t428 * t505 + t429 * t506) * t1591 + (t688 * t762 + t689 * t761) * t1587 + (-t421 * t768 + t422 * t770) * t1588 + (t402 * t543 + t1717) * m(9) + (t267 * t453 - t268 * t454) * t1590 + t1680 + t1703;
t920 = 0.4e1 * t919;
t903 = t969 * t1489;
t902 = t970 * t1491;
t767 = -t920 * t1501 - t899;
t766 = -t920 * t969 / 0.4e1 + t901;
t717 = t903 - t1326;
t716 = t902 - t1327;
t704 = t794 * t1226;
t693 = t1581 * t1228;
t690 = t766 * t1228;
t670 = t960 * t1191;
t669 = t959 * t1191;
t648 = t1579 * t1221;
t647 = t1580 * t1222;
t640 = -0.2e1 * t1339;
t638 = 0.2e1 * t1340;
t636 = t644 - t1490;
t529 = t795 * t1142 + t794 * t1537;
t526 = 0.2e1 * t1358;
t504 = 0.4e1 * t766 * t794 - 0.4e1 * t1339 + 0.4e1 * t1340;
t485 = 0.4e1 * t1651 * t794;
t484 = 0.4e1 * t1651 * t773;
t473 = t1143 - t1013;
t457 = -t1143 * t970 + t666;
t376 = t1466 / 0.4e1;
t361 = t1468 / 0.4e1;
t339 = t342 + t903;
t338 = -t1763 + t902;
t260 = Icges(6,5) * t354 - Icges(6,6) * t355;
t259 = -Icges(6,5) * t352 - Icges(6,6) * t353;
t248 = -0.4e1 * t420 * t520 + 0.4e1 * t1767;
t199 = 0.4e1 * t342 * t524 + 0.4e1 * t1768;
t154 = -0.4e1 * t1488 * t473 - 0.4e1 * t265 * t455 + 0.4e1 * t266 * t456;
t97 = t1276 * t354 + t1277 * t355 + t271 * t1390;
t96 = -t1276 * t352 + t1277 * t353 + t271 * t1391;
t95 = -t260 * t358 + (-t1278 * t967 + t1280 * t971) * t357;
t94 = -t259 * t358 + (-t1279 * t967 + t1281 * t971) * t357;
t83 = (0.2e1 * t1271 + (t1044 - t1070) * t385) * t1529;
t23 = t199 * t1524 + t1058 + t1689;
t22 = (t1133 / 0.2e1 - t853 / 0.2e1) * t969 + (t1726 + t1729) * t691 + (-t1128 / 0.2e1 - t511 / 0.2e1) * t619 + (t1730 + t1725) * t692 + (t1572 / 0.2e1 + t512 / 0.2e1) * t618 - t484 * t1532 + t1636 + t1637;
t21 = t485 * t1532 - t1744 + t981;
t20 = t1029 * t974 + t1030 * t970 + t24;
t19 = t504 * t1529 + t1481 / 0.4e1 + t529 * t1532 + t248 * t1524 + t154 * t1526 + t1059;
t12 = (t1215 * t974 + t1216 * t970) * t357;
t11 = (-t241 / 0.4e1 + t1753) * t970 + (t240 / 0.4e1 + t1750) * t974 + t1000 + t991 + t1546;
t10 = (t335 / 0.4e1 + t1750) * t974 + (-t336 / 0.4e1 + t1753) * t970 + t991 + t986 + t1545;
t9 = t1000 + t990 + t986 - t1555;
t8 = t376 - t1556;
t7 = t376 + t1556;
t6 = -t1466 / 0.4e1 - t1556;
t5 = -t1557 + t361;
t4 = t361 + t1557;
t3 = -t1468 / 0.4e1 - t1557;
t1 = [t22 * qJD(2) + t21 * qJD(3) - qJD(4) * t1769, t22 * qJD(1) + t9 * qJD(3) + t5 * qJD(4) + (t453 * t458 + t454 * t457) * t1450 + (-t421 * t698 + t422 * t700) * t1452 + (t338 * t506 + t339 * t505 - t345 * t520 + t346 * t519) * t1430 + (-t1370 / 0.2e1 - t1369 / 0.2e1 + t569 * t1627 - t571 * t1512 - t335 / 0.2e1 + (-t492 / 0.2e1 - t500 / 0.2e1) * t618 + (t678 * t1525 - t751 / 0.2e1 - t786 / 0.2e1) * t969 + (-t490 / 0.2e1 - t502 / 0.2e1) * t619 + (-t749 / 0.2e1 + t788 / 0.2e1) * t973 + t1150 * t974 - t1029 - t583 * t1463 - t1683 * t1475 + (-t856 ^ 2 * t970 + t1570 * t857) * m(3) + t1751) * t1233 + (t1460 / 0.2e1 + t990 + (-t688 * t716 + t689 * t717) * m(4) + (t550 * t1512 - t570 * t1515 + (-t750 / 0.2e1 - t789 / 0.2e1) * t973 + (t503 / 0.2e1 + t491 / 0.2e1) * t619 + (-t752 / 0.2e1 + t787 / 0.2e1 - t677 * t1525) * t969 + t1150 * t970 + (t721 * t857 + t791 * t856) * m(3) + t584 * t1463 + t632 * t1475 - t1030 + (t501 / 0.2e1 + t493 / 0.2e1) * t618 + t336 / 0.2e1 + t1699 * t1627 + t1752) * t970) * qJD(2), t21 * qJD(1) + t9 * qJD(2) + (t1257 * t1535 + t1462 / 0.2e1 + t1260 * t1533 + t1272 * t1530 + t990) * qJD(3) + t8 * qJD(4) + (-t428 * t524 + t429 * t523 + t347 / 0.2e1 + t348 / 0.2e1) * t1429 + (-t761 * t1479 - t240 / 0.2e1 + (-t267 * t692 + t453 * t691) * t1161 + t1751) * t1231 + (t762 * t1479 + t241 / 0.2e1 + (t268 * t692 + t454 * t691) * t1161 + t1752) * t1232, -t1771 + t5 * qJD(2) + t8 * qJD(3) + (-t421 * t497 + t422 * t494) * t1451 + (-t1719 + ((t95 / 0.2e1 + t97 / 0.2e1 - t1215) * t974 + (t96 / 0.2e1 + t94 / 0.2e1 - t1216) * t970) * t357) * qJD(4); t20 * qJD(2) + t10 * qJD(3) + t4 * qJD(4) + t484 * t1453 / 0.4e1 + ((-t1133 + t853) * t969 / 0.2e1 - t1637 - t1703 + t1708 - (t1572 + t512) * t618 / 0.2e1 + (t511 + t1128) * t619 / 0.2e1) * qJD(1), t20 * qJD(1) + t19 * qJD(3) + ((t1091 * t1162 * t1092 + t1242 * (t970 * (t1092 * t970 - t952) + t974 * t1259) * t1592) * t1528 + (t1142 * t1581 - t1537 * t773) * t1532 + t1220 * t857 * t856 * m(3) + t1058 + (t338 * t519 + t339 * t520 + t382 * t538) * m(11) - m(8) * pkin(1) ^ 2 * t969 * t1151 * t1501 + ((t1541 * t974 - t1560 * t970 + t1542) * t974 + t1559 * t959) * t1503 + (((t1541 - t1559) * t974 + t1542) * t970 + t1560 * t960) * t1500 + m(10) * (-t455 * t457 + t456 * t458 + t472 * t473) + m(4) * (-t1579 * t717 - t1580 * t716 + t621 * t636) + m(6) * (t697 * t700 - t698 * t699 + t766 * t767)) * qJD(2), t10 * qJD(1) + t19 * qJD(2) - t1328 + (t226 / 0.2e1 + t1747 - t1624) * t1231 + (t1748 - t227 / 0.2e1 + t1623) * t1232 + (-t199 / 0.4e1 + (t538 + t540) * t1141 - (t520 + t524) * t420 - (t519 + t523) * t1763) * t1429 + ((pkin(5) * t1186 * t1223 + t669 + t670 + t693 + t704) * t1533 + (-0.2e1 * t1332 + t638 + t640 + t690 + t704 + 0.2e1 * t1333) * t1530 + t1461 / 0.2e1 - t1459 / 0.4e1 + t1273 + (t526 / 0.2e1 + t647 / 0.2e1 - t648 / 0.2e1 - t468 / 0.4e1 + t1567) * m(4) + (t433 / 0.4e1 + ((-t265 * t970 - t266 * t974 - t473) * t692 + (t1488 + t1100) * t691) * pkin(2)) * m(10) - t1605) * qJD(3), t4 * qJD(1) - t82 * qJD(3) + (t766 * t1756 + t494 * t697 - t497 * t699) * t1451; (-t981 + t1708) * qJD(1) + t11 * qJD(2) + t24 * qJD(3) + t7 * qJD(4) - t485 * t1453 / 0.4e1, t11 * qJD(1) + t23 * qJD(3) + t1328 + (-t504 / 0.4e1 + t767 * t1485 - t698 * t770 + t700 * t768 + t638 / 0.2e1 + t640 / 0.2e1 + t690 / 0.2e1) * t1452 + (t148 / 0.2e1 + t1747 - t717 * t1479 - t1622) * t1233 + (-t154 / 0.4e1 + ((-t457 * t970 - t458 * t974 - t473) * t692 + (t1100 - t472) * t691) * pkin(2)) * t1450 + (-t248 / 0.4e1 + t338 * t523 + t339 * t524 + t342 * t520 + t382 * t540 + t1767) * t1430 + (-t1481 / 0.4e1 - t1461 / 0.4e1 + (0.2e1 * t622 * t636 + t526 + t647 - t648) * t1535 + t1459 / 0.2e1 + t1273 + (-t529 / 0.4e1 + t669 / 0.2e1 + t670 / 0.2e1 + t693 / 0.2e1 + (t1186 * t1581 - t773 * t835) * pkin(5)) * m(5) + (t1748 - t149 / 0.2e1 - t716 * t1479 + t1621) * t970) * qJD(2), t24 * qJD(1) + t23 * qJD(2) + ((-t420 * t524 + t1768) * m(11) + t1059 + t1689) * qJD(3), t7 * qJD(1) + t82 * qJD(2) + (t1485 * t1756 + t494 * t768 - t497 * t770) * t1451; t3 * qJD(2) + t6 * qJD(3) + t12 * qJD(4) + t1771, t3 * qJD(1) + (-t357 * t767 * t385 + t698 * t1539 + 0.2e1 * t496 * t700) * t1452 / 0.2e1 + t83 * qJD(3), t6 * qJD(1) + t83 * qJD(2) + (t385 * t1044 + t1271) * qJD(3) * t1530, t12 * qJD(1) + (-(-t1719 + (t94 * t970 + t95 * t974) * t357) * t358 / 0.2e1 + ((-t1278 * t352 + t1280 * t353 + t1391 * t260) * t1390 + (-t1279 * t352 + t1281 * t353 + t1391 * t259) * t1391 - t96 * t358) * t1391 / 0.2e1 + (t357 ^ 2 * t385 * t411 + 0.4e1 * t494 * t496 + 0.4e1 * t495 * t497) * t1529 + ((t1278 * t354 + t1280 * t355 + t1390 * t260) * t1390 + (t1279 * t354 + t1281 * t355 + t1390 * t259) * t1391 - t97 * t358) * t1390 / 0.2e1) * qJD(4);];
Cq = t1;