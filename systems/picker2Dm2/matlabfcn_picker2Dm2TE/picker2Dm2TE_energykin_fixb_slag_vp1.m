% Calculate kinetic energy for
% picker2Dm2TE
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm2TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2TE_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:35:40
% EndTime: 2020-05-09 10:36:07
% DurationCPUTime: 24.99s
% Computational Cost: add. (253847->703), mult. (789006->1166), div. (3626->23), fcn. (134952->10), ass. (0->495)
t1967 = 4 * pkin(1);
t1736 = (pkin(3) ^ 2);
t2027 = -2 * t1736;
t1743 = (pkin(7) ^ 2);
t1726 = -2 * t1743;
t1693 = cos(pkin(8));
t1696 = sin(qJ(1));
t1698 = cos(qJ(1));
t1974 = sin(pkin(8));
t1577 = t1693 * t1696 - t1974 * t1698;
t1566 = t1577 * qJD(1);
t1579 = t1693 * t1698 + t1696 * t1974;
t1567 = t1579 * qJD(1);
t1995 = 0.1e1 / t1579 ^ 2;
t1491 = (t1566 * t1577 * t1995 + t1567 / t1579) / (t1577 ^ 2 * t1995 + 0.1e1);
t1695 = sin(qJ(2));
t1994 = 0.2e1 * t1695;
t1916 = pkin(7) * t1994;
t1632 = pkin(3) * t1916;
t1895 = -t1736 + t1743;
t1660 = t1695 ^ 2;
t1936 = t1660 * t1736;
t1581 = t1632 + t1895 + 0.2e1 * t1936;
t1741 = pkin(1) ^ 2;
t1731 = (pkin(4) ^ 2);
t1891 = t1743 - t1731;
t1638 = t1741 + t1891;
t1590 = t1632 + t1638;
t1712 = 2 * t1736;
t1589 = t1712 + t1590;
t1649 = pkin(3) * t1695;
t1630 = t1649 + pkin(7);
t1651 = pkin(1) * t1698;
t1889 = 0.2e1 * t1651;
t1516 = t1581 * t1889 + t1589 * t1630;
t1664 = t1698 ^ 2;
t1522 = t1589 * t1698 + (0.4e1 * t1664 - 0.2e1) * t1630 * pkin(1);
t1697 = cos(qJ(2));
t1928 = t1696 * t1697;
t1875 = pkin(3) * t1928;
t1611 = -pkin(1) + t1875;
t1942 = t1630 * t1698;
t1560 = -t1611 + t1942;
t1650 = pkin(3) * t1697;
t1887 = 0.2e1 * t1650;
t1866 = pkin(7) * t1887;
t1623 = qJD(2) * t1866;
t1959 = qJD(2) * t1736;
t1811 = t1695 * t1697 * t1959;
t1798 = 0.4e1 * t1811;
t1568 = t1623 + t1798;
t1827 = pkin(1) * t1875;
t1612 = -0.2e1 * t1827;
t1556 = t1612 + t1590;
t2015 = 4 * t1741;
t1665 = t1736 * t2015;
t1920 = t1736 * t1743;
t1636 = t1665 - 4 * t1920;
t1739 = t1741 ^ 2;
t1998 = -0.4e1 * pkin(7);
t1867 = pkin(3) * t1638 * t1998;
t1944 = t1581 * t1664;
t2018 = 0.2e1 * pkin(3);
t1495 = t1636 * t1660 + t1695 * t1867 - t1739 - (t1743 - (t2018 + pkin(4)) * pkin(4)) * (t1743 + (t2018 - pkin(4)) * pkin(4)) + (t1726 + (2 * t1731) - (4 * t1736) - 0.4e1 * t1944) * t1741 + (-t1556 * t1942 + t1590 * t1875) * t1967;
t1744 = sqrt(t1495);
t1962 = qJD(1) * t1698;
t1862 = t1696 * t1962;
t1812 = t1581 * t1862;
t1930 = t1695 * t1696;
t1989 = pkin(1) * pkin(3);
t1828 = t1930 * t1989;
t1800 = qJD(2) * t1828;
t1926 = t1697 * t1698;
t1861 = qJD(1) * t1926;
t1801 = t1861 * t1989;
t2024 = -t1801 + t1800;
t1908 = 0.2e1 * t2024;
t1844 = t1623 + t1908;
t1945 = t1568 * t1664;
t1960 = qJD(2) * t1697;
t1964 = qJD(1) * t1696;
t1969 = 0.2e1 * pkin(7);
t1951 = ((0.8e1 * t1812 - 0.4e1 * t1945) * t1741 + (t1636 * t1994 + t1867) * t1960 + (t1696 * t1697 ^ 2 * t1959 * t1969 + (t1556 * t1964 - t1698 * t1844) * t1630 + (t1590 * t1861 + (-t1556 * t1926 - t1590 * t1930) * qJD(2)) * pkin(3)) * t1967) / t1744;
t1847 = t1951 / 0.2e1;
t1878 = pkin(1) * t1942;
t1925 = t1698 * t1744;
t1927 = t1696 * t1744;
t1963 = qJD(1) * t1697;
t1968 = 2 * pkin(1);
t1973 = pkin(1) * t1664;
t1987 = pkin(7) * t1698;
t1464 = t1560 * t1847 + (t1516 * t1698 - t1630 * t1927) * qJD(1) + (t1568 * t1698 - t1581 * t1964) * t1696 * t1968 + ((-t1925 + (-t1589 - 0.8e1 * t1878) * t1696) * t1963 + ((-t1522 + t1927) * t1695 + (t1925 + (t1630 * t1969 + t1589) * t1696 + (-pkin(1) + t1987 + 0.2e1 * t1973) * t1697 * t2018) * t1697) * qJD(2)) * pkin(3);
t1718 = 3 * t1741;
t1625 = t1712 + t1718 + t1891;
t1816 = -0.4e1 * t1827;
t1552 = t1625 + t1632 + t1816;
t1874 = pkin(3) * t1926;
t1943 = t1630 * t1696;
t1562 = t1874 + t1943;
t1859 = qJD(2) * t1930;
t2003 = t1859 - t1861;
t1765 = -0.4e1 * t2003;
t1763 = t1765 * t1651;
t1929 = t1695 * t1698;
t1786 = -t1928 + t1929;
t1465 = t1562 * t1847 + (-t1623 * t1698 + (t1552 * t1696 + t1925) * qJD(1)) * t1630 + (t1798 + 0.4e1 * t1812 - 0.2e1 * t1945) * pkin(1) + ((t1625 * t1698 - t1927) * t1963 + (-t1552 * t1926 - t1625 * t1930 - t1744 * t1786) * qJD(2) + (t1611 * t1960 - t1649 * t2003) * t1969 + t1630 * t1763) * pkin(3);
t1484 = t1516 * t1696 + t1522 * t1650 + t1560 * t1744;
t1655 = t1741 + t1743;
t1835 = t1736 + t1655;
t1533 = t1612 + t1632 + t1835 + 0.2e1 * t1878;
t1531 = 0.1e1 / t1533;
t1737 = 0.1e1 / pkin(3);
t1480 = -t1552 * t1942 + t1562 * t1744 + (t1611 * t1916 + t1625 * t1928) * pkin(3) + (-0.2e1 * t1944 + (0.2e1 * t1660 - 0.4e1) * t1736 - t1638) * pkin(1);
t1979 = t1579 / 0.2e1;
t1980 = t1577 / 0.2e1;
t1777 = t1480 * t1980 + t1484 * t1979;
t1950 = ((-qJD(1) * t1943 + qJD(2) * t1874) * t1968 + t1844) / t1533 ^ 2;
t1982 = -t1566 / 0.2e1;
t2013 = t1480 / 0.2e1;
t1445 = (-t1777 * t1950 + (t1464 * t1979 + t1465 * t1980 + t1484 * t1982 + t1567 * t2013) * t1531) * t1737;
t1981 = -t1577 / 0.2e1;
t1776 = t1480 * t1979 + t1484 * t1981;
t1984 = t1484 / 0.2e1;
t1446 = (-t1776 * t1950 + (t1464 * t1981 + t1465 * t1979 + t1480 * t1982 - t1567 * t1984) * t1531) * t1737;
t1949 = t1531 * t1737;
t1472 = t1777 * t1949;
t1699 = cos(pkin(9));
t1766 = t1776 * t1949;
t1988 = sin(pkin(9));
t1461 = t1472 * t1988 + t1699 * t1766;
t1459 = t1472 * t1699 - t1766 * t1988;
t1996 = 0.1e1 / t1461 ^ 2;
t2002 = 0.1e1 / (t1459 ^ 2 * t1996 + 0.1e1);
t2026 = t2002 * (-t1699 * t1445 + t1446 * t1988) / t1461;
t2025 = -0.6e1 * t1801 + 0.6e1 * t1800;
t1608 = t1743 + t1736 / 0.4e1 + t1741 / 0.4e1 - t1731 / 0.8e1;
t1730 = t1731 ^ 2;
t1742 = t1743 ^ 2;
t1749 = t1736 ^ 2;
t1904 = 0.4e1 / 0.7e1 * t1743 - t1731 / 0.7e1;
t1918 = t1743 * t1731;
t1508 = -0.32e2 / 0.21e2 * t1608 * t1827 + 0.5e1 / 0.42e2 * t1749 + (0.16e2 / 0.21e2 * t1741 + t1904) * t1736 + t1739 / 0.7e1 + t1904 * t1741 + t1742 - 0.3e1 / 0.7e1 * t1918 + t1730 / 0.42e2;
t1906 = t1741 / 0.3e1 + t1743;
t1676 = -t1731 / 0.4e1;
t2006 = t1676 + t1736 / 0.2e1;
t1610 = t1906 + t2006;
t1677 = -t1731 / 0.3e1;
t1679 = -0.2e1 / 0.3e1 * t1731;
t1683 = 0.4e1 / 0.3e1 * t1736;
t1689 = 0.4e1 / 0.3e1 * t1741;
t1510 = -0.8e1 / 0.3e1 * t1610 * t1827 + 0.5e1 / 0.18e2 * t1749 + (t1689 + t1677) * t1736 + t1742 - t1739 / 0.3e1 + t1730 / 0.18e2 + (t1683 + 0.2e1 / 0.3e1 * t1741 + t1679) * t1743;
t1894 = t1739 + t1742;
t1725 = 2 * t1743;
t1898 = t1725 - t1731;
t1774 = (t1898 * t1741) + t1730 / 0.6e1 + t1894 - t1918;
t1565 = -t1749 / 0.6e1 + t1774;
t1691 = t1741 / 0.2e1;
t1905 = t1691 + t1743;
t1569 = -0.2e1 / 0.3e1 * t1827 + t1676 + t1905;
t1723 = 4 * t1743;
t1642 = (t1723 + t1731) * t1741;
t1688 = -0.2e1 / 0.3e1 * t1736;
t1646 = t1688 + t1743;
t1647 = -t1741 / 0.3e1 + t1743;
t1656 = -t1741 + t1743;
t1678 = -t1731 / 0.2e1;
t1619 = t1678 + t1835;
t1833 = -0.4e1 * t1875;
t1807 = t1619 * t1833;
t1663 = t1698 * t1664;
t1745 = pkin(1) * t1741;
t1933 = t1663 * t1745;
t1886 = pkin(7) * t1933;
t1831 = 0.16e2 * t1886;
t1662 = t1664 ^ 2;
t1935 = t1662 * t1739;
t1883 = 0.8e1 * t1935;
t1914 = 0.6e1 * t1987;
t1919 = t1741 * t1664;
t1923 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t1486 = t1646 * t1883 + t1569 * t1831 + 0.14e2 * t1508 * t1919 + (t1656 * t1749) + (t1642 - 0.10e2 / 0.3e1 * t1739 + (2 * t1742) - t1918) * t1736 + t1565 * t1923 + (t1510 * t1914 + t1647 * t1807) * pkin(1);
t1666 = 0.10e2 / 0.3e1 * t1741;
t1769 = 0.5e1 / 0.6e1 * t1749 + t1774;
t1544 = (t1666 + t1898) * t1736 + t1769;
t1680 = -0.3e1 / 0.2e1 * t1731;
t1700 = 15 * t1739;
t1707 = 18 * t1743;
t1720 = 3 * t1742;
t1907 = t1730 / 0.2e1 - t1749 / 0.2e1;
t1808 = -(3 * t1918) + t1720 + t1907;
t1724 = 3 * t1743;
t1903 = 15 * t1741 + t1724;
t1748 = pkin(3) * t1736;
t1733 = t1748 ^ 2;
t1911 = t1655 * ((t1680 + t1725) * t1741 - 0.3e1 / 0.2e1 * t1918 + t1894 + t1907) + t1733;
t1496 = -0.6e1 * t1544 * t1827 + (t1700 + ((t1707 - 9 * t1731) * t1741) + t1808) * t1736 + (t1680 + t1903) * t1749 + t1911;
t1708 = -2 * t1731;
t1721 = 8 * t1743;
t1991 = 4 * t1739;
t1557 = t1745 * t1833 + t1665 + t1991 + ((t1708 + t1721) * t1741);
t1790 = t1743 - t1827;
t1563 = -t1741 + t1790 + t2006;
t1657 = -3 * t1741 + t1743;
t1832 = -0.2e1 * t1875;
t1845 = 0.8e1 * t1886;
t1915 = 0.4e1 * t1987;
t1498 = t1845 + t1557 * t1664 + t1619 * t1657 + (t1563 * t1915 + t1832 * t1923) * pkin(1);
t1687 = -t1736 / 0.3e1;
t1645 = t1687 + t1743;
t1791 = pkin(1) * t1807;
t1893 = t1742 - t1739;
t1502 = t1645 * t1791 - t1733 + (-t1666 - t1891) * t1749 + (t1642 + t1749 / 0.6e1 - t1730 / 0.6e1 + t1893) * t1736 + t1565 * t1743;
t1709 = -5 * t1731;
t1715 = 7 * t1739;
t1506 = (t1680 + t1724 + (7 * t1741)) * t1749 + (t1715 + ((t1709 + 10 * t1743) * t1741) + t1808) * t1736 + t1911;
t1753 = pkin(7) * t1743;
t1637 = -0.12e2 * pkin(7) * t1745 + t1753 * t1967;
t1917 = t1743 * t1741;
t1644 = -8 * t1739 + 12 * t1917;
t1524 = t1637 * t1698 + t1644 * t1664 + t1831 + t1883 + t1894 - (6 * t1917);
t1654 = -3 * t1736 + t1743;
t1924 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t1539 = t1612 * t1924 + t1619 * t1654;
t1892 = t1742 + t1749;
t1603 = 16 * (t1892 - 6 * t1920) * t1739;
t1626 = -t1731 + t1835;
t1631 = t1651 + pkin(7);
t1652 = -30 * t1731 + 60 * t1743;
t1658 = t1660 ^ 2;
t1710 = -6 * t1731;
t1717 = 6 * t1741;
t1897 = t1730 - t1749;
t1780 = 6 * t1742 + t1897 - 6 * t1918;
t1830 = 0.32e2 * t1886;
t1659 = t1695 * t1660;
t1937 = t1659 * t1748;
t1648 = t1655 ^ 2;
t1939 = t1648 * (-t1736 + t1638);
t1515 = t1791 + ((t1717 + t1898) * t1736) + t1769;
t1583 = t1645 * t1612;
t1537 = t1619 * t1924 + t1583;
t1602 = t1654 * t1845;
t1890 = pkin(7) * t1651;
t1865 = 0.6e1 * t1890;
t1868 = 0.12e2 * t1919;
t1487 = t1515 * t1865 + t1537 * t1868 + t1496 + t1602;
t2017 = 0.8e1 * t1487;
t1463 = t1603 * t1662 + t1539 * t1830 + 0.24e2 * t1502 * t1919 + (t1708 + t1723 + 28 * t1741) * t1733 + (t1626 * t1939) + (0.24e2 * t1486 * t1660 + (t1652 * t1739) + (t1710 * t1742) + (t1780 * t1717) + (t1897 * t1725) + (28 * t1745 ^ 2) + 0.4e1 * t1753 ^ 2) * t1736 + (0.32e2 * t1498 * t1937 + t1649 * t2017) * t1631 + 0.8e1 * (t1496 * t1987 - t1506 * t1875) * pkin(1) + (0.16e2 * t1524 * t1658 + (t1652 * t1741) + (70 * t1739) + t1749 + t1780) * t1749;
t1675 = -t1731 / 0.6e1;
t1606 = 0.7e1 / 0.6e1 * t1736 + t1675 + t1905;
t1685 = t1736 / 0.3e1;
t1840 = t1675 + t1685 + t1743;
t1609 = t1689 + t1840;
t1970 = t1696 * pkin(1);
t1534 = -t1606 * t1970 + t1609 * t1650;
t1841 = t1731 / 0.3e1 + t1685 + t1725;
t1540 = (t1736 * t1656) - 0.5e1 / 0.3e1 * t1739 + t1841 * t1741 + t1743 * (t1677 + t1645);
t1667 = -0.20e2 / 0.3e1 * t1741;
t1684 = 0.2e1 / 0.3e1 * t1736;
t1842 = 0.2e1 / 0.3e1 * t1731 + t1684 + t1723;
t1843 = 0.4e1 / 0.3e1 * t1731 + t1683 + t1726;
t1542 = -t1749 + (t1667 + t1842) * t1736 - (3 * t1739) + t1843 * t1741 + t1742;
t1616 = t1741 + t1840;
t1641 = t1712 + t1656;
t1546 = t1616 * t1650 - t1641 * t1970 / 0.2e1;
t1934 = t1663 * t1739;
t1854 = t1696 * t1934;
t1880 = 0.4e1 * t1919;
t1978 = -t1696 / 0.2e1;
t1488 = t1854 * t1998 + t1534 * t1880 + (-0.8e1 / 0.3e1 * t1935 + t1540) * t1650 + (t1542 * t1978 + t1546 * t1915) * pkin(1);
t1716 = 5 * t1739;
t1722 = 6 * t1743;
t1901 = t1708 + t2027;
t1836 = t1722 + t1901;
t1702 = 10 * t1741;
t1902 = t1702 + t1725;
t1528 = t1749 + (t1679 + t1688 + t1902) * t1736 + t1716 + (t1836 * t1741) + t1743 * (t1679 + t1646);
t1521 = t1528 * t1650;
t1839 = t1677 + t1655;
t1613 = 0.8e1 / 0.3e1 * t1736 + t1839;
t1899 = t1718 + t1743;
t1615 = t1677 + t1684 + t1899;
t1538 = -t1613 * t1970 + t1615 * t1650;
t1620 = 0.5e1 / 0.6e1 * t1736 + t1691 + t1675;
t1853 = t1696 * t1924;
t1550 = pkin(1) * t1853 + t1620 * t1887;
t1872 = t1745 * t1650;
t1881 = -0.4e1 * t1919;
t1711 = 5 * t1749;
t1837 = t1679 + t1655;
t1545 = t1711 + ((t1702 + t1836) * t1736) + (t1688 + t1837) * t1655;
t1948 = t1545 * t1696;
t1999 = -0.8e1 * pkin(7);
t1492 = t1663 * t1872 * t1999 + t1550 * t1881 + t1521 + (t1538 * t1915 - t1948) * pkin(1);
t1503 = -pkin(1) * t1948 + t1521;
t1543 = -(3 * t1749) + (t1667 + t1843) * t1736 + t1842 * t1741 + t1893;
t2016 = 0.10e2 / 0.3e1 * t1749 + (-t1741 + t1841) * t2027 + (t1687 + t1839) * t1726;
t1504 = t1543 * t1650 + t1970 * t2016;
t2005 = t1724 - t1731 - t1736;
t1622 = t2005 * t1702;
t1900 = t1709 - 5 * t1736;
t1505 = t1733 + ((21 * t1741 + t2005) * t1749) + ((t1743 * t1901 + t1622 + t1720 + 35 * t1739) * t1736) + ((t1715 + (t1721 + t1900) * t1741 + t1743 * (-t1736 + t1891)) * t1655);
t1873 = t1741 * t1650;
t1778 = -t1696 * t1745 + t1873;
t1592 = 0.4e1 * t1778;
t1604 = t1650 + 0.2e1 * t1970;
t1639 = t1736 + t1899;
t1618 = t1678 + t1639;
t1507 = t1657 * t1650 + t1592 * t1664 + (t1604 * t1987 + t1618 * t1696) * t1968;
t1513 = 0.7e1 * t1733 + ((35 * t1741 + 15 * t1743 + t1900) * t1749) + ((21 * t1739 + t1622 + 9 * t1742 + (t1710 - 6 * t1736) * t1743) * t1736) + t1939;
t1607 = t1743 + 0.5e1 / 0.2e1 * t1736 + 0.3e1 / 0.2e1 * t1741 + t1678;
t1871 = (pkin(1) * t1654) / 0.2e1;
t1548 = t1607 * t1650 + t1696 * t1871;
t1633 = pkin(7) * t1889;
t1571 = 0.4e1 / 0.3e1 * t1919 + t1633 + t1647;
t1938 = t1658 * t1749;
t1818 = -0.24e2 * t1571 * t1938;
t1856 = t1631 * t1937;
t1825 = -0.8e1 * t1856;
t1870 = -0.12e2 * t1936;
t1882 = -0.6e1 * t1919;
t1941 = t1631 * t1695;
t1470 = t1548 * t1831 + t1507 * t1825 + t1488 * t1870 + t1504 * t1882 + (-0.6e1 * t1492 * t1941 + (0.24e2 * t1645 * t1935 - t1505) * t1697) * pkin(3) + (-0.6e1 * t1503 * t1987 + (t1513 + t1818) * t1696) * pkin(1);
t1593 = t1649 + t1631;
t1454 = t1463 * t1593 + t1470 * t1744;
t1559 = -0.4e1 / 0.9e1 * t1827 + 0.4e1 / 0.9e1 * t1736 - t1731 / 0.9e1 + t1906;
t1570 = t1675 + t1684 + t1790;
t1614 = t1683 + t1839;
t1497 = 0.4e1 * t1886 + 0.6e1 * t1559 * t1919 + t1614 * t1923 + (t1570 * t1915 + t1647 * t1832) * pkin(1);
t1838 = t1679 + t1684 + t1725;
t1909 = (t1684 + t1837) * t1655 + t1749;
t1512 = t1614 * t1816 + (t1717 + t1838) * t1736 + t1909;
t1564 = t1612 + t1614;
t1635 = t1895 * t2015;
t1888 = 0.4e1 * t1651;
t1501 = pkin(7) * t1564 * t1888 + t1635 * t1664 + t1512;
t1535 = t1614 * t1924 + t1583;
t1536 = (t1666 + t1838) * t1736 + t1909;
t1582 = t1633 + t1880 + t1657;
t1869 = 0.12e2 * t1936;
t1966 = 6 * pkin(1);
t1476 = t1497 * t1869 + t1602 + t1535 * t1868 + t1733 + ((-t1731 + t1736 + t1903) * t1749) + ((t1700 + (t1707 + t1710 + 6 * t1736) * t1741 + t1720 + (t1708 + t1712) * t1743) * t1736) + (t1648 * t1626) + (0.6e1 * t1501 * t1649 + 0.8e1 * t1582 * t1937) * t1631 + (t1512 * t1987 - t1536 * t1875) * t1966;
t1601 = t1639 * t1650;
t1605 = t1650 - t1970;
t1640 = 3 * t1736 + t1655;
t1940 = t1640 * t1696;
t1511 = -0.2e1 * t1664 * t1873 + t1601 + (0.2e1 * t1605 * t1987 - t1940) * pkin(1);
t1591 = 0.2e1 * t1778;
t1653 = t1736 + 2 * t1741;
t1514 = t1656 * t1650 + t1591 * t1664 + (t1653 * t1696 + t1698 * t1866) * pkin(1);
t1551 = -pkin(1) * t1940 + t1601;
t2014 = 8 * t1741;
t1629 = pkin(3) * t1991 + t1748 * t2014;
t1990 = 4 * t1745;
t1553 = t1629 * t1697 + t1853 * t1990;
t1576 = t1736 * t1902 + t1716 + t1892 + 6 * t1917;
t1587 = t1711 + (t1702 + t1722) * t1736 + t1648;
t1884 = -0.4e1 * t1936;
t1485 = t1514 * t1884 + t1553 * t1664 + (-0.4e1 * t1511 * t1941 + (-t1576 + t1845) * t1697) * pkin(3) + (-0.4e1 * t1551 * t1987 + (t1587 + t1825) * t1696) * pkin(1);
t1469 = t1476 * t1593 + t1485 * t1744;
t1467 = 0.1e1 / t1469;
t1985 = t1467 / 0.4e1;
t1849 = t1480 * t1985;
t1976 = -t1744 / 0.4e1;
t1771 = t1454 * t1849 + t1484 * t1976;
t1732 = 0.1e1 / pkin(4);
t1921 = t1732 / pkin(3) ^ 2;
t1857 = t1531 * t1921;
t1442 = t1771 * t1857;
t1848 = t1484 * t1985;
t1975 = t1744 / 0.4e1;
t1773 = t1454 * t1848 + t1480 * t1975;
t1443 = t1773 * t1857;
t1785 = t1926 + t1930;
t1421 = -t1442 * t1785 - t1443 * t1786;
t1422 = -t1442 * t1786 + t1443 * t1785;
t1977 = -t1698 / 0.2e1;
t1473 = (t1480 * t1978 + t1484 * t1977) * t1949;
t1474 = (t1480 * t1977 + t1696 * t1984) * t1949;
t1986 = t1467 / 0.2e1;
t1850 = t1454 * t1986;
t1922 = t1732 * t1737;
t1852 = t1744 * t1922;
t1775 = -t1473 * t1852 / 0.2e1 + t1474 * t1850 * t1922;
t2008 = (t1473 * t1850 + t1474 * t1744 / 0.2e1) * t1922;
t2023 = t1421 * t2008 + t1422 * t1775;
t2022 = -t1421 * t1775 + t1422 * t2008;
t1438 = t1459 * t1474 - t1461 * t1473;
t1824 = -t1459 * t1473 - t1474 * t1461;
t2021 = -t1438 * t1459 + t1461 * t1824;
t2020 = t1438 * t1461 + t1459 * t1824;
t2019 = t1469 ^ 2;
t2007 = qJD(1) - qJD(2);
t1802 = qJD(1) * t1828;
t1860 = t1631 * t1960;
t1822 = pkin(3) * t1860;
t2004 = t1802 - t1822;
t1997 = -4 * pkin(1);
t1993 = -0.8e1 * t1698;
t1992 = -0.2e1 * t1698;
t1983 = t1531 / 0.2e1;
t1768 = t2003 * pkin(3);
t1541 = t1741 * t1768 * t1915;
t1634 = pkin(1) * t1964;
t1588 = pkin(3) * t1960 - t1634;
t1834 = pkin(7) * t1634;
t1624 = -0.2e1 * t1834;
t1661 = t1696 ^ 2;
t1762 = pkin(1) * t1645 * t1768 * t1919;
t1764 = t2003 * t1968;
t1961 = qJD(2) * t1695;
t1877 = pkin(3) * t1961;
t1879 = pkin(1) * t1962;
t1767 = -t1528 * t1877 - t1545 * t1879;
t1789 = t1660 * t1748 * t1860;
t1779 = -0.24e2 * t1789;
t1781 = t1647 * t1800;
t1782 = t1647 * t1801;
t1932 = t1664 * t1745;
t1885 = pkin(7) * t1932;
t1799 = t1885 * t1964;
t1784 = -0.48e2 * t1799;
t1810 = t1741 * t1862;
t1794 = -0.24e2 * t1810;
t1795 = t1877 * t1933;
t1796 = t1877 * t1935;
t1863 = qJD(1) * t1937;
t1823 = pkin(1) * t1863;
t1797 = t1696 * t1823;
t1813 = t1659 * t1749 * t1960;
t1814 = t1924 * t1962;
t1815 = qJD(1) * t1854;
t1817 = -0.6e1 * t1834;
t1858 = qJD(2) * t1929;
t1821 = pkin(3) * t1858;
t1829 = t1919 * t1989;
t1855 = t1696 * t1938;
t1876 = pkin(3) * t1941;
t1783 = -0.24e2 * t1799;
t1912 = t1654 * t1783 + 0.24e2 * t1762;
t1913 = t2025 * t1544;
t1424 = (t1818 * t1879 - 0.24e2 * pkin(1) * (-0.8e1 / 0.3e1 * t1810 + t1624) * t1855 - 0.96e2 * t1571 * t1813 * t1970 + ((-t1657 + t1881) * t1877 + 0.2e1 * (pkin(1) * t1618 - t1592 * t1696 - 0.2e1 * t1932) * t1962 + (-0.2e1 * t1821 + (-0.2e1 * t1604 * t1696 + 0.4e1 * t1973) * qJD(1)) * pkin(1) * pkin(7)) * t1825 + (0.8e1 / 0.3e1 * t1796 + (-t1606 * t1879 - t1609 * t1877) * t1880 - t1540 * t1877 - t1542 * t1879 / 0.2e1 + (0.32e2 / 0.3e1 * t1934 * t1650 + t1741 * t1534 * t1993) * t1964 + (t1616 * t1821 * t1997 + ((0.12e2 * t1661 * t1664 - 0.4e1 * t1662) * t1739 + (-0.4e1 * t1546 * t1696 - 0.2e1 * t1641 * t1973) * pkin(1)) * qJD(1)) * pkin(7)) * t1870 - 0.24e2 * t1488 * t1811 - 0.6e1 * ((-0.4e1 * (pkin(1) * t1814 - 0.2e1 * t1620 * t1877) * t1664 + 0.8e1 * t1550 * t1862) * t1741 + (0.8e1 * t1795 + (-t1613 * t1879 - t1615 * t1877) * t1888 + (t1538 * t1997 + 0.24e2 * t1664 * t1872) * t1964) * pkin(7) + t1767) * t1876 + (-t1607 * t1877 + t1871 * t1962) * t1831 + t1548 * t1784 + (-t1543 * t1877 + t1879 * t2016) * t1882 + 0.12e2 * t1504 * t1810 - 0.6e1 * t1767 * t1890 + 0.6e1 * t1503 * t1834 + t1505 * t1877 + t1513 * t1879 + (-0.96e2 * t1815 * t1650 - 0.24e2 * t1796) * t1645 + (0.8e1 * t1797 + t1779) * t1507 + 0.6e1 * t2004 * t1492) * t1744 + t1470 * t1847 + (0.16e2 * (t1644 * t1992 - t1637 - 0.48e2 * t1885 - 0.32e2 * t1934) * qJD(1) * t1855 + 0.64e2 * t1524 * t1813 + 0.32e2 * (t1541 + (t1557 * t1992 + (t1563 * t1997 - 0.24e2 * t1932) * pkin(7)) * t1964 + (t1764 * t1923 - t1765 * t1932) * pkin(3)) * t1856 + 0.24e2 * (-0.32e2 * t1646 * t1815 + (-0.2e1 / 0.3e1 * t1861 + 0.2e1 / 0.3e1 * t1859) * t1831 * t1989 + t1569 * t1784 - 0.28e2 * t1508 * t1810 + (-0.8e1 / 0.3e1 * t1861 + 0.8e1 / 0.3e1 * t1859) * t1610 * pkin(3) * t1741 * t1914 + t1510 * t1817 + 0.4e1 * (-t1782 + t1781) * t1619 + 0.64e2 / 0.3e1 * t2003 * t1608 * t1829) * t1936 + 0.48e2 * t1486 * t1811 - 0.8e1 * t1487 * t1802 + 0.8e1 * (t1537 * t1794 + (-pkin(3) * t1619 * t1763 - t1515 * t1964) * pkin(7) * t1966 + t1912 + t1913) * t1876 + t1822 * t2017 - 0.4e1 * t1603 * t1663 * t1964 + pkin(3) * t1764 * t1830 * t1924 - 0.96e2 * t1539 * t1799 + 0.96e2 * t1619 * t1762 - 0.48e2 * t1502 * t1810 + 0.8e1 * t1913 * t1890 - 0.8e1 * t1496 * t1834 + 0.8e1 * t2024 * t1506 + (-0.32e2 * t1797 + 0.96e2 * t1789) * t1498) * t1593 + t1463 * t1588;
t1479 = 0.1e1 / t1480 ^ 2;
t1846 = -t1950 / 0.2e1;
t1432 = qJD(1) + ((t1464 * t1983 + t1484 * t1846) / t2013 - 0.2e1 * (t1465 * t1983 + t1480 * t1846) * t1484 * t1479) * pkin(3) / (t1479 * t1484 ^ 2 + 0.1e1) * t1533 * t1737;
t1952 = 0.1e1 / t1454 ^ 2 * t2019;
t1819 = t1640 * t1879;
t1826 = -0.2e1 * t1862;
t1910 = 0.4e1 * t2024 * t1614;
t1954 = ((t1631 * t1823 * t1993 + t1661 * t1863 * t2014 + t1779 * t1970 + (0.2e1 * (-t1741 * t1877 - t1745 * t1962) * t1664 + t1591 * t1826 - t1656 * t1877 + (t1653 * t1962 + (-qJD(1) * t1928 - t1858) * pkin(7) * t2018) * pkin(1)) * t1884 - 0.8e1 * t1514 * t1811 - 0.4e1 * (-t1819 + (0.4e1 * t1697 * t1810 + (-t1639 + 0.2e1 * t1919) * t1961) * pkin(3) + ((-t1877 - t1879) * t1651 - t1605 * t1634) * t1969) * t1876 + t1795 * t1999 + t1783 * t1650 + (-t1629 * t1961 + t1814 * t1990) * t1664 + t1553 * t1826 - 0.4e1 * (-t1639 * t1877 - t1819) * t1890 + 0.4e1 * t1551 * t1834 + t1576 * t1877 + t1587 * t1879 + 0.4e1 * t2004 * t1511) * t1744 + t1485 * t1847 + (0.8e1 * (t1624 - 0.8e1 * t1810) * t1856 + (-0.12e2 * t1799 + 0.6e1 * (-0.4e1 / 0.9e1 * t1861 + 0.4e1 / 0.9e1 * t1859) * t1829 - 0.12e2 * t1559 * t1810 + t1541 - 0.4e1 * t1570 * t1834 - 0.2e1 * t1782 + 0.2e1 * t1781) * t1869 + 0.24e2 * t1497 * t1811 + 0.6e1 * (t1635 * t1826 + (-t1564 * t1964 + t1698 * t1908) * pkin(7) * t1967 + t1910) * t1876 + t1535 * t1794 + t1910 * t1865 + t1512 * t1817 + t1912 + (-0.8e1 * t1797 + 0.24e2 * t1789) * t1582 + t2025 * t1536 - 0.6e1 * t2004 * t1501) * t1593 + t1476 * t1588) / t2019;
t1412 = t1432 + (0.1e1 / t1454 * t1469 * t1847 - 0.2e1 * (t1424 * t1986 - t1454 * t1954 / 0.2e1) * pkin(4) * pkin(3) * t1852 * t1952) / (t1495 * t1952 + 0.1e1);
t1972 = pkin(4) * t1412;
t1809 = t1432 + t2026;
t1414 = (t1445 * t1988 + t1699 * t1446) * t1459 * t1996 * t2002 + t1809;
t1971 = pkin(6) * t1414;
t1956 = t1432 * t1473;
t1955 = t1432 * t1474;
t1851 = -t1954 / 0.4e1;
t1806 = pkin(3) * t1955 - t1879;
t1805 = pkin(2) * t1955 - t1879;
t1804 = -pkin(2) * t1956 + t1634;
t1803 = -pkin(3) * t1956 + t1634;
t1788 = t1577 * t1698 - t1579 * t1696;
t1787 = t1577 * t1696 + t1579 * t1698;
t1600 = -rSges(2,1) * t1698 + rSges(2,2) * t1696;
t1599 = -rSges(8,1) * t1697 + rSges(8,2) * t1695;
t1598 = -rSges(2,1) * t1696 - rSges(2,2) * t1698;
t1597 = rSges(8,1) * t1695 + rSges(8,2) * t1697;
t1527 = t2007 * t1785;
t1526 = t2007 * t1786;
t1518 = t1577 * t1974 - t1693 * t1579;
t1517 = -t1693 * t1577 - t1579 * t1974;
t1500 = rSges(6,1) * t1518 - rSges(6,2) * t1517;
t1499 = rSges(6,1) * t1517 + rSges(6,2) * t1518;
t1489 = qJD(1) - t1491;
t1482 = -t1879 + t1489 * (rSges(9,1) * t1787 + rSges(9,2) * t1788);
t1481 = t1634 - t1489 * (-rSges(9,1) * t1788 + rSges(9,2) * t1787);
t1429 = -t1879 + t1432 * (rSges(3,1) * t1474 - rSges(3,2) * t1473);
t1428 = t1634 - t1432 * (rSges(3,1) * t1473 + rSges(3,2) * t1474);
t1420 = 0.1e1 / t1422 ^ 2;
t1416 = (-t1773 * t1950 + (t1465 * t1975 + t1480 * t1951 / 0.8e1 + t1424 * t1848 + (t1464 * t1985 + t1484 * t1851) * t1454) * t1531) * t1921;
t1415 = (-t1771 * t1950 + (t1424 * t1849 + t1464 * t1976 - t1484 * t1951 / 0.8e1 + (t1465 * t1985 + t1480 * t1851) * t1454) * t1531) * t1921;
t1411 = -t1879 + t1414 * (rSges(7,1) * t1824 - rSges(7,2) * t1438);
t1410 = t1634 - t1414 * (rSges(7,1) * t1438 + rSges(7,2) * t1824);
t1409 = t1414 * (rSges(4,1) * t1824 - rSges(4,2) * t1438) + t1805;
t1408 = -t1414 * (rSges(4,1) * t1438 + rSges(4,2) * t1824) + t1804;
t1407 = t1412 * (rSges(5,1) * t1775 - rSges(5,2) * t2008) + t1806;
t1406 = -t1412 * (rSges(5,1) * t2008 + rSges(5,2) * t1775) + t1803;
t1405 = t1809 - t2026;
t1404 = t1824 * t1971 + t1405 * (rSges(10,1) * t2021 - rSges(10,2) * t2020) + t1805;
t1403 = -t1438 * t1971 - t1405 * (rSges(10,1) * t2020 + rSges(10,2) * t2021) + t1804;
t1402 = (-(-t1415 * t1785 - t1416 * t1786 - t1442 * t1526 + t1443 * t1527) / t1422 - (t1415 * t1786 - t1416 * t1785 - t1442 * t1527 - t1443 * t1526) * t1421 * t1420) / (t1420 * t1421 ^ 2 + 0.1e1) + t1412;
t1401 = t1775 * t1972 + t1402 * (rSges(11,1) * t2023 - rSges(11,2) * t2022) + t1806;
t1400 = -t2008 * t1972 - t1402 * (rSges(11,1) * t2022 + rSges(11,2) * t2023) + t1803;
t1 = t1402 ^ 2 * Icges(11,3) / 0.2e1 + t1432 ^ 2 * Icges(3,3) / 0.2e1 + t1412 ^ 2 * Icges(5,3) / 0.2e1 + t1489 ^ 2 * Icges(9,3) / 0.2e1 + t1405 ^ 2 * Icges(10,3) / 0.2e1 + m(3) * (t1428 ^ 2 + t1429 ^ 2) / 0.2e1 + m(9) * (t1481 ^ 2 + t1482 ^ 2) / 0.2e1 + m(7) * (t1410 ^ 2 + t1411 ^ 2) / 0.2e1 + m(11) * (t1400 ^ 2 + t1401 ^ 2) / 0.2e1 + m(10) * (t1403 ^ 2 + t1404 ^ 2) / 0.2e1 + m(5) * (t1406 ^ 2 + t1407 ^ 2) / 0.2e1 + m(4) * (t1408 ^ 2 + t1409 ^ 2) / 0.2e1 + (Icges(7,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t1414 ^ 2 + (m(2) * (t1598 ^ 2 + t1600 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (Icges(8,3) / 0.2e1 + m(8) * (t1597 ^ 2 + t1599 ^ 2) / 0.2e1) * qJD(2) ^ 2 + (m(6) * (t1499 ^ 2 + t1500 ^ 2) / 0.2e1 + Icges(6,3) / 0.2e1) * t1491 ^ 2;
T = t1;
