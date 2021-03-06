% Calculate joint inertia matrix with Newton Euler for
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_inertiaJ_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_inertiaJ_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_inertiaJ_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2IC_inertiaJ_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2IC_inertiaJ_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:48:55
% EndTime: 2020-05-02 23:49:07
% DurationCPUTime: 12.31s
% Computational Cost: add. (3441->558), mult. (3904->673), div. (144->13), fcn. (2258->40), ass. (0->330)
t1854 = pkin(11) * m(6) + mrSges(6,3);
t2016 = -t1854 * pkin(9) - Ifges(5,4);
t1718 = sin(qJ(7));
t1726 = cos(qJ(7));
t1710 = sin(pkin(19));
t1711 = cos(pkin(19));
t1615 = Ifges(11,5) * t1710 - Ifges(11,6) * t1711;
t1616 = Ifges(11,5) * t1711 + Ifges(11,6) * t1710;
t1712 = sin(qJ(10));
t1713 = cos(qJ(10));
t1872 = t1712 * t1615 + t1616 * t1713;
t1935 = pkin(4) * t1711;
t1970 = mrSges(11,3) * t1935 - Ifges(8,5) + t1872;
t1534 = -t1615 * t1713 + t1616 * t1712;
t1936 = pkin(4) * t1710;
t1999 = -mrSges(11,3) * t1936 - Ifges(8,6) + t1534;
t2015 = t1970 * t1718 + t1999 * t1726;
t1716 = sin(qJ(9));
t1724 = cos(qJ(9));
t1631 = Ifges(10,5) * t1724 - Ifges(10,6) * t1716;
t1598 = mrSges(10,3) * pkin(2) - Ifges(9,5) + t1631;
t1630 = Ifges(10,5) * t1716 + Ifges(10,6) * t1724;
t1614 = -Ifges(9,6) + t1630;
t1717 = sin(qJ(8));
t1725 = cos(qJ(8));
t1969 = t1598 * t1717 + t1614 * t1725;
t1722 = sin(qJ(3));
t1728 = cos(qJ(5));
t1702 = t1728 ^ 2;
t1715 = Ifges(6,1) - Ifges(6,2);
t1720 = sin(qJ(5));
t1928 = mrSges(6,1) * t1720;
t1946 = pkin(9) * mrSges(6,2);
t1545 = pkin(9) * t1928 - 0.2e1 * Ifges(6,4) * t1702 + t1728 * (-t1715 * t1720 + t1946) + Ifges(6,4) - Ifges(5,5);
t1943 = pkin(11) * mrSges(6,2);
t1679 = -Ifges(6,6) + t1943;
t1944 = pkin(11) * mrSges(6,1);
t1680 = -Ifges(6,5) + t1944;
t1978 = -t1679 * t1728 - t1680 * t1720;
t1569 = Ifges(5,6) - t1978;
t1729 = cos(qJ(4));
t1557 = t1569 * t1729;
t1721 = sin(qJ(4));
t1505 = -t1721 * t1545 + t1557;
t2009 = -Ifges(4,6) - t1505;
t2011 = t2009 * t1722;
t2014 = -t2011 + Ifges(3,6) - t1969 - t2015;
t1576 = -t1679 * t1720 + t1680 * t1728;
t2012 = t1576 - t2016;
t1556 = t1721 * t1569;
t1927 = mrSges(6,2) * t1720;
t1977 = mrSges(6,1) * t1728 - t1927;
t1609 = pkin(9) * t1977 + Ifges(6,3);
t1900 = t1576 * t1721;
t2008 = pkin(5) * t1977 + t1609 * t1729 + t1900;
t1762 = pkin(9) ^ 2;
t1736 = t1762 * m(6);
t1966 = Ifges(5,2) + Ifges(6,3) + t1736;
t1659 = t1715 * t1702;
t1942 = pkin(11) * mrSges(6,3);
t1693 = 0.2e1 * t1942;
t1760 = pkin(11) ^ 2;
t1735 = t1760 * m(6);
t1855 = -Ifges(5,1) - Ifges(6,2) - t1735;
t1976 = t1693 - t1855 + t1659;
t2007 = t1976 - t1966;
t1687 = t1722 * pkin(1);
t1925 = mrSges(10,2) * t1716;
t1803 = -mrSges(10,1) * t1724 + t1925;
t1610 = pkin(2) * t1803 + Ifges(10,3);
t1837 = -qJ(8) + pkin(17) + qJ(3);
t1652 = cos(qJ(9) - t1837);
t1845 = pkin(6) / t1716 / pkin(2);
t1903 = (-pkin(12) * t1652 + pkin(2) * cos(t1837)) / pkin(12);
t1801 = t1845 * t1903;
t1740 = pkin(2) ^ 2 * m(10);
t1863 = pkin(2) * t1925;
t1869 = 0.2e1 * t1863 + t1740;
t1948 = pkin(2) * mrSges(10,1);
t1915 = -0.2e1 * t1948;
t1805 = t1724 * t1915 + Ifges(9,3) + Ifges(10,3) + t1869;
t1875 = t1652 * t1805 * t1845 + t1610 * t1801;
t1599 = -pkin(9) * m(6) - mrSges(5,1) - t1977;
t1896 = t1599 * t1729;
t1672 = -mrSges(5,2) + t1854;
t1639 = t1672 * t1721;
t1974 = -mrSges(4,1) - (m(5) + m(6)) * pkin(5) - t1639;
t1898 = t1722 * t1974;
t1730 = cos(qJ(3));
t1640 = t1672 * t1729;
t1884 = t1721 * t1599;
t1871 = t1640 + t1884;
t1971 = -mrSges(4,2) + t1871;
t2001 = t1971 * t1730;
t2006 = t1875 - (t1687 + 0.2e1 * pkin(5)) * t1896 + (-t2001 - t1898) * pkin(1);
t1955 = -0.2e1 * pkin(1);
t1841 = -m(11) - m(4) - m(5) - m(8);
t1673 = -m(6) + t1841;
t1770 = pkin(1) ^ 2;
t1618 = mrSges(11,1) * t1711 + mrSges(11,2) * t1710;
t1554 = m(11) * t1936 + t1618 * t1712 - mrSges(8,2);
t1939 = pkin(1) * t1718;
t1843 = t1554 * t1939;
t1961 = -t1673 * t1770 + t1898 * t1955 + 0.2e1 * t1843;
t2005 = -Ifges(3,2) - t1961;
t1634 = pkin(5) * t1639;
t1765 = pkin(5) ^ 2;
t1975 = t1765 * m(5) + 0.2e1 * t1634;
t2004 = -Ifges(4,2) - t1975;
t2003 = t1545 * t1729;
t1705 = Ifges(11,1) - Ifges(11,2);
t1690 = t1711 ^ 2;
t1918 = Ifges(11,4) * t1690;
t1827 = -Ifges(11,4) + 0.2e1 * t1918;
t1887 = t1710 * t1711;
t1571 = -t1705 * t1887 + t1827;
t1901 = t1571 * t1712;
t1870 = mrSges(11,1) * t1690 + mrSges(11,2) * t1887;
t1929 = (-mrSges(11,1) / 0.2e1 + t1870) * pkin(4);
t1873 = 0.4e1 * t1901 - 0.4e1 * t1929;
t2002 = t1873 * t1713;
t1861 = pkin(9) * t1927;
t1924 = Ifges(6,4) * t1720;
t1891 = (pkin(9) * mrSges(6,1) + t1924) * t1728;
t1775 = -0.4e1 * t1891 + 0.4e1 * t1861 + 0.2e1 * t2007;
t1738 = t1765 * m(6);
t1957 = Ifges(4,1) - t1738 + t2004;
t2000 = 0.2e1 * t1957 - t1775;
t1732 = pkin(4) ^ 2 * m(11);
t1646 = mrSges(11,1) * t1887;
t1579 = t1646 + (-t1690 + 0.1e1 / 0.2e1) * mrSges(11,2);
t1914 = Ifges(11,1) + t1732;
t1666 = -Ifges(11,2) + t1914;
t1890 = t1666 * t1690;
t1934 = pkin(4) * t1712;
t1960 = -0.4e1 * t1579 * t1934 - Ifges(8,1) + Ifges(8,2) + 0.2e1 * t1890;
t1838 = Ifges(11,4) * t1887;
t1637 = 0.4e1 * t1838;
t1963 = t1637 - t1705;
t1994 = 0.2e1 * t1963;
t1997 = 0.2e1 * t1732 - 0.2e1 * t1960 - 0.2e1 * t2002 - t1994;
t1888 = t1705 * t1690;
t1996 = t1994 + 0.4e1 * t1888;
t1995 = t1576 * t1729 - t1609 * t1721;
t1926 = mrSges(6,2) * t1728;
t1633 = t1926 + t1928;
t1794 = mrSges(5,3) + t1633;
t1993 = t1794 * pkin(5) - Ifges(4,5);
t1691 = t1713 ^ 2;
t1781 = 0.2e1 * t1888 + t1963;
t1992 = t1691 * t1781;
t1991 = t1871 * t1730;
t1989 = t1718 * t1999 - t1970 * t1726;
t1723 = sin(qJ(2));
t1731 = cos(qJ(2));
t1821 = -t1598 * t1725 + t1614 * t1717;
t1986 = (t1723 * t1969 + t1821 * t1731) * t1652;
t1912 = 0.2e1 * pkin(1);
t1668 = t1687 + pkin(5);
t1714 = Ifges(10,1) - Ifges(10,2);
t1947 = pkin(2) * mrSges(10,2);
t1984 = (t1714 * t1716 + t1947) * t1724;
t1562 = -t1884 + mrSges(4,2);
t1983 = (-t1562 + t1640) * t1912;
t1617 = mrSges(11,1) * t1710 - mrSges(11,2) * t1711;
t1590 = t1617 * t1712;
t1816 = -m(11) * t1935 - mrSges(8,1);
t1553 = t1590 + t1816;
t1592 = t1618 * t1713;
t1824 = t1553 + t1592;
t1982 = t1824 * t1912;
t1981 = -t1617 * t1713 + t1554;
t1880 = qJ(4) + pkin(18);
t1881 = pkin(19) + qJ(3);
t1812 = t1880 + t1881;
t1920 = qJ(6) - qJ(2);
t1791 = t1812 + t1920;
t1780 = -0.2e1 * qJ(7) - pkin(20) + t1791;
t1792 = t1812 - t1920;
t1785 = pkin(20) + t1792;
t1980 = cos(qJ(10) - t1780) + cos(qJ(10) - t1785);
t1979 = t1592 + t1590;
t1973 = -0.2e1 * t1861 + 0.2e1 * t1891;
t1967 = -t1993 - t1556;
t1964 = Ifges(4,3) + t1975;
t1595 = t1710 * t1713 - t1711 * t1712;
t1596 = t1710 * t1712 + t1711 * t1713;
t1628 = t1936 + t1939;
t1938 = pkin(1) * t1726;
t1629 = t1935 + t1938;
t1507 = -t1595 * t1628 - t1596 * t1629;
t1508 = t1595 * t1629 - t1596 * t1628;
t1485 = mrSges(11,1) * t1507 - mrSges(11,2) * t1508 + Ifges(11,3);
t1518 = (-t1595 * t1710 - t1596 * t1711) * pkin(4);
t1519 = (t1595 * t1711 - t1596 * t1710) * pkin(4);
t1488 = mrSges(11,1) * t1518 - mrSges(11,2) * t1519 + Ifges(11,3);
t1767 = 0.1e1 / pkin(3);
t1879 = qJ(7) + pkin(20);
t1815 = t1879 - t1920;
t1653 = cos(t1815);
t1644 = 0.1e1 / t1653;
t1899 = (-pkin(3) * t1653 - pkin(1) * cos(t1920)) * t1644;
t1833 = t1767 * t1899;
t1959 = t1488 * t1833 + t1485;
t1953 = 0.2e1 * t1722;
t1952 = 0.2e1 * t1723;
t1941 = mrSges(11,1) * pkin(4);
t1940 = -t1722 / 0.2e1;
t1937 = pkin(4) * t1579;
t1933 = pkin(5) * t1599;
t1932 = pkin(5) * t1672;
t1675 = sin(t1880);
t1931 = pkin(5) * t1675;
t1930 = pkin(5) * t1721;
t1923 = Ifges(10,4) * t1716;
t1681 = t1716 * mrSges(10,1);
t1917 = Ifges(11,3) / pkin(8) ^ 2;
t1913 = qJ(10) + qJ(7);
t1787 = -qJ(7) + t1791;
t1788 = -qJ(7) + t1792;
t1483 = (-(cos(t1788) + cos(t1787)) * pkin(1) + (-cos(t1780) - cos(t1785)) * pkin(3)) * pkin(4) + ((cos(qJ(10) - t1788) + cos(qJ(10) - t1787)) * pkin(1) + t1980 * pkin(3)) * pkin(8);
t1763 = 0.1e1 / pkin(8);
t1911 = t1483 * t1763;
t1830 = -t1659 + t1973;
t1520 = -0.2e1 * t1942 + t1830 + t1855 + t1966;
t1703 = t1729 ^ 2;
t1910 = t1520 * t1703;
t1840 = -t1928 / 0.2e1;
t1532 = ((Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t1720 - t1946 / 0.2e1) * t1728 + pkin(9) * t1840 + Ifges(5,5) / 0.2e1 + (t1702 - 0.1e1 / 0.2e1) * Ifges(6,4);
t1527 = t1532 * t1729;
t1533 = 0.1e1 / t1980;
t1909 = t1533 * t1767;
t1906 = t2012 * t1703;
t1619 = cos(t1812 - t1913);
t1761 = 0.1e1 / pkin(10);
t1905 = (-pkin(10) * t1619 - pkin(5) * cos(t1881 - t1913)) * t1761;
t1852 = t1617 * t1939;
t1904 = (-t1852 - t1941) * t1713;
t1902 = t1571 * t1691;
t1699 = t1724 ^ 2;
t1658 = t1714 * t1699;
t1885 = t1717 * t1725;
t1883 = t1721 * t1722;
t1882 = t1729 * t1730;
t1498 = t2009 * t1730;
t1674 = sin(t1879);
t1719 = sin(qJ(6));
t1727 = cos(qJ(6));
t1878 = (-t1498 + Ifges(3,5) + (-mrSges(11,3) - mrSges(4,3) - mrSges(8,3) - t1794) * pkin(1) + t1821 + (0.2e1 * t1527 + t1967) * t1722 + t1989) * t1731 + pkin(1) * t1674 / pkin(7) * t1644 * (Ifges(7,5) * t1719 + Ifges(7,6) * t1727);
t1605 = -t1718 * t1723 + t1726 * t1731;
t1606 = t1718 * t1731 + t1723 * t1726;
t1496 = -t1595 * t1605 + t1596 * t1606;
t1497 = -t1595 * t1606 - t1596 * t1605;
t1484 = Ifges(11,5) * t1497 + Ifges(11,6) * t1496;
t1877 = -t1721 * t1775 - 0.2e1 * t1932;
t1540 = 0.4e1 * t2012;
t1876 = -t1540 * t1721 + 0.2e1 * t1933;
t1504 = t1556 + t2003;
t1874 = t1996 * t1712 - 0.4e1 * t1937;
t1868 = Ifges(6,5) * t1720 + Ifges(6,6) * t1728;
t1867 = t1760 + t1762;
t1866 = 0.2e1 * t1687;
t1865 = -0.2e1 * t1923;
t1862 = pkin(2) * t1681;
t1860 = mrSges(11,2) * t1934;
t1859 = 0.4e1 * t1902;
t1856 = Ifges(10,1) / 0.2e1 - Ifges(10,2) / 0.2e1;
t1851 = t1618 * t1939;
t1850 = pkin(4) * (pkin(1) * (sin(qJ(10) - t1920) + sin(qJ(10) + t1920)) + (sin(qJ(10) - t1815) + sin(qJ(10) + t1815)) * pkin(3)) * t1761;
t1848 = pkin(5) * t1884;
t1847 = pkin(5) * t1896;
t1846 = t1763 * t1931;
t1844 = t1824 * t1938;
t1842 = t1712 * t1929;
t1836 = t1483 * t1909;
t1835 = ((t1856 * t1716 + t1947 / 0.2e1) * t1724 + t1862 / 0.2e1 + Ifges(9,4) / 0.2e1 + (t1699 - 0.1e1 / 0.2e1) * Ifges(10,4)) * t1885;
t1613 = 0.1e1 / t1619;
t1834 = t1613 * t1905;
t1832 = t1668 * t1896;
t1831 = t1666 * t1887;
t1820 = -t1974 - t1896;
t1817 = Ifges(11,3) + Ifges(8,3) + 0.2e1 * t1860 + t1732;
t1813 = t1763 * t1836;
t1811 = mrSges(10,2) * t1724 - mrSges(9,2) + t1681;
t1802 = t1850 * t1909;
t1798 = 0.2e1 * t2012;
t1797 = (-0.4e1 * t1712 * t1781 + 0.8e1 * t1937) * t1713 - 0.8e1 * t1902;
t1796 = t1639 - t1896;
t1607 = -t1882 + t1883;
t1793 = Ifges(6,1) + Ifges(5,3) + t1693 + t1830;
t1784 = -0.2e1 * t1848 - (2 * Ifges(4,4)) + t1798;
t1783 = -Ifges(9,1) + Ifges(9,2) + (t1915 + 0.4e1 * t1923) * t1724 - 0.2e1 * t1658 + t1869 + t1714;
t1782 = m(10) * pkin(2) + mrSges(9,1) + t1803;
t1525 = t1735 + t1736 + t1793;
t1778 = -t1996 * t1691 + t1637 - t1666 + t1960;
t1777 = (t1765 + t1867) * m(6) + t1793 + t1964;
t1771 = -t1973 + t2007;
t1776 = 0.4e1 * (t1721 * t1771 + t1932) * t1729 - t1540 + 0.4e1 * t1848 + (4 * Ifges(4,4));
t1773 = t1738 + t1525 + t1964;
t1772 = t1775 * t1703 + t1520 + t1957;
t1758 = pkin(15) ^ 2;
t1704 = t1730 ^ 2;
t1701 = t1726 ^ 2;
t1700 = t1725 ^ 2;
t1608 = t1721 * t1730 + t1722 * t1729;
t1604 = t1716 * t1725 + t1717 * t1724;
t1603 = t1716 * t1717 - t1724 * t1725;
t1560 = -pkin(11) * t1633 + t1868;
t1539 = 0.8e1 * t2012;
t1529 = -t1607 * t1723 + t1608 * t1731;
t1528 = -t1607 * t1731 - t1608 * t1723;
t1514 = Ifges(8,4) - t1827 + t1831 + 0.2e1 * t1842;
t1513 = m(11) * t1519 - mrSges(11,2);
t1512 = m(11) * t1518 + mrSges(11,1);
t1501 = m(11) * t1508 - mrSges(11,2);
t1500 = m(11) * t1507 + mrSges(11,1);
t1499 = 0.2e1 * t1532 * t1721 + t1557;
t1495 = t1527 - t1556 / 0.2e1;
t1492 = t1525 + t1634 - t1847;
t1491 = t1504 + t1993;
t1490 = t1491 * t1730;
t1489 = t1527 + ((-t1943 / 0.2e1 + Ifges(6,6) / 0.2e1) * t1728 + (-t1944 / 0.2e1 + Ifges(6,5) / 0.2e1) * t1720 - Ifges(5,6) / 0.2e1) * t1721 + Ifges(4,5) / 0.2e1 + (-t1926 / 0.2e1 + t1840 - mrSges(5,3) / 0.2e1) * pkin(5);
t1486 = -pkin(1) * t1991 + t1668 * t1796 + t1525;
t1482 = pkin(5) * t1796 + t1525 * t1834 + t1525;
t1481 = (t1718 * t1534 - t1726 * t1872) * t1731 + t1723 * (t1534 * t1726 + t1718 * t1872);
t1480 = -t1832 + t1634 + t1867 * m(6) - t1525 * t1802 + (t1672 * t1883 - t1991) * pkin(1) + t1793;
t1479 = (-t1504 * t1722 + t1505 * t1730) * t1731 - (t1504 * t1730 + t1505 * t1722) * t1723;
t1478 = (t1495 * t1953 + t1499 * t1730) * t1731 + (t1495 * t1730 + t1499 * t1940) * t1952;
t1476 = -t1979 * t1938 + t1904 + t1712 * t1851 + t1860 + Ifges(11,3) + (Ifges(11,3) + (t1711 * (-mrSges(11,1) * t1596 - mrSges(11,2) * t1595) + t1710 * (-mrSges(11,1) * t1595 + mrSges(11,2) * t1596)) * pkin(4)) * t1833 + Ifges(11,3) * t1813;
t1 = [t1914 + (t1876 * t1729 + t1772) * t1704 - 0.2e1 * (-t2001 + (-t1816 - t1979) * t1726 + t1782 * t1725 - t1673 * pkin(1) + mrSges(3,1) + t1820 * t1722 + t1981 * t1718 + t1811 * t1717) * pkin(15) * t1723 - 0.2e1 * t1838 + 0.2e1 * (-t1852 - t1901 + pkin(4) * (-mrSges(11,1) + t1870)) * t1713 - 0.2e1 * t1728 * t1924 + (-t1982 + (t1713 * t1874 + 0.2e1 * t1514 + t1859) * t1718) * t1726 + (-t1983 + (t1877 * t1729 + t1784 - 0.4e1 * t1906) * t1722) * t1730 + 0.2e1 * (-mrSges(7,1) * pkin(14) + Ifges(7,4) * t1719) * t1727 + (-0.2e1 * t1668 * t1599 + t1798 * t1721) * t1729 + 0.4e1 * t1835 + (t1758 + t1765) * m(6) + 0.2e1 * (-mrSges(11,2) * t1690 + mrSges(11,2) + t1646) * t1934 + (((t1896 * t1912 - t1974 * t1955 + (-0.4e1 * t1910 + (0.8e1 * t2016 * t1721 - 0.8e1 * t1900 + 0.4e1 * t1933) * t1729 + t2000) * t1722) * t1730 + (t1539 * t1703 + t1776) * t1704 + (-0.8e1 * t1842 - 0.4e1 * t1831 - 0.4e1 * Ifges(8,4) + (0.8e1 * t1690 - 0.4e1) * Ifges(11,4) + t1797) * t1701 - 0.4e1 * (t1658 + (t1865 + t1948) * t1724 - t1863 - t1740 / 0.2e1 + Ifges(9,1) / 0.2e1 - Ifges(9,2) / 0.2e1 - t1856) * t1885 + t1562 * t1866 + (-t1981 * t1912 + (-0.8e1 * (t1888 + 0.2e1 * t1838 - Ifges(11,1) / 0.2e1 + Ifges(11,2) / 0.2e1) * t1691 - t1997) * t1718) * t1726 + (-0.2e1 * t1851 + t1874) * t1713 + (-0.2e1 * t1672 * t1687 + t1877) * t1729 + 0.2e1 * t1984 + t1859 - 0.2e1 * t1553 * t1939 + 0.2e1 * t1831 + 0.4e1 * t1842 + 0.2e1 * t1862 + t1784 - 0.4e1 * t1918 + 0.4e1 * (-t1984 - t1862 - Ifges(9,4)) * t1700 - (2 * Ifges(3,4)) + 0.2e1 * Ifges(8,4) + 0.2e1 * Ifges(9,4) + 0.2e1 * Ifges(11,4) - t1540 * t1703 + ((-0.8e1 * t1699 + 0.4e1) * t1700 - 0.2e1 + 0.4e1 * t1699) * Ifges(10,4)) * t1723 + 0.2e1 * (-t1782 * t1717 + t1824 * t1718 + t1971 * t1722 + t1811 * t1725 + t1981 * t1726 + t1820 * t1730 - mrSges(3,2)) * pkin(15) + (((t1539 * t1721 - 0.4e1 * t1933) * t1729 - 0.4e1 * t1771 * t1703 - t2000) * t1704 + (t1983 + (t1776 + 0.8e1 * t1906) * t1722) * t1730 + (0.4e1 * t1658 + (-0.8e1 * t1923 + 0.4e1 * t1948) * t1724 - 0.4e1 * t1863 - 0.2e1 * t1740 + 0.2e1 * Ifges(9,1) - 0.2e1 * Ifges(10,1) - 0.2e1 * Ifges(9,2) + 0.2e1 * Ifges(10,2)) * t1700 + (t1982 + (-0.4e1 * t1514 + t1797) * t1718) * t1726 - 0.8e1 * t1835 + (0.2e1 * t1852 + t1873) * t1713 + (t1599 * t1866 + t1876) * t1729 + (0.4e1 * t1992 + t1997) * t1701 + t1772 + t1778 + t1783 + Ifges(3,1) + t2005) * t1731) * t1731 + t1724 * t1865 + (-Ifges(7,1) + Ifges(7,2)) * t1727 ^ 2 + t1976 - t1890 + t1783 * t1700 + 0.2e1 * mrSges(7,2) * t1719 * pkin(14) - t2004 - t2005 + t1910 + t1658 + t1992 + (t1778 + t2002) * t1701 + Ifges(7,1) + Ifges(8,1) + Ifges(9,1) + Ifges(10,2) + Ifges(2,3) + pkin(14) ^ 2 * m(7) + (m(9) + m(10) + m(3) - t1841) * t1758, t1723 * ((t1967 - t2003) * t1730 - t2014) + ((t1723 * t2015 + t1989 * t1731) * t1899 + (-t1478 * t1850 + t1481 * t1911) * t1533) * t1767 + t1878, (t1489 * t1953 - t1498) * t1731 + (t1489 * t1730 - t1940 * t2009) * t1952 + (t1478 * t1905 + t1481 * t1846) * t1613 + (t1986 + ((t1630 * t1717 - t1631 * t1725) * t1731 + t1723 * (t1630 * t1725 + t1631 * t1717)) * t1903) * t1845, (-t1995 * t1722 - t2008 * t1730) * t1731 + (pkin(1) * t1977 + t2008 * t1722 - t1995 * t1730) * t1723 - pkin(15) * t1977; -(t1490 + t2014) * t1723 + ((-t1479 * t1850 + t1484 * t1911) * t1533 + (Ifges(8,5) * t1605 - Ifges(8,6) * t1606 + t1484 + (t1711 * (t1496 * t1595 + t1497 * t1596) + t1710 * (-t1496 * t1596 + t1497 * t1595)) * pkin(4) * mrSges(11,3)) * t1899) * t1767 + t1878, -0.2e1 * t1832 + t2001 * t1955 - 0.2e1 * t1844 + 0.2e1 * t1904 + t1770 * t1674 ^ 2 / pkin(7) ^ 2 / t1653 ^ 2 * Ifges(7,3) + t1817 + t1805 + t1777 + Ifges(3,3) + (Ifges(8,3) + (t1711 * (-t1500 * t1596 + t1501 * t1595) + t1710 * (-t1500 * t1595 - t1501 * t1596)) * pkin(4) + (mrSges(8,1) * t1726 - mrSges(8,2) * t1718) * pkin(1) + t1485 - t1844 + (-t1852 - 0.2e1 * t1941) * t1713 + t1843 + t1817 + (Ifges(8,3) + (t1711 * (-t1512 * t1596 + t1513 * t1595) + t1710 * (-t1512 * t1595 - t1513 * t1596)) * pkin(4) + t1488) * t1833) * t1833 + (t1476 + t1959) * t1813 + (-t1486 - t1480) * t1802 + t1961, (t1476 * t1846 + t1480 * t1905) * t1613 - t1492 * t1802 + t1777 + t2006, -t1978 * t1802 + (-pkin(1) * t1607 - t1930) * t1633 + t1978; (-t1491 * t1722 - t1498) * t1731 - (t1490 - t2011) * t1723 + (t1479 * t1905 + t1484 * t1846) * t1613 + (t1986 + (Ifges(10,5) * (t1603 * t1731 + t1604 * t1723) + Ifges(10,6) * (-t1603 * t1723 + t1604 * t1731)) * t1903) * t1845, -t1482 * t1802 + t1773 + (t1486 * t1905 + (t1763 * t1959 + t1836 * t1917) * t1931) * t1613 + t2006, -0.2e1 * t1847 + (t1875 * t1652 + (Ifges(10,3) * t1903 + t1610 * t1652) * t1801) * t1845 + (t1482 + t1492) * t1834 + t1773 + t1765 * t1675 ^ 2 / t1619 ^ 2 * t1917, -t1633 * t1930 + t1834 * t1978 + t1978; -Ifges(6,3) * t1528 + (Ifges(6,5) * t1728 - Ifges(6,6) * t1720) * t1529 + t1977 * (-pkin(15) + t1723 * pkin(1) - (-t1722 * t1723 + t1730 * t1731) * pkin(5) - t1529 * pkin(11) - t1528 * pkin(9)), -t1560 * t1802 - t1633 * (-pkin(1) * t1882 + t1668 * t1721 + pkin(11)) + t1868, t1560 * t1834 - t1633 * (pkin(11) + t1930) + t1868, Ifges(6,3);];
Mq = t1;
