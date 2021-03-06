% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fivebar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x2]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = fivebar1DE1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_jacobiR_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fivebar1DE1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
JR_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0; t9, 0; 0, 0; -t9, 0; -t8, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:53:14
	% EndTime: 2020-04-27 03:53:47
	% DurationCPUTime: 33.81s
	% Computational Cost: add. (604000->629), mult. (1814472->979), div. (4384->9), fcn. (217384->9), ass. (0->452)
	t1794 = pkin(5) ^ 2;
	t1793 = t1794 ^ 2;
	t1797 = pkin(4) ^ 2;
	t1796 = t1797 ^ 2;
	t1961 = t1793 - t1796;
	t2052 = 4 * pkin(3);
	t2041 = 2 * pkin(1);
	t2008 = 2 * pkin(3);
	t1801 = pkin(3) ^ 2;
	t1732 = 0.10e2 / 0.3e1 * t1801;
	t1805 = pkin(2) ^ 2;
	t1744 = -0.2e1 / 0.3e1 * t1794;
	t1749 = 0.2e1 / 0.3e1 * t1797;
	t1807 = pkin(1) ^ 2;
	t1791 = 2 * t1807;
	t1875 = t1744 + t1749 + t1791;
	t1719 = t1801 + t1807;
	t1813 = t1805 ^ 2;
	t1874 = t1744 + t1719;
	t1974 = t1719 * (t1749 + t1874) + t1813;
	t1613 = (t1732 + t1875) * t1805 + t1974;
	t2051 = -0.6e1 * t1613;
	t1799 = t1801 ^ 2;
	t1806 = t1807 ^ 2;
	t1959 = t1799 + t1806;
	t1962 = t1791 - t1794;
	t1981 = t1807 * t1794;
	t2042 = -t1793 / 0.6e1 + t1796 / 0.6e1;
	t1639 = t1962 * t1801 + t1959 - t1981 - t2042;
	t1825 = t1639 + t1813;
	t1615 = (t1732 + t1962) * t1805 + t1825;
	t2050 = -0.6e1 * t1615;
	t1762 = sin(qJ(1));
	t1717 = pkin(2) * t1762;
	t1951 = 0.2e1 * t1717;
	t2031 = 4 * t1801;
	t1781 = 6 * t1801;
	t1742 = -t1794 / 0.3e1;
	t1750 = t1797 / 0.3e1;
	t1957 = t1805 + t1807;
	t1872 = t1801 + t1957;
	t1665 = t1742 + t1750 + t1872;
	t2049 = -0.24e2 * t1665;
	t1771 = 10 * t1801;
	t1676 = t1807 + t1805 / 0.4e1 + t1801 / 0.4e1 - t1794 / 0.8e1;
	t1761 = sin(qJ(2));
	t1995 = t1761 * t1762;
	t1691 = pkin(2) * t1995;
	t1864 = pkin(3) * t1691;
	t1967 = 0.4e1 / 0.7e1 * t1807 - t1794 / 0.7e1;
	t1593 = -0.32e2 / 0.21e2 * t1676 * t1864 + t1813 / 0.7e1 + (0.16e2 / 0.21e2 * t1801 + t1967) * t1805 + t1799 / 0.7e1 + t1967 * t1801 + t1806 - 0.3e1 / 0.7e1 * t1981 + t1961 / 0.42e2;
	t1741 = -t1794 / 0.4e1;
	t1756 = t1801 / 0.3e1;
	t1759 = t1805 / 0.2e1;
	t1968 = t1759 + t1807;
	t1677 = t1741 + t1756 + t1968;
	t1755 = 0.4e1 / 0.3e1 * t1801;
	t2023 = 0.4e1 / 0.3e1 * t1805;
	t1594 = -0.8e1 / 0.3e1 * t1677 * t1864 + t1813 / 0.3e1 + (t1755 + t1742) * t1805 + t1806 - t1799 / 0.3e1 + (t2023 + 0.2e1 / 0.3e1 * t1801 + t1744) * t1807 + t1961 / 0.18e2;
	t1757 = t1801 / 0.2e1;
	t1642 = -0.2e1 / 0.3e1 * t1864 + t1807 + t1757 + t1741;
	t1789 = 4 * t1807;
	t1708 = (t1789 + t1794) * t1801;
	t1711 = -t1801 / 0.3e1 + t1807;
	t1713 = t1807 - 0.2e1 / 0.3e1 * t1805;
	t1763 = cos(qJ(2));
	t1716 = pkin(3) * t1763;
	t1720 = -t1801 + t1807;
	t1743 = -t1794 / 0.2e1;
	t1687 = t1743 + t1872;
	t1849 = -0.4e1 * t1864;
	t1836 = t1687 * t1849;
	t1726 = t1763 ^ 2;
	t1725 = t1763 * t1726;
	t1816 = pkin(3) * t1801;
	t2001 = t1725 * t1816;
	t1909 = 0.16e2 * t2001;
	t1724 = t1726 ^ 2;
	t1985 = t1799 * t1724;
	t1928 = 0.8e1 * t1985;
	t1983 = t1801 * t1726;
	t1988 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
	t1574 = t1713 * t1928 + 0.14e2 * t1593 * t1983 + t1711 * t1836 + t1720 * t1813 + (t1708 - 0.10e2 / 0.3e1 * t1799 + (2 * t1806) - t1981) * t1805 + t1639 * t1988 + (0.6e1 * t1594 * t1716 + t1642 * t1909) * pkin(1);
	t1745 = -0.3e1 / 0.2e1 * t1794;
	t1769 = 15 * t1799;
	t1772 = 18 * t1807;
	t1786 = 3 * t1806;
	t1971 = t1793 / 0.2e1 - t1796 / 0.2e1;
	t1844 = -0.3e1 * t1981 + t1786 + t1971;
	t1790 = 3 * t1807;
	t1966 = 15 * t1801 + t1790;
	t1812 = pkin(2) * t1805;
	t1802 = t1812 ^ 2;
	t1976 = t1719 * ((t1745 + t1791) * t1801 - 0.3e1 / 0.2e1 * t1981 + t1959 + t1971) + t1802;
	t1583 = t1864 * t2050 + (t1769 + (t1772 - 0.9e1 * t1794) * t1801 + t1844) * t1805 + (t1745 + t1966) * t1813 + t1976;
	t1598 = t1836 + (t1781 + t1962) * t1805 + t1825;
	t1684 = -0.2e1 * t1864;
	t1712 = t1807 - t1805 / 0.3e1;
	t1650 = t1712 * t1684;
	t1987 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
	t1612 = t1687 * t1987 + t1650;
	t1722 = -0.3e1 * t1805 + t1807;
	t1923 = pkin(1) * t2001;
	t1869 = 0.8e1 * t1923;
	t1670 = t1722 * t1869;
	t1938 = pkin(1) * t1716;
	t1899 = 0.6e1 * t1938;
	t1902 = 0.12e2 * t1983;
	t1575 = t1598 * t1899 + t1612 * t1902 + t1583 + t1670;
	t1730 = t1805 * t2031;
	t1773 = -0.2e1 * t1794;
	t1787 = 8 * t1807;
	t2032 = 4 * t1799;
	t1634 = -0.4e1 * t1816 * t1691 + t1730 + t2032 + (t1773 + t1787) * t1801;
	t1830 = -t1864 + t1968;
	t1640 = t1741 - t1801 + t1830;
	t1721 = -3 * t1801 + t1807;
	t1936 = 0.8e1 * t2001;
	t1941 = 0.4e1 * t1716;
	t1586 = t1684 * t1988 + t1634 * t1726 + t1687 * t1721 + (t1640 * t1941 + t1936) * pkin(1);
	t1954 = t1807 - t1794;
	t1956 = t1806 - t1799;
	t1587 = t1712 * t1836 - t1802 + (-t1732 - t1954) * t1813 + (t1708 + t1956 + t2042) * t1805 + t1807 * t1639;
	t1774 = -0.5e1 * t1794;
	t1779 = 7 * t1799;
	t1592 = (t1745 + t1790 + (7 * t1801)) * t1813 + (t1779 + (t1774 + (10 * t1807)) * t1801 + t1844) * t1805 + t1976;
	t1809 = pkin(1) * t1807;
	t1695 = -12 * pkin(1) * t1816 + t1809 * t2052;
	t1980 = t1807 * t1801;
	t1710 = -8 * t1799 + 12 * t1980;
	t1858 = pkin(1) * t1909;
	t1604 = t1695 * t1763 + t1710 * t1726 + t1858 + t1928 + t1959 - (6 * t1980);
	t1616 = t1684 * t1987 + t1687 * t1722;
	t1955 = t1806 + t1813;
	t1982 = t1805 * t1807;
	t1672 = 0.16e2 * (t1955 - 0.6e1 * t1982) * t1799;
	t1871 = t1801 + t1954;
	t1689 = t1797 + t1871;
	t1715 = -0.30e2 * t1794 + (60 * t1807);
	t1764 = cos(qJ(1));
	t1729 = t1764 ^ 2;
	t1727 = t1729 ^ 2;
	t1775 = -0.6e1 * t1794;
	t1832 = (6 * t1806) + t1961 - 0.6e1 * t1981;
	t1697 = pkin(1) + t1716;
	t1728 = t1764 * t1729;
	t1999 = t1728 * t1812;
	t1896 = t1697 * t1999;
	t1851 = -0.32e2 * t1896;
	t1911 = pkin(3) * t1995;
	t1714 = t1719 ^ 2;
	t2002 = t1714 * (-t1797 + t1871);
	t2004 = t1697 * t1764;
	t2013 = t1583 * pkin(3);
	t2034 = 0.8e1 * t1763;
	t1561 = t1586 * t1851 + t1672 * t1724 + 0.24e2 * t1587 * t1983 + (t1773 + t1789 + (28 * t1801)) * t1802 + t1689 * t2002 + (0.24e2 * t1574 * t1729 + t1715 * t1799 + t1806 * t1775 + t1832 * t1781 + t1961 * t1791 + (4 * t1809 ^ 2) + (28 * t1816 ^ 2)) * t1805 + 0.8e1 * (-t1575 * t2004 - t1592 * t1911) * pkin(2) + (0.32e2 * t1616 * t2001 + t2013 * t2034) * pkin(1) + (0.16e2 * t1604 * t1727 + t1715 * t1801 + (70 * t1799) + t1813 + t1832) * t1813;
	t1740 = -t1794 / 0.6e1;
	t1970 = t1740 - t1797 / 0.6e1;
	t1878 = t1807 + t1970;
	t1658 = t2023 + t1757 + t1878;
	t1843 = t1759 + t1878;
	t1659 = t1755 + t1843;
	t2014 = pkin(3) * t1761;
	t1605 = -t1658 * t2014 + t1659 * t1717;
	t1661 = t1801 + t1843;
	t1706 = 0.2e1 * t1805 + t1720;
	t1910 = -t2014 / 0.2e1;
	t1614 = t1661 * t1717 + t1706 * t1910;
	t1733 = -0.20e2 / 0.3e1 * t1801;
	t1880 = 0.2e1 / 0.3e1 * t1794 + t1749 + t1789;
	t1792 = -2 * t1807;
	t1881 = 0.4e1 / 0.3e1 * t1794 + 0.4e1 / 0.3e1 * t1797 + t1792;
	t1618 = -t1813 + (t1733 + t1880) * t1805 - (3 * t1799) + t1881 * t1801 + t1806;
	t2044 = t1742 - t1797 / 0.3e1;
	t1877 = t1807 + t2044;
	t1879 = t1794 / 0.3e1 + t1750 + t1791;
	t1845 = -0.8e1 / 0.3e1 * t1985 + t1805 * t1720 - 0.5e1 / 0.3e1 * t1799 + t1879 * t1801 + t1807 * t1877;
	t1925 = 0.4e1 * t1983;
	t1984 = t1799 * t1725;
	t2009 = 4 * pkin(1);
	t1576 = t1605 * t1925 + t1618 * t1910 + t1845 * t1717 + (t1614 * t1716 - t1761 * t1984) * t2009;
	t1780 = 5 * t1799;
	t1788 = 6 * t1807;
	t1964 = t1773 - 0.2e1 * t1797;
	t1873 = t1788 + t1964;
	t1965 = t1771 + t1791;
	t1754 = -0.2e1 / 0.3e1 * t1797;
	t1969 = t1744 + t1754;
	t1606 = t1813 + (t1965 + t1969) * t1805 + t1780 + t1873 * t1801 + t1807 * (t1807 + t1969);
	t1784 = 0.5e1 * t1813;
	t1621 = t1784 + (t1771 + t1873) * t1805 + t1719 * (t1754 + t1874);
	t1588 = t1606 * t1717 - t1621 * t2014;
	t1785 = 0.3e1 * t1805;
	t1842 = t1801 + t1877;
	t1663 = t1785 + t1842;
	t1782 = 3 * t1801;
	t1707 = t1782 + t1957;
	t1664 = t1707 + t2044;
	t1607 = -t1663 * t2014 + t1664 * t1717;
	t1675 = t1757 + t1805 + t1970;
	t1622 = t1675 * t1951 + t1987 * t2014;
	t1916 = t1816 * t1717;
	t1854 = t1725 * t1916;
	t1926 = -0.4e1 * t1983;
	t1577 = t1622 * t1926 + (t1607 * t1941 - 0.8e1 * t1854) * pkin(1) + t1588;
	t1619 = -0.3e1 * t1813 + (t1733 + t1881) * t1805 + t1880 * t1801 + t1956;
	t1623 = -0.5e1 / 0.3e1 * t1813 + (-t1801 + t1879) * t1805 + t1807 * t1842;
	t1943 = -0.2e1 * t2014;
	t1589 = t1619 * t1717 + t1623 * t1943;
	t1876 = t1743 - t1797 / 0.2e1 + t1807;
	t1660 = 0.3e1 / 0.2e1 * t1805 + t1782 + t1876;
	t1917 = t1801 * t1717;
	t1831 = -t1761 * t1816 + t1917;
	t1667 = 0.4e1 * t1831;
	t1679 = t1717 + 0.2e1 * t2014;
	t2020 = pkin(1) * t1763;
	t1590 = t1721 * t1717 + t1667 * t1726 + (t1660 * t1761 + t1679 * t2020) * t2008;
	t1662 = t1785 + 0.3e1 / 0.2e1 * t1801 + t1876;
	t1620 = t1662 * t1717 + t1722 * t2014 / 0.2e1;
	t2043 = t1790 - t1794 - t1797;
	t1688 = t2043 * t1771;
	t1963 = t1774 - 0.5e1 * t1797;
	t1837 = 0.24e2 * t1712 * t1985 - t1802 - ((21 * t1801) + t2043) * t1813 - (t1964 * t1807 + t1688 + t1786 + (35 * t1799)) * t1805 - (t1779 + (t1787 + t1963) * t1801 + t1807 * (-t1797 + t1954)) * t1719;
	t1856 = 0.8e1 * t1896;
	t1997 = t1729 * t1805;
	t1907 = -0.12e2 * t1997;
	t1927 = -0.6e1 * t1983;
	t1699 = 0.2e1 * t1938;
	t1643 = 0.4e1 / 0.3e1 * t1983 + t1699 + t1711;
	t2000 = t1727 * t1813;
	t2046 = 0.7e1 * t1802 + ((35 * t1801) + (15 * t1807) + t1963) * t1813 + ((21 * t1799) + t1688 + (9 * t1806) + (t1775 - 0.6e1 * t1797) * t1807) * t1805 + t2002 - 0.24e2 * t1643 * t2000;
	t1565 = t1620 * t1858 + t1590 * t1856 + t1576 * t1907 + t1589 * t1927 + (-0.6e1 * t1588 * t2020 + t2046 * t1761) * pkin(3) + (0.6e1 * t1577 * t2004 + t1837 * t1762) * pkin(2);
	t1979 = t1816 * t1726;
	t1901 = -0.24e2 * t1979;
	t1850 = t1761 * t1901;
	t2017 = pkin(2) * t1764;
	t1641 = t1712 * t1850 * t2017;
	t1992 = t1763 * t1764;
	t1915 = pkin(2) * t1992;
	t1865 = pkin(1) * t1915;
	t1833 = t1801 * t1761 * t1865;
	t1657 = -0.4e1 * t1833;
	t1696 = pkin(1) - t2017;
	t1666 = t1696 + t1716;
	t1674 = t1865 * t2008;
	t1870 = t1805 + t1954;
	t1841 = t1801 + t1870;
	t1682 = -t1797 + t1841;
	t1939 = pkin(1) * t2017;
	t1701 = -0.2e1 * t1939;
	t1645 = t1701 + t1682;
	t1625 = t1684 + t1645;
	t1680 = t1696 * t1763;
	t1702 = t1730 - 0.4e1 * t1982;
	t1776 = 0.2e1 * t1797;
	t1932 = 0.2e1 * t1997;
	t1958 = -t1805 + t1807;
	t1648 = t1701 + t1932 + t1958;
	t2005 = t1648 * t1726;
	t2027 = -pkin(4) + pkin(5);
	t2028 = -pkin(4) - pkin(5);
	t1808 = sqrt(t1702 * t1729 + 0.4e1 * t1682 * t1939 - t1799 - (t1807 + (pkin(2) - t2027) * (pkin(2) + t2027)) * (t1807 + (pkin(2) - t2028) * (pkin(2) + t2028)) + (t1776 + t1792 + 0.2e1 * t1794 - 0.6e1 * t1805 - 0.4e1 * t2005) * t1801 + (-t1625 * t1680 + t1645 * t1691) * t2052);
	t2019 = pkin(1) * t1799;
	t1838 = -0.64e2 * t1987 * t2019;
	t1898 = 0.32e2 / 0.3e1 * t1799;
	t1847 = t1725 * t1898;
	t1848 = 0.64e2 / 0.3e1 * t1676 * t1816;
	t2030 = -0.96e2 * t1712;
	t1852 = t1687 * t1816 * t2030;
	t2018 = pkin(1) * t1801;
	t1859 = -0.24e2 * t1687 * t2018;
	t1860 = -0.16e2 * t1677 * t2018;
	t1861 = -0.48e2 * t1615 * t2018;
	t2012 = t1711 * pkin(3);
	t1866 = -0.4e1 * t1687 * t2012;
	t2038 = -2 * pkin(3);
	t1867 = t1988 * t2038;
	t1700 = pkin(1) * t1951;
	t1989 = t1764 * t1805;
	t1894 = t1762 * t1989;
	t1652 = t1700 - 0.4e1 * t1894;
	t2047 = t1761 * t1805;
	t1920 = pkin(1) * t2047;
	t1994 = t1761 * t1764;
	t1919 = pkin(2) * t1994;
	t1863 = pkin(3) * t1919;
	t1973 = -0.2e1 * t1863 + t1700;
	t1993 = t1763 * t1762;
	t2033 = -0.2e1 * t1764;
	t2040 = -0.4e1 * pkin(2);
	t1578 = t1652 * t1926 + (pkin(1) * t1682 * t2040 + t1702 * t2033) * t1762 + (-t1973 * t1680 + 0.2e1 * t1762 ^ 2 * t1920 + (-t1625 * t1993 + t1645 * t1994) * pkin(2)) * t2052;
	t1580 = 0.1e1 / t1808;
	t2024 = t1580 / 0.2e1;
	t1885 = t1578 * t2024;
	t1996 = t1729 * t1812;
	t1895 = t1697 * t1996;
	t1900 = 0.4e1 * t1938;
	t1905 = 0.24e2 * t1997;
	t1908 = -0.32e2 * t1999;
	t1934 = 0.8e1 * t1999;
	t1940 = -0.8e1 * t1592 * pkin(3);
	t1946 = pkin(3) * t2050;
	t1947 = 0.6e1 * t2017;
	t1950 = -0.8e1 * t2017;
	t1998 = t1728 * t1813;
	t2036 = -0.8e1 * t1666;
	t2037 = 0.8e1 * t1575;
	t2045 = t1721 + t1925;
	t1553 = t1565 * t1885 + (-0.32e2 * t1657 * t1666 + 0.8e1 * t1674 * t1808) * t1896 + (0.24e2 * (0.4e1 * t1643 * t1998 * t2014 + t1576 * t1989 - t1590 * t1895) * t1808 + t1666 * (-0.48e2 * t1574 * t1989 + 0.96e2 * t1586 * t1895 - 0.64e2 * t1604 * t1998)) * t1762 + ((t1561 + (-0.6e1 * t1577 * t1808 + t1666 * t2037) * t1697) * t1762 + (t1697 * t1641 * t2036 + ((t1659 * t1925 + t1661 * t1900 + t1845) * t1907 + t1662 * t1858 + t1619 * t1927 - 0.6e1 * t1606 * t1938 + (t2045 * t1934 + (t1664 * t1900 - 0.8e1 * t1675 * t1983 + t1606 - 0.8e1 * t1923) * t1947) * t1697 + t1837) * t1808 + ((-pkin(1) * t1847 - t1726 * t1848 + t1763 * t1860 + t1866) * t1905 + t1725 * t1838 + t1726 * t1852 + t1763 * t1861 + t1940 + ((t1867 - 0.4e1 * t1979) * t1908 + (t1763 * t1859 + t1946) * t1950) * t1697) * t1666 * t1761) * t1764) * pkin(2);
	t1559 = t1561 * t1666 + t1565 * t1808;
	t1629 = -0.4e1 / 0.9e1 * t1864 + t1807 + t1805 / 0.3e1 + t1756 + t1797 / 0.9e1 - t1794 / 0.9e1;
	t1637 = t1797 / 0.6e1 + t1740 + t1830;
	t1584 = t1711 * t1684 + 0.6e1 * t1629 * t1983 + t1665 * t1988 + (t1637 * t1716 + t2001) * t2009;
	t1595 = t1665 * t1849 + (t1781 + t1875) * t1805 + t1974;
	t1635 = t1684 + t1665;
	t1703 = t1958 * t2031;
	t1585 = t1635 * t1900 + t1703 * t1726 + t1595;
	t1610 = t1665 * t1987 + t1650;
	t1647 = t1699 + t2045;
	t1570 = -0.8e1 * t1647 * t1896 + t1670 + t1610 * t1902 + t1595 * t1899 + t1802 + (-t1794 + t1797 + t1966) * t1813 + t1714 * t1689 + (0.12e2 * t1584 * t1729 + t1769 + (t1772 + t1775 + 0.6e1 * t1797) * t1801 + t1786 + (t1773 + t1776) * t1807) * t1805 + 0.6e1 * (-t1585 * t2004 - t1613 * t1911) * pkin(2);
	t1671 = t1707 * t1717;
	t1678 = t1717 - t2014;
	t1855 = t1726 * t1917;
	t1705 = t1785 + t1719;
	t2003 = t1705 * t1761;
	t1596 = -0.2e1 * t1855 + t1671 + (0.2e1 * t1678 * t2020 - t2003) * pkin(3);
	t1668 = 0.2e1 * t1831;
	t1718 = (2 * t1801) + t1805;
	t1599 = t1718 * t2014 + t1668 * t1726 + (t1720 + t1699) * t1717;
	t1628 = -pkin(3) * t2003 + t1671;
	t1704 = pkin(2) * t2032 + 0.8e1 * t1801 * t1812;
	t1893 = t1816 * t1987;
	t1630 = t1704 * t1762 + 0.4e1 * t1761 * t1893;
	t1646 = t1965 * t1805 + t1780 + t1955 + (6 * t1980);
	t1655 = t1784 + (t1771 + t1788) * t1805 + t1714;
	t1933 = -0.4e1 * t1997;
	t1573 = t1599 * t1933 + t1630 * t1726 + (-0.4e1 * t1628 * t2020 + (t1655 + t1856) * t1761) * pkin(3) + (0.4e1 * t1596 * t2004 + (-t1646 + t1869) * t1762) * pkin(2);
	t1906 = 0.12e2 * t1997;
	t1991 = t1763 * t1801;
	t1930 = -0.8e1 * t1991;
	t1931 = 0.8e1 * t1762 * t1805;
	t2039 = -4 * pkin(3);
	t1945 = t1665 * t2039;
	t1948 = 0.4e1 * t2017;
	t1952 = t1596 * t2040;
	t1560 = (t1674 * t1933 + (t1599 * t1931 + t1704 * t1726 + ((t1720 + 0.2e1 * t1983) * t1933 - t1646 + (-0.4e1 * t1707 * t1716 + t1936) * pkin(1)) * pkin(2)) * t1764 + ((t1674 + (t1707 - 0.2e1 * t1983) * t2017) * t1948 + (-0.24e2 * t1996 * t2014 + t1952) * t1762) * t1697) * t1808 + t1573 * t1885 + t1570 * t1717 + t1666 * (0.24e2 * t1647 * t1762 * t1895 + (t1657 + (-0.8e1 / 0.3e1 * t1979 - 0.2e1 * t2012) * t1919) * t1906 - 0.24e2 * t1584 * t1894 + t1641 + t1833 * t2049 + t1863 * t2051 + 0.6e1 * (-(pkin(1) * t1930 + t1945) * t1729 * t2047 + t1585 * t1717) * t1697);
	t1644 = t1701 + t1797 + t1841;
	t2016 = pkin(3) * t1696;
	t1601 = t1644 * t1763 + (0.4e1 * t1726 - 0.2e1) * t2016;
	t1972 = -t1691 + t1680;
	t1636 = pkin(3) + t1972;
	t1883 = t1636 * t2024;
	t1942 = 0.2e1 * t1716;
	t2015 = pkin(3) * t1726;
	t1566 = t1578 * t1883 + t1652 * t1761 * t1942 + ((-t1761 * t1808 + t1601) * t1764 + (t1763 * t1808 + (t1696 * t2041 + t1644) * t1761 + (-pkin(3) + 0.2e1 * t2015 + t2020) * t1951) * t1762) * pkin(2);
	t1681 = t1782 + t1797 + t1870;
	t1626 = t1681 + t1701 + t1849;
	t1683 = t1691 - pkin(3);
	t1918 = pkin(2) * t1993;
	t1638 = t1696 * t1761 + t1918;
	t1882 = t1638 * t2024;
	t1567 = (t1691 + t1915) * t1808 + t1578 * t1882 - 0.2e1 * t1652 * t2015 - (t1700 - 0.4e1 * t1863) * t1680 - 0.2e1 * t1729 * t1920 + t1681 * t1919 + (t1989 * t2039 + (-t1626 * t1763 + t1683 * t2041) * pkin(2)) * t1762;
	t1913 = pkin(3) * t1680;
	t1611 = t1684 + t1701 + t1872 + 0.2e1 * t1913;
	t1608 = 0.1e1 / t1611;
	t1862 = pkin(3) * t1918;
	t1627 = 0.2e1 * t1862 + t1973;
	t1609 = 0.1e1 / t1611 ^ 2;
	t1953 = pkin(1) * t2033;
	t2006 = t1638 * t1808;
	t1571 = -t1626 * t1680 + t2006 + (t1681 * t1995 + t1683 * t1953) * pkin(2) + (-t1689 - t1785 + t1932 - 0.2e1 * t2005) * pkin(3);
	t1914 = t1648 * t1716;
	t1600 = t1644 * t1696 + 0.2e1 * t1914;
	t1572 = t1600 * t1761 + t1601 * t1717 + t1636 * t1808;
	t1564 = t1570 * t1666 + t1573 * t1808;
	t2026 = 0.1e1 / t1564 / 0.4e1;
	t1890 = t1572 * t2026;
	t2021 = t1808 / 0.4e1;
	t1829 = t1559 * t1890 + t1571 * t2021;
	t1827 = t1609 * t1829;
	t1887 = t1571 * t1580 / 0.8e1;
	t2025 = -0.1e1 / t1564 ^ 2 / 0.4e1;
	t1888 = t1572 * t2025;
	t1986 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
	t1543 = (-t1627 * t1827 + (t1567 * t2021 + t1578 * t1887 + t1553 * t1890 + (t1560 * t1888 + t1566 * t2026) * t1559) * t1608) * t1986;
	t1891 = t1571 * t2026;
	t2022 = -t1808 / 0.4e1;
	t1828 = t1559 * t1891 + t1572 * t2022;
	t1897 = t1608 * t1986;
	t1556 = t1828 * t1897;
	t1653 = -t1992 - t1995;
	t1555 = t1556 * t1653;
	t1557 = t1829 * t1897;
	t1654 = -t1993 + t1994;
	t1551 = -t1557 * t1654 - t1555;
	t1549 = 0.1e1 / t1551 ^ 2;
	t1554 = t1556 * t1654;
	t1550 = -t1557 * t1653 + t1554;
	t1546 = 0.1e1 / (t1549 * t1550 ^ 2 + 0.1e1);
	t1548 = 0.1e1 / t1551;
	t1826 = t1609 * t1828;
	t1886 = -t1572 * t1580 / 0.8e1;
	t1889 = t1571 * t2025;
	t1977 = (-t1627 * t1826 + (t1553 * t1891 + t1566 * t2022 + t1578 * t1886 + (t1560 * t1889 + t1567 * t2026) * t1559) * t1608) * t1986 + t1557;
	t2007 = t1549 * t1550;
	t2048 = ((t1977 * t1654 + (-t1543 + t1556) * t1653) * t1548 - (-t1543 * t1654 - t1977 * t1653 + t1554) * t2007) * t1546 + 0.1e1;
	t2035 = -0.2e1 * t1763;
	t2029 = pkin(1) * pkin(3);
	t2011 = -6 * t2029;
	t2010 = -4 * t2029;
	t1990 = t1763 * t1805;
	t1839 = pkin(1) * t1855;
	t1673 = -0.4e1 * t1839;
	t1698 = pkin(1) * t1943;
	t1723 = t1761 ^ 2;
	t1840 = pkin(1) * t1726 * t1916;
	t1846 = t1725 * t1893;
	t1921 = pkin(1) * t1979;
	t1857 = -0.48e2 * t1921;
	t1922 = pkin(1) * t1983;
	t1868 = 0.4e1 * t1922;
	t1929 = 0.8e1 * t1991;
	t1944 = 0.4e1 * t2016;
	t1582 = (t1625 * t1944 + t1648 * t1929) * t1761 + (t1645 * t1941 + 0.8e1 * t1696 * t1983) * t1717;
	t1884 = t1582 * t2024;
	t1892 = t1761 * t1991;
	t1903 = -0.32e2 * t1984;
	t1904 = -0.24e2 * t1991;
	t1912 = pkin(3) * t1999;
	t1924 = pkin(1) * t1985;
	t1935 = -0.8e1 * t1999;
	t1937 = -0.4e1 * t2001;
	t1949 = -0.6e1 * t2017;
	t1975 = pkin(1) * t1722 * t1850 - 0.24e2 * t1712 * t1854;
	t1552 = ((t1660 * t1942 + t1868 + t1937) * t1856 + (0.12e2 * t1723 * t1726 * t2019 + t1658 * t1937 - 0.2e1 * t1706 * t1922 - 0.4e1 * t1924) * t1907 + 0.8e1 * t1722 * t1924 + 0.12e2 * t1623 * t2001 + 0.6e1 * t1621 * t1922 + (-t1618 * t1907 / 0.2e1 + t2046) * t1716) * t1808 + t1565 * t1884 + t1666 * t1673 * t1851 + ((0.6e1 * (-t1621 * t1716 - 0.4e1 * t1663 * t1922 - 0.4e1 * t1846) * t1808 + t1975 * t2036) * t2004 + ((-pkin(1) * t1724 * t1898 - t1725 * t1848 + t1726 * t1860 + t1763 * t1866) * t1905 + t1724 * t1838 + t1725 * t1852 + t1726 * t1861 + t1763 * t1940 + ((t1763 * t1867 + t1937) * t1908 + (t1726 * t1859 + t1763 * t1946) * t1950) * t1697) * t1666 * t1762) * pkin(2) + (-pkin(3) * t1561 + t1666 * (0.16e2 * (t1710 * t2035 - t1695 + t1857 + t1903) * t2000 + 0.32e2 * t1586 * t1912 + (pkin(1) * t1901 + t1634 * t2035 + t1640 * t2010) * t1851 + (-0.28e2 * t1593 * t1991 + t1594 * t2011 + t1642 * t1857 + t1713 * t1903) * t1905 + pkin(3) * t2017 * t2037 + t1697 * (t1598 * t2011 + t1612 * t1904) * t1950 - 0.4e1 * t1672 * t1725 - 0.96e2 * t1616 * t1921 - 0.48e2 * t1587 * t1991 - 0.8e1 * pkin(1) * t2013) + ((t1605 * t1930 + t1847 * t1717) * t1907 + t1984 * t1717 * t2030 + t1620 * t1857 + 0.12e2 * t1589 * t1991 + (t1577 * t1949 + t1590 * t1935 + (-0.24e2 * t1698 + 0.64e2 * t1892) * t2000 + (0.48e2 * t1614 * t1997 + 0.6e1 * t1588) * pkin(1)) * pkin(3) + (0.2e1 * (-t1667 * t1763 - t1679 * t2029) * t1934 + (t1607 * t2010 + t1622 * t1929 + 0.24e2 * t1840) * t1947) * t1697) * t1808) * t1761;
	t1558 = (t1801 * t1723 * t1935 + (t1718 * t1716 - 0.2e1 * t2001) * t1933 + 0.4e1 * t1846 + t1705 * t1868 + t1655 * t1716) * t1808 + t1573 * t1884 + t1666 * ((-0.8e1 / 0.3e1 * t1854 + t1673 - 0.2e1 * t1711 * t1862) * t1906 + t1839 * t2049 + t1862 * t2051 + t1975) + ((0.8e1 * t1668 * t1729 * t1990 + t1630 * t2035 - 0.24e2 * t1840) * t1808 + t1666 * (0.12e2 * (-t1629 * t1991 - t1921) * t1906 + t1610 * t1904) + (t1764 * t1808 * t1952 - t1570 + t1666 * (t1585 * t1947 + t1647 * t1934) + ((pkin(2) * t1729 * t1931 + 0.4e1 * t1628) * t1808 + t1666 * (-0.48e2 * t1637 * t1997 - 0.6e1 * t1595)) * pkin(1)) * pkin(3)) * t1761 + ((t1912 * t2034 + ((-pkin(3) * t1705 + t1691 * t2031) * t1763 + (-t1678 * t2014 - t1983) * t2041) * t1948) * t1808 + t1666 * ((t1698 - 0.8e1 * t1892) * t1935 + ((-0.2e1 * t1703 * t1761 + t1945 * t1717) * t1763 + (-0.4e1 * t1635 * t2014 - 0.8e1 * t1855) * pkin(1)) * t1949)) * t1697;
	t1568 = -t2006 + t1582 * t1883 + t1648 * t1723 * t2038 + t1600 * t1763 + (-t1644 - 0.8e1 * t1913) * t1691;
	t1569 = t1972 * t1808 + t1582 * t1882 + (t1626 * t1696 + 0.4e1 * t1914) * t1761 + (t1953 * t1990 + (t1681 * t1763 + t1726 * t1944) * pkin(2)) * t1762;
	t1631 = t1638 * t2008;
	t1978 = (t1631 * t1826 + (t1552 * t1891 + t1568 * t2022 + t1582 * t1886 + (t1558 * t1889 + t1569 * t2026) * t1559) * t1608) * t1986 - t1557;
	t1547 = atan2(t1550, t1551);
	t1544 = sin(t1547);
	t1545 = cos(t1547);
	t1835 = t1544 * t1764 + t1545 * t1762;
	t1834 = t1544 * t1762 - t1545 * t1764;
	t1541 = (t1631 * t1827 + (t1569 * t2021 + t1582 * t1887 + t1552 * t1890 + (t1558 * t1888 + t1568 * t2026) * t1559) * t1608) * t1986;
	t1538 = ((-t1541 * t1653 + t1978 * t1654 - t1555) * t1548 - ((-t1541 - t1556) * t1654 - t1978 * t1653) * t2007) * t1546;
	t1537 = t1835 * t1538;
	t1536 = t1834 * t1538;
	t1535 = t2048 * t1835;
	t1534 = t2048 * t1834;
	t1 = [t1535, t1537; t1534, t1536; 0, 0; -t1534, -t1536; t1535, t1537; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:51:59
	% EndTime: 2020-04-27 03:51:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t1 = [0, -t8; 0, t9; 0, 0; 0, -t9; 0, -t8; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 03:52:03
	% EndTime: 2020-04-27 03:52:04
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (17744->109), mult. (38512->194), div. (320->7), fcn. (9512->9), ass. (0->94)
	t368 = (pkin(4) ^ 2);
	t370 = (pkin(3) ^ 2);
	t371 = pkin(2) ^ 2;
	t372 = pkin(1) ^ 2;
	t395 = -pkin(5) ^ 2 + t372;
	t387 = (t371 + t395);
	t350 = 3 * t370 + t368 + t387;
	t365 = cos(qJ(1));
	t410 = pkin(2) * t365;
	t393 = pkin(1) * t410;
	t357 = -0.2e1 * t393;
	t362 = sin(qJ(2));
	t363 = sin(qJ(1));
	t401 = t362 * t363;
	t354 = pkin(2) * t401;
	t383 = pkin(3) * t354;
	t340 = t350 + t357 - 0.4e1 * t383;
	t355 = pkin(1) - t410;
	t364 = cos(qJ(2));
	t349 = t355 * t364;
	t360 = t364 ^ 2;
	t396 = t357 + t372;
	t361 = t365 ^ 2;
	t416 = 0.2e1 * t361;
	t347 = t371 * t416 - t371 + t396;
	t418 = -0.2e1 * t347;
	t379 = t360 * t418 - t395;
	t399 = t363 * t364;
	t344 = pkin(2) * t399 + t355 * t362;
	t380 = t370 + t387;
	t351 = -t368 + t380;
	t346 = t357 + t351;
	t353 = -0.2e1 * t383;
	t339 = t353 + t346;
	t358 = (t370 - t372) * t371;
	t414 = -pkin(4) + pkin(5);
	t415 = -pkin(4) - pkin(5);
	t373 = sqrt(0.4e1 * t358 * t361 + 0.4e1 * t351 * t393 - (t372 + (pkin(2) - t414) * (pkin(2) + t414)) * (t372 + (pkin(2) - t415) * (pkin(2) + t415)) + 0.4e1 * (-t339 * t349 + t346 * t354) * pkin(3) + ((2 * t368) - 0.6e1 * t371 + 0.2e1 * t379 - t370) * t370);
	t404 = t344 * t373;
	t412 = pkin(1) * (t354 - pkin(3));
	t327 = -t340 * t349 + t404 + (t350 * t401 - 0.2e1 * t365 * t412) * pkin(2) + (-t368 - t370 + (t416 - 0.3e1) * t371 + t379) * pkin(3);
	t345 = t357 + t368 + t380;
	t403 = t347 * t364;
	t390 = pkin(3) * t403;
	t334 = t345 * t355 + 0.2e1 * t390;
	t409 = pkin(3) * t355;
	t335 = t345 * t364 + (0.4e1 * t360 - 0.2e1) * t409;
	t397 = -t354 + t349;
	t343 = pkin(3) + t397;
	t411 = pkin(2) * t363;
	t328 = t334 * t362 + t335 * t411 + t343 * t373;
	t326 = 0.1e1 / t327 ^ 2;
	t389 = pkin(3) * t349;
	t338 = t353 + t370 + t371 + 0.2e1 * t389 + t396;
	t369 = 0.1e1 / pkin(4);
	t381 = pkin(4) / (t326 * t328 ^ 2 + 0.1e1) * t338 * t369;
	t375 = t326 * t328 * t381;
	t378 = 0.1e1 / t327 * t381;
	t384 = -0.8e1 * t389;
	t331 = 0.1e1 / t373;
	t402 = t360 * t370;
	t407 = pkin(3) * t364;
	t417 = 0.2e1 * t355;
	t392 = 0.2e1 * t331 * ((t339 * t409 + 0.2e1 * t370 * t403) * t362 + (t346 * t407 + t402 * t417) * t411);
	t394 = -0.2e1 * pkin(1) * t371;
	t398 = t364 * t365;
	t337 = 0.1e1 / t338 ^ 2;
	t405 = t337 * t344 * pkin(3);
	t408 = pkin(3) * t360;
	t336 = 0.1e1 / t338;
	t413 = t336 / 0.2e1;
	t420 = 0.2e1 * ((-t404 + t343 * t392 + t334 * t364 + ((-t345 + t384) * t411 + pkin(3) * t362 * t418) * t362) * t413 + t328 * t405) * t378 - 0.2e1 * ((t397 * t373 + t344 * t392 + (t340 * t355 + 0.4e1 * t390) * t362 + (t394 * t398 + (t350 * t364 + 0.4e1 * t355 * t408) * pkin(2)) * t363) * t413 + t327 * t405) * t375 + 0.1e1;
	t419 = -0.4e1 * pkin(2);
	t359 = pkin(1) * t411;
	t400 = t362 * t365;
	t406 = t337 * (t359 + (t399 - t400) * pkin(3) * pkin(2));
	t391 = -0.4e1 * t365 * t371;
	t388 = pkin(2) * t400;
	t356 = 0.2e1 * t359;
	t348 = t363 * t391 + t356;
	t382 = pkin(3) * t388;
	t386 = (-0.4e1 * t348 * t402 + (t359 - t382) * t384 + 0.4e1 * t346 * t382 + (t339 * t407 * t419 - 0.8e1 * t358 * t365 + (0.8e1 * pkin(3) * t371 * t401 + t351 * t419) * pkin(1)) * t363) * t331 / 0.2e1;
	t385 = t369 * t413;
	t323 = atan2(t328 * t385, t327 * t385);
	t321 = sin(t323);
	t322 = cos(t323);
	t377 = -t364 * t321 - t362 * t322;
	t376 = t362 * t321 - t364 * t322;
	t319 = 0.2e1 * ((0.2e1 * t348 * t362 * t407 + t343 * t386) * t413 - t328 * t406 + ((-t362 * t373 + t335) * t365 * t413 + (t364 * t373 / 0.2e1 + (pkin(1) * t417 + t345) * t362 / 0.2e1 + (pkin(1) * t364 - pkin(3) + 0.2e1 * t408) * t411) * t336 * t363) * pkin(2)) * t378 - 0.2e1 * (((pkin(2) * t398 + t354) * t373 + t344 * t386 - 0.2e1 * t348 * t408 - (t356 - 0.4e1 * t382) * t349 + t361 * t362 * t394 + t350 * t388 + (pkin(3) * t391 + (-t340 * t364 + 0.2e1 * t412) * pkin(2)) * t363) * t413 - t327 * t406) * t375;
	t318 = t377 * t319;
	t317 = t376 * t319;
	t316 = t420 * t377;
	t315 = t420 * t376;
	t1 = [t318, t316; -t317, -t315; 0, 0; t317, t315; t318, t316; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
end