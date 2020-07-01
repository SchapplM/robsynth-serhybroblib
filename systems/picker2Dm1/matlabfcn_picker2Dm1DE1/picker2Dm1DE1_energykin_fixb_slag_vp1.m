% Calculate kinetic energy for
% picker2Dm1DE1
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
% Datum: 2020-05-10 19:54
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm1DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm1DE1_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 09:15:05
% EndTime: 2020-05-10 09:16:01
% DurationCPUTime: 52.50s
% Computational Cost: add. (705908->768), mult. (1964211->1319), div. (25120->42), fcn. (527352->38), ass. (0->588)
t2236 = 4 * pkin(1);
t1975 = (pkin(3) ^ 2);
t2291 = -2 * t1975;
t1984 = (pkin(7) ^ 2);
t1959 = -2 * t1984;
t1926 = sin(qJ(2));
t1927 = sin(qJ(1));
t2191 = t1926 * t1927;
t2256 = pkin(1) * pkin(3);
t2074 = t2191 * t2256;
t2045 = qJD(2) * t2074;
t1929 = cos(qJ(2));
t1930 = cos(qJ(1));
t2187 = t1929 * t1930;
t2105 = qJD(1) * t2187;
t2046 = t2105 * t2256;
t2290 = -0.6e1 * t2046 + 0.6e1 * t2045;
t2289 = -t2046 + t2045;
t1963 = 2 * pkin(2);
t2262 = 0.2e1 * t1926;
t2176 = pkin(7) * t2262;
t1862 = pkin(3) * t2176;
t2153 = -t1975 + t1984;
t1890 = t1926 ^ 2;
t2196 = t1890 * t1975;
t1813 = t1862 + t2153 + 0.2e1 * t2196;
t1981 = pkin(1) ^ 2;
t1970 = (pkin(4) ^ 2);
t2149 = t1984 - t1970;
t1868 = t1981 + t2149;
t1822 = t1862 + t1868;
t1944 = 2 * t1975;
t1821 = t1944 + t1822;
t1879 = pkin(3) * t1926;
t1860 = t1879 + pkin(7);
t1881 = pkin(1) * t1930;
t2147 = 0.2e1 * t1881;
t1747 = t1813 * t2147 + t1821 * t1860;
t1894 = t1930 ^ 2;
t1751 = t1821 * t1930 + (0.4e1 * t1894 - 0.2e1) * t1860 * pkin(1);
t2189 = t1927 * t1929;
t2129 = pkin(3) * t2189;
t1843 = -pkin(1) + t2129;
t2202 = t1860 * t1930;
t1794 = -t1843 + t2202;
t1880 = pkin(3) * t1929;
t2073 = pkin(1) * t2129;
t1844 = -0.2e1 * t2073;
t1788 = t1844 + t1822;
t2279 = 4 * t1981;
t1895 = t1975 * t2279;
t2180 = t1975 * t1984;
t1866 = t1895 - 4 * t2180;
t1979 = t1981 ^ 2;
t2267 = -0.4e1 * pkin(7);
t2119 = pkin(3) * t1868 * t2267;
t2204 = t1813 * t1894;
t2282 = 2 * pkin(3);
t1719 = t1866 * t1890 + t1926 * t2119 - t1979 - ((t1984 - (t2282 + pkin(4)) * pkin(4)) * (t1984 + (t2282 - pkin(4)) * pkin(4))) + (t1959 + (2 * t1970) - (4 * t1975) - 0.4e1 * t2204) * t1981 + (-t1788 * t2202 + t1822 * t2129) * t2236;
t1985 = sqrt(t1719);
t1709 = t1747 * t1927 + t1751 * t1880 + t1794 * t1985;
t1923 = sin(pkin(8));
t1924 = cos(pkin(8));
t1812 = -t1923 * t1927 - t1924 * t1930;
t2243 = pkin(5) * t1812;
t1799 = -0.2e1 * pkin(1) * t2243;
t1967 = pkin(5) ^ 2;
t1790 = t1799 + 0.2e1 * t1967;
t1797 = -pkin(1) + t2243;
t1811 = t1923 * t1930 - t1924 * t1927;
t2168 = t1799 + t1967;
t2283 = 2 * pkin(1);
t1776 = -(t2283 + pkin(5)) * pkin(5) + t2168;
t1777 = pkin(5) * (t2283 - pkin(5)) + t2168;
t1987 = sqrt(-t1776 * t1777);
t1724 = pkin(5) * t1790 * t1811 - t1797 * t1987;
t2216 = t1709 * t1724;
t1950 = 3 * t1981;
t1856 = t1944 + t1950 + t2149;
t2062 = -0.4e1 * t2073;
t1783 = t1856 + t1862 + t2062;
t2128 = pkin(3) * t2187;
t2203 = t1860 * t1927;
t1796 = t2128 + t2203;
t1707 = -t1783 * t2202 + t1796 * t1985 + (t1843 * t2176 + t1856 * t2189) * pkin(3) + (-0.2e1 * t2204 + (0.2e1 * t1890 - 0.4e1) * t1975 - t1868) * pkin(1);
t2205 = t1811 * t1987;
t1722 = -pkin(5) * t2205 - t1790 * t1797;
t2219 = t1707 * t1722;
t2021 = t2219 / 0.4e1 + t2216 / 0.4e1;
t1784 = t1981 + t2168;
t1780 = 0.1e1 / t1784;
t2210 = t1780 / pkin(5);
t1885 = t1981 + t1984;
t2081 = t1975 + t1885;
t2133 = pkin(1) * t2202;
t1760 = t1844 + t1862 + t2081 + 0.2e1 * t2133;
t1758 = 0.1e1 / t1760;
t1976 = 1 / pkin(3);
t2212 = t1758 * t1976;
t2057 = t2210 * t2212;
t1681 = t2021 * t2057;
t2215 = t1722 * t1709;
t2218 = t1707 * t1724;
t2020 = t2218 / 0.4e1 - t2215 / 0.4e1;
t1682 = t2020 * t2057;
t1928 = sin(pkin(9));
t1931 = cos(pkin(9));
t1673 = t1681 * t1931 - t1682 * t1928;
t2241 = pkin(6) * t1673;
t1670 = t2241 * t1963;
t1965 = pkin(6) ^ 2;
t2259 = 0.2e1 * t1965;
t1668 = t1670 + t2259;
t1669 = -pkin(2) - t2241;
t2173 = t1670 + t1965;
t1663 = -(t1963 + pkin(6)) * pkin(6) + t2173;
t1664 = pkin(6) * (t1963 - pkin(6)) + t2173;
t2220 = t1663 * t1664;
t1986 = sqrt(-t2220);
t1675 = t1681 * t1928 + t1682 * t1931;
t2240 = pkin(6) * t1675;
t1618 = -t1668 * t1669 - t1986 * t2240;
t1619 = t1668 * t2240 - t1669 * t1986;
t1991 = pkin(2) ^ 2;
t1667 = t1991 + t2173;
t1966 = 0.1e1 / pkin(6);
t2102 = 0.1e1 / t1667 * t1966 / 0.2e1;
t1608 = atan2(t1619 * t2102, t1618 * t2102);
t1606 = sin(t1608);
t1607 = cos(t1608);
t2096 = t2212 / 0.2e1;
t2009 = atan2(t1709 * t2096, t1707 * t2096);
t2006 = sin(t2009);
t2007 = cos(t2009);
t1678 = -t1927 * t2007 - t1930 * t2006;
t1679 = t1927 * t2006 - t1930 * t2007;
t1599 = t1606 * t1679 + t1607 * t1678;
t2183 = t1966 / pkin(2);
t2277 = t1986 / 0.2e1;
t1635 = atan2(t2183 * t2277, -t1673);
t1633 = sin(t1635);
t1634 = cos(t1635);
t2033 = t1606 * t1678 - t1607 * t1679;
t2288 = -t1599 * t1633 - t1634 * t2033;
t2287 = t1599 * t1634 - t1633 * t2033;
t1840 = t1984 + t1975 / 0.4e1 + t1981 / 0.4e1 - t1970 / 0.8e1;
t1969 = t1970 ^ 2;
t1983 = t1984 ^ 2;
t1993 = t1975 ^ 2;
t2162 = 0.4e1 / 0.7e1 * t1984 - t1970 / 0.7e1;
t2178 = t1984 * t1970;
t1739 = -0.32e2 / 0.21e2 * t1840 * t2073 + 0.5e1 / 0.42e2 * t1993 + (0.16e2 / 0.21e2 * t1981 + t2162) * t1975 + t1979 / 0.7e1 + t2162 * t1981 + t1983 - 0.3e1 / 0.7e1 * t2178 + t1969 / 0.42e2;
t2164 = t1981 / 0.3e1 + t1984;
t1906 = -t1970 / 0.4e1;
t2271 = t1906 + t1975 / 0.2e1;
t1842 = t2164 + t2271;
t1907 = -t1970 / 0.3e1;
t1909 = -0.2e1 / 0.3e1 * t1970;
t1913 = 0.4e1 / 0.3e1 * t1975;
t1919 = 0.4e1 / 0.3e1 * t1981;
t1741 = -0.8e1 / 0.3e1 * t1842 * t2073 + 0.5e1 / 0.18e2 * t1993 + (t1919 + t1907) * t1975 + t1983 - t1979 / 0.3e1 + t1969 / 0.18e2 + (t1913 + 0.2e1 / 0.3e1 * t1981 + t1909) * t1984;
t2152 = t1979 + t1983;
t1958 = 2 * t1984;
t2156 = t1958 - t1970;
t2018 = (t2156 * t1981) + t1969 / 0.6e1 + t2152 - t2178;
t1802 = -t1993 / 0.6e1 + t2018;
t1921 = t1981 / 0.2e1;
t2163 = t1921 + t1984;
t1806 = -0.2e1 / 0.3e1 * t2073 + t1906 + t2163;
t1956 = 4 * t1984;
t1872 = (t1956 + t1970) * t1981;
t1918 = -0.2e1 / 0.3e1 * t1975;
t1876 = t1918 + t1984;
t1877 = -t1981 / 0.3e1 + t1984;
t1886 = -t1981 + t1984;
t1908 = -t1970 / 0.2e1;
t1851 = t1908 + t2081;
t2079 = -0.4e1 * t2129;
t2052 = t1851 * t2079;
t1893 = t1930 * t1894;
t1988 = pkin(1) * t1981;
t2193 = t1893 * t1988;
t2143 = pkin(7) * t2193;
t2077 = 0.16e2 * t2143;
t1892 = t1894 ^ 2;
t2195 = t1892 * t1979;
t2140 = 0.8e1 * t2195;
t2255 = pkin(7) * t1930;
t2174 = 0.6e1 * t2255;
t2179 = t1981 * t1894;
t2184 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t1713 = t1876 * t2140 + t1806 * t2077 + 0.14e2 * t1739 * t2179 + (t1886 * t1993) + (t1872 - 0.10e2 / 0.3e1 * t1979 + (2 * t1983) - t2178) * t1975 + t1802 * t2184 + (t1741 * t2174 + t1877 * t2052) * pkin(1);
t1896 = 0.10e2 / 0.3e1 * t1981;
t2014 = 0.5e1 / 0.6e1 * t1993 + t2018;
t1771 = (t1896 + t2156) * t1975 + t2014;
t1910 = -0.3e1 / 0.2e1 * t1970;
t1932 = 15 * t1979;
t1939 = 18 * t1984;
t1953 = 3 * t1983;
t2165 = t1969 / 0.2e1 - t1993 / 0.2e1;
t2053 = -(3 * t2178) + t1953 + t2165;
t1957 = 3 * t1984;
t2161 = 15 * t1981 + t1957;
t1992 = pkin(3) * t1975;
t1972 = t1992 ^ 2;
t2170 = t1885 * ((t1910 + t1958) * t1981 - 0.3e1 / 0.2e1 * t2178 + t2152 + t2165) + t1972;
t1727 = -0.6e1 * t1771 * t2073 + (t1932 + ((t1939 - 9 * t1970) * t1981) + t2053) * t1975 + (t1910 + t2161) * t1993 + t2170;
t1940 = -2 * t1970;
t1954 = 8 * t1984;
t2258 = 4 * t1979;
t1789 = t1988 * t2079 + t1895 + t2258 + ((t1940 + t1954) * t1981);
t2035 = t1984 - t2073;
t1800 = -t1981 + t2035 + t2271;
t1887 = -3 * t1981 + t1984;
t2078 = -0.2e1 * t2129;
t2091 = 0.8e1 * t2143;
t2175 = 0.4e1 * t2255;
t1729 = t2091 + t1789 * t1894 + t1851 * t1887 + (t1800 * t2175 + t2078 * t2184) * pkin(1);
t1917 = -t1975 / 0.3e1;
t1875 = t1917 + t1984;
t2036 = pkin(1) * t2052;
t2151 = t1983 - t1979;
t1731 = t1875 * t2036 - t1972 + (-t1896 - t2149) * t1993 + (t1872 + t1993 / 0.6e1 - t1969 / 0.6e1 + t2151) * t1975 + t1802 * t1984;
t1941 = -5 * t1970;
t1947 = 7 * t1979;
t1735 = (t1910 + t1957 + (7 * t1981)) * t1993 + (t1947 + ((t1941 + 10 * t1984) * t1981) + t2053) * t1975 + t2170;
t1997 = pkin(7) * t1984;
t1867 = -0.12e2 * pkin(7) * t1988 + t1997 * t2236;
t2177 = t1984 * t1981;
t1874 = -8 * t1979 + 12 * t2177;
t1753 = t1867 * t1930 + t1874 * t1894 + t2077 + t2140 + t2152 - (6 * t2177);
t1884 = -3 * t1975 + t1984;
t2185 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t1766 = t1844 * t2185 + t1851 * t1884;
t2150 = t1983 + t1993;
t1835 = 16 * (t2150 - 6 * t2180) * t1979;
t1857 = -t1970 + t2081;
t1861 = t1881 + pkin(7);
t1882 = -30 * t1970 + 60 * t1984;
t1888 = t1890 ^ 2;
t1942 = -6 * t1970;
t1949 = 6 * t1981;
t2155 = t1969 - t1993;
t2024 = 6 * t1983 + t2155 - 6 * t2178;
t2076 = 0.32e2 * t2143;
t1889 = t1926 * t1890;
t2197 = t1889 * t1992;
t1878 = t1885 ^ 2;
t2199 = t1878 * (-t1975 + t1868);
t1746 = t2036 + ((t1949 + t2156) * t1975) + t2014;
t1815 = t1875 * t1844;
t1764 = t1851 * t2185 + t1815;
t1834 = t1884 * t2091;
t2148 = pkin(7) * t1881;
t2117 = 0.6e1 * t2148;
t2120 = 0.12e2 * t2179;
t1714 = t1746 * t2117 + t1764 * t2120 + t1727 + t1834;
t2281 = 0.8e1 * t1714;
t1680 = t1835 * t1892 + t1766 * t2076 + 0.24e2 * t1731 * t2179 + ((t1940 + t1956 + 28 * t1981) * t1972) + (t1857 * t2199) + (0.24e2 * t1713 * t1890 + (t1882 * t1979) + (t1942 * t1983) + (t2024 * t1949) + (t2155 * t1958) + (28 * t1988 ^ 2) + 0.4e1 * t1997 ^ 2) * t1975 + (0.32e2 * t1729 * t2197 + t1879 * t2281) * t1861 + 0.8e1 * (t1727 * t2255 - t1735 * t2129) * pkin(1) + (0.16e2 * t1753 * t1888 + (t1882 * t1981) + (70 * t1979) + t1993 + t2024) * t1993;
t1905 = -t1970 / 0.6e1;
t1838 = 0.7e1 / 0.6e1 * t1975 + t1905 + t2163;
t1915 = t1975 / 0.3e1;
t2086 = t1905 + t1915 + t1984;
t1841 = t1919 + t2086;
t2239 = t1927 * pkin(1);
t1761 = -t1838 * t2239 + t1841 * t1880;
t2087 = t1970 / 0.3e1 + t1915 + t1958;
t1767 = (t1975 * t1886) - 0.5e1 / 0.3e1 * t1979 + t2087 * t1981 + t1984 * (t1907 + t1875);
t1848 = t1981 + t2086;
t1871 = t1944 + t1886;
t1773 = t1848 * t1880 - t1871 * t2239 / 0.2e1;
t2194 = t1893 * t1979;
t2107 = t1927 * t2194;
t2137 = 0.4e1 * t2179;
t1897 = -0.20e2 / 0.3e1 * t1981;
t1914 = 0.2e1 / 0.3e1 * t1975;
t2088 = 0.2e1 / 0.3e1 * t1970 + t1914 + t1956;
t2089 = 0.4e1 / 0.3e1 * t1970 + t1913 + t1959;
t2249 = t1993 / 0.2e1 - (t1897 + t2088) * t1975 / 0.2e1 + 0.3e1 / 0.2e1 * t1979 - t2089 * t1981 / 0.2e1 - t1983 / 0.2e1;
t1715 = t2107 * t2267 + t1761 * t2137 + (-0.8e1 / 0.3e1 * t2195 + t1767) * t1880 + (t1773 * t2175 + t1927 * t2249) * pkin(1);
t1948 = 5 * t1979;
t1955 = 6 * t1984;
t2159 = t1940 + t2291;
t2082 = t1955 + t2159;
t1934 = 10 * t1981;
t2160 = t1934 + t1958;
t1757 = t1993 + (t1909 + t1918 + t2160) * t1975 + t1948 + (t2082 * t1981) + t1984 * (t1909 + t1876);
t1750 = t1757 * t1880;
t2085 = t1907 + t1885;
t1845 = 0.8e1 / 0.3e1 * t1975 + t2085;
t2157 = t1950 + t1984;
t1847 = t1907 + t1914 + t2157;
t1765 = -t1845 * t2239 + t1847 * t1880;
t1852 = 0.5e1 / 0.6e1 * t1975 + t1921 + t1905;
t2106 = t1927 * t2185;
t2145 = 0.2e1 * t1880;
t1779 = pkin(1) * t2106 + t1852 * t2145;
t2126 = t1988 * t1880;
t2138 = -0.4e1 * t2179;
t1943 = 5 * t1993;
t2083 = t1909 + t1885;
t1772 = t1943 + ((t1934 + t2082) * t1975) + (t1918 + t2083) * t1885;
t2211 = t1772 * t1927;
t2268 = -0.8e1 * pkin(7);
t1716 = t1893 * t2126 * t2268 + t1779 * t2138 + t1750 + (t1765 * t2175 - t2211) * pkin(1);
t1732 = -pkin(1) * t2211 + t1750;
t1770 = -(3 * t1993) + (t1897 + t2089) * t1975 + t2088 * t1981 + t2151;
t2280 = 0.10e2 / 0.3e1 * t1993 + (-t1981 + t2087) * t2291 + (t1917 + t2085) * t1959;
t1733 = t1770 * t1880 + t2239 * t2280;
t2270 = t1957 - t1970 - t1975;
t1853 = t2270 * t1934;
t2158 = t1941 - 5 * t1975;
t1734 = t1972 + (21 * t1981 + t2270) * t1993 + (t1984 * t2159 + t1853 + t1953 + 35 * t1979) * t1975 + (t1947 + (t1954 + t2158) * t1981 + t1984 * (-t1975 + t2149)) * t1885;
t2127 = t1981 * t1880;
t2022 = -t1927 * t1988 + t2127;
t1824 = 0.4e1 * t2022;
t1836 = t1880 + 0.2e1 * t2239;
t1869 = t1975 + t2157;
t1850 = t1908 + t1869;
t1736 = t1887 * t1880 + t1824 * t1894 + (t1836 * t2255 + t1850 * t1927) * t2283;
t1744 = 7 * t1972 + (35 * t1981 + 15 * t1984 + t2158) * t1993 + (21 * t1979 + t1853 + 9 * t1983 + (t1942 - 6 * t1975) * t1984) * t1975 + t2199;
t1839 = t1984 + 0.5e1 / 0.2e1 * t1975 + 0.3e1 / 0.2e1 * t1981 + t1908;
t2123 = (pkin(1) * t1884) / 0.2e1;
t1775 = t1839 * t1880 + t1927 * t2123;
t1863 = pkin(7) * t2147;
t1808 = 0.4e1 / 0.3e1 * t2179 + t1863 + t1877;
t2198 = t1888 * t1993;
t2064 = -0.24e2 * t1808 * t2198;
t2109 = t1861 * t2197;
t2071 = -0.8e1 * t2109;
t2122 = -0.12e2 * t2196;
t2139 = -0.6e1 * t2179;
t2201 = t1861 * t1926;
t1689 = t1775 * t2077 + t1736 * t2071 + t1715 * t2122 + t1733 * t2139 + (-0.6e1 * t1716 * t2201 + (0.24e2 * t1875 * t2195 - t1734) * t1929) * pkin(3) + (-0.6e1 * t1732 * t2255 + (t1744 + t2064) * t1927) * pkin(1);
t1825 = t1879 + t1861;
t1662 = t1680 * t1825 + t1689 * t1985;
t1793 = -0.4e1 / 0.9e1 * t2073 + 0.4e1 / 0.9e1 * t1975 - t1970 / 0.9e1 + t2164;
t1807 = t1905 + t1914 + t2035;
t1846 = t1913 + t2085;
t1728 = 0.4e1 * t2143 + 0.6e1 * t1793 * t2179 + t1846 * t2184 + (t1807 * t2175 + t1877 * t2078) * pkin(1);
t2084 = t1909 + t1914 + t1958;
t2167 = (t1914 + t2083) * t1885 + t1993;
t1743 = t1846 * t2062 + (t1949 + t2084) * t1975 + t2167;
t1801 = t1844 + t1846;
t1865 = t2153 * t2279;
t2146 = 0.4e1 * t1881;
t1730 = pkin(7) * t1801 * t2146 + t1865 * t1894 + t1743;
t1762 = t1846 * t2185 + t1815;
t1763 = (t1896 + t2084) * t1975 + t2167;
t1814 = t1863 + t2137 + t1887;
t2121 = 0.12e2 * t2196;
t2235 = 6 * pkin(1);
t1697 = t1728 * t2121 + t1834 + t1762 * t2120 + t1972 + ((-t1970 + t1975 + t2161) * t1993) + ((t1932 + (t1939 + t1942 + 6 * t1975) * t1981 + t1953 + (t1940 + t1944) * t1984) * t1975) + (t1878 * t1857) + (0.6e1 * t1730 * t1879 + 0.8e1 * t1814 * t2197) * t1861 + (t1743 * t2255 - t1763 * t2129) * t2235;
t1833 = t1869 * t1880;
t1837 = t1880 - t2239;
t1870 = 3 * t1975 + t1885;
t2200 = t1870 * t1927;
t1742 = -0.2e1 * t1894 * t2127 + t1833 + (0.2e1 * t1837 * t2255 - t2200) * pkin(1);
t1823 = 0.2e1 * t2022;
t1951 = 2 * t1981;
t1883 = t1951 + t1975;
t2118 = pkin(7) * t2145;
t1745 = t1886 * t1880 + t1823 * t1894 + (t1883 * t1927 + t1930 * t2118) * pkin(1);
t1782 = -pkin(1) * t2200 + t1833;
t2278 = 8 * t1981;
t1859 = pkin(3) * t2258 + t1992 * t2278;
t2257 = 4 * t1988;
t1785 = t1859 * t1929 + t2106 * t2257;
t1810 = t1975 * t2160 + t1948 + t2150 + 6 * t2177;
t1819 = t1943 + (t1934 + t1955) * t1975 + t1878;
t2141 = -0.4e1 * t2196;
t1710 = t1745 * t2141 + t1785 * t1894 + (-0.4e1 * t1742 * t2201 + (-t1810 + t2091) * t1929) * pkin(3) + (-0.4e1 * t1782 * t2255 + (t1819 + t2071) * t1927) * pkin(1);
t1688 = t1697 * t1825 + t1710 * t1985;
t1686 = 0.1e1 / t1688;
t2253 = t1707 / 0.4e1;
t2101 = t1686 * t2253;
t2248 = -t1985 / 0.4e1;
t2016 = t1662 * t2101 + t1709 * t2248;
t1971 = 1 / pkin(4);
t2181 = t1971 / pkin(3) ^ 2;
t2110 = t1758 * t2181;
t1647 = t2016 * t2110;
t2252 = t1709 / 0.4e1;
t2100 = t1686 * t2252;
t2247 = t1985 / 0.4e1;
t2017 = t1662 * t2100 + t1707 * t2247;
t1648 = t2017 * t2110;
t2029 = t2187 + t2191;
t2190 = t1926 * t1930;
t2030 = -t2189 + t2190;
t1625 = -t1647 * t2029 - t1648 * t2030;
t1626 = t1647 * t2030 - t1648 * t2029;
t1614 = atan2(t1625, t1626);
t1611 = sin(t1614);
t1612 = cos(t1614);
t2182 = t1971 * t1976;
t2092 = t2182 / 0.2e1;
t1653 = atan2(t1985 * t2092, t1662 * t1686 * t2092);
t1651 = sin(t1653);
t1652 = cos(t1653);
t1621 = t1651 * t1679 + t1652 * t1678;
t2070 = -t1651 * t1678 + t1679 * t1652;
t2286 = t1611 * t1621 - t1612 * t2070;
t2285 = t1611 * t2070 + t1612 * t1621;
t2284 = t1688 ^ 2;
t2238 = 0.2e1 * pkin(7);
t2272 = qJD(1) - qJD(2);
t2047 = qJD(1) * t2074;
t2229 = qJD(2) * t1929;
t2113 = t1861 * t2229;
t2068 = pkin(3) * t2113;
t2269 = t2047 - t2068;
t2112 = qJD(2) * t2191;
t2015 = t2105 - t2112;
t2266 = -4 * pkin(1);
t2265 = 0.1e1 / t1673;
t2264 = 0.1e1 / t1722;
t2263 = 0.1e1 / t1673 ^ 2;
t2261 = -0.8e1 * t1930;
t2260 = -0.2e1 * t1930;
t2254 = t1686 / 0.4e1;
t2251 = t1722 / 0.4e1;
t2250 = t1758 / 0.2e1;
t2246 = pkin(1) * t1894;
t1854 = qJD(2) * t2118;
t2228 = qJD(2) * t1975;
t2055 = t1926 * t1929 * t2228;
t2043 = 0.4e1 * t2055;
t1805 = t1854 + t2043;
t2231 = qJD(1) * t1930;
t2114 = t1927 * t2231;
t2056 = t1813 * t2114;
t2166 = 0.2e1 * t2289;
t2090 = t1854 + t2166;
t2206 = t1805 * t1894;
t2233 = qJD(1) * t1927;
t2217 = ((0.8e1 * t2056 - 0.4e1 * t2206) * t1981 + (t1866 * t2262 + t2119) * t2229 + (t1927 * t1929 ^ 2 * t2228 * t2238 + (t1788 * t2233 - t1930 * t2090) * t1860 + (t1822 * t2105 + (-t1788 * t2187 - t1822 * t2191) * qJD(2)) * pkin(3)) * t2236) / t1985;
t2099 = t2217 / 0.2e1;
t2186 = t1930 * t1985;
t2188 = t1927 * t1985;
t2232 = qJD(1) * t1929;
t1683 = t1794 * t2099 + (t1747 * t1930 - t1860 * t2188) * qJD(1) + (t1805 * t1930 - t1813 * t2233) * t1927 * t2283 + ((-t2186 + (-t1821 - 0.8e1 * t2133) * t1927) * t2232 + ((-t1751 + t2188) * t1926 + (t2186 + (t1860 * t2238 + t1821) * t1927 + (t1894 * t2283 - pkin(1) + t2255) * t1929 * t2282) * t1929) * qJD(2)) * pkin(3);
t2011 = 0.4e1 * t2015;
t2008 = t2011 * t1881;
t1684 = t1796 * t2099 + (-t1854 * t1930 + (t1783 * t1927 + t2186) * qJD(1)) * t1860 + (t2043 + 0.4e1 * t2056 - 0.2e1 * t2206) * pkin(1) + ((t1856 * t1930 - t2188) * t2232 + (-t1783 * t2187 - t1856 * t2191 - t1985 * t2030) * qJD(2) + (t1843 * t2229 + t1879 * t2015) * t2238 + t1860 * t2008) * pkin(3);
t1803 = t1811 * qJD(1);
t1804 = t1812 * qJD(1);
t2116 = t1803 * t2283;
t2214 = (t1776 + t1777) * pkin(5) * t2116 / t1987;
t2098 = t2214 / 0.2e1;
t2019 = -t1804 * t1987 + t1811 * t2098;
t1698 = (-(t1797 * t2283 - t1790) * t1803 + t2019) * pkin(5);
t2142 = 0.2e1 * t1803 * t1811;
t2207 = t1803 * t1987;
t1699 = t1797 * t2098 + t1967 * pkin(1) * t2142 + (t1790 * t1804 + t2207) * pkin(5);
t2213 = ((-qJD(1) * t2203 + qJD(2) * t2128) * t2283 + t2090) / t1760 ^ 2;
t2058 = t2210 * t2213;
t1781 = 0.1e1 / t1784 ^ 2;
t2208 = t1781 * t1803;
t2134 = pkin(1) * t2208;
t1637 = (-t2021 * t2058 + (-(t2219 / 0.2e1 + t2216 / 0.2e1) * t2134 + (t1684 * t2251 + t1698 * t2253 + t1683 * t1724 / 0.4e1 + t1699 * t2252) * t2210) * t1758) * t1976;
t1638 = (t2020 * t2058 + (-(-t2218 / 0.2e1 + t2215 / 0.2e1) * t2134 + (-t1684 * t1724 / 0.4e1 - t1707 * t1699 / 0.4e1 + t1698 * t2252 + t1683 * t2251) * t2210) * t1758) * t1976;
t1615 = t1637 * t1931 + t1638 * t1928;
t2245 = pkin(2) * t1615;
t2013 = t2015 * pkin(3);
t1768 = t1981 * t2013 * t2175;
t1864 = pkin(1) * t2233;
t1820 = pkin(3) * t2229 - t1864;
t2080 = pkin(7) * t1864;
t1855 = -0.2e1 * t2080;
t1891 = t1927 ^ 2;
t2005 = pkin(1) * t1875 * t2013 * t2179;
t2010 = t2015 * t2283;
t2230 = qJD(2) * t1926;
t2131 = pkin(3) * t2230;
t2136 = pkin(1) * t2231;
t2012 = -t1757 * t2131 - t1772 * t2136;
t2034 = t1890 * t1992 * t2113;
t2023 = -0.24e2 * t2034;
t2025 = t1877 * t2045;
t2026 = t1877 * t2046;
t2192 = t1894 * t1988;
t2144 = pkin(7) * t2192;
t2044 = t2144 * t2233;
t2028 = -0.48e2 * t2044;
t2054 = t1981 * t2114;
t2039 = -0.24e2 * t2054;
t2040 = t2131 * t2193;
t2041 = t2131 * t2195;
t2115 = qJD(1) * t2197;
t2069 = pkin(1) * t2115;
t2042 = t1927 * t2069;
t2059 = t1889 * t1993 * t2229;
t2060 = t2185 * t2231;
t2061 = qJD(1) * t2107;
t2063 = -0.6e1 * t2080;
t2111 = qJD(2) * t2190;
t2067 = pkin(3) * t2111;
t2075 = t2179 * t2256;
t2108 = t1927 * t2198;
t2130 = pkin(3) * t2201;
t2027 = -0.24e2 * t2044;
t2171 = t1884 * t2027 - 0.24e2 * t2005;
t2172 = t2290 * t1771;
t1627 = (t2064 * t2136 - 0.24e2 * pkin(1) * (-0.8e1 / 0.3e1 * t2054 + t1855) * t2108 - 0.96e2 * t1808 * t2059 * t2239 + ((-t1887 + t2138) * t2131 + 0.2e1 * (pkin(1) * t1850 - t1824 * t1927 - 0.2e1 * t2192) * t2231 + (-0.2e1 * t2067 + (-0.2e1 * t1836 * t1927 + 0.4e1 * t2246) * qJD(1)) * pkin(1) * pkin(7)) * t2071 + (0.8e1 / 0.3e1 * t2041 + (-t1838 * t2136 - t1841 * t2131) * t2137 - t1767 * t2131 + t2136 * t2249 + (0.32e2 / 0.3e1 * t2194 * t1880 + t1981 * t1761 * t2261) * t2233 + (t1848 * t2067 * t2266 + ((0.12e2 * t1891 * t1894 - 0.4e1 * t1892) * t1979 + (-0.4e1 * t1773 * t1927 - 0.2e1 * t1871 * t2246) * pkin(1)) * qJD(1)) * pkin(7)) * t2122 - 0.24e2 * t1715 * t2055 - 0.6e1 * ((-0.4e1 * (pkin(1) * t2060 - 0.2e1 * t1852 * t2131) * t1894 + 0.8e1 * t1779 * t2114) * t1981 + (0.8e1 * t2040 + (-t1845 * t2136 - t1847 * t2131) * t2146 + (t1765 * t2266 + 0.24e2 * t1894 * t2126) * t2233) * pkin(7) + t2012) * t2130 + (-t1839 * t2131 + t2123 * t2231) * t2077 + t1775 * t2028 + (-t1770 * t2131 + t2136 * t2280) * t2139 + 0.12e2 * t1733 * t2054 - 0.6e1 * t2012 * t2148 + 0.6e1 * t1732 * t2080 + t1734 * t2131 + t1744 * t2136 + (-0.96e2 * t2061 * t1880 - 0.24e2 * t2041) * t1875 + (0.8e1 * t2042 + t2023) * t1736 + 0.6e1 * t2269 * t1716) * t1985 + t1689 * t2099 + (0.16e2 * (t1874 * t2260 - t1867 - 0.48e2 * t2144 - 0.32e2 * t2194) * qJD(1) * t2108 + 0.64e2 * t1753 * t2059 + 0.32e2 * (-t1768 + (t1789 * t2260 + (t1800 * t2266 - 0.24e2 * t2192) * pkin(7)) * t2233 + (-t2010 * t2184 - t2011 * t2192) * pkin(3)) * t2109 + 0.24e2 * (-0.32e2 * t1876 * t2061 + (-0.2e1 / 0.3e1 * t2105 + 0.2e1 / 0.3e1 * t2112) * t2077 * t2256 + t1806 * t2028 - 0.28e2 * t1739 * t2054 + (-0.8e1 / 0.3e1 * t2105 + 0.8e1 / 0.3e1 * t2112) * t1842 * pkin(3) * t1981 * t2174 + t1741 * t2063 + 0.4e1 * (-t2026 + t2025) * t1851 - 0.64e2 / 0.3e1 * t2015 * t1840 * t2075) * t2196 + 0.48e2 * t1713 * t2055 - 0.8e1 * t1714 * t2047 + 0.8e1 * (t1764 * t2039 + (-pkin(3) * t1851 * t2008 - t1746 * t2233) * pkin(7) * t2235 + t2171 + t2172) * t2130 + t2068 * t2281 - 0.4e1 * t1835 * t1893 * t2233 - pkin(3) * t2010 * t2076 * t2185 - 0.96e2 * t1766 * t2044 - 0.96e2 * t1851 * t2005 - 0.48e2 * t1731 * t2054 + 0.8e1 * t2172 * t2148 - 0.8e1 * t1727 * t2080 + 0.8e1 * t2289 * t1735 + (-0.32e2 * t2042 + 0.96e2 * t2034) * t1729) * t1825 + t1680 * t1820;
t1706 = 0.1e1 / t1707 ^ 2;
t2097 = -t2213 / 0.2e1;
t1639 = qJD(1) + 0.2e1 * ((t1683 * t2250 + t1709 * t2097) / t1707 - (t1684 * t2250 + t1707 * t2097) * t1709 * t1706) * pkin(3) / (t1706 * t1709 ^ 2 + 0.1e1) * t1760 * t1976;
t2221 = 0.1e1 / t1662 ^ 2 * t2284;
t2065 = t1870 * t2136;
t2072 = -0.2e1 * t2114;
t2169 = 0.4e1 * t2289 * t1846;
t2222 = ((t1861 * t2069 * t2261 + t1891 * t2115 * t2278 + t2023 * t2239 + (0.2e1 * (-t1981 * t2131 - t1988 * t2231) * t1894 + t1823 * t2072 - t1886 * t2131 + (t1883 * t2231 + (-qJD(1) * t2189 - t2111) * pkin(7) * t2282) * pkin(1)) * t2141 - 0.8e1 * t1745 * t2055 - 0.4e1 * (-t2065 + (0.4e1 * t1929 * t2054 + (-t1869 + 0.2e1 * t2179) * t2230) * pkin(3) + ((-t2131 - t2136) * t1881 - t1837 * t1864) * t2238) * t2130 + t2040 * t2268 + t2027 * t1880 + (-t1859 * t2230 + t2060 * t2257) * t1894 + t1785 * t2072 - 0.4e1 * (-t1869 * t2131 - t2065) * t2148 + 0.4e1 * t1782 * t2080 + t1810 * t2131 + t1819 * t2136 + 0.4e1 * t2269 * t1742) * t1985 + t1710 * t2099 + (0.8e1 * (t1855 - 0.8e1 * t2054) * t2109 + (-0.12e2 * t2044 + 0.6e1 * (-0.4e1 / 0.9e1 * t2105 + 0.4e1 / 0.9e1 * t2112) * t2075 - 0.12e2 * t1793 * t2054 - t1768 - 0.4e1 * t1807 * t2080 - 0.2e1 * t2026 + 0.2e1 * t2025) * t2121 + 0.24e2 * t1728 * t2055 + 0.6e1 * (t1865 * t2072 + (-t1801 * t2233 + t1930 * t2166) * pkin(7) * t2236 + t2169) * t2130 + t1762 * t2039 + t2169 * t2117 + t1743 * t2063 + t2171 + (-0.8e1 * t2042 + 0.24e2 * t2034) * t1814 + t2290 * t1763 - 0.6e1 * t2269 * t1730) * t1825 + t1697 * t1820) / t2284;
t1602 = t1639 + (0.1e1 / t1662 * t1688 * t2099 - 0.2e1 * (t1627 * t1686 / 0.2e1 - t1662 * t2222 / 0.2e1) * pkin(4) * pkin(3) * t1985 * t2182 * t2221) / (t1719 * t2221 + 0.1e1);
t2244 = pkin(4) * t1602;
t1616 = t1637 * t1928 - t1638 * t1931;
t1617 = 0.1e1 / t1618 ^ 2;
t2227 = (-t1663 - t1664) * pkin(6) * t2245 / t2277;
t2104 = -t2227 / 0.2e1;
t2132 = 0.1e1 / t1667 ^ 2 * t2245;
t2225 = t1615 * t1986;
t2226 = t1615 * t1675;
t1589 = t1639 + 0.2e1 * (((t1669 * t2104 + pkin(2) * t2226 * t2259 + (t1616 * t1668 + t2225) * pkin(6)) * t2102 - t1619 * t2132) / t1618 - (-t1618 * t2132 + (t1675 * t2104 - t1986 * t1616 + (-0.2e1 * pkin(2) * t1669 + t1668) * t1615) * pkin(6) * t2102) * t1619 * t1617) * pkin(6) / (t1617 * t1619 ^ 2 + 0.1e1) * t1667;
t2242 = pkin(6) * t1589;
t2224 = t1639 * t1678;
t2223 = t1639 * t1679;
t2209 = t1780 / pkin(1);
t2125 = pkin(5) * t2208;
t2103 = -t2222 / 0.4e1;
t2095 = t2210 / 0.2e1;
t2094 = -t2209 / 0.2e1;
t2093 = t2209 / 0.2e1;
t2051 = pkin(3) * t2223 - t2136;
t2050 = pkin(2) * t2223 - t2136;
t2049 = -pkin(2) * t2224 + t1864;
t2048 = -pkin(3) * t2224 + t1864;
t1646 = atan2(t1675, t1673);
t1643 = sin(t1646);
t1644 = cos(t1646);
t2032 = t1643 * t1679 + t1644 * t1678;
t2031 = t1643 * t1678 - t1644 * t1679;
t1832 = -rSges(2,1) * t1930 + rSges(2,2) * t1927;
t1831 = -rSges(8,1) * t1929 + rSges(8,2) * t1926;
t1830 = -rSges(2,1) * t1927 - rSges(2,2) * t1930;
t1829 = rSges(8,1) * t1926 + rSges(8,2) * t1929;
t1798 = -pkin(1) * t1812 + pkin(5);
t1791 = t1799 + t1951;
t1756 = t2272 * t2029;
t1755 = t2272 * t2030;
t1725 = pkin(1) * t1791 * t1811 + t1798 * t1987;
t1723 = -pkin(1) * t2205 + t1791 * t1798;
t1721 = 0.1e1 / t1723 ^ 2;
t1720 = 0.1e1 / t1722 ^ 2;
t1705 = atan2(t1725 * t2093, t1723 * t2094);
t1704 = atan2(t1724 * t2095, t1722 * t2095);
t1703 = cos(t1705);
t1702 = cos(t1704);
t1701 = sin(t1705);
t1700 = sin(t1704);
t1693 = -t1700 * t1927 + t1702 * t1930;
t1692 = t1700 * t1930 + t1702 * t1927;
t1691 = -t1701 * t1923 + t1703 * t1924;
t1690 = t1701 * t1924 + t1703 * t1923;
t1677 = rSges(6,1) * t1691 - rSges(6,2) * t1690;
t1676 = rSges(6,1) * t1690 + rSges(6,2) * t1691;
t1660 = (-((-t1798 * t2214 / 0.2e1 + t1981 * pkin(5) * t2142 + (t1791 * t1804 + t2207) * pkin(1)) * t2093 - t1725 * t2125) / t1723 - (t1723 * t2125 + (-(-0.2e1 * pkin(5) * t1798 - t1791) * t1803 + t2019) * pkin(1) * t2094) * t1725 * t1721) / (t1721 * t1725 ^ 2 + 0.1e1) * t1784 * t2283;
t1658 = qJD(1) + (t1699 * t2264 * t2210 + (-t1698 * t1720 * t2210 - (-t1722 * t1720 + t2264) * t1781 * t2116) * t1724) * t1784 / (t1720 * t1724 ^ 2 + 0.1e1) * pkin(5);
t1650 = -t2136 + t1658 * (rSges(9,1) * t1693 - rSges(9,2) * t1692);
t1649 = t1864 - t1658 * (rSges(9,1) * t1692 + rSges(9,2) * t1693);
t1629 = -t2136 + t1639 * (rSges(3,1) * t1679 - rSges(3,2) * t1678);
t1628 = t1864 - t1639 * (rSges(3,1) * t1678 + rSges(3,2) * t1679);
t1624 = 0.1e1 / t1626 ^ 2;
t1605 = (-t2017 * t2213 + (t1684 * t2247 + t1707 * t2217 / 0.8e1 + t1627 * t2100 + (t1683 * t2254 + t1709 * t2103) * t1662) * t1758) * t2181;
t1604 = (-t2016 * t2213 + (t1627 * t2101 + t1683 * t2248 - t1709 * t2217 / 0.8e1 + (t1684 * t2254 + t1707 * t2103) * t1662) * t1758) * t2181;
t1601 = (t1616 * t2265 - t2226 * t2263) / (t1675 ^ 2 * t2263 + 0.1e1) + t1639;
t1596 = t1602 * (rSges(5,1) * t2070 - rSges(5,2) * t1621) + t2051;
t1595 = -t1602 * (rSges(5,1) * t1621 + rSges(5,2) * t2070) + t2048;
t1594 = -t2136 + t1601 * (rSges(7,1) * t2031 + rSges(7,2) * t2032);
t1593 = t1864 - t1601 * (-rSges(7,1) * t2032 + rSges(7,2) * t2031);
t1590 = ((-t1604 * t2029 - t1605 * t2030 - t1647 * t1755 + t1648 * t1756) / t1626 - (t1604 * t2030 - t1605 * t2029 - t1647 * t1756 - t1648 * t1755) * t1625 * t1624) / (t1624 * t1625 ^ 2 + 0.1e1) + t1602;
t1588 = (-t2265 * t2227 / 0.4e1 + t2263 * t2225 / 0.2e1) / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t1991 * t2263 * t2220 / 0.4e1) * t2183 + t1589;
t1587 = t2070 * t2244 + t1590 * (rSges(11,1) * t2286 + rSges(11,2) * t2285) + t2051;
t1586 = -t1621 * t2244 - t1590 * (-rSges(11,1) * t2285 + rSges(11,2) * t2286) + t2048;
t1585 = t1589 * (rSges(4,1) * t2033 + rSges(4,2) * t1599) + t2050;
t1584 = -t1589 * (-rSges(4,1) * t1599 + rSges(4,2) * t2033) + t2049;
t1583 = t2033 * t2242 + t1588 * (rSges(10,1) * t2288 - t2287 * rSges(10,2)) + t2050;
t1582 = t1599 * t2242 - t1588 * (t2287 * rSges(10,1) + rSges(10,2) * t2288) + t2049;
t1 = t1590 ^ 2 * Icges(11,3) / 0.2e1 + m(10) * (t1582 ^ 2 + t1583 ^ 2) / 0.2e1 + m(9) * (t1649 ^ 2 + t1650 ^ 2) / 0.2e1 + m(5) * (t1595 ^ 2 + t1596 ^ 2) / 0.2e1 + m(4) * (t1584 ^ 2 + t1585 ^ 2) / 0.2e1 + m(7) * (t1593 ^ 2 + t1594 ^ 2) / 0.2e1 + m(3) * (t1628 ^ 2 + t1629 ^ 2) / 0.2e1 + m(11) * (t1586 ^ 2 + t1587 ^ 2) / 0.2e1 + t1589 ^ 2 * Icges(4,3) / 0.2e1 + t1639 ^ 2 * Icges(3,3) / 0.2e1 + t1602 ^ 2 * Icges(5,3) / 0.2e1 + t1588 ^ 2 * Icges(10,3) / 0.2e1 + t1658 ^ 2 * Icges(9,3) / 0.2e1 + t1601 ^ 2 * Icges(7,3) / 0.2e1 + (m(2) * (t1830 ^ 2 + t1832 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (m(8) * (t1829 ^ 2 + t1831 ^ 2) / 0.2e1 + Icges(8,3) / 0.2e1) * qJD(2) ^ 2 + (m(6) * (t1676 ^ 2 + t1677 ^ 2) / 0.2e1 + Icges(6,3) / 0.2e1) * t1660 ^ 2;
T = t1;