% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2TE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2TE_coriolismatJ_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:44:37
% EndTime: 2020-05-07 01:44:40
% DurationCPUTime: 2.28s
% Computational Cost: add. (2575->257), mult. (3311->363), div. (0->0), fcn. (1310->80), ass. (0->177)
t1676 = pkin(17) + pkin(18);
t1740 = pkin(15) + t1676;
t1725 = pkin(16) + t1740;
t1719 = qJ(2) + t1725;
t1640 = qJ(3) + t1719;
t1627 = qJ(4) + t1640;
t1720 = -qJ(2) + t1725;
t1641 = -qJ(3) + t1720;
t1630 = -qJ(4) + t1641;
t1813 = cos(t1630) + cos(t1627);
t1684 = sin(qJ(3));
t1689 = cos(qJ(2));
t1756 = t1684 * t1689;
t1652 = pkin(4) * t1756;
t1688 = cos(qJ(3));
t1661 = -pkin(4) * t1688 + pkin(1);
t1685 = sin(qJ(2));
t1576 = t1661 * t1685 - t1652;
t1757 = t1684 * t1685;
t1651 = pkin(4) * t1757;
t1577 = t1661 * t1689 + t1651;
t1686 = sin(pkin(15));
t1690 = cos(pkin(15));
t1543 = t1576 * t1690 + t1686 * t1577;
t1681 = sin(pkin(16));
t1682 = cos(pkin(16));
t1730 = -t1686 * t1576 + t1577 * t1690;
t1812 = -t1543 * t1681 + t1682 * t1730;
t1811 = t1543 * t1682 + t1681 * t1730;
t1754 = t1688 * t1689;
t1601 = t1754 - t1757;
t1755 = t1688 * t1685;
t1602 = t1755 + t1756;
t1551 = t1686 * t1601 + t1602 * t1690;
t1729 = t1601 * t1690 - t1686 * t1602;
t1810 = -t1681 * t1551 + t1682 * t1729;
t1809 = t1551 * t1682 + t1681 * t1729;
t1795 = m(5) + m(6);
t1611 = pkin(4) * t1795 + m(4) * rSges(4,1) + m(9) * rSges(9,1);
t1808 = t1611 * t1684;
t1677 = qJ(4) + qJ(2);
t1672 = qJ(3) + t1677;
t1655 = sin(t1672);
t1768 = qJ(4) - qJ(2);
t1673 = qJ(3) - t1768;
t1656 = sin(t1673);
t1658 = cos(t1672);
t1659 = cos(t1673);
t1679 = qJ(3) + qJ(2);
t1668 = sin(t1679);
t1671 = cos(t1679);
t1784 = rSges(5,3) * m(5);
t1785 = rSges(4,3) * m(4);
t1788 = m(9) * rSges(9,3);
t1793 = m(6) * rSges(6,1);
t1798 = m(6) / 0.2e1;
t1538 = (rSges(4,1) * t1785 + rSges(9,1) * t1788 - Icges(4,5) - Icges(9,5)) * t1671 - t1668 * (rSges(4,2) * t1785 + rSges(9,2) * t1788 - Icges(4,6) - Icges(9,6)) + (t1784 * t1671 - t1656 * t1793 / 0.2e1 + (rSges(6,1) * t1655 + (t1658 + t1659) * rSges(6,2)) * t1798) * pkin(4);
t1646 = qJ(4) + t1725;
t1716 = -qJ(2) + t1646;
t1647 = -qJ(4) + t1725;
t1717 = qJ(2) + t1647;
t1803 = (sin(t1716) / 0.4e1 - sin(t1717) / 0.4e1) * pkin(1);
t1802 = (cos(t1717) / 0.4e1 + cos(t1716) / 0.4e1) * pkin(1);
t1801 = 0.2e1 * qJ(2);
t1674 = pkin(14) - qJ(2) - pkin(15);
t1800 = 0.2e1 * t1674;
t1799 = -pkin(8) / 0.2e1;
t1797 = pkin(1) * m(8);
t1796 = pkin(3) * m(9);
t1794 = m(4) * rSges(4,2);
t1792 = m(6) * rSges(6,2);
t1791 = m(7) * rSges(7,2);
t1790 = m(7) * rSges(7,3);
t1789 = m(9) * rSges(9,2);
t1787 = rSges(3,2) * m(3);
t1786 = rSges(6,2) * pkin(8);
t1683 = sin(qJ(4));
t1687 = cos(qJ(4));
t1631 = t1687 * rSges(6,1) - t1683 * rSges(6,2);
t1783 = t1631 / 0.2e1;
t1628 = -qJ(4) + t1640;
t1782 = cos(t1628);
t1629 = qJ(4) + t1641;
t1781 = cos(t1629);
t1780 = cos(t1674);
t1779 = sin(t1628);
t1778 = sin(t1629);
t1664 = sin(t1676);
t1665 = cos(t1676);
t1775 = pkin(4) * (-t1809 * t1664 + t1810 * t1665);
t1772 = m(6) * t1631;
t1771 = m(6) * (t1683 * rSges(6,1) + t1687 * rSges(6,2));
t1691 = -rSges(6,3) - pkin(10);
t1770 = rSges(6,1) * t1691;
t1769 = rSges(6,2) * t1691;
t1767 = pkin(18) + pkin(15);
t1680 = t1801 + qJ(3);
t1575 = pkin(1) * t1611 * sin(t1680);
t1654 = t1789 + t1794;
t1606 = pkin(1) * t1654 * cos(t1680);
t1607 = sin(t1627);
t1608 = sin(t1630);
t1648 = -m(5) * rSges(5,2) - t1691 * m(6);
t1657 = sin(t1674);
t1660 = pkin(8) * m(6) + m(5) * rSges(5,1);
t1727 = qJ(2) + t1740;
t1649 = qJ(3) + t1727;
t1728 = -qJ(2) + t1740;
t1650 = -qJ(3) + t1728;
t1662 = 0.2e1 * t1679;
t1744 = t1792 / 0.4e1;
t1733 = pkin(4) * t1744;
t1745 = t1793 / 0.4e1;
t1736 = pkin(4) * t1745;
t1746 = -t1793 / 0.4e1;
t1737 = pkin(4) * t1746;
t1748 = t1796 / 0.2e1;
t1749 = -t1796 / 0.2e1;
t1702 = (-rSges(4,1) * t1794 - rSges(9,1) * t1789 + Icges(4,4) + Icges(9,4)) * cos(t1662) - (-Icges(4,1) - Icges(9,1) + Icges(4,2) + Icges(9,2) + t1795 * pkin(4) ^ 2 + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t1662) / 0.2e1 + rSges(9,1) * sin(t1650) * t1748 + t1779 * t1737 + t1778 * t1736 + (-cos(t1640) / 0.2e1 + cos(t1641) / 0.2e1) * pkin(4) * t1648 + (-sin(t1640) / 0.2e1 + sin(t1641) / 0.2e1) * pkin(4) * t1660 + (t1782 + t1781) * t1733 + (rSges(9,1) * sin(t1649) + (cos(t1649) + cos(t1650)) * rSges(9,2)) * t1749 + (t1611 * t1668 + t1654 * t1671) * pkin(12);
t1642 = qJ(2) + t1646;
t1643 = -qJ(2) + t1647;
t1706 = ((cos(t1642) + cos(t1643)) * t1744 + sin(t1642) * t1745 + sin(t1643) * t1746) * pkin(1);
t1732 = m(4) + m(8) + m(9) + t1795;
t1741 = qJ(2) + t1767;
t1742 = -qJ(2) + t1767;
t1750 = t1797 / 0.2e1;
t1751 = -t1797 / 0.2e1;
t1528 = -t1606 - t1702 - ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t1800) / 0.2e1 + pkin(1) * sin(t1728) * t1748 + pkin(1) * sin(t1727) * t1749 + t1607 * t1736 + t1608 * t1737 - t1706 - t1575 + (rSges(7,1) * t1791 - Icges(7,4)) * cos(t1800) - (-rSges(3,1) * t1787 + Icges(3,4)) * cos(t1801) + ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + t1732 * pkin(1) ^ 2) * sin(t1801) / 0.2e1 + (cos(t1720) / 0.2e1 - cos(t1719) / 0.2e1) * pkin(1) * t1648 + (sin(t1720) / 0.2e1 - sin(t1719) / 0.2e1) * pkin(1) * t1660 + t1793 * t1803 + t1792 * t1802 + t1813 * t1733 + ((t1732 * pkin(1) + rSges(3,1) * m(3)) * t1685 + t1689 * t1787) * pkin(12) + (cos(t1741) * t1750 + cos(t1742) * t1751) * rSges(8,2) + (sin(t1742) * t1750 + sin(t1741) * t1751) * rSges(8,1) + (m(7) * rSges(7,1) * t1657 - t1780 * t1791) * pkin(6);
t1766 = t1528 * qJD(1);
t1625 = 0.2e1 * t1646;
t1626 = 0.2e1 * t1647;
t1723 = 0.2e1 * t1725;
t1644 = qJ(4) + t1723;
t1645 = -qJ(4) + t1723;
t1653 = 0.4e1 * t1770;
t1692 = rSges(6,1) * pkin(8);
t1698 = 0.2e1 * qJ(4);
t1700 = 4 * Icges(6,5);
t1724 = t1607 * t1737 + t1608 * t1736 - t1813 * pkin(4) * t1792 / 0.4e1;
t1529 = (cos(t1698) / 0.2e1 - cos(t1626) / 0.4e1 - cos(t1625) / 0.4e1) * (rSges(6,1) * t1792 - Icges(6,4)) + (t1687 * t1799 + (cos(t1646) / 0.2e1 + cos(t1647) / 0.2e1) * pkin(12) + (-t1782 / 0.4e1 - t1781 / 0.4e1) * pkin(4) + t1802) * t1792 + (t1683 * t1799 + (sin(t1646) / 0.2e1 - sin(t1647) / 0.2e1) * pkin(12) + (t1779 / 0.4e1 - t1778 / 0.4e1) * pkin(4) + t1803) * t1793 + t1706 + t1724 + (sin(t1698) / 0.4e1 - sin(t1625) / 0.8e1 + sin(t1626) / 0.8e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2)) - ((t1692 + t1769) * m(6) + Icges(6,6)) * sin(t1644) / 0.4e1 + ((t1692 - t1769) * m(6) - Icges(6,6)) * sin(t1645) / 0.4e1 + ((t1653 - 0.4e1 * t1786) * m(6) + t1700) * cos(t1644) / 0.16e2 - ((t1653 + 0.4e1 * t1786) * m(6) + t1700) * cos(t1645) / 0.16e2;
t1765 = t1529 * qJD(1);
t1759 = t1654 * t1688;
t1530 = (t1808 / 0.2e1 + t1759 / 0.2e1) * pkin(1) + t1606 / 0.2e1 + t1702 + t1575 / 0.2e1 + t1724;
t1764 = t1530 * qJD(1);
t1666 = sin(t1677);
t1667 = -sin(t1768);
t1669 = cos(t1677);
t1670 = cos(t1768);
t1710 = (-t1659 / 0.4e1 + t1658 / 0.4e1) * rSges(6,2) + (t1656 / 0.4e1 + t1655 / 0.4e1) * rSges(6,1);
t1709 = t1710 * pkin(4);
t1705 = t1709 + ((t1670 / 0.4e1 - t1669 / 0.4e1) * rSges(6,2) + (-t1667 / 0.4e1 - t1666 / 0.4e1) * rSges(6,1)) * pkin(1);
t1760 = t1631 * t1576;
t1536 = (-t1760 / 0.2e1 + t1705) * m(6);
t1763 = t1536 * qJD(1);
t1762 = t1538 * qJD(3);
t1539 = (t1652 * t1783 + (t1755 * t1783 + t1710) * pkin(4)) * m(6);
t1761 = t1539 * qJD(1);
t1743 = qJD(4) * t1772;
t1731 = t1771 * t1775;
t1722 = -pkin(4) * t1754 + t1651;
t1721 = pkin(4) * t1755 + t1652;
t1707 = -t1686 * t1722 + t1721 * t1690;
t1708 = t1721 * t1686 + t1722 * t1690;
t1703 = (t1681 * t1707 + t1708 * t1682) * t1665;
t1704 = t1664 * (-t1681 * t1708 + t1707 * t1682);
t1531 = (t1775 / 0.2e1 + t1703 / 0.2e1 + t1704 / 0.2e1) * t1771;
t1557 = (t1759 + t1808) * pkin(1);
t1715 = -t1557 * qJD(2) + t1531 * qJD(4);
t1599 = pkin(8) * t1771;
t1559 = -t1683 * (m(6) * t1770 + Icges(6,5)) - (m(6) * t1769 + Icges(6,6)) * t1687;
t1554 = t1557 * qJD(3);
t1542 = t1559 * t1690 + t1686 * t1599;
t1541 = t1686 * t1559 - t1599 * t1690;
t1540 = -t1721 * t1772 / 0.2e1 + m(6) * t1709;
t1537 = t1705 * m(6) + t1760 * t1798;
t1532 = -(t1703 + t1704) * t1771 / 0.2e1 + t1731 / 0.2e1;
t1 = [-t1528 * qJD(2) + t1530 * qJD(3) + t1529 * qJD(4), t1537 * qJD(4) + t1762 - t1766 + (Icges(3,5) * t1689 - t1685 * Icges(3,6) + (-rSges(3,1) * t1689 + rSges(3,2) * t1685) * rSges(3,3) * m(3) + (-rSges(7,2) * t1790 + Icges(7,6)) * t1657 - (rSges(7,1) * t1790 - Icges(7,5)) * t1780 + ((-m(8) * rSges(8,3) - t1784 - t1785 - t1788) * t1689 + ((-t1670 / 0.2e1 - t1669 / 0.2e1) * rSges(6,2) + (t1667 / 0.2e1 - t1666 / 0.2e1) * rSges(6,1)) * m(6)) * pkin(1) + t1538) * qJD(2), t1538 * qJD(2) + t1540 * qJD(4) + t1762 + t1764, t1765 + t1537 * qJD(2) + t1540 * qJD(3) + ((t1541 * t1682 + t1681 * t1542) * t1665 + (-t1681 * t1541 + t1542 * t1682) * t1664 + (pkin(12) + t1577) * t1771) * qJD(4); t1536 * qJD(4) + t1766, t1554, t1554 - t1715, t1763 - t1531 * qJD(3) - (t1664 * t1812 + t1811 * t1665) * t1743; t1539 * qJD(4) - t1764, t1715, 0, t1761 + t1531 * qJD(2) + pkin(4) * (t1810 * t1664 + t1809 * t1665) * t1743; -t1536 * qJD(2) - t1539 * qJD(3) - t1765, -t1763 - (-t1664 * t1811 + t1812 * t1665) * qJD(2) * t1771 + t1532 * qJD(3), t1532 * qJD(2) + qJD(3) * t1731 - t1761, 0;];
Cq = t1;
