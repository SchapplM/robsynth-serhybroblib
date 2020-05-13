% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2DE1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE1_coriolismatJ_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:00:27
% EndTime: 2020-05-07 02:00:30
% DurationCPUTime: 2.24s
% Computational Cost: add. (2575->258), mult. (3311->364), div. (0->0), fcn. (1310->80), ass. (0->179)
t1792 = m(5) + m(6);
t1678 = pkin(17) + pkin(18);
t1736 = pkin(15) + t1678;
t1720 = pkin(16) + t1736;
t1686 = sin(qJ(3));
t1691 = cos(qJ(2));
t1750 = t1686 * t1691;
t1648 = pkin(4) * t1750;
t1690 = cos(qJ(3));
t1658 = -pkin(4) * t1690 + pkin(1);
t1687 = sin(qJ(2));
t1577 = t1658 * t1687 - t1648;
t1751 = t1686 * t1687;
t1647 = pkin(4) * t1751;
t1578 = t1658 * t1691 + t1647;
t1688 = sin(pkin(15));
t1692 = cos(pkin(15));
t1553 = t1577 * t1692 + t1578 * t1688;
t1683 = sin(pkin(16));
t1684 = cos(pkin(16));
t1725 = -t1577 * t1688 + t1578 * t1692;
t1805 = -t1553 * t1683 + t1684 * t1725;
t1804 = t1553 * t1684 + t1683 * t1725;
t1748 = t1690 * t1691;
t1597 = t1748 - t1751;
t1749 = t1690 * t1687;
t1598 = t1749 + t1750;
t1561 = t1597 * t1688 + t1598 * t1692;
t1724 = t1597 * t1692 - t1598 * t1688;
t1803 = -t1561 * t1683 + t1684 * t1724;
t1802 = t1561 * t1684 + t1683 * t1724;
t1801 = m(4) + m(8) + m(9) + t1792;
t1679 = qJ(3) + qJ(2);
t1673 = qJ(4) + t1679;
t1651 = sin(t1673);
t1674 = -qJ(4) + t1679;
t1652 = sin(t1674);
t1654 = cos(t1673);
t1655 = cos(t1674);
t1663 = sin(t1679);
t1667 = cos(t1679);
t1781 = rSges(5,3) * m(5);
t1782 = rSges(4,3) * m(4);
t1785 = m(9) * rSges(9,3);
t1790 = m(6) * rSges(6,1);
t1795 = m(6) / 0.2e1;
t1548 = (rSges(4,1) * t1782 + rSges(9,1) * t1785 - Icges(4,5) - Icges(9,5)) * t1667 - t1663 * (rSges(4,2) * t1782 + rSges(9,2) * t1785 - Icges(4,6) - Icges(9,6)) + (t1781 * t1667 - t1652 * t1790 / 0.2e1 + (rSges(6,1) * t1651 + (t1654 + t1655) * rSges(6,2)) * t1795) * pkin(4);
t1796 = -pkin(8) / 0.2e1;
t1794 = pkin(1) * m(8);
t1793 = pkin(3) * m(9);
t1791 = m(4) * rSges(4,2);
t1789 = m(6) * rSges(6,2);
t1788 = m(7) * rSges(7,2);
t1787 = m(7) * rSges(7,3);
t1786 = m(9) * rSges(9,2);
t1784 = rSges(3,2) * m(3);
t1783 = rSges(6,2) * pkin(8);
t1685 = sin(qJ(4));
t1689 = cos(qJ(4));
t1623 = rSges(6,1) * t1689 - rSges(6,2) * t1685;
t1780 = t1623 / 0.2e1;
t1637 = qJ(4) + t1720;
t1631 = qJ(2) + t1637;
t1619 = qJ(3) + t1631;
t1779 = cos(t1619);
t1638 = -qJ(4) + t1720;
t1634 = -qJ(2) + t1638;
t1622 = -qJ(3) + t1634;
t1778 = cos(t1622);
t1777 = cos(t1631);
t1776 = cos(t1634);
t1775 = sin(t1619);
t1774 = sin(t1622);
t1773 = sin(t1631);
t1772 = sin(t1634);
t1603 = pkin(4) * t1792 + m(4) * rSges(4,1) + m(9) * rSges(9,1);
t1771 = pkin(1) * t1603;
t1650 = t1786 + t1791;
t1769 = pkin(1) * t1650;
t1765 = m(6) * t1623;
t1764 = m(6) * (rSges(6,1) * t1685 + rSges(6,2) * t1689);
t1661 = sin(t1678);
t1662 = cos(t1678);
t1763 = (-t1661 * t1802 + t1662 * t1803) * pkin(4);
t1693 = -rSges(6,3) - pkin(10);
t1762 = rSges(6,1) * t1693;
t1761 = rSges(6,2) * t1693;
t1760 = pkin(18) + pkin(15);
t1632 = -qJ(2) + t1637;
t1613 = sin(t1632);
t1633 = qJ(2) + t1638;
t1614 = sin(t1633);
t1615 = cos(t1632);
t1616 = cos(t1633);
t1641 = qJ(2) + t1720;
t1642 = -qJ(2) + t1720;
t1643 = -m(5) * rSges(5,2) - m(6) * t1693;
t1675 = pkin(14) - qJ(2) - pkin(15);
t1646 = 0.2e1 * t1675;
t1653 = sin(t1675);
t1656 = cos(t1675);
t1657 = pkin(8) * m(6) + m(5) * rSges(5,1);
t1701 = 0.2e1 * qJ(2);
t1680 = t1701 + qJ(3);
t1664 = sin(t1680);
t1668 = cos(t1680);
t1671 = qJ(2) + t1760;
t1672 = -qJ(2) + t1760;
t1635 = qJ(3) + t1641;
t1636 = -qJ(3) + t1642;
t1722 = qJ(2) + t1736;
t1644 = qJ(3) + t1722;
t1723 = -qJ(2) + t1736;
t1645 = -qJ(3) + t1723;
t1659 = 0.2e1 * t1679;
t1620 = qJ(3) + t1633;
t1621 = -qJ(3) + t1632;
t1740 = -t1789 / 0.4e1;
t1728 = pkin(4) * t1740;
t1741 = t1790 / 0.4e1;
t1729 = pkin(4) * t1741;
t1742 = -t1790 / 0.4e1;
t1730 = pkin(4) * t1742;
t1710 = sin(t1620) * t1729 + sin(t1621) * t1730 + (cos(t1620) + cos(t1621)) * t1728;
t1744 = t1793 / 0.2e1;
t1745 = -t1793 / 0.2e1;
t1706 = -t1710 + rSges(9,1) * sin(t1645) * t1744 + t1775 * t1730 + t1774 * t1729 + (-rSges(4,1) * t1791 - rSges(9,1) * t1786 + Icges(4,4) + Icges(9,4)) * cos(t1659) - (-Icges(4,1) - Icges(9,1) + Icges(4,2) + Icges(9,2) + t1792 * pkin(4) ^ 2 + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t1659) / 0.2e1 + (-cos(t1635) / 0.2e1 + cos(t1636) / 0.2e1) * pkin(4) * t1643 + (-sin(t1635) / 0.2e1 + sin(t1636) / 0.2e1) * pkin(4) * t1657 + (t1779 + t1778) * t1728 + (t1603 * t1663 + t1650 * t1667) * pkin(12) + (rSges(9,1) * sin(t1644) + (cos(t1644) + cos(t1645)) * rSges(9,2)) * t1745;
t1746 = t1794 / 0.2e1;
t1747 = -t1794 / 0.2e1;
t1538 = (-rSges(3,1) * t1784 + Icges(3,4)) * cos(t1701) - (rSges(7,1) * t1788 - Icges(7,4)) * cos(t1646) + pkin(1) * sin(t1722) * t1744 + pkin(1) * sin(t1723) * t1745 + ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t1646) / 0.2e1 + t1668 * t1769 + t1664 * t1771 + t1706 - (-Icges(3,1) + Icges(3,2) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) + t1801 * pkin(1) ^ 2) * sin(t1701) / 0.2e1 + (cos(t1641) / 0.2e1 - cos(t1642) / 0.2e1) * pkin(1) * t1643 + (sin(t1641) / 0.2e1 - sin(t1642) / 0.2e1) * pkin(1) * t1657 + (t1772 + t1613) * pkin(1) * t1742 + (t1773 + t1614) * pkin(1) * t1741 + (t1616 + t1615) * pkin(1) * t1740 + (t1776 + t1777) * pkin(1) * t1789 / 0.4e1 + (-(pkin(1) * t1801 + rSges(3,1) * m(3)) * t1687 - t1691 * t1784) * pkin(12) + (cos(t1671) * t1747 + cos(t1672) * t1746) * rSges(8,2) + (sin(t1671) * t1746 + sin(t1672) * t1747) * rSges(8,1) + (-m(7) * rSges(7,1) * t1653 + t1656 * t1788) * pkin(6);
t1759 = t1538 * qJD(1);
t1617 = 0.2e1 * t1637;
t1618 = 0.2e1 * t1638;
t1719 = 0.2e1 * t1720;
t1639 = qJ(4) + t1719;
t1640 = -qJ(4) + t1719;
t1649 = 0.4e1 * t1762;
t1694 = rSges(6,1) * pkin(8);
t1700 = 0.2e1 * qJ(4);
t1703 = 4 * Icges(6,5);
t1539 = (cos(t1700) / 0.2e1 - cos(t1618) / 0.4e1 - cos(t1617) / 0.4e1) * (rSges(6,1) * t1789 - Icges(6,4)) + (t1689 * t1796 + (cos(t1637) / 0.2e1 + cos(t1638) / 0.2e1) * pkin(12) + (-t1778 / 0.4e1 - t1779 / 0.4e1) * pkin(4) + (t1777 / 0.4e1 + t1776 / 0.4e1 + t1615 / 0.4e1 + t1616 / 0.4e1) * pkin(1)) * t1789 + ((t1649 - 0.4e1 * t1783) * m(6) + t1703) * cos(t1639) / 0.16e2 - ((t1649 + 0.4e1 * t1783) * m(6) + t1703) * cos(t1640) / 0.16e2 + (t1685 * t1796 + (sin(t1637) / 0.2e1 - sin(t1638) / 0.2e1) * pkin(12) + (-t1775 / 0.4e1 + t1774 / 0.4e1) * pkin(4) + (t1773 / 0.4e1 - t1772 / 0.4e1 + t1613 / 0.4e1 - t1614 / 0.4e1) * pkin(1)) * t1790 + t1710 + (sin(t1700) / 0.4e1 - sin(t1617) / 0.8e1 + sin(t1618) / 0.8e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2)) - ((t1694 + t1761) * m(6) + Icges(6,6)) * sin(t1639) / 0.4e1 + ((t1694 - t1761) * m(6) - Icges(6,6)) * sin(t1640) / 0.4e1;
t1758 = t1539 * qJD(1);
t1540 = t1706 + (t1664 + t1686) * t1771 / 0.2e1 + (t1668 + t1690) * t1769 / 0.2e1;
t1757 = t1540 * qJD(1);
t1681 = qJ(2) + qJ(4);
t1665 = sin(t1681);
t1682 = qJ(2) - qJ(4);
t1666 = sin(t1682);
t1669 = cos(t1681);
t1670 = cos(t1682);
t1714 = (-t1655 / 0.4e1 + t1654 / 0.4e1) * rSges(6,2) + (t1652 / 0.4e1 + t1651 / 0.4e1) * rSges(6,1);
t1713 = t1714 * pkin(4);
t1709 = t1713 + ((t1670 / 0.4e1 - t1669 / 0.4e1) * rSges(6,2) + (-t1666 / 0.4e1 - t1665 / 0.4e1) * rSges(6,1)) * pkin(1);
t1753 = t1623 * t1577;
t1546 = (-t1753 / 0.2e1 + t1709) * m(6);
t1756 = t1546 * qJD(1);
t1755 = t1548 * qJD(3);
t1549 = (t1648 * t1780 + (t1749 * t1780 + t1714) * pkin(4)) * m(6);
t1754 = t1549 * qJD(1);
t1739 = qJD(4) * t1765;
t1726 = t1763 * t1764;
t1718 = -pkin(4) * t1748 + t1647;
t1717 = pkin(4) * t1749 + t1648;
t1711 = -t1688 * t1718 + t1692 * t1717;
t1712 = t1688 * t1717 + t1692 * t1718;
t1707 = (-t1683 * t1712 + t1684 * t1711) * t1661;
t1708 = (t1683 * t1711 + t1684 * t1712) * t1662;
t1541 = (t1763 / 0.2e1 + t1708 / 0.2e1 + t1707 / 0.2e1) * t1764;
t1566 = (t1603 * t1686 + t1650 * t1690) * pkin(1);
t1715 = -qJD(2) * t1566 + qJD(4) * t1541;
t1596 = pkin(8) * t1764;
t1568 = -t1685 * (m(6) * t1762 + Icges(6,5)) - (m(6) * t1761 + Icges(6,6)) * t1689;
t1564 = t1566 * qJD(3);
t1552 = t1568 * t1692 + t1596 * t1688;
t1551 = t1568 * t1688 - t1596 * t1692;
t1550 = -t1717 * t1765 / 0.2e1 + m(6) * t1713;
t1547 = m(6) * t1709 + t1753 * t1795;
t1542 = -(t1708 + t1707) * t1764 / 0.2e1 + t1726 / 0.2e1;
t1 = [qJD(2) * t1538 + qJD(3) * t1540 + qJD(4) * t1539, t1547 * qJD(4) + t1755 + t1759 + (Icges(3,5) * t1691 - t1687 * Icges(3,6) + (-rSges(3,1) * t1691 + rSges(3,2) * t1687) * rSges(3,3) * m(3) + (-rSges(7,2) * t1787 + Icges(7,6)) * t1653 - (rSges(7,1) * t1787 - Icges(7,5)) * t1656 + ((-m(8) * rSges(8,3) - t1781 - t1782 - t1785) * t1691 + ((-t1670 / 0.2e1 - t1669 / 0.2e1) * rSges(6,2) + (t1666 / 0.2e1 - t1665 / 0.2e1) * rSges(6,1)) * m(6)) * pkin(1) + t1548) * qJD(2), qJD(2) * t1548 + qJD(4) * t1550 + t1755 + t1757, t1758 + t1547 * qJD(2) + t1550 * qJD(3) + ((t1551 * t1684 + t1552 * t1683) * t1662 + (-t1551 * t1683 + t1552 * t1684) * t1661 + (pkin(12) + t1578) * t1764) * qJD(4); qJD(4) * t1546 - t1759, t1564, t1564 - t1715, t1756 - t1541 * qJD(3) - (t1661 * t1805 + t1662 * t1804) * t1739; qJD(4) * t1549 - t1757, t1715, 0, t1754 + t1541 * qJD(2) + (t1661 * t1803 + t1662 * t1802) * pkin(4) * t1739; -qJD(2) * t1546 - qJD(3) * t1549 - t1758, -t1756 - (-t1661 * t1804 + t1662 * t1805) * qJD(2) * t1764 + t1542 * qJD(3), qJD(2) * t1542 + qJD(3) * t1726 - t1754, 0;];
Cq = t1;
