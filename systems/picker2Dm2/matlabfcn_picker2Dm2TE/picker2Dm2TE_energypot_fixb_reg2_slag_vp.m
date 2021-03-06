% Calculate inertial parameters regressor of potential energy for
% picker2Dm2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% U_reg [1x(2*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = picker2Dm2TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 11:52:40
% EndTime: 2020-05-09 11:52:46
% DurationCPUTime: 6.28s
% Computational Cost: add. (40019->430), mult. (117731->572), div. (838->5), fcn. (20852->10), ass. (0->231)
t707 = sin(qJ(1));
t708 = cos(qJ(2));
t872 = t707 * t708;
t827 = pkin(3) * t872;
t799 = pkin(1) * t827;
t644 = -0.2e1 * t799;
t734 = 0.2e1 * pkin(3);
t746 = pkin(1) ^ 2;
t744 = t746 ^ 2;
t748 = pkin(7) ^ 2;
t709 = cos(qJ(1));
t676 = t709 ^ 2;
t706 = sin(qJ(2));
t662 = pkin(3) * t706;
t888 = 2 * pkin(7);
t652 = t662 * t888;
t673 = t706 ^ 2;
t741 = pkin(3) ^ 2;
t873 = t673 * t741;
t831 = 0.2e1 * t873;
t845 = -t741 + t748;
t778 = t652 + t831 + t845;
t775 = t778 * t676;
t668 = t746 + t748;
t736 = pkin(4) ^ 2;
t809 = -t736 + t668;
t792 = t652 + t809;
t883 = t662 + pkin(7);
t806 = t883 * t709;
t846 = t736 - t748;
t808 = -0.2e1 * t741 + t846;
t834 = -0.4e1 * t662;
t843 = t746 - t748;
t884 = 0.2e1 * t746;
t895 = 0.4e1 * pkin(1);
t621 = sqrt(-0.4e1 * t746 * t775 + 0.4e1 * t843 * t873 + pkin(7) * t809 * t834 - t744 + t808 * t884 - (t748 - (t734 + pkin(4)) * pkin(4)) * (t748 + (t734 - pkin(4)) * pkin(4)) + (-(t644 + t792) * t806 + t792 * t827) * t895);
t807 = t741 + t668;
t772 = 0.1e1 / (0.2e1 * pkin(1) * t806 + t644 + t652 + t807);
t788 = -pkin(1) + t827;
t725 = 0.3e1 * t746;
t789 = t725 - t808;
t794 = -0.4e1 * t799;
t871 = t708 * t709;
t826 = pkin(3) * t871;
t886 = 0.4e1 * t741;
t768 = t772 * ((t707 * t883 + t826) * t621 - (t652 + t789 + t794) * t806 + t788 * t652 + t789 * t827 + (-0.2e1 * t775 + t831 - t886 - t809) * pkin(1));
t765 = t768 / 0.4e1;
t663 = pkin(3) * t708;
t887 = 0.2e1 * t741;
t782 = t887 + t792;
t664 = pkin(1) * t709;
t837 = 0.2e1 * t664;
t771 = t772 * ((t806 - t788) * t621 + (t778 * t837 + t883 * t782) * t707 + (t782 * t709 + (0.4e1 * t676 - 0.2e1) * pkin(1) * t883) * t663);
t737 = 0.1e1 / pkin(4);
t867 = t737 / pkin(3) ^ 2;
t753 = t741 ^ 2;
t735 = t736 ^ 2;
t747 = t748 ^ 2;
t844 = t744 + t747;
t732 = 0.2e1 * t748;
t849 = t732 - t736;
t864 = t748 * t736;
t776 = t849 * t746 + t735 / 0.6e1 + t844 - t864;
t630 = -t753 / 0.6e1 + t776;
t697 = -t741 / 0.3e1;
t658 = t697 + t748;
t633 = t658 * t644;
t651 = t664 + pkin(7);
t639 = t662 + t651;
t667 = -0.3e1 * t741 + t748;
t675 = t709 * t676;
t749 = pkin(1) * t746;
t862 = t749 * t675;
t833 = pkin(7) * t862;
t804 = 0.8e1 * t833;
t642 = t667 * t804;
t666 = -t736 - t741;
t731 = 0.3e1 * t748;
t654 = t731 + t666;
t878 = t654 * t746;
t643 = 0.10e2 * t878;
t693 = 0.4e1 / 0.3e1 * t741;
t687 = -t736 / 0.3e1;
t814 = t687 + t668;
t645 = t693 + t814;
t688 = -t736 / 0.2e1;
t647 = t688 + t807;
t648 = -t736 + t807;
t653 = pkin(7) * t837;
t730 = 0.4e1 * t748;
t656 = (t730 + t736) * t746;
t659 = -t746 / 0.3e1 + t748;
t660 = 0.10e2 / 0.3e1 * t746;
t661 = t668 ^ 2;
t665 = -0.30e2 * t736 + 0.60e2 * t748;
t670 = -0.3e1 * t746 + t748;
t685 = -t736 / 0.6e1;
t686 = -t736 / 0.4e1;
t694 = 0.2e1 / 0.3e1 * t741;
t699 = 0.4e1 / 0.3e1 * t746;
t701 = t746 / 0.2e1;
t711 = 0.15e2 * t744;
t712 = 0.15e2 * t746;
t713 = 0.10e2 * t746;
t718 = -0.2e1 * t736;
t719 = -0.5e1 * t736;
t720 = 0.5e1 * t753;
t721 = 0.7e1 * t744;
t722 = 0.5e1 * t744;
t723 = 0.7e1 * t746;
t724 = 0.6e1 * t746;
t728 = 0.3e1 * t747;
t729 = 0.8e1 * t748;
t752 = pkin(3) * t741;
t738 = t752 ^ 2;
t757 = pkin(7) * t748;
t773 = 0.5e1 / 0.6e1 * t753 + t776;
t781 = t748 - t799;
t857 = t735 / 0.2e1 - t753 / 0.2e1;
t791 = -0.3e1 * t864 + t728 + t857;
t795 = -0.6e1 * t799;
t690 = -0.3e1 / 0.2e1 * t736;
t856 = t690 + t731;
t859 = t668 * ((t690 + t732) * t746 - 0.3e1 / 0.2e1 * t864 + t844 + t857) + t738;
t783 = ((t660 + t849) * t741 + t773) * t795 + (t711 + (-0.9e1 * t736 + 0.18e2 * t748) * t746 + t791) * t741 + (t712 + t856) * t753 + t859;
t784 = t647 * t794;
t850 = t725 + t748;
t811 = t741 + t850;
t881 = pkin(1) * t707;
t786 = -(0.3e1 * t741 + t668) * t881 + t811 * t663;
t689 = -0.2e1 / 0.3e1 * t736;
t698 = -0.2e1 / 0.3e1 * t741;
t812 = t689 + t668;
t851 = t713 + t732;
t855 = t698 + t748;
t787 = -(t720 + (0.5e1 * t746 + t654) * t887 + (t698 + t812) * t668) * t881 + (t753 + (t689 + t698 + t851) * t741 + t722 + 0.2e1 * t878 + t748 * (t689 + t855)) * t663;
t848 = t735 - t753;
t790 = -0.6e1 * t864 + 0.6e1 * t747 + t848;
t813 = t689 + t694 + t732;
t858 = (t694 + t812) * t668 + t753;
t793 = t645 * t794 + (t724 + t813) * t741 + t858;
t796 = t862 * t663;
t866 = t744 * t676 ^ 2;
t797 = t866 * t663;
t824 = 0.16e2 * t862;
t802 = pkin(7) * t824;
t803 = 0.20e2 / 0.3e1 * t746;
t847 = -t736 + t741;
t810 = t731 + t847;
t695 = t741 / 0.3e1;
t815 = t685 + t695 + t748;
t816 = t736 / 0.3e1 + t695 + t732;
t817 = 0.2e1 / 0.3e1 * t736 + t694 + t730;
t818 = 0.4e1 / 0.3e1 * t736 + t693 - 0.2e1 * t748;
t870 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t819 = t707 * t870;
t835 = 0.6e1 * t664;
t820 = pkin(7) * t835;
t836 = 0.4e1 * t664;
t821 = pkin(7) * t836;
t822 = -t881 / 0.2e1;
t865 = t746 * t676;
t825 = 0.12e2 * t865;
t828 = t746 * t663;
t829 = 0.4e1 * t865;
t830 = 0.8e1 * t866;
t874 = t706 * t673 * t752;
t832 = -0.8e1 * t874;
t838 = 0.2e1 * t881;
t839 = pkin(7) * t664;
t840 = 4 * pkin(7);
t841 = t747 + t753;
t842 = t747 - t744;
t852 = 0.4e1 / 0.7e1 * t748 - t736 / 0.7e1;
t853 = t701 + t748;
t854 = t746 / 0.3e1 + t748;
t863 = t748 * t746;
t869 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t875 = t673 ^ 2 * t753;
t876 = t666 * t748;
t877 = t661 * (-t741 + t809);
t879 = (-t707 * t749 + t828) * t676;
t885 = 0.4e1 * t744;
t889 = t686 + t741 / 0.2e1;
t880 = ((-0.24e2 * (0.4e1 / 0.3e1 * t865 + t653 + t659) * t875 * t881 - 0.12e2 * (-0.8e1 / 0.3e1 * t797 + ((t699 + t815) * t663 - (0.7e1 / 0.6e1 * t741 + t685 + t853) * t881) * t829 + (-t741 * t843 - 0.5e1 / 0.3e1 * t744 + t816 * t746 + t748 * (t687 + t658)) * t663 + (-t753 + (-t803 + t817) * t741 - 0.3e1 * t744 + t818 * t746 + t747) * t822 + (-t707 * t744 * t675 + ((t746 + t815) * t663 + (t887 - t843) * t822) * t664) * t840) * t873 + 0.24e2 * t658 * t797 + ((t748 + 0.5e1 / 0.2e1 * t741 + 0.3e1 / 0.2e1 * t746 + t688) * t663 + t667 * t881 / 0.2e1) * t802 - 0.6e1 * ((-0.3e1 * t753 + (-t803 + t818) * t741 + t817 * t746 + t842) * t663 - 0.2e1 * (-0.5e1 / 0.3e1 * t753 + (-t746 + t816) * t741 + t748 * (t697 + t814)) * t881) * t865 - 0.6e1 * t787 * t839 - (t738 + (0.21e2 * t746 + t654) * t753 + (t643 + t728 + 0.35e2 * t744 + 0.2e1 * t876) * t741 + (t721 + (t719 + t729 - 0.5e1 * t741) * t746 + t748 * (-t736 + t845)) * t668) * t663 + (0.7e1 * t738 + (t723 + t654) * t720 + (t643 + 0.21e2 * t744 + 0.9e1 * t747 + 0.6e1 * t876) * t741 + t877) * t881) * t621 + (0.16e2 * (t830 + t802 + (-0.8e1 * t744 + 0.12e2 * t863) * t676 + (-0.12e2 * pkin(7) * t749 + t757 * t895) * t709 - 0.6e1 * t863 + t844) * t875 + 0.24e2 * (t855 * t830 + 0.14e2 * (-0.32e2 / 0.21e2 * (t748 + t741 / 0.4e1 + t746 / 0.4e1 - t736 / 0.8e1) * t799 + 0.5e1 / 0.42e2 * t753 + (0.16e2 / 0.21e2 * t746 + t852) * t741 + t744 / 0.7e1 + t852 * t746 + t747 - 0.3e1 / 0.7e1 * t864 + t735 / 0.42e2) * t865 + t659 * t784 - t843 * t753 + (t656 - 0.10e2 / 0.3e1 * t744 + 0.2e1 * t747 - t864) * t741 + t630 * t869 + ((-0.2e1 / 0.3e1 * t799 + t686 + t853) * t824 + (-0.8e1 / 0.3e1 * (t854 + t889) * t799 + 0.5e1 / 0.18e2 * t753 + (0.4e1 / 0.3e1 * t748 + t699 + t687) * t741 + t747 + 0.2e1 / 0.3e1 * t863 - 0.2e1 / 0.3e1 * t864 - t744 / 0.3e1 + t735 / 0.18e2) * t835) * pkin(7)) * t873 + 0.16e2 * (-0.6e1 * t748 * t741 + t841) * t866 + 0.32e2 * (t644 * t870 + t647 * t667) * t833 + 0.24e2 * (t658 * t784 - t738 + (-t660 + t846) * t753 + (t656 + t753 / 0.6e1 - t735 / 0.6e1 + t842) * t741 + t630 * t748) * t865 + 0.8e1 * t783 * t839 - 0.8e1 * ((t723 + t856) * t753 + (t721 + (t719 + 0.10e2 * t748) * t746 + t791) * t741 + t859) * t799 + t753 ^ 2 + (t718 + t730 + 0.28e2 * t746) * t738 + (t665 * t746 + 0.70e2 * t744 + t790) * t753 + (t665 * t744 + t790 * t724 + t848 * t732 - 0.6e1 * t747 * t736 + 0.28e2 * t749 ^ 2 + 0.4e1 * t757 ^ 2) * t741 + t648 * t877) * t639 + (((0.4e1 * t879 + (t663 + t838) * t653 + t670 * t663 + (t688 + t811) * t838) * t832 - 0.6e1 * (-0.4e1 * (0.2e1 * (0.5e1 / 0.6e1 * t741 + t701 + t685) * t663 + pkin(1) * t819) * t865 + (-0.8e1 * t796 + ((t687 + t694 + t850) * t663 - (0.8e1 / 0.3e1 * t741 + t814) * t881) * t836) * pkin(7) + t787) * t662) * t621 + (0.32e2 * (t804 + (-0.4e1 * t749 * t827 + t885 + (t886 + t718 + t729) * t746) * t676 + (-t746 + t781 + t889) * t821 + t644 * t869 + t670 * t647) * t874 + 0.8e1 * (t642 + (t647 * t870 + t633) * t825 + (t784 + (t724 + t849) * t741 + t773) * t820 + t783) * t662) * t639) * t651) / ((-0.4e1 * (-t843 * t663 + 0.2e1 * t879 + (t826 * t888 + t707 * (t741 + t884)) * pkin(1)) * t873 + 0.8e1 * pkin(7) * t796 + ((pkin(3) * t885 + 0.8e1 * t746 * t752) * t708 + 0.4e1 * t749 * t819) * t676 - 0.4e1 * t786 * t839 - (t851 * t741 + t722 + t841 + 0.6e1 * t863) * t663 + (t720 + (t713 + 0.6e1 * t748) * t741 + t661) * t881) * t621 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t799 + 0.4e1 / 0.9e1 * t741 - t736 / 0.9e1 + t854) * t865 + t659 * t644 + t645 * t869 + (t862 + (t685 + t694 + t781) * t664) * t840) * t873 + t642 + (t645 * t870 + t633) * t825 + t793 * t820 + ((t660 + t813) * t741 + t858) * t795 + t738 + (t712 + t810) * t753 + (t810 * t724 + t847 * t732 + t711 + t728) * t741 + t661 * t648) * t639 + ((t832 * t881 + (-0.2e1 * t676 * t828 + (t663 - t881) * t653 + t786) * t834) * t621 + (0.8e1 * (t653 + t829 + t670) * t874 + 0.6e1 * (t845 * t829 + (t644 + t645) * t821 + t793) * t662) * t639) * t651);
t603 = (t765 * t880 - t621 * t771 / 0.4e1) * t867;
t604 = (t621 * t765 + t771 * t880 / 0.4e1) * t867;
t634 = -t706 * t707 - t871;
t635 = -t706 * t709 + t872;
t593 = t603 * t634 + t604 * t635;
t594 = t603 * t635 - t604 * t634;
t742 = 0.1e1 / pkin(3);
t767 = -t768 / 0.2e1;
t770 = -t771 / 0.2e1;
t617 = (t707 * t767 + t709 * t770) * t742;
t769 = t771 / 0.2e1;
t618 = (t707 * t769 + t709 * t767) * t742;
t805 = t880 / 0.2e1;
t868 = t737 * t742;
t777 = (-t617 * t621 / 0.2e1 + t618 * t805) * t868;
t890 = (t617 * t805 + t618 * t621 / 0.2e1) * t868;
t899 = t593 * t890 + t594 * t777;
t898 = -t593 * t777 + t594 * t890;
t703 = sin(pkin(8));
t704 = cos(pkin(8));
t631 = -t703 * t709 + t704 * t707;
t632 = t703 * t707 + t704 * t709;
t766 = t768 / 0.2e1;
t612 = (t631 * t766 + t632 * t769) * t742;
t710 = cos(pkin(9));
t764 = t742 * (t631 * t770 + t632 * t766);
t882 = sin(pkin(9));
t608 = t612 * t710 - t882 * t764;
t609 = t882 * t612 + t710 * t764;
t600 = t608 * t618 - t609 * t617;
t798 = -t608 * t617 - t618 * t609;
t897 = t600 * t608 - t609 * t798;
t896 = t600 * t609 + t608 * t798;
t861 = t618 * pkin(3) - t664;
t860 = t618 * pkin(2) - t664;
t801 = t617 * pkin(3) - t881;
t800 = t617 * pkin(2) - t881;
t785 = g(1) * t709 + g(2) * t707;
t780 = t631 * t709 - t632 * t707;
t779 = t631 * t707 + t632 * t709;
t637 = t785 * pkin(1);
t625 = t631 * t703 - t632 * t704;
t624 = -t631 * t704 - t632 * t703;
t592 = -g(1) * t798 - g(2) * t600;
t591 = g(1) * t600 - g(2) * t798;
t1 = [0, 0, 0, 0, 0, 0, t785, -g(1) * t707 + g(2) * t709, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(1) * t618 - g(2) * t617, g(1) * t617 - g(2) * t618, -g(3), t637, 0, 0, 0, 0, 0, 0, t592, t591, -g(3), -g(1) * t860 - g(2) * t800, 0, 0, 0, 0, 0, 0, -g(1) * t777 - g(2) * t890, g(1) * t890 - g(2) * t777, -g(3), -g(1) * t861 - g(2) * t801, 0, 0, 0, 0, 0, 0, -g(1) * t625 - g(2) * t624, g(1) * t624 - g(2) * t625, -g(3), (-g(1) * t704 - g(2) * t703) * pkin(5), 0, 0, 0, 0, 0, 0, t592, t591, -g(3), t637, 0, 0, 0, 0, 0, 0, -g(1) * t706 + g(2) * t708, -g(1) * t708 - g(2) * t706, -g(3), -g(1) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t779 + g(2) * t780, -g(1) * t780 - g(2) * t779, -g(3), t637, 0, 0, 0, 0, 0, 0, g(1) * t897 - g(2) * t896, g(1) * t896 + g(2) * t897, -g(3), -g(2) * (pkin(6) * t600 + t800) - g(1) * (pkin(6) * t798 + t860), 0, 0, 0, 0, 0, 0, -g(1) * t899 - g(2) * t898, g(1) * t898 - g(2) * t899, -g(3), -g(2) * (pkin(4) * t890 + t801) - g(1) * (pkin(4) * t777 + t861);];
U_reg = t1;
