% Calculate inertial parameters regressor of potential energy for
% picker2Dm2DE2
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
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = picker2Dm2DE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 20:38:17
% EndTime: 2020-05-09 20:38:23
% DurationCPUTime: 5.32s
% Computational Cost: add. (31067->425), mult. (92567->539), div. (594->5), fcn. (15660->31), ass. (0->235)
t884 = 4 * pkin(1);
t746 = pkin(4) ^ 2;
t694 = -t746 / 0.4e1;
t751 = pkin(3) ^ 2;
t883 = t694 + t751 / 0.2e1;
t882 = 2 * pkin(7);
t718 = cos(qJ(1));
t684 = t718 ^ 2;
t881 = -0.2e1 * t684;
t730 = 0.2e1 * t751;
t880 = 0.4e1 * t751;
t756 = pkin(1) ^ 2;
t754 = t756 ^ 2;
t879 = 4 * t754;
t878 = 2 * t756;
t734 = 6 * t756;
t758 = pkin(7) ^ 2;
t742 = 2 * t758;
t763 = t751 ^ 2;
t729 = 0.5e1 * t763;
t711 = sin(pkin(8));
t712 = cos(pkin(8));
t715 = sin(qJ(1));
t638 = t711 * t715 + t712 * t718;
t877 = t638 / 0.2e1;
t672 = pkin(1) * t718;
t659 = t672 + pkin(7);
t876 = pkin(1) * t715;
t714 = sin(qJ(2));
t670 = pkin(3) * t714;
t717 = cos(qJ(2));
t671 = pkin(3) * t717;
t660 = t670 * t882;
t681 = t714 ^ 2;
t866 = t681 * t751;
t828 = 0.2e1 * t866;
t842 = -t751 + t758;
t639 = t660 + t828 + t842;
t865 = t715 * t717;
t824 = pkin(3) * t865;
t794 = pkin(1) * t824;
t652 = -0.2e1 * t794;
t744 = 0.2e1 * pkin(3);
t676 = t756 + t758;
t804 = -t746 + t676;
t788 = t660 + t804;
t860 = t756 * t684;
t826 = -0.4e1 * t860;
t831 = -0.4e1 * t670;
t840 = t756 - t758;
t843 = t746 - t758;
t658 = t670 + pkin(7);
t872 = t658 * t718;
t622 = sqrt(t639 * t826 + 0.4e1 * t840 * t866 + pkin(7) * t804 * t831 - t754 + (-0.2e1 * t751 + t843) * t878 - (t758 - (t744 + pkin(4)) * pkin(4)) * (t758 + (t744 - pkin(4)) * pkin(4)) + (-(t652 + t788) * t872 + t788 * t824) * t884);
t745 = t746 ^ 2;
t757 = t758 ^ 2;
t841 = t754 + t757;
t846 = t742 - t746;
t859 = t758 * t746;
t775 = t846 * t756 + t745 / 0.6e1 + t841 - t859;
t636 = -t763 / 0.6e1 + t775;
t705 = -t751 / 0.3e1;
t666 = t705 + t758;
t640 = t666 * t652;
t647 = t670 + t659;
t675 = -0.3e1 * t751 + t758;
t683 = t718 * t684;
t759 = pkin(1) * t756;
t857 = t759 * t683;
t830 = pkin(7) * t857;
t799 = 0.8e1 * t830;
t650 = t675 * t799;
t674 = -t746 - t751;
t741 = 3 * t758;
t662 = t741 + t674;
t871 = t662 * t756;
t651 = 0.10e2 * t871;
t701 = 0.4e1 / 0.3e1 * t751;
t695 = -t746 / 0.3e1;
t809 = t695 + t676;
t653 = t701 + t809;
t696 = -t746 / 0.2e1;
t803 = t751 + t676;
t655 = t696 + t803;
t656 = -t746 + t803;
t834 = 0.2e1 * t672;
t661 = pkin(7) * t834;
t740 = 4 * t758;
t664 = (t740 + t746) * t756;
t667 = -t756 / 0.3e1 + t758;
t668 = 0.10e2 / 0.3e1 * t756;
t669 = t676 ^ 2;
t673 = -0.30e2 * t746 + (60 * t758);
t678 = -3 * t756 + t758;
t693 = -t746 / 0.6e1;
t702 = 0.2e1 / 0.3e1 * t751;
t707 = 0.4e1 / 0.3e1 * t756;
t709 = t756 / 0.2e1;
t720 = 15 * t754;
t721 = 15 * t756;
t722 = 10 * t756;
t727 = -0.2e1 * t746;
t728 = -0.5e1 * t746;
t731 = 7 * t754;
t732 = 5 * t754;
t733 = 7 * t756;
t738 = 3 * t757;
t739 = 8 * t758;
t762 = pkin(3) * t751;
t748 = t762 ^ 2;
t767 = pkin(7) * t758;
t774 = 0.5e1 / 0.6e1 * t763 + t775;
t776 = t758 - t794;
t854 = t745 / 0.2e1 - t763 / 0.2e1;
t787 = -0.3e1 * t859 + t738 + t854;
t791 = -0.6e1 * t794;
t698 = -0.3e1 / 0.2e1 * t746;
t853 = t698 + t741;
t856 = t676 * ((t698 + t742) * t756 - 0.3e1 / 0.2e1 * t859 + t841 + t854) + t748;
t777 = ((t668 + t846) * t751 + t774) * t791 + (t720 + (-0.9e1 * t746 + (18 * t758)) * t756 + t787) * t751 + (t721 + t853) * t763 + t856;
t790 = -0.4e1 * t794;
t778 = t655 * t790;
t735 = 3 * t756;
t847 = t735 + t758;
t806 = t751 + t847;
t782 = -(0.3e1 * t751 + t676) * t876 + t806 * t671;
t697 = -0.2e1 / 0.3e1 * t746;
t706 = -0.2e1 / 0.3e1 * t751;
t807 = t697 + t676;
t848 = t722 + t742;
t852 = t706 + t758;
t783 = -(t729 + ((5 * t756) + t662) * t730 + (t706 + t807) * t676) * t876 + (t763 + (t697 + t706 + t848) * t751 + t732 + 0.2e1 * t871 + t758 * (t697 + t852)) * t671;
t845 = t745 - t763;
t786 = -0.6e1 * t859 + (6 * t757) + t845;
t808 = t697 + t702 + t742;
t855 = (t702 + t807) * t676 + t763;
t789 = t653 * t790 + (t734 + t808) * t751 + t855;
t821 = t759 * t671;
t792 = t683 * t821;
t861 = t754 * t684 ^ 2;
t793 = t861 * t671;
t819 = 0.16e2 * t857;
t797 = pkin(7) * t819;
t798 = 0.20e2 / 0.3e1 * t756;
t844 = -t746 + t751;
t805 = t741 + t844;
t703 = t751 / 0.3e1;
t810 = t693 + t703 + t758;
t811 = t746 / 0.3e1 + t703 + t742;
t812 = 0.2e1 / 0.3e1 * t746 + t702 + t740;
t813 = 0.4e1 / 0.3e1 * t746 + t701 - (2 * t758);
t863 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t814 = t715 * t863;
t832 = 0.6e1 * t672;
t816 = pkin(7) * t832;
t833 = 0.4e1 * t672;
t817 = pkin(7) * t833;
t818 = -t876 / 0.2e1;
t820 = 0.12e2 * t860;
t822 = t756 * t671;
t864 = t717 * t718;
t823 = pkin(3) * t864;
t825 = 0.4e1 * t860;
t827 = 0.8e1 * t861;
t867 = t714 * t681 * t762;
t829 = -0.8e1 * t867;
t835 = 0.2e1 * t876;
t836 = pkin(7) * t672;
t837 = 4 * pkin(7);
t838 = t757 + t763;
t839 = t757 - t754;
t849 = 0.4e1 / 0.7e1 * t758 - t746 / 0.7e1;
t850 = t709 + t758;
t851 = t756 / 0.3e1 + t758;
t858 = t758 * t756;
t862 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t868 = t681 ^ 2 * t763;
t869 = t674 * t758;
t870 = t669 * (-t751 + t804);
t873 = (-t715 * t759 + t822) * t684;
t875 = ((-0.24e2 * (0.4e1 / 0.3e1 * t860 + t661 + t667) * t868 * t876 - 0.12e2 * (-0.8e1 / 0.3e1 * t793 + ((t707 + t810) * t671 - (0.7e1 / 0.6e1 * t751 + t693 + t850) * t876) * t825 + (-t751 * t840 - 0.5e1 / 0.3e1 * t754 + t811 * t756 + t758 * (t695 + t666)) * t671 + (-t763 + (-t798 + t812) * t751 - (3 * t754) + t813 * t756 + t757) * t818 + (-t715 * t754 * t683 + ((t756 + t810) * t671 + (t730 - t840) * t818) * t672) * t837) * t866 + 0.24e2 * t666 * t793 + ((t758 + 0.5e1 / 0.2e1 * t751 + 0.3e1 / 0.2e1 * t756 + t696) * t671 + t675 * t876 / 0.2e1) * t797 - 0.6e1 * ((-0.3e1 * t763 + (-t798 + t813) * t751 + t812 * t756 + t839) * t671 - 0.2e1 * (-0.5e1 / 0.3e1 * t763 + (-t756 + t811) * t751 + t758 * (t705 + t809)) * t876) * t860 - 0.6e1 * t783 * t836 - (t748 + ((21 * t756) + t662) * t763 + (t651 + t738 + (35 * t754) + 0.2e1 * t869) * t751 + (t731 + (t728 + t739 - 0.5e1 * t751) * t756 + t758 * (-t746 + t842)) * t676) * t671 + (0.7e1 * t748 + (t733 + t662) * t729 + (t651 + (21 * t754) + (9 * t757) + 0.6e1 * t869) * t751 + t870) * t876) * t622 + (0.16e2 * (t827 + t797 + (-8 * t754 + 12 * t858) * t684 + (-12 * pkin(7) * t759 + t767 * t884) * t718 - (6 * t858) + t841) * t868 + 0.24e2 * (t852 * t827 + 0.14e2 * (-0.32e2 / 0.21e2 * (t758 + t751 / 0.4e1 + t756 / 0.4e1 - t746 / 0.8e1) * t794 + 0.5e1 / 0.42e2 * t763 + (0.16e2 / 0.21e2 * t756 + t849) * t751 + t754 / 0.7e1 + t849 * t756 + t757 - 0.3e1 / 0.7e1 * t859 + t745 / 0.42e2) * t860 + t667 * t778 - t840 * t763 + (t664 - 0.10e2 / 0.3e1 * t754 + (2 * t757) - t859) * t751 + t636 * t862 + ((-0.2e1 / 0.3e1 * t794 + t694 + t850) * t819 + (-0.8e1 / 0.3e1 * (t851 + t883) * t794 + 0.5e1 / 0.18e2 * t763 + (0.4e1 / 0.3e1 * t758 + t707 + t695) * t751 + t757 + 0.2e1 / 0.3e1 * t858 - 0.2e1 / 0.3e1 * t859 - t754 / 0.3e1 + t745 / 0.18e2) * t832) * pkin(7)) * t866 + 0.16e2 * (-0.6e1 * t758 * t751 + t838) * t861 + 0.32e2 * (t652 * t863 + t655 * t675) * t830 + 0.24e2 * (t666 * t778 - t748 + (-t668 + t843) * t763 + (t664 + t763 / 0.6e1 - t745 / 0.6e1 + t839) * t751 + t636 * t758) * t860 + 0.8e1 * t777 * t836 - 0.8e1 * ((t733 + t853) * t763 + (t731 + (t728 + (10 * t758)) * t756 + t787) * t751 + t856) * t794 + t763 ^ 2 + (t727 + t740 + (28 * t756)) * t748 + (t673 * t756 + (70 * t754) + t786) * t763 + (t673 * t754 + t786 * t734 + t845 * t742 - 0.6e1 * t757 * t746 + (28 * t759 ^ 2) + (4 * t767 ^ 2)) * t751 + t656 * t870) * t647 + (((0.4e1 * t873 + (t671 + t835) * t661 + t678 * t671 + (t696 + t806) * t835) * t829 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t751 + t709 + t693) * t671 + pkin(1) * t814) * t826 + (-0.8e1 * t792 + ((t695 + t702 + t847) * t671 - (0.8e1 / 0.3e1 * t751 + t809) * t876) * t833) * pkin(7) + t783) * t670) * t622 + (0.32e2 * (t799 + (-0.4e1 * t715 * t821 + t879 + (t880 + t727 + t739) * t756) * t684 + (-t756 + t776 + t883) * t817 + t652 * t862 + t678 * t655) * t867 + 0.8e1 * (t650 + (t655 * t863 + t640) * t820 + (t778 + (t734 + t846) * t751 + t774) * t816 + t777) * t670) * t647) * t659) / ((-0.4e1 * (-t840 * t671 + 0.2e1 * t873 + (t823 * t882 + t715 * (t751 + t878)) * pkin(1)) * t866 + 0.8e1 * pkin(7) * t792 + ((pkin(3) * t879 + 0.8e1 * t756 * t762) * t717 + 0.4e1 * t759 * t814) * t684 - 0.4e1 * t782 * t836 - (t848 * t751 + t732 + t838 + (6 * t858)) * t671 + (t729 + (t722 + 6 * t758) * t751 + t669) * t876) * t622 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t794 + 0.4e1 / 0.9e1 * t751 - t746 / 0.9e1 + t851) * t860 + t667 * t652 + t653 * t862 + (t857 + (t693 + t702 + t776) * t672) * t837) * t866 + t650 + (t653 * t863 + t640) * t820 + t789 * t816 + ((t668 + t808) * t751 + t855) * t791 + t748 + (t721 + t805) * t763 + (t805 * t734 + t844 * t742 + t720 + t738) * t751 + t669 * t656) * t647 + ((t829 * t876 + (t822 * t881 + (t671 - t876) * t661 + t782) * t831) * t622 + (0.8e1 * (t661 + t825 + t678) * t867 + 0.6e1 * (t842 * t825 + (t652 + t653) * t817 + t789) * t670) * t647) * t659);
t631 = 0.1e1 / (t658 * t834 + t652 + t660 + t803);
t752 = 0.1e1 / pkin(3);
t874 = t631 * t752;
t784 = -pkin(1) + t824;
t785 = t730 + t735 - t843;
t620 = (t658 * t715 + t823) * t622 - (t660 + t785 + t790) * t872 + t784 * t660 + t785 * t824 + (t639 * t881 - t804 + t828 - t880) * pkin(1);
t644 = t730 + t788;
t621 = (-t784 + t872) * t622 + (t639 * t834 + t644 * t658) * t715 + (t644 * t718 + (0.4e1 * t684 - 0.2e1) * t658 * pkin(1)) * t671;
t801 = t874 / 0.2e1;
t618 = qJ(1) + atan2(t621 * t801, t620 * t801);
t747 = 0.1e1 / pkin(4);
t815 = t631 * t747 / pkin(3) ^ 2;
t800 = t747 * t752 / 0.2e1;
t602 = atan2(t622 * t800, t800 * t875) + t618;
t637 = t711 * t718 - t712 * t715;
t614 = (-t620 * t637 / 0.2e1 + t621 * t877) * t874;
t615 = (t620 * t877 + t621 * t637 / 0.2e1) * t874;
t716 = sin(pkin(9));
t719 = cos(pkin(9));
t611 = -t614 * t719 + t615 * t716;
t612 = t614 * t716 + t615 * t719;
t607 = atan2(t611, t612) + t618;
t802 = t875 / 0.4e1;
t617 = cos(t618);
t796 = -pkin(2) * t617 - t672;
t795 = -pkin(3) * t617 - t672;
t616 = sin(t618);
t781 = -pkin(2) * t616 - t876;
t780 = -pkin(3) * t616 - t876;
t779 = g(1) * t718 + g(2) * t715;
t645 = t779 * pkin(1);
t642 = -t714 * t718 + t865;
t641 = -t714 * t715 - t864;
t630 = qJ(1) + atan2(t637, t638);
t629 = pkin(8) + atan2(t637, -t638);
t628 = cos(t630);
t627 = sin(t630);
t626 = cos(t629);
t625 = sin(t629);
t606 = cos(t607);
t605 = sin(t607);
t604 = (t620 * t622 / 0.4e1 + t621 * t802) * t815;
t603 = (t620 * t802 - t621 * t622 / 0.4e1) * t815;
t601 = cos(t602);
t600 = sin(t602);
t599 = atan2(t611, -t612) + t607;
t598 = cos(t599);
t597 = sin(t599);
t596 = g(1) * t605 - g(2) * t606;
t595 = -g(1) * t606 - g(2) * t605;
t594 = atan2(t603 * t641 + t604 * t642, -t603 * t642 + t604 * t641) + t602;
t593 = cos(t594);
t592 = sin(t594);
t1 = [0, 0, 0, 0, 0, 0, t779, -g(1) * t715 + g(2) * t718, -g(3), 0, 0, 0, 0, 0, 0, 0, g(1) * t617 + g(2) * t616, -g(1) * t616 + g(2) * t617, -g(3), t645, 0, 0, 0, 0, 0, 0, t595, t596, -g(3), -g(1) * t796 - g(2) * t781, 0, 0, 0, 0, 0, 0, g(1) * t601 + g(2) * t600, -g(1) * t600 + g(2) * t601, -g(3), -g(1) * t795 - g(2) * t780, 0, 0, 0, 0, 0, 0, -g(1) * t626 - g(2) * t625, g(1) * t625 - g(2) * t626, -g(3), (-g(1) * t712 - g(2) * t711) * pkin(5), 0, 0, 0, 0, 0, 0, t595, t596, -g(3), t645, 0, 0, 0, 0, 0, 0, -g(1) * t714 + g(2) * t717, -g(1) * t717 - g(2) * t714, -g(3), -g(1) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t628 - g(2) * t627, g(1) * t627 - g(2) * t628, -g(3), t645, 0, 0, 0, 0, 0, 0, g(1) * t598 + g(2) * t597, -g(1) * t597 + g(2) * t598, -g(3), -g(2) * (pkin(6) * t605 + t781) - g(1) * (pkin(6) * t606 + t796), 0, 0, 0, 0, 0, 0, -g(1) * t593 - g(2) * t592, g(1) * t592 - g(2) * t593, -g(3), -g(2) * (-pkin(4) * t600 + t780) - g(1) * (-pkin(4) * t601 + t795);];
U_reg = t1;
