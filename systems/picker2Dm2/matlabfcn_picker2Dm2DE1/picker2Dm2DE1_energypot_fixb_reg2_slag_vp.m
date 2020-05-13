% Calculate inertial parameters regressor of potential energy for
% picker2Dm2DE1
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
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = picker2Dm2DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_energypot_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 16:03:04
% EndTime: 2020-05-09 16:03:13
% DurationCPUTime: 8.54s
% Computational Cost: add. (79863->432), mult. (235129->565), div. (1676->6), fcn. (42008->31), ass. (0->248)
t724 = sin(qJ(1));
t726 = cos(qJ(2));
t886 = t724 * t726;
t842 = pkin(3) * t886;
t814 = pkin(1) * t842;
t662 = -0.2e1 * t814;
t902 = 0.1e1 / pkin(3);
t895 = t902 / 0.2e1;
t752 = 0.1e1 / t895;
t763 = pkin(1) ^ 2;
t761 = t763 ^ 2;
t765 = pkin(7) ^ 2;
t727 = cos(qJ(1));
t693 = t727 ^ 2;
t723 = sin(qJ(2));
t669 = pkin(7) * t723 * t752;
t690 = t723 ^ 2;
t759 = pkin(3) ^ 2;
t887 = t690 * t759;
t846 = 0.2e1 * t887;
t860 = -t759 + t765;
t794 = t669 + t846 + t860;
t792 = t794 * t693;
t685 = t763 + t765;
t754 = pkin(4) ^ 2;
t823 = -t754 + t685;
t807 = t669 + t823;
t679 = pkin(3) * t723;
t897 = t679 + pkin(7);
t820 = t897 * t727;
t861 = t754 - t765;
t822 = -0.2e1 * t759 + t861;
t849 = -0.4e1 * t679;
t858 = t763 - t765;
t898 = 0.2e1 * t763;
t908 = 0.4e1 * pkin(1);
t631 = sqrt(-0.4e1 * t763 * t792 + 0.4e1 * t858 * t887 + pkin(7) * t823 * t849 - t761 + t822 * t898 - (t765 - (t752 + pkin(4)) * pkin(4)) * (t765 + (t752 - pkin(4)) * pkin(4)) + (-(t662 + t807) * t820 + t807 * t842) * t908);
t821 = t759 + t685;
t790 = 0.1e1 / (0.2e1 * pkin(1) * t820 + t662 + t669 + t821);
t803 = -pkin(1) + t842;
t743 = 0.3e1 * t763;
t804 = t743 - t822;
t809 = -0.4e1 * t814;
t885 = t726 * t727;
t841 = pkin(3) * t885;
t900 = 0.4e1 * t759;
t787 = t790 * ((t724 * t897 + t841) * t631 - (t669 + t804 + t809) * t820 + t803 * t669 + t804 * t842 + (-0.2e1 * t792 + t846 - t900 - t823) * pkin(1));
t785 = t787 / 0.4e1;
t680 = pkin(3) * t726;
t901 = 0.2e1 * t759;
t797 = t901 + t807;
t681 = pkin(1) * t727;
t852 = 0.2e1 * t681;
t789 = t790 * ((t820 - t803) * t631 + (t794 * t852 + t897 * t797) * t724 + (t797 * t727 + (0.4e1 * t693 - 0.2e1) * pkin(1) * t897) * t680);
t755 = 0.1e1 / pkin(4);
t882 = t755 / pkin(3) ^ 2;
t770 = t759 ^ 2;
t753 = t754 ^ 2;
t764 = t765 ^ 2;
t859 = t761 + t764;
t750 = 0.2e1 * t765;
t864 = t750 - t754;
t879 = t765 * t754;
t793 = t864 * t763 + t753 / 0.6e1 + t859 - t879;
t648 = -t770 / 0.6e1 + t793;
t714 = -t759 / 0.3e1;
t675 = t714 + t765;
t651 = t675 * t662;
t668 = t681 + pkin(7);
t657 = t679 + t668;
t684 = -0.3e1 * t759 + t765;
t692 = t727 * t693;
t766 = pkin(1) * t763;
t877 = t766 * t692;
t848 = pkin(7) * t877;
t819 = 0.8e1 * t848;
t660 = t684 * t819;
t683 = -t754 - t759;
t749 = 0.3e1 * t765;
t671 = t749 + t683;
t892 = t671 * t763;
t661 = 0.10e2 * t892;
t710 = 0.4e1 / 0.3e1 * t759;
t704 = -t754 / 0.3e1;
t828 = t704 + t685;
t663 = t710 + t828;
t705 = -t754 / 0.2e1;
t665 = t705 + t821;
t666 = -t754 + t821;
t670 = pkin(7) * t852;
t748 = 0.4e1 * t765;
t673 = (t748 + t754) * t763;
t676 = -t763 / 0.3e1 + t765;
t677 = 0.10e2 / 0.3e1 * t763;
t678 = t685 ^ 2;
t682 = -0.30e2 * t754 + 0.60e2 * t765;
t687 = -0.3e1 * t763 + t765;
t702 = -t754 / 0.6e1;
t703 = -t754 / 0.4e1;
t711 = 0.2e1 / 0.3e1 * t759;
t716 = 0.4e1 / 0.3e1 * t763;
t718 = t763 / 0.2e1;
t729 = 0.15e2 * t761;
t730 = 0.15e2 * t763;
t731 = 0.10e2 * t763;
t736 = -0.2e1 * t754;
t737 = -0.5e1 * t754;
t738 = 0.5e1 * t770;
t739 = 0.7e1 * t761;
t740 = 0.5e1 * t761;
t741 = 0.7e1 * t763;
t742 = 0.6e1 * t763;
t746 = 0.3e1 * t764;
t747 = 0.8e1 * t765;
t769 = pkin(3) * t759;
t756 = t769 ^ 2;
t774 = pkin(7) * t765;
t791 = 0.5e1 / 0.6e1 * t770 + t793;
t796 = t765 - t814;
t872 = t753 / 0.2e1 - t770 / 0.2e1;
t806 = -0.3e1 * t879 + t746 + t872;
t810 = -0.6e1 * t814;
t707 = -0.3e1 / 0.2e1 * t754;
t871 = t707 + t749;
t874 = t685 * ((t707 + t750) * t763 - 0.3e1 / 0.2e1 * t879 + t859 + t872) + t756;
t798 = ((t677 + t864) * t759 + t791) * t810 + (t729 + (-0.9e1 * t754 + 0.18e2 * t765) * t763 + t806) * t759 + (t730 + t871) * t770 + t874;
t799 = t665 * t809;
t865 = t743 + t765;
t825 = t759 + t865;
t896 = pkin(1) * t724;
t801 = -(0.3e1 * t759 + t685) * t896 + t825 * t680;
t706 = -0.2e1 / 0.3e1 * t754;
t715 = -0.2e1 / 0.3e1 * t759;
t826 = t706 + t685;
t866 = t731 + t750;
t870 = t715 + t765;
t802 = -(t738 + (0.5e1 * t763 + t671) * t901 + (t715 + t826) * t685) * t896 + (t770 + (t706 + t715 + t866) * t759 + t740 + 0.2e1 * t892 + t765 * (t706 + t870)) * t680;
t863 = t753 - t770;
t805 = -0.6e1 * t879 + 0.6e1 * t764 + t863;
t827 = t706 + t711 + t750;
t873 = (t711 + t826) * t685 + t770;
t808 = t663 * t809 + (t742 + t827) * t759 + t873;
t811 = t877 * t680;
t881 = t761 * t693 ^ 2;
t812 = t881 * t680;
t839 = 0.16e2 * t877;
t817 = pkin(7) * t839;
t818 = 0.20e2 / 0.3e1 * t763;
t862 = -t754 + t759;
t824 = t749 + t862;
t712 = t759 / 0.3e1;
t829 = t702 + t712 + t765;
t830 = t754 / 0.3e1 + t712 + t750;
t831 = 0.2e1 / 0.3e1 * t754 + t711 + t748;
t832 = 0.4e1 / 0.3e1 * t754 + t710 - 0.2e1 * t765;
t884 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t833 = t724 * t884;
t850 = 0.6e1 * t681;
t834 = pkin(7) * t850;
t851 = 0.4e1 * t681;
t835 = pkin(7) * t851;
t837 = -t896 / 0.2e1;
t880 = t763 * t693;
t840 = 0.12e2 * t880;
t843 = t763 * t680;
t844 = 0.4e1 * t880;
t845 = 0.8e1 * t881;
t888 = t723 * t690 * t769;
t847 = -0.8e1 * t888;
t853 = 0.2e1 * t896;
t854 = pkin(7) * t681;
t855 = 0.4e1 * pkin(7);
t856 = t764 + t770;
t857 = t764 - t761;
t867 = 0.4e1 / 0.7e1 * t765 - t754 / 0.7e1;
t868 = t718 + t765;
t869 = t763 / 0.3e1 + t765;
t878 = t765 * t763;
t883 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t889 = t690 ^ 2 * t770;
t890 = t683 * t765;
t891 = t678 * (-t759 + t823);
t893 = (-t724 * t766 + t843) * t693;
t899 = 0.4e1 * t761;
t903 = t703 + t759 / 0.2e1;
t894 = ((-0.24e2 * (0.4e1 / 0.3e1 * t880 + t670 + t676) * t889 * t896 - 0.12e2 * (-0.8e1 / 0.3e1 * t812 + ((t716 + t829) * t680 - (0.7e1 / 0.6e1 * t759 + t702 + t868) * t896) * t844 + (-t759 * t858 - 0.5e1 / 0.3e1 * t761 + t830 * t763 + t765 * (t704 + t675)) * t680 + (-t770 + (-t818 + t831) * t759 - 0.3e1 * t761 + t832 * t763 + t764) * t837 + (-t724 * t761 * t692 + ((t763 + t829) * t680 + (t901 - t858) * t837) * t681) * t855) * t887 + 0.24e2 * t675 * t812 + ((t765 + 0.5e1 / 0.2e1 * t759 + 0.3e1 / 0.2e1 * t763 + t705) * t680 + t684 * t896 / 0.2e1) * t817 - 0.6e1 * ((-0.3e1 * t770 + (-t818 + t832) * t759 + t831 * t763 + t857) * t680 - 0.2e1 * (-0.5e1 / 0.3e1 * t770 + (-t763 + t830) * t759 + t765 * (t714 + t828)) * t896) * t880 - 0.6e1 * t802 * t854 - (t756 + (0.21e2 * t763 + t671) * t770 + (t661 + t746 + 0.35e2 * t761 + 0.2e1 * t890) * t759 + (t739 + (t737 + t747 - 0.5e1 * t759) * t763 + t765 * (-t754 + t860)) * t685) * t680 + (0.7e1 * t756 + (t741 + t671) * t738 + (t661 + 0.21e2 * t761 + 0.9e1 * t764 + 0.6e1 * t890) * t759 + t891) * t896) * t631 + (0.16e2 * (t845 + t817 + (-0.8e1 * t761 + 0.12e2 * t878) * t693 + (-0.12e2 * pkin(7) * t766 + t774 * t908) * t727 - 0.6e1 * t878 + t859) * t889 + 0.24e2 * (t870 * t845 + 0.14e2 * (-0.32e2 / 0.21e2 * (t765 + t759 / 0.4e1 + t763 / 0.4e1 - t754 / 0.8e1) * t814 + 0.5e1 / 0.42e2 * t770 + (0.16e2 / 0.21e2 * t763 + t867) * t759 + t761 / 0.7e1 + t867 * t763 + t764 - 0.3e1 / 0.7e1 * t879 + t753 / 0.42e2) * t880 + t676 * t799 - t858 * t770 + (t673 - 0.10e2 / 0.3e1 * t761 + 0.2e1 * t764 - t879) * t759 + t648 * t883 + ((-0.2e1 / 0.3e1 * t814 + t703 + t868) * t839 + (-0.8e1 / 0.3e1 * (t869 + t903) * t814 + 0.5e1 / 0.18e2 * t770 + (0.4e1 / 0.3e1 * t765 + t716 + t704) * t759 + t764 + 0.2e1 / 0.3e1 * t878 - 0.2e1 / 0.3e1 * t879 - t761 / 0.3e1 + t753 / 0.18e2) * t850) * pkin(7)) * t887 + 0.16e2 * (-0.6e1 * t759 * t765 + t856) * t881 + 0.32e2 * (t662 * t884 + t665 * t684) * t848 + 0.24e2 * (t675 * t799 - t756 + (-t677 + t861) * t770 + (t673 + t770 / 0.6e1 - t753 / 0.6e1 + t857) * t759 + t648 * t765) * t880 + 0.8e1 * t798 * t854 - 0.8e1 * ((t741 + t871) * t770 + (t739 + (t737 + 0.10e2 * t765) * t763 + t806) * t759 + t874) * t814 + t770 ^ 2 + (t736 + t748 + 0.28e2 * t763) * t756 + (t682 * t763 + 0.70e2 * t761 + t805) * t770 + (t682 * t761 + t805 * t742 + t863 * t750 - 0.6e1 * t764 * t754 + 0.28e2 * t766 ^ 2 + 0.4e1 * t774 ^ 2) * t759 + t666 * t891) * t657 + (((0.4e1 * t893 + (t680 + t853) * t670 + t687 * t680 + (t705 + t825) * t853) * t847 - 0.6e1 * (-0.4e1 * ((0.5e1 / 0.6e1 * t759 + t718 + t702) * t726 * t752 + pkin(1) * t833) * t880 + (-0.8e1 * t811 + ((t704 + t711 + t865) * t680 - (0.8e1 / 0.3e1 * t759 + t828) * t896) * t851) * pkin(7) + t802) * t679) * t631 + (0.32e2 * (t819 + (-0.4e1 * t766 * t842 + t899 + (t900 + t736 + t747) * t763) * t693 + (-t763 + t796 + t903) * t835 + t662 * t883 + t687 * t665) * t888 + 0.8e1 * (t660 + (t665 * t884 + t651) * t840 + (t799 + (t742 + t864) * t759 + t791) * t834 + t798) * t679) * t657) * t668) / ((-0.4e1 * (-t858 * t680 + 0.2e1 * t893 + (0.2e1 * pkin(7) * t841 + t724 * (t759 + t898)) * pkin(1)) * t887 + 0.8e1 * pkin(7) * t811 + ((pkin(3) * t899 + 0.8e1 * t763 * t769) * t726 + 0.4e1 * t766 * t833) * t693 - 0.4e1 * t801 * t854 - (t866 * t759 + t740 + t856 + 0.6e1 * t878) * t680 + (t738 + (t731 + 0.6e1 * t765) * t759 + t678) * t896) * t631 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t814 + 0.4e1 / 0.9e1 * t759 - t754 / 0.9e1 + t869) * t880 + t676 * t662 + t663 * t883 + (t877 + (t702 + t711 + t796) * t681) * t855) * t887 + t660 + (t663 * t884 + t651) * t840 + t808 * t834 + ((t677 + t827) * t759 + t873) * t810 + t756 + (t730 + t824) * t770 + (t824 * t742 + t862 * t750 + t729 + t746) * t759 + t678 * t666) * t657 + ((t847 * t896 + (-0.2e1 * t693 * t843 + (t680 - t896) * t670 + t801) * t849) * t631 + (0.8e1 * (t670 + t844 + t687) * t888 + 0.6e1 * (t860 * t844 + (t662 + t663) * t835 + t808) * t679) * t657) * t668);
t606 = (t785 * t894 - t631 * t789 / 0.4e1) * t882;
t607 = (t631 * t785 + t789 * t894 / 0.4e1) * t882;
t652 = -t723 * t724 - t885;
t653 = -t723 * t727 + t886;
t596 = atan2(t606 * t652 + t607 * t653, -t606 * t653 + t607 * t652);
t594 = sin(t596);
t595 = cos(t596);
t836 = t755 * t895;
t610 = atan2(t631 * t836, t836 * t894);
t608 = sin(t610);
t609 = cos(t610);
t786 = t902 * t787;
t784 = t786 / 0.2e1;
t788 = t789 * t895;
t783 = atan2(t788, t784);
t781 = sin(t783);
t782 = cos(t783);
t624 = -t724 * t782 - t727 * t781;
t625 = t724 * t781 - t727 * t782;
t598 = t608 * t625 + t609 * t624;
t813 = -t608 * t624 + t625 * t609;
t912 = t594 * t598 - t595 * t813;
t911 = t594 * t813 + t595 * t598;
t720 = sin(pkin(8));
t721 = cos(pkin(8));
t649 = t720 * t727 - t721 * t724;
t650 = t720 * t724 + t721 * t727;
t627 = -t649 * t786 / 0.2e1 + t650 * t788;
t628 = t649 * t788 + t650 * t784;
t725 = sin(pkin(9));
t728 = cos(pkin(9));
t618 = -t627 * t728 + t628 * t725;
t619 = t627 * t725 + t628 * t728;
t615 = atan2(t618, t619);
t611 = sin(t615);
t613 = cos(t615);
t603 = t611 * t625 + t613 * t624;
t616 = atan2(t618, -t619);
t612 = sin(t616);
t614 = cos(t616);
t795 = t611 * t624 - t613 * t625;
t910 = t603 * t612 + t614 * t795;
t909 = t603 * t614 - t612 * t795;
t876 = t625 * pkin(3) - t681;
t875 = t625 * pkin(2) - t681;
t816 = t624 * pkin(3) - t896;
t815 = t624 * pkin(2) - t896;
t800 = g(1) * t727 + g(2) * t724;
t655 = t800 * pkin(1);
t643 = atan2(t649, -t650);
t642 = atan2(t649, t650);
t641 = cos(t643);
t640 = cos(t642);
t639 = sin(t643);
t638 = sin(t642);
t635 = -t638 * t724 + t640 * t727;
t634 = t638 * t727 + t640 * t724;
t633 = -t639 * t720 + t641 * t721;
t632 = t639 * t721 + t641 * t720;
t593 = -g(1) * t795 + g(2) * t603;
t592 = -g(1) * t603 - g(2) * t795;
t1 = [0, 0, 0, 0, 0, 0, t800, -g(1) * t724 + g(2) * t727, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(1) * t625 - g(2) * t624, g(1) * t624 - g(2) * t625, -g(3), t655, 0, 0, 0, 0, 0, 0, t593, t592, -g(3), -g(1) * t875 - g(2) * t815, 0, 0, 0, 0, 0, 0, -g(1) * t813 - g(2) * t598, g(1) * t598 - g(2) * t813, -g(3), -g(1) * t876 - g(2) * t816, 0, 0, 0, 0, 0, 0, -g(1) * t633 - g(2) * t632, g(1) * t632 - g(2) * t633, -g(3), (-g(1) * t721 - g(2) * t720) * pkin(5), 0, 0, 0, 0, 0, 0, t593, t592, -g(3), t655, 0, 0, 0, 0, 0, 0, -g(1) * t723 + g(2) * t726, -g(1) * t726 - g(2) * t723, -g(3), -g(1) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t635 - g(2) * t634, g(1) * t634 - g(2) * t635, -g(3), t655, 0, 0, 0, 0, 0, 0, g(1) * t910 - g(2) * t909, g(1) * t909 + g(2) * t910, -g(3), -g(2) * (-pkin(6) * t603 + t815) - g(1) * (pkin(6) * t795 + t875), 0, 0, 0, 0, 0, 0, -g(1) * t912 + g(2) * t911, -g(1) * t911 - g(2) * t912, -g(3), -g(2) * (pkin(4) * t598 + t816) - g(1) * (pkin(4) * t813 + t876);];
U_reg = t1;
