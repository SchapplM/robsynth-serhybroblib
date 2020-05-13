% Calculate kinetic energy for
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:13:29
% EndTime: 2020-04-27 09:13:54
% DurationCPUTime: 25.20s
% Computational Cost: add. (362011->739), mult. (1084238->1119), div. (2756->13), fcn. (134114->6), ass. (0->476)
t844 = (pkin(5) ^ 2);
t843 = t844 ^ 2;
t847 = pkin(4) ^ 2;
t846 = t847 ^ 2;
t1010 = t843 - t846;
t857 = (pkin(1) ^ 2);
t842 = -2 * t857;
t814 = cos(qJ(1));
t1062 = pkin(2) * t814;
t742 = pkin(1) - t1062;
t813 = cos(qJ(2));
t1049 = t742 * t813;
t855 = pkin(2) ^ 2;
t1007 = -t855 + t857;
t991 = pkin(1) * t1062;
t745 = -0.2e1 * t991;
t779 = t814 ^ 2;
t1039 = t779 * t855;
t985 = 0.2e1 * t1039;
t681 = t745 + t985 + t1007;
t776 = t813 ^ 2;
t1050 = t681 * t776;
t1072 = -pkin(4) + pkin(5);
t1073 = -pkin(4) - pkin(5);
t851 = pkin(3) ^ 2;
t1003 = t857 - t844;
t946 = t855 + t1003;
t915 = t851 + t946;
t731 = -t847 + t915;
t678 = t745 + t731;
t811 = sin(qJ(2));
t812 = sin(qJ(1));
t1038 = t811 * t812;
t974 = pkin(3) * t1038;
t932 = pkin(2) * t974;
t733 = -0.2e1 * t932;
t657 = t733 + t678;
t1029 = t855 * t857;
t1098 = 0.4e1 * t851;
t780 = t855 * t1098;
t746 = t780 - 0.4e1 * t1029;
t826 = 0.2e1 * t847;
t849 = t851 ^ 2;
t978 = pkin(2) * t1038;
t993 = 0.4e1 * pkin(3);
t858 = sqrt(t746 * t779 + 0.4e1 * t731 * t991 - t849 - (t857 + (pkin(2) - t1073) * (pkin(2) + t1073)) * (t857 + (pkin(2) - t1072) * (pkin(2) + t1072)) + (t826 + t842 + (2 * t844) - 0.6e1 * t855 - 0.4e1 * t1050) * t851 + (-t1049 * t657 + t678 * t978) * t993);
t1065 = -t858 / 0.4e1;
t1030 = t851 * t776;
t775 = t813 * t776;
t866 = pkin(3) * t851;
t1043 = t775 * t866;
t769 = t851 + t857;
t763 = t769 ^ 2;
t947 = t851 + t1003;
t1046 = t763 * (-t847 + t947);
t766 = pkin(3) * t813;
t743 = pkin(1) + t766;
t1048 = t743 * t814;
t1028 = t857 * t844;
t1033 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t1016 = 0.4e1 / 0.7e1 * t857 - t844 / 0.7e1;
t726 = t857 + t855 / 0.4e1 + t851 / 0.4e1 - t844 / 0.8e1;
t856 = t857 ^ 2;
t863 = t855 ^ 2;
t611 = -0.32e2 / 0.21e2 * t726 * t932 + t863 / 0.7e1 + (0.16e2 / 0.21e2 * t851 + t1016) * t855 + t849 / 0.7e1 + t1016 * t851 + t856 - 0.3e1 / 0.7e1 * t1028 + t1010 / 0.42e2;
t1066 = 0.4e1 / 0.3e1 * t855;
t809 = t855 / 0.2e1;
t1017 = t809 + t857;
t791 = -t844 / 0.4e1;
t806 = t851 / 0.3e1;
t727 = t791 + t806 + t1017;
t792 = -t844 / 0.3e1;
t794 = -0.2e1 / 0.3e1 * t844;
t805 = 0.4e1 / 0.3e1 * t851;
t613 = -0.8e1 / 0.3e1 * t727 * t932 + t863 / 0.3e1 + (t805 + t792) * t855 + t856 - t849 / 0.3e1 + (t1066 + 0.2e1 / 0.3e1 * t851 + t794) * t857 + t1010 / 0.18e2;
t1008 = t849 + t856;
t841 = 2 * t857;
t1011 = t841 - t844;
t1090 = -t843 / 0.6e1 + t846 / 0.6e1;
t671 = t1011 * t851 + t1008 - t1028 - t1090;
t807 = t851 / 0.2e1;
t674 = -0.2e1 / 0.3e1 * t932 + t857 + t807 + t791;
t839 = 4 * t857;
t753 = (t839 + t844) * t851;
t758 = -t851 / 0.3e1 + t857;
t760 = t857 - 0.2e1 / 0.3e1 * t855;
t770 = -t851 + t857;
t793 = -t844 / 0.2e1;
t1006 = t855 + t857;
t948 = t851 + t1006;
t735 = t793 + t948;
t919 = -0.4e1 * t932;
t897 = t735 * t919;
t972 = 0.16e2 * t1043;
t774 = t776 ^ 2;
t1045 = t774 * t849;
t989 = 0.8e1 * t1045;
t594 = t760 * t989 + 0.14e2 * t611 * t1030 + t758 * t897 + t770 * t863 + (-t1028 + t753 - 0.10e2 / 0.3e1 * t849 + (2 * t856)) * t855 + t671 * t1033 + (0.6e1 * t613 * t766 + t674 * t972) * pkin(1);
t840 = 3 * t857;
t1015 = 0.15e2 * t851 + t840;
t1020 = t843 / 0.2e1 - t846 / 0.2e1;
t795 = -0.3e1 / 0.2e1 * t844;
t862 = pkin(2) * t855;
t852 = t862 ^ 2;
t1024 = t769 * ((t795 + t841) * t851 - 0.3e1 / 0.2e1 * t1028 + t1008 + t1020) + t852;
t782 = 0.10e2 / 0.3e1 * t851;
t880 = t671 + t863;
t643 = (t782 + t1011) * t855 + t880;
t819 = 0.15e2 * t849;
t822 = 18 * t857;
t836 = 3 * t856;
t918 = -(3 * t1028) + t836 + t1020;
t601 = -0.6e1 * t643 * t932 + t1024 + (t819 + (t822 - 9 * t844) * t851 + t918) * t855 + (t795 + t1015) * t863;
t831 = 0.6e1 * t851;
t617 = t897 + (t831 + t1011) * t855 + t880;
t1032 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t759 = t857 - t855 / 0.3e1;
t682 = t759 * t733;
t639 = t1032 * t735 + t682;
t772 = -0.3e1 * t855 + t857;
t988 = 0.8e1 * t1043;
t934 = pkin(1) * t988;
t718 = t772 * t934;
t990 = pkin(1) * t766;
t959 = 0.6e1 * t990;
t970 = 0.12e2 * t1030;
t595 = t617 * t959 + t639 * t970 + t601 + t718;
t1079 = 0.4e1 * t849;
t1083 = -0.4e1 * t811;
t823 = -2 * t844;
t837 = 8 * t857;
t767 = pkin(2) * t812;
t976 = t866 * t767;
t666 = t976 * t1083 + t780 + t1079 + (t823 + t837) * t851;
t884 = -t932 + t1017;
t672 = t791 - t851 + t884;
t771 = -0.3e1 * t851 + t857;
t992 = 0.4e1 * t766;
t604 = t733 * t1033 + t666 * t776 + t735 * t771 + (t672 * t992 + t988) * pkin(1);
t1005 = t856 - t849;
t605 = t759 * t897 - t852 + (-t782 - t1003) * t863 + (t753 + t1005 + t1090) * t855 + t671 * t857;
t824 = -5 * t844;
t829 = 0.7e1 * t849;
t610 = (t795 + t840 + 0.7e1 * t851) * t863 + (t829 + (t824 + 10 * t857) * t851 + t918) * t855 + t1024;
t1027 = t857 * t851;
t859 = pkin(1) * t857;
t741 = -0.12e2 * pkin(1) * t866 + t859 * t993;
t755 = -0.8e1 * t849 + 0.12e2 * t1027;
t931 = pkin(1) * t972;
t626 = t741 * t813 + t755 * t776 + t1008 - 0.6e1 * t1027 + t931 + t989;
t644 = t1032 * t733 + t735 * t772;
t1004 = t856 + t863;
t720 = 0.16e2 * (t1004 - 0.6e1 * t1029) * t849;
t739 = t847 + t947;
t764 = -30 * t844 + 60 * t857;
t777 = t779 ^ 2;
t825 = -6 * t844;
t887 = t1010 + (6 * t856) - (6 * t1028);
t778 = t814 * t779;
t1040 = t778 * t862;
t962 = t743 * t1040;
t922 = -0.32e2 * t962;
t575 = t604 * t922 + t720 * t774 + 0.24e2 * t605 * t1030 + (t823 + t839 + 0.28e2 * t851) * t852 + t739 * t1046 + (t1010 * t841 + 0.24e2 * t594 * t779 + t764 * t849 + (t825 * t856) + t887 * t831 + 0.4e1 * t859 ^ 2 + 0.28e2 * t866 ^ 2) * t855 + 0.8e1 * (-t1048 * t595 - t610 * t974) * pkin(2) + (0.32e2 * t1043 * t644 + 0.8e1 * t601 * t766) * pkin(1) + (0.16e2 * t626 * t777 + t764 * t851 + 0.70e2 * t849 + t863 + t887) * t863;
t1063 = pkin(1) * t813;
t1077 = 0.24e2 * t759;
t1044 = t775 * t849;
t1059 = t811 * pkin(3);
t790 = -t844 / 0.6e1;
t1019 = t790 - t847 / 0.6e1;
t954 = t857 + t1019;
t691 = t1066 + t807 + t954;
t917 = t809 + t954;
t692 = t805 + t917;
t629 = -t1059 * t691 + t692 * t767;
t695 = t851 + t917;
t751 = 0.2e1 * t855 + t770;
t965 = -t1059 / 0.2e1;
t642 = t695 * t767 + t751 * t965;
t1094 = t792 - t847 / 0.3e1;
t953 = t857 + t1094;
t800 = t847 / 0.3e1;
t955 = t844 / 0.3e1 + t800 + t841;
t645 = t855 * t770 - 0.5e1 / 0.3e1 * t849 + t955 * t851 + t857 * t953;
t783 = -0.20e2 / 0.3e1 * t851;
t799 = 0.2e1 / 0.3e1 * t847;
t956 = 0.2e1 / 0.3e1 * t844 + t799 + t839;
t957 = 0.4e1 / 0.3e1 * t844 + 0.4e1 / 0.3e1 * t847 + t842;
t646 = -t863 + (t783 + t956) * t855 - 0.3e1 * t849 + t957 * t851 + t856;
t983 = 0.4e1 * t1030;
t995 = 0.4e1 * pkin(1);
t596 = t629 * t983 + t646 * t965 + (-0.8e1 / 0.3e1 * t1045 + t645) * t767 + (-t1044 * t811 + t642 * t766) * t995;
t821 = 0.10e2 * t851;
t1014 = t821 + t841;
t804 = -0.2e1 / 0.3e1 * t847;
t1018 = t794 + t804;
t830 = 0.5e1 * t849;
t1013 = t823 - 0.2e1 * t847;
t838 = 6 * t857;
t949 = t838 + t1013;
t632 = t863 + (t1014 + t1018) * t855 + t830 + t949 * t851 + t857 * (t857 + t1018);
t834 = 0.5e1 * t863;
t950 = t794 + t769;
t650 = t834 + (t821 + t949) * t855 + (t804 + t950) * t769;
t606 = -t1059 * t650 + t632 * t767;
t835 = 0.3e1 * t855;
t916 = t851 + t953;
t697 = t835 + t916;
t832 = 0.3e1 * t851;
t752 = t832 + t1006;
t698 = t752 + t1094;
t633 = -t1059 * t697 + t698 * t767;
t1099 = 0.2e1 * t767;
t725 = t807 + t855 + t1019;
t651 = t1032 * t1059 + t725 * t1099;
t597 = -0.4e1 * t651 * t1030 + (t633 * t992 - 0.8e1 * t775 * t976) * pkin(1) + t606;
t1101 = 0.10e2 / 0.3e1 * t863 - 0.2e1 * (-t851 + t955) * t855 + t916 * t842;
t647 = -0.3e1 * t863 + (t783 + t957) * t855 + t956 * t851 + t1005;
t607 = t1059 * t1101 + t647 * t767;
t952 = t793 - t847 / 0.2e1 + t857;
t694 = 0.3e1 / 0.2e1 * t855 + t832 + t952;
t1061 = pkin(2) * t851;
t977 = t812 * t1061;
t885 = -t811 * t866 + t977;
t701 = 0.4e1 * t885;
t729 = t767 + 0.2e1 * t1059;
t994 = 0.2e1 * pkin(3);
t608 = t771 * t767 + t701 * t776 + (t1063 * t729 + t694 * t811) * t994;
t1012 = t824 - 0.5e1 * t847;
t1092 = t840 - t844 - t847;
t736 = t1092 * t821;
t609 = t852 + (0.21e2 * t851 + t1092) * t863 + (t1013 * t857 + t736 + t836 + 0.35e2 * t849) * t855 + (t829 + (t837 + t1012) * t851 + t857 * (-t847 + t1003)) * t769;
t616 = 0.7e1 * t852 + (0.35e2 * t851 + (15 * t857) + t1012) * t863 + (0.21e2 * t849 + t736 + (9 * t856) + (t825 - 0.6e1 * t847) * t857) * t855 + t1046;
t696 = t835 + 0.3e1 / 0.2e1 * t851 + t952;
t966 = pkin(3) * t772 / 0.2e1;
t649 = t696 * t767 + t811 * t966;
t1041 = t777 * t863;
t744 = 0.2e1 * t990;
t675 = 0.4e1 / 0.3e1 * t1030 + t744 + t758;
t923 = -0.24e2 * t675 * t1041;
t927 = 0.8e1 * t962;
t971 = -0.12e2 * t1039;
t984 = -0.6e1 * t1030;
t581 = t649 * t931 + t608 * t927 + t596 * t971 + t607 * t984 + (-0.6e1 * t606 * t1063 + (t616 + t923) * t811) * pkin(3) + (0.6e1 * t597 * t1048 + (t1045 * t1077 - t609) * t812) * pkin(2);
t700 = t742 + t766;
t566 = t575 * t700 + t581 * t858;
t677 = t745 + t847 + t915;
t619 = t677 * t742 + 0.2e1 * t681 * t766;
t622 = t677 * t813 + (0.4e1 * t776 - 0.2e1) * t742 * pkin(3);
t732 = -pkin(3) + t978;
t668 = -t732 + t1049;
t592 = t619 * t811 + t622 * t767 + t668 * t858;
t660 = -0.4e1 / 0.9e1 * t932 + t857 + t855 / 0.3e1 + t806 + t847 / 0.9e1 - t844 / 0.9e1;
t669 = t847 / 0.6e1 + t790 + t884;
t699 = t792 + t800 + t948;
t602 = t758 * t733 + 0.6e1 * t660 * t1030 + t699 * t1033 + (t669 * t766 + t1043) * t995;
t1022 = (t799 + t950) * t769 + t863;
t951 = t794 + t799 + t841;
t614 = t699 * t919 + t1022 + (t831 + t951) * t855;
t667 = t733 + t699;
t747 = t1007 * t1098;
t603 = 0.4e1 * t667 * t990 + t747 * t776 + t614;
t636 = t1032 * t699 + t682;
t641 = (t782 + t951) * t855 + t1022;
t1095 = t771 + t983;
t680 = t744 + t1095;
t928 = -0.8e1 * t962;
t586 = t680 * t928 + t718 + t636 * t970 + t614 * t959 + t852 + (-t844 + t847 + t1015) * t863 + t763 * t739 + (0.12e2 * t602 * t779 + t819 + (t822 + t825 + 0.6e1 * t847) * t851 + t836 + (t823 + t826) * t857) * t855 + 0.6e1 * (-t1048 * t603 - t641 * t974) * pkin(2);
t750 = t835 + t769;
t1047 = t750 * t811;
t1084 = -0.2e1 * t776;
t719 = t752 * t767;
t728 = t767 - t1059;
t615 = t977 * t1084 + t719 + (0.2e1 * t1063 * t728 - t1047) * pkin(3);
t702 = 0.2e1 * t885;
t768 = 0.2e1 * t851 + t855;
t914 = t770 + t744;
t618 = t1059 * t768 + t702 * t776 + t767 * t914;
t659 = -pkin(3) * t1047 + t719;
t748 = pkin(2) * t1079 + 0.8e1 * t851 * t862;
t960 = t866 * t1032;
t663 = t748 * t812 + 0.4e1 * t811 * t960;
t679 = t1014 * t855 + t1004 + 0.6e1 * t1027 + t830;
t686 = t834 + (t821 + t838) * t855 + t763;
t986 = -0.4e1 * t1039;
t593 = t618 * t986 + t663 * t776 + (-0.4e1 * t659 * t1063 + (t686 + t927) * t811) * pkin(3) + (0.4e1 * t615 * t1048 + (-t679 + t934) * t812) * pkin(2);
t580 = t586 * t700 + t593 * t858;
t1071 = 0.1e1 / t580 / 0.4e1;
t1080 = -0.2e1 * t814;
t730 = t832 + t847 + t946;
t658 = t730 + t745 + t919;
t1035 = t813 * t812;
t670 = pkin(2) * t1035 + t742 * t811;
t590 = -t658 * t1049 + t670 * t858 + (pkin(1) * t1080 * t732 + t1038 * t730) * pkin(2) + (-t739 - t835 + t985 - 0.2e1 * t1050) * pkin(3);
t937 = t590 * t1071;
t882 = t1065 * t592 + t566 * t937;
t1031 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t975 = pkin(3) * t1049;
t637 = t733 + t745 + t948 + 0.2e1 * t975;
t634 = 0.1e1 / t637;
t963 = t634 * t1031;
t562 = t882 * t963;
t1064 = t858 / 0.4e1;
t936 = t592 * t1071;
t883 = t1064 * t590 + t566 * t936;
t563 = t883 * t963;
t1037 = t811 * t814;
t684 = -t1035 + t1037;
t685 = t813 * t814 + t1038;
t555 = t562 * t684 + t563 * t685;
t556 = t562 * t685 - t563 * t684;
t550 = t555 * t814 + t556 * t812;
t1107 = Icges(3,4) * t550;
t1001 = qJD(1) * t814;
t749 = pkin(2) * t1001;
t929 = pkin(3) * t749;
t904 = t811 * t929;
t1074 = pkin(2) * pkin(3);
t997 = qJD(2) * t813;
t941 = t812 * t997;
t905 = t941 * t1074;
t1106 = -0.6e1 * t905 - 0.6e1 * t904;
t1105 = t905 + t904;
t1002 = qJD(1) * t812;
t945 = t743 * t1002;
t999 = qJD(2) * t811;
t968 = pkin(3) * t999;
t1093 = pkin(2) * t945 + t968 * t1062;
t1104 = 0.2e1 * pkin(1);
t1100 = 0.2e1 * t749;
t926 = t555 * t812 - t814 * t556;
t1097 = t926 / 0.2e1;
t1060 = pkin(3) * t776;
t1096 = Icges(3,4) * t926;
t940 = t811 * t1001;
t881 = t940 + t941;
t1091 = qJD(1) - qJD(2);
t1089 = 0.1e1 / pkin(4);
t1087 = 0.8e1 * pkin(1);
t1086 = -0.4e1 * pkin(3);
t1085 = 0.1e1 / t590;
t1082 = -0.2e1 * t811;
t1081 = -0.2e1 * t813;
t1078 = -0.8e1 * t851;
t1076 = pkin(1) * pkin(2);
t1075 = pkin(1) * pkin(3);
t1069 = -t811 / 0.2e1;
t1058 = t1089 / 0.2e1;
t964 = t813 * t1058;
t584 = (t1069 * t1089 * t592 + t590 * t964) * t634;
t1070 = t584 / 0.2e1;
t1068 = t813 / 0.2e1;
t1067 = t814 / 0.2e1;
t1057 = Icges(2,4) * t812;
t1056 = Icges(4,4) * t811;
t583 = (t1058 * t590 * t811 + t592 * t964) * t634;
t1055 = Icges(5,4) * t583;
t1000 = qJD(1) * t855;
t969 = pkin(2) * t1002;
t738 = t969 * t1104;
t912 = t812 * t814 * t1000;
t903 = -0.4e1 * t912;
t673 = t738 + t903;
t1051 = t673 * t776;
t942 = t811 * t997;
t911 = t681 * t942;
t924 = 0.4e1 * t941;
t1021 = 0.2e1 * t1105;
t958 = t738 - t1021;
t1053 = ((0.8e1 * t911 - 0.4e1 * t1051) * t851 + (-0.4e1 * t1076 * t731 + t1080 * t746) * t1002 + (t811 * t812 ^ 2 * t1000 * t1087 + 0.4e1 * (t657 * t999 - t813 * t958) * t742 + (t678 * t924 + 0.4e1 * (-t1035 * t657 + t1037 * t678) * qJD(1)) * pkin(2)) * pkin(3)) / t858;
t1052 = ((-t742 * t999 + t813 * t969) * t994 + t958) / t637 ^ 2;
t1042 = t776 * t866;
t1036 = t811 * t858;
t1034 = t813 * t858;
t1026 = t1106 * t643;
t879 = t881 * pkin(3);
t979 = pkin(2) * t1030;
t875 = t759 * t879 * t979;
t981 = pkin(1) * t1042;
t902 = t981 * t999;
t888 = -0.24e2 * t902;
t1025 = t772 * t888 - 0.24e2 * t875;
t1023 = -0.4e1 * t1105 * t699;
t998 = qJD(2) * t812;
t987 = -0.2e1 * t1042;
t982 = V_base(6) * pkin(1) + V_base(2);
t980 = pkin(2) * t1048;
t967 = pkin(3) * t997;
t961 = t811 * t1041;
t944 = t775 * t999;
t943 = qJD(2) * t1040;
t939 = t1063 * t1061;
t648 = -0.4e1 * t881 * t939;
t687 = -t968 + t969;
t930 = pkin(1) * t968;
t737 = -0.2e1 * t930;
t773 = t811 ^ 2;
t896 = t779 * t862 * t945;
t886 = -0.24e2 * t896;
t890 = t758 * t904;
t891 = t758 * t905;
t892 = t752 * t749 - t750 * t967;
t908 = t851 * t942;
t898 = -0.24e2 * t908;
t899 = t749 * t1043;
t920 = pkin(3) * t943;
t901 = t811 * t920;
t907 = -0.6e1 * t930;
t909 = t997 * t1032;
t925 = -0.2e1 * t942;
t933 = pkin(3) * t979;
t935 = t1053 / 0.2e1;
t938 = -((0.8e1 * t813 * t743 * t920 + t773 * t943 * t1078 + t886 * t1059 + ((pkin(3) * t768 + t1082 * t702 + t987) * t997 + (t812 * t737 + (t914 + 0.2e1 * t1030) * t1001) * pkin(2)) * t986 + 0.8e1 * t618 * t912 + 0.4e1 * ((t1001 * t1084 + t811 * t924) * t1061 + ((t749 - t967) * t766 - t728 * t968) * t1104 + t892) * t980 + t899 * t1087 + t888 * t767 + (t1001 * t748 + 0.4e1 * t866 * t909) * t776 + t663 * t925 - 0.4e1 * t892 * t990 + 0.4e1 * t659 * t930 - t679 * t749 + t686 * t967 - 0.4e1 * t1093 * t615) * t858 + t593 * t935 + ((t737 - 0.8e1 * t908) * t928 + 0.12e2 * (-0.12e2 * t902 + 0.6e1 * (-0.4e1 / 0.9e1 * t941 - 0.4e1 / 0.9e1 * t940) * t933 - 0.12e2 * t660 * t908 + t648 - 0.4e1 * t669 * t930 - 0.2e1 * t891 - 0.2e1 * t890) * t1039 - 0.24e2 * t602 * t912 - 0.6e1 * (t747 * t925 + (-t1021 * t813 - t667 * t999) * pkin(1) * t993 + t1023) * t980 + t636 * t898 + t1023 * t959 + t614 * t907 + t1025 + (0.8e1 * t901 + 0.24e2 * t896) * t680 + t1106 * t641 + 0.6e1 * t1093 * t603) * t700 + t586 * t687) / t580 ^ 2 / 0.4e1;
t757 = V_base(6) + qJD(1);
t756 = V_base(6) + qJD(2);
t913 = t778 * t863 * t1002;
t910 = t849 * t944;
t900 = t749 * t1045;
t893 = t632 * t749 - t650 * t967;
t889 = -0.48e2 * t902;
t878 = 0.4e1 * t881;
t877 = -0.2e1 * t879;
t876 = t878 * t766;
t762 = Icges(2,4) * t814;
t761 = Icges(4,4) * t813;
t724 = rSges(2,1) * t814 - rSges(2,2) * t812;
t723 = rSges(4,1) * t813 - rSges(4,2) * t811;
t722 = rSges(2,1) * t812 + rSges(2,2) * t814;
t721 = rSges(4,1) * t811 + rSges(4,2) * t813;
t715 = Icges(2,1) * t814 - t1057;
t714 = Icges(2,1) * t812 + t762;
t713 = Icges(4,1) * t813 - t1056;
t712 = Icges(4,1) * t811 + t761;
t711 = -Icges(2,2) * t812 + t762;
t710 = Icges(2,2) * t814 + t1057;
t709 = -Icges(4,2) * t811 + t761;
t708 = Icges(4,2) * t813 + t1056;
t690 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t689 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t688 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t656 = V_base(5) * rSges(2,3) - t722 * t757 + V_base(1);
t655 = V_base(5) * rSges(4,3) - t721 * t756 + V_base(1);
t654 = -V_base(4) * rSges(2,3) + t724 * t757 + V_base(2);
t640 = -V_base(4) * rSges(4,3) + t723 * t756 + t982;
t638 = t722 * V_base(4) - t724 * V_base(5) + V_base(3);
t631 = t1091 * t685;
t630 = t1091 * t684;
t628 = t721 * V_base(4) + V_base(3) + (-pkin(1) - t723) * V_base(5);
t589 = 0.1e1 / t590 ^ 2;
t582 = Icges(5,4) * t584;
t577 = t670 * t935 + (-t738 * t813 + (t658 * t811 + t1034) * qJD(2)) * t742 + (t903 + 0.4e1 * t911 - 0.2e1 * t1051) * pkin(3) + ((t730 * t813 - t1036) * t998 + (-t1035 * t658 + t1037 * t730 + t685 * t858) * qJD(1) + t742 * t876 + (t1002 * t732 - t1062 * t881) * t1104) * pkin(2);
t576 = t668 * t935 + (-t1036 * t742 + t619 * t813) * qJD(2) + (t673 * t813 - t681 * t999) * t811 * t994 + ((-t1034 + (-t677 - 0.8e1 * t975) * t811) * t998 + ((t622 - t1036) * t814 + (t1034 + (t1104 * t742 + t677) * t811 + (-pkin(3) + t1063 + 0.2e1 * t1060) * t1099) * t812) * qJD(1)) * pkin(2);
t574 = rSges(5,1) * t584 - rSges(5,2) * t583;
t573 = rSges(5,1) * t583 + rSges(5,2) * t584;
t572 = Icges(5,1) * t584 - t1055;
t571 = Icges(5,1) * t583 + t582;
t570 = -Icges(5,2) * t583 + t582;
t569 = Icges(5,2) * t584 + t1055;
t565 = V_base(3) + (t573 + t1059) * V_base(4) + (-t574 - t743) * V_base(5);
t560 = t756 + ((-t577 * t592 * t589 + t1085 * t576) * t634 + (t590 * t589 - t1085) * t592 * t1052) / (t589 * t592 ^ 2 + 0.1e1) * t637;
t559 = V_base(5) * rSges(5,3) - t1059 * t756 - t560 * t573 + V_base(1);
t558 = -V_base(4) * rSges(5,3) + t560 * t574 + t756 * t766 + t982;
t557 = (t923 * t967 - 0.24e2 * pkin(3) * (-0.8e1 / 0.3e1 * t908 + t737) * t961 + 0.96e2 * t675 * t913 * t1059 + (t1095 * t749 + 0.2e1 * (pkin(3) * t694 - t701 * t811 + t987) * t997 + (t813 * t1100 + (t1082 * t729 + 0.4e1 * t1060) * qJD(2)) * t1075) * t927 + (-0.8e1 / 0.3e1 * t900 + (-t691 * t967 + t692 * t749) * t983 + t645 * t749 - t646 * t967 / 0.2e1 + (0.32e2 / 0.3e1 * t1044 * t767 + t629 * t813 * t1078) * t999 + (0.4e1 * t695 * t813 * t929 + ((0.12e2 * t773 * t776 - 0.4e1 * t774) * t849 + (-0.2e1 * t1060 * t751 + t1083 * t642) * pkin(3)) * qJD(2)) * pkin(1)) * t971 + 0.24e2 * t596 * t912 + 0.6e1 * ((-0.4e1 * (pkin(3) * t909 + t725 * t1100) * t776 + 0.8e1 * t651 * t942) * t851 + (-0.8e1 * t899 + (-t697 * t967 + t698 * t749) * t992 + (t1086 * t633 + 0.24e2 * t776 * t976) * t999) * pkin(1) + t893) * t980 + t900 * t1077 - 0.96e2 * t759 * t910 * t767 + (t696 * t749 + t966 * t997) * t931 + t649 * t889 + (t967 * t1101 + t647 * t749) * t984 + 0.12e2 * t607 * t908 - 0.6e1 * t893 * t990 + 0.6e1 * t606 * t930 - t609 * t749 + t616 * t967 + (-0.8e1 * t901 + t886) * t608 - 0.6e1 * t1093 * t597) * t858 + t581 * t935 + (0.16e2 * (t1081 * t755 - 0.32e2 * t1044 - t741 - 0.48e2 * t981) * qJD(2) * t961 - 0.64e2 * t626 * t913 + (t648 + (t666 * t1081 + (t1086 * t672 - 0.24e2 * t1042) * pkin(1)) * t999 + (t1033 * t877 - t1042 * t878) * pkin(2)) * t922 + 0.24e2 * (-0.32e2 * t760 * t910 + (-0.2e1 / 0.3e1 * t941 - 0.2e1 / 0.3e1 * t940) * t931 * t1074 + t674 * t889 - 0.28e2 * t611 * t908 + 0.6e1 * (-0.8e1 / 0.3e1 * t941 - 0.8e1 / 0.3e1 * t940) * t727 * t939 + t613 * t907 + 0.4e1 * (-t891 - t890) * t735 - 0.64e2 / 0.3e1 * t881 * t726 * t933) * t1039 - 0.48e2 * t594 * t912 - 0.8e1 * (t639 * t898 + 0.6e1 * (-pkin(2) * t735 * t876 - t617 * t999) * t1075 + t1025 + t1026) * t980 - 0.4e1 * t720 * t944 + 0.32e2 * t775 * t877 * t960 * t1076 - 0.96e2 * t644 * t902 - 0.96e2 * t735 * t875 - 0.48e2 * t605 * t908 + 0.8e1 * t1026 * t990 - 0.8e1 * t601 * t930 + (0.32e2 * t901 + 0.96e2 * t896) * t604 + 0.8e1 * t1093 * t595 - 0.8e1 * t1105 * t610) * t700 + t575 * t687;
t554 = 0.1e1 / t556 ^ 2;
t547 = (-t883 * t1052 + (t577 * t1064 + t590 * t1053 / 0.8e1 + t557 * t936 + (t1071 * t576 + t592 * t938) * t566) * t634) * t1031;
t546 = (-t882 * t1052 + (t557 * t937 + t576 * t1065 - t592 * t1053 / 0.8e1 + (t1071 * t577 + t590 * t938) * t566) * t634) * t1031;
t545 = rSges(3,1) * t926 + rSges(3,2) * t550;
t544 = -rSges(3,1) * t550 + rSges(3,2) * t926;
t543 = Icges(3,1) * t926 + t1107;
t542 = -Icges(3,1) * t550 + t1096;
t541 = Icges(3,2) * t550 + t1096;
t540 = Icges(3,2) * t926 - t1107;
t537 = t544 * V_base(4) - t545 * V_base(5) + V_base(3) + (t812 * V_base(4) - t814 * V_base(5)) * pkin(2);
t536 = ((t546 * t684 + t547 * t685 - t562 * t631 + t563 * t630) / t556 - (t546 * t685 - t547 * t684 + t562 * t630 + t563 * t631) * t555 * t554) / (t554 * t555 ^ 2 + 0.1e1) + t757;
t535 = V_base(5) * rSges(3,3) - t536 * t544 - t757 * t767 + V_base(1);
t534 = -V_base(4) * rSges(3,3) + t1062 * t757 + t536 * t545 + V_base(2);
t1 = m(1) * (t688 ^ 2 + t689 ^ 2 + t690 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t638 ^ 2 + t654 ^ 2 + t656 ^ 2) / 0.2e1 + Icges(2,3) * t757 ^ 2 / 0.2e1 + m(3) * (t534 ^ 2 + t535 ^ 2 + t537 ^ 2) / 0.2e1 + Icges(3,3) * t536 ^ 2 / 0.2e1 + m(4) * (t628 ^ 2 + t640 ^ 2 + t655 ^ 2) / 0.2e1 + Icges(4,3) * t756 ^ 2 / 0.2e1 + m(5) * (t558 ^ 2 + t559 ^ 2 + t565 ^ 2) / 0.2e1 + Icges(5,3) * t560 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (-Icges(3,5) * t550 + Icges(3,6) * t926) * t536 + (Icges(5,5) * t583 + Icges(5,6) * t584) * t560 + (Icges(4,5) * t811 + Icges(4,6) * t813) * t756 + (Icges(2,5) * t812 + Icges(2,6) * t814) * t757 + (Icges(1,2) / 0.2e1 + t710 * t1067 + t714 * t812 / 0.2e1 + t540 * t1097 - t542 * t550 / 0.2e1 + t708 * t1068 + t712 * t811 / 0.2e1 + t569 * t1070 + t571 * t583 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t814 - Icges(2,6) * t812) * t757 + (Icges(3,5) * t926 + Icges(3,6) * t550) * t536 + (Icges(4,5) * t813 - Icges(4,6) * t811) * t756 + (Icges(5,5) * t584 - Icges(5,6) * t583) * t560 + ((t711 + t714) * t814 + (t709 + t712) * t813 + (-t710 + t715) * t812 + (-t708 + t713) * t811 + (t570 + t571) * t584 + (-t569 + t572) * t583 + (t541 + t542) * t926 + (t540 - t543) * t550) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t711 * t812 / 0.2e1 + t715 * t1067 + t541 * t550 / 0.2e1 + t543 * t1097 + t709 * t1069 + t713 * t1068 - t570 * t583 / 0.2e1 + t572 * t1070) * V_base(4)) * V_base(4);
T = t1;
