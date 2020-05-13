% Calculate kinetic energy for
% fivebar1DE2
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
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1DE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1DE2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 04:27:15
% EndTime: 2020-04-27 04:27:45
% DurationCPUTime: 26.01s
% Computational Cost: add. (362011->737), mult. (1084078->1112), div. (2756->13), fcn. (134114->12), ass. (0->477)
t836 = pkin(5) ^ 2;
t835 = t836 ^ 2;
t839 = pkin(4) ^ 2;
t838 = t839 ^ 2;
t1004 = t835 - t838;
t850 = (pkin(1) ^ 2);
t834 = -2 * t850;
t803 = sin(qJ(2));
t806 = cos(qJ(1));
t995 = qJD(1) * t806;
t741 = pkin(2) * t995;
t922 = pkin(3) * t741;
t897 = t803 * t922;
t1069 = pkin(2) * pkin(3);
t804 = sin(qJ(1));
t805 = cos(qJ(2));
t991 = qJD(2) * t805;
t936 = t804 * t991;
t898 = t936 * t1069;
t1097 = -0.6e1 * t898 - 0.6e1 * t897;
t1096 = t898 + t897;
t1056 = pkin(2) * t806;
t758 = pkin(3) * t805;
t735 = pkin(1) + t758;
t996 = qJD(1) * t804;
t940 = t735 * t996;
t993 = qJD(2) * t803;
t962 = pkin(3) * t993;
t1086 = pkin(2) * t940 + t962 * t1056;
t1095 = 0.2e1 * pkin(1);
t987 = 0.4e1 * pkin(3);
t844 = pkin(3) ^ 2;
t848 = pkin(2) ^ 2;
t856 = t848 ^ 2;
t784 = -t836 / 0.3e1;
t1087 = t784 - t839 / 0.3e1;
t948 = t850 + t1087;
t909 = t844 + t948;
t792 = t839 / 0.3e1;
t833 = 2 * t850;
t950 = t836 / 0.3e1 + t792 + t833;
t1092 = 0.10e2 / 0.3e1 * t856 - 0.2e1 * (-t844 + t950) * t848 + t909 * t834;
t1091 = 0.2e1 * t741;
t759 = pkin(2) * t804;
t1090 = 0.2e1 * t759;
t1089 = 0.4e1 * t844;
t823 = 0.6e1 * t844;
t813 = 0.10e2 * t844;
t1055 = pkin(2) * t844;
t768 = t805 ^ 2;
t1053 = pkin(3) * t768;
t763 = -0.3e1 * t844 + t850;
t1024 = t844 * t768;
t977 = 0.4e1 * t1024;
t1088 = t763 + t977;
t935 = t803 * t995;
t874 = t935 + t936;
t832 = 3 * t850;
t1085 = t832 - t836 - t839;
t1084 = qJD(1) - qJD(2);
t1083 = -t835 / 0.6e1 + t838 / 0.6e1;
t1081 = 0.8e1 * pkin(1);
t1080 = -0.4e1 * pkin(3);
t1079 = -0.2e1 * t768;
t1078 = -0.4e1 * t803;
t1077 = -0.2e1 * t803;
t1076 = -0.2e1 * t805;
t1075 = -0.2e1 * t806;
t817 = -0.6e1 * t836;
t842 = t844 ^ 2;
t1074 = 0.4e1 * t842;
t1073 = -0.8e1 * t844;
t751 = t850 - t848 / 0.3e1;
t1072 = 0.24e2 * t751;
t1071 = pkin(1) * pkin(2);
t1070 = pkin(1) * pkin(3);
t1068 = -pkin(4) - pkin(5);
t1067 = -pkin(4) + pkin(5);
t734 = pkin(1) - t1056;
t1043 = t734 * t805;
t1001 = -t848 + t850;
t985 = pkin(1) * t1056;
t737 = -0.2e1 * t985;
t771 = t806 ^ 2;
t1033 = t771 * t848;
t979 = 0.2e1 * t1033;
t673 = t737 + t979 + t1001;
t1044 = t673 * t768;
t997 = t850 - t836;
t941 = t848 + t997;
t908 = t844 + t941;
t723 = -t839 + t908;
t670 = t737 + t723;
t1032 = t803 * t804;
t968 = pkin(3) * t1032;
t925 = pkin(2) * t968;
t725 = -0.2e1 * t925;
t649 = t725 + t670;
t1023 = t848 * t850;
t772 = t848 * t1089;
t738 = t772 - 0.4e1 * t1023;
t818 = 0.2e1 * t839;
t972 = pkin(2) * t1032;
t851 = sqrt(t738 * t771 + 0.4e1 * t723 * t985 - t842 - (t850 + (pkin(2) - t1067) * (pkin(2) + t1067)) * (t850 + (pkin(2) - t1068) * (pkin(2) + t1068)) + (t818 + t834 + 0.2e1 * t836 - 0.6e1 * t848 - 0.4e1 * t1044) * t844 + (-t1043 * t649 + t670 * t972) * t987);
t1059 = -t851 / 0.4e1;
t767 = t805 * t768;
t859 = pkin(3) * t844;
t1037 = t767 * t859;
t761 = t844 + t850;
t755 = t761 ^ 2;
t942 = t844 + t997;
t1040 = t755 * (-t839 + t942);
t1042 = t735 * t806;
t1022 = t850 * t836;
t1027 = (pkin(1) - pkin(3)) * (pkin(1) + pkin(3));
t1010 = 0.4e1 / 0.7e1 * t850 - t836 / 0.7e1;
t718 = t850 + t848 / 0.4e1 + t844 / 0.4e1 - t836 / 0.8e1;
t849 = t850 ^ 2;
t603 = -0.32e2 / 0.21e2 * t718 * t925 + t856 / 0.7e1 + (0.16e2 / 0.21e2 * t844 + t1010) * t848 + t842 / 0.7e1 + t1010 * t844 + t849 - 0.3e1 / 0.7e1 * t1022 + t1004 / 0.42e2;
t1060 = 0.4e1 / 0.3e1 * t848;
t801 = t848 / 0.2e1;
t1011 = t801 + t850;
t783 = -t836 / 0.4e1;
t798 = t844 / 0.3e1;
t719 = t783 + t798 + t1011;
t786 = -0.2e1 / 0.3e1 * t836;
t797 = 0.4e1 / 0.3e1 * t844;
t605 = -0.8e1 / 0.3e1 * t719 * t925 + t856 / 0.3e1 + (t797 + t784) * t848 + t849 - t842 / 0.3e1 + (t1060 + 0.2e1 / 0.3e1 * t844 + t786) * t850 + t1004 / 0.18e2;
t1002 = t842 + t849;
t1005 = t833 - t836;
t663 = t1005 * t844 + t1002 - t1022 - t1083;
t799 = t844 / 0.2e1;
t666 = -0.2e1 / 0.3e1 * t925 + t850 + t799 + t783;
t831 = 4 * t850;
t745 = (t831 + t836) * t844;
t750 = -t844 / 0.3e1 + t850;
t752 = t850 - 0.2e1 / 0.3e1 * t848;
t762 = -t844 + t850;
t785 = -t836 / 0.2e1;
t1000 = t848 + t850;
t943 = t844 + t1000;
t727 = t785 + t943;
t912 = -0.4e1 * t925;
t890 = t727 * t912;
t966 = 0.16e2 * t1037;
t766 = t768 ^ 2;
t1039 = t766 * t842;
t983 = 0.8e1 * t1039;
t586 = t752 * t983 + 0.14e2 * t603 * t1024 + t750 * t890 + t762 * t856 + (-t1022 + t745 - 0.10e2 / 0.3e1 * t842 + (2 * t849)) * t848 + t663 * t1027 + (0.6e1 * t605 * t758 + t666 * t966) * pkin(1);
t1009 = 0.15e2 * t844 + t832;
t1014 = t835 / 0.2e1 - t838 / 0.2e1;
t787 = -0.3e1 / 0.2e1 * t836;
t855 = pkin(2) * t848;
t845 = t855 ^ 2;
t1018 = t761 * ((t787 + t833) * t844 - 0.3e1 / 0.2e1 * t1022 + t1002 + t1014) + t845;
t774 = 0.10e2 / 0.3e1 * t844;
t873 = t663 + t856;
t635 = (t774 + t1005) * t848 + t873;
t811 = 0.15e2 * t842;
t814 = 18 * t850;
t828 = 3 * t849;
t911 = -0.3e1 * t1022 + t828 + t1014;
t593 = -0.6e1 * t635 * t925 + t1018 + (t811 + (t814 - 0.9e1 * t836) * t844 + t911) * t848 + (t787 + t1009) * t856;
t609 = t890 + (t823 + t1005) * t848 + t873;
t1026 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t674 = t751 * t725;
t631 = t1026 * t727 + t674;
t764 = -0.3e1 * t848 + t850;
t982 = 0.8e1 * t1037;
t927 = pkin(1) * t982;
t710 = t764 * t927;
t984 = pkin(1) * t758;
t954 = 0.6e1 * t984;
t964 = 0.12e2 * t1024;
t587 = t609 * t954 + t631 * t964 + t593 + t710;
t815 = -0.2e1 * t836;
t829 = 8 * t850;
t970 = t859 * t759;
t658 = t970 * t1078 + t772 + t1074 + (t815 + t829) * t844;
t877 = -t925 + t1011;
t664 = t783 - t844 + t877;
t986 = 0.4e1 * t758;
t596 = t725 * t1027 + t658 * t768 + t727 * t763 + (t664 * t986 + t982) * pkin(1);
t999 = t849 - t842;
t597 = t751 * t890 - t845 + (-t774 - t997) * t856 + (t745 + t999 + t1083) * t848 + t850 * t663;
t816 = -0.5e1 * t836;
t821 = 0.7e1 * t842;
t602 = (t787 + t832 + 0.7e1 * t844) * t856 + (t821 + (t816 + (10 * t850)) * t844 + t911) * t848 + t1018;
t1021 = t850 * t844;
t852 = pkin(1) * t850;
t733 = -0.12e2 * pkin(1) * t859 + t852 * t987;
t747 = -0.8e1 * t842 + 0.12e2 * t1021;
t924 = pkin(1) * t966;
t618 = t733 * t805 + t747 * t768 + t1002 - 0.6e1 * t1021 + t924 + t983;
t636 = t1026 * t725 + t727 * t764;
t998 = t849 + t856;
t712 = 0.16e2 * (t998 - 0.6e1 * t1023) * t842;
t731 = t839 + t942;
t756 = -0.30e2 * t836 + (60 * t850);
t769 = t771 ^ 2;
t880 = t1004 + (6 * t849) - 0.6e1 * t1022;
t770 = t806 * t771;
t1034 = t770 * t855;
t957 = t735 * t1034;
t915 = -0.32e2 * t957;
t568 = t596 * t915 + t712 * t766 + 0.24e2 * t597 * t1024 + (t815 + t831 + 0.28e2 * t844) * t845 + t731 * t1040 + (t1004 * t833 + 0.24e2 * t586 * t771 + t756 * t842 + t817 * t849 + t880 * t823 + 0.4e1 * t852 ^ 2 + 0.28e2 * t859 ^ 2) * t848 + 0.8e1 * (-t1042 * t587 - t602 * t968) * pkin(2) + (0.32e2 * t1037 * t636 + 0.8e1 * t593 * t758) * pkin(1) + (0.16e2 * t618 * t769 + t756 * t844 + 0.70e2 * t842 + t856 + t880) * t856;
t1057 = pkin(1) * t805;
t1038 = t767 * t842;
t1054 = pkin(3) * t803;
t782 = -t836 / 0.6e1;
t1013 = t782 - t839 / 0.6e1;
t949 = t850 + t1013;
t683 = t1060 + t799 + t949;
t910 = t801 + t949;
t684 = t797 + t910;
t621 = -t1054 * t683 + t684 * t759;
t687 = t844 + t910;
t743 = 0.2e1 * t848 + t762;
t959 = -t1054 / 0.2e1;
t634 = t687 * t759 + t743 * t959;
t637 = t848 * t762 - 0.5e1 / 0.3e1 * t842 + t950 * t844 + t850 * t948;
t775 = -0.20e2 / 0.3e1 * t844;
t791 = 0.2e1 / 0.3e1 * t839;
t951 = 0.2e1 / 0.3e1 * t836 + t791 + t831;
t952 = 0.4e1 / 0.3e1 * t836 + 0.4e1 / 0.3e1 * t839 + t834;
t638 = -t856 + (t775 + t951) * t848 - 0.3e1 * t842 + t952 * t844 + t849;
t989 = 0.4e1 * pkin(1);
t588 = t621 * t977 + t638 * t959 + (-0.8e1 / 0.3e1 * t1039 + t637) * t759 + (-t1038 * t803 + t634 * t758) * t989;
t1008 = t813 + t833;
t796 = -0.2e1 / 0.3e1 * t839;
t1012 = t786 + t796;
t822 = 0.5e1 * t842;
t1007 = t815 - 0.2e1 * t839;
t830 = 6 * t850;
t944 = t830 + t1007;
t624 = t856 + (t1008 + t1012) * t848 + t822 + t944 * t844 + t850 * (t850 + t1012);
t826 = 0.5e1 * t856;
t945 = t786 + t761;
t642 = t826 + (t813 + t944) * t848 + t761 * (t796 + t945);
t598 = -t1054 * t642 + t624 * t759;
t827 = 0.3e1 * t848;
t689 = t827 + t909;
t824 = 0.3e1 * t844;
t744 = t824 + t1000;
t690 = t744 + t1087;
t625 = -t1054 * t689 + t690 * t759;
t717 = t799 + t848 + t1013;
t643 = t1026 * t1054 + t717 * t1090;
t589 = -0.4e1 * t643 * t1024 + (t625 * t986 - 0.8e1 * t767 * t970) * pkin(1) + t598;
t639 = -0.3e1 * t856 + (t775 + t952) * t848 + t951 * t844 + t999;
t599 = t1054 * t1092 + t639 * t759;
t947 = t785 - t839 / 0.2e1 + t850;
t686 = 0.3e1 / 0.2e1 * t848 + t824 + t947;
t971 = t804 * t1055;
t878 = -t803 * t859 + t971;
t693 = 0.4e1 * t878;
t721 = t759 + 0.2e1 * t1054;
t988 = 0.2e1 * pkin(3);
t600 = t763 * t759 + t693 * t768 + (t1057 * t721 + t686 * t803) * t988;
t1006 = t816 - 0.5e1 * t839;
t728 = t1085 * t813;
t601 = t845 + (0.21e2 * t844 + t1085) * t856 + (t1007 * t850 + t728 + t828 + 0.35e2 * t842) * t848 + (t821 + (t829 + t1006) * t844 + t850 * (-t839 + t997)) * t761;
t608 = 0.7e1 * t845 + (0.35e2 * t844 + (15 * t850) + t1006) * t856 + (0.21e2 * t842 + t728 + (9 * t849) + (t817 - 0.6e1 * t839) * t850) * t848 + t1040;
t688 = t827 + 0.3e1 / 0.2e1 * t844 + t947;
t960 = pkin(3) * t764 / 0.2e1;
t641 = t688 * t759 + t803 * t960;
t1035 = t769 * t856;
t736 = 0.2e1 * t984;
t667 = 0.4e1 / 0.3e1 * t1024 + t736 + t750;
t916 = -0.24e2 * t667 * t1035;
t920 = 0.8e1 * t957;
t965 = -0.12e2 * t1033;
t978 = -0.6e1 * t1024;
t574 = t641 * t924 + t600 * t920 + t588 * t965 + t599 * t978 + (-0.6e1 * t598 * t1057 + (t608 + t916) * t803) * pkin(3) + (0.6e1 * t589 * t1042 + (t1039 * t1072 - t601) * t804) * pkin(2);
t692 = t734 + t758;
t559 = t568 * t692 + t574 * t851;
t669 = t737 + t839 + t908;
t611 = t669 * t734 + 0.2e1 * t673 * t758;
t614 = t669 * t805 + (0.4e1 * t768 - 0.2e1) * t734 * pkin(3);
t724 = -pkin(3) + t972;
t660 = -t724 + t1043;
t584 = t611 * t803 + t614 * t759 + t660 * t851;
t652 = -0.4e1 / 0.9e1 * t925 + t850 + t848 / 0.3e1 + t798 + t839 / 0.9e1 - t836 / 0.9e1;
t661 = t839 / 0.6e1 + t782 + t877;
t691 = t784 + t792 + t943;
t594 = t750 * t725 + 0.6e1 * t652 * t1024 + t691 * t1027 + (t661 * t758 + t1037) * t989;
t1016 = t761 * (t791 + t945) + t856;
t946 = t786 + t791 + t833;
t606 = t691 * t912 + t1016 + (t823 + t946) * t848;
t659 = t725 + t691;
t739 = t1001 * t1089;
t595 = 0.4e1 * t659 * t984 + t739 * t768 + t606;
t628 = t1026 * t691 + t674;
t633 = (t774 + t946) * t848 + t1016;
t672 = t736 + t1088;
t921 = -0.8e1 * t957;
t580 = t672 * t921 + t710 + t628 * t964 + t606 * t954 + t845 + (-t836 + t839 + t1009) * t856 + t755 * t731 + (0.12e2 * t594 * t771 + t811 + (t814 + t817 + 0.6e1 * t839) * t844 + t828 + (t815 + t818) * t850) * t848 + 0.6e1 * (-t1042 * t595 - t633 * t968) * pkin(2);
t742 = t827 + t761;
t1041 = t742 * t803;
t711 = t744 * t759;
t720 = t759 - t1054;
t607 = t971 * t1079 + t711 + (0.2e1 * t1057 * t720 - t1041) * pkin(3);
t694 = 0.2e1 * t878;
t760 = 0.2e1 * t844 + t848;
t907 = t762 + t736;
t610 = t1054 * t760 + t694 * t768 + t759 * t907;
t651 = -pkin(3) * t1041 + t711;
t740 = pkin(2) * t1074 + 0.8e1 * t844 * t855;
t955 = t859 * t1026;
t655 = t740 * t804 + 0.4e1 * t803 * t955;
t671 = t1008 * t848 + 0.6e1 * t1021 + t822 + t998;
t678 = t826 + (t813 + t830) * t848 + t755;
t980 = -0.4e1 * t1033;
t585 = t610 * t980 + t655 * t768 + (-0.4e1 * t651 * t1057 + (t678 + t920) * t803) * pkin(3) + (0.4e1 * t607 * t1042 + (-t671 + t927) * t804) * pkin(2);
t573 = t580 * t692 + t585 * t851;
t1065 = 0.1e1 / t573 / 0.4e1;
t722 = t824 + t839 + t941;
t650 = t722 + t737 + t912;
t1029 = t805 * t804;
t662 = pkin(2) * t1029 + t734 * t803;
t582 = -t650 * t1043 + t662 * t851 + (pkin(1) * t1075 * t724 + t1032 * t722) * pkin(2) + (-t731 - t827 + t979 - 0.2e1 * t1044) * pkin(3);
t932 = t582 * t1065;
t875 = t1059 * t584 + t559 * t932;
t1025 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t969 = pkin(3) * t1043;
t629 = t725 + t737 + t943 + 0.2e1 * t969;
t626 = 0.1e1 / t629;
t958 = t626 * t1025;
t555 = t875 * t958;
t1058 = t851 / 0.4e1;
t931 = t584 * t1065;
t876 = t1058 * t582 + t559 * t931;
t556 = t876 * t958;
t1031 = t803 * t806;
t676 = -t1029 + t1031;
t677 = t805 * t806 + t1032;
t549 = t555 * t676 + t556 * t677;
t550 = t555 * t677 - t556 * t676;
t546 = qJ(1) + atan2(t549, t550);
t545 = cos(t546);
t1066 = -t545 / 0.2e1;
t1063 = t626 / 0.2e1;
t840 = 0.1e1 / pkin(4);
t928 = t840 * t1063;
t578 = qJ(2) + atan2(t584 * t928, t582 * t928);
t577 = cos(t578);
t1064 = t577 / 0.2e1;
t1062 = t805 / 0.2e1;
t1061 = t806 / 0.2e1;
t1052 = Icges(2,4) * t804;
t544 = sin(t546);
t1051 = Icges(3,4) * t544;
t1050 = Icges(3,4) * t545;
t1049 = Icges(4,4) * t803;
t576 = sin(t578);
t1048 = Icges(5,4) * t576;
t963 = pkin(2) * t996;
t730 = t963 * t1095;
t994 = qJD(1) * t848;
t901 = t806 * t804 * t994;
t896 = -0.4e1 * t901;
t665 = t730 + t896;
t1045 = t665 * t768;
t937 = t803 * t991;
t903 = t673 * t937;
t917 = 0.4e1 * t936;
t1015 = 0.2e1 * t1096;
t953 = t730 - t1015;
t1047 = ((0.8e1 * t903 - 0.4e1 * t1045) * t844 + (-0.4e1 * t1071 * t723 + t1075 * t738) * t996 + (t803 * t804 ^ 2 * t994 * t1081 + 0.4e1 * (t649 * t993 - t805 * t953) * t734 + (t670 * t917 + 0.4e1 * (-t1029 * t649 + t1031 * t670) * qJD(1)) * pkin(2)) * pkin(3)) / t851;
t1046 = ((-t734 * t993 + t805 * t963) * t988 + t953) / t629 ^ 2;
t1036 = t768 * t859;
t1030 = t803 * t851;
t1028 = t805 * t851;
t1020 = t1097 * t635;
t872 = t874 * pkin(3);
t973 = pkin(2) * t1024;
t868 = t751 * t872 * t973;
t975 = pkin(1) * t1036;
t895 = t975 * t993;
t881 = -0.24e2 * t895;
t1019 = t764 * t881 - 0.24e2 * t868;
t1017 = -0.4e1 * t1096 * t691;
t992 = qJD(2) * t804;
t981 = -0.2e1 * t1036;
t976 = V_base(6) * pkin(1) + V_base(2);
t974 = pkin(2) * t1042;
t961 = pkin(3) * t991;
t956 = t803 * t1035;
t939 = t767 * t993;
t938 = qJD(2) * t1034;
t934 = t1057 * t1055;
t640 = -0.4e1 * t874 * t934;
t679 = -t962 + t963;
t923 = pkin(1) * t962;
t729 = -0.2e1 * t923;
t765 = t803 ^ 2;
t889 = t771 * t855 * t940;
t879 = -0.24e2 * t889;
t883 = t750 * t897;
t884 = t750 * t898;
t885 = t744 * t741 - t742 * t961;
t902 = t844 * t937;
t891 = -0.24e2 * t902;
t892 = t741 * t1037;
t913 = pkin(3) * t938;
t894 = t803 * t913;
t900 = -0.6e1 * t923;
t904 = t991 * t1026;
t918 = -0.2e1 * t937;
t926 = pkin(3) * t973;
t930 = t1047 / 0.2e1;
t933 = -((0.8e1 * t805 * t735 * t913 + t765 * t938 * t1073 + t879 * t1054 + ((pkin(3) * t760 + t1077 * t694 + t981) * t991 + (t804 * t729 + (t907 + 0.2e1 * t1024) * t995) * pkin(2)) * t980 + 0.8e1 * t610 * t901 + 0.4e1 * ((t1079 * t995 + t803 * t917) * t1055 + ((t741 - t961) * t758 - t720 * t962) * t1095 + t885) * t974 + t892 * t1081 + t881 * t759 + (t740 * t995 + 0.4e1 * t859 * t904) * t768 + t655 * t918 - 0.4e1 * t885 * t984 + 0.4e1 * t651 * t923 - t671 * t741 + t678 * t961 - 0.4e1 * t1086 * t607) * t851 + t585 * t930 + t679 * t580 + t692 * ((t729 - 0.8e1 * t902) * t921 + 0.12e2 * (-0.12e2 * t895 + 0.6e1 * (-0.4e1 / 0.9e1 * t936 - 0.4e1 / 0.9e1 * t935) * t926 - 0.12e2 * t652 * t902 + t640 - 0.4e1 * t661 * t923 - 0.2e1 * t884 - 0.2e1 * t883) * t1033 - 0.24e2 * t594 * t901 - 0.6e1 * (t739 * t918 + (-t1015 * t805 - t659 * t993) * pkin(1) * t987 + t1017) * t974 + t628 * t891 + t1017 * t954 + t606 * t900 + t1019 + (0.8e1 * t894 + 0.24e2 * t889) * t672 + t1097 * t633 + 0.6e1 * t1086 * t595)) / t573 ^ 2 / 0.4e1;
t929 = -t1046 / 0.2e1;
t749 = V_base(6) + qJD(1);
t748 = V_base(6) + qJD(2);
t906 = t770 * t856 * t996;
t905 = t842 * t939;
t893 = t741 * t1039;
t886 = t624 * t741 - t642 * t961;
t882 = -0.48e2 * t895;
t871 = 0.4e1 * t874;
t870 = -0.2e1 * t872;
t869 = t871 * t758;
t754 = Icges(2,4) * t806;
t753 = Icges(4,4) * t805;
t716 = rSges(2,1) * t806 - rSges(2,2) * t804;
t715 = rSges(4,1) * t805 - rSges(4,2) * t803;
t714 = rSges(2,1) * t804 + rSges(2,2) * t806;
t713 = rSges(4,1) * t803 + rSges(4,2) * t805;
t707 = Icges(2,1) * t806 - t1052;
t706 = Icges(2,1) * t804 + t754;
t705 = Icges(4,1) * t805 - t1049;
t704 = Icges(4,1) * t803 + t753;
t703 = -Icges(2,2) * t804 + t754;
t702 = Icges(2,2) * t806 + t1052;
t701 = -Icges(4,2) * t803 + t753;
t700 = Icges(4,2) * t805 + t1049;
t682 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t681 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t680 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t648 = V_base(5) * rSges(2,3) - t714 * t749 + V_base(1);
t647 = V_base(5) * rSges(4,3) - t713 * t748 + V_base(1);
t646 = -V_base(4) * rSges(2,3) + t716 * t749 + V_base(2);
t632 = -V_base(4) * rSges(4,3) + t715 * t748 + t976;
t630 = t714 * V_base(4) - t716 * V_base(5) + V_base(3);
t623 = t1084 * t677;
t622 = t1084 * t676;
t620 = t713 * V_base(4) + V_base(3) + (-pkin(1) - t715) * V_base(5);
t581 = 0.1e1 / t582 ^ 2;
t575 = Icges(5,4) * t577;
t570 = t662 * t930 + (-t730 * t805 + (t650 * t803 + t1028) * qJD(2)) * t734 + (t896 + 0.4e1 * t903 - 0.2e1 * t1045) * pkin(3) + ((t722 * t805 - t1030) * t992 + (-t1029 * t650 + t1031 * t722 + t677 * t851) * qJD(1) + t734 * t869 + (-t1056 * t874 + t724 * t996) * t1095) * pkin(2);
t569 = t660 * t930 + (-t1030 * t734 + t611 * t805) * qJD(2) + (t665 * t805 - t673 * t993) * t803 * t988 + ((-t1028 + (-t669 - 0.8e1 * t969) * t803) * t992 + ((t614 - t1030) * t806 + (t1028 + (t1095 * t734 + t669) * t803 + (-pkin(3) + t1057 + 0.2e1 * t1053) * t1090) * t804) * qJD(1)) * pkin(2);
t567 = rSges(5,1) * t577 - rSges(5,2) * t576;
t566 = rSges(5,1) * t576 + rSges(5,2) * t577;
t565 = Icges(5,1) * t577 - t1048;
t564 = Icges(5,1) * t576 + t575;
t563 = -Icges(5,2) * t576 + t575;
t562 = Icges(5,2) * t577 + t1048;
t558 = V_base(3) + (t566 + t1054) * V_base(4) + (-t567 - t735) * V_base(5);
t554 = t748 + 0.2e1 * ((t1063 * t569 + t584 * t929) / t582 - (t1063 * t570 + t582 * t929) * t584 * t581) * pkin(4) / (t581 * t584 ^ 2 + 0.1e1) * t629 * t840;
t553 = V_base(5) * rSges(5,3) - t1054 * t748 - t554 * t566 + V_base(1);
t552 = -V_base(4) * rSges(5,3) + t554 * t567 + t748 * t758 + t976;
t551 = (t916 * t961 - 0.24e2 * (-0.8e1 / 0.3e1 * t902 + t729) * pkin(3) * t956 + 0.96e2 * t667 * t906 * t1054 + (t1088 * t741 + 0.2e1 * (pkin(3) * t686 - t693 * t803 + t981) * t991 + (t805 * t1091 + (t1077 * t721 + 0.4e1 * t1053) * qJD(2)) * t1070) * t920 + (-0.8e1 / 0.3e1 * t893 + (-t683 * t961 + t684 * t741) * t977 + t637 * t741 - t638 * t961 / 0.2e1 + (0.32e2 / 0.3e1 * t1038 * t759 + t621 * t805 * t1073) * t993 + (0.4e1 * t687 * t805 * t922 + ((0.12e2 * t765 * t768 - 0.4e1 * t766) * t842 + (-0.2e1 * t1053 * t743 + t1078 * t634) * pkin(3)) * qJD(2)) * pkin(1)) * t965 + 0.24e2 * t588 * t901 + 0.6e1 * ((-0.4e1 * (pkin(3) * t904 + t717 * t1091) * t768 + 0.8e1 * t643 * t937) * t844 + (-0.8e1 * t892 + (-t689 * t961 + t690 * t741) * t986 + (t1080 * t625 + 0.24e2 * t768 * t970) * t993) * pkin(1) + t886) * t974 + t893 * t1072 - 0.96e2 * t751 * t905 * t759 + (t688 * t741 + t960 * t991) * t924 + t641 * t882 + (t961 * t1092 + t639 * t741) * t978 + 0.12e2 * t599 * t902 - 0.6e1 * t886 * t984 + 0.6e1 * t598 * t923 - t601 * t741 + t608 * t961 + (-0.8e1 * t894 + t879) * t600 - 0.6e1 * t1086 * t589) * t851 + t574 * t930 + t679 * t568 + t692 * (0.16e2 * (t1076 * t747 - 0.32e2 * t1038 - t733 - 0.48e2 * t975) * qJD(2) * t956 - 0.64e2 * t618 * t906 + (t640 + (t658 * t1076 + (t1080 * t664 - 0.24e2 * t1036) * pkin(1)) * t993 + (t1027 * t870 - t1036 * t871) * pkin(2)) * t915 + 0.24e2 * (-0.32e2 * t752 * t905 + (-0.2e1 / 0.3e1 * t936 - 0.2e1 / 0.3e1 * t935) * t924 * t1069 + t666 * t882 - 0.28e2 * t603 * t902 + 0.6e1 * (-0.8e1 / 0.3e1 * t936 - 0.8e1 / 0.3e1 * t935) * t719 * t934 + t605 * t900 + 0.4e1 * (-t884 - t883) * t727 - 0.64e2 / 0.3e1 * t874 * t718 * t926) * t1033 - 0.48e2 * t586 * t901 - 0.8e1 * (t631 * t891 + 0.6e1 * (-pkin(2) * t727 * t869 - t609 * t993) * t1070 + t1019 + t1020) * t974 - 0.4e1 * t712 * t939 + 0.32e2 * t767 * t870 * t955 * t1071 - 0.96e2 * t636 * t895 - 0.96e2 * t727 * t868 - 0.48e2 * t597 * t902 + 0.8e1 * t1020 * t984 - 0.8e1 * t593 * t923 + (0.32e2 * t894 + 0.96e2 * t889) * t596 + 0.8e1 * t1086 * t587 - 0.8e1 * t1096 * t602);
t548 = 0.1e1 / t550 ^ 2;
t543 = (-t876 * t1046 + (t570 * t1058 + t582 * t1047 / 0.8e1 + t551 * t931 + (t1065 * t569 + t584 * t933) * t559) * t626) * t1025;
t542 = (-t875 * t1046 + (t551 * t932 + t569 * t1059 - t584 * t1047 / 0.8e1 + (t1065 * t570 + t582 * t933) * t559) * t626) * t1025;
t541 = -rSges(3,1) * t545 + rSges(3,2) * t544;
t540 = -rSges(3,1) * t544 - rSges(3,2) * t545;
t539 = -Icges(3,1) * t545 + t1051;
t538 = -Icges(3,1) * t544 - t1050;
t537 = Icges(3,2) * t544 - t1050;
t536 = -Icges(3,2) * t545 - t1051;
t533 = t540 * V_base(4) - t541 * V_base(5) + V_base(3) + (t804 * V_base(4) - t806 * V_base(5)) * pkin(2);
t532 = ((t542 * t676 + t543 * t677 - t555 * t623 + t556 * t622) / t550 - (t542 * t677 - t543 * t676 + t555 * t622 + t556 * t623) * t549 * t548) / (t548 * t549 ^ 2 + 0.1e1) + t749;
t531 = V_base(5) * rSges(3,3) - t532 * t540 - t749 * t759 + V_base(1);
t530 = -V_base(4) * rSges(3,3) + t1056 * t749 + t532 * t541 + V_base(2);
t1 = m(1) * (t680 ^ 2 + t681 ^ 2 + t682 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t630 ^ 2 + t646 ^ 2 + t648 ^ 2) / 0.2e1 + Icges(2,3) * t749 ^ 2 / 0.2e1 + m(3) * (t530 ^ 2 + t531 ^ 2 + t533 ^ 2) / 0.2e1 + Icges(3,3) * t532 ^ 2 / 0.2e1 + m(4) * (t620 ^ 2 + t632 ^ 2 + t647 ^ 2) / 0.2e1 + Icges(4,3) * t748 ^ 2 / 0.2e1 + m(5) * (t552 ^ 2 + t553 ^ 2 + t558 ^ 2) / 0.2e1 + Icges(5,3) * t554 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (-Icges(3,5) * t544 - Icges(3,6) * t545) * t532 + (Icges(5,5) * t576 + Icges(5,6) * t577) * t554 + (Icges(4,5) * t803 + Icges(4,6) * t805) * t748 + (Icges(2,5) * t804 + Icges(2,6) * t806) * t749 + (Icges(1,2) / 0.2e1 + t702 * t1061 + t706 * t804 / 0.2e1 + t536 * t1066 - t538 * t544 / 0.2e1 + t700 * t1062 + t704 * t803 / 0.2e1 + t562 * t1064 + t564 * t576 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t806 - Icges(2,6) * t804) * t749 + (-Icges(3,5) * t545 + Icges(3,6) * t544) * t532 + (Icges(4,5) * t805 - Icges(4,6) * t803) * t748 + (Icges(5,5) * t577 - Icges(5,6) * t576) * t554 + ((t703 + t706) * t806 + (t701 + t704) * t805 + (-t702 + t707) * t804 + (-t700 + t705) * t803 + (t563 + t564) * t577 + (-t562 + t565) * t576 + (-t537 - t538) * t545 + (t536 - t539) * t544) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t703 * t804 / 0.2e1 + t707 * t1061 + t537 * t544 / 0.2e1 + t539 * t1066 - t701 * t803 / 0.2e1 + t705 * t1062 - t563 * t576 / 0.2e1 + t565 * t1064) * V_base(4)) * V_base(4);
T = t1;
