% Calculate kinetic energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1TE_energykin_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:20:48
% EndTime: 2020-04-17 15:22:14
% DurationCPUTime: 84.54s
% Computational Cost: add. (1917539->548), mult. (2943133->941), div. (127608->22), fcn. (1845728->24), ass. (0->340)
t1019 = pkin(5) + pkin(6);
t1020 = pkin(5) - pkin(6);
t1018 = -pkin(6) - pkin(2);
t1006 = sin(qJ(2));
t1007 = sin(pkin(16));
t1009 = cos(qJ(2));
t1010 = cos(pkin(16));
t931 = t1006 * t1007 - t1009 * t1010;
t929 = pkin(5) * t931;
t923 = (-0.2e1 * t929 + pkin(1)) * pkin(1);
t910 = sqrt(-((-pkin(2) + t1019) * (pkin(2) + t1020) + t923) * ((pkin(5) - t1018) * (pkin(5) + t1018) + t923));
t911 = pkin(2) ^ 2;
t1021 = pkin(1) * pkin(5);
t932 = t1006 * t1010 + t1007 * t1009;
t868 = t932 * qJD(2);
t976 = t868 * t1021;
t988 = 0.4e1 * (t1019 * t1020 - t911 + t923) * t976 / t910;
t1028 = -t988 / 0.2e1;
t1027 = pkin(5) ^ 2;
t1026 = 0.1e1 / pkin(6);
t1025 = -0.2e1 * pkin(1);
t862 = t923 + t1027;
t1024 = 0.1e1 / t862;
t897 = sin(qJ(1));
t1023 = t897 ^ 2;
t900 = cos(qJ(1));
t1022 = t900 ^ 2;
t1017 = -pkin(8) - pkin(10);
t1016 = -pkin(8) + pkin(10);
t907 = pkin(3) ^ 2;
t983 = pkin(6) ^ 2 - t911;
t851 = t862 - t983;
t867 = -t929 + pkin(1);
t928 = t910 * t932;
t830 = -pkin(5) * t928 + t867 * t851;
t930 = pkin(5) * t932;
t831 = t851 * t930 + t867 * t910;
t899 = cos(qJ(3));
t908 = 0.1e1 / pkin(2);
t896 = sin(qJ(3));
t979 = t1024 / 0.2e1;
t966 = t896 * t979;
t980 = -t1024 / 0.2e1;
t820 = (t830 * t899 * t980 + t831 * t966) * t908;
t821 = (t831 * t899 * t979 + t830 * t966) * t908;
t888 = pkin(18) + pkin(19);
t884 = sin(t888);
t885 = cos(t888);
t802 = -t820 * t885 - t821 * t884;
t1002 = pkin(4) * t802;
t906 = pkin(4) ^ 2;
t986 = -0.2e1 * pkin(3) * t1002 + t906;
t785 = t907 + t986;
t783 = 0.1e1 / t785;
t1015 = t783 / 0.2e1;
t889 = sin(pkin(17));
t1014 = t889 / 0.2e1;
t894 = cos(pkin(18));
t1013 = -t894 / 0.2e1;
t1011 = cos(pkin(15));
t1008 = sin(pkin(15));
t869 = t931 * qJD(2);
t921 = t1028 * t932 + t869 * t910;
t812 = ((t1025 * t867 - t851) * t868 + t921) * pkin(5);
t859 = 0.1e1 / t862 ^ 2;
t965 = t859 * t976;
t957 = t830 * t965;
t924 = qJD(3) * t831 * t980 + t812 * t979 + t957;
t922 = (t1025 * t930 - t910) * t868;
t813 = t867 * t988 / 0.2e1 + (-t869 * t851 + t922) * pkin(5);
t956 = t831 * t965;
t925 = t956 + (qJD(3) * t830 + t813) * t979;
t760 = (t896 * t924 + t899 * t925) * t908;
t761 = (t896 * t925 - t899 * t924) * t908;
t759 = -t760 * t884 - t761 * t885;
t1005 = pkin(3) * t759;
t801 = t820 * t884 - t821 * t885;
t1003 = pkin(4) * t801;
t984 = pkin(8) ^ 2 - pkin(10) ^ 2;
t781 = t785 + t984;
t786 = -pkin(3) + t1002;
t779 = (pkin(3) - t1017) * (pkin(3) + t1017) + t986;
t780 = (pkin(3) - t1016) * (pkin(3) + t1016) + t986;
t909 = sqrt(-t780 * t779);
t756 = t1003 * t781 - t786 * t909;
t1004 = pkin(4) * t756;
t887 = qJD(2) * t897;
t872 = qJD(3) * t897 + t887;
t1001 = pkin(4) * t872;
t977 = -qJD(2) - qJD(3);
t873 = t977 * t900;
t1000 = pkin(4) * t873;
t920 = t862 + t983;
t919 = pkin(1) * t920;
t926 = pkin(1) * t931 - pkin(5);
t918 = t1026 * (-t910 * t926 + t919 * t932);
t913 = t918 * t979;
t917 = t1026 * (-pkin(1) * t928 - t920 * t926);
t915 = t1011 * t917;
t819 = t1008 * t913 + t915 * t979;
t998 = Icges(7,4) * t819;
t914 = t1008 * t917;
t822 = t1011 * t913 + t914 * t980;
t997 = Icges(7,4) * t822;
t975 = pkin(4) * t1005;
t996 = 0.2e1 * (t779 + t780) * t975 / t909;
t811 = ((-0.3e1 * t1027 + (0.4e1 * t929 - pkin(1)) * pkin(1) - t983) * t868 + t921) * pkin(1);
t814 = pkin(1) * t922 + t1028 * t926 - t869 * t919;
t816 = 0.1e1 / t819 ^ 2;
t916 = t918 * t1021;
t953 = t1026 * t1011 * t979;
t962 = t1026 * t1024 * t1008;
t987 = t859 * t868;
t748 = ((t814 * t953 - t811 * t962 / 0.2e1 + (t1011 * t916 - t1021 * t914) * t987) / t819 - (t811 * t953 + t814 * t962 / 0.2e1 + (t1008 * t916 + t1021 * t915) * t987) * t822 * t816) / (t816 * t822 ^ 2 + 0.1e1);
t995 = t748 * t897;
t994 = t748 * t900;
t993 = t759 * t909;
t890 = cos(pkin(17));
t992 = t783 * t890;
t903 = 0.1e1 / pkin(10);
t991 = t783 * t903;
t905 = 0.1e1 / pkin(8);
t990 = t783 * t905;
t989 = t801 * t909;
t974 = pkin(1) * t1009;
t985 = (t1022 + t1023) * qJD(2) * t974;
t893 = cos(pkin(19));
t968 = t893 * t980;
t891 = sin(pkin(19));
t969 = t891 * t979;
t817 = (t830 * t968 + t831 * t969) * t908;
t815 = 0.1e1 / t817 ^ 2;
t967 = t893 * t979;
t818 = (t830 * t969 + t831 * t967) * t908;
t747 = ((t812 * t969 + t813 * t967 + t891 * t957 + t893 * t956) / t817 - (t812 * t968 + t813 * t969 + t891 * t956 - t893 * t957) * t818 * t815) / (t815 * t818 ^ 2 + 0.1e1) * t908;
t744 = t747 * t897 + t887;
t982 = qJD(1) * t900;
t981 = qJD(2) * t900;
t978 = -qJD(2) - t747;
t782 = t785 - t984;
t787 = -pkin(3) * t802 + pkin(4);
t758 = -t760 * t885 + t761 * t884;
t971 = -t996 / 0.2e1;
t943 = -t758 * t909 + t801 * t971;
t711 = ((-0.2e1 * pkin(4) * t787 - t782) * t759 + t943) * pkin(3);
t712 = t787 * t996 / 0.2e1 - 0.2e1 * t907 * t759 * t1003 + (t758 * t782 - t993) * pkin(3);
t755 = -pkin(3) * t989 + t782 * t787;
t757 = pkin(3) * t782 * t801 + t787 * t909;
t734 = (-t755 * t890 / 0.2e1 + t757 * t1014) * t991;
t733 = 0.1e1 / t734 ^ 2;
t735 = (t757 * t890 / 0.2e1 + t755 * t1014) * t991;
t784 = 0.1e1 / t785 ^ 2;
t964 = t784 * t975;
t958 = t890 * t964;
t959 = t889 * t964;
t970 = t783 * t1014;
t666 = ((t712 * t992 / 0.2e1 + t757 * t958 + t711 * t970 + t755 * t959) / t734 - (-t711 * t992 / 0.2e1 - t755 * t958 + t712 * t970 + t757 * t959) * t735 * t733) / (t733 * t735 ^ 2 + 0.1e1) * t903;
t663 = t666 * t897 + t872;
t973 = Icges(3,4) * t1009;
t972 = Icges(3,4) * t1006;
t960 = pkin(1) * qJD(2) * t1006;
t664 = (-t666 + t977) * t900;
t955 = rSges(7,1) * t819 - rSges(7,2) * t822;
t954 = t900 * t960;
t952 = Icges(7,1) * t819 - t997;
t951 = -Icges(7,2) * t822 + t998;
t950 = Icges(7,5) * t819 - Icges(7,6) * t822;
t791 = -Icges(7,6) * t900 + t897 * t951;
t793 = -Icges(7,5) * t900 + t897 * t952;
t949 = t791 * t822 - t793 * t819;
t792 = Icges(7,6) * t897 + t900 * t951;
t794 = Icges(7,5) * t897 + t900 * t952;
t948 = -t792 * t822 + t794 * t819;
t806 = Icges(7,2) * t819 + t997;
t807 = Icges(7,1) * t822 + t998;
t947 = -t806 * t822 + t807 * t819;
t946 = (-t974 - pkin(13)) * t897;
t870 = t1006 * t896 - t1009 * t899;
t864 = t870 * t897;
t866 = t870 * t900;
t945 = -t1000 * t866 + t864 * t1001 + t985;
t935 = t1006 * t899 + t1009 * t896;
t944 = -t1000 * t935 - t954;
t942 = rSges(3,1) * t1009 - rSges(3,2) * t1006;
t941 = Icges(3,1) * t1009 - t972;
t940 = -Icges(3,2) * t1006 + t973;
t939 = Icges(3,5) * t1009 - Icges(3,6) * t1006;
t803 = t1006 * t817 + t1009 * t818;
t804 = -t1006 * t818 + t1009 * t817;
t855 = -Icges(3,6) * t900 + t897 * t940;
t857 = -Icges(3,5) * t900 + t897 * t941;
t938 = -t1006 * t855 + t1009 * t857;
t856 = Icges(3,6) * t897 + t900 * t940;
t858 = Icges(3,5) * t897 + t900 * t941;
t937 = -t1006 * t856 + t1009 * t858;
t875 = Icges(3,2) * t1009 + t972;
t876 = Icges(3,1) * t1006 + t973;
t936 = -t1006 * t875 + t1009 * t876;
t883 = pkin(13) * t982;
t934 = -t897 * t960 + t974 * t982 + t883;
t933 = -pkin(4) * t864 + t946;
t927 = qJD(1) * t866 * pkin(4) + t1001 * t935 + t934;
t898 = cos(qJ(4));
t895 = sin(qJ(4));
t892 = sin(pkin(18));
t879 = rSges(2,1) * t900 - rSges(2,2) * t897;
t878 = rSges(2,1) * t897 + rSges(2,2) * t900;
t877 = rSges(3,1) * t1006 + rSges(3,2) * t1009;
t874 = Icges(3,5) * t1006 + Icges(3,6) * t1009;
t865 = t935 * t900;
t863 = t935 * t897;
t861 = t897 * rSges(3,3) + t900 * t942;
t860 = -t900 * rSges(3,3) + t897 * t942;
t854 = Icges(3,3) * t897 + t900 * t939;
t853 = -Icges(3,3) * t900 + t897 * t939;
t849 = -rSges(4,1) * t935 + rSges(4,2) * t870;
t847 = -Icges(4,1) * t935 + Icges(4,4) * t870;
t846 = -Icges(4,4) * t935 + Icges(4,2) * t870;
t845 = -Icges(4,5) * t935 + Icges(4,6) * t870;
t844 = qJD(1) * t861 - t877 * t887 + t883;
t843 = -t877 * t981 + (-pkin(13) * t897 - t860) * qJD(1);
t842 = rSges(4,1) * t866 + rSges(4,2) * t865 + rSges(4,3) * t897;
t841 = rSges(4,1) * t864 + rSges(4,2) * t863 - rSges(4,3) * t900;
t840 = Icges(4,1) * t866 + Icges(4,4) * t865 + Icges(4,5) * t897;
t839 = Icges(4,1) * t864 + Icges(4,4) * t863 - Icges(4,5) * t900;
t838 = Icges(4,4) * t866 + Icges(4,2) * t865 + Icges(4,6) * t897;
t837 = Icges(4,4) * t864 + Icges(4,2) * t863 - Icges(4,6) * t900;
t836 = Icges(4,5) * t866 + Icges(4,6) * t865 + Icges(4,3) * t897;
t835 = Icges(4,5) * t864 + Icges(4,6) * t863 - Icges(4,3) * t900;
t834 = (t860 * t897 + t861 * t900) * qJD(2);
t829 = qJD(1) * t842 - t872 * t849 + t934;
t828 = -t954 + t873 * t849 + (-t841 + t946) * qJD(1);
t827 = t841 * t872 - t842 * t873 + t985;
t808 = rSges(7,1) * t822 + rSges(7,2) * t819;
t805 = Icges(7,5) * t822 + Icges(7,6) * t819;
t800 = t803 * t900;
t799 = t804 * t900;
t798 = t803 * t897;
t797 = t804 * t897;
t796 = rSges(7,3) * t897 + t900 * t955;
t795 = -rSges(7,3) * t900 + t897 * t955;
t790 = Icges(7,3) * t897 + t900 * t950;
t789 = -Icges(7,3) * t900 + t897 * t950;
t778 = rSges(8,1) * t803 + rSges(8,2) * t804;
t777 = Icges(8,1) * t803 + Icges(8,4) * t804;
t776 = Icges(8,4) * t803 + Icges(8,2) * t804;
t775 = Icges(8,5) * t803 + Icges(8,6) * t804;
t774 = (t803 * t894 - t804 * t892) * pkin(3);
t773 = rSges(8,1) * t799 - rSges(8,2) * t800 + rSges(8,3) * t897;
t772 = rSges(8,1) * t797 - rSges(8,2) * t798 - rSges(8,3) * t900;
t771 = Icges(8,1) * t799 - Icges(8,4) * t800 + Icges(8,5) * t897;
t770 = Icges(8,1) * t797 - Icges(8,4) * t798 - Icges(8,5) * t900;
t769 = Icges(8,4) * t799 - Icges(8,2) * t800 + Icges(8,6) * t897;
t768 = Icges(8,4) * t797 - Icges(8,2) * t798 - Icges(8,6) * t900;
t767 = Icges(8,5) * t799 - Icges(8,6) * t800 + Icges(8,3) * t897;
t766 = Icges(8,5) * t797 - Icges(8,6) * t798 - Icges(8,3) * t900;
t765 = (t799 * t894 + t800 * t892) * pkin(3);
t764 = (t797 * t894 + t798 * t892) * pkin(3);
t754 = -pkin(4) * t989 - t781 * t786;
t753 = 0.1e1 / t754 ^ 2;
t745 = t978 * t900;
t743 = -t808 * t995 + (-pkin(7) * t900 + t796) * qJD(1);
t742 = -t808 * t994 + (pkin(7) * t897 - t795) * qJD(1);
t740 = qJD(1) * t773 - t744 * t778 + t934;
t739 = -t954 + t745 * t778 + (-t772 + t946) * qJD(1);
t737 = (t754 * t1013 - t892 * t756 / 0.2e1) * t990;
t736 = (t892 * t754 / 0.2e1 + t756 * t1013) * t990;
t732 = (t795 * t897 + t796 * t900) * t748;
t731 = t744 * t772 - t745 * t773 + t985;
t729 = t734 * t870 + t735 * t935;
t728 = -t734 * t935 + t735 * t870;
t727 = -qJD(4) * t729 + qJD(1);
t726 = t734 * t865 - t735 * t866;
t725 = t734 * t866 + t735 * t865;
t724 = t734 * t863 - t735 * t864;
t723 = t734 * t864 + t735 * t863;
t722 = t725 * t898 + t895 * t897;
t721 = -t725 * t895 + t897 * t898;
t720 = t723 * t898 - t895 * t900;
t719 = -t723 * t895 - t898 * t900;
t718 = -t736 * t803 + t737 * t804;
t717 = t736 * t804 + t737 * t803;
t716 = -t736 * t799 - t737 * t800;
t715 = -t736 * t800 + t737 * t799;
t714 = -t736 * t797 - t737 * t798;
t713 = -t736 * t798 + t737 * t797;
t710 = pkin(9) * t728 - pkin(11) * t729;
t709 = rSges(5,1) * t728 + rSges(5,2) * t729;
t708 = Icges(5,1) * t728 + Icges(5,4) * t729;
t707 = Icges(5,4) * t728 + Icges(5,2) * t729;
t706 = Icges(5,5) * t728 + Icges(5,6) * t729;
t705 = pkin(9) * t725 - pkin(11) * t726;
t704 = pkin(9) * t723 - pkin(11) * t724;
t703 = rSges(5,1) * t725 + rSges(5,2) * t726 + rSges(5,3) * t897;
t702 = rSges(5,1) * t723 + rSges(5,2) * t724 - rSges(5,3) * t900;
t701 = Icges(5,1) * t725 + Icges(5,4) * t726 + Icges(5,5) * t897;
t700 = Icges(5,1) * t723 + Icges(5,4) * t724 - Icges(5,5) * t900;
t699 = Icges(5,4) * t725 + Icges(5,2) * t726 + Icges(5,6) * t897;
t698 = Icges(5,4) * t723 + Icges(5,2) * t724 - Icges(5,6) * t900;
t697 = Icges(5,5) * t725 + Icges(5,6) * t726 + Icges(5,3) * t897;
t696 = Icges(5,5) * t723 + Icges(5,6) * t724 - Icges(5,3) * t900;
t695 = rSges(9,1) * t717 + rSges(9,2) * t718;
t694 = Icges(9,1) * t717 + Icges(9,4) * t718;
t693 = Icges(9,4) * t717 + Icges(9,2) * t718;
t692 = Icges(9,5) * t717 + Icges(9,6) * t718;
t691 = rSges(9,1) * t715 + rSges(9,2) * t716 + rSges(9,3) * t897;
t690 = rSges(9,1) * t713 + rSges(9,2) * t714 - rSges(9,3) * t900;
t689 = Icges(9,1) * t715 + Icges(9,4) * t716 + Icges(9,5) * t897;
t688 = Icges(9,1) * t713 + Icges(9,4) * t714 - Icges(9,5) * t900;
t687 = Icges(9,4) * t715 + Icges(9,2) * t716 + Icges(9,6) * t897;
t686 = Icges(9,4) * t713 + Icges(9,2) * t714 - Icges(9,6) * t900;
t685 = Icges(9,5) * t715 + Icges(9,6) * t716 + Icges(9,3) * t897;
t684 = Icges(9,5) * t713 + Icges(9,6) * t714 - Icges(9,3) * t900;
t683 = -rSges(6,3) * t729 + (rSges(6,1) * t898 - rSges(6,2) * t895) * t728;
t682 = -Icges(6,5) * t729 + (Icges(6,1) * t898 - Icges(6,4) * t895) * t728;
t681 = -Icges(6,6) * t729 + (Icges(6,4) * t898 - Icges(6,2) * t895) * t728;
t680 = -Icges(6,3) * t729 + (Icges(6,5) * t898 - Icges(6,6) * t895) * t728;
t679 = rSges(6,1) * t722 + rSges(6,2) * t721 - rSges(6,3) * t726;
t678 = rSges(6,1) * t720 + rSges(6,2) * t719 - rSges(6,3) * t724;
t677 = Icges(6,1) * t722 + Icges(6,4) * t721 - Icges(6,5) * t726;
t676 = Icges(6,1) * t720 + Icges(6,4) * t719 - Icges(6,5) * t724;
t675 = Icges(6,4) * t722 + Icges(6,2) * t721 - Icges(6,6) * t726;
t674 = Icges(6,4) * t720 + Icges(6,2) * t719 - Icges(6,6) * t724;
t673 = Icges(6,5) * t722 + Icges(6,6) * t721 - Icges(6,3) * t726;
t672 = Icges(6,5) * t720 + Icges(6,6) * t719 - Icges(6,3) * t724;
t671 = 0.2e1 * (((t786 * t971 + (t758 * t781 - t993) * pkin(4)) * t1015 + (-t783 * t801 * t906 + t1004 * t784) * t1005) / t754 - ((-t759 * t781 + t943) * t1015 + (t754 * t784 + t783 * t786) * t1005) * t753 * t1004) * pkin(8) / (t753 * t756 ^ 2 + 0.1e1) * t785 * t905;
t670 = (-t671 + t978) * t900;
t669 = t671 * t897 + t744;
t668 = -t669 * t695 - t744 * t774 + (t691 + t765) * qJD(1) + t934;
t667 = -t954 + t670 * t695 + t745 * t774 + (-t690 - t764 + t946) * qJD(1);
t662 = -qJD(4) * t724 + t664;
t661 = -qJD(4) * t726 + t663;
t660 = qJD(1) * t703 - t663 * t709 + t927;
t659 = t664 * t709 + (-t702 + t933) * qJD(1) + t944;
t658 = t669 * t690 - t670 * t691 + t744 * t764 - t745 * t765 + t985;
t657 = t663 * t702 - t664 * t703 + t945;
t656 = qJD(1) * t705 - t661 * t683 - t663 * t710 + t727 * t679 + t927;
t655 = t662 * t683 + t664 * t710 - t727 * t678 + (-t704 + t933) * qJD(1) + t944;
t654 = t661 * t678 - t662 * t679 + t663 * t704 - t664 * t705 + t945;
t1 = t664 * ((-t697 * t900 + t699 * t724 + t701 * t723) * t663 + (-t696 * t900 + t698 * t724 + t700 * t723) * t664 + (-t706 * t900 + t707 * t724 + t708 * t723) * qJD(1)) / 0.2e1 + t670 * ((-t685 * t900 + t687 * t714 + t689 * t713) * t669 + (-t684 * t900 + t686 * t714 + t688 * t713) * t670 + (-t692 * t900 + t693 * t714 + t694 * t713) * qJD(1)) / 0.2e1 + t745 * ((-t767 * t900 - t769 * t798 + t771 * t797) * t744 + (-t766 * t900 - t768 * t798 + t770 * t797) * t745 + (-t775 * t900 - t776 * t798 + t777 * t797) * qJD(1)) / 0.2e1 + t663 * ((t697 * t897 + t699 * t726 + t701 * t725) * t663 + (t696 * t897 + t698 * t726 + t700 * t725) * t664 + (t706 * t897 + t707 * t726 + t708 * t725) * qJD(1)) / 0.2e1 + t669 * ((t685 * t897 + t687 * t716 + t689 * t715) * t669 + (t684 * t897 + t686 * t716 + t688 * t715) * t670 + (t692 * t897 + t693 * t716 + t694 * t715) * qJD(1)) / 0.2e1 + t873 * ((-t836 * t900 + t838 * t863 + t840 * t864) * t872 + (-t835 * t900 + t837 * t863 + t839 * t864) * t873 + (-t845 * t900 + t846 * t863 + t847 * t864) * qJD(1)) / 0.2e1 + t872 * ((t836 * t897 + t838 * t865 + t840 * t866) * t872 + (t835 * t897 + t837 * t865 + t839 * t866) * t873 + (t845 * t897 + t846 * t865 + t847 * t866) * qJD(1)) / 0.2e1 + t744 * ((t767 * t897 - t769 * t800 + t771 * t799) * t744 + (t766 * t897 - t768 * t800 + t770 * t799) * t745 + (t775 * t897 - t776 * t800 + t777 * t799) * qJD(1)) / 0.2e1 + m(3) * (t834 ^ 2 + t843 ^ 2 + t844 ^ 2) / 0.2e1 + m(4) * (t827 ^ 2 + t828 ^ 2 + t829 ^ 2) / 0.2e1 + t727 * ((-t673 * t661 - t672 * t662 - t680 * t727) * t729 + ((-t675 * t895 + t677 * t898) * t661 + (-t674 * t895 + t676 * t898) * t662 + (-t681 * t895 + t682 * t898) * t727) * t728) / 0.2e1 - ((-t805 * t900 + t897 * t947) * qJD(1) + (t789 * t1022 + (t948 * t897 + (-t790 + t949) * t900) * t897) * t748) * t994 / 0.2e1 + ((t897 * t805 + t900 * t947) * qJD(1) + (t790 * t1023 + (t949 * t900 + (-t789 + t948) * t897) * t900) * t748) * t995 / 0.2e1 + m(6) * (t654 ^ 2 + t655 ^ 2 + t656 ^ 2) / 0.2e1 + ((t897 * t854 + t900 * t937) * t887 - (t897 * t853 + t900 * t938) * t981 + (t897 * t874 + t900 * t936) * qJD(1)) * t887 / 0.2e1 - ((-t900 * t854 + t897 * t937) * t887 - (-t900 * t853 + t897 * t938) * t981 + (-t900 * t874 + t897 * t936) * qJD(1)) * t981 / 0.2e1 + m(5) * (t657 ^ 2 + t659 ^ 2 + t660 ^ 2) / 0.2e1 + m(7) * (t732 ^ 2 + t742 ^ 2 + t743 ^ 2) / 0.2e1 + t662 * ((-t673 * t724 + t675 * t719 + t677 * t720) * t661 + (-t672 * t724 + t674 * t719 + t676 * t720) * t662 + (-t680 * t724 + t681 * t719 + t682 * t720) * t727) / 0.2e1 + t661 * ((-t673 * t726 + t675 * t721 + t677 * t722) * t661 + (-t672 * t726 + t674 * t721 + t676 * t722) * t662 + (-t680 * t726 + t681 * t721 + t682 * t722) * t727) / 0.2e1 + m(8) * (t731 ^ 2 + t739 ^ 2 + t740 ^ 2) / 0.2e1 + m(9) * (t658 ^ 2 + t667 ^ 2 + t668 ^ 2) / 0.2e1 + (m(2) * (t878 ^ 2 + t879 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t1006 * t858 + t1009 * t856) * t897 - (t1006 * t857 + t1009 * t855) * t900) * qJD(2) + ((t792 * t819 + t794 * t822) * t897 - (t791 * t819 + t793 * t822) * t900) * t748 + (t838 * t870 - t840 * t935) * t872 + (t837 * t870 - t839 * t935) * t873 + (t769 * t804 + t771 * t803) * t744 + (t768 * t804 + t770 * t803) * t745 + (t699 * t729 + t701 * t728) * t663 + (t698 * t729 + t700 * t728) * t664 + (t687 * t718 + t689 * t717) * t669 + (t686 * t718 + t688 * t717) * t670 + (t1006 * t876 + t1009 * t875 + t693 * t718 + t694 * t717 + t707 * t729 + t708 * t728 + t776 * t804 + t777 * t803 + t806 * t819 + t807 * t822 + t846 * t870 - t847 * t935) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
