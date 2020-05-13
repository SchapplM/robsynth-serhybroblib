% Calculate kinetic energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE2_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:10
% EndTime: 2020-05-02 20:58:19
% DurationCPUTime: 8.46s
% Computational Cost: add. (8565->448), mult. (13515->664), div. (0->0), fcn. (20472->74), ass. (0->278)
t862 = sin(qJ(1));
t978 = t862 ^ 2;
t868 = cos(qJ(1));
t977 = t868 ^ 2;
t1001 = Icges(4,3) + Icges(9,3);
t855 = sin(pkin(19));
t858 = cos(pkin(19));
t860 = sin(qJ(3));
t866 = cos(qJ(3));
t763 = t855 * t866 + t858 * t860;
t765 = t855 * t860 - t858 * t866;
t736 = qJ(2) + atan2(t765, t763);
t734 = sin(t736);
t735 = cos(t736);
t848 = qJ(2) + qJ(3);
t829 = sin(t848);
t832 = cos(t848);
t1000 = Icges(4,5) * t832 - Icges(9,5) * t734 - Icges(4,6) * t829 - Icges(9,6) * t735;
t999 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t715 = atan2(t765, -t763) + t736;
t713 = sin(t715);
t714 = cos(t715);
t863 = sin(pkin(18));
t864 = sin(pkin(17));
t869 = cos(pkin(18));
t870 = cos(pkin(17));
t770 = t863 * t870 - t864 * t869;
t773 = t863 * t864 + t869 * t870;
t861 = sin(qJ(2));
t867 = cos(qJ(2));
t732 = t770 * t861 + t773 * t867;
t931 = t770 * t867 - t773 * t861;
t998 = -Icges(3,5) * t861 + Icges(7,5) * t931 + Icges(10,5) * t713 - Icges(3,6) * t867 - Icges(7,6) * t732 + Icges(10,6) * t714;
t771 = -t860 * t869 + t863 * t866;
t772 = t860 * t863 + t866 * t869;
t730 = -t771 * t867 + t772 * t861;
t731 = t771 * t861 + t772 * t867;
t853 = sin(pkin(22));
t856 = cos(pkin(22));
t766 = t853 * t863 + t856 * t869;
t767 = t853 * t869 - t856 * t863;
t841 = pkin(22) + pkin(21);
t822 = sin(t841);
t823 = cos(t841);
t668 = -qJ(2) - atan2(t766 * t867 - t767 * t861, t766 * t861 + t767 * t867) + pkin(21) - atan2(t730 * t823 + t731 * t822, -t730 * t822 + t731 * t823);
t666 = sin(t668);
t667 = cos(t668);
t997 = Icges(11,5) * t667 + Icges(11,6) * t666;
t818 = pkin(5) * t860 + pkin(1);
t948 = t866 * t867;
t994 = pkin(5) * t948 - t818 * t861;
t993 = qJD(1) / 0.2e1;
t992 = -qJD(2) / 0.2e1;
t991 = Icges(7,4) * t931;
t842 = qJD(2) + qJD(3);
t787 = t842 * t862;
t788 = t842 * t868;
t932 = pkin(18) - t841;
t816 = -pkin(20) + t932;
t885 = qJ(4) + t816;
t805 = qJ(1) + t885;
t783 = cos(t805);
t886 = -qJ(4) + t816;
t808 = -qJ(1) + t886;
t786 = cos(t808);
t849 = qJ(1) + qJ(4);
t833 = cos(t849);
t843 = qJD(1) + qJD(4);
t990 = t843 * (-t833 / 0.2e1 - t783 / 0.4e1 - t786 / 0.4e1);
t806 = qJ(1) + t886;
t780 = sin(t806);
t807 = -qJ(1) + t885;
t781 = sin(t807);
t850 = qJ(1) - qJ(4);
t831 = sin(t850);
t844 = qJD(1) - qJD(4);
t989 = t844 * (t831 / 0.2e1 + t781 / 0.4e1 - t780 / 0.4e1);
t988 = t862 * t868;
t958 = Icges(9,4) * t734;
t914 = -Icges(9,2) * t735 - t958;
t695 = -Icges(9,6) * t868 + t862 * t914;
t696 = Icges(9,6) * t862 + t868 * t914;
t957 = Icges(9,4) * t735;
t919 = -Icges(9,1) * t734 - t957;
t697 = -Icges(9,5) * t868 + t862 * t919;
t698 = Icges(9,5) * t862 + t868 * t919;
t709 = -Icges(9,2) * t734 + t957;
t710 = Icges(9,1) * t735 - t958;
t960 = Icges(4,4) * t832;
t916 = -Icges(4,2) * t829 + t960;
t741 = -Icges(4,6) * t868 + t862 * t916;
t742 = Icges(4,6) * t862 + t868 * t916;
t961 = Icges(4,4) * t829;
t921 = Icges(4,1) * t832 - t961;
t743 = -Icges(4,5) * t868 + t862 * t921;
t744 = Icges(4,5) * t862 + t868 * t921;
t760 = Icges(4,2) * t832 + t961;
t761 = Icges(4,1) * t829 + t960;
t987 = (t695 * t735 + t697 * t734 + t741 * t829 - t743 * t832) * t788 + (-t696 * t735 - t698 * t734 - t742 * t829 + t744 * t832) * t787 + (-t709 * t735 - t710 * t734 - t760 * t829 + t761 * t832) * qJD(1);
t986 = (-t1000 * t862 + t1001 * t868) * t788 + (t1000 * t868 + t1001 * t862) * t787 + (Icges(4,5) * t829 + Icges(9,5) * t735 + Icges(4,6) * t832 - Icges(9,6) * t734) * qJD(1);
t985 = rSges(5,1) * cos(t816) - rSges(5,2) * sin(t816) - pkin(15) - t994;
t941 = t862 * qJD(1);
t940 = t868 * qJD(1);
t984 = t862 * t999 + t868 * t998;
t983 = Icges(3,5) * t867 + Icges(7,5) * t732 - Icges(10,5) * t714 - Icges(3,6) * t861 + Icges(7,6) * t931 + Icges(10,6) * t713;
t982 = t862 * t998 - t868 * t999;
t955 = Icges(10,4) * t714;
t680 = Icges(10,2) * t713 - t955;
t956 = Icges(10,4) * t713;
t681 = -Icges(10,1) * t714 + t956;
t959 = Icges(7,4) * t732;
t702 = Icges(7,2) * t931 + t959;
t703 = Icges(7,1) * t732 + t991;
t962 = Icges(3,4) * t867;
t791 = -Icges(3,2) * t861 + t962;
t963 = Icges(3,4) * t861;
t792 = Icges(3,1) * t867 - t963;
t981 = t680 * t714 + t681 * t713 - t702 * t732 + t703 * t931 - t791 * t867 - t792 * t861;
t913 = Icges(10,2) * t714 + t956;
t676 = Icges(10,6) * t862 + t868 * t913;
t918 = Icges(10,1) * t713 + t955;
t678 = Icges(10,5) * t862 + t868 * t918;
t915 = -Icges(7,2) * t732 + t991;
t689 = Icges(7,6) * t862 + t868 * t915;
t920 = Icges(7,1) * t931 - t959;
t691 = Icges(7,5) * t862 + t868 * t920;
t917 = -Icges(3,2) * t867 - t963;
t751 = Icges(3,6) * t862 + t868 * t917;
t922 = -Icges(3,1) * t861 - t962;
t753 = Icges(3,5) * t862 + t868 * t922;
t980 = t676 * t714 + t678 * t713 - t689 * t732 + t691 * t931 - t751 * t867 - t753 * t861;
t675 = -Icges(10,6) * t868 + t862 * t913;
t677 = -Icges(10,5) * t868 + t862 * t918;
t690 = -Icges(7,6) * t868 + t862 * t915;
t692 = -Icges(7,5) * t868 + t862 * t920;
t750 = -Icges(3,6) * t868 + t862 * t917;
t752 = -Icges(3,5) * t868 + t862 * t922;
t979 = -t675 * t714 - t677 * t713 + t690 * t732 - t692 * t931 + t750 * t867 + t752 * t861;
t967 = pkin(5) * qJD(3);
t726 = t994 * qJD(2) + (-t860 * t861 + t948) * t967;
t813 = -qJ(1) + t816;
t802 = sin(t813);
t974 = -t802 / 0.2e1;
t812 = qJ(1) + t816;
t803 = cos(t812);
t973 = -t803 / 0.2e1;
t972 = -qJD(1) / 0.2e1 + t992;
t971 = t993 + t992;
t968 = pkin(1) * qJD(2);
t966 = t862 * rSges(6,1);
t965 = t862 * pkin(15);
t964 = t868 * rSges(6,1);
t859 = sin(qJ(4));
t954 = t859 * t862;
t953 = t859 * t868;
t952 = t861 * t866;
t854 = sin(pkin(20));
t857 = cos(pkin(20));
t764 = -t854 * t869 + t857 * t863;
t768 = t854 * t863 + t857 * t869;
t727 = t764 * t866 - t768 * t860;
t728 = t764 * t860 + t768 * t866;
t895 = t861 * t727 + t728 * t867;
t896 = t727 * t867 - t861 * t728;
t659 = atan2(t822 * t895 - t823 * t896, -t822 * t896 - t823 * t895) + t848;
t657 = sin(t659);
t951 = t862 * t657;
t865 = cos(qJ(4));
t950 = t862 * t865;
t949 = t865 * t868;
t947 = t868 * t657;
t944 = qJD(2) * t862;
t943 = qJD(2) * t868;
t942 = qJD(4) * t657;
t847 = pkin(18) - pkin(22);
t828 = -qJ(2) + t847;
t939 = pkin(21) - atan2(cos(t828), -sin(t828));
t938 = t861 * t968;
t937 = t867 * t968;
t936 = t964 / 0.2e1;
t814 = t932 - t848;
t725 = -atan2(-sin(t814), cos(t814)) + t939;
t930 = t862 * t937;
t929 = t868 * t937;
t794 = rSges(3,1) * t861 + rSges(3,2) * t867;
t926 = rSges(4,1) * t832 - rSges(4,2) * t829;
t824 = sin(t847);
t825 = cos(t847);
t925 = rSges(8,1) * t825 - rSges(8,2) * t824;
t924 = -rSges(9,1) * t734 - rSges(9,2) * t735;
t724 = -qJ(2) + t725;
t923 = -rSges(11,1) * sin(t724) + rSges(11,2) * cos(t724);
t658 = cos(t659);
t652 = -t658 * t954 - t949;
t653 = t658 * t950 - t953;
t654 = -t658 * t953 + t950;
t655 = t658 * t949 + t954;
t907 = (Icges(6,5) * t653 + Icges(6,6) * t652 + Icges(6,3) * t951) * t862 + (Icges(6,5) * t655 + Icges(6,6) * t654 + Icges(6,3) * t947) * t868;
t897 = sin(t725) * (rSges(11,1) * t861 + rSges(11,2) * t867) + cos(t725) * (rSges(11,1) * t867 - rSges(11,2) * t861);
t799 = -rSges(7,1) * t863 + rSges(7,2) * t869;
t800 = rSges(7,1) * t869 + rSges(7,2) * t863;
t891 = t799 * t864 - t800 * t870;
t888 = -Icges(11,5) * t666 + Icges(11,6) * t667;
t817 = pkin(1) * t861 - pkin(15);
t887 = t817 * t941 - t929;
t884 = -qJD(2) * (pkin(5) * t952 + t818 * t867) - (t860 * t867 + t952) * t967;
t779 = sin(t805);
t782 = sin(t808);
t830 = sin(t849);
t883 = (t830 / 0.2e1 + t779 / 0.4e1 - t782 / 0.4e1) * t843;
t784 = cos(t806);
t785 = cos(t807);
t834 = cos(t850);
t882 = (t834 / 0.2e1 - t784 / 0.4e1 - t785 / 0.4e1) * t844;
t881 = t907 * t657;
t877 = -pkin(4) * sin(qJ(2) - t939) - t817 + t923;
t876 = rSges(9,3) * qJD(1) + (-rSges(9,1) * t735 + rSges(9,2) * t734) * t842;
t872 = qJD(2) ^ 2;
t871 = pkin(11) + rSges(6,3);
t852 = qJ(1) - qJ(2);
t851 = qJ(1) + qJ(2);
t839 = t868 * rSges(6,2);
t838 = t868 * pkin(15);
t837 = t862 * rSges(6,2);
t836 = qJ(1) - t848;
t835 = qJ(1) + t848;
t827 = qJD(1) - t842;
t826 = qJD(1) + t842;
t821 = pkin(15) * t940;
t819 = -t966 / 0.2e1;
t804 = cos(t813);
t801 = sin(t812);
t798 = rSges(2,1) * t868 - rSges(2,2) * t862;
t797 = rSges(3,1) * t867 - rSges(3,2) * t861;
t795 = rSges(2,1) * t862 + rSges(2,2) * t868;
t762 = rSges(4,1) * t829 + rSges(4,2) * t832;
t746 = rSges(4,3) * t862 + t868 * t926;
t745 = -rSges(4,3) * t868 + t862 * t926;
t738 = -t799 * t870 - t800 * t864;
t719 = -t797 * t943 + (t868 * rSges(3,3) + (-pkin(15) + t794) * t862) * qJD(1);
t718 = t821 - t797 * t944 + qJD(1) * (rSges(3,3) * t862 - t794 * t868);
t717 = -t930 + (t862 * rSges(8,3) + (-t817 - t925) * t868) * qJD(1);
t716 = qJD(1) * (rSges(8,3) * t868 + t862 * t925) + t887;
t712 = t738 * t861 - t867 * t891;
t711 = t738 * t867 + t861 * t891;
t707 = rSges(7,3) * t862 + t711 * t868;
t706 = -rSges(7,3) * t868 + t711 * t862;
t705 = -t930 - t762 * t787 + (-t817 * t868 + t746) * qJD(1);
t704 = -qJD(1) * t745 - t762 * t788 + t887;
t700 = t745 * t787 + t746 * t788 - t938;
t699 = t924 * t842;
t686 = t884 * t862 + (t862 * rSges(5,3) - t985 * t868) * qJD(1);
t685 = t884 * t868 + (t868 * rSges(5,3) + t985 * t862) * qJD(1);
t684 = t842 * t923 - t938;
t683 = (t868 * t924 + t838) * qJD(1) + t876 * t862;
t682 = (-pkin(15) - t924) * t941 + t876 * t868;
t672 = -t712 * t944 + (-pkin(14) * t868 + t707) * qJD(1);
t671 = -t712 * t943 + (pkin(14) * t862 - t706) * qJD(1);
t670 = (t706 * t862 + t707 * t868) * qJD(2);
t669 = -pkin(2) * t734 * t842 + (rSges(10,1) * t713 + rSges(10,2) * t714) * qJD(2);
t665 = (rSges(10,1) * t944 + rSges(10,2) * t940) * t714 + (rSges(10,1) * t940 - rSges(10,2) * t944) * t713 + (rSges(10,3) * t862 + t838) * qJD(1) + (-t734 * t940 - t735 * t787) * pkin(2);
t664 = (rSges(10,1) * t943 - rSges(10,2) * t941) * t714 + (-rSges(10,1) * t941 - rSges(10,2) * t943) * t713 + (rSges(10,3) * t868 - t965) * qJD(1) + (t734 * t941 - t735 * t788) * pkin(2);
t663 = -t930 + t897 * t787 + (t862 * rSges(11,3) + t868 * t877) * qJD(1);
t662 = -t929 + t897 * t788 + (t868 * rSges(11,3) - t862 * t877) * qJD(1);
t661 = t821 + (t827 * cos(t836) / 0.2e1 + t826 * cos(t835) / 0.2e1) * pkin(5) + (sin(t852) * t971 + sin(t851) * t972) * pkin(1) + ((t974 - t801 / 0.2e1) * t871 + (t973 - t804 / 0.2e1) * pkin(9)) * qJD(1) + (t883 + t989) * rSges(6,2) + (t882 + t990) * rSges(6,1);
t660 = (-t827 * sin(t836) / 0.2e1 - t826 * sin(t835) / 0.2e1) * pkin(5) + (cos(t852) * t971 + cos(t851) * t972) * pkin(1) + (-t965 + (t804 / 0.2e1 + t973) * t871 + (t801 / 0.2e1 + t974) * pkin(9)) * qJD(1) + (t882 - t990) * rSges(6,2) + (t883 - t989) * rSges(6,1);
t656 = -t658 * qJD(4) + qJD(1);
t644 = Icges(11,3) * t862 + t868 * t888;
t643 = -Icges(11,3) * t868 + t862 * t888;
t642 = -((t839 + t966) * t834 + (t837 - t964) * t831 + (t839 - t966) * t833 + t830 * (t837 + t964) + (-t781 + t780) * (t936 - t837 / 0.2e1) + (-t782 + t779) * (t936 + t837 / 0.2e1) + (t785 + t784) * (-t839 / 0.2e1 + t819) + (t786 + t783) * (t839 / 0.2e1 + t819) + ((-t803 + t804) * t868 + (-t801 - t802) * t862) * rSges(6,3)) * t942 / 0.2e1 + t726;
t641 = -Icges(6,5) * t658 + (Icges(6,1) * t865 - Icges(6,4) * t859) * t657;
t640 = -Icges(6,6) * t658 + (Icges(6,4) * t865 - Icges(6,2) * t859) * t657;
t639 = -Icges(6,3) * t658 + (Icges(6,5) * t865 - Icges(6,6) * t859) * t657;
t638 = Icges(6,1) * t655 + Icges(6,4) * t654 + Icges(6,5) * t947;
t637 = Icges(6,1) * t653 + Icges(6,4) * t652 + Icges(6,5) * t951;
t636 = Icges(6,4) * t655 + Icges(6,2) * t654 + Icges(6,6) * t947;
t635 = Icges(6,4) * t653 + Icges(6,2) * t652 + Icges(6,6) * t951;
t1 = m(6) * (t642 ^ 2 + t660 ^ 2 + t661 ^ 2) / 0.2e1 + m(10) * (t664 ^ 2 + t665 ^ 2 + t669 ^ 2) / 0.2e1 + m(7) * (t670 ^ 2 + t671 ^ 2 + t672 ^ 2) / 0.2e1 + m(11) * (t662 ^ 2 + t663 ^ 2 + t684 ^ 2) / 0.2e1 + m(9) * (t682 ^ 2 + t683 ^ 2 + t699 ^ 2) / 0.2e1 + m(4) * (t700 ^ 2 + t704 ^ 2 + t705 ^ 2) / 0.2e1 + m(5) * (t685 ^ 2 + t686 ^ 2 + t726 ^ 2) / 0.2e1 + m(8) * (pkin(1) ^ 2 * t861 ^ 2 * t872 + t716 ^ 2 + t717 ^ 2) / 0.2e1 + m(3) * (t794 ^ 2 * t872 + t718 ^ 2 + t719 ^ 2) / 0.2e1 + t656 * ((-t658 * t639 + (-t640 * t859 + t641 * t865) * t657) * t656 + (((-t636 * t859 + t638 * t865) * t868 + (-t635 * t859 + t637 * t865) * t862) * t657 - t907 * t658) * t942) / 0.2e1 + (t868 * ((t639 * t947 + t640 * t654 + t641 * t655) * t656 + ((t635 * t654 + t637 * t655) * t862 + (t654 * t636 + t655 * t638 + t881) * t868) * t942) + t862 * ((t639 * t951 + t640 * t652 + t641 * t653) * t656 + ((t636 * t652 + t638 * t653) * t868 + (t652 * t635 + t653 * t637 + t881) * t862) * t942)) * t942 / 0.2e1 + ((t984 * t978 + (t979 * t868 + (t980 - t982) * t862) * t868) * qJD(2) + (t983 * t862 + t981 * t868) * qJD(1)) * t944 / 0.2e1 - ((t982 * t977 + (t980 * t862 + (t979 - t984) * t868) * t862) * qJD(2) + (t981 * t862 - t983 * t868) * qJD(1)) * t943 / 0.2e1 + (m(2) * (t795 ^ 2 + t798 ^ 2) + t825 ^ 2 * Icges(8,2) + (Icges(8,1) * t824 + 0.2e1 * Icges(8,4) * t825) * t824 + Icges(5,2) * t658 ^ 2 + (Icges(5,1) * t657 + 0.2e1 * Icges(5,4) * t658) * t657 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t695 * t734 - t697 * t735 - t741 * t832 - t743 * t829) * t788 + (-t696 * t734 + t698 * t735 + t742 * t832 + t744 * t829) * t787 + ((-t675 * t713 + t677 * t714 - t690 * t931 - t692 * t732 + t750 * t861 - t752 * t867) * t868 + (t676 * t713 - t678 * t714 + t689 * t931 + t691 * t732 - t751 * t861 + t753 * t867) * t862) * qJD(2) + (t667 ^ 2 * Icges(11,1) + t713 * t680 - t714 * t681 + t702 * t931 + t732 * t703 - t734 * t709 + t735 * t710 + t832 * t760 + t829 * t761 - t861 * t791 + t867 * t792 + (0.2e1 * Icges(11,4) * t667 + Icges(11,2) * t666) * t666) * qJD(1) + (-t977 - t978) * t842 * t997) * t993 + (-t997 * t941 + (-t643 * t988 + t978 * t644) * t842 + t986 * t862 + t987 * t868) * t787 / 0.2e1 - (t997 * t940 + (t977 * t643 - t644 * t988) * t842 + t987 * t862 - t986 * t868) * t788 / 0.2e1;
T = t1;
