% Calculate kinetic energy for
% fivebar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1DE1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1DE1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:23:45
% EndTime: 2020-04-27 02:24:03
% DurationCPUTime: 15.56s
% Computational Cost: add. (124280->600), mult. (375760->953), div. (866->13), fcn. (47077->12), ass. (0->414)
t792 = (pkin(2) ^ 2);
t1022 = -6 * t792;
t748 = sin(qJ(2));
t751 = cos(qJ(1));
t933 = qJD(1) * t751;
t690 = pkin(2) * t933;
t867 = pkin(3) * t690;
t749 = sin(qJ(1));
t750 = cos(qJ(2));
t930 = qJD(2) * t750;
t879 = t749 * t930;
t997 = pkin(2) * pkin(3);
t1020 = -t748 * t867 - t879 * t997;
t788 = pkin(3) ^ 2;
t702 = 0.10e2 / 0.3e1 * t788;
t780 = (pkin(5) ^ 2);
t779 = t780 ^ 2;
t783 = pkin(4) ^ 2;
t782 = t783 ^ 2;
t1008 = -t779 / 0.6e1 + t782 / 0.6e1;
t786 = t788 ^ 2;
t794 = (pkin(1) ^ 2);
t793 = t794 ^ 2;
t941 = t786 + t793;
t775 = 2 * t794;
t944 = t775 - t780;
t956 = t794 * t780;
t632 = t944 * t788 - t1008 + t941 - t956;
t800 = t792 ^ 2;
t819 = t632 + t800;
t1019 = (t702 + t944) * t1022 - 0.6e1 * t819;
t1021 = -2 * t794;
t940 = t788 - t794;
t687 = t940 * t792;
t704 = pkin(3) * t750;
t684 = pkin(1) + t704;
t934 = qJD(1) * t749;
t883 = t684 * t934;
t932 = qJD(2) * t748;
t904 = pkin(3) * t932;
t988 = pkin(2) * t751;
t1010 = pkin(2) * t883 + t904 * t988;
t927 = 0.2e1 * pkin(1);
t924 = 0.4e1 * pkin(3);
t729 = -t780 / 0.3e1;
t1012 = t729 - t783 / 0.3e1;
t891 = t794 + t1012;
t855 = t788 + t891;
t737 = t783 / 0.3e1;
t893 = t780 / 0.3e1 + t737 + t775;
t1018 = 0.10e2 / 0.3e1 * t800 - 0.2e1 * (-t788 + t893) * t792 + t855 * t1021;
t1017 = 0.2e1 * t690;
t705 = pkin(2) * t749;
t1016 = 0.2e1 * t705;
t1015 = 0.2e1 * t788;
t928 = t748 * qJD(1);
t878 = t751 * t928;
t820 = t878 + t879;
t1014 = -0.4e1 * t820;
t966 = t748 * t749;
t909 = pkin(3) * t966;
t870 = pkin(2) * t909;
t674 = -0.2e1 * t870;
t922 = pkin(1) * t988;
t686 = -0.2e1 * t922;
t938 = t792 + t794;
t886 = t788 + t938;
t683 = pkin(1) - t988;
t978 = t683 * t750;
t910 = pkin(3) * t978;
t605 = t674 + t686 + t886 + 0.2e1 * t910;
t963 = t750 * t749;
t965 = t748 * t751;
t646 = -t963 + t965;
t905 = pkin(2) * t934;
t681 = pkin(1) * t905;
t983 = (t681 + (-t683 * t932 + (-t646 * qJD(1) - t879) * pkin(2)) * pkin(3)) / t605 ^ 2;
t1013 = -0.2e1 * t983;
t987 = pkin(2) * t788;
t715 = t750 ^ 2;
t985 = pkin(3) * t715;
t698 = -t788 / 0.3e1 + t794;
t1011 = t1020 * t698;
t1009 = qJD(1) - qJD(2);
t1007 = 0.8e1 * pkin(1);
t1006 = -0.2e1 * t715;
t1005 = 0.4e1 * t715;
t1004 = -0.4e1 * t748;
t1003 = 0.4e1 * t786;
t1002 = -0.8e1 * t788;
t764 = 0.6e1 * t788;
t1001 = 2 * t792;
t767 = 5 * t800;
t803 = pkin(3) * t788;
t1000 = 0.4e1 * t803;
t699 = t794 - t792 / 0.3e1;
t999 = 0.24e2 * t699;
t998 = pkin(1) * pkin(2);
t996 = -pkin(4) - pkin(5);
t995 = pkin(5) - pkin(4);
t743 = t788 / 0.3e1;
t622 = -0.4e1 / 0.9e1 * t870 + t794 + t792 / 0.3e1 + t743 + t783 / 0.9e1 - t780 / 0.9e1;
t727 = -t780 / 0.6e1;
t746 = t792 / 0.2e1;
t947 = t746 + t794;
t823 = -t870 + t947;
t630 = t783 / 0.6e1 + t727 + t823;
t658 = t729 + t737 + t886;
t926 = 0.4e1 * pkin(1);
t957 = t788 * t715;
t961 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t714 = t750 * t715;
t971 = t714 * t803;
t573 = t698 * t674 + 0.6e1 * t622 * t957 + t658 * t961 + (t630 * t704 + t971) * t926;
t859 = -0.4e1 * t870;
t731 = -0.2e1 / 0.3e1 * t780;
t736 = 0.2e1 / 0.3e1 * t783;
t889 = t731 + t736 + t775;
t708 = t788 + t794;
t888 = t731 + t708;
t952 = t708 * (t736 + t888) + t800;
t586 = t658 * t859 + (t764 + t889) * t792 + t952;
t628 = t674 + t658;
t939 = -t792 + t794;
t688 = t939 * t788;
t921 = pkin(1) * t704;
t574 = t688 * t1005 + 0.4e1 * t628 * t921 + t586;
t643 = t699 * t674;
t960 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t604 = t658 * t960 + t643;
t606 = (t702 + t889) * t792 + t952;
t685 = 0.2e1 * t921;
t710 = -0.3e1 * t788 + t794;
t915 = 0.4e1 * t957;
t641 = t685 + t915 + t710;
t769 = -3 * t792;
t711 = t769 + t794;
t920 = 0.8e1 * t971;
t871 = pkin(1) * t920;
t662 = t711 * t871;
t935 = t794 - t780;
t885 = t788 + t935;
t679 = t783 + t885;
t701 = t708 ^ 2;
t718 = t751 ^ 2;
t756 = 0.15e2 * t786;
t757 = 0.15e2 * t788;
t771 = 0.3e1 * t793;
t799 = pkin(2) * t792;
t789 = t799 ^ 2;
t717 = t751 * t718;
t968 = t717 * t799;
t899 = t684 * t968;
t866 = -0.8e1 * t899;
t774 = 3 * t794;
t942 = -t780 + t783;
t887 = t774 + t942;
t896 = 0.6e1 * t921;
t906 = 0.12e2 * t957;
t977 = t684 * t751;
t559 = t641 * t866 + t662 + t604 * t906 + t586 * t896 + t789 + (t757 + t887) * t800 + t701 * t679 + (0.12e2 * t573 * t718 + t887 * t764 + t942 * t775 + t756 + t771) * t792 + 0.6e1 * (-t574 * t977 - t606 * t909) * pkin(2);
t765 = 0.3e1 * t788;
t694 = t765 + t938;
t663 = t694 * t705;
t986 = pkin(3) * t748;
t667 = t705 - t986;
t912 = t749 * t987;
t768 = 3 * t792;
t692 = t768 + t708;
t975 = t692 * t748;
t989 = pkin(1) * t750;
t587 = t912 * t1006 + t663 + (0.2e1 * t667 * t989 - t975) * pkin(3);
t707 = t1015 + t792;
t660 = -t748 * t803 + t912;
t979 = t660 * t715;
t588 = t707 * t986 + 0.2e1 * t979 + (-t940 + t685) * t705;
t621 = -pkin(3) * t975 + t663;
t689 = pkin(2) * t1003 + 0.8e1 * t788 * t799;
t897 = t748 * t960;
t624 = t897 * t1000 + t689 * t749;
t762 = 0.5e1 * t786;
t936 = t793 + t800;
t758 = 0.10e2 * t788;
t945 = t758 + t775;
t955 = t794 * t788;
t640 = t945 * t792 + t762 + t936 + 0.6e1 * t955;
t648 = t767 + (t758 + (6 * t794)) * t792 + t701;
t865 = 0.8e1 * t899;
t967 = t718 * t792;
t919 = -0.4e1 * t967;
t565 = t588 * t919 + t624 * t715 + (-0.4e1 * t621 * t989 + (t648 + t865) * t748) * pkin(3) + (0.4e1 * t587 * t977 + (-t640 + t871) * t749) * pkin(2);
t659 = t683 + t704;
t884 = t792 + t935;
t854 = t788 + t884;
t672 = -t783 + t854;
t639 = t686 + t672;
t619 = t674 + t639;
t918 = 0.2e1 * t967;
t642 = t686 + t918 + t939;
t841 = t642 * t1006 - t935;
t913 = pkin(2) * t966;
t795 = sqrt(0.4e1 * t687 * t718 + 0.4e1 * t672 * t922 - t786 - (t794 + (pkin(2) - t995) * (pkin(2) + t995)) * (t794 + (pkin(2) - t996) * (pkin(2) + t996)) + (t769 + t783 + t841) * t1015 + (-t619 * t978 + t639 * t913) * t924);
t553 = t559 * t659 + t565 * t795;
t994 = 0.1e1 / t553 / 0.4e1;
t602 = 0.1e1 / t605;
t993 = t602 / 0.2e1;
t744 = t788 / 0.2e1;
t992 = 0.4e1 / 0.3e1 * t792;
t991 = -t795 / 0.4e1;
t990 = t795 / 0.4e1;
t817 = t750 * t820;
t816 = pkin(3) * t817;
t815 = t683 * t816;
t880 = t748 * t930;
t851 = t642 * t880;
t678 = 0.2e1 * t681;
t852 = t749 * t792 * t933;
t840 = -0.4e1 * t852;
t634 = t678 + t840;
t981 = t634 * t715;
t984 = ((0.8e1 * t851 - 0.4e1 * t981) * t788 + (-0.4e1 * t672 * t998 - 0.8e1 * t687 * t751) * t934 + (t792 * t749 ^ 2 * t928 * t1007 + (0.4e1 * t619 * t932 - 0.8e1 * t681 * t750) * t683 + 0.4e1 * (-t619 * t750 * t934 + t639 * t820 + 0.2e1 * t815) * pkin(2)) * pkin(3)) / t795;
t818 = t820 * pkin(3);
t623 = pkin(2) * t818;
t982 = t623 * t750;
t929 = qJD(2) * t803;
t980 = (t788 * t690 - t750 * t929) * t715;
t706 = -t780 - t783;
t691 = t774 + t706;
t976 = t691 * t788;
t974 = t701 * (-t783 + t885);
t973 = t706 * t794;
t972 = t714 * t786;
t970 = t715 * t803;
t716 = t718 ^ 2;
t969 = t716 * t800;
t964 = t748 * t795;
t962 = t750 * t795;
t959 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t713 = t715 ^ 2;
t958 = t786 * t713;
t813 = t623 * t957;
t812 = t699 * t813;
t846 = t715 * t748 * t929;
t833 = -0.24e2 * t846;
t825 = pkin(1) * t833;
t954 = t711 * t825 - 0.24e2 * t812;
t732 = -0.3e1 / 0.2e1 * t780;
t951 = t779 / 0.2e1 - t782 / 0.2e1;
t953 = t708 * ((t732 + t775) * t788 - 0.3e1 / 0.2e1 * t956 + t941 + t951) + t789;
t950 = t727 - t783 / 0.6e1;
t741 = -0.2e1 / 0.3e1 * t783;
t949 = t731 + t741;
t948 = t732 + t774;
t946 = 0.4e1 / 0.7e1 * t794 - t780 / 0.7e1;
t943 = t779 - t782;
t937 = t793 - t786;
t931 = qJD(2) * t749;
t925 = 0.2e1 * pkin(3);
t923 = 0.4e1 * t704;
t917 = 0.8e1 * t958;
t916 = -0.6e1 * t957;
t914 = pkin(2) * t977;
t911 = t803 * t705;
t908 = 0.16e2 * t971;
t907 = -0.12e2 * t967;
t903 = pkin(3) * t930;
t902 = pkin(3) * t711 / 0.2e1;
t901 = -t986 / 0.2e1;
t900 = t602 * t959;
t898 = t748 * t969;
t895 = 0.4e1 / 0.3e1 * t780 + 0.4e1 / 0.3e1 * t783 + t1021;
t773 = 4 * t794;
t894 = 0.2e1 / 0.3e1 * t780 + t736 + t773;
t892 = t794 + t950;
t730 = -t780 / 0.2e1;
t890 = t730 - t783 / 0.2e1 + t794;
t882 = t714 * t932;
t881 = qJD(2) * t968;
t649 = -t904 + t905;
t868 = pkin(1) * t904;
t677 = -0.2e1 * t868;
t712 = t748 ^ 2;
t814 = t788 * t817 * t998;
t831 = t718 * t799 * t883;
t824 = -0.24e2 * t831;
t829 = t694 * t690 - t692 * t903;
t847 = t788 * t880;
t834 = -0.24e2 * t847;
t835 = t690 * t971;
t860 = pkin(3) * t881;
t837 = t748 * t860;
t838 = pkin(1) * t846;
t839 = -0.6e1 * t847;
t845 = -0.6e1 * t868;
t848 = t660 * t880;
t849 = t930 * t960;
t863 = -0.2e1 * t880;
t874 = t984 / 0.2e1;
t877 = -((0.8e1 * t750 * t684 * t860 + t712 * t881 * t1002 + t824 * t986 + (t707 * t903 - 0.4e1 * t848 + 0.2e1 * t980 + (-t940 * t933 + (-t748 * t931 + t750 * t933) * pkin(1) * t925) * pkin(2)) * t919 + 0.8e1 * t588 * t852 + 0.4e1 * ((t933 * t1006 + 0.4e1 * t748 * t879) * t987 + ((t690 - t903) * t704 - t667 * t904) * t927 + t829) * t914 + t835 * t1007 + t825 * t705 + (t849 * t1000 + t689 * t933) * t715 + t624 * t863 - 0.4e1 * t829 * t921 + 0.4e1 * t621 * t868 - t640 * t690 + t648 * t903 - 0.4e1 * t1010 * t587) * t795 + t565 * t874 + t649 * t559 + t659 * ((t677 - 0.8e1 * t847) * t866 + 0.24e2 * (-0.6e1 * t838 + 0.3e1 * (-0.4e1 / 0.9e1 * t879 - 0.4e1 / 0.9e1 * t878) * t957 * t997 + t622 * t839 - t623 * t685 + t630 * t677 + t1011) * t967 - 0.24e2 * t573 * t852 - 0.6e1 * (-0.8e1 * t688 * t880 + (t658 * pkin(2) * t1014 + (-0.4e1 * t628 * t932 - 0.8e1 * t982) * pkin(1)) * pkin(3)) * t914 + t604 * t834 - 0.24e2 * t658 * t814 + t586 * t845 + t954 + (0.8e1 * t837 + 0.24e2 * t831) * t641) + 0.6e1 * t659 * (t1010 * t574 + t1020 * t606)) / t553 ^ 2 / 0.4e1;
t671 = t765 + t783 + t884;
t620 = t671 + t686 + t859;
t631 = pkin(2) * t963 + t683 * t748;
t673 = -pkin(3) + t913;
t561 = -t620 * t978 + t631 * t795 + (-0.2e1 * pkin(1) * t673 * t751 + t671 * t966) * pkin(2) + (-t768 - t783 - t788 + t841 + t918) * pkin(3);
t876 = t561 * t994;
t638 = t686 + t783 + t854;
t590 = t638 * t683 + 0.2e1 * t642 * t704;
t594 = t638 * t750 + (t1005 - 0.2e1) * t683 * pkin(3);
t629 = -t673 + t978;
t563 = t590 * t748 + t594 * t705 + t629 * t795;
t875 = t563 * t994;
t784 = 0.1e1 / pkin(4);
t873 = t784 * t993;
t872 = 0.20e2 / 0.3e1 * t788;
t869 = pkin(1) * t908;
t636 = 0.4e1 / 0.3e1 * t957 + t685 + t698;
t862 = -0.24e2 * t636 * t969;
t858 = -(3 * t956) + t771 + t951;
t857 = -(6 * t956) + 0.6e1 * t793 + t943;
t856 = t746 + t892;
t853 = t717 * t800 * t934;
t850 = t786 * t882;
t598 = t800 + (t945 + t949) * t792 + t762 + 0.2e1 * t976 + t794 * (t794 + t949);
t612 = t767 + (0.5e1 * t788 + t691) * t1001 + t708 * (t741 + t888);
t576 = t598 * t705 - t612 * t986;
t836 = t690 * t958;
t676 = t730 + t886;
t832 = t676 * t859;
t572 = t870 * t1019 + (t756 + (-9 * t780 + 18 * t794) * t788 + t858) * t792 + (t757 + t948) * t800 + t953;
t647 = t750 * t751 + t966;
t830 = t598 * t690 - t612 * t903;
t826 = -0.48e2 * t838;
t665 = t794 + t792 / 0.4e1 + t788 / 0.4e1 - t780 / 0.8e1;
t579 = -0.32e2 / 0.21e2 * t665 * t870 + t800 / 0.7e1 + (0.16e2 / 0.21e2 * t788 + t946) * t792 + t786 / 0.7e1 + t946 * t788 + t793 - 0.3e1 / 0.7e1 * t956 + t779 / 0.42e2 - t782 / 0.42e2;
t728 = -t780 / 0.4e1;
t666 = t728 + t743 + t947;
t742 = 0.4e1 / 0.3e1 * t788;
t583 = -0.8e1 / 0.3e1 * t666 * t870 + t800 / 0.3e1 + (t742 + t729) * t792 + t793 - t786 / 0.3e1 + (t992 + 0.2e1 / 0.3e1 * t788 + t731) * t794 + t779 / 0.18e2 - t782 / 0.18e2;
t635 = -0.2e1 / 0.3e1 * t870 + t794 + t744 + t728;
t695 = (t773 + t780) * t788;
t700 = t794 - 0.2e1 / 0.3e1 * t792;
t564 = t700 * t917 + 0.14e2 * t579 * t957 + t698 * t832 - t940 * t800 + (t695 - 0.10e2 / 0.3e1 * t786 + 0.2e1 * t793 - t956) * t792 + t632 * t961 + (0.6e1 * t583 * t704 + t635 * t908) * pkin(1);
t589 = t832 + (t764 + t944) * t792 + t819;
t610 = t676 * t960 + t643;
t566 = t589 * t896 + t610 * t906 + t572 + t662;
t759 = -2 * t780;
t772 = 8 * t794;
t627 = t911 * t1004 + t1003 + (4 * t792 + t759 + t772) * t788;
t633 = t728 - t788 + t823;
t575 = t674 * t961 + t627 * t715 + t676 * t710 + (t633 * t923 + t920) * pkin(1);
t577 = t699 * t832 - t789 + (-t702 - t935) * t800 + (t695 + t937 + t1008) * t792 + t794 * t632;
t760 = -5 * t780;
t761 = 0.7e1 * t786;
t763 = 0.7e1 * t788;
t582 = (t763 + t948) * t800 + (t761 + (t760 + 10 * t794) * t788 + t858) * t792 + t953;
t796 = pkin(1) * t794;
t682 = -0.12e2 * pkin(1) * t803 + t796 * t924;
t697 = -0.8e1 * t786 + 0.12e2 * t955;
t596 = t682 * t750 + t697 * t715 + t869 + t917 + t941 - 0.6e1 * t955;
t613 = t674 * t960 + t676 * t711;
t670 = ((t794 * t1022) + t936) * t786;
t703 = -30 * t780 + 60 * t794;
t548 = -0.32e2 * t575 * t899 + 0.16e2 * t670 * t713 + 0.24e2 * t577 * t957 + (t759 + t773 + 0.28e2 * t788) * t789 + t679 * t974 + (0.24e2 * t564 * t718 + t703 * t786 + t857 * t764 + t943 * t775 - 0.6e1 * t793 * t780 + 0.4e1 * t796 ^ 2 + 0.28e2 * t803 ^ 2) * t792 + 0.8e1 * (-t566 * t977 - t582 * t909) * pkin(2) + (0.8e1 * t572 * t704 + 0.32e2 * t613 * t971) * pkin(1) + (0.16e2 * t596 * t716 + t703 * t788 + 0.70e2 * t786 + t800 + t857) * t800;
t650 = t992 + t744 + t892;
t651 = t742 + t856;
t597 = -t650 * t986 + t651 * t705;
t608 = -t800 + (-t872 + t894) * t792 - 0.3e1 * t786 + t895 * t788 + t793;
t654 = t788 + t856;
t693 = t1001 - t940;
t611 = t654 * t705 + t693 * t901;
t614 = -t687 - 0.5e1 / 0.3e1 * t786 + t893 * t788 + t794 * t891;
t567 = t597 * t915 + t608 * t901 + (-0.8e1 / 0.3e1 * t958 + t614) * t705 + (t611 * t704 - t748 * t972) * t926;
t656 = t768 + t855;
t657 = t694 + t1012;
t601 = -t656 * t986 + t657 * t705;
t664 = t744 + t792 + t950;
t616 = pkin(3) * t897 + t664 * t1016;
t568 = -0.4e1 * t616 * t957 + (t601 * t923 - 0.8e1 * t714 * t911) * pkin(1) + t576;
t609 = -(3 * t800) + (-t872 + t895) * t792 + t894 * t788 + t937;
t578 = t986 * t1018 + t609 * t705;
t669 = 0.10e2 * t976;
t580 = t789 + (0.21e2 * t788 + t691) * t800 + (t669 + t771 + 0.35e2 * t786 + 0.2e1 * t973) * t792 + (t761 + (t760 + t772 - 0.5e1 * t783) * t788 + t794 * (-t783 + t935)) * t708;
t653 = 0.3e1 / 0.2e1 * t792 + t765 + t890;
t668 = t705 + 0.2e1 * t986;
t581 = t710 * t705 + 0.4e1 * t979 + (t653 * t748 + t668 * t989) * t925;
t584 = 0.7e1 * t789 + (t763 + t691) * t767 + (t669 + 0.21e2 * t786 + 0.9e1 * t793 + 0.6e1 * t973) * t792 + t974;
t655 = t768 + 0.3e1 / 0.2e1 * t788 + t890;
t615 = t655 * t705 + t748 * t902;
t554 = t615 * t869 + t581 * t865 + t567 * t907 + t578 * t916 + (-0.6e1 * t576 * t989 + (t584 + t862) * t748) * pkin(3) + (0.6e1 * t568 * t977 + (t958 * t999 - t580) * t749) * pkin(2);
t547 = t548 * t659 + t554 * t795;
t822 = t547 * t875 + t561 * t990;
t821 = t547 * t876 + t563 * t991;
t600 = t1009 * t647;
t599 = t1009 * t646;
t560 = 0.1e1 / t561 ^ 2;
t557 = atan2(t563 * t873, t561 * t873);
t556 = cos(t557);
t555 = sin(t557);
t550 = t631 * t874 + (-t678 * t750 + (t620 * t748 + t962) * qJD(2)) * t683 + (t840 + 0.4e1 * t851 - 0.2e1 * t981) * pkin(3) + ((t671 * t750 - t964) * t931 + (-t620 * t963 + t647 * t795 + t671 * t965) * qJD(1) + 0.4e1 * t815 + (t673 * t934 - t820 * t988) * t927) * pkin(2);
t549 = t629 * t874 + (t590 * t750 - t683 * t964) * qJD(2) + (t634 * t750 - t642 * t932) * t748 * t925 + ((-t962 + (-t638 - 0.8e1 * t910) * t748) * t931 + ((t594 - t964) * t751 + (t962 + (t683 * t927 + t638) * t748 + (-pkin(3) + t989 + 0.2e1 * t985) * t1016) * t749) * qJD(1)) * pkin(2);
t545 = t822 * t900;
t544 = t821 * t900;
t543 = qJD(2) + 0.2e1 * ((t549 * t993 - t563 * t983) / t561 - (t550 * t993 - t561 * t983) * t563 * t560) * pkin(4) / (t560 * t563 ^ 2 + 0.1e1) * t605 * t784;
t542 = (t862 * t903 - 0.24e2 * (-0.8e1 / 0.3e1 * t847 + t677) * pkin(3) * t898 + 0.96e2 * t636 * t853 * t986 + (t710 * t690 + 0.2e1 * t653 * t903 + (t750 * t1017 + (-0.2e1 * t668 * t748 + 0.4e1 * t985) * qJD(2)) * pkin(3) * pkin(1) - 0.8e1 * t848 + 0.4e1 * t980) * t865 + (-0.8e1 / 0.3e1 * t836 + (-t650 * t903 + t651 * t690) * t915 + t614 * t690 - t608 * t903 / 0.2e1 + (0.32e2 / 0.3e1 * t972 * t705 + t597 * t750 * t1002) * t932 + (0.4e1 * t654 * t750 * t867 + ((0.12e2 * t712 * t715 - 0.4e1 * t713) * t786 + (t611 * t1004 - 0.2e1 * t693 * t985) * pkin(3)) * qJD(2)) * pkin(1)) * t907 + 0.24e2 * t567 * t852 + 0.6e1 * ((-0.4e1 * (pkin(3) * t849 + t664 * t1017) * t715 + 0.8e1 * t616 * t880) * t788 + (-0.8e1 * t835 + (-t656 * t903 + t657 * t690) * t923 + (-0.4e1 * pkin(3) * t601 + 0.24e2 * t715 * t911) * t932) * pkin(1) + t830) * t914 + t836 * t999 - 0.96e2 * t699 * t850 * t705 + (t655 * t690 + t902 * t930) * t869 + t615 * t826 + (t903 * t1018 + t609 * t690) * t916 + 0.12e2 * t578 * t847 - 0.6e1 * t830 * t921 + 0.6e1 * t576 * t868 - t580 * t690 + t584 * t903 + (-0.8e1 * t837 + t824) * t581 - 0.6e1 * t1010 * t568) * t795 + t554 * t874 + t649 * t548 + 0.8e1 * t659 * (0.2e1 * (-0.48e2 * pkin(1) * t970 - 0.2e1 * t697 * t750 - t682 - 0.32e2 * t972) * qJD(2) * t898 - 0.8e1 * t596 * t853 - 0.4e1 * (t627 * t863 + (t833 + (-t633 * t932 - t982) * t924) * pkin(1) + (t970 * t1014 - 0.2e1 * t818 * t961) * pkin(2)) * t899 + 0.3e1 * (-0.32e2 * t700 * t850 + (-0.2e1 / 0.3e1 * t879 - 0.2e1 / 0.3e1 * t878) * t869 * t997 + t635 * t826 - 0.64e2 / 0.3e1 * t665 * t813 - 0.28e2 * t579 * t847 + 0.6e1 * (-0.8e1 / 0.3e1 * t879 - 0.8e1 / 0.3e1 * t878) * t666 * t989 * t987 + t583 * t845 + 0.4e1 * t1011 * t676) * t967 - 0.6e1 * t564 * t852 - (t610 * t834 + (-0.6e1 * pkin(1) * t589 * t932 + (-0.24e2 * pkin(1) * t676 * t816 + t820 * t1019) * pkin(2)) * pkin(3) + t954) * t914 - 0.8e1 * t670 * t882 - t623 * t871 * t960 - 0.12e2 * t613 * t838 - 0.12e2 * t676 * t812 + t577 * t839 + t814 * t1019 - t572 * t868 + t1020 * t582 + (0.4e1 * t837 + 0.12e2 * t831) * t575 + t1010 * t566);
t541 = t544 * t647 - t545 * t646;
t540 = t544 * t646 + t545 * t647;
t539 = 0.1e1 / t541 ^ 2;
t538 = atan2(t540, t541);
t536 = cos(t538);
t535 = sin(t538);
t534 = (t822 * t1013 + (t550 * t990 + t561 * t984 / 0.8e1 + t542 * t875 + (t549 * t994 + t563 * t877) * t547) * t602) * t959;
t533 = (t821 * t1013 + (t542 * t876 + t549 * t991 - t563 * t984 / 0.8e1 + (t550 * t994 + t561 * t877) * t547) * t602) * t959;
t532 = qJD(1) + ((t533 * t646 + t534 * t647 - t544 * t600 + t545 * t599) / t541 - (t533 * t647 - t534 * t646 + t544 * t599 + t545 * t600) * t540 * t539) / (t539 * t540 ^ 2 + 0.1e1);
t1 = (Ifges(2,3) / 0.2e1 + m(3) * (t535 ^ 2 + t536 ^ 2) * t746) * qJD(1) ^ 2 + (Ifges(4,3) / 0.2e1 + m(5) * (t555 ^ 2 + t556 ^ 2) * t744) * qJD(2) ^ 2 + (Ifges(5,3) * t543 / 0.2e1 + (mrSges(5,1) * t556 - mrSges(5,2) * t555) * qJD(2) * pkin(3)) * t543 + (Ifges(3,3) * t532 / 0.2e1 + (-mrSges(3,1) * t536 + mrSges(3,2) * t535) * qJD(1) * pkin(2)) * t532;
T = t1;
