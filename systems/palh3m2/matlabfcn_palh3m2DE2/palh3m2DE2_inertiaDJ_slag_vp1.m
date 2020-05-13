% Calculate time derivative of joint inertia matrix for
% palh3m2DE2
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE2_inertiaDJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:11:20
% EndTime: 2020-05-07 04:11:35
% DurationCPUTime: 16.01s
% Computational Cost: add. (11732->590), mult. (5184->847), div. (947->10), fcn. (4599->139), ass. (0->450)
t739 = sin(pkin(16));
t740 = cos(pkin(16));
t745 = sin(pkin(15));
t750 = cos(pkin(15));
t600 = t739 * t750 + t740 * t745;
t601 = -t739 * t745 + t740 * t750;
t742 = sin(qJ(3));
t747 = cos(qJ(3));
t517 = t600 * t747 + t742 * t601;
t743 = sin(qJ(2));
t748 = cos(qJ(2));
t785 = t742 * t600 - t601 * t747;
t1014 = t743 * t517 + t785 * t748;
t516 = t517 * qJD(3);
t1019 = -t1014 * qJD(2) - t743 * t516;
t1018 = -m(6) / 0.2e1;
t752 = rSges(6,1) * pkin(8);
t751 = pkin(10) + rSges(6,3);
t943 = rSges(6,2) * t751;
t1017 = (t752 - t943) * t1018 - Icges(6,6) / 0.2e1;
t1016 = (t752 + t943) * t1018 + Icges(6,6) / 0.2e1;
t916 = -qJ(2) - pkin(15);
t700 = pkin(18) - t916;
t648 = pkin(17) + qJ(3) + t700;
t574 = atan2(-sin(t648), -cos(t648));
t658 = sin(t700);
t659 = cos(t700);
t591 = atan2(t658, -t659);
t886 = pkin(17) - t591;
t514 = -t574 + t886;
t510 = cos(t514);
t765 = rSges(9,2) ^ 2;
t767 = rSges(9,1) ^ 2;
t679 = t767 / 0.2e1 - t765 / 0.2e1;
t709 = -t765 + t767;
t726 = t748 ^ 2;
t917 = t743 * t748;
t855 = rSges(9,2) * t917;
t811 = rSges(9,1) * t855;
t802 = -0.2e1 * t811;
t932 = t658 ^ 2 / t659 ^ 2;
t839 = 0.1e1 + t932;
t590 = 0.1e1 / t839;
t897 = qJD(2) * t590;
t504 = t839 * t897;
t717 = qJD(2) + qJD(3);
t471 = -t717 + t504;
t509 = sin(t514);
t954 = t509 * m(9);
t860 = t471 * t954;
t1015 = (-t709 * t726 + t679 - t802) * t510 * t860;
t787 = t517 * t748 - t743 * t785;
t1013 = -0.2e1 * t917;
t1012 = pkin(1) * qJD(2);
t724 = t743 ^ 2;
t903 = t724 - t726;
t1011 = t903 * t709;
t832 = t903 * rSges(9,2);
t851 = t709 * t917;
t508 = t510 ^ 2;
t858 = m(9) * qJD(2) * t508;
t1010 = -0.4e1 * t1015 - 0.4e1 * (rSges(9,1) * t832 - t851) * t858;
t755 = m(5) + m(6);
t1008 = pkin(4) * t755 + m(4) * rSges(4,1);
t1007 = m(6) * qJD(4) * (-cos(qJ(4)) * rSges(6,2) - sin(qJ(4)) * rSges(6,1));
t1006 = -2 * pkin(12);
t628 = pkin(16) + t648;
t615 = sin(t628);
t616 = cos(t628);
t933 = t615 ^ 2 / t616 ^ 2;
t841 = 0.1e1 + t933;
t922 = t717 / t841;
t495 = t841 * t922;
t491 = qJD(4) - t495;
t486 = qJD(3) + t491;
t482 = qJD(2) + t486;
t473 = -qJD(1) + t482;
t492 = -qJD(4) - t495;
t487 = qJD(3) + t492;
t483 = qJD(2) + t487;
t474 = qJD(1) + t483;
t568 = atan2(-t615, t616);
t563 = t568 + qJ(4);
t558 = qJ(3) + t563;
t554 = qJ(2) + t558;
t543 = -qJ(1) + t554;
t528 = cos(t543);
t719 = qJD(1) - qJD(4);
t732 = qJ(1) - qJ(4);
t697 = cos(t732);
t976 = -t697 / 0.2e1;
t564 = t568 - qJ(4);
t559 = qJ(3) + t564;
t555 = qJ(2) + t559;
t544 = qJ(1) + t555;
t529 = cos(t544);
t980 = t529 / 0.4e1;
t772 = -t473 * t528 / 0.4e1 + t474 * t980 + t719 * t976;
t525 = sin(t544);
t690 = sin(t732);
t524 = sin(t543);
t981 = -t524 / 0.4e1;
t773 = t473 * t981 - t474 * t525 / 0.4e1 + t719 * t690 / 0.2e1;
t472 = qJD(1) + t482;
t475 = -qJD(1) + t483;
t542 = qJ(1) + t554;
t527 = cos(t542);
t731 = qJ(1) + qJ(4);
t696 = cos(t731);
t718 = qJD(1) + qJD(4);
t974 = t718 / 0.2e1;
t545 = -qJ(1) + t555;
t530 = cos(t545);
t979 = -t530 / 0.4e1;
t774 = t472 * t527 / 0.4e1 + t475 * t979 + t696 * t974;
t526 = sin(t545);
t689 = sin(t731);
t523 = sin(t542);
t982 = t523 / 0.4e1;
t775 = t472 * t982 + t475 * t526 / 0.4e1 + t689 * t974;
t1004 = 0.2e1 * (t772 + t774) * rSges(6,2) + 0.2e1 * (t773 + t775) * rSges(6,1);
t1003 = 0.2e1 * (t773 - t775) * rSges(6,2) + 0.2e1 * (-t772 + t774) * rSges(6,1);
t793 = t528 / 0.4e1 + t980 + t976;
t794 = -t527 / 0.4e1 + t979 - t696 / 0.2e1;
t795 = t981 + t525 / 0.4e1 - t690 / 0.2e1;
t796 = t982 - t526 / 0.4e1 + t689 / 0.2e1;
t1002 = 0.2e1 * (t793 - t794) * rSges(6,2) + 0.2e1 * (-t795 + t796) * rSges(6,1);
t1001 = 0.2e1 * (t795 + t796) * rSges(6,2) + 0.2e1 * (t793 + t794) * rSges(6,1);
t1000 = -0.2e1 * t726;
t999 = 0.4e1 * t726;
t758 = -0.2e1 * qJD(2);
t757 = 0.2e1 * qJD(2);
t998 = m(6) / 0.2e1;
t997 = pkin(1) * m(6);
t996 = pkin(1) * m(8);
t995 = pkin(1) * m(9);
t994 = pkin(3) * m(9);
t993 = pkin(4) * m(6);
t992 = m(4) * rSges(4,2);
t991 = m(4) * rSges(4,3);
t990 = m(5) * rSges(5,2);
t989 = m(5) * rSges(5,3);
t988 = m(6) * pkin(12);
t987 = m(7) * rSges(7,3);
t986 = m(9) * rSges(9,1);
t985 = m(9) * rSges(9,3);
t984 = rSges(3,2) * m(3);
t983 = rSges(6,2) * pkin(8);
t766 = rSges(6,2) ^ 2;
t768 = rSges(6,1) ^ 2;
t952 = (-t766 + t768) * m(6);
t978 = Icges(6,1) / 0.4e1 - Icges(6,2) / 0.4e1 - t952 / 0.4e1;
t749 = cos(qJ(1));
t707 = t749 * rSges(6,2);
t675 = qJD(1) * t707;
t977 = t675 / 0.2e1;
t744 = sin(qJ(1));
t706 = t744 * rSges(6,2);
t975 = t706 / 0.2e1;
t857 = rSges(9,1) * t917;
t973 = (t832 - 0.2e1 * t857) * t1012;
t972 = pkin(1) * (-rSges(9,1) * t726 + rSges(9,1) + t855);
t971 = pkin(1) * t1008;
t951 = t751 * m(6);
t641 = t951 - t990;
t970 = pkin(1) * t641;
t671 = pkin(8) * m(6) + m(5) * rSges(5,1);
t969 = pkin(1) * t671;
t733 = qJ(1) + qJ(2);
t691 = sin(t733);
t968 = pkin(1) * t691;
t698 = cos(t733);
t967 = pkin(1) * t698;
t720 = qJD(1) + qJD(2);
t966 = pkin(1) * t720;
t721 = qJD(1) - qJD(2);
t965 = pkin(1) * t721;
t964 = pkin(4) * t641;
t730 = qJ(2) + qJ(3);
t704 = qJ(1) + t730;
t664 = sin(t704);
t963 = pkin(4) * t664;
t669 = cos(t704);
t962 = pkin(4) * t669;
t961 = pkin(4) * t671;
t684 = qJD(1) + t717;
t960 = pkin(4) * t684;
t685 = qJD(1) - t717;
t959 = pkin(4) * t685;
t695 = cos(t730);
t958 = pkin(4) * t695;
t770 = pkin(1) ^ 2;
t957 = m(6) * t770;
t956 = m(9) * t510;
t955 = rSges(9,1) * rSges(9,2);
t953 = t709 * m(9);
t712 = t748 * pkin(1);
t769 = pkin(4) ^ 2;
t950 = t769 * m(5);
t949 = m(5) * qJD(3);
t946 = rSges(4,2) * t743;
t945 = rSges(4,2) * t748;
t512 = -qJ(2) + t514;
t506 = sin(t512);
t942 = rSges(9,2) * t506;
t502 = (t590 * t932 + t590 - 0.1e1) * qJD(2);
t941 = pkin(12) * t502;
t940 = t744 * rSges(6,1);
t939 = t749 * rSges(6,1);
t708 = t751 * rSges(6,1);
t469 = -t717 + t502;
t938 = t469 * t506;
t507 = cos(t512);
t937 = t469 * t507;
t936 = t471 * t510;
t493 = qJD(3) - t495;
t490 = qJD(2) + t493;
t565 = qJ(3) + t568;
t562 = qJ(2) + t565;
t935 = t490 * sin(t562);
t715 = -qJD(4) + qJD(2);
t683 = qJD(3) + t715;
t728 = -qJ(4) + qJ(2);
t702 = qJ(3) + t728;
t931 = sin(t702) * t683;
t930 = cos(t702) * t683;
t677 = m(6) * rSges(6,1) * rSges(6,2) - Icges(6,4);
t929 = 0.4e1 * t677;
t714 = qJD(4) + qJD(2);
t727 = qJ(4) + qJ(2);
t928 = sin(t727) * t714;
t927 = sin(t728) * t715;
t688 = sin(t730);
t926 = t688 * t717;
t925 = cos(t727) * t714;
t924 = cos(t728) * t715;
t923 = t695 * t717;
t920 = t742 * t747;
t915 = t523 - t526;
t914 = t524 - t525;
t913 = t528 + t529;
t912 = t530 + t527;
t556 = qJ(1) + t562;
t535 = sin(t556);
t557 = -qJ(1) + t562;
t536 = sin(t557);
t911 = -t535 - t536;
t539 = cos(t556);
t540 = cos(t557);
t910 = t539 - t540;
t810 = t920 * t950;
t792 = qJD(3) * t810;
t909 = -0.2e1 * t769 * t755 * t688 * t923 - 0.2e1 * t792;
t705 = qJ(1) - t730;
t665 = sin(t705);
t852 = -t959 / 0.2e1;
t853 = t960 / 0.2e1;
t908 = t664 * t853 + t665 * t852;
t670 = cos(t705);
t599 = t670 * t852;
t907 = t669 * t853 + t599;
t630 = -pkin(4) * t665 / 0.2e1;
t906 = t630 + t963 / 0.2e1;
t905 = pkin(4) * t670 / 0.2e1 - t962 / 0.2e1;
t725 = t747 ^ 2;
t904 = t742 ^ 2 - t725;
t902 = 0.2e1 * qJD(3) + t757;
t761 = 0.2e1 * qJ(2);
t901 = 0.2e1 * qJ(3) + t761;
t729 = t761 + qJ(3);
t900 = t766 + t768;
t899 = qJD(1) * t744;
t898 = qJD(2) * t509;
t896 = qJD(2) * t743;
t895 = qJD(2) * t748;
t894 = qJD(3) * t747;
t893 = qJD(3) * t748;
t889 = pkin(1) * pkin(4) * m(5);
t660 = t742 * t889;
t713 = m(9) * t955;
t891 = (t660 + 0.2e1 * t713 - 0.2e1 * t810) * qJD(2);
t888 = pkin(1) * t993;
t885 = pkin(1) * t992;
t503 = t757 - t504;
t884 = t503 * t996;
t883 = t504 * t996;
t882 = (t758 + t471) * t995;
t881 = (t758 - t717 + (0.2e1 + 0.2e1 * t932) * t897) * t994;
t880 = t994 * t717;
t878 = t482 * t988;
t877 = t483 * t988;
t876 = t507 * t986;
t875 = -t997 / 0.2e1;
t874 = t997 / 0.2e1;
t873 = -t993 / 0.2e1;
t871 = t993 / 0.2e1;
t547 = -0.2e1 * t922;
t868 = (m(4) * t945 - t1008 * t743) * t695 * t1012;
t867 = t698 * t966;
t734 = qJ(1) - qJ(2);
t692 = sin(t734);
t866 = t692 * t965;
t865 = t664 * t960;
t864 = t669 * t960;
t863 = pkin(4) * t926;
t861 = m(4) * t923;
t859 = t725 * t950;
t856 = rSges(9,2) * t937;
t854 = pkin(1) * t896;
t850 = -t939 / 0.2e1;
t849 = 0.2e1 * t568 + t901;
t848 = 0.2e1 * t953 - 0.4e1 * t859 + 0.2e1 * t950;
t561 = t568 + t729;
t560 = t568 + t901;
t847 = t743 * t895;
t846 = t664 * t888;
t845 = t669 * t888;
t844 = 0.2e1 * m(8) * t941;
t762 = -0.2e1 * qJ(2);
t842 = t762 + t886;
t838 = qJD(3) * t889;
t837 = -t899 / 0.2e1;
t836 = t899 / 0.2e1;
t833 = t903 * rSges(9,1);
t831 = m(4) * (-rSges(4,1) * t748 - t946) * t688 * t1012;
t830 = t936 * t995;
t829 = (-t833 - 0.2e1 * t855) * t954 * t1012;
t828 = m(6) * t665 * t959;
t827 = rSges(6,1) * t874;
t826 = rSges(6,2) * t875;
t825 = rSges(6,1) * t873;
t824 = rSges(6,2) * t871;
t823 = t938 * t986;
t822 = m(9) * t856;
t821 = (t757 + t486) * t874;
t820 = t486 * t874;
t488 = -t495 + t902;
t818 = (qJD(4) + t488) * t873;
t817 = t491 * t873;
t814 = -0.2e1 * t709 * t507 * t506;
t813 = pkin(1) * (m(4) * t946 + t1008 * t748) * t926;
t809 = t866 / 0.2e1;
t808 = t547 * t933 + t547 + t902;
t807 = t684 * t846;
t806 = t691 * t845;
t804 = 0.2e1 * t823;
t803 = t742 * t838;
t618 = t691 * t966;
t801 = -t618 + t865;
t800 = pkin(1) * (-rSges(4,1) * t743 + t945) * t861;
t799 = (rSges(9,2) * t726 - rSges(9,2) + t857) * t830;
t798 = -t963 + t968;
t797 = -t962 + t967;
t791 = t720 * t806;
t790 = t698 * t807;
t789 = t684 * t806;
t788 = t720 * t698 * t846;
t521 = -t743 * t894 - t747 * t896 + (-t893 - t895) * t742;
t604 = t743 * t742 - t748 * t747;
t522 = t717 * t604;
t786 = t521 * t750 + t745 * t522;
t605 = t748 * t742 + t743 * t747;
t784 = t604 * t750 + t745 * t605;
t783 = t822 * t712;
t782 = -t713 + t810;
t781 = t904 * t769 * t949;
t780 = m(9) * t854 * t942;
t778 = -t864 + t867;
t776 = m(9) * t469 * t814 + t909;
t701 = qJ(3) + t727;
t661 = sin(t701);
t666 = cos(t701);
t682 = qJD(3) + t714;
t771 = -(rSges(4,2) * t991 - Icges(4,6)) * t926 - (rSges(9,2) * t985 - Icges(9,6)) * t938 + (-rSges(9,1) * t985 + Icges(9,5)) * t937 + (pkin(4) * t989 + rSges(4,1) * t991 - Icges(4,5)) * t923 + t824 * t930 + t825 * t931 + (rSges(6,1) * t661 + rSges(6,2) * t666) * t682 * t871;
t764 = -4 * Icges(6,5);
t763 = -0.2e1 * Icges(6,2);
t759 = 0.2e1 * qJ(4);
t735 = 0.4e1 * t950;
t722 = pkin(17) + pkin(18);
t711 = m(4) + m(5) + m(8) + m(9);
t703 = pkin(14) + t916;
t699 = cos(t734);
t681 = cos(t722);
t680 = sin(t722);
t676 = 0.2e1 * t730;
t674 = -t940 / 0.2e1;
t668 = cos(t703);
t663 = sin(t703);
t657 = pkin(1) * t699;
t656 = pkin(1) * t692;
t655 = 0.4e1 * t708;
t652 = 0.4e1 * t747 * t889;
t651 = 0.2e1 * t703;
t650 = qJD(1) * t850;
t643 = t747 * t838;
t639 = 0.4e1 * t953;
t636 = -0.4e1 * t803;
t623 = t706 - t939;
t622 = t706 + t939;
t621 = t707 - t940;
t620 = t707 + t940;
t619 = t699 * t965;
t612 = t674 - t707 / 0.2e1;
t611 = t674 + t707 / 0.2e1;
t610 = t975 + t850;
t609 = t975 + t939 / 0.2e1;
t608 = t712 - t958;
t585 = qJ(2) + t591;
t584 = t761 + t591;
t581 = -qJ(2) + t886;
t578 = 0.2e1 * t585;
t577 = -t851 + (t1000 + 0.1e1) * t955;
t571 = 0.2e1 * t851 + (t999 - 0.2e1) * t955;
t553 = -qJ(4) + t561;
t552 = -qJ(4) + t560;
t551 = qJ(4) + t560;
t550 = -qJ(4) + t849;
t549 = qJ(4) + t849;
t548 = qJ(4) + t561;
t541 = 0.2e1 * t562;
t538 = cos(t555);
t537 = cos(t554);
t534 = sin(t555);
t533 = sin(t554);
t532 = 0.2e1 * t555;
t531 = 0.2e1 * t554;
t520 = t745 * t604 - t605 * t750;
t515 = t785 * qJD(3);
t513 = t762 - t574 - 0.2e1 * t591 + 0.2e1 * pkin(17);
t511 = -t574 + t842;
t505 = 0.2e1 * t512;
t499 = -t745 * t521 + t522 * t750;
t497 = t520 * t739 - t784 * t740;
t496 = t520 * t740 + t784 * t739;
t467 = t499 * t740 - t786 * t739;
t466 = t499 * t739 + t786 * t740;
t465 = t496 * t680 + t497 * t681;
t464 = t496 * t681 - t497 * t680;
t463 = -t787 * qJD(2) + t743 * t515 - t516 * t748;
t462 = 0.1e1 / t465 ^ 2;
t459 = -t1014 * t681 - t680 * t787;
t458 = t1014 * t680 - t787 * t681;
t457 = 0.1e1 / t459 ^ 2;
t452 = atan2(t464, t465) + t730;
t451 = sin(t452);
t449 = atan2(t458, t459) + t730;
t448 = sin(t449);
t447 = t620 * t697 + t621 * t696 + t689 * t622 + t623 * t690 + t913 * t612 + t912 * t611 + t914 * t610 + t915 * t609 + (t911 * t744 - t910 * t749) * rSges(6,3);
t444 = (-t620 * t690 + t623 * t697) * t719 + (-t621 * t689 + t696 * t622) * t718 + (t689 + t690) * t675 + t912 * (rSges(6,2) * t837 + t650) + t913 * (rSges(6,2) * t836 + t650) + t915 * (rSges(6,1) * t837 + t977) + t914 * (rSges(6,1) * t836 + t977) + (-t611 * t526 - t609 * t530) * t475 + (-t612 * t525 - t610 * t529) * t474 + (-t612 * t524 + t610 * t528) * t473 + (-t611 * t523 + t609 * t527) * t472 + ((-t696 - t697) * t706 + ((-t696 + t697) * t749 + (-t689 + t690) * t744) * rSges(6,1)) * qJD(1) + ((-t536 * t749 - t540 * t744) * (-qJD(1) + t490) + (t535 * t749 - t539 * t744) * (qJD(1) + t490) + (t910 * t744 + t911 * t749) * qJD(1)) * rSges(6,3);
t1 = [(-(t639 - (4 * Icges(9,1)) + (4 * Icges(9,2))) * sin(t505) / 0.4e1 + 0.2e1 * (-Icges(9,4) + t713) * cos(t505)) * t469 + (((t655 - 0.4e1 * t983) * m(6) + t764) * cos(t549) / 0.8e1 + sin(t549) * t1016) * (qJD(4) + t808) + (((t655 + 0.4e1 * t983) * m(6) + t764) * cos(t550) / 0.8e1 + sin(t550) * t1017) * (-qJD(4) + t808) + qJD(3) * t742 * t971 + (sin(t564) * t825 + cos(t564) * t824) * t492 + (-0.2e1 * (-m(8) * rSges(8,1) * rSges(8,2) + Icges(8,4)) * cos(t578) - pkin(3) ^ 2 * m(9) * sin(0.2e1 * t581) + ((rSges(8,1) ^ 2 - rSges(8,2) ^ 2) * m(8) - Icges(8,1) + Icges(8,2)) * sin(t578)) * t502 + (cos(t552) * t824 + sin(t552) * t825) * (-qJD(4) + t488) + (cos(t553) * t826 + sin(t553) * t827) * (t757 + t487) + (sin(t561) * t969 - cos(t561) * t970) * (t757 + t493) + (sin(t729) * t971 + cos(t729) * t885) * (t757 + qJD(3)) + ((-t504 * sin(t886) + t503 * sin(t842)) * pkin(1) - 0.2e1 * sin(t581) * t941) * t994 - (-sin(t568) * t961 + cos(t568) * t964) * t495 + (-(0.4e1 * m(6) * t769 + t735 + 0.4e1 * (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - (4 * Icges(4,1)) + (4 * Icges(4,2))) * sin(t676) / 0.4e1 + 0.2e1 * (-rSges(4,1) * t992 + Icges(4,4)) * cos(t676)) * t717 + pkin(8) * t1007 + (cos(t563) * t817 + cos(t548) * t821 + cos(t551) * t818 + cos(t558) * t820 + t537 * t878 - t538 * t877) * rSges(6,2) + (-rSges(3,1) * t984 + Icges(3,4)) * cos(t761) * t757 + ((m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t651) + ((rSges(3,1) * m(3) + (m(6) + t711) * pkin(1)) * t743 + t748 * t984) * pkin(12)) * t758 + (cos(t591) * t883 - cos(t584) * t884 + cos(t585) * t844) * rSges(8,2) + (sin(t585) * t844 + sin(t591) * t883 - sin(t584) * t884) * rSges(8,1) + (sin(t558) * t820 + sin(t563) * t817 + sin(t548) * t821 + sin(t551) * t818 + t533 * t878 + t534 * t877) * rSges(6,1) + (-(0.2e1 * (-0.2e1 * (pkin(8) + t751) * (-pkin(8) + t751) + t900) * m(6) + 0.4e1 * (rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - (4 * Icges(5,1)) - 0.2e1 * Icges(6,1) + (4 * Icges(5,2)) + t763 + (4 * Icges(6,3))) * sin(t541) / 0.4e1 + 0.2e1 * (pkin(8) * t951 - rSges(5,1) * t990 + Icges(5,4)) * cos(t541) + t641 * cos(t562) * t1006) * t490 + t885 * t894 + (pkin(1) * t860 + sin(t574) * t880 + sin(t511) * t882 + sin(t513) * t881) * rSges(9,1) + pkin(12) * t804 + (cos(t532) * t929 / 0.8e1 + sin(t532) * t978) * t483 + (-cos(t531) * t929 / 0.8e1 + sin(t531) * t978) * t482 + (sin(t565) * t969 - cos(t565) * t970) * t493 + (-(0.2e1 * Icges(6,1) + t763 - 0.2e1 * t952) * sin(t759) / 0.4e1 + t677 * cos(t759)) * qJD(4) + (-sin(t560) * t961 + cos(t560) * t964) * t488 + (cos(t559) * t826 + sin(t559) * t827) * t487 + t822 * t1006 + (cos(t574) * t880 - cos(t511) * t882 - cos(t513) * t881 - t830) * rSges(9,2) + (0.2e1 * (-rSges(7,1) * t663 + rSges(7,2) * t668) * pkin(6) * m(7) - (t957 - Icges(3,1) + Icges(3,2) + t711 * t770 + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3)) * sin(t761) + ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t651)) * qJD(2) + 0.2e1 * (rSges(4,2) * t861 + t1008 * t926 + t671 * t935) * pkin(12); (Icges(3,5) * t748 - t743 * Icges(3,6) + (-rSges(3,1) * t748 + rSges(3,2) * t743) * rSges(3,3) * m(3) + (-rSges(7,2) * t987 + Icges(7,6)) * t663 - (rSges(7,1) * t987 - Icges(7,5)) * t668 + (-m(8) * rSges(8,3) - t985 - t989 - t991) * t712) * qJD(2) + ((-t924 / 0.2e1 - t925 / 0.2e1) * rSges(6,2) + (t927 / 0.2e1 - t928 / 0.2e1) * rSges(6,1)) * t997 + t771; t1010 + (0.2e1 * t724 + t1000) * (t660 - t782) * qJD(2) + 0.2e1 * (t577 * t936 + t973 + (0.4e1 * t811 + t1011) * t898) * t956 + (t778 * t692 + (t798 * t721 + t801) * t699) * t875 + t791 / 0.2e1 - t789 / 0.2e1 + t790 / 0.2e1 - t788 / 0.2e1 + (t636 + 0.8e1 * t792) * t726 / 0.2e1 - 0.2e1 * t783 - 0.2e1 * t868 - 0.2e1 * (t577 * t509 - t972) * t860 - 0.2e1 * (-t854 - t856) * t876 + t776 + ((t619 - t778) * t630 + (t656 - t798) * t599 + t797 * t809) * m(6) + (t643 + t781) * t1013 + (-t801 - t866) * t670 * t873 + 0.2e1 * t780 + 0.2e1 * t803 + 0.2e1 * t800 + 0.2e1 * t799 + (t657 - t797) * t828 / 0.2e1 - (t652 + t848 + 0.2e1 * t957) * t847 + (t712 - t942) * t804 + 0.2e1 * t831 + 0.2e1 * t829 + 0.2e1 * t813; t771; -(t571 * t936 - t973 + (-0.8e1 * t811 - 0.2e1 * t1011) * t898) * t956 + (-t891 + t636 / 0.4e1 + 0.4e1 * t792) * t726 + (t721 * t846 / 0.4e1 - t807 / 0.4e1) * t699 - (t643 + 0.2e1 * t781) * t917 + 0.4e1 * (-rSges(9,2) * t833 + t851) * t858 + t791 / 0.4e1 - t789 / 0.4e1 + t790 / 0.4e1 - t788 / 0.4e1 - t783 - (t639 + t652 + t735 - 0.8e1 * t859) * t847 / 0.2e1 - t868 + (t571 * t509 + t972) * t860 - (-t854 - 0.2e1 * t856) * t876 + t776 + (t684 / 0.4e1 - t721 / 0.4e1) * t692 * t845 + t780 + t724 * t891 + t803 + t800 + t799 + (t657 + 0.2e1 * t962 - t967) * t828 / 0.4e1 + (t712 - 0.2e1 * t942) * t823 + t831 + t829 - 0.4e1 * t1015 + t813 - ((t619 + 0.2e1 * t864 - t867) * t665 + ((t656 + 0.2e1 * t963 - t968) * t685 + t618 - 0.2e1 * t865 - t866) * t670) * t993 / 0.4e1; (-0.2e1 * t782 * t903 - t848 * t917) * qJD(2) + (0.4e1 * (t509 ^ 2 - t508) * (t679 * t917 + (t726 - 0.1e1 / 0.2e1) * t955) * t471 + (t814 + 0.2e1 * (-t506 ^ 2 + t507 ^ 2) * t955) * t469 - 0.4e1 * (-t903 * t679 + t802) * t510 * t898) * m(9) + ((t904 * t1013 + t920 * t999) * t949 + (-t664 * t670 + t665 * t669) * m(6) * (t685 / 0.2e1 - t684 / 0.2e1)) * t769 + t909 + t1010; -(t900 * m(6) + Icges(6,3)) * t935 + (t534 * t1017 + ((t708 + t983) * m(6) - Icges(6,5)) * t538 / 0.2e1) * t483 + (t533 * t1016 + ((t708 - t983) * m(6) - Icges(6,5)) * t537 / 0.2e1) * t482 - pkin(12) * t1007 + ((t666 * t682 - t930) * rSges(6,2) + (t661 * t682 + t931) * rSges(6,1)) * t873 + ((t924 - t925) * rSges(6,2) + (-t927 - t928) * rSges(6,1)) * t875; ((t656 / 0.2e1 - t968 / 0.2e1 + t906) * t1003 + (t619 / 0.2e1 - t867 / 0.2e1 + t907) * t1002 + (-t657 / 0.2e1 + t967 / 0.2e1 + t905) * t1004 + (t809 - t618 / 0.2e1 + t908) * t1001 + t444 * t451 * t608 + ((((-t466 * t680 + t467 * t681) / t465 - (t466 * t681 + t467 * t680) * t464 * t462) / (t464 ^ 2 * t462 + 0.1e1) + t717) * cos(t452) * t608 + t451 * (-t854 + t863)) * t447) * t998; (-t444 * t448 * t958 + t908 * t1001 + t907 * t1002 + t906 * t1003 + t905 * t1004 + (-((((t785 * t893 - t1019) * t681 - t463 * t680) / t459 - (t463 * t681 - t680 * (-t515 * t748 + t1019)) * t458 * t457) / (t458 ^ 2 * t457 + 0.1e1) + t717) * cos(t449) * t958 + t448 * t863) * t447) * t998; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
