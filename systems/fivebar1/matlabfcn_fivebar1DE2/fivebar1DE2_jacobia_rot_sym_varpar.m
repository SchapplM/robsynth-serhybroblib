% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% fivebar1DE2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in fivebar1DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% Ja_rot [3x2]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = fivebar1DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_jacobia_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fivebar1DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_jacobia_rot_sym_varpar: pkin has to be [5x1] (double)');
Ja_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 05:51:21
	% EndTime: 2020-04-27 05:51:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 05:51:21
	% EndTime: 2020-04-27 05:51:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 1, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 05:51:40
	% EndTime: 2020-04-27 05:51:57
	% DurationCPUTime: 16.43s
	% Computational Cost: add. (171380->621), mult. (514126->964), div. (1274->9), fcn. (61360->6), ass. (0->422)
	t774 = (pkin(2) ^ 2);
	t770 = (pkin(3) ^ 2);
	t776 = pkin(1) ^ 2;
	t921 = t770 - t776;
	t672 = t921 * t774;
	t1009 = 2 * pkin(1);
	t1008 = 2 * pkin(3);
	t1007 = 4 * pkin(3);
	t733 = sin(qJ(2));
	t734 = sin(qJ(1));
	t953 = t733 * t734;
	t661 = pkin(2) * t953;
	t825 = pkin(3) * t661;
	t656 = -0.2e1 * t825;
	t736 = cos(qJ(1));
	t978 = pkin(2) * t736;
	t902 = pkin(1) * t978;
	t671 = -0.2e1 * t902;
	t919 = (t774 + t776);
	t844 = t770 + t919;
	t666 = pkin(1) - t978;
	t735 = cos(qJ(2));
	t650 = t666 * t735;
	t873 = pkin(3) * t650;
	t586 = t656 + t671 + t844 + 0.2e1 * t873;
	t1006 = -0.2e1 / t586 ^ 2;
	t689 = pkin(3) * t735;
	t638 = t666 + t689;
	t996 = 0.8e1 * t638;
	t1005 = 2 * t770;
	t1004 = 8 * t774;
	t703 = t736 ^ 2;
	t955 = t703 * t774;
	t901 = pkin(1) * t689;
	t669 = 0.2e1 * t901;
	t683 = -t770 / 0.3e1 + t776;
	t700 = t735 ^ 2;
	t942 = t770 * t700;
	t617 = 0.4e1 / 0.3e1 * t942 + t669 + t683;
	t763 = pkin(5) ^ 2;
	t766 = pkin(4) ^ 2;
	t691 = -t763 - t766;
	t759 = 0.3e1 * t776;
	t675 = t759 + t691;
	t964 = t675 * t770;
	t649 = 0.10e2 * t964;
	t748 = 7 * t770;
	t782 = t774 ^ 2;
	t752 = 5 * t782;
	t768 = t770 ^ 2;
	t781 = pkin(2) * t774;
	t771 = t781 ^ 2;
	t775 = t776 ^ 2;
	t701 = t703 ^ 2;
	t958 = t701 * t782;
	t961 = t691 * t776;
	t693 = t770 + t776;
	t686 = t693 ^ 2;
	t916 = t776 - t763;
	t843 = t770 + t916;
	t962 = t686 * (-t766 + t843);
	t1003 = 0.7e1 * t771 + (t748 + t675) * t752 + (t649 + (21 * t768) + 0.9e1 * t775 + 0.6e1 * t961) * t774 + t962 - 0.24e2 * t617 * t958;
	t695 = -(3 * t770) + t776;
	t885 = 0.4e1 * t942;
	t1002 = t695 + t885;
	t714 = -t763 / 0.3e1;
	t1001 = t714 - t766 / 0.3e1;
	t762 = t763 ^ 2;
	t765 = t766 ^ 2;
	t1000 = -t762 / 0.6e1 + t765 / 0.6e1;
	t999 = -0.4e1 * pkin(2);
	t687 = 0.10e2 / 0.3e1 * t770;
	t922 = t768 + t775;
	t760 = 0.2e1 * t776;
	t925 = t760 - t763;
	t941 = t776 * t763;
	t614 = t925 * t770 - t1000 + t922 - t941;
	t796 = t614 + t782;
	t588 = (t687 + t925) * t774 + t796;
	t998 = -0.6e1 * t588;
	t893 = 0.2e1 * t955;
	t920 = -t774 + t776;
	t622 = t671 + t893 + t920;
	t997 = -0.2e1 * t622;
	t699 = t735 * t700;
	t995 = -0.8e1 * t699;
	t994 = 0.4e1 * t700;
	t993 = -0.8e1 * t735;
	t992 = 0.8e1 * t735;
	t991 = 4 * t768;
	t749 = 6 * t770;
	t990 = 2 * t774;
	t989 = pkin(1) * pkin(3);
	t988 = -pkin(4) - pkin(5);
	t987 = -pkin(4) + pkin(5);
	t728 = t770 / 0.3e1;
	t604 = -0.4e1 / 0.9e1 * t825 + t776 + t774 / 0.3e1 + t728 + t766 / 0.9e1 - t763 / 0.9e1;
	t712 = -t763 / 0.6e1;
	t731 = t774 / 0.2e1;
	t928 = t731 + t776;
	t799 = -t825 + t928;
	t612 = t766 / 0.6e1 + t712 + t799;
	t722 = t766 / 0.3e1;
	t637 = t714 + t722 + t844;
	t915 = 4 * pkin(1);
	t946 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
	t785 = pkin(3) * t770;
	t959 = t699 * t785;
	t559 = t683 * t656 + 0.6e1 * t604 * t942 + t637 * t946 + (t612 * t689 + t959) * t915;
	t812 = -0.4e1 * t825;
	t716 = -0.2e1 / 0.3e1 * t763;
	t721 = 0.2e1 / 0.3e1 * t766;
	t847 = t716 + t721 + t760;
	t846 = t716 + t693;
	t934 = t693 * (t721 + t846) + t782;
	t571 = t637 * t812 + (t749 + t847) * t774 + t934;
	t610 = t656 + t637;
	t673 = t920 * t770;
	t855 = 0.4e1 * t901;
	t560 = t610 * t855 + t673 * t994 + t571;
	t684 = t776 - t774 / 0.3e1;
	t623 = t684 * t656;
	t945 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
	t585 = t637 * t945 + t623;
	t587 = (t687 + t847) * t774 + t934;
	t621 = t669 + t1002;
	t754 = -3 * t774;
	t696 = t754 + t776;
	t883 = pkin(1) * t959;
	t833 = 0.8e1 * t883;
	t641 = t696 * t833;
	t659 = t766 + t843;
	t741 = 15 * t768;
	t742 = 15 * t770;
	t756 = 0.3e1 * t775;
	t667 = pkin(1) + t689;
	t702 = t736 * t703;
	t957 = t702 * t781;
	t861 = t667 * t957;
	t819 = -0.8e1 * t861;
	t923 = -t763 + t766;
	t845 = t759 + t923;
	t854 = 0.6e1 * t901;
	t864 = 0.12e2 * t942;
	t871 = pkin(3) * t953;
	t965 = t667 * t736;
	t545 = t621 * t819 + t641 + t585 * t864 + t571 * t854 + t771 + (t742 + t845) * t782 + t686 * t659 + (0.12e2 * t559 * t703 + t845 * t749 + t923 * t760 + t741 + t756) * t774 + 0.6e1 * (-t560 * t965 - t587 * t871) * pkin(2);
	t750 = 3 * t770;
	t678 = t750 + t919;
	t690 = pkin(2) * t734;
	t642 = t678 * t690;
	t975 = pkin(3) * t733;
	t647 = t690 - t975;
	t877 = t770 * t690;
	t815 = t700 * t877;
	t753 = 3 * t774;
	t676 = t753 + t693;
	t963 = t676 * t733;
	t980 = pkin(1) * t735;
	t572 = -0.2e1 * t815 + t642 + (0.2e1 * t647 * t980 - t963) * pkin(3);
	t692 = t1005 + t774;
	t951 = t733 * t785;
	t639 = t877 - t951;
	t966 = t639 * t700;
	t573 = t692 * t975 + 0.2e1 * t966 + (-t921 + t669) * t690;
	t603 = -pkin(3) * t963 + t642;
	t674 = pkin(2) * t991 + 0.8e1 * t770 * t781;
	t857 = t785 * t945;
	t605 = t674 * t734 + 0.4e1 * t733 * t857;
	t747 = 5 * t768;
	t917 = t775 + t782;
	t743 = 10 * t770;
	t926 = t743 + t760;
	t940 = t776 * t770;
	t620 = t774 * t926 + t747 + t917 + 0.6e1 * t940;
	t628 = t752 + (t743 + 0.6e1 * t776) * t774 + t686;
	t818 = 0.8e1 * t861;
	t894 = -0.4e1 * t955;
	t549 = t573 * t894 + t605 * t700 + (-0.4e1 * t603 * t980 + (t628 + t818) * t733) * pkin(3) + (0.4e1 * t572 * t965 + (-t620 + t833) * t734) * pkin(2);
	t842 = t774 + t916;
	t807 = t770 + t842;
	t653 = -t766 + t807;
	t619 = t671 + t653;
	t600 = t656 + t619;
	t803 = t700 * t997 - t916;
	t777 = sqrt(0.4e1 * t672 * t703 + 0.4e1 * t653 * t902 - t768 - (t776 + (pkin(2) - t987) * (pkin(2) + t987)) * (t776 + (pkin(2) - t988) * (pkin(2) + t988)) + (t754 + t766 + t803) * t1005 + (-t600 * t650 + t619 * t661) * t1007);
	t539 = t545 * t638 + t549 * t777;
	t986 = 0.1e1 / t539 / 0.4e1;
	t985 = -0.1e1 / t539 ^ 2 / 0.4e1;
	t984 = 0.4e1 / 0.3e1 * t774;
	t983 = -t777 / 0.4e1;
	t982 = t777 / 0.4e1;
	t981 = pkin(1) * (t661 - pkin(3));
	t979 = pkin(1) * t768;
	t977 = pkin(3) * t666;
	t976 = pkin(3) * t700;
	t932 = t762 / 0.2e1 - t765 / 0.2e1;
	t811 = -0.3e1 * t941 + t756 + t932;
	t717 = -0.3e1 / 0.2e1 * t763;
	t929 = t717 + t759;
	t935 = t693 * ((t717 + t760) * t770 - 0.3e1 / 0.2e1 * t941 + t922 + t932) + t771;
	t558 = t825 * t998 + (t741 + (-0.9e1 * t763 + 0.18e2 * t776) * t770 + t811) * t774 + (t742 + t929) * t782 + t935;
	t974 = t558 * pkin(3);
	t745 = -0.5e1 * t763;
	t746 = 7 * t768;
	t568 = (t748 + t929) * t782 + (t746 + (t745 + 0.10e2 * t776) * t770 + t811) * t774 + t935;
	t973 = t568 * pkin(3);
	t972 = t683 * pkin(3);
	t645 = t776 + t774 / 0.4e1 + t770 / 0.4e1 - t763 / 0.8e1;
	t927 = 0.4e1 / 0.7e1 * t776 - t763 / 0.7e1;
	t565 = -0.32e2 / 0.21e2 * t645 * t825 + t782 / 0.7e1 + (0.16e2 / 0.21e2 * t770 + t927) * t774 + t768 / 0.7e1 + t927 * t770 + t775 - 0.3e1 / 0.7e1 * t941 + t762 / 0.42e2 - t765 / 0.42e2;
	t713 = -t763 / 0.4e1;
	t646 = t713 + t728 + t928;
	t727 = 0.4e1 / 0.3e1 * t770;
	t569 = -0.8e1 / 0.3e1 * t646 * t825 + t782 / 0.3e1 + (t727 + t714) * t774 + t775 - t768 / 0.3e1 + (t984 + 0.2e1 / 0.3e1 * t770 + t716) * t776 + t762 / 0.18e2 - t765 / 0.18e2;
	t729 = t770 / 0.2e1;
	t616 = -0.2e1 / 0.3e1 * t825 + t776 + t729 + t713;
	t758 = 0.4e1 * t776;
	t679 = (t758 + t763) * t770;
	t685 = t776 - 0.2e1 / 0.3e1 * t774;
	t715 = -t763 / 0.2e1;
	t658 = t715 + t844;
	t800 = t658 * t812;
	t868 = 0.16e2 * t959;
	t698 = t700 ^ 2;
	t943 = t768 * t698;
	t889 = 0.8e1 * t943;
	t548 = t685 * t889 + 0.14e2 * t565 * t942 + t683 * t800 - t921 * t782 + (t679 - 0.10e2 / 0.3e1 * t768 + 0.2e1 * t775 - t941) * t774 + t614 * t946 + (0.6e1 * t569 * t689 + t616 * t868) * pkin(1);
	t574 = t800 + (t749 + t925) * t774 + t796;
	t591 = t658 * t945 + t623;
	t550 = t574 * t854 + t591 * t864 + t558 + t641;
	t744 = -0.2e1 * t763;
	t757 = 0.8e1 * t776;
	t829 = t951 * t999;
	t609 = t734 * t829 + t991 + ((4 * t774) + t744 + t757) * t770;
	t615 = t713 - t770 + t799;
	t897 = 0.8e1 * t959;
	t904 = 0.4e1 * t689;
	t561 = t656 * t946 + t609 * t700 + t658 * t695 + (t615 * t904 + t897) * pkin(1);
	t918 = t775 - t768;
	t563 = t684 * t800 - t771 + (-t687 - t916) * t782 + (t679 + t918 + t1000) * t774 + t776 * t614;
	t778 = pkin(1) * t776;
	t665 = -(12 * pkin(1) * t785) + t778 * t1007;
	t682 = -(8 * t768) + 0.12e2 * t940;
	t822 = pkin(1) * t868;
	t579 = t665 * t735 + t682 * t700 + t822 + t889 + t922 - 0.6e1 * t940;
	t594 = t656 * t945 + t658 * t696;
	t651 = (-0.6e1 * t774 * t776 + t917) * t768;
	t688 = -0.30e2 * t763 + 0.60e2 * t776;
	t924 = t762 - t765;
	t810 = -0.6e1 * t941 + 0.6e1 * t775 + t924;
	t536 = -0.32e2 * t561 * t861 + 0.16e2 * t651 * t698 + 0.24e2 * t563 * t942 + (t744 + t758 + (28 * t770)) * t771 + t659 * t962 + (0.24e2 * t548 * t703 + t688 * t768 + t810 * t749 + t924 * t760 - 0.6e1 * t775 * t763 + 0.4e1 * t778 ^ 2 + (28 * t785 ^ 2)) * t774 + 0.8e1 * (-t550 * t965 - t568 * t871) * pkin(2) + (0.32e2 * t594 * t959 + t974 * t992) * pkin(1) + (0.16e2 * t579 * t701 + t688 * t770 + (70 * t768) + t782 + t810) * t782;
	t931 = t712 - t766 / 0.6e1;
	t850 = t776 + t931;
	t630 = t984 + t729 + t850;
	t809 = t731 + t850;
	t631 = t727 + t809;
	t580 = -t630 * t975 + t631 * t690;
	t834 = 0.20e2 / 0.3e1 * t770;
	t852 = 0.2e1 / 0.3e1 * t763 + t721 + t758;
	t853 = 0.4e1 / 0.3e1 * t763 + 0.4e1 / 0.3e1 * t766 - 0.2e1 * t776;
	t589 = -t782 + (-t834 + t852) * t774 - (3 * t768) + t853 * t770 + t775;
	t633 = t770 + t809;
	t677 = t990 - t921;
	t863 = -t975 / 0.2e1;
	t592 = t633 * t690 + t677 * t863;
	t849 = t776 + t1001;
	t851 = t763 / 0.3e1 + t722 + t760;
	t806 = -0.8e1 / 0.3e1 * t943 - t672 - 0.5e1 / 0.3e1 * t768 + t851 * t770 + t776 * t849;
	t960 = t699 * t768;
	t551 = t580 * t885 + t589 * t863 + t806 * t690 + (t592 * t689 - t733 * t960) * t915;
	t726 = -0.2e1 / 0.3e1 * t766;
	t930 = t716 + t726;
	t581 = t782 + (t926 + t930) * t774 + t747 + 0.2e1 * t964 + t776 * (t776 + t930);
	t593 = t752 + ((5 * t770) + t675) * t990 + t693 * (t726 + t846);
	t562 = t581 * t690 - t593 * t975;
	t808 = t770 + t849;
	t635 = t753 + t808;
	t636 = t678 + t1001;
	t582 = -t635 * t975 + t636 * t690;
	t644 = t729 + t774 + t931;
	t912 = 0.2e1 * t690;
	t597 = t644 * t912 + t945 * t975;
	t876 = t785 * t690;
	t887 = -0.4e1 * t942;
	t552 = t597 * t887 + (t582 * t904 + t876 * t995) * pkin(1) + t562;
	t590 = -(3 * t782) + (-t834 + t853) * t774 + t852 * t770 + t918;
	t598 = -0.5e1 / 0.3e1 * t782 + (-t770 + t851) * t774 + t776 * t808;
	t907 = -0.2e1 * t975;
	t564 = t590 * t690 + t598 * t907;
	t848 = t715 - t766 / 0.2e1 + t776;
	t632 = 0.3e1 / 0.2e1 * t774 + t750 + t848;
	t648 = t690 + 0.2e1 * t975;
	t567 = t695 * t690 + 0.4e1 * t966 + (t632 * t733 + t648 * t980) * t1008;
	t634 = t753 + 0.3e1 / 0.2e1 * t770 + t848;
	t596 = t634 * t690 + t696 * t975 / 0.2e1;
	t801 = 0.24e2 * t684 * t943 - t771 - ((21 * t770) + t675) * t782 - (t649 + t756 + (35 * t768) + 0.2e1 * t961) * t774 - (t746 + (t745 + t757 - 0.5e1 * t766) * t770 + t776 * (-t766 + t916)) * t693;
	t866 = -0.12e2 * t955;
	t888 = -0.6e1 * t942;
	t540 = t596 * t822 + t567 * t818 + t551 * t866 + t564 * t888 + (t1003 * t733 - 0.6e1 * t562 * t980) * pkin(3) + (0.6e1 * t552 * t965 + t734 * t801) * pkin(2);
	t533 = t536 * t638 + t540 * t777;
	t618 = t671 + t766 + t807;
	t874 = t622 * t689;
	t575 = t618 * t666 + 0.2e1 * t874;
	t578 = t618 * t735 + (t994 - 0.2e1) * t977;
	t933 = -t661 + t650;
	t611 = pkin(3) + t933;
	t547 = t575 * t733 + t578 * t690 + t611 * t777;
	t652 = t750 + t766 + t842;
	t601 = t652 + t671 + t812;
	t950 = t735 * t734;
	t613 = pkin(2) * t950 + t666 * t733;
	t968 = t613 * t777;
	t546 = -t601 * t650 + t968 + (t652 * t953 - 0.2e1 * t736 * t981) * pkin(2) + (-t753 - t766 - t770 + t803 + t893) * pkin(3);
	t841 = t546 * t986;
	t797 = t533 * t841 + t547 * t983;
	t583 = 0.1e1 / t586;
	t944 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
	t862 = t583 * t944;
	t531 = t797 * t862;
	t949 = t735 * t736;
	t626 = -t949 - t953;
	t530 = t531 * t626;
	t840 = t547 * t986;
	t798 = t533 * t840 + t546 * t982;
	t532 = t798 * t862;
	t952 = t733 * t736;
	t627 = -t950 + t952;
	t526 = -t532 * t627 - t530;
	t524 = 0.1e1 / t526 ^ 2;
	t529 = t531 * t627;
	t525 = -t532 * t626 + t529;
	t971 = t524 * t525;
	t555 = 0.1e1 / t777;
	t681 = pkin(1) * t690;
	t670 = 0.2e1 * t681;
	t947 = t736 * t774;
	t890 = -0.4e1 * t947;
	t817 = t734 * t890;
	t625 = t670 + t817;
	t878 = pkin(2) * t952;
	t824 = pkin(3) * t878;
	t827 = -0.8e1 * t873;
	t906 = -0.4e1 * t689;
	t970 = (t625 * t887 + (t681 - t824) * t827 + 0.4e1 * t619 * t824 + (pkin(2) * t600 * t906 - 0.8e1 * t672 * t736 + (t871 * t1004 + t653 * t999) * pkin(1)) * t734) * t555;
	t886 = 0.2e1 * t942;
	t948 = t735 * t770;
	t969 = t555 * ((t600 * t977 + 0.2e1 * t622 * t948) * t733 + (t619 * t689 + t666 * t886) * t690);
	t939 = t785 * t700;
	t882 = pkin(1) * t939;
	t820 = -0.12e2 * t882;
	t967 = t638 * t733 * t820;
	t956 = t702 * t782;
	t954 = t703 * t781;
	t938 = -2 * t989;
	t668 = pkin(1) * t907;
	t697 = t733 ^ 2;
	t835 = 0.32e2 / 0.3e1 * t768;
	t804 = t699 * t835;
	t805 = 0.64e2 / 0.3e1 * t645 * t785;
	t813 = t699 * t857;
	t816 = t945 * t979;
	t821 = -0.48e2 * t882;
	t823 = -0.16e2 * pkin(1) * t646 * t770;
	t828 = -0.4e1 * t658 * t972;
	t880 = pkin(1) * t942;
	t830 = 0.4e1 * t880;
	t831 = -0.2e1 * t880;
	t832 = -0.4e1 * t880;
	t856 = t733 * t948;
	t858 = t684 * t959;
	t867 = -0.64e2 * t957;
	t869 = -0.32e2 * t960;
	t870 = pkin(3) * t946;
	t872 = pkin(3) * t957;
	t879 = t638 * t690;
	t884 = pkin(1) * t943;
	t892 = 0.3e1 * t955;
	t895 = 0.8e1 * t957;
	t896 = -0.8e1 * t957;
	t898 = -0.2e1 * t959;
	t899 = -0.4e1 * t959;
	t900 = 0.2e1 * t969;
	t905 = 0.2e1 * t689;
	t909 = 0.6e1 * t978;
	t911 = -0.6e1 * t978;
	t527 = ((0.12e2 * t697 * t700 * t979 + t630 * t899 + t677 * t831 - 0.4e1 * t884) * t866 + 0.8e1 * t696 * t884 + 0.12e2 * t598 * t959 + 0.6e1 * t593 * t880 + (-t589 * t866 / 0.2e1 + t1003) * t689) * t777 + t540 * t900 + (((t632 * t905 + t830 + t899) * t895 + (-t593 * t689 + t635 * t832 - 0.4e1 * t813) * t909) * t777 + t867 * t967) * t667 + (0.24e2 * (-pkin(1) * t698 * t835 - t699 * t805 + t700 * t823 + t735 * t828) * t955 - 0.64e2 * t698 * t816 - 0.96e2 * t658 * t858 - 0.48e2 * t588 * t880 + t973 * t993 + ((-t735 * t870 + t831 + t898) * t867 - 0.48e2 * (-t588 * t689 + t658 * t832 - 0.4e1 * t858) * t978) * t667) * t879 + (-pkin(3) * t536 + (0.2e1 * (-0.2e1 * t682 * t735 - t665 + t821 + t869) * t958 + 0.4e1 * t561 * t872 + (-t609 * t735 + t615 * t938) * t819 + (-0.28e2 * t565 * t948 - 0.6e1 * t569 * t989 + t616 * t821 + t685 * t869) * t892 + pkin(3) * t550 * t978 + t667 * (-t574 * t989 - 0.4e1 * t591 * t948 - 0.4e1 * t696 * t882) * t911 + t651 * t995 + t594 * t820 - 0.6e1 * t563 * t948 - pkin(1) * t974) * t996 + ((-0.8e1 * t580 * t948 + t690 * t804) * t866 - 0.96e2 * t684 * t960 * t690 + t596 * t821 + 0.12e2 * t564 * t948 + (t552 * t911 + t567 * t896 + (-0.24e2 * t668 + 0.64e2 * t856) * t958 + (0.48e2 * t592 * t955 + 0.6e1 * t562) * pkin(1)) * pkin(3) + ((t639 * t993 + t648 * t938) * t895 + (0.24e2 * pkin(1) * t700 * t876 - 0.4e1 * t582 * t989 + 0.8e1 * t597 * t948) * t909) * t667) * t777) * t733;
	t865 = 0.12e2 * t955;
	t891 = t734 * t1004;
	t903 = -0.2e1 * t972;
	t908 = -0.4e1 * pkin(3) * t637;
	t910 = 0.4e1 * t978;
	t913 = t572 * t999;
	t534 = (t770 * t697 * t896 + (t689 * t692 + t898) * t894 + 0.4e1 * t813 + t676 * t830 + t628 * t689) * t777 + t549 * t900 + t865 * t967 + ((-0.8e1 / 0.3e1 * t959 + t832 + t735 * t903) * t865 - 0.24e2 * t858 - 0.24e2 * t637 * t880 - 0.6e1 * t587 * t689) * t879 + (0.24e2 * (-t638 * t696 - t690 * t777) * t882 + ((0.16e2 * t639 * t955 - 0.2e1 * t605) * t777 + t638 * (-0.144e3 * t604 * t955 - 0.24e2 * t585) * t770) * t735 + (t736 * t777 * t913 - t545 + t638 * (t560 * t909 + t621 * t895) + ((pkin(2) * t703 * t891 + 0.4e1 * t603) * t777 + t638 * (-0.48e2 * t612 * t955 - 0.6e1 * t571)) * pkin(1)) * pkin(3)) * t733 + ((t872 * t992 + ((-pkin(3) * t676 + 0.4e1 * t661 * t770) * t735 + (-t647 * t975 - t942) * t1009) * t910) * t777 + t638 * ((t668 - 0.8e1 * t856) * t896 + ((-0.8e1 * t673 * t733 + t690 * t908) * t735 + (-0.4e1 * t610 * t975 - 0.8e1 * t815) * pkin(1)) * t911)) * t667;
	t543 = -t968 + t611 * t900 + pkin(3) * t697 * t997 + t575 * t735 + (-t618 + t827) * t661;
	t914 = -2 * pkin(1) * t774;
	t544 = t933 * t777 + t613 * t900 + (t601 * t666 + 0.4e1 * t874) * t733 + (t914 * t949 + (t652 * t735 + 0.4e1 * t666 * t976) * pkin(2)) * t734;
	t606 = t613 * pkin(3);
	t794 = t797 * t1006;
	t839 = t546 * t985;
	t937 = (-t606 * t794 + (t527 * t841 + t543 * t983 - t547 * t969 / 0.2e1 + (t534 * t839 + t544 * t986) * t533) * t583) * t944 - t532;
	t875 = pkin(2) * t949;
	t826 = pkin(1) * t875;
	t629 = -0.4e1 * t770 * t733 * t826;
	t643 = t826 * t1008;
	t837 = t970 / 0.2e1;
	t859 = t684 * t939;
	t860 = t667 * t954;
	t881 = pkin(1) * t948;
	t528 = t540 * t837 + (-0.32e2 * t629 * t638 + 0.8e1 * t643 * t777) * t861 + (0.24e2 * (0.4e1 * t617 * t956 * t975 + t551 * t947 - t567 * t860) * t777 + (-0.6e1 * t548 * t947 + 0.12e2 * t561 * t860 - 0.8e1 * t579 * t956) * t996) * t734 + ((t536 + (t550 * t996 - 0.6e1 * t552 * t777) * t667) * t734 + (((t631 * t885 + t633 * t855 + t806) * t866 + t634 * t822 + t590 * t888 - 0.6e1 * t581 * t901 + (t1002 * t895 + (t636 * t855 - 0.8e1 * t644 * t942 + t581 - 0.8e1 * t883) * t909) * t667 + t801) * t777 + ((-pkin(1) * t804 - t700 * t805 + t735 * t823 + t828) * t892 + t816 * t995 - 0.12e2 * t658 * t859 + t881 * t998 - t973 + (-0.4e1 * (-(2 * t870) - 0.4e1 * t939) * t957 - (pkin(3) * t998 - 0.24e2 * t658 * t881 - 0.24e2 * t859) * t978) * t667) * t733 * t996) * t736) * pkin(2);
	t535 = (t643 * t894 + (t573 * t891 + t674 * t700 + ((-t921 + t886) * t894 - t620 + (t678 * t906 + t897) * pkin(1)) * pkin(2)) * t736 + ((t643 + (t678 - 0.2e1 * t942) * t978) * t910 + (-0.24e2 * t954 * t975 + t913) * t734) * t667) * t777 + t549 * t837 + t545 * t690 + 0.6e1 * t638 * (0.4e1 * t621 * t734 * t860 + (t629 + (-0.8e1 / 0.3e1 * t939 + t903) * t878) * t893 + t559 * t817 + t684 * t736 * t700 * t829 + t637 * t629 - t587 * t824 + (-(-0.8e1 * t881 + t908) * t733 * t955 + t560 * t690) * t667);
	t541 = t611 * t837 + t625 * t733 * t905 + ((-t733 * t777 + t578) * t736 + (t735 * t777 + (t666 * t1009 + t618) * t733 + (-pkin(3) + 0.2e1 * t976 + t980) * t912) * t734) * pkin(2);
	t542 = (t661 + t875) * t777 + t613 * t837 - 0.2e1 * t625 * t976 - (t670 - 0.4e1 * t824) * t650 + t703 * t733 * t914 + t652 * t878 + (pkin(3) * t890 + (-t601 * t735 + 0.2e1 * t981) * pkin(2)) * t734;
	t602 = -pkin(2) * pkin(3) * t627 + t681;
	t936 = (t602 * t794 + (t528 * t841 + t541 * t983 - t547 * t970 / 0.8e1 + (t535 * t839 + t542 * t986) * t533) * t583) * t944 + t532;
	t838 = t547 * t985;
	t795 = t798 * t1006;
	t523 = 0.1e1 / t526;
	t522 = 0.1e1 / (t524 * t525 ^ 2 + 0.1e1);
	t516 = (t602 * t795 + (t542 * t982 + t546 * t970 / 0.8e1 + t528 * t840 + (t535 * t838 + t541 * t986) * t533) * t583) * t944;
	t514 = (-t606 * t795 + (t544 * t982 + t546 * t969 / 0.2e1 + t527 * t840 + (t534 * t838 + t543 * t986) * t533) * t583) * t944;
	t1 = [0, 0; 0, 0; 0.1e1 + ((t936 * t627 + (-t516 + t531) * t626) * t523 - (-t516 * t627 - t626 * t936 + t529) * t971) * t522, ((-t514 * t626 + t627 * t937 - t530) * t523 - ((-t514 - t531) * t627 - t937 * t626) * t971) * t522;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 05:51:21
	% EndTime: 2020-04-27 05:51:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 1;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-27 05:51:23
	% EndTime: 2020-04-27 05:51:24
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (5310->105), mult. (11444->181), div. (114->6), fcn. (2828->6), ass. (0->81)
	t206 = -4 * pkin(2);
	t158 = pkin(2) ^ 2;
	t152 = cos(qJ(1));
	t197 = pkin(2) * t152;
	t180 = pkin(1) * t197;
	t144 = -0.2e1 * t180;
	t159 = pkin(1) ^ 2;
	t183 = t144 + t159;
	t148 = t152 ^ 2;
	t203 = 0.2e1 * t148;
	t134 = t158 * t203 - t158 + t183;
	t205 = -0.2e1 * t134;
	t142 = pkin(1) - t197;
	t204 = 0.2e1 * t142;
	t202 = (-pkin(4) - pkin(5));
	t201 = (-pkin(4) + pkin(5));
	t149 = sin(qJ(2));
	t150 = sin(qJ(1));
	t188 = t149 * t150;
	t141 = pkin(2) * t188;
	t169 = pkin(3) * t141;
	t140 = -0.2e1 * t169;
	t157 = pkin(3) ^ 2;
	t151 = cos(qJ(2));
	t136 = t142 * t151;
	t175 = pkin(3) * t136;
	t125 = t140 + t157 + t158 + 0.2e1 * t175 + t183;
	t123 = 0.1e1 / t125;
	t200 = t123 / 0.2e1;
	t199 = pkin(1) * (t141 - pkin(3));
	t198 = pkin(2) * t150;
	t196 = pkin(3) * t142;
	t147 = t151 ^ 2;
	t195 = pkin(3) * t147;
	t194 = pkin(3) * t151;
	t124 = 0.1e1 / t125 ^ 2;
	t146 = pkin(1) * t198;
	t186 = t150 * t151;
	t187 = t149 * t152;
	t193 = t124 * (t146 + (t186 - t187) * pkin(3) * pkin(2));
	t131 = pkin(2) * t186 + t149 * t142;
	t192 = t124 * t131 * pkin(3);
	t155 = pkin(4) ^ 2;
	t182 = -pkin(5) ^ 2 + t159;
	t174 = t158 + t182;
	t166 = t157 + t174;
	t138 = -t155 + t166;
	t133 = t144 + t138;
	t126 = t140 + t133;
	t145 = (t157 - t159) * t158;
	t165 = t147 * t205 - t182;
	t160 = sqrt(0.4e1 * t145 * t148 + 0.4e1 * t138 * t180 - (t159 + ((pkin(2) - t201) * (pkin(2) + t201))) * (t159 + ((pkin(2) - t202) * (pkin(2) + t202))) + 0.4e1 * (-t126 * t136 + t133 * t141) * pkin(3) + (0.2e1 * t155 - (6 * t158) + 0.2e1 * t165 - t157) * t157);
	t191 = t131 * t160;
	t190 = t134 * t151;
	t189 = t147 * t157;
	t185 = t151 * t152;
	t184 = -t141 + t136;
	t181 = -0.2e1 * pkin(1) * t158;
	t118 = 0.1e1 / t160;
	t179 = 0.2e1 * t118 * ((t126 * t196 + 0.2e1 * t157 * t190) * t149 + (t133 * t194 + t189 * t204) * t198);
	t178 = -0.4e1 * t152 * t158;
	t177 = pkin(2) * t187;
	t176 = pkin(3) * t190;
	t143 = 0.2e1 * t146;
	t135 = t150 * t178 + t143;
	t168 = pkin(3) * t177;
	t170 = -0.8e1 * t175;
	t173 = (-0.4e1 * t135 * t189 + (t146 - t168) * t170 + 0.4e1 * t133 * t168 + (t126 * t194 * t206 - 0.8e1 * t145 * t152 + (0.8e1 * pkin(3) * t158 * t188 + t138 * t206) * pkin(1)) * t150) * t118 / 0.2e1;
	t137 = 0.3e1 * t157 + t155 + t174;
	t127 = t137 + t144 - 0.4e1 * t169;
	t114 = -t127 * t136 + t191 + (t137 * t188 - 0.2e1 * t152 * t199) * pkin(2) + (-t155 - t157 + (t203 - 0.3e1) * t158 + t165) * pkin(3);
	t113 = 0.1e1 / t114 ^ 2;
	t132 = t144 + t155 + t166;
	t121 = t132 * t142 + 0.2e1 * t176;
	t122 = t132 * t151 + (0.4e1 * t147 - 0.2e1) * t196;
	t130 = pkin(3) + t184;
	t115 = t121 * t149 + t122 * t198 + t130 * t160;
	t167 = 0.1e1 / (t113 * t115 ^ 2 + 0.1e1) * t125;
	t163 = 0.1e1 / t114 * t167;
	t162 = t113 * t115 * t167;
	t1 = [0, 0; 0, 0; 0.2e1 * ((0.2e1 * t135 * t149 * t194 + t130 * t173) * t200 - t115 * t193 + ((-t149 * t160 + t122) * t152 * t200 + (t151 * t160 / 0.2e1 + (pkin(1) * t204 + t132) * t149 / 0.2e1 + (pkin(1) * t151 - pkin(3) + 0.2e1 * t195) * t198) * t123 * t150) * pkin(2)) * t163 - 0.2e1 * (((pkin(2) * t185 + t141) * t160 + t131 * t173 - 0.2e1 * t135 * t195 - (t143 - 0.4e1 * t168) * t136 + t148 * t149 * t181 + t137 * t177 + (pkin(3) * t178 + (-t127 * t151 + 0.2e1 * t199) * pkin(2)) * t150) * t200 - t114 * t193) * t162, 0.1e1 + 0.2e1 * ((-t191 + t130 * t179 + t121 * t151 + ((-t132 + t170) * t198 + pkin(3) * t149 * t205) * t149) * t200 + t115 * t192) * t163 - 0.2e1 * ((t184 * t160 + t131 * t179 + (t142 * t127 + 0.4e1 * t176) * t149 + (t181 * t185 + (t137 * t151 + 0.4e1 * t142 * t195) * pkin(2)) * t150) * t200 + t114 * t192) * t162;];
	Ja_rot = t1;
end