% Calculate kinetic energy for
% palh1m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1DE2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1DE2_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:23
% EndTime: 2020-04-14 20:03:53
% DurationCPUTime: 82.97s
% Computational Cost: add. (1825719->677), mult. (2802284->1115), div. (120872->33), fcn. (1755854->48), ass. (0->466)
t661 = sin(qJ(1));
t667 = cos(qJ(1));
t870 = -t661 * V_base(4) + V_base(5) * t667;
t869 = -2 * pkin(1);
t868 = (-pkin(2) - pkin(13));
t867 = (-pkin(2) + pkin(13));
t866 = -pkin(8) - pkin(3);
t865 = -pkin(8) + pkin(3);
t864 = -pkin(9) - pkin(11);
t863 = pkin(11) - pkin(9);
t681 = pkin(4) ^ 2;
t680 = pkin(5) ^ 2;
t678 = pkin(7) ^ 2;
t686 = pkin(1) ^ 2;
t660 = sin(qJ(2));
t662 = sin(pkin(19));
t666 = cos(qJ(2));
t668 = cos(pkin(19));
t617 = t660 * t668 - t662 * t666;
t840 = pkin(7) * t617;
t794 = t840 * t869 + t686;
t596 = t678 + t794;
t790 = pkin(3) ^ 2 - pkin(8) ^ 2;
t581 = t596 + t790;
t602 = pkin(1) - t840;
t618 = t660 * t662 + t666 * t668;
t566 = (pkin(7) - t866) * (pkin(7) + t866) + t794;
t567 = (pkin(7) - t865) * (pkin(7) + t865) + t794;
t689 = sqrt(-t567 * t566);
t543 = pkin(7) * t581 * t618 + t602 * t689;
t665 = cos(qJ(3));
t798 = t665 * t543;
t806 = t618 * t689;
t542 = -pkin(7) * t806 + t581 * t602;
t659 = sin(qJ(3));
t803 = t659 * t542;
t714 = t803 / 0.2e1 + t798 / 0.2e1;
t592 = 0.1e1 / t596;
t683 = 0.1e1 / pkin(3);
t808 = t592 * t683;
t532 = t714 * t808;
t799 = t665 * t542;
t802 = t659 * t543;
t713 = -t799 / 0.2e1 + t802 / 0.2e1;
t533 = t713 * t808;
t650 = pkin(23) + pkin(22);
t641 = sin(t650);
t642 = cos(t650);
t500 = -t532 * t642 + t533 * t641;
t846 = pkin(5) * t500;
t797 = -0.2e1 * pkin(4) * t846 + t680;
t495 = t681 + t797;
t493 = 0.1e1 / t495;
t862 = t493 / 0.2e1;
t606 = t618 * qJD(2);
t605 = t617 * qJD(2);
t782 = pkin(1) * pkin(7) * t606;
t813 = 0.2e1 * (t566 + t567) * t782 / t689;
t768 = -t813 / 0.2e1;
t711 = t605 * t689 + t618 * t768;
t518 = ((t602 * t869 - t581) * t606 + t711) * pkin(7);
t861 = -t518 / 0.2e1;
t779 = -0.2e1 * t606 * t618;
t807 = t606 * t689;
t519 = t602 * t813 / 0.2e1 + t678 * pkin(1) * t779 + (-t581 * t605 - t807) * pkin(7);
t860 = t519 / 0.2e1;
t654 = sin(pkin(20));
t657 = cos(pkin(20));
t719 = t654 * t665 + t657 * t659;
t842 = pkin(6) * t719;
t780 = pkin(1) * t842;
t599 = 0.2e1 * t780;
t679 = pkin(6) ^ 2;
t791 = t679 + t686;
t585 = t599 + t791;
t582 = 0.1e1 / t585;
t859 = t582 / 0.2e1;
t652 = sin(pkin(23));
t857 = t652 / 0.2e1;
t653 = sin(pkin(21));
t856 = t653 / 0.2e1;
t855 = -t665 / 0.2e1;
t669 = cos(pkin(18));
t854 = t669 / 0.2e1;
t685 = 0.1e1 / pkin(2);
t853 = t685 / 0.2e1;
t614 = t654 * t659 - t657 * t665;
t601 = t614 * qJD(3);
t852 = pkin(1) * t601;
t851 = pkin(1) * t666;
t684 = pkin(2) ^ 2;
t775 = -pkin(13) ^ 2 + t791;
t579 = t599 + t684 + t775;
t598 = -pkin(1) - t842;
t795 = t599 + t679;
t562 = ((pkin(1) - t868) * (pkin(1) + t868)) + t795;
t563 = ((pkin(1) - t867) * (pkin(1) + t867)) + t795;
t812 = t563 * t562;
t688 = sqrt(-t812);
t841 = pkin(6) * t614;
t539 = -t579 * t598 - t688 * t841;
t540 = t579 * t841 - t598 * t688;
t767 = t582 * t853;
t526 = qJ(2) + atan2(t540 * t767, t539 * t767);
t525 = cos(t526);
t850 = pkin(2) * t525;
t763 = 0.1e1 / t596 ^ 2 * t782;
t435 = ((-t799 + t802) * t763 + (qJD(3) * t714 + t518 * t855 + t659 * t860) * t592) * t683;
t436 = ((-t798 - t803) * t763 + (qJD(3) * t713 + t519 * t855 + t659 * t861) * t592) * t683;
t428 = t435 * t641 + t436 * t642;
t849 = pkin(4) * t428;
t655 = cos(pkin(23));
t529 = (-t655 * t542 / 0.2e1 + t543 * t857) * t808;
t530 = (t655 * t543 / 0.2e1 + t542 * t857) * t808;
t506 = qJ(2) + atan2(t530, t529);
t503 = pkin(22) - t506;
t848 = pkin(4) * cos(t503);
t792 = pkin(9) ^ 2 - pkin(11) ^ 2;
t490 = t495 + t792;
t496 = -pkin(4) + t846;
t485 = (pkin(4) - t864) * (pkin(4) + t864) + t797;
t486 = (pkin(4) - t863) * (pkin(4) + t863) + t797;
t687 = sqrt(-t486 * t485);
t501 = t532 * t641 + t533 * t642;
t845 = pkin(5) * t501;
t425 = t490 * t845 - t496 * t687;
t847 = pkin(5) * t425;
t651 = qJ(2) + qJ(3);
t646 = sin(t651);
t844 = pkin(5) * t646;
t843 = pkin(6) * t540;
t839 = pkin(16) * t661;
t524 = sin(t526);
t838 = pkin(2) * t524;
t837 = t660 * pkin(1);
t836 = Icges(2,4) * t661;
t835 = Icges(3,4) * t660;
t834 = Icges(3,4) * t666;
t833 = Icges(4,4) * t646;
t647 = cos(t651);
t832 = Icges(4,4) * t647;
t491 = t495 - t792;
t497 = -pkin(4) * t500 + pkin(5);
t815 = t501 * t687;
t424 = -pkin(4) * t815 + t491 * t497;
t426 = pkin(4) * t491 * t501 + t497 * t687;
t656 = cos(pkin(21));
t673 = 0.1e1 / pkin(11);
t816 = t493 * t673;
t402 = (t424 * t856 + t426 * t656 / 0.2e1) * t816;
t403 = (-t424 * t656 / 0.2e1 + t426 * t856) * t816;
t393 = atan2(t402, t403) + t651;
t391 = sin(t393);
t831 = Icges(5,4) * t391;
t392 = cos(t393);
t830 = Icges(5,4) * t392;
t580 = t596 - t790;
t603 = pkin(1) * t617 - pkin(7);
t541 = -pkin(1) * t806 - t580 * t603;
t544 = pkin(1) * t580 * t618 - t603 * t689;
t663 = sin(pkin(18));
t677 = 0.1e1 / pkin(8);
t809 = t592 * t677;
t534 = (t541 * t854 - t663 * t544 / 0.2e1) * t809;
t535 = (t544 * t854 + t541 * t663 / 0.2e1) * t809;
t512 = atan2(t535, t534);
t508 = sin(t512);
t829 = Icges(7,4) * t508;
t509 = cos(t512);
t828 = Icges(7,4) * t509;
t504 = sin(t506);
t827 = Icges(8,4) * t504;
t505 = cos(t506);
t826 = Icges(8,4) * t505;
t825 = Icges(9,4) * t524;
t824 = Icges(9,4) * t525;
t578 = t684 - t775 - 0.2e1 * t780;
t764 = 0.1e1 / pkin(13) * t853;
t516 = atan2(t688 * t764, t578 * t764) + t526;
t514 = sin(t516);
t823 = Icges(10,4) * t514;
t515 = cos(t516);
t822 = Icges(10,4) * t515;
t821 = t391 * t661;
t820 = t391 * t667;
t781 = pkin(5) * t849;
t819 = 0.2e1 * (t485 + t486) * t781 / t687;
t818 = t428 * t687;
t817 = t493 * t656;
t783 = pkin(6) * t852;
t814 = 0.2e1 * (t562 + t563) * t783 / t688;
t811 = t592 * t655;
t810 = t592 * t663;
t658 = sin(qJ(4));
t805 = t658 * t661;
t804 = t658 * t667;
t664 = cos(qJ(4));
t801 = t661 * t664;
t800 = t664 * t667;
t796 = pkin(4) * sin(t503);
t793 = pkin(5) * t647;
t423 = -pkin(5) * t815 - t490 * t496;
t675 = 0.1e1 / pkin(9);
t770 = t675 * t862;
t400 = -atan2(t425 * t770, t423 * t770) + t503;
t398 = sin(t400);
t789 = Icges(11,4) * t398;
t399 = cos(t400);
t788 = Icges(11,4) * t399;
t787 = qJD(4) * t391;
t528 = 0.1e1 / t529 ^ 2;
t752 = t543 * t763;
t757 = t542 * t763;
t766 = t592 * t857;
t416 = ((t518 * t766 + t652 * t757 + t655 * t752 + t811 * t860) / t529 - (t519 * t766 + t652 * t752 - t655 * t757 + t811 * t861) * t530 * t528) / (t528 * t530 ^ 2 + 0.1e1) * t683;
t786 = -qJD(2) - t416;
t538 = 0.1e1 / t539 ^ 2;
t583 = 0.1e1 / t585 ^ 2;
t600 = t719 * qJD(3);
t769 = -t814 / 0.2e1;
t434 = 0.2e1 * (((t598 * t769 + (t600 * t579 - t601 * t688) * pkin(6)) * t859 + (-t582 * t614 * t679 + t583 * t843) * t852) / t539 - ((-t601 * t579 - t600 * t688 + t614 * t769) * t859 + (t539 * t583 + t582 * t598) * t852) * t538 * t843) * pkin(2) / (t538 * t540 ^ 2 + 0.1e1) * t585 * t685;
t785 = -qJD(2) - t434;
t784 = -qJD(2) - qJD(3);
t778 = V_base(5) * pkin(14) + V_base(1);
t772 = -t819 / 0.2e1;
t771 = t493 * t856;
t765 = t592 * t854;
t638 = qJD(2) * t661 + V_base(4);
t643 = V_base(6) + qJD(1);
t494 = 0.1e1 / t495 ^ 2;
t762 = t494 * t781;
t607 = t837 * t661;
t761 = t607 - t839;
t637 = -qJD(2) * t667 + V_base(5);
t760 = t637 * t851 + t778;
t414 = t416 * t661 + t638;
t432 = t434 * t661 + t638;
t616 = qJD(3) * t661 + t638;
t756 = t544 * t763;
t755 = t424 * t762;
t754 = t426 * t762;
t753 = t541 * t763;
t564 = t793 * t661;
t751 = -t564 + t761;
t615 = t667 * t784 + V_base(5);
t750 = t615 * t844 + t760;
t749 = pkin(10) * t392 + pkin(12) * t391;
t748 = -rSges(3,1) * t660 - rSges(3,2) * t666;
t747 = rSges(4,1) * t647 - rSges(4,2) * t646;
t746 = rSges(5,1) * t392 - rSges(5,2) * t391;
t745 = rSges(7,1) * t509 - rSges(7,2) * t508;
t744 = -rSges(8,1) * t504 - rSges(8,2) * t505;
t743 = -rSges(9,1) * t524 - rSges(9,2) * t525;
t742 = rSges(10,1) * t514 + rSges(10,2) * t515;
t427 = t435 * t642 - t436 * t641;
t712 = -t427 * t687 + t501 * t772;
t372 = ((-0.2e1 * pkin(5) * t497 - t491) * t428 + t712) * pkin(4);
t373 = t497 * t819 / 0.2e1 - 0.2e1 * t681 * t428 * t845 + (t427 * t491 - t818) * pkin(4);
t401 = 0.1e1 / t403 ^ 2;
t339 = ((t372 * t771 + t653 * t755 + t373 * t817 / 0.2e1 + t656 * t754) / t403 - (-t372 * t817 / 0.2e1 - t656 * t755 + t373 * t771 + t653 * t754) * t402 * t401) / (t401 * t402 ^ 2 + 0.1e1) * t673;
t337 = t339 * t661 + t616;
t741 = -rSges(11,1) * t398 + rSges(11,2) * t399;
t740 = -Icges(3,1) * t660 - t834;
t739 = Icges(4,1) * t647 - t833;
t738 = Icges(5,1) * t392 - t831;
t737 = Icges(7,1) * t509 - t829;
t736 = -Icges(8,1) * t504 - t826;
t735 = -Icges(9,1) * t524 - t824;
t734 = Icges(10,1) * t514 + t822;
t733 = -Icges(3,2) * t666 - t835;
t732 = -Icges(4,2) * t646 + t832;
t731 = -Icges(5,2) * t391 + t830;
t730 = -Icges(7,2) * t508 + t828;
t729 = -Icges(8,2) * t505 - t827;
t728 = -Icges(9,2) * t525 - t825;
t727 = Icges(10,2) * t515 + t823;
t726 = -Icges(3,5) * t660 - Icges(3,6) * t666;
t725 = Icges(4,5) * t647 - Icges(4,6) * t646;
t724 = Icges(5,5) * t392 - Icges(5,6) * t391;
t723 = Icges(7,5) * t509 - Icges(7,6) * t508;
t722 = -Icges(8,5) * t504 - Icges(8,6) * t505;
t721 = -Icges(9,5) * t524 - Icges(9,6) * t525;
t720 = Icges(10,5) * t514 + Icges(10,6) * t515;
t718 = -Icges(11,1) * t398 + t788;
t717 = Icges(11,2) * t399 - t789;
t716 = -Icges(11,5) * t398 + Icges(11,6) * t399;
t715 = t643 * t667 * pkin(16) - V_base(4) * pkin(14) + V_base(2);
t710 = -t870 * pkin(16) + V_base(3);
t336 = V_base(5) + (-t339 + t784) * t667;
t709 = (-Icges(5,3) * t667 + t661 * t724) * t336 + (Icges(5,3) * t661 + t667 * t724) * t337 + (Icges(5,5) * t391 + Icges(5,6) * t392) * t643;
t422 = 0.1e1 / t423 ^ 2;
t344 = 0.2e1 * (((t496 * t772 + (t427 * t490 - t818) * pkin(5)) * t862 + (-t493 * t501 * t680 + t494 * t847) * t849) / t423 - ((-t428 * t490 + t712) * t862 + (t423 * t494 + t493 * t496) * t849) * t422 * t847) * pkin(9) / (t422 * t425 ^ 2 + 0.1e1) * t495 * t675;
t342 = V_base(5) + (-t344 + t786) * t667;
t343 = t344 * t661 + t414;
t708 = (-Icges(11,3) * t667 + t661 * t716) * t342 + (Icges(11,3) * t661 + t667 * t716) * t343 + (-Icges(11,5) * t399 - Icges(11,6) * t398) * t643;
t413 = t667 * t786 + V_base(5);
t707 = (-Icges(8,3) * t667 + t661 * t722) * t413 + (Icges(8,3) * t661 + t667 * t722) * t414 + (Icges(8,5) * t505 - Icges(8,6) * t504) * t643;
t517 = ((0.2e1 * pkin(7) * t603 - t580) * t606 + t711) * pkin(1);
t520 = t603 * t768 + t686 * pkin(7) * t779 + (-t580 * t605 - t807) * pkin(1);
t531 = 0.1e1 / t534 ^ 2;
t419 = ((t520 * t765 + t669 * t756 + t517 * t810 / 0.2e1 + t663 * t753) / t534 - (t517 * t765 + t669 * t753 - t520 * t810 / 0.2e1 - t663 * t756) * t535 * t531) / (t531 * t535 ^ 2 + 0.1e1) * t677;
t417 = -t419 * t667 + V_base(5);
t418 = t419 * t661 + V_base(4);
t706 = (-Icges(7,3) * t667 + t661 * t723) * t417 + (Icges(7,3) * t661 + t667 * t723) * t418 + (Icges(7,5) * t508 + Icges(7,6) * t509) * t643;
t574 = 0.1e1 / t578 ^ 2;
t513 = (0.1e1 / t578 * t814 / 0.2e1 - 0.2e1 * t688 * t574 * t783) / (-t574 * t812 + 0.1e1);
t429 = V_base(5) + (-t513 + t785) * t667;
t430 = t513 * t661 + t432;
t705 = (-Icges(10,3) * t667 + t661 * t720) * t429 + (Icges(10,3) * t661 + t667 * t720) * t430 + (-Icges(10,5) * t515 + Icges(10,6) * t514) * t643;
t431 = t667 * t785 + V_base(5);
t704 = (-Icges(9,3) * t667 + t661 * t721) * t431 + (Icges(9,3) * t661 + t667 * t721) * t432 + (Icges(9,5) * t525 - Icges(9,6) * t524) * t643;
t703 = (-Icges(4,3) * t667 + t661 * t725) * t615 + (Icges(4,3) * t661 + t667 * t725) * t616 + (Icges(4,5) * t646 + Icges(4,6) * t647) * t643;
t702 = (-Icges(3,3) * t667 + t661 * t726) * t637 + (Icges(3,3) * t661 + t667 * t726) * t638 + (Icges(3,5) * t666 - Icges(3,6) * t660) * t643;
t608 = t837 * t667;
t701 = -t643 * t608 - t638 * t851 + t715;
t700 = -t638 * t607 + t608 * t637 + t710;
t565 = t793 * t667;
t699 = t616 * t564 - t565 * t615 + t700;
t698 = t643 * t565 - t616 * t844 + t701;
t359 = -Icges(5,6) * t667 + t661 * t731;
t360 = Icges(5,6) * t661 + t667 * t731;
t361 = -Icges(5,5) * t667 + t661 * t738;
t362 = Icges(5,5) * t661 + t667 * t738;
t368 = Icges(5,2) * t392 + t831;
t369 = Icges(5,1) * t391 + t830;
t697 = (-t360 * t391 + t362 * t392) * t337 + (-t359 * t391 + t361 * t392) * t336 + (-t368 * t391 + t369 * t392) * t643;
t376 = -Icges(11,6) * t667 + t661 * t717;
t377 = Icges(11,6) * t661 + t667 * t717;
t378 = -Icges(11,5) * t667 + t661 * t718;
t379 = Icges(11,5) * t661 + t667 * t718;
t383 = -Icges(11,2) * t398 - t788;
t384 = -Icges(11,1) * t399 - t789;
t696 = (t377 * t399 - t379 * t398) * t343 + (t376 * t399 - t378 * t398) * t342 + (t383 * t399 - t384 * t398) * t643;
t441 = -Icges(8,6) * t667 + t661 * t729;
t442 = Icges(8,6) * t661 + t667 * t729;
t443 = -Icges(8,5) * t667 + t661 * t736;
t444 = Icges(8,5) * t661 + t667 * t736;
t456 = -Icges(8,2) * t504 + t826;
t457 = Icges(8,1) * t505 - t827;
t695 = (-t442 * t505 - t444 * t504) * t414 + (-t441 * t505 - t443 * t504) * t413 + (-t456 * t505 - t457 * t504) * t643;
t449 = -Icges(7,6) * t667 + t661 * t730;
t450 = Icges(7,6) * t661 + t667 * t730;
t451 = -Icges(7,5) * t667 + t661 * t737;
t452 = Icges(7,5) * t661 + t667 * t737;
t460 = Icges(7,2) * t509 + t829;
t461 = Icges(7,1) * t508 + t828;
t694 = (-t450 * t508 + t452 * t509) * t418 + (-t449 * t508 + t451 * t509) * t417 + (-t460 * t508 + t461 * t509) * t643;
t465 = -Icges(10,6) * t667 + t661 * t727;
t466 = Icges(10,6) * t661 + t667 * t727;
t467 = -Icges(10,5) * t667 + t661 * t734;
t468 = Icges(10,5) * t661 + t667 * t734;
t472 = Icges(10,2) * t514 - t822;
t473 = -Icges(10,1) * t515 + t823;
t693 = (t466 * t515 + t468 * t514) * t430 + (t465 * t515 + t467 * t514) * t429 + (t472 * t515 + t473 * t514) * t643;
t477 = -Icges(9,6) * t667 + t661 * t728;
t478 = Icges(9,6) * t661 + t667 * t728;
t479 = -Icges(9,5) * t667 + t661 * t735;
t480 = Icges(9,5) * t661 + t667 * t735;
t488 = -Icges(9,2) * t524 + t824;
t489 = Icges(9,1) * t525 - t825;
t692 = (-t478 * t525 - t480 * t524) * t432 + (-t477 * t525 - t479 * t524) * t431 + (-t488 * t525 - t489 * t524) * t643;
t570 = -Icges(4,6) * t667 + t661 * t732;
t571 = Icges(4,6) * t661 + t667 * t732;
t572 = -Icges(4,5) * t667 + t661 * t739;
t573 = Icges(4,5) * t661 + t667 * t739;
t610 = Icges(4,2) * t647 + t833;
t611 = Icges(4,1) * t646 + t832;
t691 = (-t571 * t646 + t573 * t647) * t616 + (-t570 * t646 + t572 * t647) * t615 + (-t610 * t646 + t611 * t647) * t643;
t588 = -Icges(3,6) * t667 + t661 * t733;
t589 = Icges(3,6) * t661 + t667 * t733;
t590 = -Icges(3,5) * t667 + t661 * t740;
t591 = Icges(3,5) * t661 + t667 * t740;
t627 = -Icges(3,2) * t660 + t834;
t630 = Icges(3,1) * t666 - t835;
t690 = (-t589 * t666 - t591 * t660) * t638 + (-t588 * t666 - t590 * t660) * t637 + (-t627 * t666 - t630 * t660) * t643;
t648 = Icges(2,4) * t667;
t635 = rSges(2,1) * t667 - rSges(2,2) * t661;
t634 = rSges(3,1) * t666 - rSges(3,2) * t660;
t633 = rSges(2,1) * t661 + rSges(2,2) * t667;
t632 = Icges(2,1) * t667 - t836;
t631 = Icges(2,1) * t661 + t648;
t629 = -Icges(2,2) * t661 + t648;
t628 = Icges(2,2) * t667 + t836;
t623 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t622 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t621 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t612 = rSges(4,1) * t646 + rSges(4,2) * t647;
t595 = rSges(3,3) * t661 + t667 * t748;
t594 = -rSges(3,3) * t667 + t661 * t748;
t577 = rSges(4,3) * t661 + t667 * t747;
t576 = -rSges(4,3) * t667 + t661 * t747;
t561 = V_base(5) * rSges(2,3) - t633 * t643 + t778;
t560 = t635 * t643 + V_base(2) + (-pkin(14) - rSges(2,3)) * V_base(4);
t559 = t633 * V_base(4) - t635 * V_base(5) + V_base(3);
t552 = t634 * t637 + (-t594 - t839) * t643 + t778;
t551 = t595 * t643 - t634 * t638 + t715;
t550 = t594 * t638 - t595 * t637 + t710;
t546 = t612 * t615 + (-t576 + t761) * t643 + t760;
t545 = t577 * t643 - t612 * t616 + t701;
t537 = t576 * t616 - t577 * t615 + t700;
t522 = t838 * t667;
t521 = t838 * t661;
t492 = rSges(9,1) * t525 - rSges(9,2) * t524;
t484 = t796 * t667;
t483 = t796 * t661;
t482 = rSges(9,3) * t661 + t667 * t743;
t481 = -rSges(9,3) * t667 + t661 * t743;
t474 = -rSges(10,1) * t515 + rSges(10,2) * t514;
t470 = rSges(10,3) * t661 + t667 * t742;
t469 = -rSges(10,3) * t667 + t661 * t742;
t462 = rSges(7,1) * t508 + rSges(7,2) * t509;
t458 = rSges(8,1) * t505 - rSges(8,2) * t504;
t454 = rSges(7,3) * t661 + t667 * t745;
t453 = -rSges(7,3) * t667 + t661 * t745;
t446 = rSges(8,3) * t661 + t667 * t744;
t445 = -rSges(8,3) * t667 + t661 * t744;
t421 = t431 * t492 + (-t481 - t839) * t643 + t778;
t420 = -t432 * t492 + t482 * t643 + t715;
t412 = -t431 * t482 + t432 * t481 + t710;
t410 = -V_base(5) * pkin(17) + t417 * t462 + (pkin(15) * t661 - t453) * t643 + t778;
t409 = -t418 * t462 + V_base(2) + (-pkin(14) + pkin(17)) * V_base(4) + (-pkin(15) * t667 + t454) * t643;
t408 = t413 * t458 + (-t445 + t761) * t643 + t760;
t407 = -t414 * t458 + t446 * t643 + t701;
t405 = t431 * t850 + t429 * t474 + (-t469 + t521 - t839) * t643 + t778;
t404 = -t432 * t850 - t430 * t474 + (t470 - t522) * t643 + t715;
t397 = t870 * pkin(15) - t417 * t454 + t418 * t453 + V_base(3);
t396 = -t413 * t446 + t414 * t445 + t700;
t395 = -t429 * t470 + t430 * t469 + t431 * t522 - t432 * t521 + t710;
t390 = -qJD(4) * t392 + t643;
t389 = t392 * t800 + t805;
t388 = -t392 * t804 + t801;
t387 = t392 * t801 - t804;
t386 = -t392 * t805 - t800;
t385 = -rSges(11,1) * t399 - rSges(11,2) * t398;
t381 = rSges(11,3) * t661 + t667 * t741;
t380 = -rSges(11,3) * t667 + t661 * t741;
t371 = pkin(10) * t391 - pkin(12) * t392;
t370 = rSges(5,1) * t391 + rSges(5,2) * t392;
t366 = t749 * t667;
t365 = t749 * t661;
t364 = rSges(5,3) * t661 + t667 * t746;
t363 = -rSges(5,3) * t667 + t661 * t746;
t356 = -rSges(6,3) * t392 + (rSges(6,1) * t664 - rSges(6,2) * t658) * t391;
t355 = -Icges(6,5) * t392 + (Icges(6,1) * t664 - Icges(6,4) * t658) * t391;
t354 = -Icges(6,6) * t392 + (Icges(6,4) * t664 - Icges(6,2) * t658) * t391;
t353 = -Icges(6,3) * t392 + (Icges(6,5) * t664 - Icges(6,6) * t658) * t391;
t352 = rSges(6,1) * t389 + rSges(6,2) * t388 + rSges(6,3) * t820;
t351 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t821;
t350 = Icges(6,1) * t389 + Icges(6,4) * t388 + Icges(6,5) * t820;
t349 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t821;
t348 = Icges(6,4) * t389 + Icges(6,2) * t388 + Icges(6,6) * t820;
t347 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t821;
t346 = Icges(6,5) * t389 + Icges(6,6) * t388 + Icges(6,3) * t820;
t345 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t821;
t341 = t413 * t848 + t342 * t385 + (-t380 - t483 + t761) * t643 + t760;
t340 = -t414 * t848 - t343 * t385 + (t381 + t484) * t643 + t701;
t335 = t667 * t787 + t337;
t334 = t661 * t787 + t336;
t333 = -t342 * t381 + t343 * t380 - t413 * t484 + t414 * t483 + t700;
t332 = t336 * t370 + (-t363 + t751) * t643 + t750;
t331 = -t337 * t370 + t364 * t643 + t698;
t330 = -t336 * t364 + t337 * t363 + t699;
t329 = t334 * t356 + t336 * t371 - t351 * t390 + (-t365 + t751) * t643 + t750;
t328 = -t335 * t356 - t337 * t371 + t352 * t390 + t366 * t643 + t698;
t327 = -t334 * t352 + t335 * t351 - t336 * t366 + t337 * t365 + t699;
t1 = t334 * ((t346 * t821 + t348 * t386 + t350 * t387) * t335 + (t345 * t821 + t347 * t386 + t349 * t387) * t334 + (t353 * t821 + t354 * t386 + t355 * t387) * t390) / 0.2e1 + t335 * ((t346 * t820 + t348 * t388 + t350 * t389) * t335 + (t345 * t820 + t347 * t388 + t349 * t389) * t334 + (t353 * t820 + t354 * t388 + t355 * t389) * t390) / 0.2e1 + V_base(5) * t643 * (Icges(2,5) * t661 + Icges(2,6) * t667) + t643 * V_base(4) * (Icges(2,5) * t667 - Icges(2,6) * t661) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((t628 * t667 + t631 * t661 + Icges(1,2)) * V_base(5) + (t629 * t667 + t632 * t661 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t628 * t661 + t631 * t667 + Icges(1,4)) * V_base(5) + (-t629 * t661 + t632 * t667 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t478 * t524 + t480 * t525) * t432 + (-t477 * t524 + t479 * t525) * t431 + (t571 * t647 + t573 * t646) * t616 + (t570 * t647 + t572 * t646) * t615 + (t466 * t514 - t468 * t515) * t430 + (t465 * t514 - t467 * t515) * t429 + (-t589 * t660 + t591 * t666) * t638 + (-t588 * t660 + t590 * t666) * t637 + (t450 * t509 + t452 * t508) * t418 + (t449 * t509 + t451 * t508) * t417 + (-t377 * t398 - t379 * t399) * t343 + (-t376 * t398 - t378 * t399) * t342 + (t360 * t392 + t362 * t391) * t337 + (t359 * t392 + t361 * t391) * t336 + (-t442 * t504 + t444 * t505) * t414 + (-t441 * t504 + t443 * t505) * t413 + (t610 * t647 + t611 * t646 - t627 * t660 + t630 * t666 + t472 * t514 - t473 * t515 + t460 * t509 + t461 * t508 - t383 * t398 - t384 * t399 + t368 * t392 + t369 * t391 - t456 * t504 + t457 * t505 - t488 * t524 + t489 * t525 + Icges(2,3)) * t643) * t643 / 0.2e1 + m(11) * (t333 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(2) * (t559 ^ 2 + t560 ^ 2 + t561 ^ 2) / 0.2e1 + t638 * (t702 * t661 + t667 * t690) / 0.2e1 + t637 * (t690 * t661 - t667 * t702) / 0.2e1 + t616 * (t703 * t661 + t667 * t691) / 0.2e1 + t615 * (t691 * t661 - t667 * t703) / 0.2e1 + t432 * (t704 * t661 + t667 * t692) / 0.2e1 + t431 * (t692 * t661 - t667 * t704) / 0.2e1 + t430 * (t705 * t661 + t667 * t693) / 0.2e1 + t429 * (t693 * t661 - t667 * t705) / 0.2e1 + t418 * (t706 * t661 + t667 * t694) / 0.2e1 + t417 * (t694 * t661 - t667 * t706) / 0.2e1 + t414 * (t707 * t661 + t667 * t695) / 0.2e1 + t413 * (t695 * t661 - t667 * t707) / 0.2e1 + t343 * (t708 * t661 + t667 * t696) / 0.2e1 + t342 * (t696 * t661 - t667 * t708) / 0.2e1 + t337 * (t709 * t661 + t667 * t697) / 0.2e1 + t336 * (t697 * t661 - t667 * t709) / 0.2e1 + m(9) * (t412 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(1) * (t621 ^ 2 + t622 ^ 2 + t623 ^ 2) / 0.2e1 + m(7) * (t397 ^ 2 + t409 ^ 2 + t410 ^ 2) / 0.2e1 + m(3) * (t550 ^ 2 + t551 ^ 2 + t552 ^ 2) / 0.2e1 + m(4) * (t537 ^ 2 + t545 ^ 2 + t546 ^ 2) / 0.2e1 + m(8) * (t396 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + t390 * ((-t334 * t345 - t335 * t346 - t353 * t390) * t392 + ((-t348 * t658 + t350 * t664) * t335 + (-t347 * t658 + t349 * t664) * t334 + (-t354 * t658 + t355 * t664) * t390) * t391) / 0.2e1 + m(10) * (t395 ^ 2 + t404 ^ 2 + t405 ^ 2) / 0.2e1 + m(6) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1;
T = t1;
