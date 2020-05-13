% Calculate kinetic energy for
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1TE_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:36
% EndTime: 2020-04-12 20:12:15
% DurationCPUTime: 96.77s
% Computational Cost: add. (1944303->769), mult. (2983982->1250), div. (129344->33), fcn. (1870952->28), ass. (0->439)
t674 = sin(qJ(1));
t677 = cos(qJ(1));
t838 = -t674 * V_base(4) + V_base(5) * t677;
t684 = pkin(6) ^ 2;
t667 = sin(pkin(20));
t671 = cos(pkin(20));
t673 = sin(qJ(3));
t676 = cos(qJ(3));
t734 = t667 * t676 + t671 * t673;
t800 = pkin(6) * t734;
t833 = -2 * pkin(1);
t781 = -t800 * t833 + t684;
t830 = (-pkin(2) - pkin(13));
t597 = ((pkin(1) - t830) * (pkin(1) + t830)) + t781;
t824 = (pkin(13) - pkin(2));
t598 = ((pkin(1) - t824) * (pkin(1) + t824)) + t781;
t786 = t598 * t597;
t691 = sqrt(-t786);
t816 = t691 / 0.2e1;
t810 = sin(qJ(2));
t811 = sin(pkin(19));
t813 = cos(qJ(2));
t814 = cos(pkin(19));
t719 = t810 * t814 - t811 * t813;
t717 = pkin(7) * t719;
t711 = (-0.2e1 * t717 + pkin(1)) * pkin(1);
t827 = -pkin(8) - pkin(3);
t828 = pkin(7) + pkin(8);
t829 = pkin(7) - pkin(8);
t692 = sqrt(-((-pkin(3) + t828) * (pkin(3) + t829) + t711) * ((pkin(7) - t827) * (pkin(7) + t827) + t711));
t693 = pkin(3) ^ 2;
t720 = t810 * t811 + t813 * t814;
t627 = t720 * qJD(2);
t831 = pkin(1) * pkin(7);
t772 = t627 * t831;
t787 = 0.4e1 * (t828 * t829 - t693 + t711) * t772 / t692;
t837 = -t787 / 0.2e1;
t606 = (pkin(1) ^ 2) + t781;
t678 = pkin(13) ^ 2;
t688 = pkin(2) ^ 2;
t602 = t606 - t678 + t688;
t621 = -pkin(1) - t800;
t630 = t667 * t673 - t671 * t676;
t799 = pkin(6) * t630;
t567 = -t602 * t621 - t691 * t799;
t568 = t602 * t799 - t621 * t691;
t767 = -t813 / 0.2e1;
t604 = 0.1e1 / t606;
t689 = 0.1e1 / pkin(2);
t785 = t604 * t689;
t836 = (t810 * t568 / 0.2e1 + t567 * t767) * t785;
t835 = pkin(7) ^ 2;
t834 = 0.1e1 / pkin(8);
t616 = t711 + t835;
t832 = 0.1e1 / t616;
t826 = -pkin(9) - pkin(11);
t825 = -pkin(9) + pkin(11);
t686 = pkin(4) ^ 2;
t685 = pkin(5) ^ 2;
t779 = pkin(8) ^ 2 - t693;
t603 = t616 - t779;
t625 = -t717 + pkin(1);
t715 = t692 * t720;
t569 = -pkin(7) * t715 + t625 * t603;
t718 = pkin(7) * t720;
t570 = t603 * t718 + t625 * t692;
t687 = 0.1e1 / pkin(3);
t778 = -t832 / 0.2e1;
t754 = t676 * t778;
t554 = (t569 * t673 * t778 + t570 * t754) * t687;
t777 = t832 / 0.2e1;
t555 = (t570 * t673 * t777 + t569 * t754) * t687;
t663 = pkin(23) + pkin(22);
t656 = sin(t663);
t657 = cos(t663);
t529 = t554 * t657 + t555 * t656;
t806 = pkin(5) * t529;
t782 = -0.2e1 * pkin(4) * t806 + t685;
t513 = t686 + t782;
t511 = 0.1e1 / t513;
t823 = t511 / 0.2e1;
t601 = t678 - t684 + t688 + (-pkin(1) - 0.2e1 * t800) * pkin(1);
t822 = -t601 / 0.2e1;
t821 = t604 / 0.2e1;
t666 = sin(pkin(21));
t819 = t666 / 0.2e1;
t669 = cos(pkin(22));
t818 = -t669 / 0.2e1;
t817 = -t691 / 0.2e1;
t815 = cos(pkin(18));
t812 = sin(pkin(18));
t624 = t630 * qJD(3);
t809 = pkin(1) * t624;
t626 = t719 * qJD(2);
t709 = (t718 * t833 - t692) * t627;
t542 = t625 * t787 / 0.2e1 + (-t626 * t603 + t709) * pkin(7);
t613 = 0.1e1 / t616 ^ 2;
t752 = t613 * t772;
t742 = t570 * t752;
t753 = qJD(3) * t777;
t712 = t542 * t777 + t569 * t753 + t742;
t708 = t626 * t692 + t720 * t837;
t541 = ((t625 * t833 - t603) * t627 + t708) * pkin(7);
t743 = t569 * t752;
t713 = t541 * t778 + t570 * t753 - t743;
t470 = (t673 * t712 + t676 * t713) * t687;
t471 = (t673 * t713 - t676 * t712) * t687;
t451 = t470 * t656 + t471 * t657;
t808 = pkin(4) * t451;
t780 = pkin(9) ^ 2 - pkin(11) ^ 2;
t505 = t513 + t780;
t514 = -pkin(4) + t806;
t503 = (pkin(4) - t826) * (pkin(4) + t826) + t782;
t504 = (pkin(4) - t825) * (pkin(4) + t825) + t782;
t690 = sqrt(-t504 * t503);
t530 = -t554 * t656 + t555 * t657;
t805 = pkin(5) * t530;
t448 = t505 * t805 - t514 * t690;
t807 = pkin(5) * t448;
t634 = -t673 * t810 + t676 * t813;
t617 = t634 * t674;
t804 = pkin(5) * t617;
t619 = t634 * t677;
t803 = pkin(5) * t619;
t724 = -t673 * t813 - t676 * t810;
t802 = pkin(5) * t724;
t801 = pkin(6) * t568;
t798 = pkin(16) * t674;
t797 = Icges(2,4) * t674;
t706 = t616 + t779;
t714 = pkin(1) * t719 - pkin(7);
t699 = t834 * (-pkin(1) * t715 - t706 * t714);
t695 = t699 * t777;
t704 = pkin(1) * t706;
t700 = t834 * (-t692 * t714 + t704 * t720);
t696 = t812 * t700;
t556 = t695 * t815 + t696 * t778;
t796 = Icges(7,4) * t556;
t697 = t815 * t700;
t557 = t695 * t812 + t697 * t777;
t795 = Icges(7,4) * t557;
t771 = pkin(5) * t808;
t794 = 0.2e1 * (t503 + t504) * t771 / t690;
t793 = t451 * t690;
t670 = cos(pkin(21));
t792 = t511 * t670;
t681 = 0.1e1 / pkin(11);
t791 = t511 * t681;
t683 = 0.1e1 / pkin(9);
t790 = t511 * t683;
t789 = t530 * t690;
t773 = pkin(6) * t809;
t788 = (t597 + t598) * t773 / t816;
t784 = t613 * t627;
t783 = 0.1e1 / pkin(13) * t689;
t668 = cos(pkin(23));
t756 = t668 * t778;
t664 = sin(pkin(23));
t757 = t664 * t777;
t551 = (t569 * t756 + t570 * t757) * t687;
t550 = 0.1e1 / t551 ^ 2;
t755 = t668 * t777;
t552 = (t569 * t757 + t570 * t755) * t687;
t435 = ((t541 * t757 + t542 * t755 + t664 * t743 + t668 * t742) / t551 - (t541 * t756 + t542 * t757 + t664 * t742 - t668 * t743) * t552 * t550) / (t550 * t552 ^ 2 + 0.1e1) * t687;
t776 = -qJD(2) - t435;
t566 = 0.1e1 / t567 ^ 2;
t605 = 0.1e1 / t606 ^ 2;
t623 = t734 * qJD(3);
t758 = -t788 / 0.2e1;
t457 = 0.2e1 * (((t621 * t758 + (t602 * t623 - t624 * t691) * pkin(6)) * t821 + (-t604 * t630 * t684 + t605 * t801) * t809) / t567 - ((-t624 * t602 - t623 * t691 + t630 * t758) * t821 + (t567 * t605 + t604 * t621) * t809) * t566 * t801) * pkin(2) / (t566 * t568 ^ 2 + 0.1e1) * t606 * t689;
t775 = -qJD(2) - t457;
t774 = -qJD(2) - qJD(3);
t770 = V_base(5) * pkin(14) + V_base(1);
t766 = Icges(3,4) * t813;
t765 = Icges(3,4) * t810;
t763 = t677 * t810;
t654 = qJD(2) * t674 + V_base(4);
t761 = t813 * t654;
t760 = -t794 / 0.2e1;
t759 = t511 * t819;
t658 = V_base(6) + qJD(1);
t512 = 0.1e1 / t513 ^ 2;
t751 = t512 * t771;
t748 = pkin(1) * t674 * t810;
t747 = pkin(1) * t763;
t746 = t834 * t832 * t812;
t433 = t435 * t674 + t654;
t455 = t457 * t674 + t654;
t632 = qJD(3) * t674 + t654;
t506 = t513 - t780;
t515 = -pkin(4) * t529 + pkin(5);
t447 = -pkin(4) * t789 + t506 * t515;
t745 = t447 * t751;
t449 = pkin(4) * t506 * t530 + t515 * t690;
t744 = t449 * t751;
t653 = -qJD(2) * t677 + V_base(5);
t741 = t653 * t813 * pkin(1) + t658 * t748 + t770;
t740 = -t798 - t804;
t739 = rSges(7,1) * t556 - rSges(7,2) * t557;
t738 = t834 * t815 * t777;
t450 = t470 * t657 - t471 * t656;
t731 = -t450 * t690 + t530 * t760;
t395 = ((-0.2e1 * pkin(5) * t515 - t506) * t451 + t731) * pkin(4);
t396 = t515 * t794 / 0.2e1 - 0.2e1 * t686 * t451 * t805 + (t450 * t506 - t793) * pkin(4);
t424 = (-t447 * t670 / 0.2e1 + t449 * t819) * t791;
t420 = 0.1e1 / t424 ^ 2;
t423 = (t447 * t819 + t449 * t670 / 0.2e1) * t791;
t350 = ((t395 * t759 + t666 * t745 + t396 * t792 / 0.2e1 + t670 * t744) / t424 - (-t395 * t792 / 0.2e1 - t670 * t745 + t396 * t759 + t666 * t744) * t423 * t420) / (t420 * t423 ^ 2 + 0.1e1) * t681;
t348 = t350 * t674 + t632;
t737 = Icges(7,1) * t556 - t795;
t736 = -Icges(7,2) * t557 + t796;
t735 = Icges(7,5) * t556 - Icges(7,6) * t557;
t631 = t677 * t774 + V_base(5);
t733 = -t631 * t802 + t741;
t732 = t658 * t677 * pkin(16) - V_base(4) * pkin(14) + V_base(2);
t730 = -rSges(3,1) * t810 - rSges(3,2) * t813;
t729 = -pkin(16) * t838 + V_base(3);
t728 = -Icges(3,2) * t813 - t765;
t532 = -t551 * t810 - t552 * t813;
t727 = -Icges(3,1) * t810 - t766;
t726 = -Icges(3,5) * t810 - Icges(3,6) * t813;
t725 = -t551 * t813 + t552 * t810;
t347 = V_base(5) + (-t350 + t774) * t677;
t540 = ((-0.3e1 * t835 + (0.4e1 * t717 - pkin(1)) * pkin(1) - t779) * t627 + t708) * pkin(1);
t543 = pkin(1) * t709 - t626 * t704 + t714 * t837;
t553 = 0.1e1 / t556 ^ 2;
t698 = t699 * t831;
t438 = ((t543 * t738 + t540 * t746 / 0.2e1 + (t697 * t831 + t698 * t812) * t784) / t556 - (t540 * t738 - t543 * t746 / 0.2e1 + (-t696 * t831 + t698 * t815) * t784) * t557 * t553) / (t553 * t557 ^ 2 + 0.1e1);
t436 = -t438 * t677 + V_base(5);
t437 = t438 * t674 + V_base(4);
t723 = (-Icges(7,3) * t677 + t674 * t735) * t436 + (Icges(7,3) * t674 + t677 * t735) * t437 + (Icges(7,5) * t557 + Icges(7,6) * t556) * t658;
t722 = (-Icges(3,3) * t677 + t674 * t726) * t653 + (Icges(3,3) * t674 + t677 * t726) * t654 + (Icges(3,5) * t813 - Icges(3,6) * t810) * t658;
t549 = (t568 * t767 - t810 * t567 / 0.2e1) * t785;
t710 = t653 * t747 - t654 * t748 + t729;
t707 = (-t658 * t763 - t761) * pkin(1) + t732;
t705 = -t631 * t803 + t632 * t804 + t710;
t703 = t632 * t802 + t658 * t803 + t707;
t519 = -Icges(7,6) * t677 + t674 * t736;
t520 = Icges(7,6) * t674 + t677 * t736;
t521 = -Icges(7,5) * t677 + t674 * t737;
t522 = Icges(7,5) * t674 + t677 * t737;
t534 = Icges(7,2) * t556 + t795;
t535 = Icges(7,1) * t557 + t796;
t702 = (-t520 * t557 + t522 * t556) * t437 + (-t519 * t557 + t521 * t556) * t436 + (-t534 * t557 + t535 * t556) * t658;
t609 = -Icges(3,6) * t677 + t674 * t728;
t610 = Icges(3,6) * t674 + t677 * t728;
t611 = -Icges(3,5) * t677 + t674 * t727;
t612 = Icges(3,5) * t674 + t677 * t727;
t643 = -Icges(3,2) * t810 + t766;
t646 = Icges(3,1) * t813 - t765;
t701 = (-t610 * t813 - t612 * t810) * t654 + (-t609 * t813 - t611 * t810) * t653 + (-t643 * t813 - t646 * t810) * t658;
t675 = cos(qJ(4));
t672 = sin(qJ(4));
t665 = sin(pkin(22));
t661 = Icges(2,4) * t677;
t651 = rSges(2,1) * t677 - rSges(2,2) * t674;
t650 = rSges(3,1) * t813 - rSges(3,2) * t810;
t649 = rSges(2,1) * t674 + rSges(2,2) * t677;
t648 = Icges(2,1) * t677 - t797;
t647 = Icges(2,1) * t674 + t661;
t645 = -Icges(2,2) * t674 + t661;
t644 = Icges(2,2) * t677 + t797;
t639 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t638 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t637 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t620 = t724 * t677;
t618 = t724 * t674;
t615 = t674 * rSges(3,3) + t677 * t730;
t614 = -t677 * rSges(3,3) + t674 * t730;
t600 = 0.1e1 / t601 ^ 2;
t596 = V_base(5) * rSges(2,3) - t649 * t658 + t770;
t595 = t651 * t658 + V_base(2) + (-pkin(14) - rSges(2,3)) * V_base(4);
t593 = t649 * V_base(4) - t651 * V_base(5) + V_base(3);
t592 = -rSges(4,1) * t724 + rSges(4,2) * t634;
t591 = -Icges(4,1) * t724 + Icges(4,4) * t634;
t590 = -Icges(4,4) * t724 + Icges(4,2) * t634;
t589 = -Icges(4,5) * t724 + Icges(4,6) * t634;
t587 = rSges(4,1) * t619 + rSges(4,2) * t620 + rSges(4,3) * t674;
t586 = rSges(4,1) * t617 + rSges(4,2) * t618 - rSges(4,3) * t677;
t585 = Icges(4,1) * t619 + Icges(4,4) * t620 + Icges(4,5) * t674;
t584 = Icges(4,1) * t617 + Icges(4,4) * t618 - Icges(4,5) * t677;
t583 = Icges(4,4) * t619 + Icges(4,2) * t620 + Icges(4,6) * t674;
t582 = Icges(4,4) * t617 + Icges(4,2) * t618 - Icges(4,6) * t677;
t581 = Icges(4,5) * t619 + Icges(4,6) * t620 + Icges(4,3) * t674;
t580 = Icges(4,5) * t617 + Icges(4,6) * t618 - Icges(4,3) * t677;
t576 = t650 * t653 + (-t614 - t798) * t658 + t770;
t575 = t615 * t658 - t650 * t654 + t732;
t574 = t614 * t654 - t615 * t653 + t729;
t561 = t592 * t631 + (-t586 - t798) * t658 + t741;
t560 = t658 * t587 - t632 * t592 + t707;
t559 = t632 * t586 - t631 * t587 + t710;
t547 = t677 * t836;
t546 = t677 * t549;
t545 = t674 * t836;
t544 = t674 * t549;
t539 = (0.1e1 / t601 * t788 / 0.2e1 - 0.2e1 * t691 * t600 * t773) / (-t600 * t786 + 0.1e1);
t536 = rSges(7,1) * t557 + rSges(7,2) * t556;
t528 = rSges(7,3) * t674 + t677 * t739;
t527 = -rSges(7,3) * t677 + t674 * t739;
t526 = t725 * t677;
t525 = t532 * t677;
t524 = t725 * t674;
t523 = t532 * t674;
t510 = -rSges(9,1) * t836 + rSges(9,2) * t549;
t509 = -Icges(9,1) * t836 + Icges(9,4) * t549;
t508 = -Icges(9,4) * t836 + Icges(9,2) * t549;
t507 = -Icges(9,5) * t836 + Icges(9,6) * t549;
t502 = rSges(9,1) * t546 + rSges(9,2) * t547 + rSges(9,3) * t674;
t501 = rSges(9,1) * t544 + rSges(9,2) * t545 - rSges(9,3) * t677;
t500 = Icges(9,1) * t546 + Icges(9,4) * t547 + Icges(9,5) * t674;
t499 = Icges(9,1) * t544 + Icges(9,4) * t545 - Icges(9,5) * t677;
t498 = Icges(9,4) * t546 + Icges(9,2) * t547 + Icges(9,6) * t674;
t497 = Icges(9,4) * t544 + Icges(9,2) * t545 - Icges(9,6) * t677;
t496 = Icges(9,5) * t546 + Icges(9,6) * t547 + Icges(9,3) * t674;
t495 = Icges(9,5) * t544 + Icges(9,6) * t545 - Icges(9,3) * t677;
t494 = (t549 * t822 - t816 * t836) * t783;
t493 = (t549 * t817 - t822 * t836) * t783;
t492 = (t546 * t816 + t547 * t822) * t783;
t491 = (t546 * t822 + t547 * t817) * t783;
t490 = (t544 * t816 + t545 * t822) * t783;
t489 = (t544 * t822 + t545 * t817) * t783;
t488 = -rSges(8,1) * t725 + rSges(8,2) * t532;
t487 = -Icges(8,1) * t725 + Icges(8,4) * t532;
t486 = -Icges(8,4) * t725 + Icges(8,2) * t532;
t485 = -Icges(8,5) * t725 + Icges(8,6) * t532;
t484 = (-t532 * t665 - t669 * t725) * pkin(4);
t483 = rSges(8,1) * t525 + rSges(8,2) * t526 + rSges(8,3) * t674;
t482 = rSges(8,1) * t523 + rSges(8,2) * t524 - rSges(8,3) * t677;
t481 = Icges(8,1) * t525 + Icges(8,4) * t526 + Icges(8,5) * t674;
t480 = Icges(8,1) * t523 + Icges(8,4) * t524 - Icges(8,5) * t677;
t479 = Icges(8,4) * t525 + Icges(8,2) * t526 + Icges(8,6) * t674;
t478 = Icges(8,4) * t523 + Icges(8,2) * t524 - Icges(8,6) * t677;
t477 = Icges(8,5) * t525 + Icges(8,6) * t526 + Icges(8,3) * t674;
t476 = Icges(8,5) * t523 + Icges(8,6) * t524 - Icges(8,3) * t677;
t475 = (t525 * t669 - t526 * t665) * pkin(4);
t474 = (t523 * t669 - t524 * t665) * pkin(4);
t469 = rSges(10,1) * t493 + rSges(10,2) * t494;
t468 = Icges(10,1) * t493 + Icges(10,4) * t494;
t467 = Icges(10,4) * t493 + Icges(10,2) * t494;
t466 = Icges(10,5) * t493 + Icges(10,6) * t494;
t465 = rSges(10,1) * t491 + rSges(10,2) * t492 + rSges(10,3) * t674;
t464 = rSges(10,1) * t489 + rSges(10,2) * t490 - rSges(10,3) * t677;
t463 = Icges(10,1) * t491 + Icges(10,4) * t492 + Icges(10,5) * t674;
t462 = Icges(10,1) * t489 + Icges(10,4) * t490 - Icges(10,5) * t677;
t461 = Icges(10,4) * t491 + Icges(10,2) * t492 + Icges(10,6) * t674;
t460 = Icges(10,4) * t489 + Icges(10,2) * t490 - Icges(10,6) * t677;
t459 = Icges(10,5) * t491 + Icges(10,6) * t492 + Icges(10,3) * t674;
t458 = Icges(10,5) * t489 + Icges(10,6) * t490 - Icges(10,3) * t677;
t454 = t677 * t775 + V_base(5);
t453 = t539 * t674 + t455;
t452 = V_base(5) + (-t539 + t775) * t677;
t446 = -pkin(5) * t789 - t505 * t514;
t445 = 0.1e1 / t446 ^ 2;
t444 = t454 * t510 + (-t501 - t798) * t658 + t770;
t443 = -t455 * t510 + t502 * t658 + t732;
t432 = t677 * t776 + V_base(5);
t431 = -t454 * t502 + t455 * t501 + t729;
t430 = -V_base(5) * pkin(17) + t436 * t536 + (pkin(15) * t674 - t527) * t658 + t770;
t429 = -t437 * t536 + V_base(2) + (-pkin(14) + pkin(17)) * V_base(4) + (-pkin(15) * t677 + t528) * t658;
t427 = t432 * t488 + (-t482 - t798) * t658 + t741;
t426 = -t433 * t488 + t658 * t483 + t707;
t422 = (t446 * t818 - t665 * t448 / 0.2e1) * t790;
t421 = (t665 * t446 / 0.2e1 + t448 * t818) * t790;
t419 = -pkin(2) * t454 * t836 + t452 * t469 + (-pkin(2) * t544 - t464 - t798) * t658 + t770;
t418 = -t453 * t469 + t465 * t658 + (t455 * t836 + t546 * t658) * pkin(2) + t732;
t417 = pkin(15) * t838 - t436 * t528 + t437 * t527 + V_base(3);
t416 = -t432 * t483 + t433 * t482 + t710;
t415 = -t452 * t465 + t453 * t464 + (-t454 * t546 + t455 * t544) * pkin(2) + t729;
t413 = t423 * t634 - t424 * t724;
t412 = t423 * t724 + t424 * t634;
t411 = -qJD(4) * t412 + t658;
t410 = t423 * t620 + t424 * t619;
t409 = -t423 * t619 + t424 * t620;
t408 = t423 * t618 + t424 * t617;
t407 = -t423 * t617 + t424 * t618;
t406 = t410 * t675 + t672 * t674;
t405 = -t410 * t672 + t674 * t675;
t404 = t408 * t675 - t672 * t677;
t403 = -t408 * t672 - t675 * t677;
t402 = t421 * t725 + t422 * t532;
t401 = t421 * t532 - t422 * t725;
t400 = -t421 * t525 + t422 * t526;
t399 = t421 * t526 + t422 * t525;
t398 = -t421 * t523 + t422 * t524;
t397 = t421 * t524 + t422 * t523;
t394 = pkin(10) * t413 - pkin(12) * t412;
t393 = rSges(5,1) * t413 + rSges(5,2) * t412;
t392 = Icges(5,1) * t413 + Icges(5,4) * t412;
t391 = Icges(5,4) * t413 + Icges(5,2) * t412;
t390 = Icges(5,5) * t413 + Icges(5,6) * t412;
t389 = pkin(10) * t410 - pkin(12) * t409;
t388 = pkin(10) * t408 - pkin(12) * t407;
t387 = rSges(5,1) * t410 + rSges(5,2) * t409 + rSges(5,3) * t674;
t386 = rSges(5,1) * t408 + rSges(5,2) * t407 - rSges(5,3) * t677;
t385 = Icges(5,1) * t410 + Icges(5,4) * t409 + Icges(5,5) * t674;
t384 = Icges(5,1) * t408 + Icges(5,4) * t407 - Icges(5,5) * t677;
t383 = Icges(5,4) * t410 + Icges(5,2) * t409 + Icges(5,6) * t674;
t382 = Icges(5,4) * t408 + Icges(5,2) * t407 - Icges(5,6) * t677;
t381 = Icges(5,5) * t410 + Icges(5,6) * t409 + Icges(5,3) * t674;
t380 = Icges(5,5) * t408 + Icges(5,6) * t407 - Icges(5,3) * t677;
t379 = rSges(11,1) * t401 + rSges(11,2) * t402;
t378 = Icges(11,1) * t401 + Icges(11,4) * t402;
t377 = Icges(11,4) * t401 + Icges(11,2) * t402;
t376 = Icges(11,5) * t401 + Icges(11,6) * t402;
t375 = rSges(11,1) * t399 + rSges(11,2) * t400 + rSges(11,3) * t674;
t374 = rSges(11,1) * t397 + rSges(11,2) * t398 - rSges(11,3) * t677;
t373 = Icges(11,1) * t399 + Icges(11,4) * t400 + Icges(11,5) * t674;
t372 = Icges(11,1) * t397 + Icges(11,4) * t398 - Icges(11,5) * t677;
t371 = Icges(11,4) * t399 + Icges(11,2) * t400 + Icges(11,6) * t674;
t370 = Icges(11,4) * t397 + Icges(11,2) * t398 - Icges(11,6) * t677;
t369 = Icges(11,5) * t399 + Icges(11,6) * t400 + Icges(11,3) * t674;
t368 = Icges(11,5) * t397 + Icges(11,6) * t398 - Icges(11,3) * t677;
t367 = -rSges(6,3) * t412 + (rSges(6,1) * t675 - rSges(6,2) * t672) * t413;
t366 = -Icges(6,5) * t412 + (Icges(6,1) * t675 - Icges(6,4) * t672) * t413;
t365 = -Icges(6,6) * t412 + (Icges(6,4) * t675 - Icges(6,2) * t672) * t413;
t364 = -Icges(6,3) * t412 + (Icges(6,5) * t675 - Icges(6,6) * t672) * t413;
t363 = rSges(6,1) * t406 + rSges(6,2) * t405 - rSges(6,3) * t409;
t362 = rSges(6,1) * t404 + rSges(6,2) * t403 - rSges(6,3) * t407;
t361 = Icges(6,1) * t406 + Icges(6,4) * t405 - Icges(6,5) * t409;
t360 = Icges(6,1) * t404 + Icges(6,4) * t403 - Icges(6,5) * t407;
t359 = Icges(6,4) * t406 + Icges(6,2) * t405 - Icges(6,6) * t409;
t358 = Icges(6,4) * t404 + Icges(6,2) * t403 - Icges(6,6) * t407;
t357 = Icges(6,5) * t406 + Icges(6,6) * t405 - Icges(6,3) * t409;
t356 = Icges(6,5) * t404 + Icges(6,6) * t403 - Icges(6,3) * t407;
t355 = 0.2e1 * (((t514 * t760 + (t450 * t505 - t793) * pkin(5)) * t823 + (-t511 * t530 * t685 + t512 * t807) * t808) / t446 - ((-t451 * t505 + t731) * t823 + (t446 * t512 + t511 * t514) * t808) * t445 * t807) * pkin(9) / (t445 * t448 ^ 2 + 0.1e1) * t513 * t683;
t354 = t355 * t674 + t433;
t353 = V_base(5) + (-t355 + t776) * t677;
t352 = t353 * t379 + t432 * t484 + (-t374 - t474 - t798) * t658 + t741;
t351 = -pkin(1) * t761 - t354 * t379 - t433 * t484 + (-t747 + t375 + t475) * t658 + t732;
t346 = -qJD(4) * t409 + t348;
t345 = -qJD(4) * t407 + t347;
t344 = t347 * t393 + (-t386 + t740) * t658 + t733;
t343 = -t348 * t393 + t658 * t387 + t703;
t342 = -t353 * t375 + t354 * t374 - t432 * t475 + t433 * t474 + t710;
t341 = -t347 * t387 + t348 * t386 + t705;
t340 = t345 * t367 + t347 * t394 - t362 * t411 + (-t388 + t740) * t658 + t733;
t339 = -t346 * t367 - t348 * t394 + t411 * t363 + t658 * t389 + t703;
t338 = -t345 * t363 + t346 * t362 - t347 * t389 + t348 * t388 + t705;
t1 = ((t644 * t677 + t647 * t674 + Icges(1,2)) * V_base(5) + (t645 * t677 + t648 * t674 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t644 * t674 + t647 * t677 + Icges(1,4)) * V_base(5) + (-t645 * t674 + t648 * t677 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t411 * ((-t345 * t356 - t346 * t357 - t364 * t411) * t412 + ((-t359 * t672 + t361 * t675) * t346 + (-t358 * t672 + t360 * t675) * t345 + (-t365 * t672 + t366 * t675) * t411) * t413) / 0.2e1 + m(11) * (t342 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + ((t383 * t412 + t385 * t413) * t348 + (t382 * t412 + t384 * t413) * t347 + (t461 * t494 + t463 * t493) * t453 + (t460 * t494 + t462 * t493) * t452 + (-t610 * t810 + t612 * t813) * t654 + (-t609 * t810 + t611 * t813) * t653 + (t583 * t634 - t585 * t724) * t632 + (t582 * t634 - t584 * t724) * t631 + (t520 * t556 + t522 * t557) * t437 + (t519 * t556 + t521 * t557) * t436 + (t498 * t549 - t500 * t836) * t455 + (t497 * t549 - t499 * t836) * t454 + (t371 * t402 + t373 * t401) * t354 + (t370 * t402 + t372 * t401) * t353 + (t479 * t532 - t481 * t725) * t433 + (t478 * t532 - t480 * t725) * t432 + (Icges(2,3) - t643 * t810 + t646 * t813 + t391 * t412 + t392 * t413 + t467 * t494 + t468 * t493 + t590 * t634 - t591 * t724 + t534 * t556 + t535 * t557 + t508 * t549 - t509 * t836 + t377 * t402 + t378 * t401 + t486 * t532 - t487 * t725) * t658) * t658 / 0.2e1 + m(5) * (t341 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(7) * (t417 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(8) * (t416 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + t345 * ((-t357 * t407 + t359 * t403 + t361 * t404) * t346 + (-t356 * t407 + t358 * t403 + t360 * t404) * t345 + (-t364 * t407 + t365 * t403 + t366 * t404) * t411) / 0.2e1 + t346 * ((-t357 * t409 + t359 * t405 + t361 * t406) * t346 + (-t356 * t409 + t358 * t405 + t360 * t406) * t345 + (-t364 * t409 + t365 * t405 + t366 * t406) * t411) / 0.2e1 + t658 * V_base(5) * (Icges(2,5) * t674 + Icges(2,6) * t677) + t658 * V_base(4) * (Icges(2,5) * t677 - Icges(2,6) * t674) + m(10) * (t415 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(9) * (t431 ^ 2 + t443 ^ 2 + t444 ^ 2) / 0.2e1 + t347 * ((-t381 * t677 + t383 * t407 + t385 * t408) * t348 + (-t380 * t677 + t382 * t407 + t384 * t408) * t347 + (-t390 * t677 + t391 * t407 + t392 * t408) * t658) / 0.2e1 + t452 * ((-t459 * t677 + t461 * t490 + t463 * t489) * t453 + (-t458 * t677 + t460 * t490 + t462 * t489) * t452 + (-t466 * t677 + t467 * t490 + t468 * t489) * t658) / 0.2e1 + t631 * ((-t581 * t677 + t583 * t618 + t585 * t617) * t632 + (-t580 * t677 + t582 * t618 + t584 * t617) * t631 + (-t589 * t677 + t590 * t618 + t591 * t617) * t658) / 0.2e1 + t454 * ((-t496 * t677 + t498 * t545 + t500 * t544) * t455 + (-t495 * t677 + t497 * t545 + t499 * t544) * t454 + (-t507 * t677 + t508 * t545 + t509 * t544) * t658) / 0.2e1 + t432 * ((-t477 * t677 + t479 * t524 + t481 * t523) * t433 + (-t476 * t677 + t478 * t524 + t480 * t523) * t432 + (-t485 * t677 + t486 * t524 + t487 * t523) * t658) / 0.2e1 + t353 * ((-t369 * t677 + t371 * t398 + t373 * t397) * t354 + (-t368 * t677 + t370 * t398 + t372 * t397) * t353 + (-t376 * t677 + t377 * t398 + t378 * t397) * t658) / 0.2e1 + t354 * ((t369 * t674 + t371 * t400 + t373 * t399) * t354 + (t368 * t674 + t370 * t400 + t372 * t399) * t353 + (t376 * t674 + t377 * t400 + t378 * t399) * t658) / 0.2e1 + t453 * ((t459 * t674 + t461 * t492 + t463 * t491) * t453 + (t458 * t674 + t460 * t492 + t462 * t491) * t452 + (t466 * t674 + t467 * t492 + t468 * t491) * t658) / 0.2e1 + t632 * ((t581 * t674 + t583 * t620 + t585 * t619) * t632 + (t580 * t674 + t582 * t620 + t584 * t619) * t631 + (t589 * t674 + t590 * t620 + t591 * t619) * t658) / 0.2e1 + t455 * ((t496 * t674 + t498 * t547 + t500 * t546) * t455 + (t495 * t674 + t497 * t547 + t499 * t546) * t454 + (t507 * t674 + t508 * t547 + t509 * t546) * t658) / 0.2e1 + t433 * ((t477 * t674 + t479 * t526 + t481 * t525) * t433 + (t476 * t674 + t478 * t526 + t480 * t525) * t432 + (t485 * t674 + t486 * t526 + t487 * t525) * t658) / 0.2e1 + t348 * ((t381 * t674 + t383 * t409 + t385 * t410) * t348 + (t380 * t674 + t382 * t409 + t384 * t410) * t347 + (t390 * t674 + t391 * t409 + t392 * t410) * t658) / 0.2e1 + t654 * (t674 * t722 + t677 * t701) / 0.2e1 + t653 * (t674 * t701 - t677 * t722) / 0.2e1 + t437 * (t723 * t674 + t702 * t677) / 0.2e1 + t436 * (t702 * t674 - t723 * t677) / 0.2e1 + m(1) * (t637 ^ 2 + t638 ^ 2 + t639 ^ 2) / 0.2e1 + m(2) * (t593 ^ 2 + t595 ^ 2 + t596 ^ 2) / 0.2e1 + m(3) * (t574 ^ 2 + t575 ^ 2 + t576 ^ 2) / 0.2e1 + m(4) * (t559 ^ 2 + t560 ^ 2 + t561 ^ 2) / 0.2e1;
T = t1;
