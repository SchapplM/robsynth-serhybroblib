% Calculate kinetic energy for
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1DE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1DE1_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:47:56
% EndTime: 2020-04-13 14:50:06
% DurationCPUTime: 130.08s
% Computational Cost: add. (2862993->768), mult. (4363736->1244), div. (195296->33), fcn. (2746388->46), ass. (0->457)
t708 = sin(qJ(1));
t714 = cos(qJ(1));
t865 = -t708 * V_base(4) + V_base(5) * t714;
t864 = -2 * pkin(1);
t863 = (-pkin(2) - pkin(13));
t862 = (-pkin(2) + pkin(13));
t861 = -pkin(3) - pkin(8);
t860 = -pkin(8) + pkin(3);
t859 = -pkin(9) - pkin(11);
t858 = -pkin(9) + pkin(11);
t728 = pkin(4) ^ 2;
t727 = pkin(5) ^ 2;
t725 = pkin(7) ^ 2;
t733 = pkin(1) ^ 2;
t707 = sin(qJ(2));
t709 = sin(pkin(19));
t713 = cos(qJ(2));
t715 = cos(pkin(19));
t664 = t707 * t715 - t709 * t713;
t836 = pkin(7) * t664;
t807 = t836 * t864 + t733;
t645 = t725 + t807;
t804 = pkin(3) ^ 2 - pkin(8) ^ 2;
t631 = t645 + t804;
t654 = pkin(1) - t836;
t667 = t707 * t709 + t713 * t715;
t624 = (pkin(7) - t861) * (pkin(7) + t861) + t807;
t625 = (pkin(7) - t860) * (pkin(7) + t860) + t807;
t736 = sqrt(-t625 * t624);
t590 = pkin(7) * t631 * t667 + t654 * t736;
t712 = cos(qJ(3));
t811 = t712 * t590;
t816 = t667 * t736;
t589 = -pkin(7) * t816 + t631 * t654;
t706 = sin(qJ(3));
t815 = t706 * t589;
t751 = t815 / 0.2e1 + t811 / 0.2e1;
t641 = 0.1e1 / t645;
t730 = 0.1e1 / pkin(3);
t818 = t641 * t730;
t577 = t751 * t818;
t812 = t712 * t589;
t814 = t706 * t590;
t750 = -t812 / 0.2e1 + t814 / 0.2e1;
t578 = t750 * t818;
t696 = pkin(23) + pkin(22);
t689 = sin(t696);
t690 = cos(t696);
t555 = -t577 * t690 + t578 * t689;
t844 = pkin(5) * t555;
t809 = -0.2e1 * pkin(4) * t844 + t727;
t551 = t728 + t809;
t549 = 0.1e1 / t551;
t857 = t549 / 0.2e1;
t658 = t667 * qJD(2);
t657 = t664 * qJD(2);
t799 = pkin(1) * pkin(7) * t658;
t823 = 0.2e1 * (t624 + t625) * t799 / t736;
t782 = -t823 / 0.2e1;
t748 = t657 * t736 + t667 * t782;
t567 = ((t654 * t864 - t631) * t658 + t748) * pkin(7);
t856 = -t567 / 0.2e1;
t796 = -0.2e1 * t658 * t667;
t817 = t658 * t736;
t568 = t654 * t823 / 0.2e1 + t725 * pkin(1) * t796 + (-t631 * t657 - t817) * pkin(7);
t855 = t568 / 0.2e1;
t700 = sin(pkin(20));
t704 = cos(pkin(20));
t754 = t700 * t712 + t704 * t706;
t838 = pkin(6) * t754;
t797 = pkin(1) * t838;
t651 = 0.2e1 * t797;
t726 = pkin(6) ^ 2;
t805 = t726 + t733;
t634 = t651 + t805;
t632 = 0.1e1 / t634;
t854 = t632 / 0.2e1;
t697 = sin(pkin(23));
t852 = t697 / 0.2e1;
t699 = sin(pkin(21));
t851 = t699 / 0.2e1;
t850 = -t712 / 0.2e1;
t716 = cos(pkin(18));
t849 = t716 / 0.2e1;
t732 = 0.1e1 / pkin(2);
t848 = t732 / 0.2e1;
t661 = t700 * t706 - t704 * t712;
t653 = t661 * qJD(3);
t847 = pkin(1) * t653;
t777 = 0.1e1 / t645 ^ 2 * t799;
t505 = ((-t812 + t814) * t777 + (qJD(3) * t751 + t567 * t850 + t706 * t855) * t641) * t730;
t506 = ((-t811 - t815) * t777 + (qJD(3) * t750 + t568 * t850 + t706 * t856) * t641) * t730;
t477 = t505 * t689 + t506 * t690;
t846 = pkin(4) * t477;
t806 = pkin(9) ^ 2 - pkin(11) ^ 2;
t545 = t551 + t806;
t552 = -pkin(4) + t844;
t543 = (pkin(4) - t859) * (pkin(4) + t859) + t809;
t544 = (pkin(4) - t858) * (pkin(4) + t858) + t809;
t734 = sqrt(-t544 * t543);
t556 = t577 * t689 + t578 * t690;
t843 = pkin(5) * t556;
t462 = t545 * t843 - t552 * t734;
t845 = pkin(5) * t462;
t666 = -t706 * t707 + t712 * t713;
t646 = t666 * t708;
t842 = pkin(5) * t646;
t648 = t666 * t714;
t841 = pkin(5) * t648;
t665 = t706 * t713 + t707 * t712;
t840 = pkin(5) * t665;
t731 = pkin(2) ^ 2;
t789 = -pkin(13) ^ 2 + t805;
t629 = t651 + t731 + t789;
t650 = -pkin(1) - t838;
t808 = t651 + t726;
t622 = ((pkin(1) - t863) * (pkin(1) + t863)) + t808;
t623 = ((pkin(1) - t862) * (pkin(1) + t862)) + t808;
t822 = t623 * t622;
t735 = sqrt(-t822);
t837 = pkin(6) * t661;
t587 = t629 * t837 - t650 * t735;
t839 = pkin(6) * t587;
t835 = pkin(16) * t708;
t834 = Icges(2,4) * t708;
t833 = Icges(3,4) * t707;
t832 = Icges(3,4) * t713;
t630 = t645 - t804;
t655 = pkin(1) * t664 - pkin(7);
t588 = -pkin(1) * t816 - t630 * t655;
t591 = pkin(1) * t630 * t667 - t655 * t736;
t710 = sin(pkin(18));
t724 = 0.1e1 / pkin(8);
t819 = t641 * t724;
t579 = (t588 * t849 - t710 * t591 / 0.2e1) * t819;
t580 = (t591 * t849 + t588 * t710 / 0.2e1) * t819;
t564 = atan2(t580, t579);
t560 = sin(t564);
t831 = Icges(7,4) * t560;
t561 = cos(t564);
t830 = Icges(7,4) * t561;
t798 = pkin(5) * t846;
t829 = 0.2e1 * (t543 + t544) * t798 / t734;
t828 = t477 * t734;
t703 = cos(pkin(21));
t827 = t549 * t703;
t720 = 0.1e1 / pkin(11);
t826 = t549 * t720;
t825 = t556 * t734;
t800 = pkin(6) * t847;
t824 = 0.2e1 * (t622 + t623) * t800 / t735;
t701 = cos(pkin(23));
t821 = t641 * t701;
t820 = t641 * t710;
t813 = t707 * t714;
t687 = qJD(2) * t708 + V_base(4);
t810 = t713 * t687;
t574 = (-t701 * t589 / 0.2e1 + t590 * t852) * t818;
t573 = 0.1e1 / t574 ^ 2;
t575 = (t701 * t590 / 0.2e1 + t589 * t852) * t818;
t767 = t590 * t777;
t768 = t589 * t777;
t780 = t641 * t852;
t455 = ((t567 * t780 + t697 * t768 + t701 * t767 + t821 * t855) / t574 - (t568 * t780 + t697 * t767 - t701 * t768 + t821 * t856) * t575 * t573) / (t573 * t575 ^ 2 + 0.1e1) * t730;
t803 = -qJD(2) - t455;
t586 = -t629 * t650 - t735 * t837;
t585 = 0.1e1 / t586 ^ 2;
t633 = 0.1e1 / t634 ^ 2;
t652 = t754 * qJD(3);
t783 = -t824 / 0.2e1;
t498 = 0.2e1 * (((t650 * t783 + (t652 * t629 - t653 * t735) * pkin(6)) * t854 + (-t632 * t661 * t726 + t633 * t839) * t847) / t586 - ((-t653 * t629 - t652 * t735 + t661 * t783) * t854 + (t586 * t633 + t632 * t650) * t847) * t585 * t839) * pkin(2) / (t585 * t587 ^ 2 + 0.1e1) * t634 * t732;
t802 = -qJD(2) - t498;
t801 = -qJD(2) - qJD(3);
t795 = V_base(5) * pkin(14) + V_base(1);
t794 = pkin(1) * t707 * t708;
t793 = pkin(1) * t813;
t546 = t551 - t806;
t553 = -pkin(4) * t555 + pkin(5);
t461 = -pkin(4) * t825 + t546 * t553;
t463 = pkin(4) * t546 * t556 + t553 * t734;
t442 = (t461 * t851 + t463 * t703 / 0.2e1) * t826;
t443 = (-t461 * t703 / 0.2e1 + t463 * t851) * t826;
t790 = atan2(t442, t443);
t786 = -t829 / 0.2e1;
t785 = t549 * t851;
t722 = 0.1e1 / pkin(9);
t784 = t722 * t857;
t781 = t632 * t848;
t779 = t641 * t849;
t778 = 0.1e1 / pkin(13) * t848;
t691 = V_base(6) + qJD(1);
t550 = 0.1e1 / t551 ^ 2;
t776 = t550 * t798;
t773 = cos(t790);
t453 = t455 * t708 + t687;
t496 = t498 * t708 + t687;
t663 = qJD(3) * t708 + t687;
t772 = t710 * t777;
t771 = t716 * t777;
t770 = t463 * t776;
t769 = t461 * t776;
t686 = -qJD(2) * t714 + V_base(5);
t766 = t686 * t713 * pkin(1) + t691 * t794 + t795;
t765 = -t835 - t842;
t764 = -rSges(3,1) * t707 - rSges(3,2) * t713;
t763 = rSges(7,1) * t561 - rSges(7,2) * t560;
t476 = t505 * t690 - t506 * t689;
t749 = -t476 * t734 + t556 * t786;
t424 = ((-0.2e1 * pkin(5) * t553 - t546) * t477 + t749) * pkin(4);
t425 = t553 * t829 / 0.2e1 - 0.2e1 * t728 * t477 * t843 + (t476 * t546 - t828) * pkin(4);
t441 = 0.1e1 / t443 ^ 2;
t364 = ((t424 * t785 + t699 * t769 + t425 * t827 / 0.2e1 + t703 * t770) / t443 - (-t424 * t827 / 0.2e1 - t703 * t769 + t425 * t785 + t699 * t770) * t442 * t441) / (t441 * t442 ^ 2 + 0.1e1) * t720;
t362 = t364 * t708 + t663;
t762 = -Icges(3,1) * t707 - t832;
t761 = Icges(7,1) * t561 - t831;
t760 = -Icges(3,2) * t713 - t833;
t759 = -Icges(7,2) * t560 + t830;
t758 = -Icges(3,5) * t707 - Icges(3,6) * t713;
t757 = Icges(7,5) * t561 - Icges(7,6) * t560;
t562 = atan2(t575, t574);
t557 = sin(t562);
t558 = cos(t562);
t533 = -t557 * t713 - t558 * t707;
t756 = t557 * t707 - t558 * t713;
t572 = atan2(t587 * t781, t586 * t781);
t570 = sin(t572);
t571 = cos(t572);
t547 = -t570 * t713 - t571 * t707;
t755 = t570 * t707 - t571 * t713;
t662 = t714 * t801 + V_base(5);
t753 = t662 * t840 + t766;
t752 = t691 * t714 * pkin(16) - V_base(4) * pkin(14) + V_base(2);
t747 = -pkin(16) * t865 + V_base(3);
t361 = V_base(5) + (-t364 + t801) * t714;
t566 = ((0.2e1 * pkin(7) * t655 - t630) * t658 + t748) * pkin(1);
t569 = t655 * t782 + t733 * pkin(7) * t796 + (-t630 * t657 - t817) * pkin(1);
t576 = 0.1e1 / t579 ^ 2;
t458 = ((t569 * t779 + t591 * t771 + t566 * t820 / 0.2e1 + t588 * t772) / t579 - (t566 * t779 + t588 * t771 - t569 * t820 / 0.2e1 - t591 * t772) * t580 * t576) / (t576 * t580 ^ 2 + 0.1e1) * t724;
t456 = -t458 * t714 + V_base(5);
t457 = t458 * t708 + V_base(4);
t746 = (-Icges(7,3) * t714 + t708 * t757) * t456 + (Icges(7,3) * t708 + t714 * t757) * t457 + (Icges(7,5) * t560 + Icges(7,6) * t561) * t691;
t745 = (-Icges(3,3) * t714 + t708 * t758) * t686 + (Icges(3,3) * t708 + t714 * t758) * t687 + (Icges(3,5) * t713 - Icges(3,6) * t707) * t691;
t460 = -pkin(5) * t825 - t545 * t552;
t744 = atan2(t462 * t784, t460 * t784);
t743 = sin(t744);
t742 = t686 * t793 - t687 * t794 + t747;
t741 = (-t691 * t813 - t810) * pkin(1) + t752;
t740 = -t662 * t841 + t663 * t842 + t742;
t739 = -t663 * t840 + t691 * t841 + t741;
t527 = -Icges(7,6) * t714 + t708 * t759;
t528 = Icges(7,6) * t708 + t714 * t759;
t529 = -Icges(7,5) * t714 + t708 * t761;
t530 = Icges(7,5) * t708 + t714 * t761;
t536 = Icges(7,2) * t561 + t831;
t537 = Icges(7,1) * t560 + t830;
t738 = (-t528 * t560 + t530 * t561) * t457 + (-t527 * t560 + t529 * t561) * t456 + (-t536 * t560 + t537 * t561) * t691;
t637 = -Icges(3,6) * t714 + t708 * t760;
t638 = Icges(3,6) * t708 + t714 * t760;
t639 = -Icges(3,5) * t714 + t708 * t762;
t640 = Icges(3,5) * t708 + t714 * t762;
t676 = -Icges(3,2) * t707 + t832;
t679 = Icges(3,1) * t713 - t833;
t737 = (-t638 * t713 - t640 * t707) * t687 + (-t637 * t713 - t639 * t707) * t686 + (-t676 * t713 - t679 * t707) * t691;
t711 = cos(qJ(4));
t705 = sin(qJ(4));
t702 = cos(pkin(22));
t698 = sin(pkin(22));
t694 = Icges(2,4) * t714;
t684 = rSges(2,1) * t714 - rSges(2,2) * t708;
t683 = rSges(3,1) * t713 - rSges(3,2) * t707;
t682 = rSges(2,1) * t708 + rSges(2,2) * t714;
t681 = Icges(2,1) * t714 - t834;
t680 = Icges(2,1) * t708 + t694;
t678 = -Icges(2,2) * t708 + t694;
t677 = Icges(2,2) * t714 + t834;
t672 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t671 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t670 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t649 = t665 * t714;
t647 = t665 * t708;
t644 = rSges(3,3) * t708 + t714 * t764;
t643 = -rSges(3,3) * t714 + t708 * t764;
t628 = t731 - t789 - 0.2e1 * t797;
t627 = 0.1e1 / t628 ^ 2;
t621 = V_base(5) * rSges(2,3) - t682 * t691 + t795;
t620 = t684 * t691 + V_base(2) + (-pkin(14) - rSges(2,3)) * V_base(4);
t618 = t682 * V_base(4) - t684 * V_base(5) + V_base(3);
t617 = rSges(4,1) * t665 + rSges(4,2) * t666;
t616 = Icges(4,1) * t665 + Icges(4,4) * t666;
t615 = Icges(4,4) * t665 + Icges(4,2) * t666;
t614 = Icges(4,5) * t665 + Icges(4,6) * t666;
t612 = rSges(4,1) * t648 - rSges(4,2) * t649 + rSges(4,3) * t708;
t611 = rSges(4,1) * t646 - rSges(4,2) * t647 - rSges(4,3) * t714;
t610 = Icges(4,1) * t648 - Icges(4,4) * t649 + Icges(4,5) * t708;
t609 = Icges(4,1) * t646 - Icges(4,4) * t647 - Icges(4,5) * t714;
t608 = Icges(4,4) * t648 - Icges(4,2) * t649 + Icges(4,6) * t708;
t607 = Icges(4,4) * t646 - Icges(4,2) * t647 - Icges(4,6) * t714;
t606 = Icges(4,5) * t648 - Icges(4,6) * t649 + Icges(4,3) * t708;
t605 = Icges(4,5) * t646 - Icges(4,6) * t647 - Icges(4,3) * t714;
t600 = t683 * t686 + (-t643 - t835) * t691 + t795;
t599 = t644 * t691 - t683 * t687 + t752;
t598 = t643 * t687 - t644 * t686 + t747;
t596 = atan2(t735 * t778, t628 * t778);
t594 = cos(t596);
t593 = sin(t596);
t584 = t617 * t662 + (-t611 - t835) * t691 + t766;
t583 = t612 * t691 - t617 * t663 + t741;
t582 = t611 * t663 - t612 * t662 + t742;
t565 = (0.1e1 / t628 * t824 / 0.2e1 - 0.2e1 * t735 * t627 * t800) / (-t627 * t822 + 0.1e1);
t542 = t547 * t714;
t541 = t755 * t714;
t540 = t547 * t708;
t539 = t755 * t708;
t538 = rSges(7,1) * t560 + rSges(7,2) * t561;
t532 = rSges(7,3) * t708 + t714 * t763;
t531 = -rSges(7,3) * t714 + t708 * t763;
t524 = t533 * t714;
t523 = t756 * t714;
t522 = t533 * t708;
t521 = t756 * t708;
t520 = -rSges(9,1) * t755 + rSges(9,2) * t547;
t519 = -Icges(9,1) * t755 + Icges(9,4) * t547;
t518 = -Icges(9,4) * t755 + Icges(9,2) * t547;
t517 = -Icges(9,5) * t755 + Icges(9,6) * t547;
t514 = rSges(9,1) * t542 + rSges(9,2) * t541 + rSges(9,3) * t708;
t513 = rSges(9,1) * t540 + rSges(9,2) * t539 - rSges(9,3) * t714;
t512 = Icges(9,1) * t542 + Icges(9,4) * t541 + Icges(9,5) * t708;
t511 = Icges(9,1) * t540 + Icges(9,4) * t539 - Icges(9,5) * t714;
t510 = Icges(9,4) * t542 + Icges(9,2) * t541 + Icges(9,6) * t708;
t509 = Icges(9,4) * t540 + Icges(9,2) * t539 - Icges(9,6) * t714;
t508 = Icges(9,5) * t542 + Icges(9,6) * t541 + Icges(9,3) * t708;
t507 = Icges(9,5) * t540 + Icges(9,6) * t539 - Icges(9,3) * t714;
t504 = -t547 * t593 + t594 * t755;
t503 = -t547 * t594 - t593 * t755;
t502 = -t541 * t593 - t542 * t594;
t501 = -t541 * t594 + t542 * t593;
t500 = -t539 * t593 - t540 * t594;
t499 = -t539 * t594 + t540 * t593;
t495 = t714 * t802 + V_base(5);
t494 = -rSges(8,1) * t756 + rSges(8,2) * t533;
t493 = -Icges(8,1) * t756 + Icges(8,4) * t533;
t492 = -Icges(8,4) * t756 + Icges(8,2) * t533;
t491 = -Icges(8,5) * t756 + Icges(8,6) * t533;
t490 = (-t533 * t698 - t702 * t756) * pkin(4);
t489 = rSges(8,1) * t524 + rSges(8,2) * t523 + rSges(8,3) * t708;
t488 = rSges(8,1) * t522 + rSges(8,2) * t521 - rSges(8,3) * t714;
t487 = Icges(8,1) * t524 + Icges(8,4) * t523 + Icges(8,5) * t708;
t486 = Icges(8,1) * t522 + Icges(8,4) * t521 - Icges(8,5) * t714;
t485 = Icges(8,4) * t524 + Icges(8,2) * t523 + Icges(8,6) * t708;
t484 = Icges(8,4) * t522 + Icges(8,2) * t521 - Icges(8,6) * t714;
t483 = Icges(8,5) * t524 + Icges(8,6) * t523 + Icges(8,3) * t708;
t482 = Icges(8,5) * t522 + Icges(8,6) * t521 - Icges(8,3) * t714;
t481 = (-t523 * t698 + t524 * t702) * pkin(4);
t480 = (-t521 * t698 + t522 * t702) * pkin(4);
t479 = t565 * t708 + t496;
t478 = V_base(5) + (-t565 + t802) * t714;
t475 = rSges(10,1) * t504 + rSges(10,2) * t503;
t474 = Icges(10,1) * t504 + Icges(10,4) * t503;
t473 = Icges(10,4) * t504 + Icges(10,2) * t503;
t472 = Icges(10,5) * t504 + Icges(10,6) * t503;
t471 = rSges(10,1) * t502 + rSges(10,2) * t501 + rSges(10,3) * t708;
t470 = rSges(10,1) * t500 + rSges(10,2) * t499 - rSges(10,3) * t714;
t469 = Icges(10,1) * t502 + Icges(10,4) * t501 + Icges(10,5) * t708;
t468 = Icges(10,1) * t500 + Icges(10,4) * t499 - Icges(10,5) * t714;
t467 = Icges(10,4) * t502 + Icges(10,2) * t501 + Icges(10,6) * t708;
t466 = Icges(10,4) * t500 + Icges(10,2) * t499 - Icges(10,6) * t714;
t465 = Icges(10,5) * t502 + Icges(10,6) * t501 + Icges(10,3) * t708;
t464 = Icges(10,5) * t500 + Icges(10,6) * t499 - Icges(10,3) * t714;
t459 = 0.1e1 / t460 ^ 2;
t452 = t714 * t803 + V_base(5);
t451 = t495 * t520 + (-t513 - t835) * t691 + t795;
t450 = -t496 * t520 + t514 * t691 + t752;
t448 = -t495 * t514 + t496 * t513 + t747;
t447 = -V_base(5) * pkin(17) + t456 * t538 + (pkin(15) * t708 - t531) * t691 + t795;
t446 = -t457 * t538 + V_base(2) + (-pkin(14) + pkin(17)) * V_base(4) + (-pkin(15) * t714 + t532) * t691;
t444 = cos(t744);
t439 = t452 * t494 + (-t488 - t835) * t691 + t766;
t438 = -t453 * t494 + t489 * t691 + t741;
t437 = pkin(15) * t865 - t456 * t532 + t457 * t531 + V_base(3);
t436 = -pkin(2) * t495 * t755 + t475 * t478 + (-pkin(2) * t540 - t470 - t835) * t691 + t795;
t435 = t471 * t691 - t475 * t479 + (t496 * t755 + t542 * t691) * pkin(2) + t752;
t434 = -t452 * t489 + t453 * t488 + t742;
t432 = sin(t790);
t431 = -t702 * t444 - t698 * t743;
t430 = t444 * t698 - t702 * t743;
t426 = t470 * t479 - t471 * t478 + (-t495 * t542 + t496 * t540) * pkin(2) + t747;
t423 = t666 * t432 + t665 * t773;
t422 = t432 * t665 - t666 * t773;
t421 = qJD(4) * t422 + t691;
t420 = -t649 * t432 + t648 * t773;
t419 = t432 * t648 + t649 * t773;
t418 = -t647 * t432 + t646 * t773;
t417 = t432 * t646 + t647 * t773;
t416 = t420 * t711 + t705 * t708;
t415 = -t420 * t705 + t708 * t711;
t414 = t418 * t711 - t705 * t714;
t413 = -t418 * t705 - t711 * t714;
t412 = t430 * t533 - t431 * t756;
t411 = t430 * t756 + t431 * t533;
t410 = t430 * t523 + t431 * t524;
t409 = -t430 * t524 + t431 * t523;
t408 = t430 * t521 + t431 * t522;
t407 = -t430 * t522 + t431 * t521;
t406 = pkin(10) * t423 + pkin(12) * t422;
t405 = rSges(5,1) * t423 - rSges(5,2) * t422;
t404 = Icges(5,1) * t423 - Icges(5,4) * t422;
t403 = Icges(5,4) * t423 - Icges(5,2) * t422;
t402 = Icges(5,5) * t423 - Icges(5,6) * t422;
t401 = pkin(10) * t420 + pkin(12) * t419;
t400 = pkin(10) * t418 + pkin(12) * t417;
t399 = rSges(5,1) * t420 - rSges(5,2) * t419 + rSges(5,3) * t708;
t398 = rSges(5,1) * t418 - rSges(5,2) * t417 - rSges(5,3) * t714;
t397 = Icges(5,1) * t420 - Icges(5,4) * t419 + Icges(5,5) * t708;
t396 = Icges(5,1) * t418 - Icges(5,4) * t417 - Icges(5,5) * t714;
t395 = Icges(5,4) * t420 - Icges(5,2) * t419 + Icges(5,6) * t708;
t394 = Icges(5,4) * t418 - Icges(5,2) * t417 - Icges(5,6) * t714;
t393 = Icges(5,5) * t420 - Icges(5,6) * t419 + Icges(5,3) * t708;
t392 = Icges(5,5) * t418 - Icges(5,6) * t417 - Icges(5,3) * t714;
t391 = rSges(11,1) * t412 + rSges(11,2) * t411;
t390 = Icges(11,1) * t412 + Icges(11,4) * t411;
t389 = Icges(11,4) * t412 + Icges(11,2) * t411;
t388 = Icges(11,5) * t412 + Icges(11,6) * t411;
t387 = rSges(11,1) * t410 + rSges(11,2) * t409 + rSges(11,3) * t708;
t386 = rSges(11,1) * t408 + rSges(11,2) * t407 - rSges(11,3) * t714;
t385 = Icges(11,1) * t410 + Icges(11,4) * t409 + Icges(11,5) * t708;
t384 = Icges(11,1) * t408 + Icges(11,4) * t407 - Icges(11,5) * t714;
t383 = Icges(11,4) * t410 + Icges(11,2) * t409 + Icges(11,6) * t708;
t382 = Icges(11,4) * t408 + Icges(11,2) * t407 - Icges(11,6) * t714;
t381 = Icges(11,5) * t410 + Icges(11,6) * t409 + Icges(11,3) * t708;
t380 = Icges(11,5) * t408 + Icges(11,6) * t407 - Icges(11,3) * t714;
t379 = 0.2e1 * (((t552 * t786 + (t476 * t545 - t828) * pkin(5)) * t857 + (-t549 * t556 * t727 + t550 * t845) * t846) / t460 - ((-t477 * t545 + t749) * t857 + (t460 * t550 + t549 * t552) * t846) * t459 * t845) * pkin(9) / (t459 * t462 ^ 2 + 0.1e1) * t551 * t722;
t378 = t379 * t708 + t453;
t377 = V_base(5) + (-t379 + t803) * t714;
t376 = rSges(6,3) * t422 + (rSges(6,1) * t711 - rSges(6,2) * t705) * t423;
t375 = Icges(6,5) * t422 + (Icges(6,1) * t711 - Icges(6,4) * t705) * t423;
t374 = Icges(6,6) * t422 + (Icges(6,4) * t711 - Icges(6,2) * t705) * t423;
t373 = Icges(6,3) * t422 + (Icges(6,5) * t711 - Icges(6,6) * t705) * t423;
t372 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t419;
t371 = rSges(6,1) * t414 + rSges(6,2) * t413 + rSges(6,3) * t417;
t370 = Icges(6,1) * t416 + Icges(6,4) * t415 + Icges(6,5) * t419;
t369 = Icges(6,1) * t414 + Icges(6,4) * t413 + Icges(6,5) * t417;
t368 = Icges(6,4) * t416 + Icges(6,2) * t415 + Icges(6,6) * t419;
t367 = Icges(6,4) * t414 + Icges(6,2) * t413 + Icges(6,6) * t417;
t366 = Icges(6,5) * t416 + Icges(6,6) * t415 + Icges(6,3) * t419;
t365 = Icges(6,5) * t414 + Icges(6,6) * t413 + Icges(6,3) * t417;
t360 = qJD(4) * t419 + t362;
t359 = qJD(4) * t417 + t361;
t358 = t377 * t391 + t452 * t490 + (-t386 - t480 - t835) * t691 + t766;
t357 = -pkin(1) * t810 - t378 * t391 - t453 * t490 + (t387 + t481 - t793) * t691 + t752;
t356 = t361 * t405 + (-t398 + t765) * t691 + t753;
t355 = -t362 * t405 + t399 * t691 + t739;
t354 = -t377 * t387 + t378 * t386 - t452 * t481 + t453 * t480 + t742;
t353 = -t361 * t399 + t362 * t398 + t740;
t352 = t359 * t376 + t361 * t406 - t371 * t421 + (-t400 + t765) * t691 + t753;
t351 = -t360 * t376 - t362 * t406 + t372 * t421 + t401 * t691 + t739;
t350 = -t359 * t372 + t360 * t371 - t361 * t401 + t362 * t400 + t740;
t1 = t662 * ((-t606 * t714 - t608 * t647 + t610 * t646) * t663 + (-t605 * t714 - t607 * t647 + t609 * t646) * t662 + (-t614 * t714 - t615 * t647 + t616 * t646) * t691) / 0.2e1 + t663 * ((t606 * t708 - t608 * t649 + t610 * t648) * t663 + (t605 * t708 - t607 * t649 + t609 * t648) * t662 + (t614 * t708 - t615 * t649 + t616 * t648) * t691) / 0.2e1 + t421 * ((t359 * t365 + t360 * t366 + t373 * t421) * t422 + ((-t368 * t705 + t370 * t711) * t360 + (-t367 * t705 + t369 * t711) * t359 + (-t374 * t705 + t375 * t711) * t421) * t423) / 0.2e1 + V_base(5) * t691 * (Icges(2,5) * t708 + Icges(2,6) * t714) + V_base(4) * t691 * (Icges(2,5) * t714 - Icges(2,6) * t708) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t359 * ((t366 * t417 + t368 * t413 + t370 * t414) * t360 + (t365 * t417 + t367 * t413 + t369 * t414) * t359 + (t373 * t417 + t374 * t413 + t375 * t414) * t421) / 0.2e1 + t360 * ((t366 * t419 + t368 * t415 + t370 * t416) * t360 + (t365 * t419 + t367 * t415 + t369 * t416) * t359 + (t373 * t419 + t374 * t415 + t375 * t416) * t421) / 0.2e1 + m(2) * (t618 ^ 2 + t620 ^ 2 + t621 ^ 2) / 0.2e1 + m(8) * (t434 ^ 2 + t438 ^ 2 + t439 ^ 2) / 0.2e1 + m(1) * (t670 ^ 2 + t671 ^ 2 + t672 ^ 2) / 0.2e1 + m(6) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + t479 * ((t465 * t708 + t467 * t501 + t469 * t502) * t479 + (t464 * t708 + t466 * t501 + t468 * t502) * t478 + (t472 * t708 + t473 * t501 + t474 * t502) * t691) / 0.2e1 + t378 * ((t708 * t381 + t409 * t383 + t410 * t385) * t378 + (t380 * t708 + t382 * t409 + t384 * t410) * t377 + (t388 * t708 + t389 * t409 + t390 * t410) * t691) / 0.2e1 + t687 * (t708 * t745 + t714 * t737) / 0.2e1 + t686 * (t708 * t737 - t714 * t745) / 0.2e1 + t457 * (t746 * t708 + t738 * t714) / 0.2e1 + t456 * (t738 * t708 - t746 * t714) / 0.2e1 + ((-t677 * t708 + t680 * t714 + Icges(1,4)) * V_base(5) + (-t678 * t708 + t681 * t714 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t395 * t422 + t397 * t423) * t362 + (-t394 * t422 + t396 * t423) * t361 + (t510 * t547 - t512 * t755) * t496 + (t509 * t547 - t511 * t755) * t495 + (t485 * t533 - t487 * t756) * t453 + (t484 * t533 - t486 * t756) * t452 + (t467 * t503 + t469 * t504) * t479 + (t466 * t503 + t468 * t504) * t478 + (t528 * t561 + t530 * t560) * t457 + (t527 * t561 + t529 * t560) * t456 + (t382 * t411 + t384 * t412) * t377 + (t383 * t411 + t385 * t412) * t378 + (t608 * t666 + t610 * t665) * t663 + (t607 * t666 + t609 * t665) * t662 + (-t638 * t707 + t640 * t713) * t687 + (-t637 * t707 + t639 * t713) * t686 + (Icges(2,3) - t403 * t422 + t404 * t423 + t518 * t547 - t519 * t755 + t492 * t533 - t493 * t756 + t473 * t503 + t474 * t504 + t536 * t561 + t537 * t560 + t389 * t411 + t390 * t412 + t615 * t666 + t616 * t665 - t676 * t707 + t679 * t713) * t691) * t691 / 0.2e1 + m(5) * (t353 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(3) * (t598 ^ 2 + t599 ^ 2 + t600 ^ 2) / 0.2e1 + m(11) * (t354 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + t496 * ((t508 * t708 + t510 * t541 + t512 * t542) * t496 + (t507 * t708 + t509 * t541 + t511 * t542) * t495 + (t517 * t708 + t518 * t541 + t519 * t542) * t691) / 0.2e1 + t362 * ((t708 * t393 - t419 * t395 + t420 * t397) * t362 + (t392 * t708 - t394 * t419 + t396 * t420) * t361 + (t402 * t708 - t403 * t419 + t404 * t420) * t691) / 0.2e1 + t453 * ((t483 * t708 + t485 * t523 + t487 * t524) * t453 + (t482 * t708 + t484 * t523 + t486 * t524) * t452 + (t491 * t708 + t492 * t523 + t493 * t524) * t691) / 0.2e1 + m(9) * (t448 ^ 2 + t450 ^ 2 + t451 ^ 2) / 0.2e1 + m(10) * (t426 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(4) * (t582 ^ 2 + t583 ^ 2 + t584 ^ 2) / 0.2e1 + m(7) * (t437 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + ((t677 * t714 + t680 * t708 + Icges(1,2)) * V_base(5) + (t678 * t714 + t681 * t708 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t361 * ((-t393 * t714 - t395 * t417 + t397 * t418) * t362 + (-t714 * t392 - t417 * t394 + t418 * t396) * t361 + (-t402 * t714 - t403 * t417 + t404 * t418) * t691) / 0.2e1 + t452 * ((-t483 * t714 + t485 * t521 + t487 * t522) * t453 + (-t482 * t714 + t484 * t521 + t486 * t522) * t452 + (-t491 * t714 + t492 * t521 + t493 * t522) * t691) / 0.2e1 + t377 * ((-t381 * t714 + t383 * t407 + t385 * t408) * t378 + (-t714 * t380 + t407 * t382 + t408 * t384) * t377 + (-t388 * t714 + t389 * t407 + t390 * t408) * t691) / 0.2e1 + t495 * ((-t508 * t714 + t510 * t539 + t512 * t540) * t496 + (-t507 * t714 + t509 * t539 + t511 * t540) * t495 + (-t517 * t714 + t518 * t539 + t519 * t540) * t691) / 0.2e1 + t478 * ((-t465 * t714 + t467 * t499 + t469 * t500) * t479 + (-t464 * t714 + t466 * t499 + t468 * t500) * t478 + (-t472 * t714 + t473 * t499 + t474 * t500) * t691) / 0.2e1;
T = t1;
