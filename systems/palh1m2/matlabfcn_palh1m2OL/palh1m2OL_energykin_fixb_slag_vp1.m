% Calculate kinetic energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2OL_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:30
% EndTime: 2020-05-02 21:17:36
% DurationCPUTime: 6.37s
% Computational Cost: add. (2129->464), mult. (2280->770), div. (0->0), fcn. (2036->27), ass. (0->276)
t691 = sin(qJ(1));
t797 = t691 ^ 2;
t696 = cos(qJ(1));
t796 = t696 ^ 2;
t684 = qJ(2) + qJ(8);
t674 = sin(t684);
t794 = pkin(2) * t674;
t792 = pkin(1) * qJD(2);
t690 = sin(qJ(2));
t791 = Icges(3,4) * t690;
t695 = cos(qJ(2));
t790 = Icges(3,4) * t695;
t686 = qJ(2) + qJ(3);
t676 = sin(t686);
t789 = Icges(4,4) * t676;
t679 = cos(t686);
t788 = Icges(4,4) * t679;
t681 = qJ(4) + t686;
t659 = sin(t681);
t787 = Icges(5,4) * t659;
t661 = cos(t681);
t786 = Icges(5,4) * t661;
t687 = sin(qJ(6));
t785 = Icges(7,4) * t687;
t692 = cos(qJ(6));
t784 = Icges(7,4) * t692;
t685 = qJ(2) + qJ(7);
t675 = sin(t685);
t783 = Icges(8,4) * t675;
t678 = cos(t685);
t782 = Icges(8,4) * t678;
t781 = Icges(9,4) * t674;
t677 = cos(t684);
t780 = Icges(9,4) * t677;
t680 = qJ(9) + t684;
t658 = sin(t680);
t779 = Icges(10,4) * t658;
t660 = cos(t680);
t778 = Icges(10,4) * t660;
t689 = sin(qJ(3));
t662 = pkin(5) * t689 + pkin(1);
t694 = cos(qJ(3));
t594 = pkin(5) * t694 * t695 - t662 * t690 + pkin(15);
t777 = t594 * t691;
t657 = pkin(1) * t690 - pkin(15);
t776 = t657 * t696;
t775 = t659 * t691;
t774 = t659 * t696;
t688 = sin(qJ(5));
t773 = t688 * t691;
t772 = t688 * t696;
t771 = t690 * t694;
t693 = cos(qJ(5));
t770 = t691 * t693;
t769 = t693 * t696;
t673 = qJD(2) * t691;
t627 = qJD(8) * t691 + t673;
t628 = qJD(7) * t691 + t673;
t629 = qJD(3) * t691 + t673;
t683 = -qJ(7) + pkin(19);
t668 = -qJ(10) + t683;
t656 = -qJ(2) + t668;
t652 = sin(t656);
t768 = Icges(11,4) * t652;
t653 = cos(t656);
t767 = Icges(11,4) * t653;
t766 = qJD(1) * t691;
t765 = qJD(1) * t696;
t764 = qJD(2) * t696;
t763 = qJD(5) * t659;
t762 = qJD(6) * t691;
t761 = qJD(6) * t696;
t760 = -qJD(2) - qJD(7);
t682 = -qJD(2) - qJD(8);
t759 = qJD(2) + qJD(3);
t758 = t695 * t792;
t757 = t690 * t792;
t610 = qJD(4) * t691 + t629;
t756 = -qJD(9) + t682;
t755 = -qJD(10) + t760;
t754 = t691 * t758;
t753 = pkin(5) * t759 * t679 - t757;
t752 = pkin(9) * t661 + pkin(11) * t659;
t644 = rSges(3,1) * t690 + rSges(3,2) * t695;
t751 = rSges(4,1) * t679 - rSges(4,2) * t676;
t750 = rSges(5,1) * t661 - rSges(5,2) * t659;
t648 = rSges(6,1) * t693 - rSges(6,2) * t688;
t646 = rSges(7,1) * t692 - rSges(7,2) * t687;
t749 = -rSges(8,1) * t675 - rSges(8,2) * t678;
t748 = rSges(9,1) * t674 + rSges(9,2) * t677;
t747 = rSges(10,1) * t658 + rSges(10,2) * t660;
t612 = (-qJD(4) - t759) * t696;
t746 = rSges(6,3) * t659 + t648 * t661;
t745 = -Icges(3,1) * t690 - t790;
t744 = Icges(4,1) * t679 - t789;
t743 = Icges(5,1) * t661 - t787;
t742 = Icges(7,1) * t692 - t785;
t741 = -Icges(8,1) * t675 - t782;
t740 = -Icges(9,1) * t674 - t780;
t739 = Icges(10,1) * t658 + t778;
t738 = -Icges(3,2) * t695 - t791;
t737 = -Icges(4,2) * t676 + t788;
t736 = -Icges(5,2) * t659 + t786;
t735 = -Icges(7,2) * t687 + t784;
t734 = -Icges(8,2) * t678 - t783;
t733 = -Icges(9,2) * t677 - t781;
t732 = Icges(10,2) * t660 + t779;
t731 = -Icges(3,5) * t690 - Icges(3,6) * t695;
t730 = Icges(4,5) * t679 - Icges(4,6) * t676;
t729 = Icges(5,5) * t661 - Icges(5,6) * t659;
t728 = Icges(7,5) * t692 - Icges(7,6) * t687;
t727 = -Icges(8,5) * t675 - Icges(8,6) * t678;
t726 = -Icges(9,5) * t674 - Icges(9,6) * t677;
t725 = Icges(10,5) * t658 + Icges(10,6) * t660;
t582 = -Icges(7,6) * t696 + t691 * t735;
t586 = -Icges(7,5) * t696 + t691 * t742;
t724 = t582 * t687 - t586 * t692;
t583 = Icges(7,6) * t691 + t696 * t735;
t587 = Icges(7,5) * t691 + t696 * t742;
t723 = -t583 * t687 + t587 * t692;
t584 = -Icges(3,6) * t696 + t691 * t738;
t588 = -Icges(3,5) * t696 + t691 * t745;
t722 = t584 * t695 + t588 * t690;
t585 = Icges(3,6) * t691 + t696 * t738;
t589 = Icges(3,5) * t691 + t696 * t745;
t721 = -t585 * t695 - t589 * t690;
t634 = rSges(11,1) * t690 + rSges(11,2) * t695;
t635 = rSges(11,1) * t695 - rSges(11,2) * t690;
t654 = sin(t668);
t655 = cos(t668);
t720 = -t634 * t655 + t635 * t654;
t638 = Icges(7,2) * t692 + t785;
t640 = Icges(7,1) * t687 + t784;
t719 = -t638 * t687 + t640 * t692;
t639 = -Icges(3,2) * t690 + t790;
t641 = Icges(3,1) * t695 - t791;
t718 = -t639 * t695 - t641 * t690;
t717 = -Icges(11,1) * t652 + t767;
t716 = Icges(11,2) * t653 - t768;
t715 = -Icges(11,5) * t652 + Icges(11,6) * t653;
t714 = t657 * t766 - t696 * t758;
t713 = -pkin(14) + t646;
t712 = -pkin(5) * qJD(3) * (t689 * t695 + t771) - qJD(2) * (pkin(5) * t771 + t662 * t695);
t711 = t712 * t696;
t607 = qJD(10) * t691 + t628;
t608 = t755 * t696;
t710 = (-Icges(11,5) * t653 - Icges(11,6) * t652) * qJD(1) + (-Icges(11,3) * t696 + t691 * t715) * t608 + (Icges(11,3) * t691 + t696 * t715) * t607;
t609 = qJD(9) * t691 + t627;
t611 = t756 * t696;
t709 = (-Icges(10,5) * t660 + Icges(10,6) * t658) * qJD(1) + (-Icges(10,3) * t696 + t691 * t725) * t611 + (Icges(10,3) * t691 + t696 * t725) * t609;
t708 = (Icges(5,5) * t659 + Icges(5,6) * t661) * qJD(1) + (-Icges(5,3) * t696 + t691 * t729) * t612 + (Icges(5,3) * t691 + t696 * t729) * t610;
t630 = t682 * t696;
t707 = (Icges(9,5) * t677 - Icges(9,6) * t674) * qJD(1) + (-Icges(9,3) * t696 + t691 * t726) * t630 + (Icges(9,3) * t691 + t696 * t726) * t627;
t631 = t760 * t696;
t706 = (Icges(8,5) * t678 - Icges(8,6) * t675) * qJD(1) + (-Icges(8,3) * t696 + t691 * t727) * t631 + (Icges(8,3) * t691 + t696 * t727) * t628;
t632 = t759 * t696;
t705 = (Icges(4,5) * t676 + Icges(4,6) * t679) * qJD(1) - (-Icges(4,3) * t696 + t691 * t730) * t632 + (Icges(4,3) * t691 + t696 * t730) * t629;
t704 = t594 * t765 + t691 * t712;
t525 = -Icges(11,6) * t696 + t691 * t716;
t526 = Icges(11,6) * t691 + t696 * t716;
t527 = -Icges(11,5) * t696 + t691 * t717;
t528 = Icges(11,5) * t691 + t696 * t717;
t576 = -Icges(11,2) * t652 - t767;
t577 = -Icges(11,1) * t653 - t768;
t703 = (t526 * t653 - t528 * t652) * t607 + (t525 * t653 - t527 * t652) * t608 + (t576 * t653 - t577 * t652) * qJD(1);
t537 = -Icges(10,6) * t696 + t691 * t732;
t538 = Icges(10,6) * t691 + t696 * t732;
t541 = -Icges(10,5) * t696 + t691 * t739;
t542 = Icges(10,5) * t691 + t696 * t739;
t600 = Icges(10,2) * t658 - t778;
t602 = -Icges(10,1) * t660 + t779;
t702 = (t538 * t660 + t542 * t658) * t609 + (t537 * t660 + t541 * t658) * t611 + (t600 * t660 + t602 * t658) * qJD(1);
t539 = -Icges(5,6) * t696 + t691 * t736;
t540 = Icges(5,6) * t691 + t696 * t736;
t543 = -Icges(5,5) * t696 + t691 * t743;
t544 = Icges(5,5) * t691 + t696 * t743;
t601 = Icges(5,2) * t661 + t787;
t603 = Icges(5,1) * t659 + t786;
t701 = (-t540 * t659 + t544 * t661) * t610 + (-t539 * t659 + t543 * t661) * t612 + (-t601 * t659 + t603 * t661) * qJD(1);
t555 = -Icges(9,6) * t696 + t691 * t733;
t556 = Icges(9,6) * t691 + t696 * t733;
t561 = -Icges(9,5) * t696 + t691 * t740;
t562 = Icges(9,5) * t691 + t696 * t740;
t616 = -Icges(9,2) * t674 + t780;
t619 = Icges(9,1) * t677 - t781;
t700 = (-t556 * t677 - t562 * t674) * t627 + (-t555 * t677 - t561 * t674) * t630 + (-t616 * t677 - t619 * t674) * qJD(1);
t557 = -Icges(8,6) * t696 + t691 * t734;
t558 = Icges(8,6) * t691 + t696 * t734;
t563 = -Icges(8,5) * t696 + t691 * t741;
t564 = Icges(8,5) * t691 + t696 * t741;
t617 = -Icges(8,2) * t675 + t782;
t620 = Icges(8,1) * t678 - t783;
t699 = (-t558 * t678 - t564 * t675) * t628 + (-t557 * t678 - t563 * t675) * t631 + (-t617 * t678 - t620 * t675) * qJD(1);
t559 = -Icges(4,6) * t696 + t691 * t737;
t560 = Icges(4,6) * t691 + t696 * t737;
t565 = -Icges(4,5) * t696 + t691 * t744;
t566 = Icges(4,5) * t691 + t696 * t744;
t618 = Icges(4,2) * t679 + t789;
t621 = Icges(4,1) * t676 + t788;
t698 = (-t560 * t676 + t566 * t679) * t629 - (-t559 * t676 + t565 * t679) * t632 + (-t618 * t676 + t621 * t679) * qJD(1);
t666 = cos(t683);
t665 = sin(t683);
t663 = pkin(15) * t765;
t651 = rSges(2,1) * t696 - rSges(2,2) * t691;
t650 = rSges(3,1) * t695 - rSges(3,2) * t690;
t647 = rSges(6,1) * t688 + rSges(6,2) * t693;
t645 = rSges(2,1) * t691 + rSges(2,2) * t696;
t643 = rSges(7,1) * t687 + rSges(7,2) * t692;
t637 = Icges(3,5) * t695 - Icges(3,6) * t690;
t636 = Icges(7,5) * t687 + Icges(7,6) * t692;
t633 = -qJD(5) * t661 + qJD(1);
t624 = rSges(8,1) * t678 - rSges(8,2) * t675;
t623 = rSges(9,1) * t677 - rSges(9,2) * t674;
t622 = rSges(4,1) * t676 + rSges(4,2) * t679;
t606 = pkin(9) * t659 - pkin(11) * t661;
t605 = rSges(10,1) * t660 - rSges(10,2) * t658;
t604 = rSges(5,1) * t659 + rSges(5,2) * t661;
t596 = t665 * t690 + t666 * t695;
t595 = t665 * t695 - t666 * t690;
t593 = t661 * t769 + t773;
t592 = -t661 * t772 + t770;
t591 = t661 * t770 - t772;
t590 = -t661 * t773 - t769;
t581 = Icges(3,3) * t691 + t696 * t731;
t580 = -Icges(3,3) * t696 + t691 * t731;
t579 = Icges(7,3) * t691 + t696 * t728;
t578 = -Icges(7,3) * t696 + t691 * t728;
t574 = t748 * t682;
t573 = t752 * t696;
t572 = t752 * t691;
t570 = rSges(4,3) * t691 + t696 * t751;
t569 = rSges(8,3) * t691 + t696 * t749;
t568 = -rSges(4,3) * t696 + t691 * t751;
t567 = -rSges(8,3) * t696 + t691 * t749;
t548 = rSges(5,3) * t691 + t696 * t750;
t547 = -rSges(5,3) * t696 + t691 * t750;
t546 = t691 * t763 + t612;
t545 = t696 * t763 + t610;
t532 = -rSges(6,3) * t661 + t648 * t659;
t531 = -Icges(6,5) * t661 + (Icges(6,1) * t693 - Icges(6,4) * t688) * t659;
t530 = -Icges(6,6) * t661 + (Icges(6,4) * t693 - Icges(6,2) * t688) * t659;
t529 = -Icges(6,3) * t661 + (Icges(6,5) * t693 - Icges(6,6) * t688) * t659;
t522 = t634 * t654 + t635 * t655;
t521 = t682 * t794 - t747 * t756;
t520 = -t650 * t764 + (t696 * rSges(3,3) + (-pkin(15) + t644) * t691) * qJD(1);
t519 = -t643 * t762 + (t691 * rSges(7,3) + t696 * t713) * qJD(1);
t518 = -t643 * t761 + (t696 * rSges(7,3) - t691 * t713) * qJD(1);
t517 = t663 - t650 * t673 + qJD(1) * (rSges(3,3) * t691 - t644 * t696);
t516 = t647 * t691 + t696 * t746;
t515 = -t647 * t696 + t691 * t746;
t514 = pkin(4) * t760 * sin(qJ(2) - t683) - t757 + (rSges(11,1) * t652 - rSges(11,2) * t653) * t755;
t513 = Icges(6,1) * t593 + Icges(6,4) * t592 + Icges(6,5) * t774;
t512 = Icges(6,1) * t591 + Icges(6,4) * t590 + Icges(6,5) * t775;
t511 = Icges(6,4) * t593 + Icges(6,2) * t592 + Icges(6,6) * t774;
t510 = Icges(6,4) * t591 + Icges(6,2) * t590 + Icges(6,6) * t775;
t509 = Icges(6,5) * t593 + Icges(6,6) * t592 + Icges(6,3) * t774;
t508 = Icges(6,5) * t591 + Icges(6,6) * t590 + Icges(6,3) * t775;
t507 = t663 - t627 * t623 + qJD(1) * (rSges(9,3) * t691 - t696 * t748);
t506 = t630 * t623 + (t696 * rSges(9,3) + (-pkin(15) + t748) * t691) * qJD(1);
t505 = -t754 - t622 * t629 + (t570 - t776) * qJD(1);
t504 = -t754 - t624 * t628 + (t569 - t776) * qJD(1);
t503 = -qJD(1) * t568 - t622 * t632 + t714;
t502 = -qJD(1) * t567 + t624 * t631 + t714;
t501 = t568 * t629 + t570 * t632 - t757;
t500 = t567 * t628 - t569 * t631 - t757;
t499 = t630 * pkin(2) * t677 - t611 * t605 + (t696 * rSges(10,3) + (-pkin(15) - t747 + t794) * t691) * qJD(1);
t498 = t663 + t609 * t605 + qJD(1) * (rSges(10,3) * t691 + t696 * t747) + (-t627 * t677 - t674 * t765) * pkin(2);
t497 = t547 * t610 - t548 * t612 + t753;
t496 = qJD(1) * t548 - t604 * t610 + t704;
t495 = t604 * t612 + t711 + (-t547 - t777) * qJD(1);
t494 = -t754 - t628 * pkin(4) * t596 + t522 * t607 + (t691 * rSges(11,3) + (pkin(4) * t595 - t657 - t720) * t696) * qJD(1);
t493 = -t522 * t608 + qJD(1) * (rSges(11,3) * t696 + t691 * t720) + (-t595 * t766 + t596 * t631) * pkin(4) + t714;
t492 = qJD(1) * t573 + t516 * t633 - t532 * t545 - t606 * t610 + t704;
t491 = -t515 * t633 + t532 * t546 + t606 * t612 + t711 + (-t572 - t777) * qJD(1);
t490 = t515 * t545 - t516 * t546 + t572 * t610 - t573 * t612 + t753;
t1 = -((-t696 * t636 + t691 * t719) * qJD(1) + (t796 * t578 + (t723 * t691 + (-t579 + t724) * t696) * t691) * qJD(6)) * t761 / 0.2e1 - ((-t696 * t637 + t691 * t718) * qJD(1) + (t796 * t580 + (t721 * t691 + (-t581 + t722) * t696) * t691) * qJD(2)) * t764 / 0.2e1 + ((t691 * t636 + t696 * t719) * qJD(1) + (t797 * t579 + (t724 * t696 + (-t578 + t723) * t691) * t696) * qJD(6)) * t762 / 0.2e1 + ((t691 * t637 + t696 * t718) * qJD(1) + (t797 * t581 + (t722 * t696 + (-t580 + t721) * t691) * t696) * qJD(2)) * t673 / 0.2e1 + (Icges(2,3) + m(2) * (t645 ^ 2 + t651 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t585 * t690 + t589 * t695) * t691 - (-t584 * t690 + t588 * t695) * t696) * qJD(2) + ((t583 * t692 + t587 * t687) * t691 - (t582 * t692 + t586 * t687) * t696) * qJD(6) + (t560 * t679 + t566 * t676) * t629 - (t559 * t679 + t565 * t676) * t632 + (-t558 * t675 + t564 * t678) * t628 + (-t557 * t675 + t563 * t678) * t631 + (-t556 * t674 + t562 * t677) * t627 + (-t555 * t674 + t561 * t677) * t630 + (t540 * t661 + t544 * t659) * t610 + (t539 * t661 + t543 * t659) * t612 + (t538 * t658 - t542 * t660) * t609 + (t537 * t658 - t541 * t660) * t611 + (-t526 * t652 - t528 * t653) * t607 + (-t525 * t652 - t527 * t653) * t608 + (-t639 * t690 + t641 * t695 + t638 * t692 + t640 * t687 + t679 * t618 + t676 * t621 - t675 * t617 + t678 * t620 - t674 * t616 + t677 * t619 + t661 * t601 + t659 * t603 + t658 * t600 - t660 * t602 - t652 * t576 - t653 * t577) * qJD(1)) * qJD(1) / 0.2e1 + m(4) * (t501 ^ 2 + t503 ^ 2 + t505 ^ 2) / 0.2e1 + t633 * ((-t508 * t546 - t509 * t545 - t529 * t633) * t661 + ((-t511 * t688 + t513 * t693) * t545 + (-t510 * t688 + t512 * t693) * t546 + (-t530 * t688 + t531 * t693) * t633) * t659) / 0.2e1 + m(7) * (qJD(6) ^ 2 * t646 ^ 2 + t518 ^ 2 + t519 ^ 2) / 0.2e1 + m(8) * (t500 ^ 2 + t502 ^ 2 + t504 ^ 2) / 0.2e1 + m(11) * (t493 ^ 2 + t494 ^ 2 + t514 ^ 2) / 0.2e1 + m(10) * (t498 ^ 2 + t499 ^ 2 + t521 ^ 2) / 0.2e1 + t629 * (t691 * t705 + t696 * t698) / 0.2e1 - t632 * (t691 * t698 - t696 * t705) / 0.2e1 + t628 * (t691 * t706 + t696 * t699) / 0.2e1 + t631 * (t691 * t699 - t696 * t706) / 0.2e1 + t627 * (t691 * t707 + t696 * t700) / 0.2e1 + t630 * (t691 * t700 - t696 * t707) / 0.2e1 + t610 * (t691 * t708 + t696 * t701) / 0.2e1 + t612 * (t691 * t701 - t696 * t708) / 0.2e1 + t609 * (t691 * t709 + t696 * t702) / 0.2e1 + t611 * (t691 * t702 - t696 * t709) / 0.2e1 + t607 * (t691 * t710 + t696 * t703) / 0.2e1 + t608 * (t691 * t703 - t696 * t710) / 0.2e1 + t545 * ((t509 * t774 + t592 * t511 + t593 * t513) * t545 + (t508 * t774 + t510 * t592 + t512 * t593) * t546 + (t529 * t774 + t530 * t592 + t531 * t593) * t633) / 0.2e1 + t546 * ((t509 * t775 + t511 * t590 + t513 * t591) * t545 + (t508 * t775 + t590 * t510 + t591 * t512) * t546 + (t529 * t775 + t530 * t590 + t531 * t591) * t633) / 0.2e1 + m(5) * (t495 ^ 2 + t496 ^ 2 + t497 ^ 2) / 0.2e1 + m(6) * (t490 ^ 2 + t491 ^ 2 + t492 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t644 ^ 2 + t517 ^ 2 + t520 ^ 2) / 0.2e1 + m(9) * (t506 ^ 2 + t507 ^ 2 + t574 ^ 2) / 0.2e1;
T = t1;
