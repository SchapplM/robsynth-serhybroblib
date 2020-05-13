% Calculate kinetic energy for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1OL_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:27
% EndTime: 2020-04-15 19:28:35
% DurationCPUTime: 6.95s
% Computational Cost: add. (2275->465), mult. (2447->758), div. (0->0), fcn. (2180->22), ass. (0->278)
t660 = sin(qJ(1));
t766 = t660 ^ 2;
t664 = cos(qJ(1));
t765 = t664 ^ 2;
t654 = qJ(2) + qJ(8);
t649 = cos(t654);
t763 = pkin(2) * t649;
t655 = qJ(2) + qJ(7);
t640 = pkin(19) - t655;
t762 = pkin(4) * cos(t640);
t656 = qJ(2) + qJ(3);
t648 = sin(t656);
t761 = pkin(5) * t648;
t759 = pkin(15) * t660;
t646 = sin(t654);
t758 = pkin(2) * t646;
t659 = sin(qJ(2));
t757 = t659 * pkin(1);
t756 = Icges(3,4) * t659;
t663 = cos(qJ(2));
t755 = Icges(3,4) * t663;
t754 = Icges(4,4) * t648;
t651 = cos(t656);
t753 = Icges(4,4) * t651;
t653 = qJ(4) + t656;
t636 = sin(t653);
t752 = Icges(5,4) * t636;
t638 = cos(t653);
t751 = Icges(5,4) * t638;
t657 = sin(qJ(6));
t750 = Icges(7,4) * t657;
t661 = cos(qJ(6));
t749 = Icges(7,4) * t661;
t647 = sin(t655);
t748 = Icges(8,4) * t647;
t650 = cos(t655);
t747 = Icges(8,4) * t650;
t746 = Icges(9,4) * t646;
t745 = Icges(9,4) * t649;
t652 = qJ(9) + t654;
t635 = sin(t652);
t744 = Icges(10,4) * t635;
t637 = cos(t652);
t743 = Icges(10,4) * t637;
t742 = t636 * t660;
t741 = t636 * t664;
t658 = sin(qJ(5));
t740 = t658 * t660;
t739 = t658 * t664;
t662 = cos(qJ(5));
t738 = t660 * t662;
t737 = t662 * t664;
t596 = t757 * t660;
t597 = t757 * t664;
t645 = qJD(2) * t660;
t731 = qJD(2) * t664;
t736 = -t596 * t645 - t597 * t731;
t735 = pkin(4) * sin(t640);
t734 = pkin(5) * t651;
t612 = qJD(8) * t660 + t645;
t613 = qJD(7) * t660 + t645;
t614 = qJD(3) * t660 + t645;
t633 = -qJ(10) + t640;
t630 = sin(t633);
t733 = Icges(11,4) * t630;
t631 = cos(t633);
t732 = Icges(11,4) * t631;
t730 = qJD(5) * t636;
t729 = qJD(6) * t660;
t728 = qJD(6) * t664;
t727 = -qJD(2) - qJD(3);
t726 = -qJD(2) - qJD(7);
t725 = -qJD(2) - qJD(8);
t724 = pkin(1) * qJD(2) * t663;
t593 = qJD(4) * t660 + t614;
t723 = t596 - t759;
t722 = t664 * t724;
t523 = t734 * t660;
t721 = -t523 + t723;
t720 = pkin(9) * t638 + pkin(11) * t636;
t719 = -rSges(3,1) * t659 - rSges(3,2) * t663;
t718 = rSges(4,1) * t651 - rSges(4,2) * t648;
t717 = rSges(5,1) * t638 - rSges(5,2) * t636;
t716 = rSges(7,1) * t661 - rSges(7,2) * t657;
t715 = -rSges(8,1) * t647 - rSges(8,2) * t650;
t714 = -rSges(9,1) * t646 - rSges(9,2) * t649;
t713 = rSges(10,1) * t635 + rSges(10,2) * t637;
t595 = (-qJD(4) + t727) * t664;
t712 = -rSges(11,1) * t630 + rSges(11,2) * t631;
t711 = -Icges(3,1) * t659 - t755;
t710 = Icges(4,1) * t651 - t754;
t709 = Icges(5,1) * t638 - t752;
t708 = Icges(7,1) * t661 - t750;
t707 = -Icges(8,1) * t647 - t747;
t706 = -Icges(9,1) * t646 - t745;
t705 = Icges(10,1) * t635 + t743;
t704 = -Icges(3,2) * t663 - t756;
t703 = -Icges(4,2) * t648 + t753;
t702 = -Icges(5,2) * t636 + t751;
t701 = -Icges(7,2) * t657 + t749;
t700 = -Icges(8,2) * t650 - t748;
t699 = -Icges(9,2) * t649 - t746;
t698 = Icges(10,2) * t637 + t744;
t697 = -Icges(3,5) * t659 - Icges(3,6) * t663;
t696 = Icges(4,5) * t651 - Icges(4,6) * t648;
t695 = Icges(5,5) * t638 - Icges(5,6) * t636;
t694 = Icges(7,5) * t661 - Icges(7,6) * t657;
t693 = -Icges(8,5) * t647 - Icges(8,6) * t650;
t692 = -Icges(9,5) * t646 - Icges(9,6) * t649;
t691 = Icges(10,5) * t635 + Icges(10,6) * t637;
t558 = -Icges(7,6) * t664 + t660 * t701;
t562 = -Icges(7,5) * t664 + t660 * t708;
t690 = t558 * t657 - t562 * t661;
t559 = Icges(7,6) * t660 + t664 * t701;
t563 = Icges(7,5) * t660 + t664 * t708;
t689 = -t559 * t657 + t563 * t661;
t560 = -Icges(3,6) * t664 + t660 * t704;
t564 = -Icges(3,5) * t664 + t660 * t711;
t688 = t560 * t663 + t564 * t659;
t561 = Icges(3,6) * t660 + t664 * t704;
t565 = Icges(3,5) * t660 + t664 * t711;
t687 = -t561 * t663 - t565 * t659;
t621 = Icges(7,2) * t661 + t750;
t623 = Icges(7,1) * t657 + t749;
t686 = -t621 * t657 + t623 * t661;
t622 = -Icges(3,2) * t659 + t755;
t624 = Icges(3,1) * t663 - t756;
t685 = -t622 * t663 - t624 * t659;
t684 = -Icges(11,1) * t630 + t732;
t683 = Icges(11,2) * t631 - t733;
t682 = -Icges(11,5) * t630 + Icges(11,6) * t631;
t617 = t727 * t664;
t681 = t617 * t761 - t722;
t524 = t734 * t664;
t680 = t614 * t523 - t524 * t617 + t736;
t639 = qJD(1) * t664 * pkin(15);
t679 = -qJD(1) * t597 - t660 * t724 + t639;
t590 = qJD(10) * t660 + t613;
t591 = (-qJD(10) + t726) * t664;
t678 = qJD(1) * (-Icges(11,5) * t631 - Icges(11,6) * t630) + (-Icges(11,3) * t664 + t660 * t682) * t591 + (Icges(11,3) * t660 + t664 * t682) * t590;
t592 = qJD(9) * t660 + t612;
t594 = (-qJD(9) + t725) * t664;
t677 = qJD(1) * (-Icges(10,5) * t637 + Icges(10,6) * t635) + (-Icges(10,3) * t664 + t660 * t691) * t594 + (Icges(10,3) * t660 + t664 * t691) * t592;
t676 = qJD(1) * (Icges(5,5) * t636 + Icges(5,6) * t638) + (-Icges(5,3) * t664 + t660 * t695) * t595 + (Icges(5,3) * t660 + t664 * t695) * t593;
t615 = t725 * t664;
t675 = qJD(1) * (Icges(9,5) * t649 - Icges(9,6) * t646) + (-Icges(9,3) * t664 + t660 * t692) * t615 + (Icges(9,3) * t660 + t664 * t692) * t612;
t616 = t726 * t664;
t674 = qJD(1) * (Icges(8,5) * t650 - Icges(8,6) * t647) + (-Icges(8,3) * t664 + t660 * t693) * t616 + (Icges(8,3) * t660 + t664 * t693) * t613;
t673 = qJD(1) * (Icges(4,5) * t648 + Icges(4,6) * t651) + (-Icges(4,3) * t664 + t660 * t696) * t617 + (Icges(4,3) * t660 + t664 * t696) * t614;
t672 = qJD(1) * t524 - t614 * t761 + t679;
t492 = -Icges(11,6) * t664 + t660 * t683;
t493 = Icges(11,6) * t660 + t664 * t683;
t494 = -Icges(11,5) * t664 + t660 * t684;
t495 = Icges(11,5) * t660 + t664 * t684;
t552 = -Icges(11,2) * t630 - t732;
t553 = -Icges(11,1) * t631 - t733;
t671 = (t493 * t631 - t495 * t630) * t590 + (t492 * t631 - t494 * t630) * t591 + (t552 * t631 - t553 * t630) * qJD(1);
t509 = -Icges(10,6) * t664 + t660 * t698;
t510 = Icges(10,6) * t660 + t664 * t698;
t513 = -Icges(10,5) * t664 + t660 * t705;
t514 = Icges(10,5) * t660 + t664 * t705;
t583 = Icges(10,2) * t635 - t743;
t585 = -Icges(10,1) * t637 + t744;
t670 = (t510 * t637 + t514 * t635) * t592 + (t509 * t637 + t513 * t635) * t594 + (t583 * t637 + t585 * t635) * qJD(1);
t511 = -Icges(5,6) * t664 + t660 * t702;
t512 = Icges(5,6) * t660 + t664 * t702;
t515 = -Icges(5,5) * t664 + t660 * t709;
t516 = Icges(5,5) * t660 + t664 * t709;
t584 = Icges(5,2) * t638 + t752;
t586 = Icges(5,1) * t636 + t751;
t669 = (-t512 * t636 + t516 * t638) * t593 + (-t511 * t636 + t515 * t638) * t595 + (-t584 * t636 + t586 * t638) * qJD(1);
t531 = -Icges(9,6) * t664 + t660 * t699;
t532 = Icges(9,6) * t660 + t664 * t699;
t537 = -Icges(9,5) * t664 + t660 * t706;
t538 = Icges(9,5) * t660 + t664 * t706;
t602 = -Icges(9,2) * t646 + t745;
t605 = Icges(9,1) * t649 - t746;
t668 = (-t532 * t649 - t538 * t646) * t612 + (-t531 * t649 - t537 * t646) * t615 + (-t602 * t649 - t605 * t646) * qJD(1);
t533 = -Icges(8,6) * t664 + t660 * t700;
t534 = Icges(8,6) * t660 + t664 * t700;
t539 = -Icges(8,5) * t664 + t660 * t707;
t540 = Icges(8,5) * t660 + t664 * t707;
t603 = -Icges(8,2) * t647 + t747;
t606 = Icges(8,1) * t650 - t748;
t667 = (-t534 * t650 - t540 * t647) * t613 + (-t533 * t650 - t539 * t647) * t616 + (-t603 * t650 - t606 * t647) * qJD(1);
t535 = -Icges(4,6) * t664 + t660 * t703;
t536 = Icges(4,6) * t660 + t664 * t703;
t541 = -Icges(4,5) * t664 + t660 * t710;
t542 = Icges(4,5) * t660 + t664 * t710;
t604 = Icges(4,2) * t651 + t754;
t607 = Icges(4,1) * t648 + t753;
t666 = (-t536 * t648 + t542 * t651) * t614 + (-t535 * t648 + t541 * t651) * t617 + (-t604 * t648 + t607 * t651) * qJD(1);
t628 = rSges(2,1) * t664 - rSges(2,2) * t660;
t627 = rSges(3,1) * t663 - rSges(3,2) * t659;
t626 = rSges(2,1) * t660 + rSges(2,2) * t664;
t625 = rSges(7,1) * t657 + rSges(7,2) * t661;
t620 = Icges(3,5) * t663 - Icges(3,6) * t659;
t619 = Icges(7,5) * t657 + Icges(7,6) * t661;
t618 = -qJD(5) * t638 + qJD(1);
t610 = rSges(8,1) * t650 - rSges(8,2) * t647;
t609 = rSges(9,1) * t649 - rSges(9,2) * t646;
t608 = rSges(4,1) * t648 + rSges(4,2) * t651;
t589 = pkin(9) * t636 - pkin(11) * t638;
t588 = -rSges(10,1) * t637 + rSges(10,2) * t635;
t587 = rSges(5,1) * t636 + rSges(5,2) * t638;
t579 = t758 * t664;
t578 = t758 * t660;
t576 = t638 * t737 + t740;
t575 = -t638 * t739 + t738;
t574 = t638 * t738 - t739;
t573 = -t638 * t740 - t737;
t570 = rSges(7,3) * t660 + t664 * t716;
t569 = rSges(3,3) * t660 + t664 * t719;
t568 = -rSges(7,3) * t664 + t660 * t716;
t567 = -rSges(3,3) * t664 + t660 * t719;
t566 = -rSges(11,1) * t631 - rSges(11,2) * t630;
t557 = Icges(3,3) * t660 + t664 * t697;
t556 = -Icges(3,3) * t664 + t660 * t697;
t555 = Icges(7,3) * t660 + t664 * t694;
t554 = -Icges(7,3) * t664 + t660 * t694;
t550 = t720 * t664;
t549 = t720 * t660;
t548 = rSges(4,3) * t660 + t664 * t718;
t547 = rSges(8,3) * t660 + t664 * t715;
t546 = rSges(9,3) * t660 + t664 * t714;
t545 = -rSges(4,3) * t664 + t660 * t718;
t544 = -rSges(8,3) * t664 + t660 * t715;
t543 = -rSges(9,3) * t664 + t660 * t714;
t522 = rSges(5,3) * t660 + t664 * t717;
t521 = rSges(10,3) * t660 + t664 * t713;
t520 = -rSges(5,3) * t664 + t660 * t717;
t519 = -rSges(10,3) * t664 + t660 * t713;
t518 = t660 * t730 + t595;
t517 = t664 * t730 + t593;
t504 = t735 * t664;
t503 = t735 * t660;
t501 = -rSges(6,3) * t638 + (rSges(6,1) * t662 - rSges(6,2) * t658) * t636;
t500 = -Icges(6,5) * t638 + (Icges(6,1) * t662 - Icges(6,4) * t658) * t636;
t499 = -Icges(6,6) * t638 + (Icges(6,4) * t662 - Icges(6,2) * t658) * t636;
t498 = -Icges(6,3) * t638 + (Icges(6,5) * t662 - Icges(6,6) * t658) * t636;
t497 = rSges(11,3) * t660 + t664 * t712;
t496 = -rSges(11,3) * t664 + t660 * t712;
t488 = -t625 * t729 + (-pkin(14) * t664 + t570) * qJD(1);
t487 = qJD(1) * t569 - t627 * t645 + t639;
t486 = -t625 * t728 + (pkin(14) * t660 - t568) * qJD(1);
t485 = -t627 * t731 + (-t567 - t759) * qJD(1);
t484 = (t567 * t660 + t569 * t664) * qJD(2);
t483 = (t568 * t660 + t570 * t664) * qJD(6);
t482 = rSges(6,1) * t576 + rSges(6,2) * t575 + rSges(6,3) * t741;
t481 = rSges(6,1) * t574 + rSges(6,2) * t573 + rSges(6,3) * t742;
t480 = Icges(6,1) * t576 + Icges(6,4) * t575 + Icges(6,5) * t741;
t479 = Icges(6,1) * t574 + Icges(6,4) * t573 + Icges(6,5) * t742;
t478 = Icges(6,4) * t576 + Icges(6,2) * t575 + Icges(6,6) * t741;
t477 = Icges(6,4) * t574 + Icges(6,2) * t573 + Icges(6,6) * t742;
t476 = Icges(6,5) * t576 + Icges(6,6) * t575 + Icges(6,3) * t741;
t475 = Icges(6,5) * t574 + Icges(6,6) * t573 + Icges(6,3) * t742;
t474 = qJD(1) * t546 - t609 * t612 + t639;
t473 = t609 * t615 + (-t543 - t759) * qJD(1);
t472 = t543 * t612 - t546 * t615;
t471 = qJD(1) * t548 - t608 * t614 + t679;
t470 = qJD(1) * t547 - t610 * t613 + t679;
t469 = -t722 + t608 * t617 + (-t545 + t723) * qJD(1);
t468 = -t722 + t610 * t616 + (-t544 + t723) * qJD(1);
t467 = -t612 * t763 - t588 * t592 + t639 + (t521 - t579) * qJD(1);
t466 = t615 * t763 + t588 * t594 + (-t519 + t578 - t759) * qJD(1);
t465 = t545 * t614 - t548 * t617 + t736;
t464 = t544 * t613 - t547 * t616 + t736;
t463 = qJD(1) * t522 - t587 * t593 + t672;
t462 = t587 * t595 + (-t520 + t721) * qJD(1) + t681;
t461 = t519 * t592 - t521 * t594 - t578 * t612 + t579 * t615;
t460 = -t613 * t762 - t566 * t590 + (t497 + t504) * qJD(1) + t679;
t459 = -t722 + t616 * t762 + t566 * t591 + (-t496 - t503 + t723) * qJD(1);
t458 = t520 * t593 - t522 * t595 + t680;
t457 = t496 * t590 - t497 * t591 + t503 * t613 - t504 * t616 + t736;
t456 = qJD(1) * t550 + t482 * t618 - t501 * t517 - t589 * t593 + t672;
t455 = -t481 * t618 + t501 * t518 + t589 * t595 + (-t549 + t721) * qJD(1) + t681;
t454 = t481 * t517 - t482 * t518 + t549 * t593 - t550 * t595 + t680;
t1 = t614 * (t673 * t660 + t666 * t664) / 0.2e1 + m(10) * (t461 ^ 2 + t466 ^ 2 + t467 ^ 2) / 0.2e1 + m(8) * (t464 ^ 2 + t468 ^ 2 + t470 ^ 2) / 0.2e1 + m(4) * (t465 ^ 2 + t469 ^ 2 + t471 ^ 2) / 0.2e1 + m(9) * (t472 ^ 2 + t473 ^ 2 + t474 ^ 2) / 0.2e1 + t618 * ((-t475 * t518 - t476 * t517 - t498 * t618) * t638 + ((-t478 * t658 + t480 * t662) * t517 + (-t477 * t658 + t479 * t662) * t518 + (-t499 * t658 + t500 * t662) * t618) * t636) / 0.2e1 + m(11) * (t457 ^ 2 + t459 ^ 2 + t460 ^ 2) / 0.2e1 + m(5) * (t458 ^ 2 + t462 ^ 2 + t463 ^ 2) / 0.2e1 + m(6) * (t454 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + m(3) * (t484 ^ 2 + t485 ^ 2 + t487 ^ 2) / 0.2e1 + m(7) * (t483 ^ 2 + t486 ^ 2 + t488 ^ 2) / 0.2e1 + (((-t561 * t659 + t565 * t663) * t660 - (-t560 * t659 + t564 * t663) * t664) * qJD(2) + ((t559 * t661 + t563 * t657) * t660 - (t558 * t661 + t562 * t657) * t664) * qJD(6) + (t536 * t651 + t542 * t648) * t614 + (t535 * t651 + t541 * t648) * t617 + (-t534 * t647 + t540 * t650) * t613 + (-t533 * t647 + t539 * t650) * t616 + (-t532 * t646 + t538 * t649) * t612 + (-t531 * t646 + t537 * t649) * t615 + (t512 * t638 + t516 * t636) * t593 + (t511 * t638 + t515 * t636) * t595 + (t510 * t635 - t514 * t637) * t592 + (t509 * t635 - t513 * t637) * t594 + (-t493 * t630 - t495 * t631) * t590 + (-t492 * t630 - t494 * t631) * t591 + (-t622 * t659 + t624 * t663 + t621 * t661 + t623 * t657 + t604 * t651 + t607 * t648 - t603 * t647 + t606 * t650 - t602 * t646 + t605 * t649 + t584 * t638 + t586 * t636 + t583 * t635 - t585 * t637 - t552 * t630 - t553 * t631) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t626 ^ 2 + t628 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + t617 * (t666 * t660 - t673 * t664) / 0.2e1 + t613 * (t674 * t660 + t667 * t664) / 0.2e1 + t616 * (t667 * t660 - t674 * t664) / 0.2e1 + t612 * (t675 * t660 + t668 * t664) / 0.2e1 + t615 * (t668 * t660 - t675 * t664) / 0.2e1 + t593 * (t676 * t660 + t669 * t664) / 0.2e1 + t595 * (t669 * t660 - t676 * t664) / 0.2e1 + t592 * (t677 * t660 + t670 * t664) / 0.2e1 + t594 * (t670 * t660 - t677 * t664) / 0.2e1 + t590 * (t678 * t660 + t671 * t664) / 0.2e1 + t591 * (t671 * t660 - t678 * t664) / 0.2e1 + t517 * ((t476 * t741 + t478 * t575 + t480 * t576) * t517 + (t475 * t741 + t477 * t575 + t479 * t576) * t518 + (t498 * t741 + t499 * t575 + t500 * t576) * t618) / 0.2e1 + t518 * ((t476 * t742 + t478 * t573 + t480 * t574) * t517 + (t475 * t742 + t477 * t573 + t479 * t574) * t518 + (t498 * t742 + t499 * t573 + t500 * t574) * t618) / 0.2e1 - ((-t664 * t619 + t660 * t686) * qJD(1) + (t765 * t554 + (t689 * t660 + (-t555 + t690) * t664) * t660) * qJD(6)) * t728 / 0.2e1 - ((-t620 * t664 + t685 * t660) * qJD(1) + (t556 * t765 + (t687 * t660 + (-t557 + t688) * t664) * t660) * qJD(2)) * t731 / 0.2e1 + ((t660 * t619 + t664 * t686) * qJD(1) + (t766 * t555 + (t690 * t664 + (-t554 + t689) * t660) * t664) * qJD(6)) * t729 / 0.2e1 + ((t660 * t620 + t664 * t685) * qJD(1) + (t766 * t557 + (t688 * t664 + (-t556 + t687) * t660) * t664) * qJD(2)) * t645 / 0.2e1;
T = t1;
