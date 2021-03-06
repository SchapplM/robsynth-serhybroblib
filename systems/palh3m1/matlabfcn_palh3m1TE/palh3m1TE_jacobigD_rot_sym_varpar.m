% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh3m1TE
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JgD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = palh3m1TE_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobigD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_jacobigD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobigD_rot_sym_varpar: pkin has to be [19x1] (double)');
JgD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:03
	% EndTime: 2020-04-18 09:52:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = qJD(1) * cos(qJ(1));
	t39 = qJD(1) * sin(qJ(1));
	t1 = [0, t40, t40, 0; 0, t39, t39, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:42
	% EndTime: 2020-04-18 09:53:57
	% DurationCPUTime: 45.92s
	% Computational Cost: add. (1067500->212), mult. (1666316->433), div. (62252->14), fcn. (1041364->16), ass. (0->213)
	t572 = pkin(5) ^ 2;
	t576 = pkin(1) ^ 2;
	t680 = sin(qJ(2));
	t681 = sin(pkin(16));
	t682 = cos(qJ(2));
	t683 = cos(pkin(16));
	t555 = t680 * t681 - t682 * t683;
	t676 = pkin(5) * t555;
	t642 = -0.2e1 * pkin(1) * t676 + t576;
	t550 = t572 + t642;
	t547 = 0.1e1 / t550;
	t548 = 0.1e1 / t550 ^ 2;
	t706 = 0.4e1 * t547 * t548;
	t556 = t680 * t683 + t682 * t681;
	t553 = t556 * qJD(2);
	t575 = 0.1e1 / pkin(2);
	t546 = pkin(2) ^ 2 - pkin(6) ^ 2 + t550;
	t551 = pkin(1) - t676;
	t697 = -pkin(6) - pkin(2);
	t544 = (pkin(5) - t697) * (pkin(5) + t697) + t642;
	t696 = -pkin(6) + pkin(2);
	t545 = (pkin(5) - t696) * (pkin(5) + t696) + t642;
	t649 = t545 * t544;
	t578 = sqrt(-t649);
	t675 = pkin(5) * t556;
	t538 = t546 * t675 + t551 * t578;
	t569 = cos(qJ(3));
	t650 = t538 * t569;
	t645 = t556 * t578;
	t537 = -pkin(5) * t645 + t546 * t551;
	t567 = sin(qJ(3));
	t653 = t537 * t567;
	t600 = t650 / 0.2e1 + t653 / 0.2e1;
	t698 = pkin(1) * pkin(5);
	t640 = t548 * t698;
	t606 = 0.2e1 * (t544 + t545) * t698;
	t539 = t553 * t606;
	t554 = t555 * qJD(2);
	t702 = t556 * t572;
	t635 = t553 * t702;
	t618 = pkin(1) * t635;
	t542 = 0.1e1 / t578;
	t688 = t542 / 0.2e1;
	t628 = t551 * t688;
	t647 = t553 * t578;
	t525 = t539 * t628 - 0.2e1 * t618 + (-t554 * t546 - t647) * pkin(5);
	t659 = t525 * t567;
	t703 = t539 * t542;
	t530 = -t675 * t703 / 0.2e1;
	t674 = t551 * pkin(1);
	t624 = t546 + 0.2e1 * t674;
	t646 = t554 * t578;
	t524 = t530 + (-t624 * t553 + t646) * pkin(5);
	t660 = t524 * t569;
	t651 = t538 * t567;
	t652 = t537 * t569;
	t700 = t651 - t652;
	t500 = ((-t659 / 0.2e1 + t660 / 0.2e1 - t600 * qJD(3)) * t547 - t700 * t553 * t640) * t575;
	t607 = t650 + t653;
	t595 = t607 * t553;
	t601 = -t651 / 0.2e1 + t652 / 0.2e1;
	t658 = t525 * t569;
	t661 = t524 * t567;
	t501 = (t595 * t640 + (t658 / 0.2e1 + t661 / 0.2e1 + t601 * qJD(3)) * t547) * t575;
	t564 = pkin(18) + pkin(19);
	t560 = sin(t564);
	t561 = cos(t564);
	t496 = t500 * t561 - t501 * t560;
	t573 = pkin(4) ^ 2;
	t648 = t547 * t575;
	t528 = t601 * t648;
	t529 = t600 * t648;
	t519 = -t528 * t561 + t529 * t560;
	t677 = pkin(3) * t519;
	t699 = -2 * pkin(4);
	t644 = -t677 * t699 + t573;
	t695 = -pkin(8) - pkin(10);
	t509 = (pkin(3) - t695) * (pkin(3) + t695) + t644;
	t694 = -pkin(8) + pkin(10);
	t510 = (pkin(3) - t694) * (pkin(3) + t694) + t644;
	t664 = t510 * t509;
	t577 = sqrt(-t664);
	t705 = t496 * t577;
	t641 = 2 * pkin(4);
	t605 = pkin(3) * (t509 + t510) * t641;
	t487 = t496 * t605;
	t574 = pkin(3) ^ 2;
	t515 = t574 + t644;
	t511 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t515;
	t516 = pkin(4) + t677;
	t623 = t516 * t699 - t511;
	t503 = 0.1e1 / t577;
	t609 = -t528 * t560 - t529 * t561;
	t689 = -t609 / 0.2e1;
	t632 = t503 * t689;
	t610 = -t500 * t560 - t501 * t561;
	t701 = t610 * t577;
	t469 = (t487 * t632 + t623 * t496 - t701) * pkin(3);
	t625 = t609 * t574 * t699;
	t690 = t503 / 0.2e1;
	t633 = t516 * t690;
	t470 = t487 * t633 + t496 * t625 + (t511 * t610 - t705) * pkin(3);
	t512 = 0.1e1 / t515;
	t571 = 0.1e1 / pkin(10);
	t662 = t512 * t571;
	t678 = pkin(3) * t609;
	t492 = t511 * t678 + t516 * t577;
	t565 = sin(pkin(17));
	t667 = t492 * t565;
	t491 = t511 * t516 - t577 * t678;
	t566 = cos(pkin(17));
	t668 = t491 * t566;
	t485 = (-t668 / 0.2e1 + t667 / 0.2e1) * t662;
	t482 = 0.1e1 / t485 ^ 2;
	t513 = 0.1e1 / t515 ^ 2;
	t639 = pkin(3) * pkin(4) * t513;
	t589 = (t667 - t668) * t639;
	t686 = -t566 / 0.2e1;
	t687 = t565 / 0.2e1;
	t704 = ((t469 * t686 + t470 * t687) * t512 + t496 * t589) * t571 * t482;
	t481 = 0.1e1 / t485;
	t693 = -t487 / 0.2e1;
	t540 = t556 * t606;
	t679 = pkin(1) * t572;
	t527 = t540 * t628 - 0.2e1 * t556 ^ 2 * t679 + (-t546 * t555 - t645) * pkin(5);
	t654 = t527 * t569;
	t614 = t540 * t688 + t546;
	t526 = (t555 * t578 + (-t614 - 0.2e1 * t674) * t556) * pkin(5);
	t657 = t526 * t567;
	t602 = t654 / 0.2e1 + t657 / 0.2e1;
	t621 = t556 * t640;
	t507 = (t602 * t547 + t607 * t621) * t575;
	t655 = t527 * t567;
	t656 = t526 * t569;
	t603 = -t655 / 0.2e1 + t656 / 0.2e1;
	t508 = (-t603 * t547 + t621 * t700) * t575;
	t499 = -t507 * t560 - t508 * t561;
	t488 = t499 * t605;
	t692 = -t488 / 0.2e1;
	t497 = t609 * t605;
	t691 = -t497 / 0.2e1;
	t685 = t566 / 0.2e1;
	t684 = t567 / 0.2e1;
	t666 = t492 * t566;
	t669 = t491 * t565;
	t588 = (t666 + t669) * t639;
	t462 = ((t469 * t687 + t470 * t685) * t512 + t496 * t588) * t571;
	t486 = (t666 / 0.2e1 + t669 / 0.2e1) * t662;
	t484 = t486 ^ 2;
	t475 = t482 * t484 + 0.1e1;
	t671 = t482 * t486;
	t672 = t481 * t704;
	t673 = 0.2e1 * (t462 * t671 - t484 * t672) / t475 ^ 2;
	t670 = t487 * t503 / t664;
	t665 = t496 * t574;
	t663 = t512 * t566;
	t568 = sin(qJ(1));
	t562 = qJD(1) * t568;
	t570 = cos(qJ(1));
	t563 = qJD(1) * t570;
	t637 = t573 * t665;
	t636 = 0.1e1 / t649 * t540 * t703;
	t634 = t670 / 0.4e1;
	t631 = t512 * t687;
	t630 = -t663 / 0.2e1;
	t629 = t663 / 0.2e1;
	t627 = t574 * t641;
	t626 = 0.4e1 * pkin(4) * t665;
	t473 = 0.1e1 / t475;
	t620 = t473 * t639;
	t619 = -0.8e1 * t637;
	t617 = t512 * t513 * t637;
	t616 = t576 * t635;
	t615 = -t609 * t670 / 0.4e1;
	t611 = t566 * t617;
	t604 = 0.4e1 * t565 * t617;
	t599 = -0.4e1 * t491 * t611;
	t598 = 0.4e1 * t492 * t611;
	t597 = t499 * t604;
	t596 = t609 * t604;
	t594 = -t473 * t704 - t481 * t673;
	t498 = -t507 * t561 + t508 * t560;
	t471 = (t488 * t632 - t498 * t577 + t623 * t499) * pkin(3);
	t536 = -t554 * t606 - 0.8e1 * t616;
	t505 = t530 + (t636 / 0.4e1 + t536 * t688) * t551 + (0.2e1 * t553 * t555 + 0.4e1 * t554 * t556) * t679 + (-t614 * t553 + t646) * pkin(5);
	t506 = 0.4e1 * t618 + (t647 - t556 * t636 / 0.4e1 + t624 * t554 + (t555 * t539 / 0.2e1 + t554 * t540 / 0.2e1 - t556 * t536 / 0.2e1) * t542) * pkin(5);
	t489 = (t700 * t616 * t706 + (-t506 * t569 / 0.2e1 + t505 * t684 + t602 * qJD(3)) * t547 + (-t700 * t554 - (-t655 + t656) * t553 + (t607 * qJD(3) + t659 - t660) * t556) * t640) * t575;
	t490 = (t576 * t595 * t702 * t706 + (t505 * t569 / 0.2e1 + t506 * t684 + t603 * qJD(3)) * t547 + (-t607 * t554 - (-t654 - t657) * t553 + (-qJD(3) * t700 + t658 + t661) * t556) * t640) * t575;
	t477 = -t489 * t561 - t490 * t560;
	t593 = t469 * t499 + t471 * t496 + t477 * t491;
	t479 = (t497 * t632 - t519 * t577 + t609 * t623) * pkin(3);
	t592 = t469 * t609 + t479 * t496 + t491 * t610;
	t472 = t488 * t633 + t499 * t625 + (t498 * t511 - t499 * t577) * pkin(3);
	t591 = t470 * t499 + t472 * t496 + t477 * t492;
	t480 = t497 * t633 + t609 * t625 + (t511 * t519 - t577 * t609) * pkin(3);
	t590 = t470 * t609 + t480 * t496 + t492 * t610;
	t587 = t671 * t673 + (-t462 * t482 + 0.2e1 * t486 * t672) * t473;
	t478 = t605 * t610 + t609 * t619;
	t476 = t489 * t560 - t490 * t561;
	t468 = ((t479 * t687 + t480 * t685) * t512 + t609 * t588) * t571;
	t467 = ((t479 * t686 + t480 * t687) * t512 + t609 * t589) * t571;
	t466 = t477 * t605 + t499 * t619;
	t465 = ((t471 * t687 + t472 * t685) * t512 + t499 * t588) * t571;
	t464 = ((t471 * t686 + t472 * t687) * t512 + t499 * t589) * t571;
	t463 = t609 * t626 + (t705 + t497 * t615 + t623 * t610 + (t478 * t689 + t519 * t693 + t610 * t691) * t503) * pkin(3);
	t460 = (t478 * t690 + t497 * t634) * t516 + (-t496 * t519 - 0.2e1 * t609 * t610) * t627 + (-t496 * t511 - t701 + (t496 * t691 + t609 * t693) * t503) * pkin(3);
	t459 = t499 * t626 + (-t476 * t577 + t488 * t615 + t623 * t477 + (t466 * t689 + t498 * t693 + t610 * t692) * t503) * pkin(3);
	t458 = (-t467 * t671 + t468 * t481) * t473;
	t457 = (t466 * t690 + t488 * t634) * t516 + (-t477 * t609 - t496 * t498 - t499 * t610) * t627 + (t476 * t511 - t477 * t577 + (t496 * t692 + t499 * t693) * t503) * pkin(3);
	t455 = (-t464 * t671 + t465 * t481) * t473;
	t454 = t594 * t468 + t587 * t467 + (((t460 * t629 + t463 * t631 + t491 * t596 + t598 * t609) * t481 - (t460 * t631 + t463 * t630 + t492 * t596 + t599 * t609) * t671) * t473 + ((t590 * t481 + t592 * t671) * t566 + (t592 * t481 - t590 * t671) * t565) * t620) * t571;
	t453 = t594 * t465 + t587 * t464 + (((t457 * t629 + t459 * t631 + t491 * t597 + t499 * t598) * t481 - (t457 * t631 + t459 * t630 + t492 * t597 + t499 * t599) * t671) * t473 + ((t591 * t481 + t593 * t671) * t566 + (t593 * t481 - t591 * t671) * t565) * t620) * t571;
	t1 = [0, t453 * t568 + t455 * t563 + t563, t454 * t568 + t458 * t563 + t563, 0; 0, -t453 * t570 + t455 * t562 + t562, -t454 * t570 + t458 * t562 + t562, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:54:37
	% EndTime: 2020-04-18 09:55:55
	% DurationCPUTime: 48.34s
	% Computational Cost: add. (1133697->223), mult. (1770932->451), div. (66164->14), fcn. (1106224->16), ass. (0->219)
	t790 = pkin(5) ^ 2;
	t794 = pkin(1) ^ 2;
	t784 = sin(qJ(2));
	t787 = cos(qJ(2));
	t903 = sin(pkin(16));
	t904 = cos(pkin(16));
	t771 = t784 * t903 - t787 * t904;
	t900 = pkin(5) * t771;
	t865 = -0.2e1 * pkin(1) * t900 + t794;
	t766 = t790 + t865;
	t763 = 0.1e1 / t766;
	t764 = 0.1e1 / t766 ^ 2;
	t929 = 0.4e1 * t763 * t764;
	t772 = t784 * t904 + t787 * t903;
	t769 = t772 * qJD(2);
	t793 = 0.1e1 / pkin(2);
	t762 = pkin(2) ^ 2 - pkin(6) ^ 2 + t766;
	t767 = pkin(1) - t900;
	t918 = -pkin(6) - pkin(2);
	t760 = (pkin(5) - t918) * (pkin(5) + t918) + t865;
	t917 = -pkin(6) + pkin(2);
	t761 = (pkin(5) - t917) * (pkin(5) + t917) + t865;
	t872 = t761 * t760;
	t796 = sqrt(-t872);
	t899 = pkin(5) * t772;
	t754 = t762 * t899 + t767 * t796;
	t786 = cos(qJ(3));
	t873 = t754 * t786;
	t868 = t772 * t796;
	t753 = -pkin(5) * t868 + t767 * t762;
	t783 = sin(qJ(3));
	t876 = t753 * t783;
	t821 = t873 / 0.2e1 + t876 / 0.2e1;
	t919 = pkin(1) * pkin(5);
	t863 = t764 * t919;
	t827 = 0.2e1 * (t760 + t761) * t919;
	t755 = t769 * t827;
	t770 = t771 * qJD(2);
	t925 = t772 * t790;
	t858 = t769 * t925;
	t839 = pkin(1) * t858;
	t758 = 0.1e1 / t796;
	t909 = t758 / 0.2e1;
	t849 = t767 * t909;
	t870 = t769 * t796;
	t741 = t755 * t849 - 0.2e1 * t839 + (-t770 * t762 - t870) * pkin(5);
	t882 = t741 * t783;
	t926 = t755 * t758;
	t746 = -t899 * t926 / 0.2e1;
	t898 = t767 * pkin(1);
	t845 = t762 + 0.2e1 * t898;
	t869 = t770 * t796;
	t740 = t746 + (-t845 * t769 + t869) * pkin(5);
	t883 = t740 * t786;
	t874 = t754 * t783;
	t875 = t753 * t786;
	t923 = t874 - t875;
	t716 = ((-t882 / 0.2e1 + t883 / 0.2e1 - t821 * qJD(3)) * t763 - t923 * t769 * t863) * t793;
	t830 = t873 + t876;
	t816 = t769 * t830;
	t822 = -t874 / 0.2e1 + t875 / 0.2e1;
	t881 = t741 * t786;
	t884 = t740 * t783;
	t717 = (t816 * t863 + (t881 / 0.2e1 + t884 / 0.2e1 + t822 * qJD(3)) * t763) * t793;
	t780 = pkin(18) + pkin(19);
	t776 = sin(t780);
	t777 = cos(t780);
	t712 = t777 * t716 - t776 * t717;
	t791 = pkin(4) ^ 2;
	t871 = t763 * t793;
	t744 = t822 * t871;
	t745 = t821 * t871;
	t735 = -t777 * t744 + t776 * t745;
	t901 = pkin(3) * t735;
	t920 = -2 * pkin(4);
	t867 = -t901 * t920 + t791;
	t916 = -pkin(8) - pkin(10);
	t725 = (pkin(3) - t916) * (pkin(3) + t916) + t867;
	t915 = -pkin(8) + pkin(10);
	t726 = (pkin(3) - t915) * (pkin(3) + t915) + t867;
	t887 = t726 * t725;
	t795 = sqrt(-t887);
	t928 = t712 * t795;
	t864 = 2 * pkin(4);
	t826 = pkin(3) * (t725 + t726) * t864;
	t701 = t712 * t826;
	t792 = pkin(3) ^ 2;
	t731 = t792 + t867;
	t727 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t731;
	t732 = pkin(4) + t901;
	t844 = t732 * t920 - t727;
	t719 = 0.1e1 / t795;
	t832 = -t776 * t744 - t777 * t745;
	t910 = -t832 / 0.2e1;
	t853 = t719 * t910;
	t833 = -t776 * t716 - t777 * t717;
	t924 = t833 * t795;
	t683 = (t701 * t853 + t844 * t712 - t924) * pkin(3);
	t846 = t832 * t792 * t920;
	t911 = t719 / 0.2e1;
	t854 = t732 * t911;
	t684 = t701 * t854 + t712 * t846 + (t727 * t833 - t928) * pkin(3);
	t728 = 0.1e1 / t731;
	t789 = 0.1e1 / pkin(10);
	t729 = 0.1e1 / t731 ^ 2;
	t862 = pkin(3) * pkin(4) * t729;
	t902 = pkin(3) * t832;
	t708 = t727 * t902 + t732 * t795;
	t781 = sin(pkin(17));
	t890 = t708 * t781;
	t707 = t727 * t732 - t795 * t902;
	t782 = cos(pkin(17));
	t891 = t707 * t782;
	t810 = (t890 - t891) * t862;
	t907 = -t782 / 0.2e1;
	t908 = t781 / 0.2e1;
	t673 = ((t683 * t907 + t684 * t908) * t728 + t712 * t810) * t789;
	t885 = t728 * t789;
	t699 = (-t891 / 0.2e1 + t890 / 0.2e1) * t885;
	t696 = 0.1e1 / t699 ^ 2;
	t927 = t673 * t696;
	t922 = qJD(2) + qJD(3);
	t889 = t708 * t782;
	t892 = t707 * t781;
	t700 = (t889 / 0.2e1 + t892 / 0.2e1) * t885;
	t828 = t784 * t783 - t787 * t786;
	t829 = t787 * t783 + t784 * t786;
	t921 = t829 * t699 - t828 * t700;
	t695 = 0.1e1 / t699;
	t914 = -t701 / 0.2e1;
	t756 = t772 * t827;
	t897 = t790 * pkin(1);
	t743 = t756 * t849 - 0.2e1 * t772 ^ 2 * t897 + (-t771 * t762 - t868) * pkin(5);
	t877 = t743 * t786;
	t835 = t756 * t909 + t762;
	t742 = (t771 * t796 + (-t835 - 0.2e1 * t898) * t772) * pkin(5);
	t880 = t742 * t783;
	t823 = t877 / 0.2e1 + t880 / 0.2e1;
	t842 = t772 * t863;
	t723 = (t823 * t763 + t830 * t842) * t793;
	t878 = t743 * t783;
	t879 = t742 * t786;
	t824 = -t878 / 0.2e1 + t879 / 0.2e1;
	t724 = (-t824 * t763 + t842 * t923) * t793;
	t715 = -t776 * t723 - t777 * t724;
	t702 = t715 * t826;
	t913 = -t702 / 0.2e1;
	t713 = t832 * t826;
	t912 = -t713 / 0.2e1;
	t906 = t782 / 0.2e1;
	t905 = t783 / 0.2e1;
	t809 = (t889 + t892) * t862;
	t674 = ((t683 * t908 + t684 * t906) * t728 + t712 * t809) * t789;
	t698 = t700 ^ 2;
	t689 = t696 * t698 + 0.1e1;
	t894 = t696 * t700;
	t895 = t695 * t927;
	t896 = 0.2e1 * (t674 * t894 - t698 * t895) / t689 ^ 2;
	t893 = t701 * t719 / t887;
	t888 = t712 * t792;
	t886 = t728 * t782;
	t785 = sin(qJ(1));
	t778 = qJD(1) * t785;
	t788 = cos(qJ(1));
	t779 = qJD(1) * t788;
	t860 = t791 * t888;
	t859 = 0.1e1 / t872 * t756 * t926;
	t855 = t893 / 0.4e1;
	t852 = t728 * t908;
	t851 = -t886 / 0.2e1;
	t850 = t886 / 0.2e1;
	t848 = t792 * t864;
	t847 = 0.4e1 * pkin(4) * t888;
	t687 = 0.1e1 / t689;
	t841 = t687 * t862;
	t840 = -0.8e1 * t860;
	t838 = t728 * t729 * t860;
	t837 = t794 * t858;
	t836 = -t832 * t893 / 0.4e1;
	t834 = t782 * t838;
	t825 = 0.4e1 * t781 * t838;
	t820 = -0.4e1 * t707 * t834;
	t819 = 0.4e1 * t708 * t834;
	t818 = t715 * t825;
	t817 = t832 * t825;
	t815 = -t687 * t927 - t695 * t896;
	t714 = -t777 * t723 + t776 * t724;
	t685 = (t702 * t853 - t714 * t795 + t844 * t715) * pkin(3);
	t752 = -t770 * t827 - 0.8e1 * t837;
	t721 = t746 + (t859 / 0.4e1 + t752 * t909) * t767 + (0.2e1 * t769 * t771 + 0.4e1 * t770 * t772) * t897 + (-t835 * t769 + t869) * pkin(5);
	t722 = 0.4e1 * t839 + (t870 - t772 * t859 / 0.4e1 + t845 * t770 + (t771 * t755 / 0.2e1 + t770 * t756 / 0.2e1 - t772 * t752 / 0.2e1) * t758) * pkin(5);
	t703 = (t923 * t837 * t929 + (-t722 * t786 / 0.2e1 + t721 * t905 + t823 * qJD(3)) * t763 + (-t923 * t770 - (-t878 + t879) * t769 + (t830 * qJD(3) + t882 - t883) * t772) * t863) * t793;
	t704 = (t794 * t816 * t925 * t929 + (t721 * t786 / 0.2e1 + t722 * t905 + t824 * qJD(3)) * t763 + (-t830 * t770 - (-t877 - t880) * t769 + (-qJD(3) * t923 + t881 + t884) * t772) * t863) * t793;
	t691 = -t777 * t703 - t776 * t704;
	t814 = t683 * t715 + t685 * t712 + t691 * t707;
	t693 = (t713 * t853 - t735 * t795 + t832 * t844) * pkin(3);
	t813 = t683 * t832 + t693 * t712 + t707 * t833;
	t686 = t702 * t854 + t715 * t846 + (t714 * t727 - t715 * t795) * pkin(3);
	t812 = t684 * t715 + t686 * t712 + t691 * t708;
	t694 = t713 * t854 + t832 * t846 + (t735 * t727 - t795 * t832) * pkin(3);
	t811 = t684 * t832 + t694 * t712 + t708 * t833;
	t808 = t894 * t896 + (-t674 * t696 + 0.2e1 * t700 * t895) * t687;
	t805 = -t829 * t673 + t828 * t674 + t922 * (t828 * t699 + t829 * t700);
	t692 = t826 * t833 + t832 * t840;
	t690 = t776 * t703 - t777 * t704;
	t680 = ((t693 * t908 + t694 * t906) * t728 + t832 * t809) * t789;
	t679 = ((t693 * t907 + t694 * t908) * t728 + t832 * t810) * t789;
	t678 = t691 * t826 + t715 * t840;
	t677 = ((t685 * t908 + t686 * t906) * t728 + t715 * t809) * t789;
	t676 = ((t685 * t907 + t686 * t908) * t728 + t715 * t810) * t789;
	t675 = t832 * t847 + (t928 + t713 * t836 + t844 * t833 + (t692 * t910 + t735 * t914 + t833 * t912) * t719) * pkin(3);
	t672 = (t692 * t911 + t713 * t855) * t732 + (-t712 * t735 - 0.2e1 * t832 * t833) * t848 + (-t712 * t727 - t924 + (t712 * t912 + t832 * t914) * t719) * pkin(3);
	t671 = t715 * t847 + (-t690 * t795 + t702 * t836 + t844 * t691 + (t678 * t910 + t714 * t914 + t833 * t913) * t719) * pkin(3);
	t670 = (-t679 * t894 + t680 * t695) * t687;
	t669 = (t678 * t911 + t702 * t855) * t732 + (-t691 * t832 - t712 * t714 - t715 * t833) * t848 + (t690 * t727 - t691 * t795 + (t712 * t913 + t715 * t914) * t719) * pkin(3);
	t667 = (-t676 * t894 + t677 * t695) * t687;
	t666 = t815 * t680 + t808 * t679 + (((t672 * t850 + t675 * t852 + t707 * t817 + t819 * t832) * t695 - (t672 * t852 + t675 * t851 + t708 * t817 + t820 * t832) * t894) * t687 + ((t811 * t695 + t813 * t894) * t782 + (t813 * t695 - t811 * t894) * t781) * t841) * t789;
	t665 = t815 * t677 + t808 * t676 + (((t669 * t850 + t671 * t852 + t707 * t818 + t715 * t819) * t695 - (t669 * t852 + t671 * t851 + t708 * t818 + t715 * t820) * t894) * t687 + ((t812 * t695 + t814 * t894) * t782 + (t814 * t695 - t812 * t894) * t781) * t841) * t789;
	t1 = [0, t665 * t785 + t667 * t779 + t779, t666 * t785 + t670 * t779 + t779, t921 * t778 + t805 * t788; 0, -t665 * t788 + t667 * t778 + t778, -t666 * t788 + t670 * t778 + t778, -t779 * t921 + t805 * t785; 0, 0, 0, -t828 * t673 - t829 * t674 - t921 * t922;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:04
	% EndTime: 2020-04-18 09:52:07
	% DurationCPUTime: 1.31s
	% Computational Cost: add. (12594->80), mult. (20244->171), div. (678->9), fcn. (12648->10), ass. (0->96)
	t260 = sin(qJ(2));
	t261 = sin(pkin(16));
	t262 = cos(qJ(2));
	t263 = cos(pkin(16));
	t200 = t260 * t263 + t262 * t261;
	t258 = pkin(1) * t200;
	t197 = t200 * qJD(2);
	t269 = pkin(1) * pkin(5);
	t237 = t197 * t269;
	t211 = pkin(1) ^ 2;
	t199 = t260 * t261 - t262 * t263;
	t259 = pkin(1) * t199;
	t241 = -0.2e1 * pkin(5) * t259 + t211;
	t268 = -pkin(6) - pkin(2);
	t188 = (pkin(5) - t268) * (pkin(5) + t268) + t241;
	t267 = -pkin(6) + pkin(2);
	t189 = (pkin(5) - t267) * (pkin(5) + t267) + t241;
	t242 = t188 + t189;
	t183 = t242 * t237;
	t248 = t189 * t188;
	t212 = sqrt(-t248);
	t186 = 0.1e1 / t212;
	t272 = t183 * t186;
	t179 = t258 * t272;
	t210 = pkin(5) ^ 2;
	t194 = t210 + t241;
	t190 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t194;
	t195 = -pkin(5) + t259;
	t239 = 0.2e1 * t195 * pkin(5);
	t234 = -t190 + t239;
	t198 = t199 * qJD(2);
	t244 = t198 * t212;
	t169 = -t179 + (t197 * t234 + t244) * pkin(1);
	t236 = t211 * t200 * t197;
	t230 = pkin(5) * t236;
	t245 = t197 * t212;
	t249 = t186 * t195;
	t170 = -t183 * t249 - 0.2e1 * t230 + (-t198 * t190 - t245) * pkin(1);
	t191 = 0.1e1 / t194;
	t209 = 0.1e1 / pkin(6);
	t246 = t191 * t209;
	t182 = t190 * t258 - t195 * t212;
	t205 = sin(pkin(15));
	t251 = t182 * t205;
	t243 = t200 * t212;
	t181 = -pkin(1) * t243 - t190 * t195;
	t207 = cos(pkin(15));
	t252 = t181 * t207;
	t177 = (t252 / 0.2e1 + t251 / 0.2e1) * t246;
	t174 = 0.1e1 / t177 ^ 2;
	t223 = t251 + t252;
	t192 = 0.1e1 / t194 ^ 2;
	t233 = t192 * t237;
	t264 = t207 / 0.2e1;
	t265 = t205 / 0.2e1;
	t273 = ((t169 * t264 + t170 * t265) * t191 + t223 * t233) * t209 * t174;
	t229 = t242 * t269;
	t184 = t200 * t229;
	t271 = t186 * t184;
	t228 = t210 * t236;
	t270 = -t186 * (-t198 * t229 - 0.4e1 * t228) - 0.1e1 / t248 * t183 * t271;
	t173 = 0.1e1 / t177;
	t266 = -t205 / 0.2e1;
	t257 = pkin(5) * t211;
	t250 = t182 * t207;
	t253 = t181 * t205;
	t222 = -t250 + t253;
	t161 = ((t169 * t266 + t170 * t264) * t191 - t222 * t233) * t209;
	t178 = (t250 / 0.2e1 - t253 / 0.2e1) * t246;
	t176 = t178 ^ 2;
	t168 = t174 * t176 + 0.1e1;
	t254 = t174 * t178;
	t255 = t173 * t273;
	t256 = (t161 * t254 - t176 * t255) / t168 ^ 2;
	t247 = t191 * t205;
	t231 = t190 + t271;
	t171 = (t199 * t212 + (-t231 + t239) * t200) * pkin(1);
	t172 = -t184 * t249 - 0.2e1 * t200 ^ 2 * t257 + (-t190 * t199 - t243) * pkin(1);
	t238 = t192 * t269;
	t232 = t200 * t238;
	t164 = ((t171 * t266 + t172 * t264) * t191 - t222 * t232) * t209;
	t165 = ((t171 * t264 + t172 * t265) * t191 + t223 * t232) * t209;
	t166 = 0.1e1 / t168;
	t240 = qJD(1) * (t164 * t173 - t254 * t165) * t166;
	t235 = t191 * t264;
	t225 = t191 * t192 * t228;
	t221 = t205 * t225;
	t220 = 0.4e1 * t207 * t225;
	t219 = t169 * t200 + t171 * t197 - t181 * t198;
	t218 = t170 * t200 + t172 * t197 - t182 * t198;
	t206 = cos(qJ(1));
	t204 = sin(qJ(1));
	t163 = 0.4e1 * t230 + (t199 * t272 + t245 + t270 * t200 - (t234 - t271) * t198) * pkin(1);
	t162 = -t179 + t270 * t195 + (0.2e1 * t197 * t199 + 0.4e1 * t198 * t200) * t257 + (-t197 * t231 + t244) * pkin(1);
	t157 = (-t166 * t273 - 0.2e1 * t173 * t256) * t164 + (0.2e1 * t254 * t256 + (-t161 * t174 + 0.2e1 * t178 * t255) * t166) * t165 + ((t162 * t235 + t182 * t220 - t163 * t247 / 0.2e1 - 0.4e1 * t181 * t221) * t173 - (t163 * t235 + t181 * t220 + t162 * t247 / 0.2e1 + 0.4e1 * t182 * t221) * t254 + ((t173 * t218 - t219 * t254) * t207 + (-t173 * t219 - t218 * t254) * t205) * t238) * t166 * t209;
	t1 = [0, t157 * t204 + t206 * t240, 0, 0; 0, -t157 * t206 + t204 * t240, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobigD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:06
	% EndTime: 2020-04-18 09:52:09
	% DurationCPUTime: 1.27s
	% Computational Cost: add. (12596->81), mult. (20246->173), div. (678->9), fcn. (12650->10), ass. (0->96)
	t276 = sin(qJ(2));
	t277 = sin(pkin(16));
	t278 = cos(qJ(2));
	t279 = cos(pkin(16));
	t216 = t276 * t279 + t278 * t277;
	t274 = pkin(5) * t216;
	t213 = t216 * qJD(2);
	t285 = pkin(1) * pkin(5);
	t254 = t213 * t285;
	t227 = pkin(1) ^ 2;
	t215 = t276 * t277 - t278 * t279;
	t275 = pkin(5) * t215;
	t256 = -0.2e1 * pkin(1) * t275 + t227;
	t284 = -pkin(6) - pkin(2);
	t204 = (pkin(5) - t284) * (pkin(5) + t284) + t256;
	t283 = -pkin(6) + pkin(2);
	t205 = (pkin(5) - t283) * (pkin(5) + t283) + t256;
	t257 = t204 + t205;
	t199 = t257 * t254;
	t263 = t205 * t204;
	t228 = sqrt(-t263);
	t202 = 0.1e1 / t228;
	t288 = t199 * t202;
	t195 = t274 * t288;
	t225 = pkin(5) ^ 2;
	t210 = t225 + t256;
	t206 = pkin(2) ^ 2 - pkin(6) ^ 2 + t210;
	t211 = pkin(1) - t275;
	t273 = t211 * pkin(1);
	t250 = t206 + 0.2e1 * t273;
	t214 = t215 * qJD(2);
	t259 = t214 * t228;
	t185 = -t195 + (-t250 * t213 + t259) * pkin(5);
	t253 = t225 * t216 * t213;
	t246 = pkin(1) * t253;
	t260 = t213 * t228;
	t264 = t202 * t211;
	t186 = t199 * t264 - 0.2e1 * t246 + (-t214 * t206 - t260) * pkin(5);
	t207 = 0.1e1 / t210;
	t226 = 0.1e1 / pkin(2);
	t261 = t207 * t226;
	t198 = t206 * t274 + t211 * t228;
	t220 = sin(pkin(19));
	t266 = t198 * t220;
	t258 = t216 * t228;
	t197 = -pkin(5) * t258 + t211 * t206;
	t221 = cos(pkin(19));
	t267 = t197 * t221;
	t193 = (-t267 / 0.2e1 + t266 / 0.2e1) * t261;
	t190 = 0.1e1 / t193 ^ 2;
	t239 = -t266 + t267;
	t208 = 0.1e1 / t210 ^ 2;
	t249 = t208 * t254;
	t281 = -t221 / 0.2e1;
	t282 = t220 / 0.2e1;
	t289 = ((t185 * t281 + t186 * t282) * t207 - t239 * t249) * t226 * t190;
	t245 = t257 * t285;
	t200 = t216 * t245;
	t287 = t202 * t200;
	t243 = t227 * t253;
	t286 = 0.1e1 / t263 * t199 * t287 + t202 * (-t214 * t245 - 0.4e1 * t243);
	t189 = 0.1e1 / t193;
	t280 = t221 / 0.2e1;
	t272 = t225 * pkin(1);
	t265 = t198 * t221;
	t268 = t197 * t220;
	t238 = t265 + t268;
	t176 = ((t185 * t282 + t186 * t280) * t207 + t238 * t249) * t226;
	t194 = (t265 / 0.2e1 + t268 / 0.2e1) * t261;
	t192 = t194 ^ 2;
	t184 = t192 * t190 + 0.1e1;
	t269 = t190 * t194;
	t270 = t189 * t289;
	t271 = (t176 * t269 - t192 * t270) / t184 ^ 2;
	t262 = t207 * t221;
	t255 = t208 * t285;
	t252 = t207 * t282;
	t247 = t206 + t287;
	t187 = (t215 * t228 + (-t247 - 0.2e1 * t273) * t216) * pkin(5);
	t188 = t200 * t264 - 0.2e1 * t216 ^ 2 * t272 + (-t215 * t206 - t258) * pkin(5);
	t248 = t216 * t255;
	t180 = ((t187 * t282 + t188 * t280) * t207 + t238 * t248) * t226;
	t181 = ((t187 * t281 + t188 * t282) * t207 - t239 * t248) * t226;
	t182 = 0.1e1 / t184;
	t251 = qJD(1) * ((t180 * t189 - t181 * t269) * t182 + 0.1e1);
	t240 = t207 * t208 * t243;
	t237 = t221 * t240;
	t236 = 0.4e1 * t220 * t240;
	t235 = t185 * t216 + t187 * t213 - t197 * t214;
	t234 = t186 * t216 + t188 * t213 - t198 * t214;
	t223 = cos(qJ(1));
	t222 = sin(qJ(1));
	t179 = 0.4e1 * t246 + (t215 * t288 + t260 - t286 * t216 - (-t250 - t287) * t214) * pkin(5);
	t178 = -t195 + t286 * t211 + (0.2e1 * t213 * t215 + 0.4e1 * t214 * t216) * t272 + (-t247 * t213 + t259) * pkin(5);
	t173 = (-t182 * t289 - 0.2e1 * t189 * t271) * t180 + (0.2e1 * t269 * t271 + (-t176 * t190 + 0.2e1 * t194 * t270) * t182) * t181 + ((t178 * t262 / 0.2e1 + 0.4e1 * t198 * t237 + t179 * t252 + t197 * t236) * t189 - (-t179 * t262 / 0.2e1 - 0.4e1 * t197 * t237 + t178 * t252 + t198 * t236) * t269 + ((t234 * t189 + t235 * t269) * t221 + (t235 * t189 - t234 * t269) * t220) * t255) * t182 * t226;
	t1 = [0, t173 * t222 + t223 * t251, 0, 0; 0, -t173 * t223 + t222 * t251, 0, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobigD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:51
	% EndTime: 2020-04-18 09:53:32
	% DurationCPUTime: 25.02s
	% Computational Cost: add. (542986->224), mult. (846766->490), div. (31210->18), fcn. (529626->16), ass. (0->228)
	t572 = pkin(5) ^ 2;
	t576 = pkin(1) ^ 2;
	t698 = sin(qJ(2));
	t699 = sin(pkin(16));
	t700 = cos(qJ(2));
	t701 = cos(pkin(16));
	t557 = t698 * t699 - t700 * t701;
	t691 = pkin(5) * t557;
	t654 = -0.2e1 * pkin(1) * t691 + t576;
	t552 = t572 + t654;
	t549 = 0.1e1 / t552;
	t550 = 0.1e1 / t552 ^ 2;
	t551 = t549 * t550;
	t558 = t698 * t701 + t700 * t699;
	t555 = t558 * qJD(2);
	t726 = t558 * t572;
	t643 = t555 * t726;
	t624 = t576 * t643;
	t611 = t551 * t624;
	t731 = 0.4e1 * t611;
	t575 = 0.1e1 / pkin(2);
	t548 = pkin(2) ^ 2 - pkin(6) ^ 2 + t552;
	t553 = pkin(1) - t691;
	t718 = -pkin(6) - pkin(2);
	t546 = (pkin(5) - t718) * (pkin(5) + t718) + t654;
	t717 = -pkin(6) + pkin(2);
	t547 = (pkin(5) - t717) * (pkin(5) + t717) + t654;
	t663 = t547 * t546;
	t578 = sqrt(-t663);
	t690 = pkin(5) * t558;
	t540 = t548 * t690 + t553 * t578;
	t569 = cos(qJ(3));
	t664 = t540 * t569;
	t658 = t558 * t578;
	t539 = -pkin(5) * t658 + t548 * t553;
	t567 = sin(qJ(3));
	t669 = t539 * t567;
	t595 = t664 / 0.2e1 + t669 / 0.2e1;
	t719 = pkin(1) * pkin(5);
	t651 = t550 * t719;
	t633 = t555 * t651;
	t604 = 0.2e1 * (t546 + t547) * t719;
	t541 = t555 * t604;
	t556 = t557 * qJD(2);
	t625 = pkin(1) * t643;
	t544 = 0.1e1 / t578;
	t706 = t544 / 0.2e1;
	t639 = t553 * t706;
	t660 = t555 * t578;
	t521 = t541 * t639 - 0.2e1 * t625 + (-t556 * t548 - t660) * pkin(5);
	t678 = t521 * t567;
	t727 = t541 * t544;
	t532 = -t690 * t727 / 0.2e1;
	t688 = t553 * pkin(1);
	t635 = t548 + 0.2e1 * t688;
	t659 = t556 * t578;
	t520 = t532 + (-t555 * t635 + t659) * pkin(5);
	t679 = t520 * t569;
	t665 = t540 * t567;
	t668 = t539 * t569;
	t723 = t665 - t668;
	t489 = ((-t678 / 0.2e1 + t679 / 0.2e1 - t595 * qJD(3)) * t549 - t723 * t633) * t575;
	t605 = t664 + t669;
	t593 = t605 * t555;
	t596 = -t665 / 0.2e1 + t668 / 0.2e1;
	t677 = t521 * t569;
	t680 = t520 * t567;
	t490 = (t593 * t651 + (t677 / 0.2e1 + t680 / 0.2e1 + t596 * qJD(3)) * t549) * t575;
	t564 = pkin(18) + pkin(19);
	t562 = sin(t564);
	t563 = cos(t564);
	t485 = t489 * t563 - t490 * t562;
	t574 = pkin(3) ^ 2;
	t573 = pkin(4) ^ 2;
	t661 = t549 * t575;
	t530 = t596 * t661;
	t531 = t595 * t661;
	t512 = -t530 * t563 + t531 * t562;
	t693 = pkin(4) * t512;
	t720 = -2 * pkin(3);
	t656 = -t693 * t720 + t573;
	t508 = t574 + t656;
	t504 = pkin(8) ^ 2 - pkin(10) ^ 2 + t508;
	t716 = -pkin(8) - pkin(10);
	t502 = (pkin(3) - t716) * (pkin(3) + t716) + t656;
	t715 = -pkin(8) + pkin(10);
	t503 = (pkin(3) - t715) * (pkin(3) + t715) + t656;
	t682 = t503 * t502;
	t577 = sqrt(-t682);
	t610 = -t489 * t562 - t490 * t563;
	t730 = -t485 * t577 + t610 * t504;
	t652 = 2 * pkin(3);
	t603 = pkin(4) * (t502 + t503) * t652;
	t467 = t485 * t603;
	t609 = -t530 * t562 - t531 * t563;
	t636 = t609 * t573 * t720;
	t509 = -pkin(3) - t693;
	t492 = 0.1e1 / t577;
	t711 = -t492 / 0.2e1;
	t641 = t509 * t711;
	t457 = t730 * pkin(4) + t467 * t641 + t485 * t636;
	t694 = pkin(4) * t609;
	t480 = -t504 * t509 - t577 * t694;
	t477 = 0.1e1 / t480 ^ 2;
	t481 = t504 * t694 - t509 * t577;
	t479 = t481 ^ 2;
	t470 = t477 * t479 + 0.1e1;
	t685 = t477 * t481;
	t634 = t509 * t652 - t504;
	t707 = -t609 / 0.2e1;
	t640 = t492 * t707;
	t724 = t610 * t577;
	t456 = (t467 * t640 + t485 * t634 - t724) * pkin(4);
	t476 = 0.1e1 / t480;
	t687 = t456 * t476 * t477;
	t729 = 0.2e1 * (t457 * t685 - t479 * t687) / t470 ^ 2;
	t565 = sin(pkin(19));
	t667 = t540 * t565;
	t566 = cos(pkin(19));
	t670 = t539 * t566;
	t528 = (-t670 / 0.2e1 + t667 / 0.2e1) * t661;
	t525 = 0.1e1 / t528 ^ 2;
	t608 = -t667 + t670;
	t704 = -t566 / 0.2e1;
	t705 = t565 / 0.2e1;
	t728 = ((t520 * t704 + t521 * t705) * t549 - t608 * t633) * t575 * t525;
	t468 = 0.1e1 / t470;
	t623 = pkin(3) * pkin(4) * pkin(8) * t468 * t485;
	t689 = pkin(8) * t508;
	t648 = t468 * t689;
	t627 = t477 * t648;
	t647 = t476 * t689;
	t722 = -t456 * t627 - 0.2e1 * t476 * t623 - t647 * t729;
	t646 = t481 * t689;
	t626 = t477 * t646;
	t721 = 0.2e1 * t468 * t646 * t687 - t457 * t627 + 0.2e1 * t623 * t685 + t626 * t729;
	t505 = 0.1e1 / t508;
	t524 = 0.1e1 / t528;
	t506 = 0.1e1 / t508 ^ 2;
	t714 = -t467 / 0.2e1;
	t542 = t558 * t604;
	t697 = pkin(1) * t572;
	t523 = t542 * t639 - 0.2e1 * t558 ^ 2 * t697 + (-t548 * t557 - t658) * pkin(5);
	t673 = t523 * t569;
	t621 = t542 * t706 + t548;
	t522 = (t557 * t578 + (-t621 - 0.2e1 * t688) * t558) * pkin(5);
	t676 = t522 * t567;
	t597 = t673 / 0.2e1 + t676 / 0.2e1;
	t632 = t558 * t651;
	t500 = (t549 * t597 + t605 * t632) * t575;
	t674 = t523 * t567;
	t675 = t522 * t569;
	t598 = -t674 / 0.2e1 + t675 / 0.2e1;
	t501 = (-t549 * t598 + t632 * t723) * t575;
	t488 = -t500 * t562 - t501 * t563;
	t471 = t488 * t603;
	t713 = -t471 / 0.2e1;
	t486 = t609 * t603;
	t712 = -t486 / 0.2e1;
	t538 = -t556 * t604 - 0.8e1 * t624;
	t644 = 0.1e1 / t663 * t542 * t727;
	t496 = t532 + (t644 / 0.4e1 + t538 * t706) * t553 + (0.2e1 * t555 * t557 + 0.4e1 * t556 * t558) * t697 + (-t555 * t621 + t659) * pkin(5);
	t710 = t496 / 0.2e1;
	t497 = 0.4e1 * t625 + (t660 - t558 * t644 / 0.4e1 + t635 * t556 + (t557 * t541 / 0.2e1 + t556 * t542 / 0.2e1 - t558 * t538 / 0.2e1) * t544) * pkin(5);
	t709 = -t497 / 0.2e1;
	t708 = t505 / 0.2e1;
	t703 = t566 / 0.2e1;
	t702 = t567 / 0.2e1;
	t696 = pkin(3) * t505;
	t695 = pkin(3) * t506;
	t571 = 0.1e1 / pkin(8);
	t692 = pkin(4) * t571;
	t666 = t540 * t566;
	t671 = t539 * t565;
	t607 = t666 + t671;
	t494 = ((t520 * t705 + t521 * t703) * t549 + t607 * t633) * t575;
	t529 = (t666 / 0.2e1 + t671 / 0.2e1) * t661;
	t527 = t529 ^ 2;
	t517 = t525 * t527 + 0.1e1;
	t672 = t525 * t529;
	t683 = t524 * t728;
	t686 = 0.2e1 * (t494 * t672 - t527 * t683) / t517 ^ 2;
	t684 = t485 * t573;
	t681 = t505 * t509;
	t662 = t549 * t566;
	t472 = (t723 * t731 + (qJD(3) * t597 + t496 * t702 + t569 * t709) * t549 + (-t723 * t556 - (-t674 + t675) * t555 + (qJD(3) * t605 + t678 - t679) * t558) * t651) * t575;
	t473 = (0.4e1 * t576 * t551 * t593 * t726 + (qJD(3) * t598 + t497 * t702 + t569 * t710) * t549 + (-t605 * t556 - (-t673 - t676) * t555 + (-qJD(3) * t723 + t677 + t680) * t558) * t651) * t575;
	t463 = -t472 * t563 - t473 * t562;
	t628 = -0.8e1 * t574 * t684;
	t453 = t463 * t603 + t488 * t628;
	t487 = -t500 * t563 + t501 * t562;
	t458 = (t471 * t640 - t487 * t577 + t488 * t634) * pkin(4);
	t650 = pkin(4) * t695;
	t631 = t480 * t650;
	t454 = (t458 * t708 + t488 * t631) * t571;
	t459 = t471 * t641 + t488 * t636 + (t487 * t504 - t488 * t577) * pkin(4);
	t630 = t481 * t650;
	t455 = (t459 * t708 + t488 * t630) * t571;
	t462 = t472 * t562 - t473 * t563;
	t498 = ((t522 * t705 + t523 * t703) * t549 + t607 * t632) * t575;
	t499 = ((t522 * t704 + t523 * t705) * t549 - t608 * t632) * t575;
	t515 = 0.1e1 / t517;
	t649 = 0.4e1 * t505 * t506 * t574;
	t589 = t571 * (t480 * t649 + 0.2e1 * t696) * t684;
	t590 = t521 * t558 + t523 * t555 - t540 * t556;
	t591 = t520 * t558 + t522 * t555 - t539 * t556;
	t594 = t565 * t731;
	t600 = t566 * t611;
	t612 = t468 * t626;
	t613 = t468 * t571 * t647;
	t642 = -t467 * t492 / t682 / 0.4e1;
	t622 = t609 * t642;
	t629 = t481 * t649;
	t637 = t573 * t652;
	t638 = t549 * t705;
	t657 = 0.2e1 * (((t453 * t711 + t471 * t642) * t509 + (-t463 * t609 - t488 * t610) * t637) * t708 + (-t487 * t696 + t488 * t629) * t684 + ((t462 * t504 - t463 * t577 + (t485 * t713 + t488 * t714) * t492) * t708 + (t457 * t488 + t459 * t485 + t463 * t481) * t695) * pkin(4)) * t613 - 0.2e1 * (t488 * t589 + ((-t462 * t577 + t471 * t622 - t463 * t504 + (t453 * t707 + t487 * t714 + t610 * t713) * t492) * t708 + (t463 * t681 + (t456 * t488 + t458 * t485 + t463 * t480) * t506) * pkin(3)) * t692) * t612 + 0.2e1 * t722 * t455 + 0.2e1 * t721 * t454 + (-t515 * t728 - t524 * t686) * t498 + (t672 * t686 + (-t494 * t525 + 0.2e1 * t529 * t683) * t515) * t499 + ((t497 * t638 + t539 * t594 + 0.4e1 * t540 * t600 + t662 * t710) * t524 - (t496 * t638 - 0.4e1 * t539 * t600 + t540 * t594 + t662 * t709) * t672 + ((t524 * t590 + t591 * t672) * t566 + (t524 * t591 - t590 * t672) * t565) * t651) * t515 * t575;
	t465 = (t486 * t640 - t512 * t577 + t609 * t634) * pkin(4);
	t460 = (t465 * t708 + t609 * t631) * t571;
	t466 = t486 * t641 + t609 * t636 + (t504 * t512 - t577 * t609) * pkin(4);
	t461 = (t466 * t708 + t609 * t630) * t571;
	t617 = 0.2e1 * t648;
	t653 = qJD(1) * (-t460 * t685 + t461 * t476) * t617;
	t620 = qJD(1) * ((-t454 * t685 + t455 * t476) * t617 + (t498 * t524 - t499 * t672) * t515 + 0.1e1);
	t570 = cos(qJ(1));
	t568 = sin(qJ(1));
	t464 = t603 * t610 + t609 * t628;
	t448 = 0.2e1 * (((t464 * t711 + t486 * t642) * t509 - 0.2e1 * t609 * t610 * t637) * t708 + (-t512 * t696 + t609 * t629) * t684 + ((-t485 * t504 - t724 + (t485 * t712 + t609 * t714) * t492) * t708 + (t457 * t609 + t466 * t485 + t481 * t610) * t695) * pkin(4)) * t613 - 0.2e1 * (t609 * t589 + ((t486 * t622 + (t464 * t707 + t512 * t714 + t610 * t712) * t492 - t730) * t708 + (t610 * t681 + (t456 * t609 + t465 * t485 + t480 * t610) * t506) * pkin(3)) * t692) * t612 + 0.2e1 * t722 * t461 + 0.2e1 * t721 * t460;
	t1 = [0, t568 * t657 + t570 * t620, t448 * t568 + t570 * t653, 0; 0, t568 * t620 - t570 * t657, -t448 * t570 + t568 * t653, 0; 0, 0, 0, 0;];
	JgD_rot = t1;
end