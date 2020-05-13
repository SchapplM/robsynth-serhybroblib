% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m1DE1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_jacobiR_rot_sym_varpar: pkin has to be [23x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:57
	% EndTime: 2020-04-14 18:42:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t12 = cos(qJ(1));
	t9 = sin(qJ(2));
	t15 = t12 * t9;
	t10 = sin(qJ(1));
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t8 = t10 * t9;
	t1 = [t8, -t13, 0, 0; -t15, -t14, 0, 0; 0, -t9, 0, 0; t14, t15, 0, 0; -t13, t8, 0, 0; 0, -t11, 0, 0; t12, 0, 0, 0; t10, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:42:58
	% EndTime: 2020-04-14 18:42:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (16->14), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->13)
	t34 = sin(qJ(3));
	t35 = sin(qJ(2));
	t37 = cos(qJ(3));
	t38 = cos(qJ(2));
	t44 = -t35 * t34 + t38 * t37;
	t41 = t38 * t34 + t35 * t37;
	t39 = cos(qJ(1));
	t36 = sin(qJ(1));
	t30 = t41 * t39;
	t29 = t44 * t39;
	t28 = t41 * t36;
	t27 = t44 * t36;
	t1 = [-t27, -t30, -t30, 0; t29, -t28, -t28, 0; 0, t44, t44, 0; t28, -t29, -t29, 0; -t30, -t27, -t27, 0; 0, -t41, -t41, 0; t39, 0, 0, 0; t36, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:44:14
	% EndTime: 2020-04-14 18:44:57
	% DurationCPUTime: 24.52s
	% Computational Cost: add. (712864->96), mult. (1075904->191), div. (47728->11), fcn. (681970->19), ass. (0->119)
	t383 = pkin(7) ^ 2;
	t377 = sin(qJ(2));
	t381 = cos(pkin(19));
	t443 = sin(pkin(19));
	t444 = cos(qJ(2));
	t363 = t377 * t381 - t444 * t443;
	t440 = pkin(7) * t363;
	t453 = -2 * pkin(1);
	t419 = (pkin(1) ^ 2) + t440 * t453;
	t355 = t383 + t419;
	t352 = pkin(3) ^ 2 - pkin(8) ^ 2 + t355;
	t360 = pkin(1) - t440;
	t365 = t377 * t443 + t444 * t381;
	t450 = -pkin(8) - pkin(3);
	t350 = (pkin(7) - t450) * (pkin(7) + t450) + t419;
	t449 = -pkin(8) + pkin(3);
	t351 = (pkin(7) - t449) * (pkin(7) + t449) + t419;
	t389 = sqrt(-t351 * t350);
	t439 = pkin(7) * t365;
	t418 = pkin(1) * t439;
	t434 = 0.2e1 / t389 * (t350 + t351) * t418;
	t338 = (t363 * t389 + (-t434 / 0.2e1 - t352 + t360 * t453) * t365) * pkin(7);
	t429 = t365 * t389;
	t339 = t360 * t434 / 0.2e1 + t383 * t365 ^ 2 * t453 + (-t363 * t352 - t429) * pkin(7);
	t353 = 0.1e1 / t355;
	t376 = sin(qJ(3));
	t386 = 0.1e1 / pkin(3);
	t408 = 0.1e1 / t355 ^ 2 * t418;
	t346 = t352 * t439 + t360 * t389;
	t379 = cos(qJ(3));
	t424 = t379 * t346;
	t345 = -pkin(7) * t429 + t360 * t352;
	t428 = t376 * t345;
	t445 = -t379 / 0.2e1;
	t324 = ((-t376 * t338 / 0.2e1 + t339 * t445) * t353 + (-t424 - t428) * t408) * t386;
	t425 = t379 * t345;
	t427 = t376 * t346;
	t325 = ((t338 * t445 + t376 * t339 / 0.2e1) * t353 + (-t425 + t427) * t408) * t386;
	t373 = pkin(23) + pkin(22);
	t371 = sin(t373);
	t372 = cos(t373);
	t320 = t372 * t324 + t371 * t325;
	t433 = t353 * t386;
	t341 = (t428 / 0.2e1 + t424 / 0.2e1) * t433;
	t342 = (-t425 / 0.2e1 + t427 / 0.2e1) * t433;
	t334 = t372 * t341 - t371 * t342;
	t442 = pkin(4) * t334;
	t452 = -2 * pkin(5);
	t420 = (pkin(5) ^ 2) - t442 * t452;
	t448 = -pkin(9) - pkin(11);
	t326 = (pkin(4) - t448) * (pkin(4) + t448) + t420;
	t447 = -pkin(9) + pkin(11);
	t327 = (pkin(4) - t447) * (pkin(4) + t447) + t420;
	t451 = pkin(4) * pkin(5);
	t391 = 0.2e1 * (t326 + t327) * t451;
	t316 = t320 * t391;
	t321 = -t371 * t324 + t372 * t325;
	t388 = sqrt(-t327 * t326);
	t385 = pkin(4) ^ 2;
	t331 = t385 + t420;
	t328 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t331;
	t332 = pkin(5) + t442;
	t409 = t332 * t452 - t328;
	t323 = 0.1e1 / t388;
	t402 = t371 * t341 + t372 * t342;
	t414 = -t323 * t402 / 0.2e1;
	t295 = (t316 * t414 + t409 * t320 - t321 * t388) * pkin(4);
	t410 = t385 * t402 * t452;
	t415 = t323 * t332 / 0.2e1;
	t296 = t316 * t415 + t320 * t410 + (-t320 * t388 + t321 * t328) * pkin(4);
	t441 = pkin(4) * t402;
	t317 = t332 * t328 - t388 * t441;
	t318 = t328 * t441 + t332 * t388;
	t375 = cos(pkin(21));
	t329 = 0.1e1 / t331;
	t382 = 0.1e1 / pkin(11);
	t435 = t329 * t382;
	t374 = sin(pkin(21));
	t446 = t374 / 0.2e1;
	t315 = (-t317 * t375 / 0.2e1 + t318 * t446) * t435;
	t312 = 0.1e1 / t315;
	t417 = 0.1e1 / t331 ^ 2 * t451;
	t406 = t375 * t417;
	t399 = t320 * t406;
	t407 = t374 * t417;
	t401 = t320 * t407;
	t436 = t329 * t375;
	t411 = t436 / 0.2e1;
	t412 = -t436 / 0.2e1;
	t413 = t329 * t446;
	t313 = 0.1e1 / t315 ^ 2;
	t314 = (t317 * t446 + t318 * t375 / 0.2e1) * t435;
	t437 = t313 * t314;
	t438 = 0.1e1 / (t313 * t314 ^ 2 + 0.1e1) * t382;
	t460 = ((t295 * t413 + t296 * t411 + t317 * t401 + t318 * t399) * t312 - (t295 * t412 + t296 * t413 - t317 * t399 + t318 * t401) * t437) * t438 + 0.1e1;
	t319 = t402 * t391;
	t310 = (t319 * t414 - t334 * t388 + t402 * t409) * pkin(4);
	t311 = t319 * t415 + t402 * t410 + (t334 * t328 - t388 * t402) * pkin(4);
	t398 = t402 * t406;
	t400 = t402 * t407;
	t459 = ((t310 * t413 + t311 * t411 + t317 * t400 + t318 * t398) * t312 - (t310 * t412 + t311 * t413 - t317 * t398 + t318 * t400) * t437) * t438 + 0.1e1;
	t309 = atan2(t314, t315);
	t306 = sin(t309);
	t307 = cos(t309);
	t364 = -t377 * t376 + t444 * t379;
	t378 = sin(qJ(1));
	t356 = t364 * t378;
	t390 = t444 * t376 + t377 * t379;
	t357 = t390 * t378;
	t458 = -t356 * t306 - t357 * t307;
	t457 = t357 * t306 - t356 * t307;
	t380 = cos(qJ(1));
	t358 = t364 * t380;
	t359 = t390 * t380;
	t456 = -t358 * t306 - t359 * t307;
	t455 = -t306 * t390 + t364 * t307;
	t454 = -t364 * t306 - t307 * t390;
	t395 = -t359 * t306 + t358 * t307;
	t1 = [t457, t460 * t456, t459 * t456, 0; t395, t460 * t458, t459 * t458, 0; 0, t460 * t455, t459 * t455, 0; -t458, -t460 * t395, -t459 * t395, 0; t456, t460 * t457, t459 * t457, 0; 0, t460 * t454, t459 * t454, 0; t380, 0, 0, 0; t378, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:47:52
	% EndTime: 2020-04-14 18:49:02
	% DurationCPUTime: 42.19s
	% Computational Cost: add. (1102906->110), mult. (1664338->213), div. (73992->11), fcn. (1054902->21), ass. (0->131)
	t638 = pkin(7) ^ 2;
	t631 = sin(qJ(2));
	t636 = cos(pkin(19));
	t695 = sin(pkin(19));
	t696 = cos(qJ(2));
	t615 = t631 * t636 - t696 * t695;
	t692 = pkin(7) * t615;
	t705 = -2 * pkin(1);
	t670 = (pkin(1) ^ 2) + t692 * t705;
	t607 = t638 + t670;
	t604 = pkin(3) ^ 2 - pkin(8) ^ 2 + t607;
	t612 = pkin(1) - t692;
	t617 = t631 * t695 + t696 * t636;
	t702 = -pkin(8) - pkin(3);
	t602 = (pkin(7) - t702) * (pkin(7) + t702) + t670;
	t701 = -pkin(8) + pkin(3);
	t603 = (pkin(7) - t701) * (pkin(7) + t701) + t670;
	t644 = sqrt(-t603 * t602);
	t691 = pkin(7) * t617;
	t669 = pkin(1) * t691;
	t686 = 0.2e1 / t644 * (t602 + t603) * t669;
	t590 = (t615 * t644 + (-t686 / 0.2e1 - t604 + t612 * t705) * t617) * pkin(7);
	t681 = t617 * t644;
	t591 = t612 * t686 / 0.2e1 + t638 * t617 ^ 2 * t705 + (-t615 * t604 - t681) * pkin(7);
	t605 = 0.1e1 / t607;
	t630 = sin(qJ(3));
	t641 = 0.1e1 / pkin(3);
	t659 = 0.1e1 / t607 ^ 2 * t669;
	t598 = t604 * t691 + t612 * t644;
	t634 = cos(qJ(3));
	t676 = t634 * t598;
	t597 = -pkin(7) * t681 + t612 * t604;
	t680 = t630 * t597;
	t697 = -t634 / 0.2e1;
	t576 = ((-t630 * t590 / 0.2e1 + t591 * t697) * t605 + (-t676 - t680) * t659) * t641;
	t677 = t634 * t597;
	t679 = t630 * t598;
	t577 = ((t590 * t697 + t630 * t591 / 0.2e1) * t605 + (-t677 + t679) * t659) * t641;
	t626 = pkin(23) + pkin(22);
	t624 = sin(t626);
	t625 = cos(t626);
	t572 = t625 * t576 + t624 * t577;
	t685 = t605 * t641;
	t593 = (t680 / 0.2e1 + t676 / 0.2e1) * t685;
	t594 = (-t677 / 0.2e1 + t679 / 0.2e1) * t685;
	t586 = t625 * t593 - t624 * t594;
	t694 = pkin(4) * t586;
	t704 = -2 * pkin(5);
	t671 = (pkin(5) ^ 2) - t694 * t704;
	t700 = -pkin(9) - pkin(11);
	t578 = (pkin(4) - t700) * (pkin(4) + t700) + t671;
	t699 = -pkin(9) + pkin(11);
	t579 = (pkin(4) - t699) * (pkin(4) + t699) + t671;
	t703 = pkin(4) * pkin(5);
	t646 = 0.2e1 * (t578 + t579) * t703;
	t568 = t572 * t646;
	t573 = -t624 * t576 + t625 * t577;
	t643 = sqrt(-t579 * t578);
	t640 = pkin(4) ^ 2;
	t583 = t640 + t671;
	t580 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t583;
	t584 = pkin(5) + t694;
	t660 = t584 * t704 - t580;
	t575 = 0.1e1 / t643;
	t656 = t624 * t593 + t625 * t594;
	t665 = -t575 * t656 / 0.2e1;
	t547 = (t568 * t665 + t660 * t572 - t573 * t643) * pkin(4);
	t661 = t656 * t640 * t704;
	t666 = t575 * t584 / 0.2e1;
	t548 = t568 * t666 + t572 * t661 + (-t572 * t643 + t573 * t580) * pkin(4);
	t693 = pkin(4) * t656;
	t569 = t584 * t580 - t643 * t693;
	t570 = t580 * t693 + t584 * t643;
	t628 = cos(pkin(21));
	t581 = 0.1e1 / t583;
	t637 = 0.1e1 / pkin(11);
	t687 = t581 * t637;
	t627 = sin(pkin(21));
	t698 = t627 / 0.2e1;
	t567 = (-t569 * t628 / 0.2e1 + t570 * t698) * t687;
	t564 = 0.1e1 / t567;
	t668 = 0.1e1 / t583 ^ 2 * t703;
	t657 = t628 * t668;
	t653 = t572 * t657;
	t658 = t627 * t668;
	t655 = t572 * t658;
	t688 = t581 * t628;
	t662 = t688 / 0.2e1;
	t663 = -t688 / 0.2e1;
	t664 = t581 * t698;
	t565 = 0.1e1 / t567 ^ 2;
	t566 = (t569 * t698 + t570 * t628 / 0.2e1) * t687;
	t689 = t565 * t566;
	t690 = 0.1e1 / (t566 ^ 2 * t565 + 0.1e1) * t637;
	t714 = ((t547 * t664 + t548 * t662 + t569 * t655 + t570 * t653) * t564 - (t547 * t663 + t548 * t664 - t569 * t653 + t570 * t655) * t689) * t690 + 0.1e1;
	t571 = t656 * t646;
	t562 = (t571 * t665 - t586 * t643 + t656 * t660) * pkin(4);
	t563 = t571 * t666 + t656 * t661 + (t586 * t580 - t643 * t656) * pkin(4);
	t652 = t656 * t657;
	t654 = t656 * t658;
	t713 = ((t562 * t664 + t563 * t662 + t569 * t654 + t570 * t652) * t564 - (t562 * t663 + t563 * t664 - t569 * t652 + t570 * t654) * t689) * t690 + 0.1e1;
	t561 = atan2(t566, t567);
	t558 = sin(t561);
	t559 = cos(t561);
	t616 = -t631 * t630 + t696 * t634;
	t632 = sin(qJ(1));
	t608 = t616 * t632;
	t645 = t696 * t630 + t631 * t634;
	t609 = t645 * t632;
	t543 = t609 * t558 - t608 * t559;
	t629 = sin(qJ(4));
	t633 = cos(qJ(4));
	t635 = cos(qJ(1));
	t712 = t543 * t633 + t635 * t629;
	t711 = t543 * t629 - t635 * t633;
	t710 = -t608 * t558 - t609 * t559;
	t610 = t616 * t635;
	t611 = t645 * t635;
	t709 = -t611 * t558 + t610 * t559;
	t708 = -t610 * t558 - t611 * t559;
	t707 = -t558 * t645 + t616 * t559;
	t706 = t616 * t558 + t559 * t645;
	t542 = t632 * t629 + t633 * t709;
	t541 = -t629 * t709 + t632 * t633;
	t538 = t713 * t707;
	t537 = t713 * t708;
	t536 = t713 * t710;
	t535 = t714 * t707;
	t534 = t714 * t708;
	t533 = t714 * t710;
	t1 = [t712, t534 * t633, t537 * t633, t541; t542, t533 * t633, t536 * t633, t711; 0, t535 * t633, t538 * t633, -t706 * t629; -t711, -t534 * t629, -t537 * t629, -t542; t541, -t533 * t629, -t536 * t629, t712; 0, -t535 * t629, -t538 * t629, -t706 * t633; t710, t714 * t709, t713 * t709, 0; -t708, -t714 * t543, -t713 * t543, 0; 0, t714 * t706, t713 * t706, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:00
	% EndTime: 2020-04-14 18:43:02
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (5667->42), mult. (8738->75), div. (380->7), fcn. (5634->13), ass. (0->51)
	t164 = -2 * pkin(7);
	t163 = (-pkin(8) - pkin(3));
	t162 = (-pkin(8) + pkin(3));
	t138 = cos(pkin(18));
	t161 = t138 / 0.2e1;
	t160 = cos(qJ(2));
	t159 = sin(pkin(19));
	t133 = sin(qJ(2));
	t137 = cos(pkin(19));
	t130 = t133 * t137 - t160 * t159;
	t158 = pkin(1) * t130;
	t131 = t133 * t159 + t160 * t137;
	t157 = pkin(1) * t131;
	t141 = pkin(1) ^ 2;
	t148 = t158 * t164 + t141;
	t122 = ((pkin(7) - t163) * (pkin(7) + t163)) + t148;
	t123 = ((pkin(7) - t162) * (pkin(7) + t162)) + t148;
	t142 = sqrt(-t123 * t122);
	t147 = pkin(7) * t157;
	t156 = 0.1e1 / t142 * (t122 + t123) * t147;
	t127 = (pkin(7) ^ 2) + t148;
	t125 = 0.1e1 / t127;
	t135 = sin(pkin(18));
	t155 = t125 * t135;
	t139 = 0.1e1 / pkin(8);
	t154 = t125 * t139;
	t153 = t131 * t142;
	t124 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t127;
	t128 = -pkin(7) + t158;
	t117 = -pkin(1) * t153 - t128 * t124;
	t118 = t124 * t157 - t128 * t142;
	t115 = (t117 * t161 - t135 * t118 / 0.2e1) * t154;
	t116 = (t118 * t161 + t117 * t135 / 0.2e1) * t154;
	t111 = atan2(t116, t115);
	t108 = sin(t111);
	t134 = sin(qJ(1));
	t152 = t134 * t108;
	t109 = cos(t111);
	t151 = t134 * t109;
	t136 = cos(qJ(1));
	t150 = t136 * t108;
	t149 = t136 * t109;
	t146 = t125 * t161;
	t145 = 0.1e1 / t127 ^ 2 * t147;
	t144 = t135 * t145;
	t143 = t138 * t145;
	t114 = 0.1e1 / t115 ^ 2;
	t113 = -t128 * t156 + t141 * t131 ^ 2 * t164 + (-t130 * t124 - t153) * pkin(1);
	t112 = (t130 * t142 + (0.2e1 * t128 * pkin(7) - t124 - t156) * t131) * pkin(1);
	t107 = ((t113 * t146 + t118 * t143 + t112 * t155 / 0.2e1 + t117 * t144) / t115 - (t112 * t146 + t117 * t143 - t113 * t155 / 0.2e1 - t118 * t144) * t116 * t114) / (t116 ^ 2 * t114 + 0.1e1) * t139;
	t1 = [-t151, -t107 * t150, 0, 0; t149, -t107 * t152, 0, 0; 0, t107 * t109, 0, 0; t152, -t107 * t149, 0, 0; -t150, -t107 * t151, 0, 0; 0, -t107 * t108, 0, 0; t136, 0, 0, 0; t134, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:01
	% EndTime: 2020-04-14 18:43:04
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (12662->46), mult. (19532->80), div. (856->7), fcn. (12626->13), ass. (0->56)
	t160 = pkin(7) ^ 2;
	t155 = sin(qJ(2));
	t157 = cos(qJ(2));
	t159 = cos(pkin(19));
	t186 = sin(pkin(19));
	t150 = t155 * t159 - t157 * t186;
	t185 = pkin(7) * t150;
	t190 = -2 * pkin(1);
	t174 = (pkin(1) ^ 2) + t185 * t190;
	t147 = t160 + t174;
	t144 = pkin(3) ^ 2 - pkin(8) ^ 2 + t147;
	t148 = pkin(1) - t185;
	t151 = t155 * t186 + t157 * t159;
	t189 = -pkin(8) - pkin(3);
	t142 = (pkin(7) - t189) * (pkin(7) + t189) + t174;
	t188 = -pkin(8) + pkin(3);
	t143 = (pkin(7) - t188) * (pkin(7) + t188) + t174;
	t163 = sqrt(-t143 * t142);
	t181 = t151 * t163;
	t137 = -pkin(7) * t181 + t148 * t144;
	t184 = pkin(7) * t151;
	t138 = t144 * t184 + t148 * t163;
	t154 = cos(pkin(23));
	t145 = 0.1e1 / t147;
	t161 = 0.1e1 / pkin(3);
	t182 = t145 * t161;
	t153 = sin(pkin(23));
	t187 = t153 / 0.2e1;
	t135 = (-t154 * t137 / 0.2e1 + t138 * t187) * t182;
	t136 = (t154 * t138 / 0.2e1 + t137 * t187) * t182;
	t131 = atan2(t136, t135);
	t128 = sin(t131);
	t129 = cos(t131);
	t168 = t157 * t128 + t155 * t129;
	t173 = pkin(1) * t184;
	t183 = 0.1e1 / t163 * (t142 + t143) * t173;
	t180 = t154 * t145;
	t179 = t155 * t128;
	t176 = t157 * t129;
	t156 = sin(qJ(1));
	t175 = t168 * t156;
	t172 = t145 * t187;
	t171 = 0.1e1 / t147 ^ 2 * t173;
	t170 = t153 * t171;
	t169 = t154 * t171;
	t167 = -t176 + t179;
	t132 = (t150 * t163 + (t148 * t190 - t144 - t183) * t151) * pkin(7);
	t133 = t148 * t183 + t160 * t151 ^ 2 * t190 + (-t150 * t144 - t181) * pkin(7);
	t134 = 0.1e1 / t135 ^ 2;
	t124 = ((t133 * t180 / 0.2e1 + t138 * t169 + t132 * t172 + t137 * t170) / t135 - (-t132 * t180 / 0.2e1 - t137 * t169 + t133 * t172 + t138 * t170) * t136 * t134) / (t136 ^ 2 * t134 + 0.1e1) * t161;
	t166 = t168 * t124;
	t165 = t167 * t124 - t176;
	t164 = t165 + t179;
	t158 = cos(qJ(1));
	t127 = t158 * t179;
	t1 = [t175, t165 * t158 + t127, 0, 0; -t168 * t158, t164 * t156, 0, 0; 0, (-t124 - 0.1e1) * t168, 0, 0; -t167 * t156, (t166 + t168) * t158, 0, 0; -t158 * t176 + t127, t156 * t166 + t175, 0, 0; 0, t164, 0, 0; t158, 0, 0, 0; t156, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:00
	% EndTime: 2020-04-14 18:43:02
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (6088->35), mult. (9208->67), div. (296->7), fcn. (5938->13), ass. (0->45)
	t125 = cos(qJ(1));
	t121 = sin(qJ(2));
	t124 = cos(qJ(2));
	t126 = pkin(6) ^ 2;
	t118 = sin(pkin(20));
	t119 = cos(pkin(20));
	t120 = sin(qJ(3));
	t123 = cos(qJ(3));
	t116 = t123 * t118 + t120 * t119;
	t138 = pkin(6) * t116;
	t134 = 0.2e1 * pkin(1) * t138 + t126;
	t113 = pkin(1) ^ 2 + t134;
	t110 = pkin(2) ^ 2 - pkin(13) ^ 2 + t113;
	t114 = -pkin(1) - t138;
	t117 = t120 * t118 - t123 * t119;
	t143 = -pkin(2) - pkin(13);
	t108 = (pkin(1) - t143) * (pkin(1) + t143) + t134;
	t142 = -pkin(2) + pkin(13);
	t109 = (pkin(1) - t142) * (pkin(1) + t142) + t134;
	t129 = sqrt(-t109 * t108);
	t135 = t117 * t129;
	t103 = -pkin(6) * t135 - t114 * t110;
	t137 = t117 * pkin(6);
	t104 = t110 * t137 - t114 * t129;
	t127 = 0.1e1 / pkin(2);
	t111 = 0.1e1 / t113;
	t141 = t111 / 0.2e1;
	t133 = t127 * t141;
	t100 = atan2(t104 * t133, t103 * t133);
	t98 = sin(t100);
	t99 = cos(t100);
	t131 = t121 * t99 + t124 * t98;
	t97 = t131 * t125;
	t140 = pkin(1) * t117;
	t139 = pkin(6) * t104;
	t136 = 0.1e1 / t129 * (t108 + t109) * pkin(1) * t137;
	t130 = t121 * t98 - t124 * t99;
	t122 = sin(qJ(1));
	t95 = t131 * t122;
	t94 = t130 * t122;
	t96 = t130 * t125;
	t112 = 0.1e1 / t113 ^ 2;
	t102 = 0.1e1 / t103 ^ 2;
	t93 = 0.2e1 * (((-t114 * t136 + (t116 * t110 - t135) * pkin(6)) * t141 + (-t111 * t126 * t117 + t112 * t139) * t140) / t103 - ((-t116 * t129 + (-t110 - t136) * t117) * t141 + (t103 * t112 + t111 * t114) * t140) * t102 * t139) * pkin(2) / (t104 ^ 2 * t102 + 0.1e1) * t113 * t127;
	t1 = [t95, t96, t93 * t96, 0; -t97, t94, t93 * t94, 0; 0, -t131, -t131 * t93, 0; -t94, t97, t93 * t97, 0; t96, t95, t93 * t95, 0; 0, t130, t130 * t93, 0; t125, 0, 0, 0; t122, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:43:03
	% EndTime: 2020-04-14 18:43:07
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (15580->56), mult. (23288->106), div. (936->11), fcn. (14778->16), ass. (0->72)
	t186 = sin(pkin(20));
	t187 = cos(pkin(20));
	t188 = sin(qJ(3));
	t191 = cos(qJ(3));
	t184 = t191 * t186 + t188 * t187;
	t224 = pkin(6) * t184;
	t213 = pkin(1) * t224;
	t183 = 0.2e1 * t213;
	t197 = pkin(2) ^ 2;
	t196 = pkin(6) ^ 2;
	t214 = pkin(1) ^ 2 + t196;
	t211 = -pkin(13) ^ 2 + t214;
	t178 = t183 + t197 + t211;
	t182 = -pkin(1) - t224;
	t185 = t188 * t186 - t191 * t187;
	t215 = t183 + t196;
	t230 = -pkin(2) - pkin(13);
	t174 = (pkin(1) - t230) * (pkin(1) + t230) + t215;
	t229 = -pkin(2) + pkin(13);
	t175 = (pkin(1) - t229) * (pkin(1) + t229) + t215;
	t221 = t175 * t174;
	t200 = sqrt(-t221);
	t220 = t185 * t200;
	t212 = pkin(6) * t220;
	t165 = -t182 * t178 - t212;
	t164 = 0.1e1 / t165 ^ 2;
	t223 = t185 * pkin(6);
	t166 = t178 * t223 - t182 * t200;
	t181 = t183 + t214;
	t179 = 0.1e1 / t181;
	t180 = 0.1e1 / t181 ^ 2;
	t198 = 0.1e1 / pkin(2);
	t222 = 0.1e1 / t200 * (t174 + t175) * pkin(1) * t223;
	t225 = pkin(6) * t166;
	t226 = pkin(1) * t185;
	t228 = t179 / 0.2e1;
	t148 = 0.2e1 * (((-t182 * t222 + (t184 * t178 - t220) * pkin(6)) * t228 + (-t179 * t196 * t185 + t180 * t225) * t226) / t165 - ((-t184 * t200 + (-t178 - t222) * t185) * t228 + (t165 * t180 + t179 * t182) * t226) * t164 * t225) * pkin(2) / (t166 ^ 2 * t164 + 0.1e1) * t181 * t198;
	t227 = t198 / 0.2e1;
	t210 = t179 * t227;
	t162 = atan2(t166 * t210, t165 * t210);
	t160 = sin(t162);
	t161 = cos(t162);
	t189 = sin(qJ(2));
	t192 = cos(qJ(2));
	t201 = t192 * t160 + t189 * t161;
	t146 = t201 * t148;
	t153 = t189 * t160 - t192 * t161;
	t209 = 0.1e1 / pkin(13) * t227;
	t190 = sin(qJ(1));
	t149 = t153 * t190;
	t150 = t201 * t190;
	t177 = t197 - t211 - 0.2e1 * t213;
	t169 = atan2(t200 * t209, t177 * t209);
	t167 = sin(t169);
	t168 = cos(t169);
	t207 = t149 * t168 + t150 * t167;
	t206 = t149 * t167 - t150 * t168;
	t193 = cos(qJ(1));
	t151 = t153 * t193;
	t152 = t201 * t193;
	t205 = -t151 * t168 - t152 * t167;
	t204 = t151 * t167 - t152 * t168;
	t203 = -t153 * t168 - t167 * t201;
	t202 = -t153 * t167 + t168 * t201;
	t147 = t153 * t148;
	t176 = 0.1e1 / t177 ^ 2;
	t155 = (0.1e1 / t177 * t222 - 0.2e1 * pkin(1) * t176 * t212) / (-t176 * t221 + 0.1e1);
	t145 = t193 * t146;
	t144 = t148 * t151;
	t143 = t148 * t150;
	t142 = t190 * t147;
	t1 = [t206, t205, -t144 * t168 - t145 * t167 + t205 * t155, 0; -t204, -t207, -t142 * t168 - t143 * t167 - t207 * t155, 0; 0, t202, t146 * t168 - t147 * t167 + t202 * t155, 0; t207, t204, t144 * t167 - t145 * t168 + t204 * t155, 0; t205, t206, t142 * t167 - t143 * t168 + t206 * t155, 0; 0, t203, -t146 * t167 - t147 * t168 + t203 * t155, 0; t193, 0, 0, 0; t190, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-14 18:44:18
	% EndTime: 2020-04-14 18:45:11
	% DurationCPUTime: 30.07s
	% Computational Cost: add. (682060->138), mult. (1030012->248), div. (44976->15), fcn. (653498->24), ass. (0->145)
	t422 = sin(qJ(2));
	t423 = sin(pkin(19));
	t425 = cos(qJ(2));
	t426 = cos(pkin(19));
	t404 = t422 * t426 - t425 * t423;
	t403 = pkin(7) * t404;
	t399 = (-0.2e1 * t403 + pkin(1)) * pkin(1);
	t429 = -pkin(8) + pkin(3);
	t430 = -pkin(8) - pkin(3);
	t334 = sqrt(-((pkin(7) - t429) * (pkin(7) + t429) + t399) * ((pkin(7) - t430) * (pkin(7) + t430) + t399));
	t437 = pkin(7) ^ 2;
	t325 = t399 + t437;
	t335 = pkin(3) ^ 2;
	t441 = -pkin(8) ^ 2 + t335;
	t397 = t325 + t441;
	t401 = -t403 + pkin(1);
	t326 = t422 * t423 + t425 * t426;
	t419 = pkin(7) * t326;
	t439 = 0.1e1 / pkin(3);
	t393 = t439 * (-t334 * t419 + t401 * t397);
	t433 = 0.1e1 / t325;
	t391 = t433 * t393;
	t389 = t391 / 0.2e1;
	t396 = pkin(7) * t397;
	t394 = t439 * (t326 * t396 + t401 * t334);
	t392 = t433 * t394;
	t390 = t392 / 0.2e1;
	t421 = sin(qJ(3));
	t424 = cos(qJ(3));
	t385 = t421 * t389 + t424 * t390;
	t386 = t424 * t389 - t421 * t392 / 0.2e1;
	t415 = pkin(23) + pkin(22);
	t409 = sin(t415);
	t410 = cos(t415);
	t377 = t410 * t385 + t409 * t386;
	t374 = pkin(5) * t377;
	t434 = -2 * pkin(4);
	t370 = (-t377 * t434 + pkin(5)) * pkin(5);
	t428 = (-pkin(9) - pkin(11));
	t361 = ((pkin(4) - t428) * (pkin(4) + t428)) + t370;
	t427 = (-pkin(9) + pkin(11));
	t362 = ((pkin(4) - t427) * (pkin(4) + t427)) + t370;
	t333 = sqrt(-t362 * t361);
	t440 = t409 * t385 - t410 * t386;
	t450 = pkin(5) * t440;
	t451 = t333 * t450;
	t336 = pkin(9) ^ 2;
	t438 = pkin(5) ^ 2;
	t369 = t438 + (0.2e1 * t374 + pkin(4)) * pkin(4);
	t365 = -pkin(11) ^ 2 + t336 + t369;
	t363 = pkin(5) * t365;
	t449 = t363 * t440;
	t416 = pkin(1) * t419;
	t417 = 0.4e1 / t334 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t399 - t335) * t416;
	t317 = (t404 * t334 + (-t417 / 0.2e1 - t437 - t441) * t326) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t403) * t416;
	t318 = t401 * t417 / 0.2e1 - t404 * t396 + (-t334 - 0.2e1 * t416) * t419;
	t328 = sin(pkin(23));
	t330 = cos(pkin(23));
	t320 = -t330 * t391 / 0.2e1 + t328 * t390;
	t319 = 0.1e1 / t320 ^ 2;
	t321 = t328 * t389 + t330 * t390;
	t324 = 0.1e1 / t325 ^ 2;
	t387 = t393 * t416;
	t388 = t394 * t416;
	t420 = t439 * t433;
	t413 = t420 / 0.2e1;
	t406 = t328 * t413;
	t407 = t318 * t413;
	t414 = -t420 / 0.2e1;
	t408 = t317 * t414;
	t304 = ((t330 * t407 + t317 * t406 + (t328 * t387 + t330 * t388) * t324) / t320 - (t330 * t408 + t318 * t406 + (t328 * t388 - t330 * t387) * t324) * t321 * t319) / (t319 * t321 ^ 2 + 0.1e1);
	t316 = atan2(t321, t320);
	t313 = sin(t316);
	t314 = cos(t316);
	t402 = t425 * t313 + t422 * t314;
	t448 = (t304 + 0.1e1) * t402;
	t371 = -t374 - pkin(4);
	t447 = 0.2e1 * t371;
	t367 = 0.1e1 / t369;
	t446 = t367 / 0.2e1;
	t445 = -t371 / 0.2e1;
	t444 = -t450 / 0.2e1;
	t431 = pkin(4) * pkin(5);
	t442 = 0.1e1 / t369 ^ 2 * t431;
	t435 = 0.1e1 / pkin(9);
	t432 = 0.1e1 / t333;
	t418 = cos(pkin(22));
	t412 = t425 * t314;
	t411 = t422 * t313;
	t405 = t424 * t414;
	t310 = t412 - t411;
	t400 = t310 * t304;
	t398 = t411 - t400;
	t382 = t317 * t405 + t421 * t407 + (-t424 * t387 + t421 * t388) * t324;
	t381 = t421 * t408 + t318 * t405 + (-t421 * t387 - t424 * t388) * t324;
	t373 = pkin(4) * t450;
	t372 = t438 * t440 * t434;
	t366 = t435 * t367;
	t364 = t366 / 0.2e1;
	t360 = t410 * t381 + t409 * t382;
	t359 = pkin(5) * (-t409 * t381 + t410 * t382);
	t358 = pkin(5) * t360;
	t357 = pkin(4) * t358;
	t356 = -t371 * t333 + t449;
	t355 = -t371 * t365 - t451;
	t354 = 0.1e1 / t355 ^ 2;
	t353 = t435 * t356;
	t352 = t435 * t355;
	t351 = 0.4e1 * t432 * (((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t370) * t373 - t440 * t336 * t431);
	t350 = t353 * t442;
	t349 = t352 * t442;
	t348 = 0.2e1 * t432 * (t361 + t362) * t357;
	t347 = 0.1e1 / (t356 ^ 2 * t354 + 0.1e1) / t366;
	t346 = atan2(t353 * t446, t352 * t446);
	t345 = sin(t346);
	t344 = 0.2e1 / t355 * t347;
	t343 = -0.2e1 * t354 * t356 * t347;
	t342 = ((t351 * t445 + t377 * t363 + t372 * t440 - t451) * t364 + t440 * t350) * t344 + ((-t333 * t374 + t351 * t444 + t373 * t447 - t449) * t364 + t440 * t349) * t343;
	t303 = cos(t346);
	t341 = t303 * t342;
	t340 = ((-t333 * t358 + t348 * t445 + t365 * t359 + t360 * t372) * t364 + t360 * t350) * t344 + ((-t333 * t359 + t348 * t444 + t357 * t447 - t365 * t358) * t364 + t360 * t349) * t343;
	t339 = t303 * t340;
	t338 = t342 * t345;
	t337 = t340 * t345;
	t332 = cos(qJ(1));
	t331 = sin(qJ(1));
	t329 = sin(pkin(22));
	t312 = t332 * t411;
	t311 = t331 * t412;
	t308 = t402 * t332;
	t307 = -t332 * t412 + t312;
	t306 = t402 * t331;
	t305 = t331 * t411 - t311;
	t300 = -t412 + t398;
	t299 = t312 + (-t412 - t400) * t332;
	t298 = t448 * t332;
	t297 = t398 * t331 - t311;
	t296 = t448 * t331;
	t295 = -t418 * t303 - t329 * t345;
	t294 = t303 * t329 - t418 * t345;
	t289 = -t329 * t341 + t418 * t338;
	t288 = -t329 * t338 - t418 * t341;
	t287 = -t329 * t339 + t418 * t337;
	t286 = -t329 * t337 - t418 * t339;
	t1 = [-t294 * t305 + t295 * t306, t286 * t307 - t287 * t308 + t294 * t298 + t295 * t299, t288 * t307 - t289 * t308, 0; t294 * t307 - t295 * t308, t286 * t305 - t287 * t306 + t294 * t296 + t295 * t297, t288 * t305 - t289 * t306, 0; 0, -t286 * t402 + t287 * t310 + t294 * t300 - t295 * t448, -t288 * t402 + t289 * t310, 0; -t294 * t306 - t295 * t305, t286 * t308 + t287 * t307 - t294 * t299 + t295 * t298, t288 * t308 + t289 * t307, 0; t294 * t308 + t295 * t307, t286 * t306 + t287 * t305 - t294 * t297 + t295 * t296, t288 * t306 + t289 * t305, 0; 0, -t286 * t310 - t287 * t402 + t294 * t448 + t295 * t300, -t288 * t310 - t289 * t402, 0; t332, 0, 0, 0; t331, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end