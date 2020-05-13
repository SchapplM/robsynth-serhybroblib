% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m1DE1
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
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m1DE1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_jacobiR_rot_sym_varpar: pkin has to be [19x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
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
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0; t12, -t15, 0, 0; 0, t10, 0, 0; t15, -t12, 0, 0; -t14, -t13, 0, 0; 0, -t8, 0, 0; t11, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:14
	% EndTime: 2020-04-19 18:18:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (16->8), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->13)
	t33 = sin(qJ(3));
	t34 = sin(qJ(2));
	t36 = cos(qJ(3));
	t37 = cos(qJ(2));
	t32 = t37 * t33 + t34 * t36;
	t31 = t34 * t33 - t37 * t36;
	t38 = cos(qJ(1));
	t35 = sin(qJ(1));
	t30 = t32 * t38;
	t29 = t31 * t38;
	t28 = t31 * t35;
	t27 = t32 * t35;
	t1 = [-t28, t30, t30, 0; t29, t27, t27, 0; 0, t31, t31, 0; -t27, -t29, -t29, 0; t30, -t28, -t28, 0; 0, t32, t32, 0; t38, 0, 0, 0; t35, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:19:11
	% EndTime: 2020-04-19 18:19:50
	% DurationCPUTime: 24.71s
	% Computational Cost: add. (712864->98), mult. (1075904->195), div. (47728->11), fcn. (681970->19), ass. (0->121)
	t378 = pkin(5) ^ 2;
	t374 = sin(pkin(16));
	t376 = cos(qJ(2));
	t440 = sin(qJ(2));
	t442 = cos(pkin(16));
	t359 = t440 * t374 - t376 * t442;
	t437 = pkin(5) * t359;
	t451 = -2 * pkin(1);
	t416 = (pkin(1) ^ 2) + t437 * t451;
	t351 = t378 + t416;
	t348 = pkin(2) ^ 2 - pkin(6) ^ 2 + t351;
	t356 = pkin(1) - t437;
	t361 = t376 * t374 + t440 * t442;
	t448 = -pkin(6) - pkin(2);
	t346 = (pkin(5) - t448) * (pkin(5) + t448) + t416;
	t447 = -pkin(6) + pkin(2);
	t347 = (pkin(5) - t447) * (pkin(5) + t447) + t416;
	t384 = sqrt(-t347 * t346);
	t436 = pkin(5) * t361;
	t415 = pkin(1) * t436;
	t427 = 0.2e1 / t384 * (t346 + t347) * t415;
	t335 = (t359 * t384 + (-t427 / 0.2e1 - t348 + t356 * t451) * t361) * pkin(5);
	t422 = t361 * t384;
	t336 = t356 * t427 / 0.2e1 + t378 * t361 ^ 2 * t451 + (-t359 * t348 - t422) * pkin(5);
	t349 = 0.1e1 / t351;
	t375 = cos(qJ(3));
	t381 = 0.1e1 / pkin(2);
	t404 = 0.1e1 / t351 ^ 2 * t415;
	t342 = t348 * t436 + t356 * t384;
	t428 = t342 * t375;
	t341 = -pkin(5) * t422 + t356 * t348;
	t372 = sin(qJ(3));
	t431 = t341 * t372;
	t443 = t372 / 0.2e1;
	t321 = ((t336 * t375 / 0.2e1 + t335 * t443) * t349 + (t428 + t431) * t404) * t381;
	t429 = t342 * t372;
	t430 = t341 * t375;
	t322 = ((-t335 * t375 / 0.2e1 + t336 * t443) * t349 + (t429 - t430) * t404) * t381;
	t369 = pkin(18) + pkin(19);
	t367 = sin(t369);
	t368 = cos(t369);
	t318 = -t367 * t321 - t368 * t322;
	t426 = t349 * t381;
	t337 = (-t430 / 0.2e1 + t429 / 0.2e1) * t426;
	t338 = (t428 / 0.2e1 + t431 / 0.2e1) * t426;
	t332 = t368 * t337 + t367 * t338;
	t438 = pkin(3) * t332;
	t450 = -2 * pkin(4);
	t417 = (pkin(4) ^ 2) - t438 * t450;
	t446 = -pkin(8) - pkin(10);
	t323 = (pkin(3) - t446) * (pkin(3) + t446) + t417;
	t445 = -pkin(8) + pkin(10);
	t324 = (pkin(3) - t445) * (pkin(3) + t445) + t417;
	t449 = pkin(3) * pkin(4);
	t386 = 0.2e1 * (t323 + t324) * t449;
	t313 = t318 * t386;
	t317 = -t368 * t321 + t367 * t322;
	t383 = sqrt(-t324 * t323);
	t380 = pkin(3) ^ 2;
	t328 = t380 + t417;
	t325 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t328;
	t329 = pkin(4) + t438;
	t405 = t329 * t450 - t325;
	t320 = 0.1e1 / t383;
	t387 = t367 * t337 - t368 * t338;
	t410 = -t320 * t387 / 0.2e1;
	t292 = (t313 * t410 - t317 * t383 + t405 * t318) * pkin(3);
	t406 = t380 * t387 * t450;
	t411 = t320 * t329 / 0.2e1;
	t293 = t313 * t411 + t318 * t406 + (t317 * t325 - t318 * t383) * pkin(3);
	t439 = pkin(3) * t387;
	t314 = t329 * t325 - t383 * t439;
	t315 = t325 * t439 + t329 * t383;
	t371 = cos(pkin(17));
	t326 = 0.1e1 / t328;
	t377 = 0.1e1 / pkin(10);
	t432 = t326 * t377;
	t370 = sin(pkin(17));
	t444 = t370 / 0.2e1;
	t311 = (-t314 * t371 / 0.2e1 + t315 * t444) * t432;
	t309 = 0.1e1 / t311;
	t414 = 0.1e1 / t328 ^ 2 * t449;
	t402 = t371 * t414;
	t395 = t318 * t402;
	t403 = t370 * t414;
	t397 = t318 * t403;
	t433 = t326 * t371;
	t407 = t433 / 0.2e1;
	t408 = -t433 / 0.2e1;
	t409 = t326 * t444;
	t310 = 0.1e1 / t311 ^ 2;
	t312 = (t315 * t371 / 0.2e1 + t314 * t444) * t432;
	t434 = t310 * t312;
	t435 = 0.1e1 / (t312 ^ 2 * t310 + 0.1e1) * t377;
	t458 = ((t292 * t409 + t293 * t407 + t314 * t397 + t315 * t395) * t309 - (t292 * t408 + t293 * t409 - t314 * t395 + t315 * t397) * t434) * t435 + 0.1e1;
	t316 = t387 * t386;
	t307 = (t316 * t410 - t332 * t383 + t387 * t405) * pkin(3);
	t308 = t316 * t411 + t387 * t406 + (t332 * t325 - t383 * t387) * pkin(3);
	t394 = t387 * t402;
	t396 = t387 * t403;
	t457 = ((t307 * t409 + t308 * t407 + t314 * t396 + t315 * t394) * t309 - (t307 * t408 + t308 * t409 - t314 * t394 + t315 * t396) * t434) * t435 + 0.1e1;
	t306 = atan2(t312, t311);
	t303 = sin(t306);
	t304 = cos(t306);
	t373 = sin(qJ(1));
	t385 = t376 * t372 + t440 * t375;
	t352 = t385 * t373;
	t358 = t440 * t372 - t376 * t375;
	t353 = t358 * t373;
	t456 = -t353 * t303 + t352 * t304;
	t455 = -t352 * t303 - t353 * t304;
	t441 = cos(qJ(1));
	t398 = t441 * t440;
	t413 = t376 * t441;
	t354 = -t372 * t398 + t375 * t413;
	t355 = t372 * t413 + t375 * t398;
	t454 = t354 * t303 + t355 * t304;
	t453 = t303 * t385 + t358 * t304;
	t452 = -t358 * t303 + t304 * t385;
	t391 = t355 * t303 - t354 * t304;
	t1 = [t455, t458 * t454, t457 * t454, 0; t391, t458 * t456, t457 * t456, 0; 0, t458 * t453, t457 * t453, 0; -t456, -t458 * t391, -t457 * t391, 0; t454, t458 * t455, t457 * t455, 0; 0, t458 * t452, t457 * t452, 0; t441, 0, 0, 0; t373, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:22:46
	% EndTime: 2020-04-19 18:23:56
	% DurationCPUTime: 42.83s
	% Computational Cost: add. (1102906->108), mult. (1664338->213), div. (73992->11), fcn. (1054902->21), ass. (0->131)
	t638 = pkin(5) ^ 2;
	t632 = sin(pkin(16));
	t635 = cos(qJ(2));
	t697 = sin(qJ(2));
	t698 = cos(pkin(16));
	t616 = t632 * t697 - t635 * t698;
	t693 = pkin(5) * t616;
	t706 = -2 * pkin(1);
	t673 = (pkin(1) ^ 2) + t693 * t706;
	t608 = t638 + t673;
	t605 = pkin(2) ^ 2 - pkin(6) ^ 2 + t608;
	t613 = pkin(1) - t693;
	t618 = t635 * t632 + t697 * t698;
	t703 = -pkin(6) - pkin(2);
	t603 = (pkin(5) - t703) * (pkin(5) + t703) + t673;
	t702 = -pkin(6) + pkin(2);
	t604 = (pkin(5) - t702) * (pkin(5) + t702) + t673;
	t644 = sqrt(-t604 * t603);
	t692 = pkin(5) * t618;
	t672 = pkin(1) * t692;
	t685 = 0.2e1 / t644 * (t603 + t604) * t672;
	t592 = (t616 * t644 + (-t685 / 0.2e1 - t605 + t613 * t706) * t618) * pkin(5);
	t680 = t618 * t644;
	t593 = t613 * t685 / 0.2e1 + t638 * t618 ^ 2 * t706 + (-t616 * t605 - t680) * pkin(5);
	t606 = 0.1e1 / t608;
	t634 = cos(qJ(3));
	t641 = 0.1e1 / pkin(2);
	t660 = 0.1e1 / t608 ^ 2 * t672;
	t598 = -pkin(5) * t680 + t613 * t605;
	t696 = sin(qJ(3));
	t669 = t598 * t696;
	t670 = t696 / 0.2e1;
	t599 = t605 * t692 + t613 * t644;
	t686 = t599 * t634;
	t578 = ((t593 * t634 / 0.2e1 + t592 * t670) * t606 + (t669 + t686) * t660) * t641;
	t668 = t599 * t696;
	t687 = t598 * t634;
	t579 = ((-t592 * t634 / 0.2e1 + t593 * t670) * t606 + (t668 - t687) * t660) * t641;
	t627 = pkin(18) + pkin(19);
	t625 = sin(t627);
	t626 = cos(t627);
	t575 = -t625 * t578 - t626 * t579;
	t684 = t606 * t641;
	t594 = (-t687 / 0.2e1 + t668 / 0.2e1) * t684;
	t595 = (t686 / 0.2e1 + t669 / 0.2e1) * t684;
	t589 = t626 * t594 + t625 * t595;
	t694 = pkin(3) * t589;
	t705 = -2 * pkin(4);
	t674 = (pkin(4) ^ 2) - t694 * t705;
	t701 = -pkin(8) - pkin(10);
	t580 = (pkin(3) - t701) * (pkin(3) + t701) + t674;
	t700 = -pkin(8) + pkin(10);
	t581 = (pkin(3) - t700) * (pkin(3) + t700) + t674;
	t704 = pkin(3) * pkin(4);
	t645 = 0.2e1 * (t580 + t581) * t704;
	t570 = t575 * t645;
	t574 = -t626 * t578 + t625 * t579;
	t643 = sqrt(-t581 * t580);
	t640 = pkin(3) ^ 2;
	t585 = t640 + t674;
	t582 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t585;
	t586 = pkin(4) + t694;
	t661 = t586 * t705 - t582;
	t577 = 0.1e1 / t643;
	t646 = t625 * t594 - t626 * t595;
	t666 = -t577 * t646 / 0.2e1;
	t549 = (t570 * t666 - t574 * t643 + t575 * t661) * pkin(3);
	t662 = t640 * t646 * t705;
	t667 = t577 * t586 / 0.2e1;
	t550 = t570 * t667 + t575 * t662 + (t574 * t582 - t575 * t643) * pkin(3);
	t695 = pkin(3) * t646;
	t571 = t586 * t582 - t643 * t695;
	t572 = t582 * t695 + t586 * t643;
	t629 = cos(pkin(17));
	t583 = 0.1e1 / t585;
	t637 = 0.1e1 / pkin(10);
	t688 = t583 * t637;
	t628 = sin(pkin(17));
	t699 = t628 / 0.2e1;
	t568 = (-t571 * t629 / 0.2e1 + t572 * t699) * t688;
	t566 = 0.1e1 / t568;
	t671 = 0.1e1 / t585 ^ 2 * t704;
	t658 = t629 * t671;
	t654 = t575 * t658;
	t659 = t628 * t671;
	t656 = t575 * t659;
	t689 = t583 * t629;
	t663 = t689 / 0.2e1;
	t664 = -t689 / 0.2e1;
	t665 = t583 * t699;
	t567 = 0.1e1 / t568 ^ 2;
	t569 = (t572 * t629 / 0.2e1 + t571 * t699) * t688;
	t690 = t567 * t569;
	t691 = 0.1e1 / (t567 * t569 ^ 2 + 0.1e1) * t637;
	t716 = ((t549 * t665 + t550 * t663 + t571 * t656 + t572 * t654) * t566 - (t549 * t664 + t550 * t665 - t571 * t654 + t572 * t656) * t690) * t691 + 0.1e1;
	t573 = t646 * t645;
	t564 = (t573 * t666 - t589 * t643 + t646 * t661) * pkin(3);
	t565 = t573 * t667 + t646 * t662 + (t589 * t582 - t643 * t646) * pkin(3);
	t653 = t646 * t658;
	t655 = t646 * t659;
	t715 = ((t564 * t665 + t565 * t663 + t571 * t655 + t572 * t653) * t566 - (t564 * t664 + t565 * t665 - t571 * t653 + t572 * t655) * t690) * t691 + 0.1e1;
	t630 = sin(qJ(4));
	t633 = cos(qJ(4));
	t636 = cos(qJ(1));
	t563 = atan2(t569, t568);
	t560 = sin(t563);
	t561 = cos(t563);
	t617 = -t634 * t697 - t635 * t696;
	t631 = sin(qJ(1));
	t609 = t617 * t631;
	t615 = -t635 * t634 + t697 * t696;
	t610 = t615 * t631;
	t712 = -t609 * t560 + t610 * t561;
	t714 = t636 * t630 - t633 * t712;
	t713 = t630 * t712 + t636 * t633;
	t711 = -t610 * t560 - t609 * t561;
	t611 = t615 * t636;
	t612 = t617 * t636;
	t710 = -t611 * t560 - t612 * t561;
	t709 = -t612 * t560 + t611 * t561;
	t708 = t615 * t560 + t617 * t561;
	t707 = -t617 * t560 + t615 * t561;
	t544 = t631 * t630 + t633 * t709;
	t543 = -t630 * t709 + t631 * t633;
	t540 = t715 * t707;
	t539 = t715 * t710;
	t538 = t715 * t711;
	t537 = t716 * t707;
	t536 = t716 * t710;
	t535 = t716 * t711;
	t1 = [t714, t536 * t633, t539 * t633, t543; t544, t535 * t633, t538 * t633, -t713; 0, t537 * t633, t540 * t633, -t708 * t630; t713, -t536 * t630, -t539 * t630, -t544; t543, -t535 * t630, -t538 * t630, t714; 0, -t537 * t630, -t540 * t630, -t708 * t633; t711, t716 * t709, t715 * t709, 0; -t710, t716 * t712, t715 * t712, 0; 0, t716 * t708, t715 * t708, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:16
	% EndTime: 2020-04-19 18:18:16
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (5667->42), mult. (8738->75), div. (380->7), fcn. (5634->13), ass. (0->51)
	t164 = -2 * pkin(5);
	t163 = (-pkin(6) - pkin(2));
	t162 = (-pkin(6) + pkin(2));
	t138 = cos(pkin(15));
	t161 = t138 / 0.2e1;
	t160 = cos(pkin(16));
	t159 = cos(qJ(2));
	t133 = sin(qJ(2));
	t135 = sin(pkin(16));
	t130 = t133 * t135 - t159 * t160;
	t158 = pkin(1) * t130;
	t131 = t133 * t160 + t159 * t135;
	t157 = pkin(1) * t131;
	t141 = pkin(1) ^ 2;
	t148 = t158 * t164 + t141;
	t122 = ((pkin(5) - t163) * (pkin(5) + t163)) + t148;
	t123 = ((pkin(5) - t162) * (pkin(5) + t162)) + t148;
	t142 = sqrt(-t123 * t122);
	t147 = pkin(5) * t157;
	t156 = 0.1e1 / t142 * (t122 + t123) * t147;
	t127 = (pkin(5) ^ 2) + t148;
	t125 = 0.1e1 / t127;
	t136 = sin(pkin(15));
	t155 = t125 * t136;
	t139 = 0.1e1 / pkin(6);
	t154 = t125 * t139;
	t153 = t131 * t142;
	t124 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t127;
	t128 = -pkin(5) + t158;
	t117 = -pkin(1) * t153 - t128 * t124;
	t118 = t124 * t157 - t128 * t142;
	t115 = (t117 * t161 + t118 * t136 / 0.2e1) * t154;
	t116 = (t118 * t161 - t117 * t136 / 0.2e1) * t154;
	t111 = atan2(t116, t115);
	t108 = sin(t111);
	t134 = sin(qJ(1));
	t152 = t134 * t108;
	t109 = cos(t111);
	t151 = t134 * t109;
	t137 = cos(qJ(1));
	t150 = t137 * t108;
	t149 = t137 * t109;
	t146 = t125 * t161;
	t145 = 0.1e1 / t127 ^ 2 * t147;
	t144 = t136 * t145;
	t143 = t138 * t145;
	t114 = 0.1e1 / t115 ^ 2;
	t113 = -t128 * t156 + t141 * t131 ^ 2 * t164 + (-t130 * t124 - t153) * pkin(1);
	t112 = (t130 * t142 + (0.2e1 * t128 * pkin(5) - t124 - t156) * t131) * pkin(1);
	t107 = ((t113 * t146 + t118 * t143 - t112 * t155 / 0.2e1 - t117 * t144) / t115 - (t112 * t146 + t117 * t143 + t113 * t155 / 0.2e1 + t118 * t144) * t116 * t114) / (t116 ^ 2 * t114 + 0.1e1) * t139;
	t1 = [-t151, -t107 * t150, 0, 0; t149, -t107 * t152, 0, 0; 0, t107 * t109, 0, 0; t152, -t107 * t149, 0, 0; -t150, -t107 * t151, 0, 0; 0, -t107 * t108, 0, 0; t137, 0, 0, 0; t134, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:18:16
	% EndTime: 2020-04-19 18:18:18
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (12662->44), mult. (19532->80), div. (856->7), fcn. (12626->13), ass. (0->55)
	t153 = pkin(5) ^ 2;
	t148 = sin(qJ(2));
	t150 = sin(pkin(16));
	t151 = cos(qJ(2));
	t175 = cos(pkin(16));
	t143 = t148 * t150 - t151 * t175;
	t174 = pkin(5) * t143;
	t179 = -2 * pkin(1);
	t166 = (pkin(1) ^ 2) + t174 * t179;
	t140 = t153 + t166;
	t137 = pkin(2) ^ 2 - pkin(6) ^ 2 + t140;
	t141 = pkin(1) - t174;
	t144 = t148 * t175 + t151 * t150;
	t178 = -pkin(6) - pkin(2);
	t135 = (pkin(5) - t178) * (pkin(5) + t178) + t166;
	t177 = -pkin(6) + pkin(2);
	t136 = (pkin(5) - t177) * (pkin(5) + t177) + t166;
	t156 = sqrt(-t136 * t135);
	t173 = pkin(5) * t144;
	t165 = pkin(1) * t173;
	t170 = 0.1e1 / t156 * (t135 + t136) * t165;
	t125 = (t143 * t156 + (t141 * t179 - t137 - t170) * t144) * pkin(5);
	t167 = t144 * t156;
	t126 = t141 * t170 + t153 * t144 ^ 2 * t179 + (-t143 * t137 - t167) * pkin(5);
	t130 = -pkin(5) * t167 + t137 * t141;
	t131 = t137 * t173 + t141 * t156;
	t147 = cos(pkin(19));
	t138 = 0.1e1 / t140;
	t154 = 0.1e1 / pkin(2);
	t168 = t138 * t154;
	t146 = sin(pkin(19));
	t176 = t146 / 0.2e1;
	t128 = (-t130 * t147 / 0.2e1 + t131 * t176) * t168;
	t127 = 0.1e1 / t128 ^ 2;
	t129 = (t131 * t147 / 0.2e1 + t130 * t176) * t168;
	t163 = 0.1e1 / t140 ^ 2 * t165;
	t161 = t147 * t163;
	t162 = t146 * t163;
	t164 = t138 * t176;
	t169 = t138 * t147;
	t119 = ((t126 * t169 / 0.2e1 + t131 * t161 + t125 * t164 + t130 * t162) / t128 - (-t125 * t169 / 0.2e1 - t130 * t161 + t126 * t164 + t131 * t162) * t129 * t127) / (t127 * t129 ^ 2 + 0.1e1) * t154;
	t180 = -t119 - 0.1e1;
	t124 = atan2(t129, t128);
	t121 = sin(t124);
	t172 = t121 * t148;
	t122 = cos(t124);
	t171 = t122 * t151;
	t160 = t121 * t151 + t122 * t148;
	t159 = -t171 + t172;
	t158 = t159 * t119 - t171;
	t157 = t180 * t160;
	t152 = cos(qJ(1));
	t149 = sin(qJ(1));
	t120 = t149 * t172;
	t1 = [-t149 * t171 + t120, t157 * t152, 0, 0; -t159 * t152, t157 * t149, 0, 0; 0, t180 * t159, 0, 0; t160 * t149, (t158 + t172) * t152, 0, 0; -t160 * t152, t158 * t149 + t120, 0, 0; 0, t157, 0, 0; t152, 0, 0, 0; t149, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-19 18:19:16
	% EndTime: 2020-04-19 18:20:04
	% DurationCPUTime: 29.71s
	% Computational Cost: add. (682060->136), mult. (1030012->244), div. (44976->15), fcn. (653498->24), ass. (0->143)
	t412 = sin(qJ(2));
	t413 = sin(pkin(16));
	t415 = cos(qJ(2));
	t416 = cos(pkin(16));
	t396 = t412 * t413 - t415 * t416;
	t395 = pkin(5) * t396;
	t392 = (-0.2e1 * t395 + pkin(1)) * pkin(1);
	t419 = -pkin(6) + pkin(2);
	t420 = -pkin(6) - pkin(2);
	t328 = sqrt(-((pkin(5) - t419) * (pkin(5) + t419) + t392) * ((pkin(5) - t420) * (pkin(5) + t420) + t392));
	t427 = pkin(5) ^ 2;
	t319 = t392 + t427;
	t329 = pkin(2) ^ 2;
	t431 = -pkin(6) ^ 2 + t329;
	t390 = t319 + t431;
	t393 = -t395 + pkin(1);
	t320 = t412 * t416 + t415 * t413;
	t409 = pkin(5) * t320;
	t429 = 0.1e1 / pkin(2);
	t387 = t429 * (-t328 * t409 + t390 * t393);
	t423 = 0.1e1 / t319;
	t386 = t423 * t387;
	t384 = -t386 / 0.2e1;
	t389 = pkin(5) * t390;
	t388 = t429 * (t320 * t389 + t328 * t393);
	t385 = t423 * t388 / 0.2e1;
	t411 = sin(qJ(3));
	t414 = cos(qJ(3));
	t379 = t414 * t384 + t411 * t385;
	t383 = t386 / 0.2e1;
	t380 = t411 * t383 + t414 * t385;
	t406 = pkin(18) + pkin(19);
	t401 = sin(t406);
	t402 = cos(t406);
	t372 = t402 * t379 + t401 * t380;
	t424 = -2 * pkin(3);
	t364 = (-t372 * t424 + pkin(4)) * pkin(4);
	t418 = (-pkin(8) - pkin(10));
	t355 = ((pkin(3) - t418) * (pkin(3) + t418)) + t364;
	t417 = (-pkin(8) + pkin(10));
	t356 = ((pkin(3) - t417) * (pkin(3) + t417)) + t364;
	t327 = sqrt(-t356 * t355);
	t371 = t401 * t379 - t402 * t380;
	t439 = pkin(4) * t371;
	t440 = t327 * t439;
	t369 = t372 * pkin(4);
	t330 = pkin(8) ^ 2;
	t428 = pkin(4) ^ 2;
	t363 = t428 + (0.2e1 * t369 + pkin(3)) * pkin(3);
	t359 = -pkin(10) ^ 2 + t330 + t363;
	t357 = pkin(4) * t359;
	t438 = t357 * t371;
	t365 = -t369 - pkin(3);
	t437 = 0.2e1 * t365;
	t361 = 0.1e1 / t363;
	t436 = t361 / 0.2e1;
	t435 = -t365 / 0.2e1;
	t434 = -t439 / 0.2e1;
	t421 = pkin(3) * pkin(4);
	t433 = 0.1e1 / t363 ^ 2 * t421;
	t405 = pkin(1) * t409;
	t407 = 0.4e1 / t328 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t392 - t329) * t405;
	t311 = (t396 * t328 + (-t407 / 0.2e1 - t427 - t431) * t320) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t395) * t405;
	t312 = t393 * t407 / 0.2e1 - t396 * t389 + (-t328 - 0.2e1 * t405) * t409;
	t322 = sin(pkin(19));
	t324 = cos(pkin(19));
	t314 = t322 * t385 + t324 * t384;
	t313 = 0.1e1 / t314 ^ 2;
	t315 = t322 * t383 + t324 * t385;
	t318 = 0.1e1 / t319 ^ 2;
	t381 = t387 * t405;
	t382 = t388 * t405;
	t410 = t429 * t423;
	t404 = t410 / 0.2e1;
	t398 = t322 * t404;
	t399 = t312 * t404;
	t400 = -t311 * t410 / 0.2e1;
	t298 = ((t324 * t399 + t311 * t398 + (t322 * t381 + t324 * t382) * t318) / t314 - (t324 * t400 + t312 * t398 + (t322 * t382 - t324 * t381) * t318) * t315 * t313) / (t313 * t315 ^ 2 + 0.1e1);
	t310 = atan2(t315, t314);
	t308 = cos(t310);
	t306 = t415 * t308;
	t307 = sin(t310);
	t403 = t412 * t307;
	t430 = -t403 + t306;
	t432 = t430 * t298;
	t425 = 0.1e1 / pkin(8);
	t422 = 0.1e1 / t327;
	t408 = cos(pkin(18));
	t397 = t411 * t404;
	t304 = t415 * t307 + t412 * t308;
	t391 = -t306 - t432;
	t294 = (-t298 - 0.1e1) * t304;
	t376 = t414 * t400 + t312 * t397 + (-t414 * t381 + t411 * t382) * t318;
	t375 = t414 * t399 + t311 * t397 + (t411 * t381 + t414 * t382) * t318;
	t367 = pkin(3) * t439;
	t366 = t428 * t371 * t424;
	t360 = t425 * t361;
	t358 = t360 / 0.2e1;
	t354 = -t401 * t375 - t402 * t376;
	t353 = pkin(4) * t354;
	t352 = pkin(4) * (-t402 * t375 + t401 * t376);
	t351 = pkin(3) * t353;
	t350 = -t327 * t365 + t438;
	t349 = -t359 * t365 - t440;
	t348 = 0.1e1 / t349 ^ 2;
	t347 = t425 * t350;
	t346 = t425 * t349;
	t345 = 0.4e1 * t422 * (((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t364) * t367 - t371 * t330 * t421);
	t344 = t347 * t433;
	t343 = t346 * t433;
	t342 = 0.2e1 * t422 * (t355 + t356) * t351;
	t341 = 0.1e1 / (t348 * t350 ^ 2 + 0.1e1) / t360;
	t340 = atan2(t347 * t436, t346 * t436);
	t339 = sin(t340);
	t338 = 0.2e1 / t349 * t341;
	t337 = -0.2e1 * t348 * t350 * t341;
	t336 = ((t345 * t435 + t372 * t357 + t366 * t371 - t440) * t358 + t371 * t344) * t338 + ((-t327 * t369 + t345 * t434 + t367 * t437 - t438) * t358 + t371 * t343) * t337;
	t297 = cos(t340);
	t335 = t297 * t336;
	t334 = ((-t327 * t353 + t342 * t435 + t359 * t352 + t354 * t366) * t358 + t354 * t344) * t338 + ((-t327 * t352 + t342 * t434 + t351 * t437 - t359 * t353) * t358 + t354 * t343) * t337;
	t333 = t297 * t334;
	t332 = t336 * t339;
	t331 = t334 * t339;
	t326 = cos(qJ(1));
	t325 = sin(qJ(1));
	t323 = sin(pkin(18));
	t305 = t325 * t403;
	t302 = t430 * t326;
	t301 = t304 * t326;
	t300 = t325 * t306 - t305;
	t299 = t304 * t325;
	t295 = t430 + t432;
	t293 = t294 * t326;
	t292 = (t403 + t391) * t326;
	t291 = t294 * t325;
	t290 = t325 * t391 + t305;
	t289 = -t408 * t297 - t323 * t339;
	t288 = t297 * t323 - t408 * t339;
	t283 = -t323 * t335 + t408 * t332;
	t282 = -t323 * t332 - t408 * t335;
	t281 = -t323 * t333 + t408 * t331;
	t280 = -t323 * t331 - t408 * t333;
	t1 = [t288 * t299 - t289 * t300, -t280 * t301 + t281 * t302 + t288 * t292 + t289 * t293, -t282 * t301 + t283 * t302, 0; -t288 * t301 + t289 * t302, -t280 * t299 + t281 * t300 + t288 * t290 + t289 * t291, -t282 * t299 + t283 * t300, 0; 0, t280 * t430 + t281 * t304 + t288 * t294 + t289 * t295, t282 * t430 + t283 * t304, 0; t288 * t300 + t289 * t299, -t280 * t302 - t281 * t301 - t288 * t293 + t289 * t292, -t282 * t302 - t283 * t301, 0; -t288 * t302 - t289 * t301, -t280 * t300 - t281 * t299 - t288 * t291 + t289 * t290, -t282 * t300 - t283 * t299, 0; 0, -t280 * t304 + t281 * t430 - t288 * t295 + t289 * t294, -t282 * t304 + t283 * t430, 0; t326, 0, 0, 0; t325, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end