% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m1TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_jacobiR_rot_sym_varpar: pkin has to be [23x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
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
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
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
	% StartTime: 2020-04-13 14:17:58
	% EndTime: 2020-04-13 14:17:58
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
	% StartTime: 2020-04-13 14:18:17
	% EndTime: 2020-04-13 14:18:25
	% DurationCPUTime: 8.00s
	% Computational Cost: add. (205352->115), mult. (311192->198), div. (12992->8), fcn. (197874->16), ass. (0->113)
	t376 = pkin(7) ^ 2;
	t370 = sin(qJ(2));
	t374 = cos(pkin(19));
	t439 = sin(pkin(19));
	t440 = cos(qJ(2));
	t356 = t370 * t374 - t440 * t439;
	t435 = pkin(7) * t356;
	t448 = -2 * pkin(1);
	t416 = (pkin(1) ^ 2) + t435 * t448;
	t348 = t376 + t416;
	t346 = 0.1e1 / t348;
	t379 = 0.1e1 / pkin(3);
	t428 = t346 * t379;
	t345 = pkin(3) ^ 2 - pkin(8) ^ 2 + t348;
	t353 = pkin(1) - t435;
	t445 = -pkin(8) - pkin(3);
	t343 = (pkin(7) - t445) * (pkin(7) + t445) + t416;
	t444 = -pkin(8) + pkin(3);
	t344 = (pkin(7) - t444) * (pkin(7) + t444) + t416;
	t381 = sqrt(-t344 * t343);
	t358 = t370 * t439 + t440 * t374;
	t434 = pkin(7) * t358;
	t339 = t345 * t434 + t353 * t381;
	t372 = cos(qJ(3));
	t430 = t339 * t372;
	t427 = t358 * t381;
	t338 = -pkin(7) * t427 + t345 * t353;
	t369 = sin(qJ(3));
	t433 = t338 * t369;
	t334 = (t433 / 0.2e1 + t430 / 0.2e1) * t428;
	t431 = t339 * t369;
	t432 = t338 * t372;
	t335 = (-t432 / 0.2e1 + t431 / 0.2e1) * t428;
	t366 = pkin(23) + pkin(22);
	t364 = sin(t366);
	t365 = cos(t366);
	t401 = t334 * t365 - t335 * t364;
	t391 = pkin(4) * t401;
	t326 = pkin(5) + t391;
	t459 = t326 / 0.2e1;
	t367 = sin(pkin(21));
	t458 = t367 / 0.2e1;
	t368 = cos(pkin(21));
	t457 = -t368 / 0.2e1;
	t404 = t334 * t364 + t365 * t335;
	t437 = pkin(4) * t404;
	t456 = -t437 / 0.2e1;
	t378 = pkin(4) ^ 2;
	t446 = pkin(4) * pkin(5);
	t424 = -0.2e1 * t446;
	t417 = pkin(5) ^ 2 - t401 * t424;
	t325 = t378 + t417;
	t449 = pkin(9) ^ 2;
	t323 = pkin(11) ^ 2 + t325 - t449;
	t442 = -pkin(9) + pkin(11);
	t443 = -pkin(9) - pkin(11);
	t382 = sqrt(-((pkin(4) - t442) * (pkin(4) + t442) + t417) * ((pkin(4) - t443) * (pkin(4) + t443) + t417));
	t436 = pkin(4) * t382;
	t319 = t323 * t326 - t404 * t436;
	t320 = t323 * t437 + t326 * t382;
	t375 = 0.1e1 / pkin(11);
	t447 = 0.1e1 / t325;
	t414 = t447 / 0.2e1;
	t409 = t367 * t414;
	t315 = (t320 * t368 * t414 + t319 * t409) * t375;
	t316 = (t319 * t447 * t457 + t320 * t409) * t375;
	t357 = -t370 * t369 + t440 * t372;
	t371 = sin(qJ(1));
	t349 = t357 * t371;
	t392 = t440 * t369 + t370 * t372;
	t350 = t392 * t371;
	t423 = -t349 * t315 - t350 * t316;
	t422 = t350 * t315 - t349 * t316;
	t373 = cos(qJ(1));
	t351 = t357 * t373;
	t352 = t392 * t373;
	t455 = -t352 * t315 + t351 * t316;
	t420 = -t351 * t315 - t352 * t316;
	t405 = 0.1e1 / t325 ^ 2 * t375 * t446;
	t402 = t320 * t405;
	t403 = t319 * t405;
	t454 = t367 * t402 - t368 * t403;
	t453 = t367 * t403 + t368 * t402;
	t452 = -0.2e1 * t378 * t404 * pkin(5) - t436;
	t451 = -pkin(4) * t323 + t326 * t424;
	t450 = 0.4e1 / t382 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t417 - t449) * t446;
	t441 = -t372 / 0.2e1;
	t413 = pkin(1) * t434;
	t429 = 0.2e1 / t381 * (t343 + t344) * t413;
	t425 = t375 * t447;
	t419 = -t315 * t392 + t357 * t316;
	t418 = -t357 * t315 - t316 * t392;
	t408 = t375 * t414;
	t406 = 0.1e1 / t348 ^ 2 * t413;
	t400 = t431 - t432;
	t399 = -t430 - t433;
	t331 = (t356 * t381 + (-t429 / 0.2e1 - t345 + t353 * t448) * t358) * pkin(7);
	t332 = t353 * t429 / 0.2e1 + t376 * t358 ^ 2 * t448 + (-t345 * t356 - t427) * pkin(7);
	t394 = -t369 * t331 / 0.2e1 + t332 * t441;
	t393 = t331 * t441 + t369 * t332 / 0.2e1;
	t321 = ((t364 * t393 + t365 * t394) * t346 + (t364 * t400 + t365 * t399) * t406) * t379;
	t389 = t321 * t450;
	t388 = t404 * t450;
	t387 = (-t382 * t391 + t388 * t456 + t451 * t404) * t425;
	t386 = (t323 * t391 + t388 * t459 + t452 * t404) * t408;
	t385 = pkin(4) * t379 * ((-t364 * t394 + t365 * t393) * t346 + (-t364 * t399 + t365 * t400) * t406);
	t384 = (t451 * t321 - t382 * t385 + t389 * t456) * t425;
	t383 = (t452 * t321 + t323 * t385 + t389 * t459) * t408;
	t298 = t367 * t386 + t387 * t457 + t454 * t404;
	t297 = t368 * t386 + t387 * t458 + t453 * t404;
	t296 = t454 * t321 + t367 * t383 + t384 * t457;
	t295 = t453 * t321 + t368 * t383 + t384 * t458;
	t1 = [t422, -t295 * t352 + t296 * t351 + t420, -t297 * t352 + t298 * t351 + t420, 0; t455, -t295 * t350 + t296 * t349 + t423, -t297 * t350 + t298 * t349 + t423, 0; 0, t295 * t357 + t296 * t392 + t419, t297 * t357 + t298 * t392 + t419, 0; -t423, -t295 * t351 - t296 * t352 - t455, -t297 * t351 - t298 * t352 - t455, 0; t420, -t295 * t349 - t296 * t350 + t422, -t297 * t349 - t298 * t350 + t422, 0; 0, -t295 * t392 + t296 * t357 + t418, -t297 * t392 + t298 * t357 + t418, 0; t373, 0, 0, 0; t371, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:19:16
	% EndTime: 2020-04-13 14:19:35
	% DurationCPUTime: 18.95s
	% Computational Cost: add. (324850->131), mult. (492070->219), div. (20688->8), fcn. (312822->18), ass. (0->128)
	t629 = sin(qJ(2));
	t630 = sin(pkin(19));
	t632 = cos(qJ(2));
	t633 = cos(pkin(19));
	t540 = t629 * t630 + t632 * t633;
	t601 = t629 * t633 - t632 * t630;
	t600 = pkin(7) * t601;
	t598 = (-0.2e1 * t600 + pkin(1)) * pkin(1);
	t636 = -pkin(8) + pkin(3);
	t637 = -pkin(8) - pkin(3);
	t552 = sqrt(-((pkin(7) - t636) * (pkin(7) + t636) + t598) * ((pkin(7) - t637) * (pkin(7) + t637) + t598));
	t645 = pkin(7) ^ 2;
	t533 = t598 + t645;
	t553 = pkin(3) ^ 2;
	t653 = -pkin(8) ^ 2 + t553;
	t597 = t533 + t653;
	t596 = pkin(7) * t597;
	t599 = -t600 + pkin(1);
	t646 = 0.1e1 / pkin(3);
	t595 = t646 * (t540 * t596 + t599 * t552);
	t639 = 0.1e1 / t533;
	t592 = t639 * t595 / 0.2e1;
	t625 = pkin(7) * t540;
	t594 = t646 * (-t552 * t625 + t599 * t597);
	t593 = t639 * t594;
	t628 = sin(qJ(3));
	t631 = cos(qJ(3));
	t588 = t628 * t593 / 0.2e1 + t631 * t592;
	t655 = -t631 / 0.2e1;
	t589 = t628 * t592 + t593 * t655;
	t609 = pkin(23) + pkin(22);
	t606 = sin(t609);
	t607 = cos(t609);
	t581 = t607 * t588 - t606 * t589;
	t579 = pkin(4) * t581;
	t576 = (0.2e1 * t579 + pkin(5)) * pkin(5);
	t634 = -pkin(9) + pkin(11);
	t635 = -pkin(9) - pkin(11);
	t551 = sqrt(-((pkin(4) - t634) * (pkin(4) + t634) + t576) * ((pkin(4) - t635) * (pkin(4) + t635) + t576));
	t641 = -0.2e1 * pkin(5);
	t526 = pkin(5) ^ 2 + (-t581 * t641 + pkin(4)) * pkin(4);
	t643 = pkin(9) ^ 2;
	t574 = pkin(11) ^ 2 + t526 - t643;
	t573 = pkin(4) * t574;
	t577 = t579 + pkin(5);
	t642 = 0.1e1 / pkin(11);
	t647 = t606 * t588 + t607 * t589;
	t568 = t642 * (t577 * t551 + t573 * t647);
	t640 = 0.1e1 / t526;
	t565 = t640 * t568 / 0.2e1;
	t580 = pkin(4) * t647;
	t567 = t642 * (-t551 * t580 + t577 * t574);
	t566 = t640 * t567;
	t623 = cos(pkin(21));
	t622 = sin(pkin(21));
	t657 = t622 / 0.2e1;
	t520 = t623 * t565 + t566 * t657;
	t656 = -t623 / 0.2e1;
	t521 = t622 * t565 + t566 * t656;
	t539 = -t629 * t628 + t632 * t631;
	t548 = sin(qJ(1));
	t534 = t539 * t548;
	t538 = -t632 * t628 - t629 * t631;
	t535 = t538 * t548;
	t508 = -t520 * t535 - t521 * t534;
	t547 = sin(qJ(4));
	t549 = cos(qJ(4));
	t550 = cos(qJ(1));
	t661 = t508 * t549 + t550 * t547;
	t660 = t508 * t547 - t549 * t550;
	t659 = t577 / 0.2e1;
	t658 = -t580 / 0.2e1;
	t638 = pkin(4) * pkin(5);
	t654 = 0.1e1 / t526 ^ 2 * t638;
	t615 = -t534 * t520 + t535 * t521;
	t536 = t539 * t550;
	t537 = t538 * t550;
	t614 = -t536 * t520 + t537 * t521;
	t563 = t567 * t654;
	t564 = t568 * t654;
	t652 = t622 * t563 + t623 * t564;
	t651 = -t623 * t563 + t622 * t564;
	t650 = -0.2e1 * t577 * t638 - t573;
	t626 = pkin(4) * t551;
	t649 = pkin(4) ^ 2 * t641 * t647 - t626;
	t648 = 0.4e1 / t551 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t576 - t643) * t638;
	t627 = t646 * t639;
	t624 = t642 * t640;
	t611 = pkin(1) * t625;
	t621 = 0.4e1 / t552 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t598 - t553) * t611;
	t528 = (t601 * t552 + (-t621 / 0.2e1 - t645 - t653) * t540) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t600) * t611;
	t529 = t599 * t621 / 0.2e1 - t601 * t596 + (-t552 - 0.2e1 * t611) * t625;
	t532 = 0.1e1 / t533 ^ 2;
	t590 = t594 * t611;
	t591 = t595 * t611;
	t602 = t627 * t655;
	t605 = t628 * t627;
	t584 = -t528 * t605 / 0.2e1 + t529 * t602 + (-t628 * t590 - t631 * t591) * t532;
	t585 = t528 * t602 + t529 * t605 / 0.2e1 + (-t631 * t590 + t628 * t591) * t532;
	t524 = t607 * t584 + t606 * t585;
	t569 = pkin(4) * (-t606 * t584 + t607 * t585);
	t571 = t524 * t648;
	t608 = t624 / 0.2e1;
	t555 = (t649 * t524 + t574 * t569 + t571 * t659) * t608;
	t556 = (t650 * t524 - t551 * t569 + t571 * t658) * t624;
	t500 = t652 * t524 + t623 * t555 + t556 * t657;
	t619 = t500 + t521;
	t501 = t651 * t524 + t622 * t555 + t556 * t656;
	t618 = t501 - t520;
	t570 = t647 * t648;
	t557 = (t570 * t659 + t581 * t573 + t647 * t649) * t608;
	t558 = (t570 * t658 - t581 * t626 + t647 * t650) * t624;
	t502 = t623 * t557 + t558 * t657 + t647 * t652;
	t617 = t502 + t521;
	t503 = t622 * t557 + t558 * t656 + t647 * t651;
	t616 = t503 - t520;
	t613 = t538 * t520 + t539 * t521;
	t511 = t520 * t539 - t521 * t538;
	t510 = t520 * t537 + t521 * t536;
	t507 = t510 * t549 + t547 * t548;
	t506 = -t510 * t547 + t548 * t549;
	t499 = t502 * t539 - t503 * t538 + t613;
	t498 = t502 * t537 + t503 * t536 + t614;
	t497 = t502 * t535 + t503 * t534 + t615;
	t496 = t500 * t539 - t501 * t538 + t613;
	t495 = t500 * t537 + t501 * t536 + t614;
	t494 = t500 * t535 + t501 * t534 + t615;
	t1 = [t661, t495 * t549, t498 * t549, t506; t507, t494 * t549, t497 * t549, t660; 0, t496 * t549, t499 * t549, -t511 * t547; -t660, -t495 * t547, -t498 * t547, -t507; t506, -t494 * t547, -t497 * t547, t661; 0, -t496 * t547, -t499 * t547, -t511 * t549; t615, t619 * t536 - t618 * t537, t617 * t536 - t616 * t537, 0; -t614, t619 * t534 - t618 * t535, t617 * t534 - t616 * t535, 0; 0, -t619 * t538 - t618 * t539, -t617 * t538 - t616 * t539, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:18:00
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (1338->36), mult. (2084->64), div. (64->4), fcn. (1354->10), ass. (0->41)
	t133 = pkin(1) ^ 2;
	t125 = sin(qJ(2));
	t129 = cos(pkin(19));
	t147 = sin(pkin(19));
	t148 = cos(qJ(2));
	t122 = t125 * t129 - t148 * t147;
	t146 = pkin(1) * t122;
	t153 = -2 * pkin(7);
	t141 = t146 * t153 + t133;
	t119 = (pkin(7) ^ 2) + t141;
	t116 = -pkin(3) ^ 2 + pkin(8) ^ 2 + t119;
	t120 = -pkin(7) + t146;
	t123 = t125 * t147 + t148 * t129;
	t152 = -pkin(8) - pkin(3);
	t114 = (pkin(7) - t152) * (pkin(7) + t152) + t141;
	t151 = -pkin(8) + pkin(3);
	t115 = (pkin(7) - t151) * (pkin(7) + t151) + t141;
	t134 = sqrt(-t115 * t114);
	t142 = t123 * t134;
	t109 = -pkin(1) * t142 - t120 * t116;
	t131 = 0.1e1 / pkin(8);
	t145 = t123 * pkin(1);
	t140 = pkin(7) * t145;
	t139 = 0.1e1 / t119 ^ 2 * t131 * t140;
	t143 = 0.1e1 / t119 * t131;
	t144 = 0.1e1 / t134 * (t114 + t115) * t140;
	t154 = pkin(1) * (t122 * t134 + (0.2e1 * t120 * pkin(7) - t116 - t144) * t123) * t143 / 0.2e1 + t109 * t139;
	t127 = sin(pkin(18));
	t150 = -t127 / 0.2e1;
	t130 = cos(pkin(18));
	t149 = t130 / 0.2e1;
	t110 = t116 * t145 - t120 * t134;
	t137 = t110 * t139;
	t135 = (-t120 * t144 + t133 * t123 ^ 2 * t153 + (-t122 * t116 - t142) * pkin(1)) * t143;
	t128 = cos(qJ(1));
	t126 = sin(qJ(1));
	t108 = (-t130 * t110 / 0.2e1 + t109 * t150) * t143;
	t107 = (t109 * t149 + t110 * t150) * t143;
	t104 = t154 * t127 + t130 * t137 + t135 * t149;
	t103 = -t127 * t137 + t154 * t130 + t135 * t150;
	t1 = [-t126 * t107, t128 * t103, 0, 0; t128 * t107, t126 * t103, 0, 0; 0, t104, 0, 0; -t126 * t108, -t128 * t104, 0, 0; t128 * t108, -t126 * t104, 0, 0; 0, t103, 0, 0; t128, 0, 0, 0; t126, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:00
	% EndTime: 2020-04-13 14:18:01
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (3354->45), mult. (5216->74), div. (176->4), fcn. (3394->10), ass. (0->51)
	t130 = sin(qJ(2));
	t132 = sin(pkin(19));
	t133 = cos(qJ(2));
	t135 = cos(pkin(19));
	t154 = t130 * t135 - t133 * t132;
	t153 = pkin(7) * t154;
	t150 = (-0.2e1 * t153 + pkin(1)) * pkin(1);
	t171 = pkin(7) ^ 2;
	t128 = t150 + t171;
	t170 = 0.1e1 / t128;
	t179 = t170 / 0.2e1;
	t172 = 0.1e1 / pkin(3);
	t178 = t172 * t179;
	t168 = -pkin(8) + pkin(3);
	t169 = -pkin(8) - pkin(3);
	t136 = sqrt(-((pkin(7) - t168) * (pkin(7) + t168) + t150) * ((pkin(7) - t169) * (pkin(7) + t169) + t150));
	t129 = t130 * t132 + t133 * t135;
	t137 = pkin(3) ^ 2;
	t149 = -pkin(8) ^ 2 + t128 + t137;
	t148 = pkin(7) * t149;
	t147 = t129 * t148;
	t152 = -t153 + pkin(1);
	t145 = t172 * (t136 * t152 + t147);
	t142 = t145 * t179;
	t165 = pkin(7) * t136;
	t157 = t129 * t165;
	t144 = t172 * (t149 * t152 - t157);
	t143 = t170 * t144;
	t163 = sin(pkin(23));
	t164 = cos(pkin(23));
	t124 = t164 * t142 + t163 * t143 / 0.2e1;
	t166 = pkin(7) * t129;
	t158 = pkin(1) * t166;
	t146 = 0.2e1 / t136 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t150 - t137) * t158;
	t176 = 0.1e1 / t128 ^ 2 * t158;
	t173 = (-t146 * t166 - 0.2e1 * t152 * t158 + t154 * t165 - t147) * t178 + t144 * t176;
	t174 = (-0.2e1 * pkin(1) * t129 ^ 2 * t171 + t146 * t152 - t148 * t154 - t157) * t178 + t145 * t176;
	t160 = t163 * t174 - t164 * t173 - t124;
	t156 = t160 * t133;
	t119 = t163 * t173 + t164 * t174;
	t123 = -t164 * t143 / 0.2e1 + t163 * t142;
	t159 = t119 + t123;
	t177 = -t159 * t130 + t156;
	t162 = t123 * t130;
	t161 = t124 * t133;
	t155 = t123 * t133 - t124 * t130;
	t151 = -t130 * t160 - t133 * t159;
	t134 = cos(qJ(1));
	t131 = sin(qJ(1));
	t122 = t131 * t162;
	t1 = [t131 * t161 + t122, t151 * t134, 0, 0; (-t161 - t162) * t134, t151 * t131, 0, 0; 0, t177, 0, 0; t155 * t131, -t177 * t134, 0, 0; -t155 * t134, t122 + (t119 * t130 - t156) * t131, 0, 0; 0, t151, 0, 0; t134, 0, 0, 0; t131, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:17:59
	% EndTime: 2020-04-13 14:18:00
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (1666->39), mult. (2612->67), div. (88->4), fcn. (1702->10), ass. (0->48)
	t127 = pkin(6) ^ 2;
	t119 = sin(pkin(20));
	t120 = cos(pkin(20));
	t121 = sin(qJ(3));
	t124 = cos(qJ(3));
	t117 = t119 * t124 + t120 * t121;
	t149 = pkin(6) * t117;
	t153 = 2 * pkin(1);
	t140 = t149 * t153 + t127;
	t114 = (pkin(1) ^ 2) + t140;
	t112 = 0.1e1 / t114;
	t125 = cos(qJ(2));
	t128 = 0.1e1 / pkin(2);
	t118 = t119 * t121 - t120 * t124;
	t148 = pkin(6) * t118;
	t139 = pkin(1) * t148;
	t138 = 0.1e1 / t114 ^ 2 * t139;
	t111 = pkin(2) ^ 2 - pkin(13) ^ 2 + t114;
	t115 = -pkin(1) - t149;
	t152 = -pkin(2) - pkin(13);
	t109 = (pkin(1) - t152) * (pkin(1) + t152) + t140;
	t151 = -pkin(2) + pkin(13);
	t110 = (pkin(1) - t151) * (pkin(1) + t151) + t140;
	t130 = sqrt(-t110 * t109);
	t145 = t118 * t130;
	t104 = -pkin(6) * t145 - t111 * t115;
	t142 = t125 * t104;
	t105 = t111 * t148 - t115 * t130;
	t122 = sin(qJ(2));
	t143 = t122 * t105;
	t150 = -t122 / 0.2e1;
	t147 = 0.1e1 / t130 * (t109 + t110) * t139;
	t98 = (-t117 * t130 + (t115 * t153 - t111 - t147) * t118) * pkin(6);
	t99 = -t115 * t147 - 0.2e1 * t127 * t118 ^ 2 * pkin(1) + (t111 * t117 - t145) * pkin(6);
	t132 = ((t125 * t98 / 0.2e1 + t99 * t150) * t112 + (t142 - t143) * t138) * t128;
	t141 = t125 * t105;
	t144 = t122 * t104;
	t146 = t112 * t128;
	t133 = (-t141 / 0.2e1 - t144 / 0.2e1) * t146;
	t134 = (t143 / 0.2e1 - t142 / 0.2e1) * t146;
	t131 = ((-t125 * t99 / 0.2e1 + t98 * t150) * t112 + (-t141 - t144) * t138) * t128;
	t126 = cos(qJ(1));
	t123 = sin(qJ(1));
	t103 = t126 * t134;
	t102 = t126 * t133;
	t101 = t123 * t134;
	t100 = t123 * t133;
	t1 = [-t100, t103, t126 * t131, 0; t102, t101, t123 * t131, 0; 0, t133, t132, 0; -t101, -t102, -t126 * t132, 0; t103, -t100, -t123 * t132, 0; 0, t134, t131, 0; t126, 0, 0, 0; t123, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:00
	% EndTime: 2020-04-13 14:18:01
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (4538->62), mult. (7012->122), div. (312->5), fcn. (4470->10), ass. (0->69)
	t185 = sin(pkin(20));
	t186 = cos(pkin(20));
	t187 = sin(qJ(3));
	t190 = cos(qJ(3));
	t183 = t185 * t190 + t186 * t187;
	t226 = pkin(6) * t183;
	t213 = pkin(1) * t226;
	t182 = 0.2e1 * t213;
	t196 = pkin(2) ^ 2;
	t195 = pkin(6) ^ 2;
	t215 = pkin(1) ^ 2 + t195;
	t212 = -pkin(13) ^ 2 + t215;
	t177 = t182 + t196 + t212;
	t181 = -pkin(1) - t226;
	t184 = t185 * t187 - t186 * t190;
	t216 = t182 + t195;
	t233 = -pkin(2) - pkin(13);
	t174 = (pkin(1) - t233) * (pkin(1) + t233) + t216;
	t232 = -pkin(2) + pkin(13);
	t175 = (pkin(1) - t232) * (pkin(1) + t232) + t216;
	t199 = sqrt(-t175 * t174);
	t225 = pkin(6) * t184;
	t214 = pkin(1) * t225;
	t224 = 0.1e1 / t199 * (t174 + t175) * t214;
	t153 = (-t183 * t199 + (0.2e1 * t181 * pkin(1) - t177 - t224) * t184) * pkin(6);
	t222 = t184 * t199;
	t154 = -t181 * t224 - 0.2e1 * t195 * t184 ^ 2 * pkin(1) + (t177 * t183 - t222) * pkin(6);
	t180 = t182 + t215;
	t178 = 0.1e1 / t180;
	t191 = cos(qJ(2));
	t197 = 0.1e1 / pkin(2);
	t210 = 0.1e1 / t180 ^ 2 * t214;
	t169 = -pkin(6) * t222 - t177 * t181;
	t219 = t191 * t169;
	t170 = t177 * t225 - t181 * t199;
	t188 = sin(qJ(2));
	t220 = t188 * t170;
	t229 = -t188 / 0.2e1;
	t237 = ((t191 * t153 / 0.2e1 + t154 * t229) * t178 + (t219 - t220) * t210) * t197;
	t218 = t191 * t170;
	t221 = t188 * t169;
	t223 = t178 * t197;
	t160 = (-t218 / 0.2e1 - t221 / 0.2e1) * t223;
	t159 = (t220 / 0.2e1 - t219 / 0.2e1) * t223;
	t192 = cos(qJ(1));
	t158 = t192 * t159;
	t176 = t196 - t212 - 0.2e1 * t213;
	t217 = 0.1e1 / pkin(13) * t197;
	t227 = t199 / 0.2e1;
	t157 = t192 * t160;
	t231 = t157 / 0.2e1;
	t236 = (t158 * t227 + t176 * t231) * t217;
	t189 = sin(qJ(1));
	t155 = t189 * t160;
	t156 = t189 * t159;
	t228 = -t199 / 0.2e1;
	t230 = -t176 / 0.2e1;
	t235 = (-t155 * t228 + t156 * t230) * t217;
	t211 = -t224 / 0.2e1;
	t209 = t217 * t230;
	t208 = t217 * t227;
	t206 = -t155 * t209 + t156 * t208;
	t205 = t157 * t208 + t158 * t209;
	t150 = ((-t191 * t154 / 0.2e1 + t153 * t229) * t178 + (-t218 - t221) * t210) * t197;
	t148 = t192 * t237;
	t147 = t192 * t150;
	t146 = t189 * t237;
	t145 = t189 * t150;
	t1 = [t206, t205, (t147 * t230 - t148 * t228 - t157 * t214 + t158 * t211) * t217, 0; -t236, t235, (t145 * t230 - t146 * t228 - t155 * t214 + t156 * t211) * t217, 0; 0, (t159 * t228 + t160 * t230) * t217, (t150 * t228 + t159 * t214 + t160 * t211 + t230 * t237) * t217, 0; -t235, t236, (t147 * t227 - t148 * t230 - t158 * t214 + t224 * t231) * t217, 0; t205, t206, (t145 * t227 - t146 * t230 - t155 * t211 - t156 * t214) * t217, 0; 0, (t159 * t230 + t160 * t227) * t217, (t150 * t230 + t159 * t211 - t160 * t214 + t227 * t237) * t217, 0; t192, 0, 0, 0; t189, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-13 14:18:19
	% EndTime: 2020-04-13 14:18:31
	% DurationCPUTime: 12.04s
	% Computational Cost: add. (197868->128), mult. (300332->229), div. (12336->8), fcn. (191114->18), ass. (0->124)
	t403 = sin(qJ(2));
	t404 = sin(pkin(19));
	t406 = cos(qJ(2));
	t407 = cos(pkin(19));
	t325 = t403 * t404 + t406 * t407;
	t379 = t403 * t407 - t406 * t404;
	t377 = pkin(7) * t379;
	t374 = (-0.2e1 * t377 + pkin(1)) * pkin(1);
	t410 = -pkin(8) + pkin(3);
	t411 = -pkin(8) - pkin(3);
	t331 = sqrt(-((pkin(7) - t410) * (pkin(7) + t410) + t374) * ((pkin(7) - t411) * (pkin(7) + t411) + t374));
	t399 = pkin(7) * t325;
	t391 = pkin(1) * t399;
	t332 = pkin(3) ^ 2;
	t393 = 0.4e1 / t331 * ((pkin(7) + pkin(8)) * (pkin(7) - pkin(8)) + t374 - t332) * t391;
	t418 = pkin(7) ^ 2;
	t427 = -pkin(8) ^ 2 + t332;
	t315 = (t379 * t331 + (-t393 / 0.2e1 - t418 - t427) * t325) * pkin(7) + (-0.3e1 * pkin(1) + 0.4e1 * t377) * t391;
	t434 = -t315 / 0.2e1;
	t324 = t374 + t418;
	t373 = t324 + t427;
	t372 = pkin(7) * t373;
	t375 = -t377 + pkin(1);
	t316 = t375 * t393 / 0.2e1 - t379 * t372 + (-t331 - 0.2e1 * t391) * t399;
	t433 = t316 / 0.2e1;
	t327 = sin(pkin(22));
	t432 = t327 / 0.2e1;
	t420 = 0.1e1 / pkin(3);
	t370 = t420 * (-t331 * t399 + t375 * t373);
	t413 = 0.1e1 / t324;
	t369 = t413 * t370;
	t366 = t369 / 0.2e1;
	t371 = t420 * (t325 * t372 + t375 * t331);
	t368 = t413 * t371 / 0.2e1;
	t402 = sin(qJ(3));
	t405 = cos(qJ(3));
	t362 = t402 * t366 + t405 * t368;
	t367 = -t369 / 0.2e1;
	t363 = t405 * t367 + t402 * t368;
	t389 = pkin(23) + pkin(22);
	t386 = sin(t389);
	t387 = cos(t389);
	t355 = t387 * t362 - t386 * t363;
	t353 = pkin(5) * t355;
	t351 = -t353 - pkin(4);
	t431 = -t351 / 0.2e1;
	t421 = t386 * t362 + t387 * t363;
	t354 = pkin(5) * t421;
	t430 = -t354 / 0.2e1;
	t396 = cos(pkin(22));
	t429 = -t396 / 0.2e1;
	t419 = pkin(5) ^ 2;
	t305 = t419 + (0.2e1 * t353 + pkin(4)) * pkin(4);
	t412 = pkin(4) * pkin(5);
	t428 = 0.1e1 / t305 ^ 2 * t412;
	t333 = pkin(9) ^ 2;
	t348 = -pkin(11) ^ 2 + t305 + t333;
	t347 = pkin(5) * t348;
	t426 = 0.2e1 * t351 * t412 - t347;
	t415 = -0.2e1 * pkin(4);
	t350 = (-t355 * t415 + pkin(5)) * pkin(5);
	t408 = -pkin(9) + pkin(11);
	t409 = -pkin(9) - pkin(11);
	t330 = sqrt(-((pkin(4) - t408) * (pkin(4) + t408) + t350) * ((pkin(4) - t409) * (pkin(4) + t409) + t350));
	t416 = 0.1e1 / pkin(9);
	t342 = t416 * (-t351 * t330 + t347 * t421);
	t340 = t342 * t428;
	t294 = -t330 * t354 - t351 * t348;
	t397 = t294 * t416;
	t383 = t397 * t428;
	t425 = -t327 * t340 - t396 * t383;
	t424 = t327 * t383 - t396 * t340;
	t400 = pkin(5) * t330;
	t423 = t415 * t419 * t421 - t400;
	t422 = 0.4e1 / t330 * ((pkin(4) + pkin(11)) * (pkin(4) - pkin(11)) + t350 - t333) * t412;
	t414 = 0.1e1 / t305;
	t401 = t420 * t413;
	t398 = t416 * t414;
	t395 = cos(pkin(23));
	t394 = sin(pkin(23));
	t388 = -t398 / 0.2e1;
	t385 = t402 * t401;
	t384 = t395 * t401;
	t382 = -t405 * t401 / 0.2e1;
	t380 = t394 * t401 / 0.2e1;
	t317 = t395 * t367 + t394 * t368;
	t318 = t394 * t366 + t395 * t368;
	t311 = t406 * t317 - t403 * t318;
	t376 = t403 * t317 + t406 * t318;
	t323 = 0.1e1 / t324 ^ 2;
	t364 = t370 * t391;
	t365 = t371 * t391;
	t302 = t384 * t434 + t316 * t380 + (-t395 * t364 + t394 * t365) * t323;
	t303 = t384 * t433 + t315 * t380 + (t394 * t364 + t395 * t365) * t323;
	t299 = t406 * t302 - t403 * t303 - t376;
	t300 = -t403 * t302 - t406 * t303 - t311;
	t359 = t315 * t382 + t385 * t433 + (-t405 * t364 + t402 * t365) * t323;
	t358 = t385 * t434 + t316 * t382 + (-t402 * t364 - t405 * t365) * t323;
	t301 = t387 * t358 + t386 * t359;
	t345 = t301 * t422;
	t344 = t421 * t422;
	t343 = pkin(5) * (-t386 * t358 + t387 * t359);
	t341 = -t414 * t342 / 0.2e1;
	t337 = (t344 * t430 - t355 * t400 + t421 * t426) * t398;
	t336 = (t344 * t431 + t355 * t347 + t421 * t423) * t388;
	t335 = (t426 * t301 - t330 * t343 + t345 * t430) * t398;
	t334 = (t423 * t301 + t348 * t343 + t345 * t431) * t388;
	t329 = cos(qJ(1));
	t328 = sin(qJ(1));
	t309 = t311 * t329;
	t308 = t376 * t329;
	t307 = t311 * t328;
	t306 = t376 * t328;
	t298 = t299 * t329;
	t297 = t300 * t329;
	t296 = t299 * t328;
	t295 = t300 * t328;
	t289 = t396 * t294 * t388 + t327 * t341;
	t288 = t414 * t397 * t432 + t396 * t341;
	t283 = t327 * t336 + t337 * t429 + t421 * t425;
	t282 = t396 * t336 + t337 * t432 + t421 * t424;
	t281 = t425 * t301 + t327 * t334 + t335 * t429;
	t280 = t424 * t301 + t396 * t334 + t335 * t432;
	t1 = [t288 * t307 + t289 * t306, -t280 * t309 - t281 * t308 - t288 * t298 + t289 * t297, -t282 * t309 - t283 * t308, 0; -t288 * t309 - t289 * t308, -t280 * t307 - t281 * t306 - t288 * t296 + t289 * t295, -t282 * t307 - t283 * t306, 0; 0, -t280 * t376 + t281 * t311 + t288 * t300 + t289 * t299, -t282 * t376 + t283 * t311, 0; -t288 * t306 + t289 * t307, t280 * t308 - t281 * t309 - t288 * t297 - t289 * t298, t282 * t308 - t283 * t309, 0; t288 * t308 - t289 * t309, t280 * t306 - t281 * t307 - t288 * t295 - t289 * t296, t282 * t306 - t283 * t307, 0; 0, -t280 * t311 - t281 * t376 - t288 * t299 + t289 * t300, -t282 * t311 - t283 * t376, 0; t329, 0, 0, 0; t328, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end