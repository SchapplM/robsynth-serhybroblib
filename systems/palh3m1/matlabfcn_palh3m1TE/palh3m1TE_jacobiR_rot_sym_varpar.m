% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m1TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_jacobiR_rot_sym_varpar: pkin has to be [19x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
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
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:02
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
	% StartTime: 2020-04-18 09:52:02
	% EndTime: 2020-04-18 09:52:03
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
	% StartTime: 2020-04-18 09:52:29
	% EndTime: 2020-04-18 09:52:41
	% DurationCPUTime: 8.25s
	% Computational Cost: add. (205352->118), mult. (311192->202), div. (12992->8), fcn. (197874->16), ass. (0->115)
	t367 = pkin(5) ^ 2;
	t363 = sin(pkin(16));
	t365 = cos(qJ(2));
	t432 = sin(qJ(2));
	t434 = cos(pkin(16));
	t348 = t432 * t363 - t365 * t434;
	t428 = pkin(5) * t348;
	t442 = -2 * pkin(1);
	t409 = (pkin(1) ^ 2) + t428 * t442;
	t340 = t367 + t409;
	t338 = 0.1e1 / t340;
	t370 = 0.1e1 / pkin(2);
	t421 = t338 * t370;
	t337 = pkin(2) ^ 2 - pkin(6) ^ 2 + t340;
	t345 = pkin(1) - t428;
	t439 = -pkin(6) - pkin(2);
	t335 = (pkin(5) - t439) * (pkin(5) + t439) + t409;
	t438 = -pkin(6) + pkin(2);
	t336 = (pkin(5) - t438) * (pkin(5) + t438) + t409;
	t372 = sqrt(-t336 * t335);
	t350 = t365 * t363 + t432 * t434;
	t427 = pkin(5) * t350;
	t331 = t337 * t427 + t345 * t372;
	t361 = sin(qJ(3));
	t424 = t331 * t361;
	t420 = t350 * t372;
	t330 = -pkin(5) * t420 + t337 * t345;
	t364 = cos(qJ(3));
	t425 = t330 * t364;
	t328 = (-t425 / 0.2e1 + t424 / 0.2e1) * t421;
	t423 = t331 * t364;
	t426 = t330 * t361;
	t329 = (t423 / 0.2e1 + t426 / 0.2e1) * t421;
	t358 = pkin(18) + pkin(19);
	t356 = sin(t358);
	t357 = cos(t358);
	t393 = t328 * t357 + t329 * t356;
	t382 = pkin(3) * t393;
	t321 = pkin(4) + t382;
	t452 = t321 / 0.2e1;
	t359 = sin(pkin(17));
	t451 = t359 / 0.2e1;
	t360 = cos(pkin(17));
	t450 = -t360 / 0.2e1;
	t392 = t328 * t356 - t329 * t357;
	t430 = pkin(3) * t392;
	t449 = -t430 / 0.2e1;
	t369 = pkin(3) ^ 2;
	t440 = pkin(3) * pkin(4);
	t417 = -0.2e1 * t440;
	t410 = pkin(4) ^ 2 - t393 * t417;
	t320 = t369 + t410;
	t443 = pkin(8) ^ 2;
	t318 = pkin(10) ^ 2 + t320 - t443;
	t436 = -pkin(8) + pkin(10);
	t437 = -pkin(8) - pkin(10);
	t373 = sqrt(-((pkin(3) - t436) * (pkin(3) + t436) + t410) * ((pkin(3) - t437) * (pkin(3) + t437) + t410));
	t429 = pkin(3) * t373;
	t314 = t318 * t321 - t392 * t429;
	t315 = t318 * t430 + t321 * t373;
	t366 = 0.1e1 / pkin(10);
	t441 = 0.1e1 / t320;
	t407 = t441 / 0.2e1;
	t401 = t359 * t407;
	t310 = (t314 * t441 * t450 + t315 * t401) * t366;
	t311 = (t315 * t360 * t407 + t314 * t401) * t366;
	t362 = sin(qJ(1));
	t383 = t365 * t361 + t432 * t364;
	t341 = t383 * t362;
	t347 = t432 * t361 - t365 * t364;
	t342 = t347 * t362;
	t416 = t341 * t310 - t342 * t311;
	t415 = -t342 * t310 - t341 * t311;
	t433 = cos(qJ(1));
	t396 = t433 * t432;
	t403 = t365 * t433;
	t343 = -t361 * t396 + t364 * t403;
	t344 = t361 * t403 + t364 * t396;
	t414 = t343 * t310 - t344 * t311;
	t413 = t344 * t310 + t343 * t311;
	t397 = 0.1e1 / t320 ^ 2 * t366 * t440;
	t394 = t315 * t397;
	t395 = t314 * t397;
	t448 = t359 * t395 + t360 * t394;
	t447 = t359 * t394 - t360 * t395;
	t446 = -0.2e1 * t369 * t392 * pkin(4) - t429;
	t445 = -pkin(3) * t318 + t321 * t417;
	t444 = 0.4e1 / t373 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t410 - t443) * t440;
	t435 = t361 / 0.2e1;
	t406 = pkin(1) * t427;
	t422 = 0.2e1 / t372 * (t335 + t336) * t406;
	t418 = t366 * t441;
	t412 = t347 * t310 + t311 * t383;
	t411 = t310 * t383 - t347 * t311;
	t400 = t366 * t407;
	t398 = 0.1e1 / t340 ^ 2 * t406;
	t391 = t424 - t425;
	t390 = t423 + t426;
	t326 = (t348 * t372 + (-t422 / 0.2e1 - t337 + t345 * t442) * t350) * pkin(5);
	t327 = t345 * t422 / 0.2e1 + t367 * t350 ^ 2 * t442 + (-t337 * t348 - t420) * pkin(5);
	t385 = -t326 * t364 / 0.2e1 + t327 * t435;
	t384 = t327 * t364 / 0.2e1 + t326 * t435;
	t316 = ((-t356 * t384 - t357 * t385) * t338 + (-t356 * t390 - t357 * t391) * t398) * t370;
	t380 = t316 * t444;
	t379 = t392 * t444;
	t378 = (-t373 * t382 + t379 * t449 + t445 * t392) * t418;
	t377 = (t318 * t382 + t379 * t452 + t446 * t392) * t400;
	t376 = pkin(3) * t370 * ((t356 * t385 - t357 * t384) * t338 + (t356 * t391 - t357 * t390) * t398);
	t375 = (t445 * t316 - t373 * t376 + t380 * t449) * t418;
	t374 = (t446 * t316 + t318 * t376 + t380 * t452) * t400;
	t293 = t360 * t377 + t378 * t451 + t448 * t392;
	t292 = t359 * t377 + t378 * t450 + t447 * t392;
	t291 = t448 * t316 + t360 * t374 + t375 * t451;
	t290 = t447 * t316 + t359 * t374 + t375 * t450;
	t1 = [t415, -t290 * t343 + t291 * t344 + t413, -t292 * t343 + t293 * t344 + t413, 0; -t414, t290 * t342 + t291 * t341 + t416, t292 * t342 + t293 * t341 + t416, 0; 0, -t290 * t383 + t291 * t347 + t412, -t292 * t383 + t293 * t347 + t412, 0; -t416, t290 * t344 + t291 * t343 + t414, t292 * t344 + t293 * t343 + t414, 0; t413, t290 * t341 - t291 * t342 + t415, t292 * t341 - t293 * t342 + t415, 0; 0, t290 * t347 + t291 * t383 + t411, t292 * t347 + t293 * t383 + t411, 0; t433, 0, 0, 0; t362, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:54:06
	% EndTime: 2020-04-18 09:54:37
	% DurationCPUTime: 18.92s
	% Computational Cost: add. (324850->131), mult. (492070->219), div. (20688->8), fcn. (312822->18), ass. (0->128)
	t632 = sin(qJ(2));
	t633 = sin(pkin(16));
	t635 = cos(qJ(2));
	t636 = cos(pkin(16));
	t543 = t632 * t636 + t635 * t633;
	t604 = t632 * t633 - t635 * t636;
	t603 = pkin(5) * t604;
	t601 = (-0.2e1 * t603 + pkin(1)) * pkin(1);
	t639 = -pkin(6) + pkin(2);
	t640 = -pkin(6) - pkin(2);
	t555 = sqrt(-((pkin(5) - t639) * (pkin(5) + t639) + t601) * ((pkin(5) - t640) * (pkin(5) + t640) + t601));
	t648 = pkin(5) ^ 2;
	t536 = t601 + t648;
	t556 = pkin(2) ^ 2;
	t656 = -pkin(6) ^ 2 + t556;
	t600 = t536 + t656;
	t599 = pkin(5) * t600;
	t602 = -t603 + pkin(1);
	t649 = 0.1e1 / pkin(2);
	t598 = t649 * (t543 * t599 + t602 * t555);
	t642 = 0.1e1 / t536;
	t595 = t642 * t598 / 0.2e1;
	t628 = pkin(5) * t543;
	t597 = t649 * (-t555 * t628 + t602 * t600);
	t596 = t642 * t597;
	t631 = sin(qJ(3));
	t634 = cos(qJ(3));
	t591 = -t634 * t596 / 0.2e1 + t631 * t595;
	t658 = t631 / 0.2e1;
	t592 = t634 * t595 + t596 * t658;
	t614 = pkin(18) + pkin(19);
	t609 = sin(t614);
	t610 = cos(t614);
	t585 = t610 * t591 + t609 * t592;
	t583 = pkin(3) * t585;
	t579 = (0.2e1 * t583 + pkin(4)) * pkin(4);
	t637 = -pkin(8) + pkin(10);
	t638 = -pkin(8) - pkin(10);
	t554 = sqrt(-((pkin(3) - t637) * (pkin(3) + t637) + t579) * ((pkin(3) - t638) * (pkin(3) + t638) + t579));
	t644 = -0.2e1 * pkin(4);
	t529 = pkin(4) ^ 2 + (-t585 * t644 + pkin(3)) * pkin(3);
	t646 = pkin(8) ^ 2;
	t577 = pkin(10) ^ 2 + t529 - t646;
	t576 = pkin(3) * t577;
	t580 = t583 + pkin(4);
	t645 = 0.1e1 / pkin(10);
	t650 = t609 * t591 - t610 * t592;
	t571 = t645 * (t580 * t554 + t576 * t650);
	t643 = 0.1e1 / t529;
	t568 = t643 * t571 / 0.2e1;
	t582 = pkin(3) * t650;
	t570 = t645 * (-t554 * t582 + t580 * t577);
	t569 = t643 * t570;
	t625 = sin(pkin(17));
	t626 = cos(pkin(17));
	t659 = -t626 / 0.2e1;
	t523 = t625 * t568 + t569 * t659;
	t660 = t625 / 0.2e1;
	t524 = t626 * t568 + t569 * t660;
	t542 = -t635 * t631 - t632 * t634;
	t551 = sin(qJ(1));
	t537 = t542 * t551;
	t541 = t632 * t631 - t635 * t634;
	t538 = t541 * t551;
	t511 = -t523 * t538 + t524 * t537;
	t550 = sin(qJ(4));
	t552 = cos(qJ(4));
	t553 = cos(qJ(1));
	t664 = t511 * t552 + t553 * t550;
	t663 = t511 * t550 - t552 * t553;
	t662 = t580 / 0.2e1;
	t661 = -t582 / 0.2e1;
	t641 = pkin(3) * pkin(4);
	t657 = 0.1e1 / t529 ^ 2 * t641;
	t618 = -t537 * t523 - t538 * t524;
	t539 = t541 * t553;
	t540 = t542 * t553;
	t617 = -t540 * t523 - t539 * t524;
	t566 = t570 * t657;
	t567 = t571 * t657;
	t655 = -t626 * t566 + t625 * t567;
	t654 = t625 * t566 + t626 * t567;
	t653 = -0.2e1 * t580 * t641 - t576;
	t629 = pkin(3) * t554;
	t652 = pkin(3) ^ 2 * t644 * t650 - t629;
	t651 = 0.4e1 / t554 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t579 - t646) * t641;
	t630 = t649 * t642;
	t627 = t645 * t643;
	t613 = pkin(1) * t628;
	t624 = 0.4e1 / t555 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t601 - t556) * t613;
	t531 = (t604 * t555 + (-t624 / 0.2e1 - t648 - t656) * t543) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t603) * t613;
	t532 = t602 * t624 / 0.2e1 - t604 * t599 + (-t555 - 0.2e1 * t613) * t628;
	t535 = 0.1e1 / t536 ^ 2;
	t593 = t597 * t613;
	t594 = t598 * t613;
	t605 = t630 * t658;
	t608 = t634 * t630;
	t587 = t532 * t608 / 0.2e1 + t531 * t605 + (t631 * t593 + t634 * t594) * t535;
	t588 = -t531 * t608 / 0.2e1 + t532 * t605 + (-t634 * t593 + t631 * t594) * t535;
	t527 = -t609 * t587 - t610 * t588;
	t572 = pkin(3) * (-t610 * t587 + t609 * t588);
	t574 = t527 * t651;
	t611 = t627 / 0.2e1;
	t558 = (t652 * t527 + t577 * t572 + t574 * t662) * t611;
	t559 = (t653 * t527 - t554 * t572 + t574 * t661) * t627;
	t503 = t655 * t527 + t625 * t558 + t559 * t659;
	t622 = -t503 + t524;
	t504 = t654 * t527 + t626 * t558 + t559 * t660;
	t621 = t504 + t523;
	t573 = t650 * t651;
	t560 = (t573 * t662 + t585 * t576 + t650 * t652) * t611;
	t561 = (t573 * t661 - t585 * t629 + t650 * t653) * t627;
	t505 = t625 * t560 + t561 * t659 + t650 * t655;
	t620 = -t505 + t524;
	t506 = t626 * t560 + t561 * t660 + t650 * t654;
	t619 = t506 + t523;
	t616 = t541 * t523 - t542 * t524;
	t514 = t523 * t542 + t524 * t541;
	t513 = t523 * t539 - t524 * t540;
	t510 = t513 * t552 + t550 * t551;
	t509 = -t513 * t550 + t551 * t552;
	t502 = t505 * t542 + t506 * t541 + t616;
	t501 = t505 * t539 - t506 * t540 + t617;
	t500 = t505 * t538 - t506 * t537 + t618;
	t499 = t503 * t542 + t504 * t541 + t616;
	t498 = t503 * t539 - t504 * t540 + t617;
	t497 = t503 * t538 - t504 * t537 + t618;
	t1 = [t664, t498 * t552, t501 * t552, t509; t510, t497 * t552, t500 * t552, t663; 0, t499 * t552, t502 * t552, -t514 * t550; -t663, -t498 * t550, -t501 * t550, -t510; t509, -t497 * t550, -t500 * t550, t664; 0, -t499 * t550, -t502 * t550, -t514 * t552; t618, t621 * t539 - t622 * t540, t619 * t539 - t620 * t540, 0; -t617, -t622 * t537 + t621 * t538, -t620 * t537 + t619 * t538, 0; 0, t622 * t541 + t621 * t542, t620 * t541 + t619 * t542, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:04
	% EndTime: 2020-04-18 09:52:04
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (1338->36), mult. (2084->64), div. (64->4), fcn. (1354->10), ass. (0->41)
	t135 = pkin(1) ^ 2;
	t127 = sin(qJ(2));
	t129 = sin(pkin(16));
	t149 = cos(qJ(2));
	t150 = cos(pkin(16));
	t124 = t127 * t129 - t149 * t150;
	t148 = pkin(1) * t124;
	t155 = -2 * pkin(5);
	t143 = t148 * t155 + t135;
	t121 = (pkin(5) ^ 2) + t143;
	t118 = -pkin(2) ^ 2 + pkin(6) ^ 2 + t121;
	t122 = -pkin(5) + t148;
	t154 = -pkin(6) - pkin(2);
	t116 = (pkin(5) - t154) * (pkin(5) + t154) + t143;
	t153 = -pkin(6) + pkin(2);
	t117 = (pkin(5) - t153) * (pkin(5) + t153) + t143;
	t136 = sqrt(-t117 * t116);
	t125 = t127 * t150 + t149 * t129;
	t147 = t125 * pkin(1);
	t112 = t118 * t147 - t122 * t136;
	t133 = 0.1e1 / pkin(6);
	t142 = pkin(5) * t147;
	t141 = 0.1e1 / t121 ^ 2 * t133 * t142;
	t144 = t125 * t136;
	t145 = 0.1e1 / t121 * t133;
	t146 = 0.1e1 / t136 * (t116 + t117) * t142;
	t156 = (-t122 * t146 + t135 * t125 ^ 2 * t155 + (-t124 * t118 - t144) * pkin(1)) * t145 / 0.2e1 + t112 * t141;
	t130 = sin(pkin(15));
	t152 = t130 / 0.2e1;
	t132 = cos(pkin(15));
	t151 = t132 / 0.2e1;
	t111 = -pkin(1) * t144 - t122 * t118;
	t140 = t111 * t141;
	t138 = pkin(1) * (t124 * t136 + (0.2e1 * t122 * pkin(5) - t118 - t146) * t125) * t145;
	t131 = cos(qJ(1));
	t128 = sin(qJ(1));
	t109 = (-t112 * t132 / 0.2e1 + t111 * t152) * t145;
	t108 = (t111 * t151 + t112 * t152) * t145;
	t105 = t156 * t130 + t132 * t140 + t138 * t151;
	t104 = t156 * t132 + (-t138 / 0.2e1 - t140) * t130;
	t1 = [-t128 * t108, t131 * t105, 0, 0; t131 * t108, t128 * t105, 0, 0; 0, t104, 0, 0; -t128 * t109, -t131 * t104, 0, 0; t131 * t109, -t128 * t104, 0, 0; 0, t105, 0, 0; t131, 0, 0, 0; t128, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:04
	% EndTime: 2020-04-18 09:52:06
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (3354->44), mult. (5216->72), div. (176->4), fcn. (3394->10), ass. (0->47)
	t125 = sin(qJ(2));
	t127 = sin(pkin(16));
	t128 = cos(qJ(2));
	t130 = cos(pkin(16));
	t150 = t125 * t127 - t128 * t130;
	t149 = pkin(5) * t150;
	t145 = (-0.2e1 * t149 + pkin(1)) * pkin(1);
	t165 = pkin(5) ^ 2;
	t123 = t145 + t165;
	t164 = 0.1e1 / t123;
	t172 = t164 / 0.2e1;
	t166 = 0.1e1 / pkin(2);
	t171 = t166 * t172;
	t124 = t125 * t130 + t128 * t127;
	t160 = pkin(5) * t124;
	t154 = pkin(1) * t160;
	t170 = 0.1e1 / t123 ^ 2 * t154;
	t162 = -pkin(6) + pkin(2);
	t163 = -pkin(6) - pkin(2);
	t131 = sqrt(-((pkin(5) - t162) * (pkin(5) + t162) + t145) * ((pkin(5) - t163) * (pkin(5) + t163) + t145));
	t132 = pkin(2) ^ 2;
	t144 = -pkin(6) ^ 2 + t123 + t132;
	t143 = pkin(5) * t144;
	t142 = t124 * t143;
	t148 = -t149 + pkin(1);
	t140 = t166 * (t148 * t131 + t142);
	t141 = 0.2e1 / t131 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t145 - t132) * t154;
	t159 = pkin(5) * t131;
	t153 = t124 * t159;
	t169 = (-0.2e1 * t165 * t124 ^ 2 * pkin(1) + t148 * t141 - t150 * t143 - t153) * t171 + t140 * t170;
	t139 = t166 * (t148 * t144 - t153);
	t168 = (-t141 * t160 - 0.2e1 * t148 * t154 + t150 * t159 - t142) * t171 + t139 * t170;
	t137 = t140 * t172;
	t138 = t164 * t139;
	t157 = sin(pkin(19));
	t158 = cos(pkin(19));
	t119 = t158 * t137 + t157 * t138 / 0.2e1;
	t155 = t169 * t157 - t168 * t158 - t119;
	t118 = -t158 * t138 / 0.2e1 + t157 * t137;
	t156 = t168 * t157 + t169 * t158 + t118;
	t167 = t155 * t125 + t156 * t128;
	t152 = t118 * t128 - t119 * t125;
	t151 = t118 * t125 + t119 * t128;
	t146 = -t156 * t125 + t155 * t128;
	t129 = cos(qJ(1));
	t126 = sin(qJ(1));
	t1 = [-t152 * t126, t146 * t129, 0, 0; t152 * t129, t146 * t126, 0, 0; 0, t167, 0, 0; t151 * t126, -t167 * t129, 0, 0; -t151 * t129, -t167 * t126, 0, 0; 0, t146, 0, 0; t129, 0, 0, 0; t126, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-18 09:52:31
	% EndTime: 2020-04-18 09:52:51
	% DurationCPUTime: 12.34s
	% Computational Cost: add. (197868->131), mult. (300332->231), div. (12336->8), fcn. (191114->18), ass. (0->129)
	t415 = sin(qJ(2));
	t416 = sin(pkin(16));
	t418 = cos(qJ(2));
	t419 = cos(pkin(16));
	t387 = t415 * t416 - t418 * t419;
	t385 = pkin(5) * t387;
	t383 = (-0.2e1 * t385 + pkin(1)) * pkin(1);
	t422 = -pkin(6) + pkin(2);
	t423 = -pkin(6) - pkin(2);
	t337 = sqrt(-((pkin(5) - t422) * (pkin(5) + t422) + t383) * ((pkin(5) - t423) * (pkin(5) + t423) + t383));
	t430 = pkin(5) ^ 2;
	t330 = t383 + t430;
	t338 = pkin(2) ^ 2;
	t439 = -pkin(6) ^ 2 + t338;
	t380 = t330 + t439;
	t384 = -t385 + pkin(1);
	t331 = t415 * t419 + t418 * t416;
	t411 = pkin(5) * t331;
	t432 = 0.1e1 / pkin(2);
	t376 = t432 * (-t337 * t411 + t384 * t380);
	t425 = 0.1e1 / t330;
	t375 = t425 * t376;
	t373 = -t375 / 0.2e1;
	t378 = pkin(5) * t380;
	t377 = t432 * (t331 * t378 + t384 * t337);
	t374 = t425 * t377 / 0.2e1;
	t406 = sin(pkin(19));
	t407 = cos(pkin(19));
	t323 = t407 * t373 + t406 * t374;
	t372 = t375 / 0.2e1;
	t324 = t406 * t372 + t407 * t374;
	t315 = t418 * t323 - t415 * t324;
	t402 = pkin(1) * t411;
	t405 = 0.4e1 / t337 * ((pkin(5) + pkin(6)) * (pkin(5) - pkin(6)) + t383 - t338) * t402;
	t318 = (t387 * t337 + (-t405 / 0.2e1 - t430 - t439) * t331) * pkin(5) + (-0.3e1 * pkin(1) + 0.4e1 * t385) * t402;
	t446 = -t318 / 0.2e1;
	t319 = t384 * t405 / 0.2e1 - t387 * t378 + (-t337 - 0.2e1 * t402) * t411;
	t445 = t319 / 0.2e1;
	t333 = sin(pkin(18));
	t444 = t333 / 0.2e1;
	t414 = sin(qJ(3));
	t417 = cos(qJ(3));
	t368 = t417 * t373 + t414 * t374;
	t369 = t414 * t372 + t417 * t374;
	t403 = pkin(18) + pkin(19);
	t394 = sin(t403);
	t395 = cos(t403);
	t362 = t395 * t368 + t394 * t369;
	t360 = pkin(4) * t362;
	t357 = -t360 - pkin(3);
	t443 = -t357 / 0.2e1;
	t433 = t394 * t368 - t395 * t369;
	t359 = pkin(4) * t433;
	t442 = -t359 / 0.2e1;
	t408 = cos(pkin(18));
	t441 = -t408 / 0.2e1;
	t431 = pkin(4) ^ 2;
	t308 = t431 + (0.2e1 * t360 + pkin(3)) * pkin(3);
	t424 = pkin(3) * pkin(4);
	t440 = 0.1e1 / t308 ^ 2 * t424;
	t329 = 0.1e1 / t330 ^ 2;
	t370 = t376 * t402;
	t371 = t377 * t402;
	t413 = t432 * t425;
	t400 = t413 / 0.2e1;
	t388 = t406 * t400;
	t392 = t407 * t413;
	t305 = t392 * t445 + t318 * t388 + (t406 * t370 + t407 * t371) * t329;
	t306 = t392 * t446 + t319 * t388 + (-t407 * t370 + t406 * t371) * t329;
	t302 = t418 * t305 + t415 * t306 + t315;
	t339 = pkin(8) ^ 2;
	t354 = -pkin(10) ^ 2 + t308 + t339;
	t353 = pkin(4) * t354;
	t438 = 0.2e1 * t357 * t424 - t353;
	t427 = -0.2e1 * pkin(3);
	t356 = (-t362 * t427 + pkin(4)) * pkin(4);
	t420 = -pkin(8) + pkin(10);
	t421 = -pkin(8) - pkin(10);
	t336 = sqrt(-((pkin(3) - t420) * (pkin(3) + t420) + t356) * ((pkin(3) - t421) * (pkin(3) + t421) + t356));
	t428 = 0.1e1 / pkin(8);
	t348 = t428 * (-t357 * t336 + t353 * t433);
	t346 = t348 * t440;
	t297 = -t336 * t359 - t357 * t354;
	t409 = t297 * t428;
	t391 = t409 * t440;
	t437 = -t333 * t346 - t408 * t391;
	t436 = t333 * t391 - t408 * t346;
	t412 = pkin(4) * t336;
	t435 = t427 * t431 * t433 - t412;
	t434 = 0.4e1 / t336 * ((pkin(3) + pkin(10)) * (pkin(3) - pkin(10)) + t356 - t339) * t424;
	t426 = 0.1e1 / t308;
	t410 = t428 * t426;
	t399 = -t410 / 0.2e1;
	t398 = t418 * t324;
	t397 = t415 * t323;
	t393 = t417 * t413;
	t390 = t414 * t400;
	t314 = t397 + t398;
	t381 = -t415 * t305 + t418 * t306 - t398;
	t303 = -t397 + t381;
	t365 = t393 * t446 + t319 * t390 + (-t417 * t370 + t414 * t371) * t329;
	t364 = t393 * t445 + t318 * t390 + (t414 * t370 + t417 * t371) * t329;
	t304 = -t394 * t364 - t395 * t365;
	t351 = t304 * t434;
	t350 = t433 * t434;
	t349 = pkin(4) * (-t395 * t364 + t394 * t365);
	t347 = -t426 * t348 / 0.2e1;
	t343 = (t350 * t442 - t362 * t412 + t433 * t438) * t410;
	t342 = (t350 * t443 + t362 * t353 + t433 * t435) * t399;
	t341 = (t438 * t304 - t336 * t349 + t351 * t442) * t410;
	t340 = (t435 * t304 + t354 * t349 + t351 * t443) * t399;
	t335 = cos(qJ(1));
	t334 = sin(qJ(1));
	t320 = t334 * t397;
	t312 = t314 * t335;
	t311 = t315 * t335;
	t310 = -t334 * t398 - t320;
	t309 = t315 * t334;
	t301 = t302 * t335;
	t300 = t303 * t335;
	t299 = t302 * t334;
	t298 = t381 * t334 - t320;
	t292 = t408 * t297 * t399 + t333 * t347;
	t291 = t426 * t409 * t444 + t408 * t347;
	t286 = t333 * t342 + t343 * t441 + t433 * t437;
	t285 = t408 * t342 + t343 * t444 + t433 * t436;
	t284 = t437 * t304 + t333 * t340 + t341 * t441;
	t283 = t436 * t304 + t408 * t340 + t341 * t444;
	t1 = [-t291 * t310 - t292 * t309, -t283 * t312 + t284 * t311 - t291 * t301 + t292 * t300, -t285 * t312 + t286 * t311, 0; -t291 * t312 + t292 * t311, t283 * t310 + t284 * t309 - t291 * t299 + t292 * t298, t285 * t310 + t286 * t309, 0; 0, t283 * t315 + t284 * t314 + t291 * t303 + t292 * t302, t285 * t315 + t286 * t314, 0; t291 * t309 - t292 * t310, -t283 * t311 - t284 * t312 - t291 * t300 - t292 * t301, -t285 * t311 - t286 * t312, 0; -t291 * t311 - t292 * t312, -t283 * t309 + t284 * t310 - t291 * t298 - t292 * t299, -t285 * t309 + t286 * t310, 0; 0, -t283 * t314 + t284 * t315 - t291 * t302 + t292 * t303, -t285 * t314 + t286 * t315, 0; t335, 0, 0, 0; t334, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end