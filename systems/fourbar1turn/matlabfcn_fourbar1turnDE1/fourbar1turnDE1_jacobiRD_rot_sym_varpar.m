% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% JRD_rot [9x2]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbar1turnDE1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_jacobiRD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_jacobiRD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
JRD_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0; -t31, 0; 0, 0; t31, 0; -t30, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t31 = sin(qJ(1));
	t38 = qJD(1) * t31;
	t33 = cos(qJ(1));
	t37 = qJD(1) * t33;
	t30 = sin(qJ(2));
	t36 = qJD(2) * t30;
	t32 = cos(qJ(2));
	t35 = qJD(2) * t32;
	t34 = qJD(2) * t33;
	t29 = t31 * t36 - t32 * t37;
	t28 = t30 * t37 + t31 * t35;
	t27 = t30 * t34 + t32 * t38;
	t26 = t30 * t38 - t32 * t34;
	t1 = [t29, t26; -t27, -t28; 0, -t36; t28, t27; t26, t29; 0, -t35; -t38, 0; t37, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:11
	% EndTime: 2020-04-12 19:28:14
	% DurationCPUTime: 3.90s
	% Computational Cost: add. (59267->133), mult. (83690->293), div. (2972->12), fcn. (22834->8), ass. (0->137)
	t515 = sin(qJ(2));
	t522 = pkin(2) ^ 2;
	t523 = pkin(1) ^ 2;
	t517 = cos(qJ(2));
	t627 = pkin(1) * t517;
	t601 = -0.2e1 * pkin(2) * t627 + t523;
	t511 = t522 + t601;
	t526 = pkin(3) ^ 2;
	t507 = -pkin(4) ^ 2 + t511 + t526;
	t628 = pkin(1) * t507;
	t504 = t515 * t628;
	t512 = -pkin(2) + t627;
	t633 = -0.2e1 * t512;
	t596 = pkin(2) * t633;
	t631 = -pkin(3) - pkin(4);
	t505 = (pkin(2) - t631) * (pkin(2) + t631) + t601;
	t630 = pkin(4) - pkin(3);
	t506 = (pkin(2) - t630) * (pkin(2) + t630) + t601;
	t613 = t505 * t506;
	t525 = sqrt(-t613);
	t604 = t525 * t517;
	t632 = pkin(1) * pkin(2);
	t573 = (-t505 - t506) * t632;
	t498 = t515 * t573;
	t501 = 0.1e1 / t525;
	t615 = t501 * t498;
	t481 = t504 + (-t604 + (t596 - t615) * t515) * pkin(1);
	t608 = t515 * t525;
	t592 = pkin(1) * t608;
	t494 = -t507 * t512 - t592;
	t508 = 0.1e1 / t511;
	t520 = 0.1e1 / pkin(3);
	t509 = 0.1e1 / t511 ^ 2;
	t593 = -0.2e1 * t509 * t515;
	t569 = t593 * t632;
	t471 = (t481 * t508 + t494 * t569) * t520;
	t491 = 0.1e1 / t494 ^ 2;
	t495 = -t512 * t525 + t504;
	t621 = t491 * t495;
	t587 = t471 * t621;
	t493 = t495 ^ 2;
	t488 = t491 * t493 + 0.1e1;
	t486 = 0.1e1 / t488;
	t625 = pkin(3) * t511;
	t591 = t486 * t625;
	t514 = t515 ^ 2;
	t611 = t514 * t523;
	t614 = t501 * t512;
	t482 = -t498 * t614 + 0.2e1 * pkin(2) * t611 + (t507 * t517 + t608) * pkin(1);
	t472 = (t482 * t508 + t495 * t569) * t520;
	t490 = 0.1e1 / t494;
	t624 = t472 * t490;
	t461 = (-t587 + t624) * t591;
	t629 = t461 + 0.1e1;
	t497 = qJD(2) * t498;
	t616 = t501 * t497;
	t479 = (-t515 * t616 + (-t604 + (t507 + t596) * t515) * qJD(2)) * pkin(1);
	t580 = qJD(2) * t611;
	t574 = pkin(2) * t580;
	t597 = qJD(2) * t517;
	t602 = qJD(2) * t592 + t597 * t628;
	t480 = -t497 * t614 + 0.2e1 * t574 + t602;
	t496 = (t517 * t573 - 0.4e1 * t522 * t611) * qJD(2);
	t510 = t508 * t509;
	t557 = 0.8e1 * t510 * t522 * t580;
	t563 = 0.1e1 / t613 * t497 * t615;
	t568 = 0.2e1 * t471 * t495 * t625;
	t598 = qJD(2) * t515;
	t589 = pkin(2) * t598;
	t576 = pkin(1) * t589;
	t590 = t490 * t625;
	t594 = 0.2e1 * t509;
	t606 = t517 * t508;
	t607 = t517 * t495;
	t619 = t494 * t517;
	t489 = t494 ^ 2;
	t623 = t479 * t490 / t489;
	t626 = pkin(2) * t509;
	t460 = (-0.2e1 * t472 * t590 + t491 * t568) * (t480 * t621 - t493 * t623) / t488 ^ 2 + (-t472 * t479 - ((0.4e1 * t574 + t602) * t508 + t494 * t557 + ((-0.2e1 * t517 * t616 + (-t496 * t501 - t563) * t515) * t508 + (t479 * t593 + (t606 * t633 + (-t481 * t515 - t619) * t594) * qJD(2)) * pkin(2)) * pkin(1)) * t520 * t495 - t471 * t480) * t491 * t591 + ((0.2e1 * t624 - 0.2e1 * t587) * pkin(3) * t576 + ((0.6e1 * t517 * t523 * t589 - t496 * t614 - t512 * t563) * t508 + t495 * t557 + ((-0.2e1 * t480 * t626 + t508 * t616) * t515 + ((t604 + (-t507 + t615) * t515) * t508 + (-t482 * t515 - t607) * pkin(2) * t594) * qJD(2)) * pkin(1)) * t520 * t590 + t568 * t623) * t486;
	t636 = qJD(2) * t629;
	t538 = t460 * t494 + t461 * t479 - t495 * t636;
	t537 = t479 + t538;
	t536 = t460 * t495 + t461 * t480 + t494 * t636 + t480;
	t643 = t536 * t515;
	t646 = -t517 * t537 + t643;
	t610 = t515 * t494;
	t551 = -t495 * t514 + t517 * t610;
	t645 = 0.2e1 * t551;
	t582 = t520 * t606;
	t521 = 0.1e1 / t526;
	t603 = t489 + t493;
	t485 = t603 * t521 * t509;
	t483 = t485 ^ (-0.1e1 / 0.2e1);
	t622 = t483 * t520;
	t640 = ((t479 * t494 + t480 * t495) * t594 - 0.4e1 * t603 * t510 * t576) * t521 * t483 / t485;
	t642 = t569 * t597 * t622 - t582 * t640 / 0.2e1;
	t577 = pkin(1) * qJD(2) * t626;
	t612 = t508 * t520;
	t583 = t515 * t612;
	t641 = 0.2e1 * t514 * t577 * t622 + t583 * t640 / 0.2e1;
	t516 = sin(qJ(1));
	t600 = qJD(1) * t516;
	t635 = t629 * t600;
	t609 = t515 * t495;
	t560 = -t609 + t619;
	t634 = qJD(2) * t560 + t479 * t515 + t480 * t517;
	t620 = t494 * t514;
	t617 = t495 * t516;
	t518 = cos(qJ(1));
	t599 = qJD(1) * t518;
	t595 = 0.2e1 * t461;
	t588 = t508 * t640;
	t586 = t483 * t612;
	t584 = t515 * t607;
	t581 = t495 * t598;
	t572 = t516 * t588;
	t571 = t518 * t588;
	t570 = t483 * t582;
	t567 = t516 * t577;
	t566 = t518 * t577;
	t564 = t518 * t570;
	t559 = t607 + t610;
	t556 = t494 * t564;
	t555 = qJD(1) * t556 + t641 * t617 + (t479 * t570 + t494 * t642) * t516;
	t553 = t584 + t620;
	t549 = t609 / 0.2e1 - t619 / 0.2e1;
	t548 = t607 / 0.2e1 + t610 / 0.2e1;
	t546 = qJD(1) * t559;
	t545 = qJD(2) * t556 + t480 * t564 + (t479 * t483 * t583 - t494 * t641 + t495 * t642) * t518;
	t544 = t549 * t461;
	t543 = t548 * t461;
	t540 = t553 * t595;
	t539 = t551 * t595;
	t534 = (-t543 - t548) * t588;
	t533 = (-t540 - 0.2e1 * t584 - 0.2e1 * t620) * t577;
	t532 = t515 * t537 + t517 * t536;
	t1 = [(-t597 * t617 + (-t495 * t599 + (-qJD(2) * t494 - t480) * t516) * t515) * t586 + t555, (-t543 * t571 + (-t540 * t566 + (-t559 * t635 + (t559 * t460 + t634 * t461 - t581) * t518) * t508) * t483) * t520 + t545; (-t549 * t571 + (t645 * t566 + (t560 * t600 + (qJD(2) * t559 - t479 * t517 + t480 * t515) * t518) * t508) * t483) * t520, (t516 * t534 + (t516 * t533 + (t629 * t599 * t559 + t532 * t516) * t508) * t483) * t520; 0, ((t629 * t577 * t645 + t646 * t508) * t483 - t629 * t588 * t549) * t520; (t548 * t572 + (0.2e1 * t553 * t567 + (-t634 * t516 - t518 * t546) * t508) * t483) * t520, ((t544 + t549) * t571 + ((-t539 - t645) * t566 + (-t646 * t518 - t560 * t635) * t508) * t483) * t520; (-t516 * t546 - t518 * t581) * t586 + t545, (t544 * t572 + (-t539 * t567 + ((t461 * t560 - t609) * t599 + (t517 * t538 - t643) * t516) * t508) * t483) * t520 + t555; 0, (t534 + (t532 * t508 + t533) * t483) * t520; -t600, 0; t599, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:07
	% EndTime: 2020-04-12 19:28:08
	% DurationCPUTime: 1.73s
	% Computational Cost: add. (24213->87), mult. (34040->218), div. (1126->12), fcn. (9209->8), ass. (0->102)
	t411 = pkin(2) ^ 2;
	t412 = pkin(1) ^ 2;
	t406 = cos(qJ(2));
	t475 = pkin(2) * t406;
	t456 = -0.2e1 * pkin(1) * t475 + t412;
	t400 = t411 + t456;
	t415 = pkin(4) ^ 2;
	t396 = -pkin(3) ^ 2 + t400 + t415;
	t404 = sin(qJ(2));
	t393 = pkin(2) * t404 * t396;
	t401 = pkin(1) - t475;
	t482 = 0.2e1 * pkin(1);
	t451 = t401 * t482;
	t480 = -pkin(3) - pkin(4);
	t394 = (pkin(2) - t480) * (pkin(2) + t480) + t456;
	t479 = -pkin(3) + pkin(4);
	t395 = (pkin(2) - t479) * (pkin(2) + t479) + t456;
	t466 = t394 * t395;
	t414 = sqrt(-t466);
	t459 = t406 * t414;
	t481 = pkin(1) * pkin(2);
	t437 = (-t394 - t395) * t481;
	t387 = t404 * t437;
	t390 = 0.1e1 / t414;
	t468 = t390 * t387;
	t370 = t393 + (-t459 + (t451 - t468) * t404) * pkin(2);
	t461 = t404 * t414;
	t383 = -pkin(2) * t461 + t401 * t396;
	t409 = 0.1e1 / pkin(4);
	t398 = 0.1e1 / t400 ^ 2;
	t463 = t398 * t404;
	t442 = t463 * t481;
	t397 = 0.1e1 / t400;
	t478 = -t397 / 0.2e1;
	t366 = (t370 * t478 + t383 * t442) * t409;
	t380 = 0.1e1 / t383 ^ 2;
	t384 = t401 * t414 + t393;
	t470 = t380 * t384;
	t460 = t406 * t396;
	t403 = t404 ^ 2;
	t462 = t403 * t411;
	t467 = t390 * t401;
	t371 = t387 * t467 + t462 * t482 + (t460 + t461) * pkin(2);
	t477 = t397 / 0.2e1;
	t367 = (t371 * t477 - t384 * t442) * t409;
	t379 = 0.1e1 / t383;
	t472 = t367 * t379;
	t491 = -0.2e1 * t366 * t470 - 0.2e1 * t472;
	t386 = qJD(2) * t387;
	t453 = qJD(2) * t411;
	t443 = t403 * t453;
	t438 = pkin(1) * t443;
	t473 = pkin(2) * qJD(2);
	t449 = t404 * t473;
	t457 = t414 * t449 + t460 * t473;
	t369 = t386 * t467 + 0.2e1 * t438 + t457;
	t410 = 0.1e1 / t415;
	t378 = t383 ^ 2;
	t382 = t384 ^ 2;
	t458 = t378 + t382;
	t374 = t458 * t410 * t398;
	t372 = t374 ^ (-0.1e1 / 0.2e1);
	t405 = sin(qJ(1));
	t441 = pkin(1) * t449;
	t487 = 0.2e1 * t398;
	t426 = t441 * t487;
	t407 = cos(qJ(1));
	t454 = qJD(1) * t407;
	t423 = -t397 * t454 + t405 * t426;
	t465 = t397 * t405;
	t490 = (-t369 * t465 + t423 * t384) * t372;
	t469 = t390 * t386;
	t368 = (-t404 * t469 + (-t459 + (t396 + t451) * t404) * qJD(2)) * pkin(2);
	t455 = qJD(1) * t405;
	t424 = t397 * t455 + t407 * t426;
	t464 = t397 * t407;
	t489 = (-t368 * t464 + t424 * t383) * t372;
	t377 = t382 * t380 + 0.1e1;
	t471 = t368 * t379 / t378;
	t486 = (t369 * t470 - t382 * t471) / t377 ^ 2;
	t485 = t372 * t397;
	t476 = pkin(1) * t398;
	t474 = pkin(4) * t400;
	t375 = 0.1e1 / t377;
	t450 = t375 * t474;
	t385 = (t406 * t437 - 0.4e1 * t412 * t462) * qJD(2);
	t399 = t397 * t398;
	t430 = t399 * t412 * t443;
	t436 = 0.4e1 / t466 * t386 * t468;
	t448 = 0.2e1 * (pkin(4) * t375 * t441 * t491 + 0.2e1 * (t472 * t486 + (t375 * t471 + t380 * t486) * t366 * t384) * t474 + ((-t366 * t369 + t367 * t368) * t380 + (-((t401 * t436 / 0.4e1 + t385 * t467 + 0.6e1 * t406 * t404 * pkin(1) * t453) * t477 + 0.4e1 * t384 * t430 + ((-t369 * t476 + t469 * t477) * t404 + ((t459 + (-t396 + t468) * t404) * t477 + (-t371 * t404 - t384 * t406) * t476) * qJD(2)) * pkin(2)) * t379 - ((0.4e1 * t438 + t457) * t478 - 0.4e1 * t383 * t430 + ((-0.2e1 * t406 * t469 + (-t436 / 0.4e1 - t390 * t385) * t404) * t478 + (t368 * t463 + (-t397 * t401 * t406 + (t370 * t404 + t383 * t406) * t398) * qJD(2)) * pkin(1)) * pkin(2)) * t470) * t409) * t450) * t485;
	t446 = 0.1e1 / t374 * ((t368 * t383 + t369 * t384) * t487 - 0.4e1 * t458 * t399 * t441) * t410 * t485;
	t435 = t405 * t448;
	t434 = t407 * t448;
	t433 = -t446 / 0.2e1;
	t432 = t446 / 0.2e1;
	t429 = t384 * t432;
	t428 = t405 * t433;
	t427 = t407 * t432;
	t422 = t383 * t428 + (t368 * t465 - t423 * t383) * t372;
	t421 = t384 * t427 + (-t369 * t464 + t424 * t384) * t372;
	t363 = t450 * t491;
	t1 = [t422 * t409, (t421 * t363 - t384 * t434) * t409; (t383 * t427 + t489) * t409, (-t384 * t435 + (t405 * t429 + t490) * t363) * t409; 0, (-t383 * t448 + (t383 * t432 + (-t368 * t397 + t383 * t426) * t372) * t363) * t409; (t384 * t428 - t490) * t409, (t383 * t434 + (t407 * t383 * t433 - t489) * t363) * t409; t421 * t409, (t422 * t363 + t383 * t435) * t409; 0, (-t384 * t448 + (t429 + (-t369 * t397 + t384 * t426) * t372) * t363) * t409; -t455, 0; t454, 0; 0, 0;];
	JRD_rot = t1;
end