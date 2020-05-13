% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbar1turnDE2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_jacobiRD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_jacobiRD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
JRD_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
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
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
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
	% StartTime: 2020-04-12 19:35:45
	% EndTime: 2020-04-12 19:35:46
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (19775->80), mult. (27380->172), div. (850->9), fcn. (7464->9), ass. (0->88)
	t360 = sin(qJ(2));
	t367 = pkin(1) ^ 2;
	t362 = cos(qJ(2));
	t418 = pkin(1) * t362;
	t398 = -0.2e1 * pkin(2) * t418 + t367;
	t420 = -pkin(3) - pkin(4);
	t350 = (pkin(2) - t420) * (pkin(2) + t420) + t398;
	t419 = -pkin(3) + pkin(4);
	t351 = (pkin(2) - t419) * (pkin(2) + t419) + t398;
	t421 = pkin(1) * pkin(2);
	t384 = (-t350 - t351) * t421;
	t343 = t360 * t384;
	t342 = qJD(2) * t343;
	t403 = t360 ^ 2 * t367;
	t388 = qJD(2) * t403;
	t385 = pkin(2) * t388;
	t405 = t350 * t351;
	t368 = sqrt(-t405);
	t402 = t360 * t368;
	t392 = pkin(1) * t402;
	t366 = pkin(2) ^ 2;
	t356 = t366 + t398;
	t352 = pkin(3) ^ 2 - pkin(4) ^ 2 + t356;
	t401 = t362 * t352;
	t399 = (pkin(1) * t401 + t392) * qJD(2);
	t346 = 0.1e1 / t368;
	t357 = -pkin(2) + t418;
	t406 = t346 * t357;
	t326 = -t342 * t406 + 0.2e1 * t385 + t399;
	t339 = -t357 * t352 - t392;
	t336 = 0.1e1 / t339 ^ 2;
	t349 = pkin(1) * t360 * t352;
	t340 = -t357 * t368 + t349;
	t338 = t340 ^ 2;
	t334 = t338 * t336 + 0.1e1;
	t409 = t336 * t340;
	t422 = -0.2e1 * t357;
	t395 = pkin(2) * t422;
	t400 = t368 * t362;
	t408 = t346 * t342;
	t325 = (-t360 * t408 + (-t400 + (t352 + t395) * t360) * qJD(2)) * pkin(1);
	t335 = 0.1e1 / t339;
	t414 = t325 * t335 * t336;
	t424 = (t326 * t409 - t338 * t414) / t334 ^ 2;
	t407 = t346 * t343;
	t353 = 0.1e1 / t356;
	t354 = 0.1e1 / t356 ^ 2;
	t423 = -0.2e1 * t354;
	t417 = pkin(3) * t356;
	t393 = t360 * t423;
	t383 = t393 * t421;
	t380 = qJD(2) * t383;
	t365 = 0.1e1 / pkin(3);
	t332 = 0.1e1 / t334;
	t391 = t332 * t417;
	t386 = t365 * t391;
	t320 = qJD(2) + ((t326 * t353 + t340 * t380) * t335 - (t325 * t353 + t339 * t380) * t409) * t386;
	t327 = t349 + (-t400 + (t395 - t407) * t360) * pkin(1);
	t323 = (t327 * t353 + t339 * t383) * t365;
	t389 = t323 * t409;
	t328 = -t343 * t406 + 0.2e1 * pkin(2) * t403 + (t401 + t402) * pkin(1);
	t324 = (t328 * t353 + t340 * t383) * t365;
	t415 = t324 * t335;
	t321 = 0.1e1 + (-t389 + t415) * t391;
	t416 = t320 * t321;
	t404 = t353 * t365;
	t331 = qJ(2) + atan2(t340 * t404, t339 * t404);
	t329 = sin(t331);
	t361 = sin(qJ(1));
	t413 = t329 * t361;
	t363 = cos(qJ(1));
	t412 = t329 * t363;
	t330 = cos(t331);
	t411 = t330 * t361;
	t410 = t330 * t363;
	t397 = qJD(1) * t361;
	t396 = qJD(1) * t363;
	t394 = 0.2e1 * t354;
	t390 = pkin(2) * qJD(2) * t360;
	t381 = 0.1e1 / t405 * t342 * t407;
	t378 = 0.8e1 * t353 * t354 * t366 * t388;
	t377 = -t320 * t413 + t330 * t396;
	t376 = t320 * t412 + t330 * t397;
	t375 = t320 * t411 + t329 * t396;
	t374 = t320 * t410 - t329 * t397;
	t341 = (t362 * t384 - 0.4e1 * t366 * t403) * qJD(2);
	t319 = ((0.6e1 * t367 * t362 * t390 - t341 * t406 - t357 * t381) * t353 + t340 * t378 + ((t326 * pkin(2) * t423 + t353 * t408) * t360 + ((t400 + (-t352 + t407) * t360) * t353 + (-t328 * t360 - t340 * t362) * pkin(2) * t394) * qJD(2)) * pkin(1)) * t335 * t386 + (0.2e1 * t415 - 0.2e1 * t389) * pkin(1) * pkin(3) * t332 * t390 + (-t323 * t326 - t324 * t325 - ((0.4e1 * t385 + t399) * t353 + t339 * t378 + ((-0.2e1 * t362 * t408 + (-t341 * t346 - t381) * t360) * t353 + (t325 * t393 + (t353 * t362 * t422 + (-t327 * t360 - t339 * t362) * t394) * qJD(2)) * pkin(2)) * pkin(1)) * t340 * t365) * t336 * t391 + (-0.2e1 * t415 * t424 + 0.2e1 * (t332 * t414 + t336 * t424) * t323 * t340) * t417;
	t1 = [t377, t319 * t412 + t374 * t321; t376, t319 * t413 + t375 * t321; 0, -t319 * t330 + t329 * t416; -t375, t319 * t410 - t376 * t321; t374, t319 * t411 + t377 * t321; 0, t319 * t329 + t330 * t416; -t397, 0; t396, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:46
	% EndTime: 2020-04-12 19:35:47
	% DurationCPUTime: 1.78s
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