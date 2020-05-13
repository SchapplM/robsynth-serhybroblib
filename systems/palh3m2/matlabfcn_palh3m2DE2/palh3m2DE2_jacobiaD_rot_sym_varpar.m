% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m2DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m2DE2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobiaD_rot_sym_varpar: pkin has to be [18x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:25
	% EndTime: 2020-05-07 02:13:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:32
	% EndTime: 2020-05-07 02:13:44
	% DurationCPUTime: 10.88s
	% Computational Cost: add. (154841->144), mult. (275138->292), div. (3391->16), fcn. (407303->20), ass. (0->147)
	t441 = pkin(17) + pkin(18);
	t439 = sin(t441);
	t440 = cos(t441);
	t444 = sin(pkin(16));
	t445 = cos(pkin(16));
	t450 = sin(pkin(15));
	t455 = cos(pkin(15));
	t437 = t444 * t455 + t445 * t450;
	t438 = -t444 * t450 + t445 * t455;
	t447 = sin(qJ(3));
	t452 = cos(qJ(3));
	t431 = t447 * t437 - t438 * t452;
	t448 = sin(qJ(2));
	t453 = cos(qJ(2));
	t479 = t437 * t452 + t447 * t438;
	t569 = t431 * t448 - t479 * t453;
	t570 = t431 * t453 + t448 * t479;
	t573 = t570 * t439 + t440 * t569;
	t578 = -t569 * t439 + t440 * t570;
	t403 = qJ(2) + qJ(3) + atan2(t573, -t578);
	t402 = cos(t403);
	t401 = sin(t403);
	t449 = sin(qJ(1));
	t513 = t449 * t401;
	t385 = atan2(t513, t402);
	t382 = cos(t385);
	t381 = sin(t385);
	t497 = t381 * t513;
	t373 = t382 * t402 + t497;
	t371 = 0.1e1 / t373 ^ 2;
	t409 = t573 ^ 2;
	t411 = 0.1e1 / t578 ^ 2;
	t408 = t409 * t411 + 0.1e1;
	t406 = 0.1e1 / t408;
	t410 = 0.1e1 / t578;
	t527 = t411 * t573;
	t429 = t479 * qJD(3);
	t568 = t431 * qJD(3);
	t467 = qJD(2) * t570 + t448 * t429 + t453 * t568;
	t577 = -qJD(2) * t569 + t429 * t453 - t448 * t568;
	t574 = t439 * t577 + t440 * t467;
	t583 = -t467 * t439 + t577 * t440;
	t374 = qJD(2) + qJD(3) + (-t410 * t574 + t527 * t583) * t406;
	t454 = cos(qJ(1));
	t442 = t449 ^ 2;
	t397 = t401 ^ 2;
	t399 = 0.1e1 / t402 ^ 2;
	t536 = t397 * t399;
	t386 = t442 * t536 + 0.1e1;
	t383 = 0.1e1 / t386;
	t398 = 0.1e1 / t402;
	t507 = qJD(1) * t454;
	t491 = t401 * t507;
	t531 = t402 * t449;
	t534 = t397 * t449;
	t360 = ((t374 * t531 + t491) * t398 + t374 * t399 * t534) * t383;
	t541 = t382 * t401;
	t544 = t374 * t449;
	t354 = (t360 * t449 - t374) * t541 + (t491 + (-t360 + t544) * t402) * t381;
	t370 = 0.1e1 / t373;
	t506 = -0.2e1 * t354 * t370 * t371;
	t508 = qJD(1) * t449;
	t580 = t410 * t578 + t527 * t573;
	t584 = -t406 * t580 + 0.1e1;
	t490 = 0.1e1 + t536;
	t563 = t449 * t490;
	t477 = t383 * t563;
	t588 = t584 * t477;
	t599 = (-t381 * t402 + t382 * t513) * t588 + (t381 * t531 - t541) * t584;
	t600 = t370 * t584;
	t603 = (-t371 * t508 + t454 * t506) * t599 + t454 * t374 * t600;
	t602 = (t354 * t584 + t374 * t599) * t371;
	t547 = t371 * t401;
	t443 = t454 ^ 2;
	t535 = t397 * t443;
	t369 = t371 * t535 + 0.1e1;
	t483 = t507 * t534;
	t545 = t374 * t402;
	t554 = (t506 * t535 + 0.2e1 * (t401 * t443 * t545 - t483) * t371) / t369 ^ 2;
	t601 = (t402 * t600 - t547 * t599) * t454 * t554;
	t542 = t381 * t449;
	t598 = (-t360 * t542 + t374 * t381 + t382 * t507) * t588 + (t360 * t381 - t374 * t542) * t584;
	t597 = ((t449 * t588 - t584) * t374 + (t449 * t584 - t588) * t360) * t382;
	t446 = sin(qJ(4));
	t485 = qJD(1) * t402 + qJD(4);
	t532 = t401 * t454;
	t469 = t374 * t532 + t449 * t485;
	t486 = qJD(4) * t402 + qJD(1);
	t451 = cos(qJ(4));
	t509 = t454 * t451;
	t365 = t446 * t469 - t486 * t509;
	t510 = t454 * t446;
	t366 = t451 * t469 + t486 * t510;
	t511 = t449 * t451;
	t393 = t402 * t510 + t511;
	t387 = t393 ^ 2;
	t512 = t449 * t446;
	t395 = -t402 * t509 + t512;
	t389 = 0.1e1 / t395 ^ 2;
	t380 = t387 * t389 + 0.1e1;
	t378 = 0.1e1 / t380;
	t388 = 0.1e1 / t395;
	t537 = t393 * t451;
	t539 = t388 * t446;
	t474 = t537 * t389 + t539;
	t471 = t474 * t454;
	t549 = t366 * t388 * t389;
	t505 = -0.2e1 * t549;
	t489 = t393 * t505;
	t538 = t389 * t393;
	t553 = 0.2e1 * (-t365 * t538 - t387 * t549) / t380 ^ 2;
	t594 = (-t474 * t532 * t553 + (t471 * t545 + (-t474 * t508 + ((qJD(4) * t388 + t489) * t451 + (-t365 * t451 + (-qJD(4) * t393 - t366) * t446) * t389) * t454) * t401) * t378) * t584;
	t533 = t398 * t401;
	t567 = t401 * t536;
	t470 = 0.2e1 * t374 * (t398 * t567 + t533);
	t482 = t490 * t454;
	t550 = (0.2e1 * t399 * t483 + t442 * t470) / t386 ^ 2;
	t593 = (-t550 * t563 + (qJD(1) * t482 + t449 * t470) * t383) * t584;
	t591 = t507 * t584;
	t587 = t508 * t600;
	t566 = t583 * t411;
	t530 = t410 * t566;
	t579 = t574 * t527;
	t585 = 0.2e1 * t580 * (-t409 * t530 + t579) / t408 ^ 2;
	t582 = -0.2e1 * t573;
	t581 = t583 * t410;
	t555 = t383 - 0.1e1;
	t551 = t360 * t401;
	t540 = t382 * t449;
	t504 = t401 * t554;
	t501 = t401 * t550;
	t500 = t371 * t532;
	t499 = t371 * t402 * t454;
	t496 = t398 * t534;
	t488 = t530 * t582;
	t484 = t383 * t496;
	t478 = t485 * t454;
	t468 = t401 * t378 * t471;
	t392 = t402 * t511 + t510;
	t391 = t402 * t512 - t509;
	t367 = 0.1e1 / t369;
	t364 = (-t381 * t401 * t555 + t382 * t484) * t454;
	t358 = t585 + ((-t411 * t574 - t488) * t573 - t581 + t578 * t566 - t579) * t406;
	t357 = t585 + (-t581 - t573 * t488 + (t574 * t582 + t578 * t583) * t411) * t406;
	t352 = t357 * t477 + t593;
	t351 = t358 * t477 + t593;
	t1 = [-t454 * t398 * t501 + (t374 * t482 - t508 * t533) * t383, t352, t351, 0; (-t370 * t504 + (t370 * t545 + (-qJD(1) * t364 - t354) * t547) * t367) * t449 + (-t371 * t504 * t364 + (((-t360 * t484 - t545 * t555 + t501) * t381 + (-t496 * t550 + t551 + (-t551 + (0.2e1 * t401 + t567) * t544) * t383) * t382) * t500 + (t371 * t545 + t401 * t506) * t364 + (t370 + ((-t442 + t443) * t398 * t397 * t383 * t382 + t555 * t497) * t371) * t401 * qJD(1)) * t367) * t454, t601 + ((t587 + (-t357 * t370 + t602) * t454) * t402 + ((t352 * t540 - t357 * t382 + t598) * t500 + ((t357 * t449 - t352 + t591) * t381 + t597) * t499 + t603) * t401) * t367, t601 + ((t587 + (-t358 * t370 + t602) * t454) * t402 + ((t351 * t540 - t358 * t382 + t598) * t500 + ((t358 * t449 - t351 + t591) * t381 + t597) * t499 + t603) * t401) * t367, 0; (-t388 * t391 - t392 * t538) * t553 + (t392 * t489 + t486 * t388 * t511 + (-t374 * t513 + t478) * t539 + (-t392 * t365 - t391 * t366 + t478 * t537 + (-t374 * t401 * t451 - t446 * t486) * t393 * t449) * t389) * t378, t357 * t468 + t594, t358 * t468 + t594, -t553 + (-0.2e1 * t365 * t389 * t378 + (t378 * t505 - t389 * t553) * t393) * t393;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:53
	% EndTime: 2020-05-07 02:13:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:54
	% EndTime: 2020-05-07 02:13:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end