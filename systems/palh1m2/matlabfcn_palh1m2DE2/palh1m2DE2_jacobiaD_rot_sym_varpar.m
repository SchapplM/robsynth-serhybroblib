% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m2DE2
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
%   Wie in palh1m2DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh1m2DE2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2DE2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_jacobiaD_rot_sym_varpar: pkin has to be [22x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:31
	% EndTime: 2020-05-02 21:08:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:34
	% EndTime: 2020-05-02 21:08:46
	% DurationCPUTime: 11.65s
	% Computational Cost: add. (155435->145), mult. (275138->294), div. (3391->16), fcn. (407303->20), ass. (0->147)
	t429 = pkin(22) + pkin(21);
	t427 = sin(t429);
	t428 = cos(t429);
	t432 = cos(pkin(20));
	t437 = sin(pkin(18));
	t541 = sin(pkin(20));
	t542 = cos(pkin(18));
	t422 = t437 * t432 - t542 * t541;
	t423 = t542 * t432 + t437 * t541;
	t434 = sin(qJ(3));
	t439 = cos(qJ(3));
	t414 = t422 * t439 - t434 * t423;
	t435 = sin(qJ(2));
	t440 = cos(qJ(2));
	t550 = t434 * t422 + t423 * t439;
	t466 = t414 * t440 - t435 * t550;
	t568 = t435 * t414 + t440 * t550;
	t572 = t466 * t427 + t428 * t568;
	t389 = 0.1e1 / t572;
	t390 = 0.1e1 / t572 ^ 2;
	t546 = t550 * qJD(3);
	t562 = t414 * qJD(3);
	t571 = qJD(2) * t568 + t435 * t562 + t546 * t440;
	t579 = qJD(2) * t466 - t435 * t546 + t562 * t440;
	t596 = -t571 * t427 + t579 * t428;
	t555 = t596 * t390;
	t515 = t389 * t555;
	t573 = t568 * t427 - t428 * t466;
	t578 = -0.2e1 * t573;
	t601 = -t573 * t515 * t578 - t596 * t389;
	t382 = qJ(2) + qJ(3) + atan2(t573, -t572);
	t381 = cos(t382);
	t380 = sin(t382);
	t436 = sin(qJ(1));
	t500 = t436 * t380;
	t362 = atan2(-t500, -t381);
	t361 = cos(t362);
	t360 = sin(t362);
	t482 = t360 * t500;
	t352 = -t361 * t381 - t482;
	t350 = 0.1e1 / t352 ^ 2;
	t383 = t427 * t579 + t571 * t428;
	t388 = t573 ^ 2;
	t387 = t388 * t390 + 0.1e1;
	t385 = 0.1e1 / t387;
	t570 = t390 * t573;
	t353 = qJD(2) + qJD(3) + (-t383 * t389 + t570 * t596) * t385;
	t441 = cos(qJ(1));
	t430 = t436 ^ 2;
	t376 = t380 ^ 2;
	t378 = 0.1e1 / t381 ^ 2;
	t521 = t376 * t378;
	t365 = t430 * t521 + 0.1e1;
	t363 = 0.1e1 / t365;
	t377 = 0.1e1 / t381;
	t493 = qJD(1) * t441;
	t474 = t380 * t493;
	t516 = t381 * t436;
	t519 = t376 * t436;
	t339 = (-(-t353 * t516 - t474) * t377 + t353 * t378 * t519) * t363;
	t526 = t361 * t380;
	t530 = t353 * t436;
	t333 = (-t339 * t436 + t353) * t526 + (-t474 + (t339 - t530) * t381) * t360;
	t349 = 0.1e1 / t352;
	t539 = t333 * t349 * t350;
	t491 = 0.2e1 * t539;
	t494 = qJD(1) * t436;
	t575 = t389 * t572 + t570 * t573;
	t580 = -t385 * t575 + 0.1e1;
	t473 = 0.1e1 + t521;
	t554 = t436 * t473;
	t463 = t363 * t554;
	t584 = t580 * t463;
	t595 = (t360 * t381 - t361 * t500) * t584 + (-t360 * t516 + t526) * t580;
	t597 = t349 * t580;
	t600 = (t350 * t494 + t441 * t491) * t595 - t441 * t353 * t597;
	t599 = (-t333 * t580 - t353 * t595) * t350;
	t533 = t350 * t380;
	t431 = t441 ^ 2;
	t520 = t376 * t431;
	t348 = t350 * t520 + 0.1e1;
	t468 = t493 * t519;
	t531 = t353 * t381;
	t540 = 0.2e1 * (-t520 * t539 + (t380 * t431 * t531 - t468) * t350) / t348 ^ 2;
	t598 = (-t381 * t597 + t533 * t595) * t441 * t540;
	t527 = t360 * t436;
	t594 = (t339 * t527 - t353 * t360 - t361 * t493) * t584 + (-t339 * t360 + t353 * t527) * t580;
	t593 = ((-t436 * t584 + t580) * t353 + (-t436 * t580 + t584) * t339) * t361;
	t438 = cos(qJ(4));
	t496 = t441 * t438;
	t433 = sin(qJ(4));
	t499 = t436 * t433;
	t373 = t381 * t496 + t499;
	t458 = t381 * t499 + t496;
	t517 = t380 * t441;
	t483 = t353 * t517;
	t344 = qJD(1) * t458 - t373 * qJD(4) + t433 * t483;
	t470 = -qJD(1) * t381 + qJD(4);
	t471 = qJD(4) * t381 - qJD(1);
	t497 = t441 * t433;
	t345 = -t471 * t497 + (t436 * t470 - t483) * t438;
	t498 = t436 * t438;
	t372 = t381 * t497 - t498;
	t366 = t372 ^ 2;
	t368 = 0.1e1 / t373 ^ 2;
	t359 = t366 * t368 + 0.1e1;
	t357 = 0.1e1 / t359;
	t367 = 0.1e1 / t373;
	t522 = t372 * t438;
	t524 = t367 * t433;
	t459 = t368 * t522 - t524;
	t455 = t459 * t441;
	t535 = t345 * t367 * t368;
	t481 = t372 * t535;
	t523 = t368 * t372;
	t528 = 0.2e1 / t359 ^ 2 * (-t344 * t523 - t366 * t535);
	t590 = (-t459 * t517 * t528 + (t455 * t531 + (-t459 * t494 + ((-qJD(4) * t367 - 0.2e1 * t481) * t438 + (-t344 * t438 + (-qJD(4) * t372 + t345) * t433) * t368) * t441) * t380) * t357) * t580;
	t518 = t377 * t380;
	t556 = t380 * t521;
	t454 = 0.2e1 * t353 * (t377 * t556 + t518);
	t467 = t473 * t441;
	t536 = (0.2e1 * t378 * t468 + t430 * t454) / t365 ^ 2;
	t589 = (-t536 * t554 + (qJD(1) * t467 + t436 * t454) * t363) * t580;
	t587 = t493 * t580;
	t583 = t494 * t597;
	t566 = t383 * t570;
	t581 = 0.2e1 * t575 * (-t388 * t515 + t566) / t387 ^ 2;
	t544 = -0.2e1 * t357;
	t543 = t363 - 0.1e1;
	t537 = t339 * t380;
	t525 = t361 * t436;
	t490 = t380 * t540;
	t487 = t380 * t536;
	t486 = t350 * t517;
	t485 = t350 * t381 * t441;
	t480 = t377 * t519;
	t469 = t363 * t480;
	t464 = t470 * t441;
	t453 = t380 * t357 * t455;
	t371 = -t381 * t498 + t497;
	t346 = 0.1e1 / t348;
	t343 = (t543 * t380 * t360 - t361 * t469) * t441;
	t337 = t581 + (t572 * t555 - 0.2e1 * t566 + t601) * t385;
	t336 = t581 + ((t383 * t578 + t572 * t596) * t390 + t601) * t385;
	t331 = t336 * t463 + t589;
	t330 = t337 * t463 + t589;
	t1 = [-t441 * t377 * t487 + (t353 * t467 - t494 * t518) * t363, t331, t330, 0; (t349 * t490 + (-t349 * t531 + (qJD(1) * t343 + t333) * t533) * t346) * t436 + (t350 * t490 * t343 + (-((t339 * t469 + t543 * t531 - t487) * t360 + (t480 * t536 - t537 + (t537 + (-0.2e1 * t380 - t556) * t530) * t363) * t361) * t486 + (-t350 * t531 + t380 * t491) * t343 + (-t349 + ((-t430 + t431) * t377 * t376 * t363 * t361 + t543 * t482) * t350) * t380 * qJD(1)) * t346) * t441, t598 + ((-t583 + (t336 * t349 + t599) * t441) * t381 + (-(-t331 * t525 + t336 * t361 + t594) * t486 - ((-t336 * t436 + t331 - t587) * t360 + t593) * t485 + t600) * t380) * t346, t598 + ((-t583 + (t337 * t349 + t599) * t441) * t381 + (-(-t330 * t525 + t337 * t361 + t594) * t486 - ((-t337 * t436 + t330 - t587) * t360 + t593) * t485 + t600) * t380) * t346, 0; (t367 * t458 + t371 * t523) * t528 + (0.2e1 * t371 * t481 - t471 * t367 * t498 + (t353 * t500 + t464) * t524 + (t371 * t344 + t458 * t345 - t464 * t522 - (t353 * t380 * t438 + t433 * t471) * t372 * t436) * t368) * t357, t336 * t453 + t590, t337 * t453 + t590, -t528 + (t344 * t368 * t544 + (-t368 * t528 + t535 * t544) * t372) * t372;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:28
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:29
	% EndTime: 2020-05-02 21:08:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 21:08:30
	% EndTime: 2020-05-02 21:08:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end