% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbar1turnTE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_jacobiRD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_jacobiRD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnTE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_jacobiRD_rot_sym_varpar: pkin has to be [5x1] (double)');
JRD_rot=NaN(9,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:26
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
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:26
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
	% StartTime: 2020-04-12 19:20:28
	% EndTime: 2020-04-12 19:20:29
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (6747->112), mult. (10078->208), div. (316->5), fcn. (2774->6), ass. (0->99)
	t428 = pkin(2) ^ 2;
	t429 = pkin(1) ^ 2;
	t424 = cos(qJ(2));
	t510 = -2 * pkin(2);
	t478 = t424 * t510;
	t484 = pkin(1) * t478 + t429;
	t416 = t428 + t484;
	t413 = 0.1e1 / t416;
	t422 = sin(qJ(2));
	t412 = pkin(3) ^ 2 - pkin(4) ^ 2 + t416;
	t409 = pkin(1) * t422 * t412;
	t417 = pkin(1) * t424 - pkin(2);
	t479 = t417 * t510;
	t507 = -pkin(3) - pkin(4);
	t410 = (pkin(2) - t507) * (pkin(2) + t507) + t484;
	t506 = -pkin(3) + pkin(4);
	t411 = (pkin(2) - t506) * (pkin(2) + t506) + t484;
	t498 = t410 * t411;
	t430 = sqrt(-t498);
	t487 = t430 * t424;
	t508 = pkin(1) * pkin(2);
	t465 = (-t410 - t411) * t508;
	t403 = t422 * t465;
	t406 = 0.1e1 / t430;
	t500 = t406 * t403;
	t389 = t409 + (-t487 + (t479 - t500) * t422) * pkin(1);
	t489 = t424 * t412;
	t494 = t422 * t430;
	t499 = t406 * t417;
	t502 = pkin(2) * t429;
	t420 = t422 ^ 2;
	t509 = 0.2e1 * t420;
	t390 = -t403 * t499 + t502 * t509 + (t489 + t494) * pkin(1);
	t476 = pkin(1) * t494;
	t399 = -t417 * t412 - t476;
	t400 = -t417 * t430 + t409;
	t419 = t422 * t420;
	t414 = 0.1e1 / t416 ^ 2;
	t488 = t428 * t429;
	t462 = 0.4e1 * t413 * t414 * t488;
	t477 = t414 * t508;
	t483 = -t424 ^ 2 + t509;
	t490 = t424 * t400;
	t495 = t422 * t424;
	t437 = (t483 * t400 - (t390 + 0.3e1 * t399) * t495 - t389 * t420) * t477 + (t399 * t419 + t420 * t490) * t462;
	t467 = -t400 / 0.2e1 + t389 / 0.2e1;
	t468 = t399 / 0.2e1 + t390 / 0.2e1;
	t438 = (t468 * t422 - t467 * t424) * t413;
	t402 = qJD(2) * t403;
	t474 = t406 * t422 * t402;
	t387 = (-t474 + (-t487 + (t412 + t479) * t422) * qJD(2)) * pkin(1);
	t475 = qJD(2) * t502;
	t466 = t420 * t475;
	t485 = (pkin(1) * t489 + t476) * qJD(2);
	t388 = -t402 * t499 + 0.2e1 * t466 + t485;
	t492 = t424 * t388;
	t444 = (t387 * t420 + t422 * t492) * t477;
	t401 = (-0.4e1 * t420 * t488 + t424 * t465) * qJD(2);
	t461 = 0.1e1 / t498 * t402 * t500;
	t505 = -t422 / 0.2e1;
	t381 = 0.4e1 * t466 + (-t422 * t461 + t417 * qJD(2) * t478 + 0.2e1 * (t401 * t505 - t424 * t402) * t406) * pkin(1) + t485;
	t469 = -t388 / 0.2e1 + t381 / 0.2e1;
	t382 = 0.6e1 * t475 * t495 + (-t406 * t401 - t461) * t417 + (t474 + (t487 + (-t412 + t500) * t422) * qJD(2)) * pkin(1);
	t470 = -t387 / 0.2e1 - t382 / 0.2e1;
	t523 = (t437 - t438) * qJD(2) + (t469 * t422 - t470 * t424) * t413 - t444;
	t521 = (t467 * t422 + t468 * t424) * t413;
	t425 = cos(qJ(1));
	t427 = 0.1e1 / pkin(3);
	t497 = t413 * t427;
	t472 = t425 * t497;
	t463 = t472 / 0.2e1;
	t459 = t424 * t463;
	t471 = t400 * t505;
	t520 = t399 * t459 + t471 * t472;
	t423 = sin(qJ(1));
	t473 = t423 * t497;
	t464 = -t473 / 0.2e1;
	t460 = t422 * t464;
	t519 = t399 * t460 + t464 * t490;
	t491 = t424 * t399;
	t517 = (-t400 * t420 + t422 * t491) * t477;
	t516 = (t399 * t420 + t422 * t490) * t477;
	t504 = t422 / 0.2e1;
	t503 = t424 / 0.2e1;
	t496 = t422 * t387;
	t493 = t424 * t387;
	t482 = qJD(1) * t423;
	t481 = qJD(1) * t425;
	t480 = qJD(2) * t427;
	t458 = t473 * t493 / 0.2e1 + t388 * t460 + t519 * qJD(2) + t520 * qJD(1);
	t457 = t519 * qJD(1) + t520 * qJD(2) + t388 * t459 + t463 * t496;
	t448 = t491 / 0.2e1 + t471;
	t447 = -t490 / 0.2e1 + t399 * t505;
	t446 = (t389 * t504 + t390 * t503) * t413;
	t445 = (t389 * t503 + t390 * t505) * t413;
	t443 = (-t388 * t420 + t422 * t493) * t477;
	t436 = (-t390 * t420 + (t389 - 0.3e1 * t400) * t495 - t483 * t399) * t477 + (t400 * t419 - t420 * t491) * t462;
	t435 = t443 + (-t470 * t422 - t469 * t424) * t413 + (t436 + t521) * qJD(2);
	t1 = [-t423 * t480 * t517 + t458, ((t445 - t517) * t482 + ((-t424 * t381 / 0.2e1 + t382 * t504) * t413 + t443 + (t446 + t436) * qJD(2)) * t425) * t427 + t457; (t425 * qJD(2) * t517 + (t448 * t482 + (-t493 / 0.2e1 + t388 * t504 - t447 * qJD(2)) * t425) * t413) * t427, ((t517 + t438) * t481 + t435 * t423) * t427; 0, -t523 * t427; (t423 * qJD(2) * t516 + (t447 * t481 + (-t492 / 0.2e1 - t496 / 0.2e1 - t448 * qJD(2)) * t423) * t413) * t427, ((t516 - t521) * t482 + t523 * t425) * t427; -t425 * t480 * t516 + t457, ((t446 - t516) * t481 + ((t381 * t504 + t382 * t503) * t413 - t444 + (t445 + t437) * qJD(2)) * t423) * t427 + t458; 0, t435 * t427; -t482, 0; t481, 0; 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:27
	% EndTime: 2020-04-12 19:20:27
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (2395->62), mult. (3562->124), div. (100->5), fcn. (945->6), ass. (0->68)
	t322 = sin(qJ(2));
	t329 = pkin(1) ^ 2;
	t324 = cos(qJ(2));
	t369 = pkin(2) * t324;
	t360 = -0.2e1 * pkin(1) * t369 + t329;
	t374 = -pkin(3) - pkin(4);
	t312 = (pkin(2) - t374) * (pkin(2) + t374) + t360;
	t373 = -pkin(3) + pkin(4);
	t313 = (pkin(2) - t373) * (pkin(2) + t373) + t360;
	t375 = pkin(1) * pkin(2);
	t345 = (-t312 - t313) * t375;
	t305 = t322 * t345;
	t304 = qJD(2) * t305;
	t321 = t322 ^ 2;
	t328 = pkin(2) ^ 2;
	t353 = pkin(1) * qJD(2) * t328;
	t346 = t321 * t353;
	t365 = t312 * t313;
	t330 = sqrt(-t365);
	t363 = t322 * t330;
	t354 = pkin(2) * t363;
	t357 = qJD(2) * t324;
	t318 = t328 + t360;
	t314 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t318;
	t370 = pkin(2) * t314;
	t361 = qJD(2) * t354 + t357 * t370;
	t308 = 0.1e1 / t330;
	t319 = pkin(1) - t369;
	t366 = t308 * t319;
	t298 = t304 * t366 + 0.2e1 * t346 + t361;
	t364 = t321 * t328;
	t376 = 0.2e1 * pkin(1);
	t300 = t305 * t366 + t364 * t376 + (t324 * t314 + t363) * pkin(2);
	t311 = t322 * t370;
	t302 = t319 * t330 + t311;
	t315 = 0.1e1 / t318;
	t316 = 0.1e1 / t318 ^ 2;
	t347 = -0.4e1 * t329 * t364;
	t341 = t315 * t316 * qJD(2) * t347;
	t379 = t308 * t304;
	t352 = t322 * t379;
	t355 = t316 * t375;
	t362 = t324 * t330;
	t367 = t308 * t305;
	t372 = -t315 / 0.2e1;
	t377 = 0.1e1 / t365 * t304 * t367 + t308 * (t324 * t345 + t347) * qJD(2);
	t334 = (t298 * t322 + (t300 * t322 + t302 * t324) * qJD(2)) * t355 + (0.6e1 * t324 * t322 * t353 + t377 * t319 + (t352 + (t362 + (-t314 + t367) * t322) * qJD(2)) * pkin(2)) * t372 + t302 * t341;
	t371 = t315 / 0.2e1;
	t323 = sin(qJ(1));
	t359 = qJD(1) * t323;
	t325 = cos(qJ(1));
	t358 = qJD(1) * t325;
	t356 = t319 * t376;
	t351 = t323 * t371;
	t350 = t325 * t372;
	t349 = qJD(1) * t371;
	t348 = t322 * t355;
	t343 = qJD(2) * t348;
	t299 = t311 + (-t362 + (t356 - t367) * t322) * pkin(2);
	t301 = t319 * t314 - t354;
	t340 = t299 * t371 - t301 * t348;
	t339 = t300 * t371 - t302 * t348;
	t337 = t323 * t349 + t325 * t343;
	t336 = -t323 * t343 + t325 * t349;
	t297 = (-t352 + (-t362 + (t314 + t356) * t322) * qJD(2)) * pkin(2);
	t335 = (0.4e1 * t346 + (-t377 * t322 - 0.2e1 * t324 * t379 + t356 * t357) * pkin(2) + t361) * t372 + t301 * t341 + (t297 * t322 + (t299 * t322 + t301 * t324) * qJD(2)) * t355;
	t327 = 0.1e1 / pkin(4);
	t1 = [(t297 * t351 + t336 * t301) * t327, (t335 * t325 + t340 * t359) * t327; (t297 * t350 + t337 * t301) * t327, (t335 * t323 - t340 * t358) * t327; 0, -t334 * t327; (t298 * t351 + t336 * t302) * t327, (t334 * t325 + t339 * t359) * t327; (t298 * t350 + t337 * t302) * t327, (t334 * t323 - t339 * t358) * t327; 0, t335 * t327; -t359, 0; t358, 0; 0, 0;];
	JRD_rot = t1;
end