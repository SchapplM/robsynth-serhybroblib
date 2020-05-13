% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% JRD_rot [9x1]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = fourbar1TE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_jacobiRD_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_jacobiRD_rot_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1TE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_jacobiRD_rot_sym_varpar: pkin has to be [4x1] (double)');
JRD_rot=NaN(9,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30; -t31; 0; t31; -t30; 0; 0; 0; 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:44
	% EndTime: 2020-04-24 19:49:44
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (1156->50), mult. (1712->92), div. (40->5), fcn. (434->4), ass. (0->57)
	t332 = sin(qJ(1));
	t337 = pkin(1) ^ 2;
	t333 = cos(qJ(1));
	t371 = pkin(2) * t333;
	t377 = -2 * pkin(1);
	t362 = t371 * t377 + t337;
	t374 = -pkin(3) - pkin(4);
	t322 = (pkin(2) - t374) * (pkin(2) + t374) + t362;
	t373 = -pkin(3) + pkin(4);
	t323 = (pkin(2) - t373) * (pkin(2) + t373) + t362;
	t375 = pkin(1) * pkin(2);
	t353 = (-t322 - t323) * t375;
	t313 = t332 * t353;
	t312 = qJD(1) * t313;
	t336 = pkin(2) ^ 2;
	t328 = t336 + t362;
	t324 = pkin(3) ^ 2 - pkin(4) ^ 2 + t328;
	t370 = qJD(1) * pkin(2);
	t359 = t324 * t370;
	t320 = t333 * t359;
	t372 = pkin(2) * t332;
	t321 = t324 * t372;
	t329 = -pkin(1) + t371;
	t368 = t322 * t323;
	t338 = sqrt(-t368);
	t366 = t332 * t338;
	t315 = pkin(2) * t366;
	t364 = t336 * t332 ^ 2;
	t345 = 0.2e1 * pkin(1) * t364 - t315;
	t335 = 0.1e1 / pkin(3);
	t381 = 0.1e1 / t328 ^ 2 * t335;
	t355 = t375 * t381;
	t351 = qJD(1) * t355;
	t347 = t332 * t351;
	t352 = t332 * t355;
	t317 = 0.1e1 / t338;
	t369 = t317 * t329;
	t325 = 0.1e1 / t328;
	t354 = 0.4e1 * qJD(1) * t364;
	t380 = t325 * t337 * t354 * t381 - t333 * t351;
	t383 = -(t345 * qJD(1) + t312 * t369 + t320) * t352 + t380 * (t329 * t338 + t321) - (t313 * t369 + t324 * t371 + t345) * t347;
	t382 = t317 * t313;
	t376 = 0.2e1 * t317;
	t357 = t312 / 0.2e1;
	t311 = (t333 * t353 - 0.4e1 * t337 * t364) * qJD(1);
	t367 = t329 * t311;
	t365 = t333 * t338;
	t363 = -t332 * t359 - t365 * t370;
	t361 = t329 * t377;
	t360 = pkin(1) * qJD(1) * t332;
	t358 = t325 * t335 / 0.2e1;
	t349 = 0.2e1 * t357;
	t348 = 0.1e1 / t368 * t312 * t382;
	t342 = 0.6e1 * t336 * t333 * t360 + t329 * t348 + t363;
	t307 = t317 * t312 * t372;
	t306 = (pkin(1) * t354 + t320 + (t332 * t348 + (t333 * t361 - t366) * qJD(1) + (t332 * t311 / 0.2e1 + t349 * t333) * t376) * pkin(2)) * t358 - (t321 + (t365 + (t361 + t382) * t332) * pkin(2)) * t347 - (-0.2e1 * pkin(2) * t329 * t360 + t307 - t363) * t352 + t380 * (-t324 * t329 + t315);
	t1 = [t306; (t307 + (t357 * t372 - t367 / 0.2e1) * t376 - t342) * t358 - t383; 0; ((t367 / 0.2e1 - t349 * t372) * t376 + t342) * t358 + t383; t306; 0; 0; 0; 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:44
	% EndTime: 2020-04-24 19:49:44
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (1156->54), mult. (1712->97), div. (40->5), fcn. (434->4), ass. (0->56)
	t323 = sin(qJ(1));
	t328 = pkin(1) ^ 2;
	t324 = cos(qJ(1));
	t362 = t324 * pkin(1);
	t354 = -0.2e1 * pkin(2) * t362 + t328;
	t366 = -pkin(3) - pkin(4);
	t313 = (pkin(2) - t366) * (pkin(2) + t366) + t354;
	t365 = -pkin(3) + pkin(4);
	t314 = (pkin(2) - t365) * (pkin(2) + t365) + t354;
	t360 = t313 * t314;
	t329 = sqrt(-t360);
	t356 = t329 * t323;
	t327 = pkin(2) ^ 2;
	t319 = t327 + t354;
	t315 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t319;
	t359 = t315 * t324;
	t374 = (t356 + t359) * pkin(2);
	t343 = pkin(1) * pkin(2) * (-t313 - t314);
	t307 = t323 * t343;
	t311 = 0.1e1 / t329;
	t373 = t311 * t307;
	t326 = 0.1e1 / pkin(4);
	t372 = 0.1e1 / t319 ^ 2 * t326;
	t357 = t327 * t323 ^ 2;
	t305 = (t324 * t343 - 0.4e1 * t328 * t357) * qJD(1);
	t363 = pkin(2) * t324;
	t320 = -pkin(1) + t363;
	t306 = qJD(1) * t307;
	t340 = 0.1e1 / t360 * t306 * t373;
	t371 = (t311 * t305 + t340) * t320;
	t316 = 0.1e1 / t319;
	t369 = pkin(1) * t372;
	t344 = qJD(1) * t369;
	t348 = qJD(1) * t357;
	t370 = 0.4e1 * t316 * t328 * t348 * t372 - t344 * t363;
	t309 = pkin(2) * t356;
	t341 = t323 * t344;
	t351 = pkin(1) * t357;
	t346 = 0.2e1 * t351;
	t361 = t311 * t320;
	t364 = pkin(2) * t323;
	t368 = -t370 * (t315 * t364 - t320 * t329) + (pkin(2) * t359 - t307 * t361 + t309 + t346) * pkin(2) * t341;
	t367 = 0.2e1 * t311;
	t358 = t324 * t329;
	t308 = pkin(2) * qJD(1) * t358;
	t355 = t311 * t306 * t364 + t308;
	t353 = qJD(1) * t323;
	t352 = 0.2e1 * t320 * pkin(1);
	t350 = t327 * t362;
	t349 = t306 * t361;
	t347 = t316 * t326 / 0.2e1;
	t345 = -t315 + t373;
	t342 = t364 * t369;
	t336 = t345 + t352;
	t300 = (-0.4e1 * pkin(1) * t348 + (t323 * t340 + (t324 * t306 / 0.2e1 + t323 * t305 / 0.2e1) * t367 + (t336 * t324 - t356) * qJD(1)) * pkin(2)) * t347 - (t336 * t323 + t358) * t327 * t341 - ((-t315 + t352) * pkin(2) * t353 + t355) * t342 + t370 * (t320 * t315 + t309);
	t1 = [t300; (-t371 + (t345 * pkin(2) + 0.6e1 * t350) * t353 + t355) * t347 - (-t349 + (t346 + t374) * qJD(1)) * t342 - t368; 0; (-t308 + t371 + (-0.6e1 * qJD(1) * t350 + (qJD(1) * t315 - t306 * t367) * pkin(2)) * t323) * t347 - (t349 + (-0.2e1 * t351 - t374) * qJD(1)) * t342 + t368; t300; 0; 0; 0; 0;];
	JRD_rot = t1;
end