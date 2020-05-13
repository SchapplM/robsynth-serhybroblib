% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m1DE2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m1DE2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_jacobiR_rot_sym_varpar: pkin has to be [23x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
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
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
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
	% StartTime: 2020-04-15 18:49:35
	% EndTime: 2020-04-15 18:49:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t21 = qJ(2) + qJ(3);
	t19 = sin(t21);
	t22 = sin(qJ(1));
	t27 = t22 * t19;
	t20 = cos(t21);
	t26 = t22 * t20;
	t23 = cos(qJ(1));
	t25 = t23 * t19;
	t24 = t23 * t20;
	t1 = [-t26, -t25, -t25, 0; t24, -t27, -t27, 0; 0, t20, t20, 0; t27, -t24, -t24, 0; -t25, -t26, -t26, 0; 0, -t19, -t19, 0; t23, 0, 0, 0; t22, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:18
	% EndTime: 2020-04-15 18:50:36
	% DurationCPUTime: 11.63s
	% Computational Cost: add. (322864->99), mult. (487440->175), div. (21464->11), fcn. (308998->19), ass. (0->111)
	t387 = -2 * pkin(1);
	t386 = -2 * pkin(5);
	t385 = pkin(4) * pkin(5);
	t384 = -pkin(8) - pkin(3);
	t383 = -pkin(8) + pkin(3);
	t382 = (-pkin(9) - pkin(11));
	t381 = (-pkin(9) + pkin(11));
	t322 = sin(pkin(21));
	t380 = t322 / 0.2e1;
	t327 = cos(qJ(3));
	t379 = -t327 / 0.2e1;
	t378 = cos(qJ(2));
	t377 = sin(pkin(19));
	t331 = pkin(7) ^ 2;
	t325 = sin(qJ(2));
	t329 = cos(pkin(19));
	t316 = t325 * t329 - t378 * t377;
	t374 = pkin(7) * t316;
	t356 = (pkin(1) ^ 2) + t374 * t387;
	t313 = t331 + t356;
	t310 = pkin(3) ^ 2 - pkin(8) ^ 2 + t313;
	t314 = pkin(1) - t374;
	t308 = (pkin(7) - t384) * (pkin(7) + t384) + t356;
	t309 = (pkin(7) - t383) * (pkin(7) + t383) + t356;
	t337 = sqrt(-t309 * t308);
	t317 = t325 * t377 + t378 * t329;
	t373 = pkin(7) * t317;
	t304 = t310 * t373 + t314 * t337;
	t360 = t327 * t304;
	t366 = t317 * t337;
	t303 = -pkin(7) * t366 + t314 * t310;
	t324 = sin(qJ(3));
	t365 = t324 * t303;
	t311 = 0.1e1 / t313;
	t334 = 0.1e1 / pkin(3);
	t367 = t311 * t334;
	t299 = (t365 / 0.2e1 + t360 / 0.2e1) * t367;
	t361 = t327 * t303;
	t364 = t324 * t304;
	t300 = (-t361 / 0.2e1 + t364 / 0.2e1) * t367;
	t321 = pkin(23) + pkin(22);
	t319 = sin(t321);
	t320 = cos(t321);
	t292 = t320 * t299 - t319 * t300;
	t376 = pkin(4) * t292;
	t343 = t319 * t299 + t320 * t300;
	t375 = pkin(4) * t343;
	t333 = pkin(4) ^ 2;
	t357 = (pkin(5) ^ 2) - t376 * t386;
	t289 = t333 + t357;
	t286 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t289;
	t290 = pkin(5) + t376;
	t284 = ((pkin(4) - t382) * (pkin(4) + t382)) + t357;
	t285 = ((pkin(4) - t381) * (pkin(4) + t381)) + t357;
	t336 = sqrt(-t285 * t284);
	t275 = t290 * t286 - t336 * t375;
	t276 = t286 * t375 + t290 * t336;
	t323 = cos(pkin(21));
	t287 = 0.1e1 / t289;
	t330 = 0.1e1 / pkin(11);
	t369 = t287 * t330;
	t273 = (-t275 * t323 / 0.2e1 + t276 * t380) * t369;
	t271 = 0.1e1 / t273 ^ 2;
	t272 = (t275 * t380 + t276 * t323 / 0.2e1) * t369;
	t372 = 0.1e1 / (t271 * t272 ^ 2 + 0.1e1) * t330;
	t371 = t271 * t272;
	t370 = t287 * t323;
	t355 = pkin(1) * t373;
	t368 = 0.2e1 / t337 * (t308 + t309) * t355;
	t266 = qJ(2) + qJ(3) + atan2(t272, t273);
	t264 = sin(t266);
	t326 = sin(qJ(1));
	t363 = t326 * t264;
	t265 = cos(t266);
	t362 = t326 * t265;
	t328 = cos(qJ(1));
	t359 = t328 * t264;
	t358 = t328 * t265;
	t354 = 0.1e1 / t289 ^ 2 * t385;
	t281 = 0.1e1 / t336;
	t353 = t281 * t290 / 0.2e1;
	t352 = -t281 * t343 / 0.2e1;
	t351 = t287 * t380;
	t350 = -t370 / 0.2e1;
	t349 = t370 / 0.2e1;
	t348 = t333 * t343 * t386;
	t347 = t290 * t386 - t286;
	t346 = 0.1e1 / t313 ^ 2 * t355;
	t345 = t322 * t354;
	t344 = t323 * t354;
	t296 = (t316 * t337 + (-t368 / 0.2e1 - t310 + t314 * t387) * t317) * pkin(7);
	t297 = t314 * t368 / 0.2e1 + t331 * t317 ^ 2 * t387 + (-t316 * t310 - t366) * pkin(7);
	t282 = ((-t324 * t296 / 0.2e1 + t297 * t379) * t311 + (-t360 - t365) * t346) * t334;
	t283 = ((t296 * t379 + t324 * t297 / 0.2e1) * t311 + (-t361 + t364) * t346) * t334;
	t278 = t320 * t282 + t319 * t283;
	t342 = t278 * t345;
	t341 = t343 * t345;
	t340 = t278 * t344;
	t339 = t343 * t344;
	t338 = 0.2e1 * (t284 + t285) * t385;
	t279 = -t319 * t282 + t320 * t283;
	t277 = t343 * t338;
	t274 = t278 * t338;
	t270 = 0.1e1 / t273;
	t269 = t277 * t353 + t343 * t348 + (t292 * t286 - t336 * t343) * pkin(4);
	t268 = (t277 * t352 - t292 * t336 + t343 * t347) * pkin(4);
	t263 = t274 * t353 + t278 * t348 + (-t278 * t336 + t279 * t286) * pkin(4);
	t262 = (t274 * t352 + t347 * t278 - t279 * t336) * pkin(4);
	t261 = 0.1e1 + ((t268 * t351 + t269 * t349 + t275 * t341 + t276 * t339) * t270 - (t268 * t350 + t269 * t351 - t275 * t339 + t276 * t341) * t371) * t372;
	t260 = 0.1e1 + ((t262 * t351 + t263 * t349 + t275 * t342 + t276 * t340) * t270 - (t262 * t350 + t263 * t351 - t275 * t340 + t276 * t342) * t371) * t372;
	t1 = [-t362, -t260 * t359, -t261 * t359, 0; t358, -t260 * t363, -t261 * t363, 0; 0, t260 * t265, t261 * t265, 0; t363, -t260 * t358, -t261 * t358, 0; -t359, -t260 * t362, -t261 * t362, 0; 0, -t260 * t264, -t261 * t264, 0; t328, 0, 0, 0; t326, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:51:45
	% EndTime: 2020-04-15 18:52:14
	% DurationCPUTime: 17.50s
	% Computational Cost: add. (501095->99), mult. (756394->193), div. (33396->11), fcn. (479460->21), ass. (0->124)
	t518 = pkin(4) ^ 2;
	t516 = pkin(7) ^ 2;
	t509 = sin(qJ(2));
	t514 = cos(pkin(19));
	t572 = sin(pkin(19));
	t573 = cos(qJ(2));
	t499 = t509 * t514 - t572 * t573;
	t569 = pkin(7) * t499;
	t582 = -2 * pkin(1);
	t545 = (pkin(1) ^ 2) + t569 * t582;
	t496 = t516 + t545;
	t493 = pkin(3) ^ 2 - pkin(8) ^ 2 + t496;
	t497 = pkin(1) - t569;
	t579 = -pkin(8) - pkin(3);
	t491 = (pkin(7) - t579) * (pkin(7) + t579) + t545;
	t578 = -pkin(8) + pkin(3);
	t492 = (pkin(7) - t578) * (pkin(7) + t578) + t545;
	t522 = sqrt(-t492 * t491);
	t500 = t509 * t572 + t514 * t573;
	t568 = pkin(7) * t500;
	t487 = t493 * t568 + t497 * t522;
	t512 = cos(qJ(3));
	t548 = t512 * t487;
	t557 = t500 * t522;
	t486 = -pkin(7) * t557 + t493 * t497;
	t508 = sin(qJ(3));
	t554 = t508 * t486;
	t494 = 0.1e1 / t496;
	t519 = 0.1e1 / pkin(3);
	t558 = t494 * t519;
	t482 = (t554 / 0.2e1 + t548 / 0.2e1) * t558;
	t549 = t512 * t486;
	t553 = t508 * t487;
	t483 = (-t549 / 0.2e1 + t553 / 0.2e1) * t558;
	t504 = pkin(23) + pkin(22);
	t502 = sin(t504);
	t503 = cos(t504);
	t475 = t482 * t503 - t483 * t502;
	t571 = pkin(4) * t475;
	t581 = -2 * pkin(5);
	t546 = (pkin(5) ^ 2) - t571 * t581;
	t472 = t518 + t546;
	t469 = -pkin(9) ^ 2 + pkin(11) ^ 2 + t472;
	t473 = pkin(5) + t571;
	t577 = -pkin(9) - pkin(11);
	t467 = (pkin(4) - t577) * (pkin(4) + t577) + t546;
	t576 = -pkin(9) + pkin(11);
	t468 = (pkin(4) - t576) * (pkin(4) + t576) + t546;
	t521 = sqrt(-t468 * t467);
	t528 = t482 * t502 + t503 * t483;
	t570 = pkin(4) * t528;
	t458 = t469 * t473 - t521 * t570;
	t459 = t469 * t570 + t473 * t521;
	t506 = cos(pkin(21));
	t470 = 0.1e1 / t472;
	t515 = 0.1e1 / pkin(11);
	t560 = t470 * t515;
	t505 = sin(pkin(21));
	t575 = t505 / 0.2e1;
	t456 = (-t458 * t506 / 0.2e1 + t459 * t575) * t560;
	t453 = 0.1e1 / t456;
	t580 = pkin(4) * pkin(5);
	t543 = 0.1e1 / t472 ^ 2 * t580;
	t529 = t506 * t543;
	t530 = t505 * t543;
	t454 = 0.1e1 / t456 ^ 2;
	t455 = (t458 * t575 + t459 * t506 / 0.2e1) * t560;
	t562 = t454 * t455;
	t583 = (t458 * t530 + t459 * t529) * t453 - (-t458 * t529 + t459 * t530) * t562;
	t574 = -t512 / 0.2e1;
	t449 = qJ(2) + qJ(3) + atan2(t455, t456);
	t448 = cos(t449);
	t507 = sin(qJ(4));
	t567 = t448 * t507;
	t510 = sin(qJ(1));
	t566 = t448 * t510;
	t511 = cos(qJ(4));
	t565 = t448 * t511;
	t513 = cos(qJ(1));
	t564 = t448 * t513;
	t563 = 0.1e1 / (t454 * t455 ^ 2 + 0.1e1) * t515;
	t561 = t470 * t506;
	t544 = pkin(1) * t568;
	t559 = 0.2e1 / t522 * (t491 + t492) * t544;
	t556 = t507 * t510;
	t555 = t507 * t513;
	t447 = sin(t449);
	t552 = t510 * t447;
	t551 = t510 * t511;
	t550 = t511 * t513;
	t547 = t513 * t447;
	t542 = t507 * t552;
	t541 = t507 * t547;
	t540 = t447 * t551;
	t539 = t511 * t547;
	t464 = 0.1e1 / t521;
	t538 = t464 * t473 / 0.2e1;
	t537 = -t464 * t528 / 0.2e1;
	t536 = t470 * t575;
	t535 = -t561 / 0.2e1;
	t534 = t561 / 0.2e1;
	t533 = t518 * t528 * t581;
	t532 = t473 * t581 - t469;
	t531 = 0.1e1 / t496 ^ 2 * t544;
	t523 = 0.2e1 * (t467 + t468) * t580;
	t480 = t497 * t559 / 0.2e1 + t516 * t500 ^ 2 * t582 + (-t493 * t499 - t557) * pkin(7);
	t479 = (t499 * t522 + (-t559 / 0.2e1 - t493 + t497 * t582) * t500) * pkin(7);
	t466 = ((t479 * t574 + t508 * t480 / 0.2e1) * t494 + (-t549 + t553) * t531) * t519;
	t465 = ((-t508 * t479 / 0.2e1 + t480 * t574) * t494 + (-t548 - t554) * t531) * t519;
	t462 = -t465 * t502 + t466 * t503;
	t461 = t465 * t503 + t466 * t502;
	t460 = t528 * t523;
	t457 = t461 * t523;
	t452 = t460 * t538 + t528 * t533 + (t469 * t475 - t521 * t528) * pkin(4);
	t451 = (t460 * t537 - t475 * t521 + t528 * t532) * pkin(4);
	t446 = t448 * t550 + t556;
	t445 = -t448 * t555 + t551;
	t444 = -t448 * t551 + t555;
	t443 = t448 * t556 + t550;
	t442 = t457 * t538 + t461 * t533 + (-t461 * t521 + t462 * t469) * pkin(4);
	t441 = (t457 * t537 + t461 * t532 - t462 * t521) * pkin(4);
	t440 = 0.1e1 + ((t451 * t536 + t452 * t534) * t453 - (t451 * t535 + t452 * t536) * t562 + t583 * t528) * t563;
	t439 = 0.1e1 + ((t441 * t536 + t442 * t534) * t453 - (t441 * t535 + t442 * t536) * t562 + t583 * t461) * t563;
	t1 = [t444, -t439 * t539, -t440 * t539, t445; t446, -t439 * t540, -t440 * t540, -t443; 0, t439 * t565, t440 * t565, -t447 * t507; t443, t439 * t541, t440 * t541, -t446; t445, t439 * t542, t440 * t542, t444; 0, -t439 * t567, -t440 * t567, -t447 * t511; -t552, t439 * t564, t440 * t564, 0; t547, t439 * t566, t440 * t566, 0; 0, t439 * t447, t440 * t447, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:37
	% EndTime: 2020-04-15 18:49:39
	% DurationCPUTime: 0.38s
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
	% StartTime: 2020-04-15 18:49:37
	% EndTime: 2020-04-15 18:49:38
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (5682->43), mult. (8738->74), div. (380->7), fcn. (5634->13), ass. (0->51)
	t169 = -2 * pkin(1);
	t168 = -pkin(8) - pkin(3);
	t167 = -pkin(8) + pkin(3);
	t138 = sin(pkin(23));
	t166 = t138 / 0.2e1;
	t165 = cos(qJ(2));
	t164 = sin(pkin(19));
	t140 = sin(qJ(2));
	t143 = cos(pkin(19));
	t135 = t140 * t143 - t165 * t164;
	t163 = pkin(7) * t135;
	t136 = t140 * t164 + t165 * t143;
	t162 = pkin(7) * t136;
	t153 = (pkin(1) ^ 2) + t163 * t169;
	t127 = (pkin(7) - t168) * (pkin(7) + t168) + t153;
	t128 = (pkin(7) - t167) * (pkin(7) + t167) + t153;
	t147 = sqrt(-t128 * t127);
	t152 = pkin(1) * t162;
	t161 = 0.1e1 / t147 * (t127 + t128) * t152;
	t144 = pkin(7) ^ 2;
	t132 = t144 + t153;
	t130 = 0.1e1 / t132;
	t145 = 0.1e1 / pkin(3);
	t160 = t130 * t145;
	t159 = t136 * t147;
	t139 = cos(pkin(23));
	t158 = t139 * t130;
	t129 = pkin(3) ^ 2 - pkin(8) ^ 2 + t132;
	t133 = pkin(1) - t163;
	t122 = -pkin(7) * t159 + t133 * t129;
	t123 = t129 * t162 + t133 * t147;
	t120 = (-t139 * t122 / 0.2e1 + t123 * t166) * t160;
	t121 = (t139 * t123 / 0.2e1 + t122 * t166) * t160;
	t115 = qJ(2) + atan2(t121, t120);
	t113 = sin(t115);
	t141 = sin(qJ(1));
	t157 = t141 * t113;
	t114 = cos(t115);
	t156 = t141 * t114;
	t142 = cos(qJ(1));
	t155 = t142 * t113;
	t154 = t142 * t114;
	t151 = t130 * t166;
	t150 = 0.1e1 / t132 ^ 2 * t152;
	t149 = t138 * t150;
	t148 = t139 * t150;
	t119 = 0.1e1 / t120 ^ 2;
	t118 = t133 * t161 + t144 * t136 ^ 2 * t169 + (-t135 * t129 - t159) * pkin(7);
	t117 = (t135 * t147 + (t133 * t169 - t129 - t161) * t136) * pkin(7);
	t112 = 0.1e1 + ((t118 * t158 / 0.2e1 + t123 * t148 + t117 * t151 + t122 * t149) / t120 - (-t117 * t158 / 0.2e1 - t122 * t148 + t118 * t151 + t123 * t149) * t121 * t119) / (t121 ^ 2 * t119 + 0.1e1) * t145;
	t1 = [t157, -t112 * t154, 0, 0; -t155, -t112 * t156, 0, 0; 0, -t112 * t113, 0, 0; t156, t112 * t155, 0, 0; -t154, t112 * t157, 0, 0; 0, -t112 * t114, 0, 0; t142, 0, 0, 0; t141, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:36
	% EndTime: 2020-04-15 18:49:37
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (3062->40), mult. (4588->63), div. (148->7), fcn. (2954->11), ass. (0->41)
	t134 = -pkin(2) - pkin(13);
	t133 = -pkin(2) + pkin(13);
	t116 = pkin(6) ^ 2;
	t110 = sin(pkin(20));
	t111 = cos(pkin(20));
	t112 = sin(qJ(3));
	t114 = cos(qJ(3));
	t108 = t114 * t110 + t112 * t111;
	t129 = pkin(6) * t108;
	t122 = 0.2e1 * pkin(1) * t129 + t116;
	t105 = pkin(1) ^ 2 + t122;
	t102 = pkin(2) ^ 2 - pkin(13) ^ 2 + t105;
	t106 = -pkin(1) - t129;
	t100 = (pkin(1) - t134) * (pkin(1) + t134) + t122;
	t101 = (pkin(1) - t133) * (pkin(1) + t133) + t122;
	t119 = sqrt(-t101 * t100);
	t109 = t112 * t110 - t114 * t111;
	t128 = pkin(6) * t109;
	t96 = t102 * t128 - t106 * t119;
	t132 = pkin(6) * t96;
	t103 = 0.1e1 / t105;
	t131 = t103 / 0.2e1;
	t130 = pkin(1) * t109;
	t127 = 0.1e1 / t119 * (t100 + t101) * pkin(1) * t128;
	t113 = sin(qJ(1));
	t117 = 0.1e1 / pkin(2);
	t121 = t117 * t131;
	t123 = t109 * t119;
	t95 = -pkin(6) * t123 - t106 * t102;
	t92 = qJ(2) + atan2(t96 * t121, t95 * t121);
	t90 = sin(t92);
	t89 = t113 * t90;
	t91 = cos(t92);
	t126 = t113 * t91;
	t115 = cos(qJ(1));
	t125 = t115 * t90;
	t124 = t115 * t91;
	t104 = 0.1e1 / t105 ^ 2;
	t94 = 0.1e1 / t95 ^ 2;
	t88 = 0.2e1 * (((-t106 * t127 + (t108 * t102 - t123) * pkin(6)) * t131 + (-t103 * t116 * t109 + t104 * t132) * t130) / t95 - ((-t108 * t119 + (-t102 - t127) * t109) * t131 + (t103 * t106 + t104 * t95) * t130) * t94 * t132) * pkin(2) * t105 * t117 / (t96 ^ 2 * t94 + 0.1e1);
	t1 = [t89, -t124, -t88 * t124, 0; -t125, -t126, -t88 * t126, 0; 0, -t90, -t88 * t90, 0; t126, t125, t88 * t125, 0; -t124, t89, t88 * t89, 0; 0, -t91, -t88 * t91, 0; t115, 0, 0, 0; t113, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:49:37
	% EndTime: 2020-04-15 18:49:38
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (3994->43), mult. (5830->72), div. (248->11), fcn. (3626->12), ass. (0->52)
	t183 = -pkin(2) - pkin(13);
	t182 = -pkin(2) + pkin(13);
	t149 = sin(pkin(20));
	t150 = cos(pkin(20));
	t151 = sin(qJ(3));
	t153 = cos(qJ(3));
	t147 = t153 * t149 + t151 * t150;
	t177 = pkin(6) * t147;
	t167 = pkin(1) * t177;
	t146 = 0.2e1 * t167;
	t157 = pkin(6) ^ 2;
	t168 = pkin(1) ^ 2 + t157;
	t144 = t146 + t168;
	t142 = 0.1e1 / t144;
	t181 = t142 / 0.2e1;
	t159 = 0.1e1 / pkin(2);
	t180 = t159 / 0.2e1;
	t148 = t151 * t149 - t153 * t150;
	t179 = pkin(1) * t148;
	t158 = pkin(2) ^ 2;
	t165 = -pkin(13) ^ 2 + t168;
	t141 = t146 + t158 + t165;
	t145 = -pkin(1) - t177;
	t169 = t146 + t157;
	t137 = (pkin(1) - t183) * (pkin(1) + t183) + t169;
	t138 = (pkin(1) - t182) * (pkin(1) + t182) + t169;
	t174 = t138 * t137;
	t161 = sqrt(-t174);
	t176 = t148 * pkin(6);
	t132 = t141 * t176 - t145 * t161;
	t178 = pkin(6) * t132;
	t175 = 0.1e1 / t161 * (t137 + t138) * pkin(1) * t176;
	t173 = t148 * t161;
	t166 = pkin(6) * t173;
	t131 = -t145 * t141 - t166;
	t140 = t158 - t165 - 0.2e1 * t167;
	t163 = 0.1e1 / pkin(13) * t180;
	t164 = t142 * t180;
	t128 = qJ(2) + atan2(t132 * t164, t131 * t164) + atan2(t161 * t163, t140 * t163);
	t126 = sin(t128);
	t152 = sin(qJ(1));
	t172 = t152 * t126;
	t127 = cos(t128);
	t171 = t152 * t127;
	t154 = cos(qJ(1));
	t170 = t154 * t126;
	t125 = t154 * t127;
	t143 = 0.1e1 / t144 ^ 2;
	t139 = 0.1e1 / t140 ^ 2;
	t130 = 0.1e1 / t131 ^ 2;
	t124 = 0.2e1 * (((-t145 * t175 + (t147 * t141 - t173) * pkin(6)) * t181 + (-t142 * t157 * t148 + t143 * t178) * t179) / t131 - ((-t147 * t161 + (-t141 - t175) * t148) * t181 + (t131 * t143 + t142 * t145) * t179) * t130 * t178) * pkin(2) / (t132 ^ 2 * t130 + 0.1e1) * t144 * t159 + (0.1e1 / t140 * t175 - 0.2e1 * pkin(1) * t139 * t166) / (-t139 * t174 + 0.1e1);
	t1 = [-t172, t125, t124 * t125, 0; t170, t171, t124 * t171, 0; 0, t126, t124 * t126, 0; -t171, -t170, -t124 * t170, 0; t125, -t172, -t124 * t172, 0; 0, t127, t124 * t127, 0; t154, 0, 0, 0; t152, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-15 18:50:05
	% EndTime: 2020-04-15 18:50:13
	% DurationCPUTime: 5.99s
	% Computational Cost: add. (160844->99), mult. (242824->181), div. (10548->14), fcn. (154026->20), ass. (0->105)
	t379 = 2 * pkin(4);
	t378 = -2 * pkin(1);
	t377 = -pkin(8) - pkin(3);
	t376 = -pkin(8) + pkin(3);
	t375 = (-pkin(9) - pkin(11));
	t374 = (-pkin(9) + pkin(11));
	t325 = pkin(5) ^ 2;
	t324 = pkin(7) ^ 2;
	t318 = sin(qJ(2));
	t322 = cos(pkin(19));
	t366 = sin(pkin(19));
	t367 = cos(qJ(2));
	t309 = t318 * t322 - t367 * t366;
	t362 = pkin(7) * t309;
	t343 = (pkin(1) ^ 2) + t362 * t378;
	t306 = t324 + t343;
	t303 = pkin(3) ^ 2 - pkin(8) ^ 2 + t306;
	t307 = pkin(1) - t362;
	t301 = (pkin(7) - t377) * (pkin(7) + t377) + t343;
	t302 = (pkin(7) - t376) * (pkin(7) + t376) + t343;
	t330 = sqrt(-t302 * t301);
	t310 = t318 * t366 + t367 * t322;
	t361 = pkin(7) * t310;
	t297 = t303 * t361 + t307 * t330;
	t320 = cos(qJ(3));
	t347 = t320 * t297;
	t356 = t310 * t330;
	t296 = -pkin(7) * t356 + t303 * t307;
	t317 = sin(qJ(3));
	t352 = t317 * t296;
	t304 = 0.1e1 / t306;
	t327 = 0.1e1 / pkin(3);
	t357 = t304 * t327;
	t292 = (t352 / 0.2e1 + t347 / 0.2e1) * t357;
	t348 = t320 * t296;
	t351 = t317 * t297;
	t293 = (-t348 / 0.2e1 + t351 / 0.2e1) * t357;
	t314 = pkin(23) + pkin(22);
	t312 = sin(t314);
	t313 = cos(t314);
	t281 = t292 * t313 - t293 * t312;
	t364 = t281 * pkin(5);
	t344 = t364 * t379 + t325;
	t273 = ((pkin(4) - t375) * (pkin(4) + t375)) + t344;
	t274 = ((pkin(4) - t374) * (pkin(4) + t374)) + t344;
	t329 = sqrt(-t274 * t273);
	t373 = -0.1e1 / t329 / 0.2e1;
	t278 = (pkin(4) ^ 2) + t344;
	t276 = 0.1e1 / t278;
	t372 = t276 / 0.2e1;
	t342 = pkin(1) * t361;
	t358 = 0.2e1 / t330 * (t301 + t302) * t342;
	t286 = (t309 * t330 + (-t358 / 0.2e1 - t303 + t307 * t378) * t310) * pkin(7);
	t371 = -t286 / 0.2e1;
	t287 = t307 * t358 / 0.2e1 + t324 * t310 ^ 2 * t378 + (-t309 * t303 - t356) * pkin(7);
	t370 = t287 / 0.2e1;
	t315 = sin(pkin(23));
	t369 = t315 / 0.2e1;
	t368 = -t320 / 0.2e1;
	t275 = pkin(9) ^ 2 - pkin(11) ^ 2 + t278;
	t279 = -pkin(4) - t364;
	t337 = t292 * t312 + t313 * t293;
	t363 = pkin(5) * t337;
	t265 = t275 * t363 - t279 * t329;
	t365 = pkin(5) * t265;
	t359 = t296 * t315;
	t355 = t315 * t297;
	t316 = cos(pkin(23));
	t354 = t316 * t296;
	t353 = t316 * t297;
	t264 = -t275 * t279 - t329 * t363;
	t290 = (-t354 / 0.2e1 + t355 / 0.2e1) * t357;
	t291 = (t353 / 0.2e1 + t359 / 0.2e1) * t357;
	t323 = 0.1e1 / pkin(9);
	t340 = t323 * t372;
	t259 = -qJ(2) - atan2(t291, t290) + pkin(22) - atan2(t265 * t340, t264 * t340);
	t257 = sin(t259);
	t319 = sin(qJ(1));
	t350 = t319 * t257;
	t258 = cos(t259);
	t349 = t319 * t258;
	t321 = cos(qJ(1));
	t346 = t321 * t257;
	t345 = t321 * t258;
	t341 = t337 * t373;
	t339 = t279 * t373;
	t338 = 0.1e1 / t306 ^ 2 * t342;
	t263 = 0.1e1 / t264 ^ 2;
	t336 = pkin(9) / (t263 * t265 ^ 2 + 0.1e1) * t278 * t323;
	t335 = 0.1e1 / t264 * t336;
	t334 = pkin(5) * (t273 + t274) * t379;
	t277 = 0.1e1 / t278 ^ 2;
	t333 = pkin(4) * (t264 * t277 + t276 * t279);
	t332 = t263 * t336 * t365;
	t331 = pkin(4) * (-t276 * t325 * t337 + t277 * t365);
	t289 = 0.1e1 / t290 ^ 2;
	t272 = ((t286 * t368 + t317 * t370) * t304 + (-t348 + t351) * t338) * t327;
	t271 = ((t287 * t368 + t317 * t371) * t304 + (-t347 - t352) * t338) * t327;
	t268 = -t271 * t312 + t272 * t313;
	t267 = t271 * t313 + t272 * t312;
	t266 = t337 * t334;
	t261 = t267 * t334;
	t256 = -0.2e1 * ((t266 * t339 + (t281 * t275 - t329 * t337) * pkin(5)) * t372 + t337 * t331) * t335 + 0.2e1 * ((t266 * t341 - t275 * t337 - t281 * t329) * t372 + t337 * t333) * t332;
	t255 = -0.1e1 - 0.2e1 * ((t261 * t339 + (-t267 * t329 + t268 * t275) * pkin(5)) * t372 + t267 * t331) * t335 + 0.2e1 * ((t261 * t341 - t267 * t275 - t268 * t329) * t372 + t267 * t333) * t332 + (-((t286 * t369 + t316 * t370) * t304 + (t353 + t359) * t338) / t290 + ((t287 * t369 + t316 * t371) * t304 + (-t354 + t355) * t338) * t291 * t289) / (t289 * t291 ^ 2 + 0.1e1) * t327;
	t1 = [t350, -t255 * t345, -t256 * t345, 0; -t346, -t255 * t349, -t256 * t349, 0; 0, t255 * t257, t256 * t257, 0; -t349, -t255 * t346, -t256 * t346, 0; t345, -t255 * t350, -t256 * t350, 0; 0, -t255 * t258, -t256 * t258, 0; t321, 0, 0, 0; t319, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end