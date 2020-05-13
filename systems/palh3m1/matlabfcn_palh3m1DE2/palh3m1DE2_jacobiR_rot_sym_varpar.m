% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m1DE2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m1DE2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_jacobiR_rot_sym_varpar: pkin has to be [19x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:41
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
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
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
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
	% StartTime: 2020-04-20 16:20:42
	% EndTime: 2020-04-20 16:20:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = qJ(2) + qJ(3);
	t31 = cos(t32);
	t30 = sin(t32);
	t29 = t34 * t31;
	t28 = t34 * t30;
	t27 = t33 * t31;
	t26 = t33 * t30;
	t1 = [t27, t28, t28, 0; -t29, t26, t26, 0; 0, -t31, -t31, 0; -t26, t29, t29, 0; t28, t27, t27, 0; 0, t30, t30, 0; t34, 0, 0, 0; t33, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:30
	% EndTime: 2020-04-20 16:22:01
	% DurationCPUTime: 12.00s
	% Computational Cost: add. (322856->91), mult. (487440->175), div. (21464->11), fcn. (308998->19), ass. (0->111)
	t386 = -2 * pkin(1);
	t385 = -2 * pkin(4);
	t384 = pkin(3) * pkin(4);
	t383 = -pkin(6) - pkin(2);
	t382 = -pkin(6) + pkin(2);
	t381 = (-pkin(8) - pkin(10));
	t380 = (-pkin(8) + pkin(10));
	t321 = sin(pkin(17));
	t379 = t321 / 0.2e1;
	t323 = sin(qJ(3));
	t378 = t323 / 0.2e1;
	t377 = cos(pkin(16));
	t376 = cos(qJ(2));
	t330 = pkin(5) ^ 2;
	t324 = sin(qJ(2));
	t326 = sin(pkin(16));
	t315 = t324 * t326 - t376 * t377;
	t373 = pkin(5) * t315;
	t355 = (pkin(1) ^ 2) + t373 * t386;
	t312 = t330 + t355;
	t310 = 0.1e1 / t312;
	t333 = 0.1e1 / pkin(2);
	t362 = t310 * t333;
	t309 = pkin(2) ^ 2 - pkin(6) ^ 2 + t312;
	t313 = pkin(1) - t373;
	t307 = (pkin(5) - t383) * (pkin(5) + t383) + t355;
	t308 = (pkin(5) - t382) * (pkin(5) + t382) + t355;
	t336 = sqrt(-t308 * t307);
	t316 = t324 * t377 + t376 * t326;
	t372 = pkin(5) * t316;
	t303 = t309 * t372 + t313 * t336;
	t365 = t303 * t323;
	t361 = t316 * t336;
	t302 = -pkin(5) * t361 + t313 * t309;
	t327 = cos(qJ(3));
	t366 = t302 * t327;
	t298 = (-t366 / 0.2e1 + t365 / 0.2e1) * t362;
	t364 = t303 * t327;
	t367 = t302 * t323;
	t299 = (t364 / 0.2e1 + t367 / 0.2e1) * t362;
	t320 = pkin(18) + pkin(19);
	t318 = sin(t320);
	t319 = cos(t320);
	t338 = t318 * t298 - t319 * t299;
	t375 = pkin(3) * t338;
	t293 = t319 * t298 + t318 * t299;
	t374 = pkin(3) * t293;
	t332 = pkin(3) ^ 2;
	t356 = (pkin(4) ^ 2) - t374 * t385;
	t289 = t332 + t356;
	t286 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t289;
	t290 = pkin(4) + t374;
	t284 = ((pkin(3) - t381) * (pkin(3) + t381)) + t356;
	t285 = ((pkin(3) - t380) * (pkin(3) + t380)) + t356;
	t335 = sqrt(-t285 * t284);
	t275 = t290 * t286 - t335 * t375;
	t276 = t286 * t375 + t290 * t335;
	t322 = cos(pkin(17));
	t287 = 0.1e1 / t289;
	t329 = 0.1e1 / pkin(10);
	t368 = t287 * t329;
	t272 = (-t275 * t322 / 0.2e1 + t276 * t379) * t368;
	t271 = 0.1e1 / t272 ^ 2;
	t273 = (t276 * t322 / 0.2e1 + t275 * t379) * t368;
	t371 = 0.1e1 / (t273 ^ 2 * t271 + 0.1e1) * t329;
	t370 = t271 * t273;
	t369 = t287 * t322;
	t354 = pkin(1) * t372;
	t363 = 0.2e1 / t336 * (t307 + t308) * t354;
	t266 = qJ(2) + qJ(3) + atan2(t273, t272);
	t264 = sin(t266);
	t325 = sin(qJ(1));
	t360 = t325 * t264;
	t265 = cos(t266);
	t359 = t325 * t265;
	t328 = cos(qJ(1));
	t358 = t328 * t264;
	t357 = t328 * t265;
	t353 = 0.1e1 / t289 ^ 2 * t384;
	t281 = 0.1e1 / t335;
	t352 = t281 * t290 / 0.2e1;
	t351 = -t281 * t338 / 0.2e1;
	t350 = t287 * t379;
	t349 = -t369 / 0.2e1;
	t348 = t369 / 0.2e1;
	t347 = t332 * t338 * t385;
	t346 = t290 * t385 - t286;
	t345 = 0.1e1 / t312 ^ 2 * t354;
	t344 = t321 * t353;
	t343 = t322 * t353;
	t296 = (t315 * t336 + (-t363 / 0.2e1 - t309 + t313 * t386) * t316) * pkin(5);
	t297 = t313 * t363 / 0.2e1 + t330 * t316 ^ 2 * t386 + (-t315 * t309 - t361) * pkin(5);
	t282 = ((t297 * t327 / 0.2e1 + t296 * t378) * t310 + (t364 + t367) * t345) * t333;
	t283 = ((-t296 * t327 / 0.2e1 + t297 * t378) * t310 + (t365 - t366) * t345) * t333;
	t279 = -t318 * t282 - t319 * t283;
	t342 = t279 * t344;
	t341 = t338 * t344;
	t340 = t279 * t343;
	t339 = t338 * t343;
	t337 = 0.2e1 * (t284 + t285) * t384;
	t278 = -t319 * t282 + t318 * t283;
	t277 = t338 * t337;
	t274 = t279 * t337;
	t270 = 0.1e1 / t272;
	t269 = t277 * t352 + t338 * t347 + (t293 * t286 - t335 * t338) * pkin(3);
	t268 = (t277 * t351 - t293 * t335 + t338 * t346) * pkin(3);
	t263 = t274 * t352 + t279 * t347 + (t278 * t286 - t279 * t335) * pkin(3);
	t262 = (t274 * t351 - t278 * t335 + t346 * t279) * pkin(3);
	t261 = 0.1e1 + ((t268 * t350 + t269 * t348 + t275 * t341 + t276 * t339) * t270 - (t268 * t349 + t269 * t350 - t275 * t339 + t276 * t341) * t370) * t371;
	t260 = 0.1e1 + ((t262 * t350 + t263 * t348 + t275 * t342 + t276 * t340) * t270 - (t262 * t349 + t263 * t350 - t275 * t340 + t276 * t342) * t370) * t371;
	t1 = [t359, t260 * t358, t261 * t358, 0; -t357, t260 * t360, t261 * t360, 0; 0, -t260 * t265, -t261 * t265, 0; -t360, t260 * t357, t261 * t357, 0; t358, t260 * t359, t261 * t359, 0; 0, t260 * t264, t261 * t264, 0; t328, 0, 0, 0; t325, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:23:33
	% EndTime: 2020-04-20 16:24:12
	% DurationCPUTime: 17.91s
	% Computational Cost: add. (501099->106), mult. (756394->197), div. (33396->11), fcn. (479460->21), ass. (0->125)
	t579 = -2 * pkin(1);
	t578 = -2 * pkin(4);
	t577 = pkin(3) * pkin(4);
	t576 = -pkin(6) - pkin(2);
	t575 = -pkin(6) + pkin(2);
	t574 = (-pkin(8) - pkin(10));
	t573 = (-pkin(8) + pkin(10));
	t504 = sin(pkin(17));
	t572 = t504 / 0.2e1;
	t507 = sin(qJ(3));
	t571 = t507 / 0.2e1;
	t570 = cos(pkin(16));
	t569 = cos(qJ(2));
	t515 = pkin(5) ^ 2;
	t508 = sin(qJ(2));
	t510 = sin(pkin(16));
	t498 = t508 * t510 - t569 * t570;
	t566 = pkin(5) * t498;
	t544 = (pkin(1) ^ 2) + t566 * t579;
	t495 = t515 + t544;
	t493 = 0.1e1 / t495;
	t518 = 0.1e1 / pkin(2);
	t551 = t493 * t518;
	t492 = pkin(2) ^ 2 - pkin(6) ^ 2 + t495;
	t496 = pkin(1) - t566;
	t490 = (pkin(5) - t576) * (pkin(5) + t576) + t544;
	t491 = (pkin(5) - t575) * (pkin(5) + t575) + t544;
	t521 = sqrt(-t491 * t490);
	t499 = t508 * t570 + t569 * t510;
	t565 = pkin(5) * t499;
	t486 = t492 * t565 + t496 * t521;
	t554 = t486 * t507;
	t550 = t499 * t521;
	t485 = -pkin(5) * t550 + t496 * t492;
	t512 = cos(qJ(3));
	t555 = t485 * t512;
	t481 = (-t555 / 0.2e1 + t554 / 0.2e1) * t551;
	t553 = t486 * t512;
	t556 = t485 * t507;
	t482 = (t553 / 0.2e1 + t556 / 0.2e1) * t551;
	t503 = pkin(18) + pkin(19);
	t501 = sin(t503);
	t502 = cos(t503);
	t523 = t501 * t481 - t502 * t482;
	t568 = pkin(3) * t523;
	t476 = t502 * t481 + t501 * t482;
	t567 = pkin(3) * t476;
	t517 = pkin(3) ^ 2;
	t545 = (pkin(4) ^ 2) - t567 * t578;
	t472 = t517 + t545;
	t469 = -pkin(8) ^ 2 + pkin(10) ^ 2 + t472;
	t473 = pkin(4) + t567;
	t467 = ((pkin(3) - t574) * (pkin(3) + t574)) + t545;
	t468 = ((pkin(3) - t573) * (pkin(3) + t573)) + t545;
	t520 = sqrt(-t468 * t467);
	t458 = t473 * t469 - t520 * t568;
	t459 = t469 * t568 + t473 * t520;
	t505 = cos(pkin(17));
	t470 = 0.1e1 / t472;
	t514 = 0.1e1 / pkin(10);
	t557 = t470 * t514;
	t455 = (-t458 * t505 / 0.2e1 + t459 * t572) * t557;
	t456 = (t459 * t505 / 0.2e1 + t458 * t572) * t557;
	t449 = qJ(2) + qJ(3) + atan2(t456, t455);
	t448 = cos(t449);
	t506 = sin(qJ(4));
	t564 = t448 * t506;
	t509 = sin(qJ(1));
	t563 = t448 * t509;
	t511 = cos(qJ(4));
	t562 = t448 * t511;
	t513 = cos(qJ(1));
	t561 = t448 * t513;
	t454 = 0.1e1 / t455 ^ 2;
	t560 = 0.1e1 / (t456 ^ 2 * t454 + 0.1e1) * t514;
	t559 = t454 * t456;
	t558 = t470 * t505;
	t543 = pkin(1) * t565;
	t552 = 0.2e1 / t521 * (t490 + t491) * t543;
	t549 = t509 * t506;
	t548 = t509 * t511;
	t547 = t513 * t506;
	t546 = t513 * t511;
	t542 = 0.1e1 / t472 ^ 2 * t577;
	t447 = sin(t449);
	t541 = t447 * t549;
	t540 = t447 * t547;
	t539 = t447 * t548;
	t538 = t447 * t546;
	t464 = 0.1e1 / t520;
	t537 = t464 * t473 / 0.2e1;
	t536 = -t464 * t523 / 0.2e1;
	t535 = t470 * t572;
	t534 = -t558 / 0.2e1;
	t533 = t558 / 0.2e1;
	t532 = t523 * t517 * t578;
	t531 = t473 * t578 - t469;
	t530 = 0.1e1 / t495 ^ 2 * t543;
	t479 = (t498 * t521 + (-t552 / 0.2e1 - t492 + t496 * t579) * t499) * pkin(5);
	t480 = t496 * t552 / 0.2e1 + t515 * t499 ^ 2 * t579 + (-t498 * t492 - t550) * pkin(5);
	t465 = ((t480 * t512 / 0.2e1 + t479 * t571) * t493 + (t553 + t556) * t530) * t518;
	t466 = ((-t479 * t512 / 0.2e1 + t480 * t571) * t493 + (t554 - t555) * t530) * t518;
	t462 = -t501 * t465 - t502 * t466;
	t529 = t462 * t542;
	t528 = t523 * t542;
	t527 = t458 * t529;
	t526 = t458 * t528;
	t525 = t459 * t529;
	t524 = t459 * t528;
	t522 = 0.2e1 * (t467 + t468) * t577;
	t461 = -t502 * t465 + t501 * t466;
	t460 = t523 * t522;
	t457 = t462 * t522;
	t453 = 0.1e1 / t455;
	t452 = t460 * t537 + t523 * t532 + (t476 * t469 - t520 * t523) * pkin(3);
	t451 = (t460 * t536 - t476 * t520 + t523 * t531) * pkin(3);
	t446 = t448 * t546 - t549;
	t445 = t448 * t547 + t548;
	t444 = t448 * t548 + t547;
	t443 = t448 * t549 - t546;
	t442 = t457 * t537 + t462 * t532 + (t461 * t469 - t462 * t520) * pkin(3);
	t441 = (t457 * t536 - t461 * t520 + t531 * t462) * pkin(3);
	t440 = 0.1e1 + ((t451 * t535 + t452 * t533 + t504 * t526 + t505 * t524) * t453 - (t451 * t534 + t452 * t535 + t504 * t524 - t505 * t526) * t559) * t560;
	t439 = 0.1e1 + ((t441 * t535 + t442 * t533 + t504 * t527 + t505 * t525) * t453 - (t441 * t534 + t442 * t535 + t504 * t525 - t505 * t527) * t559) * t560;
	t1 = [t444, t439 * t538, t440 * t538, t445; -t446, t439 * t539, t440 * t539, t443; 0, -t439 * t562, -t440 * t562, t447 * t506; -t443, -t439 * t540, -t440 * t540, t446; t445, -t439 * t541, -t440 * t541, t444; 0, t439 * t564, t440 * t564, t447 * t511; t509 * t447, -t439 * t561, -t440 * t561, 0; -t513 * t447, -t439 * t563, -t440 * t563, 0; 0, -t439 * t447, -t440 * t447, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:20:44
	% EndTime: 2020-04-20 16:20:45
	% DurationCPUTime: 0.38s
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
	% StartTime: 2020-04-20 16:20:44
	% EndTime: 2020-04-20 16:20:45
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (5683->44), mult. (8738->74), div. (380->7), fcn. (5634->13), ass. (0->51)
	t169 = -2 * pkin(1);
	t168 = -pkin(6) - pkin(2);
	t167 = -pkin(6) + pkin(2);
	t138 = sin(pkin(19));
	t166 = t138 / 0.2e1;
	t165 = cos(pkin(16));
	t164 = cos(qJ(2));
	t140 = sin(qJ(2));
	t142 = sin(pkin(16));
	t135 = t140 * t142 - t164 * t165;
	t163 = pkin(5) * t135;
	t136 = t140 * t165 + t164 * t142;
	t162 = pkin(5) * t136;
	t153 = (pkin(1) ^ 2) + t163 * t169;
	t127 = (pkin(5) - t168) * (pkin(5) + t168) + t153;
	t128 = (pkin(5) - t167) * (pkin(5) + t167) + t153;
	t147 = sqrt(-t128 * t127);
	t152 = pkin(1) * t162;
	t161 = 0.1e1 / t147 * (t127 + t128) * t152;
	t144 = pkin(5) ^ 2;
	t132 = t144 + t153;
	t130 = 0.1e1 / t132;
	t139 = cos(pkin(19));
	t160 = t130 * t139;
	t145 = 0.1e1 / pkin(2);
	t159 = t130 * t145;
	t158 = t136 * t147;
	t129 = pkin(2) ^ 2 - pkin(6) ^ 2 + t132;
	t133 = pkin(1) - t163;
	t122 = -pkin(5) * t158 + t133 * t129;
	t123 = t129 * t162 + t133 * t147;
	t120 = (-t122 * t139 / 0.2e1 + t123 * t166) * t159;
	t121 = (t123 * t139 / 0.2e1 + t122 * t166) * t159;
	t115 = qJ(2) + atan2(t121, t120);
	t113 = sin(t115);
	t141 = sin(qJ(1));
	t157 = t141 * t113;
	t114 = cos(t115);
	t156 = t141 * t114;
	t143 = cos(qJ(1));
	t155 = t143 * t113;
	t154 = t143 * t114;
	t151 = t130 * t166;
	t150 = 0.1e1 / t132 ^ 2 * t152;
	t149 = t138 * t150;
	t148 = t139 * t150;
	t119 = 0.1e1 / t120 ^ 2;
	t118 = t133 * t161 + t144 * t136 ^ 2 * t169 + (-t135 * t129 - t158) * pkin(5);
	t117 = (t135 * t147 + (t133 * t169 - t129 - t161) * t136) * pkin(5);
	t112 = 0.1e1 + ((t118 * t160 / 0.2e1 + t123 * t148 + t117 * t151 + t122 * t149) / t120 - (-t117 * t160 / 0.2e1 - t122 * t148 + t118 * t151 + t123 * t149) * t121 * t119) / (t121 ^ 2 * t119 + 0.1e1) * t145;
	t1 = [-t156, -t112 * t155, 0, 0; t154, -t112 * t157, 0, 0; 0, t112 * t114, 0, 0; t157, -t112 * t154, 0, 0; -t155, -t112 * t156, 0, 0; 0, -t112 * t113, 0, 0; t143, 0, 0, 0; t141, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-20 16:21:14
	% EndTime: 2020-04-20 16:21:27
	% DurationCPUTime: 6.13s
	% Computational Cost: add. (160838->93), mult. (242824->181), div. (10548->14), fcn. (154026->20), ass. (0->105)
	t375 = 2 * pkin(3);
	t374 = -2 * pkin(1);
	t373 = -pkin(6) - pkin(2);
	t372 = -pkin(6) + pkin(2);
	t371 = (-pkin(8) - pkin(10));
	t370 = (-pkin(8) + pkin(10));
	t321 = pkin(4) ^ 2;
	t320 = pkin(5) ^ 2;
	t314 = sin(qJ(2));
	t316 = sin(pkin(16));
	t362 = cos(qJ(2));
	t363 = cos(pkin(16));
	t305 = t314 * t316 - t362 * t363;
	t358 = pkin(5) * t305;
	t339 = (pkin(1) ^ 2) + t358 * t374;
	t302 = t320 + t339;
	t300 = 0.1e1 / t302;
	t323 = 0.1e1 / pkin(2);
	t346 = t300 * t323;
	t299 = pkin(2) ^ 2 - pkin(6) ^ 2 + t302;
	t303 = pkin(1) - t358;
	t297 = (pkin(5) - t373) * (pkin(5) + t373) + t339;
	t298 = (pkin(5) - t372) * (pkin(5) + t372) + t339;
	t326 = sqrt(-t298 * t297);
	t306 = t314 * t363 + t362 * t316;
	t357 = pkin(5) * t306;
	t293 = t299 * t357 + t303 * t326;
	t313 = sin(qJ(3));
	t349 = t293 * t313;
	t345 = t306 * t326;
	t292 = -pkin(5) * t345 + t303 * t299;
	t317 = cos(qJ(3));
	t352 = t292 * t317;
	t288 = (-t352 / 0.2e1 + t349 / 0.2e1) * t346;
	t348 = t293 * t317;
	t353 = t292 * t313;
	t289 = (t348 / 0.2e1 + t353 / 0.2e1) * t346;
	t310 = pkin(18) + pkin(19);
	t308 = sin(t310);
	t309 = cos(t310);
	t279 = t309 * t288 + t308 * t289;
	t359 = t279 * pkin(4);
	t340 = t359 * t375 + t321;
	t270 = ((pkin(3) - t371) * (pkin(3) + t371)) + t340;
	t271 = ((pkin(3) - t370) * (pkin(3) + t370)) + t340;
	t325 = sqrt(-t271 * t270);
	t369 = -0.1e1 / t325 / 0.2e1;
	t275 = (pkin(3) ^ 2) + t340;
	t273 = 0.1e1 / t275;
	t368 = t273 / 0.2e1;
	t338 = pkin(1) * t357;
	t347 = 0.2e1 / t326 * (t297 + t298) * t338;
	t283 = (t305 * t326 + (-t347 / 0.2e1 - t299 + t303 * t374) * t306) * pkin(5);
	t367 = -t283 / 0.2e1;
	t284 = t303 * t347 / 0.2e1 + t320 * t306 ^ 2 * t374 + (-t305 * t299 - t345) * pkin(5);
	t366 = t284 / 0.2e1;
	t311 = sin(pkin(19));
	t365 = t311 / 0.2e1;
	t364 = t313 / 0.2e1;
	t272 = pkin(8) ^ 2 - pkin(10) ^ 2 + t275;
	t276 = -pkin(3) - t359;
	t331 = t308 * t288 - t309 * t289;
	t360 = pkin(4) * t331;
	t262 = t272 * t360 - t276 * t325;
	t361 = pkin(4) * t262;
	t355 = t292 * t311;
	t312 = cos(pkin(19));
	t354 = t292 * t312;
	t351 = t293 * t311;
	t350 = t293 * t312;
	t261 = -t276 * t272 - t325 * t360;
	t286 = (-t354 / 0.2e1 + t351 / 0.2e1) * t346;
	t287 = (t350 / 0.2e1 + t355 / 0.2e1) * t346;
	t319 = 0.1e1 / pkin(8);
	t336 = t319 * t368;
	t256 = -qJ(2) - atan2(t287, t286) + pkin(18) - atan2(t262 * t336, t261 * t336);
	t254 = sin(t256);
	t315 = sin(qJ(1));
	t344 = t315 * t254;
	t255 = cos(t256);
	t343 = t315 * t255;
	t318 = cos(qJ(1));
	t342 = t318 * t254;
	t341 = t318 * t255;
	t337 = t331 * t369;
	t335 = t276 * t369;
	t334 = 0.1e1 / t302 ^ 2 * t338;
	t260 = 0.1e1 / t261 ^ 2;
	t333 = pkin(8) / (t260 * t262 ^ 2 + 0.1e1) * t275 * t319;
	t332 = 0.1e1 / t261 * t333;
	t330 = pkin(4) * (t270 + t271) * t375;
	t274 = 0.1e1 / t275 ^ 2;
	t329 = pkin(3) * (t261 * t274 + t273 * t276);
	t328 = t260 * t333 * t361;
	t327 = pkin(3) * (-t273 * t321 * t331 + t274 * t361);
	t285 = 0.1e1 / t286 ^ 2;
	t269 = ((t284 * t364 + t317 * t367) * t300 + (t349 - t352) * t334) * t323;
	t268 = ((t283 * t364 + t317 * t366) * t300 + (t348 + t353) * t334) * t323;
	t265 = -t308 * t268 - t309 * t269;
	t264 = -t309 * t268 + t308 * t269;
	t263 = t331 * t330;
	t258 = t265 * t330;
	t253 = -0.2e1 * ((t263 * t335 + (t279 * t272 - t325 * t331) * pkin(4)) * t368 + t331 * t327) * t332 + 0.2e1 * ((t263 * t337 - t272 * t331 - t279 * t325) * t368 + t331 * t329) * t328;
	t252 = -0.1e1 - 0.2e1 * ((t258 * t335 + (t264 * t272 - t265 * t325) * pkin(4)) * t368 + t265 * t327) * t332 + 0.2e1 * ((t258 * t337 - t264 * t325 - t265 * t272) * t368 + t265 * t329) * t328 + (-((t283 * t365 + t312 * t366) * t300 + (t350 + t355) * t334) / t286 + ((t284 * t365 + t312 * t367) * t300 + (t351 - t354) * t334) * t287 * t285) / (t287 ^ 2 * t285 + 0.1e1) * t323;
	t1 = [t343, t252 * t342, t253 * t342, 0; -t341, t252 * t344, t253 * t344, 0; 0, t252 * t255, t253 * t255, 0; t344, -t252 * t341, -t253 * t341, 0; -t342, -t252 * t343, -t253 * t343, 0; 0, t252 * t254, t253 * t254, 0; t318, 0, 0, 0; t315, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end