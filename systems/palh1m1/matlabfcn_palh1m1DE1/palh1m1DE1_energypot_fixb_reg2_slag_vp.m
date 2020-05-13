% Calculate inertial parameters regressor of potential energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m1DE1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energypot_fixb_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 05:12:26
% EndTime: 2020-04-14 05:12:34
% DurationCPUTime: 7.03s
% Computational Cost: add. (159769->175), mult. (240205->285), div. (11448->9), fcn. (152251->46), ass. (0->147)
t480 = -pkin(2) - pkin(13);
t479 = -pkin(2) + pkin(13);
t478 = -pkin(8) - pkin(3);
t477 = -pkin(8) + pkin(3);
t476 = -pkin(9) - pkin(11);
t475 = -pkin(9) + pkin(11);
t474 = g(3) * pkin(14);
t430 = pkin(1) ^ 2;
t404 = sin(qJ(2));
t406 = sin(pkin(19));
t410 = cos(qJ(2));
t412 = cos(pkin(19));
t383 = t404 * t412 - t410 * t406;
t462 = pkin(7) * t383;
t452 = -0.2e1 * pkin(1) * t462 + t430;
t445 = pkin(7) ^ 2 + t452;
t449 = pkin(3) ^ 2 - pkin(8) ^ 2;
t365 = t445 + t449;
t376 = pkin(1) - t462;
t361 = sqrt(-((pkin(7) - t477) * (pkin(7) + t477) + t452) * ((pkin(7) - t478) * (pkin(7) + t478) + t452));
t386 = t404 * t406 + t410 * t412;
t457 = t361 * t386;
t354 = -pkin(7) * t457 + t376 * t365;
t473 = -t354 / 0.2e1;
t355 = pkin(7) * t386 * t365 + t376 * t361;
t472 = t355 / 0.2e1;
t471 = sin(pkin(23)) / 0.2e1;
t470 = sin(pkin(21)) / 0.2e1;
t409 = cos(qJ(3));
t469 = -t409 / 0.2e1;
t468 = cos(pkin(18)) / 0.2e1;
t467 = 0.1e1 / pkin(2) / 0.2e1;
t466 = pkin(1) * t404;
t403 = sin(qJ(3));
t367 = 0.1e1 / t445;
t455 = t367 / pkin(3);
t351 = (t355 * t469 + t403 * t473) * t455;
t352 = (t354 * t469 + t403 * t472) * t455;
t393 = pkin(23) + pkin(22);
t388 = sin(t393);
t389 = cos(t393);
t340 = t389 * t351 + t388 * t352;
t465 = pkin(5) * t340;
t397 = sin(pkin(20));
t401 = cos(pkin(20));
t464 = pkin(6) * (-t409 * t397 - t403 * t401);
t463 = pkin(6) * (t403 * t397 - t409 * t401);
t454 = -0.2e1 * pkin(4) * t465 + pkin(5) ^ 2;
t446 = pkin(4) ^ 2 + t454;
t451 = pkin(9) ^ 2 - pkin(11) ^ 2;
t333 = t446 - t451;
t338 = -pkin(4) * t340 + pkin(5);
t321 = sqrt(-((pkin(4) - t475) * (pkin(4) + t475) + t454) * ((pkin(4) - t476) * (pkin(4) + t476) + t454));
t341 = -t388 * t351 + t389 * t352;
t459 = t321 * t341;
t319 = -pkin(4) * t459 + t338 * t333;
t320 = pkin(4) * t341 * t333 + t338 * t321;
t400 = cos(pkin(21));
t336 = 0.1e1 / t446;
t458 = t336 / pkin(11);
t316 = atan2((t319 * t470 + t320 * t400 / 0.2e1) * t458, (-t319 * t400 / 0.2e1 + t320 * t470) * t458);
t314 = sin(t316);
t315 = cos(t316);
t384 = t410 * t403 + t404 * t409;
t385 = -t404 * t403 + t410 * t409;
t308 = t385 * t314 + t384 * t315;
t461 = g(3) * t308;
t460 = t410 * pkin(1) + pkin(14);
t456 = t367 / pkin(8);
t448 = pkin(1) * t464;
t375 = -0.2e1 * t448;
t423 = pkin(6) ^ 2;
t453 = t375 + t423;
t450 = t423 + t430;
t447 = t384 * pkin(5) + t460;
t444 = -pkin(13) ^ 2 + t450;
t443 = t336 / pkin(9) / 0.2e1;
t442 = 0.1e1 / (t375 + t450) * t467;
t441 = 0.1e1 / pkin(13) * t467;
t405 = sin(qJ(1));
t390 = t405 * pkin(16);
t440 = -t405 * t466 + t390;
t411 = cos(qJ(1));
t391 = t411 * pkin(16);
t439 = -t411 * t466 + t391;
t438 = g(1) * t411 + g(2) * t405;
t398 = cos(pkin(23));
t346 = atan2((t354 * t471 + t398 * t472) * t455, (t355 * t471 + t398 * t473) * t455);
t342 = sin(t346);
t343 = cos(t346);
t326 = -t410 * t342 - t404 * t343;
t437 = t404 * t342 - t410 * t343;
t360 = sqrt(-((pkin(1) - t479) * (pkin(1) + t479) + t453) * ((pkin(1) - t480) * (pkin(1) + t480) + t453));
t428 = pkin(2) ^ 2;
t363 = t375 + t428 + t444;
t374 = -pkin(1) + t464;
t350 = atan2((-t374 * t360 + t363 * t463) * t442, (-t360 * t463 - t374 * t363) * t442);
t348 = sin(t350);
t349 = cos(t350);
t334 = -t410 * t348 - t404 * t349;
t436 = t404 * t348 - t410 * t349;
t370 = t385 * t405;
t435 = t370 * pkin(5) + t440;
t372 = t385 * t411;
t434 = t372 * pkin(5) + t439;
t371 = t384 * t405;
t303 = t370 * t314 + t371 * t315;
t373 = t384 * t411;
t305 = t372 * t314 + t373 * t315;
t307 = t384 * t314 - t385 * t315;
t433 = g(1) * t305 + g(2) * t303 + g(3) * t307;
t362 = -g(1) * t439 - g(2) * t440 - g(3) * t460;
t332 = t446 + t451;
t337 = -pkin(4) + t465;
t432 = atan2((pkin(5) * t341 * t332 - t337 * t321) * t443, (-pkin(5) * t459 - t337 * t332) * t443);
t431 = sin(t432);
t408 = cos(qJ(4));
t407 = sin(pkin(18));
t402 = sin(qJ(4));
t399 = cos(pkin(22));
t395 = sin(pkin(22));
t387 = -g(1) * t405 + g(2) * t411;
t380 = -t438 * pkin(16) - t474;
t377 = pkin(1) * t383 - pkin(7);
t364 = t445 - t449;
t359 = atan2(t360 * t441, (t428 - t444 + 0.2e1 * t448) * t441);
t358 = cos(t359);
t357 = sin(t359);
t356 = pkin(1) * t386 * t364 - t377 * t361;
t353 = -pkin(1) * t457 - t377 * t364;
t347 = atan2((t356 * t468 + t353 * t407 / 0.2e1) * t456, (t353 * t468 - t407 * t356 / 0.2e1) * t456);
t345 = cos(t347);
t344 = sin(t347);
t331 = t334 * t411;
t330 = t436 * t411;
t329 = t334 * t405;
t328 = t436 * t405;
t325 = t326 * t411;
t324 = t437 * t411;
t323 = t326 * t405;
t322 = t437 * t405;
t318 = cos(t432);
t313 = -t399 * t318 - t395 * t431;
t312 = t395 * t318 - t399 * t431;
t306 = -t373 * t314 + t372 * t315;
t304 = -t371 * t314 + t370 * t315;
t1 = [0, 0, 0, 0, 0, 0, -t438, -t387, -g(3), -t474, 0, 0, 0, 0, 0, 0, -g(3) * t410 + t438 * t404, g(3) * t404 + t438 * t410, t387, t380, 0, 0, 0, 0, 0, 0, -g(1) * t372 - g(2) * t370 - g(3) * t384, g(1) * t373 + g(2) * t371 - g(3) * t385, t387, t362, 0, 0, 0, 0, 0, 0, -g(1) * t306 - g(2) * t304 - t461, t433, t387, -g(1) * t434 - g(2) * t435 - g(3) * t447, 0, 0, 0, 0, 0, 0, -t408 * t461 - g(2) * (t304 * t408 - t411 * t402) - g(1) * (t306 * t408 + t405 * t402), t402 * t461 - g(2) * (-t304 * t402 - t411 * t408) - g(1) * (-t306 * t402 + t405 * t408), -t433, -g(3) * (t308 * pkin(10) + t307 * pkin(12) + t447) - g(2) * (t304 * pkin(10) + t303 * pkin(12) + t435) - g(1) * (t306 * pkin(10) + t305 * pkin(12) + t434), 0, 0, 0, 0, 0, 0, -g(3) * t344 - t438 * t345, -g(3) * t345 + t438 * t344, t387, -g(3) * (-pkin(17) + pkin(14)) + t438 * pkin(15), 0, 0, 0, 0, 0, 0, -g(1) * t325 - g(2) * t323 + g(3) * t437, -g(1) * t324 - g(2) * t322 - g(3) * t326, t387, t362, 0, 0, 0, 0, 0, 0, -g(1) * t331 - g(2) * t329 + g(3) * t436, -g(1) * t330 - g(2) * t328 - g(3) * t334, t387, t380, 0, 0, 0, 0, 0, 0, -g(3) * (-t334 * t357 + t358 * t436) - g(2) * (-t328 * t357 - t329 * t358) - g(1) * (-t330 * t357 - t331 * t358), -g(3) * (-t334 * t358 - t357 * t436) - g(2) * (-t328 * t358 + t329 * t357) - g(1) * (-t330 * t358 + t331 * t357), t387, -g(3) * (-pkin(2) * t436 + pkin(14)) - g(2) * (t329 * pkin(2) + t390) - g(1) * (t331 * pkin(2) + t391), 0, 0, 0, 0, 0, 0, -g(3) * (t326 * t312 - t313 * t437) - g(2) * (t322 * t312 + t323 * t313) - g(1) * (t324 * t312 + t325 * t313), -g(3) * (t312 * t437 + t326 * t313) - g(2) * (-t323 * t312 + t322 * t313) - g(1) * (-t325 * t312 + t324 * t313), t387, (-g(3) * (-t326 * t395 - t399 * t437) - g(2) * (-t322 * t395 + t323 * t399) - g(1) * (-t324 * t395 + t325 * t399)) * pkin(4) + t362;];
U_reg = t1;
