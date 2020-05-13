% Calculate inertial parameters regressor of potential energy for
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m1TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energypot_fixb_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 04:17:41
% EndTime: 2020-04-13 04:17:46
% DurationCPUTime: 5.24s
% Computational Cost: add. (80011->175), mult. (120409->296), div. (5724->9), fcn. (76231->28), ass. (0->145)
t370 = sin(pkin(20));
t372 = cos(pkin(20));
t374 = sin(qJ(3));
t379 = cos(qJ(3));
t459 = pkin(6) * (-t379 * t370 - t374 * t372);
t441 = pkin(1) * t459;
t354 = -0.2e1 * t441;
t387 = pkin(6) ^ 2;
t446 = t354 + t387;
t471 = -pkin(2) + pkin(13);
t472 = -pkin(2) - pkin(13);
t339 = sqrt(-((pkin(1) - t471) * (pkin(1) + t471) + t446) * ((pkin(1) - t472) * (pkin(1) + t472) + t446));
t388 = pkin(2) ^ 2;
t390 = pkin(1) ^ 2;
t444 = t387 + t390;
t433 = -pkin(13) ^ 2 + t444;
t343 = t354 + t388 + t433;
t353 = -pkin(1) + t459;
t458 = pkin(6) * (t374 * t370 - t379 * t372);
t337 = -t339 * t458 - t353 * t343;
t338 = -t353 * t339 + t343 * t458;
t380 = cos(qJ(2));
t389 = 0.1e1 / pkin(2);
t450 = 0.1e1 / (t354 + t444) * t389;
t375 = sin(qJ(2));
t462 = -t375 / 0.2e1;
t474 = (t380 * t337 / 0.2e1 + t338 * t462) * t450;
t473 = 0.1e1 / pkin(3);
t470 = -pkin(8) - pkin(3);
t469 = -pkin(8) + pkin(3);
t468 = -pkin(9) - pkin(11);
t467 = -pkin(9) + pkin(11);
t466 = g(3) * pkin(14);
t465 = -t339 / 0.2e1;
t464 = t339 / 0.2e1;
t463 = -t388 / 0.2e1 + t433 / 0.2e1 - t441;
t461 = sin(pkin(18));
t460 = pkin(1) * t375;
t377 = sin(pkin(19));
t382 = cos(pkin(19));
t361 = t375 * t382 - t380 * t377;
t457 = pkin(7) * t361;
t445 = -0.2e1 * pkin(1) * t457 + t390;
t434 = pkin(7) ^ 2 + t445;
t346 = 0.1e1 / t434;
t443 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t426 = t434 - t443;
t432 = pkin(1) - t457;
t340 = sqrt(-((pkin(7) - t469) * (pkin(7) + t469) + t445) * ((pkin(7) - t470) * (pkin(7) + t470) + t445));
t364 = t375 * t377 + t380 * t382;
t448 = t364 * t340;
t417 = t473 * (-pkin(7) * t448 + t432 * t426);
t415 = -t417 / 0.2e1;
t418 = t473 * (pkin(7) * t364 * t426 + t432 * t340);
t416 = t418 / 0.2e1;
t413 = t374 * t416 + t379 * t415;
t414 = t374 * t415 - t379 * t418 / 0.2e1;
t440 = pkin(23) + pkin(22);
t430 = sin(t440);
t431 = cos(t440);
t411 = (t430 * t413 + t431 * t414) * t346;
t410 = pkin(5) * t411;
t407 = -0.2e1 * pkin(4) * t410 + pkin(5) ^ 2;
t400 = sqrt(-((pkin(4) - t467) * (pkin(4) + t467) + t407) * ((pkin(4) - t468) * (pkin(4) + t468) + t407));
t406 = pkin(4) ^ 2 + t407;
t442 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t402 = t406 + t442;
t405 = 0.1e1 / t406;
t403 = 0.1e1 / pkin(11) * t405;
t408 = -pkin(4) * t411 + pkin(5);
t412 = t346 * (t431 * t413 - t430 * t414);
t395 = (pkin(4) * t402 * t412 + t408 * t400) * t403 / 0.2e1;
t399 = t400 * t412;
t397 = (-pkin(4) * t399 + t408 * t402) * t403;
t452 = sin(pkin(21));
t454 = cos(pkin(21));
t311 = t452 * t397 / 0.2e1 + t454 * t395;
t312 = -t454 * t397 / 0.2e1 + t452 * t395;
t362 = t380 * t374 + t375 * t379;
t363 = -t375 * t374 + t380 * t379;
t308 = t363 * t311 + t362 * t312;
t456 = g(3) * t308;
t455 = t380 * pkin(1) + pkin(14);
t453 = cos(pkin(23));
t451 = sin(pkin(23));
t449 = t346 / pkin(8);
t447 = 0.1e1 / pkin(13) * t389;
t439 = g(1) * t447;
t438 = g(2) * t447;
t437 = g(3) * t447;
t436 = cos(pkin(18)) / 0.2e1;
t435 = t362 * pkin(5) + t455;
t376 = sin(qJ(1));
t366 = t376 * pkin(16);
t429 = -t376 * t460 + t366;
t381 = cos(qJ(1));
t367 = t381 * pkin(16);
t428 = -t381 * t460 + t367;
t427 = g(1) * t381 + g(2) * t376;
t329 = (t453 * t415 + t451 * t416) * t346;
t330 = (t453 * t416 + t451 * t417 / 0.2e1) * t346;
t321 = t380 * t329 - t375 * t330;
t322 = -t375 * t329 - t380 * t330;
t349 = t363 * t376;
t425 = t349 * pkin(5) + t429;
t351 = t363 * t381;
t424 = t351 * pkin(5) + t428;
t344 = t434 + t443;
t355 = pkin(1) * t361 - pkin(7);
t423 = pkin(1) * t364 * t344 - t355 * t340;
t422 = -pkin(1) * t448 - t355 * t344;
t350 = t362 * t376;
t303 = -t349 * t311 - t350 * t312;
t352 = t362 * t381;
t305 = -t351 * t311 - t352 * t312;
t307 = -t362 * t311 + t363 * t312;
t420 = g(1) * t305 + g(2) * t303 + g(3) * t307;
t328 = (-t380 * t338 / 0.2e1 + t337 * t462) * t450;
t341 = -g(1) * t428 - g(2) * t429 - g(3) * t455;
t409 = pkin(4) - t410;
t404 = 0.1e1 / pkin(9) * t405;
t401 = t406 - t442;
t398 = (-pkin(5) * t399 + t409 * t401) * t404;
t396 = -(pkin(5) * t401 * t412 + t409 * t400) * t404 / 0.2e1;
t378 = cos(qJ(4));
t373 = sin(qJ(4));
t371 = cos(pkin(22));
t369 = sin(pkin(22));
t365 = -g(1) * t376 + g(2) * t381;
t358 = -t427 * pkin(16) - t466;
t332 = (t423 * t436 + t422 * t461 / 0.2e1) * t449;
t331 = (t422 * t436 - t461 * t423 / 0.2e1) * t449;
t326 = t381 * t474;
t325 = t381 * t328;
t324 = t376 * t474;
t323 = t376 * t328;
t320 = t321 * t381;
t319 = t322 * t381;
t318 = t321 * t376;
t317 = t322 * t376;
t310 = -t371 * t398 / 0.2e1 + t369 * t396;
t309 = t369 * t398 / 0.2e1 + t371 * t396;
t306 = -t352 * t311 + t351 * t312;
t304 = -t350 * t311 + t349 * t312;
t1 = [0, 0, 0, 0, 0, 0, -t427, -t365, -g(3), -t466, 0, 0, 0, 0, 0, 0, -g(3) * t380 + t427 * t375, g(3) * t375 + t427 * t380, t365, t358, 0, 0, 0, 0, 0, 0, -g(1) * t351 - g(2) * t349 - g(3) * t362, g(1) * t352 + g(2) * t350 - g(3) * t363, t365, t341, 0, 0, 0, 0, 0, 0, -g(1) * t306 - g(2) * t304 - t456, -t420, t365, -g(1) * t424 - g(2) * t425 - g(3) * t435, 0, 0, 0, 0, 0, 0, -t378 * t456 - g(2) * (t304 * t378 - t381 * t373) - g(1) * (t306 * t378 + t376 * t373), t373 * t456 - g(2) * (-t304 * t373 - t381 * t378) - g(1) * (-t306 * t373 + t376 * t378), t420, -g(3) * (t308 * pkin(10) - t307 * pkin(12) + t435) - g(2) * (t304 * pkin(10) - t303 * pkin(12) + t425) - g(1) * (t306 * pkin(10) - t305 * pkin(12) + t424), 0, 0, 0, 0, 0, 0, -g(3) * t332 - t427 * t331, -g(3) * t331 + t427 * t332, t365, -g(3) * (-pkin(17) + pkin(14)) + t427 * pkin(15), 0, 0, 0, 0, 0, 0, -g(1) * t319 - g(2) * t317 - g(3) * t321, g(1) * t320 + g(2) * t318 - g(3) * t322, t365, t341, 0, 0, 0, 0, 0, 0, -g(1) * t325 - g(2) * t323 - g(3) * t474, g(1) * t326 + g(2) * t324 - g(3) * t328, t365, t358, 0, 0, 0, 0, 0, 0, -(t328 * t465 + t463 * t474) * t437 - (t323 * t463 - t324 * t465) * t438 - (t325 * t463 - t326 * t465) * t439, -(t328 * t463 + t464 * t474) * t437 - (t323 * t464 - t324 * t463) * t438 - (t325 * t464 - t326 * t463) * t439, t365, -g(3) * (pkin(2) * t474 + pkin(14)) - g(2) * (t323 * pkin(2) + t366) - g(1) * (t325 * pkin(2) + t367), 0, 0, 0, 0, 0, 0, -g(3) * (t322 * t309 + t321 * t310) - g(2) * (-t318 * t309 + t317 * t310) - g(1) * (-t320 * t309 + t319 * t310), -g(3) * (-t321 * t309 + t322 * t310) - g(2) * (-t317 * t309 - t318 * t310) - g(1) * (-t319 * t309 - t320 * t310), t365, (-g(3) * (t321 * t371 - t322 * t369) - g(2) * (t317 * t371 + t318 * t369) - g(1) * (t319 * t371 + t320 * t369)) * pkin(4) + t341;];
U_reg = t1;
