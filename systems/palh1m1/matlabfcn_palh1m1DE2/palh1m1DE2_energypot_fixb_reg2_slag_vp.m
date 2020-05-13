% Calculate inertial parameters regressor of potential energy for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m1DE2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energypot_fixb_reg2_slag_vp: pkin has to be [23x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 06:37:46
% EndTime: 2020-04-15 06:37:51
% DurationCPUTime: 4.17s
% Computational Cost: add. (70597->131), mult. (105916->191), div. (5052->9), fcn. (67021->48), ass. (0->122)
t382 = sin(qJ(1));
t388 = cos(qJ(1));
t456 = -g(1) * t388 - g(2) * t382;
t455 = -pkin(2) - pkin(13);
t454 = -pkin(2) + pkin(13);
t453 = -pkin(8) - pkin(3);
t452 = -pkin(8) + pkin(3);
t451 = -pkin(9) - pkin(11);
t450 = -pkin(9) + pkin(11);
t449 = g(3) * pkin(14);
t407 = pkin(1) ^ 2;
t381 = sin(qJ(2));
t383 = sin(pkin(19));
t387 = cos(qJ(2));
t389 = cos(pkin(19));
t360 = t381 * t389 - t387 * t383;
t438 = pkin(7) * t360;
t422 = -0.2e1 * pkin(1) * t438 + t407;
t415 = pkin(7) ^ 2 + t422;
t419 = pkin(3) ^ 2 - pkin(8) ^ 2;
t349 = t415 + t419;
t354 = pkin(1) - t438;
t345 = sqrt(-((pkin(7) - t452) * (pkin(7) + t452) + t422) * ((pkin(7) - t453) * (pkin(7) + t453) + t422));
t361 = t381 * t383 + t387 * t389;
t431 = t345 * t361;
t341 = -pkin(7) * t431 + t354 * t349;
t448 = -t341 / 0.2e1;
t342 = pkin(7) * t361 * t349 + t354 * t345;
t447 = t342 / 0.2e1;
t446 = sin(pkin(23)) / 0.2e1;
t445 = sin(pkin(21)) / 0.2e1;
t386 = cos(qJ(3));
t444 = -t386 / 0.2e1;
t443 = cos(pkin(18)) / 0.2e1;
t442 = 0.1e1 / pkin(2) / 0.2e1;
t380 = sin(qJ(3));
t351 = 0.1e1 / t415;
t429 = t351 / pkin(3);
t338 = (t342 * t444 + t380 * t448) * t429;
t339 = (t341 * t444 + t380 * t447) * t429;
t371 = pkin(23) + pkin(22);
t366 = sin(t371);
t367 = cos(t371);
t320 = t367 * t338 + t366 * t339;
t441 = pkin(5) * t320;
t375 = sin(pkin(20));
t378 = cos(pkin(20));
t440 = pkin(6) * (-t386 * t375 - t380 * t378);
t439 = pkin(6) * (t380 * t375 - t386 * t378);
t424 = -0.2e1 * pkin(4) * t441 + pkin(5) ^ 2;
t416 = pkin(4) ^ 2 + t424;
t421 = pkin(9) ^ 2 - pkin(11) ^ 2;
t314 = t416 - t421;
t317 = -pkin(4) * t320 + pkin(5);
t312 = sqrt(-((pkin(4) - t450) * (pkin(4) + t450) + t424) * ((pkin(4) - t451) * (pkin(4) + t451) + t424));
t321 = -t366 * t338 + t367 * t339;
t433 = t312 * t321;
t310 = -pkin(4) * t433 + t317 * t314;
t311 = pkin(4) * t321 * t314 + t317 * t312;
t372 = qJ(2) + qJ(3);
t377 = cos(pkin(21));
t315 = 0.1e1 / t416;
t432 = t315 / pkin(11);
t306 = atan2((t310 * t445 + t311 * t377 / 0.2e1) * t432, (-t310 * t377 / 0.2e1 + t311 * t445) * t432) + t372;
t304 = sin(t306);
t435 = g(3) * t304;
t434 = t387 * pkin(1) + pkin(14);
t430 = t351 / pkin(8);
t379 = sin(qJ(4));
t428 = t382 * t379;
t385 = cos(qJ(4));
t427 = t382 * t385;
t426 = t388 * t379;
t425 = t388 * t385;
t376 = cos(pkin(23));
t325 = qJ(2) + atan2((t341 * t446 + t376 * t447) * t429, (t342 * t446 + t376 * t448) * t429);
t418 = pkin(1) * t440;
t353 = -0.2e1 * t418;
t400 = pkin(6) ^ 2;
t423 = t353 + t400;
t344 = sqrt(-((pkin(1) - t454) * (pkin(1) + t454) + t423) * ((pkin(1) - t455) * (pkin(1) + t455) + t423));
t405 = pkin(2) ^ 2;
t420 = t400 + t407;
t414 = -pkin(13) ^ 2 + t420;
t347 = t353 + t405 + t414;
t352 = -pkin(1) + t440;
t412 = 0.1e1 / (t353 + t420) * t442;
t336 = qJ(2) + atan2((-t352 * t344 + t347 * t439) * t412, (-t344 * t439 - t352 * t347) * t412);
t368 = sin(t372);
t417 = pkin(5) * t368 + t434;
t413 = t315 / pkin(9) / 0.2e1;
t411 = 0.1e1 / pkin(13) * t442;
t410 = -t381 * pkin(1) + pkin(16);
t322 = pkin(22) - t325;
t384 = sin(pkin(18));
t369 = cos(t372);
t363 = -g(1) * t382 + g(2) * t388;
t362 = pkin(5) * t369 + t410;
t357 = pkin(16) * t456 - t449;
t355 = pkin(1) * t360 - pkin(7);
t348 = t415 - t419;
t346 = -g(3) * t434 + t410 * t456;
t343 = pkin(1) * t361 * t348 - t355 * t345;
t340 = -pkin(1) * t431 - t355 * t348;
t335 = cos(t336);
t334 = sin(t336);
t332 = atan2(t344 * t411, (t405 - t414 + 0.2e1 * t418) * t411) + t336;
t331 = cos(t332);
t330 = sin(t332);
t329 = atan2((t343 * t443 + t340 * t384 / 0.2e1) * t430, (t340 * t443 - t384 * t343 / 0.2e1) * t430);
t327 = cos(t329);
t326 = sin(t329);
t324 = cos(t325);
t323 = sin(t325);
t316 = -pkin(4) + t441;
t313 = t416 + t421;
t309 = -atan2((pkin(5) * t321 * t313 - t316 * t312) * t413, (-pkin(5) * t433 - t316 * t313) * t413) + t322;
t308 = cos(t309);
t307 = sin(t309);
t305 = cos(t306);
t303 = g(3) * t305 + t456 * t304;
t1 = [0, 0, 0, 0, 0, 0, t456, -t363, -g(3), -t449, 0, 0, 0, 0, 0, 0, -g(3) * t387 - t381 * t456, g(3) * t381 - t387 * t456, t363, t357, 0, 0, 0, 0, 0, 0, -g(3) * t368 + t369 * t456, -g(3) * t369 - t368 * t456, t363, t346, 0, 0, 0, 0, 0, 0, t456 * t305 - t435, -t303, t363, -g(3) * t417 + t362 * t456, 0, 0, 0, 0, 0, 0, -t385 * t435 - g(2) * (t305 * t427 - t426) - g(1) * (t305 * t425 + t428), t379 * t435 - g(2) * (-t305 * t428 - t425) - g(1) * (-t305 * t426 + t427), t303, -g(3) * (t304 * pkin(10) - t305 * pkin(12) + t417) + t456 * (pkin(10) * t305 + pkin(12) * t304 + t362), 0, 0, 0, 0, 0, 0, -g(3) * t326 + t327 * t456, -g(3) * t327 - t326 * t456, t363, -g(3) * (-pkin(17) + pkin(14)) - t456 * pkin(15), 0, 0, 0, 0, 0, 0, -g(3) * t324 - t323 * t456, g(3) * t323 - t324 * t456, t363, t346, 0, 0, 0, 0, 0, 0, -g(3) * t335 - t334 * t456, g(3) * t334 - t335 * t456, t363, t357, 0, 0, 0, 0, 0, 0, g(3) * t331 + t330 * t456, -g(3) * t330 + t331 * t456, t363, -g(3) * (pkin(2) * t335 + pkin(14)) + t456 * (-pkin(2) * t334 + pkin(16)), 0, 0, 0, 0, 0, 0, g(3) * t308 - t307 * t456, g(3) * t307 + t308 * t456, t363, -g(3) * (pkin(4) * cos(t322) + t434) + t456 * (pkin(4) * sin(t322) + t410);];
U_reg = t1;
