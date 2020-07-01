% Calculate vector of inverse dynamics joint torques for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1OL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1OL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'palh2m1OL_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:21
% EndTime: 2020-06-30 17:45:31
% DurationCPUTime: 4.78s
% Computational Cost: add. (3686->369), mult. (8527->489), div. (0->0), fcn. (7025->14), ass. (0->169)
t404 = qJ(2) + qJ(3);
t398 = qJ(4) + t404;
t385 = sin(t398);
t386 = cos(t398);
t411 = sin(qJ(1));
t416 = cos(qJ(1));
t448 = g(1) * t416 + g(2) * t411;
t519 = g(3) * t386 + t448 * t385;
t412 = cos(qJ(5));
t480 = qJD(5) * t412;
t414 = cos(qJ(3));
t409 = sin(qJ(3));
t410 = sin(qJ(2));
t487 = qJD(1) * t410;
t468 = t409 * t487;
t415 = cos(qJ(2));
t486 = qJD(1) * t415;
t354 = -t414 * t486 + t468;
t437 = t409 * t415 + t410 * t414;
t355 = t437 * qJD(1);
t408 = sin(qJ(4));
t413 = cos(qJ(4));
t320 = t413 * t354 + t355 * t408;
t528 = t320 * t412;
t536 = t480 - t528;
t368 = qJD(1) * pkin(1) + pkin(2) * t486;
t334 = -pkin(3) * t354 + t368;
t513 = pkin(2) * t414;
t391 = qJDD(2) * t513;
t400 = qJDD(2) + qJDD(3);
t508 = pkin(2) * qJD(2);
t472 = t409 * t508;
t341 = pkin(3) * t400 - qJD(3) * t472 + t391;
t452 = t408 * t472;
t461 = -qJD(4) * t452 + t408 * t341;
t509 = g(3) * t385;
t534 = -t334 * t320 + t448 * t386 - t461 - t509;
t401 = qJD(2) + qJD(3);
t360 = pkin(3) * t401 + t414 * t508;
t485 = qJD(3) * t414;
t465 = qJD(2) * t485;
t474 = qJDD(2) * t409;
t426 = (t465 + t474) * pkin(2);
t533 = (qJD(4) * t360 + t426) * t413;
t532 = -t533 + t534;
t441 = t354 * t408 - t413 * t355;
t531 = -t334 * t441 + t519;
t476 = qJDD(1) * t410;
t477 = qJD(1) * qJD(2);
t466 = t410 * t477;
t475 = qJDD(1) * t415;
t523 = -t475 + t466;
t309 = qJD(3) * t468 + (-t401 * t486 - t476) * t414 + t523 * t409;
t331 = t401 * t437;
t356 = t409 * t410 - t414 * t415;
t310 = t331 * qJD(1) + t356 * qJDD(1);
t482 = qJD(4) * t413;
t483 = qJD(4) * t408;
t287 = t413 * t309 + t408 * t310 + t354 * t482 + t355 * t483;
t394 = qJDD(4) + t400;
t395 = qJD(4) + t401;
t407 = sin(qJ(5));
t469 = t412 * t287 + t407 * t394 + t395 * t480;
t481 = qJD(5) * t407;
t278 = -t441 * t481 + t469;
t276 = t278 * t407;
t277 = t278 * t412;
t313 = t395 * t407 + t412 * t441;
t377 = t412 * t394;
t279 = t313 * qJD(5) + t287 * t407 - t377;
t288 = t441 * qJD(4) + t309 * t408 - t413 * t310;
t286 = qJDD(5) + t288;
t284 = t412 * t286;
t311 = -t412 * t395 + t407 * t441;
t478 = -qJD(5) + t320;
t283 = t407 * t286;
t492 = -t478 * t480 + t283;
t526 = t478 * t407;
t530 = t394 * MDP(22) - t288 * MDP(21) - t320 ^ 2 * MDP(19) + (-t320 * t395 + t287) * MDP(20) + (-t320 * MDP(18) + MDP(19) * t441 + t395 * MDP(21) + MDP(29) * t478) * t441 + (t536 * t313 + t276) * MDP(25) + (-t313 * t441 + t478 * t528 + t492) * MDP(27) + (t311 * t441 - t478 * t526 + t284) * MDP(28) + (-t407 * t279 - t536 * t311 + t313 * t526 + t277) * MDP(26);
t467 = t409 * t482;
t529 = pkin(2) * ((t408 * t485 + t467) * qJD(2) + t408 * t474);
t335 = t360 * t413 - t452;
t332 = -pkin(4) * t395 - t335;
t527 = t332 * t320;
t338 = t413 * t341;
t445 = t360 * t483 - t338;
t524 = -t445 + t531;
t297 = pkin(4) * t441 - pkin(6) * t320;
t291 = -pkin(4) * t320 - pkin(6) * t441 + t334;
t336 = t360 * t408 + t413 * t472;
t333 = pkin(6) * t395 + t336;
t289 = t291 * t412 - t333 * t407;
t300 = -pkin(4) * t394 + t445 + t529;
t423 = -t289 * t441 + t332 * t481 + (-t300 + t519) * t412;
t290 = t291 * t407 + t333 * t412;
t446 = t290 * t441 + t300 * t407 + t332 * t480;
t330 = t401 * t356;
t440 = t413 * t356 + t408 * t437;
t293 = t440 * qJD(4) + t330 * t413 + t331 * t408;
t329 = t356 * t408 - t413 * t437;
t390 = t415 * pkin(2) + pkin(1);
t340 = -pkin(3) * t356 + t390;
t296 = -pkin(4) * t440 - pkin(6) * t329 + t340;
t299 = pkin(6) * t394 + t461 + t533;
t457 = qJD(5) * t291 + t299;
t516 = qJD(5) * t296 * t478 + t332 * t293 + t300 * t329 + t440 * t457;
t512 = pkin(3) * t355;
t405 = qJDD(1) * pkin(1);
t504 = t329 * t332;
t500 = t394 * t413;
t499 = t407 * t411;
t498 = t407 * t416;
t497 = t408 * t409;
t496 = t409 * t413;
t495 = t411 * t412;
t494 = t412 * t416;
t389 = pkin(3) + t513;
t490 = pkin(2) * t496 + t408 * t389;
t402 = t410 ^ 2;
t489 = -t415 ^ 2 + t402;
t473 = pkin(2) * t487;
t471 = t410 * t508;
t463 = -0.2e1 * pkin(1) * t477;
t352 = -t523 * pkin(2) + t405;
t303 = -pkin(3) * t310 + t352;
t272 = pkin(4) * t288 - pkin(6) * t287 + t303;
t458 = qJD(5) * t333 - t272;
t295 = t297 - t512;
t349 = pkin(6) + t490;
t456 = qJD(5) * t349 + t295 - t473;
t387 = pkin(3) * t408 + pkin(6);
t455 = qJD(5) * t387 + t295;
t454 = qJD(2) * (-qJD(3) + t401);
t453 = qJD(3) * (-qJD(2) - t401);
t451 = qJD(2) * t467;
t439 = t408 * t414 + t496;
t350 = t439 * t508;
t450 = pkin(3) * t483 - t350;
t438 = t413 * t414 - t497;
t351 = t438 * t508;
t449 = -pkin(3) * t482 + t351;
t447 = g(1) * t411 - g(2) * t416;
t294 = t329 * qJD(4) + t330 * t408 - t413 * t331;
t322 = -pkin(3) * t331 - t471;
t444 = -(pkin(4) * t294 - pkin(6) * t293 + t322) * t478 + t296 * t286;
t443 = -t286 * t349 - t527;
t442 = -t387 * t286 - t527;
t436 = t293 * t412 - t329 * t481;
t396 = sin(t404);
t397 = cos(t404);
t434 = g(3) * t397 + t355 * t368 + t448 * t396 + t391;
t418 = qJD(1) ^ 2;
t432 = pkin(1) * t418 + t448;
t431 = -g(3) * t396 - t354 * t368 + t448 * t397;
t430 = t447 + 0.2e1 * t405;
t421 = (-pkin(3) * t395 - t360) * qJD(4) - t426;
t419 = t355 * t354 * MDP(11) + (-t354 * t401 + t309) * MDP(13) + (-t355 * t401 + t310) * MDP(14) + (-t354 ^ 2 + t355 ^ 2) * MDP(12) + t400 * MDP(15) + t530;
t417 = qJD(2) ^ 2;
t388 = -pkin(3) * t413 - pkin(4);
t348 = pkin(2) * t497 - t389 * t413 - pkin(4);
t347 = t386 * t494 - t499;
t346 = -t386 * t498 - t495;
t345 = -t386 * t495 - t498;
t344 = t386 * t499 - t494;
t337 = -t473 - t512;
t324 = t389 * t483 + (t439 * qJD(3) + t467) * pkin(2);
t323 = t389 * t482 + (t438 * qJD(3) - t409 * t483) * pkin(2);
t269 = t412 * t272;
t1 = [(t287 * t340 + t293 * t334 + t303 * t329 + t322 * t441 - t447 * t385) * MDP(24) + (t287 * t329 + t293 * t441) * MDP(18) + (-t430 * t410 + t415 * t463) * MDP(10) + (t410 * t463 + t430 * t415) * MDP(9) + t448 * MDP(3) + t447 * MDP(2) + (-g(1) * t344 - g(2) * t346 - t290 * t294 + t516 * t412 + (-qJD(5) * t504 - t440 * t458 - t444) * t407) * MDP(31) + (-g(1) * t345 - g(2) * t347 - t269 * t440 + t289 * t294 + ((t333 * t440 + t504) * qJD(5) + t444) * t412 + t516 * t407) * MDP(30) + (-t294 * t395 + t394 * t440) * MDP(21) + (-t278 * t440 + t329 * t284 + t294 * t313 - t436 * t478) * MDP(27) + (-t329 * t283 + t279 * t440 - t294 * t311 - (-t293 * t407 - t329 * t480) * t478) * MDP(28) + (-t286 * t440 - t294 * t478) * MDP(29) + (t287 * t440 - t288 * t329 + t293 * t320 - t294 * t441) * MDP(19) + (t288 * t340 + t294 * t334 - t303 * t440 - t320 * t322 + t447 * t386) * MDP(23) + (t309 * t390 + t330 * t368 - t352 * t437 + t355 * t471 - t447 * t396) * MDP(17) + (t330 * t401 - t400 * t437) * MDP(13) + (-t309 * t437 - t330 * t355) * MDP(11) + (t309 * t356 - t310 * t437 + t330 * t354 - t331 * t355) * MDP(12) + (-t310 * t390 - t331 * t368 - t352 * t356 + t354 * t471 + t447 * t397) * MDP(16) + 0.2e1 * (t410 * t475 - t489 * t477) * MDP(5) + (t293 * t395 + t329 * t394) * MDP(20) + (t331 * t401 + t356 * t400) * MDP(14) + qJDD(1) * MDP(1) + (qJDD(1) * t402 + 0.2e1 * t415 * t466) * MDP(4) + (-qJDD(2) * t410 - t415 * t417) * MDP(6) + (-qJDD(2) * t415 + t410 * t417) * MDP(7) + (t329 * t277 + t436 * t313) * MDP(25) + ((-t311 * t412 - t313 * t407) * t293 + (-t276 - t279 * t412 + (t311 * t407 - t313 * t412) * qJD(5)) * t329) * MDP(26); (t279 * t348 + t311 * t324 + t443 * t407 - (-t323 * t407 - t456 * t412) * t478 + t423) * MDP(30) + t419 - MDP(6) * t476 + (t278 * t348 + t313 * t324 + (t323 * t478 + t443) * t412 + (-t456 * t478 - t519) * t407 + t446) * MDP(31) + ((-t354 * t487 + t400 * t414 + t409 * t453) * pkin(2) + t434) * MDP(16) + (t389 * t500 + t337 * t320 - t324 * t395 + (-t451 + (-t465 + (-qJDD(2) - t394) * t409) * t408) * pkin(2) + t524) * MDP(23) + (g(3) * t415 + t432 * t410) * MDP(9) + (-g(3) * t410 + t432 * t415) * MDP(10) + (-t323 * t395 - t337 * t441 - t490 * t394 + t532) * MDP(24) + qJDD(2) * MDP(8) - MDP(7) * t475 + ((-t355 * t487 + (-qJDD(2) - t400) * t409 + t414 * t453) * pkin(2) + t431) * MDP(17) + (-t410 * t415 * MDP(4) + t489 * MDP(5)) * t418; t419 + (t351 * t395 + (t355 * t441 - t394 * t408) * pkin(3) + t421 * t413 + t534) * MDP(24) + (t388 * t279 + t442 * t407 + t450 * t311 - (t449 * t407 - t455 * t412) * t478 + t423) * MDP(30) + (t388 * t278 + t450 * t313 + (-t449 * t478 + t442) * t412 + (-t455 * t478 - t519) * t407 + t446) * MDP(31) + (-pkin(2) * t451 + t350 * t395 + t338 + (-t320 * t355 + t500) * pkin(3) + t421 * t408 + t531) * MDP(23) + (t409 * pkin(2) * t454 + t434) * MDP(16) + ((t414 * t454 - t474) * pkin(2) + t431) * MDP(17); (t336 * t395 + t524 - t529) * MDP(23) + (t335 * t395 + t532) * MDP(24) + (-pkin(4) * t279 + (t297 * t412 - t335 * t407) * t478 - t336 * t311 - t407 * t527 - t492 * pkin(6) + t423) * MDP(30) + (-pkin(4) * t278 - t336 * t313 + (-pkin(6) * t286 - t335 * t478 - t527) * t412 + (-(pkin(6) * qJD(5) + t297) * t478 - t519) * t407 + t446) * MDP(31) + t530; t313 * t311 * MDP(25) + (-t311 ^ 2 + t313 ^ 2) * MDP(26) + (-t311 * t478 + t469) * MDP(27) + (-t313 * t478 + t377) * MDP(28) + t286 * MDP(29) + (-g(1) * t346 + g(2) * t344 - t290 * t478 - t313 * t332 + t269) * MDP(30) + (g(1) * t347 - g(2) * t345 - t289 * t478 + t311 * t332) * MDP(31) + ((-t299 - t509) * MDP(31) + (-MDP(28) * t441 - MDP(30) * t333 - MDP(31) * t291) * qJD(5)) * t412 + (-qJD(5) * t441 * MDP(27) + (-qJD(5) * t395 - t287) * MDP(28) + (-t457 - t509) * MDP(30) + t458 * MDP(31)) * t407;];
tau = t1;
