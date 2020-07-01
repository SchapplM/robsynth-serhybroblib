% Calculate minimal parameter regressor of Coriolis joint torque vector for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnTE_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnTE_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:23:01
% EndTime: 2020-06-27 16:23:19
% DurationCPUTime: 6.43s
% Computational Cost: add. (67309->278), mult. (94967->649), div. (2796->15), fcn. (25412->4), ass. (0->244)
t390 = 0.1e1 / pkin(4);
t394 = pkin(2) ^ 2;
t395 = pkin(1) ^ 2;
t385 = cos(qJ(2));
t546 = pkin(2) * t385;
t503 = -0.2e1 * pkin(1) * t546 + t395;
t378 = t394 + t503;
t579 = pkin(4) ^ 2;
t501 = pkin(3) ^ 2 - t579;
t374 = t378 - t501;
t379 = pkin(1) - t546;
t384 = sin(qJ(2));
t565 = -pkin(3) - pkin(4);
t371 = (pkin(2) - t565) * (pkin(2) + t565) + t503;
t564 = pkin(4) - pkin(3);
t372 = (pkin(2) - t564) * (pkin(2) + t564) + t503;
t524 = t371 * t372;
t396 = sqrt(-t524);
t516 = t384 * t396;
t354 = -pkin(2) * t516 + t374 * t379;
t376 = 0.1e1 / t378 ^ 2;
t497 = qJD(2) * t384;
t480 = pkin(2) * t497;
t449 = t376 * t480;
t431 = pkin(1) * t449;
t446 = pkin(1) * pkin(2) * (-t371 - t372);
t361 = t384 * t446;
t360 = qJD(2) * t361;
t365 = 0.1e1 / t396;
t529 = t360 * t365;
t474 = t384 * t529;
t493 = 0.2e1 * t379 * pkin(1);
t512 = t385 * t396;
t324 = (-t474 + (-t512 + (t374 + t493) * t384) * qJD(2)) * pkin(2);
t375 = 0.1e1 / t378;
t540 = t324 * t375;
t407 = -t540 / 0.2e1 + t354 * t431;
t315 = t407 * t390;
t547 = pkin(2) * t384;
t369 = t374 * t547;
t355 = t379 * t396 + t369;
t457 = pkin(1) * t480;
t430 = t355 * t457;
t382 = t384 ^ 2;
t519 = t382 * t394;
t482 = pkin(1) * t519;
t448 = qJD(2) * t482;
t467 = t396 * t497;
t496 = qJD(2) * t385;
t505 = (t374 * t496 + t467) * pkin(2);
t527 = t365 * t379;
t325 = t360 * t527 + 0.2e1 * t448 + t505;
t539 = t325 * t375;
t406 = t539 / 0.2e1 - t376 * t430;
t317 = t406 * t390;
t348 = 0.1e1 / t354;
t349 = 0.1e1 / t354 ^ 2;
t533 = t349 * t355;
t418 = t315 * t533 + t317 * t348;
t578 = -0.4e1 * t354;
t577 = 0.6e1 * t385;
t350 = t348 * t349;
t351 = t355 ^ 2;
t532 = t350 * t351;
t309 = -t324 * t532 + t325 * t533;
t339 = t349 * t351 + 0.1e1;
t336 = 0.1e1 / t339 ^ 2;
t576 = t309 * t336;
t575 = t360 * t361;
t335 = 0.1e1 / t339;
t544 = pkin(4) * t378;
t484 = t335 * t544;
t452 = t349 * t484;
t574 = t390 * t355 * t452;
t373 = t378 + t501;
t380 = pkin(1) * t385 - pkin(2);
t566 = -0.2e1 * t380;
t494 = pkin(2) * t566;
t460 = t373 + t494;
t323 = (-t474 + (t384 * t460 - t512) * qJD(2)) * pkin(1);
t353 = -pkin(1) * t516 - t373 * t380;
t393 = 0.1e1 / pkin(3);
t424 = -0.2e1 * t431;
t316 = (t323 * t375 + t353 * t424) * t393;
t346 = 0.1e1 / t353 ^ 2;
t550 = pkin(1) * t384;
t370 = t373 * t550;
t356 = -t380 * t396 + t370;
t535 = t346 * t356;
t518 = t382 * t395;
t481 = pkin(2) * t518;
t447 = qJD(2) * t481;
t504 = (t373 * t496 + t467) * pkin(1);
t526 = t365 * t380;
t326 = -t360 * t526 + 0.2e1 * t447 + t504;
t318 = (t326 * t375 + t356 * t424) * t393;
t345 = 0.1e1 / t353;
t541 = t318 * t345;
t573 = -t316 * t535 + t541;
t572 = -t384 * t385 * MDP(4) + (-t385 ^ 2 + t382) * MDP(5);
t528 = t365 * t361;
t328 = t369 + (-t512 + (t493 - t528) * t384) * pkin(2);
t332 = t361 * t527 + 0.2e1 * t482 + (t374 * t385 + t516) * pkin(2);
t567 = 0.2e1 * t378;
t441 = 0.2e1 * t336 * t567;
t491 = pkin(1) * t547;
t469 = -0.4e1 * t491;
t443 = t348 * t469;
t530 = t355 * t378;
t456 = 0.4e1 * t350 * t530;
t488 = t349 * t567;
t472 = t394 * t518;
t408 = t385 * t446 - 0.4e1 * t472;
t358 = t408 * qJD(2);
t525 = 0.4e1 * t365 / t524;
t465 = -t525 / 0.4e1;
t425 = t465 * t575;
t554 = -t384 / 0.2e1;
t404 = (t384 * t425 + 0.2e1 * (t358 * t554 - t360 * t385) * t365) * t375;
t377 = t375 * t376;
t444 = t377 * t472;
t421 = 0.4e1 * t355 * t444;
t426 = qJD(2) * t444;
t483 = t348 * t544;
t428 = t335 * t390 * t483;
t487 = t394 * t550;
t436 = t487 * t577;
t440 = t379 * t525 / 0.4e1;
t522 = t375 * t385;
t473 = t379 * t522;
t475 = t375 * t529;
t520 = t376 * t384;
t538 = t325 * t376;
t551 = pkin(1) * t376;
t555 = t375 / 0.2e1;
t556 = -t375 / 0.2e1;
t506 = -0.2e1 * ((0.4e1 * t448 + t505) * t556 + t426 * t578 + (-t404 / 0.2e1 + (t324 * t520 + (-t473 + (t328 * t384 + t354 * t385) * t376) * qJD(2)) * pkin(1)) * pkin(2)) * t574 - 0.2e1 * ((qJD(2) * t436 + t358 * t527 + t440 * t575) * t555 + qJD(2) * t421 + ((t475 / 0.2e1 - pkin(1) * t538) * t384 + ((t512 + (-t374 + t528) * t384) * t555 + (-t332 * t384 - t355 * t385) * t551) * qJD(2)) * pkin(2)) * t428;
t568 = -0.2e1 * t378;
t289 = (t418 * (-t328 * t532 + t332 * t533) * t441 + ((t328 * t488 + t443) * t317 + (t328 * t456 + (t332 * t568 + t355 * t469) * t349) * t315) * t335) * pkin(4) + t506;
t459 = t376 * t491;
t537 = t328 * t375;
t319 = (-t537 / 0.2e1 + t354 * t459) * t390;
t466 = t332 * t555;
t321 = (-t355 * t459 + t466) * t390;
t416 = t319 * t533 + t321 * t348;
t291 = (t416 * t309 * t441 + ((qJD(2) * t443 + t324 * t488) * t321 + (t324 * t456 + (t325 * t568 - 0.4e1 * t430) * t349) * t319) * t335) * pkin(4) + t506;
t438 = 0.2e1 * t484;
t300 = t418 * t438;
t301 = t416 * t438;
t571 = (t291 / 0.2e1 - t289 / 0.2e1) * t375 - (-qJD(2) * t301 + t300) * t459;
t570 = 0.2e1 * t356;
t569 = -0.2e1 * t375;
t563 = -t323 / 0.2e1;
t562 = -t324 / 0.2e1;
t561 = t325 / 0.2e1;
t560 = t326 / 0.2e1;
t331 = -t361 * t526 + 0.2e1 * t481 + (t373 * t385 + t516) * pkin(1);
t559 = -t331 / 0.2e1;
t558 = -t353 / 0.2e1;
t557 = t356 / 0.2e1;
t553 = t384 / 0.2e1;
t552 = -t385 / 0.2e1;
t549 = pkin(1) * t393;
t548 = pkin(2) * t376;
t545 = pkin(3) * t378;
t420 = 0.8e1 * t426;
t486 = t395 * t547;
t435 = t486 * t577;
t490 = 0.2e1 * t376;
t495 = -0.2e1 * t548;
t514 = t385 * t356;
t543 = ((qJD(2) * t435 - t358 * t526 + t380 * t425) * t375 + t356 * t420 + ((t326 * t495 + t475) * t384 + ((t512 + (-t373 + t528) * t384) * t375 + (-t331 * t384 - t514) * pkin(2) * t490) * qJD(2)) * pkin(1)) * t393 * t345;
t352 = t356 ^ 2;
t340 = t346 * t352 + 0.1e1;
t338 = 0.1e1 / t340 ^ 2;
t536 = t338 * t378;
t347 = t345 * t346;
t534 = t347 * t352;
t531 = t355 * t376;
t387 = qJD(2) ^ 2;
t357 = t408 * t387;
t359 = t360 ^ 2;
t439 = t359 * t465;
t508 = t396 * t387;
t405 = -t357 * t365 + t439 + t508;
t468 = qJD(2) * t529;
t451 = 0.2e1 * t468;
t413 = -t374 * t387 + t451;
t470 = t385 * t508;
t513 = t385 * t387;
t293 = -0.2e1 * ((t357 * t527 + t359 * t440 + t387 * t436) * t555 + t387 * t421 + ((t384 * t413 + t470) * t555 + (-0.2e1 * t325 * t497 - t355 * t513) * t551) * pkin(2)) * t428 - 0.2e1 * ((t377 * t395 * t578 + pkin(1) * t569) * t387 * t519 + ((t384 * t405 - t385 * t413) * t556 + (-t387 * t473 + (0.2e1 * t324 * t497 + t354 * t513) * t376) * pkin(1)) * pkin(2)) * t574 - 0.2e1 * t315 * t325 * t452 + 0.2e1 * (t324 * t452 + 0.2e1 * t483 * t576) * t317 + 0.2e1 * (0.2e1 * (t324 * t335 * t350 + t349 * t576) * t315 * t530 - 0.2e1 * t418 * t335 * t457) * pkin(4);
t523 = t375 * t293;
t521 = t375 * t393;
t517 = t384 * t353;
t515 = t385 * t353;
t511 = t387 * t393;
t327 = t370 + (-t512 + (t494 - t528) * t384) * pkin(1);
t337 = 0.1e1 / t340;
t485 = t337 * t545;
t453 = t346 * t485;
t429 = t356 * t453;
t296 = (((0.4e1 * t447 + t504) * t375 + t353 * t420) * t393 + (t404 + (-0.2e1 * t323 * t520 + (t522 * t566 + (-t327 * t384 - t515) * t490) * qJD(2)) * pkin(2)) * t549) * t429;
t445 = 0.2e1 * t491;
t489 = t347 * t570;
t290 = -t296 + (-0.2e1 * t573 * (-t327 * t534 + t331 * t535) * t536 + (t573 * t445 + (-t318 * t327 * t346 + t543 + (t327 * t489 - t331 * t346) * t316) * t378) * t337) * pkin(3);
t442 = -0.2e1 * t459;
t320 = (t327 * t375 + t353 * t442) * t393;
t322 = (t331 * t375 + t356 * t442) * t393;
t415 = t320 * t535 - t322 * t345;
t310 = -t323 * t534 + t326 * t535;
t478 = t310 * t536;
t292 = -t296 + (0.2e1 * t415 * t478 + (-t415 * qJD(2) * t445 + (-t322 * t323 * t346 + t543 + (t323 * t489 - t326 * t346) * t320) * t378) * t337) * pkin(3);
t507 = -t290 + t292;
t500 = qJD(1) * t390;
t498 = qJD(1) * t393;
t492 = pkin(1) * t548;
t463 = t292 / 0.2e1 - t290 / 0.2e1;
t462 = -t327 / 0.2e1 + t557;
t461 = t558 + t559;
t458 = t377 * t491;
t455 = t376 * t486;
t454 = t345 * t485;
t432 = t351 * t458;
t422 = 0.2e1 * t395 * t449;
t414 = t514 / 0.2e1 + t517 / 0.2e1;
t302 = -t320 * t429 + t322 * t454 + 0.1e1;
t412 = t302 * t376 * t487 * t511;
t411 = (t356 * t553 - t515 / 0.2e1) * t375;
t334 = t393 * t411;
t410 = (t353 * t382 + t384 * t514) * t492;
t409 = (-t356 * t382 + t384 * t515) * t492;
t305 = ((t323 * t554 + t326 * t552) * t375 + (t411 + t410) * qJD(2)) * t393;
t306 = ((t323 * t552 + t326 * t553) * t375 + (t375 * t414 + t409) * qJD(2)) * t393;
t391 = 0.1e1 / t579;
t333 = t414 * t521;
t330 = qJD(1) * t334;
t329 = (t514 + t517) * t498 * t556;
t308 = (t409 + (-t384 * t461 + t385 * t462) * t375) * t498;
t307 = (t410 + (t384 * t462 + t385 * t461) * t375) * t498;
t304 = qJD(1) * t306;
t303 = qJD(1) * t305;
t299 = -t316 * t429 + t318 * t454 + qJD(2);
t294 = ((-t357 * t526 + t380 * t439 + t387 * t435) * t375 + 0.8e1 * t356 * t387 * t444 + ((t470 + (-t373 * t387 + t451) * t384) * t375 + (-0.4e1 * t326 * t497 - 0.2e1 * t356 * t513) * t548) * pkin(1)) * t393 * t454 - ((0.8e1 * t353 * t377 * t394 + 0.4e1 * pkin(2) * t375) * t511 * t518 + ((t468 * t569 + (t353 * t495 + t375 * t460) * t387) * t385 + (-0.4e1 * qJD(2) * t323 * t548 + t375 * t405) * t384) * t549) * t429 + (-t316 * t326 - t318 * t323) * t453 + (t310 * t338 * t346 + t323 * t337 * t347) * t316 * t545 * t570 + 0.2e1 * (t337 * t457 * t573 - t478 * t541) * pkin(3);
t1 = [MDP(6) * t513 - t387 * t384 * MDP(7) + (-t303 * t333 + t305 * t329) * MDP(11) + (t303 * t334 - t304 * t333 + t305 * t330 + t306 * t329) * MDP(12) + (-t294 * t333 + t299 * t305) * MDP(13) + (t304 * t334 + t306 * t330) * MDP(14) + (t294 * t334 + t299 * t306) * MDP(15) + ((-pkin(1) * t540 + t354 * t422) * MDP(24) + (-pkin(1) * t539 + t355 * t422) * MDP(25)) * t500 + ((-qJD(2) * t432 + t531 * t561) * MDP(19) + (t531 * t562 + (-t538 / 0.2e1 + 0.2e1 * t377 * t430) * t354) * MDP(20)) * qJD(1) * t391 - 0.2e1 * t572 * qJD(1) * qJD(2) + ((t355 * t523 / 0.2e1 - t406 * t300) * MDP(21) + (-t354 * t523 / 0.2e1 - t407 * t300) * MDP(22)) * t390 + ((-t330 * t497 + t304 * t385 + (t306 * t385 - t334 * t497) * qJD(1)) * MDP(17) + (t329 * t497 - t303 * t385 + (-t305 * t385 - t333 * t497) * qJD(1)) * MDP(18)) * pkin(2); -t329 * t307 * MDP(11) + (-t307 * t330 - t308 * t329) * MDP(12) + (-t299 * t307 + t302 * t303 + t329 * t507) * MDP(13) - t330 * t308 * MDP(14) + (-t299 * t308 + t302 * t304 + t330 * t507) * MDP(15) + (t294 * t302 + t299 * t507) * MDP(16) + (t353 * t412 + ((-t308 * t385 + t330 * t384) * qJD(1) + (t299 * t563 + t294 * t558 + (t302 * t563 + t327 * t299 / 0.2e1 - t463 * t353) * qJD(2)) * t521) * pkin(2)) * MDP(17) + (-t356 * t412 + ((t307 * t385 - t329 * t384) * qJD(1) + (t299 * t560 + t294 * t557 + (t299 * t559 + t302 * t560 + t356 * t463) * qJD(2)) * t521) * pkin(2)) * MDP(18) + (-t293 * t301 - (-t289 + t291) * t300) * MDP(23) + (((-t301 * t561 + t332 * t300 / 0.2e1) * t375 + t571 * t355) * MDP(21) + ((-t301 * t562 - t328 * t300 / 0.2e1) * t375 - t571 * t354) * MDP(22)) * t500 + (((pkin(1) * t537 / 0.2e1 - t354 * t455) * MDP(24) + (pkin(1) * t466 - t355 * t455) * MDP(25)) * t390 + ((-t332 * t531 / 0.4e1 + t432 / 0.2e1) * MDP(19) + (t328 * t531 / 0.4e1 + (t332 * t376 / 0.4e1 - t355 * t458) * t354) * MDP(20)) * t391 + t572) * qJD(1) ^ 2;];
tauc = t1;
