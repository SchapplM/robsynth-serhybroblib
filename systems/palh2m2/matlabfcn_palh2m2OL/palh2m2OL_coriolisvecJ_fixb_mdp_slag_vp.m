% Calculate minimal parameter regressor of Coriolis joint torque vector for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m2OL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'palh2m2OL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:07:04
% EndTime: 2020-06-30 18:07:37
% DurationCPUTime: 7.65s
% Computational Cost: add. (7861->396), mult. (20579->550), div. (0->0), fcn. (17486->10), ass. (0->178)
t443 = sin(qJ(3));
t447 = cos(qJ(3));
t448 = cos(qJ(2));
t519 = qJD(1) * t448;
t444 = sin(qJ(2));
t520 = qJD(1) * t444;
t412 = t443 * t520 - t447 * t519;
t414 = -t443 * t519 - t447 * t520;
t442 = sin(qJ(4));
t446 = cos(qJ(4));
t416 = t443 * t444 - t447 * t448;
t472 = t416 * qJD(3);
t394 = -t416 * qJD(2) - t472;
t456 = t394 * qJD(1);
t417 = t443 * t448 + t444 * t447;
t437 = qJD(2) + qJD(3);
t575 = t437 * qJD(1);
t457 = t417 * t575;
t517 = qJD(4) * t446;
t518 = qJD(4) * t442;
t337 = -t412 * t517 + t414 * t518 - t442 * t457 + t446 * t456;
t441 = sin(qJ(5));
t481 = t416 * t442 - t417 * t446;
t452 = t481 * t575;
t560 = t442 * t412 + t414 * t446;
t474 = t560 * qJD(4);
t451 = t474 + t452;
t483 = -t412 * t446 + t442 * t414;
t550 = cos(qJ(5));
t503 = qJD(5) * t550;
t516 = qJD(5) * t441;
t312 = t550 * t337 + t441 * t451 + t483 * t503 + t516 * t560;
t445 = cos(qJ(6));
t311 = t445 * t312;
t436 = qJD(4) + t437;
t428 = qJD(5) + t436;
t440 = sin(qJ(6));
t572 = t441 * t483 - t550 * t560;
t341 = t428 * t445 + t440 * t572;
t305 = -t341 * qJD(6) + t311;
t515 = qJD(6) * t440;
t501 = -t312 * t440 + t428 * t515;
t514 = qJD(6) * t445;
t306 = t514 * t572 - t501;
t313 = qJD(5) * t572 + t441 * t337 - t550 * t451;
t343 = -t440 * t428 + t445 * t572;
t351 = t441 * t560 + t483 * t550;
t589 = qJD(6) + t351;
t500 = t313 * t445 + t515 * t589;
t545 = t313 * t440;
t578 = t341 * t589;
t583 = t589 * t445;
t591 = -t313 * MDP(28) - t351 ^ 2 * MDP(26) + (-t351 * t428 + t312) * MDP(27) + (t351 * t440 * t589 + t500) * MDP(35) + (-t351 * MDP(25) + MDP(26) * t572 + t428 * MDP(28) - t341 * MDP(35) + MDP(36) * t589) * t572 + (t343 * t572 - t583 * t589 + t545) * MDP(34) + ((-t305 + t578) * t445 + t343 * t515 + (t343 * t351 + t306) * t440) * MDP(33);
t547 = pkin(4) * qJD(2);
t506 = t447 * t547;
t418 = pkin(2) * t437 + t506;
t507 = t443 * t547;
t489 = t442 * t507;
t398 = t446 * t418 - t489;
t393 = pkin(5) * t436 + t398;
t399 = t418 * t442 + t446 * t507;
t505 = t550 * t399;
t359 = t441 * t393 + t505;
t488 = qJD(3) * t506;
t526 = (qJD(3) + qJD(4)) * t489;
t368 = t446 * (qJD(4) * t418 + t488) - t526;
t529 = t443 * t446;
t480 = -t442 * t447 - t529;
t464 = (qJD(3) * t480 - t443 * t517) * pkin(4);
t369 = qJD(2) * t464 - t418 * t518;
t497 = -t441 * t368 + t550 * t369;
t323 = -qJD(5) * t359 + t497;
t432 = -pkin(4) * t448 - pkin(1);
t420 = t432 * qJD(1);
t397 = pkin(2) * t412 + t420;
t356 = -pkin(5) * t483 + t397;
t564 = t356 * t572;
t460 = t323 - t564;
t465 = -t368 * t550 - t441 * t369 - t393 * t503 + t399 * t516;
t588 = -t356 * t351 + t465;
t532 = t441 * t399;
t358 = t550 * t393 - t532;
t355 = pkin(3) * t428 + t358;
t326 = -pkin(3) * t351 + t356;
t317 = -t326 * t440 + t359 * t445;
t543 = t317 * t572;
t587 = -t323 * t440 - t355 * t583 - t543;
t546 = t305 * t440;
t586 = (-t343 * t583 - t546) * MDP(32) + t591;
t523 = pkin(3) * t572;
t316 = -t326 * t445 - t359 * t440;
t582 = t316 * t572 + t323 * t445;
t548 = pkin(5) * t560;
t536 = t483 * t397;
t576 = -t368 - t536;
t455 = t397 * t560 + t369;
t570 = t483 * MDP(18) * t560 + (-t436 * t483 + t337) * MDP(20) + (-t483 ^ 2 + t560 ^ 2) * MDP(19) + (-t436 * t560 + t451) * MDP(21);
t508 = qJD(1) * qJD(2);
t566 = -0.2e1 * t508;
t565 = MDP(5) * (t444 ^ 2 - t448 ^ 2);
t562 = t448 * MDP(4);
t559 = MDP(10) * t448;
t549 = pkin(4) * t444;
t541 = t343 * t445;
t391 = t446 * t416 + t417 * t442;
t354 = -t441 * t391 - t481 * t550;
t538 = t354 * t445;
t537 = t355 * t440;
t535 = t412 * t420;
t534 = t414 * t420;
t531 = t441 * t442;
t530 = t442 * t443;
t408 = t480 * t547;
t479 = t446 * t447 - t530;
t409 = t479 * t547;
t430 = pkin(2) * t446 + pkin(5);
t504 = t550 * t442;
t528 = -t408 * t550 + t441 * t409 - t430 * t516 + (-t442 * t503 + (-t441 * t446 - t504) * qJD(4)) * pkin(2);
t527 = t441 * t408 + t409 * t550 - t430 * t503 - (-t442 * t516 + (t446 * t550 - t531) * qJD(4)) * pkin(2);
t450 = qJD(1) ^ 2;
t524 = pkin(1) * t450;
t512 = t414 * pkin(2);
t435 = t444 * t547;
t502 = t444 * t508;
t493 = t589 * t355;
t473 = t417 * qJD(3);
t395 = t417 * qJD(2) + t473;
t381 = t395 * pkin(2) + t435;
t360 = -t512 - t548;
t328 = t360 + t523;
t431 = pkin(4) * t447 + pkin(2);
t406 = -pkin(4) * t530 + t431 * t446 + pkin(5);
t410 = pkin(4) * t529 + t431 * t442;
t373 = t441 * t406 + t410 * t550;
t434 = pkin(4) * t520;
t492 = qJD(6) * t373 - t328 - t434;
t407 = pkin(2) * t504 + t441 * t430;
t491 = qJD(6) * t407 - t328;
t490 = qJD(3) * (-qJD(2) - t437);
t487 = pkin(1) * t566;
t486 = t589 * t503;
t485 = (-qJD(3) + t437) * t547;
t484 = t394 * t442 + t395 * t446;
t401 = pkin(2) * t416 + t432;
t477 = -t391 * t550 + t441 * t481;
t476 = t406 * t550 - t441 * t410;
t339 = -qJD(4) * t391 + t394 * t446 - t395 * t442;
t463 = -t481 * qJD(4) + t484;
t314 = t477 * qJD(5) + t339 * t550 - t441 * t463;
t475 = -t314 * t445 + t354 * t515;
t365 = pkin(5) * t391 + t401;
t330 = -pkin(3) * t477 + t365;
t470 = qJD(6) * t330 * t589 + t355 * t314 + t323 * t354;
t387 = t431 * t517 + (qJD(3) * t479 - t443 * t518) * pkin(4);
t388 = -t431 * t518 + t464;
t332 = t476 * qJD(5) + t387 * t550 + t441 * t388;
t469 = t313 * t373 - t332 * t589 - t493;
t468 = t407 * t313 + t527 * t589 - t493;
t427 = pkin(4) * t502;
t319 = -pkin(5) * t474 + t427 + (pkin(2) * t417 - pkin(5) * t481) * t575;
t304 = t313 * pkin(3) + t319;
t315 = t354 * qJD(5) + t441 * t339 + t463 * t550;
t329 = pkin(5) * t463 + t381;
t459 = -(t315 * pkin(3) + t329) * t589 + t330 * t313 + (-qJD(6) * t359 - t304) * t477 + t355 * qJD(6) * t354;
t453 = -t414 * t412 * MDP(11) + (-t541 * t589 - t546) * MDP(32) + (t412 * t437 + t456) * MDP(13) + (-t414 * t437 - t457) * MDP(14) + (-t412 ^ 2 + t414 ^ 2) * MDP(12) + t570 + t591;
t449 = qJD(2) ^ 2;
t429 = pkin(5) * t550 + pkin(3);
t403 = -pkin(2) * t531 + t430 * t550 + pkin(3);
t400 = t434 - t512;
t371 = pkin(3) + t476;
t370 = pkin(2) * t457 + t427;
t364 = t398 * t550 - t532;
t363 = -t441 * t398 - t505;
t357 = t434 + t360;
t333 = -t373 * qJD(5) - t441 * t387 + t388 * t550;
t331 = -t548 + t523;
t324 = t326 * t515;
t1 = [(t313 * t365 + t315 * t356 - t319 * t477 - t329 * t351) * MDP(30) + (t312 * t477 - t313 * t354 + t314 * t351 - t315 * t572) * MDP(26) + (MDP(27) * t314 - MDP(28) * t315) * t428 + (-t449 * MDP(7) + MDP(9) * t487) * t444 + (MDP(20) * t339 - MDP(21) * t463) * t436 + (-t313 * t477 - t315 * t589) * MDP(36) + (t354 * t545 - t306 * t477 + t315 * t341 - (t314 * t440 + t354 * t514) * t589) * MDP(35) + (t305 * t477 - t313 * t538 - t315 * t343 - t475 * t589) * MDP(34) + (-t414 * t435 + t420 * t394 + (-t432 * t472 + (-t416 * t432 + t417 * t549) * qJD(2)) * qJD(1)) * MDP(17) + (-t381 * t483 + t370 * t391 + t397 * t484 + (-t397 * t481 - t401 * t560) * qJD(4) - t401 * t452) * MDP(23) + 0.2e1 * t502 * t562 + (t317 * t315 + (-(-qJD(6) * t326 - t465) * t477 + t470) * t445 - t459 * t440) * MDP(38) + (-t316 * t315 + t324 * t477 + (t465 * t477 + t470) * t440 + t459 * t445) * MDP(37) + (t337 * t401 + t339 * t397 - t370 * t481 - t381 * t560) * MDP(24) + (-t337 * t481 - t339 * t560) * MDP(18) + (MDP(13) * t394 - MDP(14) * t395) * t437 + (-t414 * t394 + t417 * t456) * MDP(11) + (t312 * t365 + t314 * t356 + t319 * t354 + t329 * t572) * MDP(31) + (t312 * t354 + t314 * t572) * MDP(25) + (-t394 * t412 + t414 * t395 + (t416 ^ 2 - t417 ^ 2) * t575) * MDP(12) + (-t337 * t391 + t339 * t483 + t560 * t484 + (-t452 - 0.2e1 * t474) * t481) * MDP(19) + t449 * t448 * MDP(6) + t487 * t559 + t565 * t566 + ((-t341 * t445 - t343 * t440) * t314 + (-t546 - t306 * t445 + (t341 * t440 - t541) * qJD(6)) * t354) * MDP(33) + (t412 * t435 + t420 * t395 + (t432 * t473 + (t416 * t549 + t417 * t432) * qJD(2)) * qJD(1)) * MDP(16) + (t305 * t538 - t343 * t475) * MDP(32); t450 * t565 + (t306 * t371 + t333 * t341 + t440 * t469 - t492 * t583 + t582) * MDP(37) + (t388 * t436 + t400 * t483 + t455) * MDP(23) + (t305 * t371 - t543 + t333 * t343 + (t492 * t589 - t323) * t440 + t469 * t445) * MDP(38) + (-t332 * t428 - t357 * t572 + t588) * MDP(31) + t453 + (-t387 * t436 + t400 * t560 + t576) * MDP(24) + (t333 * t428 + t351 * t357 + t460) * MDP(30) + t524 * t559 + (t534 + (-t412 * t520 + t443 * t490) * pkin(4)) * MDP(16) + (t535 + (t414 * t520 + t447 * t490) * pkin(4)) * MDP(17) + (MDP(9) * t524 - t450 * t562) * t444; (t447 * t485 + t535) * MDP(17) + (t403 * t306 + t341 * t528 + t440 * t468 - t491 * t583 + t582) * MDP(37) + (t351 * t360 + t428 * t528 + t460) * MDP(30) + (-t360 * t572 + t428 * t527 + t588) * MDP(31) + (t403 * t305 - t543 + t528 * t343 + (t491 * t589 - t323) * t440 + t468 * t445) * MDP(38) + (-t408 * t436 + (-t414 * t483 - t436 * t518) * pkin(2) + t455) * MDP(23) + (-t560 * t512 - t536 + t409 * t436 + (-t488 + (-pkin(2) * t436 - t418) * qJD(4)) * t446 + t526) * MDP(24) + t453 + (t443 * t485 + t534) * MDP(16); (t399 * t436 + t455) * MDP(23) + (t398 * t436 + t576) * MDP(24) + (-t351 * t548 - t564 - t363 * t428 + (-t505 + (-pkin(5) * t428 - t393) * t441) * qJD(5) + t497) * MDP(30) + (t364 * t428 + (-t428 * t503 + t560 * t572) * pkin(5) + t588) * MDP(31) + (t429 * t306 - (-t331 * t445 - t364 * t440) * t589 - t363 * t341 - t589 * t537 + (-t440 * t486 + (-qJD(5) * t341 - t514 * t589 + t545) * t441) * pkin(5) + t582) * MDP(37) + (t429 * t305 + (-t331 * t440 + t364 * t445) * t589 - t363 * t343 + (-t445 * t486 + (-qJD(5) * t343 + t500) * t441) * pkin(5) + t587) * MDP(38) + t570 + t586; (t359 * t428 + t460) * MDP(30) + (t358 * t428 + t588) * MDP(31) + (pkin(3) * t306 + t359 * t341 + t582 + (t358 * t440 + t445 * t523 - t537) * t589) * MDP(37) + (pkin(3) * t305 + (t358 * t445 - t440 * t523) * t589 + t359 * t343 + t587) * MDP(38) + t586; t343 * t341 * MDP(32) + (-t341 ^ 2 + t343 ^ 2) * MDP(33) + (t311 + t578) * MDP(34) + (t343 * t589 + t501) * MDP(35) - t313 * MDP(36) + (-t304 * t445 + t317 * t589 - t343 * t355 + t440 * t465 + t324) * MDP(37) + (t304 * t440 + t316 * t589 + t341 * t355 + t445 * t465) * MDP(38) + ((-MDP(34) * t572 + MDP(38) * t359) * t440 + (-MDP(34) * t428 - MDP(35) * t572 - MDP(37) * t359 + MDP(38) * t326) * t445) * qJD(6);];
tauc = t1;
