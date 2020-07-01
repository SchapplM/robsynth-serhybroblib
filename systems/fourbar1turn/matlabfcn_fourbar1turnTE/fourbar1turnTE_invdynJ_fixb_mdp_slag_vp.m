% Calculate vector of inverse dynamics joint torques for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnTE_convert_par2_MPV_fixb.m
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnTE_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnTE_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:23:08
% EndTime: 2020-06-27 16:23:24
% DurationCPUTime: 9.28s
% Computational Cost: add. (71829->347), mult. (101397->757), div. (3060->15), fcn. (27231->6), ass. (0->280)
t444 = 0.1e1 / pkin(4);
t448 = pkin(2) ^ 2;
t449 = pkin(1) ^ 2;
t438 = cos(qJ(2));
t632 = pkin(2) * t438;
t584 = -0.2e1 * pkin(1) * t632 + t449;
t430 = t448 + t584;
t673 = pkin(4) ^ 2;
t582 = pkin(3) ^ 2 - t673;
t426 = t430 - t582;
t436 = sin(qJ(2));
t655 = -pkin(3) - pkin(4);
t423 = (pkin(2) - t655) * (pkin(2) + t655) + t584;
t654 = pkin(4) - pkin(3);
t424 = (pkin(2) - t654) * (pkin(2) + t654) + t584;
t514 = (-t423 - t424) * pkin(1) * pkin(2);
t413 = t436 * t514;
t412 = qJD(2) * t413;
t609 = t423 * t424;
t450 = sqrt(-t609);
t417 = 0.1e1 / t450;
t614 = t412 * t417;
t541 = t436 * t614;
t431 = pkin(1) - t632;
t564 = 0.2e1 * t431 * pkin(1);
t595 = t438 * t450;
t370 = (-t541 + (-t595 + (t426 + t564) * t436) * qJD(2)) * pkin(2);
t600 = t436 * t450;
t407 = -pkin(2) * t600 + t426 * t431;
t428 = 0.1e1 / t430 ^ 2;
t574 = qJD(2) * t436;
t549 = pkin(2) * t574;
t525 = pkin(1) * t549;
t501 = t428 * t525;
t427 = 0.1e1 / t430;
t648 = -t427 / 0.2e1;
t466 = t370 * t648 + t407 * t501;
t361 = t466 * t444;
t434 = t436 ^ 2;
t604 = t434 * t448;
t553 = pkin(1) * t604;
t518 = qJD(2) * t553;
t538 = t450 * t574;
t573 = qJD(2) * t438;
t548 = pkin(2) * t573;
t587 = pkin(2) * t538 + t426 * t548;
t612 = t417 * t431;
t371 = t412 * t612 + 0.2e1 * t518 + t587;
t633 = pkin(2) * t436;
t421 = t426 * t633;
t408 = t431 * t450 + t421;
t647 = t427 / 0.2e1;
t465 = t371 * t647 - t408 * t501;
t363 = t465 * t444;
t401 = 0.1e1 / t407;
t402 = 0.1e1 / t407 ^ 2;
t618 = t402 * t408;
t481 = t361 * t618 + t363 * t401;
t425 = t430 + t582;
t637 = pkin(1) * t436;
t422 = t425 * t637;
t432 = pkin(1) * t438 - pkin(2);
t409 = -t432 * t450 + t422;
t598 = t438 * t409;
t406 = -pkin(1) * t600 - t425 * t432;
t602 = t436 * t406;
t676 = t598 + t602;
t657 = -0.2e1 * t432;
t565 = pkin(2) * t657;
t529 = t425 + t565;
t369 = (-t541 + (t436 * t529 - t595) * qJD(2)) * pkin(1);
t447 = 0.1e1 / pkin(3);
t488 = -0.2e1 * t501;
t362 = (t369 * t427 + t406 * t488) * t447;
t399 = 0.1e1 / t406 ^ 2;
t620 = t399 * t409;
t603 = t434 * t449;
t539 = qJD(2) * t603;
t517 = pkin(2) * t539;
t586 = (t425 * t573 + t538) * pkin(1);
t611 = t417 * t432;
t372 = -t412 * t611 + 0.2e1 * t517 + t586;
t364 = (t372 * t427 + t409 * t488) * t447;
t398 = 0.1e1 / t406;
t625 = t364 * t398;
t675 = -t362 * t620 + t625;
t613 = t417 * t413;
t373 = t422 + (-t595 + (t565 - t613) * t436) * pkin(1);
t531 = -t409 / 0.2e1 + t373 / 0.2e1;
t377 = -t413 * t611 + 0.2e1 * pkin(2) * t603 + (t425 * t438 + t600) * pkin(1);
t649 = -t406 / 0.2e1;
t532 = t649 - t377 / 0.2e1;
t634 = pkin(2) * t428;
t563 = pkin(1) * t634;
t667 = (t406 * t434 + t436 * t598) * t563;
t674 = (t436 * t531 - t438 * t532) * t427 - t667;
t672 = -0.4e1 * t407;
t561 = 0.2e1 * t428;
t671 = 0.6e1 * t436;
t670 = 2 * qJDD(2);
t636 = pkin(1) * t447;
t403 = t401 * t402;
t404 = t408 ^ 2;
t617 = t403 * t404;
t355 = -t370 * t617 + t371 * t618;
t385 = t402 * t404 + 0.1e1;
t382 = 0.1e1 / t385 ^ 2;
t669 = t355 * t382;
t381 = 0.1e1 / t385;
t630 = pkin(4) * t430;
t555 = t381 * t630;
t520 = t402 * t555;
t668 = t444 * t408 * t520;
t472 = (t598 / 0.2e1 + t602 / 0.2e1) * t427;
t572 = qJD(1) * qJD(2);
t530 = t436 * t572;
t512 = pkin(2) * t530;
t571 = qJDD(1) * t438;
t666 = pkin(2) * t571 - t512;
t441 = qJD(2) ^ 2;
t596 = t438 * t441;
t485 = qJDD(2) * t436 + t596;
t607 = t427 * t447;
t535 = -t607 / 0.2e1;
t375 = t676 * qJD(1) * t535;
t665 = MDP(11) * t375;
t599 = t438 * t406;
t601 = t436 * t409;
t471 = (t601 / 0.2e1 - t599 / 0.2e1) * t427;
t380 = t447 * t471;
t376 = qJD(1) * t380;
t664 = MDP(14) * t376;
t374 = t421 + (-t595 + (t564 - t613) * t436) * pkin(2);
t378 = t413 * t612 + 0.2e1 * t553 + (t426 * t438 + t600) * pkin(2);
t658 = 0.2e1 * t430;
t509 = 0.2e1 * t382 * t658;
t615 = t408 * t430;
t524 = 0.4e1 * t403 * t615;
t562 = pkin(1) * t633;
t540 = -0.4e1 * t562;
t559 = t402 * t658;
t523 = -0.4e1 * t448 * t603;
t410 = (t438 * t514 + t523) * qJD(2);
t610 = 0.4e1 * t417 / t609;
t511 = t412 * t413 * t610;
t489 = -t511 / 0.4e1;
t645 = -t436 / 0.2e1;
t462 = (t436 * t489 + 0.2e1 * (t410 * t645 - t412 * t438) * t417) * t427;
t429 = t427 * t428;
t605 = t429 * t448;
t490 = t539 * t605;
t554 = t401 * t630;
t495 = t381 * t444 * t554;
t558 = t448 * t637;
t504 = 0.6e1 * t438 * t558;
t543 = t427 * t614;
t597 = t438 * t427;
t606 = t428 * t436;
t624 = t371 * t428;
t638 = pkin(1) * t428;
t661 = 0.4e1 * t408;
t590 = -0.2e1 * ((0.4e1 * t518 + t587) * t648 + t490 * t672 + (-t462 / 0.2e1 + (t370 * t606 + (-t431 * t597 + (t374 * t436 + t407 * t438) * t428) * qJD(2)) * pkin(1)) * pkin(2)) * t668 - 0.2e1 * ((t431 * t511 / 0.4e1 + t410 * t612 + qJD(2) * t504) * t647 + t490 * t661 + ((t543 / 0.2e1 - pkin(1) * t624) * t436 + ((t595 + (-t426 + t613) * t436) * t647 + (-t378 * t436 - t408 * t438) * t638) * qJD(2)) * pkin(2)) * t495;
t659 = -0.2e1 * t430;
t335 = (t481 * (-t374 * t617 + t378 * t618) * t509 + ((t374 * t559 + t401 * t540) * t363 + (t374 * t524 + (t378 * t659 + t408 * t540) * t402) * t361) * t381) * pkin(4) + t590;
t557 = pkin(2) * t606;
t528 = pkin(1) * t557;
t365 = (t374 * t648 + t407 * t528) * t444;
t367 = (t378 * t647 - t408 * t528) * t444;
t479 = t365 * t618 + t367 * t401;
t507 = -0.4e1 * t525;
t337 = (t479 * t355 * t509 + ((t370 * t559 + t401 * t507) * t367 + (t370 * t524 + (t371 * t659 + t408 * t507) * t402) * t365) * t381) * pkin(4) + t590;
t506 = 0.2e1 * t555;
t346 = t481 * t506;
t347 = t479 * t506;
t663 = t427 * (t337 / 0.2e1 - t335 / 0.2e1) - (-qJD(2) * t347 + t346) * t528;
t660 = 0.2e1 * t409;
t653 = -t370 / 0.2e1;
t652 = t371 / 0.2e1;
t651 = t374 / 0.2e1;
t650 = -t378 / 0.2e1;
t646 = t428 / 0.4e1;
t644 = -t438 / 0.2e1;
t437 = sin(qJ(1));
t643 = g(1) * t437;
t439 = cos(qJ(1));
t642 = g(1) * t439;
t641 = g(2) * t437;
t640 = g(2) * t439;
t639 = pkin(1) * t427;
t635 = pkin(2) * t427;
t631 = pkin(3) * t430;
t629 = pkin(2) * qJD(1);
t484 = 0.8e1 * t490;
t566 = -0.2e1 * t634;
t622 = t377 * t436;
t628 = ((t449 * t548 * t671 - t410 * t611 + t432 * t489) * t427 + t409 * t484 + ((t372 * t566 + t543) * t436 + ((t595 + (-t425 + t613) * t436) * t427 + (-t598 - t622) * pkin(2) * t561) * qJD(2)) * pkin(1)) * t447 * t398;
t510 = -0.2e1 * t528;
t366 = (t373 * t427 + t406 * t510) * t447;
t368 = (t377 * t427 + t409 * t510) * t447;
t405 = t409 ^ 2;
t386 = t399 * t405 + 0.1e1;
t383 = 0.1e1 / t386;
t556 = t383 * t631;
t521 = t399 * t556;
t496 = t409 * t521;
t522 = t398 * t556;
t348 = -t366 * t496 + t368 * t522 + 0.1e1;
t627 = t348 * t447;
t623 = t373 * t436;
t384 = 0.1e1 / t386 ^ 2;
t621 = t384 * t430;
t400 = t398 * t399;
t619 = t400 * t405;
t616 = t408 * t428;
t608 = t427 * t432;
t594 = t439 * t447;
t593 = t441 * t449;
t592 = t450 * t441;
t342 = (((0.4e1 * t517 + t586) * t427 + t406 * t484) * t447 + (t462 + (-0.2e1 * t369 * t606 + (t597 * t657 + (-t599 - t623) * t561) * qJD(2)) * pkin(2)) * t636) * t496;
t513 = 0.2e1 * t562;
t560 = t400 * t660;
t336 = -t342 + (-0.2e1 * t675 * (-t373 * t619 + t377 * t620) * t621 + (t675 * t513 + (-t364 * t373 * t399 + t628 + (t373 * t560 - t377 * t399) * t362) * t430) * t383) * pkin(3);
t478 = t366 * t620 - t368 * t398;
t356 = -t369 * t619 + t372 * t620;
t546 = t356 * t621;
t338 = -t342 + (0.2e1 * t478 * t546 + (-t478 * qJD(2) * t513 + (-t368 * t369 * t399 + t628 + (t369 * t560 - t372 * t399) * t366) * t430) * t383) * pkin(3);
t591 = -t336 + t338;
t534 = t607 / 0.2e1;
t589 = (t534 * t599 + t535 * t601) * t437;
t588 = t676 * t439 * t534;
t583 = -t438 ^ 2 + t434;
t445 = 0.1e1 / t673;
t581 = MDP(19) * t445;
t580 = MDP(20) * t445;
t579 = MDP(21) * t444;
t578 = MDP(22) * t444;
t577 = MDP(24) * t444;
t576 = MDP(25) * t444;
t575 = qJD(2) * t428;
t570 = qJDD(2) * t425;
t569 = qJDD(2) * t426;
t568 = qJDD(2) * t427;
t552 = qJD(1) * t639;
t551 = t436 * t629;
t550 = t438 * t629;
t542 = t441 * t604;
t537 = t408 * t647;
t442 = qJD(1) ^ 2;
t536 = t442 * t647;
t527 = t429 * t562;
t526 = -qJDD(1) * t407 / 0.2e1;
t516 = pkin(2) * t535;
t515 = pkin(2) * t534;
t500 = t429 * t525;
t499 = t641 + t642;
t497 = t442 * t449 * t557;
t493 = pkin(2) * t568 * t627;
t492 = qJD(2) * t516;
t491 = qJD(2) * t515;
t487 = t407 * t428 - t427 * t431;
t483 = t642 / 0.2e1 + t641 / 0.2e1;
t475 = 0.2e1 * qJD(2) * t614 + qJDD(2) * t450;
t474 = t428 * t441 * t558 * t627;
t473 = t412 ^ 2 * t610 / 0.4e1 + t417 * (t441 * t523 + t485 * t514);
t468 = (-t409 * t434 + t436 * t599) * t563;
t467 = -t426 * t441 + t475;
t464 = -t473 + t592;
t463 = t449 * t512 * t561 + (-qJDD(1) * pkin(1) - t643 / 0.2e1 + t640 / 0.2e1) * t427;
t461 = t447 * t674;
t460 = t447 * (t468 + (-t436 * t532 - t438 * t531) * t427);
t459 = (t369 * t645 + t372 * t644) * t427 + (t471 + t667) * qJD(2);
t458 = (t436 * t372 / 0.2e1 + t369 * t644) * t427 + (t468 + t472) * qJD(2);
t354 = qJD(1) * t460;
t353 = qJD(1) * t461;
t350 = (qJD(1) * t458 + qJDD(1) * t471) * t447;
t349 = (qJD(1) * t459 - qJDD(1) * t472) * t447;
t345 = -t362 * t496 + t364 * t522 + qJD(2);
t340 = qJDD(2) + ((-t473 * t608 + (0.8e1 * t409 * t429 * t542 + (t434 * t670 + t596 * t671) * t635) * t449) * t447 + (((t570 + t592) * t427 + t409 * t441 * t566) * t438 + ((-t425 * t441 + t475) * t427 + (-0.4e1 * t372 * qJD(2) - 0.2e1 * t409 * qJDD(2)) * t634) * t436) * t636) * t522 - ((0.8e1 * t406 * t605 + 0.4e1 * t635) * t447 * t434 * t593 + ((-t475 * t427 + (t406 * t566 + t427 * t529) * t441) * t438 + ((t464 + t570) * t427 + (-0.4e1 * t369 * t575 + (-t406 * t428 - t608) * t670) * pkin(2)) * t436) * t636) * t496 + (-t362 * t372 - t364 * t369) * t521 + (t356 * t384 * t399 + t369 * t383 * t400) * t362 * t631 * t660 + (0.2e1 * t383 * t525 * t675 - 0.2e1 * t546 * t625) * pkin(3);
t339 = -0.2e1 * ((t473 * t431 + t441 * t504) * t647 + (t429 * t593 * t661 + pkin(1) * t568) * t604 + (((t569 + t592) * t438 + t467 * t436) * t647 + (-0.2e1 * t371 * t574 - t408 * t485) * t638) * pkin(2)) * t495 - 0.2e1 * ((t429 * t449 * t672 - 0.2e1 * t639) * t542 + ((pkin(1) * t441 * t487 + t467 * t647) * t438 + ((t464 + t569) * t648 + (qJDD(2) * t487 + 0.2e1 * t370 * t575) * pkin(1)) * t436) * pkin(2)) * t668 - 0.2e1 * t361 * t371 * t520 + 0.2e1 * (t370 * t520 + 0.2e1 * t554 * t669) * t363 + 0.2e1 * (0.2e1 * (t370 * t381 * t403 + t402 * t669) * t361 * t615 - 0.2e1 * t481 * t381 * t525) * pkin(4);
t1 = [qJDD(1) * MDP(1) + t499 * MDP(3) + (qJDD(1) * t434 + 0.2e1 * t438 * t530) * MDP(4) + 0.2e1 * (t436 * t571 - t572 * t583) * MDP(5) + t485 * MDP(6) + (qJDD(2) * t438 - t436 * t441) * MDP(7) + (-g(1) * t589 + t350 * t632 - t376 * t549) * MDP(17) + (-g(2) * t588 - t349 * t632 + t375 * t549) * MDP(18) + (t404 * qJDD(1) * t646 + (-t404 * t500 + t616 * t652) * qJD(1)) * t581 + (t526 * t616 + (t616 * t653 + (-t624 / 0.2e1 + 0.2e1 * t408 * t500) * t407) * qJD(1)) * t580 + (t339 * t537 - t346 * t465) * t579 + (t339 * t407 * t648 - t346 * t466) * t578 + (-t370 * t552 + t407 * t463) * t577 + (-t371 * t552 + t408 * t463) * t576 + (t349 * MDP(12) + t350 * MDP(14) + t340 * MDP(15) + (-t640 + t666) * MDP(17)) * t380 + (-t436 * MDP(10) + t438 * MDP(9) + MDP(2)) * (-t640 + t643) + ((MDP(12) * t375 + MDP(15) * t345 + MDP(17) * t550 + t664) * t458 + (MDP(12) * t376 + MDP(13) * t345 - MDP(18) * t550 + t665) * t459 - (t349 * MDP(11) + t350 * MDP(12) + t340 * MDP(13) + (-t643 - t666) * MDP(18)) * t472) * t447; t436 * qJDD(1) * MDP(6) + MDP(7) * t571 + (qJDD(2) * MDP(8)) + (-g(3) * t438 + t436 * t499) * MDP(9) + (g(3) * t436 + t438 * t499) * MDP(10) + t353 * t665 + (t353 * t376 - t354 * t375) * MDP(12) + (t345 * t353 + t348 * t349 + t375 * t591) * MDP(13) - t354 * t664 + (-t345 * t354 + t348 * t350 + t376 * t591) * MDP(15) + (t340 * t348 + t345 * t591) * MDP(16) + (t369 * t348 * t492 + t493 * t649 + t376 * t551 - t354 * t550 - g(1) * (((t373 * t644 + t622 / 0.2e1) * t427 + t468) * t594 + t588) - t460 * t641 + g(3) * t461 + (t369 * t516 + t373 * t491) * t345 + (t336 * t491 + t338 * t492 + t340 * t516 + t474) * t406) * MDP(17) + (t372 * t348 * t491 - t375 * t551 - t353 * t550 - g(1) * t674 * t594 - g(2) * (((t438 * t377 / 0.2e1 + t623 / 0.2e1) * t427 - t667) * t447 * t437 + t589) - g(3) * t460 + (t372 * t515 + t377 * t492) * t345 + (t340 * t515 - t474 + t493 / 0.2e1 + t338 * t491 + t336 * t492) * t409) * MDP(18) + (-qJDD(1) * t347 * t537 + ((-t346 * t650 - t347 * t652) * t427 + t663 * t408) * qJD(1)) * t579 + (-t427 * t347 * t526 + ((-t346 * t651 - t347 * t653) * t427 - t663 * t407) * qJD(1)) * t578 + (-t339 * t347 - (-t335 + t337) * t346) * MDP(23) + (-t407 * t497 + (g(3) * t650 + t374 * t483) * t427 + (t374 * t536 + (g(3) * t408 - t407 * t499) * t557) * pkin(1)) * t577 + (-t408 * t497 + (g(3) * t651 + t378 * t483) * t427 + (t378 * t536 + (-g(3) * t407 - t408 * t499) * t557) * pkin(1)) * t576 + (-t436 * t438 * MDP(4) + t583 * MDP(5) + (-t378 * t616 / 0.4e1 + t404 * t527 / 0.2e1) * t581 + (t374 * t616 / 0.4e1 + (t378 * t646 - t408 * t527) * t407) * t580) * t442;];
tau = t1;
