% Calculate vector of inverse dynamics joint torques for
% fourbar1turnDE2
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
%   see fourbar1turnDE2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnDE2_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:49:12
% EndTime: 2020-06-27 16:49:39
% DurationCPUTime: 10.86s
% Computational Cost: add. (106593->328), mult. (152129->732), div. (5460->21), fcn. (40979->13), ass. (0->277)
t478 = pkin(2) ^ 2;
t479 = pkin(1) ^ 2;
t467 = cos(qJ(2));
t671 = pkin(2) * t467;
t695 = -2 * pkin(1);
t618 = t671 * t695 + t479;
t459 = t478 + t618;
t456 = 0.1e1 / t459;
t457 = 0.1e1 / t459 ^ 2;
t458 = t456 * t457;
t717 = 0.8e1 * t458;
t716 = pkin(1) * pkin(2);
t470 = qJD(2) ^ 2;
t715 = t470 * t478;
t465 = sin(qJ(2));
t463 = t465 ^ 2;
t634 = t463 * t479;
t714 = 0.2e1 * t634 * pkin(2);
t707 = pkin(3) ^ 2;
t708 = pkin(4) ^ 2;
t616 = t707 - t708;
t454 = t459 + t616;
t461 = pkin(1) * t467 - pkin(2);
t599 = -0.2e1 * pkin(2) * t461;
t567 = t454 + t599;
t687 = -pkin(3) - pkin(4);
t452 = (pkin(2) - t687) * (pkin(2) + t687) + t618;
t686 = pkin(4) - pkin(3);
t453 = (pkin(2) - t686) * (pkin(2) + t686) + t618;
t552 = (-t452 - t453) * t716;
t442 = t465 * t552;
t441 = qJD(2) * t442;
t642 = t452 * t453;
t482 = sqrt(-t642);
t446 = 0.1e1 / t482;
t647 = t441 * t446;
t573 = t465 * t647;
t630 = t467 * t482;
t399 = (-t573 + (t567 * t465 - t630) * qJD(2)) * pkin(1);
t607 = qJD(2) * t465;
t570 = t482 * t607;
t606 = qJD(2) * t467;
t644 = t446 * t461;
t402 = -t441 * t644 + qJD(2) * t714 + (t454 * t606 + t570) * pkin(1);
t632 = t465 * t482;
t435 = -pkin(1) * t632 - t454 * t461;
t425 = t435 ^ 2;
t426 = 0.1e1 / t435;
t428 = t426 / t425;
t676 = pkin(1) * t465;
t451 = t454 * t676;
t438 = -t461 * t482 + t451;
t434 = t438 ^ 2;
t656 = t428 * t434;
t427 = 0.1e1 / t435 ^ 2;
t657 = t427 * t438;
t380 = -t399 * t656 + t402 * t657;
t423 = t427 * t434 + 0.1e1;
t421 = 0.1e1 / t423 ^ 2;
t658 = t421 * t459;
t713 = -0.2e1 * t380 * t658;
t712 = 0.2e1 * qJD(2);
t477 = 0.1e1 / t707;
t672 = pkin(2) * t465;
t597 = pkin(1) * t672;
t565 = t458 * t597;
t529 = 0.4e1 * t565;
t623 = t425 + t434;
t503 = t623 * t529;
t596 = 0.2e1 * t457;
t376 = ((t399 * t435 + t402 * t438) * t596 - qJD(2) * t503) * t477;
t711 = t376 / 0.2e1;
t710 = t676 * t715;
t476 = 0.1e1 / pkin(3);
t586 = pkin(2) * t607;
t563 = pkin(1) * t586;
t538 = t457 * t563;
t524 = -0.2e1 * t538;
t390 = (t399 * t456 + t435 * t524) * t476;
t392 = (t402 * t456 + t438 * t524) * t476;
t420 = 0.1e1 / t423;
t670 = pkin(3) * t459;
t591 = t420 * t670;
t559 = t427 * t591;
t534 = t438 * t559;
t560 = t426 * t591;
t371 = -t390 * t534 + t392 * t560 + qJD(2);
t646 = t446 * t442;
t403 = t451 + (-t630 + (t599 - t646) * t465) * pkin(1);
t592 = t457 * t672;
t566 = pkin(1) * t592;
t546 = -0.2e1 * t566;
t394 = (t403 * t456 + t435 * t546) * t476;
t405 = -t442 * t644 + t714 + (t454 * t467 + t632) * pkin(1);
t396 = (t405 * t456 + t438 * t546) * t476;
t374 = -t394 * t534 + t396 * t560 + 0.1e1;
t378 = ((t403 * t435 + t405 * t438) * t596 - t503) * t477;
t417 = t623 * t477 * t457;
t412 = t417 ^ (-0.1e1 / 0.2e1);
t413 = t412 / t417;
t709 = (t371 * t711 + (t374 * t711 - t371 * t378 / 0.2e1) * qJD(2)) * t413;
t665 = t392 * t426;
t513 = t390 * t657 - t665;
t473 = 0.1e1 / pkin(4);
t455 = t459 - t616;
t460 = pkin(1) - t671;
t436 = -pkin(2) * t632 + t455 * t460;
t518 = t436 * t538;
t598 = 0.2e1 * t460 * pkin(1);
t400 = (-t573 + (-t630 + (t455 + t598) * t465) * qJD(2)) * pkin(2);
t664 = t400 * t456;
t389 = (-t664 / 0.2e1 + t518) * t473;
t450 = t455 * t672;
t437 = t460 * t482 + t450;
t517 = t437 * t538;
t635 = t463 * t478;
t588 = pkin(1) * t635;
t645 = t446 * t460;
t401 = t588 * t712 + t441 * t645 + (t455 * t606 + t570) * pkin(2);
t663 = t401 * t456;
t391 = (t663 / 0.2e1 - t517) * t473;
t430 = 0.1e1 / t436;
t431 = 0.1e1 / t436 ^ 2;
t655 = t431 * t437;
t514 = t389 * t655 + t391 * t430;
t705 = -0.2e1 * t401;
t702 = 2 * qJDD(2);
t404 = t450 + (-t630 + (t598 - t646) * t465) * pkin(2);
t406 = t442 * t645 + 0.2e1 * t588 + (t455 * t467 + t632) * pkin(2);
t474 = 0.1e1 / t708;
t429 = t436 ^ 2;
t433 = t437 ^ 2;
t622 = t429 + t433;
t502 = t622 * t529;
t377 = ((t404 * t436 + t406 * t437) * t596 - t502) * t474;
t701 = t377 / 0.2e1;
t700 = pkin(1) * t476;
t432 = t430 / t429;
t654 = t432 * t433;
t379 = -t400 * t654 + t401 * t655;
t422 = t431 * t433 + 0.1e1;
t419 = 0.1e1 / t422 ^ 2;
t699 = t379 * t419;
t648 = t438 * t467;
t653 = t435 * t465;
t698 = (t648 + t653) * t456;
t631 = t467 * t470;
t519 = qJDD(2) * t465 + t631;
t418 = 0.1e1 / t422;
t669 = pkin(4) * t459;
t590 = t418 * t669;
t543 = 0.2e1 * t590;
t372 = t514 * t543;
t684 = -t456 / 0.2e1;
t393 = (t404 * t684 + t436 * t566) * t473;
t683 = t456 / 0.2e1;
t395 = (t406 * t683 - t437 * t566) * t473;
t512 = t393 * t655 + t395 * t430;
t373 = t512 * t543;
t550 = 0.2e1 * t597;
t544 = -0.4e1 * t563;
t572 = -0.4e1 * t597;
t547 = t437 * t572;
t650 = t437 * t459;
t562 = 0.4e1 * t432 * t650;
t690 = 0.2e1 * t459;
t594 = t431 * t690;
t691 = -0.2e1 * t459;
t627 = (0.2e1 * (t514 * (-t404 * t654 + t406 * t655) - t512 * t379) * t419 * t690 + ((t404 * t594 + t430 * t572) * t391 + (t404 * t562 + (t406 * t691 + t547) * t431) * t389 - (t400 * t594 + t430 * t544) * t395 - (t400 * t562 + (t401 * t691 + t437 * t544) * t431) * t393) * t418) * pkin(4);
t696 = t457 * t550 * (-qJD(2) * t373 + t372) + t627 * t456;
t416 = t622 * t474 * t457;
t414 = 0.1e1 / t416;
t410 = t416 ^ (-0.1e1 / 0.2e1);
t693 = -0.2e1 * t438;
t692 = 0.2e1 * t438;
t375 = ((t400 * t436 + t401 * t437) * t596 - qJD(2) * t502) * t474;
t685 = t375 / 0.2e1;
t468 = cos(qJ(1));
t682 = g(1) * t468;
t466 = sin(qJ(1));
t681 = g(2) * t466;
t680 = g(3) * t436;
t679 = g(3) * t437;
t678 = pkin(1) * t456;
t674 = pkin(2) * t456;
t673 = pkin(2) * t457;
t661 = t410 * t456;
t411 = t410 * t414;
t660 = t411 * t456;
t659 = t413 * t456;
t652 = t435 * t467;
t651 = t437 * t457;
t649 = t438 * t465;
t641 = t456 * t461;
t471 = qJD(1) ^ 2;
t639 = t456 * t471;
t638 = t456 * t476;
t633 = t465 * t467;
t629 = t470 * t479;
t628 = t482 * t470;
t511 = t394 * t657 - t396 * t426;
t626 = (0.2e1 * t513 * (-t403 * t656 + t405 * t657) * t658 + t511 * t713 + ((t511 * qJD(2) - t513) * t550 + ((t403 * t390 - t399 * t394) * t428 * t692 + (-t405 * t390 - t392 * t403 + t402 * t394 + t396 * t399) * t427) * t459) * t420) * pkin(3);
t624 = -t403 + t438;
t617 = -t467 ^ 2 + t463;
t615 = MDP(19) * t474;
t614 = MDP(20) * t474;
t613 = MDP(21) * t473;
t612 = MDP(22) * t473;
t611 = MDP(24) * t473;
t610 = MDP(25) * t473;
t501 = t412 * t698;
t387 = t476 * t501;
t385 = qJD(1) * t387;
t609 = qJD(1) * t476;
t608 = qJD(2) * t457;
t605 = qJD(1) * qJD(2);
t604 = qJDD(1) * t467;
t603 = qJDD(2) * t454;
t602 = qJDD(2) * t455;
t600 = -0.2e1 * t673;
t589 = t430 * t669;
t587 = qJD(1) * t678;
t504 = (t441 ^ 2 / t642 + t519 * t552 - 0.4e1 * t634 * t715) * t446;
t496 = -t504 + t628;
t508 = qJDD(2) * t482 + t647 * t712;
t497 = -t455 * t470 + t508;
t522 = t436 * t457 - t456 * t460;
t558 = t431 * t590;
t574 = t470 * t635;
t359 = -0.2e1 * t389 * t401 * t558 + 0.2e1 * (t400 * t558 + 0.2e1 * t589 * t699) * t391 + 0.2e1 * (-((t504 * t460 + 0.6e1 * t467 * t710) * t683 + (0.4e1 * t437 * t458 * t629 + qJDD(2) * t678) * t635 + (((t602 + t628) * t467 + t497 * t465) * t683 + (-t519 * t437 + t607 * t705) * t457 * pkin(1)) * pkin(2)) * t418 * t589 - ((-0.4e1 * t436 * t458 * t479 - 0.2e1 * t678) * t574 + ((t522 * t470 * pkin(1) + t497 * t683) * t467 + ((t496 + t602) * t684 + (t522 * qJDD(2) + 0.2e1 * t400 * t608) * pkin(1)) * t465) * pkin(2)) * t437 * t558) * t473 + 0.2e1 * (0.2e1 * (t432 * t418 * t400 + t431 * t699) * t389 * t650 - 0.2e1 * t514 * t418 * t563) * pkin(4);
t585 = t359 * t661;
t584 = t377 * t660;
t580 = t436 * t660;
t579 = t437 * t660;
t415 = 0.1e1 / t416 ^ 2;
t578 = t415 * t433 * t457;
t577 = t435 * t638;
t576 = t438 * t638;
t568 = qJDD(1) * t414 * t457;
t564 = -0.2e1 * t587;
t360 = qJDD(2) + ((-t504 * t641 + (t438 * t574 * t717 + (t463 * t702 + 0.6e1 * t465 * t631) * t674) * t479) * t476 + (((t603 + t628) * t456 + t438 * t470 * t600) * t467 + ((-t454 * t470 + t508) * t456 + (-0.4e1 * qJD(2) * t402 + qJDD(2) * t693) * t673) * t465) * t700) * t560 - ((t435 * t478 * t717 + 0.4e1 * t674) * t476 * t463 * t629 + ((-t508 * t456 + (t435 * t600 + t567 * t456) * t470) * t467 + ((t496 + t603) * t456 + (-0.4e1 * t399 * t608 + (-t435 * t457 - t641) * t702) * pkin(2)) * t465) * t700) * t534 + (-t390 * t402 - t392 * t399) * t559 + (t427 * t421 * t380 + t428 * t420 * t399) * t390 * t670 * t692 + (-0.2e1 * t420 * t513 * t563 + t665 * t713) * pkin(3);
t553 = qJDD(2) * t374 + t360;
t551 = t596 * t716;
t549 = t415 * t436 * t651;
t542 = qJDD(1) * t373 * t661;
t537 = t458 * t563;
t536 = t681 + t682;
t535 = g(1) * t466 - g(2) * t468;
t531 = t375 * t411 * t587;
t528 = t412 * t467 * t577;
t527 = t412 * t465 * t576;
t520 = -0.2e1 * t471 * t479 * t592;
t510 = t372 * t701 - t373 * t685;
t507 = 0.2e1 * t536;
t506 = (t649 - t652) * t456;
t505 = t412 * t457 * t476 * t710;
t500 = -t471 * pkin(1) / 0.2e1 - t682 / 0.2e1 - t681 / 0.2e1;
t499 = (t653 / 0.2e1 + t648 / 0.2e1) * t659;
t498 = (-t649 / 0.2e1 + t652 / 0.2e1) * t659;
t495 = (t435 * t463 + t438 * t633) * t551;
t494 = (t435 * t633 - t438 * t463) * t551;
t493 = 0.4e1 * qJD(1) * t457 * t479 * t586 + (qJDD(1) * t695 - t535) * t456;
t491 = t376 * t499 + ((-t399 * t465 - t402 * t467) * t456 + (t506 + t495) * qJD(2)) * t412;
t490 = t376 * t498 + ((-t399 * t467 + t402 * t465) * t456 + (t494 + t698) * qJD(2)) * t412;
t409 = qJ(2) + atan2(t576, t577);
t408 = cos(t409);
t407 = sin(t409);
t397 = qJD(1) * t528;
t388 = t527 - t528;
t386 = qJD(1) * t527 - t397;
t366 = (t378 * t498 + (t494 + (t624 * t467 + (t405 + t435) * t465) * t456) * t412) * t609;
t365 = -t397 + (t378 * t499 + (t495 + (-t405 * t467 + t624 * t465) * t456) * t412) * t609;
t364 = t490 * t476;
t363 = t491 * t476;
t362 = (t412 * qJDD(1) * t506 + t490 * qJD(1)) * t476;
t361 = (t491 * qJD(1) - qJDD(1) * t501) * t476;
t1 = [qJDD(1) * MDP(1) + t536 * MDP(3) + (qJDD(1) * t463 + 0.2e1 * t605 * t633) * MDP(4) + 0.2e1 * (t465 * t604 - t617 * t605) * MDP(5) + t519 * MDP(6) + (qJDD(2) * t467 - t465 * t470) * MDP(7) + (-t361 * t387 - t363 * t385) * MDP(11) + (t361 * t388 - t362 * t387 + t363 * t386 - t364 * t385) * MDP(12) + (-t360 * t387 + t363 * t371) * MDP(13) + (t362 * t388 + t364 * t386) * MDP(14) + (t360 * t388 + t364 * t371) * MDP(15) + (t433 * t568 + (-t375 * t578 + (0.2e1 * t401 * t651 - 0.4e1 * t433 * t537) * t414) * qJD(1)) * t615 + (-0.2e1 * t436 * t437 * t568 + (0.2e1 * t375 * t549 + (-0.2e1 * t400 * t651 + (0.8e1 * t437 * t537 + t457 * t705) * t436) * t414) * qJD(1)) * t614 + (t437 * t585 - (-t375 * t579 / 0.2e1 + (-0.2e1 * t517 + t663) * t410) * t372) * t613 + (-t436 * t585 - (t580 * t685 + (0.2e1 * t518 - t664) * t410) * t372) * t612 + (t436 * t531 + (t400 * t564 + t493 * t436) * t410) * t611 + (t437 * t531 + (t401 * t564 + t493 * t437) * t410) * t610 + (((-qJD(1) * t388 - t386) * t607 + (qJD(1) * t364 + qJDD(1) * t388 + t362) * t467) * MDP(17) + (-0.2e1 * t385 * t607 + (-qJD(1) * t363 + qJDD(1) * t387 - t361) * t467) * MDP(18)) * pkin(2) + (-t465 * MDP(10) - t408 * MDP(17) + t407 * MDP(18) + t467 * MDP(9) + MDP(2)) * t535; t465 * qJDD(1) * MDP(6) + MDP(7) * t604 + (qJDD(2) * MDP(8)) + (-g(3) * t467 + t536 * t465) * MDP(9) + (g(3) * t465 + t536 * t467) * MDP(10) + t385 * t365 * MDP(11) + (-t365 * t386 + t366 * t385) * MDP(12) + (t361 * t374 - t365 * t371 + t626 * t385) * MDP(13) - t386 * t366 * MDP(14) + (t362 * t374 - t366 * t371 - t626 * t386) * MDP(15) + (t360 * t374 - t626 * t371) * MDP(16) + ((g(3) * t408 - t536 * t407 + 0.2e1 * t435 * t505) * t374 + ((-t366 * t467 + t386 * t465) * qJD(1) + (t435 * t709 + (-t399 * t371 - t553 * t435 + (t371 * t403 - t374 * t399 + t626 * t435) * qJD(2)) * t412) * t638) * pkin(2)) * MDP(17) + ((-g(3) * t407 - t536 * t408 + t505 * t693) * t374 + ((t365 * t467 + t385 * t465) * qJD(1) + (-t438 * t709 + (t402 * t371 + t553 * t438 + (-t371 * t405 + t374 * t402 - t626 * t438) * qJD(2)) * t412) * t638) * pkin(2)) * MDP(18) + (-t437 * t542 + (-t510 * t579 + ((t372 * t406 - t373 * t401) * t456 - t696 * t437) * t410) * qJD(1)) * t613 + (t436 * t542 + (t510 * t580 + ((-t372 * t404 + t373 * t400) * t456 + t696 * t436) * t410) * qJD(1)) * t612 + (-t359 * t373 + t627 * t372) * MDP(23) + ((t679 / 0.2e1 + t500 * t436) * t584 + (t436 * t520 + (-g(3) * t406 + t536 * t404) * t456 + (t404 * t639 + (-t436 * t507 + 0.2e1 * t679) * t592) * pkin(1)) * t410) * t611 + ((-t680 / 0.2e1 + t500 * t437) * t584 + (t437 * t520 + (g(3) * t404 + t536 * t406) * t456 + (t406 * t639 + (-t437 * t507 - 0.2e1 * t680) * t592) * pkin(1)) * t410) * t610 + (-t633 * MDP(4) + t617 * MDP(5) + (t578 * t701 + (-t406 * t651 + 0.2e1 * t433 * t565) * t414) * t615 + (-t377 * t549 + (t404 * t651 + (t406 * t457 + t458 * t547) * t436) * t414) * t614) * t471;];
tau = t1;
