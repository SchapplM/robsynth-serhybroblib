% Calculate vector of inverse dynamics joint torques for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2OL_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m2OL_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(5,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'palh2m2OL_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:07:11
% EndTime: 2020-06-30 18:07:52
% DurationCPUTime: 11.18s
% Computational Cost: add. (9440->493), mult. (22574->649), div. (0->0), fcn. (19272->18), ass. (0->221)
t576 = sin(qJ(3));
t581 = cos(qJ(3));
t582 = cos(qJ(2));
t681 = qJD(1) * t582;
t577 = sin(qJ(2));
t682 = qJD(1) * t577;
t510 = t576 * t682 - t581 * t681;
t512 = -t576 * t681 - t581 * t682;
t575 = sin(qJ(4));
t580 = cos(qJ(4));
t516 = t576 * t577 - t581 * t582;
t517 = t576 * t582 + t577 * t581;
t569 = qJD(2) + qJD(3);
t746 = t569 * qJD(1);
t596 = t517 * t746;
t589 = -t516 * qJDD(1) - t596;
t484 = t569 * t516;
t617 = t517 * qJDD(1);
t590 = -t484 * qJD(1) + t617;
t677 = qJD(4) * t580;
t678 = qJD(4) * t575;
t421 = -t510 * t677 + t512 * t678 + t575 * t589 + t580 * t590;
t574 = sin(qJ(5));
t628 = t516 * t575 - t517 * t580;
t629 = t580 * t516 + t517 * t575;
t731 = t575 * t510 + t512 * t580;
t586 = t731 * qJD(4) - t629 * qJDD(1) + t628 * t746;
t631 = -t510 * t580 + t575 * t512;
t718 = cos(qJ(5));
t658 = qJD(5) * t718;
t676 = qJD(5) * t574;
t401 = t718 * t421 + t574 * t586 + t631 * t658 + t676 * t731;
t579 = cos(qJ(6));
t399 = t579 * t401;
t564 = qJD(4) + t569;
t549 = qJD(5) + t564;
t573 = sin(qJ(6));
t742 = t574 * t631 - t718 * t731;
t431 = t549 * t579 + t573 * t742;
t568 = qJDD(2) + qJDD(3);
t563 = qJDD(4) + t568;
t547 = qJDD(5) + t563;
t394 = -t431 * qJD(6) - t547 * t573 + t399;
t675 = qJD(6) * t573;
t524 = t549 * t675;
t642 = qJD(6) * t742 + t547;
t395 = t401 * t573 + t642 * t579 - t524;
t402 = qJD(5) * t742 + t574 * t421 - t718 * t586;
t433 = -t573 * t549 + t579 * t742;
t441 = t574 * t731 + t631 * t718;
t400 = -qJDD(6) + t402;
t766 = qJD(6) + t441;
t655 = t400 * t579 + t675 * t766;
t708 = t433 * t579;
t710 = t400 * t573;
t711 = t394 * t573;
t751 = t766 * t431;
t771 = t547 * MDP(29) + ((-t394 + t751) * t579 + t433 * t675 + (t433 * t441 + t395) * t573) * MDP(33) - t402 * MDP(28) - t441 ^ 2 * MDP(26) + (-t441 * t549 + t401) * MDP(27) + (t441 * t573 * t766 + t655) * MDP(35) + (-t441 * MDP(25) + MDP(26) * t742 + t549 * MDP(28) - t431 * MDP(35) + MDP(36) * t766) * t742 + (-t579 * t766 ^ 2 + t433 * t742 + t710) * MDP(34) + (-t708 * t766 - t711) * MDP(32);
t770 = t563 * MDP(22) + t631 * MDP(18) * t731 + (-t564 * t631 + t421) * MDP(20) + (-t631 ^ 2 + t731 ^ 2) * MDP(19) + (-t564 * t731 + t586) * MDP(21) + t771;
t712 = pkin(4) * qJD(2);
t518 = pkin(2) * t569 + t581 * t712;
t664 = t576 * t712;
t639 = t575 * t664;
t488 = t580 * t518 - t639;
t483 = pkin(5) * t564 + t488;
t489 = t518 * t575 + t580 * t664;
t698 = t574 * t489;
t454 = t718 * t483 - t698;
t449 = pkin(3) * t549 + t454;
t769 = t766 * t449;
t717 = pkin(4) * t581;
t559 = qJDD(2) * t717;
t499 = pkin(2) * t568 - qJD(3) * t664 + t559;
t491 = t580 * t499;
t659 = t576 * t677;
t660 = t518 * t678;
t666 = qJDD(2) * t576;
t680 = qJD(3) * t581;
t591 = t491 + (-t575 * t666 + (-t575 * t680 - t659) * qJD(2)) * pkin(4) - t660;
t444 = pkin(5) * t563 + t591;
t651 = -qJD(4) * t639 + t499 * t575;
t656 = qJD(2) * t680;
t612 = (t656 + t666) * pkin(4);
t760 = t580 * (qJD(4) * t518 + t612);
t451 = t651 + t760;
t662 = t718 * t489;
t455 = t574 * t483 + t662;
t606 = -t455 * qJD(5) + t718 * t444 - t574 * t451;
t407 = t547 * pkin(3) + t606;
t572 = qJ(2) + qJ(3);
t567 = qJ(4) + t572;
t560 = qJ(5) + t567;
t545 = sin(t560);
t578 = sin(qJ(1));
t583 = cos(qJ(1));
t635 = g(1) * t583 + g(2) * t578;
t750 = t635 * t545;
t768 = t407 + t750;
t553 = -t582 * pkin(4) - pkin(1);
t529 = qJD(1) * t553;
t487 = pkin(2) * t510 + t529;
t452 = -pkin(5) * t631 + t487;
t546 = cos(t560);
t411 = t574 * t444 + t718 * t451 + t483 * t658 - t489 * t676;
t539 = g(3) * t545;
t732 = -t411 + t539;
t765 = -t452 * t441 + t546 * t635 + t732;
t416 = -pkin(3) * t441 + t452;
t408 = -t416 * t579 - t455 * t573;
t625 = t408 * t742 + t768 * t579;
t714 = g(3) * t546;
t764 = -t573 * t769 - t579 * t714 + t625;
t557 = sin(t567);
t558 = cos(t567);
t701 = t558 * t583;
t702 = t558 * t578;
t758 = g(1) * t701 + g(2) * t702 + g(3) * t557 - t631 * t487 - t651;
t686 = pkin(3) * t742;
t409 = -t416 * t573 + t455 * t579;
t757 = -t409 * t742 + t573 * t714;
t755 = -t452 * t742 + t606 - t714 + t750;
t752 = pkin(5) * t731;
t749 = t758 - t760;
t743 = -g(3) * t558 + t487 * t731 + t635 * t557;
t673 = t512 * pkin(2);
t456 = -t673 - t752;
t422 = t456 + t686;
t551 = pkin(2) * t580 + pkin(5);
t661 = t718 * t575;
t504 = pkin(2) * t661 + t574 * t551;
t734 = (qJD(6) * t504 - t422) * t766;
t552 = pkin(2) + t717;
t533 = t580 * t552;
t695 = t575 * t576;
t503 = -pkin(4) * t695 + pkin(5) + t533;
t694 = t576 * t580;
t508 = pkin(4) * t694 + t552 * t575;
t465 = t574 * t503 + t718 * t508;
t561 = pkin(4) * t682;
t733 = (qJD(6) * t465 - t422 - t561) * t766;
t652 = t718 * t503 - t574 * t508;
t493 = pkin(2) * t516 + t553;
t485 = t569 * t517;
t723 = t628 * qJD(4) + t484 * t575 - t485 * t580;
t446 = -t574 * t629 - t628 * t718;
t705 = t446 * t579;
t700 = t573 * t578;
t699 = t573 * t583;
t696 = t574 * t575;
t693 = t578 * t579;
t692 = t579 * t583;
t627 = -t575 * t581 - t694;
t506 = t627 * t712;
t626 = t580 * t581 - t695;
t507 = t626 * t712;
t691 = t574 * t506 + t718 * t507 - t551 * t658 - (-t575 * t676 + (t718 * t580 - t696) * qJD(4)) * pkin(2);
t690 = -t718 * t506 + t574 * t507 - t551 * t676 + (-t575 * t658 + (-t574 * t580 - t661) * qJD(4)) * pkin(2);
t570 = t577 ^ 2;
t688 = -t582 ^ 2 + t570;
t674 = qJD(6) * t579;
t669 = qJD(1) * qJD(2);
t668 = qJDD(1) * t553;
t667 = qJDD(1) * t582;
t665 = qJDD(1) * pkin(1);
t562 = t577 * t712;
t657 = t577 * t669;
t471 = t485 * pkin(2) + t562;
t548 = pkin(4) * t657;
t410 = -t589 * pkin(2) - pkin(4) * t667 - t586 * pkin(5) + t548 - t665;
t646 = t402 * pkin(3) + qJD(6) * t455 + t410;
t645 = -qJD(6) * t416 + t411;
t641 = qJD(2) * (-qJD(3) + t569);
t640 = qJD(3) * (-qJD(2) - t569);
t638 = -0.2e1 * pkin(1) * t669;
t637 = t766 * t658;
t636 = qJD(2) * t659;
t634 = g(1) * t578 - g(2) * t583;
t633 = -pkin(2) * t696 + t718 * t551;
t565 = sin(t572);
t566 = cos(t572);
t624 = g(3) * t565 + t510 * t529 + t566 * t635;
t622 = t574 * t628 - t629 * t718;
t429 = -t629 * qJD(4) - t484 * t580 - t485 * t575;
t404 = t622 * qJD(5) + t718 * t429 + t574 * t723;
t621 = -t404 * t579 + t446 * t675;
t585 = qJD(1) ^ 2;
t616 = pkin(1) * t585 + t635;
t615 = t634 + 0.2e1 * t665;
t613 = -g(3) * t566 + t512 * t529 + t565 * t635 + t559;
t459 = pkin(5) * t629 + t493;
t424 = -pkin(3) * t622 + t459;
t611 = qJD(6) * t424 * t766 + t449 * t404 + t407 * t446;
t477 = t552 * t677 + (t626 * qJD(3) - t576 * t678) * pkin(4);
t478 = -t552 * t678 + (t627 * qJD(3) - t659) * pkin(4);
t426 = qJD(5) * t652 + t718 * t477 + t574 * t478;
t610 = t400 * t465 - t426 * t766 - t769;
t609 = t504 * t400 + t691 * t766 - t769;
t608 = t491 + t743;
t601 = (-pkin(2) * t564 - t518) * qJD(4) - t612;
t405 = t446 * qJD(5) + t574 * t429 - t718 * t723;
t423 = -pkin(5) * t723 + t471;
t600 = t449 * qJD(6) * t446 - (t405 * pkin(3) + t423) * t766 + t424 * t400 - t646 * t622;
t588 = -t512 * t510 * MDP(11) + (t510 * t569 + t590) * MDP(13) + (-t512 * t569 + t589) * MDP(14) + (-t510 ^ 2 + t512 ^ 2) * MDP(12) + t568 * MDP(15) + t770;
t584 = qJD(2) ^ 2;
t550 = t718 * pkin(5) + pkin(3);
t505 = t548 + t668;
t500 = pkin(3) + t633;
t498 = t546 * t692 - t700;
t497 = -t546 * t699 - t693;
t496 = -t546 * t693 - t699;
t495 = t546 * t700 - t692;
t490 = t561 - t673;
t463 = pkin(3) + t652;
t458 = t718 * t488 - t698;
t457 = -t574 * t488 - t662;
t453 = t561 + t456;
t450 = pkin(2) * t596 + qJDD(1) * t493 + t548;
t427 = -t465 * qJD(5) - t574 * t477 + t718 * t478;
t425 = t686 - t752;
t414 = t416 * t675;
t1 = [t634 * MDP(2) + t635 * MDP(3) + (g(1) * t702 - g(2) * t701 + t450 * t629 - t471 * t631 - t487 * t723 - t493 * t586) * MDP(23) + (t401 * t622 - t402 * t446 + t404 * t441 - t405 * t742) * MDP(26) + (t402 * t459 + t405 * t452 - t410 * t622 - t423 * t441 + t634 * t546) * MDP(30) + 0.2e1 * (t577 * t667 - t688 * t669) * MDP(5) + (-t615 * t577 + t582 * t638) * MDP(10) + (t577 * t638 + t615 * t582) * MDP(9) + ((-t431 * t579 - t433 * t573) * t404 + (-t711 - t395 * t579 + (t431 * t573 - t708) * qJD(6)) * t446) * MDP(33) + (t446 * t710 - t395 * t622 + t405 * t431 - (t404 * t573 + t446 * t674) * t766) * MDP(35) + (qJDD(1) * t570 + 0.2e1 * t582 * t657) * MDP(4) + (-t400 * t622 - t405 * t766) * MDP(36) + (t394 * t622 - t400 * t705 - t405 * t433 - t621 * t766) * MDP(34) + qJDD(1) * MDP(1) + (-t563 * t629 + t564 * t723) * MDP(21) + (t401 * t459 + t404 * t452 + t410 * t446 + t423 * t742 - t634 * t545) * MDP(31) + (t401 * t446 + t404 * t742) * MDP(25) + (-0.2e1 * t516 * t617 + t484 * t510 + t512 * t485 + (t516 ^ 2 - t517 ^ 2) * t746) * MDP(12) + (t394 * t705 - t621 * t433) * MDP(32) + (-t405 * t549 + t547 * t622) * MDP(28) + (-g(1) * t496 - g(2) * t498 - t408 * t405 + t414 * t622 + (-t411 * t622 + t611) * t573 + t600 * t579) * MDP(37) + (-g(1) * t495 - g(2) * t497 + t409 * t405 + (-t622 * t645 + t611) * t579 - t600 * t573) * MDP(38) + (t429 * t564 - t563 * t628) * MDP(20) + (-t421 * t629 + t429 * t631 - t586 * t628 - t723 * t731) * MDP(19) + (t421 * t493 + t429 * t487 - t450 * t628 - t471 * t731 - t634 * t557) * MDP(24) + (-t421 * t628 - t429 * t731) * MDP(18) + (-t529 * t484 + t505 * t517 - t512 * t562 + t553 * t590 - t634 * t565) * MDP(17) + (t512 * t484 + t517 * t590) * MDP(11) + (-t484 * t569 + t517 * t568) * MDP(13) + (t404 * t549 + t446 * t547) * MDP(27) + (-t485 * t569 - t516 * t568) * MDP(14) + (qJDD(2) * t577 + t582 * t584) * MDP(6) + (qJDD(2) * t582 - t577 * t584) * MDP(7) + (0.2e1 * t529 * t485 + t510 * t562 + t634 * t566 + (t505 + t668) * t516) * MDP(16); t588 + (t427 * t549 + t441 * t453 + t652 * t547 + t755) * MDP(30) + (-t477 * t564 + t490 * t731 - t508 * t563 + t749) * MDP(24) + (t395 * t463 + t427 * t431 + (-t714 - t733) * t579 + t610 * t573 + t625) * MDP(37) + ((t512 * t682 + (-qJDD(2) - t568) * t576 + t581 * t640) * pkin(4) + t624) * MDP(17) + ((-t510 * t682 + t568 * t581 + t576 * t640) * pkin(4) + t613) * MDP(16) + t688 * MDP(5) * t585 + qJDD(2) * MDP(8) + MDP(7) * t667 + (t394 * t463 + t427 * t433 + t610 * t579 + (-t768 + t733) * t573 + t757) * MDP(38) + (t478 * t564 + t533 * t563 - t660 + t490 * t631 + (-t636 + (-t656 + (-qJDD(2) - t563) * t576) * t575) * pkin(4) + t608) * MDP(23) + (-t426 * t549 - t453 * t742 - t465 * t547 + t765) * MDP(31) + (t616 * MDP(10) - g(3) * MDP(9)) * t582 + (-t585 * t582 * MDP(4) + g(3) * MDP(10) + qJDD(1) * MDP(6) + t616 * MDP(9)) * t577; t588 + (t507 * t564 + (-t512 * t731 - t563 * t575) * pkin(2) + t601 * t580 + t758) * MDP(24) + (t500 * t395 + t690 * t431 + (-t714 - t734) * t579 + t609 * t573 + t625) * MDP(37) + (t500 * t394 + t690 * t433 + t609 * t579 + (-t768 + t734) * t573 + t757) * MDP(38) + (t441 * t456 + t633 * t547 + t690 * t549 + t755) * MDP(30) + (-t456 * t742 - t504 * t547 + t691 * t549 + t765) * MDP(31) + ((t581 * t641 - t666) * pkin(4) + t624) * MDP(17) + (t576 * pkin(4) * t641 + t613) * MDP(16) + (-pkin(4) * t636 - t506 * t564 + (-t512 * t631 + t563 * t580) * pkin(2) + t601 * t575 + t608) * MDP(23); (t489 * t564 + t591 + t743) * MDP(23) + (t488 * t564 + t749) * MDP(24) + (-t457 * t549 + (-t441 * t731 + t718 * t547 - t549 * t676) * pkin(5) + t755) * MDP(30) + (t458 * t549 + (-t547 * t574 - t549 * t658 + t731 * t742) * pkin(5) + t765) * MDP(31) + (t550 * t395 - (-t425 * t579 - t458 * t573) * t766 - t457 * t431 + (-t573 * t637 + (-qJD(5) * t431 - t674 * t766 + t710) * t574) * pkin(5) + t764) * MDP(37) + (t550 * t394 - t457 * t433 + (t458 * t766 - t769) * t579 + (-t425 * t766 - t768) * t573 + (-t579 * t637 + (-qJD(5) * t433 + t655) * t574) * pkin(5) + t757) * MDP(38) + t770; (t455 * t549 + t755) * MDP(30) + (t454 * t549 + t765) * MDP(31) + (pkin(3) * t395 - (-t454 * t573 - t579 * t686) * t766 + t455 * t431 + t764) * MDP(37) + (pkin(3) * t394 + t455 * t433 + (t454 * t766 - t769) * t579 + (-t686 * t766 - t768) * t573 + t757) * MDP(38) + t771; t433 * t431 * MDP(32) + (-t431 ^ 2 + t433 ^ 2) * MDP(33) + (t399 + t751) * MDP(34) + (t433 * t766 + t524) * MDP(35) - t400 * MDP(36) + (-g(1) * t497 + g(2) * t495 + t409 * t766 - t433 * t449 + t414) * MDP(37) + (g(1) * t498 - g(2) * t496 + t408 * t766 + t431 * t449) * MDP(38) + (-t642 * MDP(34) - t401 * MDP(35) + MDP(37) * t732 + t646 * MDP(38)) * t573 + (-qJD(6) * t549 * MDP(34) - t642 * MDP(35) - t646 * MDP(37) + (-t645 + t539) * MDP(38)) * t579;];
tau = t1;
