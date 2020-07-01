% Calculate vector of inverse dynamics joint torques for
% fourbar1turnDE1
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
%   see fourbar1turnDE1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(5,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'fourbar1turnDE1_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:59
% EndTime: 2020-06-27 16:36:25
% DurationCPUTime: 13.28s
% Computational Cost: add. (110603->338), mult. (158189->747), div. (5768->21), fcn. (42667->10), ass. (0->293)
t500 = sin(qJ(2));
t514 = pkin(1) ^ 2;
t502 = cos(qJ(2));
t736 = pkin(2) * t502;
t762 = -2 * pkin(1);
t676 = t736 * t762 + t514;
t754 = -pkin(3) - pkin(4);
t487 = (pkin(2) - t754) * (pkin(2) + t754) + t676;
t753 = pkin(4) - pkin(3);
t488 = (pkin(2) - t753) * (pkin(2) + t753) + t676;
t788 = pkin(1) * pkin(2);
t606 = (-t487 - t488) * t788;
t477 = t500 * t606;
t513 = pkin(2) ^ 2;
t494 = t513 + t676;
t775 = pkin(3) ^ 2;
t776 = pkin(4) ^ 2;
t674 = t775 - t776;
t489 = t494 + t674;
t706 = t487 * t488;
t517 = sqrt(-t706);
t694 = t500 * t517;
t481 = 0.1e1 / t517;
t496 = pkin(1) * t502 - pkin(2);
t708 = t481 * t496;
t498 = t500 ^ 2;
t698 = t498 * t514;
t787 = 0.2e1 * pkin(2) * t698;
t443 = -t477 * t708 + t787 + (t489 * t502 + t694) * pkin(1);
t492 = 0.1e1 / t494 ^ 2;
t512 = 0.1e1 / t775;
t470 = -pkin(1) * t694 - t489 * t496;
t460 = t470 ^ 2;
t743 = pkin(1) * t489;
t486 = t500 * t743;
t473 = -t496 * t517 + t486;
t469 = t473 ^ 2;
t681 = t460 + t469;
t452 = t681 * t512 * t492;
t447 = t452 ^ (-0.1e1 / 0.2e1);
t491 = 0.1e1 / t494;
t511 = 0.1e1 / pkin(3);
t658 = -0.2e1 * pkin(2) * t496;
t692 = t502 * t517;
t710 = t481 * t477;
t441 = t486 + (-t692 + (t658 - t710) * t500) * pkin(1);
t493 = t491 * t492;
t737 = pkin(2) * t500;
t656 = pkin(1) * t737;
t621 = t493 * t656;
t578 = 0.4e1 * t621;
t544 = t681 * t578;
t655 = 0.2e1 * t492;
t413 = ((t441 * t470 + t443 * t473) * t655 - t544) * t512;
t712 = t473 * t502;
t717 = t470 * t500;
t448 = t447 / t452;
t723 = t448 * t491;
t539 = (-t712 / 0.2e1 - t717 / 0.2e1) * t723;
t534 = t413 * t539;
t683 = t441 - t473;
t619 = t683 * t500;
t605 = t655 * t788;
t696 = t500 * t502;
t769 = (t470 * t498 + t473 * t696) * t605;
t528 = (t447 * (t491 * (t443 * t502 + t619) - t769) + t534) * t511;
t505 = qJD(2) ^ 2;
t690 = t505 * t513;
t785 = 0.8e1 * t493;
t781 = 2 * qJD(2);
t780 = 0.6e1 * t505 * t696;
t508 = 0.1e1 / pkin(4);
t490 = t494 - t674;
t495 = pkin(1) - t736;
t471 = -pkin(2) * t694 + t490 * t495;
t666 = qJD(2) * t500;
t642 = pkin(2) * t666;
t618 = pkin(1) * t642;
t591 = t492 * t618;
t567 = t471 * t591;
t476 = qJD(2) * t477;
t711 = t476 * t481;
t629 = t500 * t711;
t657 = 0.2e1 * t495 * pkin(1);
t438 = (-t629 + (-t692 + (t490 + t657) * t500) * qJD(2)) * pkin(2);
t729 = t438 * t491;
t424 = (-t729 / 0.2e1 + t567) * t508;
t485 = t490 * t737;
t472 = t495 * t517 + t485;
t566 = t472 * t591;
t626 = t517 * t666;
t699 = t498 * t513;
t646 = pkin(1) * t699;
t665 = qJD(2) * t502;
t709 = t481 * t495;
t439 = t646 * t781 + t476 * t709 + (t490 * t665 + t626) * pkin(2);
t728 = t439 * t491;
t426 = (t728 / 0.2e1 - t566) * t508;
t465 = 0.1e1 / t471;
t466 = 0.1e1 / t471 ^ 2;
t719 = t466 * t472;
t558 = t424 * t719 + t426 * t465;
t778 = t712 + t717;
t623 = t489 + t658;
t437 = (-t629 + (t500 * t623 - t692) * qJD(2)) * pkin(1);
t440 = pkin(1) * t626 + qJD(2) * t787 - t476 * t708 + t665 * t743;
t461 = 0.1e1 / t470;
t463 = t461 / t460;
t720 = t463 * t469;
t462 = 0.1e1 / t470 ^ 2;
t721 = t462 * t473;
t415 = -t437 * t720 + t440 * t721;
t574 = -0.2e1 * t591;
t425 = (t437 * t491 + t470 * t574) * t511;
t427 = (t440 * t491 + t473 * t574) * t511;
t458 = t462 * t469 + 0.1e1;
t455 = 0.1e1 / t458;
t456 = 0.1e1 / t458 ^ 2;
t568 = qJDD(2) * t500 + t502 * t505;
t548 = (t476 ^ 2 / t706 + t568 * t606 - 0.4e1 * t698 * t690) * t481;
t688 = t517 * t505;
t536 = -t548 + t688;
t551 = qJDD(2) * t517 + t711 * t781;
t735 = pkin(3) * t494;
t649 = t455 * t735;
t613 = t462 * t649;
t584 = t473 * t613;
t614 = t461 * t649;
t630 = t498 * t690;
t722 = t456 * t494;
t639 = t415 * t722;
t738 = pkin(2) * t492;
t659 = -0.2e1 * t738;
t662 = qJDD(2) * t489;
t667 = qJD(2) * t492;
t689 = t505 * t514;
t695 = t500 * t511;
t705 = t491 * t496;
t730 = t427 * t461;
t739 = pkin(2) * t491;
t759 = 0.2e1 * t473;
t760 = -0.2e1 * t473;
t766 = -t425 * t721 + t730;
t773 = 0.2e1 * qJDD(2);
t395 = qJDD(2) + (-t548 * t705 + (t473 * t630 * t785 + (t498 * t773 + t780) * t739) * t514 + (((t662 + t688) * t491 + t473 * t505 * t659) * t502 + ((-t489 * t505 + t551) * t491 + (-0.4e1 * qJD(2) * t440 + qJDD(2) * t760) * t738) * t500) * pkin(1)) * t511 * t614 - ((t470 * t513 * t785 + 0.4e1 * t739) * t511 * t498 * t689 + ((-t551 * t491 + (t470 * t659 + t491 * t623) * t505) * t511 * t502 + ((t536 + t662) * t491 + (-0.4e1 * t437 * t667 + (-t470 * t492 - t705) * t773) * pkin(2)) * t695) * pkin(1)) * t584 + (-t425 * t440 - t427 * t437) * t613 + (t415 * t456 * t462 + t437 * t455 * t463) * t425 * t735 * t759 + 0.2e1 * (t455 * t618 * t766 - t639 * t730) * pkin(3);
t406 = -t425 * t584 + t427 * t614 + qJD(2);
t650 = t492 * t737;
t622 = pkin(1) * t650;
t599 = -0.2e1 * t622;
t429 = (t441 * t491 + t470 * t599) * t511;
t431 = (t443 * t491 + t473 * t599) * t511;
t409 = -t429 * t584 + t431 * t614 + 0.1e1;
t411 = ((t437 * t470 + t440 * t473) * t655 - qJD(2) * t544) * t512;
t702 = t491 * t511;
t651 = pkin(2) * t702;
t610 = qJD(2) * t651;
t580 = t447 * t610;
t634 = t447 * t702;
t615 = pkin(2) * t634;
t555 = t429 * t721 - t431 * t461;
t604 = 0.2e1 * t656;
t768 = (-0.2e1 * t766 * (-t441 * t720 + t443 * t721) * t722 - 0.2e1 * t555 * t639 + ((qJD(2) * t555 + t766) * t604 + ((t425 * t441 - t429 * t437) * t463 * t759 + (-t425 * t443 - t427 * t441 + t429 * t440 + t431 * t437) * t462) * t494) * t455) * pkin(3);
t777 = -t768 * t580 + (qJDD(2) * t409 + t395) * t615 + (-t411 * t651 * (qJD(2) * t409 + t406) / 0.2e1 + t406 * t413 * t610 / 0.2e1) * t448;
t774 = -0.2e1 * t439;
t442 = t485 + (-t692 + (t657 - t710) * t500) * pkin(2);
t444 = t477 * t709 + 0.2e1 * t646 + (t490 * t502 + t694) * pkin(2);
t509 = 0.1e1 / t776;
t464 = t471 ^ 2;
t468 = t472 ^ 2;
t680 = t464 + t468;
t543 = t680 * t578;
t412 = ((t442 * t471 + t444 * t472) * t655 - t543) * t509;
t772 = t412 / 0.2e1;
t467 = t465 / t464;
t718 = t467 * t468;
t414 = -t438 * t718 + t439 * t719;
t457 = t466 * t468 + 0.1e1;
t454 = 0.1e1 / t457 ^ 2;
t771 = t414 * t454;
t770 = t491 * t778;
t767 = t406 * t615 + t409 * t580;
t453 = 0.1e1 / t457;
t734 = pkin(4) * t494;
t648 = t453 * t734;
t596 = 0.2e1 * t648;
t407 = t558 * t596;
t751 = -t491 / 0.2e1;
t428 = (t442 * t751 + t471 * t622) * t508;
t750 = t491 / 0.2e1;
t430 = (t444 * t750 - t472 * t622) * t508;
t556 = t428 * t719 + t430 * t465;
t408 = t556 * t596;
t597 = -0.4e1 * t618;
t628 = -0.4e1 * t656;
t600 = t472 * t628;
t714 = t472 * t494;
t617 = 0.4e1 * t467 * t714;
t757 = 0.2e1 * t494;
t653 = t466 * t757;
t758 = -0.2e1 * t494;
t687 = (0.2e1 * (t558 * (-t442 * t718 + t444 * t719) - t556 * t414) * t454 * t757 + ((t442 * t653 + t465 * t628) * t426 + (t442 * t617 + (t444 * t758 + t600) * t466) * t424 - (t438 * t653 + t465 * t597) * t430 - (t438 * t617 + (t439 * t758 + t472 * t597) * t466) * t428) * t453) * pkin(4);
t764 = t492 * t604 * (-qJD(2) * t408 + t407) + t491 * t687;
t451 = t680 * t509 * t492;
t449 = 0.1e1 / t451;
t445 = t451 ^ (-0.1e1 / 0.2e1);
t410 = ((t438 * t471 + t439 * t472) * t655 - qJD(2) * t543) * t509;
t752 = t410 / 0.2e1;
t501 = sin(qJ(1));
t749 = g(1) * t501;
t503 = cos(qJ(1));
t748 = g(1) * t503;
t747 = g(2) * t501;
t746 = g(2) * t503;
t745 = g(3) * t471;
t744 = g(3) * t472;
t742 = pkin(1) * t491;
t741 = pkin(1) * t492;
t733 = pkin(2) * qJD(1);
t725 = t445 * t491;
t446 = t445 * t449;
t724 = t446 * t491;
t716 = t470 * t502;
t715 = t472 * t492;
t713 = t473 * t500;
t506 = qJD(1) ^ 2;
t703 = t491 * t506;
t691 = t503 * t511;
t684 = t778 * t503 * t634;
t682 = t443 + t470;
t675 = -t502 ^ 2 + t498;
t673 = MDP(19) * t509;
t672 = MDP(20) * t509;
t671 = MDP(21) * t508;
t670 = MDP(22) * t508;
t669 = MDP(24) * t508;
t668 = MDP(25) * t508;
t542 = t447 * t770;
t422 = t511 * t542;
t420 = qJD(1) * t422;
t664 = qJD(1) * qJD(2);
t663 = qJDD(1) * t502;
t661 = qJDD(2) * t490;
t647 = t465 * t734;
t645 = qJD(1) * t742;
t644 = t500 * t733;
t643 = t502 * t733;
t537 = -t490 * t505 + t551;
t571 = t471 * t492 - t491 * t495;
t612 = t466 * t648;
t394 = -0.2e1 * t424 * t439 * t612 + 0.2e1 * (t438 * t612 + 0.2e1 * t647 * t771) * t426 + 0.2e1 * (-((pkin(1) * t513 * t780 + t548 * t495) * t750 + (0.4e1 * t472 * t493 * t689 + qJDD(2) * t742) * t699 + (((t661 + t688) * t502 + t537 * t500) * t750 + (-t472 * t568 + t666 * t774) * t741) * pkin(2)) * t453 * t647 - ((-0.4e1 * t471 * t493 * t514 - 0.2e1 * t742) * t630 + ((pkin(1) * t505 * t571 + t537 * t750) * t502 + ((t536 + t661) * t751 + (qJDD(2) * t571 + 0.2e1 * t438 * t667) * pkin(1)) * t500) * pkin(2)) * t472 * t612) * t508 + 0.2e1 * (0.2e1 * (t438 * t453 * t467 + t466 * t771) * t424 * t714 - 0.2e1 * t558 * t453 * t618) * pkin(4);
t641 = t394 * t725;
t640 = t412 * t724;
t636 = t471 * t724;
t635 = t472 * t724;
t633 = t447 * t695;
t450 = 0.1e1 / t451 ^ 2;
t632 = t450 * t468 * t492;
t624 = qJDD(1) * t449 * t492;
t620 = -0.2e1 * t645;
t602 = t450 * t471 * t715;
t595 = qJDD(1) * t408 * t725;
t590 = t493 * t618;
t589 = t747 + t748;
t588 = -t746 + t749;
t581 = t410 * t446 * t645;
t577 = t473 * t491 * t633;
t436 = t634 * t716;
t573 = t713 - t716;
t569 = -0.2e1 * t506 * t514 * t650;
t564 = t406 * t580;
t553 = t407 * t772 - t408 * t752;
t550 = 0.2e1 * t589;
t549 = t573 * t491;
t541 = t409 * t633 * t690 * t741;
t540 = -t506 * pkin(1) / 0.2e1 - t748 / 0.2e1 - t747 / 0.2e1;
t538 = (t716 / 0.2e1 - t713 / 0.2e1) * t723;
t535 = t413 * t538;
t531 = (t470 * t696 - t473 * t498) * t605;
t530 = 0.4e1 * qJD(1) * t492 * t514 * t642 + (qJDD(1) * t762 - t588) * t491;
t527 = t511 * (t535 + (t531 + (t500 * t682 - t502 * t683) * t491) * t447);
t526 = -t411 * t539 + ((-t437 * t500 - t440 * t502) * t491 + (t549 + t769) * qJD(2)) * t447;
t525 = t411 * t538 + ((-t437 * t502 + t440 * t500) * t491 + (t531 + t770) * qJD(2)) * t447;
t435 = qJD(1) * t436;
t432 = t501 * t436;
t423 = -t436 + t577;
t421 = qJD(1) * t577 - t435;
t401 = qJD(1) * t527;
t400 = -qJD(1) * t528 - t435;
t399 = t525 * t511;
t398 = t526 * t511;
t397 = (qJDD(1) * t447 * t549 + qJD(1) * t525) * t511;
t396 = (qJD(1) * t526 - qJDD(1) * t542) * t511;
t1 = [qJDD(1) * MDP(1) + t589 * MDP(3) + (qJDD(1) * t498 + 0.2e1 * t664 * t696) * MDP(4) + 0.2e1 * (t500 * t663 - t664 * t675) * MDP(5) + t568 * MDP(6) + (qJDD(2) * t502 - t500 * t505) * MDP(7) + (-t396 * t422 - t398 * t420) * MDP(11) + (t396 * t423 - t397 * t422 + t398 * t421 - t399 * t420) * MDP(12) + (-t395 * t422 + t398 * t406) * MDP(13) + (t397 * t423 + t399 * t421) * MDP(14) + (t395 * t423 + t399 * t406) * MDP(15) + (-g(1) * t432 + (-t573 * t746 + t713 * t749) * t634 + ((-qJD(1) * t423 - t421) * t666 + (qJD(1) * t399 + qJDD(1) * t423 + t397) * t502) * pkin(2)) * MDP(17) + (-g(2) * t684 + t422 * t749 + (-0.2e1 * t420 * t666 + (-qJD(1) * t398 + qJDD(1) * t422 - t396) * t502) * pkin(2)) * MDP(18) + (t468 * t624 + (-t410 * t632 + (0.2e1 * t439 * t715 - 0.4e1 * t468 * t590) * t449) * qJD(1)) * t673 + (-0.2e1 * t471 * t472 * t624 + (0.2e1 * t410 * t602 + (-0.2e1 * t438 * t715 + (0.8e1 * t472 * t590 + t492 * t774) * t471) * t449) * qJD(1)) * t672 + (t472 * t641 - (-t410 * t635 / 0.2e1 + (-0.2e1 * t566 + t728) * t445) * t407) * t671 + (-t471 * t641 - (t636 * t752 + (0.2e1 * t567 - t729) * t445) * t407) * t670 + (t471 * t581 + (t438 * t620 + t471 * t530) * t445) * t669 + (t472 * t581 + (t439 * t620 + t472 * t530) * t445) * t668 + (-t500 * MDP(10) + t502 * MDP(9) + MDP(2)) * t588; t500 * qJDD(1) * MDP(6) + MDP(7) * t663 + qJDD(2) * MDP(8) + (-g(3) * t502 + t500 * t589) * MDP(9) + (g(3) * t500 + t502 * t589) * MDP(10) + t420 * t400 * MDP(11) + (-t400 * t421 + t401 * t420) * MDP(12) + (t396 * t409 - t400 * t406 + t420 * t768) * MDP(13) - t421 * t401 * MDP(14) + (t397 * t409 - t401 * t406 - t421 * t768) * MDP(15) + (t395 * t409 - t406 * t768) * MDP(16) + (t421 * t644 - t401 * t643 + t441 * t564 - g(1) * ((t535 + ((-t502 * t441 + t500 * t443) * t491 + t531) * t447) * t691 + t684) - t527 * t747 - g(3) * (-t436 - t528) - t767 * t437 + (0.2e1 * t541 - t777) * t470) * MDP(17) + (t541 * t760 + t420 * t644 + t400 * t643 - t443 * t564 - g(1) * (t534 + (-t769 + (t502 * t682 + t619) * t491) * t447) * t691 - g(2) * (t501 * t528 + t432) - g(3) * t527 + t767 * t440 + t777 * t473) * MDP(18) + (-t472 * t595 + (-t553 * t635 + ((t407 * t444 - t408 * t439) * t491 - t764 * t472) * t445) * qJD(1)) * t671 + (t471 * t595 + (t553 * t636 + ((-t407 * t442 + t408 * t438) * t491 + t764 * t471) * t445) * qJD(1)) * t670 + (-t394 * t408 + t687 * t407) * MDP(23) + ((t744 / 0.2e1 + t540 * t471) * t640 + (t471 * t569 + (-g(3) * t444 + t442 * t589) * t491 + (t442 * t703 + (-t471 * t550 + 0.2e1 * t744) * t650) * pkin(1)) * t445) * t669 + ((-t745 / 0.2e1 + t540 * t472) * t640 + (t472 * t569 + (g(3) * t442 + t444 * t589) * t491 + (t444 * t703 + (-t472 * t550 - 0.2e1 * t745) * t650) * pkin(1)) * t445) * t668 + (-t696 * MDP(4) + t675 * MDP(5) + (t632 * t772 + (-t444 * t715 + 0.2e1 * t468 * t621) * t449) * t673 + (-t412 * t602 + (t442 * t715 + (t444 * t492 + t493 * t600) * t471) * t449) * t672) * t506;];
tau = t1;
