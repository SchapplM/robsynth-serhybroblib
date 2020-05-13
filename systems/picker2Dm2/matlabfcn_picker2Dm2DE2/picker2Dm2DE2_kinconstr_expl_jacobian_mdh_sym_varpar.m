% Jacobian of explicit kinematic constraints of
% picker2Dm2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% W [15x2]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = picker2Dm2DE2_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 19:11:06
% EndTime: 2020-05-09 19:11:22
% DurationCPUTime: 14.16s
% Computational Cost: add. (94437->647), mult. (293325->1080), div. (1348->20), fcn. (47053->10), ass. (0->466)
t522 = sin(qJ(2));
t795 = pkin(7) * t522;
t706 = 0.2e1 * t795;
t456 = pkin(3) * t706;
t564 = pkin(1) ^ 2;
t554 = pkin(4) ^ 2;
t566 = pkin(7) ^ 2;
t707 = t566 - t554;
t464 = t564 + t707;
t418 = t456 + t464;
t559 = pkin(3) ^ 2;
t538 = 0.2e1 * t559;
t416 = t538 + t418;
t477 = pkin(3) * t522;
t454 = t477 + pkin(7);
t526 = cos(qJ(1));
t493 = t526 ^ 2;
t800 = 0.4e1 * t493;
t358 = t416 * t526 + (t800 - 0.2e1) * t454 * pkin(1);
t433 = t454 * t526;
t523 = sin(qJ(1));
t525 = cos(qJ(2));
t742 = t523 * t525;
t450 = pkin(3) * t742;
t723 = -t450 + t433;
t391 = pkin(1) + t723;
t478 = pkin(3) * t525;
t470 = pkin(7) * t478;
t457 = 0.2e1 * t470;
t740 = t525 * t559;
t799 = 0.4e1 * t522;
t614 = t740 * t799;
t414 = t457 + t614;
t617 = pkin(1) * t450;
t437 = -0.2e1 * t617;
t387 = t437 + t418;
t710 = t564 - t566;
t461 = t710 * t559;
t489 = t522 ^ 2;
t552 = 0.2e1 * pkin(3);
t562 = t564 ^ 2;
t749 = t489 * t559;
t683 = 0.2e1 * t749;
t712 = -t559 + t566;
t408 = t456 + t683 + t712;
t801 = -0.2e1 * t493;
t601 = t408 * t801 - t707;
t658 = -0.4e1 * pkin(3) * pkin(7) * t464;
t694 = 0.2e1 * t564;
t821 = 0.4e1 * pkin(1);
t335 = 0.4e1 * t461 * t489 + t522 * t658 - t562 - (t566 - (t552 + pkin(4)) * pkin(4)) * (t566 + (t552 - pkin(4)) * pkin(4)) + (-0.2e1 * t559 + t601) * t694 + (-t387 * t433 + t418 * t450) * t821;
t567 = sqrt(t335);
t732 = t564 * t493;
t677 = -0.4e1 * t732;
t794 = pkin(7) * t523;
t690 = t559 * t794;
t745 = t522 * t523;
t667 = pkin(3) * t745;
t442 = pkin(1) * t667;
t724 = t442 + t470;
t741 = t525 * t526;
t820 = 0.4e1 * pkin(3);
t332 = t414 * t677 + (0.8e1 * t461 * t522 + t658) * t525 + (-0.8e1 * t724 * t433 + 0.8e1 * t525 ^ 2 * t690 + (-t387 * t741 - t418 * t745) * t820) * pkin(1);
t334 = 0.1e1 / t567;
t768 = t332 * t334;
t633 = t768 / 0.2e1;
t695 = 0.2e1 * t478;
t479 = pkin(1) * t526;
t700 = 0.2e1 * t479;
t737 = t526 * t567;
t793 = pkin(7) * t526;
t812 = 0.2e1 * pkin(7);
t822 = 0.2e1 * pkin(1);
t316 = t391 * t633 + t414 * t523 * t700 + ((t523 * t567 - t358) * t522 + (t737 + (t454 * t812 + t416) * t523 + (t493 * t822 - pkin(1) + t793) * t695) * t525) * pkin(3);
t543 = 0.3e1 * t564;
t447 = t538 + t543 + t707;
t608 = -0.4e1 * t617;
t384 = t447 + t456 + t608;
t666 = pkin(3) * t741;
t758 = t454 * t523;
t393 = t666 + t758;
t435 = t450 - pkin(1);
t317 = t450 * t567 + t393 * t633 - (0.4e1 * t442 + t457) * t433 - 0.2e1 * t489 * t690 + (t414 * t801 + t614) * pkin(1) + ((-t384 * t526 + t435 * t812) * t525 + (-t447 * t523 - t737) * t522) * pkin(3);
t484 = t564 + t566;
t639 = t559 + t484;
t671 = pkin(1) * t433;
t365 = t437 + t456 + t639 + 0.2e1 * t671;
t363 = 0.1e1 / t365;
t382 = pkin(1) * t666 + t724;
t560 = 0.1e1 / pkin(3);
t762 = t393 * t567;
t325 = -t384 * t433 + t762 + (t435 * t706 + t447 * t742) * pkin(3) + (-t564 + (0.2e1 * t489 - 0.4e1) * t559 + t601) * pkin(1);
t672 = t408 * t479;
t355 = t416 * t454 + 0.2e1 * t672;
t326 = t355 * t523 + t358 * t478 + t391 * t567;
t520 = cos(pkin(8));
t776 = sin(pkin(8));
t406 = t520 * t526 + t776 * t523;
t783 = t406 / 0.2e1;
t404 = t520 * t523 - t776 * t526;
t784 = t404 / 0.2e1;
t595 = t325 * t784 + t326 * t783;
t364 = 0.1e1 / t365 ^ 2;
t819 = -0.2e1 * t364;
t588 = t595 * t819;
t293 = ((t316 * t783 + t317 * t784) * t363 + t382 * t588) * t560;
t769 = t326 * t404;
t770 = t325 * t406;
t587 = 0.2e1 * t364 * (-t770 / 0.2e1 + t769 / 0.2e1);
t294 = ((t317 * t783 - t316 * t404 / 0.2e1) * t363 + t382 * t587) * t560;
t524 = sin(pkin(9));
t527 = cos(pkin(9));
t766 = t363 * t560;
t314 = t595 * t766;
t322 = -t766 * t769 / 0.2e1;
t315 = t766 * t770 / 0.2e1 + t322;
t306 = t314 * t524 + t315 * t527;
t809 = 0.1e1 / t306;
t305 = -t314 * t527 + t315 * t524;
t808 = 0.1e1 / t306 ^ 2;
t813 = 0.1e1 / (t305 ^ 2 * t808 + 0.1e1);
t816 = t808 * t305;
t827 = ((t293 * t524 + t294 * t527) * t816 - (-t293 * t527 + t294 * t524) * t809) * t813;
t739 = t526 * t523;
t759 = t454 * t493;
t338 = (t408 * t739 + t759 * t478) * t694 + (t387 * t758 + t418 * t666) * pkin(1);
t767 = t334 * t338;
t687 = 0.2e1 * t767;
t744 = t522 * t526;
t691 = pkin(7) * t744;
t319 = t723 * t567 + t393 * t687 + (t384 * t454 + 0.4e1 * t672) * t523 + (t691 * t538 + (t447 * t526 + t759 * t821) * pkin(3)) * t525;
t386 = t393 * pkin(1);
t490 = t523 ^ 2;
t811 = -0.2e1 * pkin(1);
t318 = -t762 + t391 * t687 + t408 * t490 * t811 + t355 * t526 + (-t416 - 0.8e1 * t671) * t450;
t788 = t325 / 0.2e1;
t628 = t788 + t318 / 0.2e1;
t291 = t322 + (-t386 * t588 + (t319 * t784 + t628 * t406) * t363) * t560;
t787 = -t326 / 0.2e1;
t292 = (-t386 * t587 + ((t319 / 0.2e1 + t787) * t406 - t628 * t404) * t363) * t560;
t826 = ((t291 * t524 + t292 * t527) * t816 - (-t291 * t527 + t292 * t524) * t809) * t813;
t720 = t564 / 0.3e1 + t566;
t390 = -0.4e1 / 0.9e1 * t617 + 0.4e1 / 0.9e1 * t559 - t554 / 0.9e1 + t720;
t502 = -t554 / 0.6e1;
t511 = 0.2e1 / 0.3e1 * t559;
t597 = t566 - t617;
t398 = t502 + t511 + t597;
t510 = 0.4e1 / 0.3e1 * t559;
t504 = -t554 / 0.3e1;
t643 = t504 + t484;
t439 = t510 + t643;
t474 = -t564 / 0.3e1 + t566;
t618 = -0.2e1 * t450;
t674 = 0.6e1 * t732;
t492 = t526 * t493;
t568 = pkin(1) * t564;
t746 = t492 * t568;
t693 = pkin(7) * t746;
t705 = 0.4e1 * t793;
t735 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t340 = 0.4e1 * t693 + t390 * t674 + t439 * t735 + (t398 * t705 + t474 * t618) * pkin(1);
t542 = 0.6e1 * t564;
t506 = -0.2e1 / 0.3e1 * t554;
t550 = 0.2e1 * t566;
t642 = t506 + t511 + t550;
t572 = t559 ^ 2;
t641 = t506 + t484;
t725 = (t511 + t641) * t484 + t572;
t352 = t439 * t608 + (t542 + t642) * t559 + t725;
t475 = 0.10e2 / 0.3e1 * t564;
t366 = (t475 + t642) * t559 + t725;
t514 = -t559 / 0.3e1;
t472 = t514 + t566;
t410 = t472 * t437;
t736 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t368 = t439 * t736 + t410;
t483 = -0.3e1 * t559 + t566;
t632 = 0.8e1 * t693;
t423 = t483 * t632;
t448 = -t554 + t639;
t455 = t479 + pkin(7);
t476 = t484 ^ 2;
t528 = 0.15e2 * t562;
t529 = 0.15e2 * t564;
t565 = t566 ^ 2;
t546 = 0.3e1 * t565;
t571 = pkin(3) * t559;
t556 = t571 ^ 2;
t395 = t437 + t439;
t460 = t712 * t564;
t699 = 0.4e1 * t479;
t339 = pkin(7) * t395 * t699 + t460 * t800 + t352;
t459 = pkin(7) * t700;
t486 = -0.3e1 * t564 + t566;
t675 = 0.4e1 * t732;
t409 = t459 + t675 + t486;
t488 = t522 * t489;
t751 = t488 * t571;
t685 = 0.8e1 * t751;
t696 = 0.6e1 * t477;
t593 = t339 * t696 + t409 * t685;
t549 = 0.3e1 * t566;
t713 = -t554 + t559;
t640 = t549 + t713;
t659 = 0.12e2 * t732;
t661 = 0.12e2 * t749;
t321 = t340 * t661 + t423 + t368 * t659 + t556 + (t529 + t640) * t572 + (t640 * t542 + t713 * t550 + t528 + t546) * t559 + t476 * t448 + t593 * t455 + 0.6e1 * (t352 * t793 - t366 * t450) * pkin(1);
t482 = t559 + t694;
t665 = t564 * t478;
t419 = -t523 * t568 + t665;
t761 = t419 * t493;
t353 = -t710 * t478 + 0.2e1 * t761 + (t482 * t523 + t666 * t812) * pkin(1);
t716 = t543 + t566;
t465 = t559 + t716;
t422 = t465 * t478;
t466 = 0.3e1 * t559 + t484;
t755 = t466 * t523;
t383 = -pkin(1) * t755 + t422;
t673 = 0.8e1 * t564 * t571;
t797 = 0.4e1 * t562;
t453 = pkin(3) * t797 + t673;
t649 = t568 * t736;
t385 = t453 * t525 + 0.4e1 * t523 * t649;
t540 = 0.5e1 * t562;
t708 = t565 + t572;
t530 = 0.10e2 * t564;
t717 = t530 + t550;
t730 = t566 * t564;
t403 = t717 * t559 + t540 + t708 + 0.6e1 * t730;
t537 = 0.5e1 * t572;
t415 = t537 + (t530 + 0.6e1 * t566) * t559 + t476;
t653 = t455 * t751;
t615 = -0.8e1 * t653;
t684 = -0.4e1 * t749;
t757 = t455 * t522;
t779 = pkin(1) * t523;
t425 = t478 - t779;
t612 = t493 * t665;
t351 = -0.2e1 * t612 + t422 + (0.2e1 * t425 * t793 - t755) * pkin(1);
t805 = -0.4e1 * t351;
t328 = t353 * t684 + t385 * t493 + (t757 * t805 + (-t403 + t632) * t525) * pkin(3) + (-0.4e1 * t383 * t793 + (t415 + t615) * t523) * pkin(1);
t420 = t477 + t455;
t312 = t321 * t420 + t328 * t567;
t823 = t312 ^ 2;
t803 = 0.8e1 * t420;
t399 = 0.4e1 / 0.3e1 * t732 + t459 + t474;
t481 = -t554 - t559;
t463 = t549 + t481;
t756 = t463 * t564;
t432 = 0.10e2 * t756;
t541 = 0.7e1 * t564;
t487 = t489 ^ 2;
t752 = t487 * t572;
t753 = t481 * t566;
t754 = t476 * (-t559 + t464);
t815 = 0.7e1 * t556 + (t541 + t463) * t537 + (t432 + 0.21e2 * t562 + 0.9e1 * t565 + 0.6e1 * t753) * t559 + t754 - 0.24e2 * t399 * t752;
t503 = -t554 / 0.4e1;
t814 = t503 + t559 / 0.2e1;
t515 = -0.2e1 / 0.3e1 * t559;
t473 = t515 + t566;
t360 = t572 + (t506 + t515 + t717) * t559 + t540 + 0.2e1 * t756 + t566 * (t506 + t473);
t356 = t360 * t478;
t438 = 0.8e1 / 0.3e1 * t559 + t643;
t440 = t504 + t511 + t716;
t372 = -t438 * t779 + t440 * t478;
t518 = t564 / 0.2e1;
t446 = 0.5e1 / 0.6e1 * t559 + t518 + t502;
t381 = t446 * t695 + t736 * t779;
t624 = t568 * t470;
t374 = t537 + (0.5e1 * t564 + t463) * t538 + (t515 + t641) * t484;
t763 = t374 * t523;
t802 = -0.8e1 * t492;
t331 = t624 * t802 + t381 * t677 + t356 + (t372 * t705 - t763) * pkin(1);
t806 = -0.6e1 * t331;
t553 = t554 ^ 2;
t711 = t562 + t565;
t715 = t550 - t554;
t731 = t566 * t554;
t594 = t715 * t564 + t553 / 0.6e1 + t711 - t731;
t589 = 0.5e1 / 0.6e1 * t572 + t594;
t373 = (t475 + t715) * t559 + t589;
t804 = 0.6e1 * t373;
t798 = -0.8e1 * t526;
t796 = pkin(1) * pkin(7);
t792 = pkin(7) * t562;
t310 = 0.1e1 / t312;
t791 = t310 / 0.2e1;
t790 = t310 / 0.4e1;
t311 = 0.1e1 / t823;
t789 = -t311 / 0.4e1;
t786 = t363 / 0.2e1;
t625 = 0.20e2 / 0.3e1 * t564;
t548 = 0.4e1 * t566;
t646 = 0.2e1 / 0.3e1 * t554 + t511 + t548;
t647 = 0.4e1 / 0.3e1 * t554 + t510 - 0.2e1 * t566;
t785 = t572 / 0.2e1 - (-t625 + t646) * t559 / 0.2e1 + 0.3e1 / 0.2e1 * t562 - t647 * t564 / 0.2e1 - t565 / 0.2e1;
t782 = -t567 / 0.4e1;
t781 = t567 / 0.4e1;
t780 = pkin(1) * t474;
t778 = pkin(3) * t567;
t536 = -0.5e1 * t554;
t539 = 0.7e1 * t562;
t722 = t553 / 0.2e1 - t572 / 0.2e1;
t607 = -0.3e1 * t731 + t546 + t722;
t507 = -0.3e1 / 0.2e1 * t554;
t721 = t507 + t549;
t726 = t484 * ((t507 + t550) * t564 - 0.3e1 / 0.2e1 * t731 + t711 + t722) + t556;
t348 = (t541 + t721) * t572 + (t539 + (t536 + 0.10e2 * t566) * t564 + t607) * t559 + t726;
t777 = t348 * pkin(1);
t775 = -0.2e1 * t796;
t429 = t566 + t559 / 0.4e1 + t564 / 0.4e1 - t554 / 0.8e1;
t718 = 0.4e1 / 0.7e1 * t566 - t554 / 0.7e1;
t345 = -0.32e2 / 0.21e2 * t429 * t617 + 0.5e1 / 0.42e2 * t572 + (0.16e2 / 0.21e2 * t564 + t718) * t559 + t562 / 0.7e1 + t718 * t564 + t565 - 0.3e1 / 0.7e1 * t731 + t553 / 0.42e2;
t431 = t720 + t814;
t516 = 0.4e1 / 0.3e1 * t564;
t349 = -0.8e1 / 0.3e1 * t431 * t617 + 0.5e1 / 0.18e2 * t572 + (t516 + t504) * t559 + t565 - t562 / 0.3e1 + t553 / 0.18e2 + (t510 + 0.2e1 / 0.3e1 * t564 + t506) * t566;
t396 = -t572 / 0.6e1 + t594;
t719 = t518 + t566;
t397 = -0.2e1 / 0.3e1 * t617 + t503 + t719;
t468 = (t548 + t554) * t564;
t505 = -t554 / 0.2e1;
t445 = t505 + t639;
t619 = -0.4e1 * t450;
t602 = t445 * t619;
t622 = 0.16e2 * t693;
t491 = t493 ^ 2;
t733 = t562 * t491;
t678 = 0.8e1 * t733;
t327 = t473 * t678 + t397 * t622 + 0.14e2 * t345 * t732 - t710 * t572 + (t468 - 0.10e2 / 0.3e1 * t562 + 0.2e1 * t565 - t731) * t559 + t396 * t735 + (0.6e1 * t349 * t793 + t474 * t602) * pkin(1);
t337 = -0.6e1 * t373 * t617 + (t528 + (-0.9e1 * t554 + 0.18e2 * t566) * t564 + t607) * t559 + (t529 + t721) * t572 + t726;
t598 = pkin(1) * t602;
t709 = t565 - t562;
t343 = t472 * t598 - t556 + (-t475 - t707) * t572 + (t468 + t572 / 0.6e1 - t553 / 0.6e1 + t709) * t559 + t396 * t566;
t576 = pkin(7) * t566;
t462 = -0.12e2 * pkin(7) * t568 + t576 * t821;
t471 = -0.8e1 * t562 + 0.12e2 * t730;
t359 = t462 * t526 + t471 * t493 + t622 + t678 + t711 - 0.6e1 * t730;
t375 = t437 * t736 + t445 * t483;
t434 = (-0.6e1 * t559 * t566 + t708) * t562;
t480 = -0.30e2 * t554 + 0.60e2 * t566;
t535 = -0.2e1 * t554;
t714 = t553 - t572;
t606 = -0.6e1 * t731 + 0.6e1 * t565 + t714;
t547 = 0.8e1 * t566;
t388 = t568 * t619 + t797 + (0.4e1 * t559 + t535 + t547) * t564;
t394 = -t564 + t597 + t814;
t341 = t632 + t388 * t493 + t445 * t486 + (t394 * t705 + t618 * t735) * pkin(1);
t655 = t341 * t751;
t354 = t598 + (t542 + t715) * t559 + t589;
t371 = t445 * t736 + t410;
t702 = pkin(7) * t479;
t656 = 0.6e1 * t702;
t329 = t354 * t656 + t371 * t659 + t337 + t423;
t669 = t329 * t477;
t308 = 0.16e2 * t434 * t491 + 0.32e2 * t375 * t693 + 0.24e2 * t343 * t732 + (t535 + t548 + 0.28e2 * t564) * t556 + t448 * t754 + (0.24e2 * t327 * t489 + t480 * t562 + t606 * t542 + t714 * t550 - 0.6e1 * t565 * t554 + 0.28e2 * t568 ^ 2 + 0.4e1 * t576 ^ 2) * t559 + (0.32e2 * t655 + 0.8e1 * t669) * t455 + 0.8e1 * (t337 * t793 - t348 * t450) * pkin(1) + (0.16e2 * t359 * t487 + t480 * t564 + 0.70e2 * t562 + t572 + t606) * t572;
t427 = 0.7e1 / 0.6e1 * t559 + t502 + t719;
t512 = t559 / 0.3e1;
t644 = t502 + t512 + t566;
t430 = t516 + t644;
t367 = -t427 * t779 + t430 * t478;
t441 = t564 + t644;
t467 = t538 - t710;
t377 = t441 * t478 - t467 * t779 / 0.2e1;
t645 = t554 / 0.3e1 + t512 + t550;
t605 = -0.8e1 / 0.3e1 * t733 - t461 - 0.5e1 / 0.3e1 * t562 + t645 * t564 + t566 * (t504 + t472);
t704 = -0.4e1 * t792;
t330 = t523 * t492 * t704 + t367 * t675 + t605 * t478 + (t377 * t705 + t523 * t785) * pkin(1);
t342 = -pkin(1) * t763 + t356;
t370 = -0.3e1 * t572 + (-t625 + t647) * t559 + t646 * t564 + t709;
t378 = -0.5e1 / 0.3e1 * t572 + (-t564 + t645) * t559 + t566 * (t514 + t643);
t701 = -0.2e1 * t779;
t344 = t370 * t478 + t378 * t701;
t424 = t478 + 0.2e1 * t779;
t444 = t505 + t465;
t347 = t486 * t478 + 0.4e1 * t761 + (t424 * t793 + t444 * t523) * t822;
t428 = t566 + 0.5e1 / 0.2e1 * t559 + 0.3e1 / 0.2e1 * t564 + t505;
t379 = t428 * t478 + t483 * t779 / 0.2e1;
t599 = 0.24e2 * t472 * t733 - t556 - (0.21e2 * t564 + t463) * t572 - (t432 + t546 + 0.35e2 * t562 + 0.2e1 * t753) * t559 - (t539 + (t536 + t547 - 0.5e1 * t559) * t564 + t566 * (-t559 + t707)) * t484;
t662 = -0.12e2 * t749;
t313 = t379 * t622 + t347 * t615 + t330 * t662 - 0.6e1 * t344 * t732 + (t525 * t599 + t757 * t806) * pkin(3) + (-0.6e1 * t342 * t793 + t523 * t815) * pkin(1);
t297 = t308 * t420 + t313 * t567;
t636 = t326 * t790;
t592 = t297 * t636 + t325 * t781;
t555 = 0.1e1 / pkin(4);
t734 = t555 / pkin(3) ^ 2;
t654 = t363 * t734;
t287 = t592 * t654;
t411 = t741 + t745;
t284 = t411 * t287;
t637 = t325 * t790;
t591 = t297 * t637 + t326 * t782;
t286 = t591 * t654;
t412 = t742 - t744;
t277 = -t286 * t412 - t284;
t275 = 0.1e1 / t277 ^ 2;
t285 = t412 * t287;
t276 = -t286 * t411 + t285;
t774 = t275 * t276;
t773 = 0.1e1 / t297 ^ 2 * t823;
t765 = t364 * t382;
t764 = t364 * t386;
t729 = t568 * t493;
t692 = pkin(7) * t729;
t620 = -0.12e2 * t692;
t760 = t420 * t523 * t620;
t750 = t488 * t572;
t748 = t489 * t571;
t747 = t492 * t562;
t743 = t522 * t559;
t738 = t526 * t564;
t458 = pkin(7) * t701;
t626 = 0.32e2 / 0.3e1 * t562;
t603 = t492 * t626;
t604 = 0.64e2 / 0.3e1 * t429 * t568;
t609 = t492 * t649;
t616 = t736 * t792;
t621 = -0.48e2 * t692;
t629 = pkin(7) * t675;
t676 = -0.2e1 * t732;
t630 = pkin(7) * t676;
t631 = pkin(7) * t677;
t648 = t523 * t738;
t650 = t472 * t746;
t660 = -0.32e2 * t747;
t663 = 0.64e2 * t751;
t664 = t474 * t479;
t668 = t420 * t478;
t670 = pkin(1) * t735;
t679 = 0.8e1 * t738;
t680 = -0.2e1 * t746;
t681 = -0.4e1 * t746;
t682 = 0.3e1 * t749;
t686 = -0.8e1 * t751;
t688 = pkin(7) * t732;
t698 = -0.6e1 * t477;
t280 = ((0.12e2 * t490 * t493 * t792 + t427 * t681 + t467 * t630 + t491 * t704) * t662 + 0.12e2 * t378 * t746 + (t785 * t662 + t815) * t479 + (t374 * t674 + t483 * t678) * pkin(7)) * t567 + t313 * t687 + (((t444 * t700 + t629 + t681) * t686 + (-t374 * t479 + t438 * t631 - 0.4e1 * t609) * t698) * t567 + t663 * t760) * t455 + (0.24e2 * (-pkin(7) * t491 * t626 - 0.16e2 * t431 * t688 - 0.4e1 * t445 * t664 - t492 * t604) * t749 - 0.64e2 * t491 * t616 - 0.96e2 * t445 * t650 - 0.48e2 * t373 * t688 + t777 * t798 + ((-t526 * t670 + t630 + t680) * t663 + 0.48e2 * (-t373 * t479 + t445 * t631 - 0.4e1 * t650) * t477) * t455) * t668 + ((0.2e1 * (-0.2e1 * t471 * t526 - t462 + t621 + t660) * t752 + 0.8e1 * (-t388 * t526 + t394 * t775) * t653 + (-0.28e2 * t345 * t738 - 0.6e1 * t349 * t796 + t397 * t621 + t473 * t660) * t682 + t455 * (-t354 * t796 - 0.4e1 * t371 * t738 - 0.4e1 * t483 * t692) * t696 + t434 * t802 + t375 * t620 - 0.6e1 * t343 * t738 - t337 * t796 + (-0.4e1 * t655 - t669) * pkin(1)) * t803 - t308 * pkin(1) + ((-0.8e1 * t367 * t738 + t603 * t478) * t662 - 0.96e2 * t472 * t747 * t478 + t379 * t621 + 0.12e2 * t344 * t738 + (t331 * t696 + t347 * t685 + (-0.24e2 * t458 + 0.64e2 * t648) * t752 + (0.48e2 * t377 * t749 + 0.6e1 * t342) * pkin(7)) * pkin(1) + ((t419 * t798 + t424 * t775) * t686 + (-0.4e1 * t372 * t796 + t381 * t679 + 0.24e2 * t493 * t624) * t698) * t455) * t567) * t523;
t623 = t739 * t820;
t697 = -0.4e1 * t477;
t298 = (t490 * t488 * t673 + (t482 * t479 + t680) * t684 + 0.4e1 * t609 + t466 * t629 + t415 * t479) * t567 + t328 * t687 + t661 * t760 + ((-0.8e1 / 0.3e1 * t746 + t631 - 0.2e1 * t664) * t661 - 0.24e2 * t650 - 0.24e2 * t439 * t688 - 0.6e1 * t366 * t479) * t668 + (0.24e2 * (-t420 * t483 - t525 * t778) * t692 + ((0.16e2 * t419 * t749 - 0.2e1 * t385) * t567 + (-0.144e3 * t390 * t749 - 0.24e2 * t368) * t420 * t564) * t526 + (t351 * t778 * t799 - t593 * t420 - t321 + ((0.8e1 * pkin(3) * t489 * t740 + 0.4e1 * t383) * t567 + (-0.48e2 * t398 * t749 - 0.6e1 * t352) * t420) * pkin(7)) * pkin(1)) * t523 + ((t686 * t479 + ((pkin(7) * t801 + t525 * t623) * t564 + (-0.2e1 * t425 * t794 - t466 * t526) * pkin(1)) * t697) * t567 + ((t458 - 0.8e1 * t648) * t685 + (-0.8e1 * pkin(7) * t612 - 0.8e1 * t460 * t739 + (-t395 * t794 - t439 * t666) * t821) * t696) * t420) * t455;
t586 = t592 * t819;
t634 = t326 * t789;
t728 = (-t386 * t586 + (t319 * t781 + t767 * t788 + t280 * t636 + (t298 * t634 + t318 * t790) * t297) * t363) * t734 + t286;
t417 = t564 * t623 * t795;
t426 = pkin(3) * t691 * t811;
t651 = t472 * t729;
t652 = t455 * t748;
t657 = -0.4e1 * t702;
t689 = pkin(7) * t738;
t281 = t313 * t633 + (0.32e2 * t417 * t420 - 0.8e1 * t426 * t567) * t653 + (0.24e2 * (-0.4e1 * t399 * t750 * t779 - t330 * t743 - t347 * t652) * t567 + (0.6e1 * t327 * t743 + 0.12e2 * t341 * t652 + 0.8e1 * t359 * t750) * t803) * t525 + ((t308 + (t329 * t803 + t567 * t806) * t455) * t525 + (((t430 * t677 + t441 * t657 - t605) * t662 - 0.16e2 * t428 * t693 + t370 * t674 + t360 * t656 + ((-t486 + t677) * t686 + (t440 * t657 + 0.8e1 * t446 * t732 - t360 + t632) * t698) * t455 - t599) * t567 + ((pkin(7) * t603 + 0.16e2 * t431 * t689 + 0.4e1 * t445 * t780 + t493 * t604) * t682 + 0.8e1 * t492 * t616 + 0.12e2 * t445 * t651 + t689 * t804 + t777 + (0.4e1 * (0.2e1 * t670 + 0.4e1 * t729) * t751 + (pkin(1) * t804 + 0.24e2 * t445 * t689 + 0.24e2 * t651) * t477) * t455) * t523 * t803) * t522) * pkin(3);
t299 = (t426 * t684 + (-0.8e1 * t353 * t740 - t453 * t493 + ((t710 + t676) * t684 + t403 + (t465 * t699 - 0.8e1 * t746) * pkin(7)) * pkin(3)) * t522 + ((t426 + (-t465 + 0.2e1 * t732) * t477) * t697 + (pkin(3) * t805 - 0.24e2 * t748 * t779) * t525) * t455) * t567 + t328 * t633 + 0.6e1 * (0.4e1 * t409 * t525 * t652 + (t417 + (0.8e1 / 0.3e1 * t729 + 0.2e1 * t780) * t667) * t683 + t340 * t614 + 0.4e1 * t651 * t667 + t439 * t417 + t366 * t442 + ((pkin(7) * t679 + t439 * t821) * t523 * t749 + t339 * t478) * t455) * t420 + t321 * t478;
t727 = (t382 * t586 + (t317 * t781 + t325 * t768 / 0.8e1 + t281 * t636 + (t299 * t634 + t316 * t790) * t297) * t363) * t734 - t286;
t638 = -t297 * t311 / 0.2e1;
t635 = t325 * t789;
t324 = 0.1e1 / t325 ^ 2;
t613 = pkin(3) / (t324 * t326 ^ 2 + 0.1e1) * t365 * t560;
t288 = 0.1e1 / (t335 * t773 + 0.1e1);
t610 = t288 / t297 * t312 * t334;
t600 = 0.1e1 / t325 * t613;
t596 = t324 * t326 * t613;
t590 = pkin(4) * t288 * t555 * t560 * t773 * t778;
t585 = t591 * t819;
t274 = 0.1e1 / t277;
t273 = 0.1e1 / (t275 * t276 ^ 2 + 0.1e1);
t269 = (t382 * t585 + (t281 * t637 + t316 * t782 - t326 * t768 / 0.8e1 + (t299 * t635 + t317 * t790) * t297) * t363) * t734;
t265 = (-t386 * t585 + (t280 * t637 + t318 * t782 + t767 * t787 + (t298 * t635 + t319 * t790) * t297) * t363) * t734;
t1 = [1, 0; 0.2e1 * (t318 * t786 + t326 * t764) * t600 - 0.2e1 * (t319 * t786 + t325 * t764) * t596, 0.2e1 * (t316 * t786 - t326 * t765) * t600 - 0.2e1 * (t317 * t786 - t325 * t765) * t596; -t826, -t827; 0.2e1 * t338 * t610 - 0.2e1 * (t280 * t791 + t298 * t638) * t590, t332 * t610 / 0.2e1 - 0.2e1 * (t281 * t791 + t299 * t638) * t590; 1, 0; -t826, -t827; 0, 1; -1, 0; t826, t827; ((t728 * t412 + (-t265 + t287) * t411) * t274 - (-t265 * t412 - t728 * t411 + t285) * t774) * t273, ((-t269 * t411 + t727 * t412 - t284) * t274 - ((-t269 - t287) * t412 - t727 * t411) * t774) * t273; 1, 0; t826, t827; 0, 0; 0, 0; 0, 0;];
W = t1;
