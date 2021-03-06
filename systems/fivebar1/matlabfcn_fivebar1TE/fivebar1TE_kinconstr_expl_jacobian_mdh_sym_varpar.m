% Jacobian of explicit kinematic constraints of
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% W [6x2]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = fivebar1TE_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:09:47
% EndTime: 2020-04-27 09:09:58
% DurationCPUTime: 10.99s
% Computational Cost: add. (77256->632), mult. (232770->1011), div. (516->16), fcn. (28005->6), ass. (0->440)
t526 = pkin(2) ^ 2;
t522 = pkin(3) ^ 2;
t528 = pkin(1) ^ 2;
t676 = t522 - t528;
t423 = t676 * t526;
t479 = t522 / 0.3e1;
t514 = pkin(5) ^ 2;
t517 = pkin(4) ^ 2;
t484 = sin(qJ(2));
t485 = sin(qJ(1));
t710 = t484 * t485;
t412 = pkin(2) * t710;
t581 = pkin(3) * t412;
t355 = -0.4e1 / 0.9e1 * t581 + t528 + t526 / 0.3e1 + t479 + t517 / 0.9e1 - t514 / 0.9e1;
t463 = -t514 / 0.6e1;
t482 = t526 / 0.2e1;
t684 = t482 + t528;
t552 = -t581 + t684;
t363 = t517 / 0.6e1 + t463 + t552;
t465 = -t514 / 0.3e1;
t473 = t517 / 0.3e1;
t674 = t526 + t528;
t600 = t522 + t674;
t388 = t465 + t473 + t600;
t407 = -0.2e1 * t581;
t434 = -t522 / 0.3e1 + t528;
t486 = cos(qJ(2));
t440 = pkin(3) * t486;
t671 = 4 * pkin(1);
t451 = t486 ^ 2;
t699 = t522 * t451;
t703 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t450 = t486 * t451;
t537 = pkin(3) * t522;
t715 = t450 * t537;
t310 = t434 * t407 + 0.6e1 * t355 * t699 + t388 * t703 + (t363 * t440 + t715) * t671;
t500 = 0.6e1 * t522;
t566 = -0.4e1 * t581;
t467 = -0.2e1 / 0.3e1 * t514;
t472 = 0.2e1 / 0.3e1 * t517;
t511 = 0.2e1 * t528;
t603 = t467 + t472 + t511;
t444 = t522 + t528;
t534 = t526 ^ 2;
t602 = t467 + t444;
t690 = (t472 + t602) * t444 + t534;
t322 = t388 * t566 + (t500 + t603) * t526 + t690;
t361 = t407 + t388;
t675 = -t526 + t528;
t424 = t675 * t522;
t657 = pkin(1) * t440;
t611 = 0.4e1 * t657;
t754 = 0.4e1 * t451;
t311 = t361 * t611 + t424 * t754 + t322;
t435 = t528 - t526 / 0.3e1;
t374 = t435 * t407;
t702 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t336 = t388 * t702 + t374;
t438 = 0.10e2 / 0.3e1 * t522;
t338 = (t438 + t603) * t526 + t690;
t420 = 0.2e1 * t657;
t446 = -0.3e1 * t522 + t528;
t644 = 0.4e1 * t699;
t762 = t446 + t644;
t372 = t420 + t762;
t505 = -0.3e1 * t526;
t447 = t505 + t528;
t639 = pkin(1) * t715;
t589 = 0.8e1 * t639;
t392 = t447 * t589;
t599 = -t514 + t444;
t410 = t517 + t599;
t437 = t444 ^ 2;
t487 = cos(qJ(1));
t454 = t487 ^ 2;
t520 = t522 ^ 2;
t492 = 0.15e2 * t520;
t493 = 0.15e2 * t522;
t527 = t528 ^ 2;
t507 = 0.3e1 * t527;
t533 = pkin(2) * t526;
t523 = t533 ^ 2;
t418 = pkin(1) + t440;
t453 = t487 * t454;
t713 = t453 * t533;
t617 = t418 * t713;
t575 = -0.8e1 * t617;
t510 = 0.3e1 * t528;
t679 = -t514 + t517;
t601 = t510 + t679;
t610 = 0.6e1 * t657;
t622 = 0.12e2 * t699;
t627 = pkin(3) * t710;
t721 = t418 * t487;
t294 = t372 * t575 + t392 + t336 * t622 + t322 * t610 + t523 + (t493 + t601) * t534 + t437 * t410 + (0.12e2 * t310 * t454 + t601 * t500 + t679 * t511 + t492 + t507) * t526 + 0.6e1 * (-t311 * t721 - t338 * t627) * pkin(2);
t501 = 0.3e1 * t522;
t429 = t501 + t674;
t441 = pkin(2) * t485;
t393 = t429 * t441;
t733 = pkin(3) * t484;
t398 = t441 - t733;
t633 = t522 * t441;
t570 = t451 * t633;
t504 = 0.3e1 * t526;
t427 = t504 + t444;
t719 = t427 * t484;
t739 = pkin(1) * t486;
t323 = -0.2e1 * t570 + t393 + (0.2e1 * t398 * t739 - t719) * pkin(3);
t765 = 0.2e1 * t522;
t443 = t765 + t526;
t708 = t484 * t537;
t390 = t633 - t708;
t722 = t390 * t451;
t324 = t443 * t733 + 0.2e1 * t722 + (-t676 + t420) * t441;
t354 = -pkin(3) * t719 + t393;
t752 = 0.4e1 * t520;
t425 = pkin(2) * t752 + 0.8e1 * t522 * t533;
t613 = t537 * t702;
t356 = t425 * t485 + 0.4e1 * t484 * t613;
t498 = 0.5e1 * t520;
t672 = t527 + t534;
t494 = 0.10e2 * t522;
t682 = t494 + t511;
t696 = t528 * t522;
t371 = t682 * t526 + t498 + t672 + 0.6e1 * t696;
t503 = 0.5e1 * t534;
t379 = t503 + (t494 + 0.6e1 * t528) * t526 + t437;
t574 = 0.8e1 * t617;
t698 = t526 * t454;
t643 = -0.4e1 * t698;
t300 = t324 * t643 + t356 * t451 + (-0.4e1 * t354 * t739 + (t379 + t574) * t484) * pkin(3) + (0.4e1 * t323 * t721 + (-t371 + t589) * t485) * pkin(2);
t737 = pkin(2) * t487;
t417 = pkin(1) - t737;
t389 = t417 + t440;
t598 = -t514 + t674;
t561 = t522 + t598;
t404 = -t517 + t561;
t658 = pkin(1) * t737;
t422 = -0.2e1 * t658;
t370 = t422 + t404;
t351 = t407 + t370;
t401 = t417 * t486;
t678 = t514 - t528;
t642 = 0.2e1 * t698;
t373 = t422 + t642 + t675;
t757 = -0.2e1 * t373;
t557 = t451 * t757 + t678;
t748 = -pkin(4) + pkin(5);
t749 = -pkin(4) - pkin(5);
t767 = 0.4e1 * pkin(3);
t307 = 0.4e1 * t423 * t454 + 0.4e1 * t404 * t658 - t520 - (t528 + (pkin(2) - t749) * (pkin(2) + t749)) * (t528 + (pkin(2) - t748) * (pkin(2) + t748)) + (t505 + t517 + t557) * t765 + (-t351 * t401 + t370 * t412) * t767;
t529 = sqrt(t307);
t287 = t294 * t389 + t300 * t529;
t770 = t287 ^ 2;
t769 = 2 * pkin(1);
t768 = 0.2e1 * pkin(3);
t629 = pkin(3) * t401;
t337 = t407 + t422 + t600 + 0.2e1 * t629;
t335 = 0.1e1 / t337 ^ 2;
t766 = -0.2e1 * t335;
t756 = 0.8e1 * t389;
t764 = 0.8e1 * t526;
t368 = 0.4e1 / 0.3e1 * t699 + t420 + t434;
t442 = -t514 - t517;
t426 = t510 + t442;
t720 = t426 * t522;
t400 = 0.10e2 * t720;
t499 = 0.7e1 * t522;
t452 = t454 ^ 2;
t714 = t452 * t534;
t717 = t442 * t528;
t718 = t437 * (-t517 + t599);
t763 = 0.7e1 * t523 + (t499 + t426) * t503 + (t400 + 0.21e2 * t520 + 0.9e1 * t527 + 0.6e1 * t717) * t526 + t718 - 0.24e2 * t368 * t714;
t761 = t465 - t517 / 0.3e1;
t513 = t514 ^ 2;
t516 = t517 ^ 2;
t760 = -t513 / 0.6e1 + t516 / 0.6e1;
t759 = -0.4e1 * pkin(2);
t677 = t520 + t527;
t681 = t511 - t514;
t697 = t528 * t514;
t365 = t681 * t522 + t677 - t697 - t760;
t548 = t365 + t534;
t339 = (t438 + t681) * t526 + t548;
t758 = -0.6e1 * t339;
t755 = -0.8e1 * t450;
t753 = -0.8e1 * t486;
t751 = 0.2e1 * t526;
t750 = pkin(1) * pkin(3);
t285 = 0.1e1 / t287;
t747 = t285 / 0.2e1;
t746 = t285 / 0.4e1;
t286 = 0.1e1 / t770;
t745 = -t286 / 0.4e1;
t334 = 0.1e1 / t337;
t744 = t334 / 0.2e1;
t743 = 0.4e1 / 0.3e1 * t526;
t742 = -t529 / 0.4e1;
t741 = t529 / 0.4e1;
t740 = pkin(1) * (t412 - pkin(3));
t738 = pkin(1) * t520;
t736 = pkin(3) * t417;
t735 = pkin(3) * t434;
t734 = pkin(3) * t451;
t496 = -0.5e1 * t514;
t497 = 0.7e1 * t520;
t688 = t513 / 0.2e1 - t516 / 0.2e1;
t565 = -0.3e1 * t697 + t507 + t688;
t468 = -0.3e1 / 0.2e1 * t514;
t685 = t468 + t510;
t691 = t444 * ((t468 + t511) * t522 - 0.3e1 / 0.2e1 * t697 + t677 + t688) + t523;
t319 = (t499 + t685) * t534 + (t497 + (t496 + 0.10e2 * t528) * t522 + t565) * t526 + t691;
t731 = t319 * pkin(3);
t396 = t528 + t526 / 0.4e1 + t522 / 0.4e1 - t514 / 0.8e1;
t683 = 0.4e1 / 0.7e1 * t528 - t514 / 0.7e1;
t316 = -0.32e2 / 0.21e2 * t396 * t581 + t534 / 0.7e1 + (0.16e2 / 0.21e2 * t522 + t683) * t526 + t520 / 0.7e1 + t683 * t522 + t527 - 0.3e1 / 0.7e1 * t697 + t513 / 0.42e2 - t516 / 0.42e2;
t464 = -t514 / 0.4e1;
t397 = t464 + t479 + t684;
t478 = 0.4e1 / 0.3e1 * t522;
t320 = -0.8e1 / 0.3e1 * t397 * t581 + t534 / 0.3e1 + (t478 + t465) * t526 + t527 - t520 / 0.3e1 + (t743 + 0.2e1 / 0.3e1 * t522 + t467) * t528 + t513 / 0.18e2 - t516 / 0.18e2;
t480 = t522 / 0.2e1;
t367 = -0.2e1 / 0.3e1 * t581 + t528 + t480 + t464;
t509 = 0.4e1 * t528;
t430 = (t509 + t514) * t522;
t436 = t528 - 0.2e1 / 0.3e1 * t526;
t466 = -t514 / 0.2e1;
t409 = t466 + t600;
t554 = t409 * t566;
t624 = 0.16e2 * t715;
t449 = t451 ^ 2;
t700 = t520 * t449;
t648 = 0.8e1 * t700;
t299 = t436 * t648 + 0.14e2 * t316 * t699 + t434 * t554 - t676 * t534 + (t430 - 0.10e2 / 0.3e1 * t520 + 0.2e1 * t527 - t697) * t526 + t365 * t703 + (0.6e1 * t320 * t440 + t367 * t624) * pkin(1);
t309 = t581 * t758 + (t492 + (-0.9e1 * t514 + 0.18e2 * t528) * t522 + t565) * t526 + (t493 + t685) * t534 + t691;
t325 = t554 + (t500 + t681) * t526 + t548;
t342 = t409 * t702 + t374;
t301 = t325 * t610 + t342 * t622 + t309 + t392;
t495 = -0.2e1 * t514;
t508 = 0.8e1 * t528;
t585 = t708 * t759;
t360 = t485 * t585 + t752 + (0.4e1 * t526 + t495 + t508) * t522;
t366 = t464 - t522 + t552;
t653 = 0.8e1 * t715;
t659 = 0.4e1 * t440;
t312 = t407 * t703 + t360 * t451 + t409 * t446 + (t366 * t659 + t653) * pkin(1);
t673 = t527 - t520;
t314 = t435 * t554 - t523 + (-t438 + t678) * t534 + (t430 + t673 + t760) * t526 + t365 * t528;
t530 = pkin(1) * t528;
t416 = -0.12e2 * pkin(1) * t537 + t530 * t767;
t433 = -0.8e1 * t520 + 0.12e2 * t696;
t578 = pkin(1) * t624;
t330 = t416 * t486 + t433 * t451 + t578 + t648 + t677 - 0.6e1 * t696;
t345 = t407 * t702 + t409 * t447;
t402 = (-0.6e1 * t526 * t528 + t672) * t520;
t439 = -0.30e2 * t514 + 0.60e2 * t528;
t680 = t513 - t516;
t564 = -0.6e1 * t697 + 0.6e1 * t527 + t680;
t283 = -0.32e2 * t312 * t617 + 0.16e2 * t402 * t449 + 0.24e2 * t314 * t699 + (t495 + t509 + 0.28e2 * t522) * t523 + t410 * t718 + (0.24e2 * t299 * t454 + t439 * t520 + t564 * t500 + t680 * t511 - 0.6e1 * t527 * t514 + 0.4e1 * t530 ^ 2 + 0.28e2 * t537 ^ 2) * t526 + 0.8e1 * (-t301 * t721 - t319 * t627) * pkin(2) + (0.8e1 * t309 * t440 + 0.32e2 * t345 * t715) * pkin(1) + (0.16e2 * t330 * t452 + t439 * t522 + 0.70e2 * t520 + t534 + t564) * t534;
t687 = t463 - t517 / 0.6e1;
t606 = t528 + t687;
t381 = t743 + t480 + t606;
t563 = t482 + t606;
t382 = t478 + t563;
t331 = -t381 * t733 + t382 * t441;
t384 = t522 + t563;
t428 = t751 - t676;
t343 = t384 * t441 - t428 * t733 / 0.2e1;
t605 = t528 + t761;
t607 = t514 / 0.3e1 + t473 + t511;
t560 = -0.8e1 / 0.3e1 * t700 - t423 - 0.5e1 / 0.3e1 * t520 + t607 * t522 + t528 * t605;
t590 = 0.20e2 / 0.3e1 * t522;
t608 = 0.2e1 / 0.3e1 * t514 + t472 + t509;
t609 = 0.4e1 / 0.3e1 * t514 + 0.4e1 / 0.3e1 * t517 - 0.2e1 * t528;
t619 = -(-t534 + (-t590 + t608) * t526 - 0.3e1 * t520 + t609 * t522 + t527) * pkin(3) / 0.2e1;
t716 = t450 * t520;
t302 = t331 * t644 + t484 * t619 + t560 * t441 + (t343 * t440 - t484 * t716) * t671;
t477 = -0.2e1 / 0.3e1 * t517;
t686 = t467 + t477;
t332 = t534 + (t682 + t686) * t526 + t498 + 0.2e1 * t720 + t528 * (t528 + t686);
t344 = t503 + (0.5e1 * t522 + t426) * t751 + (t477 + t602) * t444;
t313 = t332 * t441 - t344 * t733;
t562 = t522 + t605;
t386 = t504 + t562;
t387 = t429 + t761;
t333 = -t386 * t733 + t387 * t441;
t395 = t480 + t526 + t687;
t668 = 0.2e1 * t441;
t348 = t395 * t668 + t702 * t733;
t632 = t537 * t441;
t646 = -0.4e1 * t699;
t303 = t348 * t646 + (t333 * t659 + t632 * t755) * pkin(1) + t313;
t341 = -0.3e1 * t534 + (-t590 + t609) * t526 + t608 * t522 + t673;
t349 = -0.5e1 / 0.3e1 * t534 + (-t522 + t607) * t526 + t528 * t562;
t662 = -0.2e1 * t733;
t315 = t341 * t441 + t349 * t662;
t604 = t466 - t517 / 0.2e1 + t528;
t383 = 0.3e1 / 0.2e1 * t526 + t501 + t604;
t399 = t441 + 0.2e1 * t733;
t318 = t446 * t441 + 0.4e1 * t722 + (t383 * t484 + t399 * t739) * t768;
t385 = t504 + 0.3e1 / 0.2e1 * t522 + t604;
t347 = t385 * t441 + t447 * t733 / 0.2e1;
t555 = 0.24e2 * t435 * t700 - t523 - (0.21e2 * t522 + t426) * t534 - (t400 + t507 + 0.35e2 * t520 + 0.2e1 * t717) * t526 - (t497 + (t496 + t508 - 0.5e1 * t517) * t522 + t528 * (-t517 - t678)) * t444;
t621 = -0.12e2 * t698;
t647 = -0.6e1 * t699;
t288 = t347 * t578 + t318 * t574 + t302 * t621 + t315 * t647 + (-0.6e1 * t313 * t739 + t484 * t763) * pkin(3) + (0.6e1 * t303 * t721 + t555 * t485) * pkin(2);
t280 = t283 * t389 + t288 * t529;
t369 = t422 + t517 + t561;
t630 = t373 * t440;
t326 = t369 * t417 + 0.2e1 * t630;
t329 = t369 * t486 + (t754 - 0.2e1) * t736;
t689 = -t412 + t401;
t362 = pkin(3) + t689;
t298 = t326 * t484 + t329 * t441 + t362 * t529;
t403 = t501 + t517 + t598;
t352 = t403 + t422 + t566;
t707 = t486 * t485;
t364 = pkin(2) * t707 + t417 * t484;
t724 = t364 * t529;
t297 = -t352 * t401 + t724 + (t403 * t710 - 0.2e1 * t487 * t740) * pkin(2) + (-t504 - t517 - t522 + t557 + t642) * pkin(3);
t596 = t297 * t746;
t550 = t280 * t596 + t298 * t742;
t515 = 0.1e1 / pkin(5);
t701 = t515 / pkin(4) ^ 2;
t618 = t334 * t701;
t275 = t550 * t618;
t706 = t486 * t487;
t377 = -t706 - t710;
t274 = t275 * t377;
t595 = t298 * t746;
t551 = t280 * t595 + t297 * t741;
t276 = t551 * t618;
t709 = t484 * t487;
t378 = -t707 + t709;
t270 = -t276 * t378 - t274;
t268 = 0.1e1 / t270 ^ 2;
t273 = t275 * t378;
t269 = -t276 * t377 + t273;
t730 = t268 * t269;
t729 = 0.1e1 / t280 ^ 2 * t770;
t432 = pkin(1) * t441;
t421 = 0.2e1 * t432;
t704 = t487 * t526;
t649 = -0.4e1 * t704;
t573 = t485 * t649;
t376 = t421 + t573;
t634 = pkin(2) * t709;
t580 = pkin(3) * t634;
t583 = -0.8e1 * t629;
t661 = -0.4e1 * t440;
t304 = t376 * t646 + (t432 - t580) * t583 + 0.4e1 * t370 * t580 + (pkin(2) * t351 * t661 - 0.8e1 * t423 * t487 + (t404 * t759 + t627 * t764) * pkin(1)) * t485;
t306 = 0.1e1 / t529;
t728 = t304 * t306;
t645 = 0.2e1 * t699;
t705 = t486 * t522;
t308 = (t351 * t736 + 0.2e1 * t373 * t705) * t484 + (t370 * t440 + t417 * t645) * t441;
t727 = t306 * t308;
t353 = -t378 * pkin(3) * pkin(2) + t432;
t726 = t335 * t353;
t357 = t364 * pkin(3);
t725 = t335 * t357;
t695 = t537 * t451;
t637 = pkin(1) * t695;
t576 = -0.12e2 * t637;
t723 = t389 * t484 * t576;
t712 = t453 * t534;
t711 = t454 * t533;
t694 = -0.2e1 * t750;
t419 = pkin(1) * t662;
t448 = t484 ^ 2;
t591 = 0.32e2 / 0.3e1 * t520;
t558 = t450 * t591;
t559 = 0.64e2 / 0.3e1 * t396 * t537;
t567 = t450 * t613;
t572 = t702 * t738;
t577 = -0.48e2 * t637;
t579 = -0.16e2 * pkin(1) * t397 * t522;
t584 = -0.4e1 * t409 * t735;
t638 = pkin(1) * t699;
t586 = 0.4e1 * t638;
t587 = -0.2e1 * t638;
t588 = -0.4e1 * t638;
t612 = t484 * t705;
t614 = t435 * t715;
t623 = -0.64e2 * t713;
t625 = -0.32e2 * t716;
t626 = pkin(3) * t703;
t628 = pkin(3) * t713;
t635 = t389 * t441;
t640 = pkin(1) * t700;
t641 = 0.3e1 * t698;
t651 = 0.8e1 * t713;
t652 = -0.8e1 * t713;
t654 = -0.2e1 * t715;
t655 = -0.4e1 * t715;
t656 = 0.2e1 * t727;
t660 = 0.2e1 * t440;
t665 = 0.6e1 * t737;
t667 = -0.6e1 * t737;
t271 = ((0.12e2 * t448 * t451 * t738 + t381 * t655 + t428 * t587 + t486 * t619 - 0.4e1 * t640) * t621 + 0.8e1 * t447 * t640 + 0.12e2 * t349 * t715 + 0.6e1 * t344 * t638 + t763 * t440) * t529 + t288 * t656 + (((t383 * t660 + t586 + t655) * t651 + (-t344 * t440 + t386 * t588 - 0.4e1 * t567) * t665) * t529 + t623 * t723) * t418 + (0.24e2 * (-pkin(1) * t449 * t591 - t450 * t559 + t451 * t579 + t486 * t584) * t698 - 0.64e2 * t449 * t572 - 0.96e2 * t409 * t614 - 0.48e2 * t339 * t638 + t731 * t753 + ((-t486 * t626 + t587 + t654) * t623 - 0.48e2 * (-t339 * t440 + t409 * t588 - 0.4e1 * t614) * t737) * t418) * t635 + ((0.2e1 * (-0.2e1 * t433 * t486 - t416 + t577 + t625) * t714 + (-t360 * t486 + t366 * t694) * t575 + 0.4e1 * t312 * t628 + (-0.28e2 * t316 * t705 - 0.6e1 * t320 * t750 + t367 * t577 + t436 * t625) * t641 + pkin(3) * t301 * t737 + t418 * (-t325 * t750 - 0.4e1 * t342 * t705 - 0.4e1 * t447 * t637) * t667 + t402 * t755 + t345 * t576 - 0.6e1 * t314 * t705 - t309 * t750) * t756 - t283 * pkin(3) + ((-0.8e1 * t331 * t705 + t558 * t441) * t621 - 0.96e2 * t435 * t716 * t441 + t347 * t577 + 0.12e2 * t315 * t705 + (t303 * t667 + t318 * t652 + (-0.24e2 * t419 + 0.64e2 * t612) * t714 + (0.48e2 * t343 * t698 + 0.6e1 * t313) * pkin(1)) * pkin(3) + ((t390 * t753 + t399 * t694) * t651 + (0.24e2 * pkin(1) * t451 * t632 - 0.4e1 * t333 * t750 + 0.8e1 * t348 * t705) * t665) * t418) * t529) * t484;
t620 = 0.12e2 * t698;
t650 = t485 * t764;
t663 = -0.2e1 * t735;
t664 = -0.4e1 * pkin(3) * t388;
t666 = 0.4e1 * t737;
t669 = t323 * t759;
t281 = (t522 * t448 * t652 + (t443 * t440 + t654) * t643 + 0.4e1 * t567 + t427 * t586 + t379 * t440) * t529 + t300 * t656 + t620 * t723 + ((-0.8e1 / 0.3e1 * t715 + t588 + t486 * t663) * t620 - 0.24e2 * t614 - 0.24e2 * t388 * t638 - 0.6e1 * t338 * t440) * t635 + (0.24e2 * (-t389 * t447 - t529 * t441) * t637 + ((0.16e2 * t390 * t698 - 0.2e1 * t356) * t529 + (-0.144e3 * t355 * t698 - 0.24e2 * t336) * t389 * t522) * t486 + (t487 * t529 * t669 + (t311 * t665 + t372 * t651) * t389 - t294 + ((pkin(2) * t454 * t650 + 0.4e1 * t354) * t529 + (-0.48e2 * t363 * t698 - 0.6e1 * t322) * t389) * pkin(1)) * pkin(3)) * t484 + ((0.8e1 * t486 * t628 + ((-pkin(3) * t427 + 0.4e1 * t522 * t412) * t486 + (-t398 * t733 - t699) * t769) * t666) * t529 + ((t419 - 0.8e1 * t612) * t652 + ((-0.8e1 * t424 * t484 + t664 * t441) * t486 + (-0.4e1 * t361 * t733 - 0.8e1 * t570) * pkin(1)) * t667) * t389) * t418;
t291 = -t724 + t362 * t656 + pkin(3) * t448 * t757 + t326 * t486 + (-t369 + t583) * t412;
t670 = -0.2e1 * pkin(1) * t526;
t292 = t689 * t529 + t364 * t656 + (t352 * t417 + 0.4e1 * t630) * t484 + (t670 * t706 + (t403 * t486 + 0.4e1 * t417 * t734) * pkin(2)) * t485;
t546 = t550 * t766;
t594 = t297 * t745;
t693 = (-t357 * t546 + (t271 * t596 + t291 * t742 - t298 * t727 / 0.2e1 + (t281 * t594 + t292 * t746) * t280) * t334) * t701 - t276;
t631 = pkin(2) * t706;
t582 = pkin(1) * t631;
t380 = -0.4e1 * t522 * t484 * t582;
t394 = t582 * t768;
t592 = t728 / 0.2e1;
t615 = t435 * t695;
t616 = t418 * t711;
t636 = pkin(1) * t705;
t272 = t288 * t592 + (-0.32e2 * t380 * t389 + 0.8e1 * t394 * t529) * t617 + (0.24e2 * (0.4e1 * t368 * t712 * t733 + t302 * t704 - t318 * t616) * t529 + (-0.6e1 * t299 * t704 + 0.12e2 * t312 * t616 - 0.8e1 * t330 * t712) * t756) * t485 + ((t283 + (t301 * t756 - 0.6e1 * t303 * t529) * t418) * t485 + (((t382 * t644 + t384 * t611 + t560) * t621 + t385 * t578 + t341 * t647 - 0.6e1 * t332 * t657 + (t762 * t651 + (t387 * t611 - 0.8e1 * t395 * t699 + t332 - 0.8e1 * t639) * t665) * t418 + t555) * t529 + ((-pkin(1) * t558 - t451 * t559 + t486 * t579 + t584) * t641 + t572 * t755 - 0.12e2 * t409 * t615 + t636 * t758 - t731 + (-0.4e1 * (-0.2e1 * t626 - 0.4e1 * t695) * t713 - (pkin(3) * t758 - 0.24e2 * t409 * t636 - 0.24e2 * t615) * t737) * t418) * t484 * t756) * t487) * pkin(2);
t282 = (t394 * t643 + (t324 * t650 + t425 * t451 + ((-t676 + t645) * t643 - t371 + (t429 * t661 + t653) * pkin(1)) * pkin(2)) * t487 + ((t394 + (t429 - 0.2e1 * t699) * t737) * t666 + (-0.24e2 * t711 * t733 + t669) * t485) * t418) * t529 + t300 * t592 + 0.6e1 * (0.4e1 * t372 * t485 * t616 + (t380 + (-0.8e1 / 0.3e1 * t695 + t663) * t634) * t642 + t310 * t573 + t435 * t487 * t451 * t585 + t388 * t380 - t338 * t580 + (-(-0.8e1 * t636 + t664) * t484 * t698 + t311 * t441) * t418) * t389 + t294 * t441;
t289 = t362 * t592 + t376 * t484 * t660 + ((-t484 * t529 + t329) * t487 + (t486 * t529 + (t417 * t769 + t369) * t484 + (-pkin(3) + 0.2e1 * t734 + t739) * t668) * t485) * pkin(2);
t290 = (t412 + t631) * t529 + t364 * t592 - 0.2e1 * t376 * t734 - (t421 - 0.4e1 * t580) * t401 + t454 * t484 * t670 + t403 * t634 + (pkin(3) * t649 + (-t352 * t486 + 0.2e1 * t740) * pkin(2)) * t485;
t692 = (t353 * t546 + (t272 * t596 + t289 * t742 - t298 * t728 / 0.8e1 + (t282 * t594 + t290 * t746) * t280) * t334) * t701 + t276;
t597 = -t280 * t286 / 0.2e1;
t593 = t298 * t745;
t296 = 0.1e1 / t297 ^ 2;
t571 = 0.1e1 / (t296 * t298 ^ 2 + 0.1e1) * t337;
t277 = 0.1e1 / (t307 * t729 + 0.1e1);
t568 = t277 / t280 * t287 * t306;
t556 = 0.1e1 / t297 * t571;
t553 = t296 * t298 * t571;
t549 = pkin(5) * t277 * t515 * t529 * t729;
t547 = t551 * t766;
t267 = 0.1e1 / t270;
t266 = 0.1e1 / (t268 * t269 ^ 2 + 0.1e1);
t265 = (t353 * t547 + (t290 * t741 + t297 * t728 / 0.8e1 + t272 * t595 + (t282 * t593 + t289 * t746) * t280) * t334) * t701;
t263 = (-t357 * t547 + (t292 * t741 + t297 * t727 / 0.2e1 + t271 * t595 + (t281 * t593 + t291 * t746) * t280) * t334) * t701;
t1 = [1, 0; ((t692 * t378 + (-t265 + t275) * t377) * t267 - (-t265 * t378 - t692 * t377 + t273) * t730) * t266, ((-t263 * t377 + t693 * t378 - t274) * t267 - ((-t263 - t275) * t378 - t693 * t377) * t730) * t266; 0, 1; 0.2e1 * (t289 * t744 - t298 * t726) * t556 - 0.2e1 * (t290 * t744 - t297 * t726) * t553, 0.2e1 * (t291 * t744 + t298 * t725) * t556 - 0.2e1 * (t292 * t744 + t297 * t725) * t553; t304 * t568 / 0.2e1 - 0.2e1 * (t272 * t747 + t282 * t597) * t549, 0.2e1 * t308 * t568 - 0.2e1 * (t271 * t747 + t281 * t597) * t549; 0, 0;];
W = t1;
