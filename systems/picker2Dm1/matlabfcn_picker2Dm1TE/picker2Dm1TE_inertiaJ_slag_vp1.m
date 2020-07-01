% Calculate joint inertia matrix for
% picker2Dm1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1TE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1TE_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm1TE_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:15:58
% EndTime: 2020-05-10 00:18:03
% DurationCPUTime: 109.36s
% Computational Cost: add. (2177856->869), mult. (5986846->1495), div. (79592->42), fcn. (1586054->14), ass. (0->658)
t788 = 4 * pkin(1);
t380 = sin(qJ(2));
t333 = pkin(3) * t380;
t621 = 0.8e1 * t333;
t345 = t380 ^ 2;
t344 = t380 * t345;
t426 = pkin(3) ^ 2;
t443 = pkin(3) * t426;
t683 = t344 * t443;
t583 = 0.32e2 * t683;
t763 = 0.2e1 * t380;
t634 = pkin(7) * t763;
t315 = pkin(3) * t634;
t432 = pkin(1) ^ 2;
t421 = pkin(4) ^ 2;
t435 = pkin(7) ^ 2;
t636 = t435 - t421;
t322 = t432 + t636;
t272 = t315 + t322;
t397 = 0.2e1 * t426;
t270 = t397 + t272;
t313 = t333 + pkin(7);
t640 = -t426 + t435;
t681 = t345 * t426;
t261 = t315 + t640 + 0.2e1 * t681;
t383 = cos(qJ(1));
t336 = pkin(1) * t383;
t595 = t261 * t336;
t203 = t270 * t313 + 0.2e1 * t595;
t349 = t383 ^ 2;
t205 = t270 * t383 + (0.4e1 * t349 - 0.2e1) * t313 * pkin(1);
t293 = t313 * t383;
t381 = sin(qJ(1));
t382 = cos(qJ(2));
t675 = t381 * t382;
t309 = pkin(3) * t675;
t653 = -t309 + t293;
t244 = pkin(1) + t653;
t334 = pkin(3) * t382;
t516 = pkin(1) * t309;
t296 = -0.2e1 * t516;
t238 = t296 + t272;
t785 = 4 * t432;
t350 = t426 * t785;
t665 = t426 * t435;
t320 = t350 - 0.4e1 * t665;
t412 = -0.2e1 * t435;
t413 = 0.2e1 * pkin(3);
t430 = t432 ^ 2;
t575 = -0.4e1 * pkin(3) * pkin(7) * t322;
t691 = t261 * t349;
t175 = t320 * t345 + t380 * t575 - t430 - (t435 - (t413 + pkin(4)) * pkin(4)) * (t435 + (t413 - pkin(4)) * pkin(4)) + (t412 + 0.2e1 * t421 - 0.4e1 * t426 - 0.4e1 * t691) * t432 + (-t238 * t293 + t272 * t309) * t788;
t436 = sqrt(t175);
t164 = t203 * t381 + t205 * t334 + t244 * t436;
t378 = cos(pkin(8));
t721 = sin(pkin(8));
t259 = -t378 * t383 - t721 * t381;
t729 = pkin(5) * t259;
t772 = -2 * pkin(1);
t249 = t729 * t772;
t418 = pkin(5) ^ 2;
t240 = t249 + 0.2e1 * t418;
t247 = -pkin(1) + t729;
t258 = t378 * t381 - t721 * t383;
t657 = t249 + t418;
t789 = 2 * pkin(1);
t226 = -(t789 + pkin(5)) * pkin(5) + t657;
t227 = pkin(5) * (t789 - pkin(5)) + t657;
t438 = sqrt(-t226 * t227);
t693 = t258 * t438;
t178 = pkin(5) * t693 - t240 * t247;
t703 = t178 * t164;
t403 = 3 * t432;
t306 = t397 + t403 + t636;
t497 = -0.4e1 * t516;
t234 = t306 + t315 + t497;
t294 = t309 - pkin(1);
t673 = t382 * t383;
t588 = pkin(3) * t673;
t689 = t313 * t381;
t246 = t588 + t689;
t694 = t246 * t436;
t163 = -t234 * t293 + t694 + (t294 * t634 + t306 * t675) * pkin(3) + (-0.2e1 * t691 + (0.2e1 * t345 - 0.4e1) * t426 - t322) * pkin(1);
t180 = -pkin(5) * t258 * t240 - t247 * t438;
t705 = t163 * t180;
t468 = t705 / 0.4e1 - t703 / 0.4e1;
t340 = t432 + t435;
t549 = t426 + t340;
t593 = pkin(1) * t293;
t211 = t296 + t315 + t549 + 0.2e1 * t593;
t209 = 0.1e1 / t211;
t427 = 0.1e1 / pkin(3);
t235 = t432 + t657;
t231 = 0.1e1 / t235;
t697 = t231 / pkin(5);
t569 = t427 * t697;
t500 = t209 * t569;
t137 = t468 * t500;
t384 = cos(pkin(9));
t704 = t164 * t180;
t706 = t163 * t178;
t469 = t706 / 0.4e1 + t704 / 0.4e1;
t456 = t469 * t500;
t734 = sin(pkin(9));
t127 = -t734 * t137 + t384 * t456;
t701 = t209 * t427;
t739 = -t383 / 0.2e1;
t741 = -t381 / 0.2e1;
t145 = (t163 * t741 + t164 * t739) * t701;
t414 = 2 * pkin(2);
t728 = pkin(6) * t127;
t122 = t728 * t414;
t416 = pkin(6) ^ 2;
t120 = t122 + 0.2e1 * t416;
t121 = -pkin(2) - t728;
t659 = t122 + t416;
t115 = -(t414 + pkin(6)) * pkin(6) + t659;
t116 = pkin(6) * (t414 - pkin(6)) + t659;
t711 = t115 * t116;
t437 = sqrt(-t711);
t129 = -t137 * t384 - t734 * t456;
t727 = pkin(6) * t129;
t62 = -t120 * t727 - t121 * t437;
t442 = pkin(2) ^ 2;
t119 = t442 + t659;
t417 = 0.1e1 / pkin(6);
t710 = 0.1e1 / t119 * t417;
t740 = t381 / 0.2e1;
t146 = (t163 * t739 + t164 * t740) * t701;
t749 = t146 / 0.2e1;
t61 = -t120 * t121 + t437 * t727;
t783 = t61 / 0.2e1;
t778 = (t145 * t783 + t62 * t749) * t710;
t793 = t127 * t778;
t290 = t435 + t426 / 0.4e1 + t432 / 0.4e1 - t421 / 0.8e1;
t420 = t421 ^ 2;
t434 = t435 ^ 2;
t444 = t426 ^ 2;
t649 = 0.4e1 / 0.7e1 * t435 - t421 / 0.7e1;
t662 = t435 * t421;
t196 = -0.32e2 / 0.21e2 * t290 * t516 + 0.5e1 / 0.42e2 * t444 + (0.16e2 / 0.21e2 * t432 + t649) * t426 + t430 / 0.7e1 + t649 * t432 + t434 - 0.3e1 / 0.7e1 * t662 + t420 / 0.42e2;
t651 = t432 / 0.3e1 + t435;
t361 = -t421 / 0.4e1;
t776 = t361 + t426 / 0.2e1;
t292 = t651 + t776;
t362 = -t421 / 0.3e1;
t364 = -0.2e1 / 0.3e1 * t421;
t368 = 0.4e1 / 0.3e1 * t426;
t374 = 0.4e1 / 0.3e1 * t432;
t197 = -0.8e1 / 0.3e1 * t292 * t516 + 0.5e1 / 0.18e2 * t444 + (t374 + t362) * t426 + t434 - t430 / 0.3e1 + t420 / 0.18e2 + (t368 + 0.2e1 / 0.3e1 * t432 + t364) * t435;
t639 = t430 + t434;
t411 = 0.2e1 * t435;
t643 = t411 - t421;
t467 = t643 * t432 + t420 / 0.6e1 + t639 - t662;
t252 = -t444 / 0.6e1 + t467;
t376 = t432 / 0.2e1;
t650 = t376 + t435;
t254 = -0.2e1 / 0.3e1 * t516 + t361 + t650;
t409 = 0.4e1 * t435;
t326 = (t409 + t421) * t432;
t373 = -0.2e1 / 0.3e1 * t426;
t330 = t373 + t435;
t331 = -t432 / 0.3e1 + t435;
t341 = -t432 + t435;
t363 = -t421 / 0.2e1;
t303 = t363 + t549;
t521 = -0.4e1 * t309;
t485 = t303 * t521;
t348 = t383 * t349;
t439 = pkin(1) * t432;
t679 = t348 * t439;
t610 = pkin(7) * t679;
t522 = 0.16e2 * t610;
t347 = t349 ^ 2;
t664 = t430 * t347;
t603 = 0.8e1 * t664;
t663 = t432 * t349;
t669 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t753 = pkin(7) * t383;
t168 = t330 * t603 + t254 * t522 + 0.14e2 * t196 * t663 + t341 * t444 + (t326 - 0.10e2 / 0.3e1 * t430 + 0.2e1 * t434 - t662) * t426 + t252 * t669 + (0.6e1 * t197 * t753 + t331 * t485) * pkin(1);
t365 = -0.3e1 / 0.2e1 * t421;
t385 = 15 * t430;
t392 = 0.18e2 * t435;
t406 = 0.3e1 * t434;
t652 = t420 / 0.2e1 - t444 / 0.2e1;
t495 = -0.3e1 * t662 + t406 + t652;
t410 = 0.3e1 * t435;
t648 = (15 * t432) + t410;
t423 = t443 ^ 2;
t658 = t340 * ((t365 + t411) * t432 - 0.3e1 / 0.2e1 * t662 + t639 + t652) + t423;
t351 = 0.10e2 / 0.3e1 * t432;
t460 = 0.5e1 / 0.6e1 * t444 + t467;
t221 = (t351 + t643) * t426 + t460;
t786 = -0.6e1 * t221;
t182 = t516 * t786 + (t385 + (t392 - 0.9e1 * t421) * t432 + t495) * t426 + (t365 + t648) * t444 + t658;
t372 = -t426 / 0.3e1;
t329 = t372 + t435;
t481 = pkin(1) * t485;
t638 = t434 - t430;
t188 = t329 * t481 - t423 + (-t351 - t636) * t444 + (t326 + t444 / 0.6e1 - t420 / 0.6e1 + t638) * t426 + t252 * t435;
t394 = -0.5e1 * t421;
t400 = 7 * t430;
t192 = (t365 + t410 + (7 * t432)) * t444 + (t400 + (t394 + 0.10e2 * t435) * t432 + t495) * t426 + t658;
t448 = pkin(7) * t435;
t321 = -0.12e2 * pkin(7) * t439 + t448 * t788;
t661 = t435 * t432;
t328 = -(8 * t430) + 0.12e2 * t661;
t206 = t321 * t383 + t328 * t349 + t522 + t603 + t639 - 0.6e1 * t661;
t339 = -0.3e1 * t426 + t435;
t670 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t217 = t296 * t670 + t303 * t339;
t637 = t434 + t444;
t283 = 0.16e2 * (t637 - 0.6e1 * t665) * t430;
t307 = -t421 + t549;
t314 = t336 + pkin(7);
t337 = -0.30e2 * t421 + 0.60e2 * t435;
t343 = t345 ^ 2;
t393 = -0.2e1 * t421;
t395 = -0.6e1 * t421;
t402 = 6 * t432;
t642 = t420 - t444;
t475 = 0.6e1 * t434 + t642 - 0.6e1 * t662;
t332 = t340 ^ 2;
t686 = t332 * (-t426 + t322);
t202 = t481 + (t402 + t643) * t426 + t460;
t263 = t329 * t296;
t215 = t303 * t670 + t263;
t529 = 0.8e1 * t610;
t282 = t339 * t529;
t628 = pkin(7) * t336;
t573 = 0.6e1 * t628;
t577 = 0.12e2 * t663;
t169 = t202 * t573 + t215 * t577 + t182 + t282;
t407 = 0.8e1 * t435;
t239 = t439 * t521 + t350 + (4 * t430) + (t393 + t407) * t432;
t480 = t435 - t516;
t250 = -t432 + t480 + t776;
t342 = -(3 * t432) + t435;
t520 = -0.2e1 * t309;
t632 = 0.4e1 * t753;
t186 = t529 + t239 * t349 + t303 * t342 + (t250 * t632 + t520 * t669) * pkin(1);
t774 = t169 * t621 + t186 * t583;
t135 = t283 * t347 + 0.32e2 * t217 * t610 + 0.24e2 * t188 * t663 + (t393 + t409 + (28 * t432)) * t423 + t307 * t686 + (0.24e2 * t168 * t345 + t337 * t430 + t434 * t395 + t475 * t402 + t642 * t411 + (28 * t439 ^ 2) + 0.4e1 * t448 ^ 2) * t426 + t774 * t314 + 0.8e1 * (t182 * t753 - t192 * t309) * pkin(1) + (0.16e2 * t206 * t343 + t337 * t432 + (70 * t430) + t444 + t475) * t444;
t360 = -t421 / 0.6e1;
t288 = 0.7e1 / 0.6e1 * t426 + t360 + t650;
t370 = t426 / 0.3e1;
t554 = t360 + t370 + t435;
t291 = t374 + t554;
t335 = pkin(1) * t381;
t212 = -t288 * t335 + t291 * t334;
t352 = -0.20e2 / 0.3e1 * t432;
t369 = 0.2e1 / 0.3e1 * t426;
t556 = 0.2e1 / 0.3e1 * t421 + t369 + t409;
t557 = 0.4e1 / 0.3e1 * t421 + t368 + t412;
t219 = -t444 + (t352 + t556) * t426 - (3 * t430) + t557 * t432 + t434;
t300 = t432 + t554;
t325 = t397 + t341;
t223 = t300 * t334 - t325 * t335 / 0.2e1;
t555 = t421 / 0.3e1 + t370 + t411;
t492 = -0.8e1 / 0.3e1 * t664 + t426 * t341 - 0.5e1 / 0.3e1 * t430 + t555 * t432 + t435 * (t362 + t329);
t600 = 0.4e1 * t663;
t752 = pkin(7) * t430;
t631 = -0.4e1 * t752;
t170 = t381 * t348 * t631 + t212 * t600 + t492 * t334 + (t219 * t741 + t223 * t632) * pkin(1);
t401 = 5 * t430;
t408 = 0.6e1 * t435;
t646 = t393 - 0.2e1 * t426;
t550 = t408 + t646;
t387 = 10 * t432;
t647 = t387 + t411;
t208 = t444 + (t364 + t373 + t647) * t426 + t401 + t550 * t432 + t435 * (t364 + t330);
t204 = t208 * t334;
t553 = t362 + t340;
t297 = 0.8e1 / 0.3e1 * t426 + t553;
t644 = t403 + t435;
t299 = t362 + t369 + t644;
t216 = -t297 * t335 + t299 * t334;
t304 = 0.5e1 / 0.6e1 * t426 + t376 + t360;
t620 = 0.2e1 * t334;
t229 = t304 * t620 + t670 * t335;
t586 = t439 * t334;
t504 = t348 * t586;
t602 = -0.4e1 * t663;
t396 = 0.5e1 * t444;
t551 = t364 + t340;
t222 = t396 + (t387 + t550) * t426 + (t373 + t551) * t340;
t698 = t222 * t381;
t171 = -0.8e1 * pkin(7) * t504 + t229 * t602 + t204 + (t216 * t632 - t698) * pkin(1);
t189 = -pkin(1) * t698 + t204;
t220 = -0.3e1 * t444 + (t352 + t557) * t426 + t556 * t432 + t638;
t224 = -0.5e1 / 0.3e1 * t444 + (-t432 + t555) * t426 + t435 * (t372 + t553);
t627 = -0.2e1 * t335;
t190 = t220 * t334 + t224 * t627;
t587 = t432 * t334;
t473 = -t381 * t439 + t587;
t274 = 0.4e1 * t473;
t284 = 0.2e1 * t335 + t334;
t323 = t426 + t644;
t302 = t363 + t323;
t193 = t342 * t334 + t274 * t349 + (t284 * t753 + t302 * t381) * t789;
t289 = t435 + 0.5e1 / 0.2e1 * t426 + 0.3e1 / 0.2e1 * t432 + t363;
t685 = t339 * t381;
t225 = t289 * t334 + pkin(1) * t685 / 0.2e1;
t775 = t410 - t421 - t426;
t305 = t775 * t387;
t645 = t394 - 0.5e1 * t426;
t484 = 0.24e2 * t329 * t664 - t423 - ((21 * t432) + t775) * t444 - (t646 * t435 + t305 + t406 + (35 * t430)) * t426 - (t400 + (t407 + t645) * t432 + t435 * (-t426 + t636)) * t340;
t568 = t314 * t683;
t508 = -0.8e1 * t568;
t582 = -0.12e2 * t681;
t688 = t314 * t380;
t626 = 0.2e1 * t336;
t318 = pkin(7) * t626;
t256 = 0.4e1 / 0.3e1 * t663 + t318 + t331;
t684 = t343 * t444;
t777 = 0.7e1 * t423 + ((35 * t432) + 0.15e2 * t435 + t645) * t444 + ((21 * t430) + t305 + 0.9e1 * t434 + (t395 - 0.6e1 * t426) * t435) * t426 + t686 - 0.24e2 * t256 * t684;
t142 = t225 * t522 + t193 * t508 + t170 * t582 - 0.6e1 * t190 * t663 + (-0.6e1 * t171 * t688 + t382 * t484) * pkin(3) + (-0.6e1 * t189 * t753 + t777 * t381) * pkin(1);
t275 = t333 + t314;
t112 = t135 * t275 + t142 * t436;
t243 = -0.4e1 / 0.9e1 * t516 + 0.4e1 / 0.9e1 * t426 - t421 / 0.9e1 + t651;
t255 = t360 + t369 + t480;
t298 = t368 + t553;
t599 = 0.6e1 * t663;
t185 = 0.4e1 * t610 + t243 * t599 + t298 * t669 + (t255 * t632 + t331 * t520) * pkin(1);
t552 = t364 + t369 + t411;
t655 = (t369 + t551) * t340 + t444;
t199 = t298 * t497 + (t402 + t552) * t426 + t655;
t213 = t298 * t670 + t263;
t214 = (t351 + t552) * t426 + t655;
t251 = t296 + t298;
t319 = t640 * t785;
t625 = 0.4e1 * t336;
t187 = pkin(7) * t251 * t625 + t319 * t349 + t199;
t262 = t318 + t600 + t342;
t608 = 0.8e1 * t683;
t622 = 0.6e1 * t333;
t466 = t187 * t622 + t262 * t608;
t581 = 0.12e2 * t681;
t787 = 6 * pkin(1);
t152 = t185 * t581 + t282 + t213 * t577 + t423 + (-t421 + t426 + t648) * t444 + (t385 + (t392 + t395 + 0.6e1 * t426) * t432 + t406 + (t393 + t397) * t435) * t426 + t332 * t307 + t466 * t314 + (t199 * t753 - t214 * t309) * t787;
t273 = 0.2e1 * t473;
t316 = pkin(7) * t620;
t404 = 2 * t432;
t338 = t404 + t426;
t201 = t341 * t334 + t273 * t349 + (t316 * t383 + t338 * t381) * pkin(1);
t281 = t323 * t334;
t324 = 0.3e1 * t426 + t340;
t687 = t324 * t381;
t233 = -pkin(1) * t687 + t281;
t784 = 8 * t432;
t598 = t443 * t784;
t730 = pkin(3) * t430;
t312 = t598 + 0.4e1 * t730;
t564 = t439 * t670;
t236 = t312 * t382 + 0.4e1 * t381 * t564;
t257 = t647 * t426 + t401 + t637 + 0.6e1 * t661;
t268 = t396 + (t387 + t408) * t426 + t332;
t607 = -0.4e1 * t681;
t285 = -t335 + t334;
t505 = t349 * t587;
t633 = 0.2e1 * t753;
t198 = -0.2e1 * t505 + t281 + (t285 * t633 - t687) * pkin(1);
t766 = -0.4e1 * t198;
t165 = t201 * t607 + t236 * t349 + (t688 * t766 + (-t257 + t529) * t382) * pkin(3) + (-0.4e1 * t233 * t753 + (t268 + t508) * t381) * pkin(1);
t141 = t152 * t275 + t165 * t436;
t139 = 0.1e1 / t141;
t751 = t139 / 0.2e1;
t548 = t112 * t751;
t422 = 0.1e1 / pkin(4);
t667 = t422 * t427;
t562 = t436 * t667;
t471 = -t145 * t562 / 0.2e1 + t146 * t548 * t667;
t678 = t380 * t381;
t264 = t673 + t678;
t677 = t380 * t383;
t265 = t675 - t677;
t746 = t163 / 0.4e1;
t544 = t139 * t746;
t737 = -t436 / 0.4e1;
t461 = t112 * t544 + t164 * t737;
t666 = t422 / pkin(3) ^ 2;
t570 = t209 * t666;
t96 = t461 * t570;
t745 = t164 / 0.4e1;
t543 = t139 * t745;
t736 = t436 / 0.4e1;
t463 = t112 * t543 + t163 * t736;
t97 = t463 * t570;
t67 = -t264 * t96 + t265 * t97;
t95 = t264 * t97;
t68 = t265 * t96 + t95;
t779 = (t145 * t548 + t436 * t749) * t667;
t792 = t471 * t68 + t67 * t779;
t791 = -t471 * t67 + t68 * t779;
t790 = t141 ^ 2;
t546 = t710 / 0.2e1;
t474 = t145 * t546 * t62 - t146 * t61 * t710 / 0.2e1;
t780 = t127 * t474;
t773 = 0.2e1 * pkin(7);
t770 = -2 * pkin(2);
t769 = 0.1e1 / t127;
t768 = 0.1e1 / t178;
t767 = 0.1e1 / t127 ^ 2;
t765 = 0.8e1 * t275;
t762 = -0.2e1 * t383;
t761 = -0.6e1 * t436;
t760 = pkin(1) * pkin(7);
t759 = pkin(4) * t779;
t758 = pkin(4) * t471;
t757 = pkin(6) * t778;
t756 = pkin(6) * t474;
t755 = pkin(7) * t349;
t754 = pkin(7) * t381;
t750 = t139 / 0.4e1;
t674 = t381 * t383;
t690 = t313 * t349;
t183 = (t261 * t674 + t690 * t334) * t784 + (t238 * t689 + t272 * t588) * t788;
t346 = t381 ^ 2;
t174 = 0.1e1 / t436;
t744 = t174 / 0.2e1;
t536 = t244 * t744;
t149 = -t694 + t183 * t536 + t261 * t346 * t772 + t203 * t383 + (-t270 - 0.8e1 * t593) * t309;
t748 = t149 / 0.4e1;
t747 = -t163 / 0.4e1;
t743 = t178 / 0.4e1;
t742 = t209 / 0.2e1;
t738 = t383 / 0.2e1;
t735 = 0.2e1 * t258 ^ 2;
t733 = pkin(1) * t331;
t732 = pkin(1) * t349;
t731 = pkin(2) / t119 ^ 2;
t726 = t192 * pkin(1);
t66 = 0.1e1 / t68 ^ 2;
t725 = t66 * t67;
t486 = pkin(7) * t505;
t287 = -0.4e1 * t486;
t317 = pkin(7) * t627;
t487 = t586 * t755;
t498 = t348 * t564;
t518 = pkin(1) * t588;
t530 = pkin(7) * t600;
t537 = t183 * t744;
t671 = t383 * t432;
t561 = t381 * t671;
t578 = -0.24e2 * t671;
t606 = 0.8e1 * t681;
t609 = -0.8e1 * t683;
t660 = t439 * t349;
t613 = pkin(7) * t660;
t623 = -0.4e1 * t333;
t523 = -0.24e2 * t613;
t656 = -0.24e2 * t329 * t504 + t523 * t685;
t113 = (t346 * t344 * t598 + (t338 * t336 - 0.2e1 * t679) * t607 + 0.4e1 * t498 + t324 * t530 + t268 * t336) * t436 + t165 * t537 + ((-0.8e1 / 0.3e1 * t504 + t287 - 0.2e1 * t331 * t518) * t581 - 0.24e2 * t298 * t486 - 0.6e1 * t214 * t518 + t656) * t275 + ((t273 * t383 * t606 + t236 * t762 - 0.24e2 * t487) * t436 + (0.12e2 * (-t243 * t671 - t613) * t581 + t213 * t578) * t275 + (0.4e1 * t198 * t436 * t333 - t466 * t275 - t152 + ((t606 * t334 + 0.4e1 * t233) * t436 + (-0.48e2 * t255 * t681 - 0.6e1 * t199) * t275) * pkin(7)) * pkin(1)) * t381 + ((t609 * t336 + ((0.4e1 * t309 * t383 - 0.2e1 * t755) * t432 + (-0.2e1 * t285 * t754 - t324 * t383) * pkin(1)) * t623) * t436 + ((t317 - 0.8e1 * t561) * t608 + (-0.8e1 * t486 - 0.2e1 * t319 * t674 + (-t251 * t754 - t298 * t588) * t788) * t622) * t275) * t314;
t535 = t246 * t744;
t676 = t380 * t426;
t151 = t653 * t436 + t183 * t535 + (t234 * t313 + 0.4e1 * t595) * t381 + (t633 * t676 + (t306 * t383 + t690 * t788) * pkin(3)) * t382;
t237 = t246 * t789;
t210 = 0.1e1 / t211 ^ 2;
t457 = t210 * t461;
t539 = -t164 * t174 / 0.8e1;
t140 = 0.1e1 / t790;
t542 = t140 * t747;
t528 = 0.32e2 / 0.3e1 * t430;
t488 = t348 * t528;
t489 = 0.64e2 / 0.3e1 * t290 * t439;
t502 = 0.32e2 * t568;
t510 = t670 * t752;
t524 = -0.48e2 * t613;
t579 = -0.32e2 * t348 * t430;
t580 = 0.24e2 * t681;
t584 = -0.96e2 * t329 * t348;
t592 = pkin(1) * t669;
t594 = t303 * t733;
t601 = -0.2e1 * t663;
t604 = 0.8e1 * t671;
t605 = -0.4e1 * t679;
t614 = pkin(7) * t663;
t624 = -0.6e1 * t333;
t715 = -0.4e1 * t760;
t716 = -0.6e1 * t760;
t77 = ((t302 * t626 + t530 + t605) * t508 + (0.12e2 * t346 * t349 * t752 + t288 * t605 + t347 * t631) * t582 + 0.12e2 * t224 * t679 + (-t219 * t582 / 0.2e1 + t777) * t336 + (t325 * t601 * t582 + t222 * t599 + t339 * t603) * pkin(7)) * t436 + t142 * t537 + t287 * t275 * t502 + (((pkin(7) * t297 * t602 - t222 * t336 - 0.4e1 * t498) * t761 + t656 * t765) * t688 + ((-pkin(7) * t347 * t528 - 0.16e2 * t292 * t614 - t348 * t489 - 0.4e1 * t383 * t594) * t580 - 0.64e2 * t347 * t510 + t439 * t303 * t584 - 0.48e2 * t221 * t614 - 0.8e1 * t383 * t726 + ((t592 * t762 + t605) * t583 + (-0.24e2 * t303 * t614 + t336 * t786) * t621) * t314) * t275 * t382) * pkin(3) + ((0.16e2 * (t328 * t762 - t321 + t524 + t579) * t684 + (t239 * t762 + t250 * t715 + t523) * t502 + (-0.28e2 * t196 * t671 + t197 * t716 + t254 * t524 + t330 * t579) * t580 + t314 * (t202 * t716 + t215 * t578) * t621 - 0.4e1 * t283 * t348 - 0.96e2 * t217 * t613 - 0.48e2 * t188 * t671 - 0.8e1 * t182 * t760) * t275 + ((-0.8e1 * t212 * t671 + t488 * t334) * t582 + t382 * t584 * t730 + t225 * t524 + 0.12e2 * t190 * t671 + (0.2e1 * (-t274 * t383 - t284 * t760) * t609 + (t216 * t715 + t229 * t604 + 0.24e2 * t487) * t624) * t314) * t436 + (-t774 * t275 - t135 + (t171 * t622 + t193 * t608 + (-0.24e2 * t317 + 0.64e2 * t561) * t684 + (0.48e2 * t223 * t681 + 0.6e1 * t189) * pkin(7)) * t436) * pkin(1)) * t381;
t724 = -(t237 * t457 + (t77 * t544 + t149 * t737 + t183 * t539 + (t113 * t542 + t151 * t750) * t112) * t209) * t666 + t97;
t458 = t210 * t463;
t540 = t163 * t174 / 0.8e1;
t541 = -t140 * t164 / 0.4e1;
t723 = (t237 * t458 + (t151 * t736 + t183 * t540 + t77 * t543 + (t113 * t541 + t139 * t748) * t112) * t209) * t666 + t96;
t566 = t329 * t660;
t590 = pkin(3) * t678;
t253 = 0.24e2 * t566 * t590;
t589 = pkin(3) * t677;
t526 = pkin(7) * t589;
t478 = t432 * t381 * t526;
t271 = 0.4e1 * t478;
t286 = t526 * t772;
t517 = pkin(1) * t590;
t672 = t382 * t426;
t565 = t380 * t672;
t267 = t316 + 0.4e1 * t565;
t612 = t426 * t754;
t654 = 0.2e1 * t517 + t316;
t172 = t267 * t602 + (t320 * t763 + t575) * t382 + (-t654 * t293 + 0.2e1 * t382 ^ 2 * t612 + (-t238 * t673 - t272 * t678) * pkin(3)) * t788;
t538 = t172 * t744;
t680 = t345 * t443;
t567 = t314 * t680;
t114 = (t286 * t607 + (-0.8e1 * t201 * t672 - t312 * t349 + ((-t341 + t601) * t607 + t257 + (t323 * t625 - 0.8e1 * t679) * pkin(7)) * pkin(3)) * t380 + ((t286 + (-t323 + 0.2e1 * t663) * t333) * t623 + (pkin(3) * t766 - 0.24e2 * t680 * t335) * t382) * t314) * t436 + t165 * t538 + (0.24e2 * t262 * t382 * t567 + (t271 + (0.8e1 / 0.3e1 * t660 + 0.2e1 * t733) * t590) * t581 + 0.24e2 * t185 * t565 + t253 + 0.24e2 * t298 * t478 + 0.6e1 * t214 * t517 + 0.6e1 * ((pkin(7) * t604 + t298 * t788) * t381 * t681 + t187 * t334) * t314) * t275 + t152 * t334;
t147 = t172 * t536 + t267 * t381 * t626 + ((t381 * t436 - t205) * t380 + (t383 * t436 + (t313 * t773 + t270) * t381 + (-pkin(1) + 0.2e1 * t732 + t753) * t620) * t382) * pkin(3);
t148 = (t309 - t589) * t436 + t172 * t535 - 0.2e1 * t267 * t732 - (t316 + 0.4e1 * t517) * t293 - 0.2e1 * t345 * t612 - t306 * t590 + (t676 * t788 + (-t234 * t383 + t294 * t773) * pkin(3)) * t382;
t230 = 0.2e1 * t518 + t654;
t574 = -0.4e1 * t628;
t611 = pkin(7) * t671;
t682 = t344 * t444;
t84 = t142 * t538 + (0.32e2 * t271 * t275 - 0.8e1 * t286 * t436) * t568 + (0.24e2 * (-0.4e1 * t256 * t682 * t335 - t170 * t676 - t193 * t567) * t436 + (0.48e2 * t168 * t676 + 0.96e2 * t186 * t567 + 0.64e2 * t206 * t682) * t275) * t382 + ((t135 + (t169 * t765 + t171 * t761) * t314) * t382 + (t314 * t253 * t765 + ((t291 * t602 + t300 * t574 - t492) * t582 - 0.16e2 * t289 * t610 + t220 * t599 + t208 * t573 + ((-t342 + t602) * t609 + (t299 * t574 + 0.8e1 * t304 * t663 - t208 + t529) * t624) * t314 - t484) * t436 + ((pkin(7) * t488 + 0.16e2 * t292 * t611 + t349 * t489 + 0.4e1 * t594) * t580 + 0.64e2 * t348 * t510 + 0.96e2 * t303 * t566 + 0.48e2 * t221 * t611 + 0.8e1 * t726 + ((0.2e1 * t592 + 0.4e1 * t660) * t583 + (t221 * t787 + 0.24e2 * t303 * t611) * t621) * t314) * t275 * t381) * t380) * pkin(3);
t722 = (-t230 * t458 + (t148 * t736 + t172 * t540 + t84 * t543 + (t114 * t541 + t147 * t750) * t112) * t209) * t666 - t96;
t162 = 0.1e1 / t163 ^ 2;
t506 = pkin(3) / (t162 * t164 ^ 2 + 0.1e1) * t211 * t427;
t470 = -0.2e1 * t162 * t164 * t506;
t477 = 0.2e1 / t163 * t506;
t700 = t210 * t230;
t533 = -t700 / 0.2e1;
t707 = t148 * t209;
t708 = t147 * t209;
t90 = (t708 / 0.2e1 + t164 * t533) * t477 + (t707 / 0.2e1 + t163 * t533) * t470;
t720 = t145 * t90;
t699 = t210 * t237;
t531 = t699 / 0.2e1;
t91 = 0.1e1 + (t149 * t742 + t164 * t531) * t477 + (t151 * t742 + t163 * t531) * t470;
t719 = t145 * t91;
t718 = t146 * t90;
t717 = t146 * t91;
t714 = pkin(3) * t717 - t336;
t713 = pkin(2) * t717 - t336;
t712 = 0.1e1 / t112 ^ 2 * t790;
t709 = t767 * t129;
t527 = t258 * t789;
t702 = 0.1e1 / t438 * (t226 + t227) * pkin(5) * t527;
t696 = t231 / pkin(1);
t232 = 0.1e1 / t235 ^ 2;
t695 = t232 * t258;
t692 = t259 * t438;
t668 = t417 / pkin(2);
t534 = t702 / 0.2e1;
t153 = (-t692 + (t247 * t789 - t240 + t534) * t258) * pkin(5);
t154 = -t247 * t702 / 0.2e1 + t418 * pkin(1) * t735 + (t240 * t259 - t693) * pkin(5);
t499 = t697 * t699;
t596 = pkin(1) * t695;
t78 = (t469 * t499 + ((t706 / 0.2e1 + t704 / 0.2e1) * t596 + (t151 * t743 + t153 * t746 + t154 * t745 + t180 * t748) * t697) * t209) * t427;
t79 = (-t468 * t499 + ((-t705 / 0.2e1 + t703 / 0.2e1) * t596 + (-t151 * t180 / 0.4e1 + t154 * t747 + t153 * t745 + t149 * t743) * t697) * t209) * t427;
t57 = t384 * t78 + t734 * t79;
t619 = t57 * t731;
t532 = -t700 / 0.4e1;
t464 = t707 / 0.4e1 + t163 * t532;
t465 = t708 / 0.4e1 + t164 * t532;
t102 = (t178 * t464 + t180 * t465) * t569;
t103 = (t178 * t465 - t180 * t464) * t569;
t73 = t384 * t102 + t734 * t103;
t618 = t73 * t731;
t617 = pkin(2) * t720;
t616 = pkin(3) * t720;
t60 = 0.1e1 / t61 ^ 2;
t615 = pkin(6) * t119 / (t60 * t62 ^ 2 + 0.1e1);
t585 = pkin(5) * t695;
t576 = 0.1e1 / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t442 * t767 * t711 / 0.4e1) * t668;
t572 = -t721 / 0.2e1;
t563 = t437 * t668;
t93 = 0.1e1 / t437;
t560 = -t121 * t93 / 0.2e1;
t559 = -t769 * t93 / 0.4e1;
t558 = t129 * t93 / 0.2e1;
t547 = -t112 * t140 / 0.2e1;
t545 = t767 * t437 / 0.2e1;
t525 = t129 * t416 * t770;
t519 = t121 * t770 + t120;
t507 = -t146 * t127 - t129 * t145;
t483 = -0.2e1 * t60 * t62 * t615;
t496 = t615 / t783;
t501 = pkin(6) * t546;
t476 = pkin(6) * (-t115 - t116) * t414;
t53 = t73 * t476;
t74 = t734 * t102 - t384 * t103;
t16 = ((t53 * t560 + t73 * t525 + (t120 * t74 + t437 * t73) * pkin(6)) * t546 - t62 * t618) * t496 + (-t61 * t618 + (-t437 * t74 + t519 * t73 + t53 * t558) * t501) * t483 + t90;
t99 = 0.1e1 / (t175 * t712 + 0.1e1);
t459 = -0.2e1 * pkin(3) * pkin(4) * t99 * t562 * t712;
t482 = 0.1e1 / t112 * t141 * t99 * t744;
t39 = t172 * t482 + (t114 * t547 + t84 * t751) * t459 + t90;
t494 = -t563 / 0.2e1;
t493 = t563 / 0.2e1;
t491 = -pkin(2) * t719 + t335;
t490 = -pkin(3) * t719 + t335;
t44 = t57 * t476;
t58 = -t384 * t79 + t734 * t78;
t14 = ((t44 * t560 + t57 * t525 + (t120 * t58 + t437 * t57) * pkin(6)) * t546 - t62 * t619) * t496 + (-t61 * t619 + (-t437 * t58 + t44 * t558 + t519 * t57) * t501) * t483 + t91;
t38 = t183 * t482 + (t113 * t547 + t77 * t751) * t459 + t91;
t479 = t127 * t145 - t129 * t146;
t280 = -rSges(2,1) * t383 + rSges(2,2) * t381;
t279 = rSges(8,1) * t382 - rSges(8,2) * t380;
t278 = rSges(2,1) * t381 + rSges(2,2) * t383;
t277 = rSges(8,1) * t380 + rSges(8,2) * t382;
t248 = -pkin(1) * t259 + pkin(5);
t241 = t249 + t404;
t181 = -pkin(1) * t258 * t241 + t248 * t438;
t179 = pkin(1) * t693 + t241 * t248;
t177 = 0.1e1 / t179 ^ 2;
t176 = 0.1e1 / t178 ^ 2;
t159 = (t178 * t738 + t180 * t741) * t697;
t158 = (t178 * t740 + t180 * t738) * t697;
t157 = (-t378 * t179 / 0.2e1 + t181 * t572) * t696;
t156 = (t179 * t572 + t378 * t181 / 0.2e1) * t696;
t144 = rSges(6,1) * t157 - rSges(6,2) * t156;
t143 = rSges(6,1) * t156 + rSges(6,2) * t157;
t134 = rSges(3,1) * t146 - rSges(3,2) * t145;
t133 = rSges(3,1) * t145 + rSges(3,2) * t146;
t111 = (-((t248 * t534 + t432 * pkin(5) * t735 + (t241 * t259 - t693) * pkin(1)) * t696 / 0.2e1 + t181 * t585) / t179 - (-t179 * t585 - (-t692 + (-0.2e1 * t248 * pkin(5) - t241 + t534) * t258) * pkin(1) * t696 / 0.2e1) * t181 * t177) / (t177 * t181 ^ 2 + 0.1e1) * t235 * t789;
t107 = 0.1e1 + (t154 * t768 * t697 + (-t153 * t176 * t697 + (-t178 * t176 + t768) * t232 * t527) * t180) * t235 / (t176 * t180 ^ 2 + 0.1e1) * pkin(5);
t101 = -t336 + t107 * (rSges(9,1) * t159 - rSges(9,2) * t158);
t100 = t335 - t107 * (rSges(9,1) * t158 + rSges(9,2) * t159);
t94 = 0.1e1 / (t129 ^ 2 * t767 + 0.1e1);
t81 = pkin(2) * t718;
t80 = pkin(3) * t718;
t76 = t134 * t91 - t336;
t75 = -t133 * t91 + t335;
t65 = 0.1e1 / t68;
t64 = rSges(5,1) * t471 - rSges(5,2) * t779;
t63 = rSges(5,1) * t779 + rSges(5,2) * t471;
t56 = rSges(7,1) * t507 + rSges(7,2) * t479;
t55 = -rSges(7,1) * t479 + rSges(7,2) * t507;
t52 = 0.1e1 / (t66 * t67 ^ 2 + 0.1e1);
t47 = (-t230 * t457 + (t84 * t544 + t147 * t737 + t172 * t539 + (t114 * t542 + t148 * t750) * t112) * t209) * t666;
t37 = (t73 * t709 + t74 * t769) * t94 + t90;
t35 = (t57 * t709 + t58 * t769) * t94 + t91;
t34 = t39 * t64 + t80;
t33 = -t39 * t63 - t616;
t32 = t38 * t64 + t714;
t31 = -t38 * t63 + t490;
t30 = t35 * t56 - t336;
t29 = -t35 * t55 + t335;
t28 = t792 * rSges(11,1) - t791 * rSges(11,2);
t27 = t791 * rSges(11,1) + t792 * rSges(11,2);
t26 = rSges(4,1) * t474 + rSges(4,2) * t778;
t25 = -rSges(4,1) * t778 + rSges(4,2) * t474;
t22 = (t494 * t778 + t780) * rSges(10,1) + (t474 * t493 + t793) * rSges(10,2);
t21 = (t474 * t494 - t793) * rSges(10,1) + (-t493 * t778 + t780) * rSges(10,2);
t18 = (-(-t264 * t47 + t722 * t265 - t95) * t65 - ((-t47 - t97) * t265 - t722 * t264) * t725) * t52 + t39;
t17 = (-(t724 * t264 + t723 * t265) * t65 - (-t723 * t264 + t724 * t265) * t725) * t52 + t38;
t15 = (t53 * t559 + t545 * t73) * t576 + t16;
t13 = t16 * t26 + t81;
t12 = -t16 * t25 - t617;
t11 = t18 * t28 + t39 * t758 + t80;
t10 = -t18 * t27 - t39 * t759 - t616;
t9 = t17 * t28 + t38 * t758 + t714;
t8 = -t17 * t27 - t38 * t759 + t490;
t7 = t14 * t26 + t713;
t6 = -t14 * t25 + t491;
t5 = (t44 * t559 + t545 * t57) * t576 + t14;
t4 = t15 * t22 + t16 * t756 + t81;
t3 = -t15 * t21 + t16 * t757 - t617;
t2 = t14 * t756 + t22 * t5 + t713;
t1 = t14 * t757 - t21 * t5 + t491;
t19 = [(t29 ^ 2 + t30 ^ 2) * m(7) + (t75 ^ 2 + t76 ^ 2) * m(3) + (t1 ^ 2 + t2 ^ 2) * m(10) + m(9) * (t100 ^ 2 + t101 ^ 2) + (t8 ^ 2 + t9 ^ 2) * m(11) + (t6 ^ 2 + t7 ^ 2) * m(4) + (t31 ^ 2 + t32 ^ 2) * m(5) + t5 ^ 2 * Icges(10,3) + t17 ^ 2 * Icges(11,3) + t14 ^ 2 * Icges(4,3) + t38 ^ 2 * Icges(5,3) + t35 ^ 2 * Icges(7,3) + t91 ^ 2 * Icges(3,3) + t107 ^ 2 * Icges(9,3) + m(2) * (t278 ^ 2 + t280 ^ 2) + Icges(2,3) + (m(6) * (t143 ^ 2 + t144 ^ 2) + Icges(6,3)) * t111 ^ 2; (t1 * t3 + t2 * t4) * m(10) + t5 * Icges(10,3) * t15 + t17 * Icges(11,3) * t18 + (t10 * t8 + t11 * t9) * m(11) + t14 * Icges(4,3) * t16 + (t12 * t6 + t13 * t7) * m(4) + (t31 * t33 + t32 * t34) * m(5) + t38 * Icges(5,3) * t39 + ((-t133 * t75 + t134 * t76) * m(3) + t91 * Icges(3,3)) * t90 + ((-t29 * t55 + t30 * t56) * m(7) + t35 * Icges(7,3)) * t37; t16 ^ 2 * Icges(4,3) + (t12 ^ 2 + t13 ^ 2) * m(4) + t15 ^ 2 * Icges(10,3) + (t10 ^ 2 + t11 ^ 2) * m(11) + t18 ^ 2 * Icges(11,3) + t39 ^ 2 * Icges(5,3) + m(8) * (t277 ^ 2 + t279 ^ 2) + (t33 ^ 2 + t34 ^ 2) * m(5) + (t3 ^ 2 + t4 ^ 2) * m(10) + Icges(8,3) + (Icges(3,3) + (t133 ^ 2 + t134 ^ 2) * m(3)) * t90 ^ 2 + (Icges(7,3) + (t55 ^ 2 + t56 ^ 2) * m(7)) * t37 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t19(1), t19(2); t19(2), t19(3);];
Mq = res;