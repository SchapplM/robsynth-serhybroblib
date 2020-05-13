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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-10 08:43
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1TE_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1TE_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1TE_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 00:16:10
% EndTime: 2020-05-10 00:18:12
% DurationCPUTime: 108.81s
% Computational Cost: add. (2403244->845), mult. (6593214->1480), div. (87768->44), fcn. (1768036->14), ass. (0->645)
t758 = 4 * pkin(1);
t342 = sin(qJ(2));
t296 = pkin(3) * t342;
t585 = 0.8e1 * t296;
t307 = t342 ^ 2;
t306 = t342 * t307;
t391 = pkin(3) ^ 2;
t408 = pkin(3) * t391;
t647 = t306 * t408;
t546 = 0.32e2 * t647;
t385 = pkin(4) ^ 2;
t343 = sin(qJ(1));
t344 = cos(qJ(2));
t639 = t343 * t344;
t272 = pkin(3) * t639;
t481 = pkin(1) * t272;
t398 = pkin(1) ^ 2;
t401 = pkin(7) ^ 2;
t615 = t398 / 0.3e1 + t401;
t210 = -0.4e1 / 0.9e1 * t481 + 0.4e1 / 0.9e1 * t391 - t385 / 0.9e1 + t615;
t322 = -t385 / 0.6e1;
t331 = 0.2e1 / 0.3e1 * t391;
t444 = t401 - t481;
t222 = t322 + t331 + t444;
t330 = 0.4e1 / 0.3e1 * t391;
t302 = t398 + t401;
t324 = -t385 / 0.3e1;
t513 = t324 + t302;
t261 = t330 + t513;
t294 = -t398 / 0.3e1 + t401;
t485 = -0.2e1 * t272;
t345 = cos(qJ(1));
t311 = t345 ^ 2;
t626 = t398 * t311;
t563 = 0.6e1 * t626;
t310 = t345 * t311;
t405 = pkin(1) * t398;
t643 = t310 * t405;
t575 = pkin(7) * t643;
t721 = pkin(7) * t345;
t596 = 0.4e1 * t721;
t633 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t152 = 0.4e1 * t575 + t210 * t563 + t261 * t633 + (t222 * t596 + t294 * t485) * pkin(1);
t364 = 6 * t398;
t457 = -0.4e1 * t481;
t326 = -0.2e1 / 0.3e1 * t385;
t373 = 0.2e1 * t401;
t512 = t326 + t331 + t373;
t409 = t391 ^ 2;
t511 = t326 + t302;
t619 = (t331 + t511) * t302 + t409;
t166 = t261 * t457 + (t364 + t512) * t391 + t619;
t259 = -0.2e1 * t481;
t334 = -t391 / 0.3e1;
t292 = t334 + t401;
t230 = t292 * t259;
t634 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t180 = t261 * t634 + t230;
t313 = 0.10e2 / 0.3e1 * t398;
t181 = (t313 + t512) * t391 + t619;
t301 = -0.3e1 * t391 + t401;
t494 = 0.8e1 * t575;
t245 = t301 * t494;
t509 = t391 + t302;
t270 = -t385 + t509;
t298 = pkin(1) * t345;
t277 = t298 + pkin(7);
t295 = t302 ^ 2;
t396 = t398 ^ 2;
t347 = 15 * t396;
t354 = 0.18e2 * t401;
t355 = -0.2e1 * t385;
t357 = -0.6e1 * t385;
t359 = 0.2e1 * t391;
t400 = t401 ^ 2;
t368 = 0.3e1 * t400;
t388 = t408 ^ 2;
t218 = t259 + t261;
t604 = -t391 + t401;
t755 = 4 * t398;
t282 = t604 * t755;
t589 = 0.4e1 * t298;
t154 = pkin(7) * t218 * t589 + t282 * t311 + t166;
t590 = 0.2e1 * t298;
t281 = pkin(7) * t590;
t304 = -(3 * t398) + t401;
t564 = 0.4e1 * t626;
t229 = t281 + t564 + t304;
t572 = 0.8e1 * t647;
t586 = 0.6e1 * t296;
t434 = t154 * t586 + t229 * t572;
t540 = 0.12e2 * t626;
t645 = t307 * t391;
t544 = 0.12e2 * t645;
t372 = 0.3e1 * t401;
t612 = (15 * t398) + t372;
t757 = 6 * pkin(1);
t121 = t152 * t544 + t245 + t180 * t540 + t388 + (-t385 + t391 + t612) * t409 + (t347 + (t354 + t357 + 0.6e1 * t391) * t398 + t368 + (t355 + t359) * t401) * t391 + t295 * t270 + t434 * t277 + (t166 * t721 - t181 * t272) * t757;
t297 = pkin(3) * t344;
t550 = t398 * t297;
t439 = -t343 * t405 + t550;
t240 = 0.2e1 * t439;
t584 = 0.2e1 * t297;
t279 = pkin(7) * t584;
t366 = 2 * t398;
t300 = t366 + t391;
t303 = -t398 + t401;
t168 = t303 * t297 + t240 * t311 + (t279 * t345 + t300 * t343) * pkin(1);
t365 = 3 * t398;
t608 = t365 + t401;
t286 = t391 + t608;
t244 = t286 * t297;
t287 = 0.3e1 * t391 + t302;
t651 = t287 * t343;
t200 = -pkin(1) * t651 + t244;
t754 = 8 * t398;
t562 = t408 * t754;
t699 = pkin(3) * t396;
t275 = t562 + 0.4e1 * t699;
t526 = t405 * t634;
t203 = t275 * t344 + 0.4e1 * t343 * t526;
t363 = 5 * t396;
t601 = t400 + t409;
t349 = 10 * t398;
t611 = t349 + t373;
t624 = t401 * t398;
t224 = t391 * t611 + t363 + t601 + 0.6e1 * t624;
t358 = 0.5e1 * t409;
t370 = 0.6e1 * t401;
t235 = t358 + (t349 + t370) * t391 + t295;
t530 = t277 * t647;
t476 = -0.8e1 * t530;
t571 = -0.4e1 * t645;
t652 = t277 * t342;
t700 = pkin(1) * t343;
t248 = t297 - t700;
t473 = t311 * t550;
t597 = 0.2e1 * t721;
t165 = -0.2e1 * t473 + t244 + (t248 * t597 - t651) * pkin(1);
t736 = -0.4e1 * t165;
t131 = t168 * t571 + t203 * t311 + (t652 * t736 + (-t224 + t494) * t344) * pkin(3) + (-0.4e1 * t200 * t721 + (t235 + t476) * t343) * pkin(1);
t242 = t296 + t277;
t733 = 0.2e1 * t342;
t598 = pkin(7) * t733;
t278 = pkin(3) * t598;
t600 = t401 - t385;
t285 = t398 + t600;
t239 = t278 + t285;
t205 = t259 + t239;
t276 = t296 + pkin(7);
t256 = t276 * t345;
t312 = t391 * t755;
t630 = t391 * t401;
t283 = t312 - 0.4e1 * t630;
t374 = -0.2e1 * t401;
t375 = 0.2e1 * pkin(3);
t536 = -0.4e1 * pkin(3) * pkin(7) * t285;
t228 = t278 + t604 + 0.2e1 * t645;
t655 = t228 * t311;
t141 = t283 * t307 + t342 * t536 - t396 - (t401 - (t375 + pkin(4)) * pkin(4)) * (t401 + (t375 - pkin(4)) * pkin(4)) + (t374 + 0.2e1 * t385 - 0.4e1 * t391 - 0.4e1 * t655) * t398 + (-t205 * t256 + t239 * t272) * t758;
t402 = sqrt(t141);
t114 = t121 * t242 + t131 * t402;
t760 = t114 ^ 2;
t759 = 2 * pkin(1);
t376 = 2 * pkin(2);
t384 = t385 ^ 2;
t603 = t396 + t400;
t607 = t373 - t385;
t625 = t401 * t385;
t435 = t607 * t398 + t384 / 0.6e1 + t603 - t625;
t429 = 0.5e1 / 0.6e1 * t409 + t435;
t188 = (t313 + t607) * t391 + t429;
t756 = -0.6e1 * t188;
t346 = cos(pkin(9));
t637 = t344 * t345;
t551 = pkin(3) * t637;
t638 = t343 * t345;
t653 = t276 * t343;
t654 = t276 * t311;
t150 = (t228 * t638 + t297 * t654) * t754 + (t205 * t653 + t239 * t551) * t758;
t269 = t359 + t365 + t600;
t201 = t269 + t278 + t457;
t213 = t551 + t653;
t140 = 0.1e1 / t402;
t712 = t140 / 0.2e1;
t501 = t213 * t712;
t558 = t228 * t298;
t617 = -t272 + t256;
t640 = t342 * t391;
t120 = t617 * t402 + t150 * t501 + (t201 * t276 + 0.4e1 * t558) * t343 + (t597 * t640 + (t269 * t345 + t654 * t758) * pkin(3)) * t344;
t340 = cos(pkin(8));
t682 = sin(pkin(8));
t226 = -t340 * t345 - t343 * t682;
t698 = pkin(5) * t226;
t742 = -2 * pkin(1);
t216 = t698 * t742;
t381 = pkin(5) ^ 2;
t207 = t216 + 0.2e1 * t381;
t214 = -pkin(1) + t698;
t225 = t340 * t343 - t345 * t682;
t621 = t216 + t381;
t193 = -(t759 + pkin(5)) * pkin(5) + t621;
t194 = pkin(5) * (t759 - pkin(5)) + t621;
t403 = sqrt(-t193 * t194);
t492 = t225 * t759;
t665 = 0.1e1 / t403 * (t193 + t194) * pkin(5) * t492;
t500 = t665 / 0.2e1;
t656 = t226 * t403;
t123 = (-t656 + (t214 * t759 - t207 + t500) * t225) * pkin(5);
t657 = t225 * t403;
t705 = 0.2e1 * t225 ^ 2;
t124 = -t214 * t665 / 0.2e1 + t381 * pkin(1) * t705 + (t207 * t226 - t657) * pkin(5);
t147 = -pkin(5) * t225 * t207 - t214 * t403;
t556 = pkin(1) * t256;
t178 = t259 + t278 + t509 + 0.2e1 * t556;
t176 = 0.1e1 / t178;
t392 = 0.1e1 / pkin(3);
t237 = t359 + t239;
t170 = t237 * t276 + 0.2e1 * t558;
t172 = t237 * t345 + (0.4e1 * t311 - 0.2e1) * t276 * pkin(1);
t211 = pkin(1) + t617;
t130 = t170 * t343 + t172 * t297 + t211 * t402;
t667 = t130 * t147;
t257 = t272 - pkin(1);
t658 = t213 * t402;
t128 = -t201 * t256 + t658 + (t257 * t598 + t269 * t639) * pkin(3) + (-0.2e1 * t655 + (0.2e1 * t307 - 0.4e1) * t391 - t285) * pkin(1);
t145 = pkin(5) * t657 - t207 * t214;
t669 = t128 * t145;
t437 = t669 / 0.4e1 + t667 / 0.4e1;
t202 = t398 + t621;
t198 = 0.1e1 / t202;
t661 = t198 / pkin(5);
t177 = 0.1e1 / t178 ^ 2;
t204 = t213 * t759;
t663 = t177 * t204;
t464 = t661 * t663;
t199 = 0.1e1 / t202 ^ 2;
t659 = t199 * t225;
t559 = pkin(1) * t659;
t711 = t145 / 0.4e1;
t713 = t130 / 0.4e1;
t715 = t128 / 0.4e1;
t308 = t343 ^ 2;
t502 = t211 * t712;
t118 = -t658 + t150 * t502 + t228 * t308 * t742 + t170 * t345 + (-t237 - 0.8e1 * t556) * t272;
t717 = t118 / 0.4e1;
t62 = (t437 * t464 + ((t669 / 0.2e1 + t667 / 0.2e1) * t559 + (t120 * t711 + t123 * t715 + t124 * t713 + t147 * t717) * t661) * t176) * t392;
t666 = t145 * t130;
t668 = t128 * t147;
t436 = t668 / 0.4e1 - t666 / 0.4e1;
t63 = (-t436 * t464 + ((-t668 / 0.2e1 + t666 / 0.2e1) * t559 + (-t120 * t147 / 0.4e1 - t128 * t124 / 0.4e1 + t123 * t713 + t118 * t711) * t661) * t176) * t392;
t704 = sin(pkin(9));
t41 = t346 * t62 + t63 * t704;
t378 = pkin(6) ^ 2;
t531 = t392 * t661;
t465 = t176 * t531;
t110 = t436 * t465;
t425 = t437 * t465;
t102 = -t110 * t704 + t346 * t425;
t697 = pkin(6) * t102;
t97 = t697 * t376;
t672 = t378 + t97;
t90 = -(t376 + pkin(6)) * pkin(6) + t672;
t91 = pkin(6) * (t376 - pkin(6)) + t672;
t443 = pkin(6) * (-t90 - t91) * t376;
t29 = t41 * t443;
t686 = t90 * t91;
t404 = sqrt(-t686);
t42 = -t346 * t63 + t62 * t704;
t104 = -t110 * t346 - t425 * t704;
t696 = pkin(6) * t104;
t95 = t97 + 0.2e1 * t378;
t96 = -pkin(2) - t697;
t46 = t404 * t696 - t95 * t96;
t423 = t46 ^ 2;
t44 = 0.1e1 / t423;
t47 = -t404 * t96 - t696 * t95;
t45 = t47 ^ 2;
t394 = pkin(2) ^ 2;
t94 = t394 + t672;
t581 = pkin(6) / (t44 * t45 + 0.1e1) * t94;
t448 = -0.2e1 * t44 * t47 * t581;
t466 = 0.2e1 / t46 * t581;
t379 = 0.1e1 / pkin(6);
t92 = 0.1e1 / t94;
t675 = t379 * t92;
t518 = t675 / 0.2e1;
t470 = pkin(6) * t518;
t740 = -2 * pkin(2);
t490 = t104 * t378 * t740;
t508 = t740 * t96 + t95;
t73 = 0.1e1 / t404;
t523 = t104 * t73 / 0.2e1;
t537 = -t73 * t96 / 0.2e1;
t93 = 0.1e1 / t94 ^ 2;
t726 = pkin(2) * t93;
t583 = t41 * t726;
t417 = t128 ^ 2;
t127 = 0.1e1 / t417;
t129 = t130 ^ 2;
t474 = pkin(3) / (t127 * t129 + 0.1e1) * t178 * t392;
t438 = -0.2e1 * t127 * t130 * t474;
t716 = t128 / 0.2e1;
t441 = t474 / t716;
t497 = t663 / 0.2e1;
t710 = t176 / 0.2e1;
t71 = 0.1e1 + (t118 * t710 + t130 * t497) * t441 + (t120 * t710 + t128 * t497) * t438;
t5 = ((t29 * t537 + t41 * t490 + (t404 * t41 + t42 * t95) * pkin(6)) * t518 - t47 * t583) * t466 + (-t46 * t583 + (t29 * t523 - t404 * t42 + t41 * t508) * t470) * t448 + t71;
t753 = Ifges(4,3) * t5;
t484 = pkin(1) * t551;
t642 = t342 * t343;
t553 = pkin(3) * t642;
t482 = pkin(1) * t553;
t618 = 0.2e1 * t482 + t279;
t197 = 0.2e1 * t484 + t618;
t664 = t177 * t197;
t498 = -t664 / 0.4e1;
t636 = t344 * t391;
t527 = t342 * t636;
t234 = t279 + 0.4e1 * t527;
t566 = -0.4e1 * t626;
t722 = pkin(7) * t343;
t577 = t391 * t722;
t138 = t234 * t566 + (t283 * t733 + t536) * t344 + (-t618 * t256 + 0.2e1 * t344 ^ 2 * t577 + (-t205 * t637 - t239 * t642) * pkin(3)) * t758;
t641 = t342 * t345;
t552 = pkin(3) * t641;
t701 = pkin(1) * t311;
t743 = 0.2e1 * pkin(7);
t117 = (t272 - t552) * t402 + t138 * t501 - 0.2e1 * t234 * t701 - (t279 + 0.4e1 * t482) * t256 - 0.2e1 * t307 * t577 - t269 * t553 + (t640 * t758 + (-t201 * t345 + t257 * t743) * pkin(3)) * t344;
t670 = t117 * t176;
t430 = t670 / 0.4e1 + t128 * t498;
t116 = t138 * t502 + t234 * t343 * t590 + ((t343 * t402 - t172) * t342 + (t345 * t402 + (t276 * t743 + t237) * t343 + (-pkin(1) + 0.2e1 * t701 + t721) * t584) * t344) * pkin(3);
t671 = t116 * t176;
t431 = t671 / 0.4e1 + t130 * t498;
t79 = (t145 * t430 + t147 * t431) * t531;
t80 = (t145 * t431 - t147 * t430) * t531;
t59 = t346 * t79 + t704 * t80;
t39 = t59 * t443;
t582 = t59 * t726;
t60 = -t346 * t80 + t704 * t79;
t499 = -t664 / 0.2e1;
t70 = (t671 / 0.2e1 + t130 * t499) * t441 + (t670 / 0.2e1 + t128 * t499) * t438;
t10 = ((t39 * t537 + t59 * t490 + (t404 * t59 + t60 * t95) * pkin(6)) * t518 - t47 * t582) * t466 + (-t46 * t582 + (t39 * t523 - t404 * t60 + t508 * t59) * t470) * t448 + t70;
t737 = 0.1e1 / t102 ^ 2;
t507 = t737 * t404 / 0.2e1;
t739 = 0.1e1 / t102;
t524 = -t739 * t73 / 0.4e1;
t395 = 1 / pkin(2);
t380 = 0.1e1 / pkin(6) ^ 2;
t708 = -t380 / 0.4e1;
t539 = t379 * t395 / (0.1e1 + 0.1e1 / (pkin(2) ^ 2) * t737 * t686 * t708);
t6 = (t39 * t524 + t507 * t59) * t539 + t10;
t752 = Ifges(10,3) * t6;
t386 = 0.1e1 / pkin(4);
t629 = t392 * t402;
t253 = t401 + t391 / 0.4e1 + t398 / 0.4e1 - t385 / 0.8e1;
t613 = 0.4e1 / 0.7e1 * t401 - t385 / 0.7e1;
t163 = -0.32e2 / 0.21e2 * t253 * t481 + 0.5e1 / 0.42e2 * t409 + (0.16e2 / 0.21e2 * t398 + t613) * t391 + t396 / 0.7e1 + t613 * t398 + t400 - 0.3e1 / 0.7e1 * t625 + t384 / 0.42e2;
t323 = -t385 / 0.4e1;
t746 = t323 + t391 / 0.2e1;
t255 = t615 + t746;
t336 = 0.4e1 / 0.3e1 * t398;
t164 = -0.8e1 / 0.3e1 * t255 * t481 + 0.5e1 / 0.18e2 * t409 + (t336 + t324) * t391 + t400 - t396 / 0.3e1 + t384 / 0.18e2 + (t330 + 0.2e1 / 0.3e1 * t398 + t326) * t401;
t219 = -t409 / 0.6e1 + t435;
t338 = t398 / 0.2e1;
t614 = t338 + t401;
t221 = -0.2e1 / 0.3e1 * t481 + t323 + t614;
t371 = 0.4e1 * t401;
t289 = (t371 + t385) * t398;
t335 = -0.2e1 / 0.3e1 * t391;
t293 = t335 + t401;
t325 = -t385 / 0.2e1;
t266 = t325 + t509;
t486 = -0.4e1 * t272;
t449 = t266 * t486;
t487 = 0.16e2 * t575;
t309 = t311 ^ 2;
t627 = t396 * t309;
t567 = 0.8e1 * t627;
t134 = t293 * t567 + t221 * t487 + 0.14e2 * t163 * t626 + t303 * t409 + (t289 - 0.10e2 / 0.3e1 * t396 + 0.2e1 * t400 - t625) * t391 + t219 * t633 + (0.6e1 * t164 * t721 + t294 * t449) * pkin(1);
t327 = -0.3e1 / 0.2e1 * t385;
t616 = t384 / 0.2e1 - t409 / 0.2e1;
t456 = -0.3e1 * t625 + t368 + t616;
t622 = t302 * ((t327 + t373) * t398 - 0.3e1 / 0.2e1 * t625 + t603 + t616) + t388;
t149 = t481 * t756 + (t347 + (t354 - 0.9e1 * t385) * t398 + t456) * t391 + (t327 + t612) * t409 + t622;
t445 = pkin(1) * t449;
t602 = t400 - t396;
t155 = t292 * t445 - t388 + (-t313 - t600) * t409 + (t289 + t409 / 0.6e1 - t384 / 0.6e1 + t602) * t391 + t219 * t401;
t356 = -0.5e1 * t385;
t362 = 7 * t396;
t159 = (t327 + t372 + (7 * t398)) * t409 + (t362 + (t356 + 0.10e2 * t401) * t398 + t456) * t391 + t622;
t413 = pkin(7) * t401;
t284 = -0.12e2 * pkin(7) * t405 + t413 * t758;
t291 = -(8 * t396) + 0.12e2 * t624;
t173 = t284 * t345 + t291 * t311 + t487 + t567 + t603 - 0.6e1 * t624;
t184 = t259 * t634 + t266 * t301;
t246 = 0.16e2 * (t601 - 0.6e1 * t630) * t396;
t299 = -0.30e2 * t385 + 0.60e2 * t401;
t305 = t307 ^ 2;
t606 = t384 - t409;
t440 = 0.6e1 * t400 + t606 - 0.6e1 * t625;
t650 = t295 * (-t391 + t285);
t169 = t445 + (t364 + t607) * t391 + t429;
t182 = t266 * t634 + t230;
t592 = pkin(7) * t298;
t534 = 0.6e1 * t592;
t135 = t169 * t534 + t182 * t540 + t149 + t245;
t369 = 0.8e1 * t401;
t206 = t405 * t486 + t312 + (4 * t396) + (t355 + t369) * t398;
t217 = -t398 + t444 + t746;
t153 = t494 + t206 * t311 + t266 * t304 + (t217 * t596 + t485 * t633) * pkin(1);
t744 = t135 * t585 + t153 * t546;
t108 = t246 * t309 + 0.32e2 * t184 * t575 + 0.24e2 * t155 * t626 + (t355 + t371 + (28 * t398)) * t388 + t270 * t650 + (0.24e2 * t134 * t307 + t299 * t396 + t357 * t400 + t440 * t364 + t606 * t373 + (28 * t405 ^ 2) + 0.4e1 * t413 ^ 2) * t391 + t744 * t277 + 0.8e1 * (t149 * t721 - t159 * t272) * pkin(1) + (0.16e2 * t173 * t305 + t299 * t398 + (70 * t396) + t409 + t440) * t409;
t251 = 0.7e1 / 0.6e1 * t391 + t322 + t614;
t332 = t391 / 0.3e1;
t514 = t322 + t332 + t401;
t254 = t336 + t514;
t179 = -t251 * t700 + t254 * t297;
t263 = t398 + t514;
t288 = t359 + t303;
t190 = t263 * t297 - t288 * t700 / 0.2e1;
t515 = t385 / 0.3e1 + t332 + t373;
t455 = -0.8e1 / 0.3e1 * t627 + t391 * t303 - 0.5e1 / 0.3e1 * t396 + t515 * t398 + t401 * (t324 + t292);
t720 = pkin(7) * t396;
t595 = -0.4e1 * t720;
t314 = -0.20e2 / 0.3e1 * t398;
t516 = 0.2e1 / 0.3e1 * t385 + t331 + t371;
t517 = 0.4e1 / 0.3e1 * t385 + t330 + t374;
t709 = t409 / 0.2e1 - (t314 + t516) * t391 / 0.2e1 + 0.3e1 / 0.2e1 * t396 - t517 * t398 / 0.2e1 - t400 / 0.2e1;
t136 = t343 * t310 * t595 + t179 * t564 + t455 * t297 + (t190 * t596 + t343 * t709) * pkin(1);
t610 = t355 - 0.2e1 * t391;
t510 = t370 + t610;
t175 = t409 + (t326 + t335 + t611) * t391 + t363 + t510 * t398 + t401 * (t326 + t293);
t171 = t175 * t297;
t260 = 0.8e1 / 0.3e1 * t391 + t513;
t262 = t324 + t331 + t608;
t183 = -t260 * t700 + t262 * t297;
t267 = 0.5e1 / 0.6e1 * t391 + t338 + t322;
t196 = t267 * t584 + t634 * t700;
t549 = t405 * t297;
t472 = t310 * t549;
t189 = t358 + (t349 + t510) * t391 + (t335 + t511) * t302;
t662 = t189 * t343;
t137 = -0.8e1 * pkin(7) * t472 + t196 * t566 + t171 + (t183 * t596 - t662) * pkin(1);
t156 = -pkin(1) * t662 + t171;
t187 = -0.3e1 * t409 + (t314 + t517) * t391 + t516 * t398 + t602;
t191 = -0.5e1 / 0.3e1 * t409 + (-t398 + t515) * t391 + t401 * (t334 + t513);
t591 = -0.2e1 * t700;
t157 = t187 * t297 + t191 * t591;
t241 = 0.4e1 * t439;
t247 = t297 + 0.2e1 * t700;
t265 = t325 + t286;
t160 = t304 * t297 + t241 * t311 + (t247 * t721 + t265 * t343) * t759;
t252 = t401 + 0.5e1 / 0.2e1 * t391 + 0.3e1 / 0.2e1 * t398 + t325;
t649 = t301 * t343;
t192 = t252 * t297 + pkin(1) * t649 / 0.2e1;
t745 = t372 - t385 - t391;
t268 = t745 * t349;
t609 = t356 - 0.5e1 * t391;
t447 = 0.24e2 * t292 * t627 - t388 - ((21 * t398) + t745) * t409 - (t401 * t610 + t268 + t368 + (35 * t396)) * t391 - (t362 + (t369 + t609) * t398 + t401 * (-t391 + t600)) * t302;
t545 = -0.12e2 * t645;
t223 = 0.4e1 / 0.3e1 * t626 + t281 + t294;
t648 = t305 * t409;
t747 = 0.7e1 * t388 + ((35 * t398) + 0.15e2 * t401 + t609) * t409 + ((21 * t396) + t268 + 0.9e1 * t400 + (t357 - 0.6e1 * t391) * t401) * t391 + t650 - 0.24e2 * t223 * t648;
t115 = t192 * t487 + t160 * t476 + t136 * t545 - 0.6e1 * t157 * t626 + (-0.6e1 * t137 * t652 + t344 * t447) * pkin(3) + (-0.6e1 * t156 * t721 + t747 * t343) * pkin(1);
t87 = t108 * t242 + t115 * t402;
t424 = t87 ^ 2;
t678 = t760 / t424;
t78 = 0.1e1 / (t141 * t678 + 0.1e1);
t428 = -0.2e1 * pkin(3) * pkin(4) * t386 * t78 * t629 * t678;
t446 = t114 * t78 / t87 * t712;
t113 = 0.1e1 / t760;
t676 = t113 * t87;
t520 = -t676 / 0.2e1;
t451 = pkin(7) * t473;
t250 = -0.4e1 * t451;
t280 = pkin(7) * t591;
t723 = pkin(7) * t311;
t452 = t549 * t723;
t493 = 0.32e2 / 0.3e1 * t396;
t453 = t310 * t493;
t454 = 0.64e2 / 0.3e1 * t253 * t405;
t463 = t310 * t526;
t468 = 0.32e2 * t530;
t478 = t634 * t720;
t623 = t405 * t311;
t578 = pkin(7) * t623;
t488 = -0.24e2 * t578;
t489 = -0.48e2 * t578;
t495 = pkin(7) * t564;
t503 = t150 * t712;
t635 = t345 * t398;
t525 = t343 * t635;
t541 = -0.24e2 * t635;
t542 = -0.32e2 * t310 * t396;
t543 = 0.24e2 * t645;
t547 = -0.96e2 * t292 * t310;
t555 = pkin(1) * t633;
t702 = pkin(1) * t294;
t557 = t266 * t702;
t565 = -0.2e1 * t626;
t568 = 0.8e1 * t635;
t569 = -0.4e1 * t643;
t573 = -0.8e1 * t647;
t579 = pkin(7) * t626;
t588 = -0.6e1 * t296;
t620 = -0.24e2 * t292 * t472 + t488 * t649;
t730 = pkin(1) * pkin(7);
t673 = -0.4e1 * t730;
t674 = -0.6e1 * t730;
t695 = t159 * pkin(1);
t731 = -0.6e1 * t402;
t732 = -0.2e1 * t345;
t735 = 0.8e1 * t242;
t61 = ((t265 * t590 + t495 + t569) * t476 + (0.12e2 * t308 * t311 * t720 + t251 * t569 + t309 * t595) * t545 + 0.12e2 * t191 * t643 + (t545 * t709 + t747) * t298 + (t288 * t545 * t565 + t189 * t563 + t301 * t567) * pkin(7)) * t402 + t115 * t503 + t250 * t242 * t468 + (((pkin(7) * t260 * t566 - t189 * t298 - 0.4e1 * t463) * t731 + t620 * t735) * t652 + ((-pkin(7) * t309 * t493 - 0.16e2 * t255 * t579 - t310 * t454 - 0.4e1 * t345 * t557) * t543 - 0.64e2 * t309 * t478 + t405 * t266 * t547 - 0.48e2 * t188 * t579 - 0.8e1 * t345 * t695 + ((t555 * t732 + t569) * t546 + (-0.24e2 * t266 * t579 + t298 * t756) * t585) * t277) * t242 * t344) * pkin(3) + ((0.16e2 * (t291 * t732 - t284 + t489 + t542) * t648 + (t206 * t732 + t217 * t673 + t488) * t468 + (-0.28e2 * t163 * t635 + t164 * t674 + t221 * t489 + t293 * t542) * t543 + t277 * (t169 * t674 + t182 * t541) * t585 - 0.4e1 * t246 * t310 - 0.96e2 * t184 * t578 - 0.48e2 * t155 * t635 - 0.8e1 * t149 * t730) * t242 + ((-0.8e1 * t179 * t635 + t297 * t453) * t545 + t344 * t547 * t699 + t192 * t489 + 0.12e2 * t157 * t635 + (0.2e1 * (-t241 * t345 - t247 * t730) * t573 + (t183 * t673 + t196 * t568 + 0.24e2 * t452) * t588) * t277) * t402 + (-t744 * t242 - t108 + (t137 * t586 + t160 * t572 + (-0.24e2 * t280 + 0.64e2 * t525) * t648 + (0.48e2 * t190 * t645 + 0.6e1 * t156) * pkin(7)) * t402) * pkin(1)) * t343;
t112 = 0.1e1 / t114;
t718 = t112 / 0.2e1;
t570 = 0.8e1 * t645;
t587 = -0.4e1 * t296;
t88 = (t308 * t306 * t562 + (t298 * t300 - 0.2e1 * t643) * t571 + 0.4e1 * t463 + t287 * t495 + t235 * t298) * t402 + t131 * t503 + ((-0.8e1 / 0.3e1 * t472 + t250 - 0.2e1 * t294 * t484) * t544 - 0.24e2 * t261 * t451 - 0.6e1 * t181 * t484 + t620) * t242 + ((t240 * t345 * t570 + t203 * t732 - 0.24e2 * t452) * t402 + (0.12e2 * (-t210 * t635 - t578) * t544 + t180 * t541) * t242 + (0.4e1 * t165 * t402 * t296 - t434 * t242 - t121 + ((t297 * t570 + 0.4e1 * t200) * t402 + (-0.48e2 * t222 * t645 - 0.6e1 * t166) * t242) * pkin(7)) * pkin(1)) * t343 + ((t573 * t298 + ((0.4e1 * t272 * t345 - 0.2e1 * t723) * t398 + (-0.2e1 * t248 * t722 - t287 * t345) * pkin(1)) * t587) * t402 + ((t280 - 0.8e1 * t525) * t572 + (-0.8e1 * t451 - 0.2e1 * t282 * t638 + (-t218 * t722 - t261 * t551) * t758) * t586) * t242) * t277;
t27 = t150 * t446 + (t520 * t88 + t61 * t718) * t428 + t71;
t751 = Ifges(5,3) * t27;
t679 = t104 * t737;
t74 = 0.1e1 / (t104 ^ 2 * t737 + 0.1e1);
t21 = (t41 * t679 + t42 * t739) * t74 + t71;
t750 = Ifges(7,3) * t21;
t749 = t71 * Ifges(3,3);
t231 = t637 + t642;
t232 = t639 - t641;
t528 = t292 * t623;
t220 = 0.24e2 * t528 * t553;
t491 = pkin(7) * t552;
t442 = t398 * t343 * t491;
t238 = 0.4e1 * t442;
t249 = t491 * t742;
t504 = t138 * t712;
t644 = t307 * t408;
t529 = t277 * t644;
t535 = -0.4e1 * t592;
t576 = pkin(7) * t635;
t646 = t306 * t409;
t64 = t115 * t504 + (0.32e2 * t238 * t242 - 0.8e1 * t249 * t402) * t530 + (0.24e2 * (-0.4e1 * t223 * t646 * t700 - t136 * t640 - t160 * t529) * t402 + (0.48e2 * t134 * t640 + 0.96e2 * t153 * t529 + 0.64e2 * t173 * t646) * t242) * t344 + ((t108 + (t135 * t735 + t137 * t731) * t277) * t344 + (t277 * t220 * t735 + ((t254 * t566 + t263 * t535 - t455) * t545 - 0.16e2 * t252 * t575 + t187 * t563 + t175 * t534 + ((-t304 + t566) * t573 + (t262 * t535 + 0.8e1 * t267 * t626 - t175 + t494) * t588) * t277 - t447) * t402 + ((pkin(7) * t453 + 0.16e2 * t255 * t576 + t311 * t454 + 0.4e1 * t557) * t543 + 0.64e2 * t310 * t478 + 0.96e2 * t266 * t528 + 0.48e2 * t188 * t576 + 0.8e1 * t695 + ((0.2e1 * t555 + 0.4e1 * t623) * t546 + (t188 * t757 + 0.24e2 * t266 * t576) * t585) * t277) * t242 * t343) * t342) * pkin(3);
t89 = (t249 * t571 + (-0.8e1 * t168 * t636 - t275 * t311 + ((-t303 + t565) * t571 + t224 + (t286 * t589 - 0.8e1 * t643) * pkin(7)) * pkin(3)) * t342 + ((t249 + (-t286 + 0.2e1 * t626) * t296) * t587 + (pkin(3) * t736 - 0.24e2 * t644 * t700) * t344) * t277) * t402 + t131 * t504 + (0.24e2 * t229 * t344 * t529 + (t238 + (0.8e1 / 0.3e1 * t623 + 0.2e1 * t702) * t553) * t544 + 0.24e2 * t152 * t527 + t220 + 0.24e2 * t261 * t442 + 0.6e1 * t181 * t482 + 0.6e1 * ((pkin(7) * t568 + t261 * t758) * t343 * t645 + t154 * t297) * t277) * t242 + t121 * t297;
t28 = t138 * t446 + (t520 * t89 + t64 * t718) * t428 + t70;
t677 = t112 * t87;
t521 = t677 / 0.4e1;
t707 = -t402 / 0.4e1;
t432 = t128 * t521 + t130 * t707;
t426 = t177 * t432;
t519 = -t676 / 0.4e1;
t461 = t128 * t519;
t505 = -t130 * t140 / 0.8e1;
t393 = 0.1e1 / pkin(3) ^ 2;
t632 = t386 * t393;
t727 = t87 / 0.4e1;
t32 = (-t197 * t426 + (t89 * t461 + t116 * t707 + t138 * t505 + (t117 * t727 + t64 * t715) * t112) * t176) * t632;
t458 = t130 * t521;
t706 = t402 / 0.4e1;
t433 = t128 * t706 + t458;
t532 = t176 * t632;
t77 = t433 * t532;
t75 = t231 * t77;
t76 = t432 * t532;
t53 = t232 * t76 + t75;
t51 = 0.1e1 / t53 ^ 2;
t52 = -t231 * t76 + t232 * t77;
t38 = 0.1e1 / (t51 * t52 ^ 2 + 0.1e1);
t50 = 0.1e1 / t53;
t427 = t177 * t433;
t460 = t130 * t519;
t506 = t128 * t140 / 0.8e1;
t683 = (-t197 * t427 + (t117 * t706 + t138 * t506 + t89 * t460 + (t116 * t727 + t64 * t713) * t112) * t176) * t632 - t76;
t689 = t51 * t52;
t12 = (-(-t231 * t32 + t232 * t683 - t75) * t50 - ((-t32 - t77) * t232 - t683 * t231) * t689) * t38 + t28;
t748 = Ifges(11,3) * t12;
t738 = 0.1e1 / t145;
t703 = pkin(1) * t176;
t560 = t392 * t703;
t471 = t128 * t560;
t125 = t471 / 0.2e1;
t729 = -pkin(2) * t71 / 0.2e1 - t125 / 0.2e1;
t574 = t47 * t70 * t92;
t580 = pkin(2) * t675;
t483 = t70 * t580;
t450 = t46 * t483;
t9 = -t450 / 0.2e1 + t10 * pkin(6);
t8 = t404 * t574 * t708 + t102 * t9;
t725 = mrSges(10,1) * t8;
t628 = t395 * t404;
t496 = -t628 / 0.2e1;
t719 = -t102 / 0.2e1;
t7 = (pkin(2) * t574 * t719 + t496 * t9) * t379;
t724 = mrSges(10,2) * t7;
t714 = -t130 / 0.2e1;
t467 = -t130 * t703 / 0.4e1;
t25 = (t392 * t46 * t467 + t47 * t729) * t675;
t694 = t25 * mrSges(4,2);
t475 = t130 * t560;
t26 = (t46 * t729 + t47 * t475 / 0.4e1) * t675;
t693 = t26 * mrSges(4,1);
t692 = t47 * mrSges(4,2);
t67 = pkin(3) * t71 + t125;
t48 = (t393 * t458 * t703 + t67 * t629 / 0.2e1) * t386;
t691 = t48 * mrSges(5,2);
t522 = t677 / 0.2e1;
t462 = t386 * t522;
t631 = t386 * t402;
t49 = t392 * t67 * t462 + t393 * t467 * t631;
t690 = t49 * mrSges(5,1);
t65 = (t102 * t714 + t104 * t716) * t560;
t688 = t65 * mrSges(7,2);
t66 = (t104 * t714 + t128 * t719) * t560;
t687 = t66 * mrSges(7,1);
t685 = -(t204 * t426 + (t88 * t461 + t118 * t707 + t150 * t505 + (t120 * t727 + t61 * t715) * t112) * t176) * t632 + t77;
t684 = (t204 * t427 + (t120 * t706 + t150 * t506 + t88 * t460 + (t61 * t713 + t717 * t87) * t112) * t176) * t632 + t76;
t23 = t28 * pkin(4) + t462 * t70;
t538 = t70 * t631;
t459 = t538 / 0.2e1;
t19 = t53 * t23 + t459 * t52;
t681 = mrSges(11,1) * t19;
t20 = -t52 * t23 + t459 * t53;
t680 = mrSges(11,2) * t20;
t660 = t198 / pkin(1);
t548 = pkin(5) * t659;
t418 = t145 ^ 2;
t215 = -pkin(1) * t226 + pkin(5);
t208 = t216 + t366;
t148 = -pkin(1) * t225 * t208 + t215 * t403;
t146 = pkin(1) * t657 + t208 * t215;
t144 = t147 ^ 2;
t143 = 0.1e1 / t146 ^ 2;
t142 = 0.1e1 / t418;
t86 = (-((t215 * t500 + t398 * pkin(5) * t705 + (t208 * t226 - t657) * pkin(1)) * t660 / 0.2e1 + t148 * t548) / t146 - (-t146 * t548 - (-t656 + (-0.2e1 * t215 * pkin(5) - t208 + t500) * t225) * pkin(1) * t660 / 0.2e1) * t148 * t143) / (t143 * t148 ^ 2 + 0.1e1) * t202 * t759;
t83 = 0.1e1 + (t124 * t738 * t661 + (-t123 * t142 * t661 + (-t145 * t142 + t738) * t199 * t492) * t147) * t202 / (t142 * t144 + 0.1e1) * pkin(5);
t69 = t70 ^ 2;
t24 = (t59 * t679 + t60 * t739) * t74 + t70;
t22 = pkin(4) * t27 + t49;
t18 = -t22 * t52 + t48 * t53;
t17 = t22 * t53 + t48 * t52;
t11 = (-(t231 * t685 + t232 * t684) * t50 - (-t231 * t684 + t232 * t685) * t689) * t38 + t27;
t4 = pkin(6) * t5 + t26;
t3 = (t29 * t524 + t41 * t507) * t539 + t5;
t2 = t379 * t4 * t496 + t102 * t25;
t1 = t379 * t25 * t628 / 0.2e1 + t102 * t4;
t13 = [(t17 ^ 2 + t18 ^ 2) * m(11) + (t1 ^ 2 + t2 ^ 2) * m(10) + Ifges(2,3) + t86 ^ 2 * Ifges(6,3) + (m(3) * (t129 / 0.2e1 + t417 / 0.2e1) * t393 * t177 + m(9) * (t144 / 0.2e1 + t418 / 0.2e1) / pkin(5) ^ 2 * t199) * t338 + (mrSges(3,1) * t471 - mrSges(3,2) * t475 + t749) * t71 + m(7) * (t65 ^ 2 + t66 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t25 ^ 2 + t26 ^ 2) + (0.2e1 * t693 - 0.2e1 * t694 + t753) * t5 + (0.2e1 * mrSges(10,1) * t1 - 0.2e1 * mrSges(10,2) * t2 + Ifges(10,3) * t3) * t3 + (0.2e1 * t690 - 0.2e1 * t691 + t751) * t27 + (0.2e1 * t687 - 0.2e1 * t688 + t750) * t21 + (0.2e1 * mrSges(11,1) * t17 - 0.2e1 * mrSges(11,2) * t18 + Ifges(11,3) * t11) * t11 + ((-mrSges(9,1) * t145 + mrSges(9,2) * t147) * pkin(1) * t661 + Ifges(9,3) * t83) * t83; (m(10) * t7 - mrSges(10,2) * t6) * t2 + (m(11) * t20 - mrSges(11,2) * t12) * t18 + (-t724 + t725 + t752) * t3 + (t690 - t691 + t751) * t28 + (t687 - t688 + t750) * t24 + (-t680 + t681 + t748) * t11 + (t693 - t694 + t753) * t10 + (m(11) * t19 + mrSges(11,1) * t12) * t17 + (m(10) * t8 + mrSges(10,1) * t6) * t1 + (t749 + (mrSges(3,1) * t716 + mrSges(3,2) * t714) * t560 + (m(5) * (t402 * t48 + t49 * t677) / 0.2e1 + (mrSges(5,1) * t522 - mrSges(5,2) * t402 / 0.2e1) * t27) * t386 + (m(4) * (-t25 * t47 - t26 * t46) / 0.2e1 + (-t46 * mrSges(4,1) / 0.2e1 + t692 / 0.2e1) * t5) * t580) * t70; t24 ^ 2 * Ifges(7,3) + t69 * Ifges(3,3) + Ifges(8,3) + (t19 ^ 2 + t20 ^ 2) * m(11) + (t7 ^ 2 + t8 ^ 2) * m(10) + (m(5) * (t141 / 0.2e1 + t424 * t113 / 0.2e1) / pkin(4) ^ 2 + m(4) * (t45 / 0.2e1 + t423 / 0.2e1) * t93 * t394 * t380) * t69 / 0.2e1 + (-0.2e1 * t724 + 0.2e1 * t725 + t752) * t6 + (mrSges(5,1) * t386 * t677 * t70 - mrSges(5,2) * t538 + Ifges(5,3) * t28) * t28 + (-0.2e1 * t680 + 0.2e1 * t681 + t748) * t12 + (-mrSges(4,1) * t450 + Ifges(4,3) * t10 + t483 * t692) * t10;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t13(1), t13(2); t13(2), t13(3);];
Mq = res;
