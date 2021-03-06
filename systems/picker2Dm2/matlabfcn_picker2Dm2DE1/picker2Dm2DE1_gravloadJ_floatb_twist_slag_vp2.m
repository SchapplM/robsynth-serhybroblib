% Calculate Gravitation load on the joints for
% picker2Dm2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 18:54
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2DE1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 14:32:11
% EndTime: 2020-05-09 14:39:11
% DurationCPUTime: 78.09s
% Computational Cost: add. (1596568->778), mult. (4827108->1263), div. (26008->22), fcn. (841754->28), ass. (0->573)
t304 = sin(pkin(9));
t307 = cos(pkin(9));
t301 = cos(pkin(8));
t303 = sin(qJ(1));
t306 = cos(qJ(1));
t647 = sin(pkin(8));
t198 = t301 * t303 - t306 * t647;
t200 = t301 * t306 + t303 * t647;
t273 = t306 ^ 2;
t351 = pkin(4) ^ 2;
t353 = pkin(7) ^ 2;
t592 = (t353 - t351);
t346 = (pkin(3) ^ 2);
t670 = 2 * t346;
t531 = t670 + t592;
t340 = (pkin(1) ^ 2);
t673 = 3 * t340;
t490 = t673 + t531;
t305 = cos(qJ(2));
t623 = t303 * t305;
t566 = pkin(1) * t623;
t656 = sin(qJ(2));
t567 = t656 * pkin(7);
t419 = (-0.4e1 * t566 + 0.2e1 * t567) * pkin(3) + t490;
t552 = pkin(3) * t656;
t492 = t552 + pkin(7);
t412 = t492 * t419;
t514 = pkin(7) * t552;
t464 = t353 + 0.2e1 * t514;
t361 = t656 ^ 2;
t615 = t346 * t361;
t569 = 0.2e1 * t615;
t432 = -t346 + t464 + t569;
t428 = pkin(1) * t432;
t423 = -0.2e1 * t428;
t562 = pkin(3) * t623;
t685 = 0.1e1 / pkin(3);
t641 = t685 / 0.2e1;
t718 = 0.1e1 / t641;
t444 = pkin(7) * (-pkin(1) + t562) * t718;
t456 = pkin(3) * t490;
t446 = t305 * t456;
t533 = t340 + t592;
t669 = 4 * t346;
t338 = t340 ^ 2;
t413 = pkin(1) * (0.2e1 * (-t566 + t567) * pkin(3) + t533);
t411 = t492 * t413;
t448 = t340 + t464;
t440 = -t351 + t448;
t667 = pkin(1) * pkin(3);
t431 = t440 * t667;
t422 = 0.4e1 * t305 * t431;
t425 = t340 * t432;
t668 = pkin(3) * pkin(7);
t454 = -0.4e1 * t533 * t668;
t597 = t340 - t353;
t462 = t597 * t669;
t376 = -0.4e1 * t273 * t425 - 0.4e1 * t306 * t411 + t303 * t422 + t361 * t462 + t656 * t454 - t338 - (2 * t531 * t340) - (t353 - (t718 + pkin(4)) * pkin(4)) * (t353 + (t718 - pkin(4)) * pkin(4));
t342 = sqrt(t376);
t621 = t305 * t306;
t561 = pkin(3) * t621;
t429 = t303 * t492 + t561;
t698 = t429 * t342;
t403 = t273 * t423 + t303 * t446 - t306 * t412 + t698 + t656 * t444 + (t569 - t669 - t533) * pkin(1);
t452 = t492 * t306;
t430 = t452 - t562;
t589 = 0.2e1 * pkin(1);
t162 = t430 * t589 + t346 + t448;
t681 = 0.1e1 / t162;
t396 = t681 * t403;
t388 = t396 * t641;
t421 = t306 * t428;
t436 = t670 + t440;
t409 = t436 * t492 + 0.2e1 * t421;
t453 = t273 * t492;
t683 = 0.4e1 * pkin(1);
t684 = -0.2e1 * pkin(1);
t410 = pkin(3) * (t306 * t436 + t453 * t683 + t492 * t684);
t424 = pkin(1) + t430;
t408 = t303 * t409 + t305 * t410 + t342 * t424;
t407 = t685 * t408;
t404 = t407 / 0.2e1;
t399 = t681 * t404;
t107 = t198 * t388 + t200 * t399;
t405 = -t407 / 0.2e1;
t700 = t198 * t681;
t108 = t200 * t388 + t405 * t700;
t97 = t107 * t304 + t108 * t307;
t682 = 0.1e1 / t97;
t678 = 0.1e1 / t97 ^ 2;
t96 = -t107 * t307 + t108 * t304;
t688 = 0.1e1 / (t678 * t96 ^ 2 + 0.1e1);
t703 = t678 * t96;
t584 = t305 * t718;
t244 = pkin(7) * t584;
t488 = t346 * t656 * t305;
t471 = 0.4e1 * t488;
t208 = t471 + t244;
t260 = pkin(1) * t303;
t447 = pkin(1) * t452;
t527 = t656 * t303;
t616 = t340 * t273;
t574 = -0.4e1 * t616;
t477 = t527 * t667;
t609 = 0.2e1 * t477 + t244;
t699 = t346 * t305 ^ 2;
t127 = 0.8e1 * pkin(7) * t260 * t699 + t208 * t574 - 0.4e1 * t413 * t561 - 0.4e1 * t431 * t527 - 0.4e1 * t447 * t609 + (0.2e1 * t462 * t656 + t454) * t305;
t205 = t527 + t621;
t259 = pkin(3) * t305;
t129 = 0.1e1 / t342;
t659 = t129 / 0.2e1;
t418 = t424 * t659;
t511 = 0.4e1 * t552;
t261 = pkin(1) * t306;
t585 = 0.2e1 * t261;
t660 = pkin(7) * t306;
t591 = 0.2e1 * t660;
t654 = pkin(3) * t342;
t112 = t205 * t654 + t127 * t418 + (t208 * t585 + (t340 + (t511 + 0.2e1 * pkin(7)) * pkin(7) + t531) * t259) * t303 - t656 * t410 + (t591 + (0.4e1 * t273 - 0.2e1) * pkin(1)) * t699;
t526 = t656 * t306;
t206 = -t526 + t623;
t420 = t429 * t659;
t545 = t303 * t615;
t113 = t206 * t654 + t127 * t420 + t208 * t273 * t684 - t419 * t561 - (0.4e1 * t477 + t244) * t452 + pkin(1) * t471 - 0.2e1 * pkin(7) * t545 + t305 * t444 - t456 * t527;
t642 = -t685 / 0.2e1;
t394 = t403 * t642;
t400 = t200 * t405;
t551 = t681 * t641;
t496 = t200 * t551;
t497 = t198 * t551;
t159 = 0.1e1 / t162 ^ 2;
t512 = pkin(1) * t561;
t179 = 0.2e1 * t512 + t609;
t636 = t159 * t179;
t84 = t113 * t497 + t112 * t496 + (t198 * t394 + t400) * t636;
t398 = t198 * t404;
t498 = t642 * t700;
t85 = t113 * t496 + t112 * t498 + (t200 * t394 + t398) * t636;
t732 = ((t304 * t84 + t307 * t85) * t703 - (t304 * t85 - t307 * t84) * t682) * t688;
t439 = t453 * t259;
t672 = 8 * t340;
t677 = 0.4e1 * t303;
t132 = t439 * t672 + t411 * t677 + (0.8e1 * t303 * t425 + t422) * t306;
t270 = t303 ^ 2;
t114 = -t698 + t132 * t418 + t270 * t423 + t409 * t306 + (-0.8e1 * t447 - t436) * t562;
t115 = t132 * t420 + t303 * t412 + t306 * t446 + t342 * t430 + t421 * t677 + t439 * t683 + t488 * t591;
t182 = t429 * t589;
t635 = t159 * t182;
t383 = t394 * t635;
t82 = t114 * t496 + t115 * t497 - t198 * t383 - t400 * t635 + t108;
t83 = t114 * t498 + t115 * t496 - t200 * t383 - t398 * t635 - t107;
t731 = ((t304 * t82 + t307 * t83) * t703 - (t304 * t83 - t307 * t82) * t682) * t688;
t392 = t396 / 0.4e1;
t406 = t681 * t408;
t333 = 0.1e1 / pkin(4);
t617 = t333 / t346;
t513 = pkin(1) * t562;
t607 = t340 / 0.3e1 + t353;
t185 = -0.4e1 / 0.9e1 * t513 + 0.4e1 / 0.9e1 * t346 - t351 / 0.9e1 + t607;
t283 = -t351 / 0.6e1;
t292 = 0.2e1 / 0.3e1 * t346;
t463 = t353 - t513;
t192 = t283 + t292 + t463;
t291 = 0.4e1 / 0.3e1 * t346;
t265 = t340 + t353;
t285 = -t351 / 0.3e1;
t537 = t285 + t265;
t230 = t291 + t537;
t257 = -t340 / 0.3e1 + t353;
t515 = -0.2e1 * t562;
t572 = 0.6e1 * t616;
t272 = t306 * t273;
t343 = pkin(1) * t340;
t625 = t272 * t343;
t582 = pkin(7) * t625;
t590 = 0.4e1 * t660;
t618 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t133 = 0.4e1 * t582 + t185 * t572 + t230 * t618 + (t192 * t590 + t257 * t515) * pkin(1);
t324 = 6 * t340;
t287 = -0.2e1 / 0.3e1 * t351;
t331 = 0.2e1 * t353;
t536 = t287 + t292 + t331;
t348 = t346 ^ 2;
t535 = t287 + t265;
t610 = (t292 + t535) * t265 + t348;
t145 = -0.4e1 * t230 * t513 + (t324 + t536) * t346 + t610;
t295 = -t346 / 0.3e1;
t255 = t295 + t353;
t493 = -0.2e1 * t513;
t204 = t255 * t493;
t619 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t164 = t230 * t619 + t204;
t274 = 0.10e2 / 0.3e1 * t340;
t165 = (t274 + t536) * t346 + t610;
t264 = -(3 * t346) + t353;
t523 = 0.8e1 * t582;
t217 = t264 * t523;
t238 = t346 + t533;
t243 = t261 + pkin(7);
t258 = t265 ^ 2;
t308 = 15 * t338;
t315 = 0.18e2 * t353;
t316 = -0.2e1 * t351;
t318 = -0.6e1 * t351;
t341 = t353 ^ 2;
t326 = 0.3e1 * t341;
t347 = pkin(3) * t346;
t334 = t347 ^ 2;
t188 = t230 + t493;
t247 = 0.4e1 * (-t346 + t353) * t340;
t587 = pkin(7) * t261;
t135 = 0.4e1 * t188 * t587 + t247 * t273 + t145;
t246 = pkin(7) * t585;
t267 = -(3 * t340) + t353;
t573 = 0.4e1 * t616;
t203 = t246 + t573 + t267;
t510 = 0.6e1 * t552;
t269 = t656 * t361;
t628 = t269 * t347;
t578 = 0.8e1 * t628;
t437 = t135 * t510 + t203 * t578;
t555 = 0.12e2 * t615;
t557 = 0.12e2 * t616;
t588 = 0.6e1 * pkin(1);
t330 = 0.3e1 * t353;
t604 = (15 * t340) + t330;
t118 = t133 * t555 + t217 + t164 * t557 + t334 + (t346 - t351 + t604) * t348 + (t308 + (t315 + t318 + (6 * t346)) * t340 + t326 + (t316 + t670) * t353) * t346 + t258 * t238 + t437 * t243 + (t145 * t660 - t165 * t562) * t588;
t563 = t340 * t259;
t450 = -t303 * t343 + t563;
t212 = 0.2e1 * t450;
t263 = 2 * t340 + t346;
t147 = -t597 * t259 + t212 * t273 + (t244 * t306 + t263 * t303) * pkin(1);
t593 = t353 + t673;
t249 = t346 + t593;
t216 = t249 * t259;
t250 = (3 * t346) + t265;
t632 = t250 * t303;
t180 = -pkin(1) * t632 + t216;
t571 = t347 * t672;
t674 = 4 * t338;
t242 = pkin(3) * t674 + t571;
t544 = t343 * t619;
t181 = t242 * t305 + t544 * t677;
t323 = 5 * t338;
t595 = t341 + t348;
t310 = 10 * t340;
t603 = t310 + t331;
t614 = t353 * t340;
t197 = t346 * t603 + t323 + t595 + 0.6e1 * t614;
t319 = 5 * t348;
t328 = 0.6e1 * t353;
t209 = t319 + (t310 + t328) * t346 + t258;
t548 = t243 * t628;
t508 = -0.8e1 * t548;
t530 = t243 * t656;
t570 = -0.4e1 * t615;
t220 = -t260 + t259;
t507 = t273 * t563;
t144 = -0.2e1 * t507 + t216 + (t220 * t591 - t632) * pkin(1);
t716 = -0.4e1 * t144;
t122 = t147 * t570 + t181 * t273 + (t530 * t716 + (-t197 + t523) * t305) * pkin(3) + (-0.4e1 * t180 * t660 + (t209 + t508) * t303) * pkin(1);
t214 = t261 + t492;
t105 = t118 * t214 + t122 * t342;
t103 = 0.1e1 / t105;
t225 = t353 + t346 / 0.4e1 + t340 / 0.4e1 - t351 / 0.8e1;
t332 = t351 ^ 2;
t605 = 0.4e1 / 0.7e1 * t353 - t351 / 0.7e1;
t613 = t353 * t351;
t142 = -0.32e2 / 0.21e2 * t225 * t513 + 0.5e1 / 0.42e2 * t348 + (0.16e2 / 0.21e2 * t340 + t605) * t346 + t338 / 0.7e1 + t605 * t340 + t341 - 0.3e1 / 0.7e1 * t613 + t332 / 0.42e2;
t284 = -t351 / 0.4e1;
t696 = t284 + t346 / 0.2e1;
t227 = t607 + t696;
t297 = 0.4e1 / 0.3e1 * t340;
t143 = -0.8e1 / 0.3e1 * t227 * t513 + 0.5e1 / 0.18e2 * t348 + (t297 + t285) * t346 + t341 - t338 / 0.3e1 + t332 / 0.18e2 + (t291 + 0.2e1 / 0.3e1 * t340 + t287) * t353;
t598 = t338 + t341;
t600 = t331 - t351;
t443 = t600 * t340 + t332 / 0.6e1 + t598 - t613;
t189 = -t348 / 0.6e1 + t443;
t299 = t340 / 0.2e1;
t606 = t299 + t353;
t191 = -0.2e1 / 0.3e1 * t513 + t284 + t606;
t329 = 0.4e1 * t353;
t252 = (t329 + t351) * t340;
t296 = -0.2e1 / 0.3e1 * t346;
t256 = t296 + t353;
t286 = -t351 / 0.2e1;
t235 = t286 + t346 + t265;
t714 = -0.4e1 * t235;
t480 = t562 * t714;
t518 = 0.16e2 * t582;
t271 = t273 ^ 2;
t627 = t271 * t338;
t577 = 0.8e1 * t627;
t123 = t256 * t577 + t191 * t518 + 0.14e2 * t142 * t616 - t597 * t348 + (t252 - 0.10e2 / 0.3e1 * t338 + 0.2e1 * t341 - t613) * t346 + t189 * t618 + (0.6e1 * t143 * t660 + t257 * t480) * pkin(1);
t438 = 0.5e1 / 0.6e1 * t348 + t443;
t172 = (t274 + t600) * t346 + t438;
t288 = -0.3e1 / 0.2e1 * t351;
t608 = t332 / 0.2e1 - t348 / 0.2e1;
t491 = -0.3e1 * t613 + t326 + t608;
t612 = t265 * ((t288 + t331) * t340 - 0.3e1 / 0.2e1 * t613 + t598 + t608) + t334;
t131 = -0.6e1 * t172 * t513 + (t308 + (t315 - 0.9e1 * t351) * t340 + t491) * t346 + (t288 + t604) * t348 + t612;
t465 = pkin(1) * t480;
t596 = t341 - t338;
t136 = t255 * t465 - t334 + (-t274 - t592) * t348 + (t252 + t348 / 0.6e1 - t332 / 0.6e1 + t596) * t346 + t189 * t353;
t317 = -0.5e1 * t351;
t322 = 7 * t338;
t140 = (t288 + t330 + (7 * t340)) * t348 + (t322 + (t317 + 0.10e2 * t353) * t340 + t491) * t346 + t612;
t354 = pkin(7) * t353;
t248 = -0.12e2 * pkin(7) * t343 + t354 * t683;
t675 = -8 * t338;
t254 = t675 + 0.12e2 * t614;
t150 = t248 * t306 + t254 * t273 + t518 + t577 + t598 - 0.6e1 * t614;
t168 = t235 * t264 + t493 * t619;
t671 = -6 * t346;
t218 = 0.16e2 * (t353 * t671 + t595) * t338;
t262 = -0.30e2 * t351 + 0.60e2 * t353;
t268 = t361 ^ 2;
t148 = t465 + (t324 + t600) * t346 + t438;
t166 = t235 * t619 + t204;
t124 = 0.6e1 * t148 * t587 + t166 * t557 + t131 + t217;
t327 = 0.8e1 * t353;
t505 = t343 * t562;
t183 = -0.4e1 * t505 + t674 + (t669 + t316 + t327) * t340;
t187 = -t340 + t463 + t696;
t134 = t523 + t183 * t273 + t235 * t267 + (t187 * t590 + t515 * t618) * pkin(1);
t509 = 0.8e1 * t552;
t560 = 0.32e2 * t628;
t435 = t124 * t509 + t134 * t560;
t599 = t332 - t348;
t451 = 0.6e1 * t341 + t599 - 0.6e1 * t613;
t532 = -t346 + t592;
t631 = t258 * (t340 + t532);
t101 = t218 * t271 + 0.32e2 * t168 * t582 + 0.24e2 * t136 * t616 + (t316 + t329 + (28 * t340)) * t334 + t238 * t631 + (0.24e2 * t123 * t361 + t262 * t338 + t318 * t341 + t451 * t324 + t599 * t331 + 0.28e2 * t343 ^ 2 + 0.4e1 * t354 ^ 2) * t346 + t435 * t243 + 0.8e1 * (t131 * t660 - t140 * t562) * pkin(1) + (0.16e2 * t150 * t268 + t262 * t340 + (70 * t338) + t348 + t451) * t348;
t223 = 0.7e1 / 0.6e1 * t346 + t283 + t606;
t293 = t346 / 0.3e1;
t538 = t283 + t293 + t353;
t226 = t297 + t538;
t163 = -t223 * t260 + t226 * t259;
t539 = t351 / 0.3e1 + t293 + t331;
t169 = -t346 * t597 - 0.5e1 / 0.3e1 * t338 + t539 * t340 + t353 * (t285 + t255);
t232 = t340 + t538;
t251 = t670 - t597;
t174 = t232 * t259 - t251 * t260 / 0.2e1;
t626 = t272 * t338;
t275 = -0.20e2 / 0.3e1 * t340;
t540 = 0.2e1 / 0.3e1 * t351 + t292 + t329;
t541 = 0.4e1 / 0.3e1 * t351 + t291 - 0.2e1 * t353;
t658 = t348 / 0.2e1 - (t275 + t540) * t346 / 0.2e1 + 0.3e1 / 0.2e1 * t338 - t541 * t340 / 0.2e1 - t341 / 0.2e1;
t661 = pkin(7) * t303;
t125 = -0.4e1 * t626 * t661 + t163 * t573 + (-0.8e1 / 0.3e1 * t627 + t169) * t259 + (t174 * t590 + t303 * t658) * pkin(1);
t602 = t316 - (2 * t346);
t534 = t328 + t602;
t156 = t348 + (t287 + t296 + t603) * t346 + t323 + t534 * t340 + t353 * (t287 + t256);
t149 = t156 * t259;
t173 = t319 + (t310 + t534) * t346 + (t296 + t535) * t265;
t634 = t173 * t303;
t137 = -pkin(1) * t634 + t149;
t171 = -(3 * t348) + (t275 + t541) * t346 + t540 * t340 + t596;
t175 = -0.5e1 / 0.3e1 * t348 + (-t340 + t539) * t346 + t353 * (t295 + t537);
t586 = -0.2e1 * t260;
t138 = t171 * t259 + t175 * t586;
t695 = t330 - t346 - t351;
t237 = t695 * t310;
t601 = t317 - (5 * t346);
t139 = t334 + ((21 * t340) + t695) * t348 + (t353 * t602 + t237 + t326 + (35 * t338)) * t346 + (t322 + (t327 + t601) * t340 + t353 * t532) * t265;
t213 = 0.4e1 * t450;
t219 = 0.2e1 * t260 + t259;
t234 = t286 + t249;
t141 = t267 * t259 + t213 * t273 + (t219 * t660 + t234 * t303) * t589;
t224 = t353 + 0.5e1 / 0.2e1 * t346 + 0.3e1 / 0.2e1 * t340 + t286;
t630 = t264 * t303;
t176 = t224 * t259 + pkin(1) * t630 / 0.2e1;
t556 = -0.12e2 * t615;
t193 = 0.4e1 / 0.3e1 * t616 + t246 + t257;
t629 = t268 * t348;
t697 = 0.7e1 * t334 + ((35 * t340) + 0.15e2 * t353 + t601) * t348 + ((21 * t338) + t237 + 0.9e1 * t341 + (t318 + t671) * t353) * t346 + t631 - 0.24e2 * t193 * t629;
t229 = 0.8e1 / 0.3e1 * t346 + t537;
t231 = t285 + t292 + t593;
t167 = -t229 * t260 + t231 * t259;
t236 = 0.5e1 / 0.6e1 * t346 + t299 + t283;
t546 = t303 * t619;
t178 = pkin(1) * t546 + t236 * t584;
t564 = t272 * t259;
t506 = t343 * t564;
t126 = -0.8e1 * pkin(7) * t506 + t178 * t574 + t149 + (t167 * t590 - t634) * pkin(1);
t717 = -0.6e1 * t126;
t106 = t176 * t518 + t141 * t508 + t125 * t556 - 0.6e1 * t138 * t616 + (t530 * t717 + (0.24e2 * t255 * t627 - t139) * t305) * pkin(3) + (-0.6e1 * t137 * t660 + t697 * t303) * pkin(1);
t88 = t101 * t214 + t106 * t342;
t644 = t103 * t88;
t69 = (t342 * t392 + t406 * t644 / 0.4e1) * t617;
t67 = t206 * t69;
t387 = t103 * t392;
t657 = -t342 / 0.4e1;
t68 = (t387 * t88 + t406 * t657) * t617;
t55 = -t205 * t68 + t67;
t66 = t205 * t69;
t56 = -t206 * t68 - t66;
t34 = atan2(t55, t56);
t31 = sin(t34);
t32 = cos(t34);
t719 = t105 ^ 2;
t645 = t719 / t88 ^ 2;
t73 = 0.1e1 / (t376 * t645 + 0.1e1);
t442 = pkin(4) * t333 * t342 * t73 / t642 * t645;
t466 = t105 * t73 / t88 * t659;
t104 = 0.1e1 / t719;
t643 = t104 * t88;
t504 = t642 * t643;
t550 = t103 * t641;
t481 = pkin(7) * t507;
t222 = -0.4e1 * t481;
t245 = pkin(7) * t586;
t519 = 0.32e2 / 0.3e1 * t338;
t482 = t272 * t519;
t483 = 0.64e2 / 0.3e1 * t225 * t343;
t494 = t272 * t544;
t624 = t273 * t343;
t580 = pkin(7) * t624;
t517 = -0.48e2 * t580;
t522 = pkin(7) * t573;
t524 = t132 * t659;
t620 = t306 * t340;
t542 = t303 * t620;
t554 = 0.24e2 * t615;
t558 = -0.24e2 * t624;
t559 = -0.32e2 * t626;
t622 = t303 * t306;
t575 = -0.2e1 * t622;
t576 = -0.4e1 * t625;
t579 = -0.8e1 * t628;
t581 = pkin(7) * t616;
t583 = pkin(7) * t627;
t516 = pkin(7) * t558;
t713 = -0.24e2 * t255;
t611 = t506 * t713 + t516 * t630;
t662 = pkin(7) * t273;
t676 = -0.2e1 * t306;
t59 = ((0.12e2 * t270 * t338 * t662 + t223 * t576 - 0.2e1 * t251 * t581 - 0.4e1 * t583) * t556 + 0.12e2 * t175 * t625 + (t556 * t658 + t697) * t261 + (t173 * t572 + t264 * t577) * pkin(7)) * t342 + t106 * t524 + ((-pkin(7) * t271 * t519 + t257 * t261 * t714 - 0.16e2 * t227 * t581 - t272 * t483) * t554 - 0.64e2 * t583 * t619 - 0.96e2 * t255 * t235 * t625 - 0.48e2 * t172 * t581 - 0.8e1 * t140 * t261) * t214 * t259 + (0.12e2 * (-(-0.8e1 * t163 * t620 + t259 * t482) * t615 + t255 * t564 * t675 - 0.4e1 * t176 * t580 + t138 * t620) * t342 + (0.16e2 * (t254 * t676 - t248 + t517 + t559) * t629 + (-0.28e2 * t142 * t620 + t191 * t517 + t256 * t559) * t554 - 0.4e1 * t218 * t272 - 0.96e2 * t168 * t580 - 0.48e2 * t136 * t620) * t214 + (-t435 * t214 - t101 + (t126 * t510 + t141 * t578 + (-0.24e2 * t245 + 0.64e2 * t542) * t629) * t342 + ((0.48e2 * t174 * t615 + 0.6e1 * t137) * t342 + (-0.144e3 * t143 * t615 - 0.8e1 * t131) * t214) * pkin(7)) * pkin(1)) * t303 + (((t522 + t213 * t575 + t576 + (-t219 * t661 + t234 * t306) * t589) * t579 - 0.6e1 * (-0.4e1 * t494 + (t178 * t303 * t672 - pkin(1) * t173) * t306 + (-0.4e1 * t167 * t260 + (-0.4e1 * t229 * t340 + 0.24e2 * t505) * t273) * pkin(7)) * t552) * t342 + ((t222 + (-0.2e1 * t261 * t618 + t576) * t259 + (t183 * t676 + (-0.4e1 * pkin(1) * t187 + t558) * pkin(7)) * t303) * t560 + (0.24e2 * (-t235 * t259 * t662 - t166 * t622) * t340 + (-t148 * t661 - t172 * t561) * t588 + t611) * t509) * t214) * t243;
t568 = 0.8e1 * t615;
t89 = (t270 * t269 * t571 + (t261 * t263 - 0.2e1 * t625) * t570 + 0.4e1 * t494 + t250 * t522 + t209 * t261) * t342 + t122 * t524 + ((-0.8e1 / 0.3e1 * t506 + t222 - 0.2e1 * t257 * t512) * t555 - 0.24e2 * t230 * t481 - 0.6e1 * t165 * t512 + t611) * t214 + ((t212 * t306 * t568 + t181 * t676 + t516 * t259) * t342 + (0.12e2 * (-t185 * t620 - t580) * t555 - 0.24e2 * t164 * t620) * t214 + (t144 * t342 * t511 - t437 * t214 - t118 + ((t259 * t568 + 0.4e1 * t180) * t342 + (-0.48e2 * t192 * t615 - 0.6e1 * t145) * t214) * pkin(7)) * pkin(1)) * t303 + ((t579 * t261 - 0.4e1 * ((t561 * t677 - 0.2e1 * t662) * t340 + (-0.2e1 * t220 * t661 - t250 * t306) * pkin(1)) * t552) * t342 + ((t245 - 0.8e1 * t542) * t578 + (-0.8e1 * t481 + t247 * t575 + (-t188 * t661 - t230 * t561) * t683) * t510) * t214) * t243;
t20 = t132 * t466 + (t504 * t89 + t550 * t59) * t442;
t374 = atan2(t399, t388);
t372 = sin(t374);
t373 = cos(t374);
t100 = t303 * t372 - t306 * t373;
t401 = 0.1e1 / t403 ^ 2;
t375 = 0.1e1 / (t401 * t408 ^ 2 + 0.1e1);
t370 = t375 * t401 * t408;
t371 = t375 / t403;
t368 = t114 * t371 - t115 * t370;
t364 = t368 * t372;
t365 = t368 * t373;
t41 = -t303 * t364 + t306 * t365 - t100;
t99 = -t303 * t373 - t306 * t372;
t42 = t303 * t365 + t306 * t364 - t99;
t549 = t333 * t641;
t72 = atan2(t342 * t549, t549 * t644);
t70 = sin(t72);
t71 = cos(t72);
t46 = t100 * t70 + t71 * t99;
t427 = -t20 * t46 + t41 * t70 + t42 * t71;
t702 = t100 * t71 - t70 * t99;
t6 = t20 * t702 - t41 * t71 + t42 * t70;
t730 = t31 * t6 - t32 * t427;
t528 = t343 * t656;
t486 = t273 * t528;
t655 = pkin(3) * t303;
t441 = t255 * t486 * t655;
t190 = 0.24e2 * t441;
t489 = t340 * t526;
t472 = pkin(7) * t489;
t445 = t472 * t655;
t211 = 0.4e1 * t445;
t501 = pkin(1) * t526;
t461 = t501 * t668;
t221 = -0.2e1 * t461;
t479 = t272 * t514;
t455 = t343 * t479;
t633 = t243 * t305;
t495 = t347 * t361 * t633;
t469 = -0.24e2 * t495;
t529 = t273 * t656;
t487 = t340 * t529;
t470 = -0.4e1 * t487;
t478 = -0.4e1 * t501;
t485 = t656 * t627;
t499 = t346 * t530;
t500 = pkin(3) * t530;
t553 = pkin(1) * t656;
t502 = t257 * t553;
t525 = t127 * t659;
t547 = t269 * t305 * t348;
t565 = pkin(3) * t633;
t715 = 0.6e1 * t172;
t62 = (-0.96e2 * t193 * t547 * t260 + t221 * t508 + t141 * t469 - 0.24e2 * t125 * t488 - 0.6e1 * (0.8e1 * t236 * t487 - t656 * t156 + (t231 * t478 + 0.8e1 * t272 * t528) * pkin(7)) * t499 + t565 * t717 - 0.16e2 * t224 * t455 + 0.6e1 * t156 * t461 + t139 * t552 + ((-t267 * t656 + t470) * t508 + (0.8e1 / 0.3e1 * t485 + t226 * t470 + t232 * pkin(7) * t478 - t169 * t656) * t556 + t485 * t713 + 0.6e1 * t171 * t487) * pkin(3)) * t342 + t106 * t525 + 0.8e1 * (0.8e1 * t150 * t547 + 0.4e1 * (t211 + (0.2e1 * t553 * t618 + 0.4e1 * t486) * t655) * t548 + 0.12e2 * t134 * t495 + 0.3e1 * (t483 * t529 + 0.4e1 * t235 * t502 + (0.16e2 * t227 * t489 + t482 * t656) * pkin(7)) * pkin(3) * t545 + 0.6e1 * t123 * t488 + (t190 + (0.24e2 * t235 * t472 + t553 * t715) * t655) * t500 + t124 * t565 + 0.8e1 * t338 * t479 * t546 + 0.12e2 * t235 * t441 + t445 * t715 + t140 * t477) * t214 + t101 * t259;
t90 = (t469 * t260 + (t221 + (t597 * t656 - 0.2e1 * t487) * pkin(3)) * t570 - 0.8e1 * t147 * t488 - 0.4e1 * (t221 + (-t249 * t656 + 0.2e1 * t487) * pkin(3)) * t500 + t565 * t716 - 0.8e1 * t455 - t242 * t529 + 0.4e1 * t249 * t461 + t197 * t552) * t342 + t122 * t525 + (0.24e2 * t203 * t495 + (t211 + (0.8e1 / 0.3e1 * t486 + 0.2e1 * t502) * t655) * t555 + 0.24e2 * t133 * t488 + 0.6e1 * (0.4e1 * t230 * t553 + 0.8e1 * t472) * t303 * t499 + 0.6e1 * t135 * t565 + t190 + 0.24e2 * t230 * t445 + 0.6e1 * t165 * t477) * t214 + t118 * t259;
t21 = t127 * t466 + (t504 * t90 + t550 * t62) * t442;
t369 = t112 * t371 - t113 * t370;
t366 = t369 * t372;
t367 = t369 * t373;
t43 = t303 * t366 - t306 * t367;
t44 = t303 * t367 + t306 * t366;
t426 = -t21 * t46 - t43 * t70 + t44 * t71;
t8 = t21 * t702 + t43 * t71 + t44 * t70;
t729 = t31 * t8 - t32 * t426;
t728 = t31 * t427 + t32 * t6;
t727 = t31 * t426 + t32 * t8;
t712 = -m(4) - m(10);
t726 = pkin(2) * t712 - mrSges(3,1);
t725 = t31 * t46 - t32 * t702;
t724 = t31 * t702 + t32 * t46;
t80 = atan2(t96, t97);
t74 = sin(t80);
t76 = cos(t80);
t467 = t100 * t76 - t74 * t99;
t646 = t100 * t74;
t51 = t76 * t99 + t646;
t81 = atan2(t96, -t97);
t75 = sin(t81);
t77 = cos(t81);
t721 = t467 * t77 - t51 * t75;
t720 = -t467 * t75 - t51 * t77;
t701 = t129 * t681;
t393 = t403 * t617;
t386 = -t393 / 0.4e1;
t382 = t159 * t386;
t402 = t408 * t617;
t397 = t402 / 0.4e1;
t694 = t159 * t342 * t397 + t382 * t644;
t395 = -t88 * t402 / 0.4e1;
t693 = t103 * t159 * t395 + t342 * t382;
t692 = m(3) + m(7) + m(9);
t690 = -t41 * t76 + t42 * t74 - t467 * t731;
t689 = t43 * t76 + t44 * t74 - t467 * t732;
t687 = mrSges(10,1) * t75 + mrSges(10,2) * t77 - mrSges(4,2) - mrSges(7,2);
t686 = -pkin(6) * m(10) + mrSges(10,1) * t77 - mrSges(10,2) * t75 - mrSges(4,1) - mrSges(7,1);
t664 = pkin(3) * t43;
t663 = pkin(3) * t44;
t54 = 0.1e1 / t56 ^ 2;
t653 = t54 * t55;
t380 = t393 * t701 / 0.8e1;
t384 = t104 * t681 * t395;
t391 = t103 * t681 * t397;
t543 = t681 * t617;
t484 = t543 / 0.4e1;
t449 = t484 * t644;
t459 = t342 * t484;
t649 = t114 * t449 + t115 * t459 + t132 * t380 - t693 * t182 + t384 * t89 + t391 * t59 + t68;
t648 = t112 * t449 + t113 * t459 + t127 * t380 + t693 * t179 + t384 * t90 + t391 * t62 - t68;
t640 = t42 * pkin(3) + t260;
t638 = -t41 * pkin(3) - t261;
t460 = t543 * t657;
t434 = -t41 * t74 - t731 * t646 + (-t731 * t99 - t42) * t76;
t433 = t43 * t74 - t732 * t646 + (-t732 * t99 - t44) * t76;
t417 = t721 * mrSges(10,1) + t720 * mrSges(10,2);
t416 = t720 * mrSges(10,1) - t721 * mrSges(10,2);
t415 = t725 * mrSges(11,1) + t724 * mrSges(11,2);
t414 = t724 * mrSges(11,1) - t725 * mrSges(11,2);
t390 = -t402 * t701 / 0.8e1;
t381 = t387 * t617;
t377 = t681 * t386 * t643;
t161 = atan2(-t198, -t200);
t155 = cos(t161);
t153 = sin(t161);
t120 = -t153 * t301 - t155 * t647;
t119 = -t153 * t647 + t155 * t301;
t53 = 0.1e1 / t56;
t33 = 0.1e1 / (t54 * t55 ^ 2 + 0.1e1);
t27 = t112 * t460 + t113 * t449 + t127 * t390 + t694 * t179 + t377 * t90 + t381 * t62;
t23 = t114 * t460 + t115 * t449 + t132 * t390 - t694 * t182 + t377 * t89 + t381 * t59;
t2 = ((-t205 * t27 + t206 * t648 - t66) * t53 - ((-t27 - t69) * t206 - t648 * t205) * t653) * t33;
t1 = ((t649 * t206 + (-t23 + t69) * t205) * t53 - (-t205 * t649 - t206 * t23 + t67) * t653) * t33;
t3 = [(mrSges(2,1) * t306 - mrSges(2,2) * t303 - m(11) * (pkin(4) * t6 + t638) - m(5) * t638 + t728 * mrSges(11,1) - t730 * mrSges(11,2) - t427 * mrSges(5,2) - mrSges(6,1) * t119 - mrSges(6,2) * t120 - t6 * mrSges(5,1) + mrSges(3,1) * t41 - mrSges(3,2) * t42 - t415 * t1 - t417 * t731 + t712 * (-t41 * pkin(2) - t261) - t686 * t690 + t687 * t434 + t692 * t261) * g(2) + (-mrSges(3,1) * t42 - mrSges(3,2) * t41 - mrSges(2,1) * t303 - mrSges(2,2) * t306 - mrSges(6,1) * t120 + mrSges(6,2) * t119 - m(11) * (pkin(4) * t427 + t640) - m(5) * t640 - t730 * mrSges(11,1) - t728 * mrSges(11,2) + t6 * mrSges(5,2) - t427 * mrSges(5,1) - t414 * t1 - t416 * t731 + t712 * (t42 * pkin(2) + t260) + t686 * t434 + t687 * t690 - t692 * t260) * g(1), (-t656 * mrSges(8,1) - t305 * mrSges(8,2) - m(11) * (pkin(4) * t8 + t664) + t727 * mrSges(11,1) - t729 * mrSges(11,2) - t415 * t2 - m(5) * t664 - mrSges(5,1) * t8 - mrSges(5,2) * t426 - mrSges(3,2) * t44 - t417 * t732 - t686 * t689 + t687 * t433 + t726 * t43) * g(2) + (-t305 * mrSges(8,1) + t656 * mrSges(8,2) - t416 * t732 + mrSges(3,2) * t43 - m(5) * t663 - mrSges(5,1) * t426 + mrSges(5,2) * t8 - m(11) * (pkin(4) * t426 + t663) - t729 * mrSges(11,1) - t727 * mrSges(11,2) - t414 * t2 + t686 * t433 + t687 * t689 + t726 * t44) * g(1)];
taug = t3(:);
