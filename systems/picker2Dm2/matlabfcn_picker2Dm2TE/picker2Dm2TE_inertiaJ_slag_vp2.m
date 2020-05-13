% Calculate joint inertia matrix for
% picker2Dm2TE
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
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2TE_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2TE_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2TE_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 10:36:27
% EndTime: 2020-05-09 10:37:21
% DurationCPUTime: 47.41s
% Computational Cost: add. (963995->765), mult. (2955384->1263), div. (13831->23), fcn. (493872->10), ass. (0->549)
t311 = cos(pkin(9));
t619 = sin(pkin(9));
t308 = sin(qJ(1));
t309 = cos(qJ(2));
t568 = t308 * t309;
t237 = pkin(3) * t568;
t426 = pkin(1) * t237;
t224 = -0.2e1 * t426;
t307 = sin(qJ(2));
t641 = 0.2e1 * t307;
t530 = pkin(7) * t641;
t243 = pkin(3) * t530;
t353 = pkin(1) ^ 2;
t355 = pkin(7) ^ 2;
t267 = t353 + t355;
t347 = pkin(3) ^ 2;
t449 = t347 + t267;
t261 = pkin(3) * t307;
t241 = t261 + pkin(7);
t310 = cos(qJ(1));
t221 = t241 * t310;
t491 = pkin(1) * t221;
t148 = t224 + t243 + t449 + 0.2e1 * t491;
t146 = 0.1e1 / t148;
t566 = t309 * t310;
t486 = pkin(3) * t566;
t582 = t241 * t308;
t176 = t486 + t582;
t526 = 0.2e1 * pkin(1);
t169 = t176 * t526;
t305 = cos(pkin(8));
t593 = sin(pkin(8));
t188 = t305 * t308 - t593 * t310;
t190 = t305 * t310 + t593 * t308;
t348 = 0.1e1 / pkin(3);
t147 = 0.1e1 / t148 ^ 2;
t324 = 0.2e1 * t347;
t330 = 0.3e1 * t353;
t341 = pkin(4) ^ 2;
t532 = t355 - t341;
t234 = t324 + t330 + t532;
t409 = -0.4e1 * t426;
t167 = t234 + t243 + t409;
t222 = t237 - pkin(1);
t250 = t353 + t532;
t272 = t307 ^ 2;
t536 = -t347 + t355;
t574 = t272 * t347;
t193 = t243 + t536 + 0.2e1 * t574;
t276 = t310 ^ 2;
t584 = t193 * t276;
t204 = t243 + t250;
t170 = t224 + t204;
t663 = 0.4e1 * t353;
t277 = t347 * t663;
t560 = t347 * t355;
t248 = t277 - 0.4e1 * t560;
t338 = -0.2e1 * t355;
t339 = 0.2e1 * pkin(3);
t351 = t353 ^ 2;
t470 = -0.4e1 * pkin(3) * pkin(7) * t250;
t666 = 0.4e1 * pkin(1);
t118 = t248 * t272 + t307 * t470 - t351 - (t355 - (t339 + pkin(4)) * pkin(4)) * (t355 + (t339 - pkin(4)) * pkin(4)) + (t338 + 0.2e1 * t341 - 0.4e1 * t347 - 0.4e1 * t584) * t353 + (-t170 * t221 + t204 * t237) * t666;
t356 = sqrt(t118);
t586 = t176 * t356;
t106 = -t167 * t221 + t586 + (t222 * t530 + t234 * t568) * pkin(3) + (-0.2e1 * t584 + (0.2e1 * t272 - 0.4e1) * t347 - t250) * pkin(1);
t202 = t324 + t204;
t263 = pkin(1) * t310;
t493 = t193 * t263;
t138 = t202 * t241 + 0.2e1 * t493;
t140 = t202 * t310 + (0.4e1 * t276 - 0.2e1) * t241 * pkin(1);
t549 = -t237 + t221;
t174 = pkin(1) + t549;
t262 = pkin(3) * t309;
t108 = t138 * t308 + t140 * t262 + t174 * t356;
t622 = t190 / 0.2e1;
t623 = t188 / 0.2e1;
t390 = t106 * t623 + t108 * t622;
t380 = t147 * t390;
t567 = t308 * t310;
t583 = t241 * t276;
t662 = 0.8e1 * t353;
t121 = (t193 * t567 + t583 * t262) * t662 + (t170 * t582 + t204 * t486) * t666;
t117 = 0.1e1 / t356;
t628 = t117 / 0.2e1;
t442 = t176 * t628;
t633 = pkin(7) * t310;
t529 = 0.2e1 * t633;
t569 = t307 * t347;
t100 = t549 * t356 + t121 * t442 + (t167 * t241 + 0.4e1 * t493) * t308 + (t529 * t569 + (t234 * t310 + t583 * t666) * pkin(3)) * t309;
t437 = t100 / 0.2e1 - t108 / 0.2e1;
t629 = t106 / 0.2e1;
t273 = t308 ^ 2;
t443 = t174 * t628;
t649 = -0.2e1 * pkin(1);
t98 = -t586 + t121 * t443 + t193 * t273 * t649 + t138 * t310 + (-t202 - 0.8e1 * t491) * t237;
t458 = t629 + t98 / 0.2e1;
t64 = (t169 * t380 + (t188 * t437 + t190 * t458) * t146) * t348;
t624 = -t188 / 0.2e1;
t389 = t106 * t622 + t108 * t624;
t379 = t147 * t389;
t65 = (t169 * t379 + (-t188 * t458 + t190 * t437) * t146) * t348;
t588 = t146 * t348;
t376 = t389 * t588;
t95 = t390 * t588;
t86 = t311 * t376 + t619 * t95;
t375 = t86 ^ 2;
t645 = 0.1e1 / t375;
t84 = t311 * t95 - t619 * t376;
t79 = t84 ^ 2;
t651 = 0.1e1 / (t645 * t79 + 0.1e1);
t668 = t651 / t86;
t671 = (-t311 * t64 + t619 * t65) * t668;
t428 = pkin(1) * t486;
t515 = 0.2e1 * t262;
t244 = pkin(7) * t515;
t571 = t307 * t308;
t488 = pkin(3) * t571;
t427 = pkin(1) * t488;
t550 = 0.2e1 * t427 + t244;
t165 = 0.2e1 * t428 + t550;
t565 = t309 * t347;
t462 = t307 * t565;
t199 = t244 + 0.4e1 * t462;
t557 = t353 * t276;
t500 = -0.4e1 * t557;
t634 = pkin(7) * t308;
t511 = t347 * t634;
t114 = t199 * t500 + (t248 * t641 + t470) * t309 + (-t550 * t221 + 0.2e1 * t309 ^ 2 * t511 + (-t170 * t566 - t204 * t571) * pkin(3)) * t666;
t521 = 0.2e1 * t263;
t617 = pkin(1) * t276;
t650 = 0.2e1 * pkin(7);
t96 = t114 * t443 + t199 * t308 * t521 + ((t308 * t356 - t140) * t307 + (t310 * t356 + (t241 * t650 + t202) * t308 + (-pkin(1) + 0.2e1 * t617 + t633) * t515) * t309) * pkin(3);
t570 = t307 * t310;
t487 = pkin(3) * t570;
t97 = (t237 - t487) * t356 + t114 * t442 - 0.2e1 * t199 * t617 - (t244 + 0.4e1 * t427) * t221 - 0.2e1 * t272 * t511 - t234 * t488 + (t569 * t666 + (-t167 * t310 + t222 * t650) * pkin(3)) * t309;
t67 = (-t165 * t380 + (t96 * t622 + t97 * t623) * t146) * t348;
t68 = (-t165 * t379 + (t97 * t622 + t96 * t624) * t146) * t348;
t670 = (-t311 * t67 + t619 * t68) * t668;
t516 = 0.8e1 * t261;
t271 = t307 * t272;
t360 = pkin(3) * t347;
t576 = t271 * t360;
t482 = 0.32e2 * t576;
t368 = t106 ^ 2;
t105 = 0.1e1 / t368;
t107 = t108 ^ 2;
t423 = pkin(3) * t148 * t348 / (t105 * t107 + 0.1e1);
t391 = -0.2e1 * t105 * t108 * t423;
t396 = t423 / t629;
t626 = -t147 / 0.2e1;
t440 = t169 * t626;
t627 = t146 / 0.2e1;
t55 = 0.1e1 + (-t108 * t440 + t98 * t627) * t396 + (t100 * t627 - t106 * t440) * t391;
t436 = t55 + t671;
t508 = t84 * t651 * t645;
t20 = (t311 * t65 + t619 * t64) * t508 + t436;
t657 = Ifges(4,3) + Ifges(7,3);
t669 = t657 * t20;
t547 = t353 / 0.3e1 + t355;
t173 = -0.4e1 / 0.9e1 * t426 + 0.4e1 / 0.9e1 * t347 - t341 / 0.9e1 + t547;
t287 = -t341 / 0.6e1;
t296 = 0.2e1 / 0.3e1 * t347;
t397 = t355 - t426;
t182 = t287 + t296 + t397;
t295 = 0.4e1 / 0.3e1 * t347;
t289 = -t341 / 0.3e1;
t453 = t289 + t267;
t226 = t295 + t453;
t259 = -t353 / 0.3e1 + t355;
t429 = -0.2e1 * t237;
t497 = 0.6e1 * t557;
t275 = t310 * t276;
t357 = pkin(1) * t353;
t572 = t275 * t357;
t509 = pkin(7) * t572;
t528 = 0.4e1 * t633;
t562 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t122 = 0.4e1 * t509 + t173 * t497 + t226 * t562 + (t182 * t528 + t259 * t429) * pkin(1);
t329 = 0.6e1 * t353;
t291 = -0.2e1 / 0.3e1 * t341;
t337 = 0.2e1 * t355;
t452 = t291 + t296 + t337;
t361 = t347 ^ 2;
t451 = t291 + t267;
t551 = (t296 + t451) * t267 + t361;
t134 = t226 * t409 + (t329 + t452) * t347 + t551;
t299 = -t347 / 0.3e1;
t257 = t299 + t355;
t195 = t257 * t224;
t563 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t150 = t226 * t563 + t195;
t278 = 0.10e2 / 0.3e1 * t353;
t151 = (t278 + t452) * t347 + t551;
t266 = -0.3e1 * t347 + t355;
t438 = 0.8e1 * t509;
t210 = t266 * t438;
t235 = -t341 + t449;
t242 = t263 + pkin(7);
t260 = t267 ^ 2;
t312 = 0.15e2 * t351;
t319 = 0.18e2 * t355;
t320 = -0.2e1 * t341;
t322 = -0.6e1 * t341;
t354 = t355 ^ 2;
t332 = 0.3e1 * t354;
t344 = t360 ^ 2;
t178 = t224 + t226;
t247 = t536 * t663;
t520 = 0.4e1 * t263;
t124 = pkin(7) * t178 * t520 + t247 * t276 + t134;
t246 = pkin(7) * t521;
t269 = -0.3e1 * t353 + t355;
t498 = 0.4e1 * t557;
t194 = t246 + t498 + t269;
t506 = 0.8e1 * t576;
t517 = 0.6e1 * t261;
t385 = t124 * t517 + t194 * t506;
t476 = 0.12e2 * t557;
t480 = 0.12e2 * t574;
t336 = 0.3e1 * t355;
t544 = 0.15e2 * t353 + t336;
t665 = 0.6e1 * pkin(1);
t101 = t122 * t480 + t210 + t150 * t476 + t344 + (-t341 + t347 + t544) * t361 + (t312 + (t319 + t322 + 0.6e1 * t347) * t353 + t332 + (t320 + t324) * t355) * t347 + t260 * t235 + t385 * t242 + (t134 * t633 - t151 * t237) * t665;
t485 = t353 * t262;
t393 = -t308 * t357 + t485;
t205 = 0.2e1 * t393;
t265 = t347 + 0.2e1 * t353;
t268 = -t353 + t355;
t136 = t268 * t262 + t205 * t276 + (t244 * t310 + t265 * t308) * pkin(1);
t540 = t330 + t355;
t251 = t347 + t540;
t209 = t251 * t262;
t252 = 0.3e1 * t347 + t267;
t580 = t252 * t308;
t166 = -pkin(1) * t580 + t209;
t496 = t360 * t662;
t615 = pkin(3) * t351;
t240 = t496 + 0.4e1 * t615;
t461 = t357 * t563;
t168 = t240 * t309 + 0.4e1 * t308 * t461;
t328 = 0.5e1 * t351;
t533 = t354 + t361;
t314 = 0.10e2 * t353;
t543 = t314 + t337;
t555 = t355 * t353;
t187 = t347 * t543 + t328 + t533 + 0.6e1 * t555;
t323 = 0.5e1 * t361;
t334 = 0.6e1 * t355;
t200 = t323 + (t314 + t334) * t347 + t260;
t465 = t242 * t576;
t420 = -0.8e1 * t465;
t505 = -0.4e1 * t574;
t581 = t242 * t307;
t616 = pkin(1) * t308;
t213 = t262 - t616;
t419 = t276 * t485;
t133 = -0.2e1 * t419 + t209 + (t213 * t529 - t580) * pkin(1);
t644 = -0.4e1 * t133;
t109 = t136 * t505 + t168 * t276 + (t581 * t644 + (-t187 + t438) * t309) * pkin(3) + (-0.4e1 * t166 * t633 + (t200 + t420) * t308) * pkin(1);
t207 = t261 + t242;
t92 = t101 * t207 + t109 * t356;
t667 = t92 ^ 2;
t340 = t341 ^ 2;
t535 = t351 + t354;
t539 = t337 - t341;
t556 = t355 * t341;
t386 = t539 * t353 + t340 / 0.6e1 + t535 - t556;
t381 = 0.5e1 / 0.6e1 * t361 + t386;
t158 = (t278 + t539) * t347 + t381;
t664 = -0.6e1 * t158;
t7 = t436 - t671;
t661 = Ifges(10,3) * t7;
t342 = 0.1e1 / pkin(4);
t559 = t348 * t356;
t218 = t355 + t347 / 0.4e1 + t353 / 0.4e1 - t341 / 0.8e1;
t545 = 0.4e1 / 0.7e1 * t355 - t341 / 0.7e1;
t131 = -0.32e2 / 0.21e2 * t218 * t426 + 0.5e1 / 0.42e2 * t361 + (0.16e2 / 0.21e2 * t353 + t545) * t347 + t351 / 0.7e1 + t545 * t353 + t354 - 0.3e1 / 0.7e1 * t556 + t340 / 0.42e2;
t288 = -t341 / 0.4e1;
t655 = t288 + t347 / 0.2e1;
t220 = t547 + t655;
t301 = 0.4e1 / 0.3e1 * t353;
t132 = -0.8e1 / 0.3e1 * t220 * t426 + 0.5e1 / 0.18e2 * t361 + (t301 + t289) * t347 + t354 - t351 / 0.3e1 + t340 / 0.18e2 + (t295 + 0.2e1 / 0.3e1 * t353 + t291) * t355;
t179 = -t361 / 0.6e1 + t386;
t303 = t353 / 0.2e1;
t546 = t303 + t355;
t181 = -0.2e1 / 0.3e1 * t426 + t288 + t546;
t335 = 0.4e1 * t355;
t254 = (t335 + t341) * t353;
t300 = -0.2e1 / 0.3e1 * t347;
t258 = t300 + t355;
t290 = -t341 / 0.2e1;
t231 = t290 + t449;
t430 = -0.4e1 * t237;
t402 = t231 * t430;
t431 = 0.16e2 * t509;
t274 = t276 ^ 2;
t558 = t351 * t274;
t501 = 0.8e1 * t558;
t110 = t258 * t501 + t181 * t431 + 0.14e2 * t131 * t557 + t268 * t361 + (t254 - 0.10e2 / 0.3e1 * t351 + 0.2e1 * t354 - t556) * t347 + t179 * t562 + (0.6e1 * t132 * t633 + t259 * t402) * pkin(1);
t292 = -0.3e1 / 0.2e1 * t341;
t548 = t340 / 0.2e1 - t361 / 0.2e1;
t408 = -0.3e1 * t556 + t332 + t548;
t553 = t267 * ((t292 + t337) * t353 - 0.3e1 / 0.2e1 * t556 + t535 + t548) + t344;
t120 = t426 * t664 + (t312 + (t319 - 0.9e1 * t341) * t353 + t408) * t347 + (t292 + t544) * t361 + t553;
t398 = pkin(1) * t402;
t534 = t354 - t351;
t125 = t257 * t398 - t344 + (-t278 - t532) * t361 + (t254 + t361 / 0.6e1 - t340 / 0.6e1 + t534) * t347 + t179 * t355;
t321 = -0.5e1 * t341;
t327 = 0.7e1 * t351;
t129 = (t292 + t336 + 0.7e1 * t353) * t361 + (t327 + (t321 + 0.10e2 * t355) * t353 + t408) * t347 + t553;
t365 = pkin(7) * t355;
t249 = -0.12e2 * pkin(7) * t357 + t365 * t666;
t256 = -0.8e1 * t351 + 0.12e2 * t555;
t141 = t249 * t310 + t256 * t276 + t431 + t501 + t535 - 0.6e1 * t555;
t154 = t224 * t563 + t231 * t266;
t211 = 0.16e2 * (t533 - 0.6e1 * t560) * t351;
t264 = -0.30e2 * t341 + 0.60e2 * t355;
t270 = t272 ^ 2;
t538 = t340 - t361;
t394 = 0.6e1 * t354 + t538 - 0.6e1 * t556;
t579 = t260 * (-t347 + t250);
t137 = t398 + (t329 + t539) * t347 + t381;
t152 = t231 * t563 + t195;
t524 = pkin(7) * t263;
t468 = 0.6e1 * t524;
t111 = t137 * t468 + t152 * t476 + t120 + t210;
t333 = 0.8e1 * t355;
t171 = t357 * t430 + t277 + 0.4e1 * t351 + (t320 + t333) * t353;
t177 = -t353 + t397 + t655;
t123 = t438 + t171 * t276 + t231 * t269 + (t177 * t528 + t429 * t562) * pkin(1);
t652 = t111 * t516 + t123 * t482;
t88 = t211 * t274 + 0.32e2 * t154 * t509 + 0.24e2 * t125 * t557 + (t320 + t335 + 0.28e2 * t353) * t344 + t235 * t579 + (0.24e2 * t110 * t272 + t264 * t351 + t354 * t322 + t394 * t329 + t538 * t337 + 0.28e2 * t357 ^ 2 + 0.4e1 * t365 ^ 2) * t347 + t652 * t242 + 0.8e1 * (t120 * t633 - t129 * t237) * pkin(1) + (0.16e2 * t141 * t270 + t264 * t353 + 0.70e2 * t351 + t361 + t394) * t361;
t216 = 0.7e1 / 0.6e1 * t347 + t287 + t546;
t297 = t347 / 0.3e1;
t454 = t287 + t297 + t355;
t219 = t301 + t454;
t149 = -t216 * t616 + t219 * t262;
t228 = t353 + t454;
t253 = t324 + t268;
t160 = t228 * t262 - t253 * t616 / 0.2e1;
t455 = t341 / 0.3e1 + t297 + t337;
t407 = -0.8e1 / 0.3e1 * t558 + t347 * t268 - 0.5e1 / 0.3e1 * t351 + t455 * t353 + t355 * (t289 + t257);
t632 = pkin(7) * t351;
t527 = -0.4e1 * t632;
t279 = -0.20e2 / 0.3e1 * t353;
t456 = 0.2e1 / 0.3e1 * t341 + t296 + t335;
t457 = 0.4e1 / 0.3e1 * t341 + t295 + t338;
t625 = t361 / 0.2e1 - (t279 + t456) * t347 / 0.2e1 + 0.3e1 / 0.2e1 * t351 - t457 * t353 / 0.2e1 - t354 / 0.2e1;
t112 = t308 * t275 * t527 + t149 * t498 + t407 * t262 + (t160 * t528 + t308 * t625) * pkin(1);
t542 = t320 - 0.2e1 * t347;
t450 = t334 + t542;
t143 = t361 + (t291 + t300 + t543) * t347 + t328 + t450 * t353 + t355 * (t291 + t258);
t139 = t143 * t262;
t225 = 0.8e1 / 0.3e1 * t347 + t453;
t227 = t289 + t296 + t540;
t153 = -t225 * t616 + t227 * t262;
t232 = 0.5e1 / 0.6e1 * t347 + t303 + t287;
t164 = t232 * t515 + t563 * t616;
t484 = t357 * t262;
t418 = t275 * t484;
t159 = t323 + (t314 + t450) * t347 + (t300 + t451) * t267;
t587 = t159 * t308;
t113 = -0.8e1 * pkin(7) * t418 + t164 * t500 + t139 + (t153 * t528 - t587) * pkin(1);
t126 = -pkin(1) * t587 + t139;
t157 = -0.3e1 * t361 + (t279 + t457) * t347 + t456 * t353 + t534;
t161 = -0.5e1 / 0.3e1 * t361 + (-t353 + t455) * t347 + t355 * (t299 + t453);
t522 = -0.2e1 * t616;
t127 = t157 * t262 + t161 * t522;
t206 = 0.4e1 * t393;
t212 = t262 + 0.2e1 * t616;
t230 = t290 + t251;
t130 = t269 * t262 + t206 * t276 + (t212 * t633 + t230 * t308) * t526;
t217 = t355 + 0.5e1 / 0.2e1 * t347 + 0.3e1 / 0.2e1 * t353 + t290;
t578 = t266 * t308;
t162 = t217 * t262 + pkin(1) * t578 / 0.2e1;
t654 = t336 - t341 - t347;
t233 = t654 * t314;
t541 = t321 - 0.5e1 * t347;
t400 = 0.24e2 * t257 * t558 - t344 - (0.21e2 * t353 + t654) * t361 - (t355 * t542 + t233 + t332 + 0.35e2 * t351) * t347 - (t327 + (t333 + t541) * t353 + t355 * (-t347 + t532)) * t267;
t481 = -0.12e2 * t574;
t183 = 0.4e1 / 0.3e1 * t557 + t246 + t259;
t577 = t270 * t361;
t656 = 0.7e1 * t344 + (0.35e2 * t353 + 0.15e2 * t355 + t541) * t361 + (0.21e2 * t351 + t233 + 0.9e1 * t354 + (t322 - 0.6e1 * t347) * t355) * t347 + t579 - 0.24e2 * t183 * t577;
t93 = t162 * t431 + t130 * t420 + t112 * t481 - 0.6e1 * t127 * t557 + (-0.6e1 * t113 * t581 + t309 * t400) * pkin(3) + (-0.6e1 * t126 * t633 + t656 * t308) * pkin(1);
t72 = t207 * t88 + t356 * t93;
t374 = t72 ^ 2;
t599 = 0.1e1 / t374 * t667;
t61 = 0.1e1 / (t118 * t599 + 0.1e1);
t382 = -0.2e1 * pkin(3) * pkin(4) * t342 * t61 * t559 * t599;
t401 = t61 / t72 * t92 * t628;
t91 = 0.1e1 / t667;
t597 = t72 * t91;
t472 = -t597 / 0.2e1;
t441 = t165 * t626;
t54 = (t108 * t441 + t96 * t627) * t396 + (t106 * t441 + t97 * t627) * t391;
t554 = t357 * t276;
t463 = t257 * t554;
t180 = 0.24e2 * t463 * t488;
t434 = pkin(7) * t487;
t395 = t353 * t308 * t434;
t203 = 0.4e1 * t395;
t214 = t434 * t649;
t435 = 0.32e2 / 0.3e1 * t351;
t405 = t275 * t435;
t406 = 0.64e2 / 0.3e1 * t218 * t357;
t422 = t563 * t632;
t445 = t114 * t628;
t573 = t272 * t360;
t464 = t242 * t573;
t469 = -0.4e1 * t524;
t479 = 0.24e2 * t574;
t490 = pkin(1) * t562;
t618 = pkin(1) * t259;
t492 = t231 * t618;
t507 = -0.8e1 * t576;
t564 = t310 * t353;
t510 = pkin(7) * t564;
t519 = -0.6e1 * t261;
t575 = t271 * t361;
t612 = t129 * pkin(1);
t639 = -0.6e1 * t356;
t643 = 0.8e1 * t207;
t90 = 0.1e1 / t92;
t605 = (t93 * t445 + (0.32e2 * t203 * t207 - 0.8e1 * t214 * t356) * t465 + (0.24e2 * (-0.4e1 * t183 * t575 * t616 - t112 * t569 - t130 * t464) * t356 + (0.48e2 * t110 * t569 + 0.96e2 * t123 * t464 + 0.64e2 * t141 * t575) * t207) * t309 + ((t88 + (t111 * t643 + t113 * t639) * t242) * t309 + (t242 * t180 * t643 + ((t219 * t500 + t228 * t469 - t407) * t481 - 0.16e2 * t217 * t509 + t157 * t497 + t143 * t468 + ((-t269 + t500) * t507 + (t227 * t469 + 0.8e1 * t232 * t557 - t143 + t438) * t519) * t242 - t400) * t356 + ((pkin(7) * t405 + 0.16e2 * t220 * t510 + t276 * t406 + 0.4e1 * t492) * t479 + 0.64e2 * t275 * t422 + 0.96e2 * t231 * t463 + 0.48e2 * t158 * t510 + 0.8e1 * t612 + ((0.2e1 * t490 + 0.4e1 * t554) * t482 + (t158 * t665 + 0.24e2 * t231 * t510) * t516) * t242) * t207 * t308) * t307) * pkin(3)) * t90;
t499 = -0.2e1 * t557;
t502 = 0.8e1 * t564;
t518 = -0.4e1 * t261;
t74 = (t214 * t505 + (-0.8e1 * t136 * t565 - t240 * t276 + ((-t268 + t499) * t505 + t187 + (t251 * t520 - 0.8e1 * t572) * pkin(7)) * pkin(3)) * t307 + ((t214 + (-t251 + 0.2e1 * t557) * t261) * t518 + (pkin(3) * t644 - 0.24e2 * t573 * t616) * t309) * t242) * t356 + t109 * t445 + (0.24e2 * t194 * t309 * t464 + (t203 + (0.8e1 / 0.3e1 * t554 + 0.2e1 * t618) * t488) * t480 + 0.24e2 * t122 * t462 + t180 + 0.24e2 * t226 * t395 + 0.6e1 * t151 * t427 + 0.6e1 * ((pkin(7) * t502 + t226 * t666) * t308 * t574 + t124 * t262) * t242) * t207 + t101 * t262;
t18 = t114 * t401 + (t605 / 0.2e1 + t74 * t472) * t382 + t54;
t196 = t566 + t571;
t197 = t568 - t570;
t621 = -t356 / 0.4e1;
t446 = t108 * t621;
t598 = t72 * t90;
t473 = t598 / 0.4e1;
t383 = t106 * t473 + t446;
t377 = t147 * t383;
t471 = -t597 / 0.4e1;
t387 = t605 / 0.4e1 + t74 * t471;
t447 = -t108 * t117 / 0.8e1;
t349 = 0.1e1 / pkin(3) ^ 2;
t561 = t342 * t349;
t25 = (-t165 * t377 + (t387 * t106 + t114 * t447 + t97 * t473 + t96 * t621) * t146) * t561;
t412 = t108 * t473;
t620 = t356 / 0.4e1;
t384 = t106 * t620 + t412;
t466 = t146 * t561;
t58 = t384 * t466;
t56 = t196 * t58;
t57 = t383 * t466;
t37 = t197 * t57 + t56;
t35 = 0.1e1 / t37 ^ 2;
t36 = -t196 * t57 + t197 * t58;
t29 = 0.1e1 / (t35 * t36 ^ 2 + 0.1e1);
t34 = 0.1e1 / t37;
t378 = t147 * t384;
t448 = t106 * t117 / 0.8e1;
t594 = (-t165 * t378 + (t387 * t108 + t114 * t448 + t96 * t473 + t97 * t620) * t146) * t561 - t57;
t607 = t35 * t36;
t2 = (-(-t196 * t25 + t594 * t197 - t56) * t34 - ((-t25 - t58) * t197 - t594 * t196) * t607) * t29 + t18;
t660 = Ifges(11,3) * t2;
t403 = pkin(7) * t419;
t215 = -0.4e1 * t403;
t245 = pkin(7) * t522;
t635 = pkin(7) * t276;
t404 = t484 * t635;
t411 = t275 * t461;
t415 = 0.32e2 * t465;
t512 = pkin(7) * t554;
t432 = -0.24e2 * t512;
t433 = -0.48e2 * t512;
t439 = pkin(7) * t498;
t444 = t121 * t628;
t460 = t308 * t564;
t477 = -0.24e2 * t564;
t478 = -0.32e2 * t275 * t351;
t483 = -0.96e2 * t257 * t275;
t503 = -0.4e1 * t572;
t513 = pkin(7) * t557;
t552 = -0.24e2 * t257 * t418 + t432 * t578;
t638 = pkin(1) * pkin(7);
t589 = -0.4e1 * t638;
t590 = -0.6e1 * t638;
t640 = -0.2e1 * t310;
t606 = (((t230 * t521 + t439 + t503) * t420 + (0.12e2 * t273 * t276 * t632 + t216 * t503 + t274 * t527) * t481 + 0.12e2 * t161 * t572 + (t625 * t481 + t656) * t263 + (t253 * t499 * t481 + t159 * t497 + t266 * t501) * pkin(7)) * t356 + t93 * t444 + t215 * t207 * t415 + (((pkin(7) * t225 * t500 - t159 * t263 - 0.4e1 * t411) * t639 + t552 * t643) * t581 + ((-pkin(7) * t274 * t435 - 0.16e2 * t220 * t513 - t275 * t406 - 0.4e1 * t310 * t492) * t479 - 0.64e2 * t274 * t422 + t357 * t231 * t483 - 0.48e2 * t158 * t513 - 0.8e1 * t310 * t612 + ((t490 * t640 + t503) * t482 + (-0.24e2 * t231 * t513 + t263 * t664) * t516) * t242) * t207 * t309) * pkin(3) + ((0.16e2 * (t256 * t640 - t249 + t433 + t478) * t577 + (t171 * t640 + t177 * t589 + t432) * t415 + (-0.28e2 * t131 * t564 + t132 * t590 + t181 * t433 + t258 * t478) * t479 + t242 * (t137 * t590 + t152 * t477) * t516 - 0.4e1 * t211 * t275 - 0.96e2 * t154 * t512 - 0.48e2 * t125 * t564 - 0.8e1 * t120 * t638) * t207 + ((-0.8e1 * t149 * t564 + t405 * t262) * t481 + t309 * t483 * t615 + t162 * t433 + 0.12e2 * t127 * t564 + (0.2e1 * (-t206 * t310 - t212 * t638) * t507 + (t153 * t589 + t164 * t502 + 0.24e2 * t404) * t519) * t242) * t356 + (-t652 * t207 - t88 + (t113 * t517 + t130 * t506 + (-0.24e2 * t245 + 0.64e2 * t460) * t577 + (0.48e2 * t160 * t574 + 0.6e1 * t126) * pkin(7)) * t356) * pkin(1)) * t308) * t90;
t504 = 0.8e1 * t574;
t73 = (t273 * t271 * t496 + (t265 * t263 - 0.2e1 * t572) * t505 + 0.4e1 * t411 + t252 * t439 + t200 * t263) * t356 + t109 * t444 + ((-0.8e1 / 0.3e1 * t418 + t215 - 0.2e1 * t259 * t428) * t480 - 0.24e2 * t226 * t403 - 0.6e1 * t151 * t428 + t552) * t207 + ((t205 * t310 * t504 + t168 * t640 - 0.24e2 * t404) * t356 + (0.12e2 * (-t173 * t564 - t512) * t480 + t150 * t477) * t207 + (0.4e1 * t133 * t356 * t261 - t385 * t207 - t101 + ((t504 * t262 + 0.4e1 * t166) * t356 + (-0.48e2 * t182 * t574 - 0.6e1 * t134) * t207) * pkin(7)) * pkin(1)) * t308 + ((t507 * t263 + ((0.4e1 * t237 * t310 - 0.2e1 * t635) * t353 + (-0.2e1 * t213 * t634 - t252 * t310) * pkin(1)) * t518) * t356 + ((t245 - 0.8e1 * t460) * t506 + (-0.8e1 * t403 - 0.2e1 * t247 * t567 + (-t178 * t634 - t226 * t486) * t666) * t517) * t207) * t242;
t17 = t121 * t401 + (t606 / 0.2e1 + t73 * t472) * t382 + t55;
t659 = Ifges(5,3) * t17;
t658 = t55 * Ifges(3,3);
t494 = pkin(1) * t588;
t414 = -t494 / 0.2e1;
t399 = t108 * t414;
t77 = t84 * t399;
t59 = t106 * t414 * t86 + t77;
t417 = t106 * t494;
t103 = t417 / 0.2e1;
t78 = t86 * t399;
t60 = t103 * t84 + t78;
t653 = mrSges(7,1) * t59 - mrSges(7,2) * t60;
t184 = t188 ^ 2;
t637 = m(5) / 0.2e1;
t636 = pkin(2) * t54;
t50 = pkin(2) * t55 + t103;
t39 = -t86 * t50 + t77;
t15 = pkin(6) * t20 + t39;
t40 = t50 * t84 + t78;
t8 = t15 * t86 - t40 * t84;
t631 = t8 * mrSges(10,1);
t9 = t15 * t84 + t40 * t86;
t630 = t9 * mrSges(10,2);
t611 = t20 * mrSges(4,1);
t610 = t20 * mrSges(4,2);
t49 = pkin(3) * t55 + t103;
t32 = (t349 * t146 * pkin(1) * t412 + t49 * t559 / 0.2e1) * t342;
t609 = t32 * mrSges(5,2);
t474 = t598 / 0.2e1;
t413 = t342 * t474;
t33 = pkin(1) * t446 * t466 + t348 * t49 * t413;
t608 = t33 * mrSges(5,1);
t14 = t18 * pkin(4) + t413 * t54;
t591 = t342 * t54;
t475 = t356 * t591;
t410 = t475 / 0.2e1;
t5 = t37 * t14 + t36 * t410;
t604 = t5 * mrSges(11,1);
t6 = -t36 * t14 + t37 * t410;
t603 = t6 * mrSges(11,2);
t388 = t606 / 0.4e1 + t73 * t471;
t596 = -(t169 * t377 + (t100 * t473 + t388 * t106 + t121 * t447 + t98 * t621) * t146) * t561 + t58;
t595 = (t169 * t378 + (t100 * t620 + t388 * t108 + t121 * t448 + t98 * t473) * t146) * t561 + t57;
t592 = mrSges(3,2) * t108;
t369 = t190 ^ 2;
t585 = t184 / t369;
t514 = t86 * t636;
t495 = t54 + t670;
t144 = 0.1e1 / (0.1e1 + t585);
t115 = -t144 * t585 - t144 + 0.1e1;
t22 = t508 * (t311 * t68 + t619 * t67) + t495;
t16 = pkin(6) * t22 - t514;
t13 = pkin(4) * t17 + t33;
t12 = (t16 + t514) * t84;
t11 = t16 * t86 - t79 * t636;
t10 = t495 - t670;
t4 = -t13 * t36 + t32 * t37;
t3 = t13 * t37 + t32 * t36;
t1 = (-(t596 * t196 + t595 * t197) * t34 - (-t595 * t196 + t596 * t197) * t607) * t29 + t17;
t19 = [-0.2e1 * t40 * t610 + 0.2e1 * t39 * t611 + (t3 ^ 2 + t4 ^ 2) * m(11) + (t8 ^ 2 + t9 ^ 2) * m(10) + (t39 ^ 2 + t40 ^ 2) * m(4) + m(3) * (t107 / 0.2e1 + t368 / 0.2e1) * t349 * t147 * t303 + Ifges(2,3) + Ifges(6,3) + (-0.2e1 * t630 + 0.2e1 * t631 + t661) * t7 + (mrSges(3,1) * t417 - t494 * t592 + t658) * t55 + m(9) * (t369 + t184) * t353 + m(7) * (t59 ^ 2 + t60 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + (0.2e1 * t608 - 0.2e1 * t609 + t659) * t17 + (0.2e1 * mrSges(11,1) * t3 - 0.2e1 * mrSges(11,2) * t4 + Ifges(11,3) * t1) * t1 + (0.2e1 * t653 + t669) * t20 + ((t190 * mrSges(9,1) + t188 * mrSges(9,2)) * t649 + Ifges(9,3) * t115) * t115; (m(11) * t6 - t2 * mrSges(11,2)) * t4 + (m(10) * t9 - t7 * mrSges(10,2)) * t12 + (t608 - t609 + t659) * t18 + (-t630 + t631 + t661) * t10 + (-t603 + t604 + t660) * t1 + (m(11) * t5 + t2 * mrSges(11,1)) * t3 + (m(10) * t8 + t7 * mrSges(10,1)) * t11 + (mrSges(4,1) * t39 - mrSges(4,2) * t40 + t653 + t669) * t22 + (t658 + (mrSges(3,1) * t629 - t592 / 0.2e1) * t494 + ((t32 * t356 + t33 * t598) * t637 + (-t356 * mrSges(5,2) / 0.2e1 + mrSges(5,1) * t474) * t17) * t342 + ((-m(4) * t39 - t611) * t86 + (m(4) * t40 - t610) * t84) * pkin(2)) * t54; (t5 ^ 2 + t6 ^ 2) * m(11) + (t11 ^ 2 + t12 ^ 2) * m(10) + Ifges(8,3) + ((t375 + t79) * pkin(2) ^ 2 * m(4) + (t118 / 0.2e1 + t374 * t91 / 0.2e1) / pkin(4) ^ 2 * t637 + Ifges(3,3)) * t54 ^ 2 + (-0.2e1 * t603 + 0.2e1 * t604 + t660) * t2 + (mrSges(5,1) * t591 * t598 - mrSges(5,2) * t475 + Ifges(5,3) * t18) * t18 + (0.2e1 * mrSges(10,1) * t11 - 0.2e1 * mrSges(10,2) * t12 + Ifges(10,3) * t10) * t10 + (-0.2e1 * (t86 * mrSges(4,1) + t84 * mrSges(4,2)) * t636 + t657 * t22) * t22;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t19(1), t19(2); t19(2), t19(3);];
Mq = res;
