% Calculate inertial parameters regressor of joint inertia matrix for
% picker2Dm2TE
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
% MM_reg [((2+1)*2/2)x(2*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 14:06
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = picker2Dm2TE_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2TE_inertiaJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2TE_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t312 = cos(pkin(9));
t609 = sin(pkin(9));
t309 = sin(qJ(1));
t310 = cos(qJ(2));
t566 = t309 * t310;
t238 = pkin(3) * t566;
t425 = pkin(1) * t238;
t225 = -0.2e1 * t425;
t308 = sin(qJ(2));
t627 = 0.2e1 * t308;
t529 = pkin(7) * t627;
t244 = pkin(3) * t529;
t354 = pkin(1) ^ 2;
t356 = pkin(7) ^ 2;
t268 = t354 + t356;
t348 = pkin(3) ^ 2;
t448 = t348 + t268;
t262 = pkin(3) * t308;
t242 = t262 + pkin(7);
t311 = cos(qJ(1));
t222 = t242 * t311;
t490 = pkin(1) * t222;
t149 = t225 + t244 + t448 + 0.2e1 * t490;
t147 = 0.1e1 / t149;
t564 = t310 * t311;
t485 = pkin(3) * t564;
t580 = t242 * t309;
t177 = t485 + t580;
t525 = 0.2e1 * pkin(1);
t170 = t177 * t525;
t306 = cos(pkin(8));
t590 = sin(pkin(8));
t189 = t306 * t309 - t590 * t311;
t191 = t306 * t311 + t590 * t309;
t349 = 0.1e1 / pkin(3);
t148 = 0.1e1 / t149 ^ 2;
t325 = 0.2e1 * t348;
t331 = 0.3e1 * t354;
t342 = pkin(4) ^ 2;
t530 = t356 - t342;
t235 = t325 + t331 + t530;
t411 = -0.4e1 * t425;
t168 = t235 + t244 + t411;
t223 = t238 - pkin(1);
t251 = t354 + t530;
t273 = t308 ^ 2;
t534 = -t348 + t356;
t572 = t273 * t348;
t194 = t244 + t534 + 0.2e1 * t572;
t277 = t311 ^ 2;
t582 = t194 * t277;
t205 = t244 + t251;
t171 = t225 + t205;
t646 = 0.4e1 * t354;
t278 = t348 * t646;
t558 = t348 * t356;
t249 = t278 - 0.4e1 * t558;
t339 = -0.2e1 * t356;
t340 = 0.2e1 * pkin(3);
t352 = t354 ^ 2;
t469 = -0.4e1 * pkin(3) * pkin(7) * t251;
t649 = 0.4e1 * pkin(1);
t119 = t249 * t273 + t308 * t469 - t352 - (t356 - (t340 + pkin(4)) * pkin(4)) * (t356 + (t340 - pkin(4)) * pkin(4)) + (t339 + 0.2e1 * t342 - 0.4e1 * t348 - 0.4e1 * t582) * t354 + (-t171 * t222 + t205 * t238) * t649;
t357 = sqrt(t119);
t584 = t177 * t357;
t107 = -t168 * t222 + t584 + (t223 * t529 + t235 * t566) * pkin(3) + (-0.2e1 * t582 + (0.2e1 * t273 - 0.4e1) * t348 - t251) * pkin(1);
t203 = t325 + t205;
t264 = pkin(1) * t311;
t492 = t194 * t264;
t139 = t203 * t242 + 0.2e1 * t492;
t141 = t203 * t311 + (0.4e1 * t277 - 0.2e1) * t242 * pkin(1);
t547 = -t238 + t222;
t175 = pkin(1) + t547;
t263 = pkin(3) * t310;
t109 = t139 * t309 + t141 * t263 + t175 * t357;
t612 = t191 / 0.2e1;
t613 = t189 / 0.2e1;
t391 = t107 * t613 + t109 * t612;
t381 = t148 * t391;
t565 = t309 * t311;
t581 = t242 * t277;
t645 = 0.8e1 * t354;
t122 = (t194 * t565 + t581 * t263) * t645 + (t171 * t580 + t205 * t485) * t649;
t118 = 0.1e1 / t357;
t618 = t118 / 0.2e1;
t441 = t177 * t618;
t620 = pkin(7) * t311;
t528 = 0.2e1 * t620;
t567 = t308 * t348;
t101 = t547 * t357 + t122 * t441 + (t168 * t242 + 0.4e1 * t492) * t309 + (t528 * t567 + (t235 * t311 + t581 * t649) * pkin(3)) * t310;
t436 = t101 / 0.2e1 - t109 / 0.2e1;
t644 = t107 / 0.2e1;
t274 = t309 ^ 2;
t442 = t175 * t618;
t637 = -0.2e1 * pkin(1);
t99 = -t584 + t122 * t442 + t194 * t274 * t637 + t139 * t311 + (-t203 - 0.8e1 * t490) * t238;
t457 = t644 + t99 / 0.2e1;
t65 = (t170 * t381 + (t436 * t189 + t457 * t191) * t147) * t349;
t614 = -t189 / 0.2e1;
t390 = t107 * t612 + t109 * t614;
t586 = t147 * t349;
t377 = t390 * t586;
t96 = t391 * t586;
t87 = t312 * t377 + t609 * t96;
t376 = t87 ^ 2;
t631 = 0.1e1 / t376;
t85 = t312 * t96 - t609 * t377;
t80 = t85 ^ 2;
t639 = 0.1e1 / (t631 * t80 + 0.1e1);
t651 = t639 / t87;
t380 = t148 * t390;
t66 = (t170 * t380 + (-t457 * t189 + t436 * t191) * t147) * t349;
t653 = (-t312 * t65 + t609 * t66) * t651;
t427 = pkin(1) * t485;
t514 = 0.2e1 * t263;
t245 = pkin(7) * t514;
t569 = t308 * t309;
t487 = pkin(3) * t569;
t426 = pkin(1) * t487;
t548 = 0.2e1 * t426 + t245;
t166 = 0.2e1 * t427 + t548;
t563 = t310 * t348;
t461 = t308 * t563;
t200 = t245 + 0.4e1 * t461;
t555 = t354 * t277;
t499 = -0.4e1 * t555;
t621 = pkin(7) * t309;
t510 = t348 * t621;
t115 = t200 * t499 + (t249 * t627 + t469) * t310 + (-t548 * t222 + 0.2e1 * t310 ^ 2 * t510 + (-t171 * t564 - t205 * t569) * pkin(3)) * t649;
t520 = 0.2e1 * t264;
t607 = pkin(1) * t277;
t638 = 0.2e1 * pkin(7);
t97 = t115 * t442 + t200 * t309 * t520 + ((t309 * t357 - t141) * t308 + (t311 * t357 + (t242 * t638 + t203) * t309 + (-pkin(1) + 0.2e1 * t607 + t620) * t514) * t310) * pkin(3);
t568 = t308 * t311;
t486 = pkin(3) * t568;
t98 = (t238 - t486) * t357 + t115 * t441 - 0.2e1 * t200 * t607 - (t245 + 0.4e1 * t426) * t222 - 0.2e1 * t273 * t510 - t235 * t487 + (t567 * t649 + (-t168 * t311 + t223 * t638) * pkin(3)) * t310;
t68 = (-t166 * t381 + (t97 * t612 + t98 * t613) * t147) * t349;
t69 = (-t166 * t380 + (t98 * t612 + t97 * t614) * t147) * t349;
t652 = (-t312 * t68 + t609 * t69) * t651;
t515 = 0.8e1 * t262;
t272 = t308 * t273;
t361 = pkin(3) * t348;
t574 = t272 * t361;
t481 = 0.32e2 * t574;
t545 = t354 / 0.3e1 + t356;
t174 = -0.4e1 / 0.9e1 * t425 + 0.4e1 / 0.9e1 * t348 - t342 / 0.9e1 + t545;
t288 = -t342 / 0.6e1;
t297 = 0.2e1 / 0.3e1 * t348;
t398 = t356 - t425;
t183 = t288 + t297 + t398;
t296 = 0.4e1 / 0.3e1 * t348;
t290 = -t342 / 0.3e1;
t452 = t290 + t268;
t227 = t296 + t452;
t260 = -t354 / 0.3e1 + t356;
t428 = -0.2e1 * t238;
t496 = 0.6e1 * t555;
t276 = t311 * t277;
t358 = pkin(1) * t354;
t570 = t276 * t358;
t508 = pkin(7) * t570;
t527 = 0.4e1 * t620;
t560 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t123 = 0.4e1 * t508 + t174 * t496 + t227 * t560 + (t183 * t527 + t260 * t428) * pkin(1);
t330 = 0.6e1 * t354;
t292 = -0.2e1 / 0.3e1 * t342;
t338 = 0.2e1 * t356;
t451 = t292 + t297 + t338;
t362 = t348 ^ 2;
t450 = t292 + t268;
t549 = (t297 + t450) * t268 + t362;
t135 = t227 * t411 + (t330 + t451) * t348 + t549;
t300 = -t348 / 0.3e1;
t258 = t300 + t356;
t196 = t258 * t225;
t561 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t151 = t227 * t561 + t196;
t279 = 0.10e2 / 0.3e1 * t354;
t152 = (t279 + t451) * t348 + t549;
t267 = -0.3e1 * t348 + t356;
t437 = 0.8e1 * t508;
t211 = t267 * t437;
t236 = -t342 + t448;
t243 = t264 + pkin(7);
t261 = t268 ^ 2;
t313 = 0.15e2 * t352;
t320 = 0.18e2 * t356;
t321 = -0.2e1 * t342;
t323 = -0.6e1 * t342;
t355 = t356 ^ 2;
t333 = 0.3e1 * t355;
t345 = t361 ^ 2;
t179 = t225 + t227;
t248 = t534 * t646;
t519 = 0.4e1 * t264;
t125 = pkin(7) * t179 * t519 + t248 * t277 + t135;
t247 = pkin(7) * t520;
t270 = -0.3e1 * t354 + t356;
t497 = 0.4e1 * t555;
t195 = t247 + t497 + t270;
t505 = 0.8e1 * t574;
t516 = 0.6e1 * t262;
t386 = t125 * t516 + t195 * t505;
t475 = 0.12e2 * t555;
t479 = 0.12e2 * t572;
t337 = 0.3e1 * t356;
t542 = 0.15e2 * t354 + t337;
t648 = 0.6e1 * pkin(1);
t102 = t123 * t479 + t211 + t151 * t475 + t345 + (-t342 + t348 + t542) * t362 + (t313 + (t320 + t323 + 0.6e1 * t348) * t354 + t333 + (t321 + t325) * t356) * t348 + t261 * t236 + t386 * t243 + (t135 * t620 - t152 * t238) * t648;
t484 = t354 * t263;
t394 = -t309 * t358 + t484;
t206 = 0.2e1 * t394;
t266 = t348 + 0.2e1 * t354;
t269 = -t354 + t356;
t137 = t269 * t263 + t206 * t277 + (t311 * t245 + t266 * t309) * pkin(1);
t538 = t331 + t356;
t252 = t348 + t538;
t210 = t252 * t263;
t253 = 0.3e1 * t348 + t268;
t578 = t253 * t309;
t167 = -pkin(1) * t578 + t210;
t495 = t361 * t645;
t605 = pkin(3) * t352;
t241 = t495 + 0.4e1 * t605;
t460 = t358 * t561;
t169 = t241 * t310 + 0.4e1 * t309 * t460;
t329 = 0.5e1 * t352;
t531 = t355 + t362;
t315 = 0.10e2 * t354;
t541 = t315 + t338;
t553 = t356 * t354;
t188 = t541 * t348 + t329 + t531 + 0.6e1 * t553;
t324 = 0.5e1 * t362;
t335 = 0.6e1 * t356;
t201 = t324 + (t315 + t335) * t348 + t261;
t464 = t243 * t574;
t423 = -0.8e1 * t464;
t504 = -0.4e1 * t572;
t579 = t243 * t308;
t606 = pkin(1) * t309;
t214 = t263 - t606;
t421 = t277 * t484;
t134 = -0.2e1 * t421 + t210 + (t214 * t528 - t578) * pkin(1);
t630 = -0.4e1 * t134;
t110 = t137 * t504 + t169 * t277 + (t579 * t630 + (-t188 + t437) * t310) * pkin(3) + (-0.4e1 * t167 * t620 + (t201 + t423) * t309) * pkin(1);
t208 = t262 + t243;
t93 = t102 * t208 + t110 * t357;
t650 = t93 ^ 2;
t341 = t342 ^ 2;
t533 = t352 + t355;
t537 = t338 - t342;
t554 = t356 * t342;
t387 = t537 * t354 + t341 / 0.6e1 + t533 - t554;
t382 = 0.5e1 / 0.6e1 * t362 + t387;
t159 = (t279 + t537) * t348 + t382;
t647 = -0.6e1 * t159;
t184 = 0.4e1 / 0.3e1 * t555 + t247 + t260;
t641 = t337 - t342 - t348;
t234 = t641 * t315;
t322 = -0.5e1 * t342;
t539 = t322 - 0.5e1 * t348;
t271 = t273 ^ 2;
t575 = t271 * t362;
t577 = t261 * (-t348 + t251);
t643 = 0.7e1 * t345 + (0.35e2 * t354 + 0.15e2 * t356 + t539) * t362 + (0.21e2 * t352 + t234 + 0.9e1 * t355 + (t323 - 0.6e1 * t348) * t356) * t348 + t577 - 0.24e2 * t184 * t575;
t289 = -t342 / 0.4e1;
t642 = t289 + t348 / 0.2e1;
t293 = -0.3e1 / 0.2e1 * t342;
t546 = t341 / 0.2e1 - t362 / 0.2e1;
t410 = -0.3e1 * t554 + t333 + t546;
t551 = t268 * ((t293 + t338) * t354 - 0.3e1 / 0.2e1 * t554 + t533 + t546) + t345;
t121 = t425 * t647 + (t313 + (t320 - 0.9e1 * t342) * t354 + t410) * t348 + (t293 + t542) * t362 + t551;
t291 = -t342 / 0.2e1;
t232 = t291 + t448;
t429 = -0.4e1 * t238;
t404 = t232 * t429;
t399 = pkin(1) * t404;
t138 = t399 + (t330 + t537) * t348 + t382;
t153 = t232 * t561 + t196;
t523 = pkin(7) * t264;
t467 = 0.6e1 * t523;
t112 = t138 * t467 + t153 * t475 + t121 + t211;
t334 = 0.8e1 * t356;
t172 = t358 * t429 + t278 + 0.4e1 * t352 + (t321 + t334) * t354;
t178 = -t354 + t398 + t642;
t124 = t437 + t172 * t277 + t232 * t270 + (t178 * t527 + t428 * t560) * pkin(1);
t640 = t112 * t515 + t124 * t481;
t369 = t107 ^ 2;
t106 = 0.1e1 / t369;
t108 = t109 ^ 2;
t422 = pkin(3) / (t106 * t108 + 0.1e1) * t149 * t349;
t392 = -0.2e1 * t106 * t109 * t422;
t396 = t422 / t644;
t616 = -t148 / 0.2e1;
t439 = t170 * t616;
t617 = t147 / 0.2e1;
t56 = 0.1e1 + (-t109 * t439 + t99 * t617) * t396 + (t101 * t617 - t107 * t439) * t392;
t435 = t56 + t653;
t507 = t85 * t639 * t631;
t21 = (t312 * t66 + t609 * t65) * t507 + t435;
t634 = -0.2e1 * t21;
t633 = 0.2e1 * t21;
t185 = t189 ^ 2;
t629 = 0.8e1 * t208;
t626 = -0.2e1 * t311;
t625 = -0.6e1 * t357;
t624 = pkin(1) * pkin(7);
t440 = t166 * t616;
t55 = (t109 * t440 + t97 * t617) * t396 + (t107 * t440 + t98 * t617) * t392;
t623 = pkin(2) * t55;
t622 = pkin(7) * t277;
t619 = pkin(7) * t352;
t280 = -0.20e2 / 0.3e1 * t354;
t336 = 0.4e1 * t356;
t455 = 0.2e1 / 0.3e1 * t342 + t297 + t336;
t456 = 0.4e1 / 0.3e1 * t342 + t296 + t339;
t615 = t362 / 0.2e1 - (t280 + t455) * t348 / 0.2e1 + 0.3e1 / 0.2e1 * t352 - t456 * t354 / 0.2e1 - t355 / 0.2e1;
t611 = -t357 / 0.4e1;
t610 = t357 / 0.4e1;
t608 = pkin(1) * t260;
t328 = 0.7e1 * t352;
t130 = (t293 + t337 + 0.7e1 * t354) * t362 + (t328 + (t322 + 0.10e2 * t356) * t354 + t410) * t348 + t551;
t604 = t130 * pkin(1);
t198 = t566 - t568;
t197 = t564 + t569;
t219 = t356 + t348 / 0.4e1 + t354 / 0.4e1 - t342 / 0.8e1;
t543 = 0.4e1 / 0.7e1 * t356 - t342 / 0.7e1;
t132 = -0.32e2 / 0.21e2 * t219 * t425 + 0.5e1 / 0.42e2 * t362 + (0.16e2 / 0.21e2 * t354 + t543) * t348 + t352 / 0.7e1 + t543 * t354 + t355 - 0.3e1 / 0.7e1 * t554 + t341 / 0.42e2;
t221 = t545 + t642;
t302 = 0.4e1 / 0.3e1 * t354;
t133 = -0.8e1 / 0.3e1 * t221 * t425 + 0.5e1 / 0.18e2 * t362 + (t302 + t290) * t348 + t355 - t352 / 0.3e1 + t341 / 0.18e2 + (t296 + 0.2e1 / 0.3e1 * t354 + t292) * t356;
t180 = -t362 / 0.6e1 + t387;
t304 = t354 / 0.2e1;
t544 = t304 + t356;
t182 = -0.2e1 / 0.3e1 * t425 + t289 + t544;
t255 = (t336 + t342) * t354;
t301 = -0.2e1 / 0.3e1 * t348;
t259 = t301 + t356;
t430 = 0.16e2 * t508;
t275 = t277 ^ 2;
t556 = t352 * t275;
t500 = 0.8e1 * t556;
t111 = t259 * t500 + t182 * t430 + 0.14e2 * t132 * t555 + t269 * t362 + (t255 - 0.10e2 / 0.3e1 * t352 + 0.2e1 * t355 - t554) * t348 + t180 * t560 + (0.6e1 * t133 * t620 + t260 * t404) * pkin(1);
t532 = t355 - t352;
t126 = t258 * t399 - t345 + (-t279 - t530) * t362 + (t255 + t362 / 0.6e1 - t341 / 0.6e1 + t532) * t348 + t180 * t356;
t366 = pkin(7) * t356;
t250 = -0.12e2 * pkin(7) * t358 + t366 * t649;
t257 = -0.8e1 * t352 + 0.12e2 * t553;
t142 = t250 * t311 + t257 * t277 + t430 + t500 + t533 - 0.6e1 * t553;
t155 = t225 * t561 + t232 * t267;
t212 = 0.16e2 * (t531 - 0.6e1 * t558) * t352;
t265 = -0.30e2 * t342 + 0.60e2 * t356;
t536 = t341 - t362;
t395 = 0.6e1 * t355 + t536 - 0.6e1 * t554;
t89 = t212 * t275 + 0.32e2 * t155 * t508 + 0.24e2 * t126 * t555 + (t321 + t336 + 0.28e2 * t354) * t345 + t236 * t577 + (0.24e2 * t111 * t273 + t265 * t352 + t355 * t323 + t395 * t330 + t536 * t338 + 0.28e2 * t358 ^ 2 + 0.4e1 * t366 ^ 2) * t348 + t640 * t243 + 0.8e1 * (t121 * t620 - t130 * t238) * pkin(1) + (0.16e2 * t142 * t271 + t265 * t354 + 0.70e2 * t352 + t362 + t395) * t362;
t217 = 0.7e1 / 0.6e1 * t348 + t288 + t544;
t298 = t348 / 0.3e1;
t453 = t288 + t298 + t356;
t220 = t302 + t453;
t150 = -t217 * t606 + t220 * t263;
t229 = t354 + t453;
t254 = t325 + t269;
t161 = t229 * t263 - t254 * t606 / 0.2e1;
t454 = t342 / 0.3e1 + t298 + t338;
t409 = -0.8e1 / 0.3e1 * t556 + t348 * t269 - 0.5e1 / 0.3e1 * t352 + t454 * t354 + t356 * (t290 + t258);
t526 = -0.4e1 * t619;
t113 = t309 * t276 * t526 + t150 * t497 + t409 * t263 + (t161 * t527 + t309 * t615) * pkin(1);
t540 = t321 - 0.2e1 * t348;
t449 = t335 + t540;
t144 = t362 + (t292 + t301 + t541) * t348 + t329 + t449 * t354 + t356 * (t292 + t259);
t140 = t144 * t263;
t226 = 0.8e1 / 0.3e1 * t348 + t452;
t228 = t290 + t297 + t538;
t154 = -t226 * t606 + t228 * t263;
t233 = 0.5e1 / 0.6e1 * t348 + t304 + t288;
t165 = t233 * t514 + t561 * t606;
t483 = t358 * t263;
t420 = t276 * t483;
t160 = t324 + (t315 + t449) * t348 + (t301 + t450) * t268;
t585 = t160 * t309;
t114 = -0.8e1 * pkin(7) * t420 + t165 * t499 + t140 + (t154 * t527 - t585) * pkin(1);
t127 = -pkin(1) * t585 + t140;
t158 = -0.3e1 * t362 + (t280 + t456) * t348 + t455 * t354 + t532;
t162 = -0.5e1 / 0.3e1 * t362 + (-t354 + t454) * t348 + t356 * (t300 + t452);
t521 = -0.2e1 * t606;
t128 = t158 * t263 + t162 * t521;
t207 = 0.4e1 * t394;
t213 = t263 + 0.2e1 * t606;
t231 = t291 + t252;
t131 = t270 * t263 + t207 * t277 + (t213 * t620 + t231 * t309) * t525;
t218 = t356 + 0.5e1 / 0.2e1 * t348 + 0.3e1 / 0.2e1 * t354 + t291;
t576 = t267 * t309;
t163 = t218 * t263 + pkin(1) * t576 / 0.2e1;
t401 = 0.24e2 * t258 * t556 - t345 - (0.21e2 * t354 + t641) * t362 - (t540 * t356 + t234 + t333 + 0.35e2 * t352) * t348 - (t328 + (t334 + t539) * t354 + t356 * (-t348 + t530)) * t268;
t480 = -0.12e2 * t572;
t94 = t163 * t430 + t131 * t423 + t113 * t480 - 0.6e1 * t128 * t555 + (-0.6e1 * t114 * t579 + t401 * t310) * pkin(3) + (-0.6e1 * t127 * t620 + t643 * t309) * pkin(1);
t73 = t208 * t89 + t357 * t94;
t91 = 0.1e1 / t93;
t595 = t73 * t91;
t472 = t595 / 0.4e1;
t414 = t109 * t472;
t385 = t107 * t610 + t414;
t343 = 0.1e1 / pkin(4);
t350 = 0.1e1 / pkin(3) ^ 2;
t559 = t343 * t350;
t465 = t147 * t559;
t59 = t385 * t465;
t57 = t197 * t59;
t445 = t109 * t611;
t384 = t107 * t472 + t445;
t58 = t384 * t465;
t38 = t198 * t58 + t57;
t36 = 0.1e1 / t38 ^ 2;
t37 = -t197 * t58 + t198 * t59;
t603 = t36 * t37;
t493 = pkin(1) * t586;
t419 = t107 * t493;
t104 = t419 / 0.2e1;
t51 = pkin(2) * t56 + t104;
t416 = -t493 / 0.2e1;
t400 = t109 * t416;
t79 = t87 * t400;
t41 = t51 * t85 + t79;
t602 = t41 * t85;
t405 = pkin(7) * t421;
t216 = -0.4e1 * t405;
t246 = pkin(7) * t521;
t406 = t483 * t622;
t434 = 0.32e2 / 0.3e1 * t352;
t407 = t276 * t434;
t408 = 0.64e2 / 0.3e1 * t219 * t358;
t413 = t276 * t460;
t417 = 0.32e2 * t464;
t424 = t561 * t619;
t552 = t358 * t277;
t511 = pkin(7) * t552;
t431 = -0.24e2 * t511;
t432 = -0.48e2 * t511;
t438 = pkin(7) * t497;
t443 = t122 * t618;
t562 = t311 * t354;
t459 = t309 * t562;
t476 = -0.24e2 * t562;
t477 = -0.32e2 * t276 * t352;
t478 = 0.24e2 * t572;
t482 = -0.96e2 * t258 * t276;
t489 = pkin(1) * t560;
t491 = t232 * t608;
t498 = -0.2e1 * t555;
t501 = 0.8e1 * t562;
t502 = -0.4e1 * t570;
t506 = -0.8e1 * t574;
t512 = pkin(7) * t555;
t518 = -0.6e1 * t262;
t550 = -0.24e2 * t258 * t420 + t431 * t576;
t587 = -0.4e1 * t624;
t588 = -0.6e1 * t624;
t601 = (((t231 * t520 + t438 + t502) * t423 + (0.12e2 * t274 * t277 * t619 + t217 * t502 + t275 * t526) * t480 + 0.12e2 * t162 * t570 + (t615 * t480 + t643) * t264 + (t254 * t498 * t480 + t160 * t496 + t267 * t500) * pkin(7)) * t357 + t94 * t443 + t216 * t208 * t417 + (((pkin(7) * t226 * t499 - t160 * t264 - 0.4e1 * t413) * t625 + t550 * t629) * t579 + ((-pkin(7) * t275 * t434 - 0.16e2 * t221 * t512 - t276 * t408 - 0.4e1 * t311 * t491) * t478 - 0.64e2 * t275 * t424 + t358 * t232 * t482 - 0.48e2 * t159 * t512 - 0.8e1 * t311 * t604 + ((t489 * t626 + t502) * t481 + (-0.24e2 * t232 * t512 + t264 * t647) * t515) * t243) * t208 * t310) * pkin(3) + ((0.16e2 * (t257 * t626 - t250 + t432 + t477) * t575 + (t172 * t626 + t178 * t587 + t431) * t417 + (-0.28e2 * t132 * t562 + t133 * t588 + t182 * t432 + t259 * t477) * t478 + t243 * (t138 * t588 + t153 * t476) * t515 - 0.4e1 * t212 * t276 - 0.96e2 * t155 * t511 - 0.48e2 * t126 * t562 - 0.8e1 * t121 * t624) * t208 + ((-0.8e1 * t150 * t562 + t407 * t263) * t480 + t310 * t482 * t605 + t163 * t432 + 0.12e2 * t128 * t562 + (0.2e1 * (-t207 * t311 - t213 * t624) * t506 + (t154 * t587 + t165 * t501 + 0.24e2 * t406) * t518) * t243) * t357 + (-t640 * t208 - t89 + (t114 * t516 + t131 * t505 + (-0.24e2 * t246 + 0.64e2 * t459) * t575 + (0.48e2 * t161 * t572 + 0.6e1 * t127) * pkin(7)) * t357) * pkin(1)) * t309) * t91;
t462 = t258 * t552;
t181 = 0.24e2 * t462 * t487;
t433 = pkin(7) * t486;
t397 = t354 * t309 * t433;
t204 = 0.4e1 * t397;
t215 = t433 * t637;
t444 = t115 * t618;
t571 = t273 * t361;
t463 = t243 * t571;
t468 = -0.4e1 * t523;
t509 = pkin(7) * t562;
t573 = t272 * t362;
t600 = (t94 * t444 + (0.32e2 * t204 * t208 - 0.8e1 * t215 * t357) * t464 + (0.24e2 * (-0.4e1 * t184 * t573 * t606 - t113 * t567 - t131 * t463) * t357 + (0.48e2 * t111 * t567 + 0.96e2 * t124 * t463 + 0.64e2 * t142 * t573) * t208) * t310 + ((t89 + (t112 * t629 + t114 * t625) * t243) * t310 + (t243 * t181 * t629 + ((t220 * t499 + t229 * t468 - t409) * t480 - 0.16e2 * t218 * t508 + t158 * t496 + t144 * t467 + ((-t270 + t499) * t506 + (t228 * t468 + 0.8e1 * t233 * t555 - t144 + t437) * t518) * t243 - t401) * t357 + ((pkin(7) * t407 + 0.16e2 * t221 * t509 + t277 * t408 + 0.4e1 * t491) * t478 + 0.64e2 * t276 * t424 + 0.96e2 * t232 * t462 + 0.48e2 * t159 * t509 + 0.8e1 * t604 + ((0.2e1 * t489 + 0.4e1 * t552) * t481 + (t159 * t648 + 0.24e2 * t232 * t509) * t515) * t243) * t208 * t309) * t308) * pkin(3)) * t91;
t375 = t73 ^ 2;
t596 = 0.1e1 / t375 * t650;
t92 = 0.1e1 / t650;
t594 = t73 * t92;
t378 = t148 * t384;
t470 = -t594 / 0.4e1;
t503 = 0.8e1 * t572;
t517 = -0.4e1 * t262;
t74 = (t274 * t272 * t495 + (t266 * t264 - 0.2e1 * t570) * t504 + 0.4e1 * t413 + t253 * t438 + t201 * t264) * t357 + t110 * t443 + ((-0.8e1 / 0.3e1 * t420 + t216 - 0.2e1 * t260 * t427) * t479 - 0.24e2 * t227 * t405 - 0.6e1 * t152 * t427 + t550) * t208 + ((t206 * t311 * t503 + t169 * t626 - 0.24e2 * t406) * t357 + (0.12e2 * (-t174 * t562 - t511) * t479 + t151 * t476) * t208 + (0.4e1 * t134 * t357 * t262 - t386 * t208 - t102 + ((t503 * t263 + 0.4e1 * t167) * t357 + (-0.48e2 * t183 * t572 - 0.6e1 * t135) * t208) * pkin(7)) * pkin(1)) * t309 + ((t506 * t264 + ((0.4e1 * t311 * t238 - 0.2e1 * t622) * t354 + (-0.2e1 * t214 * t621 - t253 * t311) * pkin(1)) * t517) * t357 + ((t246 - 0.8e1 * t459) * t505 + (-0.8e1 * t405 - 0.2e1 * t248 * t565 + (-t179 * t621 - t227 * t485) * t649) * t516) * t208) * t243;
t389 = t601 / 0.4e1 + t74 * t470;
t446 = -t109 * t118 / 0.8e1;
t593 = -(t170 * t378 + (t101 * t472 + t389 * t107 + t122 * t446 + t99 * t611) * t147) * t559 + t59;
t379 = t148 * t385;
t447 = t107 * t118 / 0.8e1;
t592 = (t170 * t379 + (t101 * t610 + t389 * t109 + t122 * t447 + t99 * t472) * t147) * t559 + t58;
t75 = (t215 * t504 + (-0.8e1 * t137 * t563 - t241 * t277 + ((-t269 + t498) * t504 + t188 + (t252 * t519 - 0.8e1 * t570) * pkin(7)) * pkin(3)) * t308 + ((t215 + (-t252 + 0.2e1 * t555) * t262) * t517 + (pkin(3) * t630 - 0.24e2 * t571 * t606) * t310) * t243) * t357 + t110 * t444 + (0.24e2 * t195 * t310 * t463 + (t204 + (0.8e1 / 0.3e1 * t552 + 0.2e1 * t608) * t487) * t479 + 0.24e2 * t123 * t461 + t181 + 0.24e2 * t227 * t397 + 0.6e1 * t152 * t426 + 0.6e1 * ((pkin(7) * t501 + t227 * t649) * t309 * t572 + t125 * t263) * t243) * t208 + t102 * t263;
t388 = t600 / 0.4e1 + t75 * t470;
t591 = (-t166 * t379 + (t388 * t109 + t115 * t447 + t97 * t472 + t98 * t610) * t147) * t559 - t58;
t78 = t85 * t400;
t40 = -t87 * t51 + t78;
t589 = t343 * t55;
t473 = t595 / 0.2e1;
t415 = t343 * t473;
t50 = pkin(3) * t56 + t104;
t34 = pkin(1) * t445 * t465 + t349 * t50 * t415;
t370 = t191 ^ 2;
t583 = t185 / t370;
t557 = t349 * t357;
t145 = 0.1e1 / (0.1e1 + t583);
t116 = -t145 * t583 - t145 + 0.1e1;
t522 = t116 * t637;
t513 = t87 * t623;
t494 = t55 + t652;
t474 = t357 * t589;
t471 = -t594 / 0.2e1;
t23 = (t312 * t69 + t609 * t68) * t507 + t494;
t458 = -0.2e1 * t23 * t623;
t62 = 0.1e1 / (t119 * t596 + 0.1e1);
t383 = -0.2e1 * pkin(3) * pkin(4) * t343 * t62 * t557 * t596;
t403 = t62 / t73 * t93 * t618;
t19 = t115 * t403 + (t600 / 0.2e1 + t75 * t471) * t383 + t55;
t412 = t474 / 0.2e1;
t18 = t122 * t403 + (t601 / 0.2e1 + t74 * t471) * t383 + t56;
t402 = t55 * t415;
t61 = t85 * t104 + t79;
t60 = t87 * t107 * t416 + t78;
t54 = t55 ^ 2;
t35 = 0.1e1 / t38;
t33 = (t350 * t147 * pkin(1) * t414 + t50 * t557 / 0.2e1) * t343;
t30 = 0.1e1 / (t36 * t37 ^ 2 + 0.1e1);
t26 = (-t166 * t378 + (t388 * t107 + t115 * t446 + t98 * t472 + t97 * t611) * t147) * t559;
t22 = t23 ^ 2;
t20 = t21 ^ 2;
t17 = pkin(6) * t23 - t513;
t16 = pkin(6) * t21 + t40;
t15 = t19 * pkin(4) + t402;
t14 = pkin(4) * t18 + t34;
t13 = (t17 + t513) * t85;
t12 = t17 * t87 - t80 * t623;
t11 = t494 - t652;
t10 = t16 * t85 + t41 * t87;
t9 = t16 * t87 - t602;
t8 = t435 - t653;
t7 = t21 * t23;
t6 = -t37 * t15 + t38 * t412;
t5 = t38 * t15 + t37 * t412;
t4 = -t14 * t37 + t33 * t38;
t3 = t14 * t38 + t33 * t37;
t2 = (-(-t197 * t26 + t591 * t198 - t57) * t35 - ((-t26 - t59) * t198 - t591 * t197) * t603) * t30 + t19;
t1 = (-(t593 * t197 + t592 * t198) * t35 - (-t592 * t197 + t593 * t198) * t603) * t30 + t18;
t24 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 ^ 2, t56 * t419, -t109 * t56 * t493, 0, (t108 / 0.4e1 + t369 / 0.4e1) * t354 * t350 * t148, 0, 0, 0, 0, 0, t20, t40 * t633, t41 * t634, 0, t40 ^ 2 + t41 ^ 2, 0, 0, 0, 0, 0, t18 ^ 2, 0.2e1 * t34 * t18, -0.2e1 * t33 * t18, 0, t33 ^ 2 + t34 ^ 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t60 * t633, t61 * t634, 0, t60 ^ 2 + t61 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116 ^ 2, t191 * t522, t189 * t522, 0, (t370 + t185) * t354, 0, 0, 0, 0, 0, t8 ^ 2, 0.2e1 * t9 * t8, -0.2e1 * t10 * t8, 0, t10 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, t1 ^ 2, 0.2e1 * t3 * t1, -0.2e1 * t4 * t1, 0, t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t55, t55 * t104, t55 * t400, 0, 0, 0, 0, 0, 0, 0, t7, -t21 * t513 + t23 * t40, -t21 * t85 * t623 - t23 * t41, 0, (-t40 * t87 + t602) * t623, 0, 0, 0, 0, 0, t18 * t19, t18 * t402 + t34 * t19, -t18 * t474 / 0.2e1 - t33 * t19, 0, (t33 * t357 / 0.2e1 + t34 * t473) * t589, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t60 * t23, -t61 * t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t11, t11 * t9 + t12 * t8, -t10 * t11 - t13 * t8, 0, t10 * t13 + t12 * t9, 0, 0, 0, 0, 0, t1 * t2, t1 * t5 + t2 * t3, -t1 * t6 - t2 * t4, 0, t3 * t5 + t4 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t87 * t458, t85 * t458, 0, (t376 + t80) * t54 * pkin(2) ^ 2, 0, 0, 0, 0, 0, t19 ^ 2, t19 * t589 * t595, -t19 * t474, 0, (t119 / 0.4e1 + t375 * t92 / 0.4e1) * t54 / pkin(4) ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2, 0.2e1 * t12 * t11, -0.2e1 * t13 * t11, 0, t12 ^ 2 + t13 ^ 2, 0, 0, 0, 0, 0, t2 ^ 2, 0.2e1 * t5 * t2, -0.2e1 * t6 * t2, 0, t5 ^ 2 + t6 ^ 2;];
MM_reg = t24;