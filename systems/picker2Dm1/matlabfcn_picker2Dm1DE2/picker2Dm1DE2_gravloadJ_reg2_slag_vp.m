% Calculate inertial parameters regressor of gravitation load for
% picker2Dm1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = picker2Dm1DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t739 = 4 * pkin(1);
t343 = sin(qJ(2));
t296 = pkin(3) * t343;
t583 = 0.8e1 * t296;
t308 = t343 ^ 2;
t307 = t343 * t308;
t390 = pkin(3) ^ 2;
t407 = pkin(3) * t390;
t643 = t307 * t407;
t547 = 0.32e2 * t643;
t385 = pkin(4) ^ 2;
t344 = sin(qJ(1));
t346 = cos(qJ(2));
t635 = t344 * t346;
t272 = pkin(3) * t635;
t481 = pkin(1) * t272;
t396 = pkin(1) ^ 2;
t399 = pkin(7) ^ 2;
t612 = t396 / 0.3e1 + t399;
t209 = -0.4e1 / 0.9e1 * t481 + 0.4e1 / 0.9e1 * t390 - t385 / 0.9e1 + t612;
t323 = -t385 / 0.6e1;
t332 = 0.2e1 / 0.3e1 * t390;
t439 = t399 - t481;
t221 = t323 + t332 + t439;
t331 = 0.4e1 / 0.3e1 * t390;
t303 = t396 + t399;
t325 = -t385 / 0.3e1;
t515 = t325 + t303;
t261 = t331 + t515;
t294 = -t396 / 0.3e1 + t399;
t484 = -0.2e1 * t272;
t347 = cos(qJ(1));
t312 = t347 ^ 2;
t623 = t396 * t312;
t563 = 0.6e1 * t623;
t311 = t347 * t312;
t403 = pkin(1) * t396;
t639 = t311 * t403;
t574 = pkin(7) * t639;
t705 = pkin(7) * t347;
t594 = 0.4e1 * t705;
t629 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t151 = 0.4e1 * t574 + t209 * t563 + t261 * t629 + (t221 * t594 + t294 * t484) * pkin(1);
t366 = 6 * t396;
t466 = -0.4e1 * t481;
t327 = -0.2e1 / 0.3e1 * t385;
t375 = 0.2e1 * t399;
t514 = t327 + t332 + t375;
t408 = t390 ^ 2;
t513 = t327 + t303;
t616 = (t332 + t513) * t303 + t408;
t165 = t261 * t466 + (t366 + t514) * t390 + t616;
t259 = -0.2e1 * t481;
t335 = -t390 / 0.3e1;
t292 = t335 + t399;
t229 = t292 * t259;
t630 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t179 = t261 * t630 + t229;
t314 = 0.10e2 / 0.3e1 * t396;
t180 = (t314 + t514) * t390 + t616;
t302 = -0.3e1 * t390 + t399;
t492 = 0.8e1 * t574;
t245 = t302 * t492;
t511 = t390 + t303;
t270 = -t385 + t511;
t299 = pkin(1) * t347;
t277 = t299 + pkin(7);
t295 = t303 ^ 2;
t394 = t396 ^ 2;
t349 = 15 * t394;
t356 = 0.18e2 * t399;
t357 = -0.2e1 * t385;
t359 = -0.6e1 * t385;
t361 = 0.2e1 * t390;
t398 = t399 ^ 2;
t370 = 0.3e1 * t398;
t387 = t407 ^ 2;
t217 = t259 + t261;
t601 = -t390 + t399;
t736 = 4 * t396;
t282 = t601 * t736;
t587 = 0.4e1 * t299;
t153 = pkin(7) * t217 * t587 + t282 * t312 + t165;
t588 = 0.2e1 * t299;
t281 = pkin(7) * t588;
t305 = -(3 * t396) + t399;
t564 = 0.4e1 * t623;
t228 = t281 + t564 + t305;
t572 = 0.8e1 * t643;
t584 = 0.6e1 * t296;
t428 = t153 * t584 + t228 * t572;
t541 = 0.12e2 * t623;
t641 = t308 * t390;
t545 = 0.12e2 * t641;
t374 = 0.3e1 * t399;
t609 = (15 * t396) + t374;
t738 = 6 * pkin(1);
t118 = t151 * t545 + t245 + t179 * t541 + t387 + (-t385 + t390 + t609) * t408 + (t349 + (t356 + t359 + 0.6e1 * t390) * t396 + t370 + (t357 + t361) * t399) * t390 + t295 * t270 + t428 * t277 + (t165 * t705 - t180 * t272) * t738;
t297 = pkin(3) * t346;
t551 = t396 * t297;
t434 = -t344 * t403 + t551;
t240 = 0.2e1 * t434;
t582 = 0.2e1 * t297;
t279 = pkin(7) * t582;
t368 = 2 * t396;
t301 = t368 + t390;
t304 = -t396 + t399;
t167 = t304 * t297 + t240 * t312 + (t347 * t279 + t301 * t344) * pkin(1);
t367 = 3 * t396;
t605 = t367 + t399;
t286 = t390 + t605;
t244 = t286 * t297;
t287 = 0.3e1 * t390 + t303;
t647 = t287 * t344;
t199 = -pkin(1) * t647 + t244;
t735 = 8 * t396;
t562 = t407 * t735;
t689 = pkin(3) * t394;
t275 = t562 + 0.4e1 * t689;
t526 = t403 * t630;
t202 = t275 * t346 + 0.4e1 * t344 * t526;
t365 = 5 * t394;
t598 = t398 + t408;
t351 = 10 * t396;
t608 = t351 + t375;
t621 = t399 * t396;
t223 = t608 * t390 + t365 + t598 + 0.6e1 * t621;
t360 = 0.5e1 * t408;
t372 = 0.6e1 * t399;
t234 = t360 + (t351 + t372) * t390 + t295;
t530 = t277 * t643;
t479 = -0.8e1 * t530;
t571 = -0.4e1 * t641;
t648 = t277 * t343;
t298 = pkin(1) * t344;
t248 = -t298 + t297;
t477 = t312 * t551;
t595 = 0.2e1 * t705;
t164 = -0.2e1 * t477 + t244 + (t248 * t595 - t647) * pkin(1);
t724 = -0.4e1 * t164;
t131 = t167 * t571 + t202 * t312 + (t648 * t724 + (-t223 + t492) * t346) * pkin(3) + (-0.4e1 * t199 * t705 + (t234 + t479) * t344) * pkin(1);
t242 = t296 + t277;
t721 = 0.2e1 * t343;
t596 = pkin(7) * t721;
t278 = pkin(3) * t596;
t597 = t399 - t385;
t285 = t396 + t597;
t239 = t278 + t285;
t204 = t259 + t239;
t276 = t296 + pkin(7);
t256 = t276 * t347;
t313 = t390 * t736;
t625 = t390 * t399;
t283 = t313 - 0.4e1 * t625;
t376 = -0.2e1 * t399;
t377 = 0.2e1 * pkin(3);
t536 = -0.4e1 * pkin(3) * pkin(7) * t285;
t227 = t278 + t601 + 0.2e1 * t641;
t651 = t227 * t312;
t141 = t283 * t308 + t343 * t536 - t394 - (t399 - (t377 + pkin(4)) * pkin(4)) * (t399 + (t377 - pkin(4)) * pkin(4)) + (t376 + 0.2e1 * t385 - 0.4e1 * t390 - 0.4e1 * t651) * t396 + (-t204 * t256 + t239 * t272) * t739;
t400 = sqrt(t141);
t107 = t118 * t242 + t131 * t400;
t741 = t107 ^ 2;
t740 = 2 * pkin(1);
t378 = 2 * pkin(2);
t384 = t385 ^ 2;
t600 = t394 + t398;
t604 = t375 - t385;
t622 = t399 * t385;
t429 = t604 * t396 + t384 / 0.6e1 + t600 - t622;
t423 = 0.5e1 / 0.6e1 * t408 + t429;
t187 = (t314 + t604) * t390 + t423;
t737 = -0.6e1 * t187;
t222 = 0.4e1 / 0.3e1 * t623 + t281 + t294;
t732 = t374 - t385 - t390;
t268 = t732 * t351;
t358 = -0.5e1 * t385;
t606 = t358 - 0.5e1 * t390;
t306 = t308 ^ 2;
t644 = t306 * t408;
t646 = t295 * (-t390 + t285);
t734 = 0.7e1 * t387 + ((35 * t396) + 0.15e2 * t399 + t606) * t408 + ((21 * t394) + t268 + 0.9e1 * t398 + (t359 - 0.6e1 * t390) * t399) * t390 + t646 - 0.24e2 * t222 * t644;
t324 = -t385 / 0.4e1;
t733 = t324 + t390 / 0.2e1;
t328 = -0.3e1 / 0.2e1 * t385;
t613 = t384 / 0.2e1 - t408 / 0.2e1;
t465 = -0.3e1 * t622 + t370 + t613;
t619 = t303 * ((t328 + t375) * t396 - 0.3e1 / 0.2e1 * t622 + t600 + t613) + t387;
t148 = t481 * t737 + (t349 + (t356 - 0.9e1 * t385) * t396 + t465) * t390 + (t328 + t609) * t408 + t619;
t326 = -t385 / 0.2e1;
t266 = t326 + t511;
t485 = -0.4e1 * t272;
t457 = t266 * t485;
t440 = pkin(1) * t457;
t168 = t440 + (t366 + t604) * t390 + t423;
t181 = t266 * t630 + t229;
t590 = pkin(7) * t299;
t534 = 0.6e1 * t590;
t135 = t168 * t534 + t181 * t541 + t148 + t245;
t371 = 0.8e1 * t399;
t205 = t403 * t485 + t313 + (4 * t394) + (t357 + t371) * t396;
t216 = -t396 + t439 + t733;
t152 = t492 + t205 * t312 + t266 * t305 + (t216 * t594 + t484 * t629) * pkin(1);
t731 = t135 * t583 + t152 * t547;
t730 = 0.2e1 * pkin(7);
t729 = -2 * pkin(1);
t236 = t361 + t239;
t559 = t227 * t299;
t169 = t236 * t276 + 0.2e1 * t559;
t171 = t236 * t347 + (0.4e1 * t312 - 0.2e1) * t276 * pkin(1);
t614 = -t272 + t256;
t210 = pkin(1) + t614;
t130 = t169 * t344 + t171 * t297 + t210 * t400;
t341 = cos(pkin(8));
t680 = sin(pkin(8));
t225 = -t341 * t347 - t680 * t344;
t687 = pkin(5) * t225;
t215 = t687 * t729;
t382 = pkin(5) ^ 2;
t206 = t215 + 0.2e1 * t382;
t213 = -pkin(1) + t687;
t224 = t341 * t344 - t680 * t347;
t618 = t215 + t382;
t192 = -(t740 + pkin(5)) * pkin(5) + t618;
t193 = pkin(5) * (t740 - pkin(5)) + t618;
t401 = sqrt(-t192 * t193);
t146 = -pkin(5) * t224 * t206 - t213 * t401;
t663 = t130 * t146;
t269 = t361 + t367 + t597;
t200 = t269 + t278 + t466;
t257 = t272 - pkin(1);
t633 = t346 * t347;
t552 = pkin(3) * t633;
t649 = t276 * t344;
t212 = t552 + t649;
t654 = t212 * t400;
t129 = -t200 * t256 + t654 + (t257 * t596 + t269 * t635) * pkin(3) + (-0.2e1 * t651 + (0.2e1 * t308 - 0.4e1) * t390 - t285) * pkin(1);
t653 = t224 * t401;
t144 = pkin(5) * t653 - t206 * t213;
t665 = t129 * t144;
t431 = t665 / 0.4e1 + t663 / 0.4e1;
t557 = pkin(1) * t256;
t177 = t259 + t278 + t511 + 0.2e1 * t557;
t175 = 0.1e1 / t177;
t391 = 0.1e1 / pkin(3);
t201 = t396 + t618;
t197 = 0.1e1 / t201;
t657 = t197 / pkin(5);
t531 = t391 * t657;
t471 = t175 * t531;
t102 = t431 * t471;
t662 = t144 * t130;
t664 = t129 * t146;
t430 = t664 / 0.4e1 - t662 / 0.4e1;
t103 = t430 * t471;
t345 = sin(pkin(9));
t348 = cos(pkin(9));
t96 = t102 * t348 - t103 * t345;
t727 = 0.1e1 / t96;
t726 = 0.1e1 / t144;
t725 = 0.1e1 / t96 ^ 2;
t723 = 0.8e1 * t242;
t720 = -0.2e1 * t347;
t380 = (pkin(6) ^ 2);
t719 = 2 * t380;
t718 = -0.6e1 * t400;
t717 = pkin(1) * pkin(7);
t709 = pkin(6) * t96;
t91 = t709 * t378;
t670 = t380 + t91;
t84 = -(t378 + pkin(6)) * pkin(6) + t670;
t85 = pkin(6) * (t378 - pkin(6)) + t670;
t684 = t84 * t85;
t402 = sqrt(-t684);
t62 = 0.1e1 / t402;
t716 = -t62 / 0.2e1;
t253 = t399 + t390 / 0.4e1 + t396 / 0.4e1 - t385 / 0.8e1;
t610 = 0.4e1 / 0.7e1 * t399 - t385 / 0.7e1;
t162 = -0.32e2 / 0.21e2 * t253 * t481 + 0.5e1 / 0.42e2 * t408 + (0.16e2 / 0.21e2 * t396 + t610) * t390 + t394 / 0.7e1 + t610 * t396 + t398 - 0.3e1 / 0.7e1 * t622 + t384 / 0.42e2;
t255 = t612 + t733;
t337 = 0.4e1 / 0.3e1 * t396;
t163 = -0.8e1 / 0.3e1 * t255 * t481 + 0.5e1 / 0.18e2 * t408 + (t337 + t325) * t390 + t398 - t394 / 0.3e1 + t384 / 0.18e2 + (t331 + 0.2e1 / 0.3e1 * t396 + t327) * t399;
t218 = -t408 / 0.6e1 + t429;
t339 = t396 / 0.2e1;
t611 = t339 + t399;
t220 = -0.2e1 / 0.3e1 * t481 + t324 + t611;
t373 = 0.4e1 * t399;
t289 = (t373 + t385) * t396;
t336 = -0.2e1 / 0.3e1 * t390;
t293 = t336 + t399;
t486 = 0.16e2 * t574;
t310 = t312 ^ 2;
t624 = t394 * t310;
t567 = 0.8e1 * t624;
t134 = t293 * t567 + t220 * t486 + 0.14e2 * t162 * t623 + t304 * t408 + (t289 - 0.10e2 / 0.3e1 * t394 + 0.2e1 * t398 - t622) * t390 + t218 * t629 + (0.6e1 * t163 * t705 + t294 * t457) * pkin(1);
t599 = t398 - t394;
t154 = t292 * t440 - t387 + (-t314 - t597) * t408 + (t289 + t408 / 0.6e1 - t384 / 0.6e1 + t599) * t390 + t218 * t399;
t364 = 7 * t394;
t158 = (t328 + t374 + (7 * t396)) * t408 + (t364 + (t358 + 0.10e2 * t399) * t396 + t465) * t390 + t619;
t412 = pkin(7) * t399;
t284 = -0.12e2 * pkin(7) * t403 + t412 * t739;
t291 = -(8 * t394) + 0.12e2 * t621;
t172 = t284 * t347 + t291 * t312 + t486 + t567 + t600 - 0.6e1 * t621;
t183 = t259 * t630 + t266 * t302;
t246 = 0.16e2 * (t598 - 0.6e1 * t625) * t394;
t300 = -0.30e2 * t385 + 0.60e2 * t399;
t603 = t384 - t408;
t435 = 0.6e1 * t398 + t603 - 0.6e1 * t622;
t101 = t246 * t310 + 0.32e2 * t183 * t574 + 0.24e2 * t154 * t623 + (t357 + t373 + (28 * t396)) * t387 + t270 * t646 + (0.24e2 * t134 * t308 + t300 * t394 + t398 * t359 + t435 * t366 + t603 * t375 + (28 * t403 ^ 2) + 0.4e1 * t412 ^ 2) * t390 + t731 * t277 + 0.8e1 * (t148 * t705 - t158 * t272) * pkin(1) + (0.16e2 * t172 * t306 + t300 * t396 + (70 * t394) + t408 + t435) * t408;
t251 = 0.7e1 / 0.6e1 * t390 + t323 + t611;
t333 = t390 / 0.3e1;
t516 = t323 + t333 + t399;
t254 = t337 + t516;
t178 = -t251 * t298 + t254 * t297;
t263 = t396 + t516;
t288 = t361 + t304;
t189 = t263 * t297 - t288 * t298 / 0.2e1;
t517 = t385 / 0.3e1 + t333 + t375;
t464 = -0.8e1 / 0.3e1 * t624 + t390 * t304 - 0.5e1 / 0.3e1 * t394 + t517 * t396 + t399 * (t325 + t292);
t704 = pkin(7) * t394;
t593 = -0.4e1 * t704;
t315 = -0.20e2 / 0.3e1 * t396;
t518 = 0.2e1 / 0.3e1 * t385 + t332 + t373;
t519 = 0.4e1 / 0.3e1 * t385 + t331 + t376;
t696 = t408 / 0.2e1 - (t315 + t518) * t390 / 0.2e1 + 0.3e1 / 0.2e1 * t394 - t519 * t396 / 0.2e1 - t398 / 0.2e1;
t136 = t344 * t311 * t593 + t178 * t564 + t464 * t297 + (t189 * t594 + t344 * t696) * pkin(1);
t607 = t357 - 0.2e1 * t390;
t512 = t372 + t607;
t174 = t408 + (t327 + t336 + t608) * t390 + t365 + t512 * t396 + t399 * (t327 + t293);
t170 = t174 * t297;
t260 = 0.8e1 / 0.3e1 * t390 + t515;
t262 = t325 + t332 + t605;
t182 = -t260 * t298 + t262 * t297;
t267 = 0.5e1 / 0.6e1 * t390 + t339 + t323;
t195 = t267 * t582 + t630 * t298;
t550 = t403 * t297;
t476 = t311 * t550;
t566 = -0.4e1 * t623;
t188 = t360 + (t351 + t512) * t390 + (t336 + t513) * t303;
t658 = t188 * t344;
t137 = -0.8e1 * pkin(7) * t476 + t195 * t566 + t170 + (t182 * t594 - t658) * pkin(1);
t155 = -pkin(1) * t658 + t170;
t186 = -0.3e1 * t408 + (t315 + t519) * t390 + t518 * t396 + t599;
t190 = -0.5e1 / 0.3e1 * t408 + (-t396 + t517) * t390 + t399 * (t335 + t515);
t589 = -0.2e1 * t298;
t156 = t186 * t297 + t190 * t589;
t241 = 0.4e1 * t434;
t247 = 0.2e1 * t298 + t297;
t265 = t326 + t286;
t159 = t305 * t297 + t241 * t312 + (t247 * t705 + t265 * t344) * t740;
t252 = t399 + 0.5e1 / 0.2e1 * t390 + 0.3e1 / 0.2e1 * t396 + t326;
t645 = t302 * t344;
t191 = t252 * t297 + pkin(1) * t645 / 0.2e1;
t442 = 0.24e2 * t292 * t624 - t387 - ((21 * t396) + t732) * t408 - (t607 * t399 + t268 + t370 + (35 * t394)) * t390 - (t364 + (t371 + t606) * t396 + t399 * (-t390 + t597)) * t303;
t546 = -0.12e2 * t641;
t108 = t191 * t486 + t159 * t479 + t136 * t546 - 0.6e1 * t156 * t623 + (-0.6e1 * t137 * t648 + t346 * t442) * pkin(3) + (-0.6e1 * t155 * t705 + t344 * t734) * pkin(1);
t81 = t101 * t242 + t108 * t400;
t715 = t81 / 0.4e1;
t406 = pkin(2) ^ 2;
t88 = t406 + t670;
t714 = pkin(2) / t88 ^ 2;
t697 = t175 / 0.2e1;
t502 = t391 * t697;
t111 = qJ(1) + atan2(t130 * t502, t129 * t502);
t386 = 0.1e1 / pkin(4);
t627 = t386 * t391;
t495 = t627 / 0.2e1;
t105 = 0.1e1 / t107;
t678 = t105 * t81;
t68 = atan2(t400 * t495, t495 * t678) + t111;
t64 = sin(t68);
t713 = pkin(4) * t64;
t65 = cos(t68);
t712 = pkin(4) * t65;
t98 = t102 * t345 + t103 * t348;
t708 = pkin(6) * t98;
t89 = t91 + t719;
t90 = -pkin(2) - t709;
t39 = -t402 * t708 - t89 * t90;
t40 = -t402 * t90 + t89 * t708;
t381 = 0.1e1 / pkin(6);
t521 = t381 / t88 / 0.2e1;
t20 = atan2(t40 * t521, t39 * t521) + t111;
t18 = sin(t20);
t711 = pkin(6) * t18;
t19 = cos(t20);
t710 = pkin(6) * t19;
t707 = pkin(7) * t312;
t706 = pkin(7) * t344;
t703 = t105 / 0.2e1;
t634 = t344 * t347;
t650 = t276 * t312;
t149 = (t227 * t634 + t650 * t297) * t735 + (t204 * t649 + t239 * t552) * t739;
t309 = t344 ^ 2;
t140 = 0.1e1 / t400;
t699 = t140 / 0.2e1;
t505 = t210 * t699;
t115 = -t654 + t149 * t505 + t227 * t309 * t729 + t169 * t347 + (-t236 - 0.8e1 * t557) * t272;
t702 = t115 / 0.4e1;
t701 = t129 / 0.4e1;
t700 = t130 / 0.4e1;
t698 = t144 / 0.4e1;
t695 = -t400 / 0.4e1;
t694 = t400 / 0.4e1;
t693 = t402 / 0.2e1;
t692 = 0.2e1 * t224 ^ 2;
t691 = pkin(1) * t294;
t690 = pkin(1) * t312;
t688 = pkin(3) * t400;
t686 = t158 * pkin(1);
t637 = t343 * t347;
t231 = t635 - t637;
t638 = t343 * t344;
t230 = t633 + t638;
t524 = t678 / 0.4e1;
t427 = t129 * t694 + t130 * t524;
t626 = t386 / pkin(3) ^ 2;
t532 = t175 * t626;
t70 = t427 * t532;
t66 = t230 * t70;
t426 = t129 * t524 + t130 * t695;
t69 = t426 * t532;
t44 = -t231 * t69 - t66;
t42 = 0.1e1 / t44 ^ 2;
t67 = t231 * t70;
t43 = -t230 * t69 + t67;
t685 = t42 * t43;
t683 = t725 * t98;
t504 = t212 * t699;
t636 = t343 * t390;
t117 = t614 * t400 + t149 * t504 + (t200 * t276 + 0.4e1 * t559) * t344 + (t595 * t636 + (t269 * t347 + t650 * t739) * pkin(3)) * t346;
t203 = t212 * t740;
t176 = 0.1e1 / t177 ^ 2;
t421 = t176 * t427;
t677 = 0.1e1 / t741 * t81;
t522 = -t677 / 0.4e1;
t467 = t130 * t522;
t458 = pkin(7) * t477;
t250 = -0.4e1 * t458;
t280 = pkin(7) * t589;
t459 = t550 * t707;
t491 = 0.32e2 / 0.3e1 * t394;
t460 = t311 * t491;
t461 = 0.64e2 / 0.3e1 * t253 * t403;
t469 = t311 * t526;
t473 = 0.32e2 * t530;
t480 = t630 * t704;
t620 = t403 * t312;
t577 = pkin(7) * t620;
t487 = -0.24e2 * t577;
t488 = -0.48e2 * t577;
t493 = pkin(7) * t564;
t506 = t149 * t699;
t631 = t347 * t396;
t525 = t344 * t631;
t542 = -0.32e2 * t394 * t311;
t543 = -0.24e2 * t631;
t544 = 0.24e2 * t641;
t548 = -0.96e2 * t292 * t311;
t556 = pkin(1) * t629;
t558 = t266 * t691;
t565 = -0.2e1 * t623;
t568 = 0.8e1 * t631;
t569 = -0.4e1 * t639;
t573 = -0.8e1 * t643;
t578 = pkin(7) * t623;
t586 = -0.6e1 * t296;
t617 = -0.24e2 * t292 * t476 + t487 * t645;
t671 = -0.4e1 * t717;
t672 = -0.6e1 * t717;
t50 = ((t265 * t588 + t493 + t569) * t479 + (0.12e2 * t309 * t312 * t704 + t251 * t569 + t310 * t593) * t546 + 0.12e2 * t190 * t639 + (t696 * t546 + t734) * t299 + (t288 * t565 * t546 + t188 * t563 + t302 * t567) * pkin(7)) * t400 + t108 * t506 + t250 * t242 * t473 + (((pkin(7) * t260 * t566 - t188 * t299 - 0.4e1 * t469) * t718 + t617 * t723) * t648 + ((-pkin(7) * t310 * t491 - 0.16e2 * t255 * t578 - t311 * t461 - 0.4e1 * t347 * t558) * t544 - 0.64e2 * t310 * t480 + t403 * t266 * t548 - 0.48e2 * t187 * t578 - 0.8e1 * t347 * t686 + ((t556 * t720 + t569) * t547 + (-0.24e2 * t266 * t578 + t299 * t737) * t583) * t277) * t242 * t346) * pkin(3) + ((0.16e2 * (t291 * t720 - t284 + t488 + t542) * t644 + (t205 * t720 + t216 * t671 + t487) * t473 + (-0.28e2 * t162 * t631 + t163 * t672 + t220 * t488 + t293 * t542) * t544 + t277 * (t168 * t672 + t181 * t543) * t583 - 0.4e1 * t246 * t311 - 0.96e2 * t183 * t577 - 0.48e2 * t154 * t631 - 0.8e1 * t148 * t717) * t242 + ((-0.8e1 * t178 * t631 + t460 * t297) * t546 + t346 * t548 * t689 + t191 * t488 + 0.12e2 * t156 * t631 + (0.2e1 * (-t241 * t347 - t247 * t717) * t573 + (t182 * t671 + t195 * t568 + 0.24e2 * t459) * t586) * t277) * t400 + (-t731 * t242 - t101 + (t137 * t584 + t159 * t572 + (-0.24e2 * t280 + 0.64e2 * t525) * t644 + (0.48e2 * t189 * t641 + 0.6e1 * t155) * pkin(7)) * t400) * pkin(1)) * t344;
t509 = t129 * t140 / 0.8e1;
t483 = pkin(1) * t552;
t570 = 0.8e1 * t641;
t585 = -0.4e1 * t296;
t82 = (t309 * t307 * t562 + (t301 * t299 - 0.2e1 * t639) * t571 + 0.4e1 * t469 + t287 * t493 + t234 * t299) * t400 + t131 * t506 + ((-0.8e1 / 0.3e1 * t476 + t250 - 0.2e1 * t294 * t483) * t545 - 0.24e2 * t261 * t458 - 0.6e1 * t180 * t483 + t617) * t242 + ((t240 * t347 * t570 + t202 * t720 - 0.24e2 * t459) * t400 + (0.12e2 * (-t209 * t631 - t577) * t545 + t179 * t543) * t242 + (0.4e1 * t164 * t343 * t688 - t428 * t242 - t118 + ((t570 * t297 + 0.4e1 * t199) * t400 + (-0.48e2 * t221 * t641 - 0.6e1 * t165) * t242) * pkin(7)) * pkin(1)) * t344 + ((t573 * t299 + ((0.4e1 * t347 * t272 - 0.2e1 * t707) * t396 + (-0.2e1 * t248 * t706 - t287 * t347) * pkin(1)) * t585) * t400 + ((t280 - 0.8e1 * t525) * t572 + (-0.8e1 * t458 - 0.2e1 * t282 * t634 + (-t217 * t706 - t261 * t552) * t739) * t584) * t242) * t277;
t682 = (t203 * t421 + (t117 * t694 + t149 * t509 + t82 * t467 + (t50 * t700 + t81 * t702) * t105) * t175) * t626 + t69;
t632 = t346 * t390;
t527 = t343 * t632;
t233 = t279 + 0.4e1 * t527;
t576 = t390 * t706;
t554 = pkin(3) * t638;
t482 = pkin(1) * t554;
t615 = 0.2e1 * t482 + t279;
t138 = t233 * t566 + (t283 * t721 + t536) * t346 + (-t615 * t256 + 0.2e1 * t346 ^ 2 * t576 + (-t204 * t633 - t239 * t638) * pkin(3)) * t739;
t113 = t138 * t505 + t233 * t344 * t588 + ((t344 * t400 - t171) * t343 + (t347 * t400 + (t276 * t730 + t236) * t344 + (-pkin(1) + 0.2e1 * t690 + t705) * t582) * t346) * pkin(3);
t553 = pkin(3) * t637;
t114 = (t272 - t553) * t400 + t138 * t504 - 0.2e1 * t233 * t690 - (t279 + 0.4e1 * t482) * t256 - 0.2e1 * t308 * t576 - t269 * t554 + (t636 * t739 + (-t200 * t347 + t257 * t730) * pkin(3)) * t346;
t196 = 0.2e1 * t483 + t615;
t528 = t292 * t620;
t219 = 0.24e2 * t528 * t554;
t489 = pkin(7) * t553;
t437 = t396 * t344 * t489;
t237 = 0.4e1 * t437;
t249 = t489 * t729;
t507 = t138 * t699;
t640 = t308 * t407;
t529 = t277 * t640;
t535 = -0.4e1 * t590;
t575 = pkin(7) * t631;
t642 = t307 * t408;
t58 = t108 * t507 + (0.32e2 * t237 * t242 - 0.8e1 * t249 * t400) * t530 + (0.24e2 * (-0.4e1 * t222 * t642 * t298 - t136 * t636 - t159 * t529) * t400 + (0.48e2 * t134 * t636 + 0.96e2 * t152 * t529 + 0.64e2 * t172 * t642) * t242) * t346 + ((t101 + (t135 * t723 + t137 * t718) * t277) * t346 + (t277 * t219 * t723 + ((t254 * t566 + t263 * t535 - t464) * t546 - 0.16e2 * t252 * t574 + t186 * t563 + t174 * t534 + ((-t305 + t566) * t573 + (t262 * t535 + 0.8e1 * t267 * t623 - t174 + t492) * t586) * t277 - t442) * t400 + ((pkin(7) * t460 + 0.16e2 * t255 * t575 + t312 * t461 + 0.4e1 * t558) * t544 + 0.64e2 * t311 * t480 + 0.96e2 * t266 * t528 + 0.48e2 * t187 * t575 + 0.8e1 * t686 + ((0.2e1 * t556 + 0.4e1 * t620) * t547 + (t187 * t738 + 0.24e2 * t266 * t575) * t583) * t277) * t242 * t344) * t343) * pkin(3);
t83 = (t249 * t571 + (-0.8e1 * t167 * t632 - t275 * t312 + ((-t304 + t565) * t571 + t223 + (t286 * t587 - 0.8e1 * t639) * pkin(7)) * pkin(3)) * t343 + ((t249 + (-t286 + 0.2e1 * t623) * t296) * t585 + (pkin(3) * t724 - 0.24e2 * t640 * t298) * t346) * t277) * t400 + t131 * t507 + (0.24e2 * t228 * t346 * t529 + (t237 + (0.8e1 / 0.3e1 * t620 + 0.2e1 * t691) * t554) * t545 + 0.24e2 * t151 * t527 + t219 + 0.24e2 * t261 * t437 + 0.6e1 * t180 * t482 + 0.6e1 * ((pkin(7) * t568 + t261 * t739) * t344 * t641 + t153 * t297) * t277) * t242 + t118 * t297;
t681 = (-t196 * t421 + (t114 * t694 + t138 * t509 + t83 * t467 + (t113 * t715 + t58 * t700) * t105) * t175) * t626 - t69;
t128 = 0.1e1 / t129 ^ 2;
t478 = pkin(3) / (t128 * t130 ^ 2 + 0.1e1) * t177 * t391;
t432 = -0.2e1 * t128 * t130 * t478;
t436 = 0.2e1 / t129 * t478;
t660 = t176 * t196;
t501 = -t660 / 0.2e1;
t666 = t114 * t175;
t667 = t113 * t175;
t59 = (t667 / 0.2e1 + t130 * t501) * t436 + (t666 / 0.2e1 + t129 * t501) * t432;
t679 = t741 / t81 ^ 2;
t109 = sin(t111);
t676 = t109 * t59;
t659 = t176 * t203;
t499 = t659 / 0.2e1;
t60 = 0.1e1 + (t115 * t697 + t130 * t499) * t436 + (t117 * t697 + t129 * t499) * t432;
t675 = t109 * t60;
t110 = cos(t111);
t674 = t110 * t59;
t673 = t110 * t60;
t669 = pkin(3) * t675 + t298;
t668 = pkin(2) * t675 + t298;
t490 = t224 * t740;
t661 = 0.1e1 / t401 * (t192 + t193) * pkin(5) * t490;
t656 = t197 / pkin(1);
t198 = 0.1e1 / t201 ^ 2;
t655 = t198 * t224;
t652 = t225 * t401;
t628 = t381 / pkin(2);
t503 = t661 / 0.2e1;
t119 = (-t652 + (t213 * t740 - t206 + t503) * t224) * pkin(5);
t120 = -t213 * t661 / 0.2e1 + t382 * pkin(1) * t692 + (t206 * t225 - t653) * pkin(5);
t470 = t657 * t659;
t560 = pkin(1) * t655;
t51 = (t431 * t470 + ((t665 / 0.2e1 + t663 / 0.2e1) * t560 + (t117 * t698 + t119 * t701 + t120 * t700 + t146 * t702) * t657) * t175) * t391;
t52 = (-t430 * t470 + ((-t664 / 0.2e1 + t662 / 0.2e1) * t560 + (-t117 * t146 / 0.4e1 - t129 * t120 / 0.4e1 + t119 * t700 + t115 * t698) * t657) * t175) * t391;
t35 = t345 * t52 + t348 * t51;
t581 = t35 * t714;
t500 = -t660 / 0.4e1;
t424 = t666 / 0.4e1 + t129 * t500;
t425 = t667 / 0.4e1 + t130 * t500;
t73 = (t144 * t424 + t146 * t425) * t531;
t74 = (t144 * t425 - t146 * t424) * t531;
t48 = t345 * t74 + t348 * t73;
t580 = t48 * t714;
t38 = 0.1e1 / t39 ^ 2;
t579 = pkin(6) / (t38 * t40 ^ 2 + 0.1e1) * t88;
t549 = pkin(5) * t655;
t540 = 0.1e1 / (0.1e1 - 0.1e1 / pkin(6) ^ 2 / t406 * t725 * t684 / 0.4e1) * t628;
t539 = t90 * t716;
t538 = -t62 * t727 / 0.4e1;
t537 = t98 * t716;
t523 = -t677 / 0.2e1;
t520 = t725 * t693;
t510 = -0.2e1 * pkin(2) * t90 + t89;
t508 = -t130 * t140 / 0.8e1;
t498 = t657 / 0.2e1;
t497 = -t656 / 0.2e1;
t496 = t656 / 0.2e1;
t494 = pkin(2) * t98 * t719;
t438 = pkin(6) * (-t84 - t85) * t378;
t34 = t48 * t438;
t443 = -0.2e1 * t38 * t40 * t579;
t472 = 0.2e1 / t39 * t579;
t475 = pkin(6) * t521;
t49 = t345 * t73 - t348 * t74;
t4 = ((t34 * t539 + t48 * t494 + (t402 * t48 + t49 * t89) * pkin(6)) * t521 - t40 * t580) * t472 + (-t39 * t580 + (t34 * t537 - t402 * t49 + t510 * t48) * t475) * t443 + t59;
t72 = 0.1e1 / (t141 * t679 + 0.1e1);
t422 = -0.2e1 * pkin(4) * t72 * t627 * t679 * t688;
t441 = t107 * t72 / t81 * t699;
t17 = t138 * t441 + (t83 * t523 + t58 * t703) * t422 + t59;
t21 = t35 * t438;
t36 = t345 * t51 - t348 * t52;
t2 = ((t21 * t539 + t35 * t494 + (t35 * t402 + t36 * t89) * pkin(6)) * t521 - t40 * t581) * t472 + (-t39 * t581 + (t21 * t537 + t510 * t35 - t402 * t36) * t475) * t443 + t60;
t468 = t129 * t522;
t16 = t149 * t441 + (t50 * t703 + t82 * t523) * t422 + t60;
t463 = -pkin(2) * t673 - t299;
t462 = -pkin(3) * t673 - t299;
t14 = atan2(t628 * t693, -t96) + t20;
t12 = sin(t14);
t13 = cos(t14);
t456 = -g(1) * t13 - g(2) * t12;
t455 = -g(1) * t12 + g(2) * t13;
t454 = g(1) * t19 + g(2) * t18;
t453 = g(1) * t18 - g(2) * t19;
t30 = atan2(t43, t44) + t68;
t28 = sin(t30);
t29 = cos(t30);
t452 = g(1) * t29 + g(2) * t28;
t451 = g(1) * t28 - g(2) * t29;
t57 = atan2(t98, t96) + t111;
t55 = sin(t57);
t56 = cos(t57);
t450 = g(1) * t56 + g(2) * t55;
t449 = g(1) * t55 - g(2) * t56;
t448 = -g(1) * t65 - g(2) * t64;
t447 = -g(1) * t64 + g(2) * t65;
t446 = -g(1) * t110 - g(2) * t109;
t445 = -g(1) * t109 + g(2) * t110;
t444 = -g(1) * t344 + g(2) * t347;
t433 = t445 * t59;
t420 = t176 * t426;
t238 = t444 * pkin(1);
t214 = -pkin(1) * t225 + pkin(5);
t207 = t215 + t368;
t147 = -pkin(1) * t224 * t207 + t214 * t401;
t145 = pkin(1) * t653 + t207 * t214;
t143 = 0.1e1 / t145 ^ 2;
t142 = 0.1e1 / t144 ^ 2;
t126 = qJ(1) + atan2(t146 * t498, t144 * t498);
t125 = pkin(8) + atan2(t147 * t496, t145 * t497);
t124 = cos(t126);
t123 = sin(t126);
t122 = cos(t125);
t121 = sin(t125);
t80 = (-((t214 * t503 + t396 * pkin(5) * t692 + (t207 * t225 - t653) * pkin(1)) * t496 + t147 * t549) / t145 - (-t145 * t549 + (-t652 + (-0.2e1 * t214 * pkin(5) - t207 + t503) * t224) * pkin(1) * t497) * t147 * t143) / (t143 * t147 ^ 2 + 0.1e1) * t201 * t740;
t77 = 0.1e1 + (t120 * t726 * t657 + (-t119 * t142 * t657 + (-t144 * t142 + t726) * t198 * t490) * t146) * t201 / (t142 * t146 ^ 2 + 0.1e1) * pkin(5);
t63 = 0.1e1 / (t725 * t98 ^ 2 + 0.1e1);
t41 = 0.1e1 / t44;
t33 = 0.1e1 / (t42 * t43 ^ 2 + 0.1e1);
t25 = (-t196 * t420 + (t83 * t468 + t113 * t695 + t138 * t508 + (t114 * t715 + t58 * t701) * t105) * t175) * t626;
t23 = (t203 * t420 + (t82 * t468 + t115 * t695 + t149 * t508 + (t117 * t715 + t50 * t701) * t105) * t175) * t626;
t15 = (-t48 * t683 + t49 * t727) * t63 + t59;
t11 = (-t35 * t683 + t36 * t727) * t63 + t60;
t6 = ((-t230 * t25 + t681 * t231 - t66) * t41 - ((-t25 - t70) * t231 - t681 * t230) * t685) * t33 + t17;
t5 = ((t682 * t231 + (-t23 + t70) * t230) * t41 - (-t23 * t231 - t682 * t230 + t67) * t685) * t33 + t16;
t3 = (t34 * t538 + t48 * t520) * t540 + t4;
t1 = (t21 * t538 + t35 * t520) * t540 + t2;
t7 = [0, 0, 0, 0, 0, 0, t444, -g(1) * t347 - g(2) * t344, 0, 0, 0, 0, 0, 0, 0, 0, t445 * t60, t446 * t60, 0, t238, 0, 0, 0, 0, 0, 0, t453 * t2, t454 * t2, 0, -g(1) * t668 - g(2) * t463, 0, 0, 0, 0, 0, 0, t447 * t16, t448 * t16, 0, -g(1) * t669 - g(2) * t462, 0, 0, 0, 0, 0, 0, (g(1) * t121 - g(2) * t122) * t80, (g(1) * t122 + g(2) * t121) * t80, 0, 0, 0, 0, 0, 0, 0, 0, t449 * t11, t450 * t11, 0, t238, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (g(1) * t123 - g(2) * t124) * t77, (g(1) * t124 + g(2) * t123) * t77, 0, t238, 0, 0, 0, 0, 0, 0, t455 * t1, t456 * t1, 0, -g(2) * (t2 * t710 + t463) - g(1) * (-t2 * t711 + t668), 0, 0, 0, 0, 0, 0, t451 * t5, t452 * t5, 0, -g(2) * (-t16 * t712 + t462) - g(1) * (t16 * t713 + t669); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t433, t446 * t59, 0, 0, 0, 0, 0, 0, 0, 0, t453 * t4, t454 * t4, 0, pkin(2) * t433, 0, 0, 0, 0, 0, 0, t447 * t17, t448 * t17, 0, pkin(3) * t433, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449 * t15, t450 * t15, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t346 - g(2) * t343, g(1) * t343 - g(2) * t346, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455 * t3, t456 * t3, 0, -g(2) * (-pkin(2) * t674 + t4 * t710) - g(1) * (pkin(2) * t676 - t4 * t711), 0, 0, 0, 0, 0, 0, t451 * t6, t452 * t6, 0, -g(2) * (-pkin(3) * t674 - t17 * t712) - g(1) * (pkin(3) * t676 + t17 * t713);];
taug_reg = t7;