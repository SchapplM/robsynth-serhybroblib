% Calculate inertial parameters regressor of gravitation load for
% fivebar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fivebar1DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t265 = pkin(5) ^ 2;
t264 = t265 ^ 2;
t268 = pkin(4) ^ 2;
t267 = t268 ^ 2;
t444 = t264 - t267;
t530 = 4 * pkin(3);
t520 = 2 * pkin(1);
t432 = 2 * pkin(3);
t273 = pkin(3) ^ 2;
t203 = 0.10e2 / 0.3e1 * t273;
t277 = pkin(2) ^ 2;
t215 = -0.2e1 / 0.3e1 * t265;
t220 = 0.2e1 / 0.3e1 * t268;
t279 = pkin(1) ^ 2;
t262 = 2 * t279;
t355 = t215 + t220 + t262;
t190 = t273 + t279;
t285 = t277 ^ 2;
t354 = t215 + t190;
t457 = t190 * (t220 + t354) + t285;
t84 = (t203 + t355) * t277 + t457;
t529 = -0.6e1 * t84;
t271 = t273 ^ 2;
t278 = t279 ^ 2;
t442 = t271 + t278;
t445 = t262 - t265;
t462 = t279 * t265;
t521 = -t264 / 0.6e1 + t267 / 0.6e1;
t110 = t445 * t273 + t442 - t462 - t521;
t297 = t110 + t285;
t86 = (t203 + t445) * t277 + t297;
t528 = -0.6e1 * t86;
t233 = sin(qJ(1));
t188 = pkin(2) * t233;
t430 = 0.2e1 * t188;
t511 = 4 * t273;
t252 = 6 * t273;
t213 = -t265 / 0.3e1;
t221 = t268 / 0.3e1;
t440 = t277 + t279;
t352 = t273 + t440;
t136 = t213 + t221 + t352;
t527 = -0.24e2 * t136;
t242 = 10 * t273;
t234 = cos(qJ(2));
t187 = pkin(3) * t234;
t419 = pkin(1) * t187;
t170 = 0.2e1 * t419;
t182 = -t273 / 0.3e1 + t279;
t197 = t234 ^ 2;
t464 = t273 * t197;
t114 = 0.4e1 / 0.3e1 * t464 + t170 + t182;
t261 = 3 * t279;
t522 = t261 - t265 - t268;
t159 = t522 * t242;
t246 = -0.6e1 * t265;
t284 = pkin(2) * t277;
t274 = t284 ^ 2;
t245 = -0.5e1 * t265;
t446 = t245 - 0.5e1 * t268;
t235 = cos(qJ(1));
t200 = t235 ^ 2;
t198 = t200 ^ 2;
t481 = t198 * t285;
t185 = t190 ^ 2;
t437 = t279 - t265;
t351 = t273 + t437;
t483 = t185 * (-t268 + t351);
t526 = 0.7e1 * t274 + ((35 * t273) + (15 * t279) + t446) * t285 + ((21 * t271) + t159 + (9 * t278) + (t246 - 0.6e1 * t268) * t279) * t277 + t483 - 0.24e2 * t114 * t481;
t232 = sin(qJ(2));
t525 = t232 * t277;
t192 = -3 * t273 + t279;
t406 = 0.4e1 * t464;
t524 = t192 + t406;
t523 = t213 - t268 / 0.3e1;
t519 = -0.4e1 * pkin(2);
t518 = -4 * pkin(3);
t517 = -2 * pkin(3);
t193 = -0.3e1 * t277 + t279;
t196 = t234 * t197;
t288 = pkin(3) * t273;
t482 = t196 * t288;
t404 = pkin(1) * t482;
t348 = 0.8e1 * t404;
t141 = t193 * t348;
t362 = 0.6e1 * t419;
t385 = 0.12e2 * t464;
t216 = -0.3e1 / 0.2e1 * t265;
t240 = 15 * t271;
t243 = 18 * t279;
t257 = 3 * t278;
t454 = t264 / 0.2e1 - t267 / 0.2e1;
t324 = -0.3e1 * t462 + t257 + t454;
t476 = t232 * t233;
t162 = pkin(2) * t476;
t342 = pkin(3) * t162;
t449 = 15 * t273 + t261;
t488 = t274 + t190 * ((t216 + t262) * t273 - 0.3e1 / 0.2e1 * t462 + t442 + t454);
t54 = t342 * t528 + (t240 + (t243 - 0.9e1 * t265) * t273 + t324) * t277 + (t216 + t449) * t285 + t488;
t214 = -t265 / 0.2e1;
t158 = t214 + t352;
t326 = -0.4e1 * t342;
t307 = t158 * t326;
t69 = t307 + (t252 + t445) * t277 + t297;
t155 = -0.2e1 * t342;
t183 = t279 - t277 / 0.3e1;
t121 = t183 * t155;
t468 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t83 = t158 * t468 + t121;
t46 = t69 * t362 + t83 * t385 + t141 + t54;
t516 = 0.8e1 * t46;
t494 = pkin(2) * t235;
t167 = pkin(1) - t494;
t137 = t167 + t187;
t515 = -0.8e1 * t137;
t514 = -0.2e1 * t234;
t513 = -0.2e1 * t235;
t512 = 4 * t271;
t510 = pkin(1) * pkin(3);
t350 = t277 + t437;
t321 = t273 + t350;
t153 = -t268 + t321;
t420 = pkin(1) * t494;
t172 = -0.2e1 * t420;
t116 = t172 + t153;
t151 = t167 * t234;
t201 = t277 * t511;
t463 = t277 * t279;
t173 = t201 - 0.4e1 * t463;
t247 = 0.2e1 * t268;
t263 = -2 * t279;
t478 = t200 * t277;
t413 = 0.2e1 * t478;
t441 = -t277 + t279;
t119 = t172 + t413 + t441;
t486 = t119 * t197;
t503 = -pkin(4) + pkin(5);
t504 = -pkin(4) - pkin(5);
t96 = t155 + t116;
t280 = sqrt(t173 * t200 + 0.4e1 * t153 * t420 - t271 - (t279 + (pkin(2) - t503) * (pkin(2) + t503)) * (t279 + (pkin(2) - t504) * (pkin(2) + t504)) + (t247 + t263 + 0.2e1 * t265 - 0.6e1 * t277 - 0.4e1 * t486) * t273 + (t116 * t162 - t96 * t151) * t530);
t118 = t170 + t524;
t160 = t268 + t351;
t244 = -0.2e1 * t265;
t168 = pkin(1) + t187;
t199 = t235 * t200;
t480 = t199 * t284;
t372 = t168 * t480;
t393 = pkin(3) * t476;
t485 = t168 * t235;
t227 = t273 / 0.3e1;
t100 = -0.4e1 / 0.9e1 * t342 + t279 + t277 / 0.3e1 + t227 + t268 / 0.9e1 - t265 / 0.9e1;
t211 = -t265 / 0.6e1;
t230 = t277 / 0.2e1;
t451 = t230 + t279;
t302 = -t342 + t451;
t108 = t268 / 0.6e1 + t211 + t302;
t433 = 4 * pkin(1);
t469 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t55 = t182 * t155 + 0.6e1 * t100 * t464 + t136 * t469 + (t108 * t187 + t482) * t433;
t106 = t155 + t136;
t174 = t441 * t511;
t363 = 0.4e1 * t419;
t66 = t136 * t326 + (t252 + t355) * t277 + t457;
t56 = t106 * t363 + t174 * t197 + t66;
t81 = t136 * t468 + t121;
t39 = -0.8e1 * t118 * t372 + t141 + t81 * t385 + t66 * t362 + t274 + (-t265 + t268 + t449) * t285 + t185 * t160 + (0.12e2 * t55 * t200 + t240 + (t243 + t246 + 0.6e1 * t268) * t273 + t257 + (t244 + t247) * t279) * t277 + 0.6e1 * (-t84 * t393 - t56 * t485) * pkin(2);
t175 = pkin(2) * t512 + 0.8e1 * t273 * t284;
t369 = t288 * t468;
t101 = t175 * t233 + 0.4e1 * t232 * t369;
t251 = 5 * t271;
t438 = t278 + t285;
t448 = t242 + t262;
t461 = t279 * t273;
t117 = t448 * t277 + t251 + t438 + (6 * t461);
t255 = 0.5e1 * t285;
t259 = 6 * t279;
t126 = t255 + (t242 + t259) * t277 + t185;
t334 = 0.8e1 * t372;
t414 = -0.4e1 * t478;
t497 = pkin(1) * t234;
t253 = 3 * t273;
t178 = t253 + t440;
t142 = t178 * t188;
t491 = pkin(3) * t232;
t149 = t188 - t491;
t398 = t273 * t188;
t333 = t197 * t398;
t256 = 0.3e1 * t277;
t176 = t256 + t190;
t484 = t176 * t232;
t67 = -0.2e1 * t333 + t142 + (0.2e1 * t149 * t497 - t484) * pkin(3);
t303 = -t232 * t288 + t398;
t139 = 0.2e1 * t303;
t189 = (2 * t273) + t277;
t191 = -t273 + t279;
t70 = t189 * t491 + t139 * t197 + (t191 + t170) * t188;
t99 = -pkin(3) * t484 + t142;
t44 = t70 * t414 + t101 * t197 + (-0.4e1 * t99 * t497 + (t126 + t334) * t232) * pkin(3) + (0.4e1 * t67 * t485 + (-t117 + t348) * t233) * pkin(2);
t29 = t137 * t39 + t280 * t44;
t509 = 0.1e1 / t29 / 0.4e1;
t508 = -0.1e1 / t29 ^ 2 / 0.4e1;
t51 = 0.1e1 / t280;
t507 = t51 / 0.2e1;
t394 = pkin(3) * t151;
t82 = t155 + t172 + t352 + 0.2e1 * t394;
t79 = 0.1e1 / t82;
t506 = t79 / 0.2e1;
t80 = 0.1e1 / t82 ^ 2;
t505 = -t80 / 0.2e1;
t502 = 0.4e1 / 0.3e1 * t277;
t501 = -t280 / 0.4e1;
t500 = t280 / 0.4e1;
t474 = t234 * t233;
t399 = pkin(2) * t474;
t109 = t167 * t232 + t399;
t102 = t109 * t432;
t258 = 8 * t279;
t105 = -0.4e1 * t288 * t162 + t201 + t512 + (t244 + t258) * t273;
t212 = -t265 / 0.4e1;
t111 = t212 - t273 + t302;
t228 = t273 / 0.2e1;
t113 = -0.2e1 / 0.3e1 * t342 + t279 + t228 + t212;
t453 = t211 - t268 / 0.6e1;
t358 = t279 + t453;
t129 = t502 + t228 + t358;
t356 = t214 - t268 / 0.2e1 + t279;
t131 = 0.3e1 / 0.2e1 * t277 + t253 + t356;
t357 = t279 + t523;
t322 = t273 + t357;
t134 = t256 + t322;
t138 = 0.4e1 * t303;
t143 = 0.16e2 * (t438 - 0.6e1 * t463) * t271;
t312 = pkin(1) * t333;
t144 = -0.4e1 * t312;
t150 = t188 + 0.2e1 * t491;
t281 = pkin(1) * t279;
t166 = -12 * pkin(1) * t288 + t281 * t530;
t423 = -0.2e1 * t491;
t169 = pkin(1) * t423;
t177 = 0.2e1 * t277 + t191;
t181 = -8 * t271 + 12 * t461;
t184 = t279 - 0.2e1 / 0.3e1 * t277;
t194 = t232 ^ 2;
t195 = t197 ^ 2;
t186 = -0.30e2 * t265 + (60 * t279);
t260 = 4 * t279;
t304 = (6 * t278) + t444 - 0.6e1 * t462;
t329 = -0.32e2 * t372;
t179 = (t260 + t265) * t273;
t392 = 0.16e2 * t482;
t466 = t271 * t195;
t409 = 0.8e1 * t466;
t147 = t279 + t277 / 0.4e1 + t273 / 0.4e1 - t265 / 0.8e1;
t450 = 0.4e1 / 0.7e1 * t279 - t265 / 0.7e1;
t64 = -0.32e2 / 0.21e2 * t147 * t342 + t285 / 0.7e1 + (0.16e2 / 0.21e2 * t273 + t450) * t277 + t271 / 0.7e1 + t450 * t273 + t278 - 0.3e1 / 0.7e1 * t462 + t444 / 0.42e2;
t148 = t212 + t227 + t451;
t226 = 0.4e1 / 0.3e1 * t273;
t65 = -0.8e1 / 0.3e1 * t148 * t342 + t285 / 0.3e1 + (t226 + t213) * t277 + t278 - t271 / 0.3e1 + (t502 + 0.2e1 / 0.3e1 * t273 + t215) * t279 + t444 / 0.18e2;
t45 = t184 * t409 + 0.14e2 * t64 * t464 + t182 * t307 + t191 * t285 + (t179 - 0.10e2 / 0.3e1 * t271 + (2 * t278) - t462) * t277 + t110 * t469 + (t113 * t392 + 0.6e1 * t65 * t187) * pkin(1);
t417 = 0.8e1 * t482;
t421 = 0.4e1 * t187;
t57 = t155 * t469 + t105 * t197 + t158 * t192 + (t111 * t421 + t417) * pkin(1);
t439 = t278 - t271;
t58 = t183 * t307 - t274 + (-t203 - t437) * t285 + (t179 + t439 + t521) * t277 + t279 * t110;
t250 = 7 * t271;
t63 = (t216 + t261 + (7 * t273)) * t285 + (t250 + (t245 + (10 * t279)) * t273 + t324) * t277 + t488;
t336 = pkin(1) * t392;
t75 = t166 * t234 + t181 * t197 + t336 + t409 + t442 - (6 * t461);
t87 = t155 * t468 + t158 * t193;
t26 = t57 * t329 + t143 * t195 + 0.24e2 * t58 * t464 + (t244 + t260 + (28 * t273)) * t274 + t160 * t483 + (t186 * t271 + 0.24e2 * t45 * t200 + t278 * t246 + t252 * t304 + t262 * t444 + (4 * t281 ^ 2) + (28 * t288 ^ 2)) * t277 + 0.8e1 * (-t393 * t63 - t46 * t485) * pkin(2) + (0.8e1 * t54 * t187 + 0.32e2 * t87 * t482) * pkin(1) + (t186 * t273 + 0.16e2 * t198 * t75 + (70 * t271) + t285 + t304) * t285;
t447 = t244 - 0.2e1 * t268;
t308 = 0.24e2 * t183 * t466 - t274 - ((21 * t273) + t522) * t285 - (t447 * t279 + t159 + t257 + (35 * t271)) * t277 - (t250 + (t258 + t446) * t273 + t279 * (-t268 + t437)) * t190;
t390 = -0.12e2 * t478;
t408 = -0.6e1 * t464;
t359 = t265 / 0.3e1 + t221 + t262;
t325 = -0.8e1 / 0.3e1 * t466 + t277 * t191 - 0.5e1 / 0.3e1 * t271 + t359 * t273 + t279 * t357;
t382 = -t491 / 0.2e1;
t465 = t271 * t196;
t323 = t230 + t358;
t130 = t226 + t323;
t76 = -t129 * t491 + t130 * t188;
t132 = t273 + t323;
t85 = t132 * t188 + t177 * t382;
t204 = -0.20e2 / 0.3e1 * t273;
t360 = 0.2e1 / 0.3e1 * t265 + t220 + t260;
t361 = 0.4e1 / 0.3e1 * t265 + 0.4e1 / 0.3e1 * t268 + t263;
t89 = -t285 + (t204 + t360) * t277 - (3 * t271) + t361 * t273 + t278;
t47 = t76 * t406 + t89 * t382 + t325 * t188 + (t85 * t187 - t232 * t465) * t433;
t397 = t288 * t188;
t332 = t196 * t397;
t407 = -0.4e1 * t464;
t353 = t259 + t447;
t225 = -0.2e1 / 0.3e1 * t268;
t452 = t215 + t225;
t77 = t285 + (t448 + t452) * t277 + t251 + t353 * t273 + t279 * (t279 + t452);
t92 = t255 + (t242 + t353) * t277 + t190 * (t225 + t354);
t59 = t77 * t188 - t92 * t491;
t135 = t178 + t523;
t78 = -t134 * t491 + t135 * t188;
t146 = t228 + t277 + t453;
t93 = t146 * t430 + t468 * t491;
t48 = t93 * t407 + (t78 * t421 - 0.8e1 * t332) * pkin(1) + t59;
t90 = -0.3e1 * t285 + (t204 + t361) * t277 + t360 * t273 + t439;
t94 = -0.5e1 / 0.3e1 * t285 + (-t273 + t359) * t277 + t279 * t322;
t60 = t90 * t188 + t94 * t423;
t61 = t192 * t188 + t138 * t197 + (t131 * t232 + t150 * t497) * t432;
t133 = t256 + 0.3e1 / 0.2e1 * t273 + t356;
t91 = t133 * t188 + t193 * t491 / 0.2e1;
t30 = t91 * t336 + t61 * t334 + t60 * t408 + t47 * t390 + (t232 * t526 - 0.6e1 * t59 * t497) * pkin(3) + (t233 * t308 + 0.6e1 * t48 * t485) * pkin(2);
t496 = pkin(1) * t271;
t311 = -0.64e2 * t468 * t496;
t313 = pkin(1) * t197 * t397;
t349 = 0.32e2 / 0.3e1 * t271;
t317 = t196 * t349;
t318 = 0.64e2 / 0.3e1 * t147 * t288;
t327 = t196 * t369;
t330 = -0.96e2 * t158 * t183 * t288;
t460 = t288 * t197;
t402 = pkin(1) * t460;
t335 = -0.48e2 * t402;
t495 = pkin(1) * t273;
t337 = -0.24e2 * t158 * t495;
t338 = -0.16e2 * t148 * t495;
t489 = t182 * pkin(3);
t344 = -0.4e1 * t158 * t489;
t345 = t469 * t517;
t346 = -0.48e2 * t86 * t495;
t403 = pkin(1) * t464;
t347 = 0.4e1 * t403;
t472 = t234 * t273;
t368 = t232 * t472;
t410 = 0.8e1 * t472;
t493 = pkin(3) * t167;
t424 = 0.4e1 * t493;
t53 = (t119 * t410 + t96 * t424) * t232 + (t116 * t421 + 0.8e1 * t167 * t464) * t188;
t374 = t53 * t507;
t384 = -0.24e2 * t460;
t386 = -0.32e2 * t465;
t387 = -0.24e2 * t472;
t388 = 0.24e2 * t478;
t391 = -0.32e2 * t480;
t405 = pkin(1) * t466;
t411 = -0.8e1 * t472;
t415 = 0.8e1 * t480;
t416 = -0.8e1 * t480;
t418 = -0.4e1 * t482;
t422 = 0.2e1 * t187;
t426 = 0.6e1 * t494;
t428 = -0.6e1 * t494;
t429 = -0.8e1 * t494;
t434 = -0.8e1 * t63 * pkin(3);
t435 = pkin(3) * t528;
t328 = t232 * t384;
t458 = pkin(1) * t193 * t328 - 0.24e2 * t183 * t332;
t459 = -4 * t510;
t15 = ((t131 * t422 + t347 + t418) * t334 + (0.12e2 * t194 * t197 * t496 + t129 * t418 - 0.2e1 * t177 * t403 - 0.4e1 * t405) * t390 + 0.8e1 * t193 * t405 + 0.12e2 * t94 * t482 + 0.6e1 * t92 * t403 + (-t89 * t390 / 0.2e1 + t526) * t187) * t280 + t30 * t374 + t137 * t144 * t329 + ((0.6e1 * (-0.4e1 * t134 * t403 - t92 * t187 - 0.4e1 * t327) * t280 + t458 * t515) * t485 + ((-pkin(1) * t195 * t349 - t196 * t318 + t197 * t338 + t234 * t344) * t388 + t195 * t311 + t196 * t330 + t197 * t346 + t234 * t434 + ((t234 * t345 + t418) * t391 + (t197 * t337 + t234 * t435) * t429) * t168) * t137 * t233) * pkin(2) + (0.12e2 * (-(t317 * t188 + t76 * t411) * t478 - 0.8e1 * t183 * t465 * t188 - 0.4e1 * t91 * t402 + t60 * t472) * t280 + t137 * (0.16e2 * (t181 * t514 - t166 + t335 + t386) * t481 + (t113 * t335 + t184 * t386 - 0.28e2 * t64 * t472) * t388 - 0.4e1 * t143 * t196 - 0.96e2 * t87 * t402 - 0.48e2 * t58 * t472) + (-t26 + t137 * (0.32e2 * t57 * t480 + t494 * t516) + (t48 * t428 + t61 * t416 + (-0.24e2 * t169 + 0.64e2 * t368) * t481) * t280 + ((0.48e2 * t85 * t478 + 0.6e1 * t59) * t280 + t137 * (-0.144e3 * t65 * t478 - 0.8e1 * t54)) * pkin(1)) * pkin(3) + ((0.2e1 * (-t138 * t234 - t150 * t510) * t415 + (t93 * t410 + t78 * t459 + 0.24e2 * t313) * t426) * t280 + t137 * ((pkin(1) * t384 + t105 * t514 + t111 * t459) * t391 + (t387 * t83 - 0.6e1 * t69 * t510) * t429)) * t168) * t232;
t24 = t137 * t26 + t280 * t30;
t455 = -t162 + t151;
t107 = pkin(3) + t455;
t115 = t172 + t268 + t321;
t395 = t119 * t187;
t71 = t115 * t167 + 0.2e1 * t395;
t72 = t115 * t234 + (0.4e1 * t197 - 0.2e1) * t493;
t43 = t107 * t280 + t72 * t188 + t232 * t71;
t380 = t43 * t509;
t152 = t253 + t268 + t350;
t154 = t162 - pkin(3);
t431 = pkin(1) * t513;
t487 = t109 * t280;
t97 = t152 + t172 + t326;
t42 = -t97 * t151 + t487 + (t152 * t476 + t154 * t431) * pkin(2) + (-t160 - t256 + t413 - 0.2e1 * t486) * pkin(3);
t301 = t24 * t380 + t42 * t500;
t467 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t383 = t79 * t467;
t22 = t301 * t383;
t340 = pkin(3) * t399;
t389 = 0.12e2 * t478;
t412 = 0.8e1 * t233 * t277;
t425 = t136 * t518;
t427 = 0.4e1 * t494;
t436 = t67 * t519;
t471 = t234 * t277;
t23 = (t273 * t194 * t416 + (t189 * t187 - 0.2e1 * t482) * t414 + 0.4e1 * t327 + t176 * t347 + t126 * t187) * t280 + t44 * t374 + t137 * ((-0.8e1 / 0.3e1 * t332 + t144 - 0.2e1 * t182 * t340) * t389 + t312 * t527 + t340 * t529 + t458) + ((0.8e1 * t139 * t200 * t471 + t101 * t514 - 0.24e2 * t313) * t280 + t137 * (0.12e2 * (-t100 * t472 - t402) * t389 + t81 * t387) + (t235 * t280 * t436 - t39 + t137 * (t118 * t415 + t56 * t426) + ((pkin(2) * t200 * t412 + 0.4e1 * t99) * t280 + t137 * (-0.48e2 * t108 * t478 - 0.6e1 * t66)) * pkin(1)) * pkin(3)) * t232 + ((t415 * t187 + ((-pkin(3) * t176 + t162 * t511) * t234 + (-t149 * t491 - t464) * t520) * t427) * t280 + t137 * ((t169 - 0.8e1 * t368) * t416 + ((-0.2e1 * t174 * t232 + t425 * t188) * t234 + (-0.4e1 * t106 * t491 - 0.8e1 * t333) * pkin(1)) * t428)) * t168;
t381 = t42 * t509;
t300 = t24 * t381 + t43 * t501;
t298 = t80 * t300;
t366 = t107 * t507;
t36 = -t487 + t53 * t366 + t119 * t194 * t517 + t71 * t234 + (-t115 - 0.8e1 * t394) * t162;
t365 = t109 * t507;
t37 = t455 * t280 + t53 * t365 + (t167 * t97 + 0.4e1 * t395) * t232 + (t431 * t471 + (t152 * t234 + t197 * t424) * pkin(2)) * t233;
t376 = -t43 * t51 / 0.8e1;
t379 = t42 * t508;
t499 = t22 - (t102 * t298 + (t15 * t381 + t36 * t501 + t53 * t376 + (t23 * t379 + t37 * t509) * t24) * t79) * t467;
t112 = t183 * t328 * t494;
t473 = t234 * t235;
t396 = pkin(2) * t473;
t343 = pkin(1) * t396;
t305 = t273 * t232 * t343;
t128 = -0.4e1 * t305;
t145 = t343 * t432;
t477 = t200 * t284;
t371 = t168 * t477;
t171 = pkin(1) * t430;
t470 = t235 * t277;
t370 = t233 * t470;
t123 = t171 - 0.4e1 * t370;
t401 = pkin(1) * t525;
t475 = t232 * t235;
t400 = pkin(2) * t475;
t341 = pkin(3) * t400;
t456 = -0.2e1 * t341 + t171;
t49 = t123 * t407 + (pkin(1) * t153 * t519 + t173 * t513) * t233 + (-t456 * t151 + 0.2e1 * t233 ^ 2 * t401 + (t116 * t475 - t96 * t474) * pkin(2)) * t530;
t375 = t49 * t507;
t479 = t199 * t285;
t16 = t30 * t375 + (-0.32e2 * t128 * t137 + 0.8e1 * t145 * t280) * t372 + (0.24e2 * (0.4e1 * t114 * t479 * t491 - t371 * t61 + t47 * t470) * t280 + t137 * (0.96e2 * t371 * t57 - 0.48e2 * t45 * t470 - 0.64e2 * t75 * t479)) * t233 + ((t26 + (t137 * t516 - 0.6e1 * t280 * t48) * t168) * t233 + (((t130 * t406 + t132 * t363 + t325) * t390 + t133 * t336 + t90 * t408 - 0.6e1 * t77 * t419 + t308) * t280 + (t112 * t515 + (t524 * t415 + (t135 * t363 - 0.8e1 * t146 * t464 - 0.8e1 * t404 + t77) * t426) * t280) * t168 + ((-pkin(1) * t317 - t197 * t318 + t234 * t338 + t344) * t388 + t196 * t311 + t197 * t330 + t234 * t346 + t434 + ((t345 - 0.4e1 * t460) * t391 + (t234 * t337 + t435) * t429) * t168) * t137 * t232) * t235) * pkin(2);
t25 = (t145 * t414 + (t70 * t412 + t175 * t197 + ((t191 + 0.2e1 * t464) * t414 - t117 + (-0.4e1 * t178 * t187 + t417) * pkin(1)) * pkin(2)) * t235 + ((t145 + (t178 - 0.2e1 * t464) * t494) * t427 + (-0.24e2 * t477 * t491 + t436) * t233) * t168) * t280 + t44 * t375 + t39 * t188 + t137 * (0.24e2 * t118 * t233 * t371 + (t128 + (-0.8e1 / 0.3e1 * t460 - 0.2e1 * t489) * t400) * t389 - 0.24e2 * t55 * t370 + t112 + t305 * t527 + t341 * t529 + 0.6e1 * (-(pkin(1) * t411 + t425) * t200 * t525 + t56 * t188) * t168);
t492 = pkin(3) * t197;
t34 = t49 * t366 + t123 * t232 * t422 + ((-t232 * t280 + t72) * t235 + (t234 * t280 + (t167 * t520 + t115) * t232 + (-pkin(3) + 0.2e1 * t492 + t497) * t430) * t233) * pkin(2);
t35 = (t162 + t396) * t280 + t49 * t365 - 0.2e1 * t123 * t492 - (t171 - 0.4e1 * t341) * t151 - 0.2e1 * t200 * t401 + t152 * t400 + (t470 * t518 + (t154 * t520 - t234 * t97) * pkin(2)) * t233;
t98 = 0.2e1 * t340 + t456;
t498 = t22 + (-t98 * t298 + (t16 * t381 + t34 * t501 + t49 * t376 + (t25 * t379 + t35 * t509) * t24) * t79) * t467;
t125 = -t474 + t475;
t124 = -t473 - t476;
t21 = t300 * t383;
t20 = t21 * t124;
t14 = -t125 * t22 - t20;
t12 = 0.1e1 / t14 ^ 2;
t19 = t21 * t125;
t13 = -t124 * t22 + t19;
t490 = t12 * t13;
t378 = t43 * t508;
t377 = t42 * t51 / 0.8e1;
t373 = t98 * t505;
t367 = t102 * t505;
t269 = 0.1e1 / pkin(4);
t364 = t269 * t506;
t41 = 0.1e1 / t42 ^ 2;
t339 = pkin(4) * t269 / (t41 * t43 ^ 2 + 0.1e1) * t82;
t9 = qJ(1) + atan2(t13, t14);
t7 = sin(t9);
t8 = cos(t9);
t320 = -g(1) * t8 - g(2) * t7;
t319 = -g(1) * t7 + g(2) * t8;
t316 = 0.1e1 / t42 * t339;
t33 = qJ(2) + atan2(t43 * t364, t42 * t364);
t31 = sin(t33);
t32 = cos(t33);
t315 = g(1) * t32 + g(2) * t31;
t314 = g(1) * t31 - g(2) * t32;
t310 = g(1) * t232 - g(2) * t234;
t309 = g(1) * t233 - g(2) * t235;
t306 = t41 * t43 * t339;
t299 = t80 * t301;
t18 = 0.1e1 + 0.2e1 * (t36 * t506 - t43 * t367) * t316 - 0.2e1 * (-t42 * t367 + t37 * t506) * t306;
t17 = 0.2e1 * (t34 * t506 + t43 * t373) * t316 - 0.2e1 * (t35 * t506 + t42 * t373) * t306;
t11 = 0.1e1 / t14;
t10 = 0.1e1 / (t12 * t13 ^ 2 + 0.1e1);
t6 = (-t98 * t299 + (t35 * t500 + t49 * t377 + t16 * t380 + (t25 * t378 + t34 * t509) * t24) * t79) * t467;
t4 = (t102 * t299 + (t37 * t500 + t53 * t377 + t15 * t380 + (t23 * t378 + t36 * t509) * t24) * t79) * t467;
t2 = 0.1e1 + ((t498 * t125 + (t21 - t6) * t124) * t11 - (-t498 * t124 - t125 * t6 + t19) * t490) * t10;
t1 = ((-t124 * t4 - t499 * t125 - t20) * t11 - ((-t21 - t4) * t125 + t499 * t124) * t490) * t10;
t3 = [0, 0, 0, 0, 0, 0, t309, g(1) * t235 + g(2) * t233, 0, 0, 0, 0, 0, 0, 0, 0, t319 * t2, t320 * t2, 0, t309 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314 * t17, t315 * t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t319 * t1, t320 * t1, 0, 0, 0, 0, 0, 0, 0, 0, t310, g(1) * t234 + g(2) * t232, 0, 0, 0, 0, 0, 0, 0, 0, t314 * t18, t315 * t18, 0, t310 * pkin(3);];
taug_reg = t3;
