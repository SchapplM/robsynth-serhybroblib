% Calculate Gravitation load on the joints for
% fivebar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1DE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1DE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 04:28:41
% EndTime: 2020-04-27 04:30:25
% DurationCPUTime: 18.32s
% Computational Cost: add. (296758->655), mult. (884968->1033), div. (2224->13), fcn. (108452->12), ass. (0->460)
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
t353 = t215 + t220 + t262;
t190 = t273 + t279;
t285 = t277 ^ 2;
t352 = t215 + t190;
t457 = t190 * (t220 + t352) + t285;
t84 = (t203 + t353) * t277 + t457;
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
t350 = t273 + t440;
t136 = t213 + t221 + t350;
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
t480 = t198 * t285;
t185 = t190 ^ 2;
t437 = t279 - t265;
t349 = t273 + t437;
t483 = t185 * (-t268 + t349);
t526 = 0.7e1 * t274 + ((35 * t273) + (15 * t279) + t446) * t285 + ((21 * t271) + t159 + (9 * t278) + (t246 - 0.6e1 * t268) * t279) * t277 + t483 - 0.24e2 * t114 * t480;
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
t481 = t196 * t288;
t404 = pkin(1) * t481;
t346 = 0.8e1 * t404;
t141 = t193 * t346;
t360 = 0.6e1 * t419;
t385 = 0.12e2 * t464;
t216 = -0.3e1 / 0.2e1 * t265;
t240 = 15 * t271;
t243 = 18 * t279;
t257 = 3 * t278;
t454 = t264 / 0.2e1 - t267 / 0.2e1;
t322 = -0.3e1 * t462 + t257 + t454;
t475 = t232 * t233;
t162 = pkin(2) * t475;
t340 = pkin(3) * t162;
t449 = 15 * t273 + t261;
t488 = t274 + t190 * ((t216 + t262) * t273 - 0.3e1 / 0.2e1 * t462 + t442 + t454);
t54 = t340 * t528 + (t240 + (t243 - 0.9e1 * t265) * t273 + t322) * t277 + (t216 + t449) * t285 + t488;
t214 = -t265 / 0.2e1;
t158 = t214 + t350;
t324 = -0.4e1 * t340;
t307 = t158 * t324;
t69 = t307 + (t252 + t445) * t277 + t297;
t155 = -0.2e1 * t340;
t183 = t279 - t277 / 0.3e1;
t121 = t183 * t155;
t467 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t83 = t158 * t467 + t121;
t46 = t360 * t69 + t385 * t83 + t141 + t54;
t516 = 0.8e1 * t46;
t494 = pkin(2) * t235;
t167 = pkin(1) - t494;
t137 = t167 + t187;
t515 = -0.8e1 * t137;
t514 = -0.2e1 * t234;
t513 = -0.2e1 * t235;
t512 = 4 * t271;
t510 = pkin(1) * pkin(3);
t348 = t277 + t437;
t319 = t273 + t348;
t153 = -t268 + t319;
t420 = pkin(1) * t494;
t172 = -0.2e1 * t420;
t116 = t172 + t153;
t151 = t167 * t234;
t201 = t277 * t511;
t463 = t277 * t279;
t173 = t201 - 0.4e1 * t463;
t247 = 0.2e1 * t268;
t263 = -2 * t279;
t477 = t200 * t277;
t413 = 0.2e1 * t477;
t441 = -t277 + t279;
t119 = t172 + t413 + t441;
t486 = t119 * t197;
t503 = -pkin(4) + pkin(5);
t504 = -pkin(4) - pkin(5);
t96 = t155 + t116;
t280 = sqrt(t173 * t200 + 0.4e1 * t153 * t420 - t271 - (t279 + (pkin(2) - t503) * (pkin(2) + t503)) * (t279 + (pkin(2) - t504) * (pkin(2) + t504)) + (t247 + t263 + 0.2e1 * t265 - 0.6e1 * t277 - 0.4e1 * t486) * t273 + (t116 * t162 - t96 * t151) * t530);
t118 = t170 + t524;
t160 = t268 + t349;
t244 = -0.2e1 * t265;
t168 = pkin(1) + t187;
t199 = t235 * t200;
t479 = t199 * t284;
t370 = t168 * t479;
t393 = pkin(3) * t475;
t485 = t168 * t235;
t227 = t273 / 0.3e1;
t100 = -0.4e1 / 0.9e1 * t340 + t279 + t277 / 0.3e1 + t227 + t268 / 0.9e1 - t265 / 0.9e1;
t211 = -t265 / 0.6e1;
t230 = t277 / 0.2e1;
t451 = t230 + t279;
t302 = -t340 + t451;
t108 = t268 / 0.6e1 + t211 + t302;
t433 = 4 * pkin(1);
t468 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t55 = t182 * t155 + 0.6e1 * t100 * t464 + t136 * t468 + (t108 * t187 + t481) * t433;
t106 = t155 + t136;
t174 = t441 * t511;
t361 = 0.4e1 * t419;
t66 = t136 * t324 + (t252 + t353) * t277 + t457;
t56 = t106 * t361 + t174 * t197 + t66;
t81 = t136 * t467 + t121;
t39 = -0.8e1 * t118 * t370 + t141 + t81 * t385 + t66 * t360 + t274 + (-t265 + t268 + t449) * t285 + t185 * t160 + (0.12e2 * t55 * t200 + t240 + (t243 + t246 + 0.6e1 * t268) * t273 + t257 + (t244 + t247) * t279) * t277 + 0.6e1 * (-t393 * t84 - t56 * t485) * pkin(2);
t175 = pkin(2) * t512 + 0.8e1 * t273 * t284;
t367 = t288 * t467;
t101 = t175 * t233 + 0.4e1 * t232 * t367;
t251 = 5 * t271;
t438 = t278 + t285;
t448 = t242 + t262;
t461 = t279 * t273;
t117 = t448 * t277 + t251 + t438 + (6 * t461);
t255 = 0.5e1 * t285;
t259 = 6 * t279;
t126 = t255 + (t242 + t259) * t277 + t185;
t332 = 0.8e1 * t370;
t414 = -0.4e1 * t477;
t497 = pkin(1) * t234;
t253 = 3 * t273;
t178 = t253 + t440;
t142 = t178 * t188;
t491 = pkin(3) * t232;
t149 = t188 - t491;
t398 = t273 * t188;
t331 = t197 * t398;
t256 = 0.3e1 * t277;
t176 = t256 + t190;
t484 = t176 * t232;
t67 = -0.2e1 * t331 + t142 + (0.2e1 * t149 * t497 - t484) * pkin(3);
t303 = -t232 * t288 + t398;
t139 = 0.2e1 * t303;
t189 = (2 * t273) + t277;
t191 = -t273 + t279;
t70 = t189 * t491 + t139 * t197 + (t191 + t170) * t188;
t99 = -pkin(3) * t484 + t142;
t44 = t70 * t414 + t101 * t197 + (-0.4e1 * t99 * t497 + (t126 + t332) * t232) * pkin(3) + (0.4e1 * t67 * t485 + (-t117 + t346) * t233) * pkin(2);
t29 = t137 * t39 + t280 * t44;
t509 = 0.1e1 / t29 / 0.4e1;
t508 = -0.1e1 / t29 ^ 2 / 0.4e1;
t51 = 0.1e1 / t280;
t507 = t51 / 0.2e1;
t394 = pkin(3) * t151;
t82 = t155 + t172 + t350 + 0.2e1 * t394;
t79 = 0.1e1 / t82;
t506 = t79 / 0.2e1;
t80 = 0.1e1 / t82 ^ 2;
t505 = -t80 / 0.2e1;
t502 = 0.4e1 / 0.3e1 * t277;
t501 = -t280 / 0.4e1;
t500 = t280 / 0.4e1;
t473 = t234 * t233;
t399 = pkin(2) * t473;
t109 = t167 * t232 + t399;
t102 = t109 * t432;
t258 = 8 * t279;
t105 = -0.4e1 * t288 * t162 + t201 + t512 + (t244 + t258) * t273;
t212 = -t265 / 0.4e1;
t111 = t212 - t273 + t302;
t228 = t273 / 0.2e1;
t113 = -0.2e1 / 0.3e1 * t340 + t279 + t228 + t212;
t453 = t211 - t268 / 0.6e1;
t356 = t279 + t453;
t129 = t502 + t228 + t356;
t354 = t214 - t268 / 0.2e1 + t279;
t131 = 0.3e1 / 0.2e1 * t277 + t253 + t354;
t355 = t279 + t523;
t320 = t273 + t355;
t134 = t256 + t320;
t138 = 0.4e1 * t303;
t143 = 0.16e2 * (t438 - 0.6e1 * t463) * t271;
t312 = pkin(1) * t331;
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
t327 = -0.32e2 * t370;
t179 = (t260 + t265) * t273;
t391 = 0.16e2 * t481;
t465 = t271 * t195;
t409 = 0.8e1 * t465;
t147 = t279 + t277 / 0.4e1 + t273 / 0.4e1 - t265 / 0.8e1;
t450 = 0.4e1 / 0.7e1 * t279 - t265 / 0.7e1;
t64 = -0.32e2 / 0.21e2 * t147 * t340 + t285 / 0.7e1 + (0.16e2 / 0.21e2 * t273 + t450) * t277 + t271 / 0.7e1 + t450 * t273 + t278 - 0.3e1 / 0.7e1 * t462 + t444 / 0.42e2;
t148 = t212 + t227 + t451;
t226 = 0.4e1 / 0.3e1 * t273;
t65 = -0.8e1 / 0.3e1 * t148 * t340 + t285 / 0.3e1 + (t226 + t213) * t277 + t278 - t271 / 0.3e1 + (t502 + 0.2e1 / 0.3e1 * t273 + t215) * t279 + t444 / 0.18e2;
t45 = t184 * t409 + 0.14e2 * t64 * t464 + t182 * t307 + t191 * t285 + (t179 - 0.10e2 / 0.3e1 * t271 + (2 * t278) - t462) * t277 + t110 * t468 + (t113 * t391 + 0.6e1 * t65 * t187) * pkin(1);
t417 = 0.8e1 * t481;
t421 = 0.4e1 * t187;
t57 = t155 * t468 + t105 * t197 + t158 * t192 + (t111 * t421 + t417) * pkin(1);
t439 = t278 - t271;
t58 = t183 * t307 - t274 + (-t203 - t437) * t285 + (t179 + t439 + t521) * t277 + t279 * t110;
t250 = 7 * t271;
t63 = (t216 + t261 + (7 * t273)) * t285 + (t250 + (t245 + (10 * t279)) * t273 + t322) * t277 + t488;
t334 = pkin(1) * t391;
t75 = t166 * t234 + t181 * t197 + t334 + t409 + t442 - (6 * t461);
t87 = t155 * t467 + t158 * t193;
t26 = t57 * t327 + t143 * t195 + 0.24e2 * t58 * t464 + (t244 + t260 + (28 * t273)) * t274 + t160 * t483 + (t186 * t271 + 0.24e2 * t45 * t200 + t278 * t246 + t304 * t252 + t444 * t262 + (4 * t281 ^ 2) + (28 * t288 ^ 2)) * t277 + 0.8e1 * (-t393 * t63 - t46 * t485) * pkin(2) + (0.8e1 * t54 * t187 + 0.32e2 * t87 * t481) * pkin(1) + (t186 * t273 + 0.16e2 * t198 * t75 + (70 * t271) + t285 + t304) * t285;
t447 = t244 - 0.2e1 * t268;
t308 = 0.24e2 * t183 * t465 - t274 - ((21 * t273) + t522) * t285 - (t447 * t279 + t159 + t257 + (35 * t271)) * t277 - (t250 + (t258 + t446) * t273 + t279 * (-t268 + t437)) * t190;
t389 = -0.12e2 * t477;
t408 = -0.6e1 * t464;
t357 = t265 / 0.3e1 + t221 + t262;
t323 = -0.8e1 / 0.3e1 * t465 + t277 * t191 - 0.5e1 / 0.3e1 * t271 + t357 * t273 + t279 * t355;
t382 = -t491 / 0.2e1;
t482 = t196 * t271;
t321 = t230 + t356;
t130 = t226 + t321;
t76 = -t129 * t491 + t130 * t188;
t132 = t273 + t321;
t85 = t132 * t188 + t177 * t382;
t204 = -0.20e2 / 0.3e1 * t273;
t358 = 0.2e1 / 0.3e1 * t265 + t220 + t260;
t359 = 0.4e1 / 0.3e1 * t265 + 0.4e1 / 0.3e1 * t268 + t263;
t89 = -t285 + (t204 + t358) * t277 - (3 * t271) + t359 * t273 + t278;
t47 = t76 * t406 + t89 * t382 + t323 * t188 + (t85 * t187 - t232 * t482) * t433;
t397 = t288 * t188;
t330 = t196 * t397;
t407 = -0.4e1 * t464;
t351 = t259 + t447;
t225 = -0.2e1 / 0.3e1 * t268;
t452 = t215 + t225;
t77 = t285 + (t448 + t452) * t277 + t251 + t351 * t273 + t279 * (t279 + t452);
t92 = t255 + (t242 + t351) * t277 + t190 * (t225 + t352);
t59 = t77 * t188 - t92 * t491;
t135 = t178 + t523;
t78 = -t134 * t491 + t135 * t188;
t146 = t228 + t277 + t453;
t93 = t146 * t430 + t467 * t491;
t48 = t93 * t407 + (t421 * t78 - 0.8e1 * t330) * pkin(1) + t59;
t90 = -0.3e1 * t285 + (t204 + t359) * t277 + t358 * t273 + t439;
t94 = -0.5e1 / 0.3e1 * t285 + (-t273 + t357) * t277 + t279 * t320;
t60 = t90 * t188 + t423 * t94;
t61 = t192 * t188 + t138 * t197 + (t131 * t232 + t150 * t497) * t432;
t133 = t256 + 0.3e1 / 0.2e1 * t273 + t354;
t91 = t133 * t188 + t193 * t491 / 0.2e1;
t30 = t91 * t334 + t61 * t332 + t60 * t408 + t47 * t389 + (t526 * t232 - 0.6e1 * t59 * t497) * pkin(3) + (t233 * t308 + 0.6e1 * t48 * t485) * pkin(2);
t496 = pkin(1) * t271;
t311 = -0.64e2 * t467 * t496;
t313 = pkin(1) * t197 * t397;
t347 = 0.32e2 / 0.3e1 * t271;
t317 = t196 * t347;
t318 = 0.64e2 / 0.3e1 * t147 * t288;
t325 = t196 * t367;
t328 = -0.96e2 * t158 * t183 * t288;
t460 = t288 * t197;
t402 = pkin(1) * t460;
t333 = -0.48e2 * t402;
t495 = pkin(1) * t273;
t335 = -0.24e2 * t158 * t495;
t336 = -0.16e2 * t148 * t495;
t489 = t182 * pkin(3);
t342 = -0.4e1 * t158 * t489;
t343 = t468 * t517;
t344 = -0.48e2 * t86 * t495;
t403 = pkin(1) * t464;
t345 = 0.4e1 * t403;
t471 = t234 * t273;
t366 = t232 * t471;
t410 = 0.8e1 * t471;
t493 = pkin(3) * t167;
t424 = 0.4e1 * t493;
t53 = (t119 * t410 + t424 * t96) * t232 + (t116 * t421 + 0.8e1 * t167 * t464) * t188;
t374 = t53 * t507;
t384 = -0.24e2 * t460;
t386 = -0.24e2 * t471;
t387 = 0.24e2 * t477;
t390 = -0.32e2 * t479;
t392 = -0.32e2 * t482;
t405 = pkin(1) * t465;
t411 = -0.8e1 * t471;
t415 = 0.8e1 * t479;
t416 = -0.8e1 * t479;
t418 = -0.4e1 * t481;
t422 = 0.2e1 * t187;
t426 = 0.6e1 * t494;
t428 = -0.6e1 * t494;
t429 = -0.8e1 * t494;
t434 = -0.8e1 * t63 * pkin(3);
t435 = pkin(3) * t528;
t326 = t232 * t384;
t458 = pkin(1) * t193 * t326 - 0.24e2 * t183 * t330;
t459 = -4 * t510;
t15 = ((t131 * t422 + t345 + t418) * t332 + (0.12e2 * t194 * t197 * t496 + t129 * t418 - 0.2e1 * t177 * t403 - 0.4e1 * t405) * t389 + 0.8e1 * t193 * t405 + 0.12e2 * t94 * t481 + 0.6e1 * t92 * t403 + (-t89 * t389 / 0.2e1 + t526) * t187) * t280 + t30 * t374 + t137 * t144 * t327 + ((0.6e1 * (-0.4e1 * t134 * t403 - t92 * t187 - 0.4e1 * t325) * t280 + t458 * t515) * t485 + ((-pkin(1) * t195 * t347 - t196 * t318 + t197 * t336 + t234 * t342) * t387 + t195 * t311 + t196 * t328 + t197 * t344 + t234 * t434 + ((t234 * t343 + t418) * t390 + (t197 * t335 + t234 * t435) * t429) * t168) * t137 * t233) * pkin(2) + (0.12e2 * (-(t317 * t188 + t411 * t76) * t477 - 0.8e1 * t183 * t482 * t188 - 0.4e1 * t91 * t402 + t60 * t471) * t280 + t137 * (0.16e2 * (t181 * t514 - t166 + t333 + t392) * t480 + (t113 * t333 + t184 * t392 - 0.28e2 * t64 * t471) * t387 - 0.4e1 * t143 * t196 - 0.96e2 * t87 * t402 - 0.48e2 * t58 * t471) + (-t26 + t137 * (0.32e2 * t57 * t479 + t494 * t516) + (t48 * t428 + t61 * t416 + (-0.24e2 * t169 + 0.64e2 * t366) * t480) * t280 + ((0.48e2 * t85 * t477 + 0.6e1 * t59) * t280 + t137 * (-0.144e3 * t65 * t477 - 0.8e1 * t54)) * pkin(1)) * pkin(3) + ((0.2e1 * (-t138 * t234 - t150 * t510) * t415 + (t410 * t93 + t78 * t459 + 0.24e2 * t313) * t426) * t280 + t137 * ((pkin(1) * t384 + t105 * t514 + t111 * t459) * t390 + (t386 * t83 - 0.6e1 * t69 * t510) * t429)) * t168) * t232;
t24 = t137 * t26 + t280 * t30;
t455 = -t162 + t151;
t107 = pkin(3) + t455;
t115 = t172 + t268 + t319;
t395 = t119 * t187;
t71 = t115 * t167 + 0.2e1 * t395;
t72 = t115 * t234 + (0.4e1 * t197 - 0.2e1) * t493;
t43 = t107 * t280 + t72 * t188 + t232 * t71;
t380 = t43 * t509;
t152 = t253 + t268 + t348;
t154 = t162 - pkin(3);
t431 = pkin(1) * t513;
t487 = t109 * t280;
t97 = t152 + t172 + t324;
t42 = -t97 * t151 + t487 + (t152 * t475 + t154 * t431) * pkin(2) + (-t160 - t256 + t413 - 0.2e1 * t486) * pkin(3);
t301 = t24 * t380 + t42 * t500;
t466 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t383 = t79 * t466;
t22 = t301 * t383;
t338 = pkin(3) * t399;
t388 = 0.12e2 * t477;
t412 = 0.8e1 * t233 * t277;
t425 = t136 * t518;
t427 = 0.4e1 * t494;
t436 = t67 * t519;
t470 = t234 * t277;
t23 = (t273 * t194 * t416 + (t189 * t187 - 0.2e1 * t481) * t414 + 0.4e1 * t325 + t176 * t345 + t126 * t187) * t280 + t44 * t374 + t137 * ((-0.8e1 / 0.3e1 * t330 + t144 - 0.2e1 * t182 * t338) * t388 + t312 * t527 + t338 * t529 + t458) + ((0.8e1 * t139 * t200 * t470 + t101 * t514 - 0.24e2 * t313) * t280 + t137 * (0.12e2 * (-t100 * t471 - t402) * t388 + t81 * t386) + (t235 * t280 * t436 - t39 + t137 * (t118 * t415 + t426 * t56) + ((pkin(2) * t200 * t412 + 0.4e1 * t99) * t280 + t137 * (-0.48e2 * t108 * t477 - 0.6e1 * t66)) * pkin(1)) * pkin(3)) * t232 + ((t415 * t187 + ((-pkin(3) * t176 + t162 * t511) * t234 + (-t149 * t491 - t464) * t520) * t427) * t280 + t137 * ((t169 - 0.8e1 * t366) * t416 + ((-0.2e1 * t174 * t232 + t425 * t188) * t234 + (-0.4e1 * t106 * t491 - 0.8e1 * t331) * pkin(1)) * t428)) * t168;
t381 = t42 * t509;
t300 = t24 * t381 + t43 * t501;
t298 = t80 * t300;
t364 = t107 * t507;
t36 = -t487 + t53 * t364 + t119 * t194 * t517 + t71 * t234 + (-t115 - 0.8e1 * t394) * t162;
t363 = t109 * t507;
t37 = t455 * t280 + t53 * t363 + (t167 * t97 + 0.4e1 * t395) * t232 + (t431 * t470 + (t152 * t234 + t197 * t424) * pkin(2)) * t233;
t376 = -t43 * t51 / 0.8e1;
t379 = t42 * t508;
t499 = t22 - (t102 * t298 + (t15 * t381 + t36 * t501 + t53 * t376 + (t23 * t379 + t37 * t509) * t24) * t79) * t466;
t112 = t183 * t326 * t494;
t472 = t234 * t235;
t396 = pkin(2) * t472;
t341 = pkin(1) * t396;
t305 = t273 * t232 * t341;
t128 = -0.4e1 * t305;
t145 = t341 * t432;
t476 = t200 * t284;
t369 = t168 * t476;
t171 = pkin(1) * t430;
t469 = t235 * t277;
t368 = t233 * t469;
t123 = t171 - 0.4e1 * t368;
t401 = pkin(1) * t525;
t474 = t232 * t235;
t400 = pkin(2) * t474;
t339 = pkin(3) * t400;
t456 = -0.2e1 * t339 + t171;
t49 = t123 * t407 + (pkin(1) * t153 * t519 + t173 * t513) * t233 + (-t456 * t151 + 0.2e1 * t233 ^ 2 * t401 + (t116 * t474 - t96 * t473) * pkin(2)) * t530;
t375 = t49 * t507;
t478 = t199 * t285;
t16 = t30 * t375 + (-0.32e2 * t128 * t137 + 0.8e1 * t145 * t280) * t370 + (0.24e2 * (0.4e1 * t114 * t478 * t491 - t369 * t61 + t47 * t469) * t280 + t137 * (0.96e2 * t369 * t57 - 0.48e2 * t45 * t469 - 0.64e2 * t75 * t478)) * t233 + ((t26 + (t137 * t516 - 0.6e1 * t280 * t48) * t168) * t233 + (((t130 * t406 + t132 * t361 + t323) * t389 + t133 * t334 + t90 * t408 - 0.6e1 * t77 * t419 + t308) * t280 + (t112 * t515 + (t524 * t415 + (t135 * t361 - 0.8e1 * t146 * t464 - 0.8e1 * t404 + t77) * t426) * t280) * t168 + ((-pkin(1) * t317 - t197 * t318 + t234 * t336 + t342) * t387 + t196 * t311 + t197 * t328 + t234 * t344 + t434 + ((t343 - 0.4e1 * t460) * t390 + (t234 * t335 + t435) * t429) * t168) * t137 * t232) * t235) * pkin(2);
t25 = (t145 * t414 + (t70 * t412 + t175 * t197 + ((t191 + 0.2e1 * t464) * t414 - t117 + (-0.4e1 * t178 * t187 + t417) * pkin(1)) * pkin(2)) * t235 + ((t145 + (t178 - 0.2e1 * t464) * t494) * t427 + (-0.24e2 * t476 * t491 + t436) * t233) * t168) * t280 + t44 * t375 + t39 * t188 + t137 * (0.24e2 * t118 * t233 * t369 + (t128 + (-0.8e1 / 0.3e1 * t460 - 0.2e1 * t489) * t400) * t388 - 0.24e2 * t55 * t368 + t112 + t305 * t527 + t339 * t529 + 0.6e1 * (-(pkin(1) * t411 + t425) * t200 * t525 + t56 * t188) * t168);
t492 = pkin(3) * t197;
t34 = t49 * t364 + t123 * t232 * t422 + ((-t232 * t280 + t72) * t235 + (t234 * t280 + (t167 * t520 + t115) * t232 + (-pkin(3) + 0.2e1 * t492 + t497) * t430) * t233) * pkin(2);
t35 = (t162 + t396) * t280 + t49 * t363 - 0.2e1 * t123 * t492 - (t171 - 0.4e1 * t339) * t151 - 0.2e1 * t200 * t401 + t152 * t400 + (t469 * t518 + (t154 * t520 - t234 * t97) * pkin(2)) * t233;
t98 = 0.2e1 * t338 + t456;
t498 = t22 + (-t98 * t298 + (t16 * t381 + t34 * t501 + t49 * t376 + (t25 * t379 + t35 * t509) * t24) * t79) * t466;
t125 = -t473 + t474;
t124 = -t472 - t475;
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
t372 = pkin(2) * m(3) + mrSges(2,1);
t371 = pkin(3) * m(5) + mrSges(4,1);
t365 = t102 * t505;
t269 = 0.1e1 / pkin(4);
t362 = t269 * t506;
t41 = 0.1e1 / t42 ^ 2;
t337 = pkin(4) * t269 / (t41 * t43 ^ 2 + 0.1e1) * t82;
t316 = 0.1e1 / t42 * t337;
t9 = qJ(1) + atan2(t13, t14);
t7 = sin(t9);
t8 = cos(t9);
t315 = mrSges(3,1) * t8 - mrSges(3,2) * t7;
t314 = mrSges(3,1) * t7 + mrSges(3,2) * t8;
t33 = qJ(2) + atan2(t43 * t362, t42 * t362);
t31 = sin(t33);
t32 = cos(t33);
t310 = mrSges(5,1) * t32 - mrSges(5,2) * t31;
t309 = mrSges(5,1) * t31 + mrSges(5,2) * t32;
t306 = t41 * t43 * t337;
t299 = t80 * t301;
t18 = 0.1e1 + 0.2e1 * (t36 * t506 - t43 * t365) * t316 - 0.2e1 * (-t42 * t365 + t37 * t506) * t306;
t11 = 0.1e1 / t14;
t10 = 0.1e1 / (t12 * t13 ^ 2 + 0.1e1);
t6 = (-t98 * t299 + (t35 * t500 + t49 * t377 + t16 * t380 + (t25 * t378 + t34 * t509) * t24) * t79) * t466;
t4 = (t102 * t299 + (t37 * t500 + t53 * t377 + t15 * t380 + (t23 * t378 + t36 * t509) * t24) * t79) * t466;
t2 = 0.1e1 + ((t498 * t125 + (t21 - t6) * t124) * t11 - (-t498 * t124 - t125 * t6 + t19) * t490) * t10;
t1 = [(mrSges(2,2) * t233 + t315 * t2 - t235 * t372) * g(2) + (mrSges(2,2) * t235 - t314 * t2 + t233 * t372) * g(1) + 0.2e1 * (g(1) * t309 - g(2) * t310) * ((t34 * t506 + t43 * t373) * t316 - (t35 * t506 + t42 * t373) * t306), (mrSges(4,2) * t232 - t310 * t18 - t234 * t371) * g(2) + (mrSges(4,2) * t234 + t309 * t18 + t232 * t371) * g(1) + (-g(1) * t314 + g(2) * t315) * ((-t124 * t4 - t499 * t125 - t20) * t11 - ((-t21 - t4) * t125 + t499 * t124) * t490) * t10];
taug = t1(:);