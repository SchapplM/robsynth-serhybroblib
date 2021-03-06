% Calculate joint inertia matrix for
% fivebar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1DE2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1DE2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1DE2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 04:28:19
% EndTime: 2020-04-27 04:28:38
% DurationCPUTime: 18.46s
% Computational Cost: add. (289659->650), mult. (864196->1036), div. (2160->13), fcn. (105948->12), ass. (0->462)
t265 = pkin(5) ^ 2;
t264 = t265 ^ 2;
t268 = pkin(4) ^ 2;
t267 = t268 ^ 2;
t438 = t264 - t267;
t531 = 4 * pkin(3);
t521 = 2 * pkin(1);
t532 = 2 * pkin(3);
t273 = pkin(3) ^ 2;
t203 = 0.10e2 / 0.3e1 * t273;
t277 = pkin(2) ^ 2;
t215 = -0.2e1 / 0.3e1 * t265;
t220 = 0.2e1 / 0.3e1 * t268;
t279 = pkin(1) ^ 2;
t262 = 2 * t279;
t349 = t215 + t220 + t262;
t190 = t273 + t279;
t285 = t277 ^ 2;
t348 = t215 + t190;
t451 = t190 * (t220 + t348) + t285;
t84 = (t203 + t349) * t277 + t451;
t530 = -0.6e1 * t84;
t271 = t273 ^ 2;
t278 = t279 ^ 2;
t436 = t271 + t278;
t439 = t262 - t265;
t456 = t279 * t265;
t522 = -t264 / 0.6e1 + t267 / 0.6e1;
t110 = t439 * t273 + t436 - t456 - t522;
t297 = t110 + t285;
t86 = (t203 + t439) * t277 + t297;
t529 = -0.6e1 * t86;
t233 = sin(qJ(1));
t188 = pkin(2) * t233;
t424 = 0.2e1 * t188;
t511 = 4 * t273;
t252 = 6 * t273;
t213 = -t265 / 0.3e1;
t221 = t268 / 0.3e1;
t434 = t277 + t279;
t346 = t273 + t434;
t136 = t213 + t221 + t346;
t528 = -0.24e2 * t136;
t242 = 10 * t273;
t234 = cos(qJ(2));
t187 = pkin(3) * t234;
t413 = pkin(1) * t187;
t170 = 0.2e1 * t413;
t182 = -t273 / 0.3e1 + t279;
t197 = t234 ^ 2;
t458 = t273 * t197;
t114 = 0.4e1 / 0.3e1 * t458 + t170 + t182;
t261 = 3 * t279;
t523 = t261 - t265 - t268;
t159 = t523 * t242;
t246 = -0.6e1 * t265;
t284 = pkin(2) * t277;
t274 = t284 ^ 2;
t245 = -0.5e1 * t265;
t440 = t245 - 0.5e1 * t268;
t235 = cos(qJ(1));
t200 = t235 ^ 2;
t198 = t200 ^ 2;
t474 = t198 * t285;
t185 = t190 ^ 2;
t431 = t279 - t265;
t345 = t273 + t431;
t477 = t185 * (-t268 + t345);
t527 = 0.7e1 * t274 + ((35 * t273) + (15 * t279) + t440) * t285 + ((21 * t271) + t159 + (9 * t278) + (t246 - 0.6e1 * t268) * t279) * t277 + t477 - 0.24e2 * t114 * t474;
t232 = sin(qJ(2));
t526 = t232 * t277;
t192 = -3 * t273 + t279;
t400 = 0.4e1 * t458;
t525 = t192 + t400;
t524 = t213 - t268 / 0.3e1;
t520 = -0.4e1 * pkin(2);
t519 = -4 * pkin(3);
t518 = -2 * pkin(3);
t193 = -0.3e1 * t277 + t279;
t196 = t234 * t197;
t288 = pkin(3) * t273;
t475 = t196 * t288;
t398 = pkin(1) * t475;
t342 = 0.8e1 * t398;
t141 = t193 * t342;
t356 = 0.6e1 * t413;
t379 = 0.12e2 * t458;
t216 = -0.3e1 / 0.2e1 * t265;
t240 = 15 * t271;
t243 = 18 * t279;
t257 = 3 * t278;
t448 = t264 / 0.2e1 - t267 / 0.2e1;
t318 = -0.3e1 * t456 + t257 + t448;
t469 = t232 * t233;
t162 = pkin(2) * t469;
t336 = pkin(3) * t162;
t443 = 15 * t273 + t261;
t482 = t274 + t190 * ((t216 + t262) * t273 - 0.3e1 / 0.2e1 * t456 + t436 + t448);
t54 = t336 * t529 + (t240 + (t243 - 0.9e1 * t265) * t273 + t318) * t277 + (t216 + t443) * t285 + t482;
t214 = -t265 / 0.2e1;
t158 = t214 + t346;
t320 = -0.4e1 * t336;
t307 = t158 * t320;
t69 = t307 + (t252 + t439) * t277 + t297;
t155 = -0.2e1 * t336;
t183 = t279 - t277 / 0.3e1;
t121 = t183 * t155;
t461 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t83 = t158 * t461 + t121;
t46 = t356 * t69 + t379 * t83 + t141 + t54;
t516 = 0.8e1 * t46;
t492 = pkin(2) * t235;
t167 = pkin(1) - t492;
t137 = t167 + t187;
t515 = -0.8e1 * t137;
t514 = -0.2e1 * t234;
t513 = -0.2e1 * t235;
t512 = 4 * t271;
t510 = pkin(1) * pkin(3);
t344 = t277 + t431;
t315 = t273 + t344;
t153 = -t268 + t315;
t414 = pkin(1) * t492;
t172 = -0.2e1 * t414;
t116 = t172 + t153;
t151 = t167 * t234;
t201 = t277 * t511;
t457 = t277 * t279;
t173 = t201 - 0.4e1 * t457;
t247 = 0.2e1 * t268;
t263 = -2 * t279;
t471 = t200 * t277;
t407 = 0.2e1 * t471;
t435 = -t277 + t279;
t119 = t172 + t407 + t435;
t480 = t119 * t197;
t503 = pkin(5) - pkin(4);
t504 = -pkin(4) - pkin(5);
t96 = t155 + t116;
t280 = sqrt(t173 * t200 + 0.4e1 * t153 * t414 - t271 - (t279 + (pkin(2) - t503) * (pkin(2) + t503)) * (t279 + (pkin(2) - t504) * (pkin(2) + t504)) + (t247 + t263 + 0.2e1 * t265 - 0.6e1 * t277 - 0.4e1 * t480) * t273 + (t116 * t162 - t96 * t151) * t531);
t118 = t170 + t525;
t160 = t268 + t345;
t244 = -0.2e1 * t265;
t168 = pkin(1) + t187;
t199 = t235 * t200;
t473 = t199 * t284;
t366 = t168 * t473;
t387 = pkin(3) * t469;
t479 = t168 * t235;
t227 = t273 / 0.3e1;
t100 = -0.4e1 / 0.9e1 * t336 + t279 + t277 / 0.3e1 + t227 + t268 / 0.9e1 - t265 / 0.9e1;
t211 = -t265 / 0.6e1;
t230 = t277 / 0.2e1;
t445 = t230 + t279;
t302 = -t336 + t445;
t108 = t268 / 0.6e1 + t211 + t302;
t427 = 4 * pkin(1);
t462 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t55 = t182 * t155 + 0.6e1 * t100 * t458 + t136 * t462 + (t108 * t187 + t475) * t427;
t106 = t155 + t136;
t174 = t435 * t511;
t357 = 0.4e1 * t413;
t66 = t136 * t320 + (t252 + t349) * t277 + t451;
t56 = t106 * t357 + t174 * t197 + t66;
t81 = t136 * t461 + t121;
t39 = -0.8e1 * t118 * t366 + t141 + t81 * t379 + t66 * t356 + t274 + (-t265 + t268 + t443) * t285 + t185 * t160 + (0.12e2 * t55 * t200 + t240 + (t243 + t246 + 0.6e1 * t268) * t273 + t257 + (t244 + t247) * t279) * t277 + 0.6e1 * (-t387 * t84 - t56 * t479) * pkin(2);
t175 = pkin(2) * t512 + 0.8e1 * t273 * t284;
t363 = t288 * t461;
t101 = t175 * t233 + 0.4e1 * t232 * t363;
t251 = 5 * t271;
t432 = t278 + t285;
t442 = t242 + t262;
t455 = t279 * t273;
t117 = t442 * t277 + t251 + t432 + (6 * t455);
t255 = 0.5e1 * t285;
t259 = 6 * t279;
t126 = t255 + (t242 + t259) * t277 + t185;
t328 = 0.8e1 * t366;
t408 = -0.4e1 * t471;
t495 = pkin(1) * t234;
t253 = 3 * t273;
t178 = t253 + t434;
t142 = t178 * t188;
t489 = pkin(3) * t232;
t149 = t188 - t489;
t392 = t273 * t188;
t327 = t197 * t392;
t256 = 0.3e1 * t277;
t176 = t256 + t190;
t478 = t176 * t232;
t67 = -0.2e1 * t327 + t142 + (0.2e1 * t149 * t495 - t478) * pkin(3);
t303 = -t232 * t288 + t392;
t139 = 0.2e1 * t303;
t189 = (2 * t273) + t277;
t191 = -t273 + t279;
t70 = t189 * t489 + t139 * t197 + (t191 + t170) * t188;
t99 = -pkin(3) * t478 + t142;
t44 = t70 * t408 + t101 * t197 + (-0.4e1 * t99 * t495 + (t126 + t328) * t232) * pkin(3) + (0.4e1 * t67 * t479 + (-t117 + t342) * t233) * pkin(2);
t29 = t137 * t39 + t280 * t44;
t509 = 0.1e1 / t29 / 0.4e1;
t508 = -0.1e1 / t29 ^ 2 / 0.4e1;
t51 = 0.1e1 / t280;
t507 = t51 / 0.2e1;
t388 = pkin(3) * t151;
t82 = t155 + t172 + t346 + 0.2e1 * t388;
t79 = 0.1e1 / t82;
t506 = t79 / 0.2e1;
t80 = 0.1e1 / t82 ^ 2;
t505 = -t80 / 0.2e1;
t466 = t234 * t235;
t124 = -t466 - t469;
t467 = t234 * t233;
t468 = t232 * t235;
t125 = -t467 + t468;
t143 = 0.16e2 * (t432 - 0.6e1 * t457) * t271;
t186 = -0.30e2 * t265 + (60 * t279);
t195 = t197 ^ 2;
t260 = 4 * t279;
t281 = pkin(1) * t279;
t304 = (6 * t278) + t438 - 0.6e1 * t456;
t323 = -0.32e2 * t366;
t212 = -t265 / 0.4e1;
t228 = t273 / 0.2e1;
t113 = -0.2e1 / 0.3e1 * t336 + t279 + t228 + t212;
t179 = (t260 + t265) * t273;
t184 = t279 - 0.2e1 / 0.3e1 * t277;
t385 = 0.16e2 * t475;
t459 = t271 * t195;
t403 = 0.8e1 * t459;
t147 = t279 + t277 / 0.4e1 + t273 / 0.4e1 - t265 / 0.8e1;
t444 = 0.4e1 / 0.7e1 * t279 - t265 / 0.7e1;
t64 = -0.32e2 / 0.21e2 * t147 * t336 + t285 / 0.7e1 + (0.16e2 / 0.21e2 * t273 + t444) * t277 + t271 / 0.7e1 + t444 * t273 + t278 - 0.3e1 / 0.7e1 * t456 + t438 / 0.42e2;
t148 = t212 + t227 + t445;
t226 = 0.4e1 / 0.3e1 * t273;
t500 = 0.4e1 / 0.3e1 * t277;
t65 = -0.8e1 / 0.3e1 * t148 * t336 + t285 / 0.3e1 + (t226 + t213) * t277 + t278 - t271 / 0.3e1 + (t500 + 0.2e1 / 0.3e1 * t273 + t215) * t279 + t438 / 0.18e2;
t45 = t184 * t403 + 0.14e2 * t64 * t458 + t182 * t307 + t191 * t285 + (t179 - 0.10e2 / 0.3e1 * t271 + (2 * t278) - t456) * t277 + t110 * t462 + (t113 * t385 + 0.6e1 * t65 * t187) * pkin(1);
t258 = 8 * t279;
t105 = -0.4e1 * t288 * t162 + t201 + t512 + (t244 + t258) * t273;
t111 = t212 - t273 + t302;
t411 = 0.8e1 * t475;
t415 = 0.4e1 * t187;
t57 = t155 * t462 + t105 * t197 + t158 * t192 + (t111 * t415 + t411) * pkin(1);
t433 = t278 - t271;
t58 = t183 * t307 - t274 + (-t203 - t431) * t285 + (t179 + t433 + t522) * t277 + t279 * t110;
t250 = 7 * t271;
t63 = (t216 + t261 + (7 * t273)) * t285 + (t250 + (t245 + (10 * t279)) * t273 + t318) * t277 + t482;
t166 = -12 * pkin(1) * t288 + t281 * t531;
t181 = -8 * t271 + 12 * t455;
t330 = pkin(1) * t385;
t75 = t166 * t234 + t181 * t197 + t330 + t403 + t436 - (6 * t455);
t87 = t155 * t461 + t158 * t193;
t26 = t57 * t323 + t143 * t195 + 0.24e2 * t58 * t458 + (t244 + t260 + (28 * t273)) * t274 + t160 * t477 + (t186 * t271 + 0.24e2 * t45 * t200 + t278 * t246 + t304 * t252 + t438 * t262 + (4 * t281 ^ 2) + (28 * t288 ^ 2)) * t277 + 0.8e1 * (-t387 * t63 - t46 * t479) * pkin(2) + (0.8e1 * t54 * t187 + 0.32e2 * t87 * t475) * pkin(1) + (t186 * t273 + 0.16e2 * t198 * t75 + (70 * t271) + t285 + t304) * t285;
t441 = t244 - 0.2e1 * t268;
t308 = 0.24e2 * t183 * t459 - t274 - ((21 * t273) + t523) * t285 - (t441 * t279 + t159 + t257 + (35 * t271)) * t277 - (t250 + (t258 + t440) * t273 + t279 * (-t268 + t431)) * t190;
t383 = -0.12e2 * t471;
t402 = -0.6e1 * t458;
t351 = t279 + t524;
t353 = t265 / 0.3e1 + t221 + t262;
t319 = -0.8e1 / 0.3e1 * t459 + t277 * t191 - 0.5e1 / 0.3e1 * t271 + t353 * t273 + t279 * t351;
t376 = -t489 / 0.2e1;
t476 = t196 * t271;
t447 = t211 - t268 / 0.6e1;
t352 = t279 + t447;
t129 = t500 + t228 + t352;
t317 = t230 + t352;
t130 = t226 + t317;
t76 = -t129 * t489 + t130 * t188;
t132 = t273 + t317;
t177 = 0.2e1 * t277 + t191;
t85 = t132 * t188 + t177 * t376;
t204 = -0.20e2 / 0.3e1 * t273;
t354 = 0.2e1 / 0.3e1 * t265 + t220 + t260;
t355 = 0.4e1 / 0.3e1 * t265 + 0.4e1 / 0.3e1 * t268 + t263;
t89 = -t285 + (t204 + t354) * t277 - (3 * t271) + t355 * t273 + t278;
t47 = t76 * t400 + t89 * t376 + t319 * t188 + (t85 * t187 - t232 * t476) * t427;
t391 = t288 * t188;
t326 = t196 * t391;
t401 = -0.4e1 * t458;
t347 = t259 + t441;
t225 = -0.2e1 / 0.3e1 * t268;
t446 = t215 + t225;
t77 = t285 + (t442 + t446) * t277 + t251 + t347 * t273 + t279 * (t279 + t446);
t92 = t255 + (t242 + t347) * t277 + t190 * (t225 + t348);
t59 = t77 * t188 - t92 * t489;
t316 = t273 + t351;
t134 = t256 + t316;
t135 = t178 + t524;
t78 = -t134 * t489 + t135 * t188;
t146 = t228 + t277 + t447;
t93 = t146 * t424 + t461 * t489;
t48 = t93 * t401 + (t415 * t78 - 0.8e1 * t326) * pkin(1) + t59;
t417 = -0.2e1 * t489;
t90 = -0.3e1 * t285 + (t204 + t355) * t277 + t354 * t273 + t433;
t94 = -0.5e1 / 0.3e1 * t285 + (-t273 + t353) * t277 + t279 * t316;
t60 = t90 * t188 + t417 * t94;
t350 = t214 - t268 / 0.2e1 + t279;
t131 = 0.3e1 / 0.2e1 * t277 + t253 + t350;
t138 = 0.4e1 * t303;
t150 = t188 + 0.2e1 * t489;
t61 = t192 * t188 + t138 * t197 + (t131 * t232 + t150 * t495) * t532;
t133 = t256 + 0.3e1 / 0.2e1 * t273 + t350;
t91 = t133 * t188 + t193 * t489 / 0.2e1;
t30 = t91 * t330 + t61 * t328 + t60 * t402 + t47 * t383 + (t527 * t232 - 0.6e1 * t59 * t495) * pkin(3) + (t233 * t308 + 0.6e1 * t48 * t479) * pkin(2);
t24 = t137 * t26 + t280 * t30;
t152 = t253 + t268 + t344;
t154 = t162 - pkin(3);
t425 = pkin(1) * t513;
t393 = pkin(2) * t467;
t109 = t167 * t232 + t393;
t481 = t109 * t280;
t97 = t152 + t172 + t320;
t42 = -t97 * t151 + t481 + (t152 * t469 + t154 * t425) * pkin(2) + (-t160 - t256 + t407 - 0.2e1 * t480) * pkin(3);
t375 = t42 * t509;
t449 = -t162 + t151;
t107 = pkin(3) + t449;
t115 = t172 + t268 + t315;
t389 = t119 * t187;
t71 = t115 * t167 + 0.2e1 * t389;
t491 = pkin(3) * t167;
t72 = t115 * t234 + (0.4e1 * t197 - 0.2e1) * t491;
t43 = t107 * t280 + t72 * t188 + t232 * t71;
t499 = -t280 / 0.4e1;
t300 = t24 * t375 + t43 * t499;
t460 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t377 = t79 * t460;
t21 = t300 * t377;
t19 = t21 * t125;
t374 = t43 * t509;
t498 = t280 / 0.4e1;
t301 = t24 * t374 + t42 * t498;
t22 = t301 * t377;
t13 = -t124 * t22 + t19;
t20 = t21 * t124;
t14 = -t125 * t22 - t20;
t10 = atan2(t13, t14);
t8 = cos(t10);
t502 = mrSges(3,1) * t8;
t11 = 0.1e1 / t14;
t12 = 0.1e1 / t14 ^ 2;
t484 = t12 * t13;
t454 = t288 * t197;
t378 = -0.24e2 * t454;
t322 = t232 * t378;
t112 = t183 * t322 * t492;
t390 = pkin(2) * t466;
t337 = pkin(1) * t390;
t305 = t273 * t232 * t337;
t128 = -0.4e1 * t305;
t145 = t337 * t532;
t494 = pkin(1) * t271;
t309 = -0.64e2 * t461 * t494;
t343 = 0.32e2 / 0.3e1 * t271;
t313 = t196 * t343;
t314 = 0.64e2 / 0.3e1 * t147 * t288;
t324 = -0.96e2 * t158 * t183 * t288;
t493 = pkin(1) * t273;
t331 = -0.24e2 * t158 * t493;
t332 = -0.16e2 * t148 * t493;
t483 = t182 * pkin(3);
t338 = -0.4e1 * t158 * t483;
t339 = t462 * t518;
t340 = -0.48e2 * t86 * t493;
t470 = t200 * t284;
t365 = t168 * t470;
t171 = pkin(1) * t424;
t463 = t235 * t277;
t364 = t233 * t463;
t123 = t171 - 0.4e1 * t364;
t395 = pkin(1) * t526;
t394 = pkin(2) * t468;
t335 = pkin(3) * t394;
t450 = -0.2e1 * t335 + t171;
t49 = t123 * t401 + (pkin(1) * t153 * t520 + t173 * t513) * t233 + (-t450 * t151 + 0.2e1 * t233 ^ 2 * t395 + (t116 * t468 - t96 * t467) * pkin(2)) * t531;
t369 = t49 * t507;
t381 = 0.24e2 * t471;
t384 = -0.32e2 * t473;
t409 = 0.8e1 * t473;
t420 = 0.6e1 * t492;
t423 = -0.8e1 * t492;
t428 = -0.8e1 * t63 * pkin(3);
t429 = pkin(3) * t529;
t472 = t199 * t285;
t16 = t30 * t369 + (-0.32e2 * t128 * t137 + 0.8e1 * t145 * t280) * t366 + (0.24e2 * (0.4e1 * t114 * t472 * t489 - t365 * t61 + t47 * t463) * t280 + t137 * (0.96e2 * t365 * t57 - 0.48e2 * t45 * t463 - 0.64e2 * t75 * t472)) * t233 + ((t26 + (t137 * t516 - 0.6e1 * t280 * t48) * t168) * t233 + (((t130 * t400 + t132 * t357 + t319) * t383 + t133 * t330 + t90 * t402 - 0.6e1 * t77 * t413 + t308) * t280 + (t112 * t515 + (t525 * t409 + (t135 * t357 - 0.8e1 * t146 * t458 - 0.8e1 * t398 + t77) * t420) * t280) * t168 + ((-pkin(1) * t313 - t197 * t314 + t234 * t332 + t338) * t381 + t196 * t309 + t197 * t324 + t234 * t340 + t428 + ((t339 - 0.4e1 * t454) * t384 + (t234 * t331 + t429) * t423) * t168) * t137 * t232) * t235) * pkin(2);
t382 = 0.12e2 * t471;
t465 = t234 * t273;
t405 = -0.8e1 * t465;
t406 = 0.8e1 * t233 * t277;
t419 = t136 * t519;
t421 = 0.4e1 * t492;
t430 = t67 * t520;
t25 = (t145 * t408 + (t70 * t406 + t175 * t197 + ((t191 + 0.2e1 * t458) * t408 - t117 + (-0.4e1 * t178 * t187 + t411) * pkin(1)) * pkin(2)) * t235 + ((t145 + (t178 - 0.2e1 * t458) * t492) * t421 + (-0.24e2 * t470 * t489 + t430) * t233) * t168) * t280 + t44 * t369 + t39 * t188 + t137 * (0.24e2 * t118 * t233 * t365 + (t128 + (-0.8e1 / 0.3e1 * t454 - 0.2e1 * t483) * t394) * t382 - 0.24e2 * t55 * t364 + t112 + t305 * t528 + t335 * t530 + 0.6e1 * (-(pkin(1) * t405 + t419) * t200 * t526 + t56 * t188) * t168);
t298 = t80 * t300;
t360 = t107 * t507;
t416 = 0.2e1 * t187;
t490 = pkin(3) * t197;
t33 = t49 * t360 + t123 * t232 * t416 + ((-t232 * t280 + t72) * t235 + (t234 * t280 + (t167 * t521 + t115) * t232 + (-pkin(3) + 0.2e1 * t490 + t495) * t424) * t233) * pkin(2);
t359 = t109 * t507;
t35 = (t162 + t390) * t280 + t49 * t359 - 0.2e1 * t123 * t490 - (t171 - 0.4e1 * t335) * t151 - 0.2e1 * t200 * t395 + t152 * t394 + (t463 * t519 + (t154 * t521 - t234 * t97) * pkin(2)) * t233;
t370 = -t43 * t51 / 0.8e1;
t373 = t42 * t508;
t334 = pkin(3) * t393;
t98 = 0.2e1 * t334 + t450;
t496 = t22 + (-t98 * t298 + (t16 * t375 + t33 * t499 + t49 * t370 + (t25 * t373 + t35 * t509) * t24) * t79) * t460;
t299 = t80 * t301;
t371 = t42 * t51 / 0.8e1;
t372 = t43 * t508;
t6 = (-t98 * t299 + (t35 * t498 + t49 * t371 + t16 * t374 + (t25 * t372 + t33 * t509) * t24) * t79) * t460;
t9 = 0.1e1 / (t12 * t13 ^ 2 + 0.1e1);
t2 = 0.1e1 + ((t496 * t125 + (t21 - t6) * t124) * t11 - (-t496 * t124 - t125 * t6 + t19) * t484) * t9;
t501 = Ifges(3,3) * t2;
t102 = t109 * t532;
t310 = pkin(1) * t327;
t144 = -0.4e1 * t310;
t169 = pkin(1) * t417;
t194 = t232 ^ 2;
t311 = pkin(1) * t197 * t391;
t321 = t196 * t363;
t396 = pkin(1) * t454;
t329 = -0.48e2 * t396;
t397 = pkin(1) * t458;
t341 = 0.4e1 * t397;
t362 = t232 * t465;
t404 = 0.8e1 * t465;
t418 = 0.4e1 * t491;
t53 = (t119 * t404 + t418 * t96) * t232 + (t116 * t415 + 0.8e1 * t167 * t458) * t188;
t368 = t53 * t507;
t380 = -0.24e2 * t465;
t386 = -0.32e2 * t476;
t399 = pkin(1) * t459;
t410 = -0.8e1 * t473;
t412 = -0.4e1 * t475;
t422 = -0.6e1 * t492;
t452 = pkin(1) * t193 * t322 - 0.24e2 * t183 * t326;
t453 = -4 * t510;
t15 = ((t131 * t416 + t341 + t412) * t328 + (0.12e2 * t194 * t197 * t494 + t129 * t412 - 0.2e1 * t177 * t397 - 0.4e1 * t399) * t383 + 0.8e1 * t193 * t399 + 0.12e2 * t94 * t475 + 0.6e1 * t92 * t397 + (-t89 * t383 / 0.2e1 + t527) * t187) * t280 + t30 * t368 + t137 * t144 * t323 + ((0.6e1 * (-0.4e1 * t134 * t397 - t92 * t187 - 0.4e1 * t321) * t280 + t452 * t515) * t479 + ((-pkin(1) * t195 * t343 - t196 * t314 + t197 * t332 + t234 * t338) * t381 + t195 * t309 + t196 * t324 + t197 * t340 + t234 * t428 + ((t234 * t339 + t412) * t384 + (t197 * t331 + t234 * t429) * t423) * t168) * t137 * t233) * pkin(2) + (0.12e2 * (-(t313 * t188 + t405 * t76) * t471 - 0.8e1 * t183 * t476 * t188 - 0.4e1 * t91 * t396 + t60 * t465) * t280 + t137 * (0.16e2 * (t181 * t514 - t166 + t329 + t386) * t474 + (t113 * t329 + t184 * t386 - 0.28e2 * t64 * t465) * t381 - 0.4e1 * t143 * t196 - 0.96e2 * t87 * t396 - 0.48e2 * t58 * t465) + (-t26 + t137 * (0.32e2 * t57 * t473 + t492 * t516) + (t48 * t422 + t61 * t410 + (-0.24e2 * t169 + 0.64e2 * t362) * t474) * t280 + ((0.48e2 * t85 * t471 + 0.6e1 * t59) * t280 + t137 * (-0.144e3 * t65 * t471 - 0.8e1 * t54)) * pkin(1)) * pkin(3) + ((0.2e1 * (-t138 * t234 - t150 * t510) * t409 + (t404 * t93 + t78 * t453 + 0.24e2 * t311) * t420) * t280 + t137 * ((pkin(1) * t378 + t105 * t514 + t111 * t453) * t384 + (t380 * t83 - 0.6e1 * t69 * t510) * t423)) * t168) * t232;
t464 = t234 * t277;
t23 = (t273 * t194 * t410 + (t189 * t187 - 0.2e1 * t475) * t408 + 0.4e1 * t321 + t176 * t341 + t126 * t187) * t280 + t44 * t368 + t137 * ((-0.8e1 / 0.3e1 * t326 + t144 - 0.2e1 * t182 * t334) * t382 + t310 * t528 + t334 * t530 + t452) + ((0.8e1 * t139 * t200 * t464 + t101 * t514 - 0.24e2 * t311) * t280 + t137 * (0.12e2 * (-t100 * t465 - t396) * t382 + t81 * t380) + (t235 * t280 * t430 - t39 + t137 * (t118 * t409 + t420 * t56) + ((pkin(2) * t200 * t406 + 0.4e1 * t99) * t280 + t137 * (-0.48e2 * t108 * t471 - 0.6e1 * t66)) * pkin(1)) * pkin(3)) * t232 + ((t409 * t187 + ((-pkin(3) * t176 + t162 * t511) * t234 + (-t149 * t489 - t458) * t521) * t421) * t280 + t137 * ((t169 - 0.8e1 * t362) * t410 + ((-0.2e1 * t174 * t232 + t419 * t188) * t234 + (-0.4e1 * t106 * t489 - 0.8e1 * t327) * pkin(1)) * t422)) * t168;
t36 = -t481 + t53 * t360 + t119 * t194 * t518 + t71 * t234 + (-t115 - 0.8e1 * t388) * t162;
t37 = t449 * t280 + t53 * t359 + (t167 * t97 + 0.4e1 * t389) * t232 + (t425 * t464 + (t152 * t234 + t197 * t418) * pkin(2)) * t233;
t497 = t22 - (t102 * t298 + (t15 * t375 + t36 * t499 + t53 * t370 + (t23 * t373 + t37 * t509) * t24) * t79) * t460;
t488 = m(3) * t277;
t487 = m(5) * t273;
t269 = 0.1e1 / pkin(4);
t358 = t269 * t506;
t34 = atan2(t43 * t358, t42 * t358);
t32 = cos(t34);
t486 = mrSges(5,1) * t32;
t41 = 0.1e1 / t42 ^ 2;
t333 = pkin(4) * t269 / (t41 * t43 ^ 2 + 0.1e1) * t82;
t306 = t41 * t43 * t333;
t312 = 0.1e1 / t42 * t333;
t361 = t102 * t505;
t18 = 0.1e1 + 0.2e1 * (t36 * t506 - t43 * t361) * t312 - 0.2e1 * (-t42 * t361 + t37 * t506) * t306;
t485 = Ifges(5,3) * t18;
t367 = t98 * t505;
t31 = sin(t34);
t17 = 0.2e1 * (t33 * t506 + t43 * t367) * t312 - 0.2e1 * (t35 * t506 + t42 * t367) * t306;
t7 = sin(t10);
t4 = (t102 * t299 + (t37 * t498 + t53 * t371 + t15 * t374 + (t23 * t372 + t36 * t509) * t24) * t79) * t460;
t1 = ((-t124 * t4 - t497 * t125 - t20) * t11 - ((-t21 - t4) * t125 + t497 * t124) * t484) * t9;
t3 = [t8 ^ 2 * t488 + t17 ^ 2 * Ifges(5,3) + Ifges(2,3) + (-0.2e1 * pkin(2) * t502 + t501) * t2 + (0.2e1 * pkin(2) * mrSges(3,2) * t2 + t7 * t488) * t7; (t485 + (-mrSges(5,2) * t31 + t486) * pkin(3)) * t17 + (t501 + (mrSges(3,2) * t7 - t502) * pkin(2)) * t1; t32 ^ 2 * t487 + t1 ^ 2 * Ifges(3,3) + Ifges(4,3) + (t486 * t532 + t485) * t18 + (mrSges(5,2) * t18 * t518 + t31 * t487) * t31;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t3(1), t3(2); t3(2), t3(3);];
Mq = res;
