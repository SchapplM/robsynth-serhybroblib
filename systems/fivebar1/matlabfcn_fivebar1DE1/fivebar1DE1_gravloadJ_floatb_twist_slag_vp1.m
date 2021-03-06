% Calculate Gravitation load on the joints for
% fivebar1DE1
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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1DE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1DE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:25:08
% EndTime: 2020-04-27 02:27:59
% DurationCPUTime: 31.60s
% Computational Cost: add. (621770->667), mult. (1853030->1066), div. (4704->13), fcn. (226908->12), ass. (0->470)
t273 = pkin(5) ^ 2;
t272 = t273 ^ 2;
t276 = pkin(4) ^ 2;
t275 = t276 ^ 2;
t448 = t272 - t275;
t540 = 4 * pkin(3);
t530 = 2 * pkin(1);
t436 = 2 * pkin(3);
t281 = pkin(3) ^ 2;
t211 = 0.10e2 / 0.3e1 * t281;
t285 = pkin(2) ^ 2;
t223 = -0.2e1 / 0.3e1 * t273;
t228 = 0.2e1 / 0.3e1 * t276;
t287 = pkin(1) ^ 2;
t270 = 2 * t287;
t357 = t223 + t228 + t270;
t198 = t281 + t287;
t293 = t285 ^ 2;
t356 = t223 + t198;
t461 = t198 * (t228 + t356) + t293;
t92 = (t211 + t357) * t285 + t461;
t539 = -0.6e1 * t92;
t279 = t281 ^ 2;
t286 = t287 ^ 2;
t446 = t279 + t286;
t449 = t270 - t273;
t467 = t287 * t273;
t531 = -t272 / 0.6e1 + t275 / 0.6e1;
t118 = t449 * t281 + t446 - t467 - t531;
t305 = t118 + t293;
t94 = (t211 + t449) * t285 + t305;
t538 = -0.6e1 * t94;
t241 = sin(qJ(1));
t196 = pkin(2) * t241;
t434 = 0.2e1 * t196;
t521 = 4 * t281;
t260 = 6 * t281;
t221 = -t273 / 0.3e1;
t229 = t276 / 0.3e1;
t444 = t285 + t287;
t354 = t281 + t444;
t144 = t221 + t229 + t354;
t537 = -0.24e2 * t144;
t250 = 10 * t281;
t242 = cos(qJ(2));
t195 = pkin(3) * t242;
t423 = pkin(1) * t195;
t178 = 0.2e1 * t423;
t190 = -t281 / 0.3e1 + t287;
t205 = t242 ^ 2;
t469 = t281 * t205;
t122 = 0.4e1 / 0.3e1 * t469 + t178 + t190;
t269 = 3 * t287;
t532 = t269 - t273 - t276;
t167 = t532 * t250;
t254 = -0.6e1 * t273;
t292 = pkin(2) * t285;
t282 = t292 ^ 2;
t253 = -0.5e1 * t273;
t450 = t253 - 0.5e1 * t276;
t243 = cos(qJ(1));
t208 = t243 ^ 2;
t206 = t208 ^ 2;
t486 = t206 * t293;
t193 = t198 ^ 2;
t441 = t287 - t273;
t353 = t281 + t441;
t488 = t193 * (-t276 + t353);
t536 = 0.7e1 * t282 + ((35 * t281) + (15 * t287) + t450) * t293 + ((21 * t279) + t167 + (9 * t286) + (t254 - 0.6e1 * t276) * t287) * t285 + t488 - 0.24e2 * t122 * t486;
t240 = sin(qJ(2));
t535 = t240 * t285;
t200 = -3 * t281 + t287;
t410 = 0.4e1 * t469;
t534 = t200 + t410;
t533 = t221 - t276 / 0.3e1;
t529 = -0.4e1 * pkin(2);
t528 = -4 * pkin(3);
t527 = -2 * pkin(3);
t201 = -0.3e1 * t285 + t287;
t204 = t242 * t205;
t296 = pkin(3) * t281;
t487 = t204 * t296;
t408 = pkin(1) * t487;
t350 = 0.8e1 * t408;
t149 = t201 * t350;
t364 = 0.6e1 * t423;
t389 = 0.12e2 * t469;
t224 = -0.3e1 / 0.2e1 * t273;
t248 = 15 * t279;
t251 = 18 * t287;
t265 = 3 * t286;
t458 = t272 / 0.2e1 - t275 / 0.2e1;
t326 = -0.3e1 * t467 + t265 + t458;
t481 = t240 * t241;
t170 = pkin(2) * t481;
t344 = pkin(3) * t170;
t453 = 15 * t281 + t269;
t463 = t198 * ((t224 + t270) * t281 - 0.3e1 / 0.2e1 * t467 + t446 + t458) + t282;
t62 = t344 * t538 + (t248 + (t251 - 0.9e1 * t273) * t281 + t326) * t285 + (t224 + t453) * t293 + t463;
t222 = -t273 / 0.2e1;
t166 = t222 + t354;
t328 = -0.4e1 * t344;
t315 = t166 * t328;
t77 = t315 + (t260 + t449) * t285 + t305;
t163 = -0.2e1 * t344;
t191 = t287 - t285 / 0.3e1;
t129 = t191 * t163;
t473 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t91 = t166 * t473 + t129;
t54 = t364 * t77 + t389 * t91 + t149 + t62;
t526 = 0.8e1 * t54;
t502 = pkin(2) * t243;
t175 = pkin(1) - t502;
t145 = t175 + t195;
t525 = -0.8e1 * t145;
t524 = -0.2e1 * t242;
t523 = -0.2e1 * t243;
t522 = 4 * t279;
t479 = t242 * t241;
t403 = pkin(2) * t479;
t342 = pkin(3) * t403;
t179 = pkin(1) * t434;
t480 = t240 * t243;
t404 = pkin(2) * t480;
t343 = pkin(3) * t404;
t460 = -0.2e1 * t343 + t179;
t106 = 0.2e1 * t342 + t460;
t465 = t296 * t205;
t388 = -0.24e2 * t465;
t330 = t240 * t388;
t120 = t191 * t330 * t502;
t478 = t242 * t243;
t400 = pkin(2) * t478;
t345 = pkin(1) * t400;
t313 = t281 * t240 * t345;
t136 = -0.4e1 * t313;
t234 = 0.4e1 / 0.3e1 * t281;
t238 = t285 / 0.2e1;
t219 = -t273 / 0.6e1;
t457 = t219 - t276 / 0.6e1;
t360 = t287 + t457;
t325 = t238 + t360;
t138 = t234 + t325;
t140 = t281 + t325;
t264 = 0.3e1 * t285;
t358 = t222 - t276 / 0.2e1 + t287;
t141 = t264 + 0.3e1 / 0.2e1 * t281 + t358;
t261 = 3 * t281;
t186 = t261 + t444;
t143 = t186 + t533;
t153 = t345 * t436;
t236 = t281 / 0.2e1;
t154 = t236 + t285 + t457;
t176 = pkin(1) + t195;
t352 = t285 + t441;
t323 = t281 + t352;
t161 = -t276 + t323;
t424 = pkin(1) * t502;
t180 = -0.2e1 * t424;
t124 = t180 + t161;
t104 = t163 + t124;
t159 = t175 * t242;
t209 = t285 * t521;
t468 = t285 * t287;
t181 = t209 - 0.4e1 * t468;
t255 = 0.2e1 * t276;
t271 = -2 * t287;
t483 = t208 * t285;
t417 = 0.2e1 * t483;
t445 = -t285 + t287;
t127 = t180 + t417 + t445;
t491 = t127 * t205;
t512 = -pkin(4) + pkin(5);
t513 = -pkin(4) - pkin(5);
t288 = sqrt(t181 * t208 + 0.4e1 * t161 * t424 - t279 - (t287 + (pkin(2) - t512) * (pkin(2) + t512)) * (t287 + (pkin(2) - t513) * (pkin(2) + t513)) + (t255 + t271 + 0.2e1 * t273 - 0.6e1 * t285 - 0.4e1 * t491) * t281 + (-t104 * t159 + t124 * t170) * t540);
t258 = 7 * t279;
t266 = 8 * t287;
t252 = -0.2e1 * t273;
t451 = t252 - 0.2e1 * t276;
t203 = t205 ^ 2;
t471 = t279 * t203;
t316 = 0.24e2 * t191 * t471 - t282 - ((21 * t281) + t532) * t293 - (t287 * t451 + t167 + t265 + (35 * t279)) * t285 - (t258 + (t266 + t450) * t281 + t287 * (-t276 + t441)) * t198;
t504 = pkin(1) * t279;
t317 = -0.64e2 * t473 * t504;
t351 = 0.32e2 / 0.3e1 * t279;
t321 = t204 * t351;
t155 = t287 + t285 / 0.4e1 + t281 / 0.4e1 - t273 / 0.8e1;
t322 = 0.64e2 / 0.3e1 * t155 * t296;
t199 = -t281 + t287;
t359 = t287 + t533;
t361 = t273 / 0.3e1 + t229 + t270;
t327 = -0.8e1 / 0.3e1 * t471 + t285 * t199 - 0.5e1 / 0.3e1 * t279 + t361 * t281 + t287 * t359;
t332 = -0.96e2 * t166 * t191 * t296;
t396 = 0.16e2 * t487;
t338 = pkin(1) * t396;
t503 = pkin(1) * t281;
t339 = -0.24e2 * t166 * t503;
t442 = t286 + t293;
t151 = 0.16e2 * (t442 - 0.6e1 * t468) * t279;
t168 = t276 + t353;
t194 = -0.30e2 * t273 + (60 * t287);
t268 = 4 * t287;
t289 = pkin(1) * t287;
t312 = (6 * t286) + t448 - 0.6e1 * t467;
t207 = t243 * t208;
t485 = t207 * t292;
t375 = t176 * t485;
t331 = -0.32e2 * t375;
t397 = pkin(3) * t481;
t490 = t176 * t243;
t220 = -t273 / 0.4e1;
t121 = -0.2e1 / 0.3e1 * t344 + t287 + t236 + t220;
t187 = (t268 + t273) * t281;
t192 = t287 - 0.2e1 / 0.3e1 * t285;
t413 = 0.8e1 * t471;
t474 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t454 = 0.4e1 / 0.7e1 * t287 - t273 / 0.7e1;
t72 = -0.32e2 / 0.21e2 * t155 * t344 + t293 / 0.7e1 + (0.16e2 / 0.21e2 * t281 + t454) * t285 + t279 / 0.7e1 + t454 * t281 + t286 - 0.3e1 / 0.7e1 * t467 + t448 / 0.42e2;
t235 = t281 / 0.3e1;
t455 = t238 + t287;
t156 = t220 + t235 + t455;
t510 = 0.4e1 / 0.3e1 * t285;
t73 = -0.8e1 / 0.3e1 * t156 * t344 + t293 / 0.3e1 + (t234 + t221) * t285 + t286 - t279 / 0.3e1 + (t510 + 0.2e1 / 0.3e1 * t281 + t223) * t287 + t448 / 0.18e2;
t53 = t192 * t413 + 0.14e2 * t72 * t469 + t190 * t315 + t199 * t293 + (t187 - 0.10e2 / 0.3e1 * t279 + (2 * t286) - t467) * t285 + t118 * t474 + (t121 * t396 + 0.6e1 * t73 * t195) * pkin(1);
t113 = -0.4e1 * t296 * t170 + t209 + t522 + (t252 + t266) * t281;
t310 = -t344 + t455;
t119 = t220 - t281 + t310;
t421 = 0.8e1 * t487;
t425 = 0.4e1 * t195;
t65 = t163 * t474 + t113 * t205 + t166 * t200 + (t119 * t425 + t421) * pkin(1);
t443 = t286 - t279;
t66 = t191 * t315 - t282 + (-t211 - t441) * t293 + (t187 + t443 + t531) * t285 + t287 * t118;
t71 = (t224 + t269 + (7 * t281)) * t293 + (t258 + (t253 + (10 * t287)) * t281 + t326) * t285 + t463;
t174 = -12 * pkin(1) * t296 + t289 * t540;
t466 = t287 * t281;
t189 = -8 * t279 + 12 * t466;
t83 = t174 * t242 + t189 * t205 + t338 + t413 + t446 - (6 * t466);
t95 = t163 * t473 + t166 * t201;
t34 = t65 * t331 + t151 * t203 + 0.24e2 * t66 * t469 + (t252 + t268 + (28 * t281)) * t282 + t168 * t488 + (t194 * t279 + 0.24e2 * t53 * t208 + t286 * t254 + t312 * t260 + t448 * t270 + (4 * t289 ^ 2) + (28 * t296 ^ 2)) * t285 + 0.8e1 * (-t397 * t71 - t54 * t490) * pkin(2) + (0.8e1 * t62 * t195 + 0.32e2 * t95 * t487) * pkin(1) + (t194 * t281 + 0.16e2 * t83 * t206 + (70 * t279) + t293 + t312) * t293;
t340 = -0.16e2 * t156 * t503;
t498 = t190 * pkin(3);
t346 = -0.4e1 * t166 * t498;
t347 = t474 * t527;
t348 = -0.48e2 * t94 * t503;
t365 = 0.4e1 * t423;
t482 = t208 * t292;
t374 = t176 * t482;
t59 = 0.1e1 / t288;
t516 = t59 / 0.2e1;
t475 = t243 * t285;
t373 = t241 * t475;
t131 = t179 - 0.4e1 * t373;
t405 = pkin(1) * t535;
t411 = -0.4e1 * t469;
t57 = t131 * t411 + (pkin(1) * t161 * t529 + t181 * t523) * t241 + (-t460 * t159 + 0.2e1 * t241 ^ 2 * t405 + (-t104 * t479 + t124 * t480) * pkin(2)) * t540;
t378 = t57 * t516;
t336 = 0.8e1 * t375;
t394 = -0.12e2 * t483;
t412 = -0.6e1 * t469;
t505 = pkin(1) * t242;
t499 = pkin(3) * t240;
t385 = -t499 / 0.2e1;
t437 = 4 * pkin(1);
t470 = t279 * t204;
t137 = t510 + t236 + t360;
t84 = -t137 * t499 + t138 * t196;
t185 = 0.2e1 * t285 + t199;
t93 = t140 * t196 + t185 * t385;
t212 = -0.20e2 / 0.3e1 * t281;
t362 = 0.2e1 / 0.3e1 * t273 + t228 + t268;
t363 = 0.4e1 / 0.3e1 * t273 + 0.4e1 / 0.3e1 * t276 + t271;
t97 = -t293 + (t212 + t362) * t285 - (3 * t279) + t363 * t281 + t286;
t55 = t84 * t410 + t97 * t385 + t327 * t196 + (t93 * t195 - t240 * t470) * t437;
t101 = t154 * t434 + t473 * t499;
t401 = t296 * t196;
t334 = t204 * t401;
t233 = -0.2e1 / 0.3e1 * t276;
t263 = 0.5e1 * t293;
t267 = 6 * t287;
t355 = t267 + t451;
t100 = t263 + (t250 + t355) * t285 + t198 * (t233 + t356);
t259 = 5 * t279;
t452 = t250 + t270;
t456 = t223 + t233;
t85 = t293 + (t452 + t456) * t285 + t259 + t355 * t281 + t287 * (t287 + t456);
t67 = -t100 * t499 + t85 * t196;
t324 = t281 + t359;
t142 = t264 + t324;
t86 = -t142 * t499 + t143 * t196;
t56 = t101 * t411 + (t425 * t86 - 0.8e1 * t334) * pkin(1) + t67;
t102 = -0.5e1 / 0.3e1 * t293 + (-t281 + t361) * t285 + t287 * t324;
t427 = -0.2e1 * t499;
t98 = -0.3e1 * t293 + (t212 + t363) * t285 + t362 * t281 + t443;
t68 = t102 * t427 + t98 * t196;
t139 = 0.3e1 / 0.2e1 * t285 + t261 + t358;
t402 = t281 * t196;
t311 = -t240 * t296 + t402;
t146 = 0.4e1 * t311;
t158 = t196 + 0.2e1 * t499;
t69 = t200 * t196 + t146 * t205 + (t139 * t240 + t158 * t505) * t436;
t99 = t141 * t196 + t201 * t499 / 0.2e1;
t38 = t99 * t338 + t69 * t336 + t68 * t412 + t55 * t394 + (t536 * t240 - 0.6e1 * t67 * t505) * pkin(3) + (t241 * t316 + 0.6e1 * t56 * t490) * pkin(2);
t392 = 0.24e2 * t483;
t395 = -0.32e2 * t485;
t419 = 0.8e1 * t485;
t430 = 0.6e1 * t502;
t433 = -0.8e1 * t502;
t438 = -0.8e1 * t71 * pkin(3);
t439 = pkin(3) * t538;
t484 = t207 * t293;
t24 = t38 * t378 + (-0.32e2 * t136 * t145 + 0.8e1 * t153 * t288) * t375 + (0.24e2 * (0.4e1 * t122 * t484 * t499 - t374 * t69 + t55 * t475) * t288 + t145 * (0.96e2 * t374 * t65 - 0.48e2 * t53 * t475 - 0.64e2 * t83 * t484)) * t241 + ((t34 + (t145 * t526 - 0.6e1 * t288 * t56) * t176) * t241 + (((t138 * t410 + t140 * t365 + t327) * t394 + t141 * t338 + t98 * t412 - 0.6e1 * t85 * t423 + t316) * t288 + (t120 * t525 + (t534 * t419 + (t143 * t365 - 0.8e1 * t154 * t469 - 0.8e1 * t408 + t85) * t430) * t288) * t176 + ((-pkin(1) * t321 - t205 * t322 + t242 * t340 + t346) * t392 + t204 * t317 + t205 * t332 + t242 * t348 + t438 + ((t347 - 0.4e1 * t465) * t395 + (t242 * t339 + t439) * t433) * t176) * t145 * t240) * t243) * pkin(2);
t32 = t145 * t34 + t288 * t38;
t459 = -t170 + t159;
t115 = pkin(3) + t459;
t123 = t180 + t276 + t323;
t399 = t127 * t195;
t79 = t123 * t175 + 0.2e1 * t399;
t501 = pkin(3) * t175;
t80 = t123 * t242 + (0.4e1 * t205 - 0.2e1) * t501;
t51 = t115 * t288 + t80 * t196 + t240 * t79;
t126 = t178 + t534;
t108 = -0.4e1 / 0.9e1 * t344 + t287 + t285 / 0.3e1 + t235 + t276 / 0.9e1 - t273 / 0.9e1;
t116 = t276 / 0.6e1 + t219 + t310;
t63 = t190 * t163 + 0.6e1 * t108 * t469 + t144 * t474 + (t116 * t195 + t487) * t437;
t114 = t163 + t144;
t182 = t445 * t521;
t74 = t144 * t328 + (t260 + t357) * t285 + t461;
t64 = t114 * t365 + t182 * t205 + t74;
t89 = t144 * t473 + t129;
t47 = -0.8e1 * t126 * t375 + t149 + t89 * t389 + t74 * t364 + t282 + (-t273 + t276 + t453) * t293 + t193 * t168 + (0.12e2 * t63 * t208 + t248 + (t251 + t254 + 0.6e1 * t276) * t281 + t265 + (t252 + t255) * t287) * t285 + 0.6e1 * (-t397 * t92 - t64 * t490) * pkin(2);
t150 = t186 * t196;
t184 = t264 + t198;
t489 = t184 * t240;
t107 = -pkin(3) * t489 + t150;
t183 = pkin(2) * t522 + 0.8e1 * t281 * t292;
t372 = t296 * t473;
t109 = t183 * t241 + 0.4e1 * t240 * t372;
t125 = t285 * t452 + t259 + t442 + (6 * t466);
t134 = t263 + (t250 + t267) * t285 + t193;
t418 = -0.4e1 * t483;
t157 = t196 - t499;
t335 = t205 * t402;
t75 = -0.2e1 * t335 + t150 + (0.2e1 * t157 * t505 - t489) * pkin(3);
t147 = 0.2e1 * t311;
t197 = (2 * t281) + t285;
t78 = t197 * t499 + t147 * t205 + (t199 + t178) * t196;
t52 = t78 * t418 + t109 * t205 + (-0.4e1 * t107 * t505 + (t134 + t336) * t240) * pkin(3) + (0.4e1 * t75 * t490 + (-t125 + t350) * t241) * pkin(2);
t37 = t145 * t47 + t288 * t52;
t518 = 0.1e1 / t37 / 0.4e1;
t383 = t51 * t518;
t160 = t261 + t276 + t352;
t105 = t160 + t180 + t328;
t162 = t170 - pkin(3);
t435 = pkin(1) * t523;
t117 = t175 * t240 + t403;
t492 = t117 * t288;
t50 = -t105 * t159 + t492 + (t160 * t481 + t162 * t435) * pkin(2) + (-t168 - t264 + t417 - 0.2e1 * t491) * pkin(3);
t508 = t288 / 0.4e1;
t309 = t32 * t383 + t50 * t508;
t398 = pkin(3) * t159;
t90 = t163 + t180 + t354 + 0.2e1 * t398;
t88 = 0.1e1 / t90 ^ 2;
t307 = t88 * t309;
t393 = 0.12e2 * t483;
t477 = t242 * t281;
t415 = -0.8e1 * t477;
t416 = 0.8e1 * t241 * t285;
t429 = t144 * t528;
t431 = 0.4e1 * t502;
t440 = t75 * t529;
t33 = (t153 * t418 + (t78 * t416 + t183 * t205 + ((t199 + 0.2e1 * t469) * t418 - t125 + (-0.4e1 * t186 * t195 + t421) * pkin(1)) * pkin(2)) * t243 + ((t153 + (t186 - 0.2e1 * t469) * t502) * t431 + (-0.24e2 * t482 * t499 + t440) * t241) * t176) * t288 + t52 * t378 + t47 * t196 + t145 * (0.24e2 * t126 * t241 * t374 + (t136 + (-0.8e1 / 0.3e1 * t465 - 0.2e1 * t498) * t404) * t393 - 0.24e2 * t63 * t373 + t120 + t313 * t537 + t343 * t539 + 0.6e1 * (-(pkin(1) * t415 + t429) * t208 * t535 + t64 * t196) * t176);
t380 = t50 * t59 / 0.8e1;
t517 = -0.1e1 / t37 ^ 2 / 0.4e1;
t381 = t51 * t517;
t368 = t115 * t516;
t426 = 0.2e1 * t195;
t500 = pkin(3) * t205;
t41 = t57 * t368 + t131 * t240 * t426 + ((-t240 * t288 + t80) * t243 + (t242 * t288 + (t175 * t530 + t123) * t240 + (-pkin(3) + 0.2e1 * t500 + t505) * t434) * t241) * pkin(2);
t367 = t117 * t516;
t43 = (t170 + t400) * t288 + t57 * t367 - 0.2e1 * t131 * t500 - (t179 - 0.4e1 * t343) * t159 - 0.2e1 * t208 * t405 + t160 * t404 + (t475 * t528 + (-t105 * t242 + t162 * t530) * pkin(2)) * t241;
t472 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t87 = 0.1e1 / t90;
t10 = (-t106 * t307 + (t43 * t508 + t57 * t380 + t24 * t383 + (t33 * t381 + t41 * t518) * t32) * t87) * t472;
t133 = -t479 + t480;
t132 = -t478 - t481;
t384 = t50 * t518;
t509 = -t288 / 0.4e1;
t308 = t32 * t384 + t51 * t509;
t386 = t87 * t472;
t29 = t308 * t386;
t28 = t29 * t132;
t30 = t309 * t386;
t22 = -t133 * t30 - t28;
t20 = 0.1e1 / t22 ^ 2;
t27 = t29 * t133;
t21 = -t132 * t30 + t27;
t13 = 0.1e1 / (t20 * t21 ^ 2 + 0.1e1);
t19 = 0.1e1 / t22;
t497 = t20 * t21;
t306 = t88 * t308;
t379 = -t51 * t59 / 0.8e1;
t382 = t50 * t517;
t506 = t30 + (-t106 * t306 + (t24 * t384 + t41 * t509 + t57 * t379 + (t33 * t382 + t43 * t518) * t32) * t87) * t472;
t520 = ((t506 * t133 + (-t10 + t29) * t132) * t19 - (-t10 * t133 - t506 * t132 + t27) * t497) * t13 + 0.1e1;
t519 = pkin(1) * pkin(3);
t515 = t87 / 0.2e1;
t514 = -t88 / 0.2e1;
t277 = 0.1e1 / pkin(4);
t49 = 0.1e1 / t50 ^ 2;
t341 = pkin(4) * t277 / (t49 * t51 ^ 2 + 0.1e1) * t90;
t314 = t49 * t51 * t341;
t320 = 0.1e1 / t50 * t341;
t110 = t117 * t436;
t369 = t110 * t514;
t202 = t240 ^ 2;
t414 = 0.8e1 * t477;
t428 = 0.4e1 * t501;
t61 = (t104 * t428 + t127 * t414) * t240 + (t124 * t425 + 0.8e1 * t175 * t469) * t196;
t44 = -t492 + t61 * t368 + t127 * t202 * t527 + t79 * t242 + (-t123 - 0.8e1 * t398) * t170;
t476 = t242 * t285;
t45 = t459 * t288 + t61 * t367 + (t105 * t175 + 0.4e1 * t399) * t240 + (t435 * t476 + (t160 * t242 + t205 * t428) * pkin(2)) * t241;
t511 = 0.2e1 * (-t51 * t369 + t44 * t515) * t320 - 0.2e1 * (-t50 * t369 + t45 * t515) * t314 + 0.1e1;
t318 = pkin(1) * t335;
t152 = -0.4e1 * t318;
t177 = pkin(1) * t427;
t319 = pkin(1) * t205 * t401;
t329 = t204 * t372;
t406 = pkin(1) * t465;
t337 = -0.48e2 * t406;
t407 = pkin(1) * t469;
t349 = 0.4e1 * t407;
t371 = t240 * t477;
t377 = t61 * t516;
t390 = -0.32e2 * t470;
t391 = -0.24e2 * t477;
t409 = pkin(1) * t471;
t420 = -0.8e1 * t485;
t422 = -0.4e1 * t487;
t432 = -0.6e1 * t502;
t462 = pkin(1) * t201 * t330 - 0.24e2 * t191 * t334;
t464 = -4 * t519;
t23 = ((t139 * t426 + t349 + t422) * t336 + (0.12e2 * t202 * t205 * t504 + t137 * t422 - 0.2e1 * t185 * t407 - 0.4e1 * t409) * t394 + 0.8e1 * t201 * t409 + 0.12e2 * t102 * t487 + 0.6e1 * t100 * t407 + (-t97 * t394 / 0.2e1 + t536) * t195) * t288 + t38 * t377 + t145 * t152 * t331 + ((0.6e1 * (-t100 * t195 - 0.4e1 * t142 * t407 - 0.4e1 * t329) * t288 + t462 * t525) * t490 + ((-pkin(1) * t203 * t351 - t204 * t322 + t205 * t340 + t242 * t346) * t392 + t203 * t317 + t204 * t332 + t205 * t348 + t242 * t438 + ((t242 * t347 + t422) * t395 + (t205 * t339 + t242 * t439) * t433) * t176) * t145 * t241) * pkin(2) + (0.12e2 * (-(t321 * t196 + t415 * t84) * t483 - 0.8e1 * t191 * t470 * t196 - 0.4e1 * t99 * t406 + t68 * t477) * t288 + t145 * (0.16e2 * (t189 * t524 - t174 + t337 + t390) * t486 + (t121 * t337 + t192 * t390 - 0.28e2 * t72 * t477) * t392 - 0.4e1 * t151 * t204 - 0.96e2 * t95 * t406 - 0.48e2 * t66 * t477) + (-t34 + t145 * (0.32e2 * t65 * t485 + t502 * t526) + (t56 * t432 + t69 * t420 + (-0.24e2 * t177 + 0.64e2 * t371) * t486) * t288 + ((0.48e2 * t93 * t483 + 0.6e1 * t67) * t288 + t145 * (-0.144e3 * t73 * t483 - 0.8e1 * t62)) * pkin(1)) * pkin(3) + ((0.2e1 * (-t146 * t242 - t158 * t519) * t419 + (t101 * t414 + t86 * t464 + 0.24e2 * t319) * t430) * t288 + t145 * ((pkin(1) * t388 + t113 * t524 + t119 * t464) * t395 + (t391 * t91 - 0.6e1 * t77 * t519) * t433)) * t176) * t240;
t31 = (t281 * t202 * t420 + (t197 * t195 - 0.2e1 * t487) * t418 + 0.4e1 * t329 + t184 * t349 + t134 * t195) * t288 + t52 * t377 + t145 * ((-0.8e1 / 0.3e1 * t334 + t152 - 0.2e1 * t190 * t342) * t393 + t318 * t537 + t342 * t539 + t462) + ((0.8e1 * t147 * t208 * t476 + t109 * t524 - 0.24e2 * t319) * t288 + t145 * (0.12e2 * (-t108 * t477 - t406) * t393 + t89 * t391) + (t243 * t288 * t440 - t47 + t145 * (t126 * t419 + t430 * t64) + ((pkin(2) * t208 * t416 + 0.4e1 * t107) * t288 + t145 * (-0.48e2 * t116 * t483 - 0.6e1 * t74)) * pkin(1)) * pkin(3)) * t240 + ((t419 * t195 + ((-pkin(3) * t184 + t170 * t521) * t242 + (-t157 * t499 - t469) * t530) * t431) * t288 + t145 * ((t177 - 0.8e1 * t371) * t420 + ((-0.2e1 * t182 * t240 + t429 * t196) * t242 + (-0.4e1 * t114 * t499 - 0.8e1 * t335) * pkin(1)) * t432)) * t176;
t507 = t30 - (t110 * t306 + (t23 * t384 + t44 * t509 + t61 * t379 + (t31 * t382 + t45 * t518) * t32) * t87) * t472;
t14 = atan2(t21, t22);
t11 = sin(t14);
t496 = t11 * t241;
t12 = cos(t14);
t495 = t12 * t241;
t366 = t277 * t515;
t42 = atan2(t51 * t366, t50 * t366);
t39 = sin(t42);
t494 = t240 * t39;
t40 = cos(t42);
t493 = t240 * t40;
t387 = t520 * t243;
t376 = t242 * t511;
t370 = t106 * t514;
t25 = 0.2e1 * (t51 * t370 + t41 * t515) * t320 - 0.2e1 * (t50 * t370 + t43 * t515) * t314;
t18 = (-t242 * t39 - t493) * t25;
t17 = (t242 * t40 - t494) * t25;
t16 = -t39 * t376 - t511 * t493;
t15 = -t40 * t376 + t511 * t494;
t8 = (t110 * t307 + (t45 * t508 + t61 * t380 + t23 * t383 + (t31 * t381 + t44 * t518) * t32) * t87) * t472;
t5 = ((-t132 * t8 - t507 * t133 - t28) * t19 - ((-t29 - t8) * t133 + t507 * t132) * t497) * t13;
t4 = (t11 * t243 + t495) * t5;
t3 = (-t12 * t243 + t496) * t5;
t2 = t11 * t387 + t520 * t495;
t1 = t12 * t387 - t520 * t496;
t6 = [-m(2) * (g(1) * (-rSges(2,1) * t241 - rSges(2,2) * t243) + g(2) * (rSges(2,1) * t243 - rSges(2,2) * t241)) - m(3) * (g(1) * (rSges(3,1) * t2 + rSges(3,2) * t1 - t196) + g(2) * (-rSges(3,1) * t1 + rSges(3,2) * t2 + t502)) - m(5) * (g(1) * (rSges(5,1) * t18 - rSges(5,2) * t17) + g(2) * (rSges(5,1) * t17 + rSges(5,2) * t18)), -m(3) * (g(1) * (rSges(3,1) * t4 - rSges(3,2) * t3) + g(2) * (rSges(3,1) * t3 + rSges(3,2) * t4)) - m(4) * (g(1) * (-rSges(4,1) * t240 - rSges(4,2) * t242) + g(2) * (rSges(4,1) * t242 - rSges(4,2) * t240)) - m(5) * (g(1) * (rSges(5,1) * t16 + rSges(5,2) * t15 - t499) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16 + t195))];
taug = t6(:);
