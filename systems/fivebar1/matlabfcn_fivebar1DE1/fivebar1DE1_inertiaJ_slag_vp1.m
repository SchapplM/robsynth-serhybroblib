% Calculate joint inertia matrix for
% fivebar1DE1
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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1DE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1DE1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1DE1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:24:36
% EndTime: 2020-04-27 02:25:07
% DurationCPUTime: 28.00s
% Computational Cost: add. (515895->670), mult. (1535274->1066), div. (3976->13), fcn. (187660->12), ass. (0->472)
t283 = pkin(5) ^ 2;
t282 = t283 ^ 2;
t286 = pkin(4) ^ 2;
t285 = t286 ^ 2;
t459 = t282 - t285;
t545 = 4 * pkin(3);
t535 = 2 * pkin(1);
t447 = 2 * pkin(3);
t291 = pkin(3) ^ 2;
t221 = 0.10e2 / 0.3e1 * t291;
t295 = pkin(2) ^ 2;
t233 = -0.2e1 / 0.3e1 * t283;
t238 = 0.2e1 / 0.3e1 * t286;
t297 = pkin(1) ^ 2;
t280 = 2 * t297;
t369 = t233 + t238 + t280;
t208 = t291 + t297;
t303 = t295 ^ 2;
t368 = t233 + t208;
t472 = t208 * (t238 + t368) + t303;
t98 = (t221 + t369) * t295 + t472;
t544 = -0.6e1 * t98;
t289 = t291 ^ 2;
t296 = t297 ^ 2;
t457 = t289 + t296;
t460 = t280 - t283;
t478 = t297 * t283;
t536 = -t282 / 0.6e1 + t285 / 0.6e1;
t124 = t460 * t291 + t457 - t478 - t536;
t315 = t124 + t303;
t100 = (t221 + t460) * t295 + t315;
t543 = -0.6e1 * t100;
t251 = sin(qJ(1));
t206 = pkin(2) * t251;
t445 = 0.2e1 * t206;
t526 = 4 * t291;
t270 = 6 * t291;
t231 = -t283 / 0.3e1;
t239 = t286 / 0.3e1;
t455 = t295 + t297;
t366 = t291 + t455;
t150 = t231 + t239 + t366;
t542 = -0.24e2 * t150;
t260 = 10 * t291;
t252 = cos(qJ(2));
t205 = pkin(3) * t252;
t433 = pkin(1) * t205;
t188 = 0.2e1 * t433;
t200 = -t291 / 0.3e1 + t297;
t215 = t252 ^ 2;
t480 = t291 * t215;
t128 = 0.4e1 / 0.3e1 * t480 + t188 + t200;
t279 = 3 * t297;
t537 = t279 - t283 - t286;
t177 = t537 * t260;
t264 = -0.6e1 * t283;
t302 = pkin(2) * t295;
t292 = t302 ^ 2;
t263 = -0.5e1 * t283;
t461 = t263 - 0.5e1 * t286;
t253 = cos(qJ(1));
t218 = t253 ^ 2;
t216 = t218 ^ 2;
t497 = t216 * t303;
t203 = t208 ^ 2;
t452 = t297 - t283;
t365 = t291 + t452;
t499 = t203 * (-t286 + t365);
t541 = 0.7e1 * t292 + ((35 * t291) + (15 * t297) + t461) * t303 + ((21 * t289) + t177 + (9 * t296) + (t264 - 0.6e1 * t286) * t297) * t295 + t499 - 0.24e2 * t128 * t497;
t250 = sin(qJ(2));
t540 = t250 * t295;
t210 = -3 * t291 + t297;
t420 = 0.4e1 * t480;
t539 = t210 + t420;
t538 = t231 - t286 / 0.3e1;
t534 = -0.4e1 * pkin(2);
t533 = -4 * pkin(3);
t532 = -2 * pkin(3);
t211 = -0.3e1 * t295 + t297;
t214 = t252 * t215;
t306 = pkin(3) * t291;
t498 = t214 * t306;
t418 = pkin(1) * t498;
t362 = 0.8e1 * t418;
t155 = t211 * t362;
t376 = 0.6e1 * t433;
t399 = 0.12e2 * t480;
t234 = -0.3e1 / 0.2e1 * t283;
t258 = 15 * t289;
t261 = 18 * t297;
t275 = 3 * t296;
t469 = t282 / 0.2e1 - t285 / 0.2e1;
t337 = -0.3e1 * t478 + t275 + t469;
t492 = t250 * t251;
t180 = pkin(2) * t492;
t357 = pkin(3) * t180;
t464 = 15 * t291 + t279;
t474 = t208 * ((t234 + t280) * t291 - 0.3e1 / 0.2e1 * t478 + t457 + t469) + t292;
t68 = t357 * t543 + (t258 + (t261 - 0.9e1 * t283) * t291 + t337) * t295 + (t234 + t464) * t303 + t474;
t232 = -t283 / 0.2e1;
t176 = t232 + t366;
t338 = -0.4e1 * t357;
t325 = t176 * t338;
t83 = t325 + (t270 + t460) * t295 + t315;
t173 = -0.2e1 * t357;
t201 = t297 - t295 / 0.3e1;
t135 = t201 * t173;
t484 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t97 = t176 * t484 + t135;
t60 = t376 * t83 + t399 * t97 + t155 + t68;
t531 = 0.8e1 * t60;
t511 = pkin(2) * t253;
t185 = pkin(1) - t511;
t151 = t185 + t205;
t530 = -0.8e1 * t151;
t529 = -0.2e1 * t252;
t528 = -0.2e1 * t253;
t527 = 4 * t289;
t525 = pkin(1) * pkin(3);
t364 = t295 + t452;
t334 = t291 + t364;
t171 = -t286 + t334;
t434 = pkin(1) * t511;
t190 = -0.2e1 * t434;
t130 = t190 + t171;
t110 = t173 + t130;
t169 = t185 * t252;
t219 = t295 * t526;
t479 = t295 * t297;
t191 = t219 - 0.4e1 * t479;
t265 = 0.2e1 * t286;
t281 = -2 * t297;
t494 = t218 * t295;
t427 = 0.2e1 * t494;
t456 = -t295 + t297;
t133 = t190 + t427 + t456;
t502 = t133 * t215;
t518 = pkin(5) - pkin(4);
t519 = -pkin(4) - pkin(5);
t298 = sqrt(t191 * t218 + 0.4e1 * t171 * t434 - t289 - (t297 + (pkin(2) - t518) * (pkin(2) + t518)) * (t297 + (pkin(2) - t519) * (pkin(2) + t519)) + (t265 + t281 + 0.2e1 * t283 - 0.6e1 * t295 - 0.4e1 * t502) * t291 + (-t110 * t169 + t130 * t180) * t545);
t132 = t188 + t539;
t178 = t286 + t365;
t262 = -0.2e1 * t283;
t186 = pkin(1) + t205;
t217 = t253 * t218;
t496 = t217 * t302;
t387 = t186 * t496;
t407 = pkin(3) * t492;
t501 = t186 * t253;
t245 = t291 / 0.3e1;
t114 = -0.4e1 / 0.9e1 * t357 + t297 + t295 / 0.3e1 + t245 + t286 / 0.9e1 - t283 / 0.9e1;
t229 = -t283 / 0.6e1;
t248 = t295 / 0.2e1;
t466 = t248 + t297;
t320 = -t357 + t466;
t122 = t286 / 0.6e1 + t229 + t320;
t448 = 4 * pkin(1);
t485 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t69 = t200 * t173 + 0.6e1 * t114 * t480 + t150 * t485 + (t122 * t205 + t498) * t448;
t120 = t173 + t150;
t192 = t456 * t526;
t377 = 0.4e1 * t433;
t80 = t150 * t338 + (t270 + t369) * t295 + t472;
t70 = t120 * t377 + t192 * t215 + t80;
t95 = t150 * t484 + t135;
t53 = -0.8e1 * t132 * t387 + t155 + t95 * t399 + t80 * t376 + t292 + (-t283 + t286 + t464) * t303 + t203 * t178 + (0.12e2 * t69 * t218 + t258 + (t261 + t264 + 0.6e1 * t286) * t291 + t275 + (t262 + t265) * t297) * t295 + 0.6e1 * (-t407 * t98 - t70 * t501) * pkin(2);
t271 = 3 * t291;
t196 = t271 + t455;
t156 = t196 * t206;
t274 = 0.3e1 * t295;
t194 = t274 + t208;
t500 = t194 * t250;
t113 = -pkin(3) * t500 + t156;
t193 = pkin(2) * t527 + 0.8e1 * t291 * t302;
t384 = t306 * t484;
t115 = t193 * t251 + 0.4e1 * t250 * t384;
t269 = 5 * t289;
t453 = t296 + t303;
t463 = t260 + t280;
t477 = t297 * t291;
t131 = t295 * t463 + t269 + t453 + (6 * t477);
t273 = 0.5e1 * t303;
t277 = 6 * t297;
t140 = t273 + (t260 + t277) * t295 + t203;
t346 = 0.8e1 * t387;
t428 = -0.4e1 * t494;
t514 = pkin(1) * t252;
t508 = pkin(3) * t250;
t167 = t206 - t508;
t412 = t291 * t206;
t345 = t215 * t412;
t81 = -0.2e1 * t345 + t156 + (0.2e1 * t167 * t514 - t500) * pkin(3);
t321 = -t250 * t306 + t412;
t153 = 0.2e1 * t321;
t207 = (2 * t291) + t295;
t209 = -t291 + t297;
t84 = t207 * t508 + t153 * t215 + (t209 + t188) * t206;
t58 = t84 * t428 + t115 * t215 + (-0.4e1 * t113 * t514 + (t140 + t346) * t250) * pkin(3) + (0.4e1 * t81 * t501 + (-t131 + t362) * t251) * pkin(2);
t43 = t151 * t53 + t298 * t58;
t524 = 0.1e1 / t43 / 0.4e1;
t523 = -0.1e1 / t43 ^ 2 / 0.4e1;
t65 = 0.1e1 / t298;
t522 = t65 / 0.2e1;
t408 = pkin(3) * t169;
t96 = t173 + t190 + t366 + 0.2e1 * t408;
t93 = 0.1e1 / t96;
t521 = t93 / 0.2e1;
t94 = 0.1e1 / t96 ^ 2;
t520 = -t94 / 0.2e1;
t517 = 0.4e1 / 0.3e1 * t295;
t516 = -t298 / 0.4e1;
t515 = t298 / 0.4e1;
t513 = pkin(1) * t289;
t512 = pkin(1) * t291;
t510 = pkin(3) * t185;
t509 = pkin(3) * t215;
t490 = t252 * t251;
t491 = t250 * t253;
t139 = -t490 + t491;
t489 = t252 * t253;
t138 = -t489 - t492;
t101 = t173 * t484 + t176 * t211;
t157 = 0.16e2 * (t453 - 0.6e1 * t479) * t289;
t204 = -0.30e2 * t283 + (60 * t297);
t213 = t215 ^ 2;
t278 = 4 * t297;
t299 = pkin(1) * t297;
t322 = (6 * t296) + t459 - 0.6e1 * t478;
t341 = -0.32e2 * t387;
t230 = -t283 / 0.4e1;
t246 = t291 / 0.2e1;
t127 = -0.2e1 / 0.3e1 * t357 + t297 + t246 + t230;
t197 = (t278 + t283) * t291;
t202 = t297 - 0.2e1 / 0.3e1 * t295;
t406 = 0.16e2 * t498;
t482 = t289 * t213;
t423 = 0.8e1 * t482;
t165 = t297 + t295 / 0.4e1 + t291 / 0.4e1 - t283 / 0.8e1;
t465 = 0.4e1 / 0.7e1 * t297 - t283 / 0.7e1;
t78 = -0.32e2 / 0.21e2 * t165 * t357 + t303 / 0.7e1 + (0.16e2 / 0.21e2 * t291 + t465) * t295 + t289 / 0.7e1 + t465 * t291 + t296 - 0.3e1 / 0.7e1 * t478 + t459 / 0.42e2;
t166 = t230 + t245 + t466;
t244 = 0.4e1 / 0.3e1 * t291;
t79 = -0.8e1 / 0.3e1 * t166 * t357 + t303 / 0.3e1 + (t244 + t231) * t295 + t296 - t289 / 0.3e1 + (t517 + 0.2e1 / 0.3e1 * t291 + t233) * t297 + t459 / 0.18e2;
t59 = t202 * t423 + 0.14e2 * t78 * t480 + t200 * t325 + t209 * t303 + (t197 - 0.10e2 / 0.3e1 * t289 + (2 * t296) - t478) * t295 + t124 * t485 + (t127 * t406 + 0.6e1 * t79 * t205) * pkin(1);
t276 = 8 * t297;
t119 = -0.4e1 * t306 * t180 + t219 + t527 + (t262 + t276) * t291;
t125 = t230 - t291 + t320;
t431 = 0.8e1 * t498;
t435 = 0.4e1 * t205;
t71 = t173 * t485 + t119 * t215 + t176 * t210 + (t125 * t435 + t431) * pkin(1);
t454 = t296 - t289;
t72 = t201 * t325 - t292 + (-t221 - t452) * t303 + (t197 + t454 + t536) * t295 + t297 * t124;
t268 = 7 * t289;
t77 = (t234 + t279 + (7 * t291)) * t303 + (t268 + (t263 + (10 * t297)) * t291 + t337) * t295 + t474;
t184 = -12 * pkin(1) * t306 + t299 * t545;
t199 = -8 * t289 + 12 * t477;
t348 = pkin(1) * t406;
t89 = t184 * t252 + t199 * t215 + t348 + t423 + t457 - (6 * t477);
t40 = t71 * t341 + t157 * t213 + 0.24e2 * t72 * t480 + (t262 + t278 + (28 * t291)) * t292 + t178 * t499 + (t204 * t289 + 0.24e2 * t59 * t218 + t296 * t264 + t322 * t270 + t459 * t280 + (4 * t299 ^ 2) + (28 * t306 ^ 2)) * t295 + 0.8e1 * (-t407 * t77 - t60 * t501) * pkin(2) + (0.32e2 * t101 * t498 + 0.8e1 * t68 * t205) * pkin(1) + (t204 * t291 + 0.16e2 * t89 * t216 + (70 * t289) + t303 + t322) * t303;
t370 = t232 - t286 / 0.2e1 + t297;
t147 = t274 + 0.3e1 / 0.2e1 * t291 + t370;
t105 = t147 * t206 + t211 * t508 / 0.2e1;
t462 = t262 - 0.2e1 * t286;
t326 = 0.24e2 * t201 * t482 - t292 - ((21 * t291) + t537) * t303 - (t297 * t462 + t177 + t275 + (35 * t289)) * t295 - (t268 + (t276 + t461) * t291 + t297 * (-t286 + t452)) * t208;
t404 = -0.12e2 * t494;
t422 = -0.6e1 * t480;
t222 = -0.20e2 / 0.3e1 * t291;
t374 = 0.2e1 / 0.3e1 * t283 + t238 + t278;
t375 = 0.4e1 / 0.3e1 * t283 + 0.4e1 / 0.3e1 * t286 + t281;
t103 = -t303 + (t222 + t374) * t295 - (3 * t289) + t375 * t291 + t296;
t371 = t297 + t538;
t373 = t283 / 0.3e1 + t239 + t280;
t333 = -0.8e1 / 0.3e1 * t482 + t295 * t209 - 0.5e1 / 0.3e1 * t289 + t373 * t291 + t297 * t371;
t396 = -t508 / 0.2e1;
t481 = t289 * t214;
t468 = t229 - t286 / 0.6e1;
t372 = t297 + t468;
t143 = t517 + t246 + t372;
t336 = t248 + t372;
t144 = t244 + t336;
t90 = -t143 * t508 + t144 * t206;
t146 = t291 + t336;
t195 = 0.2e1 * t295 + t209;
t99 = t146 * t206 + t195 * t396;
t61 = t90 * t420 + t103 * t396 + t333 * t206 + (t99 * t205 - t250 * t481) * t448;
t164 = t246 + t295 + t468;
t107 = t164 * t445 + t484 * t508;
t411 = t306 * t206;
t344 = t214 * t411;
t421 = -0.4e1 * t480;
t243 = -0.2e1 / 0.3e1 * t286;
t367 = t277 + t462;
t106 = t273 + (t260 + t367) * t295 + t208 * (t243 + t368);
t467 = t233 + t243;
t91 = t303 + (t463 + t467) * t295 + t269 + t367 * t291 + t297 * (t297 + t467);
t73 = -t106 * t508 + t91 * t206;
t335 = t291 + t371;
t148 = t274 + t335;
t149 = t196 + t538;
t92 = -t148 * t508 + t149 * t206;
t62 = t107 * t421 + (t435 * t92 - 0.8e1 * t344) * pkin(1) + t73;
t104 = -0.3e1 * t303 + (t222 + t375) * t295 + t374 * t291 + t454;
t108 = -0.5e1 / 0.3e1 * t303 + (-t291 + t373) * t295 + t297 * t335;
t437 = -0.2e1 * t508;
t74 = t104 * t206 + t108 * t437;
t145 = 0.3e1 / 0.2e1 * t295 + t271 + t370;
t152 = 0.4e1 * t321;
t168 = t206 + 0.2e1 * t508;
t75 = t210 * t206 + t152 * t215 + (t145 * t250 + t168 * t514) * t447;
t44 = t105 * t348 + t75 * t346 + t74 * t422 + t61 * t404 + (t541 * t250 - 0.6e1 * t73 * t514) * pkin(3) + (t251 * t326 + 0.6e1 * t62 * t501) * pkin(2);
t36 = t151 * t40 + t298 * t44;
t170 = t271 + t286 + t364;
t111 = t170 + t190 + t338;
t172 = t180 - pkin(3);
t446 = pkin(1) * t528;
t413 = pkin(2) * t490;
t123 = t185 * t250 + t413;
t503 = t123 * t298;
t56 = -t111 * t169 + t503 + (t170 * t492 + t172 * t446) * pkin(2) + (-t178 - t274 + t427 - 0.2e1 * t502) * pkin(3);
t395 = t56 * t524;
t470 = -t180 + t169;
t121 = pkin(3) + t470;
t129 = t190 + t286 + t334;
t409 = t133 * t205;
t85 = t129 * t185 + 0.2e1 * t409;
t86 = t129 * t252 + (0.4e1 * t215 - 0.2e1) * t510;
t57 = t121 * t298 + t86 * t206 + t250 * t85;
t318 = t36 * t395 + t57 * t516;
t483 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t397 = t93 * t483;
t31 = t318 * t397;
t30 = t31 * t138;
t394 = t57 * t524;
t319 = t36 * t394 + t56 * t515;
t32 = t319 * t397;
t21 = -t139 * t32 - t30;
t19 = 0.1e1 / t21 ^ 2;
t29 = t31 * t139;
t20 = -t138 * t32 + t29;
t507 = t19 * t20;
t506 = t200 * pkin(3);
t116 = t123 * t447;
t328 = pkin(1) * t345;
t162 = -0.4e1 * t328;
t187 = pkin(1) * t437;
t212 = t250 ^ 2;
t327 = -0.64e2 * t484 * t513;
t329 = pkin(1) * t215 * t411;
t363 = 0.32e2 / 0.3e1 * t289;
t331 = t214 * t363;
t332 = 0.64e2 / 0.3e1 * t165 * t306;
t339 = t214 * t384;
t342 = -0.96e2 * t176 * t201 * t306;
t476 = t306 * t215;
t416 = pkin(1) * t476;
t347 = -0.48e2 * t416;
t349 = -0.24e2 * t176 * t512;
t350 = -0.16e2 * t166 * t512;
t351 = -0.48e2 * t100 * t512;
t359 = -0.4e1 * t176 * t506;
t360 = t485 * t532;
t417 = pkin(1) * t480;
t361 = 0.4e1 * t417;
t488 = t252 * t291;
t383 = t250 * t488;
t424 = 0.8e1 * t488;
t438 = 0.4e1 * t510;
t67 = (t110 * t438 + t133 * t424) * t250 + (t130 * t435 + 0.8e1 * t185 * t480) * t206;
t388 = t67 * t522;
t398 = -0.24e2 * t476;
t400 = -0.32e2 * t481;
t401 = -0.24e2 * t488;
t402 = 0.24e2 * t494;
t405 = -0.32e2 * t496;
t419 = pkin(1) * t482;
t425 = -0.8e1 * t488;
t429 = 0.8e1 * t496;
t430 = -0.8e1 * t496;
t432 = -0.4e1 * t498;
t436 = 0.2e1 * t205;
t440 = pkin(3) * t543;
t441 = 0.6e1 * t511;
t443 = -0.6e1 * t511;
t444 = -0.8e1 * t511;
t449 = -0.8e1 * t77 * pkin(3);
t340 = t250 * t398;
t473 = pkin(1) * t211 * t340 - 0.24e2 * t201 * t344;
t475 = -4 * t525;
t24 = ((t145 * t436 + t361 + t432) * t346 + (0.12e2 * t212 * t215 * t513 + t143 * t432 - 0.2e1 * t195 * t417 - 0.4e1 * t419) * t404 + 0.8e1 * t211 * t419 + 0.12e2 * t108 * t498 + 0.6e1 * t106 * t417 + (-t103 * t404 / 0.2e1 + t541) * t205) * t298 + t44 * t388 + t151 * t162 * t341 + ((0.6e1 * (-t106 * t205 - 0.4e1 * t148 * t417 - 0.4e1 * t339) * t298 + t473 * t530) * t501 + ((-pkin(1) * t213 * t363 - t214 * t332 + t215 * t350 + t252 * t359) * t402 + t213 * t327 + t214 * t342 + t215 * t351 + t252 * t449 + ((t252 * t360 + t432) * t405 + (t215 * t349 + t252 * t440) * t444) * t186) * t151 * t251) * pkin(2) + (0.12e2 * (-(t331 * t206 + t425 * t90) * t494 - 0.8e1 * t201 * t481 * t206 - 0.4e1 * t105 * t416 + t74 * t488) * t298 + t151 * (0.16e2 * (t199 * t529 - t184 + t347 + t400) * t497 + (t127 * t347 + t202 * t400 - 0.28e2 * t78 * t488) * t402 - 0.4e1 * t157 * t214 - 0.96e2 * t101 * t416 - 0.48e2 * t72 * t488) + (-t40 + t151 * (0.32e2 * t71 * t496 + t511 * t531) + (t62 * t443 + t75 * t430 + (-0.24e2 * t187 + 0.64e2 * t383) * t497) * t298 + ((0.48e2 * t99 * t494 + 0.6e1 * t73) * t298 + t151 * (-0.144e3 * t79 * t494 - 0.8e1 * t68)) * pkin(1)) * pkin(3) + ((0.2e1 * (-t152 * t252 - t168 * t525) * t429 + (t107 * t424 + t92 * t475 + 0.24e2 * t329) * t441) * t298 + t151 * ((pkin(1) * t398 + t119 * t529 + t125 * t475) * t405 + (t401 * t97 - 0.6e1 * t83 * t525) * t444)) * t186) * t250;
t316 = t94 * t318;
t355 = pkin(3) * t413;
t403 = 0.12e2 * t494;
t426 = 0.8e1 * t251 * t295;
t439 = t150 * t533;
t442 = 0.4e1 * t511;
t450 = t81 * t534;
t487 = t252 * t295;
t35 = (t291 * t212 * t430 + (t207 * t205 - 0.2e1 * t498) * t428 + 0.4e1 * t339 + t194 * t361 + t140 * t205) * t298 + t58 * t388 + t151 * ((-0.8e1 / 0.3e1 * t344 + t162 - 0.2e1 * t200 * t355) * t403 + t328 * t542 + t355 * t544 + t473) + ((0.8e1 * t153 * t218 * t487 + t115 * t529 - 0.24e2 * t329) * t298 + t151 * (0.12e2 * (-t114 * t488 - t416) * t403 + t95 * t401) + (t253 * t298 * t450 - t53 + t151 * (t132 * t429 + t441 * t70) + ((pkin(2) * t218 * t426 + 0.4e1 * t113) * t298 + t151 * (-0.48e2 * t122 * t494 - 0.6e1 * t80)) * pkin(1)) * pkin(3)) * t250 + ((t429 * t205 + ((-pkin(3) * t194 + t180 * t526) * t252 + (-t167 * t508 - t480) * t535) * t442) * t298 + t151 * ((t187 - 0.8e1 * t383) * t430 + ((-0.2e1 * t192 * t250 + t439 * t206) * t252 + (-0.4e1 * t120 * t508 - 0.8e1 * t345) * pkin(1)) * t443)) * t186;
t390 = -t57 * t65 / 0.8e1;
t393 = t56 * t523;
t380 = t121 * t522;
t50 = -t503 + t67 * t380 + t133 * t212 * t532 + t85 * t252 + (-t129 - 0.8e1 * t408) * t180;
t379 = t123 * t522;
t51 = t470 * t298 + t67 * t379 + (t111 * t185 + 0.4e1 * t409) * t250 + (t446 * t487 + (t170 * t252 + t215 * t438) * pkin(2)) * t251;
t505 = (t116 * t316 + (t24 * t395 + t50 * t516 + t67 * t390 + (t35 * t393 + t51 * t524) * t36) * t93) * t483 - t32;
t189 = pkin(1) * t445;
t414 = pkin(2) * t491;
t356 = pkin(3) * t414;
t471 = -0.2e1 * t356 + t189;
t112 = 0.2e1 * t355 + t471;
t126 = t201 * t340 * t511;
t410 = pkin(2) * t489;
t358 = pkin(1) * t410;
t323 = t291 * t250 * t358;
t142 = -0.4e1 * t323;
t163 = t358 * t447;
t493 = t218 * t302;
t386 = t186 * t493;
t486 = t253 * t295;
t385 = t251 * t486;
t137 = t189 - 0.4e1 * t385;
t415 = pkin(1) * t540;
t63 = t137 * t421 + (pkin(1) * t171 * t534 + t191 * t528) * t251 + (-t471 * t169 + 0.2e1 * t251 ^ 2 * t415 + (-t110 * t490 + t130 * t491) * pkin(2)) * t545;
t389 = t63 * t522;
t495 = t217 * t303;
t25 = t44 * t389 + (-0.32e2 * t142 * t151 + 0.8e1 * t163 * t298) * t387 + (0.24e2 * (0.4e1 * t128 * t495 * t508 - t386 * t75 + t61 * t486) * t298 + t151 * (0.96e2 * t386 * t71 - 0.48e2 * t59 * t486 - 0.64e2 * t89 * t495)) * t251 + ((t40 + (t151 * t531 - 0.6e1 * t298 * t62) * t186) * t251 + (((t144 * t420 + t146 * t377 + t333) * t404 + t147 * t348 + t104 * t422 - 0.6e1 * t91 * t433 + t326) * t298 + (t126 * t530 + (t539 * t429 + (t149 * t377 - 0.8e1 * t164 * t480 - 0.8e1 * t418 + t91) * t441) * t298) * t186 + ((-pkin(1) * t331 - t215 * t332 + t252 * t350 + t359) * t402 + t214 * t327 + t215 * t342 + t252 * t351 + t449 + ((t360 - 0.4e1 * t476) * t405 + (t252 * t349 + t440) * t444) * t186) * t151 * t250) * t253) * pkin(2);
t37 = (t163 * t428 + (t84 * t426 + t193 * t215 + ((t209 + 0.2e1 * t480) * t428 - t131 + (-0.4e1 * t196 * t205 + t431) * pkin(1)) * pkin(2)) * t253 + ((t163 + (t196 - 0.2e1 * t480) * t511) * t442 + (-0.24e2 * t493 * t508 + t450) * t251) * t186) * t298 + t58 * t389 + t53 * t206 + t151 * (0.24e2 * t132 * t251 * t386 + (t142 + (-0.8e1 / 0.3e1 * t476 - 0.2e1 * t506) * t414) * t403 - 0.24e2 * t69 * t385 + t126 + t323 * t542 + t356 * t544 + 0.6e1 * (-(pkin(1) * t425 + t439) * t218 * t540 + t70 * t206) * t186);
t47 = t63 * t380 + t137 * t250 * t436 + ((-t250 * t298 + t86) * t253 + (t252 * t298 + (t185 * t535 + t129) * t250 + (-pkin(3) + 0.2e1 * t509 + t514) * t445) * t251) * pkin(2);
t49 = (t180 + t410) * t298 + t63 * t379 - 0.2e1 * t137 * t509 - (t189 - 0.4e1 * t356) * t169 - 0.2e1 * t218 * t415 + t170 * t414 + (t486 * t533 + (-t111 * t252 + t172 * t535) * pkin(2)) * t251;
t504 = (-t112 * t316 + (t25 * t395 + t47 * t516 + t63 * t390 + (t37 * t393 + t49 * t524) * t36) * t93) * t483 + t32;
t392 = t57 * t523;
t391 = t56 * t65 / 0.8e1;
t382 = t112 * t520;
t381 = t116 * t520;
t287 = 0.1e1 / pkin(4);
t378 = t287 * t521;
t55 = 0.1e1 / t56 ^ 2;
t354 = pkin(4) * t287 / (t55 * t57 ^ 2 + 0.1e1) * t96;
t330 = 0.1e1 / t56 * t354;
t324 = t55 * t57 * t354;
t317 = t94 * t319;
t161 = rSges(2,1) * t253 - rSges(2,2) * t251;
t160 = rSges(4,1) * t252 - rSges(4,2) * t250;
t159 = -rSges(2,1) * t251 - rSges(2,2) * t253;
t158 = -rSges(4,1) * t250 - rSges(4,2) * t252;
t48 = atan2(t57 * t378, t56 * t378);
t46 = cos(t48);
t45 = sin(t48);
t39 = -t250 * t45 + t252 * t46;
t38 = t250 * t46 + t252 * t45;
t34 = rSges(5,1) * t39 - rSges(5,2) * t38;
t33 = rSges(5,1) * t38 + rSges(5,2) * t39;
t28 = 0.1e1 + 0.2e1 * (-t57 * t381 + t50 * t521) * t330 - 0.2e1 * (-t56 * t381 + t51 * t521) * t324;
t27 = 0.2e1 * (t57 * t382 + t47 * t521) * t330 - 0.2e1 * (t56 * t382 + t49 * t521) * t324;
t23 = t28 * t34 + t205;
t22 = -t28 * t33 - t508;
t18 = 0.1e1 / t21;
t17 = atan2(t20, t21);
t16 = 0.1e1 / (t19 * t20 ^ 2 + 0.1e1);
t15 = cos(t17);
t14 = sin(t17);
t13 = (-t112 * t317 + (t49 * t515 + t63 * t391 + t25 * t394 + (t37 * t392 + t47 * t524) * t36) * t93) * t483;
t11 = (t116 * t317 + (t51 * t515 + t67 * t391 + t24 * t394 + (t35 * t392 + t50 * t524) * t36) * t93) * t483;
t9 = t14 * t251 - t15 * t253;
t8 = -t14 * t253 - t15 * t251;
t7 = rSges(3,1) * t9 - rSges(3,2) * t8;
t6 = rSges(3,1) * t8 + rSges(3,2) * t9;
t5 = 0.1e1 + ((t504 * t139 + (-t13 + t31) * t138) * t18 - (-t13 * t139 - t504 * t138 + t29) * t507) * t16;
t4 = ((-t11 * t138 + t505 * t139 - t30) * t18 - ((-t11 - t31) * t139 - t505 * t138) * t507) * t16;
t2 = t5 * t7 + t511;
t1 = -t5 * t6 - t206;
t3 = [m(2) * (t159 ^ 2 + t161 ^ 2) + Icges(2,3) + (t1 ^ 2 + t2 ^ 2) * m(3) + t5 ^ 2 * Icges(3,3) + ((t33 ^ 2 + t34 ^ 2) * m(5) + Icges(5,3)) * t27 ^ 2; ((-t1 * t6 + t2 * t7) * m(3) + t5 * Icges(3,3)) * t4 + ((-t22 * t33 + t23 * t34) * m(5) + Icges(5,3) * t28) * t27; Icges(4,3) + t28 ^ 2 * Icges(5,3) + (t22 ^ 2 + t23 ^ 2) * m(5) + m(4) * (t158 ^ 2 + t160 ^ 2) + ((t6 ^ 2 + t7 ^ 2) * m(3) + Icges(3,3)) * t4 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t3(1), t3(2); t3(2), t3(3);];
Mq = res;