% Calculate joint inertia matrix for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1TE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_inertiaJ_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m1TE_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:48:55
% EndTime: 2020-04-17 16:14:56
% DurationCPUTime: 328.93s
% Computational Cost: add. (8634531->675), mult. (13015981->1116), div. (591512->22), fcn. (8246584->24), ass. (0->438)
t503 = sin(qJ(2));
t504 = sin(pkin(16));
t506 = cos(qJ(2));
t507 = cos(pkin(16));
t385 = t503 * t504 - t506 * t507;
t379 = pkin(5) * t385;
t363 = (-0.2e1 * t379 + pkin(1)) * pkin(1);
t538 = pkin(5) ^ 2;
t298 = t363 + t538;
t353 = pkin(2) ^ 2;
t464 = pkin(6) ^ 2 - t353;
t285 = t298 - t464;
t303 = -t379 + pkin(1);
t307 = t503 * t507 + t504 * t506;
t525 = pkin(5) + pkin(6);
t526 = pkin(5) - pkin(6);
t352 = sqrt(-((-pkin(2) + t525) * (pkin(2) + t526) + t363) * ((pkin(2) + t525) * (-pkin(2) + t526) + t363));
t373 = t385 * t352;
t498 = pkin(5) * t307;
t459 = pkin(1) * t498;
t471 = 0.4e1 / t352 * (t525 * t526 - t353 + t363) * t459;
t434 = -t471 / 0.2e1;
t236 = (t373 + (-0.2e1 * t303 * pkin(1) - t285 + t434) * t307) * pkin(5);
t460 = -0.2e1 * pkin(5) * t307 ^ 2;
t470 = t307 * t352;
t238 = t303 * t471 / 0.2e1 + (pkin(1) * t460 - t285 * t385 - t470) * pkin(5);
t260 = -pkin(5) * t470 + t303 * t285;
t261 = t285 * t498 + t303 * t352;
t350 = 0.1e1 / pkin(2);
t336 = cos(pkin(19));
t535 = 0.1e1 / t298;
t462 = -t535 / 0.2e1;
t432 = t336 * t462;
t334 = sin(pkin(19));
t461 = t535 / 0.2e1;
t433 = t334 * t461;
t243 = (t260 * t432 + t261 * t433) * t350;
t241 = 0.1e1 / t243 ^ 2;
t431 = t336 * t461;
t244 = (t260 * t433 + t261 * t431) * t350;
t296 = 0.1e1 / t298 ^ 2;
t425 = t296 * t459;
t407 = t336 * t425;
t408 = t334 * t425;
t549 = ((t236 * t433 + t238 * t431 + t260 * t408 + t261 * t407) / t243 - (t236 * t432 + t238 * t433 - t260 * t407 + t261 * t408) * t244 * t241) / (t241 * t244 ^ 2 + 0.1e1) * t350 + 0.1e1;
t340 = sin(qJ(1));
t330 = t340 ^ 2;
t343 = cos(qJ(1));
t331 = t343 ^ 2;
t548 = t330 + t331;
t362 = t298 + t464;
t361 = pkin(1) * t362;
t365 = pkin(1) * t385 - pkin(5);
t537 = 0.1e1 / pkin(6);
t360 = t537 * (t307 * t361 - t352 * t365);
t357 = t360 * t461;
t457 = pkin(1) * t470;
t359 = t537 * (-t362 * t365 - t457);
t358 = t535 * t359;
t505 = sin(pkin(15));
t508 = cos(pkin(15));
t245 = t508 * t358 / 0.2e1 + t505 * t357;
t248 = t508 * t357 - t505 * t358 / 0.2e1;
t399 = rSges(7,1) * t245 - rSges(7,2) * t248;
t547 = 2 * pkin(3);
t349 = pkin(3) ^ 2;
t348 = pkin(4) ^ 2;
t342 = cos(qJ(3));
t429 = t342 * t462;
t339 = sin(qJ(3));
t430 = t339 * t461;
t246 = (t260 * t429 + t261 * t430) * t350;
t428 = t342 * t461;
t247 = (t260 * t430 + t261 * t428) * t350;
t329 = pkin(18) + pkin(19);
t320 = sin(t329);
t321 = cos(t329);
t221 = t246 * t321 + t247 * t320;
t500 = pkin(4) * t221;
t469 = t500 * t547 + t348;
t203 = t349 + t469;
t465 = pkin(8) ^ 2 - pkin(10) ^ 2;
t198 = t203 - t465;
t206 = pkin(3) * t221 + pkin(4);
t524 = (-pkin(8) - pkin(10));
t195 = ((pkin(3) - t524) * (pkin(3) + t524)) + t469;
t523 = (pkin(10) - pkin(8));
t196 = ((pkin(3) - t523) * (pkin(3) + t523)) + t469;
t351 = sqrt(-t196 * t195);
t390 = t246 * t320 - t247 * t321;
t541 = t390 * t351;
t172 = -pkin(3) * t541 + t198 * t206;
t174 = pkin(3) * t198 * t390 + t206 * t351;
t333 = cos(pkin(17));
t199 = 0.1e1 / t203;
t345 = 0.1e1 / pkin(10);
t476 = t199 * t345;
t332 = sin(pkin(17));
t513 = t332 / 0.2e1;
t151 = (-t172 * t333 / 0.2e1 + t174 * t513) * t476;
t152 = (t174 * t333 / 0.2e1 + t172 * t513) * t476;
t305 = t339 * t503 - t342 * t506;
t306 = -t339 * t506 - t342 * t503;
t142 = t151 * t306 + t152 * t305;
t143 = t151 * t305 - t152 * t306;
t114 = Icges(5,4) * t142 + Icges(5,2) * t143;
t546 = t114 / 0.2e1;
t115 = Icges(5,1) * t142 + Icges(5,4) * t143;
t545 = t115 / 0.2e1;
t544 = t142 / 0.2e1;
t543 = t143 / 0.2e1;
t511 = t340 / 0.2e1;
t510 = -t343 / 0.2e1;
t409 = 0.2e1 * t548;
t389 = t409 / 0.2e1;
t542 = t340 * t343;
t301 = t306 * t343;
t302 = t305 * t343;
t140 = t151 * t302 - t152 * t301;
t141 = -t151 * t301 - t152 * t302;
t106 = Icges(5,4) * t140 + Icges(5,2) * t141 + Icges(5,6) * t340;
t108 = Icges(5,1) * t140 + Icges(5,4) * t141 + Icges(5,5) * t340;
t113 = Icges(5,5) * t142 + Icges(5,6) * t143;
t540 = t106 * t543 + t108 * t544 + t113 * t511 + t140 * t545 + t141 * t546;
t299 = t306 * t340;
t300 = t305 * t340;
t138 = t151 * t300 - t152 * t299;
t139 = -t151 * t299 - t152 * t300;
t105 = Icges(5,4) * t138 + Icges(5,2) * t139 - Icges(5,6) * t343;
t107 = Icges(5,1) * t138 + Icges(5,4) * t139 - Icges(5,5) * t343;
t539 = t105 * t543 + t107 * t544 + t113 * t510 + t138 * t545 + t139 * t546;
t193 = (t238 * t428 + t236 * t430 + (t260 * t339 + t261 * t342) * t425) * t350;
t194 = (t236 * t429 + t238 * t430 + (-t260 * t342 + t261 * t339) * t425) * t350;
t177 = -t193 * t320 - t194 * t321;
t388 = pkin(4) * (t195 + t196) * t547;
t158 = t177 * t388;
t176 = -t193 * t321 + t194 * t320;
t179 = 0.1e1 / t351;
t515 = -t179 / 0.2e1;
t438 = t390 * t515;
t381 = t158 * t438 - t176 * t351;
t426 = -0.2e1 * pkin(4) * t206 - t198;
t126 = (t177 * t426 + t381) * pkin(3);
t501 = pkin(4) * t390;
t427 = -0.2e1 * t349 * t501;
t439 = t179 * t206 / 0.2e1;
t478 = t177 * t351;
t127 = t158 * t439 + t177 * t427 + (t176 * t198 - t478) * pkin(3);
t149 = 0.1e1 / t151;
t200 = 0.1e1 / t203 ^ 2;
t502 = pkin(4) * t200;
t458 = pkin(3) * t502;
t423 = t333 * t458;
t405 = t177 * t423;
t424 = t332 * t458;
t406 = t177 * t424;
t477 = t199 * t333;
t435 = t477 / 0.2e1;
t436 = -t477 / 0.2e1;
t437 = t199 * t513;
t150 = 0.1e1 / t151 ^ 2;
t479 = t150 * t152;
t480 = 0.1e1 / (t150 * t152 ^ 2 + 0.1e1) * t345;
t522 = ((t126 * t437 + t127 * t435 + t172 * t406 + t174 * t405) * t149 - (t126 * t436 + t127 * t437 - t172 * t405 + t174 * t406) * t479) * t480 + 0.1e1;
t36 = t522 * t340;
t534 = t36 / 0.2e1;
t37 = t522 * t343;
t533 = -t37 / 0.2e1;
t159 = t549 * t340;
t197 = t203 + t465;
t205 = -pkin(3) - t500;
t173 = t197 * t501 - t205 * t351;
t366 = pkin(3) * (-t199 * t348 * t390 + t173 * t502);
t171 = -pkin(4) * t541 - t197 * t205;
t170 = 0.1e1 / t171 ^ 2;
t347 = 0.1e1 / pkin(8);
t415 = pkin(8) / (t170 * t173 ^ 2 + 0.1e1) * t203 * t347;
t383 = pkin(4) * t170 * t173 * t415;
t384 = pkin(3) * (t171 * t200 + t199 * t205);
t401 = 0.1e1 / t171 * t415;
t440 = t205 * t515;
t514 = t199 / 0.2e1;
t69 = 0.2e1 * ((t158 * t440 + (t176 * t197 - t478) * pkin(4)) * t514 + t177 * t366) * t401 - 0.2e1 * ((-t177 * t197 + t381) * t514 + t177 * t384) * t383;
t67 = t340 * t69 + t159;
t528 = t67 / 0.2e1;
t68 = (-t69 - t549) * t343;
t527 = t68 / 0.2e1;
t175 = t390 * t388;
t380 = t175 * t438 - t221 * t351;
t147 = (t390 * t426 + t380) * pkin(3);
t148 = t175 * t439 + t390 * t427 + (t198 * t221 - t541) * pkin(3);
t403 = t390 * t423;
t404 = t390 * t424;
t521 = ((t147 * t437 + t148 * t435 + t172 * t404 + t174 * t403) * t149 - (t147 * t436 + t148 * t437 - t172 * t403 + t174 * t404) * t479) * t480 + 0.1e1;
t520 = -t139 / 0.2e1;
t519 = -t141 / 0.2e1;
t518 = -t143 / 0.2e1;
t337 = cos(pkin(18));
t512 = -t337 / 0.2e1;
t499 = pkin(4) * t306;
t282 = rSges(4,1) * t306 + rSges(4,2) * t305;
t497 = m(4) * t282;
t335 = sin(pkin(18));
t475 = t199 * t347;
t153 = (t335 * t171 / 0.2e1 + t173 * t512) * t475;
t154 = (t171 * t512 - t335 * t173 / 0.2e1) * t475;
t224 = t243 * t503 + t244 * t506;
t225 = t243 * t506 - t244 * t503;
t124 = t153 * t225 + t154 * t224;
t125 = -t153 * t224 + t154 * t225;
t100 = rSges(9,1) * t124 + rSges(9,2) * t125;
t496 = m(9) * t100;
t495 = t138 * pkin(9);
t338 = sin(qJ(4));
t341 = cos(qJ(4));
t130 = -t138 * t338 - t341 * t343;
t131 = t138 * t341 - t338 * t343;
t74 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t139;
t76 = Icges(6,4) * t131 + Icges(6,2) * t130 - Icges(6,6) * t139;
t78 = Icges(6,1) * t131 + Icges(6,4) * t130 - Icges(6,5) * t139;
t29 = -t143 * t74 + (-t338 * t76 + t341 * t78) * t142;
t494 = t29 * t37;
t132 = -t140 * t338 + t340 * t341;
t133 = t140 * t341 + t338 * t340;
t75 = Icges(6,5) * t133 + Icges(6,6) * t132 - Icges(6,3) * t141;
t77 = Icges(6,4) * t133 + Icges(6,2) * t132 - Icges(6,6) * t141;
t79 = Icges(6,1) * t133 + Icges(6,4) * t132 - Icges(6,5) * t141;
t30 = -t143 * t75 + (-t338 * t77 + t341 * t79) * t142;
t493 = t30 * t36;
t492 = t300 * pkin(4);
t297 = t302 * pkin(4);
t82 = -Icges(6,3) * t143 + (Icges(6,5) * t341 - Icges(6,6) * t338) * t142;
t84 = -Icges(6,5) * t143 + (Icges(6,1) * t341 - Icges(6,4) * t338) * t142;
t491 = t142 * t341 * t84 - t143 * t82;
t83 = -Icges(6,6) * t143 + (Icges(6,4) * t341 - Icges(6,2) * t338) * t142;
t488 = t338 * t83;
t487 = t343 * rSges(3,3);
t486 = t343 * rSges(7,3);
t400 = -t131 * rSges(6,1) - t130 * rSges(6,2);
t80 = -t139 * rSges(6,3) - t400;
t485 = -t139 * pkin(11) + t495 + t80;
t81 = t133 * rSges(6,1) + t132 * rSges(6,2) - t141 * rSges(6,3);
t484 = -t140 * pkin(9) + pkin(11) * t141 - t81;
t85 = -rSges(6,3) * t143 + (rSges(6,1) * t341 - rSges(6,2) * t338) * t142;
t483 = pkin(9) * t142 - pkin(11) * t143 + t85;
t482 = Icges(7,4) * t245;
t481 = Icges(7,4) * t248;
t219 = t224 * t343;
t474 = t219 * t335;
t375 = t300 * rSges(4,1) - t299 * rSges(4,2) - t343 * rSges(4,3);
t448 = t302 * rSges(4,1) - t301 * rSges(4,2) + t340 * rSges(4,3);
t262 = t340 * t375 + t343 * t448;
t468 = t343 * t297 + t340 * t492;
t454 = t506 * pkin(1);
t467 = t548 * t454;
t328 = t343 * pkin(13);
t444 = t343 * t506;
t466 = pkin(1) * t444 + t328;
t456 = t340 * t499;
t455 = t343 * t499;
t453 = t503 * pkin(1);
t31 = t130 * t83 + t131 * t84 - t139 * t82;
t452 = -t29 / 0.2e1 - t31 / 0.2e1;
t32 = t132 * t83 + t133 * t84 - t141 * t82;
t451 = -t32 / 0.2e1 - t30 / 0.2e1;
t218 = t225 * t343;
t122 = -t153 * t219 + t154 * t218;
t123 = -t153 * t218 - t154 * t219;
t95 = t122 * rSges(9,1) + t123 * rSges(9,2) + t340 * rSges(9,3);
t110 = t140 * rSges(5,1) + t141 * rSges(5,2) + t340 * rSges(5,3);
t450 = t218 * rSges(8,1) - t219 * rSges(8,2) + t340 * rSges(8,3);
t449 = t340 * rSges(7,3) + t399 * t343;
t447 = t297 + t466;
t446 = Icges(3,4) * t506;
t445 = Icges(3,4) * t503;
t443 = t343 * t503;
t414 = t537 * t535 * t505;
t413 = t340 * t453;
t412 = pkin(1) * t443;
t272 = Icges(4,4) * t300 - Icges(4,2) * t299 - Icges(4,6) * t343;
t273 = Icges(4,4) * t302 - Icges(4,2) * t301 + Icges(4,6) * t340;
t274 = Icges(4,1) * t300 - Icges(4,4) * t299 - Icges(4,5) * t343;
t275 = Icges(4,1) * t302 - Icges(4,4) * t301 + Icges(4,5) * t340;
t279 = Icges(4,5) * t306 + Icges(4,6) * t305;
t280 = Icges(4,4) * t306 + Icges(4,2) * t305;
t281 = Icges(4,1) * t306 + Icges(4,4) * t305;
t411 = (t273 * t305 + t275 * t306 + t279 * t340 - t280 * t301 + t281 * t302) * t511 + (t272 * t305 + t274 * t306 - t279 * t343 - t280 * t299 + t281 * t300) * t510;
t410 = t467 + t468;
t402 = -t453 - t282;
t398 = t537 * t508 * t461;
t397 = Icges(7,1) * t245 - t481;
t396 = -Icges(7,2) * t248 + t482;
t395 = Icges(7,5) * t245 - Icges(7,6) * t248;
t216 = t225 * t340;
t217 = t224 * t340;
t392 = t216 * t337 + t217 * t335;
t228 = Icges(7,2) * t245 + t481;
t229 = Icges(7,1) * t248 + t482;
t391 = -t228 * t248 + t229 * t245;
t387 = -t453 - t499;
t386 = (-t454 - pkin(13)) * t340;
t382 = rSges(3,1) * t444 - rSges(3,2) * t443 + t340 * rSges(3,3);
t378 = rSges(3,1) * t506 - rSges(3,2) * t503;
t377 = t387 * t340;
t376 = t387 * t343;
t109 = t138 * rSges(5,1) + t139 * rSges(5,2) - t343 * rSges(5,3);
t374 = t216 * rSges(8,1) - t217 * rSges(8,2) - t343 * rSges(8,3);
t120 = -t153 * t217 + t154 * t216;
t121 = -t153 * t216 - t154 * t217;
t94 = t120 * rSges(9,1) + t121 * rSges(9,2) - t343 * rSges(9,3);
t372 = Icges(3,1) * t506 - t445;
t371 = -Icges(3,2) * t503 + t446;
t370 = Icges(3,5) * t506 - Icges(3,6) * t503;
t310 = Icges(3,2) * t506 + t445;
t311 = Icges(3,1) * t503 + t446;
t367 = -t310 * t503 + t311 * t506;
t364 = t386 - t492;
t356 = t360 * t459;
t355 = t359 * t459;
t314 = rSges(2,1) * t343 - rSges(2,2) * t340;
t313 = -rSges(2,1) * t340 - rSges(2,2) * t343;
t312 = rSges(3,1) * t503 + rSges(3,2) * t506;
t309 = Icges(3,5) * t503 + Icges(3,6) * t506;
t289 = Icges(3,3) * t340 + t343 * t370;
t288 = -Icges(3,3) * t343 + t340 * t370;
t284 = t328 + t382;
t283 = t487 + (-pkin(13) - t378) * t340;
t277 = t402 * t343;
t276 = t402 * t340;
t271 = Icges(4,5) * t302 - Icges(4,6) * t301 + Icges(4,3) * t340;
t270 = Icges(4,5) * t300 - Icges(4,6) * t299 - Icges(4,3) * t343;
t269 = t343 * t382 + (t340 * t378 - t487) * t340;
t265 = t448 + t466;
t264 = t386 - t375;
t259 = t467 + t262;
t242 = 0.1e1 / t245 ^ 2;
t237 = pkin(1) ^ 2 * t460 - t361 * t385 + t365 * t434 - t457;
t235 = (t373 + (t434 - 0.3e1 * t538 + (0.4e1 * t379 - pkin(1)) * pkin(1) - t464) * t307) * pkin(1);
t234 = (-t271 * t343 - t273 * t299 + t275 * t300) * t340 - (-t270 * t343 - t272 * t299 + t274 * t300) * t343;
t233 = t340 * ((t271 * t340 - t273 * t301 + t275 * t302) * t340 - (t270 * t340 - t272 * t301 + t274 * t302) * t343);
t230 = rSges(7,1) * t248 + rSges(7,2) * t245;
t227 = Icges(7,5) * t248 + Icges(7,6) * t245;
t211 = Icges(7,3) * t340 + t343 * t395;
t210 = -Icges(7,3) * t343 + t340 * t395;
t204 = t218 * t337 * pkin(3);
t202 = -pkin(7) * t343 + t449;
t201 = t486 + (pkin(7) - t399) * t340;
t192 = rSges(8,1) * t224 + rSges(8,2) * t225;
t191 = Icges(8,1) * t224 + Icges(8,4) * t225;
t190 = Icges(8,4) * t224 + Icges(8,2) * t225;
t189 = Icges(8,5) * t224 + Icges(8,6) * t225;
t188 = (t224 * t337 - t225 * t335) * pkin(3);
t187 = Icges(8,1) * t218 - Icges(8,4) * t219 + Icges(8,5) * t340;
t186 = Icges(8,1) * t216 - Icges(8,4) * t217 - Icges(8,5) * t343;
t185 = Icges(8,4) * t218 - Icges(8,2) * t219 + Icges(8,6) * t340;
t184 = Icges(8,4) * t216 - Icges(8,2) * t217 - Icges(8,6) * t343;
t183 = Icges(8,5) * t218 - Icges(8,6) * t219 + Icges(8,3) * t340;
t182 = Icges(8,5) * t216 - Icges(8,6) * t217 - Icges(8,3) * t343;
t181 = t450 + t466;
t180 = t386 - t374;
t164 = ((t237 * t398 - t235 * t414 / 0.2e1 + (-t355 * t505 + t356 * t508) * t296) / t245 - (t235 * t398 + t237 * t414 / 0.2e1 + (t355 * t508 + t356 * t505) * t296) * t248 * t242) / (t242 * t248 ^ 2 + 0.1e1);
t160 = t549 * t343;
t157 = -t160 * t192 - t412;
t156 = -t159 * t192 - t413;
t146 = (t343 * t449 + (t340 * t399 - t486) * t340) * t164;
t145 = t159 * t374 + t160 * t450 + t467;
t116 = rSges(5,1) * t142 + rSges(5,2) * t143;
t104 = Icges(5,5) * t140 + Icges(5,6) * t141 + Icges(5,3) * t340;
t103 = Icges(5,5) * t138 + Icges(5,6) * t139 - Icges(5,3) * t343;
t102 = t447 + t110;
t101 = -t109 + t364;
t99 = Icges(9,1) * t124 + Icges(9,4) * t125;
t98 = Icges(9,4) * t124 + Icges(9,2) * t125;
t97 = Icges(9,5) * t124 + Icges(9,6) * t125;
t93 = Icges(9,1) * t122 + Icges(9,4) * t123 + Icges(9,5) * t340;
t92 = Icges(9,1) * t120 + Icges(9,4) * t121 - Icges(9,5) * t343;
t91 = Icges(9,4) * t122 + Icges(9,2) * t123 + Icges(9,6) * t340;
t90 = Icges(9,4) * t120 + Icges(9,2) * t121 - Icges(9,6) * t343;
t89 = Icges(9,5) * t122 + Icges(9,6) * t123 + Icges(9,3) * t340;
t88 = Icges(9,5) * t120 + Icges(9,6) * t121 - Icges(9,3) * t343;
t87 = pkin(3) * t474 + t204 + t466 + t95;
t86 = -pkin(3) * t392 + t386 - t94;
t73 = 0.2e1 * ((t175 * t440 + (t197 * t221 - t541) * pkin(4)) * t514 + t390 * t366) * t401 - 0.2e1 * ((-t197 * t390 + t380) * t514 + t390 * t384) * t383;
t66 = t447 - t484;
t65 = -t495 + (pkin(11) + rSges(6,3)) * t139 + t364 + t400;
t62 = t124 * t93 + t125 * t91;
t61 = t124 * t92 + t125 * t90;
t59 = t521 * t343;
t58 = t521 * t340;
t57 = t68 * t100 - t160 * t188 - t412;
t56 = -t67 * t100 - t159 * t188 - t413;
t53 = t141 * t85 - t143 * t81;
t52 = -t139 * t85 + t143 * t80;
t51 = t104 * t340 + t106 * t141 + t108 * t140;
t50 = t103 * t340 + t105 * t141 + t107 * t140;
t49 = -t104 * t343 + t106 * t139 + t108 * t138;
t48 = -t103 * t343 + t105 * t139 + t107 * t138;
t47 = t139 * t81 - t141 * t80;
t46 = t122 * t99 + t123 * t98 + t340 * t97;
t45 = t120 * t99 + t121 * t98 - t343 * t97;
t44 = t122 * t93 + t123 * t91 + t340 * t89;
t43 = t122 * t92 + t123 * t90 + t340 * t88;
t42 = t120 * t93 + t121 * t91 - t343 * t89;
t41 = t120 * t92 + t121 * t90 - t343 * t88;
t40 = -t116 * t59 - t455;
t39 = -t116 * t58 - t456;
t35 = (t340 * t94 + t343 * t95) * t73;
t34 = -t37 * t116 + t376;
t33 = -t36 * t116 + t377;
t28 = t132 * t77 + t133 * t79 - t141 * t75;
t27 = t132 * t76 + t133 * t78 - t141 * t74;
t26 = t130 * t77 + t131 * t79 - t139 * t75;
t25 = t130 * t76 + t131 * t78 - t139 * t74;
t24 = (-t142 * t488 + t491) * t143;
t23 = t160 * t204 + t67 * t94 - t68 * t95 + (t159 * t392 + t160 * t474) * pkin(3) + t467;
t22 = t109 * t58 + t110 * t59 + t468;
t21 = -t483 * t59 - t455;
t20 = -t483 * t58 - t456;
t19 = t109 * t36 + t110 * t37 + t410;
t18 = -t37 * t483 + t376;
t17 = -t36 * t483 + t377;
t16 = (t340 * t44 - t343 * t43) * t73;
t15 = (t340 * t42 - t343 * t41) * t73;
t14 = t43 * t68 + t44 * t67;
t13 = t41 * t68 + t42 * t67;
t12 = -t50 * t59 + t51 * t58;
t11 = -t48 * t59 + t49 * t58;
t10 = t36 * t51 - t37 * t50;
t9 = t36 * t49 - t37 * t48;
t8 = -t27 * t59 + t28 * t58;
t7 = -t25 * t59 + t26 * t58;
t6 = -t484 * t59 + t485 * t58 + t468;
t5 = -t139 * t27 - t141 * t28 - t143 * t32;
t4 = -t139 * t25 - t141 * t26 - t143 * t31;
t3 = -t27 * t37 + t28 * t36;
t2 = -t25 * t37 + t26 * t36;
t1 = t36 * t485 - t37 * t484 + t410;
t38 = [(t86 ^ 2 + t87 ^ 2) * m(9) + (t180 ^ 2 + t181 ^ 2) * m(8) + (t201 ^ 2 + t202 ^ 2) * m(7) + t491 + (t101 ^ 2 + t102 ^ 2) * m(5) + (t115 - t488) * t142 + (t65 ^ 2 + t66 ^ 2) * m(6) + t125 * t98 + t124 * t99 + t143 * t114 + t506 * t310 + t503 * t311 + t305 * t280 + t306 * t281 + t248 * t229 + t245 * t228 + (t264 ^ 2 + t265 ^ 2) * m(4) + (t283 ^ 2 + t284 ^ 2) * m(3) + m(2) * (t313 ^ 2 + t314 ^ 2) + t224 * t191 + t225 * t190 + Icges(2,3); (-t283 * t343 - t284 * t340) * t312 * m(3) + (t56 * t87 + t57 * t86) * m(9) + (t156 * t181 + t157 * t180) * m(8) + (t101 * t34 + t102 * t33) * m(5) + (t264 * t277 + t265 * t276) * m(4) + t493 / 0.2e1 - t494 / 0.2e1 + (t17 * t66 + t18 * t65) * m(6) + t411 + t31 * t533 + t32 * t534 + (t46 + t62) * t528 + (t45 + t61) * t527 + (t185 * t225 + t187 * t224 + t189 * t340 - t190 * t219 + t191 * t218) * t159 / 0.2e1 - (t184 * t225 + t186 * t224 - t189 * t343 - t190 * t217 + t191 * t216) * t160 / 0.2e1 + ((Icges(3,6) * t340 + t343 * t371) * t506 + (Icges(3,5) * t340 + t343 * t372) * t503 + t340 * t309 + t343 * t367) * t511 + ((-Icges(3,6) * t343 + t340 * t371) * t506 + (-Icges(3,5) * t343 + t340 * t372) * t503 - t343 * t309 + t340 * t367) * t510 - t539 * t37 + t540 * t36 + ((-t201 * t343 - t202 * t340) * t230 * m(7) + ((Icges(7,6) * t340 + t343 * t396) * t245 + (Icges(7,5) * t340 + t343 * t397) * t248 + t227 * t340 + t343 * t391) * t511 + ((-Icges(7,6) * t343 + t340 * t396) * t245 + (-Icges(7,5) * t343 + t340 * t397) * t248 - t227 * t343 + t340 * t391) * t510) * t164; (t145 ^ 2 + t156 ^ 2 + t157 ^ 2) * m(8) + m(3) * (t312 ^ 2 * t409 + 0.2e1 * t269 ^ 2) / 0.2e1 - t343 * t234 + t67 * t14 - t343 * (t331 * t288 - t289 * t542) + t159 * ((t183 * t340 - t185 * t219 + t187 * t218) * t159 - (t182 * t340 - t184 * t219 + t186 * t218) * t160) + t340 * (-t288 * t542 + t330 * t289) - t160 * ((-t183 * t343 - t185 * t217 + t187 * t216) * t159 - (-t182 * t343 - t184 * t217 + t186 * t216) * t160) + m(7) * t146 ^ 2 + (t23 ^ 2 + t56 ^ 2 + t57 ^ 2) * m(9) + (t259 ^ 2 + t276 ^ 2 + t277 ^ 2) * m(4) + t68 * t13 + (t19 ^ 2 + t33 ^ 2 + t34 ^ 2) * m(5) + (t1 ^ 2 + t17 ^ 2 + t18 ^ 2) * m(6) + t233 - (t9 + t2) * t37 + (t3 + t10) * t36 + (m(7) * t230 ^ 2 * t389 - t343 * (t331 * t210 - t211 * t542) + t340 * (-t210 * t542 + t330 * t211)) * t164 ^ 2; (t101 * t40 + t102 * t39) * m(5) + (t20 * t66 + t21 * t65) * m(6) - (-t452 + t539) * t59 + (-t451 + t540) * t58 + (-t264 * t497 + (-t86 * t496 - t61 / 0.2e1 - t45 / 0.2e1) * t73) * t343 + (-t265 * t497 + (-t87 * t496 + t62 / 0.2e1 + t46 / 0.2e1) * t73) * t340 + t411; t233 + t16 * t528 + m(9) * t23 * t35 + m(4) * t259 * t262 + t15 * t527 + (t19 * t22 + t33 * t39 + t34 * t40) * m(5) + (t1 * t6 + t17 * t20 + t18 * t21) * m(6) - (t2 / 0.2e1 + t9 / 0.2e1) * t59 + (t10 / 0.2e1 + t3 / 0.2e1) * t58 - (t7 / 0.2e1 + t11 / 0.2e1) * t37 + (t12 / 0.2e1 + t8 / 0.2e1) * t36 + (-t276 * t497 + (t14 / 0.2e1 - t56 * t496) * t73) * t340 + (-t277 * t497 - t234 + (-t57 * t496 - t13 / 0.2e1) * t73) * t343; t73 * t340 * t16 + (t22 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(5) + (t20 ^ 2 + t21 ^ 2 + t6 ^ 2) * m(6) + t233 - (t7 + t11) * t59 + (t8 + t12) * t58 + (-t15 * t73 - t234) * t343 + (t100 ^ 2 * t389 * t73 ^ 2 + t35 ^ 2) * m(9) + (t282 ^ 2 * t389 + t262 ^ 2) * m(4); -t24 + (t52 * t65 + t53 * t66) * m(6) + t451 * t141 + t452 * t139; (t493 - t494) * t518 + t3 * t519 + t5 * t534 + t2 * t520 + t4 * t533 + (t1 * t47 + t17 * t53 + t18 * t52) * m(6); -t59 * t4 / 0.2e1 + t8 * t519 + (-t29 * t59 + t30 * t58) * t518 + t58 * t5 / 0.2e1 + t7 * t520 + (t20 * t53 + t21 * t52 + t47 * t6) * m(6); -t141 * t5 - t143 * (-t29 * t139 - t30 * t141 - t24) + (t47 ^ 2 + t52 ^ 2 + t53 ^ 2) * m(6) - t139 * t4;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t38(1), t38(2), t38(4), t38(7); t38(2), t38(3), t38(5), t38(8); t38(4), t38(5), t38(6), t38(9); t38(7), t38(8), t38(9), t38(10);];
Mq = res;
