% Calculate time derivative of joint inertia matrix for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
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
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:33
% EndTime: 2020-04-12 19:34:24
% DurationCPUTime: 25.78s
% Computational Cost: add. (198398->514), mult. (280722->975), div. (11296->17), fcn. (78105->11), ass. (0->371)
t205 = pkin(2) ^ 2;
t206 = pkin(1) ^ 2;
t197 = cos(qJ(2));
t448 = pkin(2) * t197;
t370 = -0.2e1 * pkin(1) * t448 + t206;
t179 = t205 + t370;
t176 = 0.1e1 / t179 ^ 2;
t195 = sin(qJ(2));
t364 = qJD(2) * t195;
t352 = pkin(2) * t364;
t328 = pkin(1) * t352;
t309 = t176 * t328;
t261 = 0.2e1 * t309;
t201 = 0.1e1 / pkin(4);
t487 = pkin(4) ^ 2;
t202 = 0.1e1 / t487;
t368 = pkin(3) ^ 2 - t487;
t171 = t179 - t368;
t184 = pkin(1) - t448;
t459 = -pkin(3) - pkin(4);
t164 = (pkin(2) - t459) * (pkin(2) + t459) + t370;
t458 = pkin(4) - pkin(3);
t165 = (pkin(2) - t458) * (pkin(2) + t458) + t370;
t382 = t164 * t165;
t208 = sqrt(-t382);
t376 = t195 * t208;
t133 = -pkin(2) * t376 + t171 * t184;
t126 = t133 ^ 2;
t449 = pkin(2) * t195;
t162 = t171 * t449;
t134 = t184 * t208 + t162;
t130 = t134 ^ 2;
t373 = t126 + t130;
t112 = t373 * t202 * t176;
t110 = t112 ^ (-0.1e1 / 0.2e1);
t175 = 0.1e1 / t179;
t177 = t175 * t176;
t358 = 0.2e1 * t176;
t462 = pkin(1) * pkin(2);
t317 = (-t164 - t165) * t462;
t144 = t195 * t317;
t137 = qJD(2) * t144;
t150 = 0.1e1 / t208;
t386 = t137 * t150;
t344 = t195 * t386;
t359 = 0.2e1 * t184 * pkin(1);
t375 = t197 * t208;
t98 = (-t344 + (-t375 + (t171 + t359) * t195) * qJD(2)) * pkin(2);
t192 = t195 ^ 2;
t450 = pkin(1) * t205;
t355 = t192 * t450;
t319 = qJD(2) * t355;
t339 = t208 * t364;
t363 = qJD(2) * t197;
t372 = (t171 * t363 + t339) * pkin(2);
t384 = t150 * t184;
t99 = t137 * t384 + 0.2e1 * t319 + t372;
t351 = t110 / t112 * t175 * ((t133 * t98 + t134 * t99) * t358 - 0.4e1 * t373 * t177 * t328) * t202;
t389 = t110 * t201;
t494 = t389 * t261 + t201 * t351 / 0.2e1;
t423 = t134 * rSges(5,2);
t431 = t133 * rSges(5,1);
t303 = t423 + t431;
t346 = t175 * t389;
t240 = t303 * t346;
t493 = qJD(1) * t195;
t492 = qJD(1) * t197;
t128 = 0.1e1 / t133 ^ 2;
t387 = t128 * t134;
t127 = 0.1e1 / t133;
t106 = t144 * t384 + 0.2e1 * t355 + (t171 * t197 + t376) * pkin(2);
t379 = t176 * t195;
t329 = t379 * t462;
t455 = t175 / 0.2e1;
t91 = (t106 * t455 - t134 * t329) * t201;
t433 = t127 * t91;
t385 = t150 * t144;
t102 = t162 + (-t375 + (t359 - t385) * t195) * pkin(2);
t456 = -t175 / 0.2e1;
t89 = (t102 * t456 + t133 * t329) * t201;
t491 = -0.2e1 * t89 * t387 - 0.2e1 * t433;
t170 = t179 + t368;
t185 = pkin(1) * t197 - pkin(2);
t132 = -pkin(1) * t376 - t170 * t185;
t452 = pkin(1) * t170;
t163 = t195 * t452;
t135 = -t185 * t208 + t163;
t204 = 0.1e1 / pkin(3);
t380 = t175 * t204;
t109 = qJ(2) + atan2(t135 * t380, t132 * t380);
t108 = cos(t109);
t107 = sin(t109);
t407 = Icges(4,4) * t107;
t281 = Icges(4,2) * t108 + t407;
t378 = t192 * t206;
t340 = qJD(2) * t378;
t318 = pkin(2) * t340;
t371 = pkin(1) * t339 + t363 * t452;
t383 = t150 * t185;
t100 = -t137 * t383 + 0.2e1 * t318 + t371;
t123 = 0.1e1 / t132;
t288 = -0.2e1 * t309;
t124 = 0.1e1 / t132 ^ 2;
t131 = t135 ^ 2;
t118 = t124 * t131 + 0.1e1;
t115 = 0.1e1 / t118;
t447 = pkin(3) * t179;
t357 = t115 * t447;
t326 = t204 * t357;
t388 = t124 * t135;
t467 = -0.2e1 * t185;
t361 = pkin(2) * t467;
t97 = (-t344 + (-t375 + (t170 + t361) * t195) * qJD(2)) * pkin(1);
t49 = qJD(2) + ((t100 * t175 + t135 * t288) * t123 - (t132 * t288 + t175 * t97) * t388) * t326;
t253 = t281 * t49;
t406 = Icges(4,4) * t108;
t285 = Icges(4,1) * t107 + t406;
t254 = t285 * t49;
t105 = -t144 * t383 + 0.2e1 * pkin(2) * t378 + (t170 * t197 + t376) * pkin(1);
t315 = -0.2e1 * t329;
t92 = (t105 * t175 + t135 * t315) * t204;
t435 = t123 * t92;
t101 = t163 + (-t375 + (t361 - t385) * t195) * pkin(1);
t90 = (t101 * t175 + t132 * t315) * t204;
t490 = -t90 * t388 + t435;
t196 = sin(qJ(1));
t190 = t196 * rSges(4,3);
t198 = cos(qJ(1));
t391 = t107 * t198;
t374 = rSges(4,2) * t391 + t190;
t390 = t108 * t198;
t84 = -rSges(4,1) * t390 + t374;
t277 = Icges(4,5) * t107 + Icges(4,6) * t108;
t255 = t49 * t277;
t282 = -Icges(4,2) * t107 + t406;
t80 = Icges(4,6) * t196 - t198 * t282;
t286 = Icges(4,1) * t108 - t407;
t82 = Icges(4,5) * t196 - t198 * t286;
t301 = t107 * t80 - t108 * t82;
t471 = Icges(4,5) * t198 + t196 * t286;
t472 = Icges(4,6) * t198 + t196 * t282;
t302 = -t107 * t472 + t108 * t471;
t278 = Icges(4,5) * t108 - Icges(4,6) * t107;
t473 = Icges(4,3) * t198 + t196 * t278;
t457 = t490 * t357 + 0.1e1;
t50 = t457 * t196;
t51 = t457 * t198;
t78 = Icges(4,3) * t196 - t198 * t278;
t488 = qJD(1) * (t301 * t50 - t302 * t51) - (qJD(1) * t473 + t198 * t255) * t50 + (qJD(1) * t78 + t196 * t255) * t51;
t486 = pkin(1) - t240;
t399 = Icges(5,2) * t133;
t404 = Icges(5,4) * t134;
t72 = (-t399 + t404) * t346;
t418 = t134 * t72;
t405 = Icges(5,4) * t133;
t410 = Icges(5,1) * t134;
t73 = (-t405 + t410) * t346;
t425 = t133 * t73;
t295 = t418 + t425;
t241 = t295 * t346;
t411 = Icges(5,1) * t133;
t284 = t404 + t411;
t237 = t284 * t346;
t66 = -Icges(5,5) * t198 - t196 * t237;
t420 = t134 * t66;
t398 = Icges(5,2) * t134;
t280 = t398 + t405;
t236 = t280 * t346;
t64 = -Icges(5,6) * t198 - t196 * t236;
t429 = t133 * t64;
t298 = -t420 + t429;
t395 = Icges(5,6) * t133;
t400 = Icges(5,5) * t134;
t71 = (-t395 + t400) * t346;
t485 = -t196 * t241 - t198 * t71 - t298 * t346;
t67 = Icges(5,5) * t196 - t198 * t237;
t419 = t134 * t67;
t65 = Icges(5,6) * t196 - t198 * t236;
t428 = t133 * t65;
t296 = -t419 + t428;
t484 = t196 * t71 - t198 * t241 - t296 * t346;
t434 = t123 * t124 * t97;
t483 = 0.1e1 / t118 ^ 2 * (t100 * t388 - t131 * t434);
t304 = rSges(4,1) * t107 + rSges(4,2) * t108;
t482 = t304 * t49;
t408 = Icges(3,4) * t197;
t283 = -Icges(3,2) * t195 + t408;
t155 = Icges(3,6) * t196 + t198 * t283;
t409 = Icges(3,4) * t195;
t287 = Icges(3,1) * t197 - t409;
t157 = Icges(3,5) * t196 + t198 * t287;
t270 = t155 * t195 - t157 * t197;
t481 = t196 * t270;
t154 = -Icges(3,6) * t198 + t196 * t283;
t156 = -Icges(3,5) * t198 + t196 * t287;
t272 = t154 * t195 - t156 * t197;
t480 = t198 * t272;
t479 = (t431 / 0.2e1 + t423 / 0.2e1) * t196 * t351;
t478 = t303 * t196 * t261;
t365 = qJD(1) * t198;
t477 = -t195 * t365 - t196 * t363;
t300 = -t107 * t281 + t108 * t285;
t279 = Icges(3,5) * t197 - Icges(3,6) * t195;
t152 = -Icges(3,3) * t198 + t196 * t279;
t475 = qJD(1) * t152;
t474 = qJD(1) * t300 + t278 * t49;
t173 = Icges(3,2) * t197 + t409;
t174 = Icges(3,1) * t195 + t408;
t267 = t173 * t195 - t174 * t197;
t470 = qJD(1) * t267 + t279 * qJD(2);
t468 = 2 * m(3);
t466 = m(4) / 0.2e1;
t465 = m(5) / 0.2e1;
t451 = pkin(1) * t176;
t446 = pkin(4) * t179;
t117 = t128 * t130 + 0.1e1;
t113 = 0.1e1 / t117;
t136 = (t197 * t317 - 0.4e1 * t205 * t378) * qJD(2);
t316 = 0.4e1 * t137 / t382 * t385;
t289 = -t316 / 0.4e1;
t219 = (-0.2e1 * t197 * t386 + (-t136 * t150 + t289) * t195) * t175;
t290 = t177 * t205 * t340;
t360 = 0.2e1 * t446;
t310 = t134 * t89 * t360;
t323 = 0.6e1 * t195 * t363;
t345 = t175 * t386;
t356 = t113 * t446;
t381 = t175 * t197;
t432 = t127 / t126 * t98;
t5 = 0.2e1 * (t128 * t310 + t360 * t433) * (-t130 * t432 + t387 * t99) / t117 ^ 2 + 0.2e1 * ((-t89 * t99 + t91 * t98) * t128 + (-((t184 * t316 / 0.4e1 + t136 * t384 + t323 * t450) * t455 + 0.4e1 * t134 * t290 + ((t345 / 0.2e1 - t99 * t451) * t195 + ((t375 + (-t171 + t385) * t195) * t455 + (-t106 * t195 - t134 * t197) * t451) * qJD(2)) * pkin(2)) * t127 - ((0.4e1 * t319 + t372) * t456 - 0.4e1 * t133 * t290 + (-t219 / 0.2e1 + (t98 * t379 + (-t184 * t381 + (t102 * t195 + t133 * t197) * t176) * qJD(2)) * pkin(1)) * pkin(2)) * t387) * t201) * t356 + 0.2e1 * (pkin(4) * t328 * t491 + t310 * t432) * t113;
t445 = t196 * t5;
t444 = t198 * t5;
t443 = rSges(3,1) * t197;
t442 = rSges(4,1) * t108;
t441 = rSges(4,2) * t107;
t440 = rSges(3,3) * t196;
t439 = rSges(4,3) * t198;
t438 = rSges(5,3) * t198;
t430 = t133 * rSges(5,2);
t427 = t133 * t66;
t426 = t133 * t67;
t424 = t134 * rSges(5,1);
t422 = t134 * t64;
t421 = t134 * t65;
t53 = t356 * t491;
t417 = t175 * t53;
t189 = t196 * rSges(5,3);
t416 = t196 * t53;
t394 = Icges(5,6) * t134;
t401 = Icges(5,5) * t133;
t276 = t394 + t401;
t235 = t276 * t346;
t63 = Icges(5,3) * t196 - t198 * t235;
t415 = t196 * t63;
t414 = t198 * rSges(3,3);
t413 = t198 * t53;
t62 = -Icges(5,3) * t198 - t196 * t235;
t412 = t198 * t62;
t377 = t195 * t196;
t369 = t196 ^ 2 + t198 ^ 2;
t153 = Icges(3,3) * t196 + t198 * t279;
t367 = qJD(1) * t153;
t366 = qJD(1) * t196;
t362 = qJD(2) * t198;
t353 = t196 * t441;
t343 = t416 / 0.2e1;
t342 = -t413 / 0.2e1;
t341 = qJD(1) * t457;
t336 = t405 / 0.2e1;
t335 = qJD(1) * t472 + t198 * t253 + t49 * t82;
t334 = qJD(1) * t80 + t196 * t253 - t471 * t49;
t333 = -qJD(1) * t471 - t198 * t254 + t49 * t80;
t332 = -qJD(1) * t82 - t196 * t254 - t472 * t49;
t322 = rSges(4,3) * t365 + t366 * t442 + (rSges(4,1) * t391 + rSges(4,2) * t390) * t49;
t314 = 0.4e1 * t369;
t308 = -rSges(5,1) * t98 - rSges(5,2) * t99;
t306 = -rSges(3,2) * t195 + t443;
t178 = rSges(3,1) * t195 + rSges(3,2) * t197;
t305 = t441 - t442;
t299 = t422 + t427;
t297 = t421 + t426;
t294 = -t412 + t415;
t68 = -t196 * t240 - t438;
t69 = -t198 * t240 + t189;
t293 = t196 * t68 + t198 * t69;
t172 = Icges(3,5) * t195 + Icges(3,6) * t197;
t273 = t154 * t197 + t156 * t195;
t271 = t155 * t197 + t157 * t195;
t158 = -rSges(3,2) * t377 + t196 * t443 - t414;
t159 = t198 * t306 + t440;
t269 = t158 * t198 - t159 * t196;
t268 = t158 * t196 + t159 * t198;
t265 = 0.8e1 * t290;
t260 = t133 * t365 + t196 * t98;
t259 = t133 * t366 - t198 * t98;
t258 = t134 * t365 + t196 * t99;
t257 = t134 * t366 - t198 * t99;
t256 = t196 * t297;
t249 = qJD(2) * t172;
t60 = -t196 * t486 + t438;
t61 = t198 * t486 + t189;
t248 = -0.2e1 * t196 * t61 - 0.2e1 * t198 * t60;
t245 = -t305 - t448;
t243 = t299 * t346;
t242 = t297 * t346;
t239 = t133 * t261 - t175 * t98;
t238 = t134 * t288 + t175 * t99;
t234 = (t425 / 0.2e1 + t418 / 0.2e1) * t351;
t233 = t198 * t243;
t224 = t280 * t261;
t228 = (t336 + t398 / 0.2e1) * t351;
t13 = Icges(5,6) * t365 + (t198 * t228 + (t198 * t224 + (Icges(5,4) * t259 + Icges(5,2) * t257) * t175) * t110) * t201;
t225 = t284 * t261;
t229 = (t411 / 0.2e1 + t404 / 0.2e1) * t351;
t15 = Icges(5,5) * t365 + (t198 * t229 + (t198 * t225 + (Icges(5,1) * t259 + Icges(5,4) * t257) * t175) * t110) * t201;
t232 = -t13 * t134 - t133 * t15 - t65 * t99 - t67 * t98;
t14 = Icges(5,6) * t366 + (t196 * t228 + (t196 * t224 + (-Icges(5,4) * t260 - Icges(5,2) * t258) * t175) * t110) * t201;
t16 = Icges(5,5) * t366 + (t196 * t229 + (t196 * t225 + (-Icges(5,1) * t260 - Icges(5,4) * t258) * t175) * t110) * t201;
t231 = t133 * t16 + t134 * t14 + t64 * t99 + t66 * t98;
t18 = ((-t404 / 0.2e1 + t399 / 0.2e1) * t351 + (Icges(5,4) * t238 + Icges(5,2) * t239) * t110) * t201;
t19 = ((-t410 / 0.2e1 + t336) * t351 + (Icges(5,1) * t238 + Icges(5,4) * t239) * t110) * t201;
t230 = -t133 * t19 - t134 * t18 - t72 * t99 - t73 * t98;
t227 = (t401 / 0.2e1 + t394 / 0.2e1) * t351;
t226 = t295 * t261;
t223 = t276 * t261;
t30 = t282 * t49;
t31 = t286 * t49;
t222 = -qJD(1) * t277 + (t31 - t253) * t108 + (-t30 - t254) * t107;
t220 = t366 * t240 + rSges(5,3) * t365 + (t303 * t494 + t308 * t346) * t198;
t218 = ((t426 / 0.2e1 + t421 / 0.2e1) * t196 - (t427 / 0.2e1 + t422 / 0.2e1) * t198) * t53 * t351;
t217 = (-t198 * t299 + t256) * t53 * t261;
t216 = (t473 * t51 + t50 * t78) * qJD(1) + (-t332 * t51 + t333 * t50) * t108 + (-t334 * t51 + t335 * t50) * t107;
t169 = t306 * qJD(2);
t146 = -t196 * rSges(3,1) * t364 + (t198 * t443 + t440) * qJD(1) + t477 * rSges(3,2);
t145 = -t178 * t362 + (-t196 * t306 + t414) * qJD(1);
t139 = -t196 * t249 + t367;
t138 = -t198 * t249 - t475;
t122 = t153 * t196 - t198 * t270;
t121 = t152 * t196 - t480;
t120 = -t153 * t198 - t481;
t119 = -t152 * t198 - t196 * t272;
t83 = t196 * t305 - t439;
t76 = (-t442 + t448) * t198 + t374;
t75 = t196 * t245 + t439;
t74 = (t424 - t430) * t346;
t46 = t196 * t78 + t198 * t301;
t45 = -t196 * t473 + t198 * t302;
t44 = t196 * t301 - t198 * t78;
t43 = t196 * t302 + t198 * t473;
t42 = -t198 * t449 + t304 * t51;
t41 = -pkin(2) * t377 + t304 * t50;
t36 = -t198 * t242 + t415;
t35 = t196 * t62 - t233;
t34 = -t196 * t242 - t198 * t63;
t33 = -t196 * t243 - t412;
t32 = t305 * t49;
t22 = (t352 - t482) * t196 + (t198 * t245 - t190) * qJD(1);
t21 = -qJD(1) * t353 + (-t195 * t362 - t197 * t366) * pkin(2) + t322;
t20 = ((-t424 / 0.2e1 + t430 / 0.2e1) * t351 + (rSges(5,1) * t238 + rSges(5,2) * t239) * t110) * t201;
t17 = ((-t400 / 0.2e1 + t395 / 0.2e1) * t351 + (Icges(5,5) * t238 + Icges(5,6) * t239) * t110) * t201;
t12 = Icges(5,3) * t366 + (t196 * t227 + (t196 * t223 + (-Icges(5,5) * t260 - Icges(5,6) * t258) * t175) * t110) * t201;
t11 = Icges(5,3) * t365 + (t198 * t227 + (t198 * t223 + (Icges(5,5) * t259 + Icges(5,6) * t257) * t175) * t110) * t201;
t10 = (-pkin(1) * t198 - t189) * qJD(1) + (-t479 + (-t478 + (rSges(5,1) * t260 + rSges(5,2) * t258) * t175) * t110) * t201;
t9 = -pkin(1) * t366 + t220;
t7 = (t196 * t34 - t198 * t33) * t53;
t6 = ((pkin(2) * t206 * t323 - t136 * t383 + t185 * t289) * t175 + t135 * t265 + ((-0.2e1 * t100 * t176 * pkin(2) + t345) * t195 + ((t375 + (-t170 + t385) * t195) * t175 + (-t105 * t195 - t135 * t197) * pkin(2) * t358) * qJD(2)) * pkin(1)) * t123 * t326 + 0.2e1 * t490 * pkin(3) * t115 * t328 + (-t92 * t97 - ((0.4e1 * t318 + t371) * t175 + t132 * t265 + (t219 + (-0.2e1 * t97 * t379 + (t381 * t467 + (-t101 * t195 - t132 * t197) * t358) * qJD(2)) * pkin(2)) * pkin(1)) * t204 * t135 - t90 * t100) * t124 * t357 + (-0.2e1 * t435 * t483 + 0.2e1 * (t115 * t434 + t124 * t483) * t135 * t90) * t447;
t4 = t196 * t6 + t198 * t341;
t3 = t196 * t341 - t198 * t6;
t2 = pkin(2) * t477 + t304 * t4 - t32 * t50;
t1 = -t3 * t304 - t32 * t51 + (t195 * t366 - t197 * t362) * pkin(2);
t8 = [0.2e1 * m(5) * (t10 * t60 + t61 * t9) + (t145 * t159 + t146 * t158) * t468 + t107 * t31 + t108 * t30 + 0.2e1 * m(4) * (t21 * t76 + t22 * t75) + t300 * t49 + (-t133 * t18 + t134 * t19 - t72 * t98 + t73 * t99) * t346 + (t287 - t173) * t364 + (t174 + t283) * t363 + (t72 * t133 - t73 * t134) * t494; (t196 * t17 + t71 * t365) * t343 + (-t198 * t17 + t71 * t366) * t342 + 0.2e1 * (t1 * t75 + t2 * t76 + t21 * t41 + t22 * t42) * t466 + (t74 * t5 * t248 + (t20 * t248 + 0.2e1 * (-t10 * t198 - t196 * t9 + (t196 * t60 - t198 * t61) * qJD(1)) * t74) * t53) * t465 + m(3) * (t269 * t169 + (-qJD(1) * t268 - t145 * t196 + t146 * t198) * t178) + (t107 * t471 + t108 * t472 + t196 * t300 + t198 * t277) * t3 / 0.2e1 + (-t107 * t82 - t108 * t80 - t196 * t277 + t198 * t300) * t4 / 0.2e1 + (t107 * t333 - t108 * t335 - t196 * t474 + t222 * t198) * t50 / 0.2e1 - (t107 * t332 - t108 * t334 + t222 * t196 + t198 * t474) * t51 / 0.2e1 + (-qJD(2) * t270 - t154 * t492 - t156 * t493 + t196 * t470) * t196 / 0.2e1 - (-qJD(2) * t272 + t155 * t492 + t157 * t493 - t198 * t470) * t198 / 0.2e1 + t484 * t445 / 0.2e1 - t485 * t444 / 0.2e1 + ((t198 * t234 + (-t419 / 0.2e1 + t428 / 0.2e1) * t351 + (t198 * t226 + t296 * t261 + (-t13 * t133 + t134 * t15 + t198 * t230 + t295 * t366 - t65 * t98 + t67 * t99) * t175) * t110) * t343 + (t196 * t234 + (-t420 / 0.2e1 + t429 / 0.2e1) * t351 + (t196 * t226 + t298 * t261 + (-t133 * t14 + t134 * t16 + t196 * t230 - t295 * t365 - t64 * t98 + t66 * t99) * t175) * t110) * t342) * t201 + (-t172 * t198 - t196 * t267 + t485 * t53 + t273) * t366 / 0.2e1 + (t196 * t172 - t198 * t267 + t484 * t53 + t271) * t365 / 0.2e1; -t7 * t444 - ((-t33 * t5 + (t198 * t12 + (t34 + t233) * qJD(1)) * t53) * t198 + (t34 * t5 + (-t198 * t11 + (t294 + t33) * qJD(1)) * t53 + (t218 + (t217 + (t232 * t196 + (-qJD(1) * t297 + t231) * t198) * t417) * t110) * t201) * t196) * t413 + ((t36 * t5 + (t196 * t11 + (t256 * t346 + t35) * qJD(1)) * t53) * t196 + (-t35 * t5 + (-t196 * t12 + (t294 + t36) * qJD(1)) * t53 + (t218 + (t217 + (t231 * t198 + (-qJD(1) * t299 + t232) * t196) * t417) * t110) * t201) * t198) * t416 + (t268 * (qJD(1) * t269 + t198 * t145 + t196 * t146) + t369 * t178 * t169) * t468 + (-t121 * t198 + t122 * t196) * t365 + t196 * ((t196 * t138 + (t121 + t481) * qJD(1)) * t196 + (t122 * qJD(1) + (t154 * t363 + t156 * t364 - t475) * t198 + (-t139 - t271 * qJD(2) + (t153 - t272) * qJD(1)) * t196) * t198) - t198 * ((t198 * t139 + (t120 + t480) * qJD(1)) * t198 + (t119 * qJD(1) + (-t155 * t363 - t157 * t364 + t367) * t196 + (-t138 + t273 * qJD(2) + (-t152 - t270) * qJD(1)) * t198) * t196) + t4 * (-t45 * t51 + t46 * t50) + t50 * (-t488 * t196 + t216 * t198 + t45 * t3 + t46 * t4) + t3 * (-t43 * t51 + t44 * t50) - t51 * (t216 * t196 + t488 * t198 + t43 * t3 + t44 * t4) + ((t314 * t74 ^ 2 + 0.4e1 * t293 ^ 2) * t5 + (t74 * t20 * t314 + 0.4e1 * t293 * ((qJD(1) * t68 + t220) * t198 + ((-t69 + t189) * qJD(1) + (t479 + (t478 + (t196 * t308 - t303 * t365) * t175) * t110) * t201) * t196)) * t53) * t53 * t465 + 0.4e1 * (t42 * t1 + t41 * t2 + (t369 * t448 + t50 * t83 + t51 * t84) * (t4 * t83 - t3 * t84 + t51 * t322 + t50 * t196 * t482 + (-t51 * t353 + t50 * t84) * qJD(1) - t369 * t352)) * t466 + (t53 * t365 + t445) * (t196 * t36 - t198 * t35) * t53 + (-t119 * t198 + t120 * t196 + t53 * t7) * t366;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t8(1), t8(2); t8(2), t8(3);];
Mq = res;
