% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [2x2]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnDE2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:35
% EndTime: 2020-04-12 19:34:18
% DurationCPUTime: 19.90s
% Computational Cost: add. (278736->440), mult. (386648->835), div. (17184->17), fcn. (111507->11), ass. (0->345)
t229 = sin(qJ(2));
t239 = pkin(1) ^ 2;
t231 = cos(qJ(2));
t430 = pkin(2) * t231;
t351 = -0.2e1 * pkin(1) * t430 + t239;
t444 = -pkin(3) - pkin(4);
t194 = (pkin(2) - t444) * (pkin(2) + t444) + t351;
t443 = pkin(4) - pkin(3);
t195 = (pkin(2) - t443) * (pkin(2) + t443) + t351;
t369 = t194 * t195;
t241 = sqrt(-t369);
t360 = t229 * t241;
t166 = pkin(1) * t360;
t238 = pkin(2) ^ 2;
t212 = t238 + t351;
t472 = pkin(4) ^ 2;
t350 = pkin(3) ^ 2 - t472;
t198 = t212 + t350;
t433 = pkin(1) * t231;
t219 = -pkin(2) + t433;
t156 = -t198 * t219 - t166;
t186 = pkin(1) * t229 * t198;
t159 = -t219 * t241 + t186;
t205 = 0.1e1 / t212;
t237 = 0.1e1 / pkin(3);
t366 = t205 * t237;
t132 = qJ(2) + atan2(t159 * t366, t156 * t366);
t130 = sin(t132);
t131 = cos(t132);
t482 = rSges(4,1) * t130 + rSges(4,2) * t131;
t471 = 0.2e1 * pkin(1);
t270 = pkin(2) * (-t194 - t195) * t471;
t165 = t229 * t270;
t169 = 0.1e1 / t241;
t372 = t169 * t165;
t319 = -t372 / 0.2e1;
t357 = t231 * t241;
t456 = -0.2e1 * t219;
t124 = t186 + (-t357 + (pkin(2) * t456 + t319) * t229) * pkin(1);
t206 = 0.1e1 / t212 ^ 2;
t431 = pkin(2) * t229;
t346 = pkin(1) * t431;
t315 = t206 * t346;
t298 = -0.2e1 * t315;
t114 = (t124 * t205 + t156 * t298) * t237;
t317 = -t169 * t219 / 0.2e1;
t226 = t229 ^ 2;
t364 = t226 * t239;
t334 = pkin(2) * t364;
t355 = t198 * t433 + t166;
t128 = t165 * t317 + 0.2e1 * t334 + t355;
t116 = (t128 * t205 + t159 * t298) * t237;
t148 = 0.1e1 / t156 ^ 2;
t155 = t159 ^ 2;
t141 = t148 * t155 + 0.1e1;
t138 = 0.1e1 / t141;
t429 = pkin(3) * t212;
t339 = t148 * t429;
t306 = t159 * t339;
t292 = t138 * t306;
t147 = 0.1e1 / t156;
t340 = t147 * t429;
t308 = t138 * t340;
t66 = -t114 * t292 + t116 * t308 + 0.1e1;
t476 = t482 * t66;
t235 = 0.1e1 / t472;
t167 = pkin(2) * t360;
t199 = t212 - t350;
t218 = pkin(1) - t430;
t157 = t199 * t218 - t167;
t150 = t157 ^ 2;
t184 = t199 * t431;
t158 = t218 * t241 + t184;
t154 = t158 ^ 2;
t356 = t150 + t154;
t135 = t356 * t235 * t206;
t133 = t135 ^ (-0.1e1 / 0.2e1);
t234 = 0.1e1 / pkin(4);
t295 = t158 * t315;
t318 = t169 * t218 / 0.2e1;
t335 = t238 * t226 * pkin(1);
t354 = t199 * t430 + t167;
t129 = t165 * t318 + 0.2e1 * t335 + t354;
t382 = t129 * t205;
t261 = 0.2e1 * t295 - t382;
t296 = t157 * t315;
t125 = t184 + (-t357 + (t218 * t471 + t319) * t229) * pkin(2);
t383 = t125 * t205;
t262 = 0.2e1 * t296 - t383;
t207 = t205 * t206;
t418 = t133 / t135 * (0.2e1 * (t125 * t157 + t129 * t158) * t206 - 0.4e1 * t356 * t207 * t346) * t235;
t331 = t205 * t418;
t413 = t158 * rSges(5,2);
t417 = t157 * rSges(5,1);
t481 = t234 * ((rSges(5,1) * t262 + rSges(5,2) * t261) * t133 + (t417 / 0.2e1 + t413 / 0.2e1) * t331);
t232 = cos(qJ(1));
t230 = sin(qJ(1));
t359 = t230 * t231;
t363 = t229 * t230;
t182 = rSges(3,1) * t359 - rSges(3,2) * t363 - t232 * rSges(3,3);
t358 = t231 * t232;
t361 = t229 * t232;
t183 = rSges(3,1) * t358 - rSges(3,2) * t361 + rSges(3,3) * t230;
t208 = rSges(3,1) * t229 + rSges(3,2) * t231;
t196 = t208 * t230;
t197 = t208 * t232;
t26 = t230 * t481;
t27 = t232 * t481;
t33 = (t431 - t476) * t230;
t341 = pkin(2) * t361;
t425 = t476 * t232;
t34 = -t341 + t425;
t367 = t205 * t234;
t329 = t133 * t367;
t474 = t413 + t417;
t258 = t474 * t329;
t419 = rSges(5,3) * t232;
t79 = t419 + (-pkin(1) + t258) * t230;
t300 = t232 * t329;
t90 = rSges(5,3) * t230 - t474 * t300;
t80 = pkin(1) * t232 + t90;
t380 = t130 * t230;
t314 = -rSges(4,2) * t380 + rSges(4,3) * t232;
t422 = rSges(4,1) * t131;
t97 = (t422 - t430) * t230 + t314;
t376 = t131 * t232;
t313 = -rSges(4,1) * t376 + rSges(4,3) * t230;
t421 = rSges(4,2) * t130;
t98 = (t421 + t430) * t232 + t313;
t480 = m(3) * (-t182 * t196 - t183 * t197) + m(5) * (-t26 * t79 + t27 * t80) + m(4) * (t33 * t97 + t34 * t98);
t402 = Icges(4,4) * t130;
t278 = Icges(4,2) * t131 + t402;
t479 = t278 * t66;
t401 = Icges(4,4) * t131;
t281 = Icges(4,1) * t130 + t401;
t478 = t281 * t66;
t178 = Icges(3,4) * t359 - Icges(3,2) * t363 - Icges(3,6) * t232;
t224 = Icges(3,4) * t231;
t202 = -Icges(3,2) * t229 + t224;
t179 = Icges(3,6) * t230 + t202 * t232;
t403 = Icges(3,4) * t229;
t204 = Icges(3,1) * t231 - t403;
t181 = Icges(3,5) * t230 + t204 * t232;
t171 = t181 * t359;
t200 = Icges(3,5) * t231 - Icges(3,6) * t229;
t177 = Icges(3,3) * t230 + t200 * t232;
t310 = t177 * t232 - t171;
t176 = Icges(3,5) * t359 - Icges(3,6) * t363 - Icges(3,3) * t232;
t215 = Icges(3,4) * t363;
t180 = Icges(3,1) * t359 - Icges(3,5) * t232 - t215;
t353 = -t230 * t176 - t180 * t358;
t477 = -t178 * t361 - t179 * t363 - t310 - t353;
t228 = t232 ^ 2;
t462 = t230 ^ 2 + t228;
t113 = (-t383 / 0.2e1 + t296) * t234;
t115 = (t382 / 0.2e1 - t295) * t234;
t151 = 0.1e1 / t157;
t152 = 0.1e1 / t157 ^ 2;
t374 = t152 * t158;
t475 = t113 * t374 + t115 * t151;
t140 = t152 * t154 + 0.1e1;
t136 = 0.1e1 / t140;
t428 = pkin(4) * t212;
t470 = 0.2e1 * t136 * t428;
t65 = t475 * t470;
t445 = -t65 / 0.2e1;
t439 = -t230 / 0.2e1;
t438 = t230 / 0.2e1;
t436 = -t232 / 0.2e1;
t435 = t232 / 0.2e1;
t384 = t125 * t151 / t150;
t469 = 0.2e1 / t140 ^ 2 * (t129 * t374 - t154 * t384);
t375 = t148 * t159;
t385 = t124 * t147 * t148;
t468 = 0.2e1 / t141 ^ 2 * (t128 * t375 - t155 * t385);
t467 = Icges(5,4) * t261;
t465 = t157 * t234;
t464 = t158 * t234;
t203 = Icges(3,1) * t229 + t224;
t257 = t262 * Icges(5,4);
t400 = Icges(5,4) * t157;
t323 = t400 / 0.2e1;
t394 = Icges(5,2) * t158;
t399 = Icges(5,4) * t158;
t405 = Icges(5,1) * t157;
t461 = ((t405 / 0.2e1 + t399 / 0.2e1) * t331 + (Icges(5,1) * t262 + t467) * t133) * t464 - ((t323 + t394 / 0.2e1) * t331 + (Icges(5,2) * t261 + t257) * t133) * t465;
t378 = t131 * t281;
t381 = t130 * t278;
t460 = -t378 / 0.2e1 + t381 / 0.2e1;
t459 = t229 * t372 + t357;
t64 = t65 ^ 2;
t455 = -0.4e1 * t229;
t454 = -0.2e1 * t231;
t328 = t238 * t364;
t162 = t231 * t270 - 0.8e1 * t328;
t373 = t165 ^ 2 * t169 / t369;
t320 = -t373 / 0.4e1;
t253 = (t229 * t320 + (-t231 * t165 - t229 * t162 / 0.2e1) * t169) * t205;
t299 = t207 * t328;
t336 = t152 * t428;
t304 = t158 * t336;
t305 = t136 * t336;
t337 = t151 * t428;
t362 = t229 * t231;
t368 = t205 * t231;
t432 = pkin(2) * t206;
t10 = 0.2e1 * (t125 * t305 + t337 * t469) * t115 + 0.2e1 * (t158 * t384 * t470 - t129 * t305 + t304 * t469) * t113 + 0.2e1 * ((-((t218 * t373 / 0.4e1 + t162 * t318 - t184 + t459 * pkin(2)) * t205 / 0.2e1 + 0.4e1 * t158 * t299 + (0.3e1 * t238 * t205 * t362 + (-0.2e1 * t129 * t229 - t158 * t231) * t432) * pkin(1)) * t337 - (-(0.4e1 * t335 + t354) * t205 / 0.2e1 - 0.4e1 * t157 * t299 + (-t253 / 0.2e1 + (-t218 * t368 + (0.2e1 * t125 * t229 + t157 * t231) * t206) * pkin(1)) * pkin(2)) * t304) * t234 - 0.2e1 * t475 * pkin(4) * t346) * t136;
t452 = -t10 / 0.2e1;
t277 = t394 + t400;
t301 = t230 * t329;
t85 = Icges(5,6) * t232 + t277 * t301;
t280 = t399 + t405;
t87 = Icges(5,5) * t232 + t280 * t301;
t289 = t157 * t87 + t158 * t85;
t260 = t289 * t329;
t392 = Icges(5,6) * t158;
t398 = Icges(5,5) * t157;
t273 = t392 + t398;
t83 = Icges(5,3) * t232 + t273 * t301;
t47 = -t230 * t83 + t232 * t260;
t86 = Icges(5,6) * t230 - t277 * t300;
t88 = Icges(5,5) * t230 - t280 * t300;
t259 = (-t157 * t88 - t158 * t86) * t329;
t84 = Icges(5,3) * t230 - t273 * t300;
t82 = t230 * t84;
t48 = t232 * t259 + t82;
t451 = (t230 * t48 - t232 * t47) * t445;
t377 = t131 * t230;
t105 = Icges(4,1) * t377 - Icges(4,4) * t380 + Icges(4,5) * t232;
t101 = Icges(4,5) * t377 - Icges(4,6) * t380 + Icges(4,3) * t232;
t103 = Icges(4,4) * t377 - Icges(4,2) * t380 + Icges(4,6) * t232;
t379 = t130 * t232;
t424 = -t230 * t101 - t103 * t379;
t53 = t105 * t376 + t424;
t282 = -Icges(4,1) * t131 + t402;
t106 = Icges(4,5) * t230 + t232 * t282;
t275 = -Icges(4,5) * t131 + Icges(4,6) * t130;
t102 = Icges(4,3) * t230 + t232 * t275;
t279 = Icges(4,2) * t130 - t401;
t104 = Icges(4,6) * t230 + t232 * t279;
t407 = t230 * t102 + t104 * t379;
t54 = -t106 * t376 + t407;
t62 = t230 * t66;
t63 = t232 * t66;
t17 = -t53 * t63 + t54 * t62;
t450 = t17 / 0.2e1;
t448 = t62 / 0.2e1;
t447 = -t63 / 0.2e1;
t446 = t63 / 0.2e1;
t107 = rSges(4,1) * t377 + t314;
t108 = rSges(4,2) * t379 + t313;
t32 = -t107 * t62 + t108 * t63 + t462 * t430;
t442 = m(4) * t32;
t89 = t230 * t258 + t419;
t441 = m(5) * (-t230 * t89 + t232 * t90) * t65;
t276 = -Icges(3,5) * t229 - Icges(3,6) * t231;
t440 = t276 * t435;
t434 = m(3) * (t182 * t230 + t183 * t232);
t427 = t10 * t65;
t416 = t157 * rSges(5,2);
t395 = Icges(5,2) * t157;
t93 = (-t395 + t399) * t329;
t415 = t157 * t93;
t414 = t158 * rSges(5,1);
t404 = Icges(5,1) * t158;
t94 = (-t400 + t404) * t329;
t412 = t158 * t94;
t411 = t232 * t65;
t410 = t232 * t84;
t397 = Icges(5,5) * t158;
t393 = Icges(5,6) * t157;
t391 = t101 * t232;
t390 = t103 * t131;
t389 = t104 * t131;
t388 = t105 * t130;
t387 = t106 * t130;
t370 = t178 * t229;
t352 = t230 * t177 + t181 * t358;
t349 = qJD(2) * t230;
t348 = qJD(2) * t232;
t45 = t230 * t260 + t232 * t83;
t46 = t230 * t259 - t410;
t13 = (t230 * t46 - t232 * t45) * t65;
t302 = t158 * t329;
t303 = t157 * t329;
t256 = t302 * t86 + t303 * t88 + t83;
t333 = -t13 / 0.2e1 + ((t232 * t256 + t48 - t82) * t232 + (t230 * t256 + t410 + t47) * t230) * t445;
t332 = t451 + ((-t289 * t301 + t45 + t82) * t230 + (t46 + (-t302 * t85 - t303 * t87 + t84) * t232) * t232) * t65 / 0.2e1;
t326 = t411 / 0.2e1;
t145 = -t179 * t361 + t352;
t309 = t179 * t229 - t176;
t325 = (-t230 * (-t180 * t231 + t370) - t176 * t232) * t436 + (t232 * t309 + t145 - t352) * t435 + (t230 * t309 + t310 + t477) * t438;
t324 = t145 * t438 + t352 * t439 + (-t171 + (t177 + t370) * t232 + t353 + t477) * t436;
t201 = Icges(3,2) * t231 + t403;
t316 = t204 / 0.2e1 - t201 / 0.2e1;
t75 = t104 * t380;
t312 = -t102 * t232 + t75;
t311 = t106 * t131 + t101;
t307 = t138 * t339;
t297 = pkin(3) * t138 * t346;
t294 = 0.8e1 * t299;
t274 = Icges(4,5) * t130 + Icges(4,6) * t131;
t293 = (t230 * t63 - t232 * t62) * t274 * t66;
t11 = (((0.6e1 * t239 * pkin(2) * t362 + t162 * t317 + t219 * t320 - t186) * t205 + t159 * t294 + (t459 * t205 + (t128 * t455 + t159 * t454) * t432) * pkin(1)) * t308 - ((0.4e1 * t334 + t355) * t205 + t156 * t294 + (t253 + (t368 * t456 + (t124 * t455 + t156 * t454) * t206) * pkin(2)) * pkin(1)) * t292) * t237 + (-t124 * t307 + 0.2e1 * t147 * t297 - t340 * t468) * t116 + (0.2e1 * t159 * t138 * t385 * t429 - t128 * t307 - 0.2e1 * t297 * t375 + t306 * t468) * t114;
t290 = t11 * t482 - t430;
t288 = t157 * t85 - t158 * t87;
t287 = t157 * t86 - t158 * t88;
t286 = t157 * t94 + t158 * t93;
t285 = -t412 + t415;
t284 = 0.2e1 * t315;
t272 = t103 * t130 - t105 * t131;
t268 = t482 * t63 - t341;
t266 = t412 / 0.4e1 - t415 / 0.4e1;
t42 = t279 * t66;
t43 = t282 * t66;
t255 = (-t479 - t43) * t131 + (-t478 + t42) * t130;
t254 = t129 * t94 / 0.2e1 + ((-t404 / 0.2e1 + t323) * t331 + (-Icges(5,1) * t261 + t257) * t133) * t464 / 0.2e1 - t125 * t93 / 0.2e1 - ((-t399 / 0.2e1 + t395 / 0.2e1) * t331 + (Icges(5,2) * t262 - t467) * t133) * t465 / 0.2e1;
t37 = t230 * t479;
t38 = t232 * t479;
t39 = t230 * t478;
t40 = t232 * t478;
t252 = (t130 * t38 - t131 * t40) * t62 - (t130 * t37 - t131 * t39) * t63 + ((t387 + t389) * t62 - (-t388 - t390) * t63) * t66;
t251 = t234 * ((t398 / 0.2e1 + t392 / 0.2e1) * t331 + (Icges(5,5) * t262 + Icges(5,6) * t261) * t133);
t210 = rSges(3,1) * t231 - rSges(3,2) * t229;
t188 = t276 * t230;
t95 = (t414 - t416) * t329;
t92 = (-t393 + t397) * t329;
t55 = pkin(2) * t363 - t482 * t62;
t52 = -t106 * t377 + t312;
t51 = -t230 * t272 + t391;
t50 = t288 * t329;
t44 = (t421 - t422) * t66;
t41 = t275 * t66;
t31 = ((-t414 / 0.2e1 + t416 / 0.2e1) * t331 + (-rSges(5,1) * t261 + rSges(5,2) * t262) * t133) * t234;
t28 = ((-t397 / 0.2e1 + t393 / 0.2e1) * t331 + (-Icges(5,5) * t261 + Icges(5,6) * t262) * t133) * t234;
t21 = t232 * t251;
t20 = t230 * t251;
t16 = -t63 * t51 + t52 * t62;
t9 = t232 * t290 - t44 * t63;
t8 = t230 * t290 - t44 * t62;
t6 = (-t312 + t52 - t424) * t63 + (t103 * t380 - t391 + t407 + t51) * t62 + ((-t105 * t230 - t106 * t232) * t62 - (t105 * t232 - t106 * t230) * t63) * t131;
t5 = (t54 - t407) * t63 + (-t75 + t53) * t62 + (-t272 * t63 + t311 * t62) * t230 + ((t102 + t272) * t62 + t311 * t63) * t232;
t2 = (t203 / 0.2e1 + t202 / 0.2e1) * t231 + t316 * t229 + (t478 / 0.2e1 - t42 / 0.2e1) * t131 + (-t43 / 0.2e1 - t479 / 0.2e1) * t130 + (-t266 * t331 + (t205 * t254 + t285 * t315) * t133) * t234 + t480;
t1 = (t450 - t6 / 0.2e1) * t63 + (t5 / 0.2e1 + t16 / 0.2e1) * t62 + (-t332 * t65 + t324) * t232 + (-t333 * t65 + t325) * t230;
t3 = [t2 * qJD(2), t2 * qJD(1) + (t255 * t448 + t41 * t446 + t200 * t435 - (t28 * t435 - t332) * t65 + (-t180 / 0.2e1 + Icges(3,2) * t359 / 0.2e1 + t215 / 0.2e1) * t231 + (t178 / 0.2e1 - t203 * t439) * t229 + (-t388 / 0.2e1 - t390 / 0.2e1 - t274 * t435 + t460 * t230) * t11 + (t92 * t435 + t286 * t301 / 0.2e1 - t50 / 0.2e1) * t10 - t324) * t348 + (t41 * t448 + t255 * t447 + (t10 * t92 / 0.2e1 + t200 / 0.2e1) * t230 - (t28 * t438 - t333) * t65 + (t181 / 0.2e1 + t201 * t436) * t231 + (-t179 / 0.2e1 - t203 * t435) * t229 + (-t389 / 0.2e1 - t387 / 0.2e1 - t274 * t438 + (t378 - t381) * t435) * t11 - t325) * t349 + (((-t106 * t66 - t38) * t131 + (t104 * t66 - t40) * t130) * t448 + t6 * t446 - (t5 + t16) * t62 / 0.2e1 + ((t105 * t66 - t37) * t131 + (-t103 * t66 - t39) * t130 + t17) * t447 + (t33 * t268 - t34 * t55 + t8 * t98 + t9 * t97) * m(4) + ((t182 * t232 - t183 * t230) * t210 + (-t196 * t232 + t197 * t230) * t208) * m(3) + ((-t230 * t80 - t232 * t79) * t95 * t10 - ((t26 * t95 - t31 * t79) * t232 + (-t27 * t95 - t31 * t80) * t230) * t65) * m(5) + (-((-t87 * t232 / 0.4e1 - t88 * t230 / 0.4e1) * t158 + (t85 * t232 / 0.4e1 + t86 * t230 / 0.4e1) * t157) * t65 * t331 + ((-t288 * t284 + (t125 * t85 - t129 * t87 + t461 * t230) * t205) * t326 + (t287 * t284 * t445 + ((-t125 * t86 + t129 * t88 + t461 * t232) * t445 + t287 * t452 + t10 * t286 * t436) * t205) * t230) * t133) * t234) * qJD(2); t1 * qJD(2) + (-t50 * t65 * t438 + t131 * t42 / 0.2e1 + t130 * t43 / 0.2e1 + (t266 * t418 + (-t288 * t439 * t65 - t254) * t133) * t367 + (-pkin(1) * t133 * t234 * t285 * t432 - t316) * t229 + t460 * t66 - (t203 + t202) * t231 / 0.2e1 - t480) * qJD(1), t1 * qJD(1) + (m(4) * (t268 * t9 - t55 * t8 + t32 * (t425 * t63 - t462 * t431)) + m(3) * t462 * t208 * t210 + m(5) * t462 * (t31 * t64 - t427 * t95) * t95) * qJD(2) + (-t197 * t434 + t252 * t448 + t293 * t447 - (t10 * t90 - t27 * t65) * t441 - t13 * t452 + (-t45 * t427 / 0.2e1 + (-t188 / 0.2e1 - t20 * t64 / 0.2e1) * t232) * t232 + (-t16 / 0.2e1 - t53 * t448 - t51 * t447 + t108 * t442) * t11) * t348 + (-t196 * t434 + t228 * t440 + (t11 * t54 - t293) * t448 + (t52 * t11 + t252) * t447 - (-t10 * t89 - t26 * t65) * t441 + (-t11 * t107 + t62 * t476) * t442 + t10 * t451 + t11 * t450 + (t188 * t436 + (t10 * t48 + t20 * t411) * t445 + (t440 + t21 * t64 / 0.2e1) * t230) * t230 + (t21 * t411 + (t46 + t47) * t10) * t326) * t349;];
Cq = t3;
