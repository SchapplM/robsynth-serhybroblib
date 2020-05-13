% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnDE1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:40
% EndTime: 2020-04-12 19:26:28
% DurationCPUTime: 23.57s
% Computational Cost: add. (404562->455), mult. (577288->872), div. (27740->20), fcn. (164811->10), ass. (0->359)
t259 = pkin(2) ^ 2;
t260 = pkin(1) ^ 2;
t251 = cos(qJ(2));
t450 = pkin(2) * t251;
t375 = -0.2e1 * pkin(1) * t450 + t260;
t231 = t259 + t375;
t225 = 0.1e1 / t231 ^ 2;
t493 = pkin(3) ^ 2;
t258 = 0.1e1 / t493;
t249 = sin(qJ(2));
t463 = -pkin(3) - pkin(4);
t213 = (pkin(2) - t463) * (pkin(2) + t463) + t375;
t462 = pkin(4) - pkin(3);
t214 = (pkin(2) - t462) * (pkin(2) + t462) + t375;
t399 = t213 * t214;
t263 = sqrt(-t399);
t390 = t249 * t263;
t185 = pkin(1) * t390;
t494 = pkin(4) ^ 2;
t373 = t493 - t494;
t217 = t231 + t373;
t453 = pkin(1) * t251;
t238 = -pkin(2) + t453;
t175 = -t217 * t238 - t185;
t165 = t175 ^ 2;
t205 = pkin(1) * t249 * t217;
t178 = -t238 * t263 + t205;
t174 = t178 ^ 2;
t381 = t165 + t174;
t153 = t381 * t258 * t225;
t150 = t153 ^ (-0.1e1 / 0.2e1);
t224 = 0.1e1 / t231;
t257 = 0.1e1 / pkin(3);
t246 = t249 ^ 2;
t492 = 0.2e1 * pkin(1);
t356 = pkin(2) * t492;
t327 = t225 * t356;
t387 = t251 * t178;
t279 = (t175 * t246 + t249 * t387) * t327;
t294 = (-t213 - t214) * t356;
t184 = t249 * t294;
t188 = 0.1e1 / t263;
t402 = t188 * t184;
t342 = -t402 / 0.2e1;
t385 = t251 * t263;
t474 = -0.2e1 * t238;
t144 = t205 + (-t385 + (pkin(2) * t474 + t342) * t249) * pkin(1);
t340 = -t188 * t238 / 0.2e1;
t396 = t246 * t260;
t360 = pkin(2) * t396;
t379 = t217 * t453 + t185;
t146 = t184 * t340 + 0.2e1 * t360 + t379;
t226 = t224 * t225;
t451 = pkin(2) * t249;
t372 = pkin(1) * t451;
t310 = 0.4e1 * t226 * t372;
t371 = 0.2e1 * t225;
t485 = t150 * t224;
t354 = ((t144 * t175 + t146 * t178) * t371 - t381 * t310) * t258 / t153 * t485;
t395 = t249 * t175;
t280 = (-t387 / 0.2e1 - t395 / 0.2e1) * t354;
t382 = t146 + t175;
t334 = t382 * t251;
t413 = t144 * t249;
t505 = t257 * (t150 * (t224 * (t334 + t413) - t279) + t280);
t350 = t257 * t485;
t499 = t387 + t395;
t125 = t499 * t350;
t394 = t249 * t178;
t141 = t350 * t394;
t406 = t175 * t251;
t309 = t350 * t406;
t126 = t141 - t309;
t99 = -Icges(4,1) * t125 + Icges(4,4) * t126;
t464 = t99 / 0.2e1;
t460 = -t125 / 0.2e1;
t459 = t126 / 0.2e1;
t252 = cos(qJ(1));
t454 = t252 / 0.2e1;
t250 = sin(qJ(1));
t389 = t250 * t251;
t136 = t250 * t141;
t121 = t250 * t309 - t136;
t122 = t250 * t125;
t95 = -t121 * rSges(4,1) + t122 * rSges(4,2) - t252 * rSges(4,3);
t82 = -pkin(2) * t389 - t95;
t255 = 0.1e1 / t494;
t186 = pkin(2) * t390;
t218 = t231 - t373;
t237 = pkin(1) - t450;
t176 = t218 * t237 - t186;
t169 = t176 ^ 2;
t203 = t218 * t451;
t177 = t237 * t263 + t203;
t173 = t177 ^ 2;
t380 = t169 + t173;
t152 = t380 * t255 * t225;
t148 = t152 ^ (-0.1e1 / 0.2e1);
t254 = 0.1e1 / pkin(4);
t336 = t225 * t372;
t316 = t177 * t336;
t341 = t188 * t237 / 0.2e1;
t361 = t259 * t246 * pkin(1);
t378 = t218 * t450 + t186;
t147 = t184 * t341 + 0.2e1 * t361 + t378;
t409 = t147 * t224;
t286 = 0.2e1 * t316 - t409;
t317 = t176 * t336;
t145 = t203 + (-t385 + (t237 * t492 + t342) * t249) * pkin(2);
t411 = t145 * t224;
t287 = 0.2e1 * t317 - t411;
t418 = ((t145 * t176 + t147 * t177) * t371 - t380 * t310) * t255 * t148 / t152;
t355 = t224 * t418;
t438 = t177 * rSges(5,2);
t441 = t176 * rSges(5,1);
t503 = t254 * (t148 * (rSges(5,1) * t287 + rSges(5,2) * t286) + (t441 / 0.2e1 + t438 / 0.2e1) * t355);
t393 = t249 * t250;
t197 = Icges(3,4) * t389 - Icges(3,2) * t393 - Icges(3,6) * t252;
t243 = Icges(3,4) * t251;
t221 = -Icges(3,2) * t249 + t243;
t198 = Icges(3,6) * t250 + t221 * t252;
t428 = Icges(3,4) * t249;
t223 = Icges(3,1) * t251 - t428;
t200 = Icges(3,5) * t250 + t223 * t252;
t190 = t200 * t389;
t219 = Icges(3,5) * t251 - Icges(3,6) * t249;
t196 = Icges(3,3) * t250 + t219 * t252;
t333 = t196 * t252 - t190;
t195 = Icges(3,5) * t389 - Icges(3,6) * t393 - Icges(3,3) * t252;
t234 = Icges(3,4) * t393;
t199 = Icges(3,1) * t389 - Icges(3,5) * t252 - t234;
t386 = t251 * t252;
t377 = -t250 * t195 - t199 * t386;
t391 = t249 * t252;
t501 = -t197 * t391 - t198 * t393 - t333 - t377;
t132 = (-t411 / 0.2e1 + t317) * t254;
t134 = (t409 / 0.2e1 - t316) * t254;
t170 = 0.1e1 / t176;
t171 = 0.1e1 / t176 ^ 2;
t407 = t171 * t177;
t498 = t132 * t407 + t134 * t170;
t497 = t438 + t441;
t322 = t252 * t350;
t123 = (t394 - t406) * t322;
t124 = t499 * t322;
t86 = Icges(4,5) * t121 - Icges(4,6) * t122 + Icges(4,3) * t252;
t89 = Icges(4,4) * t121 - Icges(4,2) * t122 + Icges(4,6) * t252;
t92 = Icges(4,1) * t121 - Icges(4,4) * t122 + Icges(4,5) * t252;
t33 = -t123 * t92 - t124 * t89 - t250 * t86;
t88 = Icges(4,5) * t123 + Icges(4,6) * t124 + Icges(4,3) * t250;
t91 = Icges(4,4) * t123 + Icges(4,2) * t124 + Icges(4,6) * t250;
t94 = Icges(4,1) * t123 + Icges(4,4) * t124 + Icges(4,5) * t250;
t34 = t123 * t94 + t124 * t91 + t250 * t88;
t167 = 0.1e1 / t175 ^ 2;
t159 = t167 * t174 + 0.1e1;
t156 = 0.1e1 / t159;
t449 = pkin(3) * t231;
t367 = t156 * t449;
t320 = -0.2e1 * t336;
t133 = (t144 * t224 + t175 * t320) * t257;
t135 = (t146 * t224 + t178 * t320) * t257;
t166 = 0.1e1 / t175;
t408 = t167 * t178;
t478 = -t133 * t408 + t135 * t166;
t461 = t367 * t478 + 0.1e1;
t76 = t461 * t250;
t77 = t461 * t252;
t496 = -t33 * t77 + t34 * t76;
t31 = t121 * t92 - t122 * t89 + t252 * t86;
t158 = t171 * t173 + 0.1e1;
t154 = 0.1e1 / t158;
t448 = pkin(4) * t231;
t491 = 0.2e1 * t154 * t448;
t79 = t498 * t491;
t466 = -t79 / 0.2e1;
t490 = -t250 / 0.2e1;
t457 = t250 / 0.2e1;
t455 = -t252 / 0.2e1;
t439 = t177 * rSges(5,1);
t440 = t176 * rSges(5,2);
t489 = ((-t439 / 0.2e1 + t440 / 0.2e1) * t355 + (-rSges(5,1) * t286 + rSges(5,2) * t287) * t148) * t254 * t79;
t488 = Icges(5,4) * t286;
t412 = t145 * t170 / t169;
t487 = 0.2e1 * (t147 * t407 - t173 * t412) / t158 ^ 2;
t414 = t144 * t166 / t165;
t486 = 0.2e1 * (t146 * t408 - t174 * t414) / t159 ^ 2;
t483 = t176 * t254;
t482 = t177 * t254;
t222 = Icges(3,1) * t249 + t243;
t282 = t287 * Icges(5,4);
t427 = Icges(5,4) * t176;
t344 = t427 / 0.2e1;
t421 = Icges(5,2) * t177;
t426 = Icges(5,4) * t177;
t430 = Icges(5,1) * t176;
t480 = ((t430 / 0.2e1 + t426 / 0.2e1) * t355 + (Icges(5,1) * t287 + t488) * t148) * t482 - ((t344 + t421 / 0.2e1) * t355 + (Icges(5,2) * t286 + t282) * t148) * t483;
t479 = t249 * t402 + t385;
t419 = Icges(5,6) * t177;
t425 = Icges(5,5) * t176;
t304 = t419 + t425;
t397 = t224 * t254;
t351 = t148 * t397;
t324 = t250 * t351;
t111 = Icges(5,3) * t252 + t304 * t324;
t306 = t421 + t427;
t113 = Icges(5,6) * t252 + t306 * t324;
t307 = t426 + t430;
t115 = Icges(5,5) * t252 + t307 * t324;
t302 = t177 * t113 + t176 * t115;
t284 = t302 * t351;
t66 = -t250 * t111 + t252 * t284;
t323 = t252 * t351;
t112 = Icges(5,3) * t250 - t304 * t323;
t110 = t250 * t112;
t114 = Icges(5,6) * t250 - t306 * t323;
t116 = Icges(5,5) * t250 - t307 * t323;
t283 = (-t177 * t114 - t176 * t116) * t351;
t67 = t252 * t283 + t110;
t477 = t250 * t67 - t252 * t66;
t64 = t111 * t252 + t250 * t284;
t417 = t112 * t252;
t65 = t250 * t283 - t417;
t476 = t250 * t65 - t252 * t64;
t32 = -t121 * t94 + t122 * t91 - t252 * t88;
t247 = t250 ^ 2;
t248 = t252 ^ 2;
t473 = m(3) / 0.4e1;
t472 = m(4) / 0.4e1;
t471 = m(5) / 0.4e1;
t469 = t76 / 0.2e1;
t468 = -t77 / 0.2e1;
t467 = t77 / 0.2e1;
t98 = -Icges(4,4) * t125 + Icges(4,2) * t126;
t465 = t98 / 0.2e1;
t458 = -t224 / 0.2e1;
t452 = pkin(2) * t225;
t392 = t249 * t251;
t278 = (t175 * t392 - t178 * t246) * t327;
t281 = (t406 / 0.2e1 - t394 / 0.2e1) * t354;
t384 = t252 * t257;
t410 = t146 * t249;
t41 = (t281 + ((-t144 * t251 + t410) * t224 + t278) * t150) * t384 + t124;
t383 = t144 - t178;
t42 = (t280 + (-t279 + (t249 * t383 + t334) * t224) * t150) * t384;
t445 = t41 * rSges(4,1) + t42 * rSges(4,2);
t443 = rSges(5,3) * t252;
t349 = t259 * t396;
t181 = t251 * t294 - 0.8e1 * t349;
t403 = t184 ^ 2 * t188 / t399;
t343 = -t403 / 0.4e1;
t275 = (t249 * t343 + (-t251 * t184 - t249 * t181 / 0.2e1) * t188) * t224;
t321 = t226 * t349;
t362 = t171 * t448;
t328 = t177 * t362;
t329 = t154 * t362;
t363 = t170 * t448;
t398 = t224 * t251;
t13 = 0.2e1 * (t145 * t329 + t363 * t487) * t134 + 0.2e1 * (t177 * t412 * t491 - t147 * t329 + t328 * t487) * t132 + 0.2e1 * ((-((t237 * t403 / 0.4e1 + t181 * t341 - t203 + t479 * pkin(2)) * t224 / 0.2e1 + 0.4e1 * t177 * t321 + (0.3e1 * t259 * t224 * t392 + (-0.2e1 * t147 * t249 - t177 * t251) * t452) * pkin(1)) * t363 - ((0.4e1 * t361 + t378) * t458 - 0.4e1 * t176 * t321 + (-t275 / 0.2e1 + (-t237 * t398 + (0.2e1 * t145 * t249 + t176 * t251) * t225) * pkin(1)) * pkin(2)) * t328) * t254 - 0.2e1 * t498 * pkin(4) * t372) * t154;
t442 = t13 * t252;
t435 = t250 * t79;
t429 = Icges(5,1) * t177;
t424 = Icges(5,5) * t177;
t422 = Icges(5,2) * t176;
t420 = Icges(5,6) * t176;
t129 = (-t422 + t426) * t351;
t405 = t176 * t129;
t130 = (-t427 + t429) * t351;
t404 = t177 * t130;
t400 = t197 * t249;
t388 = t250 * t252;
t376 = t250 * t196 + t200 * t386;
t374 = t247 + t248;
t370 = pkin(2) * t393;
t369 = pkin(2) * t391;
t366 = t166 * t449;
t365 = t167 * t449;
t325 = t177 * t351;
t326 = t176 * t351;
t277 = t114 * t325 + t116 * t326 + t111;
t28 = t476 * t79;
t359 = -t28 / 0.2e1 + ((t252 * t277 - t110 + t67) * t252 + (t250 * t277 + t417 + t66) * t250) * t466;
t358 = t79 * t388;
t29 = t477 * t79;
t357 = ((-t302 * t324 + t110 + t64) * t250 + (t65 + (-t113 * t325 - t115 * t326 + t112) * t252) * t252) * t79 / 0.2e1 - t29 / 0.2e1;
t347 = -t435 / 0.2e1;
t346 = t79 * t455;
t163 = -t198 * t391 + t376;
t332 = t198 * t249 - t195;
t339 = (t252 * t332 + t163 - t376) * t454 + (-t250 * (-t199 * t251 + t400) - t195 * t252) * t455 + (t250 * t332 + t333 + t501) * t457;
t338 = t163 * t457 + t376 * t490 + (-t190 + (t196 + t400) * t252 + t377 + t501) * t455;
t220 = Icges(3,2) * t251 + t428;
t337 = t223 / 0.2e1 - t220 / 0.2e1;
t331 = t156 * t365;
t330 = t178 * t365;
t319 = 0.4e1 * t374;
t120 = rSges(5,3) * t250 - t323 * t497;
t315 = 0.8e1 * t321;
t44 = (t281 + (t278 + (t249 * t382 - t251 * t383) * t224) * t150) * t257;
t39 = t250 * t44;
t40 = t250 * t505 - t136;
t314 = rSges(4,1) * t39 + rSges(4,2) * t40;
t100 = -rSges(4,1) * t125 + rSges(4,2) * t126;
t14 = (-t144 * t331 - t366 * t486) * t135 + (0.2e1 * t178 * t367 * t414 - t146 * t331 + t330 * t486) * t133 + (0.2e1 * t478 * pkin(3) * t372 + (((0.6e1 * pkin(2) * t260 * t392 + t181 * t340 + t238 * t343 - t205) * t224 + t178 * t315 + (t479 * t224 + (-0.2e1 * t387 - 0.4e1 * t410) * t452) * pkin(1)) * t366 - ((0.4e1 * t360 + t379) * t224 + t175 * t315 + (t275 + (t398 * t474 + (-0.2e1 * t406 - 0.4e1 * t413) * t225) * pkin(2)) * pkin(1)) * t330) * t257) * t156;
t312 = -t100 * t14 - t450;
t227 = rSges(3,1) * t249 + rSges(3,2) * t251;
t311 = t249 * t327;
t305 = -Icges(3,5) * t249 - Icges(3,6) * t251;
t285 = t497 * t351;
t107 = t443 + (-pkin(1) + t285) * t250;
t108 = pkin(1) * t252 + t120;
t303 = -t107 * t252 - t108 * t250;
t301 = t176 * t113 - t177 * t115;
t300 = t114 * t176 - t116 * t177;
t299 = -(t250 * t285 + t443) * t250 + t120 * t252;
t298 = t129 * t177 + t130 * t176;
t297 = -t404 + t405;
t293 = -t100 * t77 - t369;
t290 = t404 / 0.4e1 - t405 / 0.4e1;
t96 = rSges(4,1) * t123 + rSges(4,2) * t124 + rSges(4,3) * t250;
t276 = t147 * t130 / 0.2e1 + ((-t429 / 0.2e1 + t344) * t355 + (-Icges(5,1) * t286 + t282) * t148) * t482 / 0.2e1 - t145 * t129 / 0.2e1 - ((-t426 / 0.2e1 + t422 / 0.2e1) * t355 + (Icges(5,2) * t287 - t488) * t148) * t483 / 0.2e1;
t201 = rSges(3,1) * t389 - rSges(3,2) * t393 - t252 * rSges(3,3);
t202 = rSges(3,1) * t386 - rSges(3,2) * t391 + rSges(3,3) * t250;
t215 = t227 * t250;
t216 = t227 * t252;
t43 = t141 - t505;
t25 = Icges(4,4) * t43 + Icges(4,2) * t44;
t26 = Icges(4,1) * t43 + Icges(4,4) * t44;
t53 = t250 * t503;
t54 = t252 * t503;
t274 = 0.4e1 * (-t201 * t215 - t202 * t216) * t473 + 0.4e1 * (-t107 * t53 + t108 * t54) * t471 + t25 * t459 + t26 * t460 + t43 * t464 + t44 * t465;
t273 = t254 * ((t425 / 0.2e1 + t419 / 0.2e1) * t355 + (Icges(5,5) * t287 + Icges(5,6) * t286) * t148);
t229 = rSges(3,1) * t251 - rSges(3,2) * t249;
t208 = t305 * t252;
t207 = t305 * t250;
t131 = (t439 - t440) * t351;
t128 = (-t420 + t424) * t351;
t97 = -Icges(4,5) * t125 + Icges(4,6) * t126;
t83 = pkin(2) * t386 + t96;
t75 = t301 * t351;
t61 = t100 * t76 + t370;
t55 = ((-t424 / 0.2e1 + t420 / 0.2e1) * t355 + (-Icges(5,5) * t286 + Icges(5,6) * t287) * t148) * t254;
t48 = t252 * t273;
t47 = t250 * t273;
t27 = rSges(4,1) * t43 + rSges(4,2) * t44;
t24 = Icges(4,5) * t43 + Icges(4,6) * t44;
t22 = Icges(4,1) * t41 + Icges(4,4) * t42;
t21 = Icges(4,1) * t39 + Icges(4,4) * t40;
t20 = Icges(4,4) * t41 + Icges(4,2) * t42;
t19 = Icges(4,4) * t39 + Icges(4,2) * t40;
t18 = Icges(4,5) * t41 + Icges(4,6) * t42;
t17 = Icges(4,5) * t39 + Icges(4,6) * t40;
t16 = -t369 + t445;
t15 = -t314 + t370;
t8 = 0.4e1 * t15 * t82 + 0.4e1 * t16 * t83;
t7 = t252 * t312 - t27 * t77;
t6 = t250 * t312 - t27 * t76;
t2 = t8 * t472 + (t222 / 0.2e1 + t221 / 0.2e1) * t251 + t337 * t249 + (-t290 * t355 + (t224 * t276 + t297 * t336) * t148) * t254 + t274;
t1 = (-t357 * t79 + t338) * t252 + (-t359 * t79 + t339) * t250;
t3 = [t2 * qJD(2), t2 * qJD(1) + (t496 * t467 + (t24 * t467 + t219 * t454 - (t454 * t55 - t357) * t79 + (-t199 / 0.2e1 + Icges(3,2) * t389 / 0.2e1 + t234 / 0.2e1) * t251 + (t197 / 0.2e1 - t222 * t490) * t249 + (t92 * t460 + t89 * t459 + t121 * t464 - t122 * t98 / 0.2e1 + t97 * t454) * t14 + (-t75 / 0.2e1 + t128 * t454 + t298 * t324 / 0.2e1) * t13 - t338) * t252 + ((t13 * t128 / 0.2e1 + t219 / 0.2e1) * t250 - (t457 * t55 - t359) * t79 + (t200 / 0.2e1 + t220 * t455) * t251 + (-t198 / 0.2e1 - t222 * t454) * t249 + (t123 * t464 + t124 * t465 + t457 * t97 + t459 * t91 + t460 * t94) * t14 - t339) * t250 + (t15 * t293 - t16 * t61 + t6 * t83 + t7 * t82) * m(4) + ((t201 * t252 - t202 * t250) * t229 + (-t215 * t252 + t216 * t250) * t227) * m(3) + (-t303 * t489 + (-(-t250 * t54 + t252 * t53) * t79 + t303 * t13) * t131) * m(5) + (-((-t115 * t252 / 0.4e1 - t116 * t250 / 0.4e1) * t177 + (t113 * t252 / 0.4e1 + t114 * t250 / 0.4e1) * t176) * t79 * t355 + (-(-t301 * t311 + (t113 * t145 - t115 * t147 + t250 * t480) * t224) * t346 + (t300 * t311 * t466 + ((-t114 * t145 + t116 * t147 + t252 * t480) * t466 - t300 * t13 / 0.2e1) * t224 + t298 * t442 * t458) * t250) * t148) * t254 + (t123 * t26 + t124 * t25 - t125 * t22 + t126 * t20 + t24 * t250 + t41 * t99 + t42 * t98 + t43 * t94 + t44 * t91) * t469 + (-t121 * t26 + t122 * t25 - t125 * t21 + t126 * t19 + t39 * t99 + t40 * t98 - t43 * t92 - t44 * t89 + t496) * t468) * qJD(2); t1 * qJD(2) + (t75 * t347 - t274 - m(4) * t8 / 0.4e1 + (t290 * t418 + (t301 * t435 / 0.2e1 - t276) * t148) * t397 + (-pkin(1) * t148 * t254 * t297 * t452 - t337) * t249 - (t222 + t221) * t251 / 0.2e1) * qJD(1), t1 * qJD(1) + (0.4e1 * (t293 * t7 - t61 * t6 + (t374 * t450 + t76 * t95 + t77 * t96) * (t76 * t314 + t77 * t445 + (t250 * t95 + t252 * t96) * t14 - t374 * t451)) * t472 + ((t123 * t22 + t124 * t20 + t18 * t250 + t41 * t94 + t42 * t91) * t76 - (t123 * t21 + t124 * t19 + t17 * t250 - t41 * t92 - t42 * t89) * t77 + (t250 * t34 - t252 * t33) * t14) * t469 + (-0.4e1 * t299 * t79 * (-(t250 * t53 + t252 * t54) * t79 + t299 * t13) - (t13 * t131 - t489) * t79 * t319 * t131) * t471 + ((-t121 * t22 + t122 * t20 - t18 * t252 + t39 * t94 + t40 * t91) * t76 - (-t121 * t21 + t122 * t19 - t17 * t252 - t39 * t92 - t40 * t89) * t77 + (t250 * t32 - t252 * t31) * t14) * t468 + (0.4e1 * (t201 * t250 + t202 * t252) * (-t215 * t250 - t216 * t252) + t229 * t227 * t319) * t473 + t28 * t442 / 0.2e1 - (-t248 * t47 * t79 + t13 * t476 + t358 * t48) * t346 + (-t247 * t48 * t79 + t13 * t477 + t358 * t47) * t347 + (-t13 * t29 + t14 * t496 - t207 * t388 + t208 * t247) * t457 + (t14 * (-t31 * t77 + t32 * t76) + t207 * t248 - t208 * t388) * t455) * qJD(2);];
Cq = t3;
