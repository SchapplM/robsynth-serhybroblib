% Calculate vector of inverse dynamics joint torques for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m1DE_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:29
% DurationCPUTime: 16.00s
% Computational Cost: add. (6755->614), mult. (11203->869), div. (0->0), fcn. (8712->10), ass. (0->358)
t253 = sin(qJ(1));
t252 = sin(qJ(2));
t255 = cos(qJ(3));
t251 = sin(qJ(3));
t256 = cos(qJ(2));
t410 = t251 * t256;
t184 = t252 * t255 + t410;
t441 = pkin(3) * qJD(3);
t372 = t184 * t441;
t156 = t253 * t372;
t231 = pkin(3) * t255 + pkin(2);
t158 = pkin(3) * t410 + t231 * t252;
t384 = qJD(2) * t253;
t397 = t158 * t384 + t156;
t411 = t251 * t252;
t157 = -pkin(3) * t411 + t231 * t256;
t318 = pkin(1) + t157;
t151 = rSges(5,1) + t318;
t257 = cos(qJ(1));
t104 = -rSges(5,3) * t253 + t151 * t257;
t496 = qJD(1) * t104;
t54 = -t397 + t496;
t248 = qJD(1) + qJD(4);
t209 = rSges(6,1) * t257 - rSges(6,2) * t253;
t211 = rSges(6,1) * t253 + rSges(6,2) * t257;
t250 = sin(qJ(4));
t254 = cos(qJ(4));
t337 = -t209 * t254 + t211 * t250;
t497 = t248 * t337;
t247 = qJD(2) + qJD(3);
t249 = qJ(2) + qJ(3);
t239 = sin(t249);
t240 = cos(t249);
t326 = rSges(4,1) * t239 + rSges(4,2) * t240;
t480 = t326 * t247;
t105 = -rSges(5,3) * t257 - t151 * t253;
t495 = qJD(1) * t105;
t465 = m(5) + m(6);
t494 = m(4) + t465;
t383 = qJD(2) * t256;
t490 = pkin(3) * (qJD(3) * t184 + t251 * t383);
t328 = rSges(3,1) * t252 + rSges(3,2) * t256;
t406 = t256 * t257;
t445 = rSges(3,3) * t253;
t352 = rSges(3,1) * t406 - t445;
t408 = t252 * t257;
t144 = -rSges(3,2) * t408 + t352;
t356 = pkin(1) * t257 + t144;
t84 = qJD(1) * t356 - t328 * t384;
t438 = t84 * t257;
t409 = t252 * t253;
t353 = rSges(3,2) * t409 - rSges(3,3) * t257;
t407 = t253 * t256;
t143 = rSges(3,1) * t407 - t353;
t357 = -pkin(1) * t253 - t143;
t382 = qJD(2) * t257;
t83 = t357 * qJD(1) - t328 * t382;
t319 = -t253 * t83 + t438;
t489 = t319 * t328;
t414 = t240 * t253;
t417 = t239 * t253;
t428 = Icges(4,6) * t257;
t117 = Icges(4,4) * t414 - Icges(4,2) * t417 + t428;
t432 = Icges(4,4) * t240;
t175 = -Icges(4,1) * t239 - t432;
t488 = -t175 * t253 + t117;
t314 = -Icges(4,2) * t239 + t432;
t118 = -Icges(4,6) * t253 + t257 * t314;
t487 = -t175 * t257 + t118;
t433 = Icges(4,4) * t239;
t316 = Icges(4,1) * t240 - t433;
t120 = -Icges(4,5) * t253 + t257 * t316;
t173 = -Icges(4,2) * t240 - t433;
t486 = -t173 * t257 - t120;
t485 = t173 + t316;
t484 = t314 - t175;
t205 = -rSges(6,1) * t250 - rSges(6,2) * t254;
t212 = rSges(6,1) * t254 - rSges(6,2) * t250;
t304 = t205 * t257 - t212 * t253;
t448 = rSges(3,1) * t256;
t366 = -pkin(1) - t448;
t439 = t252 * rSges(3,2);
t452 = pkin(2) * t256;
t265 = qJD(2) ^ 2;
t453 = pkin(2) * t252;
t331 = -qJDD(2) * t452 + t265 * t453;
t377 = qJDD(2) + qJDD(3);
t95 = (t239 * t247 ^ 2 - t240 * t377) * pkin(3) + t331;
t482 = t95 * (-pkin(3) * t240 - t452);
t185 = t255 * t256 - t411;
t371 = t185 * t441;
t300 = -qJD(2) * t157 - t371;
t387 = qJD(1) * t253;
t478 = t158 * t387 + t257 * t300;
t477 = qJD(2) * t185;
t199 = -Icges(3,5) * t252 - Icges(3,6) * t256;
t289 = qJD(2) * t199;
t435 = Icges(3,4) * t252;
t317 = Icges(3,1) * t256 - t435;
t140 = -Icges(3,5) * t253 + t257 * t317;
t434 = Icges(3,4) * t256;
t315 = -Icges(3,2) * t252 + t434;
t138 = -Icges(3,6) * t253 + t257 * t315;
t424 = t138 * t252;
t308 = -t140 * t256 + t424;
t313 = Icges(3,5) * t256 - Icges(3,6) * t252;
t427 = Icges(3,3) * t257;
t476 = -t257 * t289 + (t253 * t313 + t308 + t427) * qJD(1);
t429 = Icges(3,6) * t257;
t137 = Icges(3,4) * t407 - Icges(3,2) * t409 + t429;
t226 = Icges(3,4) * t409;
t431 = Icges(3,5) * t257;
t139 = Icges(3,1) * t407 - t226 + t431;
t309 = t137 * t252 - t139 * t256;
t136 = -Icges(3,3) * t253 + t257 * t313;
t389 = qJD(1) * t136;
t475 = qJD(1) * t309 - t253 * t289 - t389;
t171 = -Icges(4,5) * t239 - Icges(4,6) * t240;
t146 = t171 * t257;
t312 = Icges(4,5) * t240 - Icges(4,6) * t239;
t425 = t118 * t239;
t426 = Icges(4,3) * t257;
t474 = -t247 * t146 + (-t120 * t240 + t253 * t312 + t425 + t426) * qJD(1);
t145 = t171 * t253;
t215 = Icges(4,4) * t417;
t430 = Icges(4,5) * t257;
t119 = Icges(4,1) * t414 - t215 + t430;
t311 = t117 * t239 - t119 * t240;
t116 = -Icges(4,3) * t253 + t257 * t312;
t390 = qJD(1) * t116;
t473 = qJD(1) * t311 - t145 * t247 - t390;
t334 = m(4) * rSges(4,1) + pkin(3) * t465;
t464 = m(4) * rSges(4,2);
t142 = rSges(3,2) * m(3) + t251 * t334 + t255 * t464;
t107 = rSges(3,1) * m(3) + pkin(2) * t494 - t251 * t464 + t255 * t334;
t307 = t173 * t239 - t175 * t240;
t472 = qJD(1) * t307 + t247 * t312;
t201 = -Icges(3,2) * t256 - t435;
t203 = -Icges(3,1) * t252 - t434;
t305 = t201 * t252 - t203 * t256;
t471 = t305 * qJD(1) + qJD(2) * t313;
t229 = rSges(6,1) * g(1) + rSges(6,2) * g(2);
t230 = rSges(6,1) * g(2) - rSges(6,2) * g(1);
t469 = t253 * (t229 * t254 + t230 * t250) + t257 * (t229 * t250 - t230 * t254);
t468 = -t137 * t257 + t138 * t253;
t198 = t247 * t257;
t412 = t247 * t253;
t467 = qJD(1) * t485 + t198 * t488 - t412 * t487;
t217 = rSges(4,2) * t417;
t374 = rSges(4,1) * t414;
t121 = rSges(4,3) * t257 - t217 + t374;
t413 = t240 * t257;
t218 = rSges(4,1) * t413;
t416 = t239 * t257;
t444 = rSges(4,3) * t253;
t122 = -rSges(4,2) * t416 + t218 - t444;
t153 = t326 * t253;
t154 = t326 * t257;
t373 = pkin(2) * t383;
t196 = qJD(1) * t217;
t75 = -qJD(1) * t374 + t196 + (-rSges(4,3) * qJD(1) - t480) * t257;
t446 = rSges(4,2) * t239;
t447 = rSges(4,1) * t240;
t327 = -t446 + t447;
t76 = -t247 * t153 + (t257 * t327 - t444) * qJD(1);
t466 = (-t412 * t153 - t154 * t198 - t253 * t76 + t122 * t387 + (-qJD(1) * t121 - t75) * t257) * (-t121 * t412 - t122 * t198 - t373);
t266 = qJD(1) ^ 2;
t335 = qJD(1) * t247;
t124 = -t253 * t377 - t257 * t335;
t463 = t124 / 0.2e1;
t238 = qJDD(2) * t257;
t125 = qJDD(3) * t257 - t253 * t335 + t238;
t462 = t125 / 0.2e1;
t380 = qJD(1) * qJD(2);
t365 = t257 * t380;
t378 = qJDD(2) * t253;
t193 = -t365 - t378;
t461 = t193 / 0.2e1;
t364 = t253 * t380;
t194 = t238 - t364;
t460 = t194 / 0.2e1;
t459 = t412 / 0.2e1;
t458 = -t412 / 0.2e1;
t457 = -t198 / 0.2e1;
t456 = t198 / 0.2e1;
t455 = -t253 / 0.2e1;
t454 = t257 / 0.2e1;
t451 = -qJD(1) / 0.2e1;
t450 = qJD(1) / 0.2e1;
t440 = t240 * t95;
t437 = qJDD(1) / 0.2e1;
t436 = t116 * t257 + t120 * t414;
t155 = pkin(4) + t318;
t423 = t155 * t253;
t422 = t155 * t257;
t162 = t199 * t253;
t163 = t199 * t257;
t232 = pkin(1) + t452;
t419 = t232 * t253;
t405 = t256 * t265;
t123 = t257 * t136;
t404 = t257 * t266;
t62 = -t253 * t307 + t146;
t403 = t62 * qJD(1);
t80 = -t253 * t305 + t163;
t402 = t80 * qJD(1);
t401 = t140 * t407 + t123;
t398 = qJD(1) * t156 + t158 * t364;
t393 = t201 + t317;
t392 = t315 - t203;
t388 = qJD(1) * t313;
t386 = qJD(1) * t257;
t385 = qJD(2) * t252;
t379 = qJDD(1) * t155;
t236 = pkin(2) * t385;
t367 = t252 * t384;
t363 = -t387 / 0.2e1;
t362 = -t386 / 0.2e1;
t361 = -t384 / 0.2e1;
t360 = t384 / 0.2e1;
t359 = -t382 / 0.2e1;
t358 = t382 / 0.2e1;
t355 = -t326 - t453;
t354 = -t232 - t447;
t351 = -qJD(1) * t120 + t247 * t488;
t350 = -(-t253 * t316 - t430) * qJD(1) + t487 * t247;
t349 = -qJD(1) * t118 - t119 * t247 - t173 * t412;
t348 = -(-t253 * t314 - t428) * qJD(1) + t486 * t247;
t115 = Icges(4,5) * t414 - Icges(4,6) * t417 + t426;
t346 = t115 + t425;
t345 = -t121 - t419;
t344 = t232 * t257 + t122;
t343 = t484 * t247;
t342 = t485 * t247;
t341 = -qJD(1) * t154 - t327 * t412;
t135 = Icges(3,5) * t407 - Icges(3,6) * t409 + t427;
t340 = -t135 - t424;
t339 = t304 * qJD(4);
t109 = t205 * t253 + t212 * t257;
t336 = qJD(1) * t153 - t198 * t327;
t133 = t327 * t247;
t333 = -t133 - t373;
t332 = -t117 * t416 + t119 * t413;
t330 = g(1) * t257 + g(2) * t253;
t329 = -t137 * t408 + t139 * t406;
t210 = rSges(2,1) * t257 - rSges(2,2) * t253;
t208 = rSges(2,1) * t253 + rSges(2,2) * t257;
t206 = -t439 + t448;
t79 = -t138 * t256 - t140 * t252;
t290 = qJD(2) * t201;
t87 = t257 * t290 + (-t253 * t315 - t429) * qJD(1);
t291 = qJD(2) * t203;
t89 = t257 * t291 + (-t253 * t317 - t431) * qJD(1);
t271 = qJD(2) * t79 - t252 * t87 + t256 * t89 - t389;
t78 = -t137 * t256 - t139 * t252;
t88 = qJD(1) * t138 + t253 * t290;
t90 = qJD(1) * t140 + t253 * t291;
t272 = -qJD(1) * t135 + qJD(2) * t78 - t252 * t88 + t256 * t90;
t325 = (t253 * t475 + t257 * t272) * t257 - (t253 * t476 + t257 * t271) * t253;
t324 = (t253 * t272 - t257 * t475) * t257 - (t253 * t271 - t257 * t476) * t253;
t108 = t209 * t250 + t211 * t254;
t299 = -qJD(2) * t158 - t372;
t285 = t299 * t257;
t43 = -t108 * t248 - t155 * t387 + t285;
t301 = -t155 * t386 + t397;
t44 = t109 * t248 - t301;
t323 = -t253 * t44 - t257 * t43;
t49 = t135 * t257 - t253 * t309;
t50 = -t138 * t409 + t401;
t322 = -t253 * t50 + t257 * t49;
t51 = -t135 * t253 + t329;
t113 = t140 * t406;
t52 = -t136 * t253 - t138 * t408 + t113;
t321 = -t253 * t52 + t257 * t51;
t55 = t285 + t495;
t320 = -t253 * t54 - t257 * t55;
t66 = -t117 * t240 - t119 * t239;
t306 = -t201 * t256 - t203 * t252;
t298 = -pkin(3) * t477 - t371;
t297 = t107 * t256 - t142 * t252 + m(2) * rSges(2,1) + rSges(5,1) * m(5) + pkin(4) * m(6) + (m(3) + t494) * pkin(1);
t296 = -t116 + t311;
t295 = t328 * t253;
t98 = t120 * t413;
t293 = -t257 * t346 + t98;
t288 = t185 * qJD(3);
t284 = -qJD(1) * t312 + t145 * t198 - t146 * t412;
t283 = -(-t201 * t257 - t140) * t253 + (Icges(3,2) * t407 - t139 + t226) * t257;
t282 = (t252 * t392 - t256 * t393) * qJD(1);
t275 = -qJD(1) * t115 + t239 * t349 - t240 * t351;
t10 = t253 * t275 - t257 * t473;
t274 = t239 * t348 - t240 * t350 - t390;
t11 = t253 * t274 - t257 * t474;
t45 = t115 * t257 - t253 * t311;
t46 = -t118 * t417 + t436;
t21 = t198 * t45 - t412 * t46 + t403;
t47 = -t115 * t253 + t332;
t48 = -t116 * t253 - t118 * t416 + t98;
t63 = -t257 * t307 - t145;
t60 = t63 * qJD(1);
t22 = t198 * t47 - t412 * t48 + t60;
t277 = -t486 * t412 + (Icges(4,2) * t414 - t119 + t215) * t198 + t484 * qJD(1);
t267 = t239 * t277 - t240 * t467;
t273 = -qJD(1) * t171 + t239 * t343 - t240 * t342;
t30 = t253 * t472 + t257 * t273;
t31 = t253 * t273 - t257 * t472;
t34 = t239 * t351 + t240 * t349;
t35 = t239 * t350 + t240 * t348;
t67 = -t118 * t240 - t120 * t239;
t8 = t253 * t473 + t257 * t275;
t9 = t253 * t474 + t257 * t274;
t281 = (qJD(1) * t30 + qJDD(1) * t63 + t124 * t48 + t125 * t47 + t198 * t8 - t412 * t9) * t455 + (t239 * t467 + t240 * t277) * t451 + t21 * t363 + t22 * t362 + (qJD(1) * t31 + qJDD(1) * t62 + t10 * t198 - t11 * t412 + t124 * t46 + t125 * t45) * t454 + (-t253 * t48 + t257 * t47) * t463 + (-t253 * t46 + t257 * t45) * t462 + (-t253 * t9 + t257 * t8 + (-t253 * t47 - t257 * t48) * qJD(1)) * t458 + (t10 * t257 - t11 * t253 + (-t253 * t45 - t257 * t46) * qJD(1)) * t456 + (-t253 * t67 + t257 * t66) * t437 + (-t253 * t35 + t257 * t34 + (-t253 * t66 - t257 * t67) * qJD(1)) * t450 + (-t253 * t284 + t257 * t267) * t459 + (t253 * t267 + t257 * t284) * t457;
t168 = t184 * pkin(3);
t280 = t168 * t387 + t257 * t298;
t246 = qJDD(1) + qJDD(4);
t93 = t231 * t383 + (-t251 * t385 + t288) * pkin(3);
t99 = t288 + t477;
t276 = -qJD(2) * t93 - qJDD(2) * t158 + (-qJD(3) * t99 - qJDD(3) * t184) * pkin(3);
t269 = -t155 * t266 + t276;
t64 = qJD(1) * t304 + t339;
t92 = -t231 * t385 - t490;
t12 = t109 * t246 + t248 * t64 + (t379 + (t299 + t92) * qJD(1)) * t257 + t269 * t253;
t13 = -t108 * t246 + t248 * t497 + (-qJD(1) * t92 - t379) * t253 + t269 * t257 + t398;
t279 = -t12 * t253 - t13 * t257 + (t253 * t43 - t257 * t44) * qJD(1);
t56 = t231 * t367 + t253 * t490 - t496;
t23 = qJD(1) * t56 + qJDD(1) * t105 + t257 * t276 + t398;
t91 = t92 * t257;
t57 = t91 + t495;
t24 = qJDD(1) * t104 + (t57 + t285) * qJD(1) + t276 * t253;
t278 = -t23 * t257 - t24 * t253 + (t253 * t55 - t257 * t54) * qJD(1);
t58 = qJD(1) * t345 - t198 * t326 - t236 * t257;
t182 = t315 * qJD(2);
t183 = t317 * qJD(2);
t270 = -qJD(1) * t199 + qJD(2) * t306 + t182 * t252 - t183 * t256;
t268 = t252 * t283 + t256 * t468;
t234 = Icges(6,3) * t246;
t221 = pkin(2) * t367;
t219 = t387 * t439;
t195 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3);
t187 = rSges(3,1) * t383 - rSges(3,2) * t385;
t127 = g(3) * t334 + t330 * t464;
t103 = -g(3) * t464 + t330 * t334;
t81 = -t257 * t305 - t162;
t77 = t81 * qJD(1);
t68 = -t168 * t386 + t253 * t298;
t61 = -t158 * t386 + t253 * t300;
t59 = qJD(1) * t344 - t326 * t412 - t221;
t42 = qJD(1) * t219 + qJDD(1) * t144 + t193 * t328 + (-qJD(2) * t187 + t266 * t366) * t253 + (qJDD(1) * pkin(1) + qJD(1) * (-rSges(3,1) * t385 - rSges(3,2) * t383 - rSges(3,3) * qJD(1))) * t257;
t41 = -t187 * t382 - qJDD(1) * t143 - t194 * t328 + (-qJDD(1) * t253 - t404) * pkin(1) + ((-t206 * t257 + t445) * qJD(1) + qJD(2) * t295) * qJD(1);
t39 = t253 * t270 - t257 * t471;
t38 = t253 * t471 + t257 * t270;
t37 = qJD(2) * t308 - t252 * t89 - t256 * t87;
t36 = qJD(2) * t309 - t252 * t90 - t256 * t88;
t33 = -t232 * t404 - qJD(1) * t76 - t125 * t326 - t133 * t198 + t345 * qJDD(1) + (-t257 * t405 + (0.2e1 * t364 - t238) * t252) * pkin(2);
t32 = -t266 * t419 + qJD(1) * t75 + t124 * t326 - t133 * t412 + t344 * qJDD(1) + (-t253 * t405 + (-0.2e1 * t365 - t378) * t252) * pkin(2);
t29 = qJD(2) * t321 + t77;
t28 = qJD(2) * t322 + t402;
t25 = t121 * t124 - t122 * t125 - t198 * t75 - t412 * t76 + t331;
t1 = [(qJD(2) * t305 + t182 * t256 + t183 * t252) * qJD(1) + (g(1) * t195 - g(2) * t297) * t257 + t253 * (g(1) * t297 + g(2) * t195) + (t239 * t342 + t240 * t343) * qJD(1) + t234 + (t77 + ((t329 - t50 + t401) * t257 + (-t113 - t49 + (t136 - t309) * t253) * t253) * qJD(2)) * t359 + m(4) * (t59 * t196 + t33 * t217 + t32 * t218 + t58 * t221 + (t33 * t354 - t32 * rSges(4,3) + t58 * t480 + (rSges(4,3) * t58 + t354 * t59) * qJD(1)) * t253 + (-t33 * rSges(4,3) + t32 * (t232 - t446) + t59 * (-t236 - t480) + (t58 * (-t232 - t327) - t59 * rSges(4,3)) * qJD(1)) * t257) + (t60 + (t332 - t46 + t436) * t198 - (t45 + t293) * t412 + (-t198 * t346 - t296 * t412) * t253) * t457 + (t67 + t63) * t463 + (t66 + t62) * t462 + (t79 + t81) * t461 + (t78 + t80) * t460 + (-t403 + (-t48 + t293) * t198 - (t257 * t311 - t436 + t47) * t412 + (t198 * t296 - t346 * t412) * t253 + t21) * t459 + (t35 + t30) * t458 + (t37 + t38) * t361 + (-t402 + ((t257 * t340 + t113 - t52) * t257 + (t253 * t340 - t123 + t401 - t51) * t253) * qJD(2) + t28) * t360 + (t13 * (-t108 - t423) + t12 * (t109 + t422) + (-(-qJD(1) * t155 - t212 * t248) * t253 - (t205 * t248 + t299) * t257 + t339 + t91 + (t304 - t423) * qJD(1)) * t44 + (-t422 * qJD(1) - t253 * t92 - t301) * t43 + t469) * m(6) + (t104 * t24 + t105 * t23 + (t56 + t54) * t55 + (t57 + t151 * t387 - (-rSges(5,3) * qJD(1) + t299) * t257) * t54) * m(5) + (t41 * (t253 * t366 + t353) + t42 * ((pkin(1) - t439) * t257 + t352) + t84 * t219 + (t295 * t83 - t328 * t438 + t489) * qJD(2) + (t356 * t83 - t357 * t84 + (rSges(3,3) * t83 + t366 * t84) * t253 + (t83 * (-pkin(1) - t206) - t84 * rSges(3,3)) * t257) * qJD(1)) * m(3) + (t34 + t31 + t22) * t456 + (t36 + t39 + t29) * t358 + (m(2) * (t208 ^ 2 + t210 ^ 2) + t306 - t173 * t240 - t175 * t239 + Icges(5,2) + Icges(2,3)) * qJDD(1); ((-t49 * t253 - t50 * t257) * qJD(1) + t324) * t358 + (-t253 * t37 + t257 * t36 + (-t253 * t78 - t257 * t79) * qJD(1)) * t450 + ((-t51 * t253 - t52 * t257) * qJD(1) + t325) * t361 + ((t162 * t382 - t388) * t257 + (t282 + (-t163 * t257 + t268) * qJD(2)) * t253) * t359 + ((t163 * t384 + t388) * t253 + (t282 + (-t162 * t253 + t268) * qJD(2)) * t257) * t360 + t321 * t461 + t322 * t460 + (-t253 * t79 + t257 * t78) * t437 + t281 + (qJD(1) * t39 + qJD(2) * t324 + qJDD(1) * t80 + t193 * t50 + t194 * t49) * t454 + (-g(3) * t142 + t107 * t330) * t252 + (g(3) * t107 + t142 * t330) * t256 + ((t252 * t393 + t256 * t392) * qJD(1) + (-t252 * t468 + t283 * t256) * qJD(2)) * t451 + (qJD(1) * t38 + qJD(2) * t325 + qJDD(1) * t81 + t193 * t52 + t194 * t51) * t455 + t28 * t363 + t29 * t362 + (t158 * t279 + t323 * t93 - t43 * t478 - t44 * t61 + t482) * m(6) + (t158 * t278 + t320 * t93 - t478 * t55 - t54 * t61 + t482) * m(5) + (-t58 * (-t257 * t373 + t336) - t59 * ((-t252 * t386 - t253 * t383) * pkin(2) + t341) - t25 * t452 + (qJD(1) * t326 * t58 - t25 * t121 + t32 * t355 + t333 * t59) * t253 + (-t25 * t122 + t58 * t333 + (qJD(1) * t59 + t33) * t355) * t257 + t466) * m(4) + ((-t253 * t84 - t257 * t83) * t187 - (qJD(1) * t319 + t253 * t42 + t257 * t41) * t328 + qJD(1) * t489 + (qJDD(2) * t206 - t265 * t328 + t83 * t382 + t84 * t384) * t206) * m(3); (t103 * t255 - t127 * t251) * t252 + (t103 * t251 + t127 * t255) * t256 + t281 + (-t336 * t58 - t341 * t59 + t25 * (-t121 * t253 - t122 * t257) - (t253 * t59 + t257 * t58) * t133 - (t32 * t253 + t33 * t257 + (-t253 * t58 + t257 * t59) * qJD(1)) * t326 + t466) * m(4) + ((t184 * t279 + t323 * t99 - t440) * pkin(3) - t280 * t43 - t44 * t68) * m(6) + ((t184 * t278 + t320 * t99 - t440) * pkin(3) - t280 * t55 - t54 * t68) * m(5); t234 + (-t108 * t13 + t109 * t12 + t43 * t497 + t44 * t64 - (t304 * t44 + t337 * t43) * t248 + t469) * m(6);];
tau = t1;
