% Calculate vector of inverse dynamics joint torques for
% palh2m2DE
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh2m2DE_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:43
% DurationCPUTime: 13.25s
% Computational Cost: add. (3459->594), mult. (8241->851), div. (0->0), fcn. (6032->8), ass. (0->341)
t228 = qJD(1) + qJD(4);
t232 = sin(qJ(1));
t236 = cos(qJ(1));
t181 = rSges(7,1) * t236 - rSges(7,2) * t232;
t183 = rSges(7,1) * t232 + rSges(7,2) * t236;
t229 = sin(qJ(4));
t233 = cos(qJ(4));
t312 = -t181 * t233 + t183 * t229;
t461 = t228 * t312;
t231 = sin(qJ(2));
t435 = pkin(4) * t231;
t334 = qJD(2) * t435;
t230 = sin(qJ(3));
t346 = qJD(3) * t230;
t336 = pkin(5) * t346;
t155 = t334 + t336;
t120 = t155 * t236;
t235 = cos(qJ(2));
t434 = pkin(4) * t235;
t319 = pkin(1) + t434;
t191 = pkin(2) + t319;
t234 = cos(qJ(3));
t433 = pkin(5) * t234;
t265 = t191 + t433;
t137 = rSges(6,1) + t265;
t95 = t236 * rSges(6,3) - t137 * t232;
t459 = qJD(1) * t95 - t120;
t189 = t232 * t336;
t190 = t232 * t334;
t363 = t190 + t189;
t94 = rSges(6,3) * t232 + t137 * t236;
t404 = qJD(1) * t94;
t69 = -t363 + t404;
t175 = -t229 * rSges(7,1) - t233 * rSges(7,2);
t184 = rSges(7,1) * t233 - rSges(7,2) * t229;
t281 = t175 * t236 - t184 * t232;
t192 = rSges(4,1) + t319;
t118 = t236 * rSges(4,3) - t192 * t232;
t458 = qJD(1) * t118 - t236 * t334;
t350 = qJD(2) * t235;
t210 = pkin(4) * t350;
t344 = qJD(3) * t234;
t156 = pkin(5) * t344 + t210;
t239 = qJD(3) ^ 2;
t208 = qJDD(2) * t434;
t240 = qJD(2) ^ 2;
t306 = -t240 * t435 + t208;
t89 = (qJDD(3) * t234 - t230 * t239) * pkin(5) + t306;
t457 = t156 * t336 + t89 * t433;
t456 = t156 * t334 + t89 * t434;
t217 = Icges(5,4) * t234;
t290 = -Icges(5,2) * t230 + t217;
t169 = Icges(5,1) * t230 + t217;
t218 = Icges(3,4) * t235;
t291 = -Icges(3,2) * t231 + t218;
t171 = Icges(3,1) * t231 + t218;
t164 = Icges(3,5) * t235 - Icges(3,6) * t231;
t163 = Icges(3,5) * t231 + Icges(3,6) * t235;
t269 = qJD(2) * t163;
t414 = Icges(3,4) * t231;
t172 = Icges(3,1) * t235 - t414;
t107 = Icges(3,5) * t232 + t172 * t236;
t103 = Icges(3,6) * t232 + t236 * t291;
t397 = t103 * t231;
t286 = -t107 * t235 + t397;
t406 = Icges(3,3) * t236;
t453 = -t236 * t269 + (-t164 * t232 + t286 + t406) * qJD(1);
t383 = t231 * t232;
t199 = Icges(3,4) * t383;
t380 = t232 * t235;
t412 = Icges(3,5) * t236;
t106 = Icges(3,1) * t380 - t199 - t412;
t408 = Icges(3,6) * t236;
t102 = Icges(3,4) * t380 - Icges(3,2) * t383 - t408;
t398 = t102 * t231;
t287 = -t106 * t235 + t398;
t99 = Icges(3,3) * t232 + t164 * t236;
t401 = qJD(1) * t99;
t452 = qJD(1) * t287 - t232 * t269 + t401;
t162 = Icges(5,5) * t234 - Icges(5,6) * t230;
t161 = Icges(5,5) * t230 + Icges(5,6) * t234;
t266 = qJD(3) * t161;
t413 = Icges(5,4) * t230;
t170 = Icges(5,1) * t234 - t413;
t105 = Icges(5,5) * t232 + t170 * t236;
t101 = Icges(5,6) * t232 + t236 * t290;
t399 = t101 * t230;
t288 = -t105 * t234 + t399;
t405 = Icges(5,3) * t236;
t451 = -t236 * t266 + (-t162 * t232 + t288 + t405) * qJD(1);
t385 = t230 * t232;
t198 = Icges(5,4) * t385;
t381 = t232 * t234;
t411 = Icges(5,5) * t236;
t104 = Icges(5,1) * t381 - t198 - t411;
t407 = Icges(5,6) * t236;
t100 = Icges(5,4) * t381 - Icges(5,2) * t385 - t407;
t400 = t100 * t230;
t289 = -t104 * t234 + t400;
t97 = Icges(5,3) * t232 + t162 * t236;
t402 = qJD(1) * t97;
t450 = qJD(1) * t289 - t232 * t266 + t402;
t98 = Icges(3,5) * t380 - Icges(3,6) * t383 - t406;
t31 = -t232 * t287 - t236 * t98;
t96 = Icges(5,5) * t381 - Icges(5,6) * t385 - t405;
t29 = -t232 * t289 - t236 * t96;
t165 = Icges(5,2) * t234 + t413;
t284 = t165 * t230 - t169 * t234;
t449 = t284 * qJD(1) + t162 * qJD(3);
t167 = Icges(3,2) * t235 + t414;
t282 = t167 * t231 - t171 * t235;
t448 = t282 * qJD(1) + t164 * qJD(2);
t204 = rSges(7,1) * g(1) + rSges(7,2) * g(2);
t205 = rSges(7,1) * g(2) - rSges(7,2) * g(1);
t447 = t232 * (t204 * t233 + t205 * t229) + (t204 * t229 - t205 * t233) * t236;
t446 = t232 * (-t167 * t236 + t107) - t236 * (-Icges(3,2) * t380 + t106 - t199);
t445 = t232 * (-t165 * t236 + t105) - t236 * (-Icges(5,2) * t381 + t104 - t198);
t241 = qJD(1) ^ 2;
t444 = m(6) + m(7);
t443 = m(3) * rSges(3,2);
t442 = m(5) * rSges(5,2);
t340 = qJD(1) * qJD(3);
t157 = qJDD(3) * t232 + t236 * t340;
t441 = t157 / 0.2e1;
t158 = -qJDD(3) * t236 + t232 * t340;
t440 = t158 / 0.2e1;
t341 = qJD(1) * qJD(2);
t159 = qJDD(2) * t232 + t236 * t341;
t439 = t159 / 0.2e1;
t160 = -qJDD(2) * t236 + t232 * t341;
t438 = t160 / 0.2e1;
t437 = t232 / 0.2e1;
t436 = -t236 / 0.2e1;
t431 = -qJD(1) / 0.2e1;
t430 = qJD(1) / 0.2e1;
t379 = t234 * t236;
t429 = -t104 * t379 - t232 * t96;
t428 = t105 * t379 + t232 * t97;
t378 = t235 * t236;
t427 = -t106 * t378 - t232 * t98;
t426 = t107 * t378 + t232 * t99;
t425 = rSges(3,1) * t235;
t424 = rSges(5,1) * t234;
t423 = rSges(3,2) * t231;
t422 = rSges(5,2) * t230;
t221 = t232 * rSges(3,3);
t335 = rSges(3,1) * t380;
t361 = rSges(3,2) * t383 + t236 * rSges(3,3);
t109 = t335 - t361;
t177 = rSges(3,1) * t231 + rSges(3,2) * t235;
t349 = qJD(2) * t236;
t55 = -t177 * t349 + (-pkin(1) * t232 - t109) * qJD(1);
t421 = t236 * t55;
t359 = rSges(3,1) * t378 + t221;
t382 = t231 * t236;
t111 = -rSges(3,2) * t382 + t359;
t351 = qJD(2) * t232;
t56 = -t177 * t351 + (pkin(1) * t236 + t111) * qJD(1);
t420 = t236 * t56;
t417 = qJDD(1) / 0.2e1;
t135 = pkin(3) + t265;
t396 = t135 * t232;
t395 = t135 * t236;
t393 = t155 * t232;
t392 = t161 * t232;
t391 = t161 * t236;
t390 = t163 * t232;
t389 = t163 * t236;
t384 = t230 * t236;
t377 = t235 * t240;
t51 = -t232 * t284 - t391;
t376 = t51 * qJD(1);
t52 = -t232 * t282 - t389;
t375 = t52 * qJD(1);
t370 = -t165 + t170;
t369 = t169 + t290;
t368 = -t167 + t172;
t367 = t171 + t291;
t174 = qJD(1) * t190;
t366 = qJD(1) * t189 + t174;
t354 = qJD(1) * t232;
t333 = t230 * t354;
t353 = qJD(1) * t236;
t365 = rSges(5,2) * t333 + rSges(5,3) * t353;
t332 = t231 * t354;
t364 = rSges(3,2) * t332 + rSges(3,3) * t353;
t362 = rSges(5,2) * t385 + t236 * rSges(5,3);
t360 = rSges(5,1) * t379 + t232 * rSges(5,3);
t357 = qJD(1) * t155;
t356 = qJD(1) * t162;
t355 = qJD(1) * t164;
t352 = qJD(2) * t156;
t348 = qJD(3) * t156;
t176 = rSges(5,1) * t230 + rSges(5,2) * t234;
t347 = qJD(3) * t176;
t179 = -t422 + t424;
t148 = qJD(3) * t179;
t345 = qJD(3) * t232;
t343 = qJD(3) * t236;
t339 = qJDD(1) * t135;
t338 = t191 * qJDD(1);
t337 = pkin(1) + pkin(2);
t331 = t354 / 0.2e1;
t330 = t353 / 0.2e1;
t329 = -t351 / 0.2e1;
t328 = t351 / 0.2e1;
t327 = -t349 / 0.2e1;
t326 = t349 / 0.2e1;
t325 = -t345 / 0.2e1;
t324 = t345 / 0.2e1;
t323 = -t343 / 0.2e1;
t322 = t343 / 0.2e1;
t83 = t105 * t381;
t321 = t236 * t97 - t83;
t84 = t107 * t380;
t320 = t236 * t99 - t84;
t318 = -t191 - t424;
t317 = -t96 + t399;
t316 = -t98 + t397;
t315 = -pkin(1) - t425;
t314 = t281 * qJD(4);
t82 = t175 * t232 + t184 * t236;
t117 = rSges(4,3) * t232 + t192 * t236;
t79 = qJD(1) * t117 - t190;
t307 = pkin(5) * t444 + m(5) * rSges(5,1);
t305 = g(1) * t236 + g(2) * t232;
t48 = t101 * t234 + t105 * t230;
t267 = qJD(3) * t165;
t61 = -t236 * t267 + (-t232 * t290 + t407) * qJD(1);
t268 = qJD(3) * t169;
t65 = -t236 * t268 + (-t170 * t232 + t411) * qJD(1);
t250 = -qJD(3) * t48 - t230 * t61 + t234 * t65 + t402;
t47 = t100 * t234 + t104 * t230;
t62 = qJD(1) * t101 - t232 * t267;
t66 = qJD(1) * t105 - t232 * t268;
t251 = qJD(1) * t96 - qJD(3) * t47 - t230 * t62 + t234 * t66;
t304 = -(t232 * t450 + t251 * t236) * t236 + (t232 * t451 + t250 * t236) * t232;
t50 = t103 * t235 + t107 * t231;
t270 = qJD(2) * t167;
t63 = -t236 * t270 + (-t232 * t291 + t408) * qJD(1);
t271 = qJD(2) * t171;
t67 = -t236 * t271 + (-t172 * t232 + t412) * qJD(1);
t248 = -qJD(2) * t50 - t231 * t63 + t235 * t67 + t401;
t49 = t102 * t235 + t106 * t231;
t64 = qJD(1) * t103 - t232 * t270;
t68 = qJD(1) * t107 - t232 * t271;
t249 = qJD(1) * t98 - qJD(2) * t49 - t231 * t64 + t235 * t68;
t303 = t232 * (t232 * t453 + t248 * t236) - t236 * (t232 * t452 + t249 * t236);
t302 = t232 * (t250 * t232 - t236 * t451) - t236 * (t251 * t232 - t236 * t450);
t301 = t232 * (t248 * t232 - t236 * t453) - t236 * (t249 * t232 - t236 * t452);
t182 = rSges(2,1) * t236 - rSges(2,2) * t232;
t178 = rSges(2,1) * t232 + rSges(2,2) * t236;
t180 = -t423 + t425;
t30 = -t101 * t385 - t321;
t300 = t232 * t30 - t236 * t29;
t32 = -t103 * t383 - t320;
t299 = t232 * t32 - t236 * t31;
t33 = -t100 * t384 - t429;
t34 = -t101 * t384 + t428;
t298 = t232 * t34 - t236 * t33;
t35 = -t102 * t382 - t427;
t36 = -t103 * t382 + t426;
t297 = t232 * t36 - t236 * t35;
t81 = t181 * t229 + t183 * t233;
t39 = -t135 * t354 - t228 * t81 - t120;
t277 = -t135 * t353 + t363;
t40 = t228 * t82 - t277;
t296 = -t232 * t40 - t236 * t39;
t108 = rSges(5,1) * t381 - t362;
t43 = (-t334 - t347) * t236 + (-t191 * t232 - t108) * qJD(1);
t110 = -rSges(5,2) * t384 + t360;
t44 = -t176 * t345 - t190 + (t191 * t236 + t110) * qJD(1);
t295 = -t232 * t44 - t236 * t43;
t294 = -t232 * t69 - t236 * t459;
t285 = t165 * t234 + t169 * t230;
t283 = t167 * t235 + t171 * t231;
t278 = (m(4) + m(5) + t444) * pkin(4) + m(3) * rSges(3,1);
t276 = t232 * t177;
t275 = qJD(1) * t177;
t274 = qJD(1) * t176;
t272 = (-qJDD(2) * t231 - t377) * pkin(4);
t264 = -rSges(5,1) * t346 - rSges(5,2) * t344;
t263 = t236 * t275;
t262 = t100 * t236 - t101 * t232;
t261 = t102 * t236 - t103 * t232;
t259 = (-t230 * t369 + t234 * t370) * qJD(1);
t258 = (-t231 * t367 + t235 * t368) * qJD(1);
t253 = -qJD(3) * t148 - t241 * t191 + t272;
t15 = -t157 * t176 + qJDD(1) * t110 + qJD(1) * t365 + (t338 + (-t347 - 0.2e1 * t334) * qJD(1)) * t236 + (-t241 * t424 + t253) * t232;
t16 = -qJDD(1) * t108 + t158 * t176 + 0.2e1 * t174 + (-t338 - qJD(1) * (rSges(5,3) * qJD(1) + t264)) * t232 + (-t179 * t241 + t253) * t236;
t257 = -t15 * t232 - t16 * t236 + (t232 * t43 - t236 * t44) * qJD(1);
t227 = qJDD(1) + qJDD(4);
t252 = (-qJDD(3) * t230 - t234 * t239) * pkin(5) + t272;
t245 = -t135 * t241 + t252;
t41 = qJD(1) * t281 + t314;
t13 = t227 * t82 + t228 * t41 + (t339 - 0.2e1 * t357) * t236 + t245 * t232;
t14 = -t227 * t81 + t228 * t461 + (-t339 + t357) * t232 + t245 * t236 + t366;
t256 = (t232 * t39 - t236 * t40) * qJD(1) - t13 * t232 - t14 * t236;
t71 = -t404 + t393;
t25 = qJD(1) * t71 + qJDD(1) * t95 + t236 * t252 + t366;
t26 = qJDD(1) * t94 + (t459 - t120) * qJD(1) + t252 * t232;
t255 = (t232 * t459 - t236 * t69) * qJD(1) - t232 * t26 - t236 * t25;
t254 = m(3) * t423 - t278 * t235 - t307 * t234 - (pkin(3) + t337) * m(7) - (rSges(6,1) + t337) * m(6) - (m(3) + m(4)) * pkin(1) - m(2) * rSges(2,1) - m(4) * rSges(4,1) + (t422 - t337) * m(5);
t142 = t290 * qJD(3);
t144 = t170 * qJD(3);
t247 = qJD(1) * t161 - qJD(3) * t285 - t142 * t230 + t144 * t234;
t143 = t291 * qJD(2);
t145 = t172 * qJD(2);
t246 = qJD(1) * t163 - qJD(2) * t283 - t143 * t231 + t145 * t235;
t244 = -t230 * t445 + t262 * t234;
t243 = -t231 * t446 + t261 * t235;
t209 = Icges(7,3) * t227;
t150 = qJD(2) * t180;
t134 = m(2) * rSges(2,2) - m(3) * rSges(3,3) - m(4) * rSges(4,3) - m(5) * rSges(5,3) - m(6) * rSges(6,3);
t116 = (-t231 * t353 - t232 * t350) * pkin(4);
t115 = (-t235 * t349 + t332) * pkin(4);
t114 = (-t230 * t353 - t232 * t344) * pkin(5);
t113 = (-t234 * t343 + t333) * pkin(5);
t112 = t210 + t148;
t73 = -qJD(3) * t347 + qJDD(3) * t179 + t306;
t54 = -t236 * t282 + t390;
t53 = -t236 * t284 + t392;
t46 = t54 * qJD(1);
t45 = t53 * qJD(1);
t38 = -t79 * qJD(1) + qJDD(1) * t118 + t236 * t272 + t174;
t37 = qJD(1) * t458 + qJDD(1) * t117 + (-t159 * t231 - t232 * t377) * pkin(4);
t28 = -t159 * t177 + qJDD(1) * t111 + qJD(1) * (-qJD(1) * t335 + t364) + (-t232 * t150 - t263) * qJD(2) + (qJDD(1) * t236 - t232 * t241) * pkin(1);
t27 = -t150 * t349 - qJDD(1) * t109 + t160 * t177 + (-qJDD(1) * t232 - t236 * t241) * pkin(1) + ((-t180 * t236 - t221) * qJD(1) + qJD(2) * t276) * qJD(1);
t24 = t246 * t232 - t236 * t448;
t23 = t247 * t232 - t236 * t449;
t22 = t232 * t448 + t246 * t236;
t21 = t232 * t449 + t247 * t236;
t20 = -qJD(2) * t286 + t231 * t67 + t235 * t63;
t19 = -t287 * qJD(2) + t231 * t68 + t235 * t64;
t18 = -qJD(3) * t288 + t230 * t65 + t234 * t61;
t17 = -t289 * qJD(3) + t230 * t66 + t234 * t62;
t12 = qJD(2) * t297 + t46;
t11 = qJD(3) * t298 + t45;
t10 = qJD(2) * t299 + t375;
t9 = qJD(3) * t300 + t376;
t1 = [(t45 + ((t30 - t83 + (t97 + t400) * t236 + t429) * t236 + t428 * t232) * qJD(3)) * t322 + t232 * (-g(1) * t254 + g(2) * t134) + (g(1) * t134 + g(2) * t254) * t236 + m(5) * (t16 * t362 + t43 * t190 + t15 * t360 + t44 * t365 + (t16 * t318 + t43 * t347 + (-t43 * rSges(5,3) + t318 * t44) * qJD(1)) * t232 + (t15 * (t191 - t422) + t44 * (t264 - t334) + t43 * (-t179 - t191) * qJD(1)) * t236) + (-qJD(2) * t282 + t143 * t235 + t145 * t231) * qJD(1) + (-qJD(3) * t284 + t142 * t234 + t144 * t230) * qJD(1) + m(3) * (t27 * (t232 * t315 + t361) + t28 * ((pkin(1) - t423) * t236 + t359) + t56 * t364 + (-t177 * t420 + t276 * t55) * qJD(2) + ((-pkin(1) - t180) * t421 + (-t55 * rSges(3,3) + t315 * t56) * t232) * qJD(1)) + (t46 + ((t32 - t84 + (t99 + t398) * t236 + t427) * t236 + t426 * t232) * qJD(2)) * t326 + t209 + (t48 + t53) * t441 + (t47 + t51) * t440 + (t50 + t54) * t439 + (t49 + t52) * t438 + (t10 - t375 + ((t236 * t316 + t36 - t426) * t236 + (t232 * t316 + t320 + t35) * t232) * qJD(2)) * t329 + (t20 + t22) * t328 + (-t376 + ((t236 * t317 + t34 - t428) * t236 + (t232 * t317 + t321 + t33) * t232) * qJD(3) + t9) * t325 + (t18 + t21) * t324 + (t14 * (-t81 - t396) + t13 * (t82 + t395) + (-(-qJD(1) * t135 - t184 * t228) * t232 - (t175 * t228 - t155) * t236 - t120 + t314 + (t281 - t396) * qJD(1)) * t40 + (-t395 * qJD(1) - t277 + t393) * t39 + t447) * m(7) + m(6) * (t25 * t95 + t26 * t94 + (t71 + t69) * t459) + (t19 + t24 + t12) * t327 + (t23 + t11 + t17) * t323 + (t117 * t37 + t118 * t38) * m(4) + (m(2) * (t178 ^ 2 + t182 ^ 2) + t283 + t285 + Icges(6,2) + Icges(2,3) + Icges(4,2)) * qJDD(1); (qJD(1) * t24 + qJD(2) * t301 + qJDD(1) * t52 + t159 * t32 + t160 * t31) * t436 + (qJD(1) * t22 + qJD(2) * t303 + qJDD(1) * t54 + t159 * t36 + t160 * t35) * t437 + t299 * t438 + t297 * t439 + (t232 * t50 - t236 * t49) * t417 + (-t19 * t236 + t20 * t232 + (t49 * t232 + t236 * t50) * qJD(1)) * t430 + ((t231 * t368 + t235 * t367) * qJD(1) + (t261 * t231 + t235 * t446) * qJD(2)) * t431 + t12 * t330 + t10 * t331 + ((-t349 * t390 - t355) * t236 + (t258 + (t236 * t389 + t243) * qJD(2)) * t232) * t326 + ((-t351 * t389 + t355) * t232 + (t258 + (t232 * t390 + t243) * qJD(2)) * t236) * t329 - t235 * (g(3) * t278 - t305 * t443) + (g(3) * t443 + t278 * t305) * t231 + ((t35 * t232 + t36 * t236) * qJD(1) + t303) * t328 + ((t31 * t232 + t32 * t236) * qJD(1) + t301) * t327 + ((t296 * t350 + (t256 - t352) * t231) * pkin(4) - t115 * t39 - t116 * t40 + t456) * m(7) + ((t294 * t350 + (t255 - t352) * t231) * pkin(4) - t115 * t459 - t116 * t69 + t456) * m(6) + (((qJD(2) * t295 + t73) * t235 + (-qJD(2) * t112 + t257) * t231) * pkin(4) + t112 * t334 - t115 * t43 - t116 * t44) * m(5) + (((t208 + (-t232 * t79 - t236 * t458) * qJD(2)) * t235 + (-pkin(4) * t377 - t232 * t37 - t236 * t38 + (t232 * t458 - t236 * t79) * qJD(1)) * t231) * pkin(4) - t115 * t458 - t116 * t79) * m(4) + (-t55 * (-t180 * t349 + t232 * t275) - t56 * (-t180 * t351 - t263) + (-t232 * t56 - t421) * t150 + qJDD(2) * t180 ^ 2 + (-t28 * t232 - t27 * t236 + (t232 * t55 - t420) * qJD(1) - t180 * t240) * t177) * m(3); t11 * t330 + (qJD(1) * t21 + qJD(3) * t304 + qJDD(1) * t53 + t157 * t34 + t158 * t33) * t437 + t298 * t441 + ((t33 * t232 + t34 * t236) * qJD(1) + t304) * t324 + t9 * t331 + (qJD(1) * t23 + qJD(3) * t302 + qJDD(1) * t51 + t157 * t30 + t158 * t29) * t436 + t300 * t440 + ((t29 * t232 + t30 * t236) * qJD(1) + t302) * t323 + (t232 * t48 - t236 * t47) * t417 + (-t17 * t236 + t18 * t232 + (t47 * t232 + t236 * t48) * qJD(1)) * t430 + ((-t345 * t391 + t356) * t232 + (t259 + (t232 * t392 + t244) * qJD(3)) * t236) * t325 + ((-t343 * t392 - t356) * t236 + (t259 + (t236 * t391 + t244) * qJD(3)) * t232) * t322 + ((t230 * t370 + t234 * t369) * qJD(1) + (t262 * t230 + t234 * t445) * qJD(3)) * t431 + (g(3) * t442 + t305 * t307) * t230 - (g(3) * t307 - t305 * t442) * t234 + ((t296 * t344 + (t256 - t348) * t230) * pkin(5) - t113 * t39 - t114 * t40 + t457) * m(7) + ((t294 * t344 + (t255 - t348) * t230) * pkin(5) - t113 * t459 - t114 * t69 + t457) * m(6) + (t295 * t148 + t257 * t176 + t179 * t73 - t43 * (-t179 * t343 + t232 * t274) - t44 * (-t179 * t345 - t236 * t274)) * m(5); t209 + (t13 * t82 - t14 * t81 + t39 * t461 + t40 * t41 - (t40 * t281 + t39 * t312) * t228 + t447) * m(7);];
tau = t1;
