% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnTE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnTE_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:38
% EndTime: 2020-04-12 19:19:16
% DurationCPUTime: 20.62s
% Computational Cost: add. (195946->428), mult. (272584->808), div. (11160->14), fcn. (79629->6), ass. (0->323)
t260 = sin(qJ(1));
t262 = cos(qJ(1));
t435 = t262 * rSges(5,3);
t268 = pkin(2) ^ 2;
t269 = pkin(1) ^ 2;
t261 = cos(qJ(2));
t448 = pkin(2) * t261;
t381 = -0.2e1 * pkin(1) * t448 + t269;
t241 = t268 + t381;
t234 = 0.1e1 / t241;
t265 = 0.1e1 / pkin(4);
t407 = t234 * t265;
t379 = pkin(3) ^ 2 - pkin(4) ^ 2;
t228 = t241 - t379;
t259 = sin(qJ(2));
t449 = pkin(2) * t259;
t213 = t228 * t449;
t247 = pkin(1) - t448;
t469 = -pkin(3) - pkin(4);
t223 = (pkin(2) - t469) * (pkin(2) + t469) + t381;
t468 = pkin(4) - pkin(3);
t224 = (pkin(2) - t468) * (pkin(2) + t468) + t381;
t409 = t223 * t224;
t270 = sqrt(-t409);
t187 = t247 * t270 + t213;
t459 = t187 / 0.2e1;
t397 = t259 * t270;
t196 = pkin(2) * t397;
t186 = t228 * t247 - t196;
t460 = t186 / 0.2e1;
t488 = (rSges(5,1) * t460 + rSges(5,2) * t459) * t407;
t500 = pkin(1) - t488;
t119 = -t260 * t500 + t435;
t439 = t260 * rSges(5,3);
t120 = t262 * t500 + t439;
t267 = 0.1e1 / pkin(3);
t227 = t241 + t379;
t215 = pkin(1) * t259 * t227;
t451 = pkin(1) * t261;
t248 = -pkin(2) + t451;
t188 = -t248 * t270 + t215;
t256 = t259 ^ 2;
t235 = 0.1e1 / t241 ^ 2;
t450 = pkin(2) * t235;
t377 = pkin(1) * t450;
t195 = pkin(1) * t397;
t185 = -t227 * t248 - t195;
t393 = t261 * t185;
t286 = (-t188 * t256 + t259 * t393) * t377;
t390 = t261 * t270;
t339 = pkin(1) * pkin(2) * (-t223 - t224);
t194 = t259 * t339;
t198 = 0.1e1 / t270;
t414 = t198 * t194;
t478 = -0.2e1 * t248;
t137 = t215 + (-t390 + (pkin(2) * t478 - t414) * t259) * pkin(1);
t501 = t188 / 0.2e1;
t349 = t501 - t137 / 0.2e1;
t405 = t256 * t269;
t369 = pkin(2) * t405;
t385 = t227 * t451 + t195;
t412 = t198 * t248;
t139 = -t194 * t412 + 0.2e1 * t369 + t385;
t350 = t139 / 0.2e1 + t185 / 0.2e1;
t56 = (t286 + (t259 * t350 + t261 * t349) * t234) * t267;
t51 = t260 * t56;
t392 = t261 * t188;
t489 = (t185 * t256 + t259 * t392) * t377;
t280 = t267 * (-t489 + (-t259 * t349 + t261 * t350) * t234);
t52 = t260 * t280;
t324 = rSges(4,1) * t51 + rSges(4,2) * t52;
t400 = t259 * t260;
t375 = pkin(2) * t400;
t20 = -t324 + t375;
t398 = t259 * t262;
t374 = pkin(2) * t398;
t386 = t262 * t267;
t360 = t234 * t386;
t402 = t259 * t185;
t136 = (t392 + t402) * t360 / 0.2e1;
t403 = t259 * t139;
t455 = -t261 / 0.2e1;
t53 = ((t137 * t455 + t403 / 0.2e1) * t234 + t286) * t386 + t136;
t54 = t262 * t280;
t443 = t53 * rSges(4,1) + t54 * rSges(4,2);
t21 = -t374 + t443;
t396 = t260 * t261;
t211 = rSges(3,1) * t396 - rSges(3,2) * t400 - t262 * rSges(3,3);
t391 = t261 * t262;
t212 = rSges(3,1) * t391 - rSges(3,2) * t398 + rSges(3,3) * t260;
t237 = rSges(3,1) * t259 + rSges(3,2) * t261;
t225 = t237 * t260;
t226 = t237 * t262;
t406 = t234 * t267;
t496 = t259 * t501 - t393 / 0.2e1;
t142 = t496 * t406;
t404 = t259 * t137;
t55 = ((-t404 / 0.2e1 + t139 * t455) * t234 + t489) * t267 + t142;
t30 = Icges(4,4) * t55 + Icges(4,2) * t56;
t31 = Icges(4,1) * t55 + Icges(4,4) * t56;
t462 = t142 / 0.2e1;
t296 = t392 / 0.2e1 + t402 / 0.2e1;
t141 = t296 * t406;
t463 = -t141 / 0.2e1;
t103 = -Icges(4,1) * t141 + Icges(4,4) * t142;
t464 = t103 / 0.2e1;
t102 = -Icges(4,4) * t141 + Icges(4,2) * t142;
t465 = t102 / 0.2e1;
t394 = t260 * t265;
t370 = t268 * t256 * pkin(1);
t384 = t228 * t448 + t196;
t413 = t198 * t247;
t140 = t194 * t413 + 0.2e1 * t370 + t384;
t376 = pkin(1) * t449;
t347 = t235 * t376;
t326 = t187 * t347;
t457 = t234 / 0.2e1;
t289 = t140 * t457 - t326;
t138 = t213 + (-t390 + (0.2e1 * t247 * pkin(1) - t414) * t259) * pkin(2);
t327 = t186 * t347;
t290 = t138 * t457 - t327;
t481 = rSges(5,1) * t290 + rSges(5,2) * t289;
t72 = t481 * t394;
t387 = t262 * t265;
t73 = t481 * t387;
t362 = t260 * t406;
t133 = t496 * t362;
t134 = t296 * t362;
t99 = t133 * rSges(4,1) + t134 * rSges(4,2) - t262 * rSges(4,3);
t86 = -pkin(2) * t396 - t99;
t135 = t496 * t360;
t100 = rSges(4,1) * t135 + rSges(4,2) * t136 + rSges(4,3) * t260;
t87 = pkin(2) * t391 + t100;
t504 = t30 * t462 + t31 * t463 + t464 * t55 + t465 * t56 + m(3) * (-t211 * t225 - t212 * t226) + m(5) * (t119 * t72 - t120 * t73) + m(4) * (t20 * t86 + t21 * t87);
t105 = t290 * t265;
t181 = 0.1e1 / t186 ^ 2;
t421 = t181 * t187;
t107 = t289 * t265;
t180 = 0.1e1 / t186;
t426 = t107 * t180;
t502 = 0.2e1 * t105 * t421 - 0.2e1 * t426;
t453 = t262 / 0.2e1;
t207 = Icges(3,4) * t396 - Icges(3,2) * t400 - Icges(3,6) * t262;
t253 = Icges(3,4) * t261;
t231 = -Icges(3,2) * t259 + t253;
t208 = Icges(3,6) * t260 + t231 * t262;
t429 = Icges(3,4) * t259;
t233 = Icges(3,1) * t261 - t429;
t210 = Icges(3,5) * t260 + t233 * t262;
t200 = t210 * t396;
t229 = Icges(3,5) * t261 - Icges(3,6) * t259;
t206 = Icges(3,3) * t260 + t229 * t262;
t345 = t206 * t262 - t200;
t205 = Icges(3,5) * t396 - Icges(3,6) * t400 - Icges(3,3) * t262;
t244 = Icges(3,4) * t400;
t209 = Icges(3,1) * t396 - Icges(3,5) * t262 - t244;
t383 = -t260 * t205 - t209 * t391;
t499 = -t207 * t398 - t208 * t400 - t345 - t383;
t353 = Icges(5,4) * t459;
t304 = Icges(5,1) * t460 + t353;
t361 = t234 * t387;
t128 = Icges(5,5) * t260 - t304 * t361;
t352 = -t407 / 0.2e1;
t333 = t128 * t352;
t428 = Icges(5,4) * t186;
t303 = t428 / 0.2e1 + Icges(5,2) * t459;
t126 = Icges(5,6) * t260 - t303 * t361;
t334 = t126 * t352;
t498 = t186 * t333 + t187 * t334;
t90 = -Icges(4,5) * t133 - Icges(4,6) * t134 + Icges(4,3) * t262;
t93 = -Icges(4,4) * t133 - Icges(4,2) * t134 + Icges(4,6) * t262;
t96 = -Icges(4,1) * t133 - Icges(4,4) * t134 + Icges(4,5) * t262;
t35 = -t135 * t96 - t136 * t93 - t260 * t90;
t92 = Icges(4,5) * t135 + Icges(4,6) * t136 + Icges(4,3) * t260;
t95 = Icges(4,4) * t135 + Icges(4,2) * t136 + Icges(4,6) * t260;
t98 = Icges(4,1) * t135 + Icges(4,4) * t136 + Icges(4,5) * t260;
t36 = t135 * t98 + t136 * t95 + t260 * t92;
t178 = 0.1e1 / t185 ^ 2;
t184 = t188 ^ 2;
t153 = t178 * t184 + 0.1e1;
t150 = 0.1e1 / t153;
t447 = pkin(3) * t241;
t372 = t150 * t447;
t337 = -0.2e1 * t347;
t106 = (t137 * t234 + t185 * t337) * t267;
t422 = t178 * t188;
t108 = (t139 * t234 + t188 * t337) * t267;
t177 = 0.1e1 / t185;
t425 = t108 * t177;
t485 = -t106 * t422 + t425;
t467 = t372 * t485 + 0.1e1;
t44 = t467 * t260;
t45 = t467 * t262;
t495 = -t35 * t45 + t44 * t36;
t33 = -t133 * t96 - t134 * t93 + t262 * t90;
t183 = t187 ^ 2;
t152 = t181 * t183 + 0.1e1;
t148 = 0.1e1 / t152;
t446 = pkin(4) * t241;
t371 = t148 * t446;
t47 = t371 * t502;
t470 = t47 / 0.2e1;
t494 = -t260 / 0.2e1;
t456 = t260 / 0.2e1;
t454 = -t262 / 0.2e1;
t493 = t47 * (rSges(5,1) * t289 - rSges(5,2) * t290) * t265;
t423 = t138 * t180 * t181;
t492 = 0.1e1 / t152 ^ 2 * (t140 * t421 - t183 * t423);
t491 = Icges(5,4) * t289;
t363 = t234 * t394;
t125 = Icges(5,6) * t262 + t303 * t363;
t417 = t187 * t125;
t127 = Icges(5,5) * t262 + t304 * t363;
t418 = t186 * t127;
t490 = (t418 / 0.2e1 + t417 / 0.2e1) * t363;
t232 = Icges(3,1) * t259 + t253;
t285 = t290 * Icges(5,4);
t461 = -t186 / 0.2e1;
t487 = ((-Icges(5,1) * t290 - t491) * t459 + (-Icges(5,2) * t289 - t285) * t461) * t265;
t477 = 0.2e1 * t259;
t486 = t414 * t477 + t390;
t302 = Icges(5,5) * t460 + Icges(5,6) * t459;
t123 = Icges(5,3) * t262 + t302 * t363;
t351 = t407 / 0.2e1;
t59 = -t260 * t123 + (t417 + t418) * t262 * t351;
t432 = t262 * t59;
t124 = Icges(5,3) * t260 - t302 * t361;
t60 = t260 * t124 + t262 * t498;
t484 = t260 * t60 - t432;
t388 = t262 * t124;
t389 = t262 * t123;
t57 = t389 + t490;
t483 = -t262 * t57 + t260 * (t260 * t498 - t388);
t34 = t133 * t98 + t134 * t95 - t262 * t92;
t434 = t262 * t47;
t438 = t260 * t47;
t480 = (t138 * t333 + t140 * t334) * t438 - (t125 * t140 + t127 * t138) * t434 * t351;
t312 = -(t260 * t488 + t435) * t260 + (-t262 * t488 + t439) * t262;
t479 = 0.4e1 * t312 * t47;
t257 = t260 ^ 2;
t258 = t262 ^ 2;
t359 = t268 * t405;
t191 = t261 * t339 - 0.4e1 * t359;
t364 = 0.4e1 * t194 ^ 2 * t198 / t409;
t332 = -t364 / 0.4e1;
t282 = (-0.2e1 * t261 * t414 + (-t191 * t198 + t332) * t259) * t234;
t338 = t234 * t235 * t359;
t399 = t259 * t261;
t408 = t234 * t261;
t458 = -t234 / 0.2e1;
t10 = 0.2e1 * pkin(4) * t148 * t376 * t502 + 0.4e1 * (t426 * t492 + (-t148 * t423 - t181 * t492) * t105 * t187) * t446 + 0.2e1 * ((t105 * t140 + t107 * t138) * t181 + (-((t247 * t364 / 0.4e1 + t191 * t413 - t213 + t486 * pkin(2)) * t457 + 0.4e1 * t187 * t338 + (0.3e1 * t268 * t234 * t399 + (-0.2e1 * t140 * t259 - t187 * t261) * t450) * pkin(1)) * t180 - ((0.4e1 * t370 + t384) * t458 - 0.4e1 * t186 * t338 + (-t282 / 0.2e1 + (-t247 * t408 + (t138 * t477 + t186 * t261) * t235) * pkin(1)) * pkin(2)) * t421) * t265) * t371;
t475 = t10 / 0.2e1;
t473 = t44 / 0.2e1;
t472 = -t45 / 0.2e1;
t471 = t45 / 0.2e1;
t440 = t10 * t262;
t424 = t137 * t177 * t178;
t420 = t186 * t125;
t419 = t186 * t126;
t416 = t187 * t127;
t415 = t187 * t128;
t410 = t207 * t259;
t395 = t260 * t262;
t382 = t260 * t206 + t210 * t391;
t380 = t257 + t258;
t16 = t484 * t47;
t308 = t60 - t389;
t368 = t16 / 0.2e1 - (-t432 + (t308 + t57 - t490) * t260) * t47 / 0.2e1;
t15 = t483 * t47;
t367 = ((-t308 + t60) * t262 + (t388 + t59 + (t123 + (t126 * t187 + t128 * t186) * t351) * t260) * t260) * t470 + t15 / 0.2e1;
t358 = t438 / 0.2e1;
t357 = -t434 / 0.2e1;
t157 = -t208 * t398 + t382;
t344 = t208 * t259 - t205;
t356 = (-t260 * (-t209 * t261 + t410) - t205 * t262) * t454 + (t262 * t344 + t157 - t382) * t453 + (t260 * t344 + t345 + t499) * t456;
t355 = t157 * t456 + t382 * t494 + (-t200 + (t206 + t410) * t262 + t383 + t499) * t454;
t230 = Icges(3,2) * t261 + t429;
t348 = t233 / 0.2e1 - t230 / 0.2e1;
t329 = 0.2e1 * t106 * t188 * t447;
t325 = 0.8e1 * t338;
t104 = -rSges(4,1) * t141 + rSges(4,2) * t142;
t11 = (t178 * t329 - 0.2e1 * t425 * t447) * (t139 * t422 - t184 * t424) / t153 ^ 2 + ((-t106 * t139 - t108 * t137) * t178 + (((0.6e1 * t269 * pkin(2) * t399 - t191 * t412 + t248 * t332 - t215) * t234 + t188 * t325 + (t486 * t234 + (-0.2e1 * t392 - 0.4e1 * t403) * t450) * pkin(1)) * t177 - ((0.4e1 * t369 + t385) * t234 + t185 * t325 + (t282 + (t408 * t478 + (-0.2e1 * t393 - 0.4e1 * t404) * t235) * pkin(2)) * pkin(1)) * t422) * t267) * t372 + (0.2e1 * pkin(3) * t376 * t485 + t329 * t424) * t150;
t322 = -t104 * t11 - t448;
t314 = -Icges(3,5) * t259 - Icges(3,6) * t261;
t313 = -t119 * t262 - t120 * t260;
t310 = t265 * t327;
t309 = t265 * t326;
t307 = -t104 * t45 - t374;
t145 = (Icges(5,2) * t461 + t353) * t407;
t146 = (Icges(5,1) * t459 - t428 / 0.2e1) * t407;
t300 = t145 * t459 + t146 * t460;
t299 = t416 / 0.2e1 - t420 / 0.2e1;
t298 = t145 * t461 + t146 * t459;
t295 = t260 * t310;
t294 = t260 * t309;
t293 = t262 * t310;
t292 = t262 * t309;
t283 = t140 * t146 / 0.4e1 - t138 * t145 / 0.4e1 + (t187 * (Icges(5,1) * t289 - t285) / 0.4e1 - t186 * (-Icges(5,2) * t290 + t491) / 0.4e1) * t265;
t281 = t265 * (-Icges(5,5) * t290 - Icges(5,6) * t289);
t239 = rSges(3,1) * t261 - rSges(3,2) * t259;
t218 = t314 * t262;
t217 = t314 * t260;
t147 = (rSges(5,1) * t459 + rSges(5,2) * t461) * t407;
t144 = (Icges(5,5) * t459 + Icges(5,6) * t461) * t407;
t101 = -Icges(4,5) * t141 + Icges(4,6) * t142;
t84 = t299 * t407;
t74 = (Icges(5,5) * t289 - Icges(5,6) * t290) * t265;
t67 = t262 * t281;
t66 = t260 * t281;
t37 = t104 * t44 + t375;
t32 = rSges(4,1) * t55 + rSges(4,2) * t56;
t29 = Icges(4,5) * t55 + Icges(4,6) * t56;
t28 = Icges(4,1) * t53 + Icges(4,4) * t54;
t27 = Icges(4,1) * t51 + Icges(4,4) * t52;
t26 = Icges(4,4) * t53 + Icges(4,2) * t54;
t25 = Icges(4,4) * t51 + Icges(4,2) * t52;
t24 = Icges(4,5) * t53 + Icges(4,6) * t54;
t23 = Icges(4,5) * t51 + Icges(4,6) * t52;
t7 = t262 * t322 - t32 * t45;
t6 = t260 * t322 - t32 * t44;
t2 = (t232 / 0.2e1 + t231 / 0.2e1) * t261 + t348 * t259 + (t234 * t283 - t298 * t347) * t265 + t504;
t1 = (t368 * t47 + t355) * t262 + (t367 * t47 + t356) * t260;
t3 = [t2 * qJD(2), t2 * qJD(1) + (t495 * t471 + (t102 * t54 + t103 * t53 + t135 * t31 + t136 * t30 - t141 * t28 + t142 * t26 + t55 * t98 + t56 * t95) * t473 + (t102 * t52 + t103 * t51 + t133 * t31 + t134 * t30 - t141 * t27 + t142 * t25 - t55 * t96 - t56 * t93 + t495) * t472 + (m(3) * (t211 * t239 - t225 * t237) + t29 * t471 + t229 * t453 + (t453 * t74 - t368) * t47 + (-t209 / 0.2e1 + Icges(3,2) * t396 / 0.2e1 + t244 / 0.2e1) * t261 + (t207 / 0.2e1 - t232 * t494) * t259 + (t101 * t453 - t102 * t134 / 0.2e1 - t133 * t464 + t96 * t463 + t93 * t462) * t11 + (t84 / 0.2e1 + t144 * t453 + t300 * t260 * t351) * t10 - t355) * t262 + (m(3) * (-t212 * t239 + t226 * t237) + t29 * t473 + (t144 * t475 + t229 / 0.2e1) * t260 + (t456 * t74 - t367) * t47 + (t210 / 0.2e1 + t230 * t454) * t261 + (-t208 / 0.2e1 - t232 * t453) * t259 + (t101 * t456 + t135 * t464 + t136 * t465 + t462 * t95 + t463 * t98) * t11 - t356) * t260 + (t20 * t307 - t21 * t37 + t6 * t87 + t7 * t86) * m(4) + (t313 * t493 + ((t260 * t73 - t262 * t72) * t47 + t313 * t10) * t147) * m(5) + (((t416 - t420) * t347 + (-t140 * t127 / 0.2e1 + t138 * t125 / 0.2e1 + t487 * t260) * t234) * t357 + ((-t415 + t419) * t347 * t470 + ((t415 / 0.2e1 - t419 / 0.2e1) * t475 + (t140 * t128 / 0.2e1 - t138 * t126 / 0.2e1 + t487 * t262) * t470) * t234 + t300 * t440 * t458) * t260) * t265) * qJD(2); t1 * qJD(2) + (-t84 * t358 + (t299 * t358 - t283) * t407 + (t265 * t298 * t377 - t348) * t259 + (t232 + t231) * t455 - t504) * qJD(1), t1 * qJD(1) + (((t126 * t292 + t128 * t293 + t260 * t67) * t438 - (-t125 * t292 - t127 * t293 + t260 * t66) * t434 + t480 * t262 + t484 * t10) * t358 - t15 * t440 / 0.2e1 + ((t126 * t294 + t128 * t295 - t262 * t67) * t438 - (-t125 * t294 - t127 * t295 - t262 * t66) * t434 + t480 * t260 + t483 * t10) * t357 + ((t135 * t28 + t136 * t26 + t24 * t260 + t53 * t98 + t54 * t95) * t44 - (t135 * t27 + t136 * t25 + t23 * t260 - t53 * t96 - t54 * t93) * t45 + (t260 * t36 - t262 * t35) * t11) * t473 + ((t133 * t28 + t134 * t26 - t24 * t262 + t51 * t98 + t52 * t95) * t44 - (t133 * t27 + t134 * t25 - t23 * t262 - t51 * t96 - t52 * t93) * t45 + (t260 * t34 - t262 * t33) * t11) * t472 + m(5) * (t312 * t479 * t10 + ((-t260 * t72 - t262 * t73) * t479 + 0.4e1 * (t10 * t147 + t493) * t380 * t147) * t47) / 0.4e1 + m(3) * ((t211 * t260 + t212 * t262) * (-t225 * t260 - t226 * t262) + t380 * t239 * t237) + (t307 * t7 - t37 * t6 + (t100 * t45 + t380 * t448 + t44 * t99) * (t44 * t324 + t45 * t443 + (t100 * t262 + t260 * t99) * t11 - t380 * t449)) * m(4) + (t10 * t16 + t11 * t495 - t217 * t395 + t218 * t257) * t456 + (t11 * (-t45 * t33 + t34 * t44) + t217 * t258 - t218 * t395) * t454) * qJD(2);];
Cq = t3;
