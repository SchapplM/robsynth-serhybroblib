% Calculate time derivative of joint inertia matrix for
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
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:37
% EndTime: 2020-04-12 19:26:31
% DurationCPUTime: 26.71s
% Computational Cost: add. (271048->567), mult. (393078->1071), div. (17124->20), fcn. (108761->10), ass. (0->396)
t218 = pkin(2) ^ 2;
t219 = pkin(1) ^ 2;
t210 = cos(qJ(2));
t458 = pkin(2) * t210;
t390 = -0.2e1 * pkin(1) * t458 + t219;
t192 = t218 + t390;
t189 = 0.1e1 / t192 ^ 2;
t208 = sin(qJ(2));
t384 = qJD(2) * t208;
t368 = pkin(2) * t384;
t340 = pkin(1) * t368;
t323 = t189 * t340;
t495 = 0.2e1 * t323;
t213 = 0.1e1 / pkin(4);
t496 = pkin(3) ^ 2;
t497 = pkin(4) ^ 2;
t388 = t496 - t497;
t184 = t192 - t388;
t470 = -pkin(3) - pkin(4);
t177 = (pkin(2) - t470) * (pkin(2) + t470) + t390;
t469 = pkin(4) - pkin(3);
t178 = (pkin(2) - t469) * (pkin(2) + t469) + t390;
t473 = pkin(1) * pkin(2);
t288 = 0.2e1 * (-t177 - t178) * t473;
t157 = t208 * t288;
t150 = qJD(2) * t157;
t405 = t177 * t178;
t222 = sqrt(-t405);
t163 = 0.1e1 / t222;
t407 = t150 * t163;
t465 = -t208 / 0.2e1;
t324 = t407 * t465;
t197 = pkin(1) - t458;
t379 = 0.2e1 * t197 * pkin(1);
t395 = t210 * t222;
t112 = (t324 + (-t395 + (t184 + t379) * t208) * qJD(2)) * pkin(2);
t205 = t208 ^ 2;
t460 = pkin(1) * t218;
t370 = t205 * t460;
t332 = qJD(2) * t370;
t347 = t163 * t197 / 0.2e1;
t352 = t222 * t384;
t383 = qJD(2) * t210;
t392 = (t184 * t383 + t352) * pkin(2);
t113 = t150 * t347 + 0.2e1 * t332 + t392;
t214 = 0.1e1 / t497;
t399 = t208 * t222;
t146 = -pkin(2) * t399 + t184 * t197;
t139 = t146 ^ 2;
t459 = pkin(2) * t208;
t175 = t184 * t459;
t147 = t197 * t222 + t175;
t143 = t147 ^ 2;
t393 = t139 + t143;
t123 = t393 * t214 * t189;
t119 = t123 ^ (-0.1e1 / 0.2e1);
t188 = 0.1e1 / t192;
t190 = t188 * t189;
t282 = 0.4e1 * t190 * t340;
t378 = 0.2e1 * t189;
t367 = t119 / t123 * t188 * ((t112 * t146 + t113 * t147) * t378 - t393 * t282) * t214;
t414 = t119 * t213;
t504 = t414 * t495 + t213 * t367 / 0.2e1;
t444 = t147 * rSges(5,2);
t451 = t146 * rSges(5,1);
t320 = t444 + t451;
t363 = t188 * t414;
t261 = t320 * t363;
t503 = qJD(1) * t208;
t502 = qJD(1) * t210;
t406 = t163 * t157;
t349 = -t406 / 0.2e1;
t116 = t175 + (-t395 + (t349 + t379) * t208) * pkin(2);
t403 = t189 * t208;
t342 = t403 * t473;
t467 = -t188 / 0.2e1;
t103 = (t116 * t467 + t146 * t342) * t213;
t118 = t157 * t347 + 0.2e1 * t370 + (t184 * t210 + t399) * pkin(2);
t466 = t188 / 0.2e1;
t105 = (t118 * t466 - t147 * t342) * t213;
t140 = 0.1e1 / t146;
t141 = 0.1e1 / t146 ^ 2;
t412 = t141 * t147;
t501 = -0.2e1 * t103 * t412 - 0.2e1 * t105 * t140;
t183 = t192 + t388;
t462 = pkin(1) * t183;
t176 = t208 * t462;
t198 = pkin(1) * t210 - pkin(2);
t148 = -t198 * t222 + t176;
t397 = t210 * t148;
t145 = -pkin(1) * t399 - t183 * t198;
t401 = t208 * t145;
t299 = t397 + t401;
t272 = t299 * t188;
t217 = 0.1e1 / t496;
t135 = t145 ^ 2;
t144 = t148 ^ 2;
t394 = t135 + t144;
t124 = t394 * t217 * t189;
t121 = t124 ^ (-0.1e1 / 0.2e1);
t216 = 0.1e1 / pkin(3);
t489 = t121 * t216;
t500 = t272 * t489;
t479 = -0.2e1 * t198;
t380 = pkin(2) * t479;
t115 = t176 + (-t395 + (t349 + t380) * t208) * pkin(1);
t329 = -0.2e1 * t342;
t104 = (t115 * t188 + t145 * t329) * t216;
t346 = -t163 * t198 / 0.2e1;
t402 = t205 * t219;
t117 = t157 * t346 + 0.2e1 * pkin(2) * t402 + (t183 * t210 + t399) * pkin(1);
t106 = (t117 * t188 + t148 * t329) * t216;
t136 = 0.1e1 / t145;
t137 = 0.1e1 / t145 ^ 2;
t413 = t137 * t148;
t499 = -t104 * t413 + t106 * t136;
t494 = pkin(1) - t261;
t209 = sin(qJ(1));
t211 = cos(qJ(1));
t424 = Icges(5,2) * t146;
t428 = Icges(5,4) * t147;
t100 = (-t424 + t428) * t363;
t409 = t147 * t100;
t429 = Icges(5,4) * t146;
t432 = Icges(5,1) * t147;
t101 = (-t429 + t432) * t363;
t410 = t146 * t101;
t303 = t409 + t410;
t257 = t303 * t363;
t433 = Icges(5,1) * t146;
t308 = t428 + t433;
t260 = t308 * t363;
t86 = -Icges(5,5) * t211 - t209 * t260;
t441 = t147 * t86;
t423 = Icges(5,2) * t147;
t306 = t423 + t429;
t259 = t306 * t363;
t84 = -Icges(5,6) * t211 - t209 * t259;
t449 = t146 * t84;
t318 = -t441 + t449;
t421 = Icges(5,6) * t146;
t425 = Icges(5,5) * t147;
t99 = (-t421 + t425) * t363;
t493 = -t209 * t257 - t211 * t99 - t318 * t363;
t87 = Icges(5,5) * t209 - t211 * t260;
t440 = t147 * t87;
t85 = Icges(5,6) * t209 - t211 * t259;
t448 = t146 * t85;
t316 = -t440 + t448;
t492 = t209 * t99 - t211 * t257 - t316 * t363;
t129 = t141 * t143 + 0.1e1;
t415 = t112 * t140 / t139;
t491 = 0.2e1 / t129 ^ 2 * (t113 * t412 - t143 * t415);
t354 = qJD(2) * t402;
t331 = pkin(2) * t354;
t391 = pkin(1) * t352 + t383 * t462;
t114 = t150 * t346 + 0.2e1 * t331 + t391;
t130 = t137 * t144 + 0.1e1;
t111 = (t324 + (-t395 + (t183 + t380) * t208) * qJD(2)) * pkin(1);
t416 = t111 * t136 / t135;
t490 = 0.2e1 / t130 ^ 2 * (t114 * t413 - t144 * t416);
t430 = Icges(3,4) * t210;
t307 = -Icges(3,2) * t208 + t430;
t168 = Icges(3,6) * t209 + t211 * t307;
t431 = Icges(3,4) * t208;
t309 = Icges(3,1) * t210 - t431;
t170 = Icges(3,5) * t209 + t211 * t309;
t293 = t168 * t208 - t170 * t210;
t487 = t209 * t293;
t167 = -Icges(3,6) * t211 + t209 * t307;
t169 = -Icges(3,5) * t211 + t209 * t309;
t295 = t167 * t208 - t169 * t210;
t486 = t211 * t295;
t366 = t121 / t124 * t188 * ((t111 * t145 + t114 * t148) * t378 - t394 * t282) * t217;
t248 = (-t397 / 0.2e1 - t401 / 0.2e1) * t366;
t485 = (t451 / 0.2e1 + t444 / 0.2e1) * t209 * t367;
t345 = pkin(2) * t378;
t330 = pkin(1) * t345;
t312 = qJD(2) * t330;
t283 = t208 * t312;
t484 = t320 * t209 * t283;
t385 = qJD(1) * t211;
t483 = -t208 * t385 - t209 * t383;
t305 = Icges(3,5) * t210 - Icges(3,6) * t208;
t165 = -Icges(3,3) * t211 + t209 * t305;
t481 = qJD(1) * t165;
t186 = Icges(3,2) * t210 + t431;
t187 = Icges(3,1) * t208 + t430;
t290 = t186 * t208 - t187 * t210;
t480 = qJD(1) * t290 + t305 * qJD(2);
t478 = -0.2e1 * t209;
t477 = m(3) / 0.2e1;
t476 = m(5) / 0.2e1;
t127 = 0.1e1 / t130;
t457 = pkin(3) * t192;
t376 = t127 * t457;
t468 = t376 * t499 + 0.1e1;
t461 = pkin(1) * t189;
t456 = pkin(4) * t192;
t455 = rSges(3,3) * t209;
t454 = rSges(5,3) * t211;
t125 = 0.1e1 / t129;
t149 = (t210 * t288 - 0.8e1 * t218 * t402) * qJD(2);
t362 = t150 / t405 * t406;
t325 = -t362 / 0.4e1;
t232 = (t208 * t325 + (t149 * t465 - t150 * t210) * t163) * t188;
t310 = t190 * t218 * t354;
t335 = 0.6e1 * t208 * t383;
t371 = t141 * t456;
t336 = t147 * t371;
t337 = t125 * t371;
t348 = t406 / 0.2e1;
t361 = t188 * t407;
t372 = t140 * t456;
t373 = t125 * t456;
t404 = t188 * t210;
t13 = 0.2e1 * (t112 * t337 + t372 * t491) * t105 + 0.2e1 * (0.2e1 * t147 * t373 * t415 - t113 * t337 + t336 * t491) * t103 + 0.2e1 * (pkin(4) * t340 * t501 + (-((t197 * t362 / 0.4e1 + t149 * t347 + t335 * t460) * t466 + 0.4e1 * t147 * t310 + ((t361 / 0.4e1 - t113 * t461) * t208 + ((t395 + (t348 - t184) * t208) * t466 + (-t118 * t208 - t147 * t210) * t461) * qJD(2)) * pkin(2)) * t372 - ((0.4e1 * t332 + t392) * t467 - 0.4e1 * t146 * t310 + (-t232 / 0.2e1 + (t112 * t403 + (-t197 * t404 + (t116 * t208 + t146 * t210) * t189) * qJD(2)) * pkin(1)) * pkin(2)) * t336) * t213) * t125;
t453 = t13 * t209;
t452 = t13 * t211;
t450 = t146 * rSges(5,2);
t447 = t146 * t86;
t446 = t146 * t87;
t445 = t147 * rSges(5,1);
t443 = t147 * t84;
t442 = t147 * t85;
t58 = t373 * t501;
t439 = t188 * t58;
t202 = t209 * rSges(5,3);
t438 = t209 * t58;
t420 = Icges(5,6) * t147;
t426 = Icges(5,5) * t146;
t304 = t420 + t426;
t258 = t304 * t363;
t83 = Icges(5,3) * t209 - t211 * t258;
t437 = t209 * t83;
t436 = t211 * rSges(3,3);
t435 = t211 * t58;
t82 = -Icges(5,3) * t211 - t209 * t258;
t434 = t211 * t82;
t411 = t145 * t210;
t408 = t148 * t208;
t400 = t208 * t209;
t398 = t209 * t210;
t396 = t210 * t211;
t389 = t209 ^ 2 + t211 ^ 2;
t166 = Icges(3,3) * t209 + t211 * t305;
t387 = qJD(1) * t166;
t386 = qJD(1) * t209;
t382 = qJD(2) * t211;
t377 = 0.2e1 * qJD(1);
t375 = t136 * t457;
t374 = t137 * t457;
t302 = -t111 * t210 + t114 * t208;
t235 = qJD(2) * t299 + t302;
t243 = (-t148 * t205 + t210 * t401) * t330;
t241 = qJD(2) * t243;
t247 = (-t408 / 0.2e1 + t411 / 0.2e1) * t366;
t300 = -t408 + t411;
t23 = (t211 * t247 + (t211 * t241 + (t211 * t235 + t300 * t386) * t188) * t121) * t216;
t301 = -t111 * t208 - t114 * t210;
t236 = qJD(2) * t300 - t301;
t281 = t145 * t205 + t208 * t397;
t242 = t281 * t312;
t24 = (t211 * t248 + (-t211 * t242 + (t211 * t236 - t299 * t386) * t188) * t121) * t216;
t369 = t23 * rSges(4,1) + t24 * rSges(4,2) + rSges(4,3) * t385;
t273 = t300 * t188;
t97 = t273 * t489;
t94 = t211 * t97;
t95 = t211 * t500;
t69 = -t94 * rSges(4,1) + t95 * rSges(4,2) + t209 * rSges(4,3);
t358 = t438 / 0.2e1;
t357 = -t435 / 0.2e1;
t356 = qJD(1) * t468;
t353 = t209 * t384;
t350 = t429 / 0.2e1;
t341 = 2 * m(4);
t339 = t127 * t374;
t338 = t148 * t374;
t328 = 0.4e1 * t389;
t322 = rSges(3,1) * t210 - rSges(3,2) * t208;
t191 = rSges(3,1) * t208 + rSges(3,2) * t210;
t321 = -rSges(5,1) * t112 - rSges(5,2) * t113;
t319 = t443 + t447;
t317 = t442 + t446;
t315 = -t434 + t437;
t90 = -t209 * t261 - t454;
t91 = -t211 * t261 + t202;
t314 = t209 * t90 + t211 * t91;
t185 = Icges(3,5) * t208 + Icges(3,6) * t210;
t296 = t167 * t210 + t169 * t208;
t294 = t168 * t210 + t170 * t208;
t171 = rSges(3,1) * t398 - rSges(3,2) * t400 - t436;
t172 = t211 * t322 + t455;
t292 = t171 * t211 - t172 * t209;
t291 = t171 * t209 + t172 * t211;
t287 = 0.8e1 * t310;
t279 = t112 * t209 + t146 * t385;
t278 = -t112 * t211 + t146 * t386;
t277 = t113 * t209 + t147 * t385;
t276 = -t113 * t211 + t147 * t386;
t275 = t209 * t317;
t269 = qJD(2) * t185;
t80 = -t209 * t494 + t454;
t81 = t211 * t494 + t202;
t268 = -0.2e1 * t209 * t81 - 0.2e1 * t211 * t80;
t92 = t209 * t97;
t93 = t209 * t500;
t68 = -rSges(4,1) * t92 + rSges(4,2) * t93 - rSges(4,3) * t211;
t25 = (t209 * t247 + (t209 * t241 + (t209 * t235 - t300 * t385) * t188) * t121) * t216;
t26 = (t209 * t248 + (-t209 * t242 + (t209 * t236 + t299 * t385) * t188) * t121) * t216;
t264 = rSges(4,1) * t25 + rSges(4,2) * t26 + rSges(4,3) * t386;
t263 = t319 * t363;
t262 = t317 * t363;
t255 = -t112 * t188 + t146 * t495;
t254 = t113 * t188 - 0.2e1 * t147 * t323;
t253 = t211 * t263;
t252 = (t433 / 0.2e1 + t428 / 0.2e1) * t367;
t251 = (t350 + t423 / 0.2e1) * t367;
t250 = (t426 / 0.2e1 + t420 / 0.2e1) * t367;
t249 = (t410 / 0.2e1 + t409 / 0.2e1) * t367;
t239 = t306 * t283;
t35 = Icges(5,6) * t385 + (t211 * t251 + (t211 * t239 + (Icges(5,4) * t278 + Icges(5,2) * t276) * t188) * t119) * t213;
t240 = t308 * t283;
t37 = Icges(5,5) * t385 + (t211 * t252 + (t211 * t240 + (Icges(5,1) * t278 + Icges(5,4) * t276) * t188) * t119) * t213;
t246 = -t112 * t87 - t113 * t85 - t146 * t37 - t147 * t35;
t36 = Icges(5,6) * t386 + (t209 * t251 + (t209 * t239 + (-Icges(5,4) * t279 - Icges(5,2) * t277) * t188) * t119) * t213;
t38 = Icges(5,5) * t386 + (t209 * t252 + (t209 * t240 + (-Icges(5,1) * t279 - Icges(5,4) * t277) * t188) * t119) * t213;
t245 = t112 * t86 + t113 * t84 + t146 * t38 + t147 * t36;
t42 = ((-t428 / 0.2e1 + t424 / 0.2e1) * t367 + (Icges(5,4) * t254 + Icges(5,2) * t255) * t119) * t213;
t43 = ((-t432 / 0.2e1 + t350) * t367 + (Icges(5,1) * t254 + Icges(5,4) * t255) * t119) * t213;
t244 = -t100 * t113 - t101 * t112 - t146 * t43 - t147 * t42;
t238 = t304 * t283;
t237 = t303 * t283;
t233 = t386 * t261 + rSges(5,3) * t385 + (t320 * t504 + t321 * t363) * t211;
t231 = ((t446 / 0.2e1 + t442 / 0.2e1) * t209 - (t447 / 0.2e1 + t443 / 0.2e1) * t211) * t58 * t367;
t230 = (-t211 * t319 + t275) * t58 * t495;
t182 = t322 * qJD(2);
t159 = -rSges(3,1) * t353 + (rSges(3,1) * t396 + t455) * qJD(1) + t483 * rSges(3,2);
t158 = -t191 * t382 + (-t209 * t322 + t436) * qJD(1);
t152 = -t209 * t269 + t387;
t151 = -t211 * t269 - t481;
t134 = t166 * t209 - t211 * t293;
t133 = t165 * t209 - t486;
t132 = -t166 * t211 - t487;
t131 = -t165 * t211 - t209 * t295;
t102 = (t445 - t450) * t363;
t73 = -rSges(4,1) * t500 - rSges(4,2) * t97;
t72 = -Icges(4,1) * t500 - Icges(4,4) * t97;
t71 = -Icges(4,4) * t500 - Icges(4,2) * t97;
t70 = -Icges(4,5) * t500 - Icges(4,6) * t97;
t67 = -Icges(4,1) * t94 + Icges(4,4) * t95 + Icges(4,5) * t209;
t66 = -Icges(4,1) * t92 + Icges(4,4) * t93 - Icges(4,5) * t211;
t65 = -Icges(4,4) * t94 + Icges(4,2) * t95 + Icges(4,6) * t209;
t64 = -Icges(4,4) * t92 + Icges(4,2) * t93 - Icges(4,6) * t211;
t63 = -Icges(4,5) * t94 + Icges(4,6) * t95 + Icges(4,3) * t209;
t62 = -Icges(4,5) * t92 + Icges(4,6) * t93 - Icges(4,3) * t211;
t61 = pkin(2) * t396 + t69;
t60 = -pkin(2) * t398 - t68;
t56 = t468 * t211;
t55 = t468 * t209;
t50 = -t211 * t262 + t437;
t49 = t209 * t82 - t253;
t48 = -t209 * t262 - t211 * t83;
t47 = -t209 * t263 - t434;
t46 = -t211 * t459 - t56 * t73;
t45 = -pkin(2) * t400 - t55 * t73;
t44 = ((-t445 / 0.2e1 + t450 / 0.2e1) * t367 + (rSges(5,1) * t254 + rSges(5,2) * t255) * t119) * t213;
t41 = ((-t425 / 0.2e1 + t421 / 0.2e1) * t367 + (Icges(5,5) * t254 + Icges(5,6) * t255) * t119) * t213;
t40 = (t247 + (t302 * t188 + (t272 + t243) * qJD(2)) * t121) * t216;
t39 = (-t248 + (t301 * t188 + (t281 * t330 - t273) * qJD(2)) * t121) * t216;
t34 = Icges(5,3) * t386 + (t209 * t250 + (t209 * t238 + (-Icges(5,5) * t279 - Icges(5,6) * t277) * t188) * t119) * t213;
t33 = Icges(5,3) * t385 + (t211 * t250 + (t211 * t238 + (Icges(5,5) * t278 + Icges(5,6) * t276) * t188) * t119) * t213;
t32 = (-pkin(1) * t211 - t202) * qJD(1) + (-t485 + (-t484 + (rSges(5,1) * t279 + rSges(5,2) * t277) * t188) * t119) * t213;
t31 = -pkin(1) * t386 + t233;
t30 = t209 * t63 + t65 * t95 - t67 * t94;
t29 = t209 * t62 + t64 * t95 - t66 * t94;
t28 = -t211 * t63 + t65 * t93 - t67 * t92;
t27 = -t211 * t62 + t64 * t93 - t66 * t92;
t18 = rSges(4,1) * t39 + rSges(4,2) * t40;
t17 = Icges(4,1) * t39 + Icges(4,4) * t40;
t16 = Icges(4,4) * t39 + Icges(4,2) * t40;
t15 = Icges(4,5) * t39 + Icges(4,6) * t40;
t14 = (-t111 * t339 - t375 * t490) * t106 + (0.2e1 * t148 * t376 * t416 - t114 * t339 + t338 * t490) * t104 + ((((t219 * pkin(2) * t335 + t149 * t346 + t198 * t325) * t188 + t148 * t287 + ((t361 / 0.2e1 - 0.2e1 * t114 * t189 * pkin(2)) * t208 + ((t395 + (t348 - t183) * t208) * t188 + (-t117 * t208 - t397) * t345) * qJD(2)) * pkin(1)) * t375 - ((0.4e1 * t331 + t391) * t188 + t145 * t287 + (t232 + (-0.2e1 * t111 * t403 + (t404 * t479 + (-t115 * t208 - t411) * t378) * qJD(2)) * pkin(2)) * pkin(1)) * t338) * t216 + 0.2e1 * t499 * pkin(3) * t340) * t127;
t12 = Icges(4,1) * t25 + Icges(4,4) * t26 + Icges(4,5) * t386;
t11 = Icges(4,1) * t23 + Icges(4,4) * t24 + Icges(4,5) * t385;
t10 = Icges(4,4) * t25 + Icges(4,2) * t26 + Icges(4,6) * t386;
t9 = Icges(4,4) * t23 + Icges(4,2) * t24 + Icges(4,6) * t385;
t8 = Icges(4,5) * t25 + Icges(4,6) * t26 + Icges(4,3) * t386;
t7 = Icges(4,5) * t23 + Icges(4,6) * t24 + Icges(4,3) * t385;
t6 = (-t210 * t385 + t353) * pkin(2) - t264;
t5 = (-t208 * t382 - t210 * t386) * pkin(2) + t369;
t4 = t14 * t209 + t211 * t356;
t3 = -t14 * t211 + t209 * t356;
t2 = pkin(2) * t483 - t18 * t55 - t4 * t73;
t1 = -t18 * t56 + t3 * t73 + (t208 * t386 - t210 * t382) * pkin(2);
t19 = [t39 * t72 - t500 * t17 + t40 * t71 - t97 * t16 + (t5 * t61 + t6 * t60) * t341 + 0.2e1 * m(5) * (t31 * t81 + t32 * t80) + 0.2e1 * m(3) * (t158 * t172 + t159 * t171) + (-t100 * t112 + t101 * t113 - t146 * t42 + t147 * t43) * t363 + (t309 - t186) * t384 + (t187 + t307) * t383 + (t146 * t100 - t147 * t101) * t504; (t58 * t44 * t268 + (t13 * t268 + (t31 * t478 - 0.2e1 * t211 * t32 + (t209 * t80 - t211 * t81) * t377) * t58) * t102) * t476 + (t209 * t41 + t99 * t385) * t358 + (-t211 * t41 + t99 * t386) * t357 + m(4) * (t1 * t60 + t2 * t61 + t45 * t5 + t46 * t6) + (-t211 * t70 - t500 * t66 - t64 * t97 + t71 * t93 - t72 * t92) * t3 / 0.2e1 + (t209 * t70 - t500 * t67 - t65 * t97 + t71 * t95 - t72 * t94) * t4 / 0.2e1 + (-t11 * t500 + t15 * t209 + t16 * t95 - t17 * t94 + t23 * t72 + t24 * t71 + t385 * t70 + t39 * t67 + t40 * t65 - t9 * t97) * t55 / 0.2e1 - (-t10 * t97 - t12 * t500 - t15 * t211 + t16 * t93 - t17 * t92 + t25 * t72 + t26 * t71 + t386 * t70 + t39 * t66 + t40 * t64) * t56 / 0.2e1 + (-qJD(2) * t293 - t167 * t502 - t169 * t503 + t209 * t480) * t209 / 0.2e1 - (-qJD(2) * t295 + t168 * t502 + t170 * t503 - t211 * t480) * t211 / 0.2e1 + t492 * t453 / 0.2e1 - t493 * t452 / 0.2e1 + (((-t440 / 0.2e1 + t448 / 0.2e1) * t367 + t211 * t249 + (t211 * t237 + t316 * t283 + (-t112 * t85 + t113 * t87 - t146 * t35 + t147 * t37 + t211 * t244 + t303 * t386) * t188) * t119) * t358 + ((-t441 / 0.2e1 + t449 / 0.2e1) * t367 + t209 * t249 + (t209 * t237 + t318 * t283 + (-t112 * t84 + t113 * t86 - t146 * t36 + t147 * t38 + t209 * t244 - t303 * t385) * t188) * t119) * t357) * t213 + ((t158 * t478 + 0.2e1 * t159 * t211 - t291 * t377) * t191 + 0.2e1 * t292 * t182) * t477 + (-t185 * t211 - t290 * t209 + t493 * t58 + t296) * t386 / 0.2e1 + (t209 * t185 - t290 * t211 + t492 * t58 + t294) * t385 / 0.2e1; -((-t47 * t13 + (t211 * t34 + (t48 + t253) * qJD(1)) * t58) * t211 + (t48 * t13 + (-t211 * t33 + (t315 + t47) * qJD(1)) * t58 + (t231 + (t230 + (t246 * t209 + (-qJD(1) * t317 + t245) * t211) * t439) * t119) * t213) * t209) * t435 + ((t50 * t13 + (t209 * t33 + (t275 * t363 + t49) * qJD(1)) * t58) * t209 + (-t49 * t13 + (-t209 * t34 + (t315 + t50) * qJD(1)) * t58 + (t231 + (t230 + (t245 * t211 + (-qJD(1) * t319 + t246) * t209) * t439) * t119) * t213) * t211) * t438 + (-t131 * t211 + t132 * t209) * t386 - t211 * ((t152 * t211 + (t132 + t486) * qJD(1)) * t211 + (t131 * qJD(1) + (-t168 * t383 - t170 * t384 + t387) * t209 + (-t151 + t296 * qJD(2) + (-t165 - t293) * qJD(1)) * t211) * t209) + (t46 * t1 + t45 * t2 + (t389 * t458 + t55 * t68 + t56 * t69) * (t264 * t55 - t3 * t69 - t368 * t389 + t369 * t56 + t4 * t68)) * t341 + (-t133 * t211 + t134 * t209) * t385 + t209 * ((t151 * t209 + (t133 + t487) * qJD(1)) * t209 + (t134 * qJD(1) + (t167 * t383 + t169 * t384 - t481) * t211 + (-t152 - t294 * qJD(2) + (t166 - t295) * qJD(1)) * t209) * t211) + t4 * (-t29 * t56 + t30 * t55) + t55 * ((-t11 * t94 + t209 * t7 + t23 * t67 + t24 * t65 + t385 * t63 + t9 * t95) * t55 + t30 * t4 - (t10 * t95 - t12 * t94 + t209 * t8 + t23 * t66 + t24 * t64 + t385 * t62) * t56 + t29 * t3) + ((t102 ^ 2 * t328 + 0.4e1 * t314 ^ 2) * t13 + (t44 * t102 * t328 + 0.4e1 * t314 * ((qJD(1) * t90 + t233) * t211 + ((-t91 + t202) * qJD(1) + (t485 + (t484 + (t209 * t321 - t320 * t385) * t188) * t119) * t213) * t209)) * t58) * t58 * t476 + t3 * (-t27 * t56 + t28 * t55) - t56 * ((-t11 * t92 - t211 * t7 + t25 * t67 + t26 * t65 + t386 * t63 + t9 * t93) * t55 + t28 * t4 - (t10 * t93 - t12 * t92 - t211 * t8 + t25 * t66 + t26 * t64 + t386 * t62) * t56 + t27 * t3) + (0.4e1 * t291 * (qJD(1) * t292 + t158 * t211 + t159 * t209) + t191 * t182 * t328) * t477 + (t385 * t58 + t453) * (t209 * t50 - t211 * t49) * t58 + (t386 * t58 - t452) * (t209 * t48 - t211 * t47) * t58;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t19(1), t19(2); t19(2), t19(3);];
Mq = res;
