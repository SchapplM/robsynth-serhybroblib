% Calculate kinetic energy for
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1TE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1TE_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:13:19
% EndTime: 2020-04-27 09:13:40
% DurationCPUTime: 19.04s
% Computational Cost: add. (277115->695), mult. (834941->1080), div. (1988->13), fcn. (104258->6), ass. (0->447)
t274 = pkin(5) ^ 2;
t273 = t274 ^ 2;
t277 = pkin(4) ^ 2;
t276 = t277 ^ 2;
t441 = t273 - t276;
t288 = (pkin(1) ^ 2);
t272 = -2 * t288;
t240 = sin(qJ(2));
t243 = cos(qJ(1));
t432 = qJD(1) * t243;
t180 = pkin(2) * t432;
t361 = pkin(3) * t180;
t337 = t240 * t361;
t241 = sin(qJ(1));
t242 = cos(qJ(2));
t428 = qJD(2) * t242;
t370 = t241 * t428;
t503 = pkin(2) * pkin(3);
t338 = t370 * t503;
t532 = -0.6e1 * t338 - 0.6e1 * t337;
t531 = t338 + t337;
t195 = pkin(3) * t242;
t174 = pkin(1) + t195;
t433 = qJD(1) * t241;
t374 = t174 * t433;
t430 = qJD(2) * t240;
t401 = pkin(3) * t430;
t489 = pkin(2) * t243;
t519 = pkin(2) * t374 + t401 * t489;
t530 = 0.2e1 * pkin(1);
t424 = 0.4e1 * pkin(3);
t282 = pkin(3) ^ 2;
t286 = pkin(2) ^ 2;
t294 = t286 ^ 2;
t221 = -t274 / 0.3e1;
t521 = t221 - t277 / 0.3e1;
t382 = t288 + t521;
t349 = t282 + t382;
t229 = t277 / 0.3e1;
t271 = 2 * t288;
t384 = t274 / 0.3e1 + t229 + t271;
t527 = 0.10e2 / 0.3e1 * t294 - 0.2e1 * (-t282 + t384) * t286 + t349 * t272;
t526 = 0.2e1 * t180;
t196 = pkin(2) * t241;
t525 = 0.2e1 * t196;
t524 = 0.4e1 * t282;
t261 = 0.6e1 * t282;
t251 = 0.10e2 * t282;
t188 = V_base(6) + qJD(2);
t523 = t188 / 0.2e1;
t488 = pkin(2) * t282;
t205 = t242 ^ 2;
t486 = pkin(3) * t205;
t200 = -0.3e1 * t282 + t288;
t457 = t282 * t205;
t414 = 0.4e1 * t457;
t522 = t200 + t414;
t270 = 3 * t288;
t520 = t270 - t274 - t277;
t369 = t240 * t432;
t313 = t369 + t370;
t518 = qJD(1) - qJD(2);
t517 = -t273 / 0.6e1 + t276 / 0.6e1;
t515 = 0.8e1 * pkin(1);
t514 = -0.4e1 * pkin(3);
t513 = -0.2e1 * t205;
t512 = -0.4e1 * t240;
t511 = -0.2e1 * t240;
t510 = -0.2e1 * t242;
t509 = -0.2e1 * t243;
t255 = -0.6e1 * t274;
t280 = t282 ^ 2;
t508 = 0.4e1 * t280;
t507 = -0.8e1 * t282;
t191 = t288 - t286 / 0.3e1;
t506 = 0.24e2 * t191;
t505 = pkin(1) * pkin(2);
t504 = pkin(1) * pkin(3);
t173 = pkin(1) - t489;
t145 = t173 + t195;
t434 = t288 - t274;
t375 = t286 + t434;
t348 = t282 + t375;
t160 = -t277 + t348;
t422 = pkin(1) * t489;
t176 = -0.2e1 * t422;
t120 = t176 + t160;
t209 = t286 * t524;
t456 = t286 * t288;
t177 = t209 - 0.4e1 * t456;
t208 = t243 ^ 2;
t256 = 0.2e1 * t277;
t466 = t240 * t241;
t410 = pkin(2) * t466;
t476 = t173 * t242;
t467 = t208 * t286;
t417 = 0.2e1 * t467;
t438 = -t286 + t288;
t123 = t176 + t417 + t438;
t477 = t123 * t205;
t496 = pkin(5) - pkin(4);
t497 = -pkin(4) - pkin(5);
t406 = pkin(3) * t466;
t365 = pkin(2) * t406;
t162 = -0.2e1 * t365;
t94 = t162 + t120;
t289 = sqrt(t177 * t208 + 0.4e1 * t160 * t422 - t280 - (t288 + (pkin(2) - t497) * (pkin(2) + t497)) * (t288 + (pkin(2) - t496) * (pkin(2) + t496)) + (t256 + t272 + 0.2e1 * t274 - 0.6e1 * t286 - 0.4e1 * t477) * t282 + (t120 * t410 - t476 * t94) * t424);
t421 = pkin(1) * t195;
t175 = 0.2e1 * t421;
t122 = t175 + t522;
t201 = -0.3e1 * t286 + t288;
t204 = t242 * t205;
t297 = pkin(3) * t282;
t471 = t204 * t297;
t420 = 0.8e1 * t471;
t367 = pkin(1) * t420;
t151 = t201 * t367;
t376 = t282 + t434;
t168 = t277 + t376;
t198 = t282 + t288;
t193 = t198 ^ 2;
t249 = 0.15e2 * t280;
t252 = 18 * t288;
t253 = -0.2e1 * t274;
t287 = t288 ^ 2;
t266 = 3 * t287;
t293 = pkin(2) * t286;
t283 = t293 ^ 2;
t207 = t243 * t208;
t468 = t207 * t293;
t391 = t174 * t468;
t360 = -0.8e1 * t391;
t388 = 0.6e1 * t421;
t403 = 0.12e2 * t457;
t219 = -t274 / 0.6e1;
t238 = t286 / 0.2e1;
t448 = t238 + t288;
t316 = -t365 + t448;
t108 = t277 / 0.6e1 + t219 + t316;
t437 = t286 + t288;
t377 = t282 + t437;
t144 = t221 + t229 + t377;
t190 = -t282 / 0.3e1 + t288;
t426 = 0.4e1 * pkin(1);
t461 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t235 = t282 / 0.3e1;
t98 = -0.4e1 / 0.9e1 * t365 + t288 + t286 / 0.3e1 + t235 + t277 / 0.9e1 - t274 / 0.9e1;
t44 = t190 * t162 + t144 * t461 + 0.6e1 * t98 * t457 + (t108 * t195 + t471) * t426;
t446 = 0.15e2 * t282 + t270;
t105 = t162 + t144;
t178 = t438 * t524;
t352 = -0.4e1 * t365;
t223 = -0.2e1 / 0.3e1 * t274;
t228 = 0.2e1 / 0.3e1 * t277;
t380 = t223 + t228 + t271;
t379 = t223 + t198;
t453 = (t228 + t379) * t198 + t294;
t56 = t144 * t352 + (t261 + t380) * t286 + t453;
t45 = 0.4e1 * t105 * t421 + t178 * t205 + t56;
t475 = t174 * t243;
t124 = t191 * t162;
t460 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t77 = t144 * t460 + t124;
t211 = 0.10e2 / 0.3e1 * t282;
t80 = (t211 + t380) * t286 + t453;
t30 = t122 * t360 + t151 + t77 * t403 + t56 * t388 + t283 + (-t274 + t277 + t446) * t294 + t193 * t168 + (0.12e2 * t44 * t208 + t249 + (t252 + t255 + 0.6e1 * t277) * t282 + t266 + (t253 + t256) * t288) * t286 + 0.6e1 * (-t406 * t80 - t45 * t475) * pkin(2);
t179 = pkin(2) * t508 + 0.8e1 * t282 * t293;
t389 = t297 * t460;
t101 = t179 * t241 + 0.4e1 * t240 * t389;
t260 = 0.5e1 * t280;
t435 = t287 + t294;
t445 = t251 + t271;
t454 = t288 * t282;
t121 = t445 * t286 + t260 + t435 + 0.6e1 * t454;
t264 = 0.5e1 * t294;
t268 = 6 * t288;
t131 = t264 + (t251 + t268) * t286 + t193;
t359 = 0.8e1 * t391;
t418 = -0.4e1 * t467;
t490 = pkin(1) * t242;
t262 = 0.3e1 * t282;
t183 = t262 + t437;
t152 = t183 * t196;
t487 = pkin(3) * t240;
t157 = t196 - t487;
t409 = t241 * t488;
t265 = 0.3e1 * t286;
t181 = t265 + t198;
t474 = t181 * t240;
t57 = t409 * t513 + t152 + (0.2e1 * t157 * t490 - t474) * pkin(3);
t317 = -t240 * t297 + t409;
t147 = 0.2e1 * t317;
t197 = 0.2e1 * t282 + t286;
t199 = -t282 + t288;
t347 = t199 + t175;
t60 = t147 * t205 + t196 * t347 + t197 * t487;
t97 = -pkin(3) * t474 + t152;
t35 = t60 * t418 + t101 * t205 + (-0.4e1 * t97 * t490 + (t131 + t359) * t240) * pkin(3) + (0.4e1 * t57 * t475 + (-t121 + t367) * t241) * pkin(2);
t23 = t145 * t30 + t289 * t35;
t502 = 0.1e1 / t23 / 0.4e1;
t463 = t241 * t242;
t110 = pkin(2) * t463 + t173 * t240;
t159 = t262 + t277 + t375;
t161 = -pkin(3) + t410;
t95 = t159 + t176 + t352;
t32 = -t95 * t476 + t110 * t289 + (pkin(1) * t161 * t509 + t159 * t466) * pkin(2) + (-t168 - t265 + t417 - 0.2e1 * t477) * pkin(3);
t501 = t32 / 0.2e1;
t106 = -t161 + t476;
t119 = t176 + t277 + t348;
t61 = t119 * t173 + 0.2e1 * t123 * t195;
t64 = t119 * t242 + (0.4e1 * t205 - 0.2e1) * t173 * pkin(3);
t34 = t106 * t289 + t196 * t64 + t240 * t61;
t500 = -t34 / 0.2e1;
t407 = pkin(3) * t476;
t78 = t162 + t176 + t377 + 0.2e1 * t407;
t75 = 0.1e1 / t78;
t499 = t75 / 0.2e1;
t170 = pkin(1) * V_base(6) + V_base(2);
t113 = t242 * t170 - t240 * V_base(1);
t498 = pkin(3) * t523 + t113 / 0.2e1;
t126 = -t240 * V_base(4) + t242 * V_base(5);
t495 = t126 / 0.2e1;
t129 = t240 * V_base(5) + t242 * V_base(4);
t494 = t129 / 0.2e1;
t493 = 0.4e1 / 0.3e1 * t286;
t492 = -t289 / 0.4e1;
t491 = t289 / 0.4e1;
t371 = t240 * t428;
t344 = t123 * t371;
t357 = 0.4e1 * t370;
t402 = pkin(2) * t433;
t167 = t402 * t530;
t452 = 0.2e1 * t531;
t387 = t167 - t452;
t431 = qJD(1) * t286;
t465 = t240 * t243;
t345 = t241 * t243 * t431;
t336 = -0.4e1 * t345;
t115 = t167 + t336;
t478 = t115 * t205;
t485 = ((0.8e1 * t344 - 0.4e1 * t478) * t282 + (-0.4e1 * t160 * t505 + t177 * t509) * t433 + (t240 * t241 ^ 2 * t431 * t515 + 0.4e1 * (-t242 * t387 + t430 * t94) * t173 + (t120 * t357 + 0.4e1 * (t120 * t465 - t463 * t94) * qJD(1)) * pkin(2)) * pkin(3)) / t289;
t425 = 0.2e1 * pkin(3);
t484 = ((-t173 * t430 + t242 * t402) * t425 + t387) / t78 ^ 2;
t439 = t280 + t287;
t442 = t271 - t274;
t455 = t288 * t274;
t111 = t442 * t282 + t439 - t455 - t517;
t312 = t111 + t294;
t82 = (t211 + t442) * t286 + t312;
t483 = t532 * t82;
t278 = 0.1e1 / pkin(4);
t482 = t278 * t75;
t481 = -0.4e1 * t531 * t144;
t311 = t313 * pkin(3);
t411 = pkin(2) * t457;
t306 = t191 * t311 * t411;
t470 = t205 * t297;
t413 = pkin(1) * t470;
t335 = t413 * t430;
t320 = -0.24e2 * t335;
t480 = t201 * t320 - 0.24e2 * t306;
t224 = -0.3e1 / 0.2e1 * t274;
t451 = t273 / 0.2e1 - t276 / 0.2e1;
t479 = t283 + t198 * ((t224 + t271) * t282 - 0.3e1 / 0.2e1 * t455 + t439 + t451);
t473 = t193 * (-t277 + t376);
t472 = t204 * t280;
t206 = t208 ^ 2;
t469 = t206 * t294;
t464 = t240 * t289;
t462 = t242 * t289;
t459 = 0.1e1 / pkin(5) / pkin(4) ^ 2;
t203 = t205 ^ 2;
t458 = t280 * t203;
t450 = t219 - t277 / 0.6e1;
t233 = -0.2e1 / 0.3e1 * t277;
t449 = t223 + t233;
t447 = 0.4e1 / 0.7e1 * t288 - t274 / 0.7e1;
t444 = t253 - 0.2e1 * t277;
t254 = -0.5e1 * t274;
t443 = t254 - 0.5e1 * t277;
t436 = t287 - t280;
t429 = qJD(2) * t241;
t423 = 0.4e1 * t195;
t419 = -0.2e1 * t470;
t416 = 0.8e1 * t458;
t415 = -0.6e1 * t457;
t412 = pkin(2) * t475;
t408 = t297 * t196;
t405 = 0.16e2 * t471;
t404 = -0.12e2 * t467;
t400 = pkin(3) * t428;
t399 = t75 * t459;
t398 = pkin(3) * t201 / 0.2e1;
t397 = -t487 / 0.2e1;
t132 = -t401 + t402;
t362 = pkin(1) * t401;
t166 = -0.2e1 * t362;
t202 = t240 ^ 2;
t329 = t208 * t293 * t374;
t318 = -0.24e2 * t329;
t322 = t190 * t337;
t323 = t190 * t338;
t324 = t183 * t180 - t181 * t400;
t341 = t282 * t371;
t331 = -0.24e2 * t341;
t332 = t180 * t471;
t372 = qJD(2) * t468;
t353 = pkin(3) * t372;
t334 = t240 * t353;
t340 = -0.6e1 * t362;
t342 = t428 * t460;
t358 = -0.2e1 * t371;
t366 = pkin(3) * t411;
t393 = t485 / 0.2e1;
t368 = t490 * t488;
t87 = -0.4e1 * t313 * t368;
t396 = -((0.8e1 * t242 * t174 * t353 + t202 * t372 * t507 + t318 * t487 + ((pkin(3) * t197 + t147 * t511 + t419) * t428 + (t241 * t166 + (t347 + 0.2e1 * t457) * t432) * pkin(2)) * t418 + 0.8e1 * t60 * t345 + 0.4e1 * ((t240 * t357 + t432 * t513) * t488 + ((t180 - t400) * t195 - t157 * t401) * t530 + t324) * t412 + t332 * t515 + t320 * t196 + (t179 * t432 + 0.4e1 * t297 * t342) * t205 + t101 * t358 - 0.4e1 * t324 * t421 + 0.4e1 * t97 * t362 - t121 * t180 + t131 * t400 - 0.4e1 * t519 * t57) * t289 + t35 * t393 + ((t166 - 0.8e1 * t341) * t360 + 0.12e2 * (-0.12e2 * t335 + 0.6e1 * (-0.4e1 / 0.9e1 * t370 - 0.4e1 / 0.9e1 * t369) * t366 - 0.12e2 * t98 * t341 + t87 - 0.4e1 * t108 * t362 - 0.2e1 * t323 - 0.2e1 * t322) * t467 - 0.24e2 * t44 * t345 - 0.6e1 * (t178 * t358 + (-t105 * t430 - t242 * t452) * pkin(1) * t424 + t481) * t412 + t77 * t331 + t481 * t388 + t56 * t340 + t480 + t532 * t80 + 0.6e1 * t519 * t45 + (0.8e1 * t334 + 0.24e2 * t329) * t122) * t145 + t30 * t132) / t23 ^ 2 / 0.4e1;
t395 = t32 * t502;
t394 = t34 * t502;
t392 = -t484 / 0.2e1;
t390 = t240 * t469;
t386 = 0.4e1 / 0.3e1 * t274 + 0.4e1 / 0.3e1 * t277 + t272;
t269 = 4 * t288;
t385 = 0.2e1 / 0.3e1 * t274 + t228 + t269;
t383 = t288 + t450;
t222 = -t274 / 0.2e1;
t381 = t222 - t277 / 0.2e1 + t288;
t378 = t268 + t444;
t373 = t204 * t430;
t189 = V_base(6) + qJD(1);
t363 = pkin(1) * t405;
t117 = 0.4e1 / 0.3e1 * t457 + t175 + t190;
t356 = -0.24e2 * t117 * t469;
t355 = -0.32e2 * t391;
t351 = -0.3e1 * t455 + t266 + t451;
t350 = t238 + t383;
t73 = t294 + (t445 + t449) * t286 + t260 + t378 * t282 + t288 * (t288 + t449);
t89 = t264 + (t251 + t378) * t286 + (t233 + t379) * t198;
t48 = t73 * t196 - t487 * t89;
t346 = t207 * t294 * t433;
t343 = t280 * t373;
t135 = -t241 * V_base(1) + t243 * V_base(2);
t43 = -0.6e1 * t82 * t365 + (t249 + (t252 - 0.9e1 * t274) * t282 + t351) * t286 + (t224 + t446) * t294 + t479;
t333 = t180 * t458;
t164 = t222 + t377;
t330 = t164 * t352;
t171 = -V_base(5) * pkin(1) + V_base(3);
t326 = -t463 + t465;
t128 = t242 * t243 + t466;
t325 = t73 * t180 - t400 * t89;
t321 = -0.48e2 * t335;
t319 = (6 * t287) + t441 - 0.6e1 * t455;
t153 = 0.16e2 * (t435 - 0.6e1 * t456) * t280;
t194 = -0.30e2 * t274 + (60 * t288);
t290 = pkin(1) * t288;
t220 = -t274 / 0.4e1;
t236 = t282 / 0.2e1;
t116 = -0.2e1 / 0.3e1 * t365 + t288 + t236 + t220;
t184 = (t269 + t274) * t282;
t192 = t288 - 0.2e1 / 0.3e1 * t286;
t155 = t288 + t286 / 0.4e1 + t282 / 0.4e1 - t274 / 0.8e1;
t53 = -0.32e2 / 0.21e2 * t155 * t365 + t294 / 0.7e1 + (0.16e2 / 0.21e2 * t282 + t447) * t286 + t280 / 0.7e1 + t447 * t282 + t287 - 0.3e1 / 0.7e1 * t455 + t441 / 0.42e2;
t156 = t220 + t235 + t448;
t234 = 0.4e1 / 0.3e1 * t282;
t55 = -0.8e1 / 0.3e1 * t156 * t365 + t294 / 0.3e1 + (t234 + t221) * t286 + t287 - t280 / 0.3e1 + (t493 + 0.2e1 / 0.3e1 * t282 + t223) * t288 + t441 / 0.18e2;
t36 = t192 * t416 + 0.14e2 * t53 * t457 + t190 * t330 + t199 * t294 + (t184 - 0.10e2 / 0.3e1 * t280 + (2 * t287) - t455) * t286 + t111 * t461 + (t116 * t405 + 0.6e1 * t195 * t55) * pkin(1);
t59 = t330 + (t261 + t442) * t286 + t312;
t79 = t164 * t460 + t124;
t37 = t388 * t59 + t403 * t79 + t151 + t43;
t267 = 8 * t288;
t104 = t408 * t512 + t209 + t508 + (t253 + t267) * t282;
t112 = t220 - t282 + t316;
t46 = t162 * t461 + t104 * t205 + t164 * t200 + (t112 * t423 + t420) * pkin(1);
t47 = t191 * t330 - t283 + (-t211 - t434) * t294 + (t184 + t436 + t517) * t286 + t111 * t288;
t259 = 0.7e1 * t280;
t52 = (t224 + t270 + 0.7e1 * t282) * t294 + (t259 + (t254 + (10 * t288)) * t282 + t351) * t286 + t479;
t172 = -0.12e2 * pkin(1) * t297 + t290 * t424;
t187 = -0.8e1 * t280 + 0.12e2 * t454;
t68 = t172 * t242 + t187 * t205 + t363 + t416 + t439 - 0.6e1 * t454;
t83 = t162 * t460 + t164 * t201;
t18 = t46 * t355 + t153 * t203 + 0.24e2 * t47 * t457 + (t253 + t269 + 0.28e2 * t282) * t283 + t168 * t473 + (t194 * t280 + 0.24e2 * t36 * t208 + t255 * t287 + t319 * t261 + t441 * t271 + 0.4e1 * t290 ^ 2 + 0.28e2 * t297 ^ 2) * t286 + 0.8e1 * (-t37 * t475 - t406 * t52) * pkin(2) + (0.8e1 * t195 * t43 + 0.32e2 * t471 * t83) * pkin(1) + (t194 * t282 + 0.16e2 * t206 * t68 + 0.70e2 * t280 + t294 + t319) * t294;
t133 = t493 + t236 + t383;
t134 = t234 + t350;
t70 = -t133 * t487 + t134 * t196;
t139 = t282 + t350;
t182 = 0.2e1 * t286 + t199;
t81 = t139 * t196 + t182 * t397;
t84 = t286 * t199 - 0.5e1 / 0.3e1 * t280 + t384 * t282 + t288 * t382;
t212 = -0.20e2 / 0.3e1 * t282;
t85 = -t294 + (t212 + t385) * t286 - 0.3e1 * t280 + t386 * t282 + t287;
t38 = t70 * t414 + t85 * t397 + (-0.8e1 / 0.3e1 * t458 + t84) * t196 + (t195 * t81 - t240 * t472) * t426;
t142 = t265 + t349;
t143 = t183 + t521;
t74 = -t142 * t487 + t143 * t196;
t154 = t236 + t286 + t450;
t90 = t154 * t525 + t460 * t487;
t39 = -0.4e1 * t90 * t457 + (-0.8e1 * t204 * t408 + t423 * t74) * pkin(1) + t48;
t86 = -0.3e1 * t294 + (t212 + t386) * t286 + t385 * t282 + t436;
t49 = t196 * t86 + t487 * t527;
t138 = 0.3e1 / 0.2e1 * t286 + t262 + t381;
t146 = 0.4e1 * t317;
t158 = t196 + 0.2e1 * t487;
t50 = t200 * t196 + t146 * t205 + (t138 * t240 + t158 * t490) * t425;
t165 = t520 * t251;
t51 = t283 + (0.21e2 * t282 + t520) * t294 + (t288 * t444 + t165 + t266 + 0.35e2 * t280) * t286 + (t259 + (t267 + t443) * t282 + t288 * (-t277 + t434)) * t198;
t58 = 0.7e1 * t283 + (0.35e2 * t282 + (15 * t288) + t443) * t294 + (0.21e2 * t280 + t165 + (9 * t287) + (t255 - 0.6e1 * t277) * t288) * t286 + t473;
t140 = t265 + 0.3e1 / 0.2e1 * t282 + t381;
t88 = t140 * t196 + t240 * t398;
t24 = t88 * t363 + t50 * t359 + t49 * t415 + t38 * t404 + (-0.6e1 * t48 * t490 + (t58 + t356) * t240) * pkin(3) + (0.6e1 * t39 * t475 + (t458 * t506 - t51) * t241) * pkin(2);
t17 = t145 * t18 + t24 * t289;
t315 = t17 * t394 + t32 * t491;
t314 = t17 * t395 + t34 * t492;
t310 = 0.4e1 * t313;
t309 = t315 * t399;
t308 = -0.2e1 * t311;
t307 = t310 * t195;
t244 = V_base(3) ^ 2;
t136 = t241 * V_base(2) + t243 * V_base(1);
t130 = t241 * V_base(5) + t243 * V_base(4);
t127 = -t241 * V_base(4) + t243 * V_base(5);
t114 = t170 * t240 + t242 * V_base(1);
t109 = -pkin(2) * t127 + V_base(3);
t107 = pkin(2) * t189 + t135;
t96 = -pkin(3) * t126 + t171;
t72 = t518 * t128;
t71 = t518 * t326;
t31 = 0.1e1 / t32 ^ 2;
t28 = (t129 * t500 + t32 * t495) * t482;
t27 = (t32 * t494 + t34 * t495) * t482;
t26 = (t114 * t500 + t32 * t498) * t482;
t25 = (t114 * t501 + t34 * t498) * t482;
t20 = t110 * t393 + (-t167 * t242 + (t240 * t95 + t462) * qJD(2)) * t173 + (t336 + 0.4e1 * t344 - 0.2e1 * t478) * pkin(3) + ((t159 * t242 - t464) * t429 + (t128 * t289 + t159 * t465 - t463 * t95) * qJD(1) + t173 * t307 + (t161 * t433 - t313 * t489) * t530) * pkin(2);
t19 = t106 * t393 + (-t173 * t464 + t242 * t61) * qJD(2) + (t115 * t242 - t123 * t430) * t240 * t425 + ((-t462 + (-t119 - 0.8e1 * t407) * t240) * t429 + ((t64 - t464) * t243 + (t462 + (t173 * t530 + t119) * t240 + (-pkin(3) + t490 + 0.2e1 * t486) * t525) * t241) * qJD(1)) * pkin(2);
t15 = t314 * t399;
t13 = t188 + ((t19 * t499 + t34 * t392) / t501 - 0.2e1 * (t20 * t499 + t32 * t392) * t34 * t31) * pkin(4) * t278 / (t31 * t34 ^ 2 + 0.1e1) * t78;
t12 = (t356 * t400 - 0.24e2 * pkin(3) * (-0.8e1 / 0.3e1 * t341 + t166) * t390 + 0.96e2 * t117 * t346 * t487 + (t522 * t180 + 0.2e1 * (pkin(3) * t138 - t146 * t240 + t419) * t428 + (t242 * t526 + (t158 * t511 + 0.4e1 * t486) * qJD(2)) * t504) * t359 + (-0.8e1 / 0.3e1 * t333 + (-t133 * t400 + t134 * t180) * t414 + t84 * t180 - t85 * t400 / 0.2e1 + (0.32e2 / 0.3e1 * t472 * t196 + t70 * t242 * t507) * t430 + (0.4e1 * t139 * t242 * t361 + ((0.12e2 * t202 * t205 - 0.4e1 * t203) * t280 + (-0.2e1 * t182 * t486 + t512 * t81) * pkin(3)) * qJD(2)) * pkin(1)) * t404 + 0.24e2 * t38 * t345 + 0.6e1 * ((-0.4e1 * (pkin(3) * t342 + t154 * t526) * t205 + 0.8e1 * t90 * t371) * t282 + (-0.8e1 * t332 + (-t142 * t400 + t143 * t180) * t423 + (0.24e2 * t205 * t408 + t514 * t74) * t430) * pkin(1) + t325) * t412 + t333 * t506 - 0.96e2 * t191 * t343 * t196 + (t140 * t180 + t398 * t428) * t363 + t88 * t321 + (t180 * t86 + t400 * t527) * t415 + 0.12e2 * t49 * t341 - 0.6e1 * t325 * t421 + 0.6e1 * t48 * t362 - t51 * t180 + t58 * t400 + (-0.8e1 * t334 + t318) * t50 - 0.6e1 * t519 * t39) * t289 + t24 * t393 + (0.16e2 * (t187 * t510 - t172 - 0.48e2 * t413 - 0.32e2 * t472) * qJD(2) * t390 - 0.64e2 * t68 * t346 + (t87 + (t104 * t510 + (t112 * t514 - 0.24e2 * t470) * pkin(1)) * t430 + (t308 * t461 - t310 * t470) * pkin(2)) * t355 + 0.24e2 * (-0.32e2 * t192 * t343 + (-0.2e1 / 0.3e1 * t370 - 0.2e1 / 0.3e1 * t369) * t363 * t503 + t116 * t321 - 0.28e2 * t53 * t341 + 0.6e1 * (-0.8e1 / 0.3e1 * t370 - 0.8e1 / 0.3e1 * t369) * t156 * t368 + t55 * t340 + 0.4e1 * (-t323 - t322) * t164 - 0.64e2 / 0.3e1 * t313 * t155 * t366) * t467 - 0.48e2 * t36 * t345 - 0.8e1 * (t79 * t331 + 0.6e1 * (-pkin(2) * t164 * t307 - t430 * t59) * t504 + t480 + t483) * t412 - 0.4e1 * t153 * t373 + 0.32e2 * t204 * t308 * t389 * t505 - 0.96e2 * t83 * t335 - 0.96e2 * t164 * t306 - 0.48e2 * t47 * t341 + 0.8e1 * t483 * t421 - 0.8e1 * t43 * t362 + 0.8e1 * t519 * t37 - 0.8e1 * t531 * t52 + (0.32e2 * t334 + 0.96e2 * t329) * t46) * t145 + t18 * t132;
t11 = t128 * t15 - t309 * t326;
t10 = t128 * t309 + t15 * t326;
t9 = 0.1e1 / t11 ^ 2;
t7 = -t10 * t127 - t11 * t130;
t6 = t10 * t130 - t11 * t127;
t5 = -t10 * t107 - t11 * t136;
t4 = t10 * t136 - t107 * t11;
t3 = (-t315 * t484 + (t20 * t491 + t32 * t485 / 0.8e1 + t12 * t394 + (t19 * t502 + t34 * t396) * t17) * t75) * t459;
t2 = (-t314 * t484 + (t12 * t395 + t19 * t492 - t34 * t485 / 0.8e1 + (t20 * t502 + t32 * t396) * t17) * t75) * t459;
t1 = t189 + ((t128 * t3 - t15 * t72 + t2 * t326 + t309 * t71) / t11 - (t2 * t128 + t15 * t71 - t3 * t326 + t309 * t72) * t10 * t9) / (t10 ^ 2 * t9 + 0.1e1);
t8 = m(2) * (t135 ^ 2 + t136 ^ 2 + t244) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t244) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t171 ^ 2) / 0.2e1 + m(5) * (t25 ^ 2 + t26 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t109 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t109 * mrSges(3,2) - t4 * mrSges(3,3) + Ifges(3,1) * t7 / 0.2e1) * t7 + (-t96 * mrSges(5,1) + t25 * mrSges(5,3) + Ifges(5,2) * t28 / 0.2e1) * t28 + (t135 * mrSges(2,1) - t136 * mrSges(2,2) + Ifges(2,3) * t189 / 0.2e1) * t189 + (t113 * mrSges(4,1) - t114 * mrSges(4,2) + Ifges(4,3) * t523) * t188 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t109 * mrSges(3,1) + t5 * mrSges(3,3) + Ifges(3,4) * t7 + Ifges(3,2) * t6 / 0.2e1) * t6 + (t96 * mrSges(5,2) - t26 * mrSges(5,3) + Ifges(5,4) * t28 + Ifges(5,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(2,2) - t135 * mrSges(2,3) + Ifges(2,5) * t189 + Ifges(2,1) * t130 / 0.2e1) * t130 + (t171 * mrSges(4,2) - t113 * mrSges(4,3) + Ifges(4,1) * t494 + Ifges(4,5) * t188) * t129 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t26 * mrSges(5,1) - t25 * mrSges(5,2) + Ifges(5,5) * t27 + Ifges(5,6) * t28 + Ifges(5,3) * t13 / 0.2e1) * t13 + (-V_base(3) * mrSges(2,1) + t136 * mrSges(2,3) + Ifges(2,4) * t130 + Ifges(2,6) * t189 + Ifges(2,2) * t127 / 0.2e1) * t127 + (-t171 * mrSges(4,1) + t114 * mrSges(4,3) + Ifges(4,4) * t129 + Ifges(4,2) * t495 + Ifges(4,6) * t188) * t126 + (t4 * mrSges(3,1) - t5 * mrSges(3,2) + Ifges(3,5) * t7 + Ifges(3,6) * t6 + Ifges(3,3) * t1 / 0.2e1) * t1;
T = t8;
