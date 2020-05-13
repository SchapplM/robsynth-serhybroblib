% Calculate kinetic energy for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2TE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:10:14
% EndTime: 2020-05-01 20:10:31
% DurationCPUTime: 14.34s
% Computational Cost: add. (3357->547), mult. (5380->799), div. (0->0), fcn. (6176->24), ass. (0->313)
t545 = -Icges(3,4) - Icges(10,4);
t544 = Icges(3,1) + Icges(10,1);
t543 = Icges(3,2) + Icges(10,2);
t404 = cos(qJ(2));
t542 = t545 * t404;
t398 = sin(qJ(2));
t541 = t545 * t398;
t540 = Icges(3,5) + Icges(10,5);
t539 = Icges(3,6) + Icges(10,6);
t538 = t543 * t404 - t541;
t537 = -t544 * t398 + t542;
t536 = Icges(9,3) + Icges(4,3);
t399 = sin(qJ(1));
t405 = cos(qJ(1));
t535 = t538 * t399 + t539 * t405;
t534 = -t539 * t399 + t538 * t405;
t533 = t537 * t399 - t540 * t405;
t532 = t540 * t399 + t537 * t405;
t531 = t543 * t398 + t542;
t530 = t544 * t404 + t541;
t392 = sin(pkin(19));
t395 = cos(pkin(19));
t319 = t392 * t398 - t395 * t404;
t323 = t392 * t404 + t395 * t398;
t397 = sin(qJ(3));
t403 = cos(qJ(3));
t254 = t319 * t403 + t323 * t397;
t481 = t403 * t404;
t328 = -t397 * t398 + t481;
t482 = t398 * t403;
t330 = t397 * t404 + t482;
t529 = Icges(4,5) * t328 - Icges(9,5) * t254 - Icges(4,6) * t330;
t528 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t400 = sin(pkin(18));
t401 = sin(pkin(17));
t406 = cos(pkin(17));
t512 = cos(pkin(18));
t329 = t400 * t406 - t401 * t512;
t331 = t401 * t400 + t406 * t512;
t258 = t329 * t398 + t331 * t404;
t259 = t329 * t404 - t331 * t398;
t527 = Icges(7,5) * t259 - Icges(7,6) * t258 - t540 * t398 - t539 * t404;
t391 = sin(pkin(20));
t394 = cos(pkin(20));
t321 = -t391 * t512 + t394 * t400;
t325 = t400 * t391 + t394 * t512;
t387 = pkin(22) + pkin(21);
t374 = sin(t387);
t375 = cos(t387);
t242 = t321 * t374 + t325 * t375;
t525 = Icges(5,4) * t242;
t524 = Icges(7,4) * t258;
t494 = Icges(9,4) * t254;
t507 = pkin(1) * t398;
t311 = t507 * t399;
t312 = t507 * t405;
t364 = -qJD(2) * t405 + V_base(5);
t365 = qJD(2) * t399 + V_base(4);
t523 = -t365 * t311 + t312 * t364;
t318 = t392 * t403 + t395 * t397;
t322 = t392 * t397 - t395 * t403;
t253 = -t318 * t404 + t322 * t398;
t256 = -t319 * t397 + t323 * t403;
t327 = qJD(3) * t399 + t365;
t376 = V_base(6) + qJD(1);
t388 = qJD(2) + qJD(3);
t518 = -t388 * t405 + V_base(5);
t522 = (-t536 * t405 + (Icges(9,6) * t253 + t529) * t399) * t518 + (Icges(4,5) * t330 + Icges(9,5) * t256 + Icges(4,6) * t328 - Icges(9,6) * t254) * t376 + ((-Icges(9,6) * t256 + t529) * t405 + t536 * t399) * t327;
t497 = Icges(7,4) * t259;
t449 = -Icges(7,2) * t258 + t497;
t187 = -Icges(7,6) * t405 + t399 * t449;
t188 = Icges(7,6) * t399 + t405 * t449;
t455 = Icges(7,1) * t259 - t524;
t189 = -Icges(7,5) * t405 + t399 * t455;
t190 = Icges(7,5) * t399 + t405 * t455;
t202 = Icges(7,2) * t259 + t524;
t203 = Icges(7,1) * t258 + t497;
t520 = (-t202 * t258 + t203 * t259 - t530 * t398 + t531 * t404) * t376 + (-t188 * t258 + t190 * t259 - t532 * t398 + t534 * t404) * t365 + (-t187 * t258 + t189 * t259 - t533 * t398 + t535 * t404) * t364;
t519 = (Icges(7,5) * t258 + Icges(7,6) * t259 - t539 * t398 + t540 * t404) * t376 + (t528 * t399 + t527 * t405) * t365 + (t527 * t399 - t528 * t405) * t364;
t511 = pkin(1) * t404;
t510 = pkin(15) * t399;
t509 = pkin(15) * t405;
t508 = (-t318 * t398 - t322 * t404) * pkin(2);
t506 = cos(pkin(22));
t505 = sin(pkin(22));
t504 = Icges(9,1) * t254;
t503 = Icges(2,4) * t399;
t500 = Icges(4,4) * t328;
t499 = Icges(4,4) * t330;
t243 = -t321 * t375 + t325 * t374;
t498 = Icges(5,4) * t243;
t419 = t400 * t505 + t506 * t512;
t496 = Icges(8,4) * t419;
t324 = -t400 * t506 + t505 * t512;
t495 = Icges(8,4) * t324;
t493 = Icges(9,4) * t256;
t183 = -Icges(9,5) * t405 + (Icges(9,4) * t253 - t504) * t399;
t489 = t183 * t254;
t184 = Icges(9,5) * t399 + (-t493 - t504) * t405;
t488 = t184 * t254;
t195 = Icges(9,1) * t256 - t494;
t487 = t195 * t254;
t486 = t243 * t399;
t485 = t243 * t405;
t367 = t388 * t399;
t390 = sin(pkin(21));
t393 = cos(pkin(21));
t480 = (t324 * t390 - t393 * t419) * pkin(4);
t373 = pkin(5) * t397 + pkin(1);
t479 = pkin(5) * t481 - t373 * t398 + t507;
t389 = qJ(3) + qJ(2);
t380 = sin(t389);
t478 = Icges(11,4) * t380;
t381 = cos(t389);
t477 = Icges(11,4) * t381;
t476 = qJD(4) * t243;
t475 = t376 * t509 + V_base(2);
t469 = t399 * V_base(4);
t474 = pkin(15) * t469 + V_base(3);
t473 = V_base(5) * pkin(13) + V_base(1);
t470 = V_base(4) * pkin(13);
t468 = t405 * V_base(5);
t467 = t311 - t510;
t466 = t364 * t511 + t473;
t336 = -rSges(9,1) * t395 + rSges(9,2) * t392;
t337 = rSges(9,1) * t392 + rSges(9,2) * t395;
t465 = -t336 * t397 + t337 * t403;
t263 = t479 * t399;
t464 = -t263 + t467;
t290 = pkin(5) * t482 + (-pkin(1) + t373) * t404;
t463 = t290 * t518 + t466;
t462 = -rSges(3,1) * t398 - rSges(3,2) * t404;
t396 = sin(qJ(4));
t402 = cos(qJ(4));
t461 = rSges(6,1) * t402 - rSges(6,2) * t396;
t460 = -rSges(10,1) * t398 - rSges(10,2) * t404;
t293 = -t367 + t327;
t459 = rSges(11,1) * t381 - rSges(11,2) * t380;
t457 = Icges(4,1) * t328 - t499;
t456 = -Icges(5,1) * t242 - t498;
t454 = -Icges(8,1) * t419 - t495;
t451 = -Icges(4,2) * t330 + t500;
t450 = -Icges(5,2) * t243 - t525;
t448 = -Icges(8,2) * t324 - t496;
t444 = -Icges(5,5) * t242 - Icges(5,6) * t243;
t442 = -Icges(8,5) * t419 - Icges(8,6) * t324;
t222 = rSges(6,3) * t325 - t321 * t461;
t223 = rSges(6,3) * t321 + t325 * t461;
t440 = t222 * t374 - t223 * t375;
t353 = -t400 * rSges(5,1) - rSges(5,2) * t512;
t360 = rSges(5,1) * t512 - t400 * rSges(5,2);
t272 = -t353 * t391 + t360 * t394;
t437 = t353 * t394 + t360 * t391;
t439 = -t272 * t375 + t374 * t437;
t361 = -t400 * pkin(9) + pkin(11) * t512;
t362 = pkin(9) * t512 + t400 * pkin(11);
t278 = t361 * t394 + t362 * t391;
t279 = -t361 * t391 + t362 * t394;
t438 = t278 * t374 - t279 * t375;
t358 = -t400 * rSges(7,1) + rSges(7,2) * t512;
t359 = rSges(7,1) * t512 + t400 * rSges(7,2);
t436 = t401 * t358 - t359 * t406;
t435 = Icges(11,1) * t381 - t478;
t434 = -Icges(11,2) * t380 + t477;
t433 = Icges(11,5) * t381 - Icges(11,6) * t380;
t432 = -t470 + t475;
t431 = t474 + t523;
t430 = t396 * t242;
t429 = t242 * t402;
t428 = -t376 * t312 - t365 * t511 + t475;
t427 = -pkin(15) * t468 + t474;
t426 = (-Icges(5,3) * t405 + t399 * t444) * V_base(5) + (Icges(5,3) * t399 + t405 * t444) * t293 + (Icges(5,5) * t243 - Icges(5,6) * t242) * t376;
t339 = V_base(4) + t367;
t422 = (-Icges(11,3) * t405 + t399 * t433) * t518 + (Icges(11,3) * t399 + t405 * t433) * t339 + (Icges(11,5) * t380 + Icges(11,6) * t381) * t376;
t418 = (-Icges(8,3) * t405 + t399 * t442) * V_base(5) + (Icges(8,3) * t399 + t405 * t442) * V_base(4) + (Icges(8,5) * t324 - Icges(8,6) * t419) * t376;
t417 = t428 - t470;
t416 = t427 + t523;
t264 = t479 * t405;
t415 = t376 * t264 - t290 * t327 + t417;
t414 = t327 * t263 - t264 * t518 + t416;
t172 = -Icges(5,6) * t405 + t399 * t450;
t173 = Icges(5,6) * t399 + t405 * t450;
t174 = -Icges(5,5) * t405 + t399 * t456;
t175 = Icges(5,5) * t399 + t405 * t456;
t177 = -Icges(5,2) * t242 + t498;
t178 = Icges(5,1) * t243 - t525;
t413 = (-t173 * t243 - t175 * t242) * t293 + (-t172 * t243 - t174 * t242) * V_base(5) + (-t177 * t243 - t178 * t242) * t376;
t237 = -Icges(4,6) * t405 + t399 * t451;
t238 = Icges(4,6) * t399 + t405 * t451;
t239 = -Icges(4,5) * t405 + t399 * t457;
t240 = Icges(4,5) * t399 + t405 * t457;
t266 = Icges(4,2) * t328 + t499;
t267 = Icges(4,1) * t330 + t500;
t411 = (-t238 * t330 + t240 * t328) * t327 + (-t237 * t330 + t239 * t328) * t518 + (-t266 * t330 + t267 * t328) * t376;
t283 = -Icges(11,6) * t405 + t399 * t434;
t284 = Icges(11,6) * t399 + t405 * t434;
t285 = -Icges(11,5) * t405 + t399 * t435;
t286 = Icges(11,5) * t399 + t405 * t435;
t314 = Icges(11,2) * t381 + t478;
t315 = Icges(11,1) * t380 + t477;
t410 = (-t284 * t380 + t286 * t381) * t339 + (-t283 * t380 + t285 * t381) * t518 + (-t314 * t380 + t315 * t381) * t376;
t231 = -Icges(8,6) * t405 + t399 * t448;
t232 = Icges(8,6) * t399 + t405 * t448;
t233 = -Icges(8,5) * t405 + t399 * t454;
t234 = Icges(8,5) * t399 + t405 * t454;
t261 = -Icges(8,2) * t419 + t495;
t262 = Icges(8,1) * t324 - t496;
t407 = (-t232 * t324 - t234 * t419) * V_base(4) + (-t231 * t324 - t233 * t419) * V_base(5) + (-t261 * t324 - t262 * t419) * t376;
t382 = Icges(2,4) * t405;
t357 = rSges(2,1) * t405 - rSges(2,2) * t399;
t356 = rSges(3,1) * t404 - rSges(3,2) * t398;
t355 = rSges(10,1) * t404 - rSges(10,2) * t398;
t354 = rSges(6,1) * t396 + rSges(6,2) * t402;
t352 = rSges(2,1) * t399 + rSges(2,2) * t405;
t351 = Icges(2,1) * t405 - t503;
t350 = Icges(2,1) * t399 + t382;
t347 = -Icges(2,2) * t399 + t382;
t346 = Icges(2,2) * t405 + t503;
t343 = Icges(2,5) * t405 - Icges(2,6) * t399;
t342 = Icges(2,5) * t399 + Icges(2,6) * t405;
t335 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t334 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t333 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t316 = rSges(11,1) * t380 + rSges(11,2) * t381;
t309 = rSges(3,3) * t399 + t405 * t462;
t308 = rSges(10,3) * t399 + t405 * t460;
t307 = -rSges(3,3) * t405 + t399 * t462;
t306 = -rSges(10,3) * t405 + t399 * t460;
t288 = rSges(11,3) * t399 + t405 * t459;
t287 = -rSges(11,3) * t405 + t399 * t459;
t280 = (rSges(4,1) * t398 + rSges(4,2) * t404) * t403 + t397 * (rSges(4,1) * t404 - rSges(4,2) * t398);
t277 = -t358 * t406 - t401 * t359;
t276 = (rSges(4,1) * t403 - rSges(4,2) * t397) * t404 - t398 * (rSges(4,1) * t397 + rSges(4,2) * t403);
t275 = V_base(5) * rSges(2,3) - t352 * t376 + t473;
t274 = t357 * t376 + V_base(2) + (-pkin(13) - rSges(2,3)) * V_base(4);
t273 = (-rSges(8,1) * t512 + rSges(8,2) * t400) * t506 - t505 * (rSges(8,1) * t400 + rSges(8,2) * t512);
t270 = (-rSges(8,1) * t506 - rSges(8,2) * t505) * t400 + (rSges(8,1) * t505 - rSges(8,2) * t506) * t512;
t269 = -t336 * t403 - t337 * t397;
t268 = t352 * V_base(4) - t357 * V_base(5) + V_base(3);
t251 = rSges(4,3) * t399 + t276 * t405;
t250 = -rSges(4,3) * t405 + t276 * t399;
t249 = rSges(8,3) * t399 + t273 * t405;
t248 = -rSges(8,3) * t405 + t273 * t399;
t246 = t253 * pkin(2);
t245 = (t324 * t393 + t390 * t419) * pkin(4);
t226 = qJD(4) * t242 + t376;
t225 = t508 * t405;
t224 = t508 * t399;
t221 = t396 * t399 - t405 * t429;
t220 = t399 * t402 + t405 * t430;
t219 = -t396 * t405 - t399 * t429;
t218 = t399 * t430 - t402 * t405;
t217 = t480 * t405;
t216 = t480 * t399;
t215 = t356 * t364 + (-t307 - t510) * t376 + t473;
t214 = t309 * t376 - t356 * t365 + t432;
t213 = t277 * t398 - t404 * t436;
t212 = t277 * t404 + t398 * t436;
t211 = t398 * t269 + t404 * t465;
t210 = t269 * t404 - t398 * t465;
t209 = t405 * t476 + t293;
t208 = t399 * t476 + V_base(5);
t207 = t278 * t375 + t279 * t374;
t206 = t374 * t272 + t375 * t437;
t205 = rSges(7,3) * t399 + t212 * t405;
t204 = -rSges(7,3) * t405 + t212 * t399;
t200 = t438 * t405;
t199 = t438 * t399;
t198 = rSges(9,3) * t399 + t210 * t405;
t197 = -rSges(9,3) * t405 + t210 * t399;
t196 = t307 * t365 - t309 * t364 + t427;
t194 = -Icges(9,2) * t254 + t493;
t192 = -rSges(5,3) * t405 + t399 * t439;
t191 = rSges(5,3) * t399 + t405 * t439;
t182 = Icges(9,6) * t399 + (-Icges(9,2) * t256 - t494) * t405;
t181 = -Icges(9,6) * t405 + (Icges(9,2) * t253 - t494) * t399;
t169 = t222 * t375 + t223 * t374;
t168 = t270 * V_base(5) + (-t248 + t467) * t376 + t466;
t167 = t249 * t376 + (-pkin(13) - t270) * V_base(4) + t428;
t166 = t280 * t518 + (-t250 + t467) * t376 + t466;
t165 = t251 * t376 - t280 * t327 + t417;
t164 = Icges(6,5) * t242 + (Icges(6,1) * t402 - Icges(6,4) * t396) * t243;
t163 = Icges(6,6) * t242 + (Icges(6,4) * t402 - Icges(6,2) * t396) * t243;
t162 = Icges(6,3) * t242 + (Icges(6,5) * t402 - Icges(6,6) * t396) * t243;
t161 = t354 * t399 + t405 * t440;
t160 = -t354 * t405 + t399 * t440;
t159 = t248 * V_base(4) + (-t249 - t509) * V_base(5) + t431;
t158 = t250 * t327 - t251 * t518 + t416;
t157 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t485;
t156 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t486;
t155 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t485;
t154 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t486;
t153 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t485;
t152 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t486;
t151 = -t246 * t518 + t355 * t364 + (-t224 - t306 - t510) * t376 + t473;
t150 = t246 * t327 - t355 * t365 + (t225 + t308) * t376 + t432;
t149 = -V_base(5) * pkin(16) + t213 * t364 + (pkin(14) * t399 - t204) * t376 + t473;
t148 = -t213 * t365 + V_base(2) + (-pkin(13) + pkin(16)) * V_base(4) + (-pkin(14) * t405 + t205) * t376;
t147 = t211 * t518 + (-t197 - t510) * t376 + t473;
t146 = t198 * t376 - t211 * t327 + t432;
t145 = t204 * t365 - t205 * t364 + V_base(3) + (t468 - t469) * pkin(14);
t144 = t197 * t327 - t198 * t518 + t427;
t143 = t224 * t327 - t225 * t518 + t306 * t365 - t308 * t364 + t427;
t142 = t245 * V_base(5) + t316 * t518 + (-t216 - t287 + t467) * t376 + t466;
t141 = -t316 * t339 + (-pkin(13) - t245) * V_base(4) + (t217 + t288) * t376 + t428;
t140 = t216 * V_base(4) + t287 * t339 - t288 * t518 + (-t217 - t509) * V_base(5) + t431;
t139 = t206 * V_base(5) + (-t192 + t464) * t376 + t463;
t138 = t191 * t376 - t206 * t293 + t415;
t137 = -t191 * V_base(5) + t192 * t293 + t414;
t136 = -t160 * t226 + t169 * t208 + t207 * V_base(5) + (-t199 + t464) * t376 + t463;
t135 = t161 * t226 - t169 * t209 + t200 * t376 - t207 * t293 + t415;
t134 = t160 * t209 - t161 * t208 + t199 * t293 - t200 * V_base(5) + t414;
t1 = m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t158 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(8) * (t159 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + (Icges(1,4) * V_base(4) + Icges(1,2) * V_base(5) + t342 * t376 + (t346 * V_base(5) + t347 * V_base(4) - t418 - t426) * t405 + (t350 * V_base(5) + t351 * V_base(4) + t407 + t413) * t399) * V_base(5) / 0.2e1 + ((-t522 - t422) * t405 + ((t182 * t253 - t488) * t327 + (t181 * t253 - t489) * t518 + (t194 * t253 - t487) * t376 + t411 + t410) * t399) * t518 / 0.2e1 + m(3) * (t196 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t209 * ((t153 * t485 + t155 * t220 + t157 * t221) * t209 + (t152 * t485 + t154 * t220 + t156 * t221) * t208 + (t162 * t485 + t163 * t220 + t164 * t221) * t226) / 0.2e1 + m(2) * (t268 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(11) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(9) * (t144 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(7) * (t145 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(10) * (t143 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + t339 * (t422 * t399 + t410 * t405) / 0.2e1 + t293 * (t426 * t399 + t413 * t405) / 0.2e1 + t208 * ((t153 * t486 + t155 * t218 + t157 * t219) * t209 + (t152 * t486 + t154 * t218 + t156 * t219) * t208 + (t162 * t486 + t163 * t218 + t164 * t219) * t226) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t226 * ((t152 * t208 + t153 * t209 + t162 * t226) * t242 + ((-t155 * t396 + t157 * t402) * t209 + (-t154 * t396 + t156 * t402) * t208 + (-t163 * t396 + t164 * t402) * t226) * t243) / 0.2e1 + ((t284 * t381 + t286 * t380) * t339 + (-t173 * t242 + t175 * t243) * t293 + (-t232 * t419 + t234 * t324 + t343) * V_base(4) + (-t182 * t254 + t184 * t256 + t238 * t328 + t240 * t330) * t327 + (t188 * t259 + t190 * t258 + t534 * t398 + t532 * t404) * t365 + (t187 * t259 + t189 * t258 + t535 * t398 + t533 * t404) * t364 + (-t177 * t242 + t178 * t243 - t194 * t254 + t195 * t256 + t202 * t259 + t203 * t258 - t261 * t419 + t262 * t324 + t266 * t328 + t267 * t330 + t314 * t381 + t315 * t380 + t531 * t398 + t530 * t404 + Icges(2,3)) * t376 + (-t172 * t242 + t174 * t243 - t231 * t419 + t233 * t324 + t342) * V_base(5) + (-t181 * t254 + t183 * t256 + t237 * t328 + t239 * t330 + t283 * t381 + t285 * t380) * t518) * t376 / 0.2e1 + (((-t182 * t256 - t488) * t327 + (-t181 * t256 - t489) * t518 + (-t194 * t256 - t487) * t376 + t411) * t405 + t522 * t399) * t327 / 0.2e1 + m(1) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + (t519 * t399 + t520 * t405) * t365 / 0.2e1 + (t520 * t399 - t519 * t405) * t364 / 0.2e1 + (t343 * t376 + t399 * t418 + t405 * t407 + (-t346 * t399 + t350 * t405 + Icges(1,4)) * V_base(5) + (-t347 * t399 + t351 * t405 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
