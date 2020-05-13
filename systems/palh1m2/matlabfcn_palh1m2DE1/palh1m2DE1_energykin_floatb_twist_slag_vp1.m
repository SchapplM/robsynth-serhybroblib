% Calculate kinetic energy for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:47
% EndTime: 2020-05-01 20:56:03
% DurationCPUTime: 14.96s
% Computational Cost: add. (3142->589), mult. (5270->850), div. (0->0), fcn. (6039->24), ass. (0->319)
t571 = -Icges(3,4) - Icges(10,4);
t570 = Icges(3,1) + Icges(10,1);
t569 = Icges(3,2) + Icges(10,2);
t405 = cos(qJ(2));
t568 = t571 * t405;
t399 = sin(qJ(2));
t567 = t571 * t399;
t566 = Icges(3,5) + Icges(10,5);
t565 = Icges(3,6) + Icges(10,6);
t564 = t405 * t569 - t567;
t563 = -t399 * t570 + t568;
t393 = sin(pkin(19));
t396 = cos(pkin(19));
t303 = t393 * t399 - t396 * t405;
t307 = t393 * t405 + t396 * t399;
t398 = sin(qJ(3));
t404 = cos(qJ(3));
t241 = t303 * t404 + t307 * t398;
t514 = Icges(9,4) * t241;
t561 = Icges(9,3) + Icges(4,3);
t400 = sin(qJ(1));
t406 = cos(qJ(1));
t560 = t564 * t400 + t565 * t406;
t559 = -t565 * t400 + t564 * t406;
t558 = t563 * t400 - t566 * t406;
t557 = t566 * t400 + t563 * t406;
t556 = t399 * t569 + t568;
t555 = t405 * t570 + t567;
t392 = sin(pkin(20));
t395 = cos(pkin(20));
t401 = sin(pkin(18));
t407 = cos(pkin(18));
t305 = -t407 * t392 + t395 * t401;
t309 = t392 * t401 + t395 * t407;
t387 = pkin(22) + pkin(21);
t373 = sin(t387);
t374 = cos(t387);
t236 = t305 * t373 + t309 * t374;
t237 = -t305 * t374 + t309 * t373;
t525 = cos(pkin(22));
t474 = t401 * t525;
t390 = sin(pkin(22));
t495 = t407 * t390;
t308 = -t474 + t495;
t434 = t401 * t390 + t407 * t525;
t554 = Icges(5,5) * t236 + Icges(8,5) * t434 + Icges(5,6) * t237 + Icges(8,6) * t308;
t496 = t404 * t405;
t313 = -t398 * t399 + t496;
t497 = t399 * t404;
t315 = t398 * t405 + t497;
t553 = Icges(4,5) * t313 - Icges(9,5) * t241 - Icges(4,6) * t315;
t552 = Icges(2,2) + Icges(5,3) + Icges(8,3);
t551 = Icges(3,3) + Icges(7,3) + Icges(10,3);
t402 = sin(pkin(17));
t408 = cos(pkin(17));
t314 = t401 * t408 - t407 * t402;
t316 = t402 * t401 + t407 * t408;
t246 = t314 * t399 + t316 * t405;
t247 = t314 * t405 - t316 * t399;
t550 = Icges(7,5) * t247 - Icges(7,6) * t246 - t566 * t399 - t565 * t405;
t549 = Icges(5,4) * t236;
t548 = Icges(7,4) * t246;
t294 = t434 * pkin(4);
t391 = sin(pkin(21));
t532 = pkin(4) * t391;
t330 = t474 * t532;
t394 = cos(pkin(21));
t484 = t391 * t495;
t228 = pkin(4) * t484 - t294 * t394 + pkin(15) - t330;
t386 = pkin(1) * t399;
t475 = t228 - t386;
t481 = t400 * V_base(4);
t317 = t406 * V_base(5) - t481;
t302 = t393 * t404 + t396 * t398;
t442 = t393 * t398 - t396 * t404;
t547 = (-t302 * t399 - t405 * t442) * pkin(2);
t546 = t386 - pkin(15);
t355 = -rSges(8,1) * t401 - rSges(8,2) * t407;
t356 = rSges(8,1) * t407 - rSges(8,2) * t401;
t469 = t355 * t390 - t356 * t525;
t240 = -t302 * t405 + t399 * t442;
t243 = -t303 * t398 + t307 * t404;
t363 = qJD(2) * t400 + V_base(4);
t312 = qJD(3) * t400 + t363;
t388 = qJD(2) + qJD(3);
t328 = -t388 * t406 + V_base(5);
t375 = V_base(6) + qJD(1);
t545 = (Icges(4,5) * t315 + Icges(9,5) * t243 + Icges(4,6) * t313 - Icges(9,6) * t241) * t375 + ((-Icges(9,6) * t243 + t553) * t406 + t561 * t400) * t312 + (-t561 * t406 + (Icges(9,6) * t240 + t553) * t400) * t328;
t517 = Icges(7,4) * t247;
t453 = -Icges(7,2) * t246 + t517;
t187 = -Icges(7,6) * t406 + t400 * t453;
t188 = Icges(7,6) * t400 + t406 * t453;
t459 = Icges(7,1) * t247 - t548;
t189 = -Icges(7,5) * t406 + t400 * t459;
t190 = Icges(7,5) * t400 + t406 * t459;
t197 = Icges(7,2) * t247 + t548;
t198 = Icges(7,1) * t246 + t517;
t362 = -qJD(2) * t406 + V_base(5);
t544 = (-t197 * t246 + t198 * t247 - t555 * t399 + t556 * t405) * t375 + (-t188 * t246 + t190 * t247 - t557 * t399 + t559 * t405) * t363 + (-t187 * t246 + t189 * t247 - t558 * t399 + t560 * t405) * t362;
t543 = (Icges(7,5) * t246 + Icges(7,6) * t247 - t565 * t399 + t566 * t405) * t375 + (t551 * t400 + t550 * t406) * t363 + (t550 * t400 - t551 * t406) * t362;
t380 = Icges(2,4) * t406;
t523 = Icges(2,4) * t400;
t540 = (-Icges(5,5) * t237 - Icges(8,5) * t308 + Icges(5,6) * t236 + Icges(8,6) * t434) * t375 + (t554 * t400 + t552 * t406 + t523) * V_base(5) + (-t552 * t400 + t554 * t406 + t380) * V_base(4);
t454 = -Icges(5,2) * t237 - t549;
t172 = -Icges(5,6) * t406 + t400 * t454;
t173 = Icges(5,6) * t400 + t406 * t454;
t518 = Icges(5,4) * t237;
t460 = -Icges(5,1) * t236 - t518;
t174 = -Icges(5,5) * t406 + t400 * t460;
t175 = Icges(5,5) * t400 + t406 * t460;
t177 = -Icges(5,2) * t236 + t518;
t178 = Icges(5,1) * t237 - t549;
t516 = Icges(8,4) * t434;
t452 = -Icges(8,2) * t308 - t516;
t223 = -Icges(8,6) * t406 + t400 * t452;
t224 = Icges(8,6) * t400 + t406 * t452;
t515 = Icges(8,4) * t308;
t458 = -Icges(8,1) * t434 - t515;
t225 = -Icges(8,5) * t406 + t400 * t458;
t226 = Icges(8,5) * t400 + t406 * t458;
t249 = -Icges(8,2) * t434 + t515;
t250 = Icges(8,1) * t308 - t516;
t539 = (-t177 * t237 - t178 * t236 - t249 * t308 - t250 * t434) * t375 + (Icges(2,1) * t400 - t172 * t237 - t174 * t236 - t223 * t308 - t225 * t434 + t380) * V_base(5) + (Icges(2,1) * t406 - t173 * t237 - t175 * t236 - t224 * t308 - t226 * t434 - t523) * V_base(4);
t533 = pkin(1) * t405;
t531 = pkin(5) * t404;
t530 = pkin(1) * qJD(2);
t529 = pkin(5) * qJD(3);
t528 = t400 * rSges(3,3);
t527 = t406 * rSges(3,3);
t345 = rSges(3,1) * t399 + rSges(3,2) * t405;
t526 = -pkin(15) + t345;
t524 = Icges(9,1) * t241;
t520 = Icges(4,4) * t313;
t519 = Icges(4,4) * t315;
t513 = Icges(9,4) * t243;
t183 = Icges(9,5) * t400 + (-t513 - t524) * t406;
t509 = t183 * t241;
t184 = -Icges(9,5) * t406 + (Icges(9,4) * t240 - t524) * t400;
t508 = t184 * t241;
t193 = Icges(9,1) * t243 - t514;
t507 = t193 * t241;
t343 = pkin(9) * t392 + pkin(11) * t395;
t344 = pkin(9) * t395 - pkin(11) * t392;
t259 = t343 * t407 - t344 * t401;
t260 = t343 * t401 + t344 * t407;
t506 = (-t259 * t373 + t260 * t374) * t375;
t505 = t237 * t400;
t504 = t237 * t406;
t372 = pkin(5) * t398 + pkin(1);
t500 = t372 * t399;
t499 = t375 * t406;
t389 = qJ(3) + qJ(2);
t378 = sin(t389);
t494 = Icges(11,4) * t378;
t379 = cos(t389);
t493 = Icges(11,4) * t379;
t492 = qJD(1) * t400;
t491 = qJD(4) * t237;
t490 = t546 * qJD(1);
t382 = V_base(5) * pkin(13);
t489 = t382 + V_base(1);
t488 = pkin(5) * t496;
t485 = t405 * t530;
t383 = V_base(4) * pkin(15);
t479 = t405 * V_base(4);
t483 = t479 * t531 - V_base(4) * t500 + t383;
t482 = t399 * V_base(5);
t480 = t400 * V_base(6);
t478 = t405 * V_base(5);
t473 = -pkin(4) * t308 * t394 - pkin(13);
t438 = t488 - t500;
t472 = t438 * qJD(2) + t313 * t529 + V_base(3);
t324 = -rSges(9,1) * t396 + rSges(9,2) * t393;
t325 = rSges(9,1) * t393 + rSges(9,2) * t396;
t470 = -t324 * t398 + t325 * t404;
t289 = -pkin(15) - t438;
t468 = t375 * t289;
t397 = sin(qJ(4));
t403 = cos(qJ(4));
t467 = rSges(6,1) * t403 - rSges(6,2) * t397;
t466 = -rSges(10,1) * t399 - t405 * rSges(10,2);
t465 = -t399 * t530 + V_base(3);
t464 = -V_base(4) * pkin(13) + V_base(2);
t463 = rSges(11,1) * t379 - rSges(11,2) * t378;
t461 = Icges(4,1) * t313 - t519;
t455 = -Icges(4,2) * t315 + t520;
t215 = rSges(6,3) * t309 - t305 * t467;
t216 = rSges(6,3) * t305 + t309 * t467;
t444 = t215 * t373 - t216 * t374;
t351 = -t401 * rSges(7,1) + rSges(7,2) * t407;
t353 = rSges(7,1) * t407 + rSges(7,2) * t401;
t443 = t402 * t351 - t353 * t408;
t441 = Icges(11,1) * t379 - t494;
t440 = -Icges(11,2) * t378 + t493;
t439 = Icges(11,5) * t379 - Icges(11,6) * t378;
t290 = pkin(5) * t497 + t372 * t405;
t437 = -t400 * t485 + V_base(2);
t436 = -qJD(2) * t290 - t315 * t529;
t435 = pkin(15) * t499 + t464;
t433 = t372 * t478 + t482 * t531 + t489 + (t480 + t492) * t289;
t326 = rSges(5,1) * t392 - rSges(5,2) * t395;
t327 = rSges(5,1) * t395 + rSges(5,2) * t392;
t256 = t326 * t407 - t327 * t401;
t257 = t326 * t401 + t327 * t407;
t432 = t375 * (-t256 * t373 + t257 * t374);
t431 = t397 * t236;
t430 = t236 * t403;
t429 = -pkin(13) - t290;
t428 = pkin(1) * t478 - t406 * t485 + V_base(1);
t427 = -t317 * pkin(15) + V_base(3);
t329 = t388 * t400 + V_base(4);
t423 = (-Icges(11,3) * t406 + t400 * t439) * t328 + (Icges(11,3) * t400 + t406 * t439) * t329 + (Icges(11,5) * t378 + Icges(11,6) * t379) * t375;
t420 = t372 * t482 + (-pkin(15) - t488) * V_base(5);
t419 = rSges(5,3) * t375 + t436;
t418 = t400 * t490 + t382 + t428;
t231 = -Icges(4,6) * t406 + t400 * t455;
t232 = Icges(4,6) * t400 + t406 * t455;
t233 = -Icges(4,5) * t406 + t400 * t461;
t234 = Icges(4,5) * t400 + t406 * t461;
t252 = Icges(4,2) * t313 + t519;
t253 = Icges(4,1) * t315 + t520;
t414 = (-t232 * t315 + t234 * t313) * t312 + (-t231 * t315 + t233 * t313) * t328 + (-t252 * t315 + t253 * t313) * t375;
t268 = -Icges(11,6) * t406 + t400 * t440;
t269 = Icges(11,6) * t400 + t406 * t440;
t270 = -Icges(11,5) * t406 + t400 * t441;
t271 = Icges(11,5) * t400 + t406 * t441;
t298 = Icges(11,2) * t379 + t494;
t299 = Icges(11,1) * t378 + t493;
t413 = (-t269 * t378 + t271 * t379) * t329 + (-t268 * t378 + t270 * t379) * t328 + (-t298 * t378 + t299 * t379) * t375;
t366 = pkin(1) * t482;
t360 = pkin(9) * t407 + pkin(11) * t401;
t359 = -pkin(9) * t401 + pkin(11) * t407;
t354 = rSges(5,1) * t407 - rSges(5,2) * t401;
t352 = rSges(5,1) * t401 + rSges(5,2) * t407;
t350 = rSges(2,1) * t406 - rSges(2,2) * t400;
t349 = rSges(3,1) * t405 - rSges(3,2) * t399;
t348 = rSges(10,1) * t405 - rSges(10,2) * t399;
t347 = rSges(6,1) * t397 + rSges(6,2) * t403;
t346 = rSges(2,1) * t400 + rSges(2,2) * t406;
t334 = Icges(2,5) * t406 - Icges(2,6) * t400;
t333 = Icges(2,5) * t400 + Icges(2,6) * t406;
t322 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t321 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t320 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t318 = t546 * V_base(6);
t300 = rSges(11,1) * t378 + rSges(11,2) * t379;
t274 = rSges(11,3) * t400 + t406 * t463;
t273 = -rSges(11,3) * t406 + t400 * t463;
t265 = (rSges(4,1) * t399 + rSges(4,2) * t405) * t404 + t398 * (rSges(4,1) * t405 - rSges(4,2) * t399);
t264 = -t351 * t408 - t402 * t353;
t263 = (rSges(4,1) * t404 - rSges(4,2) * t398) * t405 - t399 * (rSges(4,1) * t398 + rSges(4,2) * t404);
t262 = V_base(5) * rSges(2,3) - t346 * t375 + t489;
t261 = t350 * t375 + V_base(2) + (-pkin(13) - rSges(2,3)) * V_base(4);
t258 = t355 * t525 + t356 * t390;
t255 = -t324 * t404 - t325 * t398;
t254 = t346 * V_base(4) - t350 * V_base(5) + V_base(3);
t239 = rSges(4,3) * t400 + t263 * t406;
t238 = -rSges(4,3) * t406 + t263 * t400;
t220 = t330 + (t394 * t434 - t484) * pkin(4) + t546;
t219 = qJD(4) * t236 + t375;
t218 = t406 * t491 + V_base(4);
t217 = t400 * t491 + V_base(5);
t214 = t397 * t400 - t406 * t430;
t213 = t400 * t403 + t406 * t431;
t212 = -t397 * t406 - t400 * t430;
t211 = t400 * t431 - t403 * t406;
t210 = t362 * t349 + (t400 * t526 + t527) * t375 + t489;
t209 = -t363 * t349 + t375 * (-t345 * t406 + t528) + t435;
t208 = t399 * t264 - t405 * t443;
t207 = t264 * t405 + t399 * t443;
t206 = t399 * t255 + t405 * t470;
t205 = t255 * t405 - t399 * t470;
t203 = t259 * t374 + t260 * t373;
t202 = -t345 * qJD(2) + t383 * t400 + V_base(3) + (-t345 * t400 - t527) * V_base(4) + (t406 * t526 - t528) * V_base(5);
t201 = t256 * t374 + t257 * t373;
t200 = rSges(7,3) * t400 + t207 * t406;
t199 = -rSges(7,3) * t406 + t207 * t400;
t195 = rSges(9,3) * t400 + t205 * t406;
t194 = -rSges(9,3) * t406 + t205 * t400;
t192 = -Icges(9,2) * t241 + t513;
t182 = -Icges(9,6) * t406 + (Icges(9,2) * t240 - t514) * t400;
t181 = Icges(9,6) * t400 + (-Icges(9,2) * t243 - t514) * t406;
t169 = (t366 - V_base(4) * rSges(8,3) + (-pkin(15) - t469) * V_base(5)) * t406 + (-V_base(5) * rSges(8,3) + (t469 - t546) * V_base(4)) * t400 + t465;
t168 = t215 * t374 + t216 * t373;
t167 = t375 * (t400 * rSges(8,3) + t406 * t469) - t546 * t499 + (-pkin(13) - t258 - t533) * V_base(4) + t437;
t166 = t546 * t480 + V_base(5) * t258 + t375 * (t406 * rSges(8,3) + t400 * ((rSges(8,1) * t525 + rSges(8,2) * t390) * t407 + t401 * (t390 * rSges(8,1) - rSges(8,2) * t525))) + t418;
t165 = t239 * t375 - t265 * t312 + (-t318 - t490) * t406 - t363 * t533 + t464;
t164 = -t238 * t375 + t265 * t328 + t318 * t400 + t418;
t163 = Icges(6,5) * t236 + (Icges(6,1) * t403 - Icges(6,4) * t397) * t237;
t162 = Icges(6,6) * t236 + (Icges(6,4) * t403 - Icges(6,2) * t397) * t237;
t161 = Icges(6,3) * t236 + (Icges(6,5) * t403 - Icges(6,6) * t397) * t237;
t160 = -t347 * t406 + t400 * t444;
t159 = t347 * t400 + t406 * t444;
t158 = t238 * t312 - t239 * t328 + (-qJD(2) + t317) * t386 + t427;
t157 = pkin(2) * t328 * t243 + t362 * t348 + (t406 * rSges(10,3) + (pkin(2) * t241 - pkin(15) - t466) * t400) * t375 + t489;
t156 = -t363 * t348 + t375 * (rSges(10,3) * t400 + t406 * t466) + (-t241 * t499 - t243 * t312) * pkin(2) + t435;
t155 = Icges(6,1) * t214 + Icges(6,4) * t213 + Icges(6,5) * t504;
t154 = Icges(6,1) * t212 + Icges(6,4) * t211 + Icges(6,5) * t505;
t153 = Icges(6,4) * t214 + Icges(6,2) * t213 + Icges(6,6) * t504;
t152 = Icges(6,4) * t212 + Icges(6,2) * t211 + Icges(6,6) * t505;
t151 = Icges(6,5) * t214 + Icges(6,6) * t213 + Icges(6,3) * t504;
t150 = Icges(6,5) * t212 + Icges(6,6) * t211 + Icges(6,3) * t505;
t149 = -V_base(5) * pkin(16) + t208 * t362 + (pkin(14) * t400 - t199) * t375 + t489;
t148 = -t208 * t363 + V_base(2) + (-pkin(13) + pkin(16)) * V_base(4) + (-pkin(14) * t406 + t200) * t375;
t147 = t206 * t328 + (-t400 * pkin(15) - t194) * t375 + t489;
t146 = t195 * t375 - t206 * t312 + t435;
t145 = pkin(14) * t317 + t199 * t363 - t200 * t362 + V_base(3);
t144 = t194 * t312 - t195 * t328 + t427;
t143 = (-t228 * V_base(5) + t366) * t406 + t329 * t273 - t328 * t274 + t475 * t481 + t465;
t142 = (-rSges(5,3) * V_base(4) + t420) * t406 + (-V_base(5) * rSges(5,3) + t483) * t400 + ((t352 * t392 + t354 * t395) * t374 + (t352 * t395 - t354 * t392) * t373) * t317 + t472;
t141 = t220 * t492 - t475 * t480 + (t294 * t391 - t473) * V_base(5) + t328 * t300 - t375 * t273 + t428;
t140 = -pkin(1) * t479 + (-t434 * t532 + t473) * V_base(4) - t329 * t300 + t375 * t274 + t437 + (-qJD(1) * t220 + t475 * V_base(6)) * t406;
t139 = (t547 + t466) * qJD(2) + qJD(3) * t547 + V_base(3) + (-rSges(10,3) * V_base(4) - V_base(5) * pkin(15)) * t406 + (-rSges(10,3) * V_base(5) + t383) * t400 + (t399 * (pkin(2) * t302 + rSges(10,1)) + t405 * (pkin(2) * t442 + rSges(10,2))) * t317;
t138 = V_base(2) + t419 * t400 + (-t432 - t468) * t406 + (-t201 + t429) * V_base(4);
t137 = t201 * V_base(5) + t400 * t432 + t406 * t419 + t433;
t136 = t420 * t406 + t483 * t400 + t218 * t160 - t217 * t159 + ((-t359 * t392 + t360 * t395) * t374 - (t359 * t395 + t360 * t392) * t373) * t317 + t472;
t135 = t159 * t219 - t168 * t218 + V_base(2) + t436 * t400 + (-t468 - t506) * t406 + (-t203 + t429) * V_base(4);
t134 = -t160 * t219 + t168 * t217 + t203 * V_base(5) + t400 * t506 + t406 * t436 + t433;
t1 = (t544 * t400 - t543 * t406) * t362 / 0.2e1 + (t543 * t400 + t544 * t406) * t363 / 0.2e1 + (Icges(1,1) * V_base(4) + t334 * t375 - t540 * t400 + t539 * t406) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t333 * t375 + t539 * t400 + t540 * t406) * V_base(5) / 0.2e1 + t219 * ((t150 * t217 + t151 * t218 + t161 * t219) * t236 + ((-t153 * t397 + t155 * t403) * t218 + (-t152 * t397 + t154 * t403) * t217 + (-t162 * t397 + t163 * t403) * t219) * t237) / 0.2e1 + ((t269 * t379 + t271 * t378) * t329 + (-t181 * t241 + t183 * t243 + t232 * t313 + t234 * t315) * t312 + (-t172 * t236 + t174 * t237 - t223 * t434 + t225 * t308 + t333) * V_base(5) + (-t173 * t236 + t175 * t237 - t224 * t434 + t226 * t308 + t334) * V_base(4) + (t188 * t247 + t190 * t246 + t559 * t399 + t557 * t405) * t363 + (t187 * t247 + t189 * t246 + t560 * t399 + t558 * t405) * t362 + (-t177 * t236 + t178 * t237 - t192 * t241 + t193 * t243 + t197 * t247 + t198 * t246 - t249 * t434 + t250 * t308 + t252 * t313 + t253 * t315 + t298 * t379 + t299 * t378 + t556 * t399 + t555 * t405 + Icges(2,3)) * t375 + (-t182 * t241 + t184 * t243 + t231 * t313 + t233 * t315 + t268 * t379 + t270 * t378) * t328) * t375 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (((-t181 * t243 - t509) * t312 + (-t182 * t243 - t508) * t328 + (-t192 * t243 - t507) * t375 + t414) * t406 + t545 * t400) * t312 / 0.2e1 + V_base(5) * V_base(4) * Icges(1,4) + m(8) * (t166 ^ 2 + t167 ^ 2 + t169 ^ 2) / 0.2e1 + m(10) * (t139 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(4) * (t158 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t142 ^ 2) / 0.2e1 + m(11) * (t140 ^ 2 + t141 ^ 2 + t143 ^ 2) / 0.2e1 + m(9) * (t144 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(7) * (t145 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t218 * ((t151 * t504 + t213 * t153 + t214 * t155) * t218 + (t150 * t504 + t152 * t213 + t154 * t214) * t217 + (t161 * t504 + t162 * t213 + t163 * t214) * t219) / 0.2e1 + t217 * ((t151 * t505 + t153 * t211 + t155 * t212) * t218 + (t150 * t505 + t211 * t152 + t212 * t154) * t217 + (t161 * t505 + t162 * t211 + t163 * t212) * t219) / 0.2e1 + m(1) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + t329 * (t423 * t400 + t413 * t406) / 0.2e1 + m(2) * (t254 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + m(3) * (t202 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + ((-t545 - t423) * t406 + ((t181 * t240 - t509) * t312 + (t182 * t240 - t508) * t328 + (t192 * t240 - t507) * t375 + t414 + t413) * t400) * t328 / 0.2e1;
T = t1;
