% Calculate kinetic energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(6,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: qJD has to be [13x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:22
% EndTime: 2020-05-02 21:17:30
% DurationCPUTime: 7.75s
% Computational Cost: add. (2521->569), mult. (2585->849), div. (0->0), fcn. (2188->27), ass. (0->312)
t344 = V_base(6) * pkin(15);
t353 = sin(qJ(3));
t325 = pkin(5) * t353 + pkin(1);
t354 = sin(qJ(2));
t359 = cos(qJ(2));
t358 = cos(qJ(3));
t478 = pkin(5) * t358;
t488 = -t325 * t354 + t359 * t478;
t489 = t488 * V_base(6) + t344 + qJD(1) * (pkin(15) + t488);
t348 = qJ(2) + qJ(8);
t342 = qJ(9) + t348;
t321 = sin(t342);
t323 = cos(t342);
t419 = -rSges(10,1) * t321 - rSges(10,2) * t323;
t335 = sin(t348);
t481 = pkin(2) * t335;
t487 = pkin(15) - t419 - t481;
t329 = V_base(6) + qJD(1);
t355 = sin(qJ(1));
t453 = t329 * t355;
t360 = cos(qJ(1));
t452 = t329 * t360;
t334 = qJD(2) * t355;
t310 = V_base(4) + t334;
t351 = sin(qJ(6));
t356 = cos(qJ(6));
t422 = rSges(7,1) * t356 - rSges(7,2) * t351;
t390 = -pkin(14) + t422;
t485 = -t360 * rSges(7,3) + t355 * t390;
t484 = t355 * rSges(7,3) + t360 * t390;
t482 = pkin(1) * t354;
t338 = cos(t348);
t480 = pkin(2) * t338;
t347 = -qJ(7) + pkin(19);
t326 = sin(t347);
t327 = cos(t347);
t479 = pkin(4) * (t326 * t354 + t327 * t359);
t477 = pkin(1) * qJD(2);
t476 = t355 * rSges(3,3);
t474 = t360 * rSges(3,3);
t298 = rSges(3,1) * t354 + rSges(3,2) * t359;
t472 = -pkin(15) + t298;
t471 = Icges(2,4) * t355;
t470 = Icges(3,4) * t354;
t469 = Icges(3,4) * t359;
t350 = qJ(2) + qJ(3);
t337 = sin(t350);
t468 = Icges(4,4) * t337;
t340 = cos(t350);
t467 = Icges(4,4) * t340;
t343 = qJ(4) + t350;
t322 = sin(t343);
t466 = Icges(5,4) * t322;
t324 = cos(t343);
t465 = Icges(5,4) * t324;
t464 = Icges(7,4) * t351;
t463 = Icges(7,4) * t356;
t349 = qJ(2) + qJ(7);
t336 = sin(t349);
t462 = Icges(8,4) * t336;
t339 = cos(t349);
t461 = Icges(8,4) * t339;
t460 = Icges(9,4) * t335;
t459 = Icges(9,4) * t338;
t458 = Icges(10,4) * t321;
t457 = Icges(10,4) * t323;
t456 = t322 * t355;
t455 = t322 * t360;
t454 = t325 * t359;
t352 = sin(qJ(5));
t451 = t352 * t360;
t450 = t354 * t358;
t449 = t355 * t352;
t357 = cos(qJ(5));
t448 = t355 * t357;
t447 = t357 * t360;
t328 = -qJ(10) + t347;
t318 = -qJ(2) + t328;
t311 = sin(t318);
t446 = Icges(11,4) * t311;
t312 = cos(t318);
t445 = Icges(11,4) * t312;
t333 = qJD(3) * t355;
t443 = qJD(5) * t322;
t319 = -pkin(15) + t482;
t442 = t319 * qJD(1);
t441 = -qJD(2) - qJD(3);
t440 = -qJD(2) - qJD(8);
t439 = qJD(2) + qJD(7);
t438 = V_base(5) * pkin(13) + V_base(1);
t434 = t359 * t477;
t433 = V_base(5) * pkin(15);
t432 = t354 * V_base(4);
t431 = t354 * V_base(5);
t430 = t355 * V_base(4);
t429 = t359 * V_base(5);
t320 = V_base(5) * t360;
t314 = pkin(1) * t431;
t272 = qJD(8) * t355 + t310;
t273 = qJD(7) * t355 + t310;
t274 = t333 + t310;
t428 = pkin(1) * t429 + t355 * t442 + t438;
t427 = t314 - t433;
t426 = pkin(9) * t324 + pkin(11) * t322;
t425 = -qJD(2) - t430;
t424 = rSges(4,1) * t340 - rSges(4,2) * t337;
t423 = rSges(5,1) * t324 - rSges(5,2) * t322;
t301 = rSges(6,1) * t357 - rSges(6,2) * t352;
t421 = -rSges(8,1) * t336 - rSges(8,2) * t339;
t420 = rSges(9,1) * t335 + rSges(9,2) * t338;
t418 = -t354 * t477 + V_base(3);
t417 = rSges(6,3) * t322 + t301 * t324;
t235 = qJD(4) * t355 + t274;
t416 = -V_base(4) * pkin(13) + V_base(2);
t415 = -Icges(3,1) * t354 - t469;
t414 = Icges(4,1) * t340 - t468;
t413 = Icges(5,1) * t324 - t466;
t412 = Icges(7,1) * t356 - t464;
t411 = -Icges(8,1) * t336 - t461;
t410 = -Icges(9,1) * t335 - t459;
t409 = Icges(10,1) * t321 + t457;
t408 = -Icges(3,2) * t359 - t470;
t407 = -Icges(4,2) * t337 + t467;
t406 = -Icges(5,2) * t322 + t465;
t405 = -Icges(7,2) * t351 + t463;
t404 = -Icges(8,2) * t339 - t462;
t403 = -Icges(9,2) * t338 - t460;
t402 = Icges(10,2) * t323 + t458;
t401 = -Icges(3,5) * t354 - Icges(3,6) * t359;
t400 = Icges(4,5) * t340 - Icges(4,6) * t337;
t399 = Icges(5,5) * t324 - Icges(5,6) * t322;
t398 = Icges(7,5) * t356 - Icges(7,6) * t351;
t397 = -Icges(8,5) * t336 - Icges(8,6) * t339;
t396 = -Icges(9,5) * t335 - Icges(9,6) * t338;
t395 = Icges(10,5) * t321 + Icges(10,6) * t323;
t346 = V_base(4) * pkin(15);
t394 = -pkin(1) * t432 + t346;
t393 = -Icges(11,1) * t311 + t445;
t392 = Icges(11,2) * t312 - t446;
t391 = -Icges(11,5) * t311 + Icges(11,6) * t312;
t389 = -t320 - t425;
t388 = pkin(15) * t452 + t416;
t387 = -qJD(7) - t389;
t249 = -qJD(8) - t389;
t283 = rSges(11,1) * t354 + rSges(11,2) * t359;
t284 = rSges(11,1) * t359 - rSges(11,2) * t354;
t316 = sin(t328);
t317 = cos(t328);
t385 = pkin(4) * (t326 * t359 - t327 * t354) + t283 * t317 - t284 * t316;
t233 = V_base(5) + (-qJD(4) + t441) * t360;
t229 = -V_base(5) + (qJD(10) + t439) * t360;
t230 = qJD(10) * t355 + t273;
t384 = (-Icges(11,3) * t360 + t355 * t391) * t229 - (Icges(11,3) * t355 + t360 * t391) * t230 - (-Icges(11,5) * t312 - Icges(11,6) * t311) * t329;
t232 = V_base(5) + (-qJD(9) + t440) * t360;
t234 = qJD(9) * t355 + t272;
t383 = (-Icges(10,3) * t360 + t355 * t395) * t232 + (Icges(10,3) * t355 + t360 * t395) * t234 + (-Icges(10,5) * t323 + Icges(10,6) * t321) * t329;
t382 = (-Icges(5,3) * t360 + t355 * t399) * t233 + (Icges(5,3) * t355 + t360 * t399) * t235 + (Icges(5,5) * t322 + Icges(5,6) * t324) * t329;
t269 = t360 * t440 + V_base(5);
t381 = (-Icges(9,3) * t360 + t355 * t396) * t269 + (Icges(9,3) * t355 + t360 * t396) * t272 + (Icges(9,5) * t338 - Icges(9,6) * t335) * t329;
t270 = -t360 * t439 + V_base(5);
t380 = (-Icges(8,3) * t360 + t355 * t397) * t270 + (Icges(8,3) * t355 + t360 * t397) * t273 + (Icges(8,5) * t339 - Icges(8,6) * t336) * t329;
t271 = t360 * t441 + V_base(5);
t379 = (-Icges(4,3) * t360 + t355 * t400) * t271 + (Icges(4,3) * t355 + t360 * t400) * t274 + (Icges(4,5) * t337 + Icges(4,6) * t340) * t329;
t307 = -qJD(6) * t360 + V_base(5);
t309 = qJD(6) * t355 + V_base(4);
t378 = (-Icges(7,3) * t360 + t355 * t398) * t307 + (Icges(7,3) * t355 + t360 * t398) * t309 + (Icges(7,5) * t351 + Icges(7,6) * t356) * t329;
t308 = -qJD(2) * t360 + V_base(5);
t377 = (-Icges(3,3) * t360 + t355 * t401) * t308 + (Icges(3,3) * t355 + t360 * t401) * t310 + (Icges(3,5) * t359 - Icges(3,6) * t354) * t329;
t376 = pkin(5) * (qJD(3) + t389) * t340 + t427 * t360 + t394 * t355 + t418;
t277 = t482 * V_base(6) - t344;
t374 = t277 * t355 - t360 * t434 + t428;
t373 = -pkin(1) * t310 * t359 + t416;
t372 = t360 * t314 + t425 * t482 + V_base(3) + (-t320 + t430) * pkin(15);
t371 = t373 + (-t277 - t442) * t360;
t240 = pkin(5) * t450 + t454;
t276 = t353 * t359 + t450;
t370 = t325 * t429 + t431 * t478 + (-pkin(5) * qJD(3) * t276 - qJD(2) * t240) * t360 + t438 - t489 * t355;
t157 = -Icges(11,6) * t360 + t355 * t392;
t158 = Icges(11,6) * t355 + t360 * t392;
t159 = -Icges(11,5) * t360 + t355 * t393;
t160 = Icges(11,5) * t355 + t360 * t393;
t211 = -Icges(11,2) * t311 - t445;
t212 = -Icges(11,1) * t312 - t446;
t369 = (t158 * t312 - t160 * t311) * t230 - (t157 * t312 - t159 * t311) * t229 + (t211 * t312 - t212 * t311) * t329;
t174 = -Icges(10,6) * t360 + t355 * t402;
t175 = Icges(10,6) * t355 + t360 * t402;
t178 = -Icges(10,5) * t360 + t355 * t409;
t179 = Icges(10,5) * t355 + t360 * t409;
t243 = Icges(10,2) * t321 - t457;
t245 = -Icges(10,1) * t323 + t458;
t368 = (t175 * t323 + t179 * t321) * t234 + (t174 * t323 + t178 * t321) * t232 + (t243 * t323 + t245 * t321) * t329;
t176 = -Icges(5,6) * t360 + t355 * t406;
t177 = Icges(5,6) * t355 + t360 * t406;
t180 = -Icges(5,5) * t360 + t355 * t413;
t181 = Icges(5,5) * t355 + t360 * t413;
t244 = Icges(5,2) * t324 + t466;
t246 = Icges(5,1) * t322 + t465;
t367 = (-t177 * t322 + t181 * t324) * t235 + (-t176 * t322 + t180 * t324) * t233 + (-t244 * t322 + t246 * t324) * t329;
t190 = -Icges(9,6) * t360 + t355 * t403;
t191 = Icges(9,6) * t355 + t360 * t403;
t196 = -Icges(9,5) * t360 + t355 * t410;
t197 = Icges(9,5) * t355 + t360 * t410;
t260 = -Icges(9,2) * t335 + t459;
t263 = Icges(9,1) * t338 - t460;
t366 = (-t191 * t338 - t197 * t335) * t272 + (-t190 * t338 - t196 * t335) * t269 + (-t260 * t338 - t263 * t335) * t329;
t192 = -Icges(8,6) * t360 + t355 * t404;
t193 = Icges(8,6) * t355 + t360 * t404;
t198 = -Icges(8,5) * t360 + t355 * t411;
t199 = Icges(8,5) * t355 + t360 * t411;
t261 = -Icges(8,2) * t336 + t461;
t264 = Icges(8,1) * t339 - t462;
t365 = (-t193 * t339 - t199 * t336) * t273 + (-t192 * t339 - t198 * t336) * t270 + (-t261 * t339 - t264 * t336) * t329;
t194 = -Icges(4,6) * t360 + t355 * t407;
t195 = Icges(4,6) * t355 + t360 * t407;
t200 = -Icges(4,5) * t360 + t355 * t414;
t201 = Icges(4,5) * t355 + t360 * t414;
t262 = Icges(4,2) * t340 + t468;
t265 = Icges(4,1) * t337 + t467;
t364 = (-t195 * t337 + t201 * t340) * t274 + (-t194 * t337 + t200 * t340) * t271 + (-t262 * t337 + t265 * t340) * t329;
t217 = -Icges(7,6) * t360 + t355 * t405;
t218 = Icges(7,6) * t355 + t360 * t405;
t221 = -Icges(7,5) * t360 + t355 * t412;
t222 = Icges(7,5) * t355 + t360 * t412;
t289 = Icges(7,2) * t356 + t464;
t293 = Icges(7,1) * t351 + t463;
t363 = (-t218 * t351 + t222 * t356) * t309 + (-t217 * t351 + t221 * t356) * t307 + (-t289 * t351 + t293 * t356) * t329;
t219 = -Icges(3,6) * t360 + t355 * t408;
t220 = Icges(3,6) * t355 + t360 * t408;
t223 = -Icges(3,5) * t360 + t355 * t415;
t224 = Icges(3,5) * t355 + t360 * t415;
t290 = -Icges(3,2) * t354 + t469;
t294 = Icges(3,1) * t359 - t470;
t362 = (-t220 * t359 - t224 * t354) * t310 + (-t219 * t359 - t223 * t354) * t308 + (-t290 * t359 - t294 * t354) * t329;
t361 = -V_base(4) * t454 + (-t276 * t333 - t358 * t432) * pkin(5) - t240 * t334 + t416 + t489 * t360;
t341 = Icges(2,4) * t360;
t303 = rSges(2,1) * t360 - rSges(2,2) * t355;
t302 = rSges(3,1) * t359 - rSges(3,2) * t354;
t300 = rSges(6,1) * t352 + rSges(6,2) * t357;
t299 = rSges(2,1) * t355 + rSges(2,2) * t360;
t297 = rSges(7,1) * t351 + rSges(7,2) * t356;
t296 = Icges(2,1) * t360 - t471;
t295 = Icges(2,1) * t355 + t341;
t292 = -Icges(2,2) * t355 + t341;
t291 = Icges(2,2) * t360 + t471;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t275 = -qJD(5) * t324 + t329;
t268 = rSges(8,1) * t339 - rSges(8,2) * t336;
t267 = rSges(9,1) * t338 - rSges(9,2) * t335;
t266 = rSges(4,1) * t337 + rSges(4,2) * t340;
t250 = pkin(9) * t322 - pkin(11) * t324;
t248 = rSges(10,1) * t323 - rSges(10,2) * t321;
t247 = rSges(5,1) * t322 + rSges(5,2) * t324;
t228 = t324 * t447 + t449;
t227 = -t324 * t451 + t448;
t226 = t324 * t448 - t451;
t225 = -t324 * t449 - t447;
t209 = t426 * t360;
t208 = t426 * t355;
t205 = rSges(4,3) * t355 + t360 * t424;
t204 = rSges(8,3) * t355 + t360 * t421;
t203 = -rSges(4,3) * t360 + t355 * t424;
t202 = -rSges(8,3) * t360 + t355 * t421;
t183 = rSges(5,3) * t355 + t360 * t423;
t182 = -rSges(5,3) * t360 + t355 * t423;
t169 = V_base(5) * rSges(2,3) - t299 * t329 + t438;
t168 = t303 * t329 + V_base(2) + (-pkin(13) - rSges(2,3)) * V_base(4);
t167 = t360 * t443 + t235;
t166 = t355 * t443 + t233;
t165 = t299 * V_base(4) - t303 * V_base(5) + V_base(3);
t164 = -rSges(6,3) * t324 + t301 * t322;
t163 = -Icges(6,5) * t324 + (Icges(6,1) * t357 - Icges(6,4) * t352) * t322;
t162 = -Icges(6,6) * t324 + (Icges(6,4) * t357 - Icges(6,2) * t352) * t322;
t161 = -Icges(6,3) * t324 + (Icges(6,5) * t357 - Icges(6,6) * t352) * t322;
t152 = t283 * t316 + t284 * t317;
t151 = t300 * t355 + t360 * t417;
t150 = -t300 * t360 + t355 * t417;
t149 = Icges(6,1) * t228 + Icges(6,4) * t227 + Icges(6,5) * t455;
t148 = Icges(6,1) * t226 + Icges(6,4) * t225 + Icges(6,5) * t456;
t147 = Icges(6,4) * t228 + Icges(6,2) * t227 + Icges(6,6) * t455;
t146 = Icges(6,4) * t226 + Icges(6,2) * t225 + Icges(6,6) * t456;
t145 = Icges(6,5) * t228 + Icges(6,6) * t227 + Icges(6,3) * t455;
t144 = Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t456;
t143 = t308 * t302 + (t355 * t472 + t474) * t329 + t438;
t142 = -t310 * t302 + t329 * (-t298 * t360 + t476) + t388;
t141 = -V_base(5) * pkin(16) + t307 * t297 - t329 * t485 + t438;
t140 = -t309 * t297 + V_base(2) + (-pkin(13) + pkin(16)) * V_base(4) + t484 * t329;
t139 = -t298 * qJD(2) + t346 * t355 + V_base(3) + (-t298 * t355 - t474) * V_base(4) + (t360 * t472 - t476) * V_base(5);
t138 = t422 * qJD(6) - t484 * V_base(5) + t485 * V_base(4) + V_base(3);
t137 = t269 * t267 + (t360 * rSges(9,3) + (-pkin(15) + t420) * t355) * t329 + t438;
t136 = -t272 * t267 + t329 * (rSges(9,3) * t355 - t360 * t420) + t388;
t135 = V_base(3) + (-rSges(9,3) * V_base(4) - t433) * t360 + (-rSges(9,3) * V_base(5) + t346) * t355 + t420 * t249;
t134 = t249 * t481 + V_base(3) + (-rSges(10,3) * V_base(4) - t433) * t360 + (-rSges(10,3) * V_base(5) + t346) * t355 + t419 * (-qJD(9) + t249);
t133 = t205 * t329 - t266 * t274 + t371;
t132 = t329 * t204 - t273 * t268 + t371;
t131 = -t203 * t329 + t266 * t271 + t374;
t130 = -t202 * t329 + t268 * t270 + t374;
t129 = t203 * t274 - t205 * t271 + t372;
t128 = t202 * t273 - t204 * t270 + t372;
t127 = rSges(10,3) * t453 + t234 * t248 - t272 * t480 + t452 * t487 + t416;
t126 = rSges(10,3) * t452 - t232 * t248 + t269 * t480 - t453 * t487 + t438;
t125 = pkin(4) * t387 * sin(qJ(2) - t347) + (-V_base(4) * rSges(11,3) + t427) * t360 + (-V_base(5) * rSges(11,3) + t394) * t355 + (rSges(11,1) * t311 - rSges(11,2) * t312) * (-qJD(10) + t387) + t418;
t124 = t182 * t235 - t183 * t233 + t376;
t123 = t183 * t329 - t235 * t247 + t361;
t122 = -t182 * t329 + t233 * t247 + t370;
t121 = rSges(11,3) * t453 + t152 * t230 - t273 * t479 + t373 + (-t319 + t385) * t452;
t120 = t270 * t479 + t152 * t229 + (rSges(11,3) * t329 - t434) * t360 + (V_base(6) * t319 - t329 * t385) * t355 + t428;
t119 = t150 * t167 - t151 * t166 + t208 * t235 - t209 * t233 + t376;
t118 = t151 * t275 - t164 * t167 + t209 * t329 - t235 * t250 + t361;
t117 = -t275 * t150 + t166 * t164 - t329 * t208 + t233 * t250 + t370;
t1 = V_base(5) * t329 * (Icges(2,5) * t355 + Icges(2,6) * t360) + t329 * V_base(4) * (Icges(2,5) * t360 - Icges(2,6) * t355) + ((t291 * t360 + t355 * t295 + Icges(1,2)) * V_base(5) + (t292 * t360 + t355 * t296 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t355 * t291 + t295 * t360 + Icges(1,4)) * V_base(5) + (-t355 * t292 + t296 * t360 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(2) * (t165 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(9) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(7) * (t138 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t139 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(5) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(11) * (t120 ^ 2 + t121 ^ 2 + t125 ^ 2) / 0.2e1 + m(8) * (t128 ^ 2 + t130 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t129 ^ 2 + t131 ^ 2 + t133 ^ 2) / 0.2e1 + m(10) * (t126 ^ 2 + t127 ^ 2 + t134 ^ 2) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + ((t218 * t356 + t222 * t351) * t309 + (t217 * t356 + t221 * t351) * t307 + (-t220 * t354 + t224 * t359) * t310 + (-t219 * t354 + t223 * t359) * t308 + (t195 * t340 + t201 * t337) * t274 + (t194 * t340 + t200 * t337) * t271 + (-t191 * t335 + t197 * t338) * t272 + (-t190 * t335 + t196 * t338) * t269 + (-t193 * t336 + t199 * t339) * t273 + (-t192 * t336 + t198 * t339) * t270 + (t177 * t324 + t181 * t322) * t235 + (t176 * t324 + t180 * t322) * t233 + (t175 * t321 - t179 * t323) * t234 + (t174 * t321 - t178 * t323) * t232 + (-t158 * t311 - t160 * t312) * t230 - (-t157 * t311 - t159 * t312) * t229 + (Icges(2,3) + t289 * t356 + t293 * t351 - t290 * t354 + t294 * t359 + t262 * t340 + t265 * t337 - t260 * t335 + t263 * t338 - t261 * t336 + t264 * t339 - t211 * t311 - t212 * t312 + t244 * t324 + t246 * t322 + t243 * t321 - t245 * t323) * t329) * t329 / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t275 * ((-t144 * t166 - t145 * t167 - t161 * t275) * t324 + ((-t147 * t352 + t149 * t357) * t167 + (-t146 * t352 + t148 * t357) * t166 + (-t162 * t352 + t163 * t357) * t275) * t322) / 0.2e1 + t166 * ((t145 * t456 + t147 * t225 + t149 * t226) * t167 + (t144 * t456 + t146 * t225 + t148 * t226) * t166 + (t161 * t456 + t162 * t225 + t163 * t226) * t275) / 0.2e1 + t167 * ((t145 * t455 + t227 * t147 + t228 * t149) * t167 + (t144 * t455 + t146 * t227 + t148 * t228) * t166 + (t161 * t455 + t162 * t227 + t163 * t228) * t275) / 0.2e1 + t269 * (t366 * t355 - t381 * t360) / 0.2e1 + t235 * (t382 * t355 + t367 * t360) / 0.2e1 + t233 * (t367 * t355 - t382 * t360) / 0.2e1 + t234 * (t383 * t355 + t368 * t360) / 0.2e1 + t232 * (t368 * t355 - t383 * t360) / 0.2e1 - t229 * (t369 * t355 + t384 * t360) / 0.2e1 + t230 * (-t384 * t355 + t369 * t360) / 0.2e1 + t310 * (t377 * t355 + t362 * t360) / 0.2e1 + t308 * (t362 * t355 - t377 * t360) / 0.2e1 + t309 * (t378 * t355 + t363 * t360) / 0.2e1 + t307 * (t363 * t355 - t378 * t360) / 0.2e1 + t274 * (t379 * t355 + t364 * t360) / 0.2e1 + t271 * (t364 * t355 - t379 * t360) / 0.2e1 + t273 * (t380 * t355 + t365 * t360) / 0.2e1 + t270 * (t365 * t355 - t380 * t360) / 0.2e1 + t272 * (t381 * t355 + t366 * t360) / 0.2e1;
T = t1;
