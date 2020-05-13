% Calculate kinetic energy for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(6,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: qJD has to be [13x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m1OL_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:19
% EndTime: 2020-04-15 19:28:27
% DurationCPUTime: 6.94s
% Computational Cost: add. (2639->528), mult. (2582->800), div. (0->0), fcn. (2234->22), ass. (0->284)
t347 = sin(qJ(1));
t351 = cos(qJ(1));
t453 = -t347 * V_base(4) + V_base(5) * t351;
t350 = cos(qJ(2));
t451 = pkin(1) * t350;
t341 = qJ(2) + qJ(8);
t334 = cos(t341);
t450 = pkin(2) * t334;
t342 = qJ(2) + qJ(7);
t324 = pkin(19) - t342;
t449 = pkin(4) * cos(t324);
t343 = qJ(2) + qJ(3);
t333 = sin(t343);
t448 = pkin(5) * t333;
t447 = pkin(15) * t347;
t331 = sin(t341);
t446 = pkin(2) * t331;
t346 = sin(qJ(2));
t445 = t346 * pkin(1);
t444 = Icges(2,4) * t347;
t443 = Icges(3,4) * t346;
t442 = Icges(3,4) * t350;
t441 = Icges(4,4) * t333;
t336 = cos(t343);
t440 = Icges(4,4) * t336;
t339 = qJ(4) + t343;
t321 = sin(t339);
t439 = Icges(5,4) * t321;
t323 = cos(t339);
t438 = Icges(5,4) * t323;
t344 = sin(qJ(6));
t437 = Icges(7,4) * t344;
t348 = cos(qJ(6));
t436 = Icges(7,4) * t348;
t332 = sin(t342);
t435 = Icges(8,4) * t332;
t335 = cos(t342);
t434 = Icges(8,4) * t335;
t433 = Icges(9,4) * t331;
t432 = Icges(9,4) * t334;
t338 = qJ(9) + t341;
t320 = sin(t338);
t431 = Icges(10,4) * t320;
t322 = cos(t338);
t430 = Icges(10,4) * t322;
t429 = t321 * t347;
t428 = t321 * t351;
t345 = sin(qJ(5));
t427 = t345 * t351;
t426 = t347 * t345;
t349 = cos(qJ(5));
t425 = t347 * t349;
t424 = t349 * t351;
t423 = pkin(4) * sin(t324);
t422 = pkin(5) * t336;
t318 = -qJ(10) + t324;
t314 = sin(t318);
t421 = Icges(11,4) * t314;
t315 = cos(t318);
t420 = Icges(11,4) * t315;
t419 = qJD(5) * t321;
t418 = -qJD(2) - qJD(3);
t417 = -qJD(2) - qJD(7);
t416 = -qJD(2) - qJD(8);
t415 = V_base(5) * pkin(13) + V_base(1);
t312 = qJD(2) * t347 + V_base(4);
t325 = V_base(6) + qJD(1);
t265 = t445 * t347;
t410 = t265 - t447;
t310 = -qJD(2) * t351 + V_base(5);
t409 = t310 * t451 + t415;
t283 = qJD(8) * t347 + t312;
t284 = qJD(7) * t347 + t312;
t285 = qJD(3) * t347 + t312;
t193 = t422 * t347;
t408 = -t193 + t410;
t282 = t351 * t418 + V_base(5);
t407 = t282 * t448 + t409;
t406 = pkin(9) * t323 + pkin(11) * t321;
t405 = -rSges(3,1) * t346 - rSges(3,2) * t350;
t404 = rSges(4,1) * t336 - rSges(4,2) * t333;
t403 = rSges(5,1) * t323 - rSges(5,2) * t321;
t402 = rSges(7,1) * t348 - rSges(7,2) * t344;
t401 = -rSges(8,1) * t332 - rSges(8,2) * t335;
t400 = -rSges(9,1) * t331 - rSges(9,2) * t334;
t399 = rSges(10,1) * t320 + rSges(10,2) * t322;
t253 = qJD(4) * t347 + t285;
t398 = -rSges(11,1) * t314 + rSges(11,2) * t315;
t397 = -Icges(3,1) * t346 - t442;
t396 = Icges(4,1) * t336 - t441;
t395 = Icges(5,1) * t323 - t439;
t394 = Icges(7,1) * t348 - t437;
t393 = -Icges(8,1) * t332 - t434;
t392 = -Icges(9,1) * t331 - t432;
t391 = Icges(10,1) * t320 + t430;
t390 = -Icges(3,2) * t350 - t443;
t389 = -Icges(4,2) * t333 + t440;
t388 = -Icges(5,2) * t321 + t438;
t387 = -Icges(7,2) * t344 + t436;
t386 = -Icges(8,2) * t335 - t435;
t385 = -Icges(9,2) * t334 - t433;
t384 = Icges(10,2) * t322 + t431;
t383 = -Icges(3,5) * t346 - Icges(3,6) * t350;
t382 = Icges(4,5) * t336 - Icges(4,6) * t333;
t381 = Icges(5,5) * t323 - Icges(5,6) * t321;
t380 = Icges(7,5) * t348 - Icges(7,6) * t344;
t379 = -Icges(8,5) * t332 - Icges(8,6) * t335;
t378 = -Icges(9,5) * t331 - Icges(9,6) * t334;
t377 = Icges(10,5) * t320 + Icges(10,6) * t322;
t376 = -Icges(11,1) * t314 + t420;
t375 = Icges(11,2) * t315 - t421;
t374 = -Icges(11,5) * t314 + Icges(11,6) * t315;
t373 = t325 * t351 * pkin(15) - V_base(4) * pkin(13) + V_base(2);
t372 = -pkin(15) * t453 + V_base(3);
t251 = V_base(5) + (-qJD(4) + t418) * t351;
t248 = V_base(5) + (-qJD(10) + t417) * t351;
t249 = qJD(10) * t347 + t284;
t371 = (-Icges(11,3) * t351 + t347 * t374) * t248 + (Icges(11,3) * t347 + t351 * t374) * t249 + (-Icges(11,5) * t315 - Icges(11,6) * t314) * t325;
t250 = V_base(5) + (-qJD(9) + t416) * t351;
t252 = qJD(9) * t347 + t283;
t370 = (-Icges(10,3) * t351 + t347 * t377) * t250 + (Icges(10,3) * t347 + t351 * t377) * t252 + (-Icges(10,5) * t322 + Icges(10,6) * t320) * t325;
t369 = (-Icges(5,3) * t351 + t347 * t381) * t251 + (Icges(5,3) * t347 + t351 * t381) * t253 + (Icges(5,5) * t321 + Icges(5,6) * t323) * t325;
t280 = t351 * t416 + V_base(5);
t368 = (-Icges(9,3) * t351 + t347 * t378) * t280 + (Icges(9,3) * t347 + t351 * t378) * t283 + (Icges(9,5) * t334 - Icges(9,6) * t331) * t325;
t281 = t351 * t417 + V_base(5);
t367 = (-Icges(8,3) * t351 + t347 * t379) * t281 + (Icges(8,3) * t347 + t351 * t379) * t284 + (Icges(8,5) * t335 - Icges(8,6) * t332) * t325;
t366 = (-Icges(4,3) * t351 + t347 * t382) * t282 + (Icges(4,3) * t347 + t351 * t382) * t285 + (Icges(4,5) * t333 + Icges(4,6) * t336) * t325;
t309 = -qJD(6) * t351 + V_base(5);
t311 = qJD(6) * t347 + V_base(4);
t365 = (-Icges(7,3) * t351 + t347 * t380) * t309 + (Icges(7,3) * t347 + t351 * t380) * t311 + (Icges(7,5) * t344 + Icges(7,6) * t348) * t325;
t364 = (-Icges(3,3) * t351 + t347 * t383) * t310 + (Icges(3,3) * t347 + t351 * t383) * t312 + (Icges(3,5) * t350 - Icges(3,6) * t346) * t325;
t266 = t445 * t351;
t363 = -t325 * t266 - t312 * t451 + t373;
t362 = -t312 * t265 + t310 * t266 + t372;
t194 = t422 * t351;
t361 = t285 * t193 - t282 * t194 + t362;
t360 = t325 * t194 - t285 * t448 + t363;
t160 = -Icges(11,6) * t351 + t347 * t375;
t161 = Icges(11,6) * t347 + t351 * t375;
t162 = -Icges(11,5) * t351 + t347 * t376;
t163 = Icges(11,5) * t347 + t351 * t376;
t223 = -Icges(11,2) * t314 - t420;
t224 = -Icges(11,1) * t315 - t421;
t359 = (t161 * t315 - t163 * t314) * t249 + (t160 * t315 - t162 * t314) * t248 + (t223 * t315 - t224 * t314) * t325;
t181 = -Icges(10,6) * t351 + t347 * t384;
t182 = Icges(10,6) * t347 + t351 * t384;
t185 = -Icges(10,5) * t351 + t347 * t391;
t186 = Icges(10,5) * t347 + t351 * t391;
t258 = Icges(10,2) * t320 - t430;
t260 = -Icges(10,1) * t322 + t431;
t358 = (t182 * t322 + t186 * t320) * t252 + (t181 * t322 + t185 * t320) * t250 + (t258 * t322 + t260 * t320) * t325;
t183 = -Icges(5,6) * t351 + t347 * t388;
t184 = Icges(5,6) * t347 + t351 * t388;
t187 = -Icges(5,5) * t351 + t347 * t395;
t188 = Icges(5,5) * t347 + t351 * t395;
t259 = Icges(5,2) * t323 + t439;
t261 = Icges(5,1) * t321 + t438;
t357 = (-t184 * t321 + t188 * t323) * t253 + (-t183 * t321 + t187 * t323) * t251 + (-t259 * t321 + t261 * t323) * t325;
t201 = -Icges(9,6) * t351 + t347 * t385;
t202 = Icges(9,6) * t347 + t351 * t385;
t207 = -Icges(9,5) * t351 + t347 * t392;
t208 = Icges(9,5) * t347 + t351 * t392;
t271 = -Icges(9,2) * t331 + t432;
t274 = Icges(9,1) * t334 - t433;
t356 = (-t202 * t334 - t208 * t331) * t283 + (-t201 * t334 - t207 * t331) * t280 + (-t271 * t334 - t274 * t331) * t325;
t203 = -Icges(8,6) * t351 + t347 * t386;
t204 = Icges(8,6) * t347 + t351 * t386;
t209 = -Icges(8,5) * t351 + t347 * t393;
t210 = Icges(8,5) * t347 + t351 * t393;
t272 = -Icges(8,2) * t332 + t434;
t275 = Icges(8,1) * t335 - t435;
t355 = (-t204 * t335 - t210 * t332) * t284 + (-t203 * t335 - t209 * t332) * t281 + (-t272 * t335 - t275 * t332) * t325;
t205 = -Icges(4,6) * t351 + t347 * t389;
t206 = Icges(4,6) * t347 + t351 * t389;
t211 = -Icges(4,5) * t351 + t347 * t396;
t212 = Icges(4,5) * t347 + t351 * t396;
t273 = Icges(4,2) * t336 + t441;
t276 = Icges(4,1) * t333 + t440;
t354 = (-t206 * t333 + t212 * t336) * t285 + (-t205 * t333 + t211 * t336) * t282 + (-t273 * t333 + t276 * t336) * t325;
t230 = -Icges(7,6) * t351 + t347 * t387;
t231 = Icges(7,6) * t347 + t351 * t387;
t234 = -Icges(7,5) * t351 + t347 * t394;
t235 = Icges(7,5) * t347 + t351 * t394;
t296 = Icges(7,2) * t348 + t437;
t300 = Icges(7,1) * t344 + t436;
t353 = (-t231 * t344 + t235 * t348) * t311 + (-t230 * t344 + t234 * t348) * t309 + (-t296 * t344 + t300 * t348) * t325;
t232 = -Icges(3,6) * t351 + t347 * t390;
t233 = Icges(3,6) * t347 + t351 * t390;
t236 = -Icges(3,5) * t351 + t347 * t397;
t237 = Icges(3,5) * t347 + t351 * t397;
t297 = -Icges(3,2) * t346 + t442;
t301 = Icges(3,1) * t350 - t443;
t352 = (-t233 * t350 - t237 * t346) * t312 + (-t232 * t350 - t236 * t346) * t310 + (-t297 * t350 - t301 * t346) * t325;
t337 = Icges(2,4) * t351;
t307 = rSges(2,1) * t351 - t347 * rSges(2,2);
t306 = rSges(3,1) * t350 - rSges(3,2) * t346;
t305 = t347 * rSges(2,1) + rSges(2,2) * t351;
t304 = rSges(7,1) * t344 + rSges(7,2) * t348;
t303 = Icges(2,1) * t351 - t444;
t302 = Icges(2,1) * t347 + t337;
t299 = -Icges(2,2) * t347 + t337;
t298 = Icges(2,2) * t351 + t444;
t291 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t290 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t289 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t286 = -qJD(5) * t323 + t325;
t279 = rSges(8,1) * t335 - rSges(8,2) * t332;
t278 = rSges(9,1) * t334 - rSges(9,2) * t331;
t277 = rSges(4,1) * t333 + rSges(4,2) * t336;
t264 = pkin(9) * t321 - pkin(11) * t323;
t263 = -rSges(10,1) * t322 + rSges(10,2) * t320;
t262 = rSges(5,1) * t321 + rSges(5,2) * t323;
t255 = t446 * t351;
t254 = t446 * t347;
t247 = t323 * t424 + t426;
t246 = -t323 * t427 + t425;
t245 = t323 * t425 - t427;
t244 = -t323 * t426 - t424;
t242 = t347 * rSges(7,3) + t351 * t402;
t241 = t347 * rSges(3,3) + t351 * t405;
t240 = -rSges(7,3) * t351 + t347 * t402;
t239 = -rSges(3,3) * t351 + t347 * t405;
t238 = -rSges(11,1) * t315 - rSges(11,2) * t314;
t221 = t406 * t351;
t220 = t406 * t347;
t219 = t347 * rSges(4,3) + t351 * t404;
t218 = t347 * rSges(8,3) + t351 * t401;
t217 = t347 * rSges(9,3) + t351 * t400;
t216 = -rSges(4,3) * t351 + t347 * t404;
t215 = -rSges(8,3) * t351 + t347 * t401;
t214 = -rSges(9,3) * t351 + t347 * t400;
t192 = t347 * rSges(5,3) + t351 * t403;
t191 = t347 * rSges(10,3) + t351 * t399;
t190 = -rSges(5,3) * t351 + t347 * t403;
t189 = -rSges(10,3) * t351 + t347 * t399;
t176 = V_base(5) * rSges(2,3) - t305 * t325 + t415;
t175 = t307 * t325 + V_base(2) + (-pkin(13) - rSges(2,3)) * V_base(4);
t174 = t423 * t351;
t173 = t423 * t347;
t172 = t351 * t419 + t253;
t171 = t347 * t419 + t251;
t170 = t305 * V_base(4) - t307 * V_base(5) + V_base(3);
t169 = -rSges(6,3) * t323 + (rSges(6,1) * t349 - rSges(6,2) * t345) * t321;
t168 = -Icges(6,5) * t323 + (Icges(6,1) * t349 - Icges(6,4) * t345) * t321;
t167 = -Icges(6,6) * t323 + (Icges(6,4) * t349 - Icges(6,2) * t345) * t321;
t166 = -Icges(6,3) * t323 + (Icges(6,5) * t349 - Icges(6,6) * t345) * t321;
t165 = t347 * rSges(11,3) + t351 * t398;
t164 = -rSges(11,3) * t351 + t347 * t398;
t155 = t247 * rSges(6,1) + t246 * rSges(6,2) + rSges(6,3) * t428;
t154 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t429;
t153 = Icges(6,1) * t247 + Icges(6,4) * t246 + Icges(6,5) * t428;
t152 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t429;
t151 = Icges(6,4) * t247 + Icges(6,2) * t246 + Icges(6,6) * t428;
t150 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t429;
t149 = Icges(6,5) * t247 + Icges(6,6) * t246 + Icges(6,3) * t428;
t148 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t429;
t147 = t306 * t310 + (-t239 - t447) * t325 + t415;
t146 = t241 * t325 - t306 * t312 + t373;
t145 = -V_base(5) * pkin(16) + t304 * t309 + (pkin(14) * t347 - t240) * t325 + t415;
t144 = -t311 * t304 + V_base(2) + (-pkin(13) + pkin(16)) * V_base(4) + (-pkin(14) * t351 + t242) * t325;
t143 = t278 * t280 + (-t214 - t447) * t325 + t415;
t142 = t217 * t325 - t278 * t283 + t373;
t141 = pkin(14) * t453 + t311 * t240 - t309 * t242 + V_base(3);
t140 = t312 * t239 - t310 * t241 + t372;
t139 = t283 * t214 - t280 * t217 + t372;
t138 = t277 * t282 + (-t216 + t410) * t325 + t409;
t137 = t279 * t281 + (-t215 + t410) * t325 + t409;
t136 = t219 * t325 - t277 * t285 + t363;
t135 = t218 * t325 - t279 * t284 + t363;
t134 = t280 * t450 + t250 * t263 + (-t189 + t254 - t447) * t325 + t415;
t133 = -t283 * t450 - t252 * t263 + (t191 - t255) * t325 + t373;
t132 = t285 * t216 - t282 * t219 + t362;
t131 = t284 * t215 - t281 * t218 + t362;
t130 = t251 * t262 + (-t190 + t408) * t325 + t407;
t129 = t192 * t325 - t253 * t262 + t360;
t128 = t252 * t189 - t250 * t191 - t283 * t254 + t280 * t255 + t372;
t127 = t281 * t449 + t238 * t248 + (-t164 - t173 + t410) * t325 + t409;
t126 = -t284 * t449 - t238 * t249 + (t165 + t174) * t325 + t363;
t125 = t253 * t190 - t251 * t192 + t361;
t124 = t249 * t164 - t248 * t165 + t284 * t173 - t281 * t174 + t362;
t123 = -t154 * t286 + t169 * t171 + t251 * t264 + (-t220 + t408) * t325 + t407;
t122 = t155 * t286 - t169 * t172 + t221 * t325 - t253 * t264 + t360;
t121 = t172 * t154 - t171 * t155 + t253 * t220 - t251 * t221 + t361;
t1 = t325 * V_base(5) * (Icges(2,5) * t347 + Icges(2,6) * t351) + V_base(4) * t325 * (Icges(2,5) * t351 - Icges(2,6) * t347) + ((-t347 * t298 + t302 * t351 + Icges(1,4)) * V_base(5) + (-t347 * t299 + t303 * t351 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((t298 * t351 + t347 * t302 + Icges(1,2)) * V_base(5) + (t299 * t351 + t347 * t303 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(2) * (t170 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(9) * (t139 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(3) * (t140 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(11) * (t124 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t125 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(10) * (t128 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(8) * (t131 ^ 2 + t135 ^ 2 + t137 ^ 2) / 0.2e1 + m(4) * (t132 ^ 2 + t136 ^ 2 + t138 ^ 2) / 0.2e1 + m(6) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + t286 * ((-t148 * t171 - t149 * t172 - t166 * t286) * t323 + ((-t151 * t345 + t153 * t349) * t172 + (-t150 * t345 + t152 * t349) * t171 + (-t167 * t345 + t168 * t349) * t286) * t321) / 0.2e1 + t312 * (t364 * t347 + t352 * t351) / 0.2e1 + t310 * (t352 * t347 - t364 * t351) / 0.2e1 + t311 * (t365 * t347 + t353 * t351) / 0.2e1 + t309 * (t353 * t347 - t365 * t351) / 0.2e1 + t285 * (t366 * t347 + t354 * t351) / 0.2e1 + t282 * (t354 * t347 - t366 * t351) / 0.2e1 + t284 * (t367 * t347 + t355 * t351) / 0.2e1 + t281 * (t355 * t347 - t367 * t351) / 0.2e1 + t283 * (t368 * t347 + t356 * t351) / 0.2e1 + t280 * (t356 * t347 - t368 * t351) / 0.2e1 + t253 * (t369 * t347 + t357 * t351) / 0.2e1 + t251 * (t357 * t347 - t369 * t351) / 0.2e1 + t252 * (t370 * t347 + t358 * t351) / 0.2e1 + t250 * (t358 * t347 - t370 * t351) / 0.2e1 + t249 * (t371 * t347 + t359 * t351) / 0.2e1 + t248 * (t359 * t347 - t371 * t351) / 0.2e1 + ((t231 * t348 + t235 * t344) * t311 + (t230 * t348 + t234 * t344) * t309 + (-t233 * t346 + t237 * t350) * t312 + (-t232 * t346 + t236 * t350) * t310 + (-t204 * t332 + t210 * t335) * t284 + (-t203 * t332 + t209 * t335) * t281 + (t206 * t336 + t212 * t333) * t285 + (t205 * t336 + t211 * t333) * t282 + (-t202 * t331 + t208 * t334) * t283 + (-t201 * t331 + t207 * t334) * t280 + (t182 * t320 - t186 * t322) * t252 + (t181 * t320 - t185 * t322) * t250 + (t184 * t323 + t188 * t321) * t253 + (t183 * t323 + t187 * t321) * t251 + (-t161 * t314 - t163 * t315) * t249 + (-t160 * t314 - t162 * t315) * t248 + (Icges(2,3) + t296 * t348 + t300 * t344 - t297 * t346 + t301 * t350 - t272 * t332 + t275 * t335 + t273 * t336 + t276 * t333 - t271 * t331 + t274 * t334 - t223 * t314 - t224 * t315 + t258 * t320 - t260 * t322 + t259 * t323 + t261 * t321) * t325) * t325 / 0.2e1 + t172 * ((t149 * t428 + t246 * t151 + t247 * t153) * t172 + (t148 * t428 + t246 * t150 + t247 * t152) * t171 + (t166 * t428 + t246 * t167 + t247 * t168) * t286) / 0.2e1 + t171 * ((t149 * t429 + t151 * t244 + t153 * t245) * t172 + (t148 * t429 + t244 * t150 + t245 * t152) * t171 + (t166 * t429 + t167 * t244 + t168 * t245) * t286) / 0.2e1;
T = t1;
