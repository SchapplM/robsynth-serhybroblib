% Calculate kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:00:31
% EndTime: 2020-05-03 00:00:33
% DurationCPUTime: 1.51s
% Computational Cost: add. (1084->225), mult. (1193->384), div. (0->0), fcn. (1101->12), ass. (0->134)
t407 = pkin(2) * qJD(2);
t356 = sin(qJ(2));
t406 = Icges(3,4) * t356;
t360 = cos(qJ(2));
t405 = Icges(3,4) * t360;
t353 = qJ(2) + qJ(3);
t350 = sin(t353);
t404 = Icges(4,4) * t350;
t351 = cos(t353);
t403 = Icges(4,4) * t351;
t352 = qJ(4) + t353;
t343 = sin(t352);
t402 = Icges(5,4) * t343;
t344 = cos(t352);
t401 = Icges(5,4) * t344;
t359 = cos(qJ(3));
t345 = pkin(3) * t359 + pkin(2);
t355 = sin(qJ(3));
t316 = -pkin(3) * t355 * t356 + t345 * t360 + pkin(1);
t357 = sin(qJ(1));
t400 = t316 * t357;
t399 = t343 * t357;
t361 = cos(qJ(1));
t398 = t343 * t361;
t354 = sin(qJ(5));
t397 = t354 * t357;
t396 = t354 * t361;
t395 = t355 * t360;
t358 = cos(qJ(5));
t394 = t357 * t358;
t393 = t358 * t361;
t349 = qJD(2) * t361;
t331 = qJD(3) * t361 + t349;
t392 = qJD(2) * t357;
t391 = qJD(5) * t343;
t390 = -qJD(2) - qJD(3);
t389 = t356 * t407;
t388 = t360 * t407;
t324 = qJD(4) * t361 + t331;
t387 = pkin(3) * t390 * t351 - t388;
t386 = pkin(4) * t344 + pkin(6) * t343;
t341 = -rSges(3,1) * t360 + rSges(3,2) * t356;
t385 = rSges(4,1) * t351 - rSges(4,2) * t350;
t384 = rSges(5,1) * t344 - rSges(5,2) * t343;
t340 = rSges(6,1) * t358 - rSges(6,2) * t354;
t323 = (-qJD(4) + t390) * t357;
t383 = rSges(6,3) * t343 + t340 * t344;
t382 = Icges(3,1) * t360 - t406;
t381 = Icges(4,1) * t351 - t404;
t380 = Icges(5,1) * t344 - t402;
t379 = -Icges(3,2) * t356 + t405;
t378 = -Icges(4,2) * t350 + t403;
t377 = -Icges(5,2) * t343 + t401;
t376 = Icges(3,5) * t360 - Icges(3,6) * t356;
t375 = Icges(4,5) * t351 - Icges(4,6) * t350;
t374 = Icges(5,5) * t344 - Icges(5,6) * t343;
t308 = Icges(3,6) * t361 + t379 * t357;
t310 = Icges(3,5) * t361 + t382 * t357;
t373 = -t308 * t356 + t310 * t360;
t309 = -Icges(3,6) * t357 + t379 * t361;
t311 = -Icges(3,5) * t357 + t382 * t361;
t372 = t309 * t356 - t311 * t360;
t334 = -Icges(3,2) * t360 - t406;
t335 = -Icges(3,1) * t356 - t405;
t371 = -t334 * t356 + t335 * t360;
t370 = pkin(1) - t341;
t369 = -pkin(3) * qJD(3) * (t356 * t359 + t395) - qJD(2) * (pkin(3) * t395 + t345 * t356);
t368 = t369 * t361;
t367 = (-Icges(5,5) * t343 - Icges(5,6) * t344) * qJD(1) + (Icges(5,3) * t361 + t374 * t357) * t324 + (-Icges(5,3) * t357 + t374 * t361) * t323;
t330 = t390 * t357;
t366 = (-Icges(4,5) * t350 - Icges(4,6) * t351) * qJD(1) + (Icges(4,3) * t361 + t375 * t357) * t331 + (-Icges(4,3) * t357 + t375 * t361) * t330;
t365 = t316 * t361 * qJD(1) + t369 * t357;
t287 = Icges(5,6) * t361 + t377 * t357;
t288 = -Icges(5,6) * t357 + t377 * t361;
t289 = Icges(5,5) * t361 + t380 * t357;
t290 = -Icges(5,5) * t357 + t380 * t361;
t319 = -Icges(5,2) * t344 - t402;
t320 = -Icges(5,1) * t343 - t401;
t364 = (-t288 * t343 + t290 * t344) * t323 + (-t287 * t343 + t289 * t344) * t324 + (-t319 * t343 + t320 * t344) * qJD(1);
t297 = Icges(4,6) * t361 + t378 * t357;
t298 = -Icges(4,6) * t357 + t378 * t361;
t299 = Icges(4,5) * t361 + t381 * t357;
t300 = -Icges(4,5) * t357 + t381 * t361;
t326 = -Icges(4,2) * t351 - t404;
t327 = -Icges(4,1) * t350 - t403;
t363 = (-t298 * t350 + t300 * t351) * t330 + (-t297 * t350 + t299 * t351) * t331 + (-t326 * t350 + t327 * t351) * qJD(1);
t346 = pkin(2) * t360 + pkin(1);
t342 = rSges(2,1) * t361 - rSges(2,2) * t357;
t339 = rSges(6,1) * t354 + rSges(6,2) * t358;
t338 = rSges(2,1) * t357 + rSges(2,2) * t361;
t337 = -rSges(3,1) * t356 - rSges(3,2) * t360;
t333 = -Icges(3,5) * t356 - Icges(3,6) * t360;
t332 = qJD(5) * t344 + qJD(1);
t328 = -rSges(4,1) * t350 - rSges(4,2) * t351;
t322 = -pkin(4) * t343 + pkin(6) * t344;
t321 = -rSges(5,1) * t343 - rSges(5,2) * t344;
t315 = t344 * t393 - t397;
t314 = -t344 * t396 - t394;
t313 = t344 * t394 + t396;
t312 = -t344 * t397 + t393;
t307 = -Icges(3,3) * t357 + t376 * t361;
t306 = Icges(3,3) * t361 + t376 * t357;
t305 = t386 * t361;
t304 = t386 * t357;
t302 = -rSges(4,3) * t357 + t385 * t361;
t301 = rSges(4,3) * t361 + t385 * t357;
t294 = -rSges(5,3) * t357 + t384 * t361;
t293 = rSges(5,3) * t361 + t384 * t357;
t292 = t357 * t391 + t324;
t291 = t361 * t391 + t323;
t284 = rSges(6,3) * t344 - t340 * t343;
t283 = Icges(6,5) * t344 + (-Icges(6,1) * t358 + Icges(6,4) * t354) * t343;
t282 = Icges(6,6) * t344 + (-Icges(6,4) * t358 + Icges(6,2) * t354) * t343;
t281 = Icges(6,3) * t344 + (-Icges(6,5) * t358 + Icges(6,6) * t354) * t343;
t280 = t337 * t392 + (-t357 * rSges(3,3) + t370 * t361) * qJD(1);
t279 = t337 * t349 + (-t361 * rSges(3,3) - t370 * t357) * qJD(1);
t278 = -t339 * t357 + t383 * t361;
t277 = t339 * t361 + t383 * t357;
t276 = Icges(6,1) * t315 + Icges(6,4) * t314 + Icges(6,5) * t398;
t275 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t399;
t274 = Icges(6,4) * t315 + Icges(6,2) * t314 + Icges(6,6) * t398;
t273 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t399;
t272 = Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t398;
t271 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t399;
t270 = -t357 * t389 - t328 * t330 + (t346 * t361 + t302) * qJD(1);
t269 = -t361 * t389 + t328 * t331 + (-t346 * t357 - t301) * qJD(1);
t268 = t301 * t330 - t302 * t331 - t388;
t267 = t293 * t323 - t294 * t324 + t387;
t266 = qJD(1) * t294 - t321 * t323 + t365;
t265 = t321 * t324 + t368 + (-t293 - t400) * qJD(1);
t264 = qJD(1) * t305 + t278 * t332 - t284 * t291 - t322 * t323 + t365;
t263 = -t277 * t332 + t284 * t292 + t322 * t324 + t368 + (-t304 - t400) * qJD(1);
t262 = t277 * t291 - t278 * t292 + t304 * t323 - t305 * t324 + t387;
t1 = m(3) * (qJD(2) ^ 2 * t341 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 - ((-t357 * t333 + t371 * t361) * qJD(1) + (t357 ^ 2 * t307 + (t373 * t361 + (-t306 + t372) * t357) * t361) * qJD(2)) * t392 / 0.2e1 + ((t361 * t333 + t371 * t357) * qJD(1) + (t361 ^ 2 * t306 + (t372 * t357 + (-t307 + t373) * t361) * t357) * qJD(2)) * t349 / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + t330 * (-t366 * t357 + t363 * t361) / 0.2e1 + t331 * (t363 * t357 + t366 * t361) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t323 * (-t367 * t357 + t364 * t361) / 0.2e1 + t324 * (t364 * t357 + t367 * t361) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t291 * ((t272 * t398 + t314 * t274 + t315 * t276) * t291 + (t271 * t398 + t273 * t314 + t275 * t315) * t292 + (t281 * t398 + t282 * t314 + t283 * t315) * t332) / 0.2e1 + t292 * ((t272 * t399 + t274 * t312 + t276 * t313) * t291 + (t271 * t399 + t312 * t273 + t313 * t275) * t292 + (t281 * t399 + t282 * t312 + t283 * t313) * t332) / 0.2e1 + t332 * ((t271 * t292 + t272 * t291 + t281 * t332) * t344 + ((t274 * t354 - t276 * t358) * t291 + (t273 * t354 - t275 * t358) * t292 + (t282 * t354 - t283 * t358) * t332) * t343) / 0.2e1 + (m(2) * (t338 ^ 2 + t342 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((-(-t309 * t360 - t311 * t356) * t357 + (-t308 * t360 - t310 * t356) * t361) * qJD(2) + (-t298 * t351 - t300 * t350) * t330 + (-t297 * t351 - t299 * t350) * t331 + (-t288 * t344 - t290 * t343) * t323 + (-t287 * t344 - t343 * t289) * t324 + (-t344 * t319 - t343 * t320 - t351 * t326 - t350 * t327 - t360 * t334 - t356 * t335) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
