% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:28
% EndTime: 2020-04-12 19:25:31
% DurationCPUTime: 2.83s
% Computational Cost: add. (22086->193), mult. (31825->381), div. (1740->13), fcn. (9214->10), ass. (0->130)
t407 = pkin(4) ^ 2;
t406 = pkin(3) ^ 2;
t349 = pkin(2) ^ 2;
t350 = pkin(1) ^ 2;
t340 = cos(qJ(2));
t396 = pkin(2) * t340;
t383 = -0.2e1 * pkin(1) * t396 + t350;
t333 = t349 + t383;
t382 = t406 - t407;
t323 = t333 + t382;
t335 = pkin(1) * t340 - pkin(2);
t338 = sin(qJ(2));
t400 = -pkin(3) - pkin(4);
t321 = (pkin(2) - t400) * (pkin(2) + t400) + t383;
t399 = -pkin(3) + pkin(4);
t322 = (pkin(2) - t399) * (pkin(2) + t399) + t383;
t351 = sqrt(-t321 * t322);
t385 = t338 * t351;
t303 = -pkin(1) * t385 - t323 * t335;
t397 = pkin(1) * t338;
t306 = t323 * t397 - t335 * t351;
t302 = t306 ^ 2;
t328 = 0.1e1 / t333;
t329 = 0.1e1 / t333 ^ 2;
t347 = 0.1e1 / pkin(3);
t352 = t303 ^ 2;
t374 = ((t302 + t352) / t406 * t329) ^ (-0.1e1 / 0.2e1) * t328 * t347;
t405 = (t303 * t338 + t306 * t340) * t374;
t339 = sin(qJ(1));
t341 = cos(qJ(1));
t334 = pkin(1) - t396;
t324 = t333 - t382;
t387 = t324 * t338;
t305 = pkin(2) * t387 + t334 * t351;
t301 = t305 ^ 2;
t344 = 0.1e1 / pkin(4);
t304 = -pkin(2) * t385 + t324 * t334;
t353 = t304 ^ 2;
t375 = ((t301 + t353) / t407 * t329) ^ (-0.1e1 / 0.2e1) * t328 * t344;
t392 = Icges(5,4) * t304;
t357 = (-Icges(5,2) * t305 - t392) * t375;
t279 = -Icges(5,6) * t341 + t339 * t357;
t280 = Icges(5,6) * t339 + t341 * t357;
t391 = Icges(5,4) * t305;
t358 = (-Icges(5,1) * t304 - t391) * t375;
t281 = -Icges(5,5) * t341 + t339 * t358;
t282 = Icges(5,5) * t339 + t341 * t358;
t404 = ((-t280 * t305 - t282 * t304) * t339 + (t279 * t305 + t281 * t304) * t341) * t375;
t403 = -0.2e1 * pkin(2);
t402 = t339 ^ 2;
t401 = t341 ^ 2;
t394 = Icges(3,4) * t338;
t393 = Icges(3,4) * t340;
t300 = 0.1e1 / t353;
t378 = t338 * qJD(2);
t373 = (-t321 - t322) * pkin(2) * pkin(1) * t378 / t351 * t328;
t370 = t373 / 0.2e1;
t384 = t340 * t351;
t386 = t328 * t338 ^ 2;
t263 = 0.2e1 * (-(t334 * t370 + (t349 * pkin(1) * t386 + ((t324 * t340 + t385) * t328 / 0.2e1 - t305 * t329 * t397) * pkin(2)) * qJD(2)) / t304 - (t338 * t370 + (-(-t384 + t387) * t328 / 0.2e1 + (t304 * t329 - t328 * t334) * t397) * qJD(2)) * pkin(2) * t305 * t300) * pkin(4) / (t300 * t301 + 0.1e1) * t333 * t344;
t390 = t263 * t339;
t389 = t263 * t341;
t388 = t306 * t338;
t381 = qJD(1) * t340;
t380 = qJD(2) * t339;
t379 = qJD(2) * t341;
t299 = 0.1e1 / t352;
t376 = t329 * t403;
t377 = qJD(2) + ((-t335 * t373 + (0.2e1 * t350 * pkin(2) * t386 + ((t323 * t340 + t385) * t328 + t376 * t388) * pkin(1)) * qJD(2)) / t303 - (-t338 * t373 + (-t328 * t384 + ((t335 * t403 + t323) * t328 + t303 * t376) * t338) * qJD(2)) * pkin(1) * t306 * t299) * pkin(3) / (t299 * t302 + 0.1e1) * t333 * t347;
t369 = rSges(3,1) * t340 - rSges(3,2) * t338;
t368 = Icges(3,1) * t340 - t394;
t367 = -Icges(3,2) * t338 + t393;
t366 = Icges(3,5) * t340 - Icges(3,6) * t338;
t315 = -Icges(3,6) * t341 + t367 * t339;
t317 = -Icges(3,5) * t341 + t368 * t339;
t363 = t315 * t338 - t317 * t340;
t316 = Icges(3,6) * t339 + t367 * t341;
t318 = Icges(3,5) * t339 + t368 * t341;
t362 = -t316 * t338 + t318 * t340;
t326 = Icges(3,2) * t340 + t394;
t327 = Icges(3,1) * t338 + t393;
t361 = -t326 * t338 + t327 * t340;
t359 = (-rSges(5,1) * t304 - rSges(5,2) * t305) * t375;
t356 = (-Icges(5,5) * t304 - Icges(5,6) * t305) * t375;
t292 = (-Icges(5,2) * t304 + t391) * t375;
t293 = (Icges(5,1) * t305 - t392) * t375;
t355 = (-t292 * t305 - t293 * t304) * t375;
t290 = (-t303 * t340 + t388) * t374;
t332 = rSges(2,1) * t341 - rSges(2,2) * t339;
t331 = rSges(2,1) * t339 + rSges(2,2) * t341;
t330 = rSges(3,1) * t338 + rSges(3,2) * t340;
t325 = Icges(3,5) * t338 + Icges(3,6) * t340;
t320 = rSges(3,3) * t339 + t369 * t341;
t319 = -rSges(3,3) * t341 + t369 * t339;
t314 = Icges(3,3) * t339 + t366 * t341;
t313 = -Icges(3,3) * t341 + t366 * t339;
t310 = -qJD(1) * t319 - t330 * t379;
t309 = qJD(1) * t320 - t330 * t380;
t307 = (t319 * t339 + t320 * t341) * qJD(2);
t294 = (rSges(5,1) * t305 - rSges(5,2) * t304) * t375;
t291 = (Icges(5,5) * t305 - Icges(5,6) * t304) * t375;
t288 = t341 * t405;
t287 = t341 * t290;
t286 = t339 * t405;
t285 = t339 * t290;
t284 = rSges(5,3) * t339 + t341 * t359;
t283 = -rSges(5,3) * t341 + t339 * t359;
t278 = Icges(5,3) * t339 + t341 * t356;
t277 = -Icges(5,3) * t341 + t339 * t356;
t276 = -rSges(4,1) * t405 + rSges(4,2) * t290;
t275 = -Icges(4,1) * t405 + Icges(4,4) * t290;
t274 = -Icges(4,4) * t405 + Icges(4,2) * t290;
t273 = -Icges(4,5) * t405 + Icges(4,6) * t290;
t272 = rSges(4,1) * t287 + rSges(4,2) * t288 + rSges(4,3) * t339;
t271 = rSges(4,1) * t285 + rSges(4,2) * t286 - rSges(4,3) * t341;
t270 = Icges(4,1) * t287 + Icges(4,4) * t288 + Icges(4,5) * t339;
t269 = Icges(4,1) * t285 + Icges(4,4) * t286 - Icges(4,5) * t341;
t268 = Icges(4,4) * t287 + Icges(4,2) * t288 + Icges(4,6) * t339;
t267 = Icges(4,4) * t285 + Icges(4,2) * t286 - Icges(4,6) * t341;
t266 = Icges(4,5) * t287 + Icges(4,6) * t288 + Icges(4,3) * t339;
t265 = Icges(4,5) * t285 + Icges(4,6) * t286 - Icges(4,3) * t341;
t262 = t377 * t341;
t261 = t377 * t339;
t260 = -t294 * t390 + (pkin(1) * t341 + t284) * qJD(1);
t259 = -t294 * t389 + (-pkin(1) * t339 - t283) * qJD(1);
t258 = -qJD(1) * t271 - t262 * t276 + (-t339 * t381 - t341 * t378) * pkin(2);
t257 = qJD(1) * t272 - t261 * t276 + (-t339 * t378 + t341 * t381) * pkin(2);
t256 = (t283 * t339 + t284 * t341) * t263;
t255 = t261 * t271 + t262 * t272 + (t401 + t402) * qJD(2) * t396;
t1 = m(3) * (t307 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + ((t339 * t325 + t361 * t341) * qJD(1) + (t402 * t314 + (t363 * t341 + (-t313 + t362) * t339) * t341) * qJD(2)) * t380 / 0.2e1 - ((-t341 * t325 + t361 * t339) * qJD(1) + (t401 * t313 + (t362 * t339 + (-t314 + t363) * t341) * t339) * qJD(2)) * t379 / 0.2e1 + m(4) * (t255 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t261 * ((t339 * t266 + t288 * t268 + t287 * t270) * t261 - (t265 * t339 + t267 * t288 + t269 * t287) * t262 + (t273 * t339 + t274 * t288 + t275 * t287) * qJD(1)) / 0.2e1 - t262 * ((-t266 * t341 + t268 * t286 + t270 * t285) * t261 - (-t341 * t265 + t286 * t267 + t285 * t269) * t262 + (-t273 * t341 + t274 * t286 + t275 * t285) * qJD(1)) / 0.2e1 + m(5) * (t256 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + ((t339 * t291 + t341 * t355) * qJD(1) + (t402 * t278 + (-t339 * t277 + t404) * t341) * t263) * t390 / 0.2e1 - ((-t341 * t291 + t339 * t355) * qJD(1) + (t401 * t277 + (-t341 * t278 + t404) * t339) * t263) * t389 / 0.2e1 + (m(2) * (t331 ^ 2 + t332 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t316 * t340 + t318 * t338) * t339 - (t315 * t340 + t338 * t317) * t341) * qJD(2) + (t268 * t290 - t270 * t405) * t261 - (t267 * t290 - t269 * t405) * t262 + ((-t280 * t304 + t282 * t305) * t339 - (-t279 * t304 + t281 * t305) * t341) * t263 * t375 + (t340 * t326 + t338 * t327 + t290 * t274 - t405 * t275 + (-t292 * t304 + t293 * t305) * t375) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
