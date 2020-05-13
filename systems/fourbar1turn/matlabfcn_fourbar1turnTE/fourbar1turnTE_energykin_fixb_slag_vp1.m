% Calculate kinetic energy for
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnTE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnTE_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:31
% EndTime: 2020-04-12 19:18:34
% DurationCPUTime: 2.55s
% Computational Cost: add. (12024->191), mult. (17083->384), div. (804->11), fcn. (5002->6), ass. (0->127)
t329 = pkin(2) ^ 2;
t330 = pkin(1) ^ 2;
t322 = cos(qJ(2));
t375 = pkin(2) * t322;
t360 = -0.2e1 * pkin(1) * t375 + t330;
t315 = t329 + t360;
t359 = pkin(3) ^ 2 - pkin(4) ^ 2;
t305 = t315 + t359;
t317 = pkin(1) * t322 - pkin(2);
t320 = sin(qJ(2));
t382 = -pkin(3) - pkin(4);
t303 = (pkin(2) - t382) * (pkin(2) + t382) + t360;
t381 = -pkin(3) + pkin(4);
t304 = (pkin(2) - t381) * (pkin(2) + t381) + t360;
t331 = sqrt(-t303 * t304);
t362 = t320 * t331;
t285 = -pkin(1) * t362 - t305 * t317;
t376 = pkin(1) * t320;
t288 = t305 * t376 - t317 * t331;
t310 = 0.1e1 / t315;
t328 = 0.1e1 / pkin(3);
t364 = t310 * t328;
t387 = (t322 * t288 / 0.2e1 + t320 * t285 / 0.2e1) * t364;
t321 = sin(qJ(1));
t323 = cos(qJ(1));
t306 = t315 - t359;
t316 = pkin(1) - t375;
t286 = -pkin(2) * t362 + t306 * t316;
t380 = -t286 / 0.2e1;
t351 = Icges(5,4) * t380;
t326 = 0.1e1 / pkin(4);
t365 = t310 * t326;
t367 = t306 * t320;
t287 = pkin(2) * t367 + t316 * t331;
t379 = -t287 / 0.2e1;
t335 = (Icges(5,2) * t379 + t351) * t365;
t265 = -Icges(5,6) * t323 + t321 * t335;
t266 = Icges(5,6) * t321 + t323 * t335;
t370 = Icges(5,4) * t287;
t336 = (Icges(5,1) * t380 - t370 / 0.2e1) * t365;
t267 = -Icges(5,5) * t323 + t321 * t336;
t268 = Icges(5,5) * t321 + t323 * t336;
t378 = t287 / 0.2e1;
t386 = (t323 * (t286 * t267 / 0.2e1 + t265 * t378) + (t266 * t379 + t268 * t380) * t321) * t365;
t385 = -0.2e1 * pkin(2);
t384 = t321 ^ 2;
t383 = t323 ^ 2;
t373 = pkin(2) * qJD(2);
t372 = Icges(3,4) * t320;
t371 = Icges(3,4) * t322;
t284 = 0.1e1 / t286 ^ 2;
t311 = 0.1e1 / t315 ^ 2;
t352 = (-t303 - t304) * t373 * t376 / t331 * t310;
t348 = t352 / 0.2e1;
t361 = t322 * t331;
t366 = t310 * t320 ^ 2;
t249 = 0.2e1 * (-(t316 * t348 + (t329 * pkin(1) * t366 + ((t306 * t322 + t362) * t310 / 0.2e1 - t287 * t311 * t376) * pkin(2)) * qJD(2)) / t286 - (t320 * t348 + (-(-t361 + t367) * t310 / 0.2e1 + (t286 * t311 - t310 * t316) * t376) * qJD(2)) * pkin(2) * t287 * t284) * pkin(4) / (t284 * t287 ^ 2 + 0.1e1) * t315 * t326;
t369 = t249 * t321;
t368 = t249 * t323;
t363 = t320 * t288;
t358 = qJD(1) * t322;
t357 = qJD(2) * t321;
t356 = qJD(2) * t323;
t283 = 0.1e1 / t285 ^ 2;
t354 = t311 * t385;
t355 = qJD(2) + ((-t317 * t352 + (0.2e1 * t330 * pkin(2) * t366 + ((t305 * t322 + t362) * t310 + t354 * t363) * pkin(1)) * qJD(2)) / t285 - (-t320 * t352 + (-t310 * t361 + ((t317 * t385 + t305) * t310 + t285 * t354) * t320) * qJD(2)) * pkin(1) * t288 * t283) * pkin(3) / (t283 * t288 ^ 2 + 0.1e1) * t315 * t328;
t347 = rSges(3,1) * t322 - rSges(3,2) * t320;
t346 = Icges(3,1) * t322 - t372;
t345 = -Icges(3,2) * t320 + t371;
t344 = Icges(3,5) * t322 - Icges(3,6) * t320;
t297 = -Icges(3,6) * t323 + t321 * t345;
t299 = -Icges(3,5) * t323 + t321 * t346;
t343 = t297 * t320 - t299 * t322;
t298 = Icges(3,6) * t321 + t323 * t345;
t300 = Icges(3,5) * t321 + t323 * t346;
t342 = -t298 * t320 + t300 * t322;
t308 = Icges(3,2) * t322 + t372;
t309 = Icges(3,1) * t320 + t371;
t341 = -t308 * t320 + t309 * t322;
t337 = (rSges(5,1) * t380 + rSges(5,2) * t379) * t365;
t334 = (Icges(5,5) * t380 + Icges(5,6) * t379) * t365;
t278 = (t370 / 0.2e1 + Icges(5,2) * t380) * t365;
t279 = (Icges(5,1) * t378 + t351) * t365;
t333 = (t278 * t379 + t279 * t380) * t365;
t276 = (t363 / 0.2e1 - t322 * t285 / 0.2e1) * t364;
t314 = rSges(2,1) * t323 - rSges(2,2) * t321;
t313 = rSges(2,1) * t321 + rSges(2,2) * t323;
t312 = rSges(3,1) * t320 + rSges(3,2) * t322;
t307 = Icges(3,5) * t320 + Icges(3,6) * t322;
t302 = rSges(3,3) * t321 + t323 * t347;
t301 = -rSges(3,3) * t323 + t321 * t347;
t296 = Icges(3,3) * t321 + t323 * t344;
t295 = -Icges(3,3) * t323 + t321 * t344;
t292 = -qJD(1) * t301 - t312 * t356;
t291 = qJD(1) * t302 - t312 * t357;
t289 = (t301 * t321 + t302 * t323) * qJD(2);
t280 = (rSges(5,1) * t378 + rSges(5,2) * t380) * t365;
t277 = (Icges(5,5) * t378 + Icges(5,6) * t380) * t365;
t274 = t323 * t387;
t273 = t323 * t276;
t272 = t321 * t387;
t271 = t321 * t276;
t270 = t321 * rSges(5,3) + t323 * t337;
t269 = -t323 * rSges(5,3) + t321 * t337;
t264 = Icges(5,3) * t321 + t323 * t334;
t263 = -Icges(5,3) * t323 + t321 * t334;
t262 = -rSges(4,1) * t387 + rSges(4,2) * t276;
t261 = -Icges(4,1) * t387 + Icges(4,4) * t276;
t260 = -Icges(4,4) * t387 + Icges(4,2) * t276;
t259 = -Icges(4,5) * t387 + Icges(4,6) * t276;
t258 = rSges(4,1) * t273 + rSges(4,2) * t274 + rSges(4,3) * t321;
t257 = rSges(4,1) * t271 + rSges(4,2) * t272 - rSges(4,3) * t323;
t256 = Icges(4,1) * t273 + Icges(4,4) * t274 + Icges(4,5) * t321;
t255 = Icges(4,1) * t271 + Icges(4,4) * t272 - Icges(4,5) * t323;
t254 = Icges(4,4) * t273 + Icges(4,2) * t274 + Icges(4,6) * t321;
t253 = Icges(4,4) * t271 + Icges(4,2) * t272 - Icges(4,6) * t323;
t252 = Icges(4,5) * t273 + Icges(4,6) * t274 + Icges(4,3) * t321;
t251 = Icges(4,5) * t271 + Icges(4,6) * t272 - Icges(4,3) * t323;
t248 = t355 * t323;
t247 = t355 * t321;
t246 = -t280 * t369 + (pkin(1) * t323 + t270) * qJD(1);
t245 = -t280 * t368 + (-pkin(1) * t321 - t269) * qJD(1);
t244 = -qJD(1) * t257 - t248 * t262 + (-t320 * t356 - t321 * t358) * pkin(2);
t243 = qJD(1) * t258 - t247 * t262 + (-t320 * t357 + t323 * t358) * pkin(2);
t242 = (t269 * t321 + t270 * t323) * t249;
t241 = t247 * t257 + t248 * t258 + (t383 + t384) * t322 * t373;
t1 = m(3) * (t289 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + ((t307 * t321 + t323 * t341) * qJD(1) + (t296 * t384 + (t343 * t323 + (-t295 + t342) * t321) * t323) * qJD(2)) * t357 / 0.2e1 - ((-t307 * t323 + t321 * t341) * qJD(1) + (t295 * t383 + (t342 * t321 + (-t296 + t343) * t323) * t321) * qJD(2)) * t356 / 0.2e1 + m(4) * (t241 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + t247 * ((t252 * t321 + t254 * t274 + t256 * t273) * t247 - (t251 * t321 + t253 * t274 + t255 * t273) * t248 + (t259 * t321 + t260 * t274 + t261 * t273) * qJD(1)) / 0.2e1 - t248 * ((-t252 * t323 + t254 * t272 + t256 * t271) * t247 - (-t251 * t323 + t253 * t272 + t255 * t271) * t248 + (-t259 * t323 + t260 * t272 + t261 * t271) * qJD(1)) / 0.2e1 + m(5) * (t242 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + ((t321 * t277 + t323 * t333) * qJD(1) + (t384 * t264 + (-t321 * t263 + t386) * t323) * t249) * t369 / 0.2e1 - ((-t323 * t277 + t321 * t333) * qJD(1) + (t383 * t263 + (-t323 * t264 + t386) * t321) * t249) * t368 / 0.2e1 + (m(2) * (t313 ^ 2 + t314 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t298 * t322 + t300 * t320) * t321 - (t297 * t322 + t299 * t320) * t323) * qJD(2) + (t254 * t276 - t256 * t387) * t247 - (t253 * t276 - t255 * t387) * t248 + ((t266 * t380 + t268 * t378) * t321 - (t265 * t380 + t267 * t378) * t323) * t249 * t365 + (t308 * t322 + t309 * t320 + t260 * t276 - t261 * t387 + (t278 * t380 + t279 * t378) * t365) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
