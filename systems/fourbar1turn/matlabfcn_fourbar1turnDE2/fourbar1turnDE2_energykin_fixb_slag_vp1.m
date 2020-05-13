% Calculate kinetic energy for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:26
% EndTime: 2020-04-12 19:33:28
% DurationCPUTime: 2.42s
% Computational Cost: add. (15378->180), mult. (21633->361), div. (1116->12), fcn. (6354->11), ass. (0->128)
t387 = pkin(4) ^ 2;
t318 = sin(qJ(1));
t320 = cos(qJ(1));
t319 = cos(qJ(2));
t378 = pkin(2) * t319;
t313 = pkin(1) - t378;
t328 = pkin(1) ^ 2;
t363 = -0.2e1 * pkin(1) * t378 + t328;
t382 = -pkin(3) - pkin(4);
t300 = (pkin(2) - t382) * (pkin(2) + t382) + t363;
t381 = -pkin(3) + pkin(4);
t301 = (pkin(2) - t381) * (pkin(2) + t381) + t363;
t329 = sqrt(-t300 * t301);
t327 = pkin(2) ^ 2;
t312 = t327 + t363;
t362 = pkin(3) ^ 2 - t387;
t303 = t312 - t362;
t317 = sin(qJ(2));
t366 = t317 * t303;
t284 = pkin(2) * t366 + t313 * t329;
t281 = t284 ^ 2;
t307 = 0.1e1 / t312;
t308 = 0.1e1 / t312 ^ 2;
t323 = 0.1e1 / pkin(4);
t365 = t317 * t329;
t283 = -pkin(2) * t365 + t303 * t313;
t330 = t283 ^ 2;
t355 = ((t281 + t330) / t387 * t308) ^ (-0.1e1 / 0.2e1) * t307 * t323;
t372 = Icges(5,4) * t283;
t334 = (-Icges(5,2) * t284 - t372) * t355;
t251 = -Icges(5,6) * t320 + t318 * t334;
t252 = Icges(5,6) * t318 + t320 * t334;
t371 = Icges(5,4) * t284;
t335 = (-Icges(5,1) * t283 - t371) * t355;
t253 = -Icges(5,5) * t320 + t318 * t335;
t254 = Icges(5,5) * t318 + t320 * t335;
t386 = ((-t252 * t284 - t254 * t283) * t318 + (t251 * t284 + t253 * t283) * t320) * t355;
t385 = -0.2e1 * pkin(2);
t384 = t318 ^ 2;
t383 = t320 ^ 2;
t379 = pkin(1) * t317;
t376 = Icges(3,4) * t317;
t375 = Icges(3,4) * t319;
t302 = t312 + t362;
t314 = pkin(1) * t319 - pkin(2);
t282 = -pkin(1) * t365 - t302 * t314;
t285 = t302 * t379 - t314 * t329;
t326 = 0.1e1 / pkin(3);
t367 = t307 * t326;
t275 = qJ(2) + atan2(t285 * t367, t282 * t367);
t273 = sin(t275);
t374 = Icges(4,4) * t273;
t274 = cos(t275);
t373 = Icges(4,4) * t274;
t280 = 0.1e1 / t330;
t358 = t317 * qJD(2);
t354 = (-t300 - t301) * pkin(2) * pkin(1) * t358 / t329 * t307;
t351 = t354 / 0.2e1;
t364 = t319 * t329;
t368 = t307 * t317 ^ 2;
t247 = 0.2e1 * (-(t313 * t351 + (t327 * pkin(1) * t368 + ((t303 * t319 + t365) * t307 / 0.2e1 - t284 * t308 * t379) * pkin(2)) * qJD(2)) / t283 - (t317 * t351 + (-(-t364 + t366) * t307 / 0.2e1 + (t283 * t308 - t307 * t313) * t379) * qJD(2)) * pkin(2) * t284 * t280) * pkin(4) / (t280 * t281 + 0.1e1) * t312 * t323;
t370 = t247 * t318;
t369 = t247 * t320;
t361 = qJD(1) * t319;
t360 = qJD(2) * t318;
t359 = qJD(2) * t320;
t279 = 0.1e1 / t282 ^ 2;
t356 = t308 * t385;
t357 = qJD(2) + ((-t314 * t354 + (0.2e1 * t328 * pkin(2) * t368 + ((t302 * t319 + t365) * t307 + t285 * t317 * t356) * pkin(1)) * qJD(2)) / t282 - (-t317 * t354 + (-t307 * t364 + ((t314 * t385 + t302) * t307 + t282 * t356) * t317) * qJD(2)) * pkin(1) * t285 * t279) * pkin(3) / (t279 * t285 ^ 2 + 0.1e1) * t312 * t326;
t350 = rSges(3,1) * t319 - rSges(3,2) * t317;
t349 = -rSges(4,1) * t274 + rSges(4,2) * t273;
t348 = Icges(3,1) * t319 - t376;
t347 = -Icges(4,1) * t274 + t374;
t346 = -Icges(3,2) * t317 + t375;
t345 = Icges(4,2) * t273 - t373;
t344 = Icges(3,5) * t319 - Icges(3,6) * t317;
t343 = -Icges(4,5) * t274 + Icges(4,6) * t273;
t294 = -Icges(3,6) * t320 + t346 * t318;
t296 = -Icges(3,5) * t320 + t348 * t318;
t341 = t294 * t317 - t296 * t319;
t295 = Icges(3,6) * t318 + t346 * t320;
t297 = Icges(3,5) * t318 + t348 * t320;
t340 = -t295 * t317 + t297 * t319;
t305 = Icges(3,2) * t319 + t376;
t306 = Icges(3,1) * t317 + t375;
t339 = -t305 * t317 + t306 * t319;
t245 = t357 * t318;
t246 = t357 * t320;
t337 = (-Icges(4,5) * t273 - Icges(4,6) * t274) * qJD(1) + (Icges(4,3) * t318 + t343 * t320) * t245 - (-Icges(4,3) * t320 + t343 * t318) * t246;
t336 = (-rSges(5,1) * t283 - rSges(5,2) * t284) * t355;
t333 = (-Icges(5,5) * t283 - Icges(5,6) * t284) * t355;
t258 = (-Icges(5,2) * t283 + t371) * t355;
t259 = (Icges(5,1) * t284 - t372) * t355;
t332 = (-t258 * t284 - t259 * t283) * t355;
t263 = -Icges(4,6) * t320 + t345 * t318;
t264 = Icges(4,6) * t318 + t345 * t320;
t265 = -Icges(4,5) * t320 + t347 * t318;
t266 = Icges(4,5) * t318 + t347 * t320;
t270 = -Icges(4,2) * t274 - t374;
t271 = -Icges(4,1) * t273 - t373;
t331 = (t264 * t273 - t266 * t274) * t245 - (t263 * t273 - t265 * t274) * t246 + (t270 * t273 - t271 * t274) * qJD(1);
t311 = rSges(2,1) * t320 - rSges(2,2) * t318;
t310 = rSges(2,1) * t318 + rSges(2,2) * t320;
t309 = rSges(3,1) * t317 + rSges(3,2) * t319;
t304 = Icges(3,5) * t317 + Icges(3,6) * t319;
t299 = rSges(3,3) * t318 + t350 * t320;
t298 = -rSges(3,3) * t320 + t350 * t318;
t293 = Icges(3,3) * t318 + t344 * t320;
t292 = -Icges(3,3) * t320 + t344 * t318;
t289 = -qJD(1) * t298 - t309 * t359;
t288 = qJD(1) * t299 - t309 * t360;
t286 = (t298 * t318 + t299 * t320) * qJD(2);
t272 = -rSges(4,1) * t273 - rSges(4,2) * t274;
t268 = rSges(4,3) * t318 + t349 * t320;
t267 = -rSges(4,3) * t320 + t349 * t318;
t260 = (rSges(5,1) * t284 - rSges(5,2) * t283) * t355;
t257 = (Icges(5,5) * t284 - Icges(5,6) * t283) * t355;
t256 = rSges(5,3) * t318 + t320 * t336;
t255 = -rSges(5,3) * t320 + t318 * t336;
t250 = Icges(5,3) * t318 + t320 * t333;
t249 = -Icges(5,3) * t320 + t318 * t333;
t244 = -qJD(1) * t267 - t246 * t272 + (-t318 * t361 - t320 * t358) * pkin(2);
t243 = qJD(1) * t268 - t245 * t272 + (-t318 * t358 + t320 * t361) * pkin(2);
t242 = -t260 * t370 + (pkin(1) * t320 + t256) * qJD(1);
t241 = -t260 * t369 + (-pkin(1) * t318 - t255) * qJD(1);
t240 = t245 * t267 + t246 * t268 + (t383 + t384) * qJD(2) * t378;
t239 = (t255 * t318 + t256 * t320) * t247;
t1 = m(3) * (t286 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + ((t318 * t304 + t339 * t320) * qJD(1) + (t384 * t293 + (t341 * t320 + (-t292 + t340) * t318) * t320) * qJD(2)) * t360 / 0.2e1 - ((-t320 * t304 + t339 * t318) * qJD(1) + (t383 * t292 + (t340 * t318 + (-t293 + t341) * t320) * t318) * qJD(2)) * t359 / 0.2e1 + m(4) * (t240 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + t245 * (t337 * t318 + t331 * t320) / 0.2e1 - t246 * (t331 * t318 - t337 * t320) / 0.2e1 + m(5) * (t239 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + ((t318 * t257 + t320 * t332) * qJD(1) + (t384 * t250 + (-t318 * t249 + t386) * t320) * t247) * t370 / 0.2e1 - ((-t320 * t257 + t318 * t332) * qJD(1) + (t383 * t249 + (-t320 * t250 + t386) * t318) * t247) * t369 / 0.2e1 + (m(2) * (t310 ^ 2 + t311 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t295 * t319 + t297 * t317) * t318 - (t294 * t319 + t296 * t317) * t320) * qJD(2) + (-t264 * t274 - t266 * t273) * t245 - (-t263 * t274 - t265 * t273) * t246 + ((-t252 * t283 + t254 * t284) * t318 - (-t251 * t283 + t253 * t284) * t320) * t247 * t355 + (t319 * t305 + t317 * t306 - t274 * t270 - t273 * t271 + (-t258 * t283 + t259 * t284) * t355) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
