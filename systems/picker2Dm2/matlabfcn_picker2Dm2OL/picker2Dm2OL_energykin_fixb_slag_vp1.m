% Calculate kinetic energy for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
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
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm2OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2OL_energykin_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:51
% EndTime: 2020-05-09 23:18:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (209->86), mult. (188->139), div. (0->0), fcn. (66->20), ass. (0->67)
t330 = qJD(1) + qJD(2);
t318 = qJD(4) + t330;
t350 = pkin(4) * t318;
t319 = qJD(3) + t330;
t349 = pkin(6) * t319;
t348 = pkin(1) * qJD(1);
t333 = qJ(1) + qJ(2);
t323 = sin(t333);
t347 = t323 * t330;
t325 = cos(t333);
t346 = t325 * t330;
t335 = sin(qJ(1));
t316 = t335 * t348;
t345 = pkin(3) * t347 + t316;
t344 = pkin(2) * t347 + t316;
t337 = cos(qJ(1));
t343 = t337 * t348;
t328 = qJ(3) + t333;
t327 = qJ(4) + t333;
t342 = -pkin(2) * t346 - t343;
t341 = -pkin(3) * t346 - t343;
t336 = cos(qJ(7));
t334 = sin(qJ(7));
t332 = qJ(1) + qJ(8);
t331 = pkin(8) + qJ(5);
t329 = qJD(1) + qJD(8);
t326 = qJ(6) + t333;
t324 = cos(t332);
t322 = sin(t332);
t321 = cos(t331);
t320 = sin(t331);
t317 = qJD(6) + t330;
t315 = qJ(9) + t328;
t314 = qJ(10) + t327;
t313 = cos(t328);
t312 = cos(t327);
t311 = cos(t326);
t310 = sin(t328);
t309 = sin(t327);
t308 = sin(t326);
t307 = qJD(9) + t319;
t306 = qJD(10) + t318;
t305 = cos(t315);
t304 = sin(t315);
t303 = cos(t314);
t302 = sin(t314);
t301 = -t337 * rSges(2,1) + t335 * rSges(2,2);
t300 = -t336 * rSges(8,1) + t334 * rSges(8,2);
t299 = -t335 * rSges(2,1) - t337 * rSges(2,2);
t298 = t334 * rSges(8,1) + t336 * rSges(8,2);
t295 = t321 * rSges(6,1) - t320 * rSges(6,2);
t294 = t320 * rSges(6,1) + t321 * rSges(6,2);
t293 = -t343 + t330 * (-t325 * rSges(3,1) + t323 * rSges(3,2));
t292 = -t343 + t329 * (t324 * rSges(9,1) - t322 * rSges(9,2));
t291 = t316 - t330 * (-t323 * rSges(3,1) - t325 * rSges(3,2));
t290 = t316 - t329 * (t322 * rSges(9,1) + t324 * rSges(9,2));
t289 = -t343 + t317 * (t311 * rSges(7,1) - t308 * rSges(7,2));
t288 = t316 - t317 * (t308 * rSges(7,1) + t311 * rSges(7,2));
t287 = t319 * (t313 * rSges(4,1) - t310 * rSges(4,2)) + t342;
t286 = t318 * (-t312 * rSges(5,1) + t309 * rSges(5,2)) + t341;
t285 = -t319 * (t310 * rSges(4,1) + t313 * rSges(4,2)) + t344;
t284 = -t318 * (-t309 * rSges(5,1) - t312 * rSges(5,2)) + t345;
t283 = t313 * t349 + t307 * (-t305 * rSges(10,1) + t304 * rSges(10,2)) + t342;
t282 = -t310 * t349 - t307 * (-t304 * rSges(10,1) - t305 * rSges(10,2)) + t344;
t281 = -t312 * t350 + t306 * (t303 * rSges(11,1) - t302 * rSges(11,2)) + t341;
t280 = t309 * t350 - t306 * (t302 * rSges(11,1) + t303 * rSges(11,2)) + t345;
t1 = m(10) * (t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(11) * (t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(9) * (t290 ^ 2 + t292 ^ 2) / 0.2e1 + m(7) * (t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(5) * (t284 ^ 2 + t286 ^ 2) / 0.2e1 + m(4) * (t285 ^ 2 + t287 ^ 2) / 0.2e1 + m(3) * (t291 ^ 2 + t293 ^ 2) / 0.2e1 + t330 ^ 2 * Icges(3,3) / 0.2e1 + t319 ^ 2 * Icges(4,3) / 0.2e1 + t318 ^ 2 * Icges(5,3) / 0.2e1 + t306 ^ 2 * Icges(11,3) / 0.2e1 + t329 ^ 2 * Icges(9,3) / 0.2e1 + t307 ^ 2 * Icges(10,3) / 0.2e1 + t317 ^ 2 * Icges(7,3) / 0.2e1 + (m(2) * (t299 ^ 2 + t301 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (m(6) * (t294 ^ 2 + t295 ^ 2) / 0.2e1 + Icges(6,3) / 0.2e1) * qJD(5) ^ 2 + (m(8) * (t298 ^ 2 + t300 ^ 2) / 0.2e1 + Icges(8,3) / 0.2e1) * qJD(7) ^ 2;
T = t1;
