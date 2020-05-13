% Calculate time derivative of joint inertia matrix for
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2TE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2TE_inertiaDJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:43:15
% EndTime: 2020-05-07 01:43:19
% DurationCPUTime: 3.94s
% Computational Cost: add. (665->217), mult. (1013->294), div. (0->0), fcn. (392->80), ass. (0->139)
t445 = m(5) + m(6);
t346 = pkin(17) + pkin(18);
t399 = pkin(15) + t346;
t378 = (pkin(16) + t399);
t453 = m(4) + m(8) + m(9) + t445;
t293 = pkin(4) * t445 + m(4) * rSges(4,1) + m(9) * rSges(9,1);
t441 = m(9) * rSges(9,2);
t444 = m(4) * rSges(4,2);
t322 = t441 + t444;
t354 = sin(qJ(3));
t358 = cos(qJ(3));
t279 = (t293 * t354 + t322 * t358) * pkin(1) * qJD(3);
t344 = -qJD(3) - qJD(2);
t452 = -0.2e1 * t344;
t451 = -0.2e1 * qJD(2);
t450 = 0.2e1 * qJD(2);
t449 = pkin(1) * m(6);
t448 = pkin(3) * m(9);
t447 = pkin(4) * m(6);
t446 = pkin(8) * m(6);
t443 = m(6) * rSges(6,1);
t442 = m(7) * rSges(7,3);
t440 = m(9) * rSges(9,3);
t439 = rSges(3,2) * m(3);
t438 = rSges(6,2) * pkin(8);
t437 = rSges(4,3) * m(4);
t436 = rSges(5,3) * m(5);
t342 = qJD(4) + qJD(2);
t435 = -t342 / 0.2e1;
t434 = pkin(1) * (t450 + qJD(3));
t430 = m(6) * qJD(4);
t361 = -rSges(6,3) - pkin(10);
t429 = rSges(6,1) * t361;
t428 = rSges(6,2) * t361;
t355 = sin(qJ(2));
t427 = t354 * t355;
t359 = cos(qJ(2));
t426 = t354 * t359;
t425 = qJ(4) - qJ(2);
t347 = qJ(4) + qJ(2);
t326 = -pkin(4) * t358 + pkin(1);
t283 = pkin(4) * t427 + t326 * t359;
t424 = qJD(2) * t355;
t423 = qJD(2) * t359;
t421 = qJD(4) * (rSges(6,2) * t443 - Icges(6,4));
t420 = pkin(18) + pkin(15);
t343 = qJD(4) - qJD(2);
t418 = pkin(12) * t452;
t416 = t449 / 0.2e1;
t415 = t447 / 0.2e1;
t408 = rSges(6,1) * t430;
t407 = rSges(6,2) * t430;
t406 = -t443 / 0.2e1;
t393 = t342 * t416;
t392 = t343 * t416;
t331 = qJD(3) + t342;
t391 = -t331 * t447 / 0.2e1;
t390 = t331 * t415;
t332 = qJD(3) - t343;
t389 = t332 * t415;
t387 = 2 * t378;
t386 = -qJ(2) + t399;
t385 = qJ(2) + t399;
t382 = rSges(6,1) * t390;
t381 = pkin(4) * t332 * t406;
t379 = rSges(6,2) * t389;
t353 = sin(qJ(4));
t357 = cos(qJ(4));
t377 = t357 * rSges(6,1) - t353 * rSges(6,2);
t301 = t353 * rSges(6,1) + t357 * rSges(6,2);
t290 = t358 * t355 + t426;
t376 = -t358 * t359 + t427;
t315 = -qJ(2) + t378;
t314 = qJ(2) + t378;
t313 = -qJ(4) + t378;
t312 = qJ(4) + t378;
t375 = t301 * qJD(4);
t305 = -qJ(3) + t315;
t304 = qJ(3) + t314;
t349 = qJ(3) + qJ(2);
t333 = sin(t349);
t334 = cos(t349);
t337 = qJ(3) + t347;
t338 = qJ(3) - t425;
t374 = rSges(6,2) * cos(t337) * t390 + cos(t338) * t379 + sin(t337) * t382 + sin(t338) * t381 + (t333 * (rSges(4,2) * t437 + rSges(9,2) * t440 - Icges(4,6) - Icges(9,6)) - (pkin(4) * t436 + rSges(4,1) * t437 + rSges(9,1) * t440 - Icges(4,5) - Icges(9,5)) * t334) * t344;
t371 = 4 * Icges(6,5);
t369 = 0.2e1 * qJ(2);
t368 = 0.2e1 * qJ(4);
t362 = rSges(6,1) * pkin(8);
t360 = cos(pkin(15));
t356 = sin(pkin(15));
t352 = cos(pkin(16));
t351 = sin(pkin(16));
t350 = t369 + qJ(3);
t339 = pkin(14) - qJ(2) - pkin(15);
t336 = -qJ(2) + t420;
t335 = qJ(2) + t420;
t330 = cos(t346);
t329 = sin(t346);
t327 = 0.2e1 * t349;
t325 = m(5) * rSges(5,1) + t446;
t324 = cos(t339);
t323 = sin(t339);
t321 = 0.4e1 * t429;
t319 = 0.2e1 * t339;
t318 = -qJ(3) + t386;
t317 = qJ(3) + t385;
t316 = -m(5) * rSges(5,2) - t361 * m(6);
t311 = -qJ(4) + t387;
t310 = qJ(4) + t387;
t309 = -qJ(2) + t313;
t308 = qJ(2) + t313;
t307 = -qJ(2) + t312;
t306 = qJ(2) + t312;
t300 = -qJ(4) + t305;
t299 = qJ(4) + t305;
t298 = -qJ(4) + t304;
t297 = qJ(4) + t304;
t296 = 0.2e1 * t313;
t295 = 0.2e1 * t312;
t291 = t377 * qJD(4);
t284 = t375 * t446;
t282 = -pkin(4) * t426 + t326 * t355;
t281 = (-(m(6) * t429 + Icges(6,5)) * t353 - t357 * (t428 * m(6) + Icges(6,6))) * qJD(4);
t278 = t344 * t376;
t277 = t344 * t290;
t276 = t290 * t360 - t356 * t376;
t275 = -t356 * t290 - t360 * t376;
t274 = t326 * t423 + (qJD(3) * t376 + t354 * t424) * pkin(4);
t273 = -t326 * t424 + (qJD(3) * t290 + t354 * t423) * pkin(4);
t272 = -t356 * t282 + t283 * t360;
t271 = t282 * t360 + t356 * t283;
t270 = t281 * t360 + t356 * t284;
t269 = t356 * t281 - t284 * t360;
t268 = t356 * t277 + t278 * t360;
t267 = t277 * t360 - t356 * t278;
t266 = t273 * t356 + t274 * t360;
t265 = t273 * t360 - t356 * t274;
t1 = [(sin(t350) * t434 + t333 * t418) * t293 + (cos(t350) * t434 + t334 * t418) * t322 + cos(t368) * t421 + (-t353 * t408 - t357 * t407) * pkin(8) + sin(t298) * t381 + sin(t300) * t382 + (cos(t299) + cos(t298)) * t379 + t279 + (-rSges(3,1) * t439 + Icges(3,4)) * cos(t369) * t450 + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t319) * t451 + (-rSges(4,1) * t444 - rSges(9,1) * t441 + Icges(4,4) + Icges(9,4)) * cos(t327) * t452 - (cos(t296) + cos(t295)) * t421 / 0.2e1 + (sin(t297) * t391 + sin(t299) * t389 + sin(t306) * t393 + sin(t307) * t392) * rSges(6,1) + (t343 * sin(t308) + t342 * sin(t309)) * pkin(1) * t406 + (((pkin(1) * t453 + rSges(3,1) * m(3)) * t355 + t359 * t439) * t451 + (cos(t312) + cos(t313)) * t407 + (sin(t312) - sin(t313)) * t408) * pkin(12) + ((cos(t300) + cos(t297)) * t391 + (cos(t308) + cos(t307)) * t392 + (cos(t309) + cos(t306)) * t393) * rSges(6,2) + (((t362 - t428) * m(6) - Icges(6,6)) * sin(t311) / 0.2e1 - ((t362 + t428) * m(6) + Icges(6,6)) * sin(t310) / 0.2e1 - ((t321 + 0.4e1 * t438) * m(6) + t371) * cos(t311) / 0.8e1 + ((t321 - 0.4e1 * t438) * m(6) + t371) * cos(t310) / 0.8e1 + (sin(t368) / 0.2e1 + sin(t296) / 0.4e1 - sin(t295) / 0.4e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2))) * qJD(4) + ((-Icges(4,1) - Icges(9,1) + Icges(4,2) + Icges(9,2) + t445 * pkin(4) ^ 2 + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t327) + ((-sin(t305) + sin(t304)) * t325 + (-cos(t305) + cos(t304)) * t316) * pkin(4) + ((cos(t318) + cos(t317)) * rSges(9,2) + (-sin(t318) + sin(t317)) * rSges(9,1)) * t448) * t344 + (0.2e1 * (-rSges(7,1) * t323 + rSges(7,2) * t324) * pkin(6) * m(7) - (-Icges(3,1) + Icges(3,2) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) + t453 * pkin(1) ^ 2) * sin(t369) + ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t319) + (((-cos(t335) + cos(t336)) * rSges(8,2) + (sin(t335) - sin(t336)) * rSges(8,1)) * m(8) + (-sin(t386) + sin(t385)) * t448 + (sin(t314) - sin(t315)) * t325 + (cos(t314) - cos(t315)) * t316) * pkin(1)) * qJD(2); (Icges(3,5) * t359 - t355 * Icges(3,6) + (-rSges(3,1) * t359 + rSges(3,2) * t355) * rSges(3,3) * m(3) + (-rSges(7,2) * t442 + Icges(7,6)) * t323 - (rSges(7,1) * t442 - Icges(7,5)) * t324 + (-m(8) * rSges(8,3) - t436 - t437 - t440) * pkin(1) * t359) * qJD(2) + ((t343 * cos(-t425) / 0.2e1 + cos(t347) * t435) * rSges(6,2) + (-t343 * sin(-t425) / 0.2e1 + sin(t347) * t435) * rSges(6,1)) * t449 + t374; 0.2e1 * t279; t374; t279; 0; (t269 * t352 + t351 * t270) * t330 + (-t351 * t269 + t270 * t352) * t329 + (-t377 * t273 + (pkin(12) + t283) * t375) * m(6); (-((t351 * t265 + t266 * t352) * t330 + t329 * (t265 * t352 - t351 * t266)) * t301 - ((t271 * t352 + t351 * t272) * t330 + t329 * (-t271 * t351 + t272 * t352)) * t291) * m(6); (((t267 * t351 + t268 * t352) * t330 + (t267 * t352 - t351 * t268) * t329) * t301 + ((t275 * t351 + t276 * t352) * t330 + (t275 * t352 - t351 * t276) * t329) * t291) * t447; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
