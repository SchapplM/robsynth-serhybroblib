% Calculate time derivative of joint inertia matrix for
% palh3m2DE1
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
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE1_inertiaDJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:59:09
% EndTime: 2020-05-07 01:59:14
% DurationCPUTime: 4.01s
% Computational Cost: add. (665->217), mult. (1013->294), div. (0->0), fcn. (392->80), ass. (0->139)
t444 = m(5) + m(6);
t345 = pkin(17) + pkin(18);
t398 = pkin(15) + t345;
t377 = (pkin(16) + t398);
t452 = m(4) + m(8) + m(9) + t444;
t292 = pkin(4) * t444 + m(4) * rSges(4,1) + m(9) * rSges(9,1);
t440 = m(9) * rSges(9,2);
t443 = m(4) * rSges(4,2);
t321 = t440 + t443;
t353 = sin(qJ(3));
t357 = cos(qJ(3));
t278 = (t292 * t353 + t321 * t357) * pkin(1) * qJD(3);
t343 = -qJD(3) - qJD(2);
t451 = -0.2e1 * t343;
t450 = -0.2e1 * qJD(2);
t449 = 0.2e1 * qJD(2);
t448 = pkin(1) * m(6);
t447 = pkin(3) * m(9);
t446 = pkin(4) * m(6);
t445 = pkin(8) * m(6);
t442 = m(6) * rSges(6,1);
t441 = m(7) * rSges(7,3);
t439 = m(9) * rSges(9,3);
t438 = rSges(3,2) * m(3);
t437 = rSges(6,2) * pkin(8);
t436 = rSges(4,3) * m(4);
t435 = rSges(5,3) * m(5);
t341 = qJD(4) + qJD(2);
t434 = -t341 / 0.2e1;
t433 = pkin(1) * (t449 + qJD(3));
t429 = m(6) * qJD(4);
t360 = -rSges(6,3) - pkin(10);
t428 = rSges(6,2) * t360;
t427 = t360 * rSges(6,1);
t354 = sin(qJ(2));
t426 = t353 * t354;
t358 = cos(qJ(2));
t425 = t353 * t358;
t346 = qJ(3) + qJ(2);
t424 = qJ(4) - qJ(2);
t325 = -pkin(4) * t357 + pkin(1);
t282 = pkin(4) * t426 + t325 * t358;
t423 = qJD(2) * t354;
t422 = qJD(2) * t358;
t420 = qJD(4) * (rSges(6,2) * t442 - Icges(6,4));
t419 = pkin(18) + pkin(15);
t342 = qJD(4) - qJD(2);
t417 = pkin(12) * t451;
t415 = t448 / 0.2e1;
t414 = t446 / 0.2e1;
t407 = rSges(6,1) * t429;
t406 = rSges(6,2) * t429;
t405 = -t442 / 0.2e1;
t392 = t341 * t415;
t391 = t342 * t415;
t330 = qJD(3) + t341;
t390 = -t330 * t446 / 0.2e1;
t389 = t330 * t414;
t331 = qJD(3) - t342;
t388 = t331 * t414;
t386 = 2 * t377;
t385 = -qJ(2) + t398;
t384 = qJ(2) + t398;
t381 = rSges(6,1) * t389;
t380 = pkin(4) * t331 * t405;
t378 = rSges(6,2) * t388;
t352 = sin(qJ(4));
t356 = cos(qJ(4));
t376 = rSges(6,1) * t356 - rSges(6,2) * t352;
t300 = rSges(6,1) * t352 + rSges(6,2) * t356;
t289 = t354 * t357 + t425;
t375 = -t357 * t358 + t426;
t314 = -qJ(2) + t377;
t313 = qJ(2) + t377;
t310 = -qJ(4) + t377;
t309 = qJ(4) + t377;
t374 = t300 * qJD(4);
t306 = -qJ(2) + t310;
t305 = qJ(2) + t310;
t304 = -qJ(2) + t309;
t303 = qJ(2) + t309;
t332 = sin(t346);
t333 = cos(t346);
t336 = qJ(4) + t346;
t337 = qJ(3) - t424;
t373 = rSges(6,2) * cos(t336) * t389 + cos(t337) * t378 + sin(t336) * t381 + sin(t337) * t380 + (t332 * (rSges(4,2) * t436 + rSges(9,2) * t439 - Icges(4,6) - Icges(9,6)) - (pkin(4) * t435 + rSges(4,1) * t436 + rSges(9,1) * t439 - Icges(4,5) - Icges(9,5)) * t333) * t343;
t370 = 4 * Icges(6,5);
t368 = 0.2e1 * qJ(2);
t367 = 0.2e1 * qJ(4);
t361 = rSges(6,1) * pkin(8);
t359 = cos(pkin(15));
t355 = sin(pkin(15));
t351 = cos(pkin(16));
t350 = sin(pkin(16));
t348 = qJ(2) + qJ(4);
t347 = t368 + qJ(3);
t338 = pkin(14) - qJ(2) - pkin(15);
t335 = -qJ(2) + t419;
t334 = qJ(2) + t419;
t329 = cos(t345);
t328 = sin(t345);
t326 = 0.2e1 * t346;
t324 = m(5) * rSges(5,1) + t445;
t323 = cos(t338);
t322 = sin(t338);
t320 = 0.4e1 * t427;
t318 = 0.2e1 * t338;
t317 = -qJ(3) + t385;
t316 = qJ(3) + t384;
t315 = -m(5) * rSges(5,2) - m(6) * t360;
t312 = -qJ(4) + t386;
t311 = qJ(4) + t386;
t308 = -qJ(3) + t314;
t307 = qJ(3) + t313;
t299 = -qJ(3) + t306;
t298 = -qJ(3) + t304;
t297 = qJ(3) + t305;
t296 = qJ(3) + t303;
t295 = 0.2e1 * t310;
t294 = 0.2e1 * t309;
t290 = t376 * qJD(4);
t283 = t374 * t445;
t281 = -pkin(4) * t425 + t325 * t354;
t280 = (-(m(6) * t427 + Icges(6,5)) * t352 - (t428 * m(6) + Icges(6,6)) * t356) * qJD(4);
t277 = t343 * t375;
t276 = t343 * t289;
t275 = t289 * t359 - t355 * t375;
t274 = -t289 * t355 - t359 * t375;
t273 = t325 * t422 + (qJD(3) * t375 + t353 * t423) * pkin(4);
t272 = -t325 * t423 + (qJD(3) * t289 + t353 * t422) * pkin(4);
t271 = -t281 * t355 + t282 * t359;
t270 = t281 * t359 + t282 * t355;
t269 = t280 * t359 + t283 * t355;
t268 = t280 * t355 - t283 * t359;
t267 = t276 * t355 + t277 * t359;
t266 = t276 * t359 - t277 * t355;
t265 = t272 * t355 + t273 * t359;
t264 = t272 * t359 - t273 * t355;
t1 = [t278 + (-rSges(3,1) * t438 + Icges(3,4)) * cos(t368) * t449 + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t318) * t450 + (-rSges(4,1) * t443 - rSges(9,1) * t440 + Icges(4,4) + Icges(9,4)) * cos(t326) * t451 + (t341 * sin(t306) + t342 * sin(t305)) * pkin(1) * t405 - (cos(t295) + cos(t294)) * t420 / 0.2e1 + (cos(t347) * t433 + t333 * t417) * t321 + (sin(t296) * t390 + sin(t298) * t388 + sin(t303) * t392 + sin(t304) * t391) * rSges(6,1) + (sin(t347) * t433 + t332 * t417) * t292 + (cos(t298) + cos(t297)) * t378 + (-t352 * t407 - t356 * t406) * pkin(8) + cos(t367) * t420 + sin(t297) * t380 + sin(t299) * t381 + ((cos(t310) + cos(t309)) * t406 + (sin(t309) - sin(t310)) * t407 + ((pkin(1) * t452 + rSges(3,1) * m(3)) * t354 + t358 * t438) * t450) * pkin(12) + ((cos(t305) + cos(t304)) * t391 + (cos(t306) + cos(t303)) * t392 + (cos(t299) + cos(t296)) * t390) * rSges(6,2) + ((sin(t367) / 0.2e1 + sin(t295) / 0.4e1 - sin(t294) / 0.4e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2)) - ((t320 + 0.4e1 * t437) * m(6) + t370) * cos(t312) / 0.8e1 + ((t320 - 0.4e1 * t437) * m(6) + t370) * cos(t311) / 0.8e1 + ((t361 - t428) * m(6) - Icges(6,6)) * sin(t312) / 0.2e1 - ((t361 + t428) * m(6) + Icges(6,6)) * sin(t311) / 0.2e1) * qJD(4) + ((-Icges(4,1) - Icges(9,1) + Icges(4,2) + Icges(9,2) + t444 * pkin(4) ^ 2 + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t326) + ((-sin(t308) + sin(t307)) * t324 + (-cos(t308) + cos(t307)) * t315) * pkin(4) + ((cos(t317) + cos(t316)) * rSges(9,2) + (sin(t316) - sin(t317)) * rSges(9,1)) * t447) * t343 + (((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t318) - (-Icges(3,1) + Icges(3,2) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) + t452 * pkin(1) ^ 2) * sin(t368) + 0.2e1 * (-rSges(7,1) * t322 + rSges(7,2) * t323) * pkin(6) * m(7) + (((-cos(t334) + cos(t335)) * rSges(8,2) + (sin(t334) - sin(t335)) * rSges(8,1)) * m(8) + (-sin(t385) + sin(t384)) * t447 + (sin(t313) - sin(t314)) * t324 + (cos(t313) - cos(t314)) * t315) * pkin(1)) * qJD(2); (Icges(3,5) * t358 - t354 * Icges(3,6) + (-rSges(3,1) * t358 + rSges(3,2) * t354) * rSges(3,3) * m(3) + (-rSges(7,2) * t441 + Icges(7,6)) * t322 - (rSges(7,1) * t441 - Icges(7,5)) * t323 + (-m(8) * rSges(8,3) - t435 - t436 - t439) * pkin(1) * t358) * qJD(2) + ((t342 * cos(-t424) / 0.2e1 + cos(t348) * t434) * rSges(6,2) + (-t342 * sin(-t424) / 0.2e1 + sin(t348) * t434) * rSges(6,1)) * t448 + t373; 0.2e1 * t278; t373; t278; 0; (t268 * t351 + t269 * t350) * t329 + (-t268 * t350 + t269 * t351) * t328 + (-t376 * t272 + (pkin(12) + t282) * t374) * m(6); (-((t264 * t350 + t265 * t351) * t329 + (t264 * t351 - t265 * t350) * t328) * t300 - ((t270 * t351 + t271 * t350) * t329 + (-t270 * t350 + t271 * t351) * t328) * t290) * m(6); (((t266 * t350 + t267 * t351) * t329 + (t266 * t351 - t267 * t350) * t328) * t300 + ((t274 * t350 + t275 * t351) * t329 + (t274 * t351 - t275 * t350) * t328) * t290) * t446; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
