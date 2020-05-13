% Calculate kinetic energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1TE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1TE_energykin_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:19:41
% EndTime: 2020-04-17 15:20:25
% DurationCPUTime: 37.61s
% Computational Cost: add. (984330->242), mult. (1526696->467), div. (62972->22), fcn. (952374->22), ass. (0->183)
t458 = -2 * pkin(1);
t457 = pkin(5) - pkin(6);
t456 = pkin(5) + pkin(6);
t455 = -pkin(8) - pkin(10);
t454 = pkin(10) - pkin(8);
t385 = pkin(3) ^ 2;
t384 = pkin(4) ^ 2;
t383 = pkin(5) ^ 2;
t388 = pkin(1) ^ 2;
t367 = sin(qJ(2));
t368 = sin(pkin(16));
t372 = cos(qJ(2));
t373 = cos(pkin(16));
t350 = t367 * t368 - t372 * t373;
t439 = pkin(5) * t350;
t415 = t439 * t458 + t388;
t342 = t383 + t415;
t413 = pkin(2) ^ 2 - pkin(6) ^ 2;
t339 = t342 + t413;
t343 = pkin(1) - t439;
t351 = t367 * t373 + t368 * t372;
t336 = (pkin(2) + t456) * (-pkin(2) + t457) + t415;
t337 = (-pkin(2) + t456) * (pkin(2) + t457) + t415;
t390 = sqrt(-t337 * t336);
t330 = pkin(5) * t339 * t351 + t343 * t390;
t366 = sin(qJ(3));
t425 = t330 * t366;
t418 = t351 * t390;
t329 = -pkin(5) * t418 + t339 * t343;
t371 = cos(qJ(3));
t428 = t329 * t371;
t397 = -t425 / 0.2e1 + t428 / 0.2e1;
t340 = 0.1e1 / t342;
t387 = 0.1e1 / pkin(2);
t420 = t340 * t387;
t325 = t397 * t420;
t424 = t330 * t371;
t429 = t329 * t366;
t396 = t424 / 0.2e1 + t429 / 0.2e1;
t326 = t396 * t420;
t358 = pkin(18) + pkin(19);
t355 = sin(t358);
t356 = cos(t358);
t310 = t325 * t356 - t326 * t355;
t440 = pkin(4) * t310;
t416 = -0.2e1 * pkin(3) * t440 + t384;
t305 = t385 + t416;
t414 = pkin(8) ^ 2 - pkin(10) ^ 2;
t301 = t305 + t414;
t306 = -pkin(3) + t440;
t309 = -t325 * t355 - t326 * t356;
t299 = (pkin(3) - t455) * (pkin(3) + t455) + t416;
t300 = (pkin(3) - t454) * (pkin(3) + t454) + t416;
t389 = sqrt(-t300 * t299);
t432 = t309 * t389;
t288 = -pkin(4) * t432 - t301 * t306;
t453 = t288 / 0.2e1;
t303 = 0.1e1 / t305;
t452 = t303 / 0.2e1;
t348 = t351 * qJD(2);
t349 = t350 * qJD(2);
t412 = pkin(1) * pkin(5) * t348;
t423 = 0.2e1 * (t336 + t337) * t412 / t390;
t406 = -t423 / 0.2e1;
t394 = t349 * t390 + t351 * t406;
t316 = ((t343 * t458 - t339) * t348 + t394) * pkin(5);
t451 = -t316 / 0.2e1;
t410 = -0.2e1 * t348 * t351;
t419 = t348 * t390;
t317 = t343 * t423 / 0.2e1 + t383 * pkin(1) * t410 + (-t339 * t349 - t419) * pkin(5);
t450 = t317 / 0.2e1;
t359 = sin(pkin(17));
t449 = t359 / 0.2e1;
t361 = sin(pkin(19));
t448 = t361 / 0.2e1;
t364 = cos(pkin(18));
t447 = -t364 / 0.2e1;
t446 = t366 / 0.2e1;
t374 = cos(pkin(15));
t445 = t374 / 0.2e1;
t427 = t330 * t361;
t363 = cos(pkin(19));
t430 = t329 * t363;
t322 = (-t430 / 0.2e1 + t427 / 0.2e1) * t420;
t391 = t322 ^ 2;
t319 = 0.1e1 / t391;
t426 = t330 * t363;
t431 = t329 * t361;
t323 = (t426 / 0.2e1 + t431 / 0.2e1) * t420;
t320 = t323 ^ 2;
t404 = 0.1e1 / t342 ^ 2 * t412;
t281 = qJD(2) + (-((t317 * t448 + t363 * t451) * t340 + (t427 - t430) * t404) * t323 * t319 + ((t316 * t448 + t363 * t450) * t340 + (t426 + t431) * t404) / t322) / (t319 * t320 + 0.1e1) * t387;
t444 = pkin(3) * t281;
t293 = ((t424 + t429) * t404 + (t397 * qJD(3) + t316 * t446 + t371 * t450) * t340) * t387;
t294 = ((t425 - t428) * t404 + (t396 * qJD(3) + t317 * t446 + t371 * t451) * t340) * t387;
t291 = -t293 * t355 - t294 * t356;
t443 = pkin(3) * t291;
t441 = pkin(4) * t309;
t289 = t301 * t441 - t306 * t389;
t442 = pkin(4) * t289;
t438 = pkin(1) * qJD(2);
t411 = pkin(4) * t443;
t437 = 0.2e1 * (t299 + t300) * t411 / t389;
t436 = t291 * t389;
t360 = cos(pkin(17));
t435 = t303 * t360;
t378 = 0.1e1 / pkin(10);
t434 = t303 * t378;
t380 = 0.1e1 / pkin(8);
t433 = t303 * t380;
t369 = sin(pkin(15));
t422 = t340 * t369;
t382 = 0.1e1 / pkin(6);
t421 = t340 * t382;
t417 = qJD(2) ^ 2 * t388;
t357 = qJD(2) + qJD(3);
t409 = t366 * t438;
t408 = -t437 / 0.2e1;
t407 = t303 * t449;
t405 = t340 * t445;
t304 = 0.1e1 / t305 ^ 2;
t403 = t304 * t411;
t338 = t342 - t413;
t344 = pkin(1) * t350 - pkin(5);
t328 = -pkin(1) * t418 - t338 * t344;
t401 = t328 * t404;
t331 = pkin(1) * t338 * t351 - t344 * t390;
t400 = t331 * t404;
t354 = (-pkin(1) * t372 - pkin(13)) * qJD(1);
t302 = t305 - t414;
t307 = -pkin(3) * t310 + pkin(4);
t398 = pkin(3) * t302 * t309 + t307 * t389;
t399 = -pkin(3) * t432 + t302 * t307;
t273 = (-t399 * t360 / 0.2e1 + t398 * t449) * t434;
t274 = (t398 * t360 / 0.2e1 + t399 * t449) * t434;
t346 = (t366 * t367 - t371 * t372) * qJD(1);
t347 = (-t366 * t372 - t367 * t371) * qJD(1);
t267 = t273 * t346 - t274 * t347;
t352 = pkin(4) * t357 - t371 * t438;
t268 = -t273 * t409 + t274 * t352;
t290 = -t293 * t356 + t294 * t355;
t395 = -t290 * t389 + t309 * t408;
t269 = t273 * t352 + t274 * t409;
t335 = -pkin(4) * t346 + t354;
t393 = t399 * t403;
t392 = t398 * t403;
t370 = cos(qJ(4));
t365 = sin(qJ(4));
t362 = sin(pkin(18));
t353 = t354 ^ 2;
t327 = (t331 * t445 - t328 * t369 / 0.2e1) * t421;
t324 = (t328 * t445 + t331 * t369 / 0.2e1) * t421;
t321 = 0.1e1 / t324 ^ 2;
t318 = t344 * t406 + t388 * pkin(5) * t410 + (-t338 * t349 - t419) * pkin(1);
t315 = ((0.2e1 * pkin(5) * t344 - t338) * t348 + t394) * pkin(1);
t312 = (t322 * t372 - t323 * t367) * qJD(1);
t311 = (t322 * t367 + t323 * t372) * qJD(1);
t297 = t354 + (-t311 * t362 - t312 * t364) * pkin(3);
t287 = 0.1e1 / t288 ^ 2;
t282 = ((t318 * t405 + t374 * t400 - t315 * t422 / 0.2e1 - t369 * t401) / t324 - (t315 * t405 + t374 * t401 + t318 * t422 / 0.2e1 + t369 * t400) * t327 * t321) / (t321 * t327 ^ 2 + 0.1e1) * t382;
t280 = t323 * t438 + t362 * t444;
t279 = t322 * t438 + t364 * t444;
t276 = (t288 * t447 - t362 * t289 / 0.2e1) * t433;
t275 = (t289 * t447 + t362 * t453) * t433;
t272 = 0.1e1 / t273 ^ 2;
t266 = t273 * t347 + t274 * t346;
t265 = qJD(4) - t267;
t264 = -t275 * t311 + t276 * t312;
t263 = t275 * t312 + t276 * t311;
t262 = t307 * t437 / 0.2e1 - 0.2e1 * t385 * t291 * t441 + (t290 * t302 - t436) * pkin(3);
t261 = ((-0.2e1 * pkin(4) * t307 - t302) * t291 + t395) * pkin(3);
t260 = -t275 * t280 + t276 * t279;
t259 = t275 * t279 + t276 * t280;
t258 = -pkin(9) * t267 - pkin(11) * t266 + t335;
t257 = t281 + (((t306 * t408 + (t290 * t301 - t436) * pkin(4)) * t452 + (-t303 * t309 * t384 + t304 * t442) * t443) / t453 - 0.2e1 * ((-t291 * t301 + t395) * t452 + (t288 * t304 + t303 * t306) * t443) * t287 * t442) * pkin(8) / (t287 * t289 ^ 2 + 0.1e1) * t305 * t380;
t256 = ((t262 * t435 / 0.2e1 + t360 * t392 + t261 * t407 + t359 * t393) / t273 - (-t261 * t435 / 0.2e1 - t360 * t393 + t262 * t407 + t359 * t392) * t274 * t272) / (t272 * t274 ^ 2 + 0.1e1) * t378 + t357;
t255 = -pkin(9) * t256 - t269;
t254 = pkin(11) * t256 + t268;
t253 = t256 * t365 + t266 * t370;
t252 = t256 * t370 - t266 * t365;
t251 = t254 * t370 + t258 * t365;
t250 = -t254 * t365 + t258 * t370;
t1 = Ifges(7,3) * t282 ^ 2 / 0.2e1 + m(6) * (t250 ^ 2 + t251 ^ 2 + t255 ^ 2) / 0.2e1 + Ifges(4,3) * t357 ^ 2 / 0.2e1 + m(8) * (t353 + (t320 + t391) * t417) / 0.2e1 + m(4) * (t353 + (t366 ^ 2 + t371 ^ 2) * t417) / 0.2e1 + m(5) * (t268 ^ 2 + t269 ^ 2 + t335 ^ 2) / 0.2e1 + m(9) * (t259 ^ 2 + t260 ^ 2 + t297 ^ 2) / 0.2e1 + (-t354 * mrSges(8,1) + Ifges(8,2) * t312 / 0.2e1) * t312 + (t354 * mrSges(8,2) + Ifges(8,4) * t312 + Ifges(8,1) * t311 / 0.2e1) * t311 + (-t335 * mrSges(5,1) + t268 * mrSges(5,3) + Ifges(5,2) * t267 / 0.2e1) * t267 + (t250 * mrSges(6,1) - t251 * mrSges(6,2) + Ifges(6,3) * t265 / 0.2e1) * t265 + (-t297 * mrSges(9,1) + t259 * mrSges(9,3) + Ifges(9,2) * t264 / 0.2e1) * t264 + (t354 * mrSges(4,2) + Ifges(4,5) * t357 + Ifges(4,1) * t347 / 0.2e1) * t347 + (t335 * mrSges(5,2) - t269 * mrSges(5,3) + Ifges(5,4) * t267 + Ifges(5,1) * t266 / 0.2e1) * t266 + (t297 * mrSges(9,2) - t260 * mrSges(9,3) + Ifges(9,4) * t264 + Ifges(9,1) * t263 / 0.2e1) * t263 + (t255 * mrSges(6,2) - t250 * mrSges(6,3) + Ifges(6,5) * t265 + Ifges(6,1) * t253 / 0.2e1) * t253 + (-t354 * mrSges(4,1) + Ifges(4,4) * t347 + Ifges(4,6) * t357 + Ifges(4,2) * t346 / 0.2e1) * t346 + (Ifges(8,5) * t311 + Ifges(8,6) * t312 + Ifges(8,3) * t281 / 0.2e1) * t281 + (t260 * mrSges(9,1) - t259 * mrSges(9,2) + Ifges(9,5) * t263 + Ifges(9,6) * t264 + Ifges(9,3) * t257 / 0.2e1) * t257 + (t269 * mrSges(5,1) - t268 * mrSges(5,2) + Ifges(5,5) * t266 + Ifges(5,6) * t267 + Ifges(5,3) * t256 / 0.2e1) * t256 + (-t255 * mrSges(6,1) + t251 * mrSges(6,3) + Ifges(6,4) * t253 + Ifges(6,6) * t265 + Ifges(6,2) * t252 / 0.2e1) * t252 + (Ifges(3,3) * qJD(2) / 0.2e1 + (-t366 * (-mrSges(4,2) * t357 + mrSges(4,3) * t346) - t371 * (mrSges(4,1) * t357 - mrSges(4,3) * t347) + t322 * (mrSges(8,1) * t281 - mrSges(8,3) * t311) + t323 * (-mrSges(8,2) * t281 + mrSges(8,3) * t312)) * pkin(1)) * qJD(2) + (t282 * (Ifges(7,5) * t327 + Ifges(7,6) * t324) + qJD(2) * (Ifges(3,5) * t367 + Ifges(3,6) * t372) + (Ifges(2,3) / 0.2e1 + m(3) * pkin(13) ^ 2 / 0.2e1 + m(7) * pkin(7) ^ 2 / 0.2e1 + (pkin(13) * mrSges(3,1) + Ifges(3,2) * t372 / 0.2e1) * t372 + (pkin(7) * mrSges(7,2) + Ifges(7,1) * t327 / 0.2e1) * t327 + (-pkin(13) * mrSges(3,2) + Ifges(3,4) * t372 + Ifges(3,1) * t367 / 0.2e1) * t367 + (-pkin(7) * mrSges(7,1) + Ifges(7,4) * t327 + Ifges(7,2) * t324 / 0.2e1) * t324) * qJD(1)) * qJD(1);
T = t1;
