% Calculate kinetic energy for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1DE1_energykin_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:26:27
% EndTime: 2020-04-18 10:27:14
% DurationCPUTime: 44.66s
% Computational Cost: add. (1168980->242), mult. (1804000->466), div. (76180->22), fcn. (1128070->34), ass. (0->191)
t464 = -2 * pkin(1);
t463 = pkin(5) - pkin(6);
t462 = pkin(5) + pkin(6);
t461 = -pkin(8) - pkin(10);
t460 = pkin(10) - pkin(8);
t395 = pkin(3) ^ 2;
t394 = pkin(4) ^ 2;
t393 = pkin(5) ^ 2;
t398 = pkin(1) ^ 2;
t377 = sin(qJ(2));
t378 = sin(pkin(16));
t382 = cos(qJ(2));
t383 = cos(pkin(16));
t360 = t377 * t378 - t382 * t383;
t446 = pkin(5) * t360;
t423 = t446 * t464 + t398;
t352 = t393 + t423;
t421 = pkin(2) ^ 2 - pkin(6) ^ 2;
t349 = t352 + t421;
t353 = pkin(1) - t446;
t361 = t377 * t383 + t378 * t382;
t346 = (pkin(2) + t462) * (-pkin(2) + t463) + t423;
t347 = (-pkin(2) + t462) * (pkin(2) + t463) + t423;
t400 = sqrt(-t347 * t346);
t340 = pkin(5) * t349 * t361 + t353 * t400;
t376 = sin(qJ(3));
t433 = t340 * t376;
t426 = t361 * t400;
t339 = -pkin(5) * t426 + t349 * t353;
t381 = cos(qJ(3));
t436 = t339 * t381;
t404 = -t433 / 0.2e1 + t436 / 0.2e1;
t350 = 0.1e1 / t352;
t397 = 0.1e1 / pkin(2);
t428 = t350 * t397;
t335 = t404 * t428;
t432 = t340 * t381;
t437 = t339 * t376;
t403 = t432 / 0.2e1 + t437 / 0.2e1;
t336 = t403 * t428;
t368 = pkin(18) + pkin(19);
t365 = sin(t368);
t366 = cos(t368);
t317 = t335 * t366 - t336 * t365;
t447 = pkin(4) * t317;
t424 = -0.2e1 * pkin(3) * t447 + t394;
t312 = t395 + t424;
t310 = 0.1e1 / t312;
t459 = t310 / 0.2e1;
t358 = t361 * qJD(2);
t359 = t360 * qJD(2);
t420 = pkin(1) * pkin(5) * t358;
t431 = 0.2e1 * (t346 + t347) * t420 / t400;
t413 = -t431 / 0.2e1;
t401 = t359 * t400 + t361 * t413;
t327 = ((t353 * t464 - t349) * t358 + t401) * pkin(5);
t458 = -t327 / 0.2e1;
t418 = -0.2e1 * t358 * t361;
t427 = t358 * t400;
t328 = t353 * t431 / 0.2e1 + t393 * pkin(1) * t418 + (-t349 * t359 - t427) * pkin(5);
t457 = t328 / 0.2e1;
t369 = sin(pkin(17));
t456 = t369 / 0.2e1;
t371 = sin(pkin(19));
t455 = t371 / 0.2e1;
t454 = t376 / 0.2e1;
t384 = cos(pkin(15));
t453 = t384 / 0.2e1;
t422 = pkin(8) ^ 2 - pkin(10) ^ 2;
t308 = t312 + t422;
t313 = -pkin(3) + t447;
t316 = -t335 * t365 - t336 * t366;
t306 = (pkin(3) - t461) * (pkin(3) + t461) + t424;
t307 = (pkin(3) - t460) * (pkin(3) + t460) + t424;
t399 = sqrt(-t307 * t306);
t440 = t316 * t399;
t291 = -pkin(4) * t440 - t308 * t313;
t448 = pkin(4) * t316;
t293 = t308 * t448 - t313 * t399;
t390 = 0.1e1 / pkin(8);
t414 = t390 * t459;
t283 = atan2(t293 * t414, t291 * t414);
t452 = sin(t283);
t435 = t340 * t371;
t373 = cos(pkin(19));
t438 = t339 * t373;
t332 = (-t438 / 0.2e1 + t435 / 0.2e1) * t428;
t330 = 0.1e1 / t332 ^ 2;
t434 = t340 * t373;
t439 = t339 * t371;
t333 = (t434 / 0.2e1 + t439 / 0.2e1) * t428;
t411 = 0.1e1 / t352 ^ 2 * t420;
t288 = qJD(2) + (-((t328 * t455 + t373 * t458) * t350 + (t435 - t438) * t411) * t333 * t330 + ((t327 * t455 + t373 * t457) * t350 + (t434 + t439) * t411) / t332) / (t330 * t333 ^ 2 + 0.1e1) * t397;
t451 = pkin(3) * t288;
t299 = ((t432 + t437) * t411 + (t404 * qJD(3) + t327 * t454 + t381 * t457) * t350) * t397;
t300 = ((t433 - t436) * t411 + (t403 * qJD(3) + t328 * t454 + t381 * t458) * t350) * t397;
t296 = -t299 * t365 - t300 * t366;
t450 = pkin(3) * t296;
t449 = pkin(4) * t293;
t445 = pkin(1) * qJD(2);
t419 = pkin(4) * t450;
t444 = 0.2e1 * (t306 + t307) * t419 / t399;
t443 = t296 * t399;
t370 = cos(pkin(17));
t442 = t310 * t370;
t388 = 0.1e1 / pkin(10);
t441 = t310 * t388;
t379 = sin(pkin(15));
t430 = t350 * t379;
t392 = 0.1e1 / pkin(6);
t429 = t350 * t392;
t425 = qJD(2) ^ 2 * t398;
t367 = qJD(2) + qJD(3);
t417 = t376 * t445;
t416 = -t444 / 0.2e1;
t415 = t310 * t456;
t412 = t350 * t453;
t311 = 0.1e1 / t312 ^ 2;
t410 = t311 * t419;
t309 = t312 - t422;
t314 = -pkin(3) * t317 + pkin(4);
t292 = -pkin(3) * t440 + t309 * t314;
t408 = t292 * t410;
t294 = pkin(3) * t309 * t316 + t314 * t399;
t407 = t294 * t410;
t348 = t352 - t421;
t354 = pkin(1) * t360 - pkin(5);
t338 = -pkin(1) * t426 - t348 * t354;
t406 = t338 * t411;
t341 = pkin(1) * t348 * t361 - t354 * t400;
t405 = t341 * t411;
t364 = (-pkin(1) * t382 - pkin(13)) * qJD(1);
t280 = (-t292 * t370 / 0.2e1 + t294 * t456) * t441;
t281 = (t294 * t370 / 0.2e1 + t292 * t456) * t441;
t277 = atan2(t281, t280);
t274 = sin(t277);
t275 = cos(t277);
t356 = (t376 * t377 - t381 * t382) * qJD(1);
t357 = (-t376 * t382 - t377 * t381) * qJD(1);
t265 = -t274 * t357 + t275 * t356;
t362 = pkin(4) * t367 - t381 * t445;
t268 = t274 * t362 - t275 * t417;
t295 = -t299 * t366 + t300 * t365;
t402 = -t295 * t399 + t316 * t416;
t267 = t274 * t417 + t275 * t362;
t345 = -pkin(4) * t356 + t364;
t380 = cos(qJ(4));
t375 = sin(qJ(4));
t374 = cos(pkin(18));
t372 = sin(pkin(18));
t363 = t364 ^ 2;
t337 = (t341 * t453 - t338 * t379 / 0.2e1) * t429;
t334 = (t338 * t453 + t341 * t379 / 0.2e1) * t429;
t331 = 0.1e1 / t334 ^ 2;
t329 = t354 * t413 + t398 * pkin(5) * t418 + (-t348 * t359 - t427) * pkin(1);
t326 = ((0.2e1 * pkin(5) * t354 - t348) * t358 + t401) * pkin(1);
t325 = atan2(t337, t334);
t324 = atan2(t333, t332);
t322 = cos(t325);
t321 = sin(t325);
t319 = cos(t324);
t318 = sin(t324);
t304 = (t318 * t382 + t319 * t377) * qJD(1);
t303 = (-t318 * t377 + t319 * t382) * qJD(1);
t297 = t364 + (-t303 * t374 - t304 * t372) * pkin(3);
t290 = 0.1e1 / t291 ^ 2;
t289 = ((t329 * t412 + t384 * t405 - t326 * t430 / 0.2e1 - t379 * t406) / t334 - (t326 * t412 + t384 * t406 + t329 * t430 / 0.2e1 + t379 * t405) * t337 * t331) / (t331 * t337 ^ 2 + 0.1e1) * t392;
t287 = t319 * t445 + t374 * t451;
t286 = t318 * t445 + t372 * t451;
t282 = cos(t283);
t279 = 0.1e1 / t280 ^ 2;
t273 = -t374 * t282 - t372 * t452;
t272 = t282 * t372 - t374 * t452;
t270 = t314 * t444 / 0.2e1 - 0.2e1 * t395 * t296 * t448 + (t295 * t309 - t443) * pkin(3);
t269 = ((-0.2e1 * pkin(4) * t314 - t309) * t296 + t402) * pkin(3);
t266 = t274 * t356 + t275 * t357;
t264 = qJD(4) - t265;
t263 = t272 * t303 + t273 * t304;
t262 = -t272 * t304 + t273 * t303;
t261 = t272 * t287 + t273 * t286;
t260 = -t272 * t286 + t273 * t287;
t259 = -pkin(9) * t265 - pkin(11) * t266 + t345;
t258 = t288 + 0.2e1 * (((t313 * t416 + (t295 * t308 - t443) * pkin(4)) * t459 + (-t310 * t316 * t394 + t311 * t449) * t450) / t291 - ((-t296 * t308 + t402) * t459 + (t291 * t311 + t310 * t313) * t450) * t290 * t449) * pkin(8) / (t290 * t293 ^ 2 + 0.1e1) * t312 * t390;
t257 = ((t270 * t442 / 0.2e1 + t370 * t407 + t269 * t415 + t369 * t408) / t280 - (-t269 * t442 / 0.2e1 - t370 * t408 + t270 * t415 + t369 * t407) * t281 * t279) / (t279 * t281 ^ 2 + 0.1e1) * t388 + t367;
t256 = pkin(11) * t257 + t268;
t255 = -pkin(9) * t257 - t267;
t254 = t257 * t375 + t266 * t380;
t253 = t257 * t380 - t266 * t375;
t252 = t256 * t380 + t259 * t375;
t251 = -t256 * t375 + t259 * t380;
t1 = m(5) * (t267 ^ 2 + t268 ^ 2 + t345 ^ 2) / 0.2e1 + m(9) * (t260 ^ 2 + t261 ^ 2 + t297 ^ 2) / 0.2e1 + m(6) * (t251 ^ 2 + t252 ^ 2 + t255 ^ 2) / 0.2e1 + Ifges(7,3) * t289 ^ 2 / 0.2e1 + Ifges(4,3) * t367 ^ 2 / 0.2e1 + m(8) * (t363 + (t318 ^ 2 + t319 ^ 2) * t425) / 0.2e1 + m(4) * (t363 + (t376 ^ 2 + t381 ^ 2) * t425) / 0.2e1 + (t364 * mrSges(8,2) + Ifges(8,1) * t304 / 0.2e1) * t304 + (-t364 * mrSges(8,1) + Ifges(8,4) * t304 + Ifges(8,2) * t303 / 0.2e1) * t303 + (t345 * mrSges(5,2) - t267 * mrSges(5,3) + Ifges(5,1) * t266 / 0.2e1) * t266 + (t251 * mrSges(6,1) - t252 * mrSges(6,2) + Ifges(6,3) * t264 / 0.2e1) * t264 + (t297 * mrSges(9,2) - t260 * mrSges(9,3) + Ifges(9,1) * t263 / 0.2e1) * t263 + (t364 * mrSges(4,2) + Ifges(4,5) * t367 + Ifges(4,1) * t357 / 0.2e1) * t357 + (-t345 * mrSges(5,1) + t268 * mrSges(5,3) + Ifges(5,4) * t266 + Ifges(5,2) * t265 / 0.2e1) * t265 + (-t297 * mrSges(9,1) + t261 * mrSges(9,3) + Ifges(9,4) * t263 + Ifges(9,2) * t262 / 0.2e1) * t262 + (t255 * mrSges(6,2) - t251 * mrSges(6,3) + Ifges(6,5) * t264 + Ifges(6,1) * t254 / 0.2e1) * t254 + (-t364 * mrSges(4,1) + Ifges(4,4) * t357 + Ifges(4,6) * t367 + Ifges(4,2) * t356 / 0.2e1) * t356 + (Ifges(8,5) * t304 + Ifges(8,6) * t303 + Ifges(8,3) * t288 / 0.2e1) * t288 + (t260 * mrSges(9,1) - t261 * mrSges(9,2) + Ifges(9,5) * t263 + Ifges(9,6) * t262 + Ifges(9,3) * t258 / 0.2e1) * t258 + (t267 * mrSges(5,1) - t268 * mrSges(5,2) + Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * t257 / 0.2e1) * t257 + (-t255 * mrSges(6,1) + t252 * mrSges(6,3) + Ifges(6,4) * t254 + Ifges(6,6) * t264 + Ifges(6,2) * t253 / 0.2e1) * t253 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t318 * (-mrSges(8,2) * t288 + mrSges(8,3) * t303) + t319 * (mrSges(8,1) * t288 - mrSges(8,3) * t304) - t376 * (-mrSges(4,2) * t367 + mrSges(4,3) * t356) - t381 * (mrSges(4,1) * t367 - mrSges(4,3) * t357)) * pkin(1)) * qJD(2) + (t289 * (Ifges(7,5) * t321 + Ifges(7,6) * t322) + qJD(2) * (Ifges(3,5) * t377 + Ifges(3,6) * t382) + (Ifges(2,3) / 0.2e1 + m(3) * pkin(13) ^ 2 / 0.2e1 + m(7) * pkin(7) ^ 2 / 0.2e1 + (pkin(13) * mrSges(3,1) + Ifges(3,2) * t382 / 0.2e1) * t382 + (-pkin(7) * mrSges(7,1) + Ifges(7,2) * t322 / 0.2e1) * t322 + (-pkin(13) * mrSges(3,2) + Ifges(3,4) * t382 + Ifges(3,1) * t377 / 0.2e1) * t377 + (pkin(7) * mrSges(7,2) + Ifges(7,4) * t322 + Ifges(7,1) * t321 / 0.2e1) * t321) * qJD(1)) * qJD(1);
T = t1;
