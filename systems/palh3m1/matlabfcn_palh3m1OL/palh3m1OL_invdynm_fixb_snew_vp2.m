% Calculate vector of cutting torques with Newton-Euler for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% m [3x9]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = palh3m1OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1OL_invdynm_fixb_snew_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:04:34
% EndTime: 2020-04-20 17:05:20
% DurationCPUTime: 9.14s
% Computational Cost: add. (116065->491), mult. (255505->634), div. (0->0), fcn. (196591->18), ass. (0->199)
t411 = sin(qJ(1));
t419 = cos(qJ(1));
t384 = -t419 * g(1) - t411 * g(2);
t420 = qJD(1) ^ 2;
t370 = -t420 * pkin(12) + t384;
t410 = sin(qJ(2));
t418 = cos(qJ(2));
t347 = -t418 * g(3) - t410 * t370;
t349 = -t410 * g(3) + t418 * t370;
t440 = qJD(1) * qJD(2);
t375 = t410 * qJDD(1) + t418 * t440;
t438 = t410 * t440;
t377 = t418 * qJDD(1) - t438;
t340 = (t410 * t418 * t420 + qJDD(2)) * pkin(1) + t347;
t341 = (-t418 ^ 2 * t420 - qJD(2) ^ 2) * pkin(1) + t349;
t409 = sin(qJ(3));
t417 = cos(qJ(3));
t298 = -t417 * t340 + t409 * t341;
t300 = -t409 * t340 - t417 * t341;
t361 = (-t409 * t418 - t410 * t417) * qJD(1);
t316 = -t361 * qJD(3) + t409 * t375 - t417 * t377;
t359 = (t409 * t410 - t417 * t418) * qJD(1);
t318 = t359 * qJD(3) - t417 * t375 - t409 * t377;
t401 = qJD(2) + qJD(3);
t324 = Ifges(4,4) * t361 + Ifges(4,2) * t359 + Ifges(4,6) * t401;
t326 = Ifges(4,1) * t361 + Ifges(4,4) * t359 + Ifges(4,5) * t401;
t399 = qJDD(2) + qJDD(3);
t408 = sin(qJ(4));
t416 = cos(qJ(4));
t332 = t408 * t359 + t416 * t361;
t260 = -t332 * qJD(4) + t416 * t316 - t408 * t318;
t331 = t416 * t359 - t408 * t361;
t261 = t331 * qJD(4) + t408 * t316 + t416 * t318;
t383 = t411 * g(1) - t419 * g(2);
t368 = -qJDD(1) * pkin(12) - t383;
t337 = t368 + (-t377 + t438) * pkin(1);
t279 = t337 + (t361 * t401 - t316) * pkin(4);
t395 = qJD(4) + t401;
t230 = (-t331 * t395 - t261) * pkin(10) + (t332 * t395 - t260) * pkin(8) + t279;
t274 = (t359 * t361 + t399) * pkin(4) + t298;
t283 = (-t359 ^ 2 - t401 ^ 2) * pkin(4) + t300;
t241 = t408 * t274 + t416 * t283;
t290 = -t331 * pkin(8) - t332 * pkin(10);
t391 = t395 ^ 2;
t393 = qJDD(4) + t399;
t232 = -t391 * pkin(8) + t393 * pkin(10) + t331 * t290 + t241;
t407 = sin(qJ(5));
t415 = cos(qJ(5));
t222 = t407 * t230 + t415 * t232;
t240 = t416 * t274 - t408 * t283;
t231 = -t393 * pkin(8) - t391 * pkin(10) + t332 * t290 - t240;
t308 = t415 * t332 + t407 * t395;
t237 = -t308 * qJD(5) - t407 * t261 + t415 * t393;
t307 = -t407 * t332 + t415 * t395;
t238 = t307 * qJD(5) + t415 * t261 + t407 * t393;
t328 = qJD(5) - t331;
t251 = Ifges(6,5) * t308 + Ifges(6,6) * t307 + Ifges(6,3) * t328;
t253 = Ifges(6,1) * t308 + Ifges(6,4) * t307 + Ifges(6,5) * t328;
t258 = qJDD(5) - t260;
t200 = -mrSges(6,1) * t231 + mrSges(6,3) * t222 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t258 - t308 * t251 + t328 * t253;
t221 = t415 * t230 - t407 * t232;
t252 = Ifges(6,4) * t308 + Ifges(6,2) * t307 + Ifges(6,6) * t328;
t201 = mrSges(6,2) * t231 - mrSges(6,3) * t221 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t258 + t307 * t251 - t328 * t252;
t285 = Ifges(5,4) * t332 + Ifges(5,2) * t331 + Ifges(5,6) * t395;
t286 = Ifges(5,1) * t332 + Ifges(5,4) * t331 + Ifges(5,5) * t395;
t280 = -t328 * mrSges(6,2) + t307 * mrSges(6,3);
t281 = t328 * mrSges(6,1) - t308 * mrSges(6,3);
t431 = -m(6) * t231 + t237 * mrSges(6,1) - t238 * mrSges(6,2) + t307 * t280 - t308 * t281;
t277 = -t307 * mrSges(6,1) + t308 * mrSges(6,2);
t215 = m(6) * t221 + t258 * mrSges(6,1) - t238 * mrSges(6,3) - t308 * t277 + t328 * t280;
t216 = m(6) * t222 - t258 * mrSges(6,2) + t237 * mrSges(6,3) + t307 * t277 - t328 * t281;
t437 = -t407 * t215 + t415 * t216;
t427 = -pkin(8) * t431 - pkin(10) * t437 - mrSges(5,1) * t240 + mrSges(5,2) * t241 - Ifges(5,5) * t261 - Ifges(5,6) * t260 - Ifges(5,3) * t393 - t415 * t200 - t407 * t201 - t332 * t285 + t331 * t286;
t289 = -t331 * mrSges(5,1) + t332 * mrSges(5,2);
t320 = t395 * mrSges(5,1) - t332 * mrSges(5,3);
t194 = m(5) * t241 - t393 * mrSges(5,2) + t260 * mrSges(5,3) + t331 * t289 - t395 * t320 + t437;
t319 = -t395 * mrSges(5,2) + t331 * mrSges(5,3);
t207 = m(5) * t240 + t393 * mrSges(5,1) - t261 * mrSges(5,3) - t332 * t289 + t395 * t319 + t431;
t447 = t408 * t194 + t416 * t207;
t424 = -pkin(4) * t447 - mrSges(4,1) * t298 + mrSges(4,2) * t300 - Ifges(4,5) * t318 - Ifges(4,6) * t316 - Ifges(4,3) * t399 - t361 * t324 + t359 * t326 + t427;
t405 = sin(qJ(7));
t413 = cos(qJ(7));
t297 = t413 * t340 - t405 * t341;
t299 = t405 * t340 + t413 * t341;
t360 = (t405 * t418 + t410 * t413) * qJD(1);
t315 = -t360 * qJD(7) - t405 * t375 + t413 * t377;
t358 = (-t405 * t410 + t413 * t418) * qJD(1);
t317 = t358 * qJD(7) + t413 * t375 + t405 * t377;
t400 = qJD(2) + qJD(7);
t323 = Ifges(8,4) * t360 + Ifges(8,2) * t358 + Ifges(8,6) * t400;
t325 = Ifges(8,1) * t360 + Ifges(8,4) * t358 + Ifges(8,5) * t400;
t398 = qJDD(2) + qJDD(7);
t403 = sin(pkin(15));
t404 = cos(pkin(15));
t327 = (-t358 * t404 - t360 * t403) * pkin(3);
t397 = t400 ^ 2;
t262 = -t360 * t327 + (t397 * t403 + t398 * t404) * pkin(3) + t297;
t263 = t358 * t327 + (-t397 * t404 + t398 * t403) * pkin(3) + t299;
t412 = cos(qJ(8));
t451 = sin(qJ(8));
t366 = t403 * t412 - t404 * t451;
t367 = -t403 * t451 - t404 * t412;
t234 = t367 * t262 - t366 * t263;
t235 = t366 * t262 + t367 * t263;
t304 = t366 * t358 + t367 * t360;
t248 = -t304 * qJD(8) + t367 * t315 - t366 * t317;
t303 = t367 * t358 - t366 * t360;
t249 = t303 * qJD(8) + t366 * t315 + t367 * t317;
t394 = qJD(8) + t400;
t265 = Ifges(9,4) * t304 + Ifges(9,2) * t303 + Ifges(9,6) * t394;
t266 = Ifges(9,1) * t304 + Ifges(9,4) * t303 + Ifges(9,5) * t394;
t392 = qJDD(8) + t398;
t430 = -mrSges(9,1) * t234 + mrSges(9,2) * t235 - Ifges(9,5) * t249 - Ifges(9,6) * t248 - Ifges(9,3) * t392 - t304 * t265 + t303 * t266;
t270 = -t303 * mrSges(9,1) + t304 * mrSges(9,2);
t301 = -t394 * mrSges(9,2) + t303 * mrSges(9,3);
t228 = m(9) * t234 + t392 * mrSges(9,1) - t249 * mrSges(9,3) - t304 * t270 + t394 * t301;
t302 = t394 * mrSges(9,1) - t304 * mrSges(9,3);
t229 = m(9) * t235 - t392 * mrSges(9,2) + t248 * mrSges(9,3) + t303 * t270 - t394 * t302;
t445 = -t366 * t228 + t367 * t229;
t446 = t367 * t228 + t366 * t229;
t449 = pkin(3) * t404;
t450 = pkin(3) * t403;
t425 = -mrSges(8,1) * t297 + mrSges(8,2) * t299 - Ifges(8,5) * t317 - Ifges(8,6) * t315 - Ifges(8,3) * t398 - t360 * t323 + t358 * t325 - t445 * t450 - t446 * t449 + t430;
t336 = -t359 * mrSges(4,1) + t361 * mrSges(4,2);
t343 = -t401 * mrSges(4,2) + t359 * mrSges(4,3);
t184 = m(4) * t298 + t399 * mrSges(4,1) - t318 * mrSges(4,3) - t361 * t336 + t401 * t343 + t447;
t345 = t401 * mrSges(4,1) - t361 * mrSges(4,3);
t185 = m(4) * t300 - t399 * mrSges(4,2) + t316 * mrSges(4,3) + t416 * t194 - t408 * t207 + t359 * t336 - t401 * t345;
t335 = -t358 * mrSges(8,1) + t360 * mrSges(8,2);
t342 = -t400 * mrSges(8,2) + t358 * mrSges(8,3);
t204 = m(8) * t297 + t398 * mrSges(8,1) - t317 * mrSges(8,3) - t360 * t335 + t400 * t342 + t446;
t344 = t400 * mrSges(8,1) - t360 * mrSges(8,3);
t205 = m(8) * t299 - t398 * mrSges(8,2) + t315 * mrSges(8,3) + t358 * t335 - t400 * t344 + t445;
t448 = -t417 * t184 - t409 * t185 + t413 * t204 + t405 * t205;
t452 = -t448 * pkin(1) - mrSges(3,1) * t347 + mrSges(3,2) * t349 - Ifges(3,5) * t375 - Ifges(3,6) * t377 - Ifges(3,3) * qJDD(2) + t424 + t425;
t196 = t415 * t215 + t407 * t216;
t406 = sin(qJ(6));
t444 = qJD(1) * t406;
t443 = qJD(1) * t410;
t414 = cos(qJ(6));
t442 = qJD(1) * t414;
t441 = qJD(1) * t418;
t439 = qJD(1) * qJD(6);
t371 = t420 * pkin(6) + t384;
t346 = -t414 * g(3) - t406 * t371;
t372 = (-mrSges(7,1) * t414 + mrSges(7,2) * t406) * qJD(1);
t374 = t406 * qJDD(1) + t414 * t439;
t381 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t442;
t295 = m(7) * t346 + qJDD(6) * mrSges(7,1) - t374 * mrSges(7,3) + qJD(6) * t381 - t372 * t444;
t348 = -t406 * g(3) + t414 * t371;
t376 = t414 * qJDD(1) - t406 * t439;
t379 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t444;
t296 = m(7) * t348 - qJDD(6) * mrSges(7,2) + t376 * mrSges(7,3) - qJD(6) * t379 + t372 * t442;
t436 = -t406 * t295 + t414 * t296;
t250 = (-t315 * t404 - t317 * t403 + (-t358 * t403 + t360 * t404) * t400) * pkin(3) + t337;
t223 = m(9) * t250 - t248 * mrSges(9,1) + t249 * mrSges(9,2) - t303 * t301 + t304 * t302;
t354 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t406 + Ifges(7,2) * t414) * qJD(1);
t356 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t406 + Ifges(7,4) * t414) * qJD(1);
t435 = t406 * t354 - t414 * t356;
t355 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t410 + Ifges(3,2) * t418) * qJD(1);
t357 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t410 + Ifges(3,4) * t418) * qJD(1);
t434 = t410 * t355 - t418 * t357;
t369 = qJDD(1) * pkin(6) - t383;
t433 = -m(7) * t369 + t376 * mrSges(7,1) - t374 * mrSges(7,2) + t381 * t442;
t432 = m(5) * t279 - t260 * mrSges(5,1) + t261 * mrSges(5,2) - t331 * t319 + t332 * t320 + t196;
t284 = Ifges(5,5) * t332 + Ifges(5,6) * t331 + Ifges(5,3) * t395;
t182 = -pkin(10) * t196 + mrSges(5,2) * t279 - mrSges(5,3) * t240 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * t393 - t407 * t200 + t415 * t201 + t331 * t284 - t395 * t285;
t426 = mrSges(6,1) * t221 - mrSges(6,2) * t222 + Ifges(6,5) * t238 + Ifges(6,6) * t237 + Ifges(6,3) * t258 + t308 * t252 - t307 * t253;
t183 = -pkin(8) * t196 - mrSges(5,1) * t279 + mrSges(5,3) * t241 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * t393 - t332 * t284 + t395 * t286 - t426;
t322 = Ifges(4,5) * t361 + Ifges(4,6) * t359 + Ifges(4,3) * t401;
t177 = -pkin(4) * t432 - mrSges(4,1) * t337 + mrSges(4,3) * t300 + Ifges(4,4) * t318 + Ifges(4,2) * t316 + Ifges(4,6) * t399 + t408 * t182 + t416 * t183 - t361 * t322 + t401 * t326;
t178 = mrSges(4,2) * t337 - mrSges(4,3) * t298 + Ifges(4,1) * t318 + Ifges(4,4) * t316 + Ifges(4,5) * t399 + t416 * t182 - t408 * t183 + t359 * t322 - t401 * t324;
t264 = Ifges(9,5) * t304 + Ifges(9,6) * t303 + Ifges(9,3) * t394;
t217 = -mrSges(9,1) * t250 + mrSges(9,3) * t235 + Ifges(9,4) * t249 + Ifges(9,2) * t248 + Ifges(9,6) * t392 - t304 * t264 + t394 * t266;
t218 = mrSges(9,2) * t250 - mrSges(9,3) * t234 + Ifges(9,1) * t249 + Ifges(9,4) * t248 + Ifges(9,5) * t392 + t303 * t264 - t394 * t265;
t321 = Ifges(8,5) * t360 + Ifges(8,6) * t358 + Ifges(8,3) * t400;
t190 = -mrSges(8,1) * t337 + mrSges(8,3) * t299 + Ifges(8,4) * t317 + Ifges(8,2) * t315 + Ifges(8,6) * t398 + t367 * t217 + t366 * t218 - t223 * t449 - t360 * t321 + t400 * t325;
t191 = mrSges(8,2) * t337 - mrSges(8,3) * t297 + Ifges(8,1) * t317 + Ifges(8,4) * t315 + Ifges(8,5) * t398 - t366 * t217 + t367 * t218 - t223 * t450 + t358 * t321 - t400 * t323;
t353 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t410 + Ifges(3,6) * t418) * qJD(1);
t423 = t315 * mrSges(8,1) + (-m(4) - m(8)) * t337 + t316 * mrSges(4,1) - t432 - t361 * t345 - t360 * t344 - t317 * mrSges(8,2) - t318 * mrSges(4,2) + t359 * t343 + t358 * t342 - t223;
t173 = pkin(1) * t423 - mrSges(3,1) * t368 + mrSges(3,3) * t349 + Ifges(3,4) * t375 + Ifges(3,2) * t377 + Ifges(3,6) * qJDD(2) + qJD(2) * t357 - t417 * t177 - t409 * t178 + t413 * t190 + t405 * t191 - t353 * t443;
t175 = mrSges(3,2) * t368 - mrSges(3,3) * t347 + Ifges(3,1) * t375 + Ifges(3,4) * t377 + Ifges(3,5) * qJDD(2) - qJD(2) * t355 + t409 * t177 - t417 * t178 - t405 * t190 + t413 * t191 + t353 * t441;
t352 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t406 + Ifges(7,6) * t414) * qJD(1);
t275 = -mrSges(7,1) * t369 + mrSges(7,3) * t348 + Ifges(7,4) * t374 + Ifges(7,2) * t376 + Ifges(7,6) * qJDD(6) + qJD(6) * t356 - t352 * t444;
t276 = mrSges(7,2) * t369 - mrSges(7,3) * t346 + Ifges(7,1) * t374 + Ifges(7,4) * t376 + Ifges(7,5) * qJDD(6) - qJD(6) * t354 + t352 * t442;
t292 = -t379 * t444 + t433;
t380 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t443;
t382 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t441;
t422 = -m(3) * t368 + t377 * mrSges(3,1) - t375 * mrSges(3,2) + t382 * t441 + t423;
t429 = -pkin(6) * t292 - mrSges(2,2) * t384 + t418 * t173 + t410 * t175 + pkin(12) * (-t380 * t443 + t422) + t406 * t276 + t414 * t275 + mrSges(2,1) * t383 + Ifges(2,3) * qJDD(1);
t428 = mrSges(7,1) * t346 - mrSges(7,2) * t348 + Ifges(7,5) * t374 + Ifges(7,6) * t376 + Ifges(7,3) * qJDD(6);
t373 = (-mrSges(3,1) * t418 + mrSges(3,2) * t410) * qJD(1);
t188 = t422 + qJDD(1) * mrSges(2,1) - t420 * mrSges(2,2) + m(2) * t383 + (-t406 * t379 - t410 * t380) * qJD(1) + t433;
t180 = m(3) * t349 - qJDD(2) * mrSges(3,2) + t377 * mrSges(3,3) - qJD(2) * t380 + t409 * t184 - t417 * t185 - t405 * t204 + t413 * t205 + t373 * t441;
t179 = m(3) * t347 + qJDD(2) * mrSges(3,1) - t375 * mrSges(3,3) + qJD(2) * t382 - t373 * t443 + t448;
t176 = m(2) * t384 - t420 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t410 * t179 + t418 * t180 + t436;
t171 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - pkin(12) * (t418 * t179 + t410 * t180) + t420 * Ifges(2,5) + pkin(6) * (t414 * t295 + t406 * t296) + pkin(13) * t436 + mrSges(2,3) * t384 + (-t434 - t435) * qJD(1) - t428 + t452;
t170 = -mrSges(2,2) * g(3) - mrSges(2,3) * t383 - pkin(13) * t292 + Ifges(2,5) * qJDD(1) - t420 * Ifges(2,6) - t410 * t173 + t418 * t175 - t406 * t275 + t414 * t276;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t419 * t170 - t411 * t171 - pkin(11) * (t411 * t176 + t419 * t188), t170, t175, t178, t182, t201, t276, t191, t218, 0, 0, 0, 0; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t411 * t170 + t419 * t171 + pkin(11) * (t419 * t176 - t411 * t188), t171, t173, t177, t183, t200, t275, t190, t217, 0, 0, 0, 0; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t429, t429, t434 * qJD(1) - t452, -t424, -t427, t426, t435 * qJD(1) + t428, -t425, -t430, 0, 0, 0, 0;];
m_new = t1;
