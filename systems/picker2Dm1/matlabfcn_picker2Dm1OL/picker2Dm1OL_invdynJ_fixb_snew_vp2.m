% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [12x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = picker2Dm1OL_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_invdynJ_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:45:07
% EndTime: 2020-05-11 05:45:08
% DurationCPUTime: 0.47s
% Computational Cost: add. (3484->129), mult. (4456->163), div. (0->0), fcn. (2496->22), ass. (0->93)
t365 = sin(qJ(1));
t374 = cos(qJ(1));
t385 = -t365 * g(1) + t374 * g(2);
t320 = qJDD(1) * pkin(1) + t385;
t386 = t374 * g(1) + t365 * g(2);
t321 = -qJD(1) ^ 2 * pkin(1) + t386;
t364 = sin(qJ(2));
t373 = cos(qJ(2));
t312 = t373 * t320 - t364 * t321;
t350 = qJDD(1) + qJDD(2);
t305 = t350 * pkin(3) + t312;
t314 = t364 * t320 + t373 * t321;
t352 = qJD(1) + qJD(2);
t348 = t352 ^ 2;
t307 = -t348 * pkin(3) + t314;
t362 = sin(qJ(4));
t371 = cos(qJ(4));
t294 = t371 * t305 - t362 * t307;
t342 = qJD(4) + t352;
t336 = t342 ^ 2;
t339 = qJDD(4) + t350;
t286 = t339 * pkin(4) + t294;
t296 = t362 * t305 + t371 * t307;
t288 = -t336 * pkin(4) + t296;
t353 = sin(qJ(10));
t355 = cos(qJ(10));
t282 = -t355 * t286 + t353 * t288;
t324 = qJDD(10) + t339;
t331 = qJD(10) + t342;
t328 = t331 ^ 2;
t276 = m(11) * t282 + t324 * mrSges(11,1) - t328 * mrSges(11,2);
t283 = -t353 * t286 - t355 * t288;
t277 = m(11) * t283 - t328 * mrSges(11,1) - t324 * mrSges(11,2);
t379 = -t355 * t276 - t353 * t277;
t270 = m(5) * t294 + t339 * mrSges(5,1) - t336 * mrSges(5,2) + t379;
t271 = m(5) * t296 - t336 * mrSges(5,1) - t339 * mrSges(5,2) + t353 * t276 - t355 * t277;
t387 = t371 * t270 + t362 * t271;
t343 = qJD(3) + t352;
t306 = t350 * pkin(2) + t312;
t308 = -t348 * pkin(2) + t314;
t363 = sin(qJ(3));
t372 = cos(qJ(3));
t295 = -t372 * t306 + t363 * t308;
t340 = qJDD(3) + t350;
t360 = sin(qJ(6));
t369 = cos(qJ(6));
t299 = -t369 * t312 + t360 * t314;
t300 = -t360 * t312 - t369 * t314;
t338 = qJDD(6) + t350;
t384 = mrSges(7,1) * t299 - mrSges(7,2) * t300 + Ifges(7,3) * t338;
t358 = sin(qJ(8));
t367 = cos(qJ(8));
t311 = -t367 * t320 + t358 * t321;
t313 = -t358 * t320 - t367 * t321;
t349 = qJDD(1) + qJDD(8);
t383 = mrSges(9,1) * t311 - mrSges(9,2) * t313 + Ifges(9,3) * t349;
t287 = t340 * pkin(6) + t295;
t297 = -t363 * t306 - t372 * t308;
t337 = t343 ^ 2;
t289 = -t337 * pkin(6) + t297;
t357 = sin(qJ(9));
t366 = cos(qJ(9));
t284 = -t366 * t287 + t357 * t289;
t285 = -t357 * t287 - t366 * t289;
t329 = qJDD(9) + t340;
t382 = mrSges(10,1) * t284 - mrSges(10,2) * t285 + Ifges(10,3) * t329;
t381 = mrSges(11,1) * t282 - mrSges(11,2) * t283 + Ifges(11,3) * t324;
t332 = qJD(9) + t343;
t330 = t332 ^ 2;
t278 = m(10) * t284 + t329 * mrSges(10,1) - t330 * mrSges(10,2);
t279 = m(10) * t285 - t330 * mrSges(10,1) - t329 * mrSges(10,2);
t378 = -t366 * t278 - t357 * t279;
t272 = m(4) * t295 + t340 * mrSges(4,1) - t337 * mrSges(4,2) + t378;
t273 = m(4) * t297 - t337 * mrSges(4,1) - t340 * mrSges(4,2) + t357 * t278 - t366 * t279;
t380 = -t372 * t272 - t363 * t273;
t377 = pkin(6) * t378 + mrSges(4,1) * t295 - mrSges(4,2) * t297 + Ifges(4,3) * t340 + t382;
t376 = pkin(4) * t379 + mrSges(5,1) * t294 - mrSges(5,2) * t296 + Ifges(5,3) * t339 + t381;
t375 = pkin(2) * t380 + pkin(3) * t387 + mrSges(3,1) * t312 - mrSges(3,2) * t314 + Ifges(3,3) * t350 + t376 + t377 + t384;
t370 = cos(qJ(5));
t368 = cos(qJ(7));
t361 = sin(qJ(5));
t359 = sin(qJ(7));
t356 = cos(pkin(8));
t354 = sin(pkin(8));
t351 = qJD(1) + qJD(8);
t347 = t351 ^ 2;
t341 = qJD(6) + t352;
t335 = t341 ^ 2;
t319 = -t354 * t361 + t356 * t370;
t318 = t354 * t370 + t356 * t361;
t291 = m(7) * t300 - t335 * mrSges(7,1) - t338 * mrSges(7,2);
t290 = m(7) * t299 + t338 * mrSges(7,1) - t335 * mrSges(7,2);
t1 = [mrSges(2,1) * t385 + Ifges(2,3) * qJDD(1) - mrSges(2,2) * t386 + (t364 * (m(3) * t314 - t348 * mrSges(3,1) - t350 * mrSges(3,2) - t362 * t270 + t371 * t271 + t363 * t272 - t372 * t273 + t360 * t290 - t369 * t291) + t373 * (m(3) * t312 + t350 * mrSges(3,1) - t348 * mrSges(3,2) - t369 * t290 - t360 * t291 + t380 + t387) - t358 * (m(9) * t313 - t347 * mrSges(9,1) - t349 * mrSges(9,2)) - t367 * (m(9) * t311 + t349 * mrSges(9,1) - t347 * mrSges(9,2))) * pkin(1) + t375 + t383; t375; t377; t376; Ifges(6,3) * qJDD(5) + mrSges(6,1) * (t318 * g(1) - t319 * g(2)) - mrSges(6,2) * (-t319 * g(1) - t318 * g(2)); t384; Ifges(8,3) * qJDD(7) + mrSges(8,1) * (-t368 * g(1) - t359 * g(2)) - mrSges(8,2) * (-t359 * g(1) + t368 * g(2)); t383; t382; t381; 0; 0;];
tauJ = t1;
