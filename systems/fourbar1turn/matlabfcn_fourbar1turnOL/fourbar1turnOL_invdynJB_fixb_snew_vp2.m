% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = fourbar1turnOL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:53
% EndTime: 2020-04-12 19:40:54
% DurationCPUTime: 0.60s
% Computational Cost: add. (2626->210), mult. (5503->282), div. (0->0), fcn. (3364->8), ass. (0->86)
t391 = sin(qJ(1));
t395 = cos(qJ(1));
t382 = -t395 * g(1) - t391 * g(2);
t390 = sin(qJ(2));
t394 = cos(qJ(2));
t355 = -t394 * g(3) - t390 * t382;
t396 = qJD(1) ^ 2;
t349 = (t390 * t394 * t396 + qJDD(2)) * pkin(2) + t355;
t356 = -t390 * g(3) + t394 * t382;
t350 = (-t394 ^ 2 * t396 - qJD(2) ^ 2) * pkin(2) + t356;
t389 = sin(qJ(3));
t393 = cos(qJ(3));
t337 = -t393 * t349 + t389 * t350;
t363 = (t389 * t390 - t393 * t394) * qJD(1);
t408 = qJD(1) * qJD(2);
t373 = t390 * qJDD(1) + t394 * t408;
t406 = t390 * t408;
t375 = t394 * qJDD(1) - t406;
t343 = t363 * qJD(3) - t393 * t373 - t389 * t375;
t364 = (-t389 * t394 - t390 * t393) * qJD(1);
t347 = -t363 * mrSges(4,1) + t364 * mrSges(4,2);
t386 = qJD(2) + qJD(3);
t351 = -t386 * mrSges(4,2) + t363 * mrSges(4,3);
t385 = qJDD(2) + qJDD(3);
t325 = m(4) * t337 + t385 * mrSges(4,1) - t343 * mrSges(4,3) - t364 * t347 + t386 * t351;
t338 = -t389 * t349 - t393 * t350;
t342 = -t364 * qJD(3) + t389 * t373 - t393 * t375;
t352 = t386 * mrSges(4,1) - t364 * mrSges(4,3);
t326 = m(4) * t338 - t385 * mrSges(4,2) + t342 * mrSges(4,3) + t363 * t347 - t386 * t352;
t320 = -t393 * t325 - t389 * t326;
t370 = (-mrSges(3,1) * t394 + mrSges(3,2) * t390) * qJD(1);
t409 = qJD(1) * t394;
t380 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t409;
t411 = qJD(1) * t390;
t318 = m(3) * t355 + qJDD(2) * mrSges(3,1) - t373 * mrSges(3,3) + qJD(2) * t380 - t370 * t411 + t320;
t378 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t411;
t319 = m(3) * t356 - qJDD(2) * mrSges(3,2) + t375 * mrSges(3,3) - qJD(2) * t378 + t389 * t325 - t393 * t326 + t370 * t409;
t371 = -t396 * pkin(1) + t382;
t388 = sin(qJ(4));
t392 = cos(qJ(4));
t353 = -t392 * g(3) - t388 * t371;
t369 = (-mrSges(5,1) * t392 + mrSges(5,2) * t388) * qJD(1);
t407 = qJD(1) * qJD(4);
t372 = t388 * qJDD(1) + t392 * t407;
t410 = qJD(1) * t392;
t379 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t410;
t412 = qJD(1) * t388;
t334 = m(5) * t353 + qJDD(4) * mrSges(5,1) - t372 * mrSges(5,3) + qJD(4) * t379 - t369 * t412;
t354 = -t388 * g(3) + t392 * t371;
t374 = t392 * qJDD(1) - t388 * t407;
t377 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t412;
t335 = m(5) * t354 - qJDD(4) * mrSges(5,2) + t374 * mrSges(5,3) - qJD(4) * t377 + t369 * t410;
t313 = m(2) * t382 - mrSges(2,1) * t396 - qJDD(1) * mrSges(2,2) - t318 * t390 + t319 * t394 - t334 * t388 + t335 * t392;
t381 = t391 * g(1) - t395 * g(2);
t348 = (-t375 + t406) * pkin(2) - t381;
t399 = m(4) * t348 - t342 * mrSges(4,1) + t343 * mrSges(4,2) - t363 * t351 + t364 * t352;
t368 = -qJDD(1) * pkin(1) - t381;
t402 = -m(5) * t368 + t374 * mrSges(5,1) - t372 * mrSges(5,2) + t379 * t410;
t322 = qJDD(1) * mrSges(2,1) + t375 * mrSges(3,1) - t396 * mrSges(2,2) - t373 * mrSges(3,2) + (m(2) + m(3)) * t381 + (-t388 * t377 - t390 * t378 + t394 * t380) * qJD(1) - t399 + t402;
t414 = t391 * t313 + t395 * t322;
t413 = t392 * t334 + t388 * t335;
t405 = t395 * t313 - t322 * t391;
t359 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t388 + Ifges(5,2) * t392) * qJD(1);
t361 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t388 + Ifges(5,4) * t392) * qJD(1);
t404 = t359 * t388 - t361 * t392;
t360 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t390 + Ifges(3,2) * t394) * qJD(1);
t362 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t390 + Ifges(3,4) * t394) * qJD(1);
t403 = t360 * t390 - t362 * t394;
t344 = Ifges(4,5) * t364 + Ifges(4,6) * t363 + Ifges(4,3) * t386;
t346 = Ifges(4,1) * t364 + Ifges(4,4) * t363 + Ifges(4,5) * t386;
t323 = -mrSges(4,1) * t348 + mrSges(4,3) * t338 + Ifges(4,4) * t343 + Ifges(4,2) * t342 + Ifges(4,6) * t385 - t364 * t344 + t386 * t346;
t345 = Ifges(4,4) * t364 + Ifges(4,2) * t363 + Ifges(4,6) * t386;
t324 = mrSges(4,2) * t348 - mrSges(4,3) * t337 + Ifges(4,1) * t343 + Ifges(4,4) * t342 + Ifges(4,5) * t385 + t363 * t344 - t386 * t345;
t358 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t390 + Ifges(3,6) * t394) * qJD(1);
t315 = mrSges(3,1) * t381 + mrSges(3,3) * t356 + Ifges(3,4) * t373 + Ifges(3,2) * t375 + Ifges(3,6) * qJDD(2) - pkin(2) * t399 + qJD(2) * t362 - t393 * t323 - t389 * t324 - t358 * t411;
t317 = -mrSges(3,2) * t381 - mrSges(3,3) * t355 + Ifges(3,1) * t373 + Ifges(3,4) * t375 + Ifges(3,5) * qJDD(2) - qJD(2) * t360 + t389 * t323 - t393 * t324 + t358 * t409;
t357 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t388 + Ifges(5,6) * t392) * qJD(1);
t329 = -mrSges(5,1) * t368 + mrSges(5,3) * t354 + Ifges(5,4) * t372 + Ifges(5,2) * t374 + Ifges(5,6) * qJDD(4) + qJD(4) * t361 - t357 * t412;
t330 = mrSges(5,2) * t368 - mrSges(5,3) * t353 + Ifges(5,1) * t372 + Ifges(5,4) * t374 + Ifges(5,5) * qJDD(4) - qJD(4) * t359 + t357 * t410;
t401 = -mrSges(2,2) * t382 + t394 * t315 + t390 * t317 + t392 * t329 + t388 * t330 + pkin(1) * (-t377 * t412 + t402) + mrSges(2,1) * t381 + Ifges(2,3) * qJDD(1);
t400 = -mrSges(4,1) * t337 + mrSges(4,2) * t338 - Ifges(4,5) * t343 - Ifges(4,6) * t342 - Ifges(4,3) * t385 - t364 * t345 + t363 * t346;
t398 = mrSges(5,1) * t353 - mrSges(5,2) * t354 + Ifges(5,5) * t372 + Ifges(5,6) * t374 + Ifges(5,3) * qJDD(4);
t397 = mrSges(3,1) * t355 - mrSges(3,2) * t356 + Ifges(3,5) * t373 + Ifges(3,6) * t375 + Ifges(3,3) * qJDD(2) + pkin(2) * t320 - t400;
t310 = (-t403 - t404) * qJD(1) - pkin(1) * t413 + mrSges(2,3) * t382 + t396 * Ifges(2,5) - t397 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - t398;
t309 = -mrSges(2,2) * g(3) - mrSges(2,3) * t381 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t396 - t315 * t390 + t317 * t394 - t329 * t388 + t330 * t392;
t1 = [-m(1) * g(1) + t405; -m(1) * g(2) + t414; t394 * t318 + t390 * t319 + (-m(1) - m(2)) * g(3) + t413; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t414 + t395 * t309 - t391 * t310; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t405 + t391 * t309 + t395 * t310; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t401; t401; t403 * qJD(1) + t397; -t400; t404 * qJD(1) + t398; 0;];
tauJB = t1;
