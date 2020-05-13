% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% fivebar1OL
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
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = fivebar1OL_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:12
% EndTime: 2020-04-27 06:13:13
% DurationCPUTime: 0.58s
% Computational Cost: add. (418->128), mult. (636->158), div. (0->0), fcn. (404->8), ass. (0->69)
t363 = sin(qJ(1));
t372 = pkin(2) * m(3);
t362 = sin(qJ(2));
t359 = t362 * mrSges(3,2);
t366 = cos(qJ(2));
t403 = t366 * mrSges(3,1);
t407 = -t403 + t359;
t421 = -mrSges(2,1) - t372 - t407;
t423 = t363 * t421;
t360 = sin(qJ(4));
t358 = t360 * mrSges(5,2);
t364 = cos(qJ(4));
t404 = t364 * mrSges(5,1);
t345 = -t358 + t404;
t331 = pkin(3) * m(5) + mrSges(4,1) + t345;
t344 = t360 * mrSges(5,1) + t364 * mrSges(5,2);
t339 = mrSges(4,2) + t344;
t361 = sin(qJ(3));
t365 = cos(qJ(3));
t390 = -t331 * t365 + t339 * t361;
t367 = cos(qJ(1));
t399 = t421 * t367;
t422 = -t390 + mrSges(1,1) - t399;
t374 = qJD(3) ^ 2;
t380 = t361 * t344 - t345 * t365;
t389 = t361 * t331 + t339 * t365;
t420 = -t390 * qJDD(3) - t380 * qJDD(4) - t389 * t374;
t347 = t362 * mrSges(3,1) + t366 * mrSges(3,2);
t413 = -mrSges(2,2) + t347;
t419 = t363 * t413;
t417 = (2 * qJD(3) + qJD(4)) * qJD(4);
t409 = (2 * qJD(1) + qJD(2)) * qJD(2);
t328 = t413 * t367;
t416 = -mrSges(1,2) + t328 - t389;
t408 = t399 - t419;
t370 = m(4) + m(5);
t355 = pkin(3) * t404;
t401 = Ifges(5,3) + t355;
t393 = pkin(3) * t358;
t343 = t364 * Ifges(5,5) - t360 * Ifges(5,6);
t312 = t328 + t423;
t385 = (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t359 - 0.2e1 * t403 + t372) * pkin(2)) * qJDD(1) + (t407 * pkin(2) + Ifges(3,3)) * qJDD(2) + t409 * pkin(2) * t347;
t346 = t366 * Ifges(3,5) - t362 * Ifges(3,6);
t349 = Ifges(3,5) * t362 + Ifges(3,6) * t366;
t342 = t360 * Ifges(5,5) + t364 * Ifges(5,6);
t329 = -pkin(3) * mrSges(5,3) + Ifges(4,5) + t343;
t338 = Ifges(4,6) + t342;
t384 = t361 * t329 + t338 * t365;
t330 = -pkin(2) * mrSges(3,3) + Ifges(2,5) - t346;
t341 = -Ifges(2,6) + t349;
t383 = t363 * t330 - t341 * t367;
t382 = t342 * t365 + t361 * t343;
t314 = -t361 * t342 + t343 * t365;
t381 = -t365 * t344 - t361 * t345;
t315 = -t346 * t367 + t363 * t349;
t379 = t363 * t346 + t349 * t367;
t316 = t347 * t367 - t363 * t407;
t313 = t363 * t347 + t367 * t407;
t378 = (pkin(3) ^ 2 * m(5) + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t355 - 0.2e1 * t393) * qJDD(3);
t376 = qJD(1) ^ 2;
t356 = m(2) + m(3) + m(1) + t370;
t352 = mrSges(2,3) + mrSges(4,3) + mrSges(5,3) + mrSges(3,3) + mrSges(1,3);
t337 = -t376 * pkin(2) - t367 * g(1) - t363 * g(2);
t336 = -t374 * pkin(3) - t365 * g(1) - t361 * g(2);
t335 = qJDD(1) * pkin(2) + t363 * g(1) - t367 * g(2);
t334 = qJDD(3) * pkin(3) + t361 * g(1) - t365 * g(2);
t309 = t330 * t367 + t363 * t341;
t308 = t329 * t365 - t361 * t338;
t1 = [-t356 * g(1) + t312 * qJDD(1) + t316 * qJDD(2) - qJDD(3) * t389 + t381 * qJDD(4) - t409 * t313 + t390 * t374 + t408 * t376 + t380 * t417; -t356 * g(2) - qJDD(1) * t408 + t313 * qJDD(2) + t312 * t376 + t409 * t316 + t381 * t417 + t420; -t356 * g(3); -t384 * t374 + t308 * qJDD(3) + t314 * qJDD(4) - (-t416 - t423) * g(3) + t315 * qJDD(2) - t383 * t376 + t309 * qJDD(1) + t352 * g(2) - t417 * t382 + t409 * t379; t382 * qJDD(4) + t308 * t374 + t384 * qJDD(3) - (-pkin(1) * t370 - t419 - t422) * g(3) - t352 * g(1) - t379 * qJDD(2) + t309 * t376 + t383 * qJDD(1) + t409 * t315 + t417 * t314; (-t393 + t401) * qJDD(4) + t378 - t416 * g(1) - t422 * g(2) + (-g(1) * t421 - g(2) * t413) * t363 + t417 * (t381 * pkin(1) - pkin(3) * t344) + (-t370 * g(2) + t420) * pkin(1) + t385; -t312 * g(1) + g(2) * t408 + t385; Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (-t366 * t335 + t362 * t337) - mrSges(3,2) * (-t362 * t335 - t366 * t337); t389 * g(1) + t390 * g(2) + t378 + t401 * qJDD(4) + (-qJDD(4) * t358 - t344 * t417) * pkin(3); Ifges(5,3) * (qJDD(3) + qJDD(4)) + mrSges(5,1) * (t364 * t334 - t360 * t336) - mrSges(5,2) * (t360 * t334 + t364 * t336); 0;];
tauJB = t1;
