% Calculate vector of inverse dynamics base forces with Newton-Euler for
% picker2Dm2OL
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = picker2Dm2OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2OL_invdynB_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:19:09
% EndTime: 2020-05-09 23:19:10
% DurationCPUTime: 1.11s
% Computational Cost: add. (10518->282), mult. (13224->341), div. (0->0), fcn. (7222->22), ass. (0->131)
t536 = pkin(5) * m(6);
t535 = -m(4) - m(10);
t534 = -m(11) - m(5);
t512 = sin(qJ(1));
t521 = cos(qJ(1));
t473 = -t512 * g(1) + t521 * g(2);
t470 = qJDD(1) * pkin(1) + t473;
t475 = t521 * g(1) + t512 * g(2);
t524 = qJD(1) ^ 2;
t471 = -t524 * pkin(1) + t475;
t511 = sin(qJ(2));
t520 = cos(qJ(2));
t455 = t520 * t470 - t511 * t471;
t497 = qJDD(1) + qJDD(2);
t450 = t497 * pkin(3) + t455;
t457 = t511 * t470 + t520 * t471;
t499 = qJD(1) + qJD(2);
t495 = t499 ^ 2;
t452 = -t495 * pkin(3) + t457;
t509 = sin(qJ(4));
t518 = cos(qJ(4));
t428 = t518 * t450 - t509 * t452;
t489 = qJD(4) + t499;
t483 = t489 ^ 2;
t486 = qJDD(4) + t497;
t422 = t486 * pkin(4) + t428;
t430 = t509 * t450 + t518 * t452;
t424 = -t483 * pkin(4) + t430;
t500 = sin(qJ(10));
t502 = cos(qJ(10));
t416 = -t502 * t422 + t500 * t424;
t476 = qJDD(10) + t486;
t480 = qJD(10) + t489;
t477 = t480 ^ 2;
t412 = m(11) * t416 + t476 * mrSges(11,1) - t477 * mrSges(11,2);
t417 = -t500 * t422 - t502 * t424;
t413 = m(11) * t417 - t477 * mrSges(11,1) - t476 * mrSges(11,2);
t527 = -t502 * t412 - t500 * t413;
t404 = m(5) * t428 + t486 * mrSges(5,1) - t483 * mrSges(5,2) + t527;
t405 = m(5) * t430 - t483 * mrSges(5,1) - t486 * mrSges(5,2) + t500 * t412 - t502 * t413;
t533 = t518 * t404 + t509 * t405;
t501 = sin(pkin(8));
t503 = cos(pkin(8));
t508 = sin(qJ(5));
t517 = cos(qJ(5));
t468 = t501 * t517 + t503 * t508;
t469 = -t501 * t508 + t503 * t517;
t458 = t468 * g(1) - t469 * g(2);
t523 = qJD(5) ^ 2;
t448 = m(6) * t458 + qJDD(5) * mrSges(6,1) - t523 * mrSges(6,2);
t459 = -t469 * g(1) - t468 * g(2);
t449 = m(6) * t459 - t523 * mrSges(6,1) - qJDD(5) * mrSges(6,2);
t532 = t469 * t448 + t468 * t449;
t531 = -t468 * t448 + t469 * t449;
t490 = qJD(3) + t499;
t451 = t497 * pkin(2) + t455;
t453 = -t495 * pkin(2) + t457;
t510 = sin(qJ(3));
t519 = cos(qJ(3));
t429 = -t519 * t451 + t510 * t453;
t506 = sin(qJ(7));
t515 = cos(qJ(7));
t472 = -t506 * g(1) + t515 * g(2);
t522 = qJD(7) ^ 2;
t463 = m(8) * t472 - t522 * mrSges(8,1) - qJDD(7) * mrSges(8,2);
t474 = -t515 * g(1) - t506 * g(2);
t464 = m(8) * t474 + qJDD(7) * mrSges(8,1) - t522 * mrSges(8,2);
t530 = -t515 * t463 + t506 * t464;
t487 = qJDD(3) + t497;
t529 = m(3) + m(7) + m(9) - t534 - t535;
t484 = t490 ^ 2;
t423 = t487 * pkin(6) + t429;
t431 = -t510 * t451 - t519 * t453;
t425 = -t484 * pkin(6) + t431;
t504 = sin(qJ(9));
t513 = cos(qJ(9));
t418 = -t513 * t423 + t504 * t425;
t478 = qJDD(9) + t487;
t481 = qJD(9) + t490;
t479 = t481 ^ 2;
t414 = m(10) * t418 + t478 * mrSges(10,1) - t479 * mrSges(10,2);
t419 = -t504 * t423 - t513 * t425;
t415 = m(10) * t419 - t479 * mrSges(10,1) - t478 * mrSges(10,2);
t526 = -t513 * t414 - t504 * t415;
t406 = m(4) * t429 + t487 * mrSges(4,1) - t484 * mrSges(4,2) + t526;
t407 = m(4) * t431 - t484 * mrSges(4,1) - t487 * mrSges(4,2) + t504 * t414 - t513 * t415;
t528 = -t519 * t406 - t510 * t407;
t507 = sin(qJ(6));
t516 = cos(qJ(6));
t432 = -t516 * t455 + t507 * t457;
t488 = qJD(6) + t499;
t482 = t488 ^ 2;
t485 = qJDD(6) + t497;
t426 = m(7) * t432 + t485 * mrSges(7,1) - t482 * mrSges(7,2);
t433 = -t507 * t455 - t516 * t457;
t427 = m(7) * t433 - t482 * mrSges(7,1) - t485 * mrSges(7,2);
t396 = m(3) * t455 + t497 * mrSges(3,1) - t495 * mrSges(3,2) - t516 * t426 - t507 * t427 + t528 + t533;
t397 = m(3) * t457 - t495 * mrSges(3,1) - t497 * mrSges(3,2) - t509 * t404 + t518 * t405 + t510 * t406 - t519 * t407 + t507 * t426 - t516 * t427;
t505 = sin(qJ(8));
t514 = cos(qJ(8));
t454 = -t514 * t470 + t505 * t471;
t498 = qJD(1) + qJD(8);
t494 = t498 ^ 2;
t496 = qJDD(1) + qJDD(8);
t440 = m(9) * t454 + t496 * mrSges(9,1) - t494 * mrSges(9,2);
t456 = -t505 * t470 - t514 * t471;
t441 = m(9) * t456 - t494 * mrSges(9,1) - t496 * mrSges(9,2);
t525 = t520 * t396 + t511 * t397 - t514 * t440 - t505 * t441;
t462 = -mrSges(8,2) * g(3) - mrSges(8,3) * t474 + Ifges(8,5) * qJDD(7) - t522 * Ifges(8,6);
t461 = mrSges(8,1) * g(3) + mrSges(8,3) * t472 + t522 * Ifges(8,5) + Ifges(8,6) * qJDD(7);
t443 = -mrSges(6,2) * g(3) - mrSges(6,3) * t458 + Ifges(6,5) * qJDD(5) - t523 * Ifges(6,6);
t442 = mrSges(6,1) * g(3) + mrSges(6,3) * t459 + t523 * Ifges(6,5) + Ifges(6,6) * qJDD(5);
t439 = -mrSges(9,2) * g(3) - mrSges(9,3) * t454 + Ifges(9,5) * t496 - t494 * Ifges(9,6);
t438 = mrSges(9,1) * g(3) + mrSges(9,3) * t456 + t494 * Ifges(9,5) + Ifges(9,6) * t496;
t421 = -mrSges(7,2) * g(3) - mrSges(7,3) * t432 + Ifges(7,5) * t485 - t482 * Ifges(7,6);
t420 = mrSges(7,1) * g(3) + mrSges(7,3) * t433 + t482 * Ifges(7,5) + Ifges(7,6) * t485;
t411 = -mrSges(10,2) * g(3) - mrSges(10,3) * t418 + Ifges(10,5) * t478 - t479 * Ifges(10,6);
t410 = mrSges(10,1) * g(3) + mrSges(10,3) * t419 + t479 * Ifges(10,5) + Ifges(10,6) * t478;
t409 = -mrSges(11,2) * g(3) - mrSges(11,3) * t416 + Ifges(11,5) * t476 - t477 * Ifges(11,6);
t408 = mrSges(11,1) * g(3) + mrSges(11,3) * t417 + t477 * Ifges(11,5) + Ifges(11,6) * t476;
t401 = -mrSges(4,2) * g(3) - mrSges(4,3) * t429 + Ifges(4,5) * t487 - t484 * Ifges(4,6) + t504 * t410 - t513 * t411;
t400 = -mrSges(5,2) * g(3) - mrSges(5,3) * t428 + Ifges(5,5) * t486 - t483 * Ifges(5,6) + t500 * t408 - t502 * t409;
t399 = mrSges(4,3) * t431 + t484 * Ifges(4,5) + Ifges(4,6) * t487 - t513 * t410 - t504 * t411 + (pkin(6) * m(10) + mrSges(4,1)) * g(3);
t398 = mrSges(5,3) * t430 + t483 * Ifges(5,5) + Ifges(5,6) * t486 - t502 * t408 - t500 * t409 + (pkin(4) * m(11) + mrSges(5,1)) * g(3);
t393 = -mrSges(3,2) * g(3) - mrSges(3,3) * t455 + Ifges(3,5) * t497 - t495 * Ifges(3,6) - t509 * t398 + t510 * t399 + t518 * t400 - t519 * t401 + t507 * t420 - t516 * t421;
t392 = mrSges(3,3) * t457 + t495 * Ifges(3,5) + Ifges(3,6) * t497 + t518 * t398 - t519 * t399 + t509 * t400 - t510 * t401 - t516 * t420 - t507 * t421 + (-pkin(2) * t535 - pkin(3) * t534 + mrSges(3,1)) * g(3);
t391 = m(2) * t475 - t524 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t511 * t396 + t520 * t397 + t505 * t440 - t514 * t441;
t390 = m(2) * t473 + qJDD(1) * mrSges(2,1) - t524 * mrSges(2,2) + t525;
t389 = -mrSges(2,2) * g(3) - mrSges(2,3) * t473 + Ifges(2,5) * qJDD(1) - t524 * Ifges(2,6) - t511 * t392 + t520 * t393 + t505 * t438 - t514 * t439;
t388 = mrSges(2,3) * t475 + t524 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t520 * t392 + t511 * t393 - t514 * t438 - t505 * t439 + (t529 * pkin(1) + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t512 * t390 - t521 * t391 + t506 * t463 + t515 * t464 + t531; -m(1) * g(2) - t521 * t390 - t512 * t391 + t530 + t532; (-m(1) - m(2) - m(6) - m(8) - t529) * g(3); mrSges(1,3) * g(2) + t512 * t388 - t521 * t389 - t468 * t442 + t469 * t443 + t515 * t461 + t506 * t462 + (-t501 * t536 - mrSges(1,2)) * g(3); -mrSges(1,3) * g(1) - t521 * t388 - t512 * t389 + t469 * t442 + t468 * t443 + t506 * t461 - t515 * t462 + (m(8) * pkin(7) + t503 * t536 + mrSges(1,1)) * g(3); pkin(7) * t530 + mrSges(11,1) * t416 - mrSges(11,2) * t417 + mrSges(10,1) * t418 - mrSges(10,2) * t419 + mrSges(5,1) * t428 + mrSges(4,1) * t429 - mrSges(5,2) * t430 - mrSges(4,2) * t431 + mrSges(7,1) * t432 - mrSges(7,2) * t433 + mrSges(9,1) * t454 + mrSges(3,1) * t455 - mrSges(9,2) * t456 - mrSges(3,2) * t457 + mrSges(6,1) * t458 - mrSges(6,2) * t459 - mrSges(8,2) * t472 + mrSges(2,1) * t473 + mrSges(8,1) * t474 - mrSges(2,2) * t475 + Ifges(11,3) * t476 + Ifges(2,3) * qJDD(1) + Ifges(8,3) * qJDD(7) + Ifges(6,3) * qJDD(5) + Ifges(10,3) * t478 + Ifges(7,3) * t485 + Ifges(5,3) * t486 + Ifges(4,3) * t487 + Ifges(9,3) * t496 + Ifges(3,3) * t497 - mrSges(1,1) * g(2) + pkin(4) * t527 + pkin(2) * t528 + mrSges(1,2) * g(1) + (-t501 * t531 + t503 * t532) * pkin(5) + pkin(3) * t533 + t525 * pkin(1) + pkin(6) * t526;];
tauB = t1;
