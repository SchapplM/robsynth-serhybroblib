% Calculate vector of cutting torques with Newton-Euler for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% m [3x11]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = palh1m1OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynm_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:29:07
% EndTime: 2020-04-15 19:29:51
% DurationCPUTime: 10.29s
% Computational Cost: add. (128225->619), mult. (281752->790), div. (0->0), fcn. (216629->22), ass. (0->247)
t516 = sin(qJ(1));
t525 = cos(qJ(1));
t479 = -t525 * g(1) - t516 * g(2);
t526 = qJD(1) ^ 2;
t465 = -t526 * pkin(15) + t479;
t515 = sin(qJ(2));
t524 = cos(qJ(2));
t442 = -t524 * g(3) - t515 * t465;
t444 = t515 * g(3) - t524 * t465;
t551 = qJD(1) * qJD(2);
t549 = t524 * t551;
t470 = -t515 * qJDD(1) - t549;
t472 = t524 * qJDD(1) - t515 * t551;
t430 = (-t515 * t524 * t526 + qJDD(2)) * pkin(1) + t444;
t431 = (-t515 ^ 2 * t526 - qJD(2) ^ 2) * pkin(1) + t442;
t514 = sin(qJ(3));
t523 = cos(qJ(3));
t370 = -t523 * t430 + t514 * t431;
t372 = t514 * t430 + t523 * t431;
t457 = (-t514 * t515 + t523 * t524) * qJD(1);
t392 = t457 * qJD(3) - t523 * t470 + t514 * t472;
t454 = (t514 * t524 + t515 * t523) * qJD(1);
t395 = -t454 * qJD(3) + t514 * t470 + t523 * t472;
t503 = qJD(2) + qJD(3);
t408 = Ifges(4,4) * t454 + Ifges(4,2) * t457 + Ifges(4,6) * t503;
t411 = Ifges(4,1) * t454 + Ifges(4,4) * t457 + Ifges(4,5) * t503;
t500 = qJDD(2) + qJDD(3);
t513 = sin(qJ(4));
t522 = cos(qJ(4));
t419 = t522 * t454 + t513 * t457;
t320 = -t419 * qJD(4) - t513 * t392 + t522 * t395;
t417 = -t513 * t454 + t522 * t457;
t322 = t417 * qJD(4) + t522 * t392 + t513 * t395;
t478 = t516 * g(1) - t525 * g(2);
t463 = -qJDD(1) * pkin(15) - t478;
t426 = t463 + (-t470 + t549) * pkin(1);
t341 = t426 + (t454 * t503 - t395) * pkin(5);
t494 = qJD(4) + t503;
t283 = (-t417 * t494 - t322) * pkin(11) + (t419 * t494 - t320) * pkin(9) + t341;
t336 = (t454 * t457 + t500) * pkin(5) + t372;
t345 = (-t457 ^ 2 - t503 ^ 2) * pkin(5) + t370;
t294 = t513 * t336 + t522 * t345;
t359 = -t417 * pkin(9) - t419 * pkin(11);
t489 = t494 ^ 2;
t491 = qJDD(4) + t500;
t285 = -t489 * pkin(9) + t491 * pkin(11) + t417 * t359 + t294;
t512 = sin(qJ(5));
t521 = cos(qJ(5));
t272 = t512 * t283 + t521 * t285;
t293 = t522 * t336 - t513 * t345;
t284 = -t491 * pkin(9) - t489 * pkin(11) + t419 * t359 - t293;
t379 = t521 * t419 + t512 * t494;
t290 = -t379 * qJD(5) - t512 * t322 + t521 * t491;
t378 = -t512 * t419 + t521 * t494;
t291 = t378 * qJD(5) + t521 * t322 + t512 * t491;
t412 = qJD(5) - t417;
t307 = Ifges(6,5) * t379 + Ifges(6,6) * t378 + Ifges(6,3) * t412;
t309 = Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t412;
t317 = qJDD(5) - t320;
t244 = -mrSges(6,1) * t284 + mrSges(6,3) * t272 + Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t317 - t379 * t307 + t412 * t309;
t271 = t521 * t283 - t512 * t285;
t308 = Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t412;
t245 = mrSges(6,2) * t284 - mrSges(6,3) * t271 + Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t317 + t378 * t307 - t412 * t308;
t349 = Ifges(5,4) * t419 + Ifges(5,2) * t417 + Ifges(5,6) * t494;
t351 = Ifges(5,1) * t419 + Ifges(5,4) * t417 + Ifges(5,5) * t494;
t342 = -t412 * mrSges(6,2) + t378 * mrSges(6,3);
t343 = t412 * mrSges(6,1) - t379 * mrSges(6,3);
t540 = -m(6) * t284 + t290 * mrSges(6,1) - t291 * mrSges(6,2) + t378 * t342 - t379 * t343;
t339 = -t378 * mrSges(6,1) + t379 * mrSges(6,2);
t260 = m(6) * t271 + t317 * mrSges(6,1) - t291 * mrSges(6,3) - t379 * t339 + t412 * t342;
t261 = m(6) * t272 - t317 * mrSges(6,2) + t290 * mrSges(6,3) + t378 * t339 - t412 * t343;
t548 = -t512 * t260 + t521 * t261;
t534 = -pkin(9) * t540 - pkin(11) * t548 - mrSges(5,1) * t293 + mrSges(5,2) * t294 - Ifges(5,5) * t322 - Ifges(5,6) * t320 - Ifges(5,3) * t491 - t521 * t244 - t512 * t245 - t419 * t349 + t417 * t351;
t358 = -t417 * mrSges(5,1) + t419 * mrSges(5,2);
t401 = t494 * mrSges(5,1) - t419 * mrSges(5,3);
t237 = m(5) * t294 - t491 * mrSges(5,2) + t320 * mrSges(5,3) + t417 * t358 - t494 * t401 + t548;
t399 = -t494 * mrSges(5,2) + t417 * mrSges(5,3);
t252 = m(5) * t293 + t491 * mrSges(5,1) - t322 * mrSges(5,3) - t419 * t358 + t494 * t399 + t540;
t558 = t513 * t237 + t522 * t252;
t530 = -pkin(5) * t558 - mrSges(4,1) * t372 + mrSges(4,2) * t370 - Ifges(4,5) * t392 - Ifges(4,6) * t395 - Ifges(4,3) * t500 - t454 * t408 + t457 * t411 + t534;
t510 = sin(qJ(7));
t519 = cos(qJ(7));
t369 = t519 * t430 - t510 * t431;
t371 = t510 * t430 + t519 * t431;
t456 = (-t510 * t515 + t519 * t524) * qJD(1);
t391 = -t456 * qJD(7) + t519 * t470 - t510 * t472;
t453 = (-t510 * t524 - t515 * t519) * qJD(1);
t394 = t453 * qJD(7) + t510 * t470 + t519 * t472;
t502 = qJD(2) + qJD(7);
t407 = Ifges(8,4) * t456 + Ifges(8,2) * t453 + Ifges(8,6) * t502;
t410 = Ifges(8,1) * t456 + Ifges(8,4) * t453 + Ifges(8,5) * t502;
t499 = qJDD(2) + qJDD(7);
t505 = sin(pkin(19));
t506 = cos(pkin(19));
t402 = (-t453 * t506 - t456 * t505) * pkin(4);
t497 = t502 ^ 2;
t323 = -t456 * t402 + (t497 * t505 + t499 * t506) * pkin(4) + t369;
t324 = t453 * t402 + (-t497 * t506 + t499 * t505) * pkin(4) + t371;
t507 = cos(qJ(10));
t560 = sin(qJ(10));
t459 = t505 * t507 - t506 * t560;
t460 = -t505 * t560 - t506 * t507;
t287 = t460 * t323 - t459 * t324;
t288 = t459 * t323 + t460 * t324;
t374 = t459 * t453 + t460 * t456;
t300 = -t374 * qJD(10) + t460 * t391 - t459 * t394;
t373 = t460 * t453 - t459 * t456;
t301 = t373 * qJD(10) + t459 * t391 + t460 * t394;
t492 = qJD(10) + t502;
t326 = Ifges(11,4) * t374 + Ifges(11,2) * t373 + Ifges(11,6) * t492;
t327 = Ifges(11,1) * t374 + Ifges(11,4) * t373 + Ifges(11,5) * t492;
t485 = qJDD(10) + t499;
t537 = -mrSges(11,1) * t287 + mrSges(11,2) * t288 - Ifges(11,5) * t301 - Ifges(11,6) * t300 - Ifges(11,3) * t485 - t374 * t326 + t373 * t327;
t332 = -t373 * mrSges(11,1) + t374 * mrSges(11,2);
t365 = -t492 * mrSges(11,2) + t373 * mrSges(11,3);
t277 = m(11) * t287 + t485 * mrSges(11,1) - t301 * mrSges(11,3) - t374 * t332 + t492 * t365;
t366 = t492 * mrSges(11,1) - t374 * mrSges(11,3);
t278 = m(11) * t288 - t485 * mrSges(11,2) + t300 * mrSges(11,3) + t373 * t332 - t492 * t366;
t556 = -t459 * t277 + t460 * t278;
t557 = t460 * t277 + t459 * t278;
t561 = pkin(4) * t506;
t562 = pkin(4) * t505;
t531 = -mrSges(8,1) * t369 + mrSges(8,2) * t371 - Ifges(8,5) * t394 - Ifges(8,6) * t391 - Ifges(8,3) * t499 - t456 * t407 + t453 * t410 - t556 * t562 - t557 * t561 + t537;
t509 = sin(qJ(8));
t518 = cos(qJ(8));
t455 = (-t509 * t515 + t518 * t524) * qJD(1);
t390 = -t455 * qJD(8) + t518 * t470 - t509 * t472;
t452 = (-t509 * t524 - t515 * t518) * qJD(1);
t393 = t452 * qJD(8) + t509 * t470 + t518 * t472;
t396 = -t509 * t442 + t518 * t444;
t397 = t518 * t442 + t509 * t444;
t501 = qJD(2) + qJD(8);
t406 = Ifges(9,4) * t455 + Ifges(9,2) * t452 + Ifges(9,6) * t501;
t409 = Ifges(9,1) * t455 + Ifges(9,4) * t452 + Ifges(9,5) * t501;
t498 = qJDD(2) + qJDD(8);
t352 = (t452 * t455 + t498) * pkin(2) + t396;
t362 = (-t452 ^ 2 - t501 ^ 2) * pkin(2) + t397;
t508 = sin(qJ(9));
t517 = cos(qJ(9));
t305 = -t517 * t352 + t508 * t362;
t306 = -t508 * t352 - t517 * t362;
t418 = -t508 * t452 - t517 * t455;
t319 = -t418 * qJD(9) - t517 * t390 + t508 * t393;
t416 = -t517 * t452 + t508 * t455;
t321 = t416 * qJD(9) - t508 * t390 - t517 * t393;
t493 = qJD(9) + t501;
t348 = Ifges(10,4) * t418 + Ifges(10,2) * t416 + Ifges(10,6) * t493;
t350 = Ifges(10,1) * t418 + Ifges(10,4) * t416 + Ifges(10,5) * t493;
t490 = qJDD(9) + t498;
t538 = -mrSges(10,1) * t305 + mrSges(10,2) * t306 - Ifges(10,5) * t321 - Ifges(10,6) * t319 - Ifges(10,3) * t490 - t418 * t348 + t416 * t350;
t357 = -t416 * mrSges(10,1) + t418 * mrSges(10,2);
t398 = -t493 * mrSges(10,2) + t416 * mrSges(10,3);
t281 = m(10) * t305 + t490 * mrSges(10,1) - t321 * mrSges(10,3) - t418 * t357 + t493 * t398;
t400 = t493 * mrSges(10,1) - t418 * mrSges(10,3);
t282 = m(10) * t306 - t490 * mrSges(10,2) + t319 * mrSges(10,3) + t416 * t357 - t493 * t400;
t546 = -t517 * t281 - t508 * t282;
t532 = -pkin(2) * t546 - mrSges(9,1) * t396 + mrSges(9,2) * t397 - Ifges(9,5) * t393 - Ifges(9,6) * t390 - Ifges(9,3) * t498 - t455 * t406 + t452 * t409 + t538;
t423 = -t457 * mrSges(4,1) + t454 * mrSges(4,2);
t437 = t503 * mrSges(4,1) - t454 * mrSges(4,3);
t229 = m(4) * t370 - t500 * mrSges(4,2) + t395 * mrSges(4,3) + t522 * t237 - t513 * t252 + t457 * t423 - t503 * t437;
t440 = -t503 * mrSges(4,2) + t457 * mrSges(4,3);
t230 = m(4) * t372 + t500 * mrSges(4,1) - t392 * mrSges(4,3) - t454 * t423 + t503 * t440 + t558;
t425 = -t453 * mrSges(8,1) + t456 * mrSges(8,2);
t436 = -t502 * mrSges(8,2) + t453 * mrSges(8,3);
t249 = m(8) * t369 + t499 * mrSges(8,1) - t394 * mrSges(8,3) - t456 * t425 + t502 * t436 + t557;
t439 = t502 * mrSges(8,1) - t456 * mrSges(8,3);
t250 = m(8) * t371 - t499 * mrSges(8,2) + t391 * mrSges(8,3) + t453 * t425 - t502 * t439 + t556;
t559 = -t523 * t229 + t514 * t230 + t519 * t249 + t510 * t250;
t563 = -t559 * pkin(1) - mrSges(3,1) * t444 + mrSges(3,2) * t442 - Ifges(3,5) * t472 - Ifges(3,6) * t470 - Ifges(3,3) * qJDD(2) + t530 + t531 + t532;
t239 = t521 * t260 + t512 * t261;
t511 = sin(qJ(6));
t555 = qJD(1) * t511;
t554 = qJD(1) * t515;
t520 = cos(qJ(6));
t553 = qJD(1) * t520;
t552 = qJD(1) * t524;
t550 = qJD(1) * qJD(6);
t466 = t526 * pkin(14) + t479;
t441 = -t520 * g(3) - t511 * t466;
t467 = (-mrSges(7,1) * t520 + mrSges(7,2) * t511) * qJD(1);
t469 = t511 * qJDD(1) + t520 * t550;
t476 = -qJD(6) * mrSges(7,2) + mrSges(7,3) * t553;
t367 = m(7) * t441 + qJDD(6) * mrSges(7,1) - t469 * mrSges(7,3) + qJD(6) * t476 - t467 * t555;
t443 = -t511 * g(3) + t520 * t466;
t471 = t520 * qJDD(1) - t511 * t550;
t474 = qJD(6) * mrSges(7,1) - mrSges(7,3) * t555;
t368 = m(7) * t443 - qJDD(6) * mrSges(7,2) + t471 * mrSges(7,3) - qJD(6) * t474 + t467 * t553;
t547 = -t511 * t367 + t520 * t368;
t303 = (-t391 * t506 - t394 * t505 + (-t453 * t505 + t456 * t506) * t502) * pkin(4) + t426;
t270 = m(11) * t303 - t300 * mrSges(11,1) + t301 * mrSges(11,2) - t373 * t365 + t374 * t366;
t448 = Ifges(7,6) * qJD(6) + (Ifges(7,4) * t511 + Ifges(7,2) * t520) * qJD(1);
t450 = Ifges(7,5) * qJD(6) + (Ifges(7,1) * t511 + Ifges(7,4) * t520) * qJD(1);
t545 = t511 * t448 - t520 * t450;
t449 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t524 - Ifges(3,2) * t515) * qJD(1);
t451 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t524 - Ifges(3,4) * t515) * qJD(1);
t544 = t524 * t449 + t515 * t451;
t475 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t554;
t477 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t552;
t543 = -t515 * t475 - t524 * t477;
t464 = qJDD(1) * pkin(14) - t478;
t542 = -m(7) * t464 + t471 * mrSges(7,1) - t469 * mrSges(7,2) + t476 * t553;
t541 = m(5) * t341 - t320 * mrSges(5,1) + t322 * mrSges(5,2) - t417 * t399 + t419 * t401 + t239;
t356 = (t455 * t501 - t390) * pkin(2) + t463;
t539 = m(10) * t356 - t319 * mrSges(10,1) + t321 * mrSges(10,2) - t416 * t398 + t418 * t400;
t536 = mrSges(7,1) * t441 - mrSges(7,2) * t443 + Ifges(7,5) * t469 + Ifges(7,6) * t471 + Ifges(7,3) * qJDD(6);
t347 = Ifges(5,5) * t419 + Ifges(5,6) * t417 + Ifges(5,3) * t494;
t224 = -pkin(11) * t239 + mrSges(5,2) * t341 - mrSges(5,3) * t293 + Ifges(5,1) * t322 + Ifges(5,4) * t320 + Ifges(5,5) * t491 - t512 * t244 + t521 * t245 + t417 * t347 - t494 * t349;
t533 = mrSges(6,1) * t271 - mrSges(6,2) * t272 + Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t317 + t379 * t308 - t378 * t309;
t225 = -pkin(9) * t239 - mrSges(5,1) * t341 + mrSges(5,3) * t294 + Ifges(5,4) * t322 + Ifges(5,2) * t320 + Ifges(5,6) * t491 - t419 * t347 + t494 * t351 - t533;
t405 = Ifges(4,5) * t454 + Ifges(4,6) * t457 + Ifges(4,3) * t503;
t219 = -pkin(5) * t541 - mrSges(4,1) * t426 + mrSges(4,3) * t370 + Ifges(4,4) * t392 + Ifges(4,2) * t395 + Ifges(4,6) * t500 + t513 * t224 + t522 * t225 - t454 * t405 + t503 * t411;
t222 = mrSges(4,2) * t426 - mrSges(4,3) * t372 + Ifges(4,1) * t392 + Ifges(4,4) * t395 + Ifges(4,5) * t500 + t522 * t224 - t513 * t225 + t457 * t405 - t503 * t408;
t325 = Ifges(11,5) * t374 + Ifges(11,6) * t373 + Ifges(11,3) * t492;
t262 = -mrSges(11,1) * t303 + mrSges(11,3) * t288 + Ifges(11,4) * t301 + Ifges(11,2) * t300 + Ifges(11,6) * t485 - t374 * t325 + t492 * t327;
t263 = mrSges(11,2) * t303 - mrSges(11,3) * t287 + Ifges(11,1) * t301 + Ifges(11,4) * t300 + Ifges(11,5) * t485 + t373 * t325 - t492 * t326;
t404 = Ifges(8,5) * t456 + Ifges(8,6) * t453 + Ifges(8,3) * t502;
t233 = -mrSges(8,1) * t426 + mrSges(8,3) * t371 + Ifges(8,4) * t394 + Ifges(8,2) * t391 + Ifges(8,6) * t499 + t460 * t262 + t459 * t263 - t270 * t561 - t456 * t404 + t502 * t410;
t234 = mrSges(8,2) * t426 - mrSges(8,3) * t369 + Ifges(8,1) * t394 + Ifges(8,4) * t391 + Ifges(8,5) * t499 - t459 * t262 + t460 * t263 - t270 * t562 + t453 * t404 - t502 * t407;
t346 = Ifges(10,5) * t418 + Ifges(10,6) * t416 + Ifges(10,3) * t493;
t279 = -mrSges(10,1) * t356 + mrSges(10,3) * t306 + Ifges(10,4) * t321 + Ifges(10,2) * t319 + Ifges(10,6) * t490 - t418 * t346 + t493 * t350;
t280 = mrSges(10,2) * t356 - mrSges(10,3) * t305 + Ifges(10,1) * t321 + Ifges(10,4) * t319 + Ifges(10,5) * t490 + t416 * t346 - t493 * t348;
t403 = Ifges(9,5) * t455 + Ifges(9,6) * t452 + Ifges(9,3) * t501;
t241 = -pkin(2) * t539 - mrSges(9,1) * t463 + mrSges(9,3) * t397 + Ifges(9,4) * t393 + Ifges(9,2) * t390 + Ifges(9,6) * t498 - t517 * t279 - t508 * t280 - t455 * t403 + t501 * t409;
t246 = mrSges(9,2) * t463 - mrSges(9,3) * t396 + Ifges(9,1) * t393 + Ifges(9,4) * t390 + Ifges(9,5) * t498 + t508 * t279 - t517 * t280 + t452 * t403 - t501 * t406;
t447 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t524 - Ifges(3,6) * t515) * qJD(1);
t529 = -t541 + t453 * t436 + t457 * t440 + (-m(4) - m(8)) * t426 + t391 * mrSges(8,1) + t395 * mrSges(4,1) - t454 * t437 - t456 * t439 - t392 * mrSges(4,2) - t394 * mrSges(8,2) - t270;
t215 = pkin(1) * t529 - mrSges(3,1) * t463 + mrSges(3,3) * t442 + Ifges(3,4) * t472 + Ifges(3,2) * t470 + Ifges(3,6) * qJDD(2) + qJD(2) * t451 + t514 * t219 - t523 * t222 + t519 * t233 + t510 * t234 + t518 * t241 + t509 * t246 - t447 * t552;
t217 = mrSges(3,2) * t463 - mrSges(3,3) * t444 + Ifges(3,1) * t472 + Ifges(3,4) * t470 + Ifges(3,5) * qJDD(2) - qJD(2) * t449 + t523 * t219 + t514 * t222 - t510 * t233 + t519 * t234 - t509 * t241 + t518 * t246 - t447 * t554;
t446 = Ifges(7,3) * qJD(6) + (Ifges(7,5) * t511 + Ifges(7,6) * t520) * qJD(1);
t337 = -mrSges(7,1) * t464 + mrSges(7,3) * t443 + Ifges(7,4) * t469 + Ifges(7,2) * t471 + Ifges(7,6) * qJDD(6) + qJD(6) * t450 - t446 * t555;
t338 = mrSges(7,2) * t464 - mrSges(7,3) * t441 + Ifges(7,1) * t469 + Ifges(7,4) * t471 + Ifges(7,5) * qJDD(6) - qJD(6) * t448 + t446 * t553;
t361 = -t474 * t555 + t542;
t435 = -t501 * mrSges(9,2) + t452 * mrSges(9,3);
t438 = t501 * mrSges(9,1) - t455 * mrSges(9,3);
t528 = t470 * mrSges(3,1) + t529 + t390 * mrSges(9,1) + t452 * t435 + (-m(3) - m(9)) * t463 - t472 * mrSges(3,2) - t455 * t438 - t393 * mrSges(9,2) - t539;
t535 = -mrSges(2,2) * t479 - pkin(14) * t361 - t515 * t215 + t524 * t217 + pkin(15) * (qJD(1) * t543 + t528) + t511 * t338 + t520 * t337 + mrSges(2,1) * t478 + Ifges(2,3) * qJDD(1);
t468 = (mrSges(3,1) * t515 + mrSges(3,2) * t524) * qJD(1);
t424 = -t452 * mrSges(9,1) + t455 * mrSges(9,2);
t265 = m(9) * t397 - t498 * mrSges(9,2) + t390 * mrSges(9,3) + t508 * t281 - t517 * t282 + t452 * t424 - t501 * t438;
t264 = m(9) * t396 + t498 * mrSges(9,1) - t393 * mrSges(9,3) - t455 * t424 + t501 * t435 + t546;
t226 = qJDD(1) * mrSges(2,1) + (-t511 * t474 + t543) * qJD(1) - t526 * mrSges(2,2) + m(2) * t478 + t528 + t542;
t221 = m(3) * t444 + qJDD(2) * mrSges(3,1) - t472 * mrSges(3,3) + qJD(2) * t475 + t518 * t264 + t509 * t265 - t468 * t552 + t559;
t220 = m(3) * t442 - qJDD(2) * mrSges(3,2) + t470 * mrSges(3,3) - qJD(2) * t477 + t514 * t229 + t523 * t230 - t510 * t249 + t519 * t250 - t509 * t264 + t518 * t265 - t468 * t554;
t218 = m(2) * t479 - t526 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t515 * t220 - t524 * t221 + t547;
t214 = mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-t544 - t545) * qJD(1) - pkin(15) * (t524 * t220 - t515 * t221) + t526 * Ifges(2,5) + pkin(14) * (t520 * t367 + t511 * t368) - pkin(16) * t547 + mrSges(2,3) * t479 - t536 + t563;
t213 = -mrSges(2,2) * g(3) - mrSges(2,3) * t478 + pkin(16) * t361 + Ifges(2,5) * qJDD(1) - t526 * Ifges(2,6) - t524 * t215 - t515 * t217 - t511 * t337 + t520 * t338;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t525 * t213 - t516 * t214 - pkin(13) * (t516 * t218 + t525 * t226), t213, t217, t222, t224, t245, t338, t234, t246, t280, t263, 0, 0, 0, 0, 0, 0; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t516 * t213 + t525 * t214 + pkin(13) * (t525 * t218 - t516 * t226), t214, t215, t219, t225, t244, t337, t233, t241, t279, t262, 0, 0, 0, 0, 0, 0; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t535, t535, t544 * qJD(1) - t563, -t530, -t534, t533, qJD(1) * t545 + t536, -t531, -t532, -t538, -t537, 0, 0, 0, 0, 0, 0;];
m_new = t1;
