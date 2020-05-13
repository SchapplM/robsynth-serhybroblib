% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh3m2TE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2TE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:47:34
% EndTime: 2020-05-07 01:47:40
% DurationCPUTime: 5.48s
% Computational Cost: add. (9404->340), mult. (18210->544), div. (0->0), fcn. (22369->20), ass. (0->194)
t464 = sin(qJ(3));
t569 = pkin(1) * qJD(2);
t512 = t464 * t569;
t465 = sin(qJ(2));
t468 = cos(qJ(2));
t460 = cos(pkin(16));
t469 = cos(pkin(15));
t547 = sin(pkin(16));
t575 = sin(pkin(15));
t429 = t460 * t575 + t547 * t469;
t430 = t460 * t469 - t547 * t575;
t467 = cos(qJ(3));
t487 = t429 * t464 - t430 * t467;
t505 = -t429 * t467 - t464 * t430;
t359 = t465 * t487 + t468 * t505;
t456 = pkin(17) + pkin(18);
t451 = sin(t456);
t452 = cos(t456);
t608 = t487 * t468;
t595 = -t505 * t465 + t608;
t606 = t359 * t452 + t451 * t595;
t344 = t359 * t451 - t595 * t452;
t455 = qJD(3) + qJD(2);
t572 = pkin(4) * t455;
t483 = t467 * t569 - t572;
t617 = t344 * t483;
t618 = -t606 * t512 + t617;
t511 = t344 * t569;
t340 = t464 * t511;
t613 = t606 * t483;
t327 = -t340 - t613;
t382 = t505 * qJD(3);
t615 = t595 * qJD(2) - t382 * t465;
t381 = t487 * qJD(3);
t345 = -t381 * t468 - t615;
t346 = t359 * qJD(2) + t381 * t465 + t382 * t468;
t335 = -t345 * t451 + t346 * t452;
t449 = -pkin(1) * t467 + pkin(4);
t537 = t346 * t451;
t538 = t606 * t467;
t478 = (-qJD(3) * t538 + (t345 * t452 + t537) * t464) * pkin(1);
t574 = pkin(1) * t464;
t513 = t344 * t574;
t523 = -qJD(3) * t513 - t335 * t449 + t327 + t478;
t612 = mrSges(9,1) + mrSges(4,1);
t611 = mrSges(4,2) + mrSges(9,2);
t470 = cos(pkin(14));
t576 = sin(pkin(14));
t435 = -t469 * t576 + t470 * t575;
t436 = t470 * t469 + t576 * t575;
t398 = t435 * t468 + t436 * t465;
t610 = Ifges(7,4) * t398;
t462 = cos(pkin(18));
t570 = sin(pkin(18));
t431 = t462 * t575 + t570 * t469;
t432 = t462 * t469 - t570 * t575;
t459 = qJ(3) + qJ(2);
t453 = sin(t459);
t454 = cos(t459);
t607 = mrSges(8,1) * t432 + mrSges(9,1) * t454 + mrSges(8,2) * t431 - mrSges(9,2) * t453;
t486 = t453 * t464 + t454 * t467;
t375 = t429 * t452 + t430 * t451;
t370 = t375 * qJD(1);
t376 = t429 * t451 - t430 * t452;
t371 = t376 * qJD(1);
t433 = t464 * t465 - t467 * t468;
t425 = t433 * qJD(1);
t447 = -pkin(1) * t468 - pkin(12);
t405 = -pkin(4) * t425 + t447 * qJD(1);
t349 = -pkin(8) * t371 - pkin(10) * t370 + t405;
t463 = sin(qJ(4));
t466 = cos(qJ(4));
t321 = -t327 * t463 + t349 * t466;
t322 = t327 * t466 + t349 * t463;
t605 = (-t321 * t466 - t322 * t463) * qJD(4);
t504 = -t435 * t465 + t436 * t468;
t587 = -t504 / 0.2e1;
t604 = Ifges(7,4) * t504;
t368 = (cos(pkin(17)) * t432 - t431 * sin(pkin(17))) * pkin(3) + t447;
t602 = t368 * (-mrSges(9,1) * t453 - mrSges(9,2) * t454);
t485 = t464 * t468 + t465 * t467;
t426 = t485 * qJD(1);
t599 = (-mrSges(4,1) * t425 - mrSges(4,2) * t426 + t607 * qJD(1)) * t465;
t516 = t453 * qJD(1);
t567 = mrSges(4,3) * t426;
t598 = (mrSges(9,3) * t516 + t612 * t455 + t567) * t464;
t597 = (m(4) + m(8)) * t447;
t334 = (qJD(3) * t608 + t615) * t452 - t537;
t476 = -t335 * t464 + (-t344 * t467 + t464 * t606) * qJD(3);
t316 = t476 * pkin(1) + t334 * t449;
t596 = t316 - t618;
t509 = qJD(3) * t569;
t501 = t464 * t509;
t313 = qJD(2) * t478 + t335 * t483 - t344 * t501;
t330 = -t344 * t449 - t574 * t606;
t594 = t313 * t330 + t523 * t618;
t491 = -t321 * t463 + t322 * t466;
t314 = t334 * t572 + (-t334 * t467 + t476) * t569;
t402 = t455 * t485;
t386 = t402 * qJD(1);
t510 = qJD(1) * t569;
t367 = -pkin(4) * t386 + t465 * t510;
t311 = qJD(4) * t321 + t314 * t466 + t367 * t463;
t312 = -qJD(4) * t322 - t314 * t463 + t367 * t466;
t593 = t311 * t466 - t312 * t463;
t592 = -0.2e1 * pkin(12);
t591 = qJD(1) ^ 2;
t387 = t398 * qJD(2);
t589 = -0.3e1 / 0.2e1 * t387;
t388 = t504 * qJD(2);
t588 = 0.3e1 / 0.2e1 * t388;
t586 = t398 / 0.2e1;
t401 = t455 * t433;
t585 = t401 / 0.2e1;
t584 = t402 / 0.2e1;
t582 = t425 / 0.2e1;
t581 = -t426 / 0.2e1;
t579 = t463 / 0.2e1;
t578 = -t465 / 0.2e1;
t577 = t468 / 0.2e1;
t573 = pkin(4) * t426;
t571 = m(9) * t368;
t568 = mrSges(4,3) * t425;
t566 = mrSges(6,3) * t370;
t565 = Ifges(7,1) * t398;
t564 = Ifges(3,4) * t465;
t563 = Ifges(3,4) * t468;
t562 = Ifges(4,4) * t426;
t561 = Ifges(6,4) * t463;
t560 = Ifges(6,4) * t466;
t559 = Ifges(9,4) * t453;
t558 = Ifges(9,4) * t454;
t369 = qJD(4) - t371;
t557 = Ifges(6,5) * t369;
t556 = Ifges(7,5) * t388;
t555 = Ifges(9,5) * t454;
t554 = Ifges(9,5) * t455;
t553 = Ifges(3,2) * t468;
t552 = Ifges(7,2) * t504;
t551 = Ifges(6,6) * t369;
t550 = Ifges(7,6) * t387;
t549 = Ifges(9,6) * t453;
t548 = Ifges(9,6) * t455;
t320 = t613 - t340 + (t344 * t464 - t538) * t569;
t544 = t320 * t618;
t524 = t466 * t370;
t520 = qJD(1) * t454;
t522 = mrSges(9,3) * t520 + t611 * t455 - t568;
t484 = -pkin(4) * t433 + t447;
t350 = -pkin(8) * t376 - pkin(10) * t375 + t484;
t518 = qJD(4) * t350;
t517 = qJD(4) * t370;
t515 = t465 * qJD(1);
t507 = t376 * t517;
t499 = t312 * mrSges(6,1) - t311 * mrSges(6,2);
t498 = mrSges(6,1) * t466 - mrSges(6,2) * t463;
t492 = -t313 * t344 - t335 * t618;
t353 = -mrSges(6,2) * t369 - t463 * t566;
t354 = mrSges(6,1) * t369 - mrSges(6,3) * t524;
t490 = t466 * t353 - t463 * t354;
t395 = t431 * t468 + t432 * t465;
t383 = t395 * qJD(2);
t394 = -t465 * t431 + t432 * t468;
t384 = t394 * qJD(2);
t488 = -t383 * t431 - t384 * t432;
t482 = Ifges(3,5) * t577 + Ifges(3,6) * t578;
t450 = qJD(1) * t563;
t481 = (Ifges(3,6) * qJD(2) + (t553 + t564) * qJD(1)) * t578 + (Ifges(3,1) * t515 + Ifges(3,5) * qJD(2) + t450) * t577;
t479 = t486 * t455;
t477 = -qJD(1) * t602 - t447 * (-mrSges(4,1) * t426 + mrSges(4,2) * t425) - t454 * t554;
t475 = qJD(4) * (-t353 * t463 - t354 * t466 + (-t463 ^ 2 - t466 ^ 2) * t566);
t373 = Ifges(4,2) * t425 + Ifges(4,6) * t455 - t562;
t415 = Ifges(4,4) * t425;
t374 = -Ifges(4,1) * t426 + Ifges(4,5) * t455 + t415;
t385 = t401 * qJD(1);
t413 = t548 + (-Ifges(9,2) * t454 - t559) * qJD(1);
t414 = t554 + (-Ifges(9,1) * t453 - t558) * qJD(1);
t474 = t373 * t581 + t426 * (Ifges(4,1) * t425 + t562) / 0.2e1 + Ifges(4,6) * t386 + Ifges(4,5) * t385 + t414 * t520 / 0.2e1 + t512 * t567 - (Ifges(4,2) * t426 + t374 + t415) * t425 / 0.2e1 - (Ifges(4,5) * t425 + Ifges(4,6) * t426 + (t549 - t555) * qJD(1)) * t455 / 0.2e1 + (-t413 / 0.2e1 + t548) * t516 + (t453 * (-Ifges(9,1) * t454 + t559) + t454 * (Ifges(9,2) * t453 - t558)) * t591 / 0.2e1 + t486 * mrSges(9,3) * t510 + t612 * t501 + t611 * t467 * t509;
t408 = pkin(1) * t515 - t573;
t377 = -pkin(4) * t402 + t465 * t569;
t358 = Ifges(7,5) * qJD(2) + qJD(1) * (t565 + t604);
t357 = Ifges(7,6) * qJD(2) + qJD(1) * (t552 + t610);
t355 = -mrSges(5,1) * t371 + mrSges(5,2) * t370;
t352 = (mrSges(6,1) * t463 + mrSges(6,2) * t466) * t370;
t351 = t498 * t517;
t348 = t557 + (t466 * Ifges(6,1) - t561) * t370;
t347 = t551 + (-t463 * Ifges(6,2) + t560) * t370;
t331 = t449 * t606 - t513;
t324 = t408 * t463 + t466 * t618;
t323 = t408 * t466 - t463 * t618;
t319 = -t467 * t511 + t617;
t318 = t319 * t466 - t463 * t573;
t317 = -t319 * t463 - t466 * t573;
t1 = [t377 * t355 + t373 * t584 - t387 * t357 / 0.2e1 + t388 * t358 / 0.2e1 + t374 * t585 + t447 * (-mrSges(4,1) * t386 + mrSges(4,2) * t385) + (t367 * mrSges(5,2) + t313 * mrSges(5,3)) * t375 + (t433 * t386 + t402 * t582) * Ifges(4,2) + (-t385 * t485 + t401 * t581) * Ifges(4,1) + (t367 * t484 + t405 * t377) * m(5) + (Ifges(4,5) * t585 + Ifges(4,6) * t584 - t454 * t414 / 0.2e1 + t453 * t413 / 0.2e1 + (-t555 / 0.2e1 + t549 / 0.2e1) * t455) * t455 + (-t367 * mrSges(5,1) + t314 * mrSges(5,3) - t499) * t376 + (t433 * t385 - t386 * t485 + t401 * t582 + t402 * t581) * Ifges(4,4) + (t556 / 0.2e1 - t550 / 0.2e1 + t482 * qJD(2) + (t488 * mrSges(8,3) + t599 + (t486 * qJD(3) - t479) * mrSges(9,3) + (t401 * t467 - t402 * t464 + (-t433 * t467 + t464 * t485) * qJD(3)) * mrSges(4,3)) * pkin(1) + t481) * qJD(2) + (Ifges(6,6) * t507 + m(6) * (t312 * t350 + t321 * t377 + t322 * t518) + t353 * t518 + t377 * t354 + (t313 * mrSges(6,2) - t312 * mrSges(6,3) + (-t551 / 0.2e1 - t322 * mrSges(6,3) - t347 / 0.2e1 + t618 * mrSges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,4) * t524) * qJD(4)) * t375) * t466 + (-t354 * t518 + Ifges(6,5) * t507 + m(6) * (t311 * t350 - t321 * t518 + t322 * t377) + t377 * t353 + (t313 * mrSges(6,1) - t311 * mrSges(6,3) + (-t348 / 0.2e1 + t321 * mrSges(6,3) - t618 * mrSges(6,2) - t557 / 0.2e1 + (0.3e1 / 0.2e1 * t561 + (0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(6,1)) * t466) * t370) * qJD(4)) * t375) * t463 + (0.2e1 * pkin(6) * (mrSges(7,1) * t387 + mrSges(7,2) * t388) + t565 * t588 + t552 * t589 + t447 * (-mrSges(4,1) * t402 + mrSges(4,2) * t401) + (t398 * t589 + t504 * t588) * Ifges(7,4) + ((mrSges(3,2) * t592 + 0.3e1 / 0.2e1 * t563) * t468 + (mrSges(3,1) * t592 - 0.3e1 / 0.2e1 * t564 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t468 + (-mrSges(4,1) * t433 - mrSges(4,2) * t485 + 0.2e1 * t571 + 0.2e1 * t597 + t607) * pkin(1)) * t465) * qJD(2) + (0.2e1 * t602 + (0.3e1 / 0.2e1 * Ifges(9,1) - 0.3e1 / 0.2e1 * Ifges(9,2)) * t453 * t454 + (0.3e1 / 0.2e1 * t454 ^ 2 - 0.3e1 / 0.2e1 * t453 ^ 2) * Ifges(9,4)) * t455) * qJD(1); t330 * t351 + t490 * t316 - t408 * t355 - t324 * t353 - t323 * t354 + t474 + t523 * t352 + (t523 * t370 + t596 * t371) * mrSges(5,3) + t331 * t475 + (-t468 * t450 / 0.2e1 + t358 * t587 + t357 * t586 - t550 + t556 + (pkin(12) * (mrSges(3,1) * t465 + mrSges(3,2) * t468) + (Ifges(3,1) * t468 - t564) * t578 + t465 * t553 / 0.2e1 - pkin(6) * (mrSges(7,1) * t398 + mrSges(7,2) * t504) - t398 * (Ifges(7,1) * t504 - t610) / 0.2e1 + (-Ifges(7,2) * t398 + t604) * t587) * qJD(1) + (Ifges(7,5) * t587 + Ifges(7,6) * t586 + t482) * qJD(2) + t477 - t481) * qJD(1) + (t491 * t316 + (t605 + t593) * t331 - t321 * t323 - t322 * t324 + t594) * m(6) + (t314 * t331 + t596 * t327 - t405 * t408 + t594) * m(5) + ((-t571 - t597) * t591 * t465 + (-t464 * t386 + (-qJD(2) * t425 + t385) * t467) * mrSges(4,3) + (t522 * t467 + t598) * qJD(3) + (-mrSges(9,3) * t479 - t599 + ((t394 * t432 + t395 * t431) * qJD(2) + t488) * mrSges(8,3)) * qJD(1) + 0.2e1 * m(8) * (-t383 * t394 + t384 * t395) * t569) * pkin(1); -t320 * t352 - m(5) * (t319 * t327 + t544) + t477 * qJD(1) - t318 * t353 - t317 * t354 + t474 - m(6) * (t317 * t321 + t318 * t322 + t544) + (-t598 + (-t522 - t568) * t467) * t569 + (-t319 * t371 - t320 * t370) * mrSges(5,3) + (-t344 * t351 + t426 * t355 + (-t370 * mrSges(5,3) - t352) * t335 + m(6) * (t605 * t606 + t492) + (m(6) * t593 + t475) * t606 + (m(6) * t491 + t371 * mrSges(5,3) + t490) * t334 + (t314 * t606 + t327 * t334 + t405 * t426 + t492) * m(5)) * pkin(4); -t321 * t353 + t322 * t354 + (t466 * t347 / 0.2e1 + t348 * t579 - t618 * t498 + ((-Ifges(6,2) * t466 - t561) * t579 - t466 * (-Ifges(6,1) * t463 - t560) / 0.2e1) * t370 + t491 * mrSges(6,3) + (-t369 / 0.2e1 + qJD(4)) * (-Ifges(6,5) * t463 - Ifges(6,6) * t466)) * t370 + t499;];
tauc = t1(:);
