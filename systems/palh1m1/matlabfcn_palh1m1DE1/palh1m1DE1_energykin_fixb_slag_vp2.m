% Calculate kinetic energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1DE1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1DE1_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:48:29
% EndTime: 2020-04-13 14:49:17
% DurationCPUTime: 48.34s
% Computational Cost: add. (1179204->307), mult. (1819709->590), div. (76812->33), fcn. (1137810->44), ass. (0->246)
t605 = -2 * pkin(1);
t604 = (-pkin(2) - pkin(13));
t603 = pkin(7) - pkin(8);
t602 = pkin(7) + pkin(8);
t601 = -pkin(9) - pkin(11);
t600 = pkin(11) - pkin(9);
t599 = (pkin(13) - pkin(2));
t512 = pkin(4) ^ 2;
t511 = pkin(5) ^ 2;
t509 = pkin(7) ^ 2;
t517 = pkin(1) ^ 2;
t490 = sin(qJ(2));
t491 = sin(pkin(19));
t495 = cos(qJ(2));
t496 = cos(pkin(19));
t470 = t490 * t496 - t491 * t495;
t579 = pkin(7) * t470;
t552 = t579 * t605 + t517;
t456 = t509 + t552;
t549 = pkin(3) ^ 2 - pkin(8) ^ 2;
t450 = t456 + t549;
t461 = pkin(1) - t579;
t471 = t490 * t491 + t495 * t496;
t444 = (pkin(3) + t602) * (-pkin(3) + t603) + t552;
t445 = (-pkin(3) + t602) * (pkin(3) + t603) + t552;
t520 = sqrt(-t445 * t444);
t429 = pkin(7) * t450 * t471 + t461 * t520;
t494 = cos(qJ(3));
t556 = t494 * t429;
t563 = t471 * t520;
t428 = -pkin(7) * t563 + t450 * t461;
t489 = sin(qJ(3));
t559 = t489 * t428;
t524 = t559 / 0.2e1 + t556 / 0.2e1;
t454 = 0.1e1 / t456;
t514 = 0.1e1 / pkin(3);
t565 = t454 * t514;
t419 = t524 * t565;
t557 = t494 * t428;
t558 = t489 * t429;
t523 = -t557 / 0.2e1 + t558 / 0.2e1;
t420 = t523 * t565;
t478 = pkin(23) + pkin(22);
t476 = sin(t478);
t477 = cos(t478);
t398 = -t419 * t477 + t420 * t476;
t584 = pkin(5) * t398;
t554 = -0.2e1 * pkin(4) * t584 + t511;
t394 = t512 + t554;
t392 = 0.1e1 / t394;
t598 = t392 / 0.2e1;
t467 = t471 * qJD(2);
t466 = t470 * qJD(2);
t547 = pkin(1) * pkin(7) * t467;
t569 = 0.2e1 * (t444 + t445) * t547 / t520;
t537 = -t569 / 0.2e1;
t521 = t466 * t520 + t471 * t537;
t409 = ((t461 * t605 - t450) * t467 + t521) * pkin(7);
t597 = -t409 / 0.2e1;
t544 = -0.2e1 * t467 * t471;
t564 = t467 * t520;
t410 = t461 * t569 / 0.2e1 + t509 * pkin(1) * t544 + (-t450 * t466 - t564) * pkin(7);
t596 = t410 / 0.2e1;
t483 = sin(pkin(20));
t487 = cos(pkin(20));
t525 = t483 * t494 + t487 * t489;
t581 = pkin(6) * t525;
t545 = pkin(1) * t581;
t458 = 0.2e1 * t545;
t510 = pkin(6) ^ 2;
t550 = t510 + t517;
t453 = t458 + t550;
t451 = 0.1e1 / t453;
t595 = t451 / 0.2e1;
t480 = sin(pkin(23));
t594 = t480 / 0.2e1;
t482 = sin(pkin(21));
t593 = t482 / 0.2e1;
t592 = -t494 / 0.2e1;
t497 = cos(pkin(18));
t591 = t497 / 0.2e1;
t516 = 0.1e1 / pkin(2);
t590 = t516 / 0.2e1;
t551 = pkin(9) ^ 2 - pkin(11) ^ 2;
t390 = t394 + t551;
t395 = -pkin(4) + t584;
t399 = t419 * t476 + t420 * t477;
t386 = (pkin(4) - t601) * (pkin(4) + t601) + t554;
t387 = (pkin(4) - t600) * (pkin(4) + t600) + t554;
t518 = sqrt(-t387 * t386);
t572 = t399 * t518;
t363 = -pkin(5) * t572 - t390 * t395;
t583 = pkin(5) * t399;
t365 = t390 * t583 - t395 * t518;
t506 = 0.1e1 / pkin(9);
t539 = t506 * t598;
t355 = atan2(t365 * t539, t363 * t539);
t589 = sin(t355);
t469 = t483 * t489 - t487 * t494;
t460 = t469 * qJD(3);
t588 = pkin(1) * t460;
t484 = cos(pkin(23));
t561 = t484 * t428;
t562 = t480 * t429;
t416 = (-t561 / 0.2e1 + t562 / 0.2e1) * t565;
t415 = 0.1e1 / t416 ^ 2;
t560 = t484 * t429;
t571 = t428 * t480;
t417 = (t560 / 0.2e1 + t571 / 0.2e1) * t565;
t533 = 0.1e1 / t456 ^ 2 * t547;
t360 = qJD(2) + (-((t410 * t594 + t484 * t597) * t454 + (-t561 + t562) * t533) * t417 * t415 + ((t409 * t594 + t484 * t596) * t454 + (t560 + t571) * t533) / t416) / (t415 * t417 ^ 2 + 0.1e1) * t514;
t587 = pkin(4) * t360;
t376 = ((-t557 + t558) * t533 + (qJD(3) * t524 + t409 * t592 + t489 * t596) * t454) * t514;
t377 = ((-t556 - t559) * t533 + (qJD(3) * t523 + t410 * t592 + t489 * t597) * t454) * t514;
t368 = t376 * t476 + t377 * t477;
t586 = pkin(4) * t368;
t585 = pkin(5) * t365;
t515 = pkin(2) ^ 2;
t542 = -pkin(13) ^ 2 + t550;
t448 = t458 + t515 + t542;
t457 = -pkin(1) - t581;
t553 = t458 + t510;
t442 = ((pkin(1) - t604) * (pkin(1) + t604)) + t553;
t443 = ((pkin(1) - t599) * (pkin(1) + t599)) + t553;
t568 = t443 * t442;
t519 = sqrt(-t568);
t580 = pkin(6) * t469;
t426 = t448 * t580 - t457 * t519;
t582 = pkin(6) * t426;
t578 = pkin(1) * qJD(2);
t577 = qJD(1) * pkin(16);
t546 = pkin(5) * t586;
t576 = 0.2e1 * (t386 + t387) * t546 / t518;
t575 = t368 * t518;
t486 = cos(pkin(21));
t574 = t392 * t486;
t504 = 0.1e1 / pkin(11);
t573 = t392 * t504;
t548 = pkin(6) * t588;
t570 = 0.2e1 * (t442 + t443) * t548 / t519;
t492 = sin(pkin(18));
t567 = t454 * t492;
t508 = 0.1e1 / pkin(8);
t566 = t454 * t508;
t555 = qJD(2) ^ 2 * t517;
t479 = qJD(2) + qJD(3);
t543 = t494 * t578;
t541 = -t576 / 0.2e1;
t540 = t392 * t593;
t538 = -t570 / 0.2e1;
t536 = t451 * t590;
t535 = t454 * t591;
t534 = 0.1e1 / pkin(13) * t590;
t425 = -t448 * t457 - t519 * t580;
t424 = 0.1e1 / t425 ^ 2;
t452 = 0.1e1 / t453 ^ 2;
t459 = t525 * qJD(3);
t372 = qJD(2) + 0.2e1 * (((t457 * t538 + (t448 * t459 - t460 * t519) * pkin(6)) * t595 + (-t451 * t469 * t510 + t452 * t582) * t588) / t425 - ((-t460 * t448 - t459 * t519 + t469 * t538) * t595 + (t425 * t452 + t451 * t457) * t588) * t424 * t582) * pkin(2) / (t424 * t426 ^ 2 + 0.1e1) * t453 * t516;
t393 = 0.1e1 / t394 ^ 2;
t532 = t393 * t546;
t474 = t490 * qJD(1) * pkin(1) - t577;
t449 = t456 - t549;
t462 = pkin(1) * t470 - pkin(7);
t430 = pkin(1) * t449 * t471 - t462 * t520;
t529 = t430 * t533;
t528 = t482 * t532;
t527 = t486 * t532;
t427 = -pkin(1) * t563 - t449 * t462;
t526 = t427 * t533;
t391 = t394 - t551;
t396 = -pkin(4) * t398 + pkin(5);
t364 = -pkin(4) * t572 + t391 * t396;
t366 = pkin(4) * t391 * t399 + t396 * t518;
t352 = (t364 * t593 + t366 * t486 / 0.2e1) * t573;
t353 = (-t364 * t486 / 0.2e1 + t366 * t593) * t573;
t349 = atan2(t352, t353);
t346 = sin(t349);
t347 = cos(t349);
t464 = (t489 * t495 + t490 * t494) * qJD(1);
t465 = (-t489 * t490 + t494 * t495) * qJD(1);
t337 = -t346 * t464 + t347 * t465;
t472 = pkin(5) * t479 + t489 * t578;
t340 = t346 * t472 - t347 * t543;
t441 = -pkin(5) * t465 + t474;
t367 = t376 * t477 - t377 * t476;
t522 = -t367 * t518 + t399 * t541;
t339 = t346 * t543 + t347 * t472;
t493 = cos(qJ(4));
t488 = sin(qJ(4));
t485 = cos(pkin(22));
t481 = sin(pkin(22));
t473 = t474 ^ 2;
t447 = t515 - t542 - 0.2e1 * t545;
t446 = 0.1e1 / t447 ^ 2;
t435 = atan2(t519 * t534, t447 * t534);
t433 = cos(t435);
t432 = sin(t435);
t422 = (t430 * t591 + t427 * t492 / 0.2e1) * t566;
t421 = (t427 * t591 - t492 * t430 / 0.2e1) * t566;
t418 = 0.1e1 / t421 ^ 2;
t414 = atan2(t426 * t536, t425 * t536);
t413 = cos(t414);
t412 = sin(t414);
t411 = t462 * t537 + t517 * pkin(7) * t544 + (-t449 * t466 - t564) * pkin(1);
t408 = ((0.2e1 * pkin(7) * t462 - t449) * t467 + t521) * pkin(1);
t407 = atan2(t422, t421);
t405 = atan2(t417, t416);
t404 = cos(t407);
t403 = sin(t407);
t401 = cos(t405);
t400 = sin(t405);
t389 = (-t412 * t490 + t413 * t495) * qJD(1);
t388 = (-t412 * t495 - t413 * t490) * qJD(1);
t385 = -pkin(2) * t388 - t577;
t381 = (-t400 * t490 + t401 * t495) * qJD(1);
t380 = (-t400 * t495 - t401 * t490) * qJD(1);
t375 = -t388 * t432 - t389 * t433;
t374 = -t388 * t433 + t389 * t432;
t370 = (-t380 * t485 - t381 * t481) * pkin(4) + t474;
t369 = (0.1e1 / t447 * t570 / 0.2e1 - 0.2e1 * t519 * t446 * t548) / (-t446 * t568 + 0.1e1) + t372;
t362 = 0.1e1 / t363 ^ 2;
t361 = ((t411 * t535 + t497 * t529 + t408 * t567 / 0.2e1 + t492 * t526) / t421 - (t408 * t535 + t497 * t526 - t411 * t567 / 0.2e1 - t492 * t529) * t422 * t418) / (t418 * t422 ^ 2 + 0.1e1) * t508;
t359 = t401 * t578 + t485 * t587;
t358 = t400 * t578 + t481 * t587;
t354 = cos(t355);
t351 = 0.1e1 / t353 ^ 2;
t345 = -t485 * t354 - t481 * t589;
t344 = t354 * t481 - t485 * t589;
t342 = t396 * t576 / 0.2e1 - 0.2e1 * t512 * t368 * t583 + (t367 * t391 - t575) * pkin(4);
t341 = ((-0.2e1 * pkin(5) * t396 - t391) * t368 + t522) * pkin(4);
t338 = t346 * t465 + t347 * t464;
t336 = qJD(4) - t337;
t335 = t344 * t380 + t345 * t381;
t334 = -t344 * t381 + t345 * t380;
t333 = t344 * t359 + t345 * t358;
t332 = -t344 * t358 + t345 * t359;
t331 = -pkin(10) * t337 - pkin(12) * t338 + t441;
t330 = t360 + 0.2e1 * (((t395 * t541 + (t367 * t390 - t575) * pkin(5)) * t598 + (-t392 * t399 * t511 + t393 * t585) * t586) / t363 - ((-t368 * t390 + t522) * t598 + (t363 * t393 + t392 * t395) * t586) * t362 * t585) * pkin(9) / (t362 * t365 ^ 2 + 0.1e1) * t394 * t506;
t329 = ((t341 * t540 + t364 * t528 + t342 * t574 / 0.2e1 + t366 * t527) / t353 - (-t341 * t574 / 0.2e1 - t364 * t527 + t342 * t540 + t366 * t528) * t352 * t351) / (t351 * t352 ^ 2 + 0.1e1) * t504 + t479;
t328 = pkin(12) * t329 + t340;
t327 = -pkin(10) * t329 - t339;
t326 = t329 * t488 + t338 * t493;
t325 = t329 * t493 - t338 * t488;
t324 = t328 * t493 + t331 * t488;
t323 = -t328 * t488 + t331 * t493;
t1 = (-t474 * mrSges(8,1) + Ifges(8,4) * t381 + Ifges(8,2) * t380 / 0.2e1) * t380 + m(10) * (t385 ^ 2 + (t432 ^ 2 + t433 ^ 2) * t515 * t372 ^ 2) / 0.2e1 + (t332 * mrSges(11,1) - t333 * mrSges(11,2) + Ifges(11,5) * t335 + Ifges(11,6) * t334 + Ifges(11,3) * t330 / 0.2e1) * t330 + (t339 * mrSges(5,1) - t340 * mrSges(5,2) + Ifges(5,5) * t338 + Ifges(5,6) * t337 + Ifges(5,3) * t329 / 0.2e1) * t329 + (Ifges(9,4) * t389 + Ifges(9,2) * t388 / 0.2e1) * t388 + (-t327 * mrSges(6,1) + t324 * mrSges(6,3) + Ifges(6,4) * t326 + Ifges(6,6) * t336 + Ifges(6,2) * t325 / 0.2e1) * t325 + (Ifges(10,5) * t375 + Ifges(10,6) * t374 + Ifges(10,3) * t369 / 0.2e1) * t369 + Ifges(7,3) * t361 ^ 2 / 0.2e1 + (Ifges(9,5) * t389 + Ifges(9,6) * t388 + Ifges(9,3) * t372 / 0.2e1 + (-t432 * (-mrSges(10,2) * t369 + mrSges(10,3) * t374) - t433 * (mrSges(10,1) * t369 - mrSges(10,3) * t375)) * pkin(2)) * t372 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t400 * (-mrSges(8,2) * t360 + mrSges(8,3) * t380) + t401 * (mrSges(8,1) * t360 - mrSges(8,3) * t381) + t489 * (mrSges(4,1) * t479 - mrSges(4,3) * t464) - t494 * (-mrSges(4,2) * t479 + mrSges(4,3) * t465)) * pkin(1)) * qJD(2) + (t474 * mrSges(8,2) + Ifges(8,1) * t381 / 0.2e1) * t381 + (-t370 * mrSges(11,1) + t333 * mrSges(11,3) + Ifges(11,4) * t335 + Ifges(11,2) * t334 / 0.2e1) * t334 + (t385 * mrSges(10,2) + Ifges(10,1) * t375 / 0.2e1) * t375 + Ifges(9,1) * t389 ^ 2 / 0.2e1 + (t370 * mrSges(11,2) - t332 * mrSges(11,3) + Ifges(11,1) * t335 / 0.2e1) * t335 + (t441 * mrSges(5,2) - t339 * mrSges(5,3) + Ifges(5,1) * t338 / 0.2e1) * t338 + (-t474 * mrSges(4,1) + Ifges(4,6) * t479 + Ifges(4,2) * t465 / 0.2e1) * t465 + m(4) * (t473 + (t489 ^ 2 + t494 ^ 2) * t555) / 0.2e1 + m(8) * (t473 + (t400 ^ 2 + t401 ^ 2) * t555) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t441 ^ 2) / 0.2e1 + (-t441 * mrSges(5,1) + t340 * mrSges(5,3) + Ifges(5,4) * t338 + Ifges(5,2) * t337 / 0.2e1) * t337 + m(11) * (t332 ^ 2 + t333 ^ 2 + t370 ^ 2) / 0.2e1 + (t474 * mrSges(4,2) + Ifges(4,4) * t465 + Ifges(4,5) * t479 + Ifges(4,1) * t464 / 0.2e1) * t464 + (t327 * mrSges(6,2) - t323 * mrSges(6,3) + Ifges(6,5) * t336 + Ifges(6,1) * t326 / 0.2e1) * t326 + (-pkin(16) * (-mrSges(9,1) * t388 + mrSges(9,2) * t389) + t361 * (Ifges(7,5) * t403 + Ifges(7,6) * t404) + qJD(2) * (Ifges(3,5) * t495 - Ifges(3,6) * t490) + (Ifges(2,3) / 0.2e1 + m(7) * pkin(15) ^ 2 / 0.2e1 + (m(9) / 0.2e1 + m(3) / 0.2e1) * pkin(16) ^ 2 + (-pkin(16) * mrSges(3,2) + Ifges(3,1) * t495 / 0.2e1) * t495 + (-pkin(15) * mrSges(7,1) + Ifges(7,2) * t404 / 0.2e1) * t404 + (-pkin(16) * mrSges(3,1) - Ifges(3,4) * t495 + Ifges(3,2) * t490 / 0.2e1) * t490 + (pkin(15) * mrSges(7,2) + Ifges(7,4) * t404 + Ifges(7,1) * t403 / 0.2e1) * t403) * qJD(1)) * qJD(1) + (-t385 * mrSges(10,1) + Ifges(10,4) * t375 + Ifges(10,2) * t374 / 0.2e1) * t374 + (Ifges(8,5) * t381 + Ifges(8,6) * t380 + Ifges(8,3) * t360 / 0.2e1) * t360 + m(6) * (t323 ^ 2 + t324 ^ 2 + t327 ^ 2) / 0.2e1 + (t323 * mrSges(6,1) - t324 * mrSges(6,2) + Ifges(6,3) * t336 / 0.2e1) * t336 + Ifges(4,3) * t479 ^ 2 / 0.2e1;
T = t1;
