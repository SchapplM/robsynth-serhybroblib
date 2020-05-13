% Calculate kinetic energy for
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m1TE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(23,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1TE_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:11:07
% EndTime: 2020-04-12 20:11:52
% DurationCPUTime: 44.24s
% Computational Cost: add. (992329->312), mult. (1539069->596), div. (63400->34), fcn. (959844->26), ass. (0->248)
t496 = pkin(7) ^ 2;
t504 = pkin(1) ^ 2;
t476 = sin(qJ(2));
t477 = sin(pkin(19));
t481 = cos(qJ(2));
t482 = cos(pkin(19));
t456 = t476 * t482 - t477 * t481;
t577 = pkin(7) * t456;
t608 = -2 * pkin(1);
t547 = t577 * t608 + t504;
t442 = t496 + t547;
t544 = pkin(3) ^ 2 - pkin(8) ^ 2;
t435 = t442 - t544;
t448 = pkin(1) * t456 - pkin(7);
t457 = t476 * t477 + t481 * t482;
t605 = pkin(7) + pkin(8);
t606 = pkin(7) - pkin(8);
t430 = (pkin(3) + t605) * (-pkin(3) + t606) + t547;
t431 = (-pkin(3) + t605) * (pkin(3) + t606) + t547;
t507 = sqrt(-t431 * t430);
t559 = t457 * t507;
t416 = -pkin(1) * t559 - t435 * t448;
t419 = pkin(1) * t435 * t457 - t448 * t507;
t478 = sin(pkin(18));
t440 = 0.1e1 / t442;
t495 = 0.1e1 / pkin(8);
t562 = t440 * t495;
t483 = cos(pkin(18));
t588 = t483 / 0.2e1;
t411 = (t419 * t588 + t416 * t478 / 0.2e1) * t562;
t613 = t411 ^ 2;
t612 = Ifges(3,5) * t481;
t611 = pkin(13) ^ 2;
t469 = sin(pkin(20));
t473 = cos(pkin(20));
t475 = sin(qJ(3));
t480 = cos(qJ(3));
t518 = t469 * t480 + t473 * t475;
t579 = pkin(6) * t518;
t539 = pkin(1) * t579;
t444 = 0.2e1 * t539;
t502 = pkin(2) ^ 2;
t497 = pkin(6) ^ 2;
t545 = t497 + t504;
t533 = t545 - t611;
t434 = t444 + t502 + t533;
t443 = -pkin(1) - t579;
t548 = t444 + t497;
t607 = (-pkin(2) - pkin(13));
t428 = ((pkin(1) - t607) * (pkin(1) + t607)) + t548;
t602 = (pkin(13) - pkin(2));
t429 = ((pkin(1) - t602) * (pkin(1) + t602)) + t548;
t564 = t429 * t428;
t506 = sqrt(-t564);
t455 = t469 * t475 - t473 * t480;
t578 = pkin(6) * t455;
t414 = -t434 * t443 - t506 * t578;
t413 = 0.1e1 / t414 ^ 2;
t415 = t434 * t578 - t443 * t506;
t439 = t444 + t545;
t437 = 0.1e1 / t439;
t438 = 0.1e1 / t439 ^ 2;
t445 = t518 * qJD(3);
t446 = t455 * qJD(3);
t503 = 0.1e1 / pkin(2);
t586 = pkin(1) * t446;
t542 = pkin(6) * t586;
t566 = 0.2e1 * (t428 + t429) * t542 / t506;
t528 = -t566 / 0.2e1;
t580 = pkin(6) * t415;
t594 = t437 / 0.2e1;
t596 = t414 / 0.2e1;
t368 = qJD(2) + (((t443 * t528 + (t434 * t445 - t446 * t506) * pkin(6)) * t594 + (-t437 * t455 * t497 + t438 * t580) * t586) / t596 - 0.2e1 * ((-t446 * t434 - t445 * t506 + t455 * t528) * t594 + (t414 * t438 + t437 * t443) * t586) * t413 * t580) * pkin(2) / (t413 * t415 ^ 2 + 0.1e1) * t439 * t503;
t488 = 0.1e1 / pkin(13);
t530 = -t368 * t488 / 0.2e1;
t610 = mrSges(10,3) * t530;
t485 = qJD(1) ^ 2;
t604 = -pkin(9) - pkin(11);
t603 = pkin(11) - pkin(9);
t499 = pkin(4) ^ 2;
t498 = pkin(5) ^ 2;
t436 = t442 + t544;
t447 = pkin(1) - t577;
t418 = pkin(7) * t436 * t457 + t447 * t507;
t552 = t480 * t418;
t417 = -pkin(7) * t559 + t436 * t447;
t555 = t475 * t417;
t515 = t555 / 0.2e1 + t552 / 0.2e1;
t501 = 0.1e1 / pkin(3);
t561 = t440 * t501;
t408 = t515 * t561;
t553 = t480 * t417;
t554 = t475 * t418;
t514 = -t553 / 0.2e1 + t554 / 0.2e1;
t409 = t514 * t561;
t464 = pkin(23) + pkin(22);
t462 = sin(t464);
t463 = cos(t464);
t390 = -t408 * t463 + t409 * t462;
t582 = pkin(5) * t390;
t549 = -0.2e1 * pkin(4) * t582 + t498;
t386 = t499 + t549;
t546 = pkin(9) ^ 2 - pkin(11) ^ 2;
t382 = t386 + t546;
t387 = -pkin(4) + t582;
t391 = t408 * t462 + t409 * t463;
t380 = (pkin(4) - t604) * (pkin(4) + t604) + t549;
t381 = (pkin(4) - t603) * (pkin(4) + t603) + t549;
t505 = sqrt(-t381 * t380);
t568 = t391 * t505;
t362 = -pkin(5) * t568 - t382 * t387;
t601 = t362 / 0.2e1;
t384 = 0.1e1 / t386;
t600 = t384 / 0.2e1;
t453 = t457 * qJD(2);
t452 = t456 * qJD(2);
t541 = pkin(1) * pkin(7) * t453;
t565 = 0.2e1 * (t430 + t431) * t541 / t507;
t527 = -t565 / 0.2e1;
t512 = t452 * t507 + t457 * t527;
t397 = ((t447 * t608 - t436) * t453 + t512) * pkin(7);
t599 = -t397 / 0.2e1;
t538 = -0.2e1 * t453 * t457;
t560 = t453 * t507;
t398 = t447 * t565 / 0.2e1 + t496 * pkin(1) * t538 + (-t436 * t452 - t560) * pkin(7);
t598 = t398 / 0.2e1;
t532 = qJD(1) * t437 * t503;
t590 = -t476 / 0.2e1;
t401 = (t415 * t590 + t481 * t596) * t532;
t597 = t401 / 0.2e1;
t433 = t502 - t533 - 0.2e1 * t539;
t595 = -t433 / 0.2e1;
t466 = sin(pkin(23));
t593 = t466 / 0.2e1;
t468 = sin(pkin(21));
t592 = t468 / 0.2e1;
t471 = cos(pkin(22));
t591 = -t471 / 0.2e1;
t589 = -t480 / 0.2e1;
t587 = t485 / 0.2e1;
t470 = cos(pkin(23));
t557 = t470 * t417;
t558 = t466 * t418;
t405 = (-t557 / 0.2e1 + t558 / 0.2e1) * t561;
t508 = t405 ^ 2;
t403 = 0.1e1 / t508;
t556 = t470 * t418;
t567 = t417 * t466;
t406 = (t556 / 0.2e1 + t567 / 0.2e1) * t561;
t404 = t406 ^ 2;
t524 = 0.1e1 / t442 ^ 2 * t541;
t355 = qJD(2) + (-((t398 * t593 + t470 * t599) * t440 + (-t557 + t558) * t524) * t406 * t403 + ((t397 * t593 + t470 * t598) * t440 + (t556 + t567) * t524) / t405) / (t403 * t404 + 0.1e1) * t501;
t585 = pkin(4) * t355;
t370 = ((-t553 + t554) * t524 + (qJD(3) * t515 + t397 * t589 + t475 * t598) * t440) * t501;
t371 = ((-t552 - t555) * t524 + (qJD(3) * t514 + t398 * t589 + t475 * t599) * t440) * t501;
t365 = t370 * t462 + t371 * t463;
t584 = pkin(4) * t365;
t581 = pkin(5) * t391;
t363 = t382 * t581 - t387 * t505;
t583 = pkin(5) * t363;
t575 = pkin(1) * qJD(2);
t574 = qJD(1) * pkin(16);
t540 = pkin(5) * t584;
t573 = 0.2e1 * (t380 + t381) * t540 / t505;
t572 = t365 * t505;
t472 = cos(pkin(21));
t571 = t384 * t472;
t491 = 0.1e1 / pkin(11);
t570 = t384 * t491;
t493 = 0.1e1 / pkin(9);
t569 = t384 * t493;
t563 = t440 * t478;
t551 = qJD(2) ^ 2 * t504;
t550 = t488 * t503;
t543 = qJD(1) * t476;
t465 = qJD(2) + qJD(3);
t537 = t405 * t575;
t536 = t406 * t575;
t535 = t475 * t575;
t534 = t480 * t575;
t531 = -t573 / 0.2e1;
t529 = t384 * t592;
t526 = t440 * t588;
t385 = 0.1e1 / t386 ^ 2;
t523 = t385 * t540;
t460 = pkin(1) * t543 - t574;
t520 = t478 * t524;
t519 = t483 * t524;
t383 = t386 - t546;
t388 = -pkin(4) * t390 + pkin(5);
t516 = pkin(4) * t383 * t391 + t388 * t505;
t517 = -pkin(4) * t568 + t383 * t388;
t349 = (t517 * t592 + t516 * t472 / 0.2e1) * t570;
t350 = (-t517 * t472 / 0.2e1 + t516 * t592) * t570;
t450 = (t475 * t481 + t476 * t480) * qJD(1);
t451 = (-t475 * t476 + t480 * t481) * qJD(1);
t340 = -t349 * t450 + t350 * t451;
t458 = pkin(5) * t465 + t535;
t343 = t349 * t458 - t350 * t534;
t427 = -pkin(5) * t451 + t460;
t364 = t370 * t463 - t371 * t462;
t513 = -t364 * t505 + t391 * t531;
t342 = t349 * t534 + t350 * t458;
t511 = t517 * t523;
t510 = t516 * t523;
t509 = t433 ^ 2;
t479 = cos(qJ(4));
t474 = sin(qJ(4));
t467 = sin(pkin(22));
t459 = t460 ^ 2;
t432 = 0.1e1 / t509;
t410 = (t416 * t588 - t478 * t419 / 0.2e1) * t562;
t407 = 0.1e1 / t410 ^ 2;
t402 = (-t415 * t481 / 0.2e1 + t414 * t590) * t532;
t400 = -pkin(2) * t402 - t574;
t399 = t448 * t527 + t504 * pkin(7) * t538 + (-t435 * t452 - t560) * pkin(1);
t396 = ((0.2e1 * pkin(7) * t448 - t435) * t453 + t512) * pkin(1);
t393 = (-t405 * t476 - t406 * t481) * qJD(1);
t392 = (t405 * t481 - t406 * t476) * qJD(1);
t379 = (t402 * t595 + t506 * t597) * t550;
t378 = (t401 * t595 - t506 * t402 / 0.2e1) * t550;
t374 = (-t392 * t467 - t393 * t471) * pkin(4) + t460;
t366 = (0.1e1 / t433 * t566 / 0.2e1 - 0.2e1 * t506 * t432 * t542) / (-t432 * t564 + 0.1e1) + t368;
t361 = 0.1e1 / t362 ^ 2;
t356 = ((t399 * t526 + t419 * t519 + t396 * t563 / 0.2e1 + t416 * t520) / t410 - (t396 * t526 + t416 * t519 - t399 * t563 / 0.2e1 - t419 * t520) * t411 * t407) / (t407 * t613 + 0.1e1) * t495;
t354 = t467 * t585 + t536;
t353 = t471 * t585 + t537;
t348 = (t362 * t591 - t467 * t363 / 0.2e1) * t569;
t347 = (t363 * t591 + t467 * t601) * t569;
t346 = 0.1e1 / t350 ^ 2;
t341 = t349 * t451 + t350 * t450;
t339 = qJD(4) - t340;
t338 = -t347 * t392 + t348 * t393;
t337 = t347 * t393 + t348 * t392;
t336 = t388 * t573 / 0.2e1 - 0.2e1 * t499 * t365 * t581 + (t364 * t383 - t572) * pkin(4);
t335 = ((-0.2e1 * pkin(5) * t388 - t383) * t365 + t513) * pkin(4);
t334 = -t347 * t354 + t348 * t353;
t333 = t347 * t353 + t348 * t354;
t332 = -pkin(10) * t340 - pkin(12) * t341 + t427;
t331 = t355 + (((t387 * t531 + (t364 * t382 - t572) * pkin(5)) * t600 + (-t384 * t391 * t498 + t385 * t583) * t584) / t601 - 0.2e1 * ((-t365 * t382 + t513) * t600 + (t362 * t385 + t384 * t387) * t584) * t361 * t583) * pkin(9) / (t361 * t363 ^ 2 + 0.1e1) * t386 * t493;
t330 = ((t335 * t529 + t468 * t511 + t336 * t571 / 0.2e1 + t472 * t510) / t350 - (-t335 * t571 / 0.2e1 - t472 * t511 + t336 * t529 + t468 * t510) * t349 * t346) / (t346 * t349 ^ 2 + 0.1e1) * t491 + t465;
t329 = pkin(12) * t330 + t343;
t328 = -pkin(10) * t330 - t342;
t327 = t330 * t474 + t341 * t479;
t326 = t330 * t479 - t341 * t474;
t325 = t329 * t479 + t332 * t474;
t324 = -t329 * t474 + t332 * t479;
t1 = (qJD(2) * t612 + (t481 * (Ifges(3,1) * t481 - Ifges(3,4) * t476) + Ifges(7,1) * t613 + (0.2e1 * Ifges(7,4) * t411 + Ifges(7,2) * t410) * t410) * qJD(1)) * qJD(1) / 0.2e1 + qJD(2) * (Ifges(3,3) * qJD(2) + (-Ifges(3,6) * t476 + t612) * qJD(1)) / 0.2e1 + ((Ifges(7,5) * t411 + Ifges(7,6) * t410) * qJD(1) + Ifges(7,3) * t356 / 0.2e1) * t356 + t368 * (Ifges(9,5) * t401 + Ifges(9,3) * t368) / 0.2e1 + (mrSges(8,1) * t537 - mrSges(8,2) * t536 + Ifges(8,5) * t392 + Ifges(8,6) * t393 + Ifges(8,3) * t355 / 0.2e1) * t355 + (t485 * (-mrSges(7,1) * t410 + mrSges(7,2) * t411) + m(7) * t587 * pkin(15)) * pkin(15) + (t328 * mrSges(6,2) - t324 * mrSges(6,3) + Ifges(6,5) * t339 + Ifges(6,1) * t327 / 0.2e1) * t327 + (t460 * mrSges(4,2) - mrSges(4,3) * t535 + Ifges(4,4) * t451 + Ifges(4,5) * t465 + Ifges(4,1) * t450 / 0.2e1) * t450 + (-t433 * t610 + t400 * mrSges(10,2) + Ifges(10,4) * t379 + Ifges(10,5) * t366 + Ifges(10,1) * t378 / 0.2e1) * t378 + (t506 * t610 - t400 * mrSges(10,1) + Ifges(10,2) * t379 / 0.2e1) * t379 + (Ifges(10,6) * t379 + Ifges(10,3) * t366 / 0.2e1 + (mrSges(10,1) * t433 - mrSges(10,2) * t506) * t530) * t366 + (t427 * mrSges(5,2) - t342 * mrSges(5,3) + Ifges(5,1) * t341 / 0.2e1) * t341 + (-t460 * mrSges(8,1) + mrSges(8,3) * t536 + Ifges(8,2) * t393 / 0.2e1) * t393 + (t460 * mrSges(8,2) - mrSges(8,3) * t537 + Ifges(8,4) * t393 + Ifges(8,1) * t392 / 0.2e1) * t392 + (-t460 * mrSges(4,1) - mrSges(4,3) * t534 + Ifges(4,6) * t465 + Ifges(4,2) * t451 / 0.2e1) * t451 + m(10) * (t400 ^ 2 + (-t564 / 0.4e1 + t509 / 0.4e1) / t611 * t368 ^ 2) / 0.2e1 + Ifges(2,3) * t587 + (Ifges(9,1) * t401 + Ifges(9,4) * t402 + Ifges(9,5) * t368) * t597 + (t334 * mrSges(11,1) - t333 * mrSges(11,2) + Ifges(11,5) * t337 + Ifges(11,6) * t338 + Ifges(11,3) * t331 / 0.2e1) * t331 + m(4) * (t459 + (t475 ^ 2 + t480 ^ 2) * t551) / 0.2e1 + m(8) * (t459 + (t404 + t508) * t551) / 0.2e1 + (-t485 * (mrSges(3,1) * t476 + mrSges(3,2) * t481) + (m(9) + m(3)) * t587 * pkin(16)) * pkin(16) - (-mrSges(9,1) * t402 + mrSges(9,2) * t401) * t574 + (t374 * mrSges(11,2) - t334 * mrSges(11,3) + Ifges(11,4) * t338 + Ifges(11,1) * t337 / 0.2e1) * t337 + (t324 * mrSges(6,1) - t325 * mrSges(6,2) + Ifges(6,3) * t339 / 0.2e1) * t339 - (Ifges(3,6) * qJD(2) + (Ifges(3,4) * t481 - Ifges(3,2) * t476) * qJD(1)) * t543 / 0.2e1 + t402 * (Ifges(9,4) * t401 + Ifges(9,2) * t402) / 0.2e1 + (t342 * mrSges(5,1) - t343 * mrSges(5,2) + Ifges(5,5) * t341 + Ifges(5,6) * t340 + Ifges(5,3) * t330 / 0.2e1) * t330 + (-t328 * mrSges(6,1) + t325 * mrSges(6,3) + Ifges(6,4) * t327 + Ifges(6,6) * t339 + Ifges(6,2) * t326 / 0.2e1) * t326 + (-t427 * mrSges(5,1) + t343 * mrSges(5,3) + Ifges(5,4) * t341 + Ifges(5,2) * t340 / 0.2e1) * t340 + t368 * t402 * Ifges(9,6) + (-t374 * mrSges(11,1) + t333 * mrSges(11,3) + Ifges(11,2) * t338 / 0.2e1) * t338 + m(5) * (t342 ^ 2 + t343 ^ 2 + t427 ^ 2) / 0.2e1 + (mrSges(4,1) * t535 + mrSges(4,2) * t534 + Ifges(4,3) * t465 / 0.2e1) * t465 + m(11) * (t333 ^ 2 + t334 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t324 ^ 2 + t325 ^ 2 + t328 ^ 2) / 0.2e1;
T = t1;
