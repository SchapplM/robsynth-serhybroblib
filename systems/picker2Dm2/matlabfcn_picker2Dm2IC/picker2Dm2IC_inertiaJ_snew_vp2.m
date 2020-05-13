% Calculate joint inertia matrix with Newton Euler for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2IC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_inertiaJ_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_inertiaJ_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_inertiaJ_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2IC_inertiaJ_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2IC_inertiaJ_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:21:06
% EndTime: 2020-05-11 09:21:08
% DurationCPUTime: 2.10s
% Computational Cost: add. (5763->326), mult. (4596->350), div. (313->10), fcn. (3988->60), ass. (0->193)
t500 = sin(qJ(9));
t507 = cos(qJ(9));
t451 = -mrSges(10,1) * t507 + mrSges(10,2) * t500;
t583 = pkin(6) * t451 + Ifges(10,3);
t498 = sin(qJ(10));
t499 = cos(qJ(10));
t445 = -mrSges(11,1) * t499 + mrSges(11,2) * t498;
t582 = pkin(4) * t445 + Ifges(11,3);
t501 = sin(qJ(8));
t508 = cos(qJ(8));
t579 = -t508 * mrSges(9,1) + t501 * mrSges(9,2);
t497 = qJ(1) + qJ(2);
t485 = sin(t497);
t486 = cos(t497);
t502 = sin(qJ(7));
t509 = cos(qJ(7));
t578 = (t485 * t502 + t486 * t509) * pkin(3);
t504 = sin(qJ(4));
t495 = 0.1e1 / t504;
t516 = 0.1e1 / pkin(4);
t556 = t495 * t516;
t555 = qJ(2) + qJ(4);
t484 = sin(t555);
t517 = 0.1e1 / pkin(3);
t570 = pkin(3) * t504;
t542 = t495 * (-pkin(1) * t484 - t570) * t517;
t577 = pkin(1) * m(9);
t576 = pkin(6) * m(10);
t575 = pkin(4) * m(11);
t506 = sin(qJ(2));
t574 = pkin(1) * t506;
t554 = qJ(3) - qJ(6);
t541 = -qJ(9) - t554;
t394 = 0.1e1 / ((-cos(qJ(4) - t554) + cos(qJ(4) + t554)) * pkin(6) + (cos(qJ(4) + t541) - cos(qJ(4) - t541)) * pkin(2));
t573 = pkin(2) * t394;
t505 = sin(qJ(3));
t572 = pkin(2) * t505;
t512 = cos(qJ(3));
t571 = pkin(2) * t512;
t511 = cos(qJ(4));
t569 = pkin(3) * t511;
t513 = cos(qJ(2));
t494 = t513 * pkin(1);
t568 = cos(qJ(1)) * pkin(1);
t515 = 0.1e1 / pkin(5);
t496 = qJ(1) + qJ(8);
t545 = pkin(8) + qJ(5);
t534 = t545 + t554;
t518 = t534 - t496;
t535 = t545 - t554;
t519 = t535 - t496;
t561 = 0.1e1 / (pkin(6) * (cos(t519) - cos(t518)) + (-cos(qJ(9) - t519) + cos(qJ(9) + t518)) * pkin(2)) * t515;
t560 = t394 * t517;
t489 = qJ(1) + t555;
t470 = sin(t489);
t471 = cos(t489);
t522 = t470 * t486 - t471 * t485;
t403 = 0.1e1 / t522;
t521 = t470 * t502 + t471 * t509;
t559 = t403 * t521;
t557 = (pkin(3) * t506 + pkin(4) * t484) * t517;
t474 = t494 + pkin(3);
t426 = t511 * t474 - t504 * t574;
t421 = pkin(4) + t426;
t429 = t474 * t504 + t511 * t574;
t382 = -t421 * t499 + t429 * t498;
t373 = m(11) * t382 + mrSges(11,1);
t383 = -t421 * t498 - t429 * t499;
t374 = m(11) * t383 - mrSges(11,2);
t528 = -t373 * t499 - t374 * t498;
t345 = m(5) * t426 + mrSges(5,1) + t528;
t346 = m(5) * t429 + t373 * t498 - t374 * t499 - mrSges(5,2);
t553 = t511 * t345 + t504 * t346;
t472 = pkin(4) + t569;
t423 = -t472 * t499 + t498 * t570;
t407 = m(11) * t423 + mrSges(11,1);
t424 = -t472 * t498 - t499 * t570;
t408 = m(11) * t424 - mrSges(11,2);
t362 = m(5) * t570 + t407 * t498 - t408 * t499 - mrSges(5,2);
t526 = -t407 * t499 - t408 * t498;
t363 = m(5) * t569 + mrSges(5,1) + t526;
t552 = t504 * t362 + t511 * t363;
t477 = -qJ(9) + t489;
t454 = t477 - t554;
t476 = qJ(9) + t489;
t455 = t476 + t554;
t551 = sin(t455) - sin(t454);
t550 = cos(t455) - cos(t454);
t549 = -cos(0.2e1 * t496) + cos(0.2e1 * t545);
t548 = sin(t477) - sin(t476);
t547 = cos(t477) - cos(t476);
t503 = sin(qJ(6));
t510 = cos(qJ(6));
t546 = t503 * mrSges(7,1) + t510 * mrSges(7,2);
t368 = -t550 * t502 + t551 * t509;
t544 = t368 * t573;
t452 = mrSges(10,1) * t500 + mrSges(10,2) * t507;
t399 = -t451 * t512 - t452 * t505;
t386 = pkin(2) * t399 + t583;
t543 = (-t549 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t515 / t549;
t446 = mrSges(11,1) * t498 + mrSges(11,2) * t499;
t398 = t445 * t511 + t446 * t504;
t379 = pkin(3) * t398 + t582;
t540 = -qJ(1) + t545;
t539 = -mrSges(7,1) * t510 + t503 * mrSges(7,2);
t475 = t494 + pkin(2);
t427 = -t475 * t512 + t505 * t574;
t536 = qJ(4) + t540;
t531 = -qJ(8) + t536;
t537 = -qJ(4) + t540;
t532 = -qJ(8) + t537;
t538 = (pkin(3) * (cos(t537) - cos(t536)) + (cos(qJ(2) - t532) - cos(qJ(2) + t531)) * pkin(5)) / (cos(t532) - cos(t531)) * t515 * t517;
t433 = pkin(1) * (t506 * t503 - t513 * t510);
t434 = pkin(1) * (-t513 * t503 - t506 * t510);
t393 = mrSges(7,1) * t433 - mrSges(7,2) * t434 + Ifges(7,3);
t422 = pkin(6) + t427;
t430 = -t475 * t505 - t512 * t574;
t384 = -t422 * t507 + t430 * t500;
t385 = -t422 * t500 - t430 * t507;
t357 = mrSges(10,1) * t384 - mrSges(10,2) * t385 + Ifges(10,3);
t473 = pkin(6) - t571;
t425 = -t473 * t507 - t500 * t572;
t428 = -t473 * t500 + t507 * t572;
t389 = mrSges(10,1) * t425 - mrSges(10,2) * t428 + Ifges(10,3);
t355 = mrSges(11,1) * t382 - mrSges(11,2) * t383 + Ifges(11,3);
t388 = mrSges(11,1) * t423 - mrSges(11,2) * t424 + Ifges(11,3);
t533 = pkin(1) * t579 + Ifges(9,3);
t375 = m(10) * t384 + mrSges(10,1);
t376 = m(10) * t385 - mrSges(10,2);
t527 = -t375 * t507 - t376 * t500;
t347 = m(4) * t427 + mrSges(4,1) + t527;
t348 = m(4) * t430 + t375 * t500 - t376 * t507 - mrSges(4,2);
t530 = -t347 * t512 - t348 * t505;
t409 = m(10) * t425 + mrSges(10,1);
t410 = m(10) * t428 - mrSges(10,2);
t365 = -m(4) * t572 + t409 * t500 - t410 * t507 - mrSges(4,2);
t525 = -t409 * t507 - t410 * t500;
t366 = -m(4) * t571 + mrSges(4,1) + t525;
t529 = -t365 * t505 - t366 * t512;
t457 = -t498 * t575 - mrSges(11,2);
t458 = -t499 * t575 + mrSges(11,1);
t524 = -t457 * t498 - t458 * t499;
t467 = -t500 * t576 - mrSges(10,2);
t468 = -t507 * t576 + mrSges(10,1);
t523 = -t467 * t500 - t468 * t507;
t392 = pkin(6) * t523 + Ifges(4,3) + t583;
t390 = pkin(4) * t524 + Ifges(5,3) + t582;
t405 = -t467 * t507 + t468 * t500 - mrSges(4,2);
t406 = mrSges(4,1) + t523;
t369 = -t405 * t505 - t406 * t512;
t356 = pkin(2) * t369 + t392;
t401 = -t457 * t499 + t458 * t498 - mrSges(5,2);
t402 = mrSges(5,1) + t524;
t364 = t401 * t504 + t402 * t511;
t354 = pkin(3) * t364 + t390;
t339 = pkin(6) * t527 + mrSges(4,1) * t427 - mrSges(4,2) * t430 + Ifges(4,3) + t357;
t337 = pkin(4) * t528 + mrSges(5,1) * t426 - mrSges(5,2) * t429 + Ifges(5,3) + t355;
t352 = pkin(6) * t525 - mrSges(4,1) * t571 + mrSges(4,2) * t572 + Ifges(4,3) + t389;
t350 = pkin(4) * t526 + mrSges(5,1) * t569 - mrSges(5,2) * t570 + Ifges(5,3) + t388;
t326 = pkin(2) * t529 + pkin(3) * t552 + Ifges(3,3) + Ifges(7,3) + t350 + t352;
t323 = pkin(2) * t530 + pkin(3) * t553 + mrSges(3,1) * t494 - mrSges(3,2) * t574 + Ifges(3,3) + t337 + t339 + t393;
t493 = sin(qJ(1)) * pkin(1);
t488 = -qJ(9) + t545;
t487 = qJ(9) + t545;
t466 = qJ(9) + t534;
t465 = -qJ(9) + t535;
t459 = -sin(qJ(8) - t540);
t443 = -pkin(5) * sin(t496) + t493;
t442 = t568 - pkin(5) * cos(t496);
t420 = m(7) * t434 - mrSges(7,2);
t419 = m(7) * t433 + mrSges(7,1);
t414 = pkin(3) * t486 + pkin(4) * t471 + t568;
t413 = pkin(3) * t485 + pkin(4) * t470 + t493;
t391 = t521 * pkin(4) + t578;
t381 = t547 * t502 - t548 * t509;
t380 = t522 * pkin(4) - t578;
t358 = t381 * Ifges(7,3) * t573;
t341 = -Ifges(7,3) * t559 + t358;
t340 = t379 * t542 + (t506 * (-t445 * t504 + t446 * t511) + t513 * t398 + (-Ifges(11,3) * t506 + t557 * t582) * t556) * pkin(1) + t379;
t338 = t583 * t544 - (Ifges(10,3) + t386) * t559;
t336 = (-t379 * t521 + (Ifges(11,3) * t380 + t391 * t582) * t516) * t403;
t333 = ((t548 * t413 + t547 * t414) * t560 + (-(sin(t488) - sin(t487)) * t443 - (cos(t488) - cos(t487)) * t442) * t561) * pkin(2);
t332 = t333 * Ifges(7,3);
t331 = ((-t551 * t413 - t550 * t414) * t560 + ((sin(t466) - sin(t465)) * t443 + (cos(t466) - cos(t465)) * t442) * t561) * pkin(2);
t330 = t392 * t544 - (t356 + t583) * t559;
t329 = (-t354 * t521 + (t380 * t582 + t390 * t391) * t516) * t403;
t328 = pkin(1) * (t506 * t546 + t513 * t539) + Ifges(7,3) + Ifges(7,3) * t542 + t332;
t327 = t354 * t542 + (t506 * (t401 * t511 - t402 * t504) + t513 * t364 + (t390 * t557 - t506 * t582) * t556) * pkin(1) + t354;
t325 = t386 * t542 + t331 * t583 + (t506 * (t451 * t505 - t452 * t512) + t513 * t399 - Ifges(10,3) * t538) * pkin(1) + t386;
t324 = t356 * t542 + t331 * t392 + (t506 * (-t405 * t512 + t406 * t505) + t513 * t369 - t583 * t538) * pkin(1) + t356;
t322 = t352 * t544 + t358 + ((t350 * t391 + t380 * t388) * t516 - (t326 + t389) * t521) * t403;
t321 = t326 + t332 + (t506 * (t362 * t511 - t363 * t504 - t365 * t512 + t366 * t505 - mrSges(3,2) + t546) + t513 * (mrSges(3,1) + t529 + t539 + t552) - t389 * t538 + (t350 * t557 - t388 * t506) * t556) * pkin(1) + t331 * t352 + t326 * t542;
t1 = [(t321 + t323) * t542 + t323 + (t324 + t339) * t331 + t501 ^ 2 / t459 ^ 2 * Ifges(6,3) + (t328 + t393) * t333 + Ifges(2,3) + t533 + (Ifges(9,3) * t543 + Ifges(9,3) + t533) * t543 + (((-t340 - t355) * t506 + (t327 + t337) * t557) * t556 + t506 * (m(3) * t574 - t345 * t504 + t346 * t511 + t347 * t505 - t348 * t512 + t419 * t503 - t420 * t510 - mrSges(3,2)) + t513 * (m(3) * t494 - t419 * t510 - t420 * t503 + mrSges(3,1) + t530 + t553) + t579 * t543 - t501 * (-t501 * t577 - mrSges(9,2)) - t508 * (-t508 * t577 + mrSges(9,1)) + (-t325 - t357) * t538) * pkin(1), (t324 * t368 + t328 * t381) * t573 + ((t327 * t391 + t340 * t380) * t516 - (t321 + t325) * t521) * t403; t322 * t542 + t330 * t331 + t341 * t333 + (t339 * t368 + t381 * t393) * t573 + (-t338 * t538 + (t329 * t557 - t336 * t506) * t556) * pkin(1) + ((t337 * t391 + t355 * t380) * t516 - (t323 + t357) * t521) * t403, Ifges(8,3) + (t330 * t368 + t341 * t381) * t573 + ((t329 * t391 + t336 * t380) * t516 - (t322 + t338) * t521) * t403;];
Mq = t1;
