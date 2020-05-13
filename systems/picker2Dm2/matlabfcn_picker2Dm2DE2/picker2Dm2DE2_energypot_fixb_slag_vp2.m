% Calculate potential energy for
% picker2Dm2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm2DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm2DE2_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm2DE2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2DE2_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2DE2_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 19:15:51
% EndTime: 2020-05-09 19:15:59
% DurationCPUTime: 5.34s
% Computational Cost: add. (31097->453), mult. (92600->537), div. (594->5), fcn. (15660->31), ass. (0->234)
t590 = 4 * pkin(1);
t582 = m(10) + m(4);
t449 = pkin(4) ^ 2;
t397 = -t449 / 0.4e1;
t454 = pkin(3) ^ 2;
t589 = t397 + t454 / 0.2e1;
t588 = 2 * pkin(7);
t421 = cos(qJ(1));
t387 = t421 ^ 2;
t587 = -0.2e1 * t387;
t433 = 0.2e1 * t454;
t586 = 0.4e1 * t454;
t459 = pkin(1) ^ 2;
t457 = t459 ^ 2;
t585 = 4 * t457;
t584 = 2 * t459;
t437 = 6 * t459;
t461 = pkin(7) ^ 2;
t445 = 2 * t461;
t466 = t454 ^ 2;
t432 = 0.5e1 * t466;
t583 = pkin(5) * m(6);
t414 = sin(pkin(8));
t415 = cos(pkin(8));
t418 = sin(qJ(1));
t342 = t414 * t418 + t415 * t421;
t581 = t342 / 0.2e1;
t580 = m(11) + m(5);
t375 = pkin(1) * t421;
t362 = t375 + pkin(7);
t417 = sin(qJ(2));
t373 = pkin(3) * t417;
t420 = cos(qJ(2));
t374 = pkin(3) * t420;
t579 = t418 * pkin(1);
t578 = mrSges(7,2) + mrSges(4,2);
t363 = t373 * t588;
t384 = t417 ^ 2;
t568 = t384 * t454;
t529 = 0.2e1 * t568;
t544 = -t454 + t461;
t343 = t363 + t529 + t544;
t567 = t418 * t420;
t525 = pkin(3) * t567;
t495 = pkin(1) * t525;
t355 = -0.2e1 * t495;
t447 = 0.2e1 * pkin(3);
t379 = t459 + t461;
t505 = -t449 + t379;
t488 = t363 + t505;
t562 = t459 * t387;
t527 = -0.4e1 * t562;
t533 = -0.4e1 * t373;
t542 = t459 - t461;
t545 = t449 - t461;
t361 = t373 + pkin(7);
t574 = t361 * t421;
t326 = sqrt(t343 * t527 + 0.4e1 * t542 * t568 + pkin(7) * t505 * t533 - t457 + (-0.2e1 * t454 + t545) * t584 - (t461 - (t447 + pkin(4)) * pkin(4)) * (t461 + (t447 - pkin(4)) * pkin(4)) + (-(t355 + t488) * t574 + t488 * t525) * t590);
t448 = t449 ^ 2;
t460 = t461 ^ 2;
t543 = t457 + t460;
t548 = t445 - t449;
t561 = t461 * t449;
t478 = t548 * t459 + t448 / 0.6e1 + t543 - t561;
t340 = -t466 / 0.6e1 + t478;
t408 = -t454 / 0.3e1;
t369 = t408 + t461;
t344 = t369 * t355;
t350 = t373 + t362;
t378 = -0.3e1 * t454 + t461;
t386 = t421 * t387;
t462 = pkin(1) * t459;
t559 = t462 * t386;
t531 = pkin(7) * t559;
t500 = 0.8e1 * t531;
t353 = t378 * t500;
t377 = -t449 - t454;
t444 = 3 * t461;
t365 = t444 + t377;
t573 = t365 * t459;
t354 = 0.10e2 * t573;
t404 = 0.4e1 / 0.3e1 * t454;
t398 = -t449 / 0.3e1;
t510 = t398 + t379;
t356 = t404 + t510;
t399 = -t449 / 0.2e1;
t504 = t454 + t379;
t358 = t399 + t504;
t359 = -t449 + t504;
t536 = 0.2e1 * t375;
t364 = pkin(7) * t536;
t443 = 4 * t461;
t367 = (t443 + t449) * t459;
t370 = -t459 / 0.3e1 + t461;
t371 = 0.10e2 / 0.3e1 * t459;
t372 = t379 ^ 2;
t376 = -0.30e2 * t449 + (60 * t461);
t381 = -3 * t459 + t461;
t396 = -t449 / 0.6e1;
t405 = 0.2e1 / 0.3e1 * t454;
t410 = 0.4e1 / 0.3e1 * t459;
t412 = t459 / 0.2e1;
t423 = 15 * t457;
t424 = 15 * t459;
t425 = 10 * t459;
t430 = -0.2e1 * t449;
t431 = -0.5e1 * t449;
t434 = 7 * t457;
t435 = 5 * t457;
t436 = 7 * t459;
t441 = 3 * t460;
t442 = 8 * t461;
t465 = pkin(3) * t454;
t451 = t465 ^ 2;
t470 = pkin(7) * t461;
t477 = 0.5e1 / 0.6e1 * t466 + t478;
t479 = t461 - t495;
t556 = t448 / 0.2e1 - t466 / 0.2e1;
t487 = -0.3e1 * t561 + t441 + t556;
t491 = -0.6e1 * t495;
t401 = -0.3e1 / 0.2e1 * t449;
t555 = t401 + t444;
t558 = t379 * ((t401 + t445) * t459 - 0.3e1 / 0.2e1 * t561 + t543 + t556) + t451;
t480 = ((t371 + t548) * t454 + t477) * t491 + (t423 + (-0.9e1 * t449 + (18 * t461)) * t459 + t487) * t454 + (t424 + t555) * t466 + t558;
t490 = -0.4e1 * t495;
t481 = t358 * t490;
t438 = 3 * t459;
t549 = t438 + t461;
t507 = t454 + t549;
t482 = -(0.3e1 * t454 + t379) * t579 + t507 * t374;
t400 = -0.2e1 / 0.3e1 * t449;
t409 = -0.2e1 / 0.3e1 * t454;
t508 = t400 + t379;
t550 = t425 + t445;
t554 = t409 + t461;
t483 = -(t432 + ((5 * t459) + t365) * t433 + (t409 + t508) * t379) * t579 + (t466 + (t400 + t409 + t550) * t454 + t435 + 0.2e1 * t573 + t461 * (t400 + t554)) * t374;
t547 = t448 - t466;
t486 = -0.6e1 * t561 + (6 * t460) + t547;
t509 = t400 + t405 + t445;
t557 = (t405 + t508) * t379 + t466;
t489 = t356 * t490 + (t437 + t509) * t454 + t557;
t522 = t462 * t374;
t493 = t386 * t522;
t563 = t457 * t387 ^ 2;
t494 = t563 * t374;
t520 = 0.16e2 * t559;
t498 = pkin(7) * t520;
t499 = 0.20e2 / 0.3e1 * t459;
t546 = -t449 + t454;
t506 = t444 + t546;
t406 = t454 / 0.3e1;
t511 = t396 + t406 + t461;
t512 = t449 / 0.3e1 + t406 + t445;
t513 = 0.2e1 / 0.3e1 * t449 + t405 + t443;
t514 = 0.4e1 / 0.3e1 * t449 + t404 - (2 * t461);
t565 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t515 = t418 * t565;
t534 = 0.6e1 * t375;
t517 = pkin(7) * t534;
t535 = 0.4e1 * t375;
t518 = pkin(7) * t535;
t519 = -t579 / 0.2e1;
t521 = 0.12e2 * t562;
t523 = t459 * t374;
t566 = t420 * t421;
t524 = pkin(3) * t566;
t526 = 0.4e1 * t562;
t528 = 0.8e1 * t563;
t569 = t417 * t384 * t465;
t530 = -0.8e1 * t569;
t532 = 0.2e1 * t579;
t537 = pkin(7) * t375;
t539 = 4 * pkin(7);
t540 = t460 + t466;
t541 = t460 - t457;
t551 = 0.4e1 / 0.7e1 * t461 - t449 / 0.7e1;
t552 = t412 + t461;
t553 = t459 / 0.3e1 + t461;
t560 = t461 * t459;
t564 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t570 = t384 ^ 2 * t466;
t571 = t377 * t461;
t572 = t372 * (-t454 + t505);
t575 = (-t418 * t462 + t523) * t387;
t577 = ((-0.24e2 * (0.4e1 / 0.3e1 * t562 + t364 + t370) * t570 * t579 - 0.12e2 * (-0.8e1 / 0.3e1 * t494 + ((t410 + t511) * t374 - (0.7e1 / 0.6e1 * t454 + t396 + t552) * t579) * t526 + (-t454 * t542 - 0.5e1 / 0.3e1 * t457 + t512 * t459 + t461 * (t398 + t369)) * t374 + (-t466 + (-t499 + t513) * t454 - (3 * t457) + t514 * t459 + t460) * t519 + (-t418 * t457 * t386 + ((t459 + t511) * t374 + (t433 - t542) * t519) * t375) * t539) * t568 + 0.24e2 * t369 * t494 + ((t461 + 0.5e1 / 0.2e1 * t454 + 0.3e1 / 0.2e1 * t459 + t399) * t374 + t378 * t579 / 0.2e1) * t498 - 0.6e1 * ((-0.3e1 * t466 + (-t499 + t514) * t454 + t513 * t459 + t541) * t374 - 0.2e1 * (-0.5e1 / 0.3e1 * t466 + (-t459 + t512) * t454 + t461 * (t408 + t510)) * t579) * t562 - 0.6e1 * t483 * t537 - (t451 + ((21 * t459) + t365) * t466 + (t354 + t441 + (35 * t457) + 0.2e1 * t571) * t454 + (t434 + (t431 + t442 - 0.5e1 * t454) * t459 + t461 * (-t449 + t544)) * t379) * t374 + (0.7e1 * t451 + (t436 + t365) * t432 + (t354 + (21 * t457) + (9 * t460) + 0.6e1 * t571) * t454 + t572) * t579) * t326 + (0.16e2 * (t528 + t498 + (-8 * t457 + 12 * t560) * t387 + (-12 * pkin(7) * t462 + t470 * t590) * t421 - (6 * t560) + t543) * t570 + 0.24e2 * (t554 * t528 + 0.14e2 * (-0.32e2 / 0.21e2 * (t461 + t454 / 0.4e1 + t459 / 0.4e1 - t449 / 0.8e1) * t495 + 0.5e1 / 0.42e2 * t466 + (0.16e2 / 0.21e2 * t459 + t551) * t454 + t457 / 0.7e1 + t551 * t459 + t460 - 0.3e1 / 0.7e1 * t561 + t448 / 0.42e2) * t562 + t370 * t481 - t542 * t466 + (t367 - 0.10e2 / 0.3e1 * t457 + (2 * t460) - t561) * t454 + t340 * t564 + ((-0.2e1 / 0.3e1 * t495 + t397 + t552) * t520 + (-0.8e1 / 0.3e1 * (t553 + t589) * t495 + 0.5e1 / 0.18e2 * t466 + (0.4e1 / 0.3e1 * t461 + t410 + t398) * t454 + t460 + 0.2e1 / 0.3e1 * t560 - 0.2e1 / 0.3e1 * t561 - t457 / 0.3e1 + t448 / 0.18e2) * t534) * pkin(7)) * t568 + 0.16e2 * (-0.6e1 * t461 * t454 + t540) * t563 + 0.32e2 * (t355 * t565 + t358 * t378) * t531 + 0.24e2 * (t369 * t481 - t451 + (-t371 + t545) * t466 + (t367 + t466 / 0.6e1 - t448 / 0.6e1 + t541) * t454 + t340 * t461) * t562 + 0.8e1 * t480 * t537 - 0.8e1 * ((t436 + t555) * t466 + (t434 + (t431 + (10 * t461)) * t459 + t487) * t454 + t558) * t495 + t466 ^ 2 + (t430 + t443 + (28 * t459)) * t451 + (t376 * t459 + (70 * t457) + t486) * t466 + (t376 * t457 + t486 * t437 + t547 * t445 - 0.6e1 * t460 * t449 + (28 * t462 ^ 2) + (4 * t470 ^ 2)) * t454 + t359 * t572) * t350 + (((0.4e1 * t575 + (t374 + t532) * t364 + t381 * t374 + (t399 + t507) * t532) * t530 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t454 + t412 + t396) * t374 + pkin(1) * t515) * t527 + (-0.8e1 * t493 + ((t398 + t405 + t549) * t374 - (0.8e1 / 0.3e1 * t454 + t510) * t579) * t535) * pkin(7) + t483) * t373) * t326 + (0.32e2 * (t500 + (-0.4e1 * t418 * t522 + t585 + (t586 + t430 + t442) * t459) * t387 + (-t459 + t479 + t589) * t518 + t355 * t564 + t381 * t358) * t569 + 0.8e1 * (t353 + (t358 * t565 + t344) * t521 + (t481 + (t437 + t548) * t454 + t477) * t517 + t480) * t373) * t350) * t362) / ((-0.4e1 * (-t542 * t374 + 0.2e1 * t575 + (t524 * t588 + t418 * (t454 + t584)) * pkin(1)) * t568 + 0.8e1 * pkin(7) * t493 + ((pkin(3) * t585 + 0.8e1 * t459 * t465) * t420 + 0.4e1 * t462 * t515) * t387 - 0.4e1 * t482 * t537 - (t454 * t550 + t435 + t540 + (6 * t560)) * t374 + (t432 + (t425 + 6 * t461) * t454 + t372) * t579) * t326 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t495 + 0.4e1 / 0.9e1 * t454 - t449 / 0.9e1 + t553) * t562 + t370 * t355 + t356 * t564 + (t559 + (t396 + t405 + t479) * t375) * t539) * t568 + t353 + (t356 * t565 + t344) * t521 + t489 * t517 + ((t371 + t509) * t454 + t557) * t491 + t451 + (t424 + t506) * t466 + (t506 * t437 + t546 * t445 + t423 + t441) * t454 + t372 * t359) * t350 + ((t530 * t579 + (t523 * t587 + (t374 - t579) * t364 + t482) * t533) * t326 + (0.8e1 * (t364 + t526 + t381) * t569 + 0.6e1 * (t544 * t526 + (t355 + t356) * t518 + t489) * t373) * t350) * t362);
t335 = 0.1e1 / (t361 * t536 + t355 + t363 + t504);
t455 = 0.1e1 / pkin(3);
t576 = t335 * t455;
t484 = -pkin(1) + t525;
t485 = t433 + t438 - t545;
t324 = (t361 * t418 + t524) * t326 - (t363 + t485 + t490) * t574 + t484 * t363 + t485 * t525 + (t343 * t587 - t505 + t529 - t586) * pkin(1);
t348 = t433 + t488;
t325 = (-t484 + t574) * t326 + (t343 * t536 + t348 * t361) * t418 + (t348 * t421 + (0.4e1 * t387 - 0.2e1) * t361 * pkin(1)) * t374;
t502 = t576 / 0.2e1;
t322 = qJ(1) + atan2(t325 * t502, t324 * t502);
t538 = m(3) + m(7) + m(9);
t450 = 0.1e1 / pkin(4);
t516 = t335 * t450 / pkin(3) ^ 2;
t501 = t450 * t455 / 0.2e1;
t306 = atan2(t326 * t501, t501 * t577) + t322;
t341 = t414 * t421 - t415 * t418;
t318 = (-t324 * t341 / 0.2e1 + t325 * t581) * t576;
t319 = (t324 * t581 + t325 * t341 / 0.2e1) * t576;
t419 = sin(pkin(9));
t422 = cos(pkin(9));
t315 = -t318 * t422 + t319 * t419;
t316 = t318 * t419 + t319 * t422;
t311 = atan2(t315, t316) + t322;
t503 = t577 / 0.4e1;
t321 = cos(t322);
t496 = -pkin(3) * t321 - t375;
t492 = -m(10) * pkin(6) - mrSges(4,1) - mrSges(7,1);
t346 = -t417 * t421 + t567;
t345 = -t417 * t418 - t566;
t334 = qJ(1) + atan2(t341, t342);
t333 = pkin(8) + atan2(t341, -t342);
t332 = cos(t334);
t331 = sin(t334);
t330 = cos(t333);
t329 = sin(t333);
t320 = sin(t322);
t310 = cos(t311);
t309 = sin(t311);
t308 = (t324 * t326 / 0.4e1 + t325 * t503) * t516;
t307 = (t324 * t503 - t325 * t326 / 0.4e1) * t516;
t305 = cos(t306);
t304 = sin(t306);
t303 = atan2(t315, -t316) + t311;
t302 = cos(t303);
t301 = sin(t303);
t300 = atan2(t307 * t345 + t308 * t346, -t307 * t346 + t308 * t345) + t306;
t299 = cos(t300);
t298 = sin(t300);
t1 = (-mrSges(2,3) - mrSges(6,3) - mrSges(8,3) - mrSges(10,3) - mrSges(11,3) - mrSges(7,3) - mrSges(4,3) - mrSges(3,3) - mrSges(5,3) - mrSges(9,3) - mrSges(1,3)) * g(3) + (-t414 * t583 - t298 * mrSges(11,1) - mrSges(6,1) * t329 + mrSges(8,1) * t420 - mrSges(9,1) * t331 + t301 * mrSges(10,1) + mrSges(2,2) * t421 - t299 * mrSges(11,2) + mrSges(3,2) * t321 + t305 * mrSges(5,2) - mrSges(6,2) * t330 - mrSges(8,2) * t417 - mrSges(9,2) * t332 + t302 * mrSges(10,2) - mrSges(1,2) - t578 * t310 + (m(11) * pkin(4) + mrSges(5,1)) * t304 + t492 * t309 + (pkin(2) * t582 + pkin(3) * t580 + mrSges(3,1)) * t320 + (mrSges(2,1) + (t538 + t580 + t582) * pkin(1)) * t418) * g(2) + (-m(11) * (-pkin(4) * t305 + t496) - m(5) * t496 + t492 * t310 + t298 * mrSges(11,2) - t299 * mrSges(11,1) - mrSges(8,1) * t417 - mrSges(8,2) * t420 - t301 * mrSges(10,2) + t305 * mrSges(5,1) + t578 * t309 - t415 * t583 + t302 * mrSges(10,1) + (pkin(1) * t538 + mrSges(2,1)) * t421 - (m(8) * pkin(7)) - t304 * mrSges(5,2) - mrSges(6,1) * t330 + mrSges(6,2) * t329 - mrSges(2,2) * t418 + mrSges(3,1) * t321 - mrSges(3,2) * t320 - mrSges(9,1) * t332 + mrSges(9,2) * t331 - mrSges(1,1) + t582 * (pkin(2) * t321 + t375)) * g(1);
U = t1;
