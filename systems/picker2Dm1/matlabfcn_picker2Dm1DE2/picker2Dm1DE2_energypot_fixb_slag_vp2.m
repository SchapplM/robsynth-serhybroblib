% Calculate potential energy for
% picker2Dm1DE2
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
% Datum: 2020-05-11 05:26
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = picker2Dm1DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(9,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'picker2Dm1DE2_energypot_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'picker2Dm1DE2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1DE2_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1DE2_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-10 20:21:30
% EndTime: 2020-05-10 20:21:43
% DurationCPUTime: 7.99s
% Computational Cost: add. (65181->483), mult. (176712->593), div. (2626->11), fcn. (45644->38), ass. (0->266)
t661 = 4 * pkin(1);
t507 = pkin(4) ^ 2;
t449 = -t507 / 0.4e1;
t512 = pkin(3) ^ 2;
t660 = t449 + t512 / 0.2e1;
t659 = 2 * pkin(7);
t500 = 2 * pkin(2);
t466 = sin(pkin(8));
t467 = cos(pkin(8));
t470 = sin(qJ(1));
t473 = cos(qJ(1));
t393 = t466 * t473 - t467 * t470;
t658 = 0.2e1 * t393;
t439 = t473 ^ 2;
t657 = -0.2e1 * t439;
t485 = 0.2e1 * t512;
t656 = 0.4e1 * t512;
t517 = pkin(1) ^ 2;
t515 = t517 ^ 2;
t655 = 4 * t515;
t654 = 2 * t517;
t489 = 6 * t517;
t520 = pkin(7) ^ 2;
t497 = 2 * t520;
t525 = t512 ^ 2;
t484 = 0.5e1 * t525;
t653 = pkin(5) * m(6);
t652 = (m(4) + m(10));
t504 = pkin(5) ^ 2;
t394 = -t466 * t470 - t467 * t473;
t647 = pkin(5) * t394;
t593 = pkin(1) * t647;
t385 = t504 - t593;
t389 = -pkin(1) + t647;
t501 = 2 * pkin(1);
t620 = -0.2e1 * t593 + t504;
t378 = sqrt(-(-(t501 + pkin(5)) * pkin(5) + t620) * (pkin(5) * (t501 - pkin(5)) + t620));
t643 = t378 * t393;
t376 = -pkin(5) * t643 - 0.2e1 * t385 * t389;
t651 = t376 / 0.4e1;
t650 = 0.1e1 / pkin(6) / 0.2e1;
t649 = m(11) + m(5);
t427 = pkin(1) * t473;
t414 = t427 + pkin(7);
t648 = pkin(1) * t470;
t469 = sin(qJ(2));
t425 = pkin(3) * t469;
t472 = cos(qJ(2));
t426 = pkin(3) * t472;
t415 = t425 * t659;
t436 = t469 ^ 2;
t632 = t436 * t512;
t590 = 0.2e1 * t632;
t606 = -t512 + t520;
t395 = t415 + t590 + t606;
t631 = t470 * t472;
t586 = pkin(3) * t631;
t554 = pkin(1) * t586;
t407 = -0.2e1 * t554;
t499 = 0.2e1 * pkin(3);
t431 = t517 + t520;
t566 = -t507 + t431;
t547 = t415 + t566;
t626 = t517 * t439;
t588 = -0.4e1 * t626;
t594 = -0.4e1 * t425;
t604 = t517 - t520;
t607 = t507 - t520;
t413 = t425 + pkin(7);
t638 = t413 * t473;
t375 = sqrt(t395 * t588 + 0.4e1 * t604 * t632 + pkin(7) * t566 * t594 - t515 + (-0.2e1 * t512 + t607) * t654 - (t520 - (t499 + pkin(4)) * pkin(4)) * (t520 + (t499 - pkin(4)) * pkin(4)) + (-(t407 + t547) * t638 + t547 * t586) * t661);
t543 = -pkin(1) + t586;
t490 = 3 * t517;
t544 = t485 + t490 - t607;
t549 = -0.4e1 * t554;
t630 = t472 * t473;
t585 = pkin(3) * t630;
t373 = (t413 * t470 + t585) * t375 - (t415 + t544 + t549) * t638 + t543 * t415 + t544 * t586 + (t395 * t657 - t566 + t590 - t656) * pkin(1);
t400 = t485 + t547;
t597 = 0.2e1 * t427;
t374 = (-t543 + t638) * t375 + (t395 * t597 + t400 * t413) * t470 + (t400 * t473 + (0.4e1 * t439 - 0.2e1) * t413 * pkin(1)) * t426;
t377 = pkin(5) * t385 * t658 - t378 * t389;
t384 = 0.1e1 / (t517 + t620);
t641 = t384 / pkin(5);
t565 = t512 + t431;
t381 = 0.1e1 / (t413 * t597 + t407 + t415 + t565);
t513 = 0.1e1 / pkin(3);
t642 = t381 * t513;
t551 = t641 * t642;
t360 = (t373 * t651 + t374 * t377 / 0.4e1) * t551;
t361 = (-t373 * t377 / 0.4e1 + t374 * t651) * t551;
t471 = sin(pkin(9));
t474 = cos(pkin(9));
t358 = t360 * t474 + t361 * t471;
t646 = pkin(6) * t358;
t359 = t360 * t471 - t361 * t474;
t645 = pkin(6) * t359;
t506 = t507 ^ 2;
t519 = t520 ^ 2;
t605 = t515 + t519;
t610 = t497 - t507;
t625 = t520 * t507;
t537 = t610 * t517 + t506 / 0.6e1 + t605 - t625;
t392 = -t525 / 0.6e1 + t537;
t460 = -t512 / 0.3e1;
t421 = t460 + t520;
t396 = t421 * t407;
t402 = t425 + t414;
t430 = -0.3e1 * t512 + t520;
t438 = t473 * t439;
t521 = pkin(1) * t517;
t623 = t521 * t438;
t592 = pkin(7) * t623;
t559 = 0.8e1 * t592;
t405 = t430 * t559;
t429 = -t507 - t512;
t496 = 3 * t520;
t417 = t496 + t429;
t637 = t417 * t517;
t406 = 0.10e2 * t637;
t456 = 0.4e1 / 0.3e1 * t512;
t450 = -t507 / 0.3e1;
t571 = t450 + t431;
t408 = t456 + t571;
t451 = -t507 / 0.2e1;
t410 = t451 + t565;
t411 = -t507 + t565;
t416 = pkin(7) * t597;
t495 = 4 * t520;
t419 = (t495 + t507) * t517;
t422 = -t517 / 0.3e1 + t520;
t423 = 0.10e2 / 0.3e1 * t517;
t424 = t431 ^ 2;
t428 = -0.30e2 * t507 + (60 * t520);
t433 = -3 * t517 + t520;
t448 = -t507 / 0.6e1;
t457 = 0.2e1 / 0.3e1 * t512;
t462 = 0.4e1 / 0.3e1 * t517;
t464 = t517 / 0.2e1;
t475 = 15 * t515;
t476 = 15 * t517;
t477 = 10 * t517;
t482 = -0.2e1 * t507;
t483 = -0.5e1 * t507;
t486 = 7 * t515;
t487 = 5 * t515;
t488 = 7 * t517;
t493 = 3 * t519;
t494 = 8 * t520;
t524 = pkin(3) * t512;
t509 = t524 ^ 2;
t529 = pkin(7) * t520;
t536 = 0.5e1 / 0.6e1 * t525 + t537;
t538 = t520 - t554;
t618 = t506 / 0.2e1 - t525 / 0.2e1;
t546 = -0.3e1 * t625 + t493 + t618;
t550 = -0.6e1 * t554;
t453 = -0.3e1 / 0.2e1 * t507;
t617 = t453 + t496;
t621 = t431 * ((t453 + t497) * t517 - 0.3e1 / 0.2e1 * t625 + t605 + t618) + t509;
t539 = ((t423 + t610) * t512 + t536) * t550 + (t475 + (-0.9e1 * t507 + (18 * t520)) * t517 + t546) * t512 + (t476 + t617) * t525 + t621;
t540 = t410 * t549;
t611 = t490 + t520;
t568 = t512 + t611;
t541 = -(0.3e1 * t512 + t431) * t648 + t568 * t426;
t452 = -0.2e1 / 0.3e1 * t507;
t461 = -0.2e1 / 0.3e1 * t512;
t569 = t452 + t431;
t612 = t477 + t497;
t616 = t461 + t520;
t542 = -(t484 + ((5 * t517) + t417) * t485 + (t461 + t569) * t431) * t648 + (t525 + (t452 + t461 + t612) * t512 + t487 + 0.2e1 * t637 + t520 * (t452 + t616)) * t426;
t609 = t506 - t525;
t545 = -0.6e1 * t625 + (6 * t519) + t609;
t570 = t452 + t457 + t497;
t619 = (t457 + t569) * t431 + t525;
t548 = t408 * t549 + (t489 + t570) * t512 + t619;
t583 = t521 * t426;
t552 = t438 * t583;
t627 = t515 * t439 ^ 2;
t553 = t627 * t426;
t581 = 0.16e2 * t623;
t557 = pkin(7) * t581;
t558 = 0.20e2 / 0.3e1 * t517;
t608 = -t507 + t512;
t567 = t496 + t608;
t458 = t512 / 0.3e1;
t572 = t448 + t458 + t520;
t573 = t507 / 0.3e1 + t458 + t497;
t574 = 0.2e1 / 0.3e1 * t507 + t457 + t495;
t575 = 0.4e1 / 0.3e1 * t507 + t456 - (2 * t520);
t629 = (pkin(7) + pkin(3)) * (pkin(7) - pkin(3));
t576 = t470 * t629;
t595 = 0.6e1 * t427;
t578 = pkin(7) * t595;
t596 = 0.4e1 * t427;
t579 = pkin(7) * t596;
t580 = -t648 / 0.2e1;
t582 = 0.12e2 * t626;
t584 = t517 * t426;
t587 = 0.4e1 * t626;
t589 = 0.8e1 * t627;
t633 = t469 * t436 * t524;
t591 = -0.8e1 * t633;
t598 = 0.2e1 * t648;
t599 = pkin(7) * t427;
t601 = 4 * pkin(7);
t602 = t519 + t525;
t603 = t519 - t515;
t613 = 0.4e1 / 0.7e1 * t520 - t507 / 0.7e1;
t614 = t464 + t520;
t615 = t517 / 0.3e1 + t520;
t624 = t520 * t517;
t628 = (pkin(7) + pkin(1)) * (pkin(7) - pkin(1));
t634 = t436 ^ 2 * t525;
t635 = t429 * t520;
t636 = t424 * (-t512 + t566);
t639 = (-t470 * t521 + t584) * t439;
t644 = ((-0.24e2 * (0.4e1 / 0.3e1 * t626 + t416 + t422) * t634 * t648 - 0.12e2 * (-0.8e1 / 0.3e1 * t553 + ((t462 + t572) * t426 - (0.7e1 / 0.6e1 * t512 + t448 + t614) * t648) * t587 + (-t512 * t604 - 0.5e1 / 0.3e1 * t515 + t573 * t517 + t520 * (t450 + t421)) * t426 + (-t525 + (-t558 + t574) * t512 - (3 * t515) + t575 * t517 + t519) * t580 + (-t470 * t515 * t438 + ((t517 + t572) * t426 + (t485 - t604) * t580) * t427) * t601) * t632 + 0.24e2 * t421 * t553 + ((t520 + 0.5e1 / 0.2e1 * t512 + 0.3e1 / 0.2e1 * t517 + t451) * t426 + t430 * t648 / 0.2e1) * t557 - 0.6e1 * ((-0.3e1 * t525 + (-t558 + t575) * t512 + t574 * t517 + t603) * t426 - 0.2e1 * (-0.5e1 / 0.3e1 * t525 + (-t517 + t573) * t512 + t520 * (t460 + t571)) * t648) * t626 - 0.6e1 * t542 * t599 - (t509 + ((21 * t517) + t417) * t525 + (t406 + t493 + (35 * t515) + 0.2e1 * t635) * t512 + (t486 + (t483 + t494 - 0.5e1 * t512) * t517 + t520 * (-t507 + t606)) * t431) * t426 + (0.7e1 * t509 + (t488 + t417) * t484 + (t406 + (21 * t515) + (9 * t519) + 0.6e1 * t635) * t512 + t636) * t648) * t375 + (0.16e2 * (t589 + t557 + (-8 * t515 + 12 * t624) * t439 + (-12 * pkin(7) * t521 + t529 * t661) * t473 - (6 * t624) + t605) * t634 + 0.24e2 * (t616 * t589 + 0.14e2 * (-0.32e2 / 0.21e2 * (t520 + t512 / 0.4e1 + t517 / 0.4e1 - t507 / 0.8e1) * t554 + 0.5e1 / 0.42e2 * t525 + (0.16e2 / 0.21e2 * t517 + t613) * t512 + t515 / 0.7e1 + t613 * t517 + t519 - 0.3e1 / 0.7e1 * t625 + t506 / 0.42e2) * t626 + t422 * t540 - t604 * t525 + (t419 - 0.10e2 / 0.3e1 * t515 + (2 * t519) - t625) * t512 + t392 * t628 + ((-0.2e1 / 0.3e1 * t554 + t449 + t614) * t581 + (-0.8e1 / 0.3e1 * (t615 + t660) * t554 + 0.5e1 / 0.18e2 * t525 + (0.4e1 / 0.3e1 * t520 + t462 + t450) * t512 + t519 + 0.2e1 / 0.3e1 * t624 - 0.2e1 / 0.3e1 * t625 - t515 / 0.3e1 + t506 / 0.18e2) * t595) * pkin(7)) * t632 + 0.16e2 * (-0.6e1 * t512 * t520 + t602) * t627 + 0.32e2 * (t407 * t629 + t410 * t430) * t592 + 0.24e2 * (t421 * t540 - t509 + (-t423 + t607) * t525 + (t419 + t525 / 0.6e1 - t506 / 0.6e1 + t603) * t512 + t392 * t520) * t626 + 0.8e1 * t539 * t599 - 0.8e1 * ((t488 + t617) * t525 + (t486 + (t483 + (10 * t520)) * t517 + t546) * t512 + t621) * t554 + t525 ^ 2 + (t482 + t495 + (28 * t517)) * t509 + (t428 * t517 + (70 * t515) + t545) * t525 + (t428 * t515 + t545 * t489 + t609 * t497 - 0.6e1 * t519 * t507 + (28 * t521 ^ 2) + (4 * t529 ^ 2)) * t512 + t411 * t636) * t402 + (((0.4e1 * t639 + (t426 + t598) * t416 + t433 * t426 + (t451 + t568) * t598) * t591 - 0.6e1 * ((0.2e1 * (0.5e1 / 0.6e1 * t512 + t464 + t448) * t426 + pkin(1) * t576) * t588 + (-0.8e1 * t552 + ((t450 + t457 + t611) * t426 - (0.8e1 / 0.3e1 * t512 + t571) * t648) * t596) * pkin(7) + t542) * t425) * t375 + (0.32e2 * (t559 + (-0.4e1 * t470 * t583 + t655 + (t656 + t482 + t494) * t517) * t439 + (-t517 + t538 + t660) * t579 + t407 * t628 + t433 * t410) * t633 + 0.8e1 * (t405 + (t410 * t629 + t396) * t582 + (t540 + (t489 + t610) * t512 + t536) * t578 + t539) * t425) * t402) * t414) / ((-0.4e1 * (-t604 * t426 + 0.2e1 * t639 + (t585 * t659 + t470 * (t512 + t654)) * pkin(1)) * t632 + 0.8e1 * pkin(7) * t552 + ((pkin(3) * t655 + 0.8e1 * t517 * t524) * t472 + 0.4e1 * t521 * t576) * t439 - 0.4e1 * t541 * t599 - (t512 * t612 + t487 + t602 + (6 * t624)) * t426 + (t484 + (t477 + 6 * t520) * t512 + t424) * t648) * t375 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t554 + 0.4e1 / 0.9e1 * t512 - t507 / 0.9e1 + t615) * t626 + t422 * t407 + t408 * t628 + (t623 + (t448 + t457 + t538) * t427) * t601) * t632 + t405 + (t408 * t629 + t396) * t582 + t548 * t578 + ((t423 + t570) * t512 + t619) * t550 + t509 + (t476 + t567) * t525 + (t567 * t489 + t608 * t497 + t475 + t493) * t512 + t424 * t411) * t402 + ((t591 * t648 + (t584 * t657 + (t426 - t648) * t416 + t541) * t594) * t375 + (0.8e1 * (t416 + t587 + t433) * t633 + 0.6e1 * (t606 * t587 + (t407 + t408) * t579 + t548) * t425) * t402) * t414);
t640 = t384 / pkin(1);
t562 = t642 / 0.2e1;
t365 = qJ(1) + atan2(t374 * t562, t373 * t562);
t357 = t646 * t500;
t502 = pkin(6) ^ 2;
t622 = t357 + t502;
t600 = m(3) + m(7) + m(9);
t508 = 0.1e1 / pkin(4);
t577 = t381 * t508 / pkin(3) ^ 2;
t346 = sqrt(-(-(t500 + pkin(6)) * pkin(6) + t622) * (pkin(6) * (t500 - pkin(6)) + t622));
t355 = t357 + 0.2e1 * t502;
t356 = -pkin(2) - t646;
t563 = 0.1e1 / ((pkin(2) ^ 2) + t622) * t650;
t338 = atan2((-t346 * t356 + t355 * t645) * t563, (-t346 * t645 - t355 * t356) * t563) + t365;
t560 = t508 * t513 / 0.2e1;
t349 = atan2(t375 * t560, t560 * t644) + t365;
t564 = t644 / 0.4e1;
t561 = t641 / 0.2e1;
t364 = cos(t365);
t556 = -pkin(2) * t364 - t427;
t555 = -pkin(3) * t364 - t427;
t398 = -t469 * t473 + t631;
t397 = -t469 * t470 - t630;
t390 = -pkin(1) * t394 + pkin(5);
t386 = t517 - t593;
t372 = qJ(1) + atan2(t377 * t561, t376 * t561);
t371 = pkin(8) + atan2((pkin(1) * t386 * t658 + t378 * t390) * t640 / 0.2e1, -(-pkin(1) * t643 + 0.2e1 * t386 * t390) * t640 / 0.2e1);
t370 = cos(t372);
t369 = sin(t372);
t368 = cos(t371);
t367 = sin(t371);
t363 = sin(t365);
t351 = (t373 * t375 / 0.4e1 + t374 * t564) * t577;
t350 = (t373 * t564 - t374 * t375 / 0.4e1) * t577;
t348 = cos(t349);
t347 = sin(t349);
t345 = atan2(t359, t358) + t365;
t344 = cos(t345);
t343 = sin(t345);
t342 = atan2(t350 * t397 + t351 * t398, -t350 * t398 + t351 * t397) + t349;
t341 = cos(t342);
t340 = sin(t342);
t337 = cos(t338);
t336 = sin(t338);
t335 = atan2(0.1e1 / pkin(2) * t346 * t650, -t358) + t338;
t334 = cos(t335);
t333 = sin(t335);
t1 = (-mrSges(10,3) - mrSges(6,3) - mrSges(3,3) - mrSges(9,3) - mrSges(5,3) - mrSges(1,3) - mrSges(8,3) - mrSges(7,3) - mrSges(2,3) - mrSges(4,3) - mrSges(11,3)) * g(3) + (-t340 * mrSges(11,1) + mrSges(3,2) * t364 + t333 * mrSges(10,1) + mrSges(8,1) * t472 - mrSges(8,2) * t469 + mrSges(2,2) * t473 - t337 * mrSges(4,2) + t348 * mrSges(5,2) + t334 * mrSges(10,2) - mrSges(6,1) * t367 - mrSges(6,2) * t368 - t466 * t653 - mrSges(9,1) * t369 - mrSges(9,2) * t370 - mrSges(7,1) * t343 - mrSges(7,2) * t344 - t341 * mrSges(11,2) - mrSges(1,2) + (m(11) * pkin(4) + mrSges(5,1)) * t347 + (-m(10) * pkin(6) - mrSges(4,1)) * t336 + ((pkin(2) * t652) + pkin(3) * t649 + mrSges(3,1)) * t363 + (mrSges(2,1) + (t600 + t649 + t652) * pkin(1)) * t470) * g(2) + (t340 * mrSges(11,2) + (pkin(1) * t600 + mrSges(2,1)) * t473 - t467 * t653 - m(11) * (-pkin(4) * t348 + t555) - m(5) * t555 - m(10) * (pkin(6) * t337 + t556) - m(4) * t556 - mrSges(8,1) * t469 - mrSges(8,2) * t472 - t333 * mrSges(10,2) - mrSges(2,2) * t470 - t337 * mrSges(4,1) + t348 * mrSges(5,1) - mrSges(9,1) * t370 + mrSges(9,2) * t369 - mrSges(7,1) * t344 + mrSges(7,2) * t343 + mrSges(3,1) * t364 - mrSges(3,2) * t363 + t334 * mrSges(10,1) - mrSges(6,1) * t368 + mrSges(6,2) * t367 - (m(8) * pkin(7)) + t336 * mrSges(4,2) - t347 * mrSges(5,2) - t341 * mrSges(11,1) - mrSges(1,1)) * g(1);
U = t1;
