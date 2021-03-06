% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
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
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = fourbar1turnOL_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:41:15
% EndTime: 2020-04-12 19:41:18
% DurationCPUTime: 3.18s
% Computational Cost: add. (2559->306), mult. (6084->418), div. (0->0), fcn. (4403->8), ass. (0->225)
t516 = sin(qJ(3));
t509 = qJDD(2) + qJDD(3);
t520 = cos(qJ(3));
t521 = cos(qJ(2));
t555 = t520 * t521;
t517 = sin(qJ(2));
t567 = t516 * t517;
t453 = (-t555 + t567) * qJD(1);
t566 = t516 * t521;
t454 = (-t517 * t520 - t566) * qJD(1);
t579 = t453 * t454;
t587 = t509 + t579;
t590 = t516 * t587;
t589 = t520 * t587;
t518 = sin(qJ(1));
t522 = cos(qJ(1));
t488 = t522 * g(1) + t518 * g(2);
t449 = -t517 * g(3) - t521 * t488;
t514 = t521 ^ 2;
t525 = qJD(1) ^ 2;
t506 = t514 * t525;
t524 = qJD(2) ^ 2;
t584 = -t506 - t524;
t418 = t584 * pkin(2) + t449;
t448 = t521 * g(3) - t517 * t488;
t497 = t521 * t525 * t517;
t543 = qJDD(2) + t497;
t526 = pkin(2) * t543 - t448;
t388 = t516 * t418 - t520 * t526;
t556 = t520 * t418;
t389 = t516 * t526 + t556;
t366 = t520 * t388 - t516 * t389;
t588 = t517 * t366;
t487 = t518 * g(1) - t522 * g(2);
t463 = t522 * t487;
t586 = -t518 * t488 + t463;
t504 = t521 * qJDD(1);
t549 = qJD(1) * qJD(2);
t538 = t517 * t549;
t473 = t504 - t538;
t585 = t473 - t538;
t451 = t453 ^ 2;
t452 = t454 ^ 2;
t510 = qJD(2) + qJD(3);
t508 = t510 ^ 2;
t466 = -t525 * pkin(1) - t488;
t515 = sin(qJ(4));
t519 = cos(qJ(4));
t444 = t519 * g(3) + t515 * t466;
t445 = -t515 * g(3) + t519 * t466;
t395 = t519 * t444 - t515 * t445;
t583 = pkin(1) * t395;
t502 = t521 * t549;
t546 = t517 * qJDD(1);
t470 = t502 + t546;
t528 = -t453 * qJD(3) + t520 * t470 + t516 * t473;
t578 = t510 * t453;
t381 = t528 + t578;
t391 = -t454 * qJD(3) + t516 * t470 - t520 * t473;
t577 = t510 * t454;
t527 = t391 + t577;
t352 = -t520 * t381 - t516 * t527;
t582 = pkin(2) * t352;
t581 = pkin(2) * t366;
t417 = t585 * pkin(2) + t487;
t580 = pkin(2) * t417;
t576 = t510 * t516;
t575 = t510 * t520;
t511 = t515 ^ 2;
t574 = t511 * t525;
t512 = t517 ^ 2;
t573 = t512 * t525;
t465 = qJDD(1) * pkin(1) + t487;
t572 = t515 * t465;
t495 = t519 * t525 * t515;
t483 = qJDD(4) + t495;
t571 = t515 * t483;
t484 = qJDD(4) - t495;
t570 = t515 * t484;
t406 = -t579 + t509;
t569 = t516 * t406;
t568 = t516 * t417;
t565 = t517 * t543;
t461 = t517 * t487;
t564 = t518 * t417;
t563 = t518 * t465;
t562 = t518 * t487;
t560 = t519 * t465;
t559 = t519 * t484;
t558 = t520 * t406;
t557 = t520 * t417;
t486 = qJDD(2) - t497;
t554 = t521 * t486;
t462 = t521 * t487;
t553 = t522 * t417;
t552 = t522 * t465;
t513 = t519 ^ 2;
t551 = t511 + t513;
t550 = t512 + t514;
t548 = qJD(1) * qJD(4);
t547 = t515 * qJDD(1);
t545 = t518 * qJDD(1);
t503 = t519 * qJDD(1);
t544 = t522 * qJDD(1);
t542 = t521 * t580;
t541 = t518 * t579;
t540 = t522 * t579;
t539 = t515 * t548;
t501 = t519 * t548;
t396 = t515 * t444 + t519 * t445;
t537 = t517 * t448 + t521 * t449;
t536 = -t522 * t488 - t562;
t535 = t518 * t495;
t534 = t518 * t497;
t533 = t522 * t495;
t532 = t522 * t497;
t478 = -t518 * t525 + t544;
t531 = -t518 * g(3) - pkin(5) * t478;
t529 = t516 * t388 + t520 * t389;
t382 = -t528 + t578;
t523 = qJD(4) ^ 2;
t505 = t513 * t525;
t494 = t506 - t524;
t493 = -t505 - t523;
t492 = t505 - t523;
t491 = t524 - t573;
t490 = -t523 - t574;
t489 = t523 - t574;
t482 = -t506 + t573;
t481 = t506 + t573;
t480 = -t505 + t574;
t479 = t505 + t574;
t477 = t522 * t525 + t545;
t476 = t550 * qJDD(1);
t475 = t551 * qJDD(1);
t474 = t504 - 0.2e1 * t538;
t472 = t503 - 0.2e1 * t539;
t471 = t503 - t539;
t469 = 0.2e1 * t502 + t546;
t468 = t501 + t547;
t467 = 0.2e1 * t501 + t547;
t460 = t519 * t483;
t459 = t550 * t549;
t458 = t551 * t548;
t457 = pkin(1) * t465;
t450 = t522 * g(3) - pkin(5) * t477;
t442 = -t452 + t508;
t441 = t451 - t508;
t440 = t521 * t470 - t512 * t549;
t439 = t519 * t468 - t511 * t548;
t438 = -t517 * t473 - t514 * t549;
t437 = -t515 * t471 - t513 * t548;
t435 = -t452 - t508;
t434 = -t554 - t517 * (-t524 - t573);
t433 = -t517 * t491 + t521 * t543;
t432 = t521 * t584 - t565;
t431 = -t517 * t486 + t521 * t494;
t430 = -t515 * t490 - t559;
t429 = -t515 * t489 + t460;
t428 = t519 * t493 - t571;
t427 = t519 * t492 - t570;
t426 = t521 * t491 + t565;
t425 = t517 * t494 + t554;
t424 = (t470 + t502) * t517;
t423 = t519 * t489 + t571;
t422 = t515 * t492 + t559;
t421 = (t468 + t501) * t515;
t420 = t585 * t521;
t419 = (t471 - t539) * t519;
t416 = pkin(1) * t472 + t560;
t415 = -pkin(1) * t467 - t572;
t414 = -t517 * t469 + t521 * t474;
t413 = -t515 * t467 + t519 * t472;
t412 = t521 * t469 + t517 * t474;
t411 = t519 * t467 + t515 * t472;
t408 = t452 - t451;
t404 = -t508 - t451;
t402 = t521 * t448 - t517 * t449;
t400 = (-t453 * t520 - t454 * t516) * t510;
t399 = (-t453 * t516 + t454 * t520) * t510;
t398 = -pkin(1) * (t519 * t490 - t570) + t445;
t397 = -pkin(1) * (t515 * t493 + t460) + t444;
t393 = -t451 - t452;
t390 = pkin(1) * t479 + t396;
t387 = -t520 * t441 + t569;
t386 = t516 * t442 - t589;
t385 = -t516 * t441 - t558;
t384 = -t520 * t442 - t590;
t383 = -t520 * t435 + t569;
t379 = t391 - t577;
t376 = t454 * t576 + t520 * t528;
t375 = -t454 * t575 + t516 * t528;
t374 = t516 * t391 + t453 * t575;
t373 = -t520 * t391 + t453 * t576;
t372 = -t516 * t404 - t589;
t371 = -t517 * t399 + t521 * t400;
t370 = t521 * t399 + t517 * t400;
t369 = pkin(2) * t379 - t557;
t368 = -pkin(2) * t382 + t568;
t367 = t556 - t516 * t448 + (t516 * t543 + t383) * pkin(2);
t364 = -t517 * t385 + t521 * t387;
t363 = -t517 * t384 + t521 * t386;
t362 = t521 * t385 + t517 * t387;
t361 = t521 * t384 + t517 * t386;
t360 = t521 * (t516 * t435 + t558) - t517 * t383;
t359 = -t517 * t368 + t417 * t555;
t358 = -t517 * t369 + t417 * t566;
t357 = t521 * t368 + t517 * t557;
t356 = t521 * t369 + t417 * t567;
t355 = pkin(2) * t372 + t388;
t354 = -t520 * t379 + t516 * t382;
t353 = -t516 * t379 - t520 * t382;
t351 = -t517 * t375 + t521 * t376;
t350 = -t517 * t373 + t521 * t374;
t349 = t521 * t375 + t517 * t376;
t348 = t521 * t373 + t517 * t374;
t347 = t521 * (-t520 * t404 + t590) - t517 * t372;
t346 = -pkin(2) * t393 + t529;
t345 = t521 * t529 + t588;
t344 = -t517 * t353 + t521 * t354;
t343 = t521 * (t516 * t381 - t520 * t527) - t517 * t352;
t342 = t521 * t353 + t517 * t354;
t341 = -t517 * t346 + t521 * t366;
t340 = t521 * t346 + t588;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t478, 0, -t477, 0, t531, -t450, -t586, -pkin(5) * t586, t522 * t440 - t534, t522 * t414 + t518 * t482, t522 * t433 + t517 * t545, t522 * t438 + t534, t522 * t431 + t518 * t504, t518 * qJDD(2) + t522 * t459, -t517 * t463 - t518 * t448 - pkin(5) * (t518 * t432 + t522 * t474), -t521 * t463 - t518 * t449 - pkin(5) * (t518 * t434 - t522 * t469), t522 * t402 - pkin(5) * (t518 * t476 + t522 * t481), -pkin(5) * (t518 * t537 + t463), t522 * t351 - t541, t522 * t344 + t518 * t408, t522 * t363 - t518 * t381, t522 * t350 + t541, t522 * t364 + t518 * t527, t522 * t371 + t518 * t509, t522 * t358 + t518 * t355 - pkin(5) * (t518 * t347 + t522 * t379), t522 * t359 + t518 * t367 - pkin(5) * (t518 * t360 - t382 * t522), t522 * t341 + t518 * t582 - pkin(5) * (t518 * t343 - t522 * t393), -pkin(5) * (t518 * t345 + t553) + (-t366 * t518 - t517 * t553) * pkin(2), t522 * t439 - t535, t522 * t413 + t518 * t480, t522 * t429 + t515 * t545, t522 * t437 + t535, t522 * t427 + t518 * t503, t518 * qJDD(4) + t522 * t458, -t515 * t552 - t518 * t397 - pkin(5) * (t518 * t428 + t522 * t472), -t519 * t552 - t518 * t398 - pkin(5) * (t518 * t430 - t522 * t467), t522 * t395 - pkin(5) * (t518 * t475 + t522 * t479), -t518 * t583 - pkin(5) * (t518 * t396 + t552); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t477, 0, t478, 0, t450, t531, t536, pkin(5) * t536, t518 * t440 + t532, t518 * t414 - t522 * t482, t518 * t433 - t517 * t544, t518 * t438 - t532, t518 * t431 - t521 * t544, -t522 * qJDD(2) + t518 * t459, -t517 * t562 + t522 * t448 + pkin(5) * (t522 * t432 - t518 * t474), -t518 * t462 + t522 * t449 + pkin(5) * (t522 * t434 + t518 * t469), t518 * t402 + pkin(5) * (t522 * t476 - t518 * t481), pkin(5) * (t522 * t537 - t562), t518 * t351 + t540, t518 * t344 - t522 * t408, t518 * t363 + t522 * t381, t518 * t350 - t540, t518 * t364 - t522 * t527, t518 * t371 - t522 * t509, t518 * t358 - t522 * t355 + pkin(5) * (t522 * t347 - t518 * t379), t518 * t359 - t522 * t367 + pkin(5) * (t522 * t360 + t382 * t518), t518 * t341 - t522 * t582 + pkin(5) * (t522 * t343 + t518 * t393), pkin(5) * (t522 * t345 - t564) + (t366 * t522 - t517 * t564) * pkin(2), t518 * t439 + t533, t518 * t413 - t522 * t480, t518 * t429 - t515 * t544, t518 * t437 - t533, t518 * t427 - t519 * t544, -t522 * qJDD(4) + t518 * t458, -t515 * t563 + t522 * t397 + pkin(5) * (t522 * t428 - t518 * t472), -t518 * t560 + t522 * t398 + pkin(5) * (t522 * t430 + t518 * t467), t518 * t395 + pkin(5) * (t522 * t475 - t518 * t479), t522 * t583 + pkin(5) * (t522 * t396 - t563); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t487, t488, 0, 0, t424, t412, t426, t420, t425, 0, t462, -t461, t537, 0, t349, t342, t361, t348, t362, t370, t356, t357, t340, t542, t421, t411, t423, t419, t422, 0, t416, t415, t390, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t525, 0, 0, -g(3), -t487, 0, t440, t414, t433, t438, t431, t459, -t461, -t462, t402, 0, t351, t344, t363, t350, t364, t371, t358, t359, t341, -t517 * t580, t439, t413, t429, t437, t427, t458, -t572, -t560, t395, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t525, 0, qJDD(1), 0, g(3), 0, -t488, 0, t497, -t482, -t546, -t497, -t504, -qJDD(2), t448, t449, 0, 0, t579, -t408, t381, -t579, -t527, -t509, -t355, -t367, -t582, t581, t495, -t480, -t547, -t495, -t503, -qJDD(4), t397, t398, 0, t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t487, t488, 0, 0, t424, t412, t426, t420, t425, 0, t462, -t461, t537, 0, t349, t342, t361, t348, t362, t370, t356, t357, t340, t542, t421, t411, t423, t419, t422, 0, t416, t415, t390, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, t474, t543, -t502, t494, t502, 0, -t487, t448, 0, t376, t354, t386, t374, t387, t400, t568, t557, t366, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t538, t469, t491, t473, t486, -t538, t487, 0, t449, 0, t375, t353, t384, t373, t385, t399, t369, t368, t346, t580, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t497, t482, t546, t497, t504, qJDD(2), -t448, -t449, 0, 0, -t579, t408, -t381, t579, t527, t509, t355, t367, t582, -t581, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t528, t379, t587, -t578, t441, t578, 0, -t417, -t388, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t382, t442, t391, t406, -t577, t417, 0, -t389, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, t408, -t381, t579, t527, t509, t388, t389, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t468, t472, t483, -t501, t492, t501, 0, -t465, t444, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t539, t467, t489, t471, t484, -t539, t465, 0, t445, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t495, t480, t547, t495, t503, qJDD(4), -t444, -t445, 0, 0;];
m_new_reg = t1;
