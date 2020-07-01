% Calculate homogenous joint transformation matrices for
% hybBKspatial
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED,L1,L2]';
% 
% Output:
% T_mdh [4x4x10]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(10+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-07-01 10:28
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = hybBKspatial_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'hybBKspatial_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'hybBKspatial_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-07-01 10:16:02
% EndTime: 2020-07-01 10:16:06
% DurationCPUTime: 3.54s
% Computational Cost: add. (8489->391), mult. (25202->473), div. (74->5), fcn. (2936->14), ass. (0->215)
t603 = 4 * pkin(3);
t460 = pkin(5) ^ 2;
t403 = -t460 / 0.3e1;
t463 = pkin(4) ^ 2;
t602 = t403 - t463 / 0.3e1;
t459 = t460 ^ 2;
t462 = t463 ^ 2;
t601 = -t459 / 0.6e1 + t462 / 0.6e1;
t431 = cos(qJ(3));
t389 = t431 ^ 2;
t600 = -0.2e1 * t389;
t468 = pkin(3) ^ 2;
t466 = t468 ^ 2;
t599 = 4 * t466;
t598 = 2 * t468;
t446 = 6 * t468;
t472 = (pkin(2) ^ 2);
t597 = 2 * t472;
t474 = (pkin(1) ^ 2);
t457 = 2 * t474;
t479 = t472 ^ 2;
t449 = 5 * t479;
t596 = -pkin(4) - pkin(5);
t595 = -pkin(4) + pkin(5);
t426 = sin(qJ(2));
t381 = pkin(2) * t426;
t432 = cos(qJ(2));
t594 = pkin(2) * t432;
t425 = sin(qJ(3));
t593 = pkin(3) * t425;
t380 = pkin(3) * t431;
t549 = pkin(1) * t594;
t370 = -0.2e1 * t549;
t392 = t432 ^ 2;
t576 = t472 * t392;
t543 = 0.2e1 * t576;
t556 = -t472 + t474;
t350 = t370 + t543 + t556;
t583 = t425 * t426;
t541 = pkin(2) * t583;
t512 = pkin(3) * t541;
t362 = -0.2e1 * t512;
t451 = -3 * t472;
t555 = t472 + t474;
t518 = -t460 + t555;
t501 = t468 + t518;
t494 = -t463 + t501;
t493 = t370 + t494;
t559 = t460 - t474;
t521 = t463 + t559;
t577 = t468 * t389;
t545 = -0.4e1 * t577;
t557 = t468 - t474;
t367 = pkin(1) - t594;
t589 = t367 * t431;
t340 = sqrt(t350 * t545 + 0.4e1 * t557 * t576 + 0.4e1 * t494 * t549 - t466 + (t451 + t521) * t598 - (t474 + (pkin(2) - t595) * (pkin(2) + t595)) * (t474 + (pkin(2) - t596) * (pkin(2) + t596)) + (-(t362 + t493) * t589 + t493 * t541) * t603);
t473 = t474 ^ 2;
t558 = t466 + t473;
t563 = t457 - t460;
t575 = t474 * t460;
t348 = t563 * t468 + t558 - t575 - t601;
t376 = t474 - t472 / 0.3e1;
t351 = t376 * t362;
t411 = t463 / 0.3e1;
t520 = t468 + t555;
t355 = t403 + t411 + t520;
t356 = t367 + t380;
t386 = t451 + t474;
t388 = t431 * t389;
t482 = pkin(3) * t468;
t573 = t482 * t388;
t542 = pkin(1) * t573;
t513 = 0.8e1 * t542;
t359 = t386 * t513;
t456 = 3 * t474;
t560 = t460 + t463;
t371 = t456 - t560;
t588 = t371 * t468;
t361 = 0.10e2 * t588;
t404 = -t460 / 0.2e1;
t364 = t404 + t520;
t383 = t468 + t474;
t519 = -t460 + t383;
t365 = t463 + t519;
t368 = pkin(1) + t380;
t548 = pkin(1) * t380;
t369 = 0.2e1 * t548;
t455 = 4 * t474;
t373 = (t455 + t460) * t468;
t375 = -t468 / 0.3e1 + t474;
t377 = t383 ^ 2;
t378 = 0.10e2 / 0.3e1 * t468;
t379 = -0.30e2 * t460 + (60 * t474);
t385 = -3 * t468 + t474;
t401 = -t460 / 0.6e1;
t402 = -t460 / 0.4e1;
t416 = 0.4e1 / 0.3e1 * t468;
t417 = t468 / 0.3e1;
t418 = t468 / 0.2e1;
t438 = 15 * t466;
t439 = 15 * t468;
t440 = 10 * t468;
t441 = -0.2e1 * t460;
t442 = -0.5e1 * t460;
t443 = 7 * t466;
t444 = 5 * t466;
t445 = 7 * t468;
t447 = 3 * t468;
t450 = 3 * t472;
t453 = 3 * t473;
t454 = 8 * t474;
t478 = pkin(2) * t472;
t469 = t478 ^ 2;
t475 = pkin(1) * t474;
t491 = t348 + t479;
t420 = t472 / 0.2e1;
t566 = t420 + t474;
t492 = -t512 + t566;
t570 = t459 / 0.2e1 - t462 / 0.2e1;
t505 = -0.3e1 * t575 + t453 + t570;
t508 = -0.6e1 * t512;
t406 = -0.3e1 / 0.2e1 * t460;
t567 = t406 + t456;
t572 = t383 * ((t406 + t457) * t468 - 0.3e1 / 0.2e1 * t575 + t558 + t570) + t469;
t496 = ((t378 + t563) * t472 + t491) * t508 + (t438 + (-0.9e1 * t460 + (18 * t474)) * t468 + t505) * t472 + (t439 + t567) * t479 + t572;
t507 = -0.4e1 * t512;
t497 = t364 * t507;
t523 = t447 + t555;
t498 = -(t450 + t383) * t593 + t523 * t381;
t415 = -0.2e1 / 0.3e1 * t463;
t405 = -0.2e1 / 0.3e1 * t460;
t524 = t405 + t383;
t564 = t440 + t457;
t568 = t405 + t415;
t499 = -(t449 + ((5 * t468) + t371) * t597 + (t415 + t524) * t383) * t593 + (t479 + (t564 + t568) * t472 + t444 + 0.2e1 * t588 + t474 * (t474 + t568)) * t381;
t527 = t474 + t602;
t502 = t468 + t527;
t569 = t401 - t463 / 0.6e1;
t528 = t474 + t569;
t503 = t420 + t528;
t562 = t459 - t462;
t504 = -0.6e1 * t575 + (6 * t473) + t562;
t410 = 0.2e1 / 0.3e1 * t463;
t525 = t405 + t410 + t457;
t571 = (t410 + t524) * t383 + t479;
t506 = t355 * t507 + (t446 + t525) * t472 + t571;
t539 = t482 * t381;
t509 = t388 * t539;
t578 = t466 * t389 ^ 2;
t510 = t578 * t381;
t537 = 0.16e2 * t573;
t511 = pkin(1) * t537;
t514 = 0.20e2 / 0.3e1 * t468;
t561 = -t460 + t463;
t522 = t456 + t561;
t526 = t404 - t463 / 0.2e1 + t474;
t529 = t460 / 0.3e1 + t411 + t457;
t530 = 0.2e1 / 0.3e1 * t460 + t410 + t455;
t531 = 0.4e1 / 0.3e1 * t460 + 0.4e1 / 0.3e1 * t463 - (2 * t474);
t532 = 0.6e1 * t548;
t533 = 0.4e1 * t548;
t580 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t534 = t425 * t580;
t536 = -t593 / 0.2e1;
t538 = 0.12e2 * t577;
t540 = t468 * t381;
t544 = 0.4e1 * t577;
t546 = 0.8e1 * t578;
t584 = t432 * t392 * t478;
t547 = 0.8e1 * t584;
t551 = 0.2e1 * t593;
t552 = 0.4e1 * pkin(1);
t553 = t473 + t479;
t554 = t473 - t466;
t565 = 0.4e1 / 0.7e1 * t474 - t460 / 0.7e1;
t574 = t474 * t468;
t581 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t585 = t392 ^ 2 * t479;
t586 = t560 * t474;
t587 = t377 * (-t463 + t519);
t590 = (-t425 * t482 + t540) * t389;
t592 = ((-0.24e2 * (0.4e1 / 0.3e1 * t577 + t369 + t375) * t585 * t593 - 0.12e2 * (-0.8e1 / 0.3e1 * t510 + ((t416 + t503) * t381 - (0.4e1 / 0.3e1 * t472 + t418 + t528) * t593) * t544 + (-(t472 * t557) - 0.5e1 / 0.3e1 * t466 + t529 * t468 + t474 * t527) * t381 + (-t479 + (-t514 + t530) * t472 - (3 * t466) + t531 * t468 + t473) * t536 + (-t425 * t466 * t388 + ((t468 + t503) * t381 + (t597 - t557) * t536) * t380) * t552) * t576 + 0.24e2 * t376 * t510 + ((t450 + 0.3e1 / 0.2e1 * t468 + t526) * t381 + t386 * t593 / 0.2e1) * t511 - 0.6e1 * ((-(3 * t479) + (-t514 + t531) * t472 + t530 * t468 + t554) * t381 - 0.2e1 * (-0.5e1 / 0.3e1 * t479 + (-t468 + t529) * t472 + t474 * t502) * t593) * t577 - 0.6e1 * t499 * t548 - (t469 + ((21 * t468) + t371) * t479 + (t361 + t453 + (35 * t466) - 0.2e1 * t586) * t472 + (t443 + (t442 + t454 - 0.5e1 * t463) * t468 - t474 * t521) * t383) * t381 + (0.7e1 * t469 + (t445 + t371) * t449 + (t361 + (21 * t466) + (9 * t473) - 0.6e1 * t586) * t472 + t587) * t593) * t340 + (0.16e2 * (t546 + t511 + (-8 * t466 + 12 * t574) * t389 + (-0.12e2 * pkin(1) * t482 + t475 * t603) * t431 - (6 * t574) + t558) * t585 + 0.24e2 * ((t474 - 0.2e1 / 0.3e1 * t472) * t546 + 0.14e2 * (-0.32e2 / 0.21e2 * (t474 + t472 / 0.4e1 + t468 / 0.4e1 - t460 / 0.8e1) * t512 + t479 / 0.7e1 + (0.16e2 / 0.21e2 * t468 + t565) * t472 + t466 / 0.7e1 + t565 * t468 + t473 - 0.3e1 / 0.7e1 * t575 + t459 / 0.42e2 - t462 / 0.42e2) * t577 + t375 * t497 - (t557 * t479) + (t373 - 0.10e2 / 0.3e1 * t466 + (2 * t473) - t575) * t472 + t348 * t581 + ((-0.2e1 / 0.3e1 * t512 + t474 + t418 + t402) * t537 + 0.6e1 * (-0.8e1 / 0.3e1 * (t402 + t417 + t566) * t512 + t479 / 0.3e1 + (0.4e1 / 0.3e1 * t474 + t416 + t403) * t472 + t473 + 0.2e1 / 0.3e1 * t574 - 0.2e1 / 0.3e1 * t575 - t466 / 0.3e1 + t459 / 0.18e2 - t462 / 0.18e2) * t380) * pkin(1)) * t576 + 0.16e2 * (-6 * t474 * t472 + t553) * t578 + 0.32e2 * (t362 * t580 + t364 * t386) * t542 + 0.24e2 * (t376 * t497 - t469 + (-t378 + t559) * t479 + (t373 + t554 + t601) * t472 + t348 * t474) * t577 + 0.8e1 * t496 * t548 - 0.8e1 * ((t445 + t567) * t479 + (t443 + (t442 + (10 * t474)) * t468 + t505) * t472 + t572) * t512 + (t479 ^ 2) + (t441 + t455 + (28 * t468)) * t469 + (t379 * t468 + (70 * t466) + t504) * t479 + (t379 * t466 + t504 * t446 + t562 * t457 - 0.6e1 * t473 * t460 + 0.4e1 * t475 ^ 2 + (28 * t482 ^ 2)) * t472 + t365 * t587) * t356 + (((0.4e1 * t590 + (t381 + t551) * t369 + t385 * t381 + (0.3e1 / 0.2e1 * t472 + t447 + t526) * t551) * t547 + 0.6e1 * ((0.2e1 * (t418 + t472 + t569) * t381 + pkin(3) * t534) * t545 + (-0.8e1 * t509 + 0.4e1 * ((t523 + t602) * t381 - (t450 + t502) * t593) * t380) * pkin(1) + t499) * t594) * t340 + (-0.32e2 * (t513 + (-0.4e1 * t425 * t539 + t599 + ((4 * t472) + t441 + t454) * t468) * t389 + (t402 - t468 + t492) * t533 + t362 * t581 + t364 * t385) * t584 - 0.8e1 * (t359 + (t364 * t580 + t351) * t538 + (t497 + (t446 + t563) * t472 + t491) * t532 + t496) * t594) * t356) * t368) / ((-0.4e1 * (0.2e1 * t590 + (t598 + t472) * t593 + (-t557 + t369) * t381) * t576 + 0.8e1 * pkin(1) * t509 + ((pkin(2) * t599 + 0.8e1 * t468 * t478) * t426 + 0.4e1 * t482 * t534) * t389 - 0.4e1 * t498 * t548 - (t564 * t472 + t444 + t553 + 6 * t574) * t381 + (t449 + (t440 + 6 * t474) * t472 + t377) * t593) * t340 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t512 + t474 + t472 / 0.3e1 + t417 + t463 / 0.9e1 - t460 / 0.9e1) * t577 + t375 * t362 + t355 * t581 + (t573 + (t463 / 0.6e1 + t401 + t492) * t380) * t552) * t576 + t359 + (t355 * t580 + t351) * t538 + t506 * t532 + ((t378 + t525) * t472 + t571) * t508 + t469 + (t439 + t522) * t479 + (t522 * t446 + t561 * t457 + t438 + t453) * t472 + t377 * t365) * t356 + ((t547 * t593 + 0.4e1 * (t540 * t600 + (t381 - t593) * t369 + t498) * t594) * t340 + (-0.8e1 * (t369 + t544 + t385) * t584 - 0.6e1 * (t556 * t544 + (t362 + t355) * t533 + t506) * t594) * t356) * t368);
t550 = 0.2e1 * t380;
t343 = 0.1e1 / (t367 * t550 + t362 + t370 + t520);
t464 = 0.1e1 / pkin(4);
t591 = t343 * t464;
t582 = t426 * t431;
t461 = 0.1e1 / pkin(5);
t579 = t461 * t464;
t535 = t343 * t461 / pkin(4) ^ 2;
t517 = t592 / 0.4e1;
t516 = -t591 / 0.2e1;
t515 = t579 / 0.2e1;
t500 = -pkin(3) + t541;
t495 = t447 + t463 + t518;
t433 = cos(qJ(1));
t430 = cos(qJ(4));
t429 = cos(qJ(5));
t428 = cos(qJ(6));
t427 = sin(qJ(1));
t424 = sin(qJ(4));
t423 = sin(qJ(5));
t422 = sin(qJ(6));
t354 = t431 * t432 + t583;
t353 = t425 * t432 - t582;
t349 = t370 + t463 + t501;
t339 = (-t500 + t589) * t340 + (t349 * t367 + t350 * t550) * t425 + (t349 * t431 + (0.4e1 * t389 - 0.2e1) * t367 * pkin(3)) * t381;
t338 = (pkin(2) * t582 + t367 * t425) * t340 - (t370 + t495 + t507) * t589 + t500 * t370 + t495 * t541 + (t350 * t600 - t365 - t450 + t543) * pkin(3);
t337 = t338 * t516;
t334 = t515 * t592;
t333 = (t338 * t340 / 0.4e1 + t339 * t517) * t535;
t332 = (t338 * t517 - t339 * t340 / 0.4e1) * t535;
t331 = -t332 * t354 + t333 * t353;
t330 = t332 * t353 + t333 * t354;
t1 = [t433, -t427, 0, 0; t427, t433, 0, 0; 0, 0, 1, pkin(6); t432, -t426, 0, -pkin(1) / 0.2e1; 0, 0, -1, 0; t426, t432, 0, 0; t331, t330, 0, pkin(2); -t330, t331, 0, 0; 0, 0, 1, 0; t431, -t425, 0, pkin(1) / 0.2e1; 0, 0, -1, 0; t425, t431, 0, 0; t337, t339 * t591 / 0.2e1, 0, pkin(3); t339 * t516, t337, 0, 0; 0, 0, 1, 0; t430, -t424, 0, pkin(7); t424, t430, 0, 0; 0, 0, 1, 0; t429, -t423, 0, 0; 0, 0, 1, 0; -t423, -t429, 0, 0; 0, 0, -1, 0; -t428, t422, 0, 0; t422, t428, 0, 0; t334, -t340 * t579 / 0.2e1, 0, -pkin(4); t340 * t515, t334, 0, 0; 0, 0, 1, 0; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,10);             % numerisch
else,                         T_mdh = sym('xx', [4,4,10]); end % symbolisch

for i = 1:10
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end