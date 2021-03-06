% Calculate homogenous joint transformation matrices for
% fivebar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = fivebar1DE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:21:30
% EndTime: 2020-04-27 02:21:34
% DurationCPUTime: 3.95s
% Computational Cost: add. (16969->389), mult. (50400->465), div. (148->5), fcn. (5856->15), ass. (0->210)
t551 = 4 * pkin(3);
t410 = pkin(5) ^ 2;
t361 = -t410 / 0.3e1;
t413 = pkin(4) ^ 2;
t550 = t361 - t413 / 0.3e1;
t409 = t410 ^ 2;
t412 = t413 ^ 2;
t549 = -t409 / 0.6e1 + t412 / 0.6e1;
t382 = cos(qJ(2));
t347 = t382 ^ 2;
t548 = -0.2e1 * t347;
t418 = pkin(3) ^ 2;
t416 = t418 ^ 2;
t547 = 4 * t416;
t546 = 2 * t418;
t396 = 6 * t418;
t422 = (pkin(2) ^ 2);
t545 = 2 * t422;
t424 = (pkin(1) ^ 2);
t407 = 2 * t424;
t429 = t422 ^ 2;
t399 = 5 * t429;
t544 = -pkin(4) - pkin(5);
t543 = -pkin(4) + pkin(5);
t542 = 0.1e1 / pkin(4) / 0.2e1;
t381 = sin(qJ(1));
t339 = pkin(2) * t381;
t383 = cos(qJ(1));
t541 = pkin(2) * t383;
t380 = sin(qJ(2));
t540 = pkin(3) * t380;
t338 = pkin(3) * t382;
t498 = pkin(1) * t541;
t328 = -0.2e1 * t498;
t350 = t383 ^ 2;
t531 = t350 * t422;
t495 = 0.2e1 * t531;
t505 = -t422 + t424;
t308 = t328 + t495 + t505;
t530 = t380 * t381;
t490 = pkin(2) * t530;
t462 = pkin(3) * t490;
t320 = -0.2e1 * t462;
t401 = -3 * t422;
t341 = t418 + t424;
t468 = -t410 + t341;
t324 = -t413 + t468;
t444 = t422 + t324;
t443 = t328 + t444;
t508 = t410 - t424;
t470 = t413 + t508;
t525 = t418 * t347;
t493 = -0.4e1 * t525;
t506 = t418 - t424;
t325 = pkin(1) - t541;
t537 = t325 * t382;
t298 = sqrt(t308 * t493 + 0.4e1 * t506 * t531 + 0.4e1 * t444 * t498 - t416 + (t401 + t470) * t546 - (t424 + (pkin(2) - t543) * (pkin(2) + t543)) * (t424 + (pkin(2) - t544) * (pkin(2) + t544)) + (-(t320 + t443) * t537 + t443 * t490) * t551);
t423 = t424 ^ 2;
t507 = t416 + t423;
t512 = t407 - t410;
t524 = t424 * t410;
t306 = t512 * t418 + t507 - t524 - t549;
t334 = t424 - t422 / 0.3e1;
t309 = t334 * t320;
t369 = t413 / 0.3e1;
t504 = t422 + t424;
t469 = t418 + t504;
t313 = t361 + t369 + t469;
t314 = t325 + t338;
t344 = t401 + t424;
t346 = t382 * t347;
t432 = pkin(3) * t418;
t522 = t432 * t346;
t491 = pkin(1) * t522;
t463 = 0.8e1 * t491;
t317 = t344 * t463;
t406 = 3 * t424;
t509 = t410 + t413;
t329 = t406 - t509;
t536 = t329 * t418;
t319 = 0.10e2 * t536;
t362 = -t410 / 0.2e1;
t322 = t362 + t469;
t323 = t413 + t468;
t326 = pkin(1) + t338;
t497 = pkin(1) * t338;
t327 = 0.2e1 * t497;
t405 = 4 * t424;
t331 = (t405 + t410) * t418;
t333 = -t418 / 0.3e1 + t424;
t335 = t341 ^ 2;
t336 = 0.10e2 / 0.3e1 * t418;
t337 = -0.30e2 * t410 + (60 * t424);
t343 = -3 * t418 + t424;
t359 = -t410 / 0.6e1;
t360 = -t410 / 0.4e1;
t374 = 0.4e1 / 0.3e1 * t418;
t375 = t418 / 0.3e1;
t376 = t418 / 0.2e1;
t388 = 15 * t416;
t389 = 15 * t418;
t390 = 10 * t418;
t391 = -0.2e1 * t410;
t392 = -0.5e1 * t410;
t393 = 7 * t416;
t394 = 5 * t416;
t395 = 7 * t418;
t397 = 3 * t418;
t400 = 3 * t422;
t403 = 3 * t423;
t404 = 8 * t424;
t428 = pkin(2) * t422;
t419 = t428 ^ 2;
t425 = pkin(1) * t424;
t441 = t306 + t429;
t378 = t422 / 0.2e1;
t515 = t378 + t424;
t442 = -t462 + t515;
t519 = t409 / 0.2e1 - t412 / 0.2e1;
t455 = -0.3e1 * t524 + t403 + t519;
t458 = -0.6e1 * t462;
t364 = -0.3e1 / 0.2e1 * t410;
t516 = t364 + t406;
t521 = t341 * ((t364 + t407) * t418 - 0.3e1 / 0.2e1 * t524 + t507 + t519) + t419;
t446 = ((t336 + t512) * t422 + t441) * t458 + (t388 + (-0.9e1 * t410 + (18 * t424)) * t418 + t455) * t422 + (t389 + t516) * t429 + t521;
t457 = -0.4e1 * t462;
t447 = t322 * t457;
t472 = t397 + t504;
t448 = -(t400 + t341) * t540 + t472 * t339;
t373 = -0.2e1 / 0.3e1 * t413;
t363 = -0.2e1 / 0.3e1 * t410;
t473 = t363 + t341;
t513 = t390 + t407;
t517 = t363 + t373;
t449 = -(t399 + ((5 * t418) + t329) * t545 + t341 * (t373 + t473)) * t540 + (t429 + (t513 + t517) * t422 + t394 + 0.2e1 * t536 + t424 * (t424 + t517)) * t339;
t476 = t424 + t550;
t452 = t418 + t476;
t518 = t359 - t413 / 0.6e1;
t477 = t424 + t518;
t453 = t378 + t477;
t511 = t409 - t412;
t454 = -0.6e1 * t524 + (6 * t423) + t511;
t368 = 0.2e1 / 0.3e1 * t413;
t474 = t363 + t368 + t407;
t520 = t341 * (t368 + t473) + t429;
t456 = t313 * t457 + (t396 + t474) * t422 + t520;
t488 = t432 * t339;
t459 = t346 * t488;
t526 = t416 * t347 ^ 2;
t460 = t526 * t339;
t486 = 0.16e2 * t522;
t461 = pkin(1) * t486;
t464 = 0.20e2 / 0.3e1 * t418;
t510 = -t410 + t413;
t471 = t406 + t510;
t475 = t362 - t413 / 0.2e1 + t424;
t478 = t410 / 0.3e1 + t369 + t407;
t479 = 0.2e1 / 0.3e1 * t410 + t368 + t405;
t480 = 0.4e1 / 0.3e1 * t410 + 0.4e1 / 0.3e1 * t413 - (2 * t424);
t481 = 0.6e1 * t497;
t482 = 0.4e1 * t497;
t527 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t483 = t380 * t527;
t485 = -t540 / 0.2e1;
t487 = 0.12e2 * t525;
t489 = t418 * t339;
t492 = 0.4e1 * t525;
t494 = 0.8e1 * t526;
t532 = t383 * t350 * t428;
t496 = 0.8e1 * t532;
t500 = 0.2e1 * t540;
t501 = 0.4e1 * pkin(1);
t502 = t423 + t429;
t503 = t423 - t416;
t514 = 0.4e1 / 0.7e1 * t424 - t410 / 0.7e1;
t523 = t424 * t418;
t528 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t533 = t350 ^ 2 * t429;
t534 = t509 * t424;
t535 = t335 * t324;
t538 = (-t380 * t432 + t489) * t347;
t539 = ((-0.24e2 * (0.4e1 / 0.3e1 * t525 + t327 + t333) * t533 * t540 - 0.12e2 * (-0.8e1 / 0.3e1 * t460 + ((t374 + t453) * t339 - (0.4e1 / 0.3e1 * t422 + t376 + t477) * t540) * t492 + (-(t422 * t506) - 0.5e1 / 0.3e1 * t416 + t478 * t418 + t424 * t476) * t339 + (-t429 + (-t464 + t479) * t422 - (3 * t416) + t480 * t418 + t423) * t485 + (-t380 * t416 * t346 + ((t418 + t453) * t339 + (t545 - t506) * t485) * t338) * t501) * t531 + 0.24e2 * t334 * t460 + ((t400 + 0.3e1 / 0.2e1 * t418 + t475) * t339 + t344 * t540 / 0.2e1) * t461 - 0.6e1 * ((-(3 * t429) + (-t464 + t480) * t422 + t479 * t418 + t503) * t339 - 0.2e1 * (-0.5e1 / 0.3e1 * t429 + (-t418 + t478) * t422 + t424 * t452) * t540) * t525 - 0.6e1 * t449 * t497 - (t419 + ((21 * t418) + t329) * t429 + (t319 + t403 + (35 * t416) - 0.2e1 * t534) * t422 + (t393 + (t392 + t404 - 0.5e1 * t413) * t418 - t424 * t470) * t341) * t339 + (0.7e1 * t419 + (t395 + t329) * t399 + (t319 + (21 * t416) + (9 * t423) - 0.6e1 * t534) * t422 + t535) * t540) * t298 + t314 * (0.16e2 * (t494 + t461 + (-8 * t416 + 12 * t523) * t347 + (-0.12e2 * pkin(1) * t432 + t425 * t551) * t382 - (6 * t523) + t507) * t533 + 0.24e2 * ((t424 - 0.2e1 / 0.3e1 * t422) * t494 + 0.14e2 * (-0.32e2 / 0.21e2 * (t424 + t422 / 0.4e1 + t418 / 0.4e1 - t410 / 0.8e1) * t462 + t429 / 0.7e1 + (0.16e2 / 0.21e2 * t418 + t514) * t422 + t416 / 0.7e1 + t514 * t418 + t423 - 0.3e1 / 0.7e1 * t524 + t409 / 0.42e2 - t412 / 0.42e2) * t525 + t333 * t447 - (t506 * t429) + (t331 - 0.10e2 / 0.3e1 * t416 + (2 * t423) - t524) * t422 + t306 * t528 + ((-0.2e1 / 0.3e1 * t462 + t424 + t376 + t360) * t486 + 0.6e1 * (-0.8e1 / 0.3e1 * (t360 + t375 + t515) * t462 + t429 / 0.3e1 + (0.4e1 / 0.3e1 * t424 + t374 + t361) * t422 + t423 + 0.2e1 / 0.3e1 * t523 - 0.2e1 / 0.3e1 * t524 - t416 / 0.3e1 + t409 / 0.18e2 - t412 / 0.18e2) * t338) * pkin(1)) * t531 + 0.16e2 * (-6 * t424 * t422 + t502) * t526 + 0.32e2 * (t320 * t527 + t322 * t344) * t491 + 0.24e2 * (t334 * t447 - t419 + (-t336 + t508) * t429 + (t331 + t503 + t549) * t422 + t424 * t306) * t525 + 0.8e1 * t446 * t497 - 0.8e1 * ((t395 + t516) * t429 + (t393 + (t392 + (10 * t424)) * t418 + t455) * t422 + t521) * t462 + (t429 ^ 2) + (t391 + t405 + (28 * t418)) * t419 + (t337 * t418 + (70 * t416) + t454) * t429 + (t337 * t416 + t454 * t396 + t511 * t407 - 0.6e1 * t423 * t410 + 0.4e1 * t425 ^ 2 + (28 * t432 ^ 2)) * t422 + t323 * t535) + (((0.4e1 * t538 + (t339 + t500) * t327 + t343 * t339 + (0.3e1 / 0.2e1 * t422 + t397 + t475) * t500) * t496 + 0.6e1 * ((0.2e1 * (t376 + t422 + t518) * t339 + pkin(3) * t483) * t493 + (-0.8e1 * t459 + 0.4e1 * ((t472 + t550) * t339 - (t400 + t452) * t540) * t338) * pkin(1) + t449) * t541) * t298 + t314 * (-0.32e2 * (t463 + (-0.4e1 * t380 * t488 + t547 + ((4 * t422) + t391 + t404) * t418) * t347 + (t360 - t418 + t442) * t482 + t320 * t528 + t322 * t343) * t532 - 0.8e1 * (t317 + (t322 * t527 + t309) * t487 + (t447 + (t396 + t512) * t422 + t441) * t481 + t446) * t541)) * t326) / ((-0.4e1 * (0.2e1 * t538 + (t546 + t422) * t540 + (-t506 + t327) * t339) * t531 + 0.8e1 * pkin(1) * t459 + ((pkin(2) * t547 + 0.8e1 * t418 * t428) * t381 + 0.4e1 * t432 * t483) * t347 - 0.4e1 * t448 * t497 - (t513 * t422 + t394 + t502 + 6 * t523) * t339 + (t399 + (t390 + 6 * t424) * t422 + t335) * t540) * t298 + t314 * (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t462 + t424 + t422 / 0.3e1 + t375 + t413 / 0.9e1 - t410 / 0.9e1) * t525 + t333 * t320 + t313 * t528 + (t522 + (t413 / 0.6e1 + t359 + t442) * t338) * t501) * t531 + t317 + (t313 * t527 + t309) * t487 + t456 * t481 + ((t336 + t474) * t422 + t520) * t458 + t419 + (t389 + t471) * t429 + (t471 * t396 + t510 * t407 + t388 + t403) * t422 + t335 * t323) + ((t496 * t540 + 0.4e1 * (t489 * t548 + (t339 - t540) * t327 + t448) * t541) * t298 + t314 * (-0.8e1 * (t327 + t492 + t343) * t532 - 0.6e1 * (t505 * t492 + (t320 + t313) * t482 + t456) * t541)) * t326);
t529 = t381 * t382;
t499 = 0.2e1 * t338;
t301 = 0.1e1 / (t325 * t499 + t320 + t328 + t469);
t411 = 0.1e1 / pkin(5);
t484 = t301 * t411 / pkin(4) ^ 2;
t467 = t539 / 0.4e1;
t466 = t301 * t542;
t465 = t411 * t542;
t451 = t504 + t510;
t450 = -pkin(3) + t490;
t445 = t397 + t451;
t312 = t382 * t383 + t530;
t311 = t380 * t383 - t529;
t307 = t328 + t418 + t451;
t297 = (-t450 + t537) * t298 + (t307 * t325 + t308 * t499) * t380 + (t307 * t382 + (0.4e1 * t347 - 0.2e1) * t325 * pkin(3)) * t339;
t296 = (pkin(2) * t529 + t325 * t380) * t298 - (t328 + t445 + t457) * t537 + t450 * t328 + t445 * t490 + (t308 * t548 - t323 - t400 + t495) * pkin(3);
t295 = atan2(t297 * t466, t296 * t466);
t294 = cos(t295);
t293 = sin(t295);
t290 = atan2(t298 * t465, t465 * t539);
t289 = cos(t290);
t288 = sin(t290);
t287 = (t296 * t298 / 0.4e1 + t297 * t467) * t484;
t286 = (t296 * t467 - t297 * t298 / 0.4e1) * t484;
t285 = atan2(t286 * t311 + t287 * t312, t286 * t312 - t287 * t311);
t284 = cos(t285);
t283 = sin(t285);
t1 = [t383, -t381, 0, 0; t381, t383, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t284, t283, 0, pkin(2); -t283, -t284, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t382, -t380, 0, pkin(1); t380, t382, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t294, -t293, 0, pkin(3); t293, t294, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t289, t288, 0, pkin(4); -t288, -t289, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
