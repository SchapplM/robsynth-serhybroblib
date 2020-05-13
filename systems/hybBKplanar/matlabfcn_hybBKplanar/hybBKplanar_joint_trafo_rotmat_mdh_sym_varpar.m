% Calculate homogenous joint transformation matrices for
% hybBKplanar
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,CF,ED]';
% 
% Output:
% T_mdh [4x4x7]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 19:03
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = hybBKplanar_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'hybBKplanar_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'hybBKplanar_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-11 18:53:37
% EndTime: 2020-04-11 18:53:40
% DurationCPUTime: 3.30s
% Computational Cost: add. (8484->386), mult. (25200->471), div. (74->5), fcn. (2924->8), ass. (0->209)
t555 = 4 * pkin(3);
t412 = pkin(6) ^ 2;
t361 = -t412 / 0.3e1;
t415 = pkin(4) ^ 2;
t554 = t361 - t415 / 0.3e1;
t411 = t412 ^ 2;
t414 = t415 ^ 2;
t553 = -t411 / 0.6e1 + t414 / 0.6e1;
t384 = cos(qJ(2));
t347 = t384 ^ 2;
t552 = -0.2e1 * t347;
t420 = pkin(3) ^ 2;
t418 = t420 ^ 2;
t551 = 4 * t418;
t550 = 2 * t420;
t398 = 6 * t420;
t424 = (pkin(2) ^ 2);
t549 = 2 * t424;
t426 = (pkin(1) ^ 2);
t409 = 2 * t426;
t431 = t424 ^ 2;
t401 = 5 * t431;
t548 = -pkin(4) - pkin(6);
t547 = -pkin(4) + pkin(6);
t382 = sin(qJ(1));
t339 = pkin(2) * t382;
t385 = cos(qJ(1));
t546 = pkin(2) * t385;
t381 = sin(qJ(2));
t545 = pkin(3) * t381;
t338 = pkin(3) * t384;
t501 = pkin(1) * t546;
t328 = -0.2e1 * t501;
t350 = t385 ^ 2;
t535 = t350 * t424;
t498 = 0.2e1 * t535;
t508 = -t424 + t426;
t308 = t328 + t498 + t508;
t534 = t381 * t382;
t493 = pkin(2) * t534;
t464 = pkin(3) * t493;
t320 = -0.2e1 * t464;
t403 = -3 * t424;
t507 = t424 + t426;
t470 = -t412 + t507;
t453 = t420 + t470;
t446 = -t415 + t453;
t445 = t328 + t446;
t511 = t412 - t426;
t473 = t415 + t511;
t528 = t420 * t347;
t496 = -0.4e1 * t528;
t509 = t420 - t426;
t325 = pkin(1) - t546;
t541 = t325 * t384;
t298 = sqrt(t308 * t496 + 0.4e1 * t509 * t535 + 0.4e1 * t446 * t501 - t418 + (t403 + t473) * t550 - (t426 + (pkin(2) - t547) * (pkin(2) + t547)) * (t426 + (pkin(2) - t548) * (pkin(2) + t548)) + (-(t320 + t445) * t541 + t445 * t493) * t555);
t425 = t426 ^ 2;
t510 = t418 + t425;
t515 = t409 - t412;
t527 = t426 * t412;
t306 = t515 * t420 + t510 - t527 - t553;
t334 = t426 - t424 / 0.3e1;
t309 = t334 * t320;
t369 = t415 / 0.3e1;
t472 = t420 + t507;
t313 = t361 + t369 + t472;
t314 = t325 + t338;
t344 = t403 + t426;
t346 = t384 * t347;
t434 = pkin(3) * t420;
t525 = t434 * t346;
t494 = pkin(1) * t525;
t465 = 0.8e1 * t494;
t317 = t344 * t465;
t408 = 3 * t426;
t512 = t412 + t415;
t329 = t408 - t512;
t540 = t329 * t420;
t319 = 0.10e2 * t540;
t362 = -t412 / 0.2e1;
t322 = t362 + t472;
t341 = t420 + t426;
t471 = -t412 + t341;
t323 = t415 + t471;
t326 = pkin(1) + t338;
t500 = pkin(1) * t338;
t327 = 0.2e1 * t500;
t407 = 4 * t426;
t331 = (t407 + t412) * t420;
t333 = -t420 / 0.3e1 + t426;
t335 = t341 ^ 2;
t336 = 0.10e2 / 0.3e1 * t420;
t337 = -0.30e2 * t412 + (60 * t426);
t343 = -3 * t420 + t426;
t359 = -t412 / 0.6e1;
t360 = -t412 / 0.4e1;
t374 = 0.4e1 / 0.3e1 * t420;
t375 = t420 / 0.3e1;
t376 = t420 / 0.2e1;
t390 = 15 * t418;
t391 = 15 * t420;
t392 = 10 * t420;
t393 = -0.2e1 * t412;
t394 = -0.5e1 * t412;
t395 = 7 * t418;
t396 = 5 * t418;
t397 = 7 * t420;
t399 = 3 * t420;
t402 = 3 * t424;
t405 = 3 * t425;
t406 = 8 * t426;
t430 = pkin(2) * t424;
t421 = t430 ^ 2;
t427 = pkin(1) * t426;
t443 = t306 + t431;
t378 = t424 / 0.2e1;
t518 = t378 + t426;
t444 = -t464 + t518;
t522 = t411 / 0.2e1 - t414 / 0.2e1;
t457 = -0.3e1 * t527 + t405 + t522;
t460 = -0.6e1 * t464;
t364 = -0.3e1 / 0.2e1 * t412;
t519 = t364 + t408;
t524 = ((t364 + t409) * t420 - 0.3e1 / 0.2e1 * t527 + t510 + t522) * t341 + t421;
t448 = ((t336 + t515) * t424 + t443) * t460 + (t390 + (-0.9e1 * t412 + (18 * t426)) * t420 + t457) * t424 + (t391 + t519) * t431 + t524;
t459 = -0.4e1 * t464;
t449 = t322 * t459;
t475 = t399 + t507;
t450 = -(t402 + t341) * t545 + t475 * t339;
t373 = -0.2e1 / 0.3e1 * t415;
t363 = -0.2e1 / 0.3e1 * t412;
t476 = t363 + t341;
t516 = t392 + t409;
t520 = t363 + t373;
t451 = -(t401 + ((5 * t420) + t329) * t549 + (t373 + t476) * t341) * t545 + (t431 + (t516 + t520) * t424 + t396 + 0.2e1 * t540 + t426 * (t426 + t520)) * t339;
t479 = t426 + t554;
t454 = t420 + t479;
t521 = t359 - t415 / 0.6e1;
t480 = t426 + t521;
t455 = t378 + t480;
t514 = t411 - t414;
t456 = -0.6e1 * t527 + (6 * t425) + t514;
t368 = 0.2e1 / 0.3e1 * t415;
t477 = t363 + t368 + t409;
t523 = (t368 + t476) * t341 + t431;
t458 = t313 * t459 + (t398 + t477) * t424 + t523;
t491 = t434 * t339;
t461 = t346 * t491;
t529 = t418 * t347 ^ 2;
t462 = t529 * t339;
t489 = 0.16e2 * t525;
t463 = pkin(1) * t489;
t466 = 0.20e2 / 0.3e1 * t420;
t513 = -t412 + t415;
t474 = t408 + t513;
t478 = t362 - t415 / 0.2e1 + t426;
t481 = t412 / 0.3e1 + t369 + t409;
t482 = 0.2e1 / 0.3e1 * t412 + t368 + t407;
t483 = 0.4e1 / 0.3e1 * t412 + 0.4e1 / 0.3e1 * t415 - (2 * t426);
t484 = 0.6e1 * t500;
t485 = 0.4e1 * t500;
t531 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t486 = t381 * t531;
t488 = -t545 / 0.2e1;
t490 = 0.12e2 * t528;
t492 = t420 * t339;
t495 = 0.4e1 * t528;
t497 = 0.8e1 * t529;
t536 = t385 * t350 * t430;
t499 = 0.8e1 * t536;
t503 = 0.2e1 * t545;
t504 = 0.4e1 * pkin(1);
t505 = t425 + t431;
t506 = t425 - t418;
t517 = 0.4e1 / 0.7e1 * t426 - t412 / 0.7e1;
t526 = t426 * t420;
t532 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t537 = t350 ^ 2 * t431;
t538 = t512 * t426;
t539 = t335 * (-t415 + t471);
t542 = (-t381 * t434 + t492) * t347;
t544 = ((-0.24e2 * (0.4e1 / 0.3e1 * t528 + t327 + t333) * t537 * t545 - 0.12e2 * (-0.8e1 / 0.3e1 * t462 + ((t374 + t455) * t339 - (0.4e1 / 0.3e1 * t424 + t376 + t480) * t545) * t495 + (-(t424 * t509) - 0.5e1 / 0.3e1 * t418 + t481 * t420 + t426 * t479) * t339 + (-t431 + (-t466 + t482) * t424 - (3 * t418) + t483 * t420 + t425) * t488 + (-t381 * t418 * t346 + ((t420 + t455) * t339 + (t549 - t509) * t488) * t338) * t504) * t535 + 0.24e2 * t334 * t462 + ((t402 + 0.3e1 / 0.2e1 * t420 + t478) * t339 + t344 * t545 / 0.2e1) * t463 - 0.6e1 * ((-(3 * t431) + (-t466 + t483) * t424 + t482 * t420 + t506) * t339 - 0.2e1 * (-0.5e1 / 0.3e1 * t431 + (-t420 + t481) * t424 + t426 * t454) * t545) * t528 - 0.6e1 * t451 * t500 - (t421 + ((21 * t420) + t329) * t431 + (t319 + t405 + (35 * t418) - 0.2e1 * t538) * t424 + (t395 + (t394 + t406 - 0.5e1 * t415) * t420 - t426 * t473) * t341) * t339 + (0.7e1 * t421 + (t397 + t329) * t401 + (t319 + (21 * t418) + (9 * t425) - 0.6e1 * t538) * t424 + t539) * t545) * t298 + t314 * (0.16e2 * (t497 + t463 + (-8 * t418 + 12 * t526) * t347 + (-0.12e2 * pkin(1) * t434 + t427 * t555) * t384 - (6 * t526) + t510) * t537 + 0.24e2 * ((t426 - 0.2e1 / 0.3e1 * t424) * t497 + 0.14e2 * (-0.32e2 / 0.21e2 * (t426 + t424 / 0.4e1 + t420 / 0.4e1 - t412 / 0.8e1) * t464 + t431 / 0.7e1 + (0.16e2 / 0.21e2 * t420 + t517) * t424 + t418 / 0.7e1 + t517 * t420 + t425 - 0.3e1 / 0.7e1 * t527 + t411 / 0.42e2 - t414 / 0.42e2) * t528 + t333 * t449 - (t509 * t431) + (t331 - 0.10e2 / 0.3e1 * t418 + (2 * t425) - t527) * t424 + t306 * t532 + ((-0.2e1 / 0.3e1 * t464 + t426 + t376 + t360) * t489 + 0.6e1 * (-0.8e1 / 0.3e1 * (t360 + t375 + t518) * t464 + t431 / 0.3e1 + (0.4e1 / 0.3e1 * t426 + t374 + t361) * t424 + t425 + 0.2e1 / 0.3e1 * t526 - 0.2e1 / 0.3e1 * t527 - t418 / 0.3e1 + t411 / 0.18e2 - t414 / 0.18e2) * t338) * pkin(1)) * t535 + 0.16e2 * (-6 * t426 * t424 + t505) * t529 + 0.32e2 * (t320 * t531 + t322 * t344) * t494 + 0.24e2 * (t334 * t449 - t421 + (-t336 + t511) * t431 + (t331 + t506 + t553) * t424 + t426 * t306) * t528 + 0.8e1 * t448 * t500 - 0.8e1 * ((t397 + t519) * t431 + (t395 + (t394 + (10 * t426)) * t420 + t457) * t424 + t524) * t464 + (t431 ^ 2) + (t393 + t407 + (28 * t420)) * t421 + (t337 * t420 + (70 * t418) + t456) * t431 + (t337 * t418 + t456 * t398 + t514 * t409 - 0.6e1 * t425 * t412 + 0.4e1 * t427 ^ 2 + (28 * t434 ^ 2)) * t424 + t323 * t539) + (((0.4e1 * t542 + (t339 + t503) * t327 + t343 * t339 + (0.3e1 / 0.2e1 * t424 + t399 + t478) * t503) * t499 + 0.6e1 * ((0.2e1 * (t376 + t424 + t521) * t339 + pkin(3) * t486) * t496 + (-0.8e1 * t461 + 0.4e1 * ((t475 + t554) * t339 - (t402 + t454) * t545) * t338) * pkin(1) + t451) * t546) * t298 + t314 * (-0.32e2 * (t465 + (-0.4e1 * t381 * t491 + t551 + ((4 * t424) + t393 + t406) * t420) * t347 + (t360 - t420 + t444) * t485 + t320 * t532 + t322 * t343) * t536 - 0.8e1 * (t317 + (t322 * t531 + t309) * t490 + (t449 + (t398 + t515) * t424 + t443) * t484 + t448) * t546)) * t326) / ((-0.4e1 * (0.2e1 * t542 + (t550 + t424) * t545 + (-t509 + t327) * t339) * t535 + 0.8e1 * pkin(1) * t461 + ((pkin(2) * t551 + 0.8e1 * t420 * t430) * t382 + 0.4e1 * t434 * t486) * t347 - 0.4e1 * t450 * t500 - (t516 * t424 + t396 + t505 + 6 * t526) * t339 + (t401 + (t392 + 6 * t426) * t424 + t335) * t545) * t298 + t314 * (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t464 + t426 + t424 / 0.3e1 + t375 + t415 / 0.9e1 - t412 / 0.9e1) * t528 + t333 * t320 + t313 * t532 + (t525 + (t415 / 0.6e1 + t359 + t444) * t338) * t504) * t535 + t317 + (t313 * t531 + t309) * t490 + t458 * t484 + ((t336 + t477) * t424 + t523) * t460 + t421 + (t391 + t474) * t431 + (t474 * t398 + t513 * t409 + t390 + t405) * t424 + t335 * t323) + ((t499 * t545 + 0.4e1 * (t492 * t552 + (t339 - t545) * t327 + t450) * t546) * t298 + t314 * (-0.8e1 * (t327 + t495 + t343) * t536 - 0.6e1 * (t508 * t495 + (t320 + t313) * t485 + t458) * t546)) * t326);
t502 = 0.2e1 * t338;
t301 = 0.1e1 / (t325 * t502 + t320 + t328 + t472);
t416 = 0.1e1 / pkin(4);
t543 = t301 * t416;
t533 = t382 * t384;
t413 = 0.1e1 / pkin(6);
t530 = t413 * t416;
t487 = t301 * t413 / pkin(4) ^ 2;
t469 = t544 / 0.4e1;
t468 = -t543 / 0.2e1;
t467 = t530 / 0.2e1;
t452 = -pkin(3) + t493;
t447 = t399 + t415 + t470;
t383 = cos(qJ(3));
t380 = sin(qJ(3));
t312 = t384 * t385 + t534;
t311 = t381 * t385 - t533;
t307 = t328 + t415 + t453;
t297 = (-t452 + t541) * t298 + (t307 * t325 + t308 * t502) * t381 + (t307 * t384 + (0.4e1 * t347 - 0.2e1) * t325 * pkin(3)) * t339;
t296 = (pkin(2) * t533 + t325 * t381) * t298 - (t328 + t447 + t459) * t541 + t452 * t328 + t447 * t493 + (t308 * t552 - t323 - t402 + t498) * pkin(3);
t295 = t296 * t468;
t292 = t467 * t544;
t291 = (t296 * t298 / 0.4e1 + t297 * t469) * t487;
t290 = (t296 * t469 - t297 * t298 / 0.4e1) * t487;
t289 = -t290 * t312 + t291 * t311;
t288 = t290 * t311 + t291 * t312;
t1 = [t385, -t382, 0, 0; t382, t385, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t289, t288, 0, pkin(2); -t288, t289, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t384, -t381, 0, pkin(1); t381, t384, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t295, t297 * t543 / 0.2e1, 0, pkin(3); t297 * t468, t295, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t383, -t380, 0, pkin(5); t380, t383, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t292, -t298 * t530 / 0.2e1, 0, -pkin(4); t298 * t467, t292, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(6); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,7);             % numerisch
else,                         T_mdh = sym('xx', [4,4,7]); end % symbolisch

for i = 1:7
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
