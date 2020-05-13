% Calculate kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2OL_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:18:31
% EndTime: 2020-05-03 01:18:33
% DurationCPUTime: 1.99s
% Computational Cost: add. (1650->280), mult. (1551->474), div. (0->0), fcn. (1441->16), ass. (0->165)
t435 = sin(qJ(4));
t436 = sin(qJ(3));
t437 = sin(qJ(2));
t440 = cos(qJ(4));
t441 = cos(qJ(3));
t442 = cos(qJ(2));
t503 = pkin(5) * (t442 * (-t436 * t435 + t441 * t440) - t437 * (t441 * t435 + t436 * t440));
t502 = pkin(5) * ((t442 * t435 + t437 * t440) * t441 + t436 * (-t437 * t435 + t442 * t440));
t500 = pkin(4) * qJD(2);
t499 = Icges(3,4) * t437;
t498 = Icges(3,4) * t442;
t433 = qJ(2) + qJ(3);
t429 = sin(t433);
t497 = Icges(4,4) * t429;
t430 = cos(t433);
t496 = Icges(4,4) * t430;
t431 = qJ(4) + t433;
t420 = sin(t431);
t495 = Icges(5,4) * t420;
t421 = cos(t431);
t494 = Icges(5,4) * t421;
t422 = qJ(5) + t431;
t416 = sin(t422);
t493 = Icges(6,4) * t416;
t417 = cos(t422);
t492 = Icges(6,4) * t417;
t438 = sin(qJ(1));
t491 = t416 * t438;
t443 = cos(qJ(1));
t490 = t416 * t443;
t434 = sin(qJ(6));
t489 = t438 * t434;
t439 = cos(qJ(6));
t488 = t438 * t439;
t487 = t443 * t434;
t486 = t443 * t439;
t428 = qJD(2) * t438;
t405 = qJD(3) * t438 + t428;
t485 = pkin(2) * t436;
t484 = pkin(3) * t417;
t483 = qJD(2) * t443;
t482 = qJD(6) * t416;
t481 = t443 * qJD(1);
t432 = qJD(2) + qJD(3);
t480 = t437 * t500;
t397 = qJD(4) * t438 + t405;
t418 = t441 * pkin(2) + pkin(4);
t387 = t418 * t442 - t437 * t485 + pkin(1);
t479 = -t387 - t503;
t478 = qJD(4) + t432;
t375 = qJD(5) * t438 + t397;
t423 = t442 * t500;
t477 = pkin(2) * t432 * t430 + t423;
t414 = t442 * rSges(3,1) - t437 * rSges(3,2);
t476 = rSges(4,1) * t430 - rSges(4,2) * t429;
t475 = rSges(5,1) * t421 - rSges(5,2) * t420;
t474 = rSges(6,1) * t417 - rSges(6,2) * t416;
t413 = t439 * rSges(7,1) - rSges(7,2) * t434;
t473 = -rSges(7,3) * t416 + t413 * t417;
t472 = Icges(3,1) * t442 - t499;
t471 = Icges(4,1) * t430 - t497;
t470 = Icges(5,1) * t421 - t495;
t469 = Icges(6,1) * t417 - t493;
t468 = -Icges(3,2) * t437 + t498;
t467 = -Icges(4,2) * t429 + t496;
t466 = -Icges(5,2) * t420 + t494;
t465 = -Icges(6,2) * t416 + t492;
t464 = Icges(3,5) * t442 - Icges(3,6) * t437;
t463 = Icges(4,5) * t430 - Icges(4,6) * t429;
t462 = Icges(5,5) * t421 - Icges(5,6) * t420;
t461 = Icges(6,5) * t417 - Icges(6,6) * t416;
t379 = -Icges(3,6) * t443 + t468 * t438;
t381 = -Icges(3,5) * t443 + t472 * t438;
t460 = t379 * t437 - t381 * t442;
t380 = Icges(3,6) * t438 + t468 * t443;
t382 = Icges(3,5) * t438 + t472 * t443;
t459 = -t380 * t437 + t382 * t442;
t408 = Icges(3,2) * t442 + t499;
t409 = Icges(3,1) * t437 + t498;
t458 = -t408 * t437 + t409 * t442;
t457 = pkin(1) + t414;
t376 = (-qJD(5) - t478) * t443;
t456 = -pkin(2) * qJD(3) * (t442 * t436 + t437 * t441) - qJD(2) * (t437 * t418 + t442 * t485);
t455 = (Icges(6,5) * t416 + Icges(6,6) * t417) * qJD(1) + (-Icges(6,3) * t443 + t461 * t438) * t376 + (Icges(6,3) * t438 + t461 * t443) * t375;
t398 = t478 * t443;
t454 = (Icges(5,5) * t420 + Icges(5,6) * t421) * qJD(1) - (-Icges(5,3) * t443 + t462 * t438) * t398 + (Icges(5,3) * t438 + t462 * t443) * t397;
t406 = t432 * t443;
t453 = (Icges(4,5) * t429 + Icges(4,6) * t430) * qJD(1) - (-Icges(4,3) * t443 + t463 * t438) * t406 + (Icges(4,3) * t438 + t463 * t443) * t405;
t452 = pkin(5) * t421 * t478 + t477;
t451 = t456 * t443;
t450 = t387 * t481 + t456 * t438;
t449 = -t398 * t502 + t451;
t448 = -t397 * t502 + t481 * t503 + t450;
t352 = -Icges(6,6) * t443 + t465 * t438;
t353 = Icges(6,6) * t438 + t465 * t443;
t354 = -Icges(6,5) * t443 + t469 * t438;
t355 = Icges(6,5) * t438 + t469 * t443;
t389 = Icges(6,2) * t417 + t493;
t390 = Icges(6,1) * t416 + t492;
t447 = (-t353 * t416 + t355 * t417) * t375 + (-t352 * t416 + t354 * t417) * t376 + (-t389 * t416 + t390 * t417) * qJD(1);
t360 = -Icges(5,6) * t443 + t466 * t438;
t361 = Icges(5,6) * t438 + t466 * t443;
t362 = -Icges(5,5) * t443 + t470 * t438;
t363 = Icges(5,5) * t438 + t470 * t443;
t394 = Icges(5,2) * t421 + t495;
t395 = Icges(5,1) * t420 + t494;
t446 = (-t361 * t420 + t363 * t421) * t397 - (-t360 * t420 + t362 * t421) * t398 + (-t394 * t420 + t395 * t421) * qJD(1);
t369 = -Icges(4,6) * t443 + t467 * t438;
t370 = Icges(4,6) * t438 + t467 * t443;
t371 = -Icges(4,5) * t443 + t471 * t438;
t372 = Icges(4,5) * t438 + t471 * t443;
t400 = Icges(4,2) * t430 + t497;
t401 = Icges(4,1) * t429 + t496;
t445 = (-t370 * t429 + t372 * t430) * t405 - (-t369 * t429 + t371 * t430) * t406 + (-t400 * t429 + t401 * t430) * qJD(1);
t419 = t442 * pkin(4) + pkin(1);
t415 = t443 * rSges(2,1) - t438 * rSges(2,2);
t412 = t434 * rSges(7,1) + t439 * rSges(7,2);
t411 = t438 * rSges(2,1) + t443 * rSges(2,2);
t410 = t437 * rSges(3,1) + t442 * rSges(3,2);
t407 = Icges(3,5) * t437 + Icges(3,6) * t442;
t404 = qJD(6) * t417 + qJD(1);
t402 = t429 * rSges(4,1) + t430 * rSges(4,2);
t396 = t420 * rSges(5,1) + t421 * rSges(5,2);
t391 = t416 * rSges(6,1) + t417 * rSges(6,2);
t386 = t417 * t486 - t489;
t385 = -t417 * t487 - t488;
t384 = t417 * t488 + t487;
t383 = -t417 * t489 + t486;
t378 = Icges(3,3) * t438 + t464 * t443;
t377 = -Icges(3,3) * t443 + t464 * t438;
t374 = t438 * rSges(4,3) + t476 * t443;
t373 = -t443 * rSges(4,3) + t476 * t438;
t365 = t438 * rSges(5,3) + t475 * t443;
t364 = -t443 * rSges(5,3) + t475 * t438;
t357 = t438 * rSges(6,3) + t474 * t443;
t356 = -t443 * rSges(6,3) + t474 * t438;
t349 = -t438 * t482 + t376;
t348 = -t443 * t482 + t375;
t347 = t417 * rSges(7,3) + t413 * t416;
t344 = Icges(7,5) * t417 + (Icges(7,1) * t439 - Icges(7,4) * t434) * t416;
t343 = Icges(7,6) * t417 + (Icges(7,4) * t439 - Icges(7,2) * t434) * t416;
t342 = Icges(7,3) * t417 + (Icges(7,5) * t439 - Icges(7,6) * t434) * t416;
t340 = -t410 * t428 + (t438 * rSges(3,3) + t457 * t443) * qJD(1);
t339 = -t410 * t483 + (t443 * rSges(3,3) - t457 * t438) * qJD(1);
t338 = -t438 * t412 + t473 * t443;
t337 = t443 * t412 + t473 * t438;
t336 = Icges(7,1) * t386 + Icges(7,4) * t385 - Icges(7,5) * t490;
t335 = Icges(7,1) * t384 + Icges(7,4) * t383 - Icges(7,5) * t491;
t334 = Icges(7,4) * t386 + Icges(7,2) * t385 - Icges(7,6) * t490;
t333 = Icges(7,4) * t384 + Icges(7,2) * t383 - Icges(7,6) * t491;
t332 = Icges(7,5) * t386 + Icges(7,6) * t385 - Icges(7,3) * t490;
t331 = Icges(7,5) * t384 + Icges(7,6) * t383 - Icges(7,3) * t491;
t330 = -t438 * t480 - t405 * t402 + (t419 * t443 + t374) * qJD(1);
t329 = -t443 * t480 - t406 * t402 + (-t419 * t438 - t373) * qJD(1);
t328 = t405 * t373 + t406 * t374 + t423;
t327 = t397 * t364 + t398 * t365 + t477;
t326 = qJD(1) * t365 - t397 * t396 + t450;
t325 = -t398 * t396 + t451 + (-t387 * t438 - t364) * qJD(1);
t324 = t375 * t356 - t376 * t357 + t452;
t323 = qJD(1) * t357 - t375 * t391 + t448;
t322 = t376 * t391 + (t479 * t438 - t356) * qJD(1) + t449;
t321 = t348 * t337 - t349 * t338 + (t438 * t375 - t443 * t376) * t484 + t452;
t320 = t404 * t338 - t348 * t347 + (-t416 * t375 + t417 * t481) * pkin(3) + t448;
t319 = t376 * t416 * pkin(3) - t404 * t337 + t349 * t347 + (t479 - t484) * t438 * qJD(1) + t449;
t1 = -((-t443 * t407 + t458 * t438) * qJD(1) + (t443 ^ 2 * t377 + (t459 * t438 + (-t378 + t460) * t443) * t438) * qJD(2)) * t483 / 0.2e1 + ((t438 * t407 + t458 * t443) * qJD(1) + (t438 ^ 2 * t378 + (t460 * t443 + (-t377 + t459) * t438) * t443) * qJD(2)) * t428 / 0.2e1 + m(3) * (t414 ^ 2 * qJD(2) ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(7) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(6) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(4) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + t349 * ((-t332 * t491 + t383 * t334 + t384 * t336) * t348 + (-t331 * t491 + t383 * t333 + t384 * t335) * t349 + (-t342 * t491 + t383 * t343 + t384 * t344) * t404) / 0.2e1 + t404 * ((t331 * t349 + t332 * t348 + t342 * t404) * t417 + ((-t334 * t434 + t336 * t439) * t348 + (-t333 * t434 + t335 * t439) * t349 + (-t343 * t434 + t344 * t439) * t404) * t416) / 0.2e1 + t348 * ((-t332 * t490 + t385 * t334 + t386 * t336) * t348 + (-t331 * t490 + t385 * t333 + t386 * t335) * t349 + (-t342 * t490 + t385 * t343 + t386 * t344) * t404) / 0.2e1 + t375 * (t455 * t438 + t447 * t443) / 0.2e1 + t376 * (t447 * t438 - t455 * t443) / 0.2e1 + t405 * (t453 * t438 + t445 * t443) / 0.2e1 - t406 * (t445 * t438 - t453 * t443) / 0.2e1 + t397 * (t454 * t438 + t446 * t443) / 0.2e1 - t398 * (t446 * t438 - t454 * t443) / 0.2e1 + m(5) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t411 ^ 2 + t415 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t430 * t370 + t429 * t372) * t405 - (t430 * t369 + t429 * t371) * t406 + (t421 * t361 + t420 * t363) * t397 - (t421 * t360 + t420 * t362) * t398 + ((t442 * t380 + t437 * t382) * t438 - (t442 * t379 + t437 * t381) * t443) * qJD(2) + (t417 * t353 + t416 * t355) * t375 + (t417 * t352 + t416 * t354) * t376 + (t417 * t389 + t416 * t390 + t421 * t394 + t420 * t395 + t430 * t400 + t429 * t401 + t442 * t408 + t437 * t409) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
