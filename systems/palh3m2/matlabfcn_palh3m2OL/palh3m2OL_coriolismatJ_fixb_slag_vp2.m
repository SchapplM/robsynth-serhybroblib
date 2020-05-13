% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [10x10]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:34:00
% EndTime: 2020-05-07 04:34:50
% DurationCPUTime: 8.48s
% Computational Cost: add. (12665->535), mult. (27075->728), div. (0->0), fcn. (31834->16), ass. (0->337)
t331 = cos(qJ(5));
t328 = sin(qJ(5));
t527 = Ifges(6,4) * t328;
t302 = Ifges(6,1) * t331 - t527;
t320 = Ifges(6,4) * t331;
t301 = Ifges(6,1) * t328 + t320;
t473 = t331 * t301;
t299 = t331 * Ifges(6,2) + t527;
t479 = t328 * t299;
t375 = t479 / 0.2e1 - t473 / 0.2e1;
t546 = -t331 / 0.2e1;
t548 = -t328 / 0.2e1;
t597 = -Ifges(6,2) * t328 + t320;
t411 = -t302 * t548 - t597 * t546 - t375;
t325 = sin(pkin(15));
t326 = cos(pkin(15));
t330 = cos(qJ(8));
t536 = sin(qJ(8));
t290 = t325 * t330 - t326 * t536;
t291 = t325 * t536 + t326 * t330;
t327 = sin(qJ(7));
t539 = sin(qJ(2));
t540 = cos(qJ(7));
t544 = cos(qJ(2));
t292 = -t327 * t539 + t540 * t544;
t293 = -t327 * t544 - t540 * t539;
t168 = t290 * t292 + t291 * t293;
t568 = -t168 / 0.2e1;
t538 = sin(qJ(3));
t543 = cos(qJ(3));
t294 = t538 * t544 + t543 * t539;
t329 = sin(qJ(4));
t354 = t538 * t539 - t543 * t544;
t542 = cos(qJ(4));
t240 = t294 * t329 + t542 * t354;
t501 = t331 * mrSges(6,2);
t504 = t328 * mrSges(6,1);
t400 = t501 + t504;
t618 = t400 * t240;
t627 = pkin(8) * t618;
t626 = Ifges(9,4) * t168;
t323 = t328 ^ 2;
t324 = t331 ^ 2;
t470 = t323 + t324;
t625 = t470 * mrSges(6,3);
t303 = pkin(1) * t327 + pkin(3) * t325;
t304 = t540 * pkin(1) + t326 * pkin(3);
t190 = t290 * t304 - t291 * t303;
t515 = t190 * mrSges(9,1);
t188 = -t290 * t303 - t291 * t304;
t516 = t188 * mrSges(9,2);
t32 = -t515 - t516;
t624 = qJD(8) * t32;
t608 = (t290 * t326 - t291 * t325) * pkin(3);
t512 = t608 * mrSges(9,1);
t229 = (-t290 * t325 - t291 * t326) * pkin(3);
t513 = t229 * mrSges(9,2);
t35 = -t512 - t513;
t623 = t35 * qJD(8);
t342 = -t542 * t294 + t329 * t354;
t552 = -t301 / 0.4e1;
t427 = t552 - t597 / 0.4e1;
t140 = t301 * t342;
t521 = Ifges(6,6) * t240;
t70 = t342 * t597 - t521;
t448 = -t140 / 0.4e1 - t70 / 0.4e1;
t345 = t342 * t427 + t448 + 0.3e1 / 0.4e1 * t521;
t622 = t190 / 0.2e1;
t560 = -t240 / 0.2e1;
t237 = Ifges(5,4) * t240;
t620 = -mrSges(5,2) + t625;
t314 = -t544 * pkin(1) - pkin(12);
t252 = -t354 * pkin(4) + t314;
t109 = -pkin(8) * t240 - pkin(10) * t342 + t252;
t424 = m(6) * t470;
t619 = t109 * t424;
t310 = t329 * t538 * pkin(1);
t468 = t543 * pkin(1);
t275 = -t542 * t468 + t310;
t420 = t470 * t275;
t480 = t328 * t240;
t606 = mrSges(6,2) * t342;
t142 = -mrSges(6,3) * t480 - t606;
t477 = t331 * t142;
t475 = t331 * t240;
t607 = mrSges(6,1) * t342;
t145 = -mrSges(6,3) * t475 + t607;
t482 = t328 * t145;
t616 = t477 / 0.2e1 - t482 / 0.2e1;
t603 = Ifges(9,6) * t168;
t417 = t290 * t293 - t291 * t292;
t605 = Ifges(9,5) * t417;
t615 = -t603 + t605;
t586 = -t605 / 0.2e1 + t603 / 0.2e1;
t151 = pkin(8) * t342 - pkin(10) * t240;
t526 = Ifges(6,5) * t342;
t613 = t302 * t240 + t526;
t522 = Ifges(6,6) * t342;
t612 = t240 * t597 + t522;
t569 = t417 / 0.2e1;
t611 = (Ifges(9,1) * t417 - t626) * t168 / 0.2e1 + (Ifges(9,2) * t417 + t626) * t568 + (0.2e1 * Ifges(9,4) * t417 + (Ifges(9,1) - Ifges(9,2)) * t168) * t569;
t575 = -mrSges(6,1) / 0.2e1;
t559 = t342 / 0.2e1;
t322 = t539 * pkin(1);
t467 = t542 * pkin(4);
t528 = Ifges(5,4) * t342;
t235 = Ifges(5,6) * t342;
t502 = t331 * mrSges(6,1);
t503 = t328 * mrSges(6,2);
t401 = t502 - t503;
t601 = mrSges(5,1) + t401;
t298 = Ifges(6,5) * t328 + Ifges(6,6) * t331;
t493 = t342 * t298;
t600 = m(5) * t252 - mrSges(5,1) * t240 + mrSges(5,2) * t342;
t319 = Ifges(6,5) * t331;
t520 = Ifges(6,6) * t328;
t598 = t319 - t520;
t531 = mrSges(6,3) * t240;
t144 = -t331 * t531 + t607;
t490 = t342 * t331;
t146 = -mrSges(6,1) * t240 - mrSges(6,3) * t490;
t414 = -t467 / 0.2e1;
t534 = pkin(4) * t329;
t317 = pkin(10) + t534;
t551 = -t317 / 0.2e1;
t596 = t144 * t551 + t146 * t414;
t141 = -t328 * t531 - t606;
t491 = t342 * t328;
t143 = mrSges(6,2) * t240 - mrSges(6,3) * t491;
t413 = t467 / 0.2e1;
t550 = t317 / 0.2e1;
t595 = t141 * t550 + t143 * t413;
t138 = t400 * t342;
t465 = t534 / 0.2e1;
t318 = -t467 - pkin(8);
t549 = t318 / 0.2e1;
t594 = t138 * t465 + t549 * t618;
t407 = t538 * t542;
t410 = -t468 + pkin(4);
t270 = pkin(1) * t407 - t329 * t410;
t269 = t542 * t410 + t310;
t264 = -pkin(8) - t269;
t558 = t264 / 0.2e1;
t592 = -t270 * t138 / 0.2e1 + t618 * t558;
t265 = pkin(10) - t270;
t557 = -t265 / 0.2e1;
t572 = -t146 / 0.2e1;
t591 = t144 * t557 + t269 * t572;
t556 = t269 / 0.2e1;
t590 = t265 * t141 / 0.2e1 + t143 * t556;
t381 = -t501 / 0.2e1 - t504 / 0.2e1;
t589 = -t503 / 0.2e1 + t502 / 0.2e1;
t587 = t477 - t482;
t518 = Ifges(6,3) * t342;
t585 = -Ifges(6,6) * t480 / 0.2e1 + Ifges(6,5) * t475 / 0.2e1 + t518 / 0.2e1;
t583 = (t538 * mrSges(4,1) + t543 * mrSges(4,2)) * pkin(1);
t232 = (-t540 * t290 + t291 * t327) * pkin(1);
t221 = t232 * mrSges(9,1);
t233 = (-t290 * t327 - t540 * t291) * pkin(1);
t511 = t233 * mrSges(9,2);
t582 = (t327 * mrSges(8,1) + t540 * mrSges(8,2)) * pkin(1) - t221 + t511;
t581 = mrSges(6,3) * t420;
t580 = (-t319 / 0.2e1 + t520 / 0.2e1) * t240;
t64 = mrSges(9,1) * t168 + mrSges(9,2) * t417;
t236 = Ifges(5,5) * t240;
t579 = t235 / 0.2e1 - t236 / 0.2e1 - t493 / 0.4e1 + t627 / 0.2e1;
t577 = -pkin(10) / 0.2e1;
t576 = pkin(4) * m(6);
t574 = -mrSges(6,2) / 0.2e1;
t573 = mrSges(6,2) / 0.2e1;
t564 = t221 / 0.2e1;
t563 = t229 / 0.2e1;
t562 = t608 / 0.2e1;
t561 = t240 / 0.2e1;
t274 = (t543 * t329 + t407) * pkin(1);
t555 = -t274 / 0.2e1;
t554 = -t275 / 0.2e1;
t553 = t299 / 0.4e1;
t547 = t328 / 0.2e1;
t545 = t331 / 0.2e1;
t541 = cos(qJ(6));
t537 = sin(qJ(6));
t535 = pkin(4) * t294;
t525 = Ifges(6,5) * t240;
t372 = -t151 + t535;
t110 = t322 - t372;
t205 = (-t292 * t326 + t293 * t325) * pkin(3) + t314;
t234 = (-t292 * t325 - t293 * t326) * pkin(3);
t212 = t322 + t234;
t149 = Ifges(5,2) * t240 + t528;
t150 = Ifges(5,1) * t342 + t237;
t72 = t302 * t342 - t525;
t450 = t72 * t546;
t453 = t70 * t547;
t68 = -Ifges(6,3) * t240 + t342 * t598;
t332 = -t252 * (mrSges(5,1) * t342 + mrSges(5,2) * t240) - t314 * (-t294 * mrSges(4,1) + t354 * mrSges(4,2)) + t149 * t559 + (t240 * t598 + t518) * t561 - t354 ^ 2 * Ifges(4,4) + t240 * t453 + t612 * t491 / 0.2e1 + t240 * t450 - t613 * t490 / 0.2e1 + (Ifges(4,4) * t294 + (Ifges(4,1) - Ifges(4,2)) * t354) * t294 - (Ifges(5,1) * t240 - t528 + t68) * t342 / 0.2e1 + (-t328 * t142 - t331 * t145) * t109 + (-Ifges(5,2) * t342 + t150 + t237) * t560;
t335 = t314 * (-mrSges(8,1) * t293 + mrSges(8,2) * t292) - Ifges(8,4) * t293 ^ 2 + t611;
t390 = t328 * t143 + t331 * t146;
t402 = Ifges(8,4) * t292 + (-Ifges(8,1) + Ifges(8,2)) * t293;
t65 = -mrSges(9,1) * t417 + mrSges(9,2) * t168;
t1 = (t390 + t619) * t110 + t600 * (t322 - t535) + (-mrSges(8,2) * t293 - t354 * mrSges(4,1) - t294 * mrSges(4,2) + (m(4) + m(8)) * t314) * t322 + (-mrSges(8,1) * t322 + t402) * t292 - pkin(12) * (t539 * mrSges(3,1) + t544 * mrSges(3,2)) + t335 + (m(9) * t212 + t64) * t205 - t332 + t212 * t65 + (Ifges(3,1) - Ifges(3,2)) * t544 * t539 + (-t539 ^ 2 + t544 ^ 2) * Ifges(3,4);
t517 = t1 * qJD(1);
t359 = t331 * t372;
t360 = t328 * t372;
t2 = t143 * t360 + t146 * t359 + t372 * t619 + t600 * t535 + t332;
t514 = t2 * qJD(1);
t510 = t269 * mrSges(5,2);
t509 = t270 * mrSges(5,1);
t508 = t274 * mrSges(5,1);
t507 = t275 * mrSges(5,2);
t3 = (t68 / 0.2e1 - t149 / 0.2e1 + t252 * mrSges(5,1) + t613 * t545 + t612 * t548 - t528 / 0.2e1) * t342 - (-t237 / 0.2e1 - t252 * mrSges(5,2) - t150 / 0.2e1 + t450 + t453 - t580 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t342) * t240 + t390 * t151 + (t328 * t141 + t331 * t144 + t151 * t424) * t109;
t505 = t3 * qJD(1);
t11 = t234 * t65 + (m(9) * t234 + t64) * t205 + t402 * t292 + t335;
t500 = t11 * qJD(1);
t499 = t168 * t608;
t498 = t417 * t229;
t139 = t342 * t299;
t476 = t331 * t143;
t481 = t328 * t146;
t17 = (t476 - t481) * t109 + (-t140 * t545 + t298 * t561 + t70 * t546 + (-t139 + t72) * t548) * t342;
t497 = t17 * qJD(1);
t18 = t205 * t64 + t611;
t496 = t18 * qJD(1);
t495 = t188 * t417;
t494 = t190 * t168;
t492 = t342 * t329;
t489 = t264 * t618;
t488 = t269 * t240;
t487 = t270 * t342;
t486 = t270 * t401;
t485 = t274 * t401;
t483 = t318 * t618;
t152 = -pkin(6) * (t537 * mrSges(7,1) + t541 * mrSges(7,2)) + (t537 ^ 2 - t541 ^ 2) * Ifges(7,4) + (-Ifges(7,1) + Ifges(7,2)) * t541 * t537;
t469 = t152 * qJD(1);
t466 = -t534 / 0.2e1;
t464 = t143 * t577;
t463 = pkin(10) * t572;
t458 = -t518 / 0.2e1;
t449 = -t139 / 0.4e1 + t72 / 0.4e1;
t447 = t542 * t240;
t446 = -t240 * t319 / 0.4e1;
t443 = t143 * t557;
t441 = t146 * t557;
t437 = t138 * t555;
t435 = t143 * t551;
t433 = t146 * t551;
t426 = t302 / 0.4e1 - t299 / 0.4e1;
t425 = t550 + t557;
t421 = t470 * t269;
t419 = t470 * t317;
t412 = Ifges(8,5) * t292 + Ifges(8,6) * t293 + t615;
t406 = mrSges(6,3) * (-t324 / 0.2e1 - t323 / 0.2e1);
t405 = t612 / 0.4e1 - t240 * t552;
t404 = t613 / 0.4e1 - t240 * t553;
t379 = pkin(8) * t400;
t395 = -t379 / 0.2e1 + t411;
t347 = mrSges(6,3) * t421 + t486 + t509 - t510;
t31 = m(6) * (-t264 * t270 + t265 * t421) + t347;
t336 = (t553 + t527 / 0.4e1 + (Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1) * t331) * t240 + pkin(10) * t145 / 0.2e1 - t526 / 0.4e1 + t404;
t363 = (t298 / 0.4e1 - Ifges(5,6) / 0.2e1) * t342 + Ifges(5,5) * t561;
t340 = t363 + t579;
t341 = (t552 - t320 / 0.4e1) * t240 + t142 * t577 - t522 / 0.4e1 + t405;
t4 = (t341 + t590) * t331 + (t336 + t591) * t328 + t340 + t592;
t394 = t4 * qJD(1) + t31 * qJD(2);
t346 = t342 * t555 - t240 * t554 - t488 / 0.2e1 + t487 / 0.2e1;
t16 = t437 + (t558 - t318 / 0.2e1) * t618 + (t275 * t143 / 0.2e1 - t425 * t142) * t331 + (t425 * t145 + t146 * t554) * t328 + ((t492 / 0.2e1 + t447 / 0.2e1) * pkin(4) + t346) * mrSges(5,3);
t33 = t601 * t274 + t583 + t620 * t275 + m(6) * (-t264 * t274 + t265 * t420) + m(5) * (t269 * t274 - t270 * t275);
t393 = t16 * qJD(1) + t33 * qJD(2);
t378 = t232 * t568 + t233 * t569;
t20 = ((-t188 / 0.2e1 + t563) * t417 - (t622 - t608 / 0.2e1) * t168 + t378) * mrSges(9,3);
t34 = -m(9) * (t188 * t232 + t190 * t233) + t582;
t392 = t20 * qJD(1) - t34 * qJD(2);
t391 = t32 * qJD(2);
t388 = pkin(10) * t406;
t361 = t342 * t426 + t449;
t349 = Ifges(6,5) * t560 + t361;
t135 = t401 * t342;
t377 = t135 * t558 + t446;
t386 = t265 * t406;
t10 = t458 + t342 * t386 + (t110 * t573 + t345 + t443) * t328 + (t110 * t575 + t349 + t441) * t331 + t377;
t374 = t264 * t400;
t120 = t374 + t411;
t387 = t10 * qJD(1) + t120 * qJD(2);
t385 = t317 * t406;
t384 = -pkin(8) * t135 / 0.2e1 + t446;
t380 = t521 / 0.4e1 + t448;
t376 = t135 * t549 + t446;
t373 = t318 * t400;
t371 = t379 / 0.2e1;
t370 = t470 * t542;
t369 = t381 * t269;
t368 = t381 * t275;
t367 = -t374 / 0.2e1;
t337 = t620 * t467 - t601 * t534;
t121 = (t370 * t317 + t318 * t329) * t576 + t337;
t333 = m(6) * (-t318 * t270 + t269 * t419 + (t264 * t329 + t370 * t265) * pkin(4)) / 0.2e1 - t510 / 0.2e1 + t509 / 0.2e1 + t486 / 0.2e1 + mrSges(5,1) * t466 - t401 * t465 + mrSges(5,2) * t414 + (t413 + t556) * t625;
t334 = -m(6) * (pkin(8) * t274 + pkin(10) * t420) / 0.2e1 - t508 / 0.2e1 - t485 / 0.2e1 + t507 / 0.2e1 - t581 / 0.2e1;
t24 = t333 + t334;
t6 = (t341 + t595) * t331 + (t336 + t596) * t328 + t340 + t594;
t366 = t6 * qJD(1) + t24 * qJD(2) + t121 * qJD(3);
t21 = t564 + (t622 + t562) * mrSges(9,1) + (-t233 / 0.2e1 + t188 / 0.2e1 + t563) * mrSges(9,2);
t343 = (t168 * t562 + t568 * t608) * mrSges(9,3) - t586;
t28 = t343 + t586;
t365 = t28 * qJD(1) - t21 * qJD(2) + t35 * qJD(7);
t13 = t458 + t342 * t385 + (t372 * t574 + t345 + t435) * t328 + (t433 + t372 * mrSges(6,1) / 0.2e1 + t349) * t331 + t376;
t147 = t373 + t411;
t338 = -t373 / 0.2e1 - t411;
t52 = t367 + t368 + t338;
t364 = t13 * qJD(1) - t52 * qJD(2) + t147 * qJD(3);
t350 = t331 * t612;
t351 = t328 * t613;
t358 = t493 / 0.2e1 + t236 + Ifges(4,6) * t294 + Ifges(4,5) * t354 + t351 / 0.2e1 + t350 / 0.2e1 - t235 - t375 * t240;
t15 = (-Ifges(6,3) / 0.2e1 + t388) * t342 + (t463 - t525 / 0.2e1 + t151 * t575 + t361) * t331 + (t151 * t573 + t345 + t464) * t328 + t384;
t153 = -t379 + t411;
t54 = t371 + t367 + t369 - t411;
t348 = t381 * t467;
t81 = t371 + t348 + t338;
t353 = t15 * qJD(1) - t54 * qJD(2) - t81 * qJD(3) + t153 * qJD(4);
t352 = t427 * t328 + t426 * t331;
t339 = t351 / 0.4e1 + t350 / 0.4e1 + t363 - t579 + (-t479 / 0.4e1 + t473 / 0.4e1) * t240 + t616 * pkin(10);
t268 = t373 / 0.2e1;
t243 = t374 / 0.2e1;
t82 = t268 + t348 + t395;
t55 = t243 + t369 + t395;
t53 = t268 + t243 + t368 + t411;
t27 = t343 - t586;
t25 = -0.2e1 * t586;
t23 = t333 - t334;
t22 = -t515 / 0.2e1 - t516 / 0.2e1 - t513 / 0.2e1 - t512 / 0.2e1 + t564 - t511 / 0.2e1;
t19 = t412 + (-t498 / 0.2e1 - t499 / 0.2e1 - t495 / 0.2e1 - t494 / 0.2e1 + t378) * mrSges(9,3);
t14 = Ifges(6,3) * t559 + (t464 + t380) * t328 + (t463 + t449) * t331 + (t388 + t352) * t342 + t384 - t580 + t589 * t151;
t12 = -t360 * t574 + t359 * t575 + (t435 + t380) * t328 + (t433 + t449) * t331 + (t385 + t352) * t342 + t376 + t585;
t9 = (t443 + t380) * t328 + (t441 + t449) * t331 + (t386 + t352) * t342 + t377 + t589 * t110 + t585;
t8 = t358 + (t476 / 0.2e1 - t481 / 0.2e1) * t275 + t483 / 0.2e1 + t437 + t489 / 0.2e1 + (t240 * t414 + t342 * t466 + t346) * mrSges(5,3) + (t317 + t265) * t616;
t7 = (t405 + t595) * t331 + (t404 + t596) * t328 + t339 + t594;
t5 = t339 + (t405 + t590) * t331 + (t404 + t591) * t328 + t592;
t26 = [qJD(2) * t1 - qJD(3) * t2 + qJD(4) * t3 + qJD(5) * t17 - qJD(6) * t152 + qJD(7) * t11 + qJD(8) * t18, t8 * qJD(3) + t5 * qJD(4) + t9 * qJD(5) + t19 * qJD(7) + t25 * qJD(8) + t517 + (Ifges(3,5) * t544 - Ifges(3,6) * t539 + t358 + t412 + t489 + ((-t540 * t292 + t327 * t293) * mrSges(8,3) + (-t538 * t294 + t543 * t354) * mrSges(4,3)) * pkin(1) + t587 * t265 + (-t494 - t495) * mrSges(9,3) + (t487 - t488) * mrSges(5,3)) * qJD(2), -t514 + t8 * qJD(2) + (t483 + t587 * t317 + (-t447 - t492) * mrSges(5,3) * pkin(4) + t358) * qJD(3) + t7 * qJD(4) + t12 * qJD(5), t5 * qJD(2) + t7 * qJD(3) + t14 * qJD(5) + t505 + (-t627 - t235 + t298 * t559 + t612 * t545 + t613 * t547 - (-Ifges(5,5) + t375) * t240 + (t331 * t141 - t328 * t144) * pkin(10)) * qJD(4), t497 + t9 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + (-t109 * t400 - t493) * qJD(5), -t469 + (Ifges(7,5) * t541 - Ifges(7,6) * t537) * qJD(6), t500 + t19 * qJD(2) + ((-t498 - t499) * mrSges(9,3) + t412) * qJD(7) + t27 * qJD(8), t25 * qJD(2) + t27 * qJD(7) + t615 * qJD(8) + t496, 0, 0; qJD(3) * t16 + qJD(4) * t4 + qJD(5) * t10 + qJD(7) * t20 - t517, qJD(3) * t33 + qJD(4) * t31 + qJD(5) * t120 - qJD(7) * t34 + t624, (t485 + m(6) * (-t274 * t318 + t275 * t419) - t507 + t508 + m(5) * (t542 * t274 + t275 * t329) * pkin(4) + t581 + t583) * qJD(3) + t23 * qJD(4) + t53 * qJD(5) + t393, t23 * qJD(3) + (m(6) * (pkin(8) * t270 + pkin(10) * t421) + t347) * qJD(4) + t55 * qJD(5) + t394, t53 * qJD(3) + t55 * qJD(4) + (-t265 * t401 + t598) * qJD(5) + t387, 0, (m(9) * (t229 * t232 + t233 * t608) - t582) * qJD(7) + t22 * qJD(8) + t392, t22 * qJD(7) + t391 + t624, 0, 0; -qJD(2) * t16 + qJD(4) * t6 + qJD(5) * t13 + t514, qJD(4) * t24 - qJD(5) * t52 - t393, qJD(4) * t121 + qJD(5) * t147, ((-pkin(8) * t329 + pkin(10) * t370) * t576 + t337) * qJD(4) + t82 * qJD(5) + t366, t82 * qJD(4) + (-t317 * t401 + t598) * qJD(5) + t364, 0, 0, 0, 0, 0; -qJD(2) * t4 - qJD(3) * t6 + qJD(5) * t15 - t505, -qJD(3) * t24 - qJD(5) * t54 - t394, -qJD(5) * t81 - t366, t153 * qJD(5), (-pkin(10) * t401 + t598) * qJD(5) + t353, 0, 0, 0, 0, 0; -qJD(2) * t10 - qJD(3) * t13 - qJD(4) * t15 - t497, qJD(3) * t52 + qJD(4) * t54 - t387, qJD(4) * t81 - t364, -t353, 0, 0, 0, 0, 0, 0; t469, 0, 0, 0, 0, 0, 0, 0, 0, 0; -qJD(2) * t20 + qJD(8) * t28 - t500, -qJD(8) * t21 - t392, 0, 0, 0, 0, t623, t365 + t623, 0, 0; -qJD(7) * t28 - t496, qJD(7) * t21 - t391, 0, 0, 0, 0, -t365, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t26;
