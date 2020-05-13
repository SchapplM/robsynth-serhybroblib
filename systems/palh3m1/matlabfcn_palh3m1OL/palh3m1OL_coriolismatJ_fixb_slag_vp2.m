% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m1OL
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
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m1OL_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1OL_coriolismatJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:04:21
% EndTime: 2020-04-20 17:05:14
% DurationCPUTime: 9.66s
% Computational Cost: add. (12665->520), mult. (27075->711), div. (0->0), fcn. (31834->16), ass. (0->335)
t331 = sin(pkin(15));
t332 = cos(pkin(15));
t335 = cos(qJ(8));
t536 = sin(qJ(8));
t299 = t331 * t335 - t332 * t536;
t300 = t331 * t536 + t332 * t335;
t537 = sin(qJ(7));
t540 = sin(qJ(2));
t541 = cos(qJ(7));
t545 = cos(qJ(2));
t301 = -t537 * t540 + t541 * t545;
t302 = -t537 * t545 - t541 * t540;
t176 = t299 * t301 + t300 * t302;
t419 = t299 * t302 - t300 * t301;
t643 = Ifges(9,5) * t419 - Ifges(9,6) * t176;
t642 = t643 * qJD(8);
t561 = t176 / 0.2e1;
t640 = -t561 - t176 / 0.2e1;
t333 = sin(qJ(5));
t336 = cos(qJ(5));
t527 = Ifges(6,4) * t336;
t309 = Ifges(6,1) * t333 + t527;
t472 = t336 * t309;
t528 = Ifges(6,4) * t333;
t308 = Ifges(6,2) * t336 + t528;
t475 = t333 * t308;
t384 = t475 / 0.2e1 - t472 / 0.2e1;
t402 = -Ifges(6,2) * t333 + t527;
t403 = Ifges(6,1) * t336 - t528;
t622 = -t336 / 0.2e1;
t623 = -t333 / 0.2e1;
t415 = -t402 * t622 - t403 * t623 - t384;
t307 = Ifges(6,5) * t333 + Ifges(6,6) * t336;
t539 = sin(qJ(3));
t544 = cos(qJ(3));
t303 = t539 * t545 + t544 * t540;
t334 = sin(qJ(4));
t361 = t539 * t540 - t544 * t545;
t543 = cos(qJ(4));
t348 = -t543 * t303 + t334 * t361;
t597 = t348 * t307;
t600 = Ifges(5,6) * t348;
t249 = t303 * t334 + t543 * t361;
t611 = Ifges(6,6) * t348 + t402 * t249;
t630 = t336 * t611;
t610 = Ifges(6,5) * t348 + t403 * t249;
t631 = t333 * t610;
t639 = t597 / 0.2e1 - t600 + t630 / 0.2e1 + t631 / 0.2e1;
t558 = -t249 / 0.4e1;
t624 = t249 / 0.4e1;
t638 = t472 * t624 + t475 * t558 + t630 / 0.4e1 + t631 / 0.4e1;
t636 = -t610 / 0.4e1;
t635 = -t611 / 0.4e1;
t621 = Ifges(5,5) * t249;
t634 = -t621 / 0.2e1;
t502 = t336 * mrSges(6,2);
t508 = t333 * mrSges(6,1);
t404 = t502 + t508;
t616 = t404 * t249;
t633 = pkin(8) * t616;
t629 = -t597 / 0.4e1 + t600 / 0.2e1 + t634;
t626 = -t249 / 0.2e1;
t625 = t249 / 0.2e1;
t246 = Ifges(5,4) * t249;
t329 = t333 ^ 2;
t330 = t336 ^ 2;
t469 = t329 + t330;
t620 = t469 * mrSges(6,3);
t462 = t537 * pkin(1);
t310 = t331 * pkin(3) + t462;
t463 = t541 * pkin(1);
t311 = t332 * pkin(3) + t463;
t199 = t299 * t311 - t300 * t310;
t185 = t199 * mrSges(9,1);
t197 = -t299 * t310 - t300 * t311;
t32 = t197 * mrSges(9,2) + t185;
t619 = qJD(8) * t32;
t321 = -t545 * pkin(1) - pkin(12);
t261 = -t361 * pkin(4) + t321;
t119 = -pkin(8) * t249 - pkin(10) * t348 + t261;
t425 = m(6) * t469;
t618 = t119 * t425;
t606 = (t299 * t332 - t300 * t331) * pkin(3);
t230 = t606 * mrSges(9,1);
t239 = (-t299 * t331 - t300 * t332) * pkin(3);
t35 = t239 * mrSges(9,2) + t230;
t617 = t35 * qJD(8);
t490 = t249 * t336;
t605 = mrSges(6,1) * t348;
t153 = -mrSges(6,3) * t490 + t605;
t497 = t153 * t333;
t476 = t333 * t249;
t604 = mrSges(6,2) * t348;
t150 = -mrSges(6,3) * t476 - t604;
t498 = t150 * t336;
t615 = -t497 / 0.2e1 + t498 / 0.2e1;
t159 = pkin(8) * t348 - pkin(10) * t249;
t609 = t640 * Ifges(9,4) * t176 + (0.2e1 * Ifges(9,1) * t561 + Ifges(9,4) * t419 + t640 * Ifges(9,2)) * t419;
t575 = -mrSges(6,1) / 0.2e1;
t317 = t334 * t539 * pkin(1);
t465 = t544 * pkin(1);
t414 = -t465 + pkin(4);
t278 = t543 * t414 + t317;
t273 = -pkin(8) - t278;
t555 = t273 / 0.2e1;
t328 = t540 * pkin(1);
t529 = Ifges(5,4) * t348;
t408 = (t329 / 0.2e1 + t330 / 0.2e1) * mrSges(6,3);
t596 = m(5) * t261 - mrSges(5,1) * t249 + mrSges(5,2) * t348;
t326 = Ifges(6,5) * t336;
t521 = Ifges(6,6) * t333;
t593 = t326 - t521;
t531 = mrSges(6,3) * t249;
t149 = -t333 * t531 - t604;
t474 = t336 * t149;
t152 = -t336 * t531 + t605;
t478 = t333 * t152;
t592 = t474 / 0.2e1 - t478 / 0.2e1;
t488 = t348 * t333;
t151 = mrSges(6,2) * t249 - mrSges(6,3) * t488;
t473 = t336 * t151;
t487 = t348 * t336;
t154 = -mrSges(6,1) * t249 - mrSges(6,3) * t487;
t477 = t333 * t154;
t591 = -t477 / 0.2e1 + t473 / 0.2e1;
t503 = t336 * mrSges(6,1);
t507 = t333 * mrSges(6,2);
t589 = t503 / 0.2e1 - t507 / 0.2e1;
t509 = t330 * mrSges(6,3);
t510 = t329 * mrSges(6,3);
t588 = t509 / 0.2e1 + t510 / 0.2e1;
t587 = -t497 + t498;
t584 = (t539 * mrSges(4,1) + t544 * mrSges(4,2)) * pkin(1);
t494 = t419 * t239;
t495 = t176 * t606;
t583 = (-t495 / 0.2e1 - t494 / 0.2e1) * mrSges(9,3);
t464 = t543 * pkin(4);
t582 = t464 * t620;
t581 = (-t326 / 0.2e1 + t521 / 0.2e1) * t249;
t71 = mrSges(9,1) * t176 + mrSges(9,2) * t419;
t580 = t633 / 0.2e1 + t629;
t578 = -pkin(10) / 0.2e1;
t577 = m(6) / 0.2e1;
t576 = pkin(4) * m(6);
t574 = -mrSges(6,2) / 0.2e1;
t573 = mrSges(6,2) / 0.2e1;
t572 = Ifges(6,3) / 0.2e1;
t522 = Ifges(6,6) * t249;
t79 = t348 * t402 - t522;
t571 = -t79 / 0.4e1;
t525 = Ifges(6,5) * t249;
t82 = t348 * t403 - t525;
t570 = t82 / 0.4e1;
t148 = t404 * t348;
t569 = t148 / 0.2e1;
t568 = -t151 / 0.2e1;
t567 = t151 / 0.2e1;
t411 = t539 * t543;
t279 = pkin(1) * t411 - t334 * t414;
t274 = pkin(10) - t279;
t554 = -t274 / 0.2e1;
t283 = (t544 * t334 + t411) * pkin(1);
t553 = -t283 / 0.2e1;
t284 = -t543 * t465 + t317;
t552 = -t284 / 0.2e1;
t551 = -t308 / 0.4e1;
t550 = -t309 / 0.4e1;
t532 = t334 * pkin(4);
t324 = pkin(10) + t532;
t549 = -t324 / 0.2e1;
t548 = t324 / 0.2e1;
t325 = -t464 - pkin(8);
t547 = t325 / 0.2e1;
t546 = t336 / 0.2e1;
t542 = cos(qJ(6));
t538 = sin(qJ(6));
t535 = pkin(4) * t303;
t519 = Ifges(6,3) * t348;
t376 = -t159 + t535;
t120 = t328 - t376;
t214 = (-t301 * t332 + t302 * t331) * pkin(3) + t321;
t243 = (-t301 * t331 - t302 * t332) * pkin(3);
t221 = t328 + t243;
t157 = Ifges(5,2) * t249 + t529;
t158 = Ifges(5,1) * t348 + t246;
t445 = t82 * t622;
t506 = t333 * t79;
t447 = t506 / 0.2e1;
t76 = -Ifges(6,3) * t249 + t348 * t593;
t337 = -t261 * (mrSges(5,1) * t348 + mrSges(5,2) * t249) - t321 * (-t303 * mrSges(4,1) + t361 * mrSges(4,2)) + t348 * t157 / 0.2e1 + (t249 * t593 + t519) * t625 - t361 ^ 2 * Ifges(4,4) + t249 * t447 + t611 * t488 / 0.2e1 + t249 * t445 - t610 * t487 / 0.2e1 + (Ifges(4,4) * t303 + (Ifges(4,1) - Ifges(4,2)) * t361) * t303 - (Ifges(5,1) * t249 - t529 + t76) * t348 / 0.2e1 + (-t333 * t150 - t336 * t153) * t119 + (-Ifges(5,2) * t348 + t158 + t246) * t626;
t340 = t321 * (-mrSges(8,1) * t302 + mrSges(8,2) * t301) - Ifges(8,4) * t302 ^ 2 + t609;
t479 = t333 * t151;
t395 = t336 * t154 + t479;
t406 = Ifges(8,4) * t301 + (-Ifges(8,1) + Ifges(8,2)) * t302;
t72 = -mrSges(9,1) * t419 + mrSges(9,2) * t176;
t1 = -t337 + t340 + (t395 + t618) * t120 + t596 * (t328 - t535) - pkin(12) * (t540 * mrSges(3,1) + t545 * mrSges(3,2)) + (-mrSges(8,1) * t328 + t406) * t301 + (-t361 * mrSges(4,1) - t303 * mrSges(4,2) - mrSges(8,2) * t302 + (m(4) + m(8)) * t321) * t328 + t221 * t72 + (m(9) * t221 + t71) * t214 + (-Ifges(3,2) + Ifges(3,1)) * t545 * t540 + (-t540 ^ 2 + t545 ^ 2) * Ifges(3,4);
t518 = t1 * qJD(1);
t365 = t336 * t376;
t366 = t333 * t376;
t2 = t151 * t366 + t154 * t365 + t376 * t618 + t596 * t535 + t337;
t517 = t2 * qJD(1);
t242 = (-t537 * t299 - t541 * t300) * pkin(1);
t516 = t242 * mrSges(9,2);
t515 = t278 * mrSges(5,2);
t514 = t279 * mrSges(5,1);
t513 = t283 * mrSges(5,1);
t512 = t284 * mrSges(5,2);
t3 = (t611 * t623 + t610 * t546 + t76 / 0.2e1 + t261 * mrSges(5,1) - t157 / 0.2e1 - t529 / 0.2e1) * t348 - (t447 + t445 - t261 * mrSges(5,2) - t246 / 0.2e1 - t158 / 0.2e1 - t581 + (t572 + Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t348) * t249 + t395 * t159 + (t333 * t149 + t336 * t152 + t159 * t425) * t119;
t511 = t3 * qJD(1);
t11 = t243 * t72 + (m(9) * t243 + t71) * t214 + t406 * t301 + t340;
t499 = t11 * qJD(1);
t17 = (t79 * t546 + t333 * t82 / 0.2e1 + t307 * t626 - t384 * t348) * t348 + (t477 - t473) * t119;
t496 = t17 * qJD(1);
t18 = t214 * t71 + t609;
t493 = t18 * qJD(1);
t491 = t348 * t334;
t486 = t273 * t616;
t485 = t278 * t249;
t484 = t279 * t348;
t405 = t503 - t507;
t483 = t279 * t405;
t482 = t283 * t405;
t480 = t325 * t616;
t160 = -pkin(6) * (t538 * mrSges(7,1) + t542 * mrSges(7,2)) + (t538 ^ 2 - t542 ^ 2) * Ifges(7,4) + (-Ifges(7,1) + Ifges(7,2)) * t542 * t538;
t468 = t160 * qJD(1);
t467 = t405 * t532;
t466 = -t543 / 0.2e1;
t461 = t197 * t419 * mrSges(9,3);
t460 = t199 * t176 * mrSges(9,3);
t145 = t405 * t348;
t459 = -pkin(8) * t145 / 0.2e1;
t458 = t154 * t578;
t457 = -t532 / 0.2e1;
t456 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t453 = Ifges(6,5) * t626;
t444 = t543 * t249;
t439 = t326 * t558;
t438 = t145 * t555;
t437 = t148 * t553;
t436 = t145 * t547;
t435 = t616 * t547;
t434 = -t479 / 0.2e1;
t431 = t334 * t569;
t428 = -t239 / 0.2e1 - t197 / 0.2e1;
t427 = t548 + t554;
t423 = t469 * t278;
t422 = t469 * t284;
t421 = t469 * t324;
t417 = -t464 / 0.2e1;
t416 = Ifges(8,5) * t301 + Ifges(8,6) * t302 + t643;
t410 = t543 * t567;
t409 = t333 * t466;
t387 = pkin(8) * t404;
t400 = -t387 / 0.2e1 + t415;
t355 = mrSges(6,3) * t423 + t483 + t514 - t515;
t31 = m(6) * (-t273 * t279 + t274 * t423) + t355;
t360 = t615 * pkin(10) - t580 + t638;
t385 = t279 * t569 - t555 * t616;
t5 = t634 + (t636 + t278 * t154 / 0.2e1 + t274 * t152 / 0.2e1 - t249 * t551) * t333 + (t149 * t554 + t278 * t568 + t309 * t558 + t635) * t336 + (Ifges(5,6) / 0.2e1 - t307 / 0.4e1) * t348 + t360 + t385;
t399 = -t5 * qJD(1) + t31 * qJD(2);
t353 = t348 * t553 - t249 * t552 - t485 / 0.2e1 + t484 / 0.2e1;
t16 = t437 + (t555 - t325 / 0.2e1) * t616 + (-t427 * t150 + t284 * t567) * t336 + (t427 * t153 + t154 * t552) * t333 + ((t444 / 0.2e1 + t491 / 0.2e1) * pkin(4) + t353) * mrSges(5,3);
t33 = (mrSges(5,1) + t405) * t283 + t584 + (-mrSges(5,2) + t620) * t284 + m(6) * (-t273 * t283 + t274 * t422) + m(5) * (t278 * t283 - t279 * t284);
t398 = t16 * qJD(1) + t33 * qJD(2);
t364 = pkin(1) * (-t541 * t299 + t537 * t300);
t341 = t461 / 0.2e1 + t460 / 0.2e1 + (-t242 * t419 / 0.2e1 + t364 * t561) * mrSges(9,3);
t20 = t341 + t583;
t231 = mrSges(9,1) * t364;
t357 = -mrSges(8,1) * t462 - mrSges(8,2) * t463 + t231 - t516;
t34 = -m(9) * (t197 * t364 + t199 * t242) - t357;
t397 = -t20 * qJD(1) - t34 * qJD(2);
t396 = t32 * qJD(2);
t393 = pkin(10) * t408;
t379 = t456 * t336 + t551;
t354 = (t379 - t528) * t348 + t570;
t347 = t154 * t554 + t354;
t378 = -t456 * t333 + t550;
t349 = (t624 + t625) * Ifges(6,6) + t378 * t348 + t571;
t386 = -t519 / 0.2e1 + t439;
t391 = t274 * t408;
t10 = t438 - t348 * t391 + (t120 * t573 + t151 * t554 + t349) * t333 + (t120 * t575 + t347 + t453) * t336 + t386;
t383 = t273 * t404;
t130 = t383 + t415;
t392 = t10 * qJD(1) + t130 * qJD(2);
t390 = t324 * t408;
t388 = t508 / 0.2e1 + t502 / 0.2e1;
t382 = t325 * t404;
t377 = -t629 + t638;
t375 = t387 / 0.2e1;
t374 = t469 * t543;
t373 = t388 * t278;
t372 = t388 * t284;
t371 = -t383 / 0.2e1;
t343 = -mrSges(5,1) * t532 - mrSges(5,2) * t464 - t467 + t582;
t131 = (t374 * t324 + t325 * t334) * t576 + t343;
t338 = (-t325 * t279 + (t273 * t334 + t374 * t274) * pkin(4)) * t577 - t515 / 0.2e1 + t514 / 0.2e1 + t483 / 0.2e1 + mrSges(5,1) * t457 - t467 / 0.2e1 + mrSges(5,2) * t417 + t582 / 0.2e1 + (t421 * t577 + t588) * t278;
t344 = (pkin(8) * t283 + pkin(10) * t422) * t577 + t513 / 0.2e1 + t482 / 0.2e1 - t512 / 0.2e1;
t23 = t588 * t284 - t338 + t344;
t7 = t435 + pkin(4) * t431 + (pkin(4) * t410 + t149 * t548 + t150 * t578 + t249 * t550 + t635) * t336 + (t636 + t154 * t417 + t308 * t624 + t152 * t549 + pkin(10) * t153 / 0.2e1) * t333 + t377 + t580;
t370 = t7 * qJD(1) - t23 * qJD(2) + t131 * qJD(3);
t352 = -t185 / 0.2e1 - t230 / 0.2e1;
t22 = -t231 / 0.2e1 + (t242 / 0.2e1 + t428) * mrSges(9,2) + t352;
t369 = t22 * qJD(2) - t35 * qJD(7);
t346 = t154 * t549 + t354;
t13 = t436 - t348 * t390 + (t151 * t549 + t376 * t574 + t349) * t333 + (t453 + t376 * mrSges(6,1) / 0.2e1 + t346) * t336 + t386;
t155 = t382 + t415;
t339 = -t382 / 0.2e1 - t415;
t57 = t371 - t372 + t339;
t368 = t13 * qJD(1) - t57 * qJD(2) + t155 * qJD(3);
t367 = t378 * t333;
t363 = Ifges(4,5) * t361 + Ifges(4,6) * t303 - t249 * t384 + t621 + t639;
t362 = t593 * t558 - t506 / 0.4e1;
t15 = t459 + t439 + (-Ifges(6,3) / 0.2e1 - t393) * t348 + (t458 + t570 - t525 / 0.2e1 + t159 * t575 + t379 * t348) * t336 + (0.3e1 / 0.4e1 * t522 + t571 + pkin(10) * t568 + t159 * t573 + (t378 - t527) * t348) * t333;
t161 = t387 - t415;
t59 = t371 + t375 - t373 - t415;
t358 = (mrSges(6,1) * t409 + t466 * t502) * pkin(4);
t91 = t375 + t358 + t339;
t359 = t15 * qJD(1) - t59 * qJD(2) - t91 * qJD(3) - t161 * qJD(4);
t356 = -Ifges(6,6) * t476 / 0.2e1 + Ifges(6,5) * t490 / 0.2e1 + t519 / 0.2e1 + t362;
t345 = t360 + t377;
t277 = t382 / 0.2e1;
t252 = t383 / 0.2e1;
t92 = t277 + t358 + t400;
t60 = t252 - t373 + t400;
t58 = t252 + t277 - t372 + t415;
t24 = t284 * t408 + t338 + t344;
t21 = t231 / 0.2e1 - t516 / 0.2e1 + t428 * mrSges(9,2) + t352;
t19 = t583 - t341 + t416;
t14 = t459 + pkin(10) * t434 + (t458 + t354) * t336 + t362 + (t572 - t393 + t367) * t348 - t581 + t589 * t159;
t12 = t436 + t324 * t434 + t365 * t575 - t366 * t574 + (-t390 + t367) * t348 + t346 * t336 + t356;
t9 = t274 * t434 + t438 + (-t391 + t367) * t348 + t347 * t336 + t356 + t589 * t120;
t8 = t363 + t480 / 0.2e1 + t437 + t486 / 0.2e1 + t591 * t284 + (t249 * t417 + t348 * t457 + t353) * mrSges(5,3) + (t324 + t274) * t615;
t6 = t345 + t435 + (t154 * t409 + t336 * t410 + t431) * pkin(4) + t592 * t324;
t4 = t592 * t274 + t591 * t278 + t345 - t385;
t25 = [qJD(2) * t1 - qJD(3) * t2 + qJD(4) * t3 - qJD(5) * t17 - qJD(6) * t160 + qJD(7) * t11 + qJD(8) * t18, t8 * qJD(3) + t4 * qJD(4) + t9 * qJD(5) + t19 * qJD(7) + t642 + t518 + (Ifges(3,5) * t545 - Ifges(3,6) * t540 + t363 + t416 - t460 - t461 + t486 + ((-t541 * t301 + t537 * t302) * mrSges(8,3) + (-t539 * t303 + t544 * t361) * mrSges(4,3)) * pkin(1) + t587 * t274 + (t484 - t485) * mrSges(5,3)) * qJD(2), -t517 + t8 * qJD(2) + (t480 + t587 * t324 + (-t444 - t491) * mrSges(5,3) * pkin(4) + t363) * qJD(3) + t6 * qJD(4) + t12 * qJD(5), t4 * qJD(2) + t6 * qJD(3) + t14 * qJD(5) + t511 + (-t633 - (-Ifges(5,5) + t384) * t249 + (-t478 + t474) * pkin(10) + t639) * qJD(4), -t496 + t9 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + (-t119 * t404 - t597) * qJD(5), -t468 + (Ifges(7,5) * t542 - Ifges(7,6) * t538) * qJD(6), t499 + t19 * qJD(2) + ((-t494 - t495) * mrSges(9,3) + t416) * qJD(7) + t642, t642 + t493 + (qJD(2) + qJD(7)) * t643, 0, 0; qJD(3) * t16 - qJD(4) * t5 + qJD(5) * t10 - qJD(7) * t20 - t518, qJD(3) * t33 + qJD(4) * t31 + qJD(5) * t130 - qJD(7) * t34 - t619, (t482 + m(6) * (-t325 * t283 + t284 * t421) + t513 - t512 + m(5) * (t283 * t543 + t284 * t334) * pkin(4) + t284 * t510 + t284 * t509 + t584) * qJD(3) + t24 * qJD(4) + t58 * qJD(5) + t398, t24 * qJD(3) + (m(6) * (pkin(8) * t279 + pkin(10) * t423) + t355) * qJD(4) + t60 * qJD(5) + t399, t58 * qJD(3) + t60 * qJD(4) + (-t274 * t405 + t593) * qJD(5) + t392, 0, (m(9) * (t239 * t364 + t242 * t606) + t357) * qJD(7) + t21 * qJD(8) + t397, t21 * qJD(7) - t396 - t619, 0, 0; -qJD(2) * t16 + qJD(4) * t7 + qJD(5) * t13 + t517, -qJD(4) * t23 - qJD(5) * t57 - t398, qJD(4) * t131 + qJD(5) * t155, ((-pkin(8) * t334 + pkin(10) * t374) * t576 + t343) * qJD(4) + t92 * qJD(5) + t370, t92 * qJD(4) + (-t324 * t405 + t593) * qJD(5) + t368, 0, 0, 0, 0, 0; qJD(2) * t5 - qJD(3) * t7 + qJD(5) * t15 - t511, qJD(3) * t23 - qJD(5) * t59 - t399, -qJD(5) * t91 - t370, -t161 * qJD(5), (-pkin(10) * t405 + t593) * qJD(5) + t359, 0, 0, 0, 0, 0; -qJD(2) * t10 - qJD(3) * t13 - qJD(4) * t15 + t496, qJD(3) * t57 + qJD(4) * t59 - t392, qJD(4) * t91 - t368, -t359, 0, 0, 0, 0, 0, 0; t468, 0, 0, 0, 0, 0, 0, 0, 0, 0; qJD(2) * t20 - t499, qJD(8) * t22 - t397, 0, 0, 0, 0, -t617, t369 - t617, 0, 0; -t493, -qJD(7) * t22 + t396, 0, 0, 0, 0, -t369, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t25;
