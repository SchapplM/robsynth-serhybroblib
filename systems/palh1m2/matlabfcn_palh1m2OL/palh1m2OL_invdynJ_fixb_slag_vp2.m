% Calculate vector of inverse dynamics joint torques for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% qJDD [13x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% tau [13x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m2OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:44
% EndTime: 2020-05-02 21:20:17
% DurationCPUTime: 38.01s
% Computational Cost: add. (11027->932), mult. (25144->1353), div. (0->0), fcn. (20325->30), ass. (0->426)
t403 = sin(qJ(2));
t411 = cos(qJ(3));
t530 = t403 * t411;
t402 = sin(qJ(3));
t412 = cos(qJ(2));
t532 = t402 * t412;
t309 = t530 + t532;
t285 = t309 * qJD(1);
t527 = t411 * t412;
t310 = -t402 * t403 + t527;
t288 = t310 * qJD(1);
t401 = sin(qJ(4));
t410 = cos(qJ(4));
t479 = -t285 * t401 + t410 * t288;
t194 = Ifges(5,4) * t479;
t453 = t410 * t285 + t288 * t401;
t507 = qJD(3) + qJD(4);
t380 = qJD(2) + t507;
t577 = Ifges(5,5) * t380;
t101 = Ifges(5,1) * t453 + t194 + t577;
t534 = t401 * t411;
t354 = pkin(1) * t534;
t363 = pkin(1) * t402 + pkin(5);
t538 = t363 * t410;
t273 = t354 + t538;
t590 = pkin(5) * qJD(3);
t495 = t410 * t590;
t243 = qJD(2) * t273 + t495;
t227 = -pkin(9) * t380 - t243;
t400 = sin(qJ(5));
t409 = cos(qJ(5));
t324 = mrSges(6,1) * t400 + mrSges(6,2) * t409;
t436 = t227 * t324;
t155 = t380 * t409 - t400 * t453;
t195 = qJD(5) - t479;
t156 = t380 * t400 + t409 * t453;
t575 = t156 * Ifges(6,4);
t46 = t155 * Ifges(6,2) + t195 * Ifges(6,6) + t575;
t154 = Ifges(6,4) * t155;
t47 = t156 * Ifges(6,1) + t195 * Ifges(6,5) + t154;
t560 = t409 * t47;
t610 = t400 / 0.2e1;
t674 = -t560 / 0.2e1 + t46 * t610 - t101 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t453 - t194 / 0.2e1 - t436 - t577 / 0.2e1;
t381 = -qJ(7) + pkin(19) - qJ(2);
t341 = pkin(4) * sin(t381);
t386 = t403 * pkin(1);
t358 = -qJ(10) + t381;
t342 = sin(t358);
t343 = cos(t358);
t482 = -t342 * mrSges(11,1) + t343 * mrSges(11,2);
t673 = m(11) * (t341 - t386) + t482;
t632 = m(5) + m(6);
t484 = -mrSges(6,1) * t409 + t400 * mrSges(6,2);
t304 = m(6) * pkin(9) + mrSges(5,1) - t484;
t377 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t651 = -t304 * t410 - t377 * t401;
t225 = pkin(5) * t632 + mrSges(4,1) - t651;
t240 = -t304 * t401 + t377 * t410;
t238 = -mrSges(4,2) + t240;
t121 = t225 * t411 + t238 * t402;
t120 = -t225 * t402 + t238 * t411;
t573 = t195 * Ifges(6,3);
t574 = t156 * Ifges(6,5);
t576 = t155 * Ifges(6,6);
t45 = t573 + t574 + t576;
t565 = t380 * Ifges(5,6);
t528 = t410 * t411;
t353 = pkin(1) * t528;
t272 = t363 * t401 - t353;
t241 = t272 * qJD(2) + t401 * t590;
t570 = t241 * mrSges(5,3);
t582 = Ifges(5,4) * t453;
t226 = pkin(11) * t380 + t241;
t362 = pkin(5) * t402 + pkin(1);
t263 = -pkin(5) * t527 + t362 * t403 - pkin(15);
t523 = qJD(1) * t263;
t86 = -pkin(9) * t479 - pkin(11) * t453 + t523;
t58 = t226 * t409 + t400 * t86;
t657 = t58 * mrSges(6,2);
t57 = -t226 * t400 + t409 * t86;
t658 = t57 * mrSges(6,1);
t99 = Ifges(5,2) * t479 + t565 + t582;
t672 = t657 + t570 + t99 / 0.2e1 - t45 / 0.2e1 + t582 / 0.2e1 + t565 / 0.2e1 - t658;
t631 = m(8) + m(4);
t535 = t401 * t402;
t500 = pkin(1) * t535;
t517 = qJD(4) * t401;
t213 = -qJD(3) * t500 + t353 * t507 - t363 * t517;
t117 = -pkin(5) * (qJD(3) * t517 - qJDD(3) * t410) + qJD(2) * t213 + qJDD(2) * t273;
t389 = qJDD(2) + qJDD(3);
t376 = qJDD(4) + t389;
t113 = -pkin(9) * t376 - t117;
t511 = qJD(1) * qJD(2);
t318 = qJDD(1) * t412 - t403 * t511;
t487 = t412 * t511;
t319 = -qJDD(1) * t403 - t487;
t427 = t309 * qJD(3);
t169 = -qJD(1) * t427 + t318 * t411 + t319 * t402;
t426 = t310 * qJD(3);
t170 = qJD(1) * t426 + t318 * t402 - t319 * t411;
t55 = qJD(4) * t479 + t169 * t401 + t170 * t410;
t17 = qJD(5) * t155 + t376 * t400 + t409 * t55;
t18 = -qJD(5) * t156 + t376 * t409 - t400 * t55;
t8 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t669 = m(6) * t113 + t8;
t668 = t58 * t400;
t352 = pkin(5) * t530;
t293 = pkin(5) * t532 + t352;
t667 = mrSges(11,1) * t343 + mrSges(11,2) * t342;
t398 = sin(qJ(7));
t407 = cos(qJ(7));
t307 = -t398 * t412 - t403 * t407;
t284 = t307 * qJD(1);
t520 = qJD(1) * t412;
t521 = qJD(1) * t403;
t287 = -t398 * t521 + t407 * t520;
t394 = sin(pkin(19));
t395 = cos(pkin(19));
t592 = sin(qJ(10));
t593 = cos(qJ(10));
t299 = t394 * t592 + t395 * t593;
t424 = -t394 * t593 + t395 * t592;
t480 = -t299 * t284 + t287 * t424;
t125 = Ifges(11,4) * t480;
t429 = t307 * qJD(7);
t167 = qJD(1) * t429 + t318 * t407 + t319 * t398;
t450 = t398 * t403 - t407 * t412;
t428 = t450 * qJD(7);
t168 = qJD(1) * t428 - t318 * t398 + t319 * t407;
t269 = t299 * qJD(10);
t270 = t424 * qJD(10);
t25 = -t167 * t299 - t168 * t424 - t269 * t284 + t270 * t287;
t26 = t167 * t424 - t168 * t299 + t269 * t287 + t270 * t284;
t388 = qJDD(2) + qJDD(7);
t370 = qJDD(10) + t388;
t591 = pkin(1) * qJD(2);
t494 = qJD(7) * t591;
t556 = pkin(1) * qJDD(2);
t300 = -t398 * t494 + t407 * t556;
t604 = pkin(4) * t395;
t251 = t388 * t604 + t300;
t301 = t398 * t556 + t407 * t494;
t605 = pkin(4) * t394;
t252 = t388 * t605 + t301;
t391 = qJD(2) + qJD(7);
t498 = t398 * t591;
t294 = t391 * t605 + t498;
t497 = t407 * t591;
t295 = t391 * t604 + t497;
t60 = -t251 * t424 - t252 * t299 - t269 * t295 + t270 * t294;
t61 = -t251 * t299 + t252 * t424 + t269 * t294 + t270 * t295;
t628 = -t480 / 0.2e1;
t130 = -t284 * t424 - t287 * t299;
t378 = qJD(10) + t391;
t64 = t130 * Ifges(11,1) + t378 * Ifges(11,5) + t125;
t665 = t61 * mrSges(11,1) - t60 * mrSges(11,2) + Ifges(11,5) * t25 + Ifges(11,6) * t26 + Ifges(11,3) * t370 + (t64 + t125) * t628;
t613 = -t378 / 0.2e1;
t627 = -t130 / 0.2e1;
t664 = Ifges(11,1) * t627 + Ifges(11,5) * t613;
t135 = t294 * t299 + t295 * t424;
t549 = t135 * mrSges(11,3);
t63 = t130 * Ifges(11,4) + Ifges(11,2) * t480 + t378 * Ifges(11,6);
t663 = Ifges(11,4) * t627 + Ifges(11,2) * t628 + Ifges(11,6) * t613 + t549 - t63 / 0.2e1;
t397 = sin(qJ(8));
t406 = cos(qJ(8));
t305 = -t397 * t412 - t403 * t406;
t283 = t305 * qJD(1);
t286 = -t397 * t521 + t406 * t520;
t396 = sin(qJ(9));
t405 = cos(qJ(9));
t455 = -t283 * t405 + t286 * t396;
t193 = Ifges(10,4) * t455;
t390 = qJD(2) + qJD(8);
t379 = qJD(9) + t390;
t454 = t283 * t396 + t286 * t405;
t100 = -Ifges(10,1) * t454 + t379 * Ifges(10,5) + t193;
t553 = qJD(1) * pkin(15);
t253 = -pkin(2) * t283 - t553;
t387 = qJDD(2) + qJDD(8);
t514 = qJD(9) * t396;
t261 = (-t387 * t405 + t390 * t514) * pkin(2);
t262 = (-qJD(9) * t390 * t405 - t387 * t396) * pkin(2);
t393 = qJ(2) + qJ(8);
t383 = qJ(9) + t393;
t360 = sin(t383);
t361 = cos(t383);
t375 = qJDD(9) + t387;
t404 = sin(qJ(1));
t413 = cos(qJ(1));
t473 = g(1) * t413 + g(2) * t404;
t431 = t305 * qJD(8);
t165 = qJD(1) * t431 + t318 * t406 + t319 * t397;
t451 = t397 * t403 - t406 * t412;
t430 = t451 * qJD(8);
t166 = qJD(1) * t430 - t318 * t397 + t319 * t406;
t53 = qJD(9) * t455 - t165 * t405 - t166 * t396;
t54 = qJD(9) * t454 + t165 * t396 - t166 * t405;
t571 = t454 * Ifges(10,4);
t621 = -t454 / 0.2e1;
t98 = Ifges(10,2) * t455 + t379 * Ifges(10,6) - t571;
t661 = -t262 * mrSges(10,2) + Ifges(10,5) * t53 + Ifges(10,6) * t54 - (g(3) * mrSges(10,1) - mrSges(10,2) * t473) * t360 + (-mrSges(10,1) * t473 - g(3) * mrSges(10,2)) * t361 + t261 * mrSges(10,1) + Ifges(10,3) * t375 - (Ifges(10,5) * t455 + Ifges(10,6) * t454) * t379 / 0.2e1 - t253 * (-mrSges(10,1) * t454 + mrSges(10,2) * t455) + t98 * t621 + (Ifges(10,1) * t455 + t571) * t454 / 0.2e1 - (Ifges(10,2) * t454 + t100 + t193) * t455 / 0.2e1;
t660 = m(9) + m(3);
t399 = sin(qJ(6));
t611 = t399 / 0.2e1;
t609 = t412 / 0.2e1;
t659 = g(3) * t482;
t441 = pkin(1) * (t299 * t398 + t407 * t424);
t654 = -qJD(2) * t441 + (t269 * t394 + t270 * t395) * pkin(4);
t440 = pkin(1) * (-t299 * t407 + t398 * t424);
t653 = -qJD(2) * t440 + (-t269 * t395 + t270 * t394) * pkin(4);
t652 = -mrSges(4,1) * t288 - mrSges(8,1) * t284 + mrSges(4,2) * t285 + mrSges(8,2) * t287;
t224 = t309 * t410 + t310 * t401;
t515 = qJD(5) * t409;
t223 = -t401 * t309 + t310 * t410;
t234 = qJD(2) * t310 + t426;
t235 = -qJD(2) * t309 - t427;
t85 = qJD(4) * t223 + t234 * t410 + t235 * t401;
t445 = t224 * t515 + t400 * t85;
t649 = -t400 * t57 + t409 * t58;
t648 = t631 + t632;
t647 = -mrSges(11,3) - mrSges(10,3) - t324 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t408 = cos(qJ(6));
t470 = -t408 * mrSges(7,1) + t399 * mrSges(7,2);
t645 = pkin(14) * m(7) + t470;
t646 = -mrSges(2,1) + t645 + (-t648 - t660 - m(11) - m(10)) * pkin(15) - t673;
t56 = -qJD(4) * t453 + t169 * t410 - t170 * t401;
t583 = Ifges(3,4) * t412;
t584 = Ifges(3,4) * t403;
t644 = pkin(15) * (mrSges(3,1) * t412 - mrSges(3,2) * t403) - t412 * (-Ifges(3,1) * t403 - t583) / 0.2e1 + t403 * (-Ifges(3,2) * t412 - t584) / 0.2e1;
t518 = qJD(2) * t412;
t214 = qJD(2) * t352 + qJD(3) * t293 + t362 * t518;
t209 = t214 * qJD(1);
t131 = t263 * qJDD(1) + t209;
t11 = -pkin(9) * t56 - pkin(11) * t55 + t131;
t533 = t402 * t410;
t449 = t533 + t534;
t212 = qJD(4) * t538 + (qJD(3) * t449 + t411 * t517) * pkin(1);
t603 = pkin(5) * t401;
t116 = t212 * qJD(2) + qJD(4) * t495 + t272 * qJDD(2) + qJDD(3) * t603;
t112 = pkin(11) * t376 + t116;
t2 = qJD(5) * t57 + t11 * t400 + t112 * t409;
t552 = qJD(5) * t58;
t3 = t11 * t409 - t112 * t400 - t552;
t643 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t642 = pkin(9) * t453 - pkin(11) * t479;
t640 = pkin(2) * m(10);
t639 = t17 / 0.2e1;
t638 = t18 / 0.2e1;
t52 = qJDD(5) - t56;
t635 = t52 / 0.2e1;
t626 = -t155 / 0.2e1;
t625 = -t156 / 0.2e1;
t624 = t156 / 0.2e1;
t623 = -t195 / 0.2e1;
t617 = t285 / 0.2e1;
t616 = t286 / 0.2e1;
t615 = t287 / 0.2e1;
t608 = pkin(1) * t412;
t607 = pkin(2) * (-mrSges(10,1) * t455 - mrSges(10,2) * t454);
t606 = pkin(4) * cos(t381);
t598 = t2 * t409;
t597 = t3 * t400;
t596 = -qJD(2) / 0.2e1;
t594 = -qJD(6) / 0.2e1;
t589 = m(11) * (-t606 - t608);
t587 = mrSges(6,3) * t409;
t586 = mrSges(8,3) * t287;
t585 = mrSges(10,3) * t390;
t581 = Ifges(6,4) * t400;
t580 = Ifges(6,4) * t409;
t579 = Ifges(7,4) * t399;
t578 = Ifges(7,4) * t408;
t572 = t479 * mrSges(5,3);
t569 = t243 * mrSges(5,3);
t568 = t285 * Ifges(4,4);
t567 = t286 * Ifges(9,4);
t566 = t287 * Ifges(8,4);
t564 = t399 * Ifges(7,1);
t561 = t408 * Ifges(7,2);
t558 = t57 * t409;
t557 = mrSges(5,1) * t380 + mrSges(6,1) * t155 - mrSges(6,2) * t156 - mrSges(5,3) * t453;
t359 = t386 - pkin(15);
t551 = t480 * mrSges(11,3);
t134 = t294 * t424 - t295 * t299;
t550 = t134 * mrSges(11,3);
t548 = t224 * t400;
t547 = t224 * t409;
t545 = t263 * (mrSges(5,1) * t453 + mrSges(5,2) * t479);
t414 = qJD(1) ^ 2;
t544 = t263 * t414;
t464 = t561 + t579;
t536 = t399 * (Ifges(7,6) * qJD(6) + qJD(1) * t464);
t468 = t412 * Ifges(3,1) - t584;
t531 = t403 * (Ifges(3,5) * qJD(2) + qJD(1) * t468);
t366 = qJD(1) * t578;
t529 = t408 * (Ifges(7,5) * qJD(6) + qJD(1) * t564 + t366);
t526 = t667 * t404;
t525 = t667 * t413;
t349 = pkin(1) * t487;
t276 = t359 * qJDD(1) + t349;
t522 = qJD(1) * t359;
t292 = pkin(1) * t533 + t354;
t519 = qJD(2) * t292;
t516 = qJD(5) * t400;
t513 = qJDD(1) * pkin(15);
t510 = qJD(1) * qJD(6);
t509 = qJD(2) * qJD(3);
t506 = t253 * t640;
t505 = sin(t393) * t640;
t504 = pkin(2) * t585;
t503 = m(11) * t606;
t502 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t52;
t501 = mrSges(4,3) * t591;
t496 = pkin(1) * t518;
t491 = t560 / 0.2e1;
t486 = -t516 / 0.2e1;
t483 = t396 * t504;
t477 = t402 * t501;
t476 = t411 * t501;
t475 = mrSges(8,3) * t497;
t469 = mrSges(7,1) * t399 + mrSges(7,2) * t408;
t326 = -t407 * mrSges(8,1) + mrSges(8,2) * t398;
t323 = -mrSges(8,1) * t398 - mrSges(8,2) * t407;
t325 = -t406 * mrSges(9,1) + t397 * mrSges(9,2);
t322 = -mrSges(9,1) * t397 - mrSges(9,2) * t406;
t467 = Ifges(6,1) * t409 - t581;
t466 = -t403 * Ifges(3,2) + t583;
t465 = -Ifges(6,2) * t400 + t580;
t463 = -Ifges(3,5) * t403 - Ifges(3,6) * t412;
t462 = Ifges(6,5) * t409 - Ifges(6,6) * t400;
t461 = Ifges(7,5) * t408 - Ifges(7,6) * t399;
t460 = t558 + t668;
t228 = t449 * t324;
t229 = (-t528 + t535) * t324;
t457 = -t228 * t403 - t229 * t412;
t392 = qJD(2) + qJD(3);
t456 = (mrSges(4,1) * t392 - mrSges(4,3) * t285) * t411 + (-mrSges(4,2) * t392 + mrSges(4,3) * t288) * t402;
t221 = -t305 * t405 - t396 * t451;
t452 = t305 * t396 - t405 * t451;
t444 = t224 * t516 - t409 * t85;
t443 = pkin(14) * t469;
t439 = t155 * t465;
t438 = t156 * t467;
t437 = t195 * t462;
t435 = t399 * (Ifges(7,1) * t408 - t579);
t182 = (-t284 * t394 + t287 * t395) * pkin(4);
t423 = -qJD(5) * t460 - t597;
t6 = t17 * Ifges(6,4) + t18 * Ifges(6,2) + t52 * Ifges(6,6);
t7 = t17 * Ifges(6,1) + t18 * Ifges(6,4) + t52 * Ifges(6,5);
t420 = -t116 * mrSges(5,2) + t2 * t587 + t409 * t6 / 0.2e1 + t7 * t610 + t117 * mrSges(5,1) + (Ifges(6,1) * t400 + t580) * t639 + (Ifges(6,2) * t409 + t581) * t638 + Ifges(5,3) * t376 + (Ifges(6,5) * t400 + Ifges(6,6) * t409) * t635 + t46 * t486 + Ifges(5,6) * t56 + Ifges(5,5) * t55 + t113 * t484 + (t436 + t491) * qJD(5) + (t439 + t438 + t437) * qJD(5) / 0.2e1;
t10 = -mrSges(6,2) * t52 + mrSges(6,3) * t18;
t9 = mrSges(6,1) * t52 - mrSges(6,3) * t17;
t90 = -mrSges(6,2) * t195 + mrSges(6,3) * t155;
t91 = mrSges(6,1) * t195 - mrSges(6,3) * t156;
t419 = m(6) * (-t515 * t57 - t516 * t58 - t597 + t598) + t409 * t10 - t400 * t9 - t91 * t515 - t90 * t516;
t188 = t284 * Ifges(8,2) + t391 * Ifges(8,6) + t566;
t267 = Ifges(8,4) * t284;
t191 = t287 * Ifges(8,1) + t391 * Ifges(8,5) + t267;
t418 = -t301 * mrSges(8,2) + t188 * t615 - t287 * (Ifges(8,1) * t284 - t566) / 0.2e1 + Ifges(8,6) * t168 + Ifges(8,5) * t167 - t391 * (Ifges(8,5) * t284 - Ifges(8,6) * t287) / 0.2e1 + t284 * t475 + t300 * mrSges(8,1) + Ifges(8,3) * t388 - (-Ifges(8,2) * t287 + t191 + t267) * t284 / 0.2e1 + (t550 + t664) * t480 - t663 * t130 + t665;
t173 = -mrSges(10,2) * t379 + mrSges(10,3) * t455;
t175 = mrSges(10,1) * t379 + mrSges(10,3) * t454;
t187 = t283 * Ifges(9,2) + t390 * Ifges(9,6) + t567;
t266 = Ifges(9,4) * t283;
t190 = t286 * Ifges(9,1) + t390 * Ifges(9,5) + t266;
t417 = t454 * t483 - t286 * (Ifges(9,1) * t283 - t567) / 0.2e1 - t286 * t607 + t187 * t616 + (mrSges(9,1) * t286 + mrSges(9,2) * t283) * t553 + Ifges(9,5) * t165 + Ifges(9,6) * t166 - t286 * t506 + g(3) * t505 + Ifges(9,3) * t387 - t390 * (Ifges(9,5) * t283 - Ifges(9,6) * t286) / 0.2e1 + (t473 * cos(t393) - t261 * t405 - t262 * t396) * t640 - (-Ifges(9,2) * t286 + t190 + t266) * t283 / 0.2e1 + (-t396 * (-mrSges(10,2) * t375 + mrSges(10,3) * t54) + (-mrSges(10,1) * t375 + mrSges(10,3) * t53 - qJD(9) * t173 - t455 * t585) * t405 + t175 * t514) * pkin(2) + t661;
t189 = t288 * Ifges(4,2) + t392 * Ifges(4,6) + t568;
t268 = Ifges(4,4) * t288;
t192 = t285 * Ifges(4,1) + t392 * Ifges(4,5) + t268;
t302 = (-qJDD(2) * t411 + t402 * t509) * pkin(1);
t303 = (qJDD(2) * t402 + t411 * t509) * pkin(1);
t416 = t420 - t285 * (Ifges(4,1) * t288 - t568) / 0.2e1 + t189 * t617 - t302 * mrSges(4,2) + t303 * mrSges(4,1) + Ifges(4,5) * t170 + Ifges(4,6) * t169 - t285 * t476 + t288 * t477 + Ifges(4,3) * t389 - t392 * (Ifges(4,5) * t288 - Ifges(4,6) * t285) / 0.2e1 - (-Ifges(4,2) * t285 + t192 + t268) * t288 / 0.2e1 + (Ifges(6,5) * t625 + Ifges(6,6) * t626 + Ifges(6,3) * t623 + t672) * t453 + (mrSges(6,3) * t668 + t462 * t623 + t465 * t626 + t467 * t625 + t57 * t587 + t569 + t674) * t479;
t328 = pkin(1) * t407 + t604;
t327 = pkin(1) * t398 + t605;
t317 = qJDD(1) * t399 + t408 * t510;
t316 = qJDD(1) * t408 - t399 * t510;
t291 = t353 - t500;
t280 = Ifges(3,6) * qJD(2) + qJD(1) * t466;
t271 = t362 * t412 + t352;
t265 = -pkin(9) - t273;
t264 = pkin(11) + t272;
t246 = mrSges(8,1) * t391 - t586;
t244 = -mrSges(8,2) * t391 + mrSges(8,3) * t284;
t233 = qJD(2) * t450 + t428;
t232 = qJD(2) * t307 + t429;
t231 = qJD(2) * t451 + t430;
t230 = qJD(2) * t305 + t431;
t218 = mrSges(4,1) * t285 + mrSges(4,2) * t288;
t215 = mrSges(8,1) * t287 + mrSges(8,2) * t284;
t211 = (-t299 * t394 - t395 * t424) * pkin(4);
t210 = (-t299 * t395 + t394 * t424) * pkin(4);
t186 = (-t307 * t395 + t394 * t450) * pkin(4) + t359;
t178 = pkin(1) * t520 + t182;
t174 = -mrSges(5,2) * t380 + t572;
t172 = -t299 * t327 - t328 * t424;
t171 = -t299 * t328 + t327 * t424;
t164 = t522 + (-t284 * t395 - t287 * t394) * pkin(4);
t145 = -pkin(2) * t166 - t513;
t133 = t240 * t411 + t402 * t651;
t132 = t240 * t402 - t411 * t651;
t119 = mrSges(11,1) * t378 - mrSges(11,3) * t130;
t118 = -mrSges(11,2) * t378 + t551;
t114 = t496 + (-t232 * t394 - t233 * t395) * pkin(4);
t109 = -mrSges(5,1) * t479 + mrSges(5,2) * t453;
t104 = -mrSges(3,2) + t121 + t322 + t323;
t94 = pkin(1) * t648 + mrSges(3,1) - t120 - t325 - t326;
t93 = qJD(7) * t441 + t269 * t327 + t270 * t328;
t92 = qJD(7) * t440 - t269 * t328 + t270 * t327;
t89 = qJD(1) * t293 + t642;
t87 = qJD(1) * t271 + t642;
t84 = qJD(4) * t224 + t234 * t401 - t410 * t235;
t83 = qJD(9) * t452 + t230 * t396 - t231 * t405;
t82 = qJD(9) * t221 - t230 * t405 - t231 * t396;
t80 = t243 * t409 + t400 * t642;
t79 = -t243 * t400 + t409 * t642;
t78 = t400 * t89 + t409 * t519;
t77 = -t400 * t519 + t409 * t89;
t71 = (-t167 * t394 - t168 * t395) * pkin(4) + t276;
t70 = -mrSges(11,1) * t480 + mrSges(11,2) * t130;
t69 = mrSges(11,1) * t130 + mrSges(11,2) * t480;
t43 = t232 * t424 - t233 * t299 - t269 * t450 + t270 * t307;
t42 = -t232 * t299 - t233 * t424 - t269 * t307 - t270 * t450;
t39 = -mrSges(5,2) * t376 + mrSges(5,3) * t56;
t38 = mrSges(5,1) * t376 - mrSges(5,3) * t55;
t20 = -mrSges(11,2) * t370 + mrSges(11,3) * t26;
t19 = mrSges(11,1) * t370 - mrSges(11,3) * t25;
t1 = [(mrSges(10,1) * t360 + mrSges(10,2) * t361 + t104 * t412 - t403 * t94 - t505) * (g(1) * t404 - g(2) * t413) + (t408 * (-Ifges(7,2) * t399 + t578) + t435) * t510 / 0.2e1 + t480 * (Ifges(11,4) * t42 + Ifges(11,2) * t43) / 0.2e1 + t455 * (Ifges(10,4) * t82 + Ifges(10,2) * t83) / 0.2e1 + t479 * (Ifges(5,4) * t85 - Ifges(5,2) * t84) / 0.2e1 - t644 * t511 + (t529 / 0.2e1 + t461 * qJD(6) / 0.2e1) * qJD(6) + (pkin(15) ^ 2 * t660 + pkin(14) * t645 + Ifges(2,3)) * qJDD(1) + (t209 + t131) * m(5) * t263 + (-mrSges(4,1) * t235 - mrSges(8,1) * t233 + mrSges(4,2) * t234 + mrSges(8,2) * t232) * t522 + t408 * (Ifges(7,4) * t317 + Ifges(7,2) * t316) / 0.2e1 + (t631 * (t276 + t349) - mrSges(4,1) * t169 - mrSges(8,1) * t168 + mrSges(4,2) * t170 + mrSges(8,2) * t167) * t359 - t403 * (Ifges(3,4) * t318 + Ifges(3,2) * t319) / 0.2e1 + (g(1) * t647 + g(2) * t646) * t413 + (-g(1) * t646 + g(2) * t647) * t404 + (mrSges(4,2) * t276 - mrSges(4,3) * t303 + Ifges(4,1) * t170 + Ifges(4,4) * t169 + Ifges(4,5) * t389) * t309 + (-mrSges(8,1) * t276 + mrSges(8,3) * t301 + Ifges(8,4) * t167 + Ifges(8,2) * t168 + Ifges(8,6) * t388) * t307 + (mrSges(9,1) * t513 + Ifges(9,4) * t165 + Ifges(9,2) * t166 + Ifges(9,6) * t387) * t305 + (mrSges(3,1) * t319 + mrSges(9,1) * t166 - mrSges(3,2) * t318 - mrSges(9,2) * t165) * pkin(15) + (-Ifges(6,1) * t444 - Ifges(6,4) * t445 + Ifges(6,5) * t84) * t624 + (t131 * mrSges(5,2) - t117 * mrSges(5,3) + Ifges(5,1) * t55 + Ifges(5,4) * t56 + Ifges(5,5) * t376 + t113 * t324 + t462 * t635 + t465 * t638 + t467 * t639 + t47 * t486) * t224 + (-Ifges(6,3) * t635 - Ifges(6,6) * t638 - Ifges(6,5) * t639 + Ifges(5,6) * t376 + Ifges(5,4) * t55 + Ifges(5,2) * t56 - t502 / 0.2e1 - t131 * mrSges(5,1) + t116 * mrSges(5,3) + t643) * t223 + t453 * (Ifges(5,1) * t85 - Ifges(5,4) * t84) / 0.2e1 + t317 * (t564 + t578) / 0.2e1 - t85 * t569 - t84 * t570 + t7 * t547 / 0.2e1 - t6 * t548 / 0.2e1 - t280 * t518 / 0.2e1 - t231 * t607 + (-mrSges(4,1) * t276 + mrSges(4,3) * t302 + Ifges(4,4) * t170 + Ifges(4,2) * t169 + Ifges(4,6) * t389) * t310 + (-t2 * t548 - t3 * t547 + t444 * t57 - t445 * t58) * mrSges(6,3) + t443 * t510 - (-mrSges(9,1) * t231 + mrSges(9,2) * t230) * t553 - t42 * t550 + (Ifges(8,1) * t232 + Ifges(8,4) * t233) * t615 + (Ifges(9,1) * t230 + Ifges(9,4) * t231) * t616 + (Ifges(4,1) * t234 + Ifges(4,4) * t235) * t617 + (Ifges(10,1) * t82 + Ifges(10,4) * t83) * t621 + (-mrSges(10,1) * t145 + mrSges(10,3) * t262 + Ifges(10,4) * t53 + Ifges(10,2) * t54 + Ifges(10,6) * t375) * t221 - t43 * t549 + t378 * (Ifges(11,5) * t42 + Ifges(11,6) * t43) / 0.2e1 + t379 * (Ifges(10,5) * t82 + Ifges(10,6) * t83) / 0.2e1 - (mrSges(8,2) * t276 - mrSges(8,3) * t300 + Ifges(8,1) * t167 + Ifges(8,4) * t168 + Ifges(8,5) * t388) * t450 + (-mrSges(11,1) * t71 + mrSges(11,3) * t60 + Ifges(11,4) * t25 + Ifges(11,2) * t26 + Ifges(11,6) * t370) * (-t299 * t307 - t424 * t450) + (mrSges(11,2) * t71 - mrSges(11,3) * t61 + Ifges(11,1) * t25 + Ifges(11,4) * t26 + Ifges(11,5) * t370) * (t299 * t450 - t307 * t424) - (-mrSges(9,2) * t513 + Ifges(9,1) * t165 + Ifges(9,4) * t166 + Ifges(9,5) * t387) * t451 + (Ifges(7,1) * t317 + Ifges(7,4) * t316) * t611 - (mrSges(3,1) * t403 + mrSges(3,2) * t412) * t513 + (Ifges(3,1) * t318 + Ifges(3,4) * t319) * t609 + (m(6) * t460 + t400 * t90 + t409 * t91) * (pkin(9) * t84 - pkin(11) * t85 + t214) - (mrSges(10,2) * t145 - mrSges(10,3) * t261 + Ifges(10,1) * t53 + Ifges(10,4) * t54 + Ifges(10,5) * t375) * t452 + (m(10) * t145 - mrSges(10,1) * t54 + mrSges(10,2) * t53) * (-pkin(2) * t305 - pkin(15)) + pkin(14) * (-mrSges(7,1) * t316 + mrSges(7,2) * t317) - t231 * t506 + t288 * (Ifges(4,4) * t234 + Ifges(4,2) * t235) / 0.2e1 + t284 * (Ifges(8,4) * t232 + Ifges(8,2) * t233) / 0.2e1 + t283 * (Ifges(9,4) * t230 + Ifges(9,2) * t231) / 0.2e1 + t263 * (-mrSges(5,1) * t56 + mrSges(5,2) * t55) + t253 * (-mrSges(10,1) * t83 + mrSges(10,2) * t82) + t235 * t189 / 0.2e1 + t233 * t188 / 0.2e1 + t234 * t192 / 0.2e1 + t230 * t190 / 0.2e1 + t231 * t187 / 0.2e1 + t232 * t191 / 0.2e1 + t214 * t109 + t186 * (-mrSges(11,1) * t26 + mrSges(11,2) * t25) + t164 * (-mrSges(11,1) * t43 + mrSges(11,2) * t42) + t130 * (Ifges(11,1) * t42 + Ifges(11,4) * t43) / 0.2e1 - t83 * t483 + t114 * t70 + (t114 * t164 + t186 * t71) * m(11) + t83 * t98 / 0.2e1 - t84 * t99 / 0.2e1 + t82 * t100 / 0.2e1 + t85 * t101 / 0.2e1 - t234 * t477 - t235 * t476 + t84 * t45 / 0.2e1 - t232 * t475 + t43 * t63 / 0.2e1 + t42 * t64 / 0.2e1 + t84 * t658 + t318 * t468 / 0.2e1 + t319 * t466 / 0.2e1 + qJD(2) ^ 2 * t463 / 0.2e1 + t316 * t464 / 0.2e1 + t536 * t594 + t531 * t596 + (mrSges(5,1) * t84 + mrSges(5,2) * t85) * t523 + t380 * (Ifges(5,5) * t85 - Ifges(5,6) * t84) / 0.2e1 + t390 * (Ifges(9,5) * t230 + Ifges(9,6) * t231) / 0.2e1 + t391 * (Ifges(8,5) * t232 + Ifges(8,6) * t233) / 0.2e1 + t392 * (Ifges(4,5) * t234 + Ifges(4,6) * t235) / 0.2e1 + (m(6) * (qJD(5) * t649 + t2 * t400 + t3 * t409) + t90 * t515 + t409 * t9 + t400 * t10 - t91 * t516) * (-pkin(9) * t223 - pkin(11) * t224 + t263) - t445 * t46 / 0.2e1 + t652 * t496 + t85 * t491 + t233 * mrSges(8,3) * t498 + t405 * t82 * t504 - t84 * t657 + (0.2e1 * Ifges(3,5) * t609 - Ifges(3,6) * t403) * qJDD(2) + (0.2e1 * Ifges(7,5) * t611 + Ifges(7,6) * t408) * qJDD(6) + t227 * (mrSges(6,1) * t445 - mrSges(6,2) * t444) + t155 * (-Ifges(6,4) * t444 - Ifges(6,2) * t445 + Ifges(6,6) * t84) / 0.2e1 + t195 * (-Ifges(6,5) * t444 - Ifges(6,6) * t445 + Ifges(6,3) * t84) / 0.2e1; m(11) * (t134 * t93 - t135 * t92 + t171 * t61 + t172 * t60) - g(1) * (t413 * t589 + t525) - g(2) * (t404 * t589 + t526) + t557 * t213 + (-t411 * (-mrSges(4,2) * t389 + mrSges(4,3) * t169) + t402 * (mrSges(4,1) * t389 - mrSges(4,3) * t170) + (mrSges(8,1) * t388 - mrSges(8,3) * t167 + qJD(7) * t244) * t407 - t631 * t414 * t412 * t359 + t456 * qJD(3) + (-mrSges(8,2) * t388 - qJD(7) * t246 + (qJD(2) * t287 + t168) * mrSges(8,3)) * t398 + m(4) * (-t302 * t411 + t303 * t402) + m(8) * (t300 * t407 + t301 * t398)) * pkin(1) + (t264 * t10 + t212 * t90 - t87 * t91 + (-t57 * mrSges(6,3) - t264 * t91) * qJD(5)) * t409 + (-t212 * t91 - t87 * t90 + (-qJD(5) * t90 - t9) * t264 + (-t3 - t552) * mrSges(6,3)) * t400 + t416 + t417 + t418 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + t272 * t39 + t273 * t38 + t265 * t8 + t212 * t174 - t178 * t70 + t171 * t19 + t172 * t20 + t92 * t118 + t93 * t119 + Ifges(3,3) * qJDD(2) - (-g(3) * t94 - t104 * t473) * t403 + (-g(3) * t104 + t473 * t94) * t412 + (-m(11) * t178 - t69) * t164 + (t113 * t265 - t213 * t227 - t460 * t87 + (t423 + t598) * t264 + t649 * t212) * m(6) + (t280 * t609 + t531 / 0.2e1 - t545 - t271 * t109 + t463 * t596 + t644 * qJD(1) + (-t215 - t218) * t359 - t652 * t608) * qJD(1) + (t116 * t272 + t117 * t273 + t212 * t241 + t213 * t243 - t271 * t544) * m(5) - g(3) * t673; ((m(5) * t116 + t39 + (-m(5) * t243 + m(6) * t227 - t557) * qJD(4)) * t401 + (m(5) * t117 + t38 + (m(5) * t241 + m(6) * t649 - t400 * t91 + t409 * t90 + t174) * qJD(4)) * t410) * pkin(5) + t419 * (pkin(11) + t603) - m(5) * (t293 * t544 + (t241 * t292 + t243 * t291) * qJD(2)) + (-pkin(1) * t456 - t292 * t174 - t291 * t557) * qJD(2) + t423 * mrSges(6,3) - m(6) * (-qJD(2) * t227 * t291 + t57 * t77 + t58 * t78) + (-g(3) * t121 - t120 * t473) * t412 + (-g(3) * t120 + t121 * t473) * t403 + (-t109 * t293 - t359 * t218 - t545) * qJD(1) + t416 - t78 * t90 - t77 * t91 + t669 * (-pkin(5) * t410 - pkin(9)); (-t174 + t572) * t243 + (-t195 * t558 + (-t195 * t58 - t3) * t400) * mrSges(6,3) + (-g(3) * t132 - t133 * t473) * t412 + (-g(3) * t133 + t132 * t473) * t403 + (-mrSges(5,2) * t523 - t439 / 0.2e1 - t438 / 0.2e1 - t437 / 0.2e1 + t674) * t479 + t420 + (-t576 / 0.2e1 - t574 / 0.2e1 - t573 / 0.2e1 - mrSges(5,1) * t523 + t672) * t453 + t557 * t241 + t419 * pkin(11) - t80 * t90 - t79 * t91 - m(6) * (t227 * t241 + t57 * t79 + t58 * t80) - t669 * pkin(9); -t227 * (mrSges(6,1) * t156 + mrSges(6,2) * t155) + (Ifges(6,1) * t155 - t575) * t625 + t46 * t624 + (Ifges(6,5) * t155 - Ifges(6,6) * t156) * t623 - t57 * t90 + t58 * t91 + (t228 * t412 - t229 * t403) * g(3) + (t404 * t457 - t413 * t484) * g(2) + (t404 * t484 + t413 * t457) * g(1) + (t155 * t57 + t156 * t58) * mrSges(6,3) + t502 + (-Ifges(6,2) * t156 + t154 + t47) * t626 - t643; Ifges(7,5) * t317 + Ifges(7,6) * t316 + Ifges(7,3) * qJDD(6) + t470 * g(3) + (-t529 / 0.2e1 + t536 / 0.2e1 - t408 * t366 / 0.2e1 + t461 * t594 + (-t443 - t435 / 0.2e1 + t561 * t611) * qJD(1)) * qJD(1) + t473 * t469; (-g(3) * t323 - t326 * t473) * t412 + (-t244 * t407 + (t246 + t586) * t398) * t591 + (-g(3) * t326 + t323 * t473) * t403 - t659 + t654 * t119 + t653 * t118 + t418 + t210 * t19 + t211 * t20 - t182 * t70 - t164 * t69 - g(2) * (-t404 * t503 + t526) - g(1) * (-t413 * t503 + t525) - t215 * t522 + (-g(3) * t341 + t134 * t654 - t135 * t653 - t164 * t182 + t210 * t61 + t211 * t60) * m(11); (-g(3) * t325 + t322 * t473) * t403 + t417 + (-g(3) * t322 - t325 * t473) * t412; (t173 * t405 - t175 * t396 + (t396 * t454 - t405 * t455) * mrSges(10,3)) * t390 * pkin(2) + t661; -t135 * t119 - g(1) * t525 - g(2) * t526 - t659 + (-t118 + t551) * t134 + (-t164 * mrSges(11,2) + t664) * t480 - (t164 * mrSges(11,1) + t663) * t130 + t665; 0; 0; 0;];
tau = t1;
