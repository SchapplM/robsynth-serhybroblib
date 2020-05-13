% Calculate vector of inverse dynamics joint torques for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% qJDD [10x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [10x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:43
% EndTime: 2020-05-07 04:35:57
% DurationCPUTime: 25.44s
% Computational Cost: add. (9711->775), mult. (22302->1096), div. (0->0), fcn. (17808->28), ass. (0->369)
t316 = qJ(2) + qJ(3);
t309 = qJ(4) + t316;
t286 = sin(t309);
t287 = cos(t309);
t328 = cos(qJ(5));
t501 = mrSges(6,1) * t328;
t599 = -mrSges(5,2) * t287 + (-pkin(8) * m(6) - mrSges(5,1) - t501) * t286;
t594 = m(8) + m(4);
t324 = sin(qJ(2));
t331 = cos(qJ(2));
t260 = -mrSges(3,1) * t331 + mrSges(3,2) * t324;
t315 = qJ(2) + qJ(7);
t304 = sin(t315);
t306 = cos(t315);
t308 = pkin(15) - t315;
t288 = -qJ(8) + t308;
t274 = sin(t288);
t275 = cos(t288);
t386 = t275 * mrSges(9,1) + t274 * mrSges(9,2);
t572 = -t306 * mrSges(8,1) + mrSges(8,2) * t304 + t386;
t598 = t260 + t572;
t330 = cos(qJ(3));
t323 = sin(qJ(3));
t526 = pkin(1) * t323;
t425 = qJD(3) * t526;
t465 = pkin(1) * qJDD(2);
t231 = qJD(2) * t425 - t330 * t465;
t312 = qJDD(2) + qJDD(3);
t194 = pkin(4) * t312 + t231;
t440 = qJD(3) * t330;
t232 = (-qJD(2) * t440 - qJDD(2) * t323) * pkin(1);
t322 = sin(qJ(4));
t329 = cos(qJ(4));
t314 = qJD(2) + qJD(3);
t505 = pkin(1) * qJD(2);
t426 = t330 * t505;
t368 = -pkin(4) * t314 + t426;
t428 = t323 * t505;
t370 = t322 * t368 + t329 * t428;
t100 = qJD(4) * t370 + t329 * t194 - t322 * t232;
t303 = qJD(4) + t314;
t321 = sin(qJ(5));
t238 = t323 * t324 - t330 * t331;
t220 = t238 * qJD(1);
t371 = t323 * t331 + t324 * t330;
t222 = t371 * qJD(1);
t374 = t220 * t322 - t329 * t222;
t125 = t303 * t328 - t321 * t374;
t126 = t303 * t321 + t328 * t374;
t399 = t329 * t220 + t222 * t322;
t152 = Ifges(5,4) * t399;
t153 = qJD(5) - t399;
t310 = t331 * pkin(1);
t285 = t310 + pkin(12);
t261 = t285 * qJD(1);
t180 = -pkin(4) * t220 - t261;
t434 = qJD(1) * qJD(2);
t245 = qJDD(1) * t331 - t324 * t434;
t246 = qJDD(1) * t324 + t331 * t434;
t353 = t238 * qJD(3);
t134 = qJD(1) * t353 - t245 * t323 - t246 * t330;
t354 = t371 * qJD(3);
t135 = qJD(1) * t354 - t245 * t330 + t246 * t323;
t49 = qJD(4) * t399 + t134 * t329 + t135 * t322;
t50 = -qJD(4) * t374 - t134 * t322 + t135 * t329;
t435 = qJDD(1) * pkin(12);
t203 = -pkin(1) * t245 - t435;
t97 = -pkin(4) * t135 + t203;
t11 = -pkin(8) * t50 - pkin(10) * t49 + t97;
t177 = pkin(10) * t303 - t370;
t71 = -pkin(8) * t399 - pkin(10) * t374 + t180;
t51 = -t177 * t321 + t328 * t71;
t301 = qJDD(4) + t312;
t269 = t322 * t428;
t181 = -t329 * t368 + t269;
t99 = t181 * qJD(4) + t322 * t194 + t329 * t232;
t93 = pkin(10) * t301 + t99;
t2 = qJD(5) * t51 + t11 * t321 + t328 * t93;
t52 = t177 * t328 + t321 * t71;
t3 = -qJD(5) * t52 + t11 * t328 - t321 * t93;
t176 = -t303 * pkin(8) - t181;
t389 = mrSges(6,1) * t321 + mrSges(6,2) * t328;
t359 = t176 * t389;
t380 = Ifges(6,5) * t328 - Ifges(6,6) * t321;
t488 = Ifges(6,4) * t328;
t383 = -Ifges(6,2) * t321 + t488;
t489 = Ifges(6,4) * t321;
t385 = Ifges(6,1) * t328 - t489;
t437 = qJD(5) * t321;
t410 = -t437 / 0.2e1;
t124 = Ifges(6,4) * t125;
t45 = t126 * Ifges(6,1) + t153 * Ifges(6,5) + t124;
t468 = t328 * t45;
t415 = t468 / 0.2e1;
t43 = t126 * Ifges(6,5) + t125 * Ifges(6,6) + t153 * Ifges(6,3);
t436 = qJD(5) * t328;
t482 = t126 * Ifges(6,4);
t44 = t125 * Ifges(6,2) + t153 * Ifges(6,6) + t482;
t479 = t370 * mrSges(5,3);
t480 = t181 * mrSges(5,3);
t481 = t374 * Ifges(5,4);
t493 = mrSges(6,3) * t328;
t494 = mrSges(6,3) * t321;
t531 = t321 / 0.2e1;
t533 = -t303 / 0.2e1;
t539 = t374 / 0.2e1;
t540 = -t374 / 0.2e1;
t541 = -t399 / 0.2e1;
t542 = -t153 / 0.2e1;
t544 = -t126 / 0.2e1;
t545 = -t125 / 0.2e1;
t48 = qJDD(5) - t50;
t549 = t48 / 0.2e1;
t18 = -qJD(5) * t126 + t301 * t328 - t321 * t49;
t551 = t18 / 0.2e1;
t17 = qJD(5) * t125 + t301 * t321 + t328 * t49;
t552 = t17 / 0.2e1;
t581 = t52 * mrSges(6,2);
t582 = t51 * mrSges(6,1);
t497 = mrSges(6,2) * t321;
t591 = t497 - t501;
t6 = t17 * Ifges(6,4) + t18 * Ifges(6,2) + t48 * Ifges(6,6);
t7 = t17 * Ifges(6,1) + t18 * Ifges(6,4) + t48 * Ifges(6,5);
t83 = Ifges(5,2) * t399 + t303 * Ifges(5,6) + t481;
t84 = Ifges(5,1) * t374 + t303 * Ifges(5,5) + t152;
t94 = -t301 * pkin(8) - t100;
t595 = -t99 * mrSges(5,2) - t3 * t494 + t2 * t493 + (Ifges(6,1) * t321 + t488) * t552 + (Ifges(6,2) * t328 + t489) * t551 + Ifges(5,3) * t301 + (Ifges(6,5) * t321 + Ifges(6,6) * t328) * t549 + t44 * t410 + t7 * t531 + Ifges(5,6) * t50 + Ifges(5,5) * t49 + t328 * t6 / 0.2e1 + t94 * t591 + t100 * mrSges(5,1) + (-t436 * t51 - t437 * t52) * mrSges(6,3) + (t359 + t415) * qJD(5) + (t125 * t383 + t126 * t385 + t153 * t380) * qJD(5) / 0.2e1 + t83 * t539 + (t152 + t84) * t541 + (-t481 + t43) * t540 + (-t180 * mrSges(5,2) + Ifges(5,1) * t540 + Ifges(5,5) * t533 + t380 * t542 + t383 * t545 + t385 * t544 + t493 * t51 + t494 * t52 - t359 + t44 * t531 - t468 / 0.2e1 + t480) * t399 + (-t180 * mrSges(5,1) + Ifges(6,5) * t544 - Ifges(5,2) * t541 - Ifges(5,6) * t533 + Ifges(6,6) * t545 + Ifges(6,3) * t542 - t479 + t581 - t582) * t374;
t593 = -mrSges(5,1) * t303 - mrSges(6,1) * t125 + mrSges(6,2) * t126 + mrSges(5,3) * t374;
t592 = t287 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t286;
t305 = sin(t316);
t307 = cos(t316);
t590 = mrSges(4,1) * t305 + mrSges(4,2) * t307;
t319 = sin(qJ(7));
t326 = cos(qJ(7));
t237 = -t319 * t324 + t326 * t331;
t219 = t237 * qJD(1);
t239 = t319 * t331 + t324 * t326;
t221 = t239 * qJD(1);
t317 = sin(pkin(15));
t318 = cos(pkin(15));
t131 = -t261 + (-t219 * t318 - t221 * t317) * pkin(3);
t313 = qJD(2) + qJD(7);
t302 = qJD(8) + t313;
t534 = -t302 / 0.2e1;
t528 = sin(qJ(8));
t529 = cos(qJ(8));
t234 = t317 * t528 + t318 * t529;
t343 = -t317 * t529 + t318 * t528;
t109 = -t219 * t343 - t221 * t234;
t546 = -t109 / 0.2e1;
t589 = -t131 * mrSges(9,2) + Ifges(9,1) * t546 + Ifges(9,5) * t534;
t400 = -t234 * t219 + t221 * t343;
t104 = Ifges(9,4) * t400;
t355 = t237 * qJD(7);
t132 = qJD(1) * t355 + t245 * t319 + t246 * t326;
t356 = t239 * qJD(7);
t133 = -qJD(1) * t356 + t245 * t326 - t246 * t319;
t213 = t234 * qJD(8);
t214 = t343 * qJD(8);
t27 = -t132 * t234 - t133 * t343 - t213 * t219 + t214 * t221;
t28 = t132 * t343 - t133 * t234 + t213 * t221 + t214 * t219;
t311 = qJDD(2) + qJDD(7);
t300 = qJDD(8) + t311;
t429 = t319 * t505;
t521 = pkin(3) * t317;
t226 = t313 * t521 + t429;
t427 = t326 * t505;
t520 = pkin(3) * t318;
t227 = t313 * t520 + t427;
t115 = t226 * t234 + t227 * t343;
t483 = t115 * mrSges(9,3);
t504 = pkin(1) * qJD(7);
t420 = qJD(2) * t504;
t229 = -t319 * t420 + t326 * t465;
t191 = t311 * t520 + t229;
t230 = t319 * t465 + t326 * t420;
t192 = t311 * t521 + t230;
t54 = -t191 * t343 - t192 * t234 - t213 * t227 + t214 * t226;
t547 = -t400 / 0.2e1;
t55 = -t191 * t234 + t192 * t343 + t213 * t226 + t214 * t227;
t57 = t109 * Ifges(9,4) + Ifges(9,2) * t400 + t302 * Ifges(9,6);
t58 = t109 * Ifges(9,1) + t302 * Ifges(9,5) + t104;
t587 = t55 * mrSges(9,1) - t54 * mrSges(9,2) + Ifges(9,5) * t27 + Ifges(9,6) * t28 + Ifges(9,3) * t300 + (t58 + t104) * t547 - (t131 * mrSges(9,1) + Ifges(9,4) * t546 + Ifges(9,2) * t547 + Ifges(9,6) * t534 + t483 - t57 / 0.2e1) * t109;
t583 = t245 / 0.2e1;
t320 = sin(qJ(6));
t532 = t320 / 0.2e1;
t580 = t331 * Ifges(3,2);
t63 = -mrSges(9,1) * t400 + mrSges(9,2) * t109;
t579 = m(9) * t131 + t63;
t90 = -mrSges(5,1) * t399 + mrSges(5,2) * t374;
t578 = m(5) * t180 + t90;
t361 = pkin(1) * (t234 * t319 + t326 * t343);
t577 = -qJD(2) * t361 + (t213 * t317 + t214 * t318) * pkin(3);
t360 = pkin(1) * (-t234 * t326 + t319 * t343);
t576 = -qJD(2) * t360 + (-t213 * t318 + t214 * t317) * pkin(3);
t575 = -mrSges(4,1) * t220 - mrSges(8,1) * t219 - mrSges(4,2) * t222 + mrSges(8,2) * t221;
t404 = -mrSges(4,1) * t307 + t305 * mrSges(4,2);
t421 = t286 * t497;
t367 = -t287 * mrSges(6,3) - t421;
t514 = pkin(10) * t287;
t518 = pkin(4) * t305;
t574 = -m(6) * (-t514 + t518) - t367 - m(5) * t518;
t573 = -t591 * t287 + t592;
t332 = cos(qJ(1));
t571 = t599 * t332;
t325 = sin(qJ(1));
t570 = t599 * t325;
t171 = t238 * t322 - t329 * t371;
t174 = qJD(2) * t238 + t353;
t175 = qJD(2) * t371 + t354;
t373 = t329 * t238 + t322 * t371;
t74 = qJD(4) * t373 + t174 * t329 + t175 * t322;
t365 = t171 * t436 + t321 * t74;
t567 = g(1) * t332 + g(2) * t325;
t566 = Ifges(3,6) * qJDD(2);
t8 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t565 = m(6) * t94 + t8;
t564 = -t590 * t332 + t571;
t563 = -t590 * t325 + t570;
t136 = -mrSges(5,2) * t303 + mrSges(5,3) * t399;
t77 = -mrSges(6,2) * t153 + mrSges(6,3) * t125;
t78 = mrSges(6,1) * t153 - mrSges(6,3) * t126;
t562 = -t321 * t78 + t328 * t77 + t136;
t327 = cos(qJ(6));
t259 = -mrSges(7,1) * t327 + mrSges(7,2) * t320;
t350 = m(7) * pkin(6) + t259;
t525 = pkin(1) * t324;
t249 = t518 - t525;
t500 = mrSges(9,1) * t274;
t522 = pkin(3) * sin(t308);
t561 = -m(6) * (t249 - t514) - t367 - m(9) * (t522 - t525) + t500 + m(4) * t525 - m(5) * t249;
t560 = m(6) * t176 + t593;
t559 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t378 = t52 * t321 + t51 * t328;
t558 = m(6) * t378 + t321 * t77 + t328 * t78;
t91 = pkin(8) * t374 - pkin(10) * t399;
t392 = pkin(8) * t287 + pkin(10) * t286;
t517 = pkin(4) * t307;
t557 = -m(6) * (-t392 - t517) - t404 + t573;
t556 = mrSges(2,2) - mrSges(3,3) - mrSges(7,3) - mrSges(9,3) - mrSges(8,3) - mrSges(5,3) - mrSges(4,3);
t10 = -mrSges(6,2) * t48 + mrSges(6,3) * t18;
t9 = mrSges(6,1) * t48 - mrSges(6,3) * t17;
t555 = m(6) * (-qJD(5) * t378 + t2 * t328 - t3 * t321) - t78 * t436 - t77 * t437 + t328 * t10 - t321 * t9;
t409 = t310 - t517;
t242 = pkin(12) + t409;
t273 = pkin(3) * cos(t308);
t444 = t273 + t310;
t554 = m(6) * (-t242 + t392) - mrSges(2,1) - m(9) * (pkin(12) + t444) - m(5) * t242 - t404 - m(3) * pkin(12) + t350 + t592 - t594 * t285 + t598;
t543 = t126 / 0.2e1;
t536 = t221 / 0.2e1;
t535 = -t222 / 0.2e1;
t530 = t331 / 0.2e1;
t527 = pkin(1) * t319;
t524 = pkin(1) * t326;
t523 = pkin(1) * t330;
t519 = pkin(4) * t222;
t516 = pkin(4) * t322;
t515 = pkin(4) * t329;
t506 = -qJD(6) / 0.2e1;
t495 = mrSges(9,2) * t275;
t492 = mrSges(8,3) * t221;
t491 = Ifges(3,4) * t324;
t490 = Ifges(3,4) * t331;
t487 = Ifges(7,4) * t320;
t486 = Ifges(7,4) * t327;
t485 = t400 * mrSges(9,3);
t114 = t226 * t343 - t227 * t234;
t484 = t114 * mrSges(9,3);
t478 = t221 * Ifges(8,4);
t477 = t222 * Ifges(4,4);
t474 = t320 * Ifges(7,1);
t470 = t327 * Ifges(7,2);
t464 = t171 * t321;
t463 = t171 * t328;
t462 = (-mrSges(8,2) * t313 + mrSges(8,3) * t219) * t326;
t382 = t470 + t487;
t456 = t320 * (Ifges(7,6) * qJD(6) + qJD(1) * t382);
t455 = t321 * t332;
t454 = t322 * t323;
t453 = t323 * t329;
t452 = t325 * t321;
t451 = t325 * t328;
t294 = qJD(1) * t486;
t450 = t327 * (Ifges(7,5) * qJD(6) + qJD(1) * t474 + t294);
t449 = t328 * t332;
t443 = qJD(1) * t324;
t442 = qJD(1) * t331;
t441 = qJD(2) * t324;
t439 = qJD(4) * t322;
t438 = qJD(4) * t329;
t433 = qJD(1) * qJD(6);
t430 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t48;
t297 = pkin(1) * t441;
t158 = -pkin(4) * t175 + t297;
t406 = t434 / 0.2e1;
t398 = mrSges(4,3) * t428;
t397 = mrSges(4,3) * t426;
t396 = mrSges(8,3) * t429;
t395 = mrSges(8,3) * t427;
t193 = -pkin(4) * t238 - t285;
t291 = pkin(4) - t523;
t210 = -pkin(1) * t453 + t322 * t291;
t390 = mrSges(3,1) * t324 + mrSges(3,2) * t331;
t388 = mrSges(7,1) * t320 + mrSges(7,2) * t327;
t387 = -mrSges(8,1) * t304 - mrSges(8,2) * t306;
t384 = t491 + t580;
t381 = Ifges(3,5) * t331 - Ifges(3,6) * t324;
t379 = Ifges(7,5) * t327 - Ifges(7,6) * t320;
t377 = -t321 * t51 + t328 * t52;
t372 = t322 * t330 + t453;
t209 = pkin(1) * t454 + t291 * t329;
t364 = t171 * t437 - t328 * t74;
t363 = pkin(6) * t388;
t362 = pkin(12) * t390;
t358 = t320 * (Ifges(7,1) * t327 - t487);
t357 = t324 * (Ifges(3,1) * t331 - t491);
t150 = (-t219 * t317 + t221 * t318) * pkin(3);
t76 = -t519 + t91;
t341 = -t421 + (-m(6) * pkin(10) - mrSges(6,3)) * t287;
t146 = t219 * Ifges(8,2) + t313 * Ifges(8,6) + t478;
t206 = Ifges(8,4) * t219;
t148 = t221 * Ifges(8,1) + t313 * Ifges(8,5) + t206;
t334 = -t230 * mrSges(8,2) + t261 * (mrSges(8,1) * t221 + mrSges(8,2) * t219) + t146 * t536 - t221 * (Ifges(8,1) * t219 - t478) / 0.2e1 + Ifges(8,6) * t133 + Ifges(8,5) * t132 - t313 * (Ifges(8,5) * t219 - Ifges(8,6) * t221) / 0.2e1 + t219 * t395 + t229 * mrSges(8,1) + Ifges(8,3) * t311 - (-Ifges(8,2) * t221 + t148 + t206) * t219 / 0.2e1 + (t484 + t589) * t400 + t587;
t147 = t220 * Ifges(4,2) + t314 * Ifges(4,6) - t477;
t207 = Ifges(4,4) * t220;
t149 = -t222 * Ifges(4,1) + t314 * Ifges(4,5) + t207;
t333 = -t220 * t397 + t222 * (Ifges(4,1) * t220 + t477) / 0.2e1 + t261 * (-mrSges(4,1) * t222 + mrSges(4,2) * t220) + Ifges(4,3) * t312 - t314 * (Ifges(4,5) * t220 + Ifges(4,6) * t222) / 0.2e1 + t231 * mrSges(4,1) - t232 * mrSges(4,2) + Ifges(4,5) * t134 + Ifges(4,6) * t135 + t222 * t398 + t147 * t535 - (Ifges(4,2) * t222 + t149 + t207) * t220 / 0.2e1 + t595;
t295 = Ifges(3,4) * t442;
t290 = -pkin(8) - t515;
t263 = t520 + t524;
t262 = t521 + t527;
t248 = t332 * t495;
t247 = t325 * t495;
t244 = qJDD(1) * t320 + t327 * t433;
t243 = qJDD(1) * t327 - t320 * t433;
t218 = Ifges(3,1) * t443 + Ifges(3,5) * qJD(2) + t295;
t216 = Ifges(3,6) * qJD(2) + qJD(1) * t384;
t205 = -t329 * t426 + t269;
t204 = t372 * t505;
t200 = t287 * t449 - t452;
t199 = t287 * t455 + t451;
t198 = t287 * t451 + t455;
t197 = t287 * t452 - t449;
t187 = mrSges(4,1) * t314 + mrSges(4,3) * t222;
t186 = mrSges(8,1) * t313 - t492;
t185 = -mrSges(4,2) * t314 + mrSges(4,3) * t220;
t173 = -qJD(2) * t239 - t356;
t172 = qJD(2) * t237 + t355;
t164 = (-t234 * t317 - t318 * t343) * pkin(3);
t163 = (-t234 * t318 + t317 * t343) * pkin(3);
t139 = -t234 * t262 - t263 * t343;
t138 = -t234 * t263 + t262 * t343;
t96 = mrSges(9,1) * t302 - mrSges(9,3) * t109;
t95 = -mrSges(9,2) * t302 + t485;
t80 = qJD(7) * t361 + t213 * t262 + t214 * t263;
t79 = qJD(7) * t360 - t213 * t263 + t214 * t262;
t75 = qJD(4) * t171 + t174 * t322 - t329 * t175;
t70 = t181 * t328 + t321 * t91;
t69 = -t181 * t321 + t328 * t91;
t65 = t205 * t328 + t321 * t76;
t64 = -t205 * t321 + t328 * t76;
t61 = (-t132 * t317 - t133 * t318) * pkin(3) + t203;
t42 = t172 * t343 - t173 * t234 + t213 * t239 + t214 * t237;
t41 = -t172 * t234 - t173 * t343 - t213 * t237 + t214 * t239;
t37 = -mrSges(5,2) * t301 + mrSges(5,3) * t50;
t36 = mrSges(5,1) * t301 - mrSges(5,3) * t49;
t20 = -mrSges(9,2) * t300 + mrSges(9,3) * t28;
t19 = mrSges(9,1) * t300 - mrSges(9,3) * t27;
t1 = [t594 * (-t203 * t285 - t261 * t297) + (Ifges(7,1) * t244 + Ifges(7,4) * t243) * t532 - (-mrSges(4,1) * t175 - mrSges(8,1) * t173 + mrSges(4,2) * t174 + mrSges(8,2) * t172) * t261 + (t200 * mrSges(6,1) - t199 * mrSges(6,2) + t325 * t556 + t332 * t554) * g(2) + (-t198 * mrSges(6,1) + t197 * mrSges(6,2) - t325 * t554 + t332 * t556) * g(1) + (mrSges(4,1) * t135 + mrSges(8,1) * t133 - mrSges(4,2) * t134 - mrSges(8,2) * t132) * t285 + t558 * (pkin(8) * t75 - pkin(10) * t74 + t158) + (t77 * t436 + m(6) * (qJD(5) * t377 + t2 * t321 + t3 * t328) + t321 * t10 + t328 * t9 - t78 * t437) * (-pkin(8) * t373 - pkin(10) * t171 + t193) - (-t99 * mrSges(5,3) - Ifges(5,4) * t49 - Ifges(5,2) * t50 - Ifges(5,6) * t301 + t97 * mrSges(5,1) + t430 / 0.2e1 + Ifges(6,3) * t549 + Ifges(6,6) * t551 + Ifges(6,5) * t552 + t559) * t373 + (-mrSges(9,1) * t61 + mrSges(9,3) * t54 + Ifges(9,4) * t27 + Ifges(9,2) * t28 + Ifges(9,6) * t300) * (-t234 * t237 + t239 * t343) + (mrSges(9,2) * t61 - mrSges(9,3) * t55 + Ifges(9,1) * t27 + Ifges(9,4) * t28 + Ifges(9,5) * t300) * (-t234 * t239 - t237 * t343) - (mrSges(4,2) * t203 - mrSges(4,3) * t231 + Ifges(4,1) * t134 + Ifges(4,4) * t135 + Ifges(4,5) * t312) * t371 + t327 * (Ifges(7,4) * t244 + Ifges(7,2) * t243) / 0.2e1 + t75 * t582 + t384 * t583 + (t218 * t530 + t381 * qJD(2) / 0.2e1) * qJD(2) - t175 * t398 - t74 * t480 + t246 * t490 / 0.2e1 + t75 * t479 - t42 * t483 + t363 * t433 + t400 * (Ifges(9,4) * t41 + Ifges(9,2) * t42) / 0.2e1 + t399 * (Ifges(5,4) * t74 - Ifges(5,2) * t75) / 0.2e1 - t41 * t484 - t216 * t441 / 0.2e1 + (-t2 * t464 - t3 * t463 + t364 * t51 - t365 * t52) * mrSges(6,3) + t244 * (t474 + t486) / 0.2e1 + (m(5) * t97 - mrSges(5,1) * t50 + mrSges(5,2) * t49) * t193 + (-mrSges(4,1) * t203 + mrSges(4,3) * t232 + Ifges(4,4) * t134 + Ifges(4,2) * t135 + Ifges(4,6) * t312) * t238 + (mrSges(8,2) * t203 - mrSges(8,3) * t229 + Ifges(8,1) * t132 + Ifges(8,4) * t133 + Ifges(8,5) * t311) * t239 - t172 * t395 + (m(9) * t61 - mrSges(9,1) * t28 + mrSges(9,2) * t27) * ((-t237 * t318 - t239 * t317) * pkin(3) - t285) + (t358 + t327 * (-Ifges(7,2) * t320 + t486)) * t433 / 0.2e1 + t7 * t463 / 0.2e1 - t6 * t464 / 0.2e1 - t362 * t434 - t260 * t435 + (t97 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,1) * t49 + Ifges(5,4) * t50 + Ifges(5,5) * t301 + t380 * t549 + t383 * t551 + t385 * t552 + t389 * t94 + t410 * t45) * t171 + t313 * (Ifges(8,5) * t172 + Ifges(8,6) * t173) / 0.2e1 + t314 * (Ifges(4,5) * t174 + Ifges(4,6) * t175) / 0.2e1 + t302 * (Ifges(9,5) * t41 + Ifges(9,6) * t42) / 0.2e1 + t303 * (Ifges(5,5) * t74 - Ifges(5,6) * t75) / 0.2e1 + t243 * t382 / 0.2e1 + pkin(6) * (-mrSges(7,1) * t243 + mrSges(7,2) * t244) - pkin(12) * (-mrSges(3,1) * t245 + mrSges(3,2) * t246) + t153 * (-Ifges(6,5) * t364 - Ifges(6,6) * t365 + Ifges(6,3) * t75) / 0.2e1 + t125 * (-Ifges(6,4) * t364 - Ifges(6,2) * t365 + Ifges(6,6) * t75) / 0.2e1 + t219 * (Ifges(8,4) * t172 + Ifges(8,2) * t173) / 0.2e1 + t220 * (Ifges(4,4) * t174 + Ifges(4,2) * t175) / 0.2e1 + t176 * (mrSges(6,1) * t365 - mrSges(6,2) * t364) + t180 * (mrSges(5,1) * t75 + mrSges(5,2) * t74) + t172 * t148 / 0.2e1 + t173 * t146 / 0.2e1 + t174 * t149 / 0.2e1 + t175 * t147 / 0.2e1 + t131 * (-mrSges(9,1) * t42 + mrSges(9,2) * t41) + t109 * (Ifges(9,1) * t41 + Ifges(9,4) * t42) / 0.2e1 - t75 * t83 / 0.2e1 + t74 * t84 / 0.2e1 + t75 * t43 / 0.2e1 + t42 * t57 / 0.2e1 + t41 * t58 / 0.2e1 + t173 * t396 + t174 * t397 + (m(3) * pkin(12) ^ 2 + pkin(6) * t350 + Ifges(2,3)) * qJDD(1) + (t566 / 0.2e1 + t490 * t406) * t331 + (Ifges(3,4) * t246 + Ifges(3,2) * t245 + t566) * t530 - t365 * t44 / 0.2e1 + t575 * t297 + t578 * t158 + t579 * (t297 + (-t172 * t317 - t173 * t318) * pkin(3)) + t456 * t506 + (Ifges(4,1) * t174 + Ifges(4,4) * t175) * t535 + (Ifges(8,1) * t172 + Ifges(8,4) * t173) * t536 + (Ifges(5,1) * t74 - Ifges(5,4) * t75) * t539 + (-Ifges(6,1) * t364 - Ifges(6,4) * t365 + Ifges(6,5) * t75) * t543 + t357 * t406 + t74 * t415 + (t450 / 0.2e1 + t379 * qJD(6) / 0.2e1) * qJD(6) - t75 * t581 + (0.2e1 * Ifges(7,5) * t532 + Ifges(7,6) * t327) * qJDD(6) + (Ifges(3,1) * t246 + Ifges(3,4) * t583 + Ifges(3,5) * qJDD(2) - t406 * t580) * t324 + (-mrSges(8,1) * t203 + mrSges(8,3) * t230 + Ifges(8,4) * t132 + Ifges(8,2) * t133 + Ifges(8,6) * t311) * t237; (-m(5) * t409 - m(9) * t444 + (-m(6) - t594) * t310 + t557 + t598) * g(3) + (-t186 * t319 + t462) * t504 + m(9) * (t114 * t80 - t115 * t79 + t138 * t55 + t139 * t54) + (-t185 * t440 + m(8) * (t229 * t326 + t230 * t319) + m(4) * (-t231 * t330 - t232 * t323)) * pkin(1) + t555 * (pkin(10) + t210) + m(5) * (t100 * t209 + t210 * t99) + t216 * t443 / 0.2e1 + (t261 * t594 - t558 - t575 - t578 - t579) * pkin(1) * t443 - (mrSges(4,1) * t312 - mrSges(4,3) * t134) * t523 - (-Ifges(3,2) * t443 + t218 + t295) * t442 / 0.2e1 + (t362 - t357 / 0.2e1) * qJD(1) ^ 2 - (-mrSges(4,2) * t312 + mrSges(4,3) * t135) * t526 - t381 * t434 / 0.2e1 + (m(8) * t525 - t387 + t390) * t567 + t334 + t187 * t425 - t558 * t76 + t578 * t519 - t579 * t150 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t245 + Ifges(3,5) * t246 + t209 * t36 + t210 * t37 + t138 * t19 + t139 * t20 + t79 * t95 + t80 * t96 + t221 * t396 + (m(5) * t181 - t560) * (-t291 * t439 + (qJD(3) * t372 + t323 * t438) * pkin(1)) + (-m(5) * t370 + m(6) * t377 + t562) * (t291 * t438 + (t323 * t439 + (-t329 * t330 + t454) * qJD(3)) * pkin(1)) + (t325 * t561 - t247 + t563) * g(2) + (t332 * t561 - t248 + t564) * g(1) + t565 * (-pkin(8) - t209) + (mrSges(8,1) * t311 - mrSges(8,3) * t132) * t524 + (-mrSges(8,2) * t311 + mrSges(8,3) * t133) * t527 + t333; -t205 * t136 + t185 * t426 - t187 * t428 + t290 * t8 + t36 * t515 + t37 * t516 + t90 * t519 - t64 * t78 - t65 * t77 + t333 + (t290 * t94 + (t176 * t322 + t329 * t377) * qJD(4) * pkin(4) + t176 * t204 - t51 * t64 - t52 * t65) * m(6) + (t180 * t519 - t181 * t204 + t370 * t205 + (t100 * t329 + t322 * t99 + (-t181 * t322 - t329 * t370) * qJD(4)) * pkin(4)) * m(5) + t562 * pkin(4) * t438 + (m(5) * t517 + t557) * g(3) + (t325 * t574 + t563) * g(2) + (t332 * t574 + t564) * g(1) + t555 * (pkin(10) + t516) + t593 * (pkin(4) * t439 + t204); -m(6) * (t51 * t69 + t52 * t70) - t181 * t136 - t70 * t77 - t69 * t78 + t560 * t370 + (m(6) * t392 + t573) * g(3) + (-t325 * t341 + t570) * g(2) + (-t332 * t341 + t571) * g(1) - t565 * pkin(8) + t555 * pkin(10) + t595; -t176 * (mrSges(6,1) * t126 + mrSges(6,2) * t125) + (Ifges(6,1) * t125 - t482) * t544 + t44 * t543 + (Ifges(6,5) * t125 - Ifges(6,6) * t126) * t542 - t51 * t77 + t52 * t78 - g(1) * (mrSges(6,1) * t199 + mrSges(6,2) * t200) - g(2) * (mrSges(6,1) * t197 + mrSges(6,2) * t198) - g(3) * t389 * t286 + (t125 * t51 + t126 * t52) * mrSges(6,3) + t430 + (-Ifges(6,2) * t126 + t124 + t45) * t545 + t559; Ifges(7,5) * t244 + Ifges(7,6) * t243 + Ifges(7,3) * qJDD(6) + g(3) * t259 + (-t450 / 0.2e1 + t456 / 0.2e1 - t327 * t294 / 0.2e1 + t379 * t506 + (-t363 - t358 / 0.2e1 + t470 * t532) * qJD(1)) * qJD(1) + t567 * t388; -g(1) * t248 - g(2) * t247 + t572 * g(3) + t334 + (-t462 + (t186 + t492) * t319) * t505 + t577 * t96 + t576 * t95 + t163 * t19 + t164 * t20 - t150 * t63 + t567 * (-m(9) * t522 - t387 + t500) + (-t273 * g(3) + t114 * t577 - t115 * t576 - t131 * t150 + t163 * t55 + t164 * t54) * m(9); -t115 * t96 + g(3) * t386 - g(1) * (-t332 * t500 + t248) - g(2) * (-t325 * t500 + t247) + (-t95 + t485) * t114 + t589 * t400 + t587; 0; 0;];
tau = t1;
