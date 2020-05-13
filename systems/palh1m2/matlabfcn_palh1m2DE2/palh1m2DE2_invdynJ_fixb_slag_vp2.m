% Calculate vector of inverse dynamics joint torques for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m2DE2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:30
% EndTime: 2020-05-02 21:00:17
% DurationCPUTime: 38.30s
% Computational Cost: add. (11864->710), mult. (21875->1000), div. (0->0), fcn. (25870->80), ass. (0->364)
t310 = cos(qJ(3));
t419 = qJD(2) * t310;
t402 = pkin(1) * t419;
t292 = pkin(22) + pkin(21);
t270 = sin(t292);
t271 = cos(t292);
t301 = cos(pkin(20));
t307 = sin(pkin(18));
t313 = cos(pkin(18));
t445 = sin(pkin(20));
t216 = t301 * t307 - t313 * t445;
t220 = t313 * t301 + t307 * t445;
t304 = sin(qJ(3));
t138 = t216 * t310 - t220 * t304;
t311 = cos(qJ(2));
t121 = t138 * t311;
t305 = sin(qJ(2));
t350 = t304 * t216 + t220 * t310;
t385 = -t305 * t350 + t121;
t590 = t305 * t138 + t311 * t350;
t587 = t270 * t385 + t271 * t590;
t293 = qJD(2) + qJD(3);
t471 = pkin(1) * qJD(2);
t403 = t304 * t471;
t349 = pkin(5) * t293 + t403;
t359 = -t590 * t270 + t271 * t385;
t595 = t349 * t359;
t18 = t587 * t402 - t595;
t131 = t138 * qJD(3);
t126 = t350 * qJD(3);
t591 = qJD(2) * t590 + t126 * t311 + t305 * t131;
t593 = -qJD(2) * t385 + t305 * t126;
t24 = t591 * t271 + t270 * (t131 * t311 - t593);
t556 = Ifges(3,4) + Ifges(10,4);
t592 = t587 * t349;
t586 = Ifges(3,1) + Ifges(10,1);
t585 = Ifges(3,2) + Ifges(10,2);
t25 = (-qJD(3) * t121 + t593) * t271 + t591 * t270;
t573 = pkin(1) * t304;
t263 = pkin(5) + t573;
t347 = pkin(1) * t359;
t491 = pkin(1) * t310;
t408 = t587 * t491;
t492 = -t24 * t491 - t25 * t263 + (-t304 * t347 + t408) * qJD(3) + t18;
t589 = t556 * t311;
t588 = t556 * t305;
t555 = Ifges(10,5) + Ifges(3,5);
t554 = Ifges(10,6) + Ifges(3,6);
t302 = cos(pkin(19));
t446 = sin(pkin(19));
t215 = t304 * t302 + t310 * t446;
t217 = t302 * t310 - t304 * t446;
t140 = -t215 * t305 + t217 * t311;
t333 = t140 * qJD(1);
t116 = Ifges(9,4) * t333;
t584 = Ifges(9,2) * t333;
t303 = sin(qJ(4));
t309 = cos(qJ(4));
t114 = -t216 * t271 + t220 * t270;
t103 = t114 * qJD(1);
t351 = t216 * t270 + t220 * t271;
t104 = t351 * qJD(1);
t226 = -t304 * t305 + t310 * t311;
t207 = t226 * qJD(1);
t259 = t305 * pkin(1) - pkin(15);
t243 = t259 * qJD(1);
t165 = -pkin(5) * t207 + t243;
t47 = pkin(9) * t104 - pkin(11) * t103 + t165;
t11 = -t18 * t303 + t309 * t47;
t12 = t18 * t309 + t303 * t47;
t534 = t11 * t303 - t12 * t309;
t583 = t534 * mrSges(6,3);
t319 = m(5) + m(6);
t578 = m(8) + m(4);
t400 = t319 + t578;
t582 = m(11) + t400;
t581 = -t585 * t305 + t589;
t580 = t586 * t311 - t588;
t308 = sin(pkin(17));
t314 = cos(pkin(17));
t224 = t307 * t314 - t313 * t308;
t227 = t307 * t308 + t313 * t314;
t148 = t224 * t311 - t227 * t305;
t577 = -t148 / 0.2e1;
t505 = t148 / 0.2e1;
t493 = t311 / 0.2e1;
t137 = t215 * t311 + t217 * t305;
t532 = t137 * qJD(1);
t575 = t532 / 0.2e1;
t17 = t359 * t402 + t592;
t6 = t24 * t263 + (-t25 * t310 + (-t304 * t587 - t310 * t359) * qJD(3)) * pkin(1);
t574 = -t17 + t6;
t528 = pkin(5) * t319 + mrSges(4,1);
t572 = g(1) * t528;
t472 = qJD(2) / 0.2e1;
t571 = Ifges(7,4) * t148;
t570 = Ifges(9,4) * t532;
t569 = mrSges(11,1) + mrSges(4,1);
t448 = mrSges(11,2) + mrSges(4,2);
t330 = qJD(2) * t347;
t16 = t310 * t330 + t592;
t20 = t263 * t587 + t310 * t347;
t412 = qJD(2) * qJD(3);
t213 = (-qJDD(2) * t310 + t304 * t412) * pkin(1);
t410 = qJDD(2) * t304;
t230 = t310 * t412 + t410;
t214 = pkin(1) * t230;
t291 = qJDD(2) + qJDD(3);
t3 = -t24 * t402 - t25 * t349 + t587 * (pkin(5) * t291 + t214) - t359 * t213;
t567 = t492 * t16 + t20 * t3;
t294 = qJ(3) + pkin(19);
t276 = qJ(2) + t294;
t296 = qJ(2) + qJ(3);
t278 = sin(t296);
t281 = cos(t296);
t524 = pkin(2) * m(10);
t289 = mrSges(9,1) + t524;
t565 = t528 * t281 - mrSges(4,2) * t278 - mrSges(9,2) * sin(t276) + t289 * cos(t276);
t225 = t304 * t311 + t305 * t310;
t206 = t225 * qJD(1);
t465 = Ifges(4,4) * t206;
t109 = t207 * Ifges(4,2) + Ifges(4,6) * t293 + t465;
t193 = Ifges(4,4) * t207;
t110 = t206 * Ifges(4,1) + Ifges(4,5) * t293 + t193;
t386 = pkin(18) - t292;
t245 = t386 - t296;
t295 = pkin(18) - pkin(22);
t277 = -qJ(2) + t295;
t390 = pkin(21) - qJ(2) - atan2(cos(t277), -sin(t277));
t117 = -atan2(-sin(t245), cos(t245)) + t390;
t387 = pkin(2) * t293;
t145 = t215 * t387;
t147 = t217 * t387;
t442 = Ifges(11,4) * t278;
t355 = Ifges(11,2) * t281 + t442;
t174 = Ifges(11,6) * t293 + qJD(1) * t355;
t423 = qJD(1) * t281;
t244 = Ifges(11,4) * t423;
t424 = qJD(1) * t278;
t175 = Ifges(11,1) * t424 + Ifges(11,5) * t293 + t244;
t422 = qJD(1) * t293;
t191 = qJDD(1) * t278 + t281 * t422;
t192 = qJDD(1) * t281 - t278 * t422;
t198 = t215 * qJD(3);
t199 = t217 * qJD(3);
t413 = qJD(1) * qJD(2);
t231 = qJDD(1) * t311 - t305 * t413;
t397 = t311 * t413;
t232 = -qJDD(1) * t305 - t397;
t468 = mrSges(10,3) * t305;
t235 = -qJD(2) * mrSges(10,2) - qJD(1) * t468;
t421 = qJD(1) * t311;
t236 = qJD(2) * mrSges(10,1) - mrSges(10,3) * t421;
t257 = qJ(1) + t276;
t258 = -qJ(1) + t276;
t264 = t289 * g(1);
t285 = qJ(1) + t296;
t286 = qJ(1) - t296;
t315 = mrSges(9,2) * g(1);
t316 = mrSges(4,2) * g(2);
t317 = mrSges(4,2) * g(1);
t321 = qJD(1) ^ 2;
t336 = t278 * (Ifges(11,1) * t281 - t442);
t354 = Ifges(11,5) * t281 - Ifges(11,6) * t278;
t443 = mrSges(11,3) * t281;
t377 = t403 * t443;
t382 = mrSges(4,3) * t402;
t383 = mrSges(4,3) * t403;
t478 = g(2) * t289;
t479 = g(2) * t528;
t490 = pkin(2) * t215;
t497 = -t293 / 0.2e1;
t503 = t206 / 0.2e1;
t511 = mrSges(9,2) * g(2);
t112 = qJ(1) + t117;
t480 = mrSges(11,2) * g(2);
t481 = mrSges(11,2) * g(1);
t482 = mrSges(11,1) * g(2);
t483 = mrSges(11,1) * g(1);
t529 = (-t481 - t482) * sin(t112) / 0.2e1 - (-t480 + t483) * cos(t112) / 0.2e1;
t353 = -t198 * t311 - t199 * t305;
t54 = qJD(1) * t353 + t215 * t232 + t217 * t231;
t352 = -t198 * t305 + t199 * t311;
t55 = qJD(1) * t352 + t215 * t231 - t217 * t232;
t113 = -qJ(1) + t117;
t553 = -(t480 + t483) * cos(t113) / 0.2e1 + (-t481 + t482) * sin(t113) / 0.2e1;
t62 = Ifges(9,6) * t293 + t570 + t584;
t63 = Ifges(9,1) * t532 + Ifges(9,5) * t293 + t116;
t76 = (-t198 * t293 + t217 * t291) * pkin(2);
t77 = (t199 * t293 + t215 * t291) * pkin(2);
t332 = t225 * qJD(3);
t92 = -qJD(1) * t332 + t231 * t310 + t232 * t304;
t331 = t226 * qJD(3);
t93 = qJD(1) * t331 + t231 * t304 - t232 * t310;
t563 = (qJDD(2) * mrSges(10,1) - mrSges(10,3) * t231) * t490 - (-Ifges(9,2) * t532 + t116 + t63) * t333 / 0.2e1 - t206 * (Ifges(4,1) * t207 - t465) / 0.2e1 + (t315 - t478) * cos(t257) / 0.2e1 + (t315 + t478) * cos(t258) / 0.2e1 + (t317 - t479) * cos(t285) / 0.2e1 + (t317 + t479) * cos(t286) / 0.2e1 + (Ifges(4,3) + Ifges(9,3) + Ifges(11,3)) * t291 + t109 * t503 + qJD(1) * t377 + (-t198 * t235 + t199 * t236 + t217 * (-qJDD(2) * mrSges(10,2) + mrSges(10,3) * t232)) * pkin(2) + (t145 * t199 - t147 * t198 + t215 * t77 + t217 * t76) * t524 + t553 + (t264 - t511) * sin(t258) / 0.2e1 + (t264 + t511) * sin(t257) / 0.2e1 + (t316 - t572) * sin(t286) / 0.2e1 + (t316 + t572) * sin(t285) / 0.2e1 - (Ifges(9,1) * t333 - t570) * t532 / 0.2e1 + t529 - t321 * t336 / 0.2e1 + (Ifges(4,5) * t207 + Ifges(9,5) * t333 - Ifges(4,6) * t206 - Ifges(9,6) * t532) * t497 + t214 * mrSges(4,1) + Ifges(11,5) * t191 + Ifges(11,6) * t192 - (-Ifges(4,2) * t206 + t110 + t193) * t207 / 0.2e1 - (-Ifges(11,2) * t424 + t175 + t244) * t423 / 0.2e1 - t206 * t382 + t207 * t383 + Ifges(4,6) * t92 + Ifges(4,5) * t93 + Ifges(9,5) * t55 + Ifges(9,6) * t54 + (sin(t117) * g(3) + t214) * mrSges(11,1) + t62 * t575 - t354 * t422 / 0.2e1 + t174 * t424 / 0.2e1;
t562 = m(9) + m(3);
t561 = pkin(14) * m(7);
t560 = t192 / 0.2e1;
t498 = t291 / 0.2e1;
t559 = t293 / 0.2e1;
t384 = t224 * t305 + t227 * t311;
t552 = Ifges(7,4) * t384;
t558 = Ifges(7,2) * t505 + t552 / 0.2e1;
t557 = pkin(2) * t532;
t61 = mrSges(5,1) * t104 + mrSges(5,2) * t103;
t548 = -m(5) * t165 - t61;
t545 = (mrSges(11,3) * t423 + mrSges(4,3) * t207 - t293 * t448) * t304;
t376 = mrSges(6,1) * t303 + mrSges(6,2) * t309;
t544 = -mrSges(11,3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3) - mrSges(10,3) - t376;
t540 = qJD(1) * t581 + qJD(2) * t554;
t539 = qJD(1) * t580 + qJD(2) * t555;
t388 = mrSges(6,1) * t309 - mrSges(6,2) * t303;
t374 = mrSges(7,1) * t307 - mrSges(7,2) * t313;
t375 = mrSges(7,1) * t313 + mrSges(7,2) * t307;
t538 = t561 + (pkin(1) * t400 + mrSges(4,2) * t310 + t304 * t528 + t308 * t374 + t314 * t375 + mrSges(3,1)) * t305 - (-t304 * mrSges(4,2) - t308 * t375 + t310 * t528 + t314 * t374 - mrSges(3,2)) * t311 - mrSges(2,1) + (-m(10) - t562 - t582) * pkin(15);
t537 = -t305 * t555 - t311 * t554;
t469 = mrSges(6,3) * t103;
t97 = qJD(4) + t104;
t59 = -mrSges(6,2) * t97 - t303 * t469;
t60 = mrSges(6,1) * t97 - t309 * t469;
t535 = -t303 * t59 - t309 * t60;
t4 = -t587 * t213 + (t24 * t293 - t291 * t359) * pkin(5) + (-t359 * t410 + (t24 * t304 + (-qJD(3) * t359 - t25) * t310) * qJD(2)) * pkin(1);
t100 = t351 * qJDD(1);
t240 = t259 * qJDD(1);
t246 = pkin(1) * t397;
t200 = t246 + t240;
t78 = -pkin(5) * t92 + t200;
t99 = t114 * qJDD(1);
t42 = pkin(9) * t100 - pkin(11) * t99 + t78;
t1 = qJD(4) * t11 + t303 * t42 + t309 * t4;
t2 = -qJD(4) * t12 - t303 * t4 + t309 * t42;
t533 = t1 * t309 - t2 * t303;
t531 = 0.2e1 * t498;
t356 = -mrSges(11,1) * t281 + mrSges(11,2) * t278;
t299 = sin(pkin(22));
t447 = cos(pkin(22));
t218 = t307 * t299 + t313 * t447;
t219 = t313 * t299 - t307 * t447;
t373 = mrSges(8,1) * t218 + mrSges(8,2) * t219;
t530 = -mrSges(4,1) * t207 + mrSges(4,2) * t206 + (t356 + t373) * qJD(1);
t527 = -pkin(15) * (mrSges(3,1) * t311 - mrSges(3,2) * t305) - (-t585 * t311 - t588) * t305 / 0.2e1 + (-t305 * t586 - t589) * t493;
t95 = (cos(pkin(21)) * t218 - t219 * sin(pkin(21))) * pkin(4) + t259;
t450 = t95 * (mrSges(11,1) * t278 + mrSges(11,2) * t281) * qJD(1);
t526 = pkin(15) * (mrSges(9,1) * t532 + mrSges(9,2) * t333) - t259 * (mrSges(4,1) * t206 + mrSges(4,2) * t207) - t450;
t522 = Ifges(7,6) * t472 + qJD(1) * t558;
t520 = m(11) * g(1);
t519 = m(11) * g(2);
t518 = mrSges(7,1) * g(2);
t517 = mrSges(10,1) * g(1);
t516 = mrSges(10,1) * g(2);
t515 = mrSges(3,2) * g(1);
t514 = mrSges(3,2) * g(2);
t513 = mrSges(7,2) * g(1);
t512 = mrSges(7,2) * g(2);
t510 = mrSges(10,2) * g(1);
t509 = mrSges(10,2) * g(2);
t9 = -t304 * t330 + t595;
t508 = t16 * t9;
t504 = t384 / 0.2e1;
t297 = qJ(1) + qJ(2);
t499 = sin(t297) / 0.2e1;
t372 = mrSges(10,1) * t305 + mrSges(10,2) * t311;
t228 = t372 * qJD(1);
t489 = pkin(2) * t228;
t488 = pkin(5) * t206;
t486 = m(11) * t95;
t440 = qJD(1) * pkin(15);
t102 = -pkin(2) * t333 - t440;
t484 = m(10) * t102;
t476 = t114 * t3;
t464 = Ifges(6,4) * t303;
t463 = Ifges(6,4) * t309;
t370 = Ifges(6,1) * t309 - t464;
t456 = t303 * (t97 * Ifges(6,5) + t103 * t370);
t366 = -Ifges(6,2) * t303 + t463;
t454 = t309 * (t97 * Ifges(6,6) + t103 * t366);
t444 = mrSges(11,3) * t278;
t441 = Ifges(11,4) * t281;
t439 = t114 * t303;
t438 = t114 * t309;
t434 = t213 * t310;
t427 = t311 * t321;
t156 = qJ(2) + atan2(-t217, t215);
t398 = mrSges(11,3) * t424;
t426 = -mrSges(4,3) * t206 + t293 * t569 - t398;
t420 = qJD(2) * t305;
t418 = qJD(2) * t311;
t417 = qJD(4) * t114;
t416 = qJD(4) * t303;
t415 = qJD(4) * t309;
t414 = qJDD(1) * pkin(15);
t57 = -t103 * t416 + t309 * t99;
t58 = -t103 * t415 - t303 * t99;
t96 = qJDD(4) + t100;
t407 = Ifges(6,5) * t57 + Ifges(6,6) * t58 + Ifges(6,3) * t96;
t406 = -t520 / 0.2e1;
t405 = t520 / 0.2e1;
t404 = -t519 / 0.2e1;
t401 = pkin(1) * t418;
t284 = qJ(2) + pkin(17) - pkin(18);
t274 = sin(t294);
t275 = cos(t294);
t164 = qJ(2) + atan2(-t275, t274) + atan2(-t275, -t274);
t381 = t419 * t444;
t177 = -pkin(5) * t226 + t259;
t306 = sin(qJ(1));
t312 = cos(qJ(1));
t378 = g(1) * t306 - g(2) * t312;
t369 = Ifges(7,1) * t384 + t571;
t361 = t11 * t309 + t12 * t303;
t141 = t218 * t311 - t219 * t305;
t142 = -t218 * t305 - t219 * t311;
t345 = t16 * t388;
t344 = t97 * (-Ifges(6,5) * t303 - Ifges(6,6) * t309);
t343 = t102 * (mrSges(10,1) * t311 - mrSges(10,2) * t305);
t342 = t303 * (-Ifges(6,2) * t309 - t464);
t339 = t309 * (-Ifges(6,1) * t303 - t463);
t43 = mrSges(6,1) * t96 - mrSges(6,3) * t57;
t44 = -mrSges(6,2) * t96 + mrSges(6,3) * t58;
t327 = qJD(4) * t535 - t303 * t43 + t309 * t44;
t318 = mrSges(7,1) * g(1);
t298 = qJ(1) - qJ(2);
t287 = t311 * mrSges(3,2);
t283 = cos(t298);
t282 = cos(t297);
t280 = sin(t298);
t267 = -qJ(1) + t284;
t266 = qJ(1) + t284;
t248 = -pkin(20) + t386;
t242 = cos(t248);
t241 = sin(t248);
t239 = pkin(1) * t582 + mrSges(3,1);
t234 = g(2) * t239;
t233 = t239 * g(1);
t181 = -qJ(1) + t390;
t180 = qJ(1) + t390;
t168 = pkin(1) * t421 + t488;
t161 = qJ(1) - t164;
t160 = qJ(1) + t164;
t155 = -qJD(2) * t225 - t332;
t154 = qJD(2) * t226 + t331;
t146 = t293 * t490;
t133 = t148 * qJD(2);
t132 = t384 * qJD(2);
t130 = t142 * qJD(2);
t129 = t141 * qJD(2);
t111 = cos(t117);
t82 = qJDD(1) * t95 + t246;
t79 = atan2(-t217, -t215) + t156;
t75 = qJD(1) * t133 + qJDD(1) * t384;
t74 = -qJD(1) * t132 + qJDD(1) * t148;
t69 = Ifges(7,5) * qJD(2) + qJD(1) * t369;
t65 = (qJD(2) * t130 + qJDD(2) * t141) * pkin(1);
t64 = (qJD(2) * t129 - qJDD(2) * t142) * pkin(1);
t56 = t376 * t103;
t49 = -pkin(2) * t54 - t414;
t27 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t21 = -t263 * t359 + t408;
t14 = t168 * t303 + t17 * t309;
t13 = t168 * t309 - t17 * t303;
t10 = -t403 * t587 + t592;
t8 = t10 * t309 + t303 * t488;
t7 = -t10 * t303 + t309 * t488;
t5 = [t376 * t476 + t213 * t443 + t417 * t583 + (-mrSges(4,1) * t155 + mrSges(4,2) * t154) * t243 + t345 * t417 + t114 * t99 * Ifges(5,1) + t74 * t558 + t355 * t560 + t100 * Ifges(5,2) * t351 + (mrSges(9,1) * pkin(15) + Ifges(9,4) * t137 + Ifges(9,2) * t140) * t54 + (-mrSges(9,2) * pkin(15) + Ifges(9,1) * t137 + Ifges(9,4) * t140) * t55 + t580 * t231 / 0.2e1 + t581 * t232 / 0.2e1 + (t59 * t415 + m(6) * (-qJD(4) * t534 + t1 * t303 + t2 * t309) - t60 * t416 + t309 * t43 + t303 * t44) * (pkin(9) * t351 - pkin(11) * t114 + t177) + (Ifges(4,5) * t154 + Ifges(4,6) * t155 + t281 * t175) * t559 + (t240 + t200) * t373 + (Ifges(9,6) * t559 + t584 / 0.2e1 + mrSges(9,1) * t440 - pkin(2) * t484 - t489 + t62 / 0.2e1 + Ifges(9,4) * t575) * (-qJD(2) * t137 + t353) + (t578 * (t200 + t246) - mrSges(4,1) * t92 + mrSges(4,2) * t93) * t259 + t527 * t413 + (t145 * t420 - t147 * t418 - t311 * t77) * mrSges(10,3) + (m(10) * t49 - mrSges(10,1) * t232 + mrSges(10,2) * t231) * (-pkin(2) * t140 - pkin(15)) + (t103 * t339 + t344) * t417 / 0.2e1 - (t103 * t342 + t454 + t456) * t417 / 0.2e1 + (mrSges(9,1) * t140 - mrSges(9,2) * t137 - t287) * t414 - t155 * t382 + t191 * t441 / 0.2e1 + (t246 + t82) * t486 + (Ifges(4,1) * t154 + Ifges(4,4) * t155) * t503 + (Ifges(7,1) * t75 + Ifges(7,4) * t74 + Ifges(7,5) * qJDD(2)) * t504 + (Ifges(7,4) * t75 + Ifges(7,2) * t74 + Ifges(7,6) * qJDD(2)) * t505 + (cos(t181) * t405 + cos(t180) * t406 + (sin(t181) + sin(t180)) * t404) * pkin(4) + (t281 * (-Ifges(11,2) * t278 + t441) + t336) * t422 / 0.2e1 - t132 * t522 + t553 + (-t218 * t65 - t219 * t64) * mrSges(8,3) + (mrSges(10,2) * cos(t79) + mrSges(10,1) * sin(t79) + sin(t295) * mrSges(8,2) - cos(t156) * mrSges(9,2) - cos(t295) * mrSges(8,1) - sin(t156) * t289 - t242 * (pkin(9) * m(6) + mrSges(5,1) + t388) - t241 * (pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3))) * t378 - t76 * t468 + (t280 * t404 + t282 * t405 + t283 * t406 + t293 * t381 + t499 * t519) * pkin(1) + t75 * t369 / 0.2e1 + t49 * t372 + t82 * t356 + (t472 * t537 + t343) * qJD(2) + t281 * (Ifges(11,4) * t191 + Ifges(11,2) * t192 + Ifges(11,6) * t291) / 0.2e1 - pkin(15) * (-mrSges(3,1) * t232 + mrSges(3,2) * t231) - t529 + (Ifges(6,1) * t57 + Ifges(6,4) * t58 + Ifges(6,5) * t96) * t438 / 0.2e1 - (Ifges(6,4) * t57 + Ifges(6,2) * t58 + Ifges(6,6) * t96) * t439 / 0.2e1 - t154 * t383 + t207 * (Ifges(4,4) * t154 + Ifges(4,2) * t155) / 0.2e1 + t95 * (-mrSges(11,1) * t192 + mrSges(11,2) * t191) + t154 * t110 / 0.2e1 + t155 * t109 / 0.2e1 + (pkin(14) * (mrSges(7,1) * t132 + mrSges(7,2) * t133) + (Ifges(7,1) * t133 - Ifges(7,4) * t132) * t504 + (Ifges(7,4) * t133 - Ifges(7,2) * t132) * t505) * qJD(1) + t133 * t69 / 0.2e1 + (Ifges(7,5) * t133 - Ifges(7,6) * t132) * t472 - t444 * t214 + (t562 * pkin(15) ^ 2 + Ifges(8,1) * t219 ^ 2 + Ifges(2,3) + (-0.2e1 * Ifges(8,4) * t219 + Ifges(8,2) * t218) * t218 + (-mrSges(7,1) * t148 + mrSges(7,2) * t384 + t561) * pkin(14)) * qJDD(1) + (Ifges(7,5) * t384 + Ifges(7,6) * t148 + t311 * t555) * qJDD(2) / 0.2e1 + (0.2e1 * Ifges(9,5) * t137 + Ifges(11,6) * t281 + 0.2e1 * Ifges(9,6) * t140) * t498 + t351 * t407 / 0.2e1 + (-t351 * t4 + t476) * mrSges(5,3) + t96 * (Ifges(6,3) * t351 + (Ifges(6,5) * t309 - Ifges(6,6) * t303) * t114) / 0.2e1 + t57 * (Ifges(6,5) * t351 + t114 * t370) / 0.2e1 + t58 * (Ifges(6,6) * t351 + t114 * t366) / 0.2e1 + t2 * (mrSges(6,1) * t351 - mrSges(6,3) * t438) + t1 * (-mrSges(6,2) * t351 - mrSges(6,3) * t439) + t78 * (mrSges(5,1) * t351 + mrSges(5,2) * t114) + (-t100 * t114 - t351 * t99) * Ifges(5,4) + (t354 * t559 - t377 + t450) * t293 + (Ifges(11,1) * t191 + Ifges(11,4) * t560 + Ifges(11,5) * t531 + t174 * t497) * t278 + pkin(14) * (-mrSges(7,1) * t74 + mrSges(7,2) * t75) + (m(6) * t361 - t535 - t548) * (-pkin(5) * t155 + t401) + (-g(1) * t538 + g(2) * t544) * t306 + (g(1) * t544 + g(2) * t538) * t312 - t539 * t420 / 0.2e1 - t540 * t418 / 0.2e1 + (-t200 * mrSges(4,1) + t213 * mrSges(4,3) + Ifges(4,4) * t93 + Ifges(4,2) * t92 + Ifges(4,6) * t531) * t226 + (t200 * mrSges(4,2) - t214 * mrSges(4,3) + Ifges(4,1) * t93 + Ifges(4,4) * t92 + Ifges(4,5) * t531) * t225 + t530 * t401 + (Ifges(9,5) * t559 + t116 / 0.2e1 - mrSges(9,2) * t440 + t63 / 0.2e1 + Ifges(9,1) * t575) * (qJD(2) * t140 + t352) + (m(5) * t78 + mrSges(5,1) * t100 + mrSges(5,2) * t99) * t177 + (-mrSges(3,1) * t414 - t585 * t232 / 0.2e1 - t556 * t231 / 0.2e1 - t554 * qJDD(2)) * t305 + (t555 * qJDD(2) + t231 * t586 + t556 * t232) * t493; t327 * t21 + (-t11 * t13 - t12 * t14 - t534 * t6 + (-qJD(4) * t361 + t533) * t21 + t567) * m(6) + t492 * t56 + (Ifges(10,3) + Ifges(7,3) + Ifges(3,3)) * qJDD(2) + (-t303 * t6 - t13) * t60 - t448 * t213 + (t309 * t6 - t14) * t59 + (t318 - t512) * cos(t267) / 0.2e1 + (t318 + t512) * cos(t266) / 0.2e1 + (t233 - t514) * t283 / 0.2e1 + (t233 + t514) * t282 / 0.2e1 + (t234 - t515) * t499 + (t234 + t515) * t280 / 0.2e1 + (-t510 - t516) * sin(t161) / 0.2e1 + (t510 - t516) * sin(t160) / 0.2e1 + (-t509 - t517) * cos(t160) / 0.2e1 + (t509 - t517) * cos(t161) / 0.2e1 + (-t513 - t518) * sin(t267) / 0.2e1 + (-t513 + t518) * sin(t266) / 0.2e1 + (-t165 * t168 + t18 * t574 + t21 * t4 + t567) * m(5) + (-t100 * t21 + t492 * t103 - t104 * t574 + t20 * t99) * mrSges(5,3) + (-t343 + t69 * t577 + t384 * t522 + (-pkin(14) * (mrSges(7,1) * t384 + mrSges(7,2) * t148) + (-Ifges(7,2) * t384 + t571) * t577 - t384 * (Ifges(7,1) * t148 - t552) / 0.2e1 - t527) * qJD(1) + (-t145 * t305 + t147 * t311) * mrSges(10,3) + t539 * t305 / 0.2e1 + t540 * t493 - (Ifges(7,5) * t148 - Ifges(7,6) * t384 + t537) * qJD(2) / 0.2e1 + t526) * qJD(1) + (m(4) * (t214 * t304 - t434) + (-t381 - t530 * t311 + (-t129 * t219 - t130 * t218 + (t141 * t219 + t142 * t218) * qJD(2)) * mrSges(8,3)) * qJD(1) + (-m(4) * t259 - t486) * t427 + (-mrSges(11,3) * t192 - mrSges(4,3) * t92 + t291 * t448) * t310 + (-mrSges(11,3) * t191 - mrSges(4,3) * t93 + t291 * t569) * t304 + (-t141 * t218 + t142 * t219) * qJDD(1) * mrSges(8,3) + (t310 * t426 + t545) * qJD(3) + (t230 * t573 - t434) * m(11) + (-t259 * t427 + t141 * t65 - t142 * t64 + 0.2e1 * (-t129 * t142 + t130 * t141) * t472 * pkin(1)) * m(8)) * pkin(1) + (mrSges(7,1) * sin(t284) - mrSges(10,1) * sin(t164) - mrSges(11,2) * t111 + mrSges(7,2) * cos(t284) - mrSges(10,2) * cos(t164) + t239 * t305 + t287 - t565) * g(3) - t168 * t61 - (t228 + t484) * t557 + t554 * t232 + t555 * t231 - t76 * mrSges(10,2) + t77 * mrSges(10,1) + Ifges(7,6) * t74 + Ifges(7,5) * t75 + t563 + t20 * t27; (-t545 + (-t398 - t426) * t310) * t471 - t532 * t489 + (-g(3) * t111 - t213) * mrSges(11,2) + t526 * qJD(1) - m(10) * (t102 * t557 + (t145 - t146) * t147) - m(5) * (t10 * t18 + t508) - m(6) * (t11 * t7 + t12 * t8 + t508) - t565 * g(3) + t146 * t235 - t147 * t236 - t213 * mrSges(4,2) - t7 * t60 - t8 * t59 - t9 * t56 + (-t25 * t56 + t587 * t27 + (-t103 * t25 + t587 * t99) * mrSges(5,3) + t319 * (-t16 * t25 + t3 * t587) + (m(5) * t18 - m(6) * t534 - t104 * mrSges(5,3) - t303 * t60 + t309 * t59) * t24 + t548 * t206 - (t327 + m(6) * (-t11 * t415 - t12 * t416 + t533) + m(5) * t4 - t100 * mrSges(5,3)) * t359) * pkin(5) + (t10 * t104 - t103 * t9) * mrSges(5,3) + t563; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t12 * t60 - t11 * t59 - t388 * t378 + (-(g(1) * t312 + g(2) * t306) * t242 - g(3) * t241) * t376 + (-t344 / 0.2e1 + t454 / 0.2e1 + t456 / 0.2e1 - t345 + (t342 / 0.2e1 - t339 / 0.2e1) * t103 - t583) * t103 + t407;];
tau = t5;
