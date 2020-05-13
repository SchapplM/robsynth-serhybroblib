% Calculate vector of inverse dynamics joint torques with ic for
% palh3m2IC
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m2IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2IC_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:59:05
% EndTime: 2020-05-07 04:59:43
% DurationCPUTime: 34.74s
% Computational Cost: add. (11544->825), mult. (26205->1129), div. (14->8), fcn. (20838->49), ass. (0->392)
t325 = qJ(2) + qJ(3);
t318 = qJ(4) + t325;
t295 = sin(t318);
t296 = cos(t318);
t337 = cos(qJ(5));
t531 = mrSges(6,1) * t337;
t634 = -mrSges(5,2) * t296 + (-pkin(8) * m(6) - mrSges(5,1) - t531) * t295;
t341 = cos(qJ(1));
t596 = t634 * t341;
t314 = sin(t325);
t316 = cos(t325);
t612 = mrSges(4,1) * t314 + mrSges(4,2) * t316;
t631 = -t612 * t341 + t596;
t334 = sin(qJ(1));
t595 = t634 * t334;
t630 = -t612 * t334 + t595;
t333 = sin(qJ(2));
t340 = cos(qJ(2));
t265 = -mrSges(3,1) * t340 + mrSges(3,2) * t333;
t424 = -t296 * mrSges(5,1) + t295 * mrSges(5,2);
t425 = -mrSges(4,1) * t316 + t314 * mrSges(4,2);
t324 = qJ(2) + qJ(7);
t313 = sin(t324);
t315 = cos(t324);
t460 = pkin(15) - qJ(7);
t434 = -qJ(8) + t460;
t297 = -qJ(2) + t434;
t283 = sin(t297);
t284 = cos(t297);
t405 = t284 * mrSges(9,1) + t283 * mrSges(9,2);
t614 = -t315 * mrSges(8,1) + mrSges(8,2) * t313 + t405;
t629 = t614 + t265 - t424 - t425;
t331 = sin(qJ(4));
t332 = sin(qJ(3));
t535 = pkin(1) * qJD(2);
t451 = t332 * t535;
t276 = t331 * t451;
t338 = cos(qJ(4));
t323 = qJD(2) + qJD(3);
t339 = cos(qJ(3));
t449 = t339 * t535;
t387 = -pkin(4) * t323 + t449;
t184 = -t338 * t387 + t276;
t312 = qJD(4) + t323;
t178 = -t312 * pkin(8) - t184;
t330 = sin(qJ(5));
t243 = t332 * t333 - t339 * t340;
t223 = t243 * qJD(1);
t390 = t332 * t340 + t333 * t339;
t225 = t390 * qJD(1);
t393 = t223 * t331 - t338 * t225;
t127 = t312 * t337 - t330 * t393;
t128 = t312 * t330 + t337 * t393;
t615 = -mrSges(5,1) * t312 - mrSges(6,1) * t127 + mrSges(6,2) * t128 + mrSges(5,3) * t393;
t628 = m(6) * t178 + t615;
t627 = m(4) + m(8);
t328 = sin(qJ(7));
t335 = cos(qJ(7));
t242 = -t328 * t333 + t335 * t340;
t222 = t242 * qJD(1);
t244 = t328 * t340 + t333 * t335;
t224 = t244 * qJD(1);
t319 = t340 * pkin(1);
t294 = t319 + pkin(12);
t266 = t294 * qJD(1);
t326 = sin(pkin(15));
t327 = cos(pkin(15));
t133 = -t266 + (-t222 * t327 - t224 * t326) * pkin(3);
t559 = sin(qJ(8));
t560 = cos(qJ(8));
t238 = t326 * t559 + t327 * t560;
t360 = -t326 * t560 + t327 * t559;
t111 = -t222 * t360 - t224 * t238;
t421 = -t238 * t222 + t224 * t360;
t65 = -mrSges(9,1) * t421 + mrSges(9,2) * t111;
t625 = m(9) * t133 + t65;
t183 = -pkin(4) * t223 - t266;
t420 = t338 * t223 + t225 * t331;
t92 = -mrSges(5,1) * t420 + mrSges(5,2) * t393;
t624 = m(5) * t183 + t92;
t623 = -mrSges(4,1) * t223 - mrSges(8,1) * t222 - mrSges(4,2) * t225 + mrSges(8,2) * t224;
t592 = g(1) * t341 + g(2) * t334;
t389 = t331 * t387 + t338 * t451;
t179 = pkin(10) * t312 - t389;
t73 = -pkin(8) * t420 - pkin(10) * t393 + t183;
t53 = -t179 * t330 + t337 * t73;
t54 = t179 * t337 + t330 * t73;
t397 = t330 * t54 + t337 * t53;
t155 = qJD(5) - t420;
t79 = -mrSges(6,2) * t155 + mrSges(6,3) * t127;
t80 = mrSges(6,1) * t155 - mrSges(6,3) * t128;
t621 = m(6) * t397 + t330 * t79 + t337 * t80;
t138 = -mrSges(5,2) * t312 + mrSges(5,3) * t420;
t620 = -t330 * t80 + t337 * t79 + t138;
t321 = qJDD(2) + qJDD(3);
t310 = qJDD(4) + t321;
t458 = qJD(1) * qJD(2);
t250 = qJDD(1) * t340 - t333 * t458;
t251 = qJDD(1) * t333 + t340 * t458;
t370 = t243 * qJD(3);
t136 = qJD(1) * t370 - t250 * t332 - t251 * t339;
t371 = t390 * qJD(3);
t137 = qJD(1) * t371 - t250 * t339 + t251 * t332;
t51 = qJD(4) * t420 + t136 * t338 + t137 * t331;
t19 = qJD(5) * t127 + t310 * t330 + t337 * t51;
t20 = -qJD(5) * t128 + t310 * t337 - t330 * t51;
t10 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t557 = pkin(1) * t332;
t448 = qJD(3) * t557;
t494 = pkin(1) * qJDD(2);
t235 = qJD(2) * t448 - t339 * t494;
t197 = pkin(4) * t321 + t235;
t466 = qJD(3) * t339;
t236 = (-qJD(2) * t466 - qJDD(2) * t332) * pkin(1);
t102 = qJD(4) * t389 + t338 * t197 - t331 * t236;
t96 = -t310 * pkin(8) - t102;
t619 = m(6) * t96 + t10;
t469 = qJD(1) * t333;
t305 = pkin(1) * t469;
t618 = -g(3) * t319 + t266 * t305;
t549 = pkin(4) * t314;
t556 = pkin(1) * t333;
t254 = t549 - t556;
t527 = mrSges(6,2) * t330;
t444 = t295 * t527;
t386 = -t296 * mrSges(6,3) - t444;
t530 = mrSges(9,1) * t283;
t545 = pkin(10) * t296;
t317 = -qJ(2) + t460;
t553 = pkin(3) * sin(t317);
t617 = -m(6) * (t254 - t545) - t386 - m(9) * (t553 - t556) + t530 + m(4) * t556 - m(5) * t254;
t101 = qJD(4) * t184 + t331 * t197 + t338 * t236;
t154 = Ifges(5,4) * t420;
t408 = mrSges(6,1) * t330 + mrSges(6,2) * t337;
t378 = t178 * t408;
t399 = Ifges(6,5) * t337 - Ifges(6,6) * t330;
t52 = -qJD(4) * t393 - t136 * t331 + t137 * t338;
t461 = qJDD(1) * pkin(12);
t206 = -pkin(1) * t250 - t461;
t99 = -pkin(4) * t137 + t206;
t13 = -pkin(8) * t52 - pkin(10) * t51 + t99;
t95 = pkin(10) * t310 + t101;
t4 = qJD(5) * t53 + t13 * t330 + t337 * t95;
t518 = Ifges(6,4) * t337;
t402 = -Ifges(6,2) * t330 + t518;
t519 = Ifges(6,4) * t330;
t404 = Ifges(6,1) * t337 - t519;
t463 = qJD(5) * t330;
t432 = -t463 / 0.2e1;
t126 = Ifges(6,4) * t127;
t47 = t128 * Ifges(6,1) + t155 * Ifges(6,5) + t126;
t497 = t337 * t47;
t438 = t497 / 0.2e1;
t45 = t128 * Ifges(6,5) + t127 * Ifges(6,6) + t155 * Ifges(6,3);
t512 = t128 * Ifges(6,4);
t46 = t127 * Ifges(6,2) + t155 * Ifges(6,6) + t512;
t462 = qJD(5) * t337;
t5 = -qJD(5) * t54 + t13 * t337 - t330 * t95;
t509 = t389 * mrSges(5,3);
t510 = t184 * mrSges(5,3);
t511 = t393 * Ifges(5,4);
t523 = mrSges(6,3) * t337;
t524 = mrSges(6,3) * t330;
t562 = t330 / 0.2e1;
t564 = -t312 / 0.2e1;
t570 = t393 / 0.2e1;
t571 = -t393 / 0.2e1;
t572 = -t420 / 0.2e1;
t573 = -t155 / 0.2e1;
t575 = -t128 / 0.2e1;
t576 = -t127 / 0.2e1;
t50 = qJDD(5) - t52;
t580 = t50 / 0.2e1;
t582 = t20 / 0.2e1;
t583 = t19 / 0.2e1;
t603 = t54 * mrSges(6,2);
t604 = t53 * mrSges(6,1);
t613 = t527 - t531;
t8 = t19 * Ifges(6,4) + t20 * Ifges(6,2) + t50 * Ifges(6,6);
t85 = Ifges(5,2) * t420 + t312 * Ifges(5,6) + t511;
t86 = Ifges(5,1) * t393 + t312 * Ifges(5,5) + t154;
t9 = t19 * Ifges(6,1) + t20 * Ifges(6,4) + t50 * Ifges(6,5);
t616 = -t101 * mrSges(5,2) - t5 * t524 + t102 * mrSges(5,1) + (Ifges(6,1) * t330 + t518) * t583 + (Ifges(6,2) * t337 + t519) * t582 + Ifges(5,3) * t310 + t4 * t523 + (Ifges(6,5) * t330 + Ifges(6,6) * t337) * t580 + t46 * t432 + Ifges(5,6) * t52 + Ifges(5,5) * t51 + t9 * t562 + t337 * t8 / 0.2e1 + t96 * t613 + (-t462 * t53 - t463 * t54) * mrSges(6,3) + (t378 + t438) * qJD(5) + (t127 * t402 + t128 * t404 + t155 * t399) * qJD(5) / 0.2e1 + t85 * t570 + (t154 + t86) * t572 + (-t511 + t45) * t571 + (-t183 * mrSges(5,2) + Ifges(5,1) * t571 + Ifges(5,5) * t564 + t399 * t573 + t402 * t576 + t404 * t575 + t523 * t53 + t524 * t54 - t378 + t46 * t562 - t497 / 0.2e1 + t510) * t420 + (-t183 * mrSges(5,1) + Ifges(6,5) * t575 - Ifges(5,2) * t572 - Ifges(5,6) * t564 + Ifges(6,6) * t576 + Ifges(6,3) * t573 - t509 + t603 - t604) * t393;
t11 = mrSges(6,1) * t50 - mrSges(6,3) * t19;
t12 = -mrSges(6,2) * t50 + mrSges(6,3) * t20;
t588 = m(6) * (-qJD(5) * t397 - t330 * t5 + t337 * t4) - t80 * t462 - t79 * t463 + t337 * t12 - t330 * t11;
t322 = qJD(2) + qJD(7);
t311 = qJD(8) + t322;
t565 = -t311 / 0.2e1;
t577 = -t111 / 0.2e1;
t611 = -t133 * mrSges(9,2) + Ifges(9,1) * t577 + Ifges(9,5) * t565;
t106 = Ifges(9,4) * t421;
t372 = t242 * qJD(7);
t134 = qJD(1) * t372 + t250 * t328 + t251 * t335;
t373 = t244 * qJD(7);
t135 = -qJD(1) * t373 + t250 * t335 - t251 * t328;
t216 = t238 * qJD(8);
t217 = t360 * qJD(8);
t29 = -t134 * t238 - t135 * t360 - t216 * t222 + t217 * t224;
t30 = t134 * t360 - t135 * t238 + t216 * t224 + t217 * t222;
t320 = qJDD(2) + qJDD(7);
t309 = qJDD(8) + t320;
t452 = t328 * t535;
t552 = pkin(3) * t326;
t229 = t322 * t552 + t452;
t450 = t335 * t535;
t551 = pkin(3) * t327;
t230 = t322 * t551 + t450;
t117 = t229 * t238 + t230 * t360;
t513 = t117 * mrSges(9,3);
t534 = pkin(1) * qJD(7);
t443 = qJD(2) * t534;
t233 = -t328 * t443 + t335 * t494;
t194 = t320 * t551 + t233;
t234 = t328 * t494 + t335 * t443;
t195 = t320 * t552 + t234;
t56 = -t194 * t360 - t195 * t238 - t216 * t230 + t217 * t229;
t57 = -t194 * t238 + t195 * t360 + t216 * t229 + t217 * t230;
t578 = -t421 / 0.2e1;
t59 = t111 * Ifges(9,4) + Ifges(9,2) * t421 + t311 * Ifges(9,6);
t60 = t111 * Ifges(9,1) + t311 * Ifges(9,5) + t106;
t609 = t57 * mrSges(9,1) - t56 * mrSges(9,2) + Ifges(9,5) * t29 + Ifges(9,6) * t30 + Ifges(9,3) * t309 + (t60 + t106) * t578 - (t133 * mrSges(9,1) + Ifges(9,4) * t577 + Ifges(9,2) * t578 + Ifges(9,6) * t565 + t513 - t59 / 0.2e1) * t111;
t606 = t250 / 0.2e1;
t329 = sin(qJ(6));
t563 = t329 / 0.2e1;
t602 = t340 * Ifges(3,2);
t380 = pkin(1) * (t238 * t328 + t335 * t360);
t601 = -qJD(2) * t380 + (t216 * t326 + t217 * t327) * pkin(3);
t379 = pkin(1) * (-t238 * t335 + t328 * t360);
t600 = -qJD(2) * t379 + (-t216 * t327 + t217 * t326) * pkin(3);
t456 = qJ(4) + pkin(14);
t431 = qJ(3) + t456;
t412 = pkin(15) + t431;
t476 = qJ(2) - qJ(6);
t374 = t412 - t476;
t354 = -0.2e1 * qJ(7) - pkin(16) + t374;
t375 = t412 + t476;
t357 = pkin(16) + t375;
t599 = -cos(qJ(8) - t354) + cos(qJ(8) - t357);
t598 = -m(6) * (-t545 + t549) - t386 - m(5) * t549;
t506 = t295 * mrSges(6,3);
t368 = t296 * t613 - t506;
t597 = -t368 - t424;
t173 = t243 * t331 - t338 * t390;
t176 = qJD(2) * t243 + t370;
t177 = qJD(2) * t390 + t371;
t392 = t338 * t243 + t331 * t390;
t76 = qJD(4) * t392 + t176 * t338 + t177 * t331;
t384 = t173 * t462 + t330 * t76;
t591 = Ifges(3,6) * qJDD(2);
t336 = cos(qJ(6));
t264 = -mrSges(7,1) * t336 + mrSges(7,2) * t329;
t367 = m(7) * pkin(6) + t264;
t590 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t93 = pkin(8) * t393 - pkin(10) * t420;
t589 = -mrSges(9,3) - mrSges(8,3) - mrSges(5,3) - mrSges(4,3) - mrSges(3,3) - mrSges(7,3) + mrSges(2,2);
t548 = pkin(4) * t316;
t430 = t319 - t548;
t247 = pkin(12) + t430;
t411 = pkin(8) * t296 + pkin(10) * t295;
t281 = pkin(3) * cos(t317);
t470 = t281 + t319;
t587 = m(6) * (-t247 + t411) - mrSges(2,1) - m(9) * (pkin(12) + t470) - m(5) * t247 - m(3) * pkin(12) + t367 + t506 - t627 * t294 + t629;
t574 = t128 / 0.2e1;
t567 = t224 / 0.2e1;
t566 = -t225 / 0.2e1;
t561 = t340 / 0.2e1;
t558 = pkin(1) * t328;
t555 = pkin(1) * t335;
t554 = pkin(1) * t339;
t550 = pkin(4) * t225;
t547 = pkin(4) * t331;
t546 = pkin(4) * t338;
t355 = -t444 + (-pkin(10) * m(6) - mrSges(6,3)) * t296;
t71 = -t184 * t330 + t337 * t93;
t72 = t184 * t337 + t330 * t93;
t539 = (-t72 * t79 - t71 * t80 - m(6) * (t53 * t71 + t54 * t72) - t184 * t138 + t628 * t389 + (m(6) * t411 + t597) * g(3) + (-t334 * t355 + t595) * g(2) + (-t341 * t355 + t596) * g(1) - t619 * pkin(8) + t588 * pkin(10) + t616) / pkin(9);
t116 = t229 * t360 - t230 * t238;
t525 = mrSges(9,2) * t284;
t252 = t334 * t525;
t253 = t341 * t525;
t515 = t421 * mrSges(9,3);
t97 = -mrSges(9,2) * t311 + t515;
t98 = mrSges(9,1) * t311 - mrSges(9,3) * t111;
t538 = (-t117 * t98 + g(3) * t405 - g(1) * (-t341 * t530 + t253) - g(2) * (-t334 * t530 + t252) + (-t97 + t515) * t116 + t611 * t421 + t609) / pkin(7);
t536 = -qJD(6) / 0.2e1;
t522 = mrSges(8,3) * t224;
t521 = Ifges(3,4) * t333;
t520 = Ifges(3,4) * t340;
t517 = Ifges(7,4) * t329;
t516 = Ifges(7,4) * t336;
t514 = t116 * mrSges(9,3);
t508 = t224 * Ifges(8,4);
t507 = t225 * Ifges(4,4);
t504 = t329 * Ifges(7,1);
t499 = t336 * Ifges(7,2);
t493 = t173 * t330;
t492 = t173 * t337;
t490 = (-mrSges(8,2) * t322 + mrSges(8,3) * t222) * t335;
t401 = t499 + t517;
t484 = t329 * (Ifges(7,6) * qJD(6) + qJD(1) * t401);
t483 = t330 * t334;
t482 = t330 * t341;
t481 = t331 * t332;
t480 = t332 * t338;
t479 = t334 * t337;
t303 = qJD(1) * t516;
t478 = t336 * (Ifges(7,5) * qJD(6) + qJD(1) * t504 + t303);
t477 = t337 * t341;
t468 = qJD(1) * t340;
t467 = qJD(2) * t333;
t465 = qJD(4) * t331;
t464 = qJD(4) * t338;
t459 = qJ(7) + pkin(16);
t457 = qJD(1) * qJD(6);
t453 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t50;
t306 = pkin(1) * t467;
t160 = -pkin(4) * t177 + t306;
t427 = t458 / 0.2e1;
t419 = mrSges(4,3) * t451;
t418 = mrSges(4,3) * t449;
t417 = mrSges(8,3) * t452;
t416 = mrSges(8,3) * t450;
t413 = t459 + t476;
t196 = -pkin(4) * t243 - t294;
t300 = pkin(4) - t554;
t213 = -pkin(1) * t480 + t331 * t300;
t409 = mrSges(3,1) * t333 + mrSges(3,2) * t340;
t407 = mrSges(7,1) * t329 + mrSges(7,2) * t336;
t406 = -mrSges(8,1) * t313 - mrSges(8,2) * t315;
t403 = t521 + t602;
t400 = Ifges(3,5) * t340 - Ifges(3,6) * t333;
t398 = Ifges(7,5) * t336 - Ifges(7,6) * t329;
t396 = -t330 * t53 + t337 * t54;
t391 = t331 * t339 + t480;
t212 = pkin(1) * t481 + t300 * t338;
t383 = t173 * t463 - t337 * t76;
t382 = pkin(6) * t407;
t381 = pkin(12) * t409;
t377 = t329 * (Ifges(7,1) * t336 - t517);
t376 = t333 * (Ifges(3,1) * t340 - t521);
t152 = (-t222 * t326 + t224 * t327) * pkin(3);
t78 = -t550 + t93;
t369 = -t411 - t548;
t359 = -qJ(7) + t374;
t358 = -qJ(7) + t375;
t148 = t222 * Ifges(8,2) + t322 * Ifges(8,6) + t508;
t209 = Ifges(8,4) * t222;
t150 = t224 * Ifges(8,1) + t322 * Ifges(8,5) + t209;
t347 = -t234 * mrSges(8,2) + t266 * (mrSges(8,1) * t224 + mrSges(8,2) * t222) + t148 * t567 - t224 * (Ifges(8,1) * t222 - t508) / 0.2e1 + Ifges(8,6) * t135 + Ifges(8,5) * t134 - t322 * (Ifges(8,5) * t222 - Ifges(8,6) * t224) / 0.2e1 + t222 * t416 + t233 * mrSges(8,1) + Ifges(8,3) * t320 - (-Ifges(8,2) * t224 + t150 + t209) * t222 / 0.2e1 + (t514 + t611) * t421 + t609;
t149 = t223 * Ifges(4,2) + t323 * Ifges(4,6) - t507;
t210 = Ifges(4,4) * t223;
t151 = -t225 * Ifges(4,1) + t323 * Ifges(4,5) + t210;
t346 = t266 * (-mrSges(4,1) * t225 + mrSges(4,2) * t223) + Ifges(4,5) * t136 + Ifges(4,6) * t137 + t149 * t566 - t223 * t418 + t225 * (Ifges(4,1) * t223 + t507) / 0.2e1 + t225 * t419 + t235 * mrSges(4,1) - t236 * mrSges(4,2) + Ifges(4,3) * t321 - t323 * (Ifges(4,5) * t223 + Ifges(4,6) * t225) / 0.2e1 - (Ifges(4,2) * t225 + t151 + t210) * t223 / 0.2e1 + t616;
t345 = 0.1e1 / pkin(2);
t304 = Ifges(3,4) * t468;
t299 = -pkin(8) - t546;
t282 = sin(t413);
t268 = t551 + t555;
t267 = t552 + t558;
t249 = qJDD(1) * t329 + t336 * t457;
t248 = qJDD(1) * t336 - t329 * t457;
t221 = Ifges(3,1) * t469 + Ifges(3,5) * qJD(2) + t304;
t219 = Ifges(3,6) * qJD(2) + qJD(1) * t403;
t208 = -t338 * t449 + t276;
t207 = t391 * t535;
t203 = t296 * t477 - t483;
t202 = t296 * t482 + t479;
t201 = t296 * t479 + t482;
t200 = t296 * t483 - t477;
t190 = mrSges(4,1) * t323 + mrSges(4,3) * t225;
t189 = mrSges(8,1) * t322 - t522;
t188 = -mrSges(4,2) * t323 + mrSges(4,3) * t223;
t175 = -qJD(2) * t244 - t373;
t174 = qJD(2) * t242 + t372;
t166 = (-t238 * t326 - t327 * t360) * pkin(3);
t165 = (-t238 * t327 + t326 * t360) * pkin(3);
t141 = -t238 * t267 - t268 * t360;
t140 = -t238 * t268 + t267 * t360;
t82 = qJD(7) * t380 + t216 * t267 + t217 * t268;
t81 = qJD(7) * t379 - t216 * t268 + t217 * t267;
t77 = qJD(4) * t173 + t176 * t331 - t338 * t177;
t67 = t208 * t337 + t330 * t78;
t66 = -t208 * t330 + t337 * t78;
t63 = (-t134 * t326 - t135 * t327) * pkin(3) + t206;
t44 = t174 * t360 - t175 * t238 + t216 * t244 + t217 * t242;
t43 = -t174 * t238 - t175 * t360 - t216 * t242 + t217 * t244;
t39 = -mrSges(5,2) * t310 + mrSges(5,3) * t52;
t38 = mrSges(5,1) * t310 - mrSges(5,3) * t51;
t22 = -mrSges(9,2) * t309 + mrSges(9,3) * t30;
t21 = mrSges(9,1) * t309 - mrSges(9,3) * t29;
t1 = [(-mrSges(4,1) * t206 + mrSges(4,3) * t236 + Ifges(4,4) * t136 + Ifges(4,2) * t137 + Ifges(4,6) * t321) * t243 + t251 * t520 / 0.2e1 + (-mrSges(8,1) * t206 + mrSges(8,3) * t234 + Ifges(8,4) * t134 + Ifges(8,2) * t135 + Ifges(8,6) * t320) * t242 + (m(5) * t99 - mrSges(5,1) * t52 + mrSges(5,2) * t51) * t196 + (t99 * mrSges(5,2) - t102 * mrSges(5,3) + Ifges(5,1) * t51 + Ifges(5,4) * t52 + Ifges(5,5) * t310 + t399 * t580 + t402 * t582 + t404 * t583 + t408 * t96 + t432 * t47) * t173 + (m(3) * pkin(12) ^ 2 + pkin(6) * t367 + Ifges(2,3)) * qJDD(1) + t421 * (Ifges(9,4) * t43 + Ifges(9,2) * t44) / 0.2e1 + t420 * (Ifges(5,4) * t76 - Ifges(5,2) * t77) / 0.2e1 + (m(9) * t63 - mrSges(9,1) * t30 + mrSges(9,2) * t29) * ((-t242 * t327 - t244 * t326) * pkin(3) - t294) + (t377 + t336 * (-Ifges(7,2) * t329 + t516)) * t457 / 0.2e1 + t624 * t160 + t625 * (t306 + (-t174 * t326 - t175 * t327) * pkin(3)) + t623 * t306 + t249 * (t504 + t516) / 0.2e1 + (t520 * t427 + t591 / 0.2e1) * t340 + (Ifges(3,4) * t251 + Ifges(3,2) * t250 + t591) * t561 - t384 * t46 / 0.2e1 + t174 * t150 / 0.2e1 + t175 * t148 / 0.2e1 + t133 * (-mrSges(9,1) * t44 + mrSges(9,2) * t43) + t77 * t509 - t44 * t513 + t111 * (Ifges(9,1) * t43 + Ifges(9,4) * t44) / 0.2e1 - t174 * t416 - t77 * t85 / 0.2e1 + t76 * t86 / 0.2e1 + t77 * t45 / 0.2e1 + t621 * (pkin(8) * t77 - pkin(10) * t76 + t160) + (t203 * mrSges(6,1) - t202 * mrSges(6,2) + t334 * t589 + t341 * t587) * g(2) + (-t201 * mrSges(6,1) + t200 * mrSges(6,2) - t334 * t587 + t341 * t589) * g(1) + (mrSges(4,1) * t137 + mrSges(8,1) * t135 - mrSges(4,2) * t136 - mrSges(8,2) * t134) * t294 - t76 * t510 + t484 * t536 + (t478 / 0.2e1 + t398 * qJD(6) / 0.2e1) * qJD(6) + (Ifges(4,1) * t176 + Ifges(4,4) * t177) * t566 + (Ifges(8,1) * t174 + Ifges(8,4) * t175) * t567 + (Ifges(5,1) * t76 - Ifges(5,4) * t77) * t570 + (-Ifges(6,1) * t383 - Ifges(6,4) * t384 + Ifges(6,5) * t77) * t574 - (-mrSges(4,1) * t177 - mrSges(8,1) * t175 + mrSges(4,2) * t176 + mrSges(8,2) * t174) * t266 + t44 * t59 / 0.2e1 + t43 * t60 / 0.2e1 - t43 * t514 + (t383 * t53 - t384 * t54 - t4 * t493 - t492 * t5) * mrSges(6,3) - t219 * t467 / 0.2e1 - t265 * t461 + t76 * t438 + t376 * t427 + (t221 * t561 + t400 * qJD(2) / 0.2e1) * qJD(2) - t177 * t419 + t175 * t417 + t176 * t418 + t336 * (Ifges(7,4) * t249 + Ifges(7,2) * t248) / 0.2e1 - (t99 * mrSges(5,1) - Ifges(5,4) * t51 - Ifges(5,2) * t52 + t453 / 0.2e1 + Ifges(6,3) * t580 + Ifges(6,6) * t582 + Ifges(6,5) * t583 - t101 * mrSges(5,3) - Ifges(5,6) * t310 + t590) * t392 + (t79 * t462 - t80 * t463 + t330 * t12 + t337 * t11 + m(6) * (qJD(5) * t396 + t330 * t4 + t337 * t5)) * (-pkin(8) * t392 - pkin(10) * t173 + t196) + (-mrSges(9,1) * t63 + mrSges(9,3) * t56 + Ifges(9,4) * t29 + Ifges(9,2) * t30 + Ifges(9,6) * t309) * (-t238 * t242 + t244 * t360) + (mrSges(9,2) * t63 - mrSges(9,3) * t57 + Ifges(9,1) * t29 + Ifges(9,4) * t30 + Ifges(9,5) * t309) * (-t238 * t244 - t242 * t360) - (mrSges(4,2) * t206 - mrSges(4,3) * t235 + Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * t321) * t390 - t381 * t458 + t382 * t457 + t9 * t492 / 0.2e1 - t8 * t493 / 0.2e1 + t176 * t151 / 0.2e1 + t627 * (-t206 * t294 - t266 * t306) + t177 * t149 / 0.2e1 + t183 * (mrSges(5,1) * t77 + mrSges(5,2) * t76) - t77 * t603 + t222 * (Ifges(8,4) * t174 + Ifges(8,2) * t175) / 0.2e1 + t223 * (Ifges(4,4) * t176 + Ifges(4,2) * t177) / 0.2e1 + (0.2e1 * Ifges(7,5) * t563 + Ifges(7,6) * t336) * qJDD(6) + (Ifges(3,1) * t251 + Ifges(3,4) * t606 + Ifges(3,5) * qJDD(2) - t427 * t602) * t333 + t77 * t604 + t403 * t606 + (mrSges(8,2) * t206 - mrSges(8,3) * t233 + Ifges(8,1) * t134 + Ifges(8,4) * t135 + Ifges(8,5) * t320) * t244 + pkin(6) * (-mrSges(7,1) * t248 + mrSges(7,2) * t249) - pkin(12) * (-mrSges(3,1) * t250 + mrSges(3,2) * t251) + t311 * (Ifges(9,5) * t43 + Ifges(9,6) * t44) / 0.2e1 + t312 * (Ifges(5,5) * t76 - Ifges(5,6) * t77) / 0.2e1 + t322 * (Ifges(8,5) * t174 + Ifges(8,6) * t175) / 0.2e1 + t323 * (Ifges(4,5) * t176 + Ifges(4,6) * t177) / 0.2e1 + (Ifges(7,1) * t249 + Ifges(7,4) * t248) * t563 + t178 * (mrSges(6,1) * t384 - mrSges(6,2) * t383) + t155 * (-Ifges(6,5) * t383 - Ifges(6,6) * t384 + Ifges(6,3) * t77) / 0.2e1 + t127 * (-Ifges(6,4) * t383 - Ifges(6,2) * t384 + Ifges(6,6) * t77) / 0.2e1 + t248 * t401 / 0.2e1; ((-t235 * t339 - t236 * t332) * pkin(1) + t618) * m(4) + ((t233 * t335 + t234 * t328) * pkin(1) + t618) * m(8) + (-t328 * t189 + t490) * t534 + t588 * (pkin(10) + t213) + ((-pkin(2) * t282 - pkin(1) * sin(t476)) * t345 * (t166 * t22 + t165 * t21 - t152 * t65 + t601 * t98 + t600 * t97 + t614 * g(3) + (-t490 + (t189 + t522) * t328) * t535 + t347 - g(1) * t253 - g(2) * t252 + t592 * (-m(9) * t553 - t406 + t530) + (-t281 * g(3) + t116 * t601 - t117 * t600 - t133 * t152 + t165 * t57 + t166 * t56) * m(9)) + pkin(1) * sin(t459) / pkin(5) * (Ifges(7,5) * t249 + Ifges(7,6) * t248 + Ifges(7,3) * qJDD(6) + g(3) * t264 + (-t478 / 0.2e1 + t484 / 0.2e1 - t336 * t303 / 0.2e1 + t398 * t536 + (-t382 - t377 / 0.2e1 + t499 * t563) * qJD(1)) * qJD(1) + t592 * t407)) / t282 - (-Ifges(3,2) * t469 + t221 + t304) * t468 / 0.2e1 + (t381 - t376 / 0.2e1) * qJD(1) ^ 2 + m(5) * (t101 * t213 + t102 * t212) - t624 * (t305 - t550) - t625 * (t305 + t152) + (m(5) * t184 - t628) * (-t300 * t465 + (qJD(3) * t391 + t332 * t464) * pkin(1)) + (-m(9) * t470 - m(5) * t430 - m(6) * (t319 + t369) - t368 + t629) * g(3) + (t334 * t617 - t252 + t630) * g(2) + (t341 * t617 - t253 + t631) * g(1) + t592 * t409 + t592 * (m(8) * t556 - t406) - t623 * t305 + t140 * t21 + t141 * t22 + t81 * t97 + t82 * t98 + t619 * (-pkin(8) - t212) + (-m(5) * t389 + m(6) * t396 + t620) * (t300 * t464 + (t332 * t465 + (-t338 * t339 + t481) * qJD(3)) * pkin(1)) - t621 * (t305 + t78) + (pkin(3) * (pkin(1) * (cos(-qJ(8) + t476) - cos(qJ(8) + t476)) + (cos(-qJ(8) + t413) - cos(qJ(8) + t413)) * pkin(2)) * t539 + (((cos(t359) - cos(t358)) * pkin(1) + (cos(t354) - cos(t357)) * pkin(2)) * pkin(3) + ((-cos(-qJ(8) + t359) + cos(-qJ(8) + t358)) * pkin(1) + t599 * pkin(2)) * pkin(7)) * t538) / t599 * t345 + (mrSges(8,1) * t320 - mrSges(8,3) * t134) * t555 + (-mrSges(8,2) * t320 + mrSges(8,3) * t135) * t558 + t347 + t346 + t219 * t469 / 0.2e1 + t224 * t417 + m(9) * (t116 * t82 - t117 * t81 + t140 * t57 + t141 * t56) - t400 * t458 / 0.2e1 - (mrSges(4,1) * t321 - mrSges(4,3) * t136) * t554 + t190 * t448 + t212 * t38 + t213 * t39 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t250 + Ifges(3,5) * t251 - pkin(1) * t188 * t466 - (-mrSges(4,2) * t321 + mrSges(4,3) * t137) * t557; -pkin(9) * t539 - t67 * t79 - t66 * t80 + t38 * t546 + t39 * t547 + t92 * t550 + t346 - t190 * t451 - m(5) * (-t183 * t550 + t184 * t207 - t208 * t389) + t188 * t449 - t208 * t138 + t299 * t10 + (t178 * t207 + t299 * t96 - t53 * t66 - t54 * t67) * m(6) + (m(5) * t548 - m(6) * t369 - t425 + t597) * g(3) + (t334 * t598 + t630) * g(2) + (t341 * t598 + t631) * g(1) + t588 * (pkin(10) + t547) + t615 * (pkin(4) * t465 + t207) + (t620 * t464 - (cos(0.2e1 * qJ(3) + t456) - cos(0.2e1 * qJ(7) - 0.2e1 * pkin(15) + 0.2e1 * qJ(8) + t456)) / (cos(0.2e1 * t431) - cos(0.2e1 * t434)) * t539 + m(5) * (t101 * t331 + t102 * t338 + (-t184 * t331 - t338 * t389) * qJD(4)) + sin(t456) / sin(-qJ(7) - qJ(8) + t412) * t538 + m(6) * (t178 * t331 + t338 * t396) * qJD(4)) * pkin(4); -t178 * (mrSges(6,1) * t128 + mrSges(6,2) * t127) + (Ifges(6,1) * t127 - t512) * t575 + t46 * t574 + (Ifges(6,5) * t127 - Ifges(6,6) * t128) * t573 - t53 * t79 + t54 * t80 - g(1) * (mrSges(6,1) * t202 + mrSges(6,2) * t203) - g(2) * (mrSges(6,1) * t200 + mrSges(6,2) * t201) - g(3) * t408 * t295 + (t127 * t53 + t128 * t54) * mrSges(6,3) + t453 + (-Ifges(6,2) * t128 + t126 + t47) * t576 + t590;];
tau = t1(:);
