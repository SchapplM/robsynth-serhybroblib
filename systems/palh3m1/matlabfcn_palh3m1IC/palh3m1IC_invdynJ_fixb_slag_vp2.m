% Calculate vector of inverse dynamics joint torques with ic for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh3m1IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [10 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: qJDD has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1IC_invdynJ_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:31:02
% EndTime: 2020-04-20 17:31:35
% DurationCPUTime: 31.13s
% Computational Cost: add. (11618->800), mult. (26339->1135), div. (16->6), fcn. (20933->36), ass. (0->393)
t473 = qJ(3) + qJ(4);
t330 = qJ(2) + t473;
t304 = sin(t330);
t305 = cos(t330);
t350 = cos(qJ(5));
t527 = mrSges(6,1) * t350;
t628 = -mrSges(5,2) * t305 + (-pkin(8) * m(6) - mrSges(5,1) - t527) * t304;
t622 = m(4) + m(8);
t354 = cos(qJ(1));
t589 = t628 * t354;
t338 = qJ(2) + qJ(3);
t324 = sin(t338);
t326 = cos(t338);
t614 = mrSges(4,1) * t324 + mrSges(4,2) * t326;
t625 = -t614 * t354 + t589;
t347 = sin(qJ(1));
t588 = t628 * t347;
t624 = -t614 * t347 + t588;
t344 = sin(qJ(4));
t345 = sin(qJ(3));
t531 = pkin(1) * qJD(2);
t452 = t345 * t531;
t281 = t344 * t452;
t351 = cos(qJ(4));
t335 = qJD(2) + qJD(3);
t352 = cos(qJ(3));
t450 = t352 * t531;
t391 = -pkin(4) * t335 + t450;
t186 = -t351 * t391 + t281;
t552 = pkin(1) * t345;
t449 = qJD(3) * t552;
t490 = pkin(1) * qJDD(2);
t239 = qJD(2) * t449 - t352 * t490;
t333 = qJDD(2) + qJDD(3);
t200 = pkin(4) * t333 + t239;
t464 = qJD(3) * t352;
t240 = (-qJD(2) * t464 - qJDD(2) * t345) * pkin(1);
t101 = qJD(4) * t186 + t344 * t200 + t351 * t240;
t393 = t344 * t391 + t351 * t452;
t102 = qJD(4) * t393 + t351 * t200 - t344 * t240;
t321 = qJD(4) + t335;
t343 = sin(qJ(5));
t346 = sin(qJ(2));
t353 = cos(qJ(2));
t246 = t345 * t346 - t352 * t353;
t226 = t246 * qJD(1);
t394 = t345 * t353 + t346 * t352;
t228 = t394 * qJD(1);
t398 = t226 * t344 - t351 * t228;
t127 = t321 * t350 - t343 * t398;
t128 = t321 * t343 + t350 * t398;
t423 = t351 * t226 + t228 * t344;
t156 = Ifges(5,4) * t423;
t157 = qJD(5) - t423;
t331 = t353 * pkin(1);
t303 = t331 + pkin(12);
t273 = t303 * qJD(1);
t185 = -pkin(4) * t226 - t273;
t319 = qJDD(4) + t333;
t181 = -t321 * pkin(8) - t186;
t413 = mrSges(6,1) * t343 + mrSges(6,2) * t350;
t382 = t181 * t413;
t458 = qJD(1) * qJD(2);
t257 = qJDD(1) * t353 - t346 * t458;
t258 = qJDD(1) * t346 + t353 * t458;
t376 = t246 * qJD(3);
t136 = qJD(1) * t376 - t257 * t345 - t258 * t352;
t377 = t394 * qJD(3);
t137 = qJD(1) * t377 - t257 * t352 + t258 * t345;
t51 = qJD(4) * t423 + t136 * t351 + t137 * t344;
t52 = -qJD(4) * t398 - t136 * t344 + t137 * t351;
t459 = qJDD(1) * pkin(12);
t209 = -pkin(1) * t257 - t459;
t99 = -pkin(4) * t137 + t209;
t13 = -pkin(8) * t52 - pkin(10) * t51 + t99;
t182 = pkin(10) * t321 - t393;
t73 = -pkin(8) * t423 - pkin(10) * t398 + t185;
t53 = -t182 * t343 + t350 * t73;
t95 = pkin(10) * t319 + t101;
t4 = qJD(5) * t53 + t13 * t343 + t350 * t95;
t404 = Ifges(6,5) * t350 - Ifges(6,6) * t343;
t514 = Ifges(6,4) * t350;
t407 = -Ifges(6,2) * t343 + t514;
t515 = Ifges(6,4) * t343;
t409 = Ifges(6,1) * t350 - t515;
t461 = qJD(5) * t343;
t434 = -t461 / 0.2e1;
t126 = Ifges(6,4) * t127;
t47 = t128 * Ifges(6,1) + t157 * Ifges(6,5) + t126;
t493 = t350 * t47;
t439 = t493 / 0.2e1;
t45 = t128 * Ifges(6,5) + t127 * Ifges(6,6) + t157 * Ifges(6,3);
t508 = t128 * Ifges(6,4);
t46 = t127 * Ifges(6,2) + t157 * Ifges(6,6) + t508;
t460 = qJD(5) * t350;
t54 = t182 * t350 + t343 * t73;
t5 = -qJD(5) * t54 + t13 * t350 - t343 * t95;
t505 = t393 * mrSges(5,3);
t506 = t186 * mrSges(5,3);
t507 = t398 * Ifges(5,4);
t519 = mrSges(6,3) * t350;
t520 = mrSges(6,3) * t343;
t557 = t343 / 0.2e1;
t559 = -t321 / 0.2e1;
t565 = t398 / 0.2e1;
t566 = -t398 / 0.2e1;
t567 = -t423 / 0.2e1;
t568 = -t157 / 0.2e1;
t570 = -t128 / 0.2e1;
t571 = -t127 / 0.2e1;
t50 = qJDD(5) - t52;
t575 = t50 / 0.2e1;
t20 = -qJD(5) * t128 + t319 * t350 - t343 * t51;
t577 = t20 / 0.2e1;
t19 = qJD(5) * t127 + t319 * t343 + t350 * t51;
t578 = t19 / 0.2e1;
t596 = t54 * mrSges(6,2);
t597 = t53 * mrSges(6,1);
t523 = mrSges(6,2) * t343;
t615 = t523 - t527;
t8 = t19 * Ifges(6,4) + t20 * Ifges(6,2) + t50 * Ifges(6,6);
t85 = Ifges(5,2) * t423 + t321 * Ifges(5,6) + t507;
t86 = Ifges(5,1) * t398 + t321 * Ifges(5,5) + t156;
t9 = t19 * Ifges(6,1) + t20 * Ifges(6,4) + t50 * Ifges(6,5);
t96 = -t319 * pkin(8) - t102;
t623 = -t101 * mrSges(5,2) - t5 * t520 + t102 * mrSges(5,1) + (Ifges(6,1) * t343 + t514) * t578 + (Ifges(6,2) * t350 + t515) * t577 + Ifges(5,3) * t319 + t4 * t519 + (Ifges(6,5) * t343 + Ifges(6,6) * t350) * t575 + t46 * t434 + Ifges(5,6) * t52 + Ifges(5,5) * t51 + t9 * t557 + t350 * t8 / 0.2e1 + t96 * t615 + (-t460 * t53 - t461 * t54) * mrSges(6,3) + (t382 + t439) * qJD(5) + (t127 * t407 + t128 * t409 + t157 * t404) * qJD(5) / 0.2e1 + t85 * t565 + (t156 + t86) * t567 + (-t507 + t45) * t566 + (-t185 * mrSges(5,2) + Ifges(5,1) * t566 + Ifges(5,5) * t559 + t404 * t568 + t407 * t571 + t409 * t570 + t519 * t53 + t520 * t54 - t382 + t46 * t557 - t493 / 0.2e1 + t506) * t423 + (-t185 * mrSges(5,1) + Ifges(6,5) * t570 - Ifges(5,2) * t567 - Ifges(5,6) * t559 + Ifges(6,6) * t571 + Ifges(6,3) * t568 - t505 + t596 - t597) * t398;
t416 = pkin(8) * t305 + pkin(10) * t304;
t543 = pkin(4) * t326;
t621 = m(6) * (-t416 - t543);
t620 = -mrSges(5,1) * t321 - mrSges(6,1) * t127 + mrSges(6,2) * t128 + mrSges(5,3) * t398;
t341 = sin(qJ(7));
t348 = cos(qJ(7));
t245 = -t341 * t346 + t348 * t353;
t225 = t245 * qJD(1);
t247 = t341 * t353 + t346 * t348;
t227 = t247 * qJD(1);
t339 = sin(pkin(15));
t340 = cos(pkin(15));
t133 = -t273 + (-t225 * t340 - t227 * t339) * pkin(3);
t554 = sin(qJ(8));
t555 = cos(qJ(8));
t242 = t339 * t554 + t340 * t555;
t366 = -t339 * t555 + t340 * t554;
t111 = -t225 * t366 - t227 * t242;
t424 = -t242 * t225 + t227 * t366;
t65 = -mrSges(9,1) * t424 + mrSges(9,2) * t111;
t619 = m(9) * t133 + t65;
t92 = -mrSges(5,1) * t423 + mrSges(5,2) * t398;
t618 = m(5) * t185 + t92;
t617 = -mrSges(4,1) * t226 - mrSges(8,1) * t225 - mrSges(4,2) * t228 + mrSges(8,2) * t227;
t337 = qJ(2) + qJ(7);
t323 = sin(t337);
t325 = cos(t337);
t336 = -qJ(7) + pkin(15);
t327 = -qJ(8) + t336;
t306 = -qJ(2) + t327;
t286 = sin(t306);
t287 = cos(t306);
t410 = t287 * mrSges(9,1) + t286 * mrSges(9,2);
t616 = -t325 * mrSges(8,1) + mrSges(8,2) * t323 + t410;
t138 = -mrSges(5,2) * t321 + mrSges(5,3) * t423;
t79 = -mrSges(6,2) * t157 + mrSges(6,3) * t127;
t80 = mrSges(6,1) * t157 - mrSges(6,3) * t128;
t613 = -t343 * t80 + t350 * t79 + t138;
t402 = t54 * t343 + t53 * t350;
t612 = -m(6) * t402 - t343 * t79 - t350 * t80;
t10 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t611 = m(6) * t96 + t10;
t544 = pkin(4) * t324;
t551 = pkin(1) * t346;
t261 = t544 - t551;
t445 = t304 * t523;
t390 = -t305 * mrSges(6,3) - t445;
t526 = mrSges(9,1) * t286;
t539 = pkin(10) * t305;
t329 = -qJ(2) + t336;
t548 = pkin(3) * sin(t329);
t610 = -m(6) * (t261 - t539) - t390 - m(9) * (t548 - t551) + t526 + m(4) * t551 - m(5) * t261;
t609 = m(6) * t181 + t620;
t272 = -mrSges(3,1) * t353 + mrSges(3,2) * t346;
t427 = -t305 * mrSges(5,1) + t304 * mrSges(5,2);
t428 = -mrSges(4,1) * t326 + t324 * mrSges(4,2);
t608 = t272 - t427 - t428 + t616;
t11 = mrSges(6,1) * t50 - mrSges(6,3) * t19;
t12 = -mrSges(6,2) * t50 + mrSges(6,3) * t20;
t581 = m(6) * (-qJD(5) * t402 - t343 * t5 + t350 * t4) - t80 * t460 - t79 * t461 + t350 * t12 - t343 * t11;
t334 = qJD(2) + qJD(7);
t320 = qJD(8) + t334;
t560 = -t320 / 0.2e1;
t572 = -t111 / 0.2e1;
t607 = -t133 * mrSges(9,2) + Ifges(9,1) * t572 + Ifges(9,5) * t560;
t328 = pkin(16) + t337;
t298 = sin(t328);
t301 = cos(t328);
t342 = sin(qJ(6));
t349 = cos(qJ(6));
t188 = 0.1e1 / (-t298 * t349 + t301 * t342) / pkin(5) / pkin(2);
t252 = -pkin(2) * t298 - t551;
t253 = pkin(2) * t301 + t331;
t606 = pkin(2) * t188 * (-t252 * t301 - t253 * t298);
t106 = Ifges(9,4) * t424;
t378 = t245 * qJD(7);
t134 = qJD(1) * t378 + t257 * t341 + t258 * t348;
t379 = t247 * qJD(7);
t135 = -qJD(1) * t379 + t257 * t348 - t258 * t341;
t219 = t242 * qJD(8);
t220 = t366 * qJD(8);
t29 = -t134 * t242 - t135 * t366 - t219 * t225 + t220 * t227;
t30 = t134 * t366 - t135 * t242 + t219 * t227 + t220 * t225;
t332 = qJDD(2) + qJDD(7);
t318 = qJDD(8) + t332;
t453 = t341 * t531;
t547 = pkin(3) * t339;
t234 = t334 * t547 + t453;
t451 = t348 * t531;
t546 = pkin(3) * t340;
t235 = t334 * t546 + t451;
t117 = t234 * t242 + t235 * t366;
t509 = t117 * mrSges(9,3);
t530 = pkin(1) * qJD(7);
t444 = qJD(2) * t530;
t237 = -t341 * t444 + t348 * t490;
t197 = t332 * t546 + t237;
t238 = t341 * t490 + t348 * t444;
t198 = t332 * t547 + t238;
t56 = -t197 * t366 - t198 * t242 - t219 * t235 + t220 * t234;
t57 = -t197 * t242 + t198 * t366 + t219 * t234 + t220 * t235;
t573 = -t424 / 0.2e1;
t59 = t111 * Ifges(9,4) + Ifges(9,2) * t424 + t320 * Ifges(9,6);
t60 = t111 * Ifges(9,1) + t320 * Ifges(9,5) + t106;
t604 = t57 * mrSges(9,1) - t56 * mrSges(9,2) + Ifges(9,5) * t29 + Ifges(9,6) * t30 + Ifges(9,3) * t318 + (t60 + t106) * t573 - (t133 * mrSges(9,1) + Ifges(9,4) * t572 + Ifges(9,2) * t573 + Ifges(9,6) * t560 + t509 - t59 / 0.2e1) * t111;
t363 = -t445 + (-m(6) * pkin(10) - mrSges(6,3)) * t305;
t502 = t304 * mrSges(6,3);
t374 = t305 * t615 - t502;
t590 = -t374 - t427;
t93 = pkin(8) * t398 - pkin(10) * t423;
t71 = -t186 * t343 + t350 * t93;
t72 = t186 * t350 + t343 * t93;
t602 = pkin(7) * (-m(6) * (t53 * t71 + t54 * t72) - t72 * t79 - t71 * t80 - t186 * t138 + t609 * t393 + (m(6) * t416 + t590) * g(3) + (-t347 * t363 + t588) * g(2) + (-t354 * t363 + t589) * g(1) - t611 * pkin(8) + t581 * pkin(10) + t623);
t116 = t234 * t366 - t235 * t242;
t521 = mrSges(9,2) * t287;
t259 = t347 * t521;
t260 = t354 * t521;
t511 = t424 * mrSges(9,3);
t97 = -mrSges(9,2) * t320 + t511;
t98 = mrSges(9,1) * t320 - mrSges(9,3) * t111;
t600 = pkin(9) * (-t117 * t98 + g(3) * t410 - g(2) * (-t347 * t526 + t259) - g(1) * (-t354 * t526 + t260) + (-t97 + t511) * t116 + t607 * t424 + t604);
t599 = t257 / 0.2e1;
t558 = t342 / 0.2e1;
t595 = t353 * Ifges(3,2);
t594 = (t252 * t349 + t253 * t342) * pkin(5) * t188;
t384 = pkin(1) * (t242 * t341 + t348 * t366);
t593 = -qJD(2) * t384 + (t219 * t339 + t220 * t340) * pkin(3);
t383 = pkin(1) * (-t242 * t348 + t341 * t366);
t592 = -qJD(2) * t383 + (-t219 * t340 + t220 * t339) * pkin(3);
t591 = -m(6) * (-t539 + t544) - t390 - m(5) * t544;
t175 = t246 * t344 - t351 * t394;
t179 = qJD(2) * t246 + t376;
t180 = qJD(2) * t394 + t377;
t397 = t351 * t246 + t344 * t394;
t76 = qJD(4) * t397 + t179 * t351 + t180 * t344;
t388 = t175 * t460 + t343 * t76;
t584 = Ifges(3,6) * qJDD(2);
t271 = -mrSges(7,1) * t349 + mrSges(7,2) * t342;
t373 = m(7) * pkin(6) + t271;
t583 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t582 = -mrSges(9,3) - mrSges(8,3) - mrSges(5,3) - mrSges(4,3) - mrSges(7,3) - mrSges(3,3) + mrSges(2,2);
t433 = t331 - t543;
t254 = pkin(12) + t433;
t285 = pkin(3) * cos(t329);
t468 = t285 + t331;
t580 = m(6) * (-t254 + t416) - mrSges(2,1) - m(9) * (pkin(12) + t468) - m(5) * t254 - m(3) * pkin(12) + t373 + t502 + t608 - t622 * t303;
t569 = t128 / 0.2e1;
t562 = t227 / 0.2e1;
t561 = -t228 / 0.2e1;
t556 = t353 / 0.2e1;
t553 = pkin(1) * t341;
t550 = pkin(1) * t348;
t549 = pkin(1) * t352;
t545 = pkin(4) * t228;
t542 = pkin(4) * t344;
t541 = pkin(4) * t351;
t532 = -qJD(6) / 0.2e1;
t518 = mrSges(8,3) * t227;
t517 = Ifges(3,4) * t346;
t516 = Ifges(3,4) * t353;
t513 = Ifges(7,4) * t342;
t512 = Ifges(7,4) * t349;
t510 = t116 * mrSges(9,3);
t504 = t227 * Ifges(8,4);
t503 = t228 * Ifges(4,4);
t500 = t342 * Ifges(7,1);
t495 = t349 * Ifges(7,2);
t489 = t175 * t343;
t488 = t175 * t350;
t487 = (-mrSges(8,2) * t334 + mrSges(8,3) * t225) * t348;
t406 = t495 + t513;
t481 = t342 * (Ifges(7,6) * qJD(6) + qJD(1) * t406);
t480 = t343 * t354;
t479 = t344 * t345;
t478 = t345 * t351;
t477 = t347 * t343;
t476 = t347 * t350;
t312 = qJD(1) * t512;
t475 = t349 * (Ifges(7,5) * qJD(6) + qJD(1) * t500 + t312);
t474 = t350 * t354;
t467 = qJD(1) * t346;
t466 = qJD(1) * t353;
t465 = qJD(2) * t346;
t463 = qJD(4) * t344;
t462 = qJD(4) * t351;
t457 = qJD(1) * qJD(6);
t454 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t50;
t315 = pkin(1) * t465;
t162 = -pkin(4) * t180 + t315;
t430 = t458 / 0.2e1;
t422 = mrSges(4,3) * t452;
t421 = mrSges(4,3) * t450;
t420 = mrSges(8,3) * t453;
t419 = mrSges(8,3) * t451;
t199 = -pkin(4) * t246 - t303;
t309 = pkin(4) - t549;
t216 = -pkin(1) * t478 + t344 * t309;
t414 = mrSges(3,1) * t346 + mrSges(3,2) * t353;
t412 = mrSges(7,1) * t342 + mrSges(7,2) * t349;
t411 = -mrSges(8,1) * t323 - mrSges(8,2) * t325;
t408 = t517 + t595;
t405 = Ifges(3,5) * t353 - Ifges(3,6) * t346;
t403 = Ifges(7,5) * t349 - Ifges(7,6) * t342;
t401 = -t343 * t53 + t350 * t54;
t395 = t344 * t352 + t478;
t215 = pkin(1) * t479 + t309 * t351;
t387 = t175 * t461 - t350 * t76;
t386 = pkin(6) * t412;
t385 = pkin(12) * t414;
t381 = t342 * (Ifges(7,1) * t349 - t513);
t380 = t346 * (Ifges(3,1) * t353 - t517);
t154 = (-t225 * t339 + t227 * t340) * pkin(3);
t78 = -t545 + t93;
t150 = t225 * Ifges(8,2) + t334 * Ifges(8,6) + t504;
t212 = Ifges(8,4) * t225;
t152 = t227 * Ifges(8,1) + t334 * Ifges(8,5) + t212;
t356 = -t238 * mrSges(8,2) + t273 * (mrSges(8,1) * t227 + mrSges(8,2) * t225) + t150 * t562 - t227 * (Ifges(8,1) * t225 - t504) / 0.2e1 + Ifges(8,6) * t135 + Ifges(8,5) * t134 - t334 * (Ifges(8,5) * t225 - Ifges(8,6) * t227) / 0.2e1 + t225 * t419 + t237 * mrSges(8,1) + Ifges(8,3) * t332 - (-Ifges(8,2) * t227 + t152 + t212) * t225 / 0.2e1 + (t510 + t607) * t424 + t604;
t151 = t226 * Ifges(4,2) + t335 * Ifges(4,6) - t503;
t213 = Ifges(4,4) * t226;
t153 = -t228 * Ifges(4,1) + t335 * Ifges(4,5) + t213;
t355 = Ifges(4,5) * t136 + Ifges(4,6) * t137 + t228 * t422 + t228 * (Ifges(4,1) * t226 + t503) / 0.2e1 + t151 * t561 + t239 * mrSges(4,1) - t240 * mrSges(4,2) + Ifges(4,3) * t333 - t335 * (Ifges(4,5) * t226 + Ifges(4,6) * t228) / 0.2e1 - t226 * t421 + t273 * (-mrSges(4,1) * t228 + mrSges(4,2) * t226) - (Ifges(4,2) * t228 + t153 + t213) * t226 / 0.2e1 + t623;
t322 = pkin(14) + t473;
t313 = Ifges(3,4) * t466;
t308 = -pkin(8) - t541;
t300 = cos(t327);
t297 = sin(t327);
t295 = cos(t322);
t294 = sin(t322);
t275 = t546 + t550;
t274 = t547 + t553;
t256 = qJDD(1) * t342 + t349 * t457;
t255 = qJDD(1) * t349 - t342 * t457;
t250 = pkin(4) * t352 + pkin(9) * t295;
t249 = -pkin(4) * t345 - pkin(9) * t294;
t233 = -t300 * pkin(7) + pkin(3) * cos(t336);
t232 = -t297 * pkin(7) + pkin(3) * sin(t336);
t224 = Ifges(3,1) * t467 + Ifges(3,5) * qJD(2) + t313;
t222 = Ifges(3,6) * qJD(2) + qJD(1) * t408;
t211 = -t351 * t450 + t281;
t210 = t395 * t531;
t206 = t305 * t474 - t477;
t205 = t305 * t480 + t476;
t204 = t305 * t476 + t480;
t203 = t305 * t477 - t474;
t193 = mrSges(4,1) * t335 + mrSges(4,3) * t228;
t192 = mrSges(8,1) * t334 - t518;
t191 = -mrSges(4,2) * t335 + mrSges(4,3) * t226;
t178 = -qJD(2) * t247 - t379;
t177 = qJD(2) * t245 + t378;
t176 = 0.1e1 / (t294 * t300 + t295 * t297) / pkin(9) / pkin(7);
t168 = (-t242 * t339 - t340 * t366) * pkin(3);
t167 = (-t242 * t340 + t339 * t366) * pkin(3);
t141 = -t242 * t274 - t275 * t366;
t140 = -t242 * t275 + t274 * t366;
t82 = qJD(7) * t384 + t219 * t274 + t220 * t275;
t81 = qJD(7) * t383 - t219 * t275 + t220 * t274;
t77 = qJD(4) * t175 + t179 * t344 - t351 * t180;
t67 = t211 * t350 + t343 * t78;
t66 = -t211 * t343 + t350 * t78;
t63 = (-t134 * t339 - t135 * t340) * pkin(3) + t209;
t44 = t177 * t366 - t178 * t242 + t219 * t247 + t220 * t245;
t43 = -t177 * t242 - t178 * t366 - t219 * t245 + t220 * t247;
t39 = -mrSges(5,2) * t319 + mrSges(5,3) * t52;
t38 = mrSges(5,1) * t319 - mrSges(5,3) * t51;
t22 = -mrSges(9,2) * t318 + mrSges(9,3) * t30;
t21 = mrSges(9,1) * t318 - mrSges(9,3) * t29;
t1 = [t622 * (-t209 * t303 - t273 * t315) + (m(5) * t99 - mrSges(5,1) * t52 + mrSges(5,2) * t51) * t199 + (-mrSges(9,1) * t63 + mrSges(9,3) * t56 + Ifges(9,4) * t29 + Ifges(9,2) * t30 + Ifges(9,6) * t318) * (-t242 * t245 + t247 * t366) + (mrSges(9,2) * t63 - mrSges(9,3) * t57 + Ifges(9,1) * t29 + Ifges(9,4) * t30 + Ifges(9,5) * t318) * (-t242 * t247 - t245 * t366) - (mrSges(4,2) * t209 - mrSges(4,3) * t239 + Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * t333) * t394 + (t516 * t430 + t584 / 0.2e1) * t353 + (Ifges(3,4) * t258 + Ifges(3,2) * t257 + t584) * t556 + (-mrSges(4,1) * t209 + mrSges(4,3) * t240 + Ifges(4,4) * t136 + Ifges(4,2) * t137 + Ifges(4,6) * t333) * t246 + t258 * t516 / 0.2e1 + (t206 * mrSges(6,1) - t205 * mrSges(6,2) + t347 * t582 + t354 * t580) * g(2) + (-t204 * mrSges(6,1) + t203 * mrSges(6,2) - t347 * t580 + t354 * t582) * g(1) + t77 * t597 + t408 * t599 + (m(6) * (qJD(5) * t401 + t343 * t4 + t350 * t5) + t79 * t460 + t343 * t12 + t350 * t11 - t80 * t461) * (-pkin(8) * t397 - pkin(10) * t175 + t199) - (t454 / 0.2e1 + Ifges(6,3) * t575 + Ifges(6,6) * t577 + Ifges(6,5) * t578 - t101 * mrSges(5,3) + t99 * mrSges(5,1) - Ifges(5,2) * t52 - Ifges(5,4) * t51 - Ifges(5,6) * t319 + t583) * t397 + (Ifges(3,1) * t258 + Ifges(3,4) * t599 + Ifges(3,5) * qJDD(2) - t430 * t595) * t346 + (0.2e1 * Ifges(7,5) * t558 + Ifges(7,6) * t349) * qJDD(6) + t9 * t488 / 0.2e1 - t612 * (pkin(8) * t77 - pkin(10) * t76 + t162) + (mrSges(8,2) * t209 - mrSges(8,3) * t237 + Ifges(8,1) * t134 + Ifges(8,4) * t135 + Ifges(8,5) * t332) * t247 + (m(3) * pkin(12) ^ 2 + pkin(6) * t373 + Ifges(2,3)) * qJDD(1) + t380 * t430 + (t475 / 0.2e1 + t403 * qJD(6) / 0.2e1) * qJD(6) + t618 * t162 + t619 * (t315 + (-t177 * t339 - t178 * t340) * pkin(3)) + t181 * (mrSges(6,1) * t388 - mrSges(6,2) * t387) - (-mrSges(4,1) * t180 - mrSges(8,1) * t178 + mrSges(4,2) * t179 + mrSges(8,2) * t177) * t273 + t255 * t406 / 0.2e1 + t133 * (-mrSges(9,1) * t44 + mrSges(9,2) * t43) + t111 * (Ifges(9,1) * t43 + Ifges(9,4) * t44) / 0.2e1 - t77 * t85 / 0.2e1 + t76 * t86 / 0.2e1 + t77 * t45 / 0.2e1 + t44 * t59 / 0.2e1 + t43 * t60 / 0.2e1 - t76 * t506 + t617 * t315 + (t99 * mrSges(5,2) - t102 * mrSges(5,3) + Ifges(5,1) * t51 + Ifges(5,4) * t52 + Ifges(5,5) * t319 + t404 * t575 + t407 * t577 + t409 * t578 + t413 * t96 + t434 * t47) * t175 + t423 * (Ifges(5,4) * t76 - Ifges(5,2) * t77) / 0.2e1 - t43 * t510 + (t224 * t556 + t405 * qJD(2) / 0.2e1) * qJD(2) + (-mrSges(8,1) * t209 + mrSges(8,3) * t238 + Ifges(8,4) * t134 + Ifges(8,2) * t135 + Ifges(8,6) * t332) * t245 + (Ifges(7,1) * t256 + Ifges(7,4) * t255) * t558 - t385 * t458 - t272 * t459 + t349 * (Ifges(7,4) * t256 + Ifges(7,2) * t255) / 0.2e1 + (mrSges(4,1) * t137 + mrSges(8,1) * t135 - mrSges(4,2) * t136 - mrSges(8,2) * t134) * t303 - t180 * t422 - t388 * t46 / 0.2e1 + (t387 * t53 - t388 * t54 - t4 * t489 - t488 * t5) * mrSges(6,3) - t177 * t419 - t77 * t596 + t157 * (-Ifges(6,5) * t387 - Ifges(6,6) * t388 + Ifges(6,3) * t77) / 0.2e1 + t127 * (-Ifges(6,4) * t387 - Ifges(6,2) * t388 + Ifges(6,6) * t77) / 0.2e1 + t256 * (t500 + t512) / 0.2e1 - t8 * t489 / 0.2e1 + t178 * t420 + t179 * t421 - t222 * t465 / 0.2e1 + t386 * t457 + t177 * t152 / 0.2e1 + t178 * t150 / 0.2e1 + t179 * t153 / 0.2e1 + t180 * t151 / 0.2e1 + t185 * (mrSges(5,1) * t77 + mrSges(5,2) * t76) + t225 * (Ifges(8,4) * t177 + Ifges(8,2) * t178) / 0.2e1 + t226 * (Ifges(4,4) * t179 + Ifges(4,2) * t180) / 0.2e1 + t77 * t505 - t44 * t509 + t481 * t532 + (Ifges(4,1) * t179 + Ifges(4,4) * t180) * t561 + (Ifges(8,1) * t177 + Ifges(8,4) * t178) * t562 + (Ifges(5,1) * t76 - Ifges(5,4) * t77) * t565 + (-Ifges(6,1) * t387 - Ifges(6,4) * t388 + Ifges(6,5) * t77) * t569 + t76 * t439 + pkin(6) * (-mrSges(7,1) * t255 + mrSges(7,2) * t256) - pkin(12) * (-mrSges(3,1) * t257 + mrSges(3,2) * t258) + t320 * (Ifges(9,5) * t43 + Ifges(9,6) * t44) / 0.2e1 + t321 * (Ifges(5,5) * t76 - Ifges(5,6) * t77) / 0.2e1 + t334 * (Ifges(8,5) * t177 + Ifges(8,6) * t178) / 0.2e1 + t335 * (Ifges(4,5) * t179 + Ifges(4,6) * t180) / 0.2e1 + (m(9) * t63 - mrSges(9,1) * t30 + mrSges(9,2) * t29) * ((-t245 * t340 - t247 * t339) * pkin(3) - t303) + (t381 + t349 * (-Ifges(7,2) * t342 + t512)) * t457 / 0.2e1 + t424 * (Ifges(9,4) * t43 + Ifges(9,2) * t44) / 0.2e1; (-t192 * t341 + t487) * t530 + t581 * (pkin(10) + t216) + (m(8) * (t237 * t348 + t238 * t341) + m(4) * (-t239 * t352 - t240 * t345) - t191 * t464) * pkin(1) + (t622 * t273 + t612 - t617 - t618 - t619) * pkin(1) * t467 + t356 + t355 + (t414 + m(8) * t551 - t411 - (-m(9) * t548 - t411 + t526) * t594 + t412 * t606) * (g(1) * t354 + g(2) * t347) - (-Ifges(3,2) * t467 + t224 + t313) * t466 / 0.2e1 + (-t380 / 0.2e1 + t385) * qJD(1) ^ 2 + m(9) * (t116 * t82 - t117 * t81 + t140 * t57 + t141 * t56) + (t271 * t606 - t621 - t374 - m(5) * t433 - m(9) * t468 + (-m(6) - t622) * t331 + t608) * g(3) + t611 * (-pkin(8) - t215) + (-m(5) * t393 + m(6) * t401 + t613) * (t309 * t462 + (t345 * t463 + (-t351 * t352 + t479) * qJD(3)) * pkin(1)) - (mrSges(4,1) * t333 - mrSges(4,3) * t136) * t549 - (-mrSges(4,2) * t333 + mrSges(4,3) * t137) * t552 + t193 * t449 + t140 * t21 + t141 * t22 + t82 * t98 + t81 * t97 + ((m(9) * t285 - t616) * g(3) + ((-t232 * t300 + t233 * t297) * t602 + (-t232 * t295 - t233 * t294) * t600) * t176 + g(2) * t259 + g(1) * t260 - t356 - (-t487 + (t192 + t518) * t341) * t531 - t592 * t97 - t167 * t21 + t154 * t65 - t593 * t98 - t168 * t22 - (t116 * t593 - t117 * t592 - t133 * t154 + t167 * t57 + t168 * t56) * m(9)) * t594 - t405 * t458 / 0.2e1 + m(5) * (t101 * t216 + t102 * t215) + t612 * t78 + t618 * t545 - t619 * t154 + t227 * t420 + (m(5) * t186 - t609) * (-t309 * t463 + (qJD(3) * t395 + t345 * t462) * pkin(1)) + (Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * qJDD(6) + (-t475 / 0.2e1 + t481 / 0.2e1 - t349 * t312 / 0.2e1 + t403 * t532 + (-t381 / 0.2e1 - t386 + t495 * t558) * qJD(1)) * qJD(1)) * t606 + t222 * t467 / 0.2e1 + t215 * t38 + t216 * t39 + (mrSges(8,1) * t332 - mrSges(8,3) * t134) * t550 + (-mrSges(8,2) * t332 + mrSges(8,3) * t135) * t553 + (t347 * t610 - t259 + t624) * g(2) + (t354 * t610 - t260 + t625) * g(1) + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t257 + Ifges(3,5) * t258; t308 * t10 - t211 * t138 + t191 * t450 - t193 * t452 + t38 * t541 + t39 * t542 + t92 * t545 - t66 * t80 - t67 * t79 + t355 + ((t249 * t295 + t250 * t294) * t600 + (t249 * t300 - t250 * t297) * t602) * t176 + (t181 * t210 - t53 * t66 - t54 * t67 + t308 * t96 + (t181 * t344 + t351 * t401) * qJD(4) * pkin(4)) * m(6) + ((t101 * t344 + t102 * t351 + (-t186 * t344 - t351 * t393) * qJD(4)) * pkin(4) + t185 * t545 - t186 * t210 + t393 * t211) * m(5) + t613 * pkin(4) * t462 + (m(5) * t543 - t428 + t590 - t621) * g(3) + (t347 * t591 + t624) * g(2) + (t354 * t591 + t625) * g(1) + t581 * (pkin(10) + t542) + t620 * (pkin(4) * t463 + t210); -t181 * (mrSges(6,1) * t128 + mrSges(6,2) * t127) + (Ifges(6,5) * t127 - Ifges(6,6) * t128) * t568 + t54 * t80 - t53 * t79 + t46 * t569 + (Ifges(6,1) * t127 - t508) * t570 - g(1) * (mrSges(6,1) * t205 + mrSges(6,2) * t206) - g(2) * (mrSges(6,1) * t203 + mrSges(6,2) * t204) - g(3) * t413 * t304 + (t127 * t53 + t128 * t54) * mrSges(6,3) + t454 + (-Ifges(6,2) * t128 + t126 + t47) * t571 + t583;];
tau = t1(:);
