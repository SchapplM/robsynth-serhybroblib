% Calculate vector of inverse dynamics joint torques for
% palh1m2TE
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
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m2TE_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:46:41
% EndTime: 2020-05-01 20:47:04
% DurationCPUTime: 23.73s
% Computational Cost: add. (2218->693), mult. (2970->920), div. (0->0), fcn. (315->136), ass. (0->453)
t646 = 2 * qJD(2) + qJD(3);
t340 = (mrSges(4,1) + mrSges(11,1));
t359 = m(5) + m(6);
t526 = (t359 * pkin(5));
t166 = t340 + t526;
t369 = qJD(1) ^ 2;
t549 = -t369 / 0.2e1;
t645 = t646 * qJD(3);
t368 = qJD(2) ^ 2;
t235 = (2 * t368) + t369;
t480 = (pkin(15) * qJDD(1));
t300 = (pkin(18) - pkin(22));
t231 = (-pkin(21) + t300);
t174 = (-pkin(20) + t231);
t141 = (-qJ(1) + t174);
t111 = qJ(4) + t141;
t140 = (qJ(1) + t174);
t114 = -qJ(4) + t140;
t355 = g(2) * mrSges(6,2);
t556 = g(1) * mrSges(6,1);
t199 = -t355 + t556;
t554 = g(2) * mrSges(6,1);
t555 = g(1) * mrSges(6,2);
t200 = t554 + t555;
t327 = -qJ(4) + qJ(1);
t644 = -(cos(t111) + cos(t114)) * t200 / 0.4e1 - t199 * sin(t327) / 0.2e1 + t200 * cos(t327) / 0.2e1;
t339 = (mrSges(4,2) + mrSges(11,2));
t622 = 16 * t339;
t643 = g(1) * t622;
t347 = (mrSges(3,2) + mrSges(10,2));
t624 = 8 * t347;
t642 = g(1) * t624;
t641 = g(2) * t622;
t640 = g(2) * t624;
t299 = (qJ(3) + pkin(19));
t371 = 2 * qJ(2);
t239 = t371 + t299;
t189 = cos(t239);
t638 = mrSges(10,1) * t189;
t185 = sin(t239);
t637 = mrSges(10,2) * t185;
t181 = qJ(2) + t231;
t151 = qJ(3) + t181;
t182 = -qJ(2) + t231;
t152 = -qJ(3) + t182;
t636 = cos(t152) + cos(t151);
t635 = cos(t182) + cos(t181);
t225 = m(11) + m(4) + m(8) + t359;
t173 = (t225 * pkin(1) ^ 2);
t634 = t173 + Ifges(3,2) + Ifges(10,2);
t366 = (qJD(4) ^ 2);
t633 = -t366 + t645;
t290 = t359 * pkin(5) ^ 2;
t509 = Ifges(11,2) + Ifges(4,2);
t510 = Ifges(11,1) + Ifges(4,1);
t615 = -t509 + t510 - t290;
t632 = 8 * mrSges(11,2) + 8 * mrSges(4,2);
t631 = (8 * mrSges(3,1)) + 0.8e1 * mrSges(10,1);
t527 = g(1) * t526;
t625 = 8 * t340;
t630 = -g(1) * t625 - 8 * t527;
t528 = g(2) * t526;
t629 = -g(2) * t625 - 8 * t528;
t619 = (pkin(1) * t225);
t546 = g(1) * t619;
t348 = mrSges(3,1) + mrSges(10,1);
t621 = 0.16e2 * t348;
t628 = -g(1) * t621 - (16 * t546);
t547 = g(2) * t619;
t627 = -g(2) * t621 - (16 * t547);
t322 = mrSges(6,2) * qJDD(1);
t245 = t322 / 0.2e1;
t472 = qJD(1) * qJD(4);
t270 = mrSges(6,1) * t472;
t104 = t245 + t270;
t120 = sin(t151);
t121 = sin(t152);
t240 = qJ(2) + t299;
t157 = 2 * t240;
t124 = sin(t157);
t125 = cos(t157);
t328 = qJ(3) + qJ(2);
t214 = 2 * t328;
t158 = sin(t214);
t159 = cos(t214);
t186 = sin(t240);
t190 = cos(t240);
t277 = qJ(1) + t328;
t207 = sin(t277);
t278 = qJ(1) - t328;
t208 = sin(t278);
t212 = cos(t277);
t213 = cos(t278);
t444 = m(6) * pkin(11) + mrSges(6,3);
t223 = -mrSges(5,2) + t444;
t361 = m(10) * pkin(2) ^ 2;
t224 = -Ifges(9,1) + Ifges(9,2) + t361;
t259 = cos(t328);
t272 = pkin(2) * mrSges(10,3) - Ifges(9,5);
t281 = pkin(9) * m(6) + mrSges(5,1);
t282 = (pkin(2) * m(10) + mrSges(9,1));
t336 = Ifges(4,6) + Ifges(11,6);
t338 = Ifges(4,4) + Ifges(11,4);
t441 = -mrSges(6,1) * t369 / 0.4e1;
t410 = pkin(5) * t441;
t425 = t281 * t549;
t426 = t223 * t549;
t515 = mrSges(6,2) * t369;
t439 = t515 / 0.4e1;
t533 = pkin(4) * t369;
t451 = -t533 / 0.2e1;
t491 = t632 * g(1);
t493 = t632 * g(2);
t502 = t339 * t369;
t512 = pkin(15) * t369;
t440 = -t515 / 0.4e1;
t142 = qJ(2) + t174;
t108 = qJ(4) + t142;
t82 = qJ(3) + t108;
t52 = sin(t82);
t573 = pkin(5) * cos(t82);
t520 = t52 * t410 + t440 * t573;
t112 = -qJ(4) + t142;
t83 = qJ(3) + t112;
t53 = sin(t83);
t501 = qJ(4) - qJ(2);
t275 = qJ(3) - t501;
t205 = sin(t275);
t531 = pkin(5) * t205;
t324 = qJ(4) + qJ(2);
t274 = qJ(3) + t324;
t204 = sin(t274);
t532 = pkin(5) * t204;
t548 = t369 / 0.2e1;
t143 = -qJ(2) + t174;
t117 = -qJ(3) + t143;
t568 = pkin(5) * cos(t117);
t116 = qJ(3) + t142;
t569 = pkin(5) * sin(t116);
t57 = cos(t83);
t572 = pkin(5) * t57;
t203 = -qJ(1) + t240;
t356 = mrSges(9,2) * g(2);
t393 = t282 * g(1);
t558 = mrSges(9,2) * g(1);
t618 = g(2) * t282;
t613 = (t558 + t618) * cos(t203) / 0.2e1 + (-t356 + t393) * sin(t203) / 0.2e1;
t202 = qJ(1) + t240;
t620 = (t558 - t618) * cos(t202) / 0.2e1 + (t356 + t393) * sin(t202) / 0.2e1;
t68 = sin(t117);
t75 = cos(t116);
t626 = t620 + t53 * t410 + (pkin(15) * t502 - g(3) * t166 + qJDD(1) * t336) * t259 + t613 + t520 + (-t531 - t532) * t104 + (-Ifges(9,4) * t125 - t338 * t159) * t369 + t426 * t568 + t425 * t569 + t439 * t572 + (g(3) * mrSges(9,2) - qJDD(1) * t272 + t282 * t512) * t186 + (mrSges(9,2) * t512 + Ifges(9,6) * qJDD(1) - g(3) * t282) * t190 + t636 * mrSges(11,2) * t451 + (t120 * t451 + t121 * t533 / 0.2e1) * mrSges(11,1) + (t224 * t124 + (t223 * t75 + t281 * t68) * pkin(5)) * t548 + (t493 + t630) * t208 / 0.16e2 + (t491 + t629) * t212 / 0.16e2 + (t493 - t630) * t207 / 0.16e2 + (t491 - t629) * t213 / 0.16e2 + t615 * t158 * t549;
t380 = 32 * t368 + 32 * t645;
t413 = pkin(1) * t441;
t581 = pkin(1) * sin(t108);
t69 = cos(t108);
t617 = t69 * t413 + t439 * t581;
t523 = Ifges(3,1) + Ifges(10,1);
t416 = -t523 + t634;
t496 = t290 + t361;
t298 = qJD(3) + qJD(2);
t397 = qJD(1) * t298;
t283 = qJD(2) + qJD(3) / 0.2e1;
t391 = t283 * qJD(1);
t109 = qJ(4) + t143;
t84 = -qJ(3) + t109;
t54 = sin(t84);
t113 = -qJ(4) + t143;
t85 = -qJ(3) + t113;
t55 = sin(t85);
t614 = pkin(5) * (t54 + t55);
t612 = 16 * m(6);
t608 = -16 * g(1);
t607 = 16 * g(1);
t606 = -16 * g(2);
t605 = 16 * g(2);
t604 = 0.16e2 * mrSges(6,1);
t602 = -0.16e2 * mrSges(6,2);
t601 = 0.64e2 * pkin(15);
t600 = -16 * Ifges(6,5);
t599 = 16 * Ifges(6,6);
t230 = (mrSges(4,3) + mrSges(5,3) + mrSges(8,3) + mrSges(11,3));
t598 = 32 * mrSges(2,2) - 32 * mrSges(3,3) - 32 * mrSges(7,3) - 32 * mrSges(9,3) - 32 * mrSges(10,3) - 32 * t230;
t335 = pkin(11) * qJDD(1);
t443 = (pkin(9) * t472);
t597 = 0.16e2 * t335 - (16 * t443);
t332 = qJDD(1) * pkin(9);
t442 = pkin(11) * t472;
t596 = 0.16e2 * t332 - 0.16e2 * t442;
t507 = qJDD(4) * pkin(11);
t525 = t366 * pkin(9);
t595 = 0.16e2 * t507 - 0.16e2 * t525;
t508 = pkin(9) * qJDD(4);
t524 = t366 * pkin(11);
t594 = 0.16e2 * t508 - 0.16e2 * t524;
t293 = (qJDD(2) + qJDD(3));
t593 = -32 * t293;
t592 = 32 * t293;
t591 = 16 * t366;
t590 = -32 * t368;
t588 = pkin(11) * g(1);
t587 = g(1) * pkin(9);
t586 = mrSges(6,1) / 0.4e1;
t585 = 0.16e2 * qJDD(1);
t584 = -32 * qJDD(2);
t62 = sin(t109);
t580 = pkin(1) * t62;
t64 = sin(t112);
t579 = pkin(1) * t64;
t65 = sin(t113);
t578 = pkin(1) * t65;
t58 = cos(t84);
t571 = pkin(5) * t58;
t59 = cos(t85);
t570 = pkin(5) * t59;
t567 = m(11) * pkin(4);
t566 = mrSges(6,1) * pkin(9);
t565 = mrSges(7,1) * g(1);
t564 = mrSges(7,1) * g(2);
t563 = mrSges(8,1) * g(1);
t562 = mrSges(6,2) * pkin(9);
t561 = mrSges(7,2) * g(1);
t560 = mrSges(8,2) * g(1);
t559 = mrSges(8,2) * g(2);
t557 = g(1) * mrSges(5,1);
t552 = t199 / 0.4e1;
t249 = sin(t324);
t545 = pkin(1) * t249;
t250 = -sin(t501);
t544 = pkin(1) * t250;
t329 = t371 + qJ(3);
t543 = pkin(1) * sin(t329);
t256 = cos(t324);
t542 = pkin(1) * t256;
t257 = cos(t501);
t541 = pkin(1) * t257;
t540 = pkin(1) * cos(t329);
t539 = pkin(1) * sin(qJ(3));
t538 = pkin(1) * cos(qJ(3));
t536 = pkin(2) * sin(t299);
t535 = pkin(2) * cos(t299);
t534 = pkin(2) * t369;
t209 = cos(t274);
t530 = pkin(5) * t209;
t210 = cos(t275);
t529 = pkin(5) * t210;
t521 = Ifges(5,2) + Ifges(6,3);
t519 = mrSges(6,1) * t293;
t518 = mrSges(6,1) * t366;
t517 = mrSges(6,2) * t293;
t516 = mrSges(6,2) * t366;
t514 = pkin(14) * t369;
t129 = t270 + t322;
t513 = pkin(15) * t129;
t503 = t166 * t369;
t474 = qJD(1) * qJD(2);
t150 = t281 * t474;
t473 = qJD(1) * qJD(3);
t500 = (mrSges(5,1) * t473) + t150;
t229 = qJDD(1) * t340;
t468 = qJDD(1) * t359;
t499 = pkin(5) * t468 + t229;
t287 = qJDD(1) * t600;
t498 = t472 * t599 + t287;
t497 = (Ifges(6,6) * t591) + qJDD(4) * t600;
t495 = mrSges(5,2) * t606 + mrSges(6,3) * t605;
t494 = t347 * t605;
t492 = t347 * t608;
t269 = mrSges(6,2) * t472;
t323 = mrSges(6,1) * qJDD(1);
t490 = t323 - t269;
t489 = m(11) * qJDD(1);
t488 = mrSges(5,1) * qJDD(1);
t487 = mrSges(6,1) * qJDD(2);
t486 = mrSges(6,1) * qJDD(4);
t485 = mrSges(8,1) * qJDD(1);
t484 = mrSges(10,1) * qJDD(1);
t483 = mrSges(6,2) * qJDD(2);
t482 = mrSges(6,2) * qJDD(4);
t481 = mrSges(8,2) * qJDD(1);
t479 = Ifges(6,6) * qJDD(1);
t478 = qJD(1) * (qJD(3) + qJD(4));
t477 = qJD(1) * (qJD(3) - qJD(4));
t476 = qJD(3) * t283;
t475 = qJDD(1) * mrSges(10,2);
t471 = qJDD(1) * t223;
t470 = qJDD(1) * t281;
t469 = qJDD(1) * t339;
t345 = Ifges(6,1) - Ifges(6,2);
t467 = t345 * qJDD(1);
t346 = mrSges(5,2) - mrSges(6,3);
t466 = t346 * qJDD(1);
t465 = 0.2e1 * pkin(14);
t464 = 0.16e2 * pkin(9) * t369;
t463 = 0.16e2 * t512;
t292 = qJDD(2) - qJDD(4);
t291 = qJDD(2) + qJDD(4);
t458 = -t571 / 0.2e1;
t457 = t570 / 0.2e1;
t455 = t567 / 0.2e1;
t453 = (Ifges(4,3) + Ifges(9,3) + Ifges(11,3));
t452 = pkin(4) * t489;
t450 = -t529 / 0.2e1;
t448 = -64 * t480;
t446 = t631 * g(2) + (8 * t547);
t445 = t631 * g(1) + (8 * t546);
t438 = -64 * t474;
t276 = qJ(2) + pkin(17) - pkin(18);
t437 = t340 * t607 + 16 * t527;
t334 = t368 * mrSges(6,1);
t436 = t633 * mrSges(6,1) + t334;
t333 = t368 * mrSges(6,2);
t435 = t633 * mrSges(6,2) + t333;
t434 = t340 * t606 - 16 * t528;
t433 = -pkin(11) * mrSges(6,1) + Ifges(6,5);
t432 = -pkin(11) * mrSges(6,2) + Ifges(6,6);
t146 = qJ(4) + t174;
t97 = cos(t146);
t430 = t97 / 0.32e2;
t147 = -qJ(4) + t174;
t98 = cos(t147);
t429 = t98 / 0.32e2;
t428 = Ifges(6,4) * t472;
t427 = Ifges(6,5) * t472;
t424 = 0.2e1 / 0.3e1 * t369 + 0.4e1 / 0.3e1 * t366;
t422 = t334 - t518;
t421 = t333 - t516;
t139 = m(9) + m(10) + m(3) + t225;
t420 = -pkin(14) * m(7) + pkin(15) * t139 + mrSges(2,1);
t418 = t534 * t638;
t417 = t534 * t637;
t415 = t503 * t540;
t414 = t502 * t543;
t396 = 2 * t174;
t389 = mrSges(5,2) * t607 + mrSges(6,3) * t608;
t126 = t223 * t474;
t388 = -(t346 * t473) + t126;
t110 = qJ(4) + t140;
t115 = -qJ(4) + t141;
t198 = t355 + t556;
t201 = t554 - t555;
t325 = qJ(4) + qJ(1);
t386 = t201 * cos(t325) / 0.2e1 - t198 * sin(t325) / 0.2e1 + (sin(t115) - sin(t110)) * t198 / 0.4e1 + (cos(t110) + cos(t115)) * t201 / 0.4e1;
t385 = (Ifges(6,5) * t591) + qJDD(4) * t599;
t384 = -Ifges(5,1) - Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - 0.2e1 * mrSges(6,3) * pkin(11);
t217 = qJ(1) + t276;
t357 = mrSges(7,2) * g(2);
t383 = (t357 + t565) * cos(t217) / 0.2e1 + (-t561 + t564) * sin(t217) / 0.2e1;
t218 = -qJ(1) + t276;
t382 = (-t357 + t565) * cos(t218) / 0.2e1 - (t561 + t564) * sin(t218) / 0.2e1;
t376 = pkin(4) ^ 2;
t374 = pkin(9) ^ 2;
t373 = pkin(11) ^ 2;
t370 = 2 * qJ(4);
t362 = pkin(11) * g(2);
t360 = g(2) * pkin(9);
t358 = mrSges(8,1) * g(2);
t354 = cos(qJ(2));
t352 = cos(qJ(4));
t351 = sin(qJ(2));
t349 = sin(qJ(4));
t344 = Ifges(7,1) - Ifges(7,2);
t343 = (Ifges(3,4) + Ifges(10,4));
t342 = (Ifges(3,5) + Ifges(10,5));
t341 = (Ifges(3,6) + Ifges(10,6));
t337 = (Ifges(4,5) + Ifges(11,5));
t331 = -qJ(2) + qJ(1);
t330 = qJ(2) + qJ(1);
t315 = cos(t371);
t314 = cos(t370);
t313 = sin(t371);
t312 = sin(t370);
t304 = mrSges(5,1) * t605;
t288 = t298 ^ 2;
t285 = -0.16e2 * t479;
t284 = 0.16e2 * t479;
t280 = pkin(9) * t473;
t279 = pkin(11) * t473;
t268 = mrSges(6,1) * t474;
t267 = mrSges(8,1) * t474;
t266 = mrSges(6,2) * t474;
t265 = mrSges(8,2) * t474;
t263 = qJDD(2) + qJDD(3) / 0.2e1;
t262 = cos(t331);
t261 = cos(t330);
t255 = sin(t331);
t254 = sin(t330);
t252 = sin(t328);
t248 = -t323 / 0.2e1;
t247 = t323 / 0.2e1;
t246 = -t322 / 0.2e1;
t244 = -qJ(1) + t300;
t243 = qJ(1) + t300;
t242 = -qJ(2) + t300;
t241 = qJ(2) + t300;
t238 = -g(3) + t332;
t237 = g(3) + t332;
t221 = qJDD(3) + t292;
t220 = qJDD(3) + t291;
t211 = cos(t276);
t206 = sin(t276);
t193 = 2 * t300;
t192 = cos(t242);
t191 = cos(t241);
t188 = sin(t242);
t187 = sin(t241);
t184 = -qJ(1) + t231;
t183 = qJ(1) + t231;
t180 = mrSges(6,2) * t235;
t178 = t235 * mrSges(6,1);
t176 = mrSges(6,2) * t397;
t175 = mrSges(6,1) * t397;
t171 = t508 + t524;
t169 = t507 + t525;
t168 = 2 * t276;
t167 = (2 * t288) + t369;
t165 = t335 - t512;
t164 = t335 + t512;
t145 = -qJ(4) + t396;
t144 = qJ(4) + t396;
t138 = cos(t174);
t136 = t332 + t442;
t134 = t335 + t443;
t133 = cos(t168);
t132 = sin(t168);
t128 = t167 * mrSges(6,1);
t127 = t167 * mrSges(6,2);
t107 = 2 * t174;
t96 = cos(t145);
t95 = cos(t144);
t94 = cos(t143);
t93 = cos(t142);
t92 = sin(t147);
t91 = sin(t146);
t90 = sin(t145);
t89 = sin(t144);
t88 = sin(t143);
t87 = sin(t142);
t79 = 2 * t147;
t78 = 2 * t146;
t73 = cos(t113);
t72 = cos(t112);
t70 = cos(t109);
t66 = sin(t114);
t63 = sin(t111);
t51 = cos(t79);
t50 = cos(t78);
t49 = sin(t79);
t48 = sin(t78);
t47 = -0.32e2 * t490 * pkin(15);
t46 = 0.4e1 / 0.3e1 * t482 + t424 * mrSges(6,1);
t45 = -0.4e1 / 0.3e1 * t486 + t424 * mrSges(6,2);
t1 = [((-t486 + t516) * t352 + (t482 + t518) * t349 - 0.2e1 * cos(t231) * t452 + t139 * t480) * pkin(15) + ((Ifges(7,4) * t438) + t344 * t585) * t133 / 0.32e2 + ((pkin(15) * t225 * t438 + (t230 * t584)) * pkin(1) + (-t347 * qJDD(1) - t348 * t474) * t601 + (t341 * t590) + (32 * qJDD(2) * t342)) * t354 / 0.32e2 + ((t225 * t448 + (32 * t368 * t230)) * pkin(1) + (-t348 * qJDD(1) + (t347 * t474)) * t601 + (t342 * t590) + (t341 * t584)) * t351 / 0.32e2 + (mrSges(6,1) * t595 + t171 * t602 + t497 + 0.32e2 * t513) * t91 / 0.32e2 + (mrSges(6,1) * t597 + t136 * t602 + t498) * t89 / 0.32e2 + (mrSges(6,1) * t596 + t134 * t602 + t284 + 0.16e2 * t427) * t96 / 0.32e2 + (mrSges(6,2) * t594 + t169 * t604 + t497 - 0.32e2 * t513) * t92 / 0.32e2 + (mrSges(6,2) * t597 + t136 * t604 + t285 - 0.16e2 * t427) * t95 / 0.32e2 + (mrSges(6,2) * t596 + t134 * t604 + t498) * t90 / 0.32e2 + ((t362 + t587) * t612 + (16 * t557) + t495) * sin(t141) / 0.32e2 + ((t362 - t587) * t612 - (16 * t557) + t495) * sin(t140) / 0.32e2 + ((t360 + t588) * t612 + t304 - t389) * cos(t140) / 0.32e2 + ((t360 - t588) * t612 + t304 + t389) * cos(t141) / 0.32e2 + ((-mrSges(9,2) * qJDD(1) - (t282 * t397)) * t601 + t272 * t593 - Ifges(9,6) * t380) * t186 / 0.32e2 + ((-(mrSges(9,2) * t397) + t282 * qJDD(1)) * t601 + Ifges(9,6) * t592 - t272 * t380) * t190 / 0.32e2 + ((-0.64e2 * t359 * pkin(15) * t397 + (mrSges(5,3) * t593)) * pkin(5) + (-(t340 * t397) - t469) * t601 + (t337 * t592) - t336 * t380) * t252 / 0.32e2 + ((-(32 * mrSges(5,3) * t288) + t468 * t601) * pkin(5) + (-(t339 * t397) + t229) * t601 + t336 * t592 + (t337 * t380)) * t259 / 0.32e2 + (Ifges(9,4) * qJDD(1) - t224 * t397) * t124 + (0.32e2 * Ifges(6,3) * qJDD(4) + t281 * t448) * t138 / 0.32e2 + (-0.32e2 * t428 + 0.8e1 * t467) * t314 / 0.32e2 + (t437 - t641) * t208 / 0.32e2 + (t437 + t641) * t207 / 0.32e2 + (t434 - t643) * t213 / 0.32e2 + (t434 + t643) * t212 / 0.32e2 + (t492 - t627) * t254 / 0.32e2 + (t492 + t627) * t255 / 0.32e2 + t620 + ((t126 - t470) * t88 + (t126 + t470) * t87 + (t150 - t471) * t93 + (t150 + t471) * t94 + (t267 - t481) * t192 + (t267 + t481) * t191 - (t265 - t485) * t187 - (t265 + t485) * t188 + (-sin(t182) + sin(t181)) * t452 + t635 * t474 * t567) * pkin(1) + (t458 - t572 / 0.2e1) * (mrSges(6,2) * t477 + t266 + t323) + (t457 + t573 / 0.2e1) * (mrSges(6,2) * t478 + t266 - t323) + (t578 / 0.2e1 - t581 / 0.2e1) * (t266 - t490) + (-t580 / 0.2e1 + t579 / 0.2e1) * (t266 + t490) + (t48 / 0.4e1 - t49 / 0.4e1 - t312 / 0.2e1) * (Ifges(6,4) * qJDD(1) + t345 * t472) + (t50 + t51) * (t428 / 0.2e1 - t467 / 0.8e1) - 0.2e1 * (mrSges(8,1) * cos(t300) + t223 * sin(t174)) * t480 + (sin(t184) * t455 - sin(t183) * t567 / 0.2e1) * g(1) + t386 + (0.32e2 * g(1) * t420 + (g(2) * t598)) * sin(qJ(1)) / 0.32e2 + ((g(1) * t598) - 0.32e2 * g(2) * t420) * cos(qJ(1)) / 0.32e2 + (0.64e2 * Ifges(9,4) * t397 + t224 * t585) * t125 / 0.32e2 + (-((t279 + t332) * m(6) + t488 + t388) * t75 - ((t280 + t335) * m(6) - t466 + t500) * t68 + (-t54 / 0.2e1 + t53 / 0.2e1) * (mrSges(6,1) * t477 + t268 - t322) + (t52 / 0.2e1 - t55 / 0.2e1) * (mrSges(6,1) * t478 + t268 + t322)) * pkin(5) + (-t484 * t185 - t475 * t189 + 0.2e1 * (t637 - t638) * t391) * pkin(2) + (t494 - t628) * t261 / 0.32e2 + (t494 + t628) * t262 / 0.32e2 - t199 * t66 / 0.4e1 - t613 - 0.2e1 * (t469 / 0.2e1 + t166 * t391) * t540 + (t166 * t473 + t469) * t538 + qJDD(1) * (pkin(9) * t444 + Ifges(5,4)) * sin(t107) + (t220 * mrSges(6,1) + t435) * t450 + t490 * pkin(9) * t352 + ((t72 + t70) * (t268 - t129) + (t73 + t69) * (t268 + t129)) * pkin(1) / 0.2e1 - (Ifges(8,1) - Ifges(8,2)) * qJDD(1) * cos(t193) / 0.2e1 + Ifges(8,4) * qJDD(1) * sin(t193) + (-Ifges(7,4) * qJDD(1) - t344 * t474) * t132 + (-qJDD(1) * t416 - (4 * t343 * t474)) * t315 / 0.2e1 + (-(Ifges(7,5) * t368) - (Ifges(7,6) * qJDD(2)) + (mrSges(7,1) * qJDD(1) - (mrSges(7,2) * t474)) * t465) * t206 + (-t343 * qJDD(1) + t416 * t474) * t313 + ((Ifges(7,5) * qJDD(2)) - (Ifges(7,6) * t368) + ((mrSges(7,1) * t474) + mrSges(7,2) * qJDD(1)) * t465) * t211 + (((-t373 + t374) * m(6) + t384 + t521) * cos(t107) + (t373 + t374) * m(6) + m(11) * t376 + t521 + 0.2e1 * pkin(14) ^ 2 * m(7) - t384 + t496 + t509 + t510 + Ifges(7,1) + Ifges(8,1) + Ifges(9,1) + Ifges(7,2) + Ifges(8,2) + Ifges(9,2) + (2 * Ifges(2,3)) + t523 + t634) * qJDD(1) / 0.2e1 - t382 + t383 - (-(qJD(1) * t339 * t646) + t499) * t543 + t63 * t552 + ((t279 - t332) * m(6) - t488 + t388) * t568 + ((t280 - t335) * m(6) + t466 + t500) * t569 + (-t636 * (mrSges(11,1) * qJDD(1) - (mrSges(11,2) * t397)) + (t120 - t121) * (mrSges(11,1) * t397 + qJDD(1) * mrSges(11,2))) * pkin(4) - t644 + (mrSges(6,1) * t594 + t169 * t602 + t385 + t47) * t429 - t129 * pkin(9) * t349 + (mrSges(10,1) * t473 - t475) * t535 + (mrSges(10,2) * t473 + t484) * t536 + (-(t339 * t473) + t499) * t539 + (mrSges(6,2) * t595 + t171 * t604 - t385 + t47) * t430 + (-t559 - t563) * sin(t243) / 0.2e1 + (-t559 + t563) * sin(t244) / 0.2e1 + (-qJDD(1) * t615 + 0.4e1 * t338 * t397) * t159 / 0.2e1 + (qJDD(1) * t338 + t397 * t615) * t158 + (t358 - t560) * cos(t243) / 0.2e1 + (t358 + t560) * cos(t244) / 0.2e1 + (cos(t184) + cos(t183)) * g(2) * t455 - (mrSges(6,2) * t292 + t422) * t542 / 0.2e1 + (-mrSges(6,2) * t291 + t422) * t541 / 0.2e1 + (-t292 * mrSges(6,1) + t421) * t545 / 0.2e1 + (t291 * mrSges(6,1) + t421) * t544 / 0.2e1 + (-mrSges(6,2) * t220 + t436) * t531 / 0.2e1 - (mrSges(6,2) * t221 + t436) * t532 / 0.2e1 - (-t221 * mrSges(6,1) + t435) * t530 / 0.2e1 + 0.2e1 * mrSges(8,2) * sin(t300) * t480 + t376 * cos((2 * t231)) * t489 / 0.2e1; (t529 + t545 - t544 - t530) * (t248 + t269) + (t445 - t640) * t262 / 0.16e2 + (t445 + t640) * t261 / 0.16e2 + (t446 - t642) * t254 / 0.16e2 + (t446 + t642) * t255 / 0.16e2 + t416 * t313 * t549 + (((t188 + t187) * mrSges(8,2) / 0.2e1 - (t192 + t191) * mrSges(8,1) / 0.2e1) * t369 - 0.3e1 / 0.8e1 * (t73 + t70) * t46 + t635 * m(11) * t451 + (t94 + t93) * t425 + (t87 + t88) * t426) * pkin(1) + t344 * t132 * t548 + t72 * t413 + (0.16e2 * (-t230 * qJDD(1) + t225 * t512) * pkin(1) + t348 * t463 + t342 * t585 + 0.16e2 * t347 * g(3)) * t354 / 0.16e2 + t617 + (-t542 - t541) * t104 + (Ifges(7,4) * t133 + t343 * t315) * t369 + (pkin(15) * t503 + (-mrSges(5,3) * pkin(5) + t337) * qJDD(1) + t339 * g(3)) * t252 - t414 + t415 + ((t173 + Ifges(3,3) + Ifges(7,3) + Ifges(10,3) + t453) * qJDD(2)) - t417 + t418 + t626 + (-t347 * t512 - t341 * qJDD(1) + (t348 + t619) * g(3)) * t351 + t382 + t383 + t440 * t579 + 0.3e1 / 0.8e1 * t46 * t614 + 0.3e1 / 0.8e1 * (t580 + t571) * t45 - 0.3e1 / 0.8e1 * (t578 + t570) * t45 + qJDD(3) * t453 + t496 * t293 + 0.2e1 * (t166 * t476 + t339 * t263) * t538 + 0.2e1 * (t166 * t263 - t339 * t476) * t539 + 0.2e1 * (mrSges(10,1) * t476 - mrSges(10,2) * t263) * t535 + 0.2e1 * (t263 * mrSges(10,1) + mrSges(10,2) * t476) * t536 + (-mrSges(7,1) * t514 + g(3) * mrSges(7,2) + Ifges(7,5) * qJDD(1)) * t211 + (g(3) * mrSges(7,1) + mrSges(7,2) * t514 - Ifges(7,6) * qJDD(1)) * t206; (t549 - t368) * mrSges(10,2) * t536 + ((mrSges(4,2) / 0.2e1 + mrSges(11,2) / 0.2e1) * t369 + (t339 * t368)) * t539 + (0.16e2 * (-qJDD(1) * mrSges(5,3) + t359 * t512) * pkin(5) + t340 * t463 + t337 * t585 + g(3) * t622) * t252 / 0.16e2 + (t453 + t496) * t293 + t415 / 0.2e1 + t418 / 0.2e1 - t417 / 0.2e1 - t414 / 0.2e1 + (t458 + t457) * (t486 + (t549 - t366) * mrSges(6,2)) + (0.2e1 * t482 + ((2 * t366) + t369) * mrSges(6,1)) * t614 / 0.4e1 + (t530 / 0.2e1 + t450) * (-0.2e1 * t269 + t323) + (mrSges(10,1) * t536 - mrSges(10,2) * t535 + t166 * t539 + t339 * t538) * qJDD(2) + t626 + (-t166 * t538 / 0.2e1 - mrSges(10,1) * t535 / 0.2e1) * t235; (mrSges(6,1) * t164 + mrSges(6,2) * t238) * t92 / 0.2e1 + (qJDD(1) * t138 + qJDD(4)) * Ifges(6,3) + t386 + t284 * t429 + t285 * t430 + (-(t127 + 0.2e1 * t519) * t58 / 0.4e1 - (t128 - 0.2e1 * t517) * t54 / 0.4e1 + (t128 + 0.2e1 * t517) * t55 / 0.4e1 - (t127 - 0.2e1 * t519) * t59 / 0.4e1 + (t175 + t245) * t204 + (t175 + t246) * t205 - (t247 + t176) * t210 + (t248 + t176) * t209) * pkin(5) + ((t178 - 0.2e1 * t483) * t70 / 0.4e1 - (t178 + 0.2e1 * t483) * t73 / 0.4e1 - (t180 + 0.2e1 * t487) * t62 / 0.4e1 - (t180 - 0.2e1 * t487) * t65 / 0.4e1 - (t248 + t266) * t249 + (t247 + t266) * t250 + (t245 + t268) * t256 + (t246 + t268) * t257) * pkin(1) + (mrSges(6,1) * t238 - mrSges(6,2) * t164) * t98 / 0.2e1 + (mrSges(6,1) * t237 + mrSges(6,2) * t165) * t97 / 0.2e1 + (mrSges(6,1) * t165 - mrSges(6,2) * t237) * t91 / 0.2e1 + (-t63 + t66) * t552 + (-(t432 - t566) * t89 / 0.4e1 - (t432 + t566) * t90 / 0.4e1 + (t433 + t562) * t95 / 0.4e1 - (t433 - t562) * t96 / 0.4e1 + (t53 * t586 - mrSges(6,2) * t57 / 0.4e1) * pkin(5) + (mrSges(6,2) * t64 / 0.4e1 + t72 * t586) * pkin(1) + (-t48 / 0.8e1 + t49 / 0.8e1 + t312 / 0.4e1) * t345 + (-t50 / 0.4e1 - t51 / 0.4e1 + t314 / 0.2e1) * Ifges(6,4)) * t369 + t520 + t617 + (mrSges(6,1) * t464 + 0.32e2 * mrSges(6,2) * t480) * t349 / 0.32e2 + (-0.32e2 * mrSges(6,1) * t480 + mrSges(6,2) * t464) * t352 / 0.32e2 + (t91 + t92) * t287 / 0.32e2 + t644;];
tau = t1;
