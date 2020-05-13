% Calculate vector of inverse dynamics joint torques with ic for
% palh1m2IC
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m2IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2IC_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:46:01
% EndTime: 2020-05-02 23:46:48
% DurationCPUTime: 46.85s
% Computational Cost: add. (12888->967), mult. (29041->1373), div. (19->10), fcn. (23430->50), ass. (0->458)
t423 = cos(qJ(2));
t740 = t423 / 0.2e1;
t545 = qJD(3) + qJD(4);
t390 = qJD(2) + t545;
t411 = sin(qJ(5));
t420 = cos(qJ(5));
t414 = sin(qJ(2));
t422 = cos(qJ(3));
t577 = t414 * t422;
t413 = sin(qJ(3));
t578 = t413 * t423;
t314 = t577 + t578;
t288 = t314 * qJD(1);
t574 = t422 * t423;
t315 = -t413 * t414 + t574;
t291 = t315 * qJD(1);
t412 = sin(qJ(4));
t421 = cos(qJ(4));
t474 = t421 * t288 + t291 * t412;
t157 = t390 * t420 - t411 * t474;
t158 = t390 * t411 + t420 * t474;
t739 = -mrSges(5,1) * t390 - mrSges(6,1) * t157 + mrSges(6,2) * t158 + mrSges(5,3) * t474;
t409 = sin(qJ(7));
t418 = cos(qJ(7));
t312 = -t409 * t423 - t414 * t418;
t287 = t312 * qJD(1);
t564 = qJD(1) * t423;
t565 = qJD(1) * t414;
t290 = -t409 * t565 + t418 * t564;
t405 = sin(pkin(19));
t406 = cos(pkin(19));
t396 = t414 * pkin(1);
t369 = t396 - pkin(15);
t566 = qJD(1) * t369;
t166 = t566 + (-t287 * t406 - t290 * t405) * pkin(4);
t645 = sin(qJ(10));
t646 = cos(qJ(10));
t303 = t405 * t645 + t406 * t646;
t441 = -t405 * t646 + t406 * t645;
t132 = -t287 * t441 - t290 * t303;
t506 = -t303 * t287 + t290 * t441;
t72 = -mrSges(11,1) * t506 + mrSges(11,2) * t132;
t738 = m(11) * t166 + t72;
t737 = -mrSges(4,1) * t291 - mrSges(8,1) * t287 + mrSges(4,2) * t288 + mrSges(8,2) * t290;
t505 = -t288 * t412 + t421 * t291;
t196 = Ifges(5,4) * t505;
t615 = t390 * Ifges(5,5);
t103 = Ifges(5,1) * t474 + t196 + t615;
t580 = t412 * t422;
t364 = pkin(1) * t580;
t662 = pkin(1) * t413;
t373 = pkin(5) + t662;
t276 = t373 * t421 + t364;
t641 = pkin(5) * qJD(3);
t246 = qJD(2) * t276 + t421 * t641;
t229 = -pkin(9) * t390 - t246;
t331 = mrSges(6,1) * t411 + mrSges(6,2) * t420;
t457 = t229 * t331;
t197 = qJD(5) - t505;
t625 = t158 * Ifges(6,4);
t48 = t157 * Ifges(6,2) + t197 * Ifges(6,6) + t625;
t372 = pkin(5) * t413 + pkin(1);
t266 = -pkin(5) * t574 + t372 * t414 - pkin(15);
t569 = qJD(1) * t266;
t156 = Ifges(6,4) * t157;
t49 = t158 * Ifges(6,1) + t197 * Ifges(6,5) + t156;
t605 = t420 * t49;
t665 = t411 / 0.2e1;
t736 = -t605 / 0.2e1 + t48 * t665 - t103 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t474 - t196 / 0.2e1 - t457 - mrSges(5,2) * t569 - t615 / 0.2e1;
t632 = Ifges(3,4) * t423;
t735 = (-Ifges(3,1) * t414 - t632) * t740 - pkin(15) * (mrSges(3,1) * t423 - mrSges(3,2) * t414);
t622 = t505 * mrSges(5,3);
t176 = -mrSges(5,2) * t390 + t622;
t92 = -mrSges(6,2) * t197 + mrSges(6,3) * t157;
t93 = mrSges(6,1) * t197 - mrSges(6,3) * t158;
t734 = -t411 * t93 + t420 * t92 + t176;
t575 = t421 * t422;
t363 = pkin(1) * t575;
t275 = t373 * t412 - t363;
t244 = t275 * qJD(2) + t412 * t641;
t228 = pkin(11) * t390 + t244;
t88 = -pkin(9) * t505 - pkin(11) * t474 + t569;
t60 = t228 * t420 + t411 * t88;
t59 = -t228 * t411 + t420 * t88;
t604 = t420 * t59;
t481 = t411 * t60 + t604;
t733 = -m(6) * t481 - t411 * t92 - t420 * t93;
t579 = t413 * t421;
t470 = t579 + t580;
t559 = qJD(4) * t421;
t560 = qJD(4) * t412;
t214 = t373 * t559 + (qJD(3) * t470 + t422 * t560) * pkin(1);
t529 = pkin(5) * t559;
t655 = pkin(5) * t412;
t118 = t214 * qJD(2) + qJD(3) * t529 + t275 * qJDD(2) + qJDD(3) * t655;
t399 = qJDD(2) + qJDD(3);
t386 = qJDD(4) + t399;
t114 = pkin(11) * t386 + t118;
t362 = pkin(5) * t577;
t297 = pkin(5) * t578 + t362;
t561 = qJD(2) * t423;
t216 = qJD(2) * t362 + qJD(3) * t297 + t372 * t561;
t211 = t216 * qJD(1);
t133 = t266 * qJDD(1) + t211;
t552 = qJD(1) * qJD(2);
t323 = qJDD(1) * t423 - t414 * t552;
t517 = t423 * t552;
t324 = -qJDD(1) * t414 - t517;
t446 = t314 * qJD(3);
t171 = -qJD(1) * t446 + t323 * t422 + t324 * t413;
t445 = t315 * qJD(3);
t172 = qJD(1) * t445 + t323 * t413 - t324 * t422;
t57 = qJD(4) * t505 + t171 * t412 + t172 * t421;
t58 = -qJD(4) * t474 + t171 * t421 - t172 * t412;
t13 = -pkin(9) * t58 - pkin(11) * t57 + t133;
t4 = qJD(5) * t59 + t114 * t420 + t13 * t411;
t5 = -qJD(5) * t60 - t114 * t411 + t13 * t420;
t493 = t4 * t420 - t411 * t5;
t19 = qJD(5) * t157 + t386 * t411 + t420 * t57;
t54 = qJDD(5) - t58;
t11 = mrSges(6,1) * t54 - mrSges(6,3) * t19;
t20 = -qJD(5) * t158 + t386 * t420 - t411 * t57;
t12 = -mrSges(6,2) * t54 + mrSges(6,3) * t20;
t557 = qJD(5) * t420;
t558 = qJD(5) * t411;
t724 = -t411 * t11 + t420 * t12 - t93 * t557 - t92 * t558;
t732 = m(6) * (-qJD(5) * t481 + t493) + t724;
t548 = pkin(19) - qJ(7);
t391 = -qJ(2) + t548;
t349 = pkin(4) * sin(t391);
t512 = -qJ(10) + t548;
t368 = -qJ(2) + t512;
t350 = sin(t368);
t351 = cos(t368);
t508 = -t350 * mrSges(11,1) + t351 * mrSges(11,2);
t731 = m(11) * (t349 - t396) + t508;
t688 = m(5) + m(6);
t511 = -mrSges(6,1) * t420 + t411 * mrSges(6,2);
t695 = m(6) * pkin(9);
t309 = mrSges(5,1) - t511 + t695;
t387 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t706 = -t309 * t421 - t387 * t412;
t227 = pkin(5) * t688 + mrSges(4,1) - t706;
t243 = -t309 * t412 + t387 * t421;
t241 = -mrSges(4,2) + t243;
t123 = t227 * t422 + t241 * t413;
t122 = -t227 * t413 + t241 * t422;
t614 = t390 * Ifges(5,6);
t631 = Ifges(5,4) * t474;
t101 = Ifges(5,2) * t505 + t614 + t631;
t623 = t197 * Ifges(6,3);
t624 = t158 * Ifges(6,5);
t626 = t157 * Ifges(6,6);
t47 = t623 + t624 + t626;
t620 = t244 * mrSges(5,3);
t714 = t60 * mrSges(6,2);
t715 = t59 * mrSges(6,1);
t730 = t714 + t620 + t101 / 0.2e1 - t47 / 0.2e1 + t631 / 0.2e1 - mrSges(5,1) * t569 + t614 / 0.2e1 - t715;
t729 = m(8) + m(4);
t725 = mrSges(11,1) * t351 + mrSges(11,2) * t350;
t401 = qJD(2) + qJD(7);
t388 = qJD(10) + t401;
t668 = -t388 / 0.2e1;
t682 = -t132 / 0.2e1;
t722 = -t166 * mrSges(11,2) + Ifges(11,1) * t682 + Ifges(11,5) * t668;
t400 = qJD(2) + qJD(8);
t389 = qJD(9) + t400;
t408 = sin(qJ(8));
t417 = cos(qJ(8));
t310 = -t408 * t423 - t414 * t417;
t286 = t310 * qJD(1);
t289 = -t408 * t565 + t417 * t564;
t407 = sin(qJ(9));
t416 = cos(qJ(9));
t476 = -t286 * t416 + t289 * t407;
t475 = t286 * t407 + t289 * t416;
t621 = t475 * Ifges(10,4);
t100 = Ifges(10,2) * t476 + t389 * Ifges(10,6) - t621;
t195 = Ifges(10,4) * t476;
t102 = -Ifges(10,1) * t475 + t389 * Ifges(10,5) + t195;
t598 = qJD(1) * pkin(15);
t256 = -pkin(2) * t286 - t598;
t397 = qJDD(2) + qJDD(8);
t556 = qJD(9) * t407;
t264 = (-t397 * t416 + t400 * t556) * pkin(2);
t555 = qJD(9) * t416;
t265 = (-t397 * t407 - t400 * t555) * pkin(2);
t404 = qJ(2) + qJ(8);
t393 = qJ(9) + t404;
t370 = sin(t393);
t371 = cos(t393);
t385 = qJDD(9) + t397;
t415 = sin(qJ(1));
t424 = cos(qJ(1));
t495 = g(1) * t424 + g(2) * t415;
t450 = t310 * qJD(8);
t167 = qJD(1) * t450 + t323 * t417 + t324 * t408;
t472 = t408 * t414 - t417 * t423;
t449 = t472 * qJD(8);
t168 = qJD(1) * t449 - t323 * t408 + t324 * t417;
t55 = qJD(9) * t476 - t167 * t416 - t168 * t407;
t56 = qJD(9) * t475 + t167 * t407 - t168 * t416;
t676 = -t475 / 0.2e1;
t721 = -t265 * mrSges(10,2) + Ifges(10,5) * t55 + Ifges(10,6) * t56 - (g(3) * mrSges(10,1) - mrSges(10,2) * t495) * t370 + (-mrSges(10,1) * t495 - g(3) * mrSges(10,2)) * t371 + t264 * mrSges(10,1) + Ifges(10,3) * t385 - (Ifges(10,5) * t476 + Ifges(10,6) * t475) * t389 / 0.2e1 - t256 * (-mrSges(10,1) * t475 + mrSges(10,2) * t476) + t100 * t676 + (Ifges(10,1) * t476 + t621) * t475 / 0.2e1 - (Ifges(10,2) * t475 + t102 + t195) * t476 / 0.2e1;
t127 = Ifges(11,4) * t506;
t448 = t312 * qJD(7);
t169 = qJD(1) * t448 + t323 * t418 + t324 * t409;
t471 = t409 * t414 - t418 * t423;
t447 = t471 * qJD(7);
t170 = qJD(1) * t447 - t323 * t409 + t324 * t418;
t272 = t303 * qJD(10);
t273 = t441 * qJD(10);
t27 = -t169 * t303 - t170 * t441 - t272 * t287 + t273 * t290;
t28 = t169 * t441 - t170 * t303 + t272 * t290 + t273 * t287;
t398 = qJDD(2) + qJDD(7);
t380 = qJDD(10) + t398;
t644 = pkin(1) * qJD(2);
t536 = t409 * t644;
t657 = pkin(4) * t405;
t298 = t401 * t657 + t536;
t534 = t418 * t644;
t656 = pkin(4) * t406;
t299 = t401 * t656 + t534;
t137 = t298 * t303 + t299 * t441;
t595 = t137 * mrSges(11,3);
t642 = pkin(1) * qJD(7);
t528 = qJD(2) * t642;
t601 = pkin(1) * qJDD(2);
t304 = -t409 * t528 + t418 * t601;
t254 = t398 * t656 + t304;
t305 = t409 * t601 + t418 * t528;
t255 = t398 * t657 + t305;
t62 = -t254 * t441 - t255 * t303 - t272 * t299 + t273 * t298;
t63 = -t254 * t303 + t255 * t441 + t272 * t298 + t273 * t299;
t65 = t132 * Ifges(11,4) + Ifges(11,2) * t506 + t388 * Ifges(11,6);
t66 = t132 * Ifges(11,1) + t388 * Ifges(11,5) + t127;
t683 = -t506 / 0.2e1;
t719 = t63 * mrSges(11,1) - t62 * mrSges(11,2) + Ifges(11,5) * t27 + Ifges(11,6) * t28 + Ifges(11,3) * t380 + (t66 + t127) * t683 - (t166 * mrSges(11,1) + Ifges(11,4) * t682 + Ifges(11,2) * t683 + Ifges(11,6) * t668 + t595 - t65 / 0.2e1) * t132;
t718 = m(9) + m(3);
t410 = sin(qJ(6));
t666 = t410 / 0.2e1;
t664 = -t414 / 0.2e1;
t716 = g(3) * t508;
t461 = pkin(1) * (-t303 * t418 + t409 * t441);
t711 = (-t272 * t406 + t273 * t405) * pkin(4) - qJD(2) * t461;
t462 = pkin(1) * (t303 * t409 + t418 * t441);
t710 = (t272 * t405 + t273 * t406) * pkin(4) - qJD(2) * t462;
t547 = qJ(4) + pkin(18);
t497 = pkin(19) + qJ(3) + t547;
t573 = qJ(6) - qJ(2);
t452 = t497 + t573;
t439 = -0.2e1 * qJ(7) - pkin(20) + t452;
t453 = t497 - t573;
t440 = pkin(20) + t453;
t709 = cos(qJ(10) - t439) + cos(qJ(10) - t440);
t226 = t314 * t421 + t315 * t412;
t225 = -t412 * t314 + t315 * t421;
t236 = qJD(2) * t315 + t445;
t237 = -qJD(2) * t314 - t446;
t87 = qJD(4) * t225 + t236 * t421 + t237 * t412;
t466 = t226 * t557 + t411 * t87;
t707 = -t59 * t557 - t60 * t558;
t705 = t688 + t729;
t704 = -mrSges(11,3) - mrSges(10,3) - t331 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t419 = cos(qJ(6));
t491 = -t419 * mrSges(7,1) + t410 * mrSges(7,2);
t702 = pkin(14) * m(7) + t491;
t703 = -mrSges(2,1) + t702 + (-t705 - t718 - m(11) - m(10)) * pkin(15) - t731;
t643 = pkin(1) * qJD(3);
t531 = t413 * t643;
t215 = t363 * t545 - t373 * t560 - t412 * t531;
t119 = -pkin(5) * (qJD(3) * t560 - qJDD(3) * t421) + qJD(2) * t215 + qJDD(2) * t276;
t701 = -t5 * mrSges(6,1) + t4 * mrSges(6,2);
t700 = pkin(9) * t474 - pkin(11) * t505;
t425 = qJD(1) ^ 2;
t697 = pkin(2) * m(10);
t694 = t19 / 0.2e1;
t693 = t20 / 0.2e1;
t690 = t54 / 0.2e1;
t681 = -t157 / 0.2e1;
t680 = -t158 / 0.2e1;
t679 = t158 / 0.2e1;
t678 = -t197 / 0.2e1;
t672 = t288 / 0.2e1;
t671 = t289 / 0.2e1;
t670 = t290 / 0.2e1;
t663 = pkin(1) * t409;
t661 = pkin(1) * t418;
t660 = pkin(1) * t423;
t659 = pkin(2) * (-mrSges(10,1) * t476 - mrSges(10,2) * t475);
t658 = pkin(4) * cos(t391);
t654 = pkin(5) * t421;
t10 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t115 = -pkin(9) * t386 - t119;
t134 = t243 * t413 - t422 * t706;
t135 = t243 * t422 + t413 * t706;
t483 = Ifges(6,5) * t420 - Ifges(6,6) * t411;
t458 = t197 * t483;
t630 = Ifges(6,4) * t411;
t488 = Ifges(6,1) * t420 - t630;
t459 = t158 * t488;
t629 = Ifges(6,4) * t420;
t486 = -Ifges(6,2) * t411 + t629;
t460 = t157 * t486;
t516 = -t558 / 0.2e1;
t521 = t605 / 0.2e1;
t636 = mrSges(6,3) * t420;
t8 = t19 * Ifges(6,4) + t20 * Ifges(6,2) + t54 * Ifges(6,6);
t9 = t19 * Ifges(6,1) + t20 * Ifges(6,4) + t54 * Ifges(6,5);
t435 = -t118 * mrSges(5,2) + t4 * t636 + t420 * t8 / 0.2e1 + t9 * t665 + t119 * mrSges(5,1) + (Ifges(6,1) * t411 + t629) * t694 + (Ifges(6,2) * t420 + t630) * t693 + Ifges(5,3) * t386 + (Ifges(6,5) * t411 + Ifges(6,6) * t420) * t690 + t48 * t516 + Ifges(5,6) * t58 + Ifges(5,5) * t57 + t115 * t511 + (t457 + t521) * qJD(5) + (t460 + t459 + t458) * qJD(5) / 0.2e1;
t81 = -t246 * t411 + t420 * t700;
t82 = t246 * t420 + t411 * t700;
t650 = (-t739 * t244 + (t622 - t176) * t246 + (m(6) * (t493 + t707) + t724) * pkin(11) + t435 - t115 * t695 + (-t626 / 0.2e1 - t624 / 0.2e1 - t623 / 0.2e1 + t730) * t474 + (-g(3) * t134 - t135 * t495) * t423 + (-t197 * t604 + (-t197 * t60 - t5) * t411) * mrSges(6,3) + (-g(3) * t135 + t134 * t495) * t414 - m(6) * (t229 * t244 + t59 * t81 + t60 * t82) - pkin(9) * t10 - t82 * t92 - t81 * t93 + (-t460 / 0.2e1 - t459 / 0.2e1 - t458 / 0.2e1 + t736) * t505) / pkin(10);
t597 = t506 * mrSges(11,3);
t120 = -mrSges(11,2) * t388 + t597;
t121 = mrSges(11,1) * t388 - mrSges(11,3) * t132;
t136 = t298 * t441 - t299 * t303;
t571 = t725 * t424;
t572 = t725 * t415;
t649 = (-t137 * t121 - g(1) * t571 - g(2) * t572 - t716 + (-t120 + t597) * t136 + t722 * t506 + t719) / pkin(8);
t647 = -qJD(6) / 0.2e1;
t639 = m(11) * (-t658 - t660);
t637 = mrSges(6,3) * t411;
t635 = mrSges(8,3) * t290;
t634 = mrSges(10,3) * t400;
t633 = Ifges(3,4) * t414;
t628 = Ifges(7,4) * t410;
t627 = Ifges(7,4) * t419;
t619 = t246 * mrSges(5,3);
t618 = t288 * Ifges(4,4);
t617 = t289 * Ifges(9,4);
t616 = t290 * Ifges(8,4);
t613 = t407 * (-mrSges(10,2) * t385 + mrSges(10,3) * t56);
t612 = t410 * Ifges(7,1);
t607 = t419 * Ifges(7,2);
t596 = t136 * mrSges(11,3);
t594 = t226 * t411;
t593 = t226 * t420;
t590 = (-mrSges(8,2) * t401 + mrSges(8,3) * t287) * t418;
t589 = t266 * t425;
t485 = t607 + t628;
t582 = t410 * (Ifges(7,6) * qJD(6) + qJD(1) * t485);
t581 = t412 * t413;
t376 = qJD(1) * t627;
t576 = t419 * (Ifges(7,5) * qJD(6) + qJD(1) * t612 + t376);
t359 = pkin(1) * t517;
t279 = t369 * qJDD(1) + t359;
t274 = t372 * t423 + t362;
t568 = qJD(1) * t274;
t567 = qJD(1) * t297;
t295 = -pkin(1) * t581 + t363;
t563 = qJD(2) * t295;
t296 = pkin(1) * t579 + t364;
t562 = qJD(2) * t296;
t554 = qJDD(1) * pkin(15);
t551 = qJD(1) * qJD(6);
t550 = qJD(2) * qJD(3);
t546 = qJ(7) + pkin(20);
t544 = t256 * t697;
t543 = sin(t404) * t697;
t542 = pkin(2) * t634;
t541 = m(11) * t658;
t540 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t54;
t537 = pkin(1) * t564;
t535 = t413 * t644;
t533 = t422 * t644;
t532 = pkin(1) * t561;
t515 = -qJ(8) + pkin(17) + qJ(3);
t514 = -t552 / 0.2e1;
t510 = t407 * t542;
t509 = t416 * t542;
t503 = mrSges(4,3) * t535;
t502 = mrSges(4,3) * t533;
t501 = mrSges(8,3) * t536;
t500 = mrSges(8,3) * t534;
t499 = t369 * t425 * t660;
t498 = t546 - t573;
t490 = mrSges(7,1) * t410 + mrSges(7,2) * t419;
t333 = -t418 * mrSges(8,1) + mrSges(8,2) * t409;
t330 = -mrSges(8,1) * t409 - mrSges(8,2) * t418;
t332 = -t417 * mrSges(9,1) + t408 * mrSges(9,2);
t329 = -mrSges(9,1) * t408 - mrSges(9,2) * t417;
t489 = t423 * Ifges(3,1) - t633;
t487 = -t414 * Ifges(3,2) + t632;
t484 = -Ifges(3,5) * t414 - Ifges(3,6) * t423;
t482 = Ifges(7,5) * t419 - Ifges(7,6) * t410;
t480 = -t411 * t59 + t420 * t60;
t230 = t470 * t331;
t231 = (-t575 + t581) * t331;
t477 = -t230 * t414 - t231 * t423;
t223 = -t310 * t416 - t407 * t472;
t473 = t310 * t407 - t416 * t472;
t465 = t226 * t558 - t420 * t87;
t464 = pkin(14) * t490;
t456 = t410 * (Ifges(7,1) * t419 - t628);
t455 = t414 * (-Ifges(3,2) * t423 - t633);
t184 = (-t287 * t405 + t290 * t406) * pkin(4);
t444 = -qJ(7) + t453;
t443 = -qJ(7) + t452;
t177 = mrSges(10,1) * t389 + mrSges(10,3) * t475;
t189 = t286 * Ifges(9,2) + t400 * Ifges(9,6) + t617;
t269 = Ifges(9,4) * t286;
t192 = t289 * Ifges(9,1) + t400 * Ifges(9,5) + t269;
t434 = t475 * t510 - t289 * (Ifges(9,1) * t286 - t617) / 0.2e1 - t289 * t659 + t189 * t671 + pkin(2) * t177 * t556 + (mrSges(9,1) * t289 + mrSges(9,2) * t286) * t598 + g(3) * t543 - t289 * t544 + Ifges(9,5) * t167 + Ifges(9,6) * t168 + Ifges(9,3) * t397 - t400 * (Ifges(9,5) * t286 - Ifges(9,6) * t289) / 0.2e1 + (t495 * cos(t404) - t264 * t416 - t265 * t407) * t697 - (-Ifges(9,2) * t289 + t192 + t269) * t286 / 0.2e1 + t721;
t190 = t287 * Ifges(8,2) + t401 * Ifges(8,6) + t616;
t270 = Ifges(8,4) * t287;
t193 = t290 * Ifges(8,1) + t401 * Ifges(8,5) + t270;
t433 = -t305 * mrSges(8,2) - (mrSges(8,1) * t290 + mrSges(8,2) * t287) * t566 + t190 * t670 - t290 * (Ifges(8,1) * t287 - t616) / 0.2e1 + Ifges(8,6) * t170 + Ifges(8,5) * t169 - t401 * (Ifges(8,5) * t287 - Ifges(8,6) * t290) / 0.2e1 + t287 * t500 + t304 * mrSges(8,1) + Ifges(8,3) * t398 - (-Ifges(8,2) * t290 + t193 + t270) * t287 / 0.2e1 + (t596 + t722) * t506 + t719;
t402 = qJD(2) + qJD(3);
t191 = t291 * Ifges(4,2) + t402 * Ifges(4,6) + t618;
t271 = Ifges(4,4) * t291;
t194 = t288 * Ifges(4,1) + t402 * Ifges(4,5) + t271;
t306 = (-qJDD(2) * t422 + t413 * t550) * pkin(1);
t307 = (qJDD(2) * t413 + t422 * t550) * pkin(1);
t432 = -(mrSges(4,1) * t288 + mrSges(4,2) * t291) * t566 - t288 * (Ifges(4,1) * t291 - t618) / 0.2e1 + t191 * t672 + t435 - t5 * t637 - t288 * t502 + t291 * t503 + t707 * mrSges(6,3) + Ifges(4,6) * t171 + Ifges(4,5) * t172 - t306 * mrSges(4,2) + t307 * mrSges(4,1) + Ifges(4,3) * t399 - t402 * (Ifges(4,5) * t291 - Ifges(4,6) * t288) / 0.2e1 - (-Ifges(4,2) * t288 + t194 + t271) * t291 / 0.2e1 + (Ifges(6,5) * t680 + Ifges(6,6) * t681 + Ifges(6,3) * t678 + t730) * t474 + (t483 * t678 + t486 * t681 + t488 * t680 + t59 * t636 + t60 * t637 + t619 + t736) * t505;
t430 = 0.1e1 / pkin(3);
t375 = -pkin(9) - t654;
t354 = cos(t498);
t353 = cos(qJ(9) - t515);
t335 = t656 + t661;
t334 = t657 + t663;
t322 = qJDD(1) * t410 + t419 * t551;
t321 = qJDD(1) * t419 - t410 * t551;
t285 = Ifges(3,5) * qJD(2) + qJD(1) * t489;
t283 = Ifges(3,6) * qJD(2) + qJD(1) * t487;
t250 = -mrSges(4,2) * t402 + mrSges(4,3) * t291;
t249 = mrSges(8,1) * t401 - t635;
t248 = mrSges(4,1) * t402 - mrSges(4,3) * t288;
t235 = qJD(2) * t471 + t447;
t234 = qJD(2) * t312 + t448;
t233 = qJD(2) * t472 + t449;
t232 = qJD(2) * t310 + t450;
t213 = (-t303 * t405 - t406 * t441) * pkin(4);
t212 = (-t303 * t406 + t405 * t441) * pkin(4);
t175 = -mrSges(10,2) * t389 + mrSges(10,3) * t476;
t174 = -t303 * t334 - t335 * t441;
t173 = -t303 * t335 + t334 * t441;
t147 = -pkin(2) * t168 - t554;
t111 = -mrSges(5,1) * t505 + mrSges(5,2) * t474;
t106 = -mrSges(3,2) + t123 + t329 + t330;
t96 = pkin(1) * t705 + mrSges(3,1) - t122 - t332 - t333;
t95 = qJD(7) * t462 + t272 * t334 + t273 * t335;
t94 = qJD(7) * t461 - t272 * t335 + t273 * t334;
t91 = t700 + t567;
t86 = qJD(4) * t226 + t236 * t412 - t421 * t237;
t85 = qJD(9) * t473 + t232 * t407 - t233 * t416;
t84 = qJD(9) * t223 - t232 * t416 - t233 * t407;
t80 = t411 * t91 + t420 * t562;
t79 = -t411 * t562 + t420 * t91;
t73 = (-t169 * t405 - t170 * t406) * pkin(4) + t279;
t45 = t234 * t441 - t235 * t303 - t272 * t471 + t273 * t312;
t44 = -t234 * t303 - t235 * t441 - t272 * t312 - t273 * t471;
t41 = -mrSges(5,2) * t386 + mrSges(5,3) * t58;
t40 = mrSges(5,1) * t386 - mrSges(5,3) * t57;
t38 = mrSges(10,1) * t385 - mrSges(10,3) * t55;
t22 = -mrSges(11,2) * t380 + mrSges(11,3) * t28;
t21 = mrSges(11,1) * t380 - mrSges(11,3) * t27;
t1 = [(pkin(15) ^ 2 * t718 + pkin(14) * t702 + Ifges(2,3)) * qJDD(1) + t737 * t532 + (Ifges(7,1) * t322 + Ifges(7,4) * t321) * t666 + t474 * (Ifges(5,1) * t87 - Ifges(5,4) * t86) / 0.2e1 + (-mrSges(8,1) * t279 + mrSges(8,3) * t305 + Ifges(8,4) * t169 + Ifges(8,2) * t170 + Ifges(8,6) * t398) * t312 + (-mrSges(11,1) * t73 + mrSges(11,3) * t62 + Ifges(11,4) * t27 + Ifges(11,2) * t28 + Ifges(11,6) * t380) * (-t303 * t312 - t441 * t471) - (mrSges(8,2) * t279 - mrSges(8,3) * t304 + Ifges(8,1) * t169 + Ifges(8,4) * t170 + Ifges(8,5) * t398) * t471 - (-mrSges(9,2) * t554 + Ifges(9,1) * t167 + Ifges(9,4) * t168 + Ifges(9,5) * t397) * t472 + t419 * (Ifges(7,4) * t322 + Ifges(7,2) * t321) / 0.2e1 + (m(11) * t73 - mrSges(11,1) * t28 + mrSges(11,2) * t27) * ((-t312 * t406 + t405 * t471) * pkin(4) + t369) + (mrSges(11,2) * t73 - mrSges(11,3) * t63 + Ifges(11,1) * t27 + Ifges(11,4) * t28 + Ifges(11,5) * t380) * (t303 * t471 - t312 * t441) + (-mrSges(10,1) * t147 + mrSges(10,3) * t265 + Ifges(10,4) * t55 + Ifges(10,2) * t56 + Ifges(10,6) * t385) * t223 + (-Ifges(6,3) * t690 - t540 / 0.2e1 - Ifges(6,6) * t693 - Ifges(6,5) * t694 + t118 * mrSges(5,3) + Ifges(5,4) * t57 + Ifges(5,2) * t58 - t133 * mrSges(5,1) + Ifges(5,6) * t386 + t701) * t225 - (mrSges(10,2) * t147 - mrSges(10,3) * t264 + Ifges(10,1) * t55 + Ifges(10,4) * t56 + Ifges(10,5) * t385) * t473 + (m(10) * t147 - mrSges(10,1) * t56 + mrSges(10,2) * t55) * (-pkin(2) * t310 - pkin(15)) + (mrSges(10,1) * t370 + mrSges(10,2) * t371 + t106 * t423 - t414 * t96 - t543) * (g(1) * t415 - g(2) * t424) + (t419 * (-Ifges(7,2) * t410 + t627) + t456) * t551 / 0.2e1 + t506 * (Ifges(11,4) * t44 + Ifges(11,2) * t45) / 0.2e1 + (Ifges(3,4) * t323 + Ifges(3,2) * t324) * t664 + t729 * (t279 + t359) * t369 + (t576 / 0.2e1 + t482 * qJD(6) / 0.2e1) * qJD(6) + t476 * (Ifges(10,4) * t84 + Ifges(10,2) * t85) / 0.2e1 + t505 * (Ifges(5,4) * t87 - Ifges(5,2) * t86) / 0.2e1 + (mrSges(9,1) * t554 + Ifges(9,4) * t167 + Ifges(9,2) * t168 + Ifges(9,6) * t397) * t310 + (-mrSges(4,1) * t237 - mrSges(8,1) * t235 + mrSges(4,2) * t236 + mrSges(8,2) * t234) * t566 + (t133 * mrSges(5,2) - t119 * mrSges(5,3) + Ifges(5,1) * t57 + Ifges(5,4) * t58 + Ifges(5,5) * t386 + t115 * t331 + t483 * t690 + t486 * t693 + t488 * t694 + t49 * t516) * t226 - (mrSges(3,1) * t414 + mrSges(3,2) * t423) * t554 - t87 * t619 - t86 * t620 + t582 * t647 - t85 * t510 + t324 * t487 / 0.2e1 + t323 * t489 / 0.2e1 + t321 * t485 / 0.2e1 + (Ifges(8,1) * t234 + Ifges(8,4) * t235) * t670 + (Ifges(9,1) * t232 + Ifges(9,4) * t233) * t671 + t738 * (t532 + (-t234 * t405 - t235 * t406) * pkin(4)) + (Ifges(4,1) * t236 + Ifges(4,4) * t237) * t672 + (Ifges(10,1) * t84 + Ifges(10,4) * t85) * t676 + (-Ifges(6,1) * t465 - Ifges(6,4) * t466 + Ifges(6,5) * t86) * t679 + (Ifges(3,1) * t323 + Ifges(3,4) * t324) * t740 - t234 * t500 + (-mrSges(4,1) * t171 - mrSges(8,1) * t170 + mrSges(4,2) * t172 + mrSges(8,2) * t169) * t369 + (t285 * t664 + t484 * qJD(2) / 0.2e1) * qJD(2) + (-t4 * t594 + t465 * t59 - t466 * t60 - t5 * t593) * mrSges(6,3) - t45 * t595 + t735 * t552 - t233 * t544 + t86 * t715 - t233 * t659 + (mrSges(5,1) * t86 + mrSges(5,2) * t87) * t569 + t9 * t593 / 0.2e1 - t8 * t594 / 0.2e1 + (m(6) * (qJD(5) * t480 + t4 * t411 + t420 * t5) + t92 * t557 + t411 * t12 + t420 * t11 - t93 * t558) * (-pkin(9) * t225 - pkin(11) * t226 + t266) - (-mrSges(9,1) * t233 + mrSges(9,2) * t232) * t598 + t322 * (t612 + t627) / 0.2e1 - t237 * t502 - t733 * (pkin(9) * t86 - pkin(11) * t87 + t216) + (t211 + t133) * m(5) * t266 - t236 * t503 - t283 * t561 / 0.2e1 + t464 * t551 - t44 * t596 + t235 * t501 + t84 * t509 + (g(1) * t704 + g(2) * t703) * t424 + (-g(1) * t703 + g(2) * t704) * t415 - t466 * t48 / 0.2e1 + t45 * t65 / 0.2e1 + t44 * t66 / 0.2e1 + t86 * t47 / 0.2e1 - t86 * t714 + t85 * t100 / 0.2e1 - t86 * t101 / 0.2e1 + t84 * t102 / 0.2e1 + t87 * t103 / 0.2e1 + (Ifges(3,5) * t423 + 0.2e1 * Ifges(3,6) * t664) * qJDD(2) + (0.2e1 * Ifges(7,5) * t666 + Ifges(7,6) * t419) * qJDD(6) + t132 * (Ifges(11,1) * t44 + Ifges(11,4) * t45) / 0.2e1 + t166 * (-mrSges(11,1) * t45 + mrSges(11,2) * t44) + (mrSges(4,2) * t279 - mrSges(4,3) * t307 + Ifges(4,1) * t172 + Ifges(4,4) * t171 + Ifges(4,5) * t399) * t314 + t216 * t111 + t232 * t192 / 0.2e1 + t233 * t189 / 0.2e1 + t234 * t193 / 0.2e1 + t235 * t190 / 0.2e1 + t236 * t194 / 0.2e1 + t237 * t191 / 0.2e1 + t256 * (-mrSges(10,1) * t85 + mrSges(10,2) * t84) + t266 * (-mrSges(5,1) * t58 + mrSges(5,2) * t57) + t286 * (Ifges(9,4) * t232 + Ifges(9,2) * t233) / 0.2e1 + t287 * (Ifges(8,4) * t234 + Ifges(8,2) * t235) / 0.2e1 + t291 * (Ifges(4,4) * t236 + Ifges(4,2) * t237) / 0.2e1 + pkin(14) * (-mrSges(7,1) * t321 + mrSges(7,2) * t322) + t455 * t514 + t87 * t521 + (-mrSges(4,1) * t279 + mrSges(4,3) * t306 + Ifges(4,4) * t172 + Ifges(4,2) * t171 + Ifges(4,6) * t399) * t315 + t388 * (Ifges(11,5) * t44 + Ifges(11,6) * t45) / 0.2e1 + t389 * (Ifges(10,5) * t84 + Ifges(10,6) * t85) / 0.2e1 + t390 * (Ifges(5,5) * t87 - Ifges(5,6) * t86) / 0.2e1 + t400 * (Ifges(9,5) * t232 + Ifges(9,6) * t233) / 0.2e1 + t401 * (Ifges(8,5) * t234 + Ifges(8,6) * t235) / 0.2e1 + t402 * (Ifges(4,5) * t236 + Ifges(4,6) * t237) / 0.2e1 + t229 * (mrSges(6,1) * t466 - mrSges(6,2) * t465) + t157 * (-Ifges(6,4) * t465 - Ifges(6,2) * t466 + Ifges(6,6) * t86) / 0.2e1 + t197 * (-Ifges(6,5) * t465 - Ifges(6,6) * t466 + Ifges(6,3) * t86) / 0.2e1 + (mrSges(3,1) * t324 + mrSges(9,1) * t168 - mrSges(3,2) * t323 - mrSges(9,2) * t167) * pkin(15); (-t175 * t555 - t38 * t416 - t613) * pkin(2) - t737 * t537 + (m(6) * t115 + t10) * (-pkin(9) - t276) + ((-pkin(3) * t354 - pkin(1) * cos(-t573)) * t430 * ((-g(3) * t333 + t330 * t495) * t414 + (-g(3) * t330 - t333 * t495) * t423 + t433 + t710 * t121 + t711 * t120 + (-t590 + (t249 + t635) * t409) * t644 - t184 * t72 + t212 * t21 + t213 * t22 - t716 - g(2) * (-t415 * t541 + t572) - g(1) * (-t424 * t541 + t571) + (-g(3) * t349 + t136 * t710 - t137 * t711 - t166 * t184 + t212 * t63 + t213 * t62) * m(11)) + pkin(1) * sin(t546) / pkin(7) * (Ifges(7,5) * t322 + Ifges(7,6) * t321 + Ifges(7,3) * qJDD(6) + t491 * g(3) + (-t576 / 0.2e1 + t582 / 0.2e1 - t419 * t376 / 0.2e1 + t482 * t647 + (-t464 - t456 / 0.2e1 + t607 * t666) * qJD(1)) * qJD(1) + t495 * t490)) / t354 + (t118 * t275 + t119 * t276 + t215 * t246 - t274 * t589) * m(5) + m(11) * (t136 * t95 - t137 * t94 + t173 * t63 + t174 * t62) - t476 * t509 + (-m(6) * t229 - t739) * t215 + (mrSges(8,1) * t398 - mrSges(8,3) * t169) * t661 + (mrSges(4,1) * t399 - mrSges(4,3) * t172) * t662 + (-mrSges(8,2) * t398 + mrSges(8,3) * t170) * t663 - g(1) * (t424 * t639 + t571) - g(2) * (t415 * t639 + t572) + ((t304 * t418 + t305 * t409) * pkin(1) - t499) * m(8) + (t248 * t643 - pkin(1) * (-mrSges(4,2) * t399 + mrSges(4,3) * t171)) * t422 + t432 - t738 * (t184 + t537) - g(3) * t731 + t732 * (pkin(11) + t275) + t290 * t501 - (-g(3) * t96 - t106 * t495) * t414 + t433 + t434 + (((-(cos(t444) + cos(t443)) * pkin(1) + (-cos(t439) - cos(t440)) * pkin(3)) * pkin(4) + ((cos(-qJ(10) + t444) + cos(-qJ(10) + t443)) * pkin(1) + t709 * pkin(3)) * pkin(8)) * t649 - (pkin(1) * (sin(qJ(10) - t573) + sin(qJ(10) + t573)) + (-sin(-qJ(10) + t498) + sin(qJ(10) + t498)) * pkin(3)) * pkin(4) * t650) / t709 * t430 + (-g(3) * t106 + t495 * t96) * t423 + (t455 / 0.2e1 - t735) * t425 + ((-t306 * t422 + t307 * t413) * pkin(1) - t499) * m(4) + t733 * (t700 + t568) + (m(5) * t244 + m(6) * t480 + t734) * t214 + t283 * t564 / 0.2e1 + t285 * t565 / 0.2e1 + (-t249 * t409 + t590) * t642 + t250 * t531 + Ifges(3,3) * qJDD(2) + t94 * t120 + t95 * t121 - t111 * t568 + t173 * t21 + t174 * t22 + t275 * t41 + t276 * t40 + Ifges(3,5) * t323 + Ifges(3,6) * t324 + t484 * t514; t40 * t654 + t41 * t655 + t432 - t250 * t535 + (-g(3) * t123 - t122 * t495) * t423 + (-g(3) * t122 + t123 * t495) * t414 - t176 * t562 - t248 * t533 - t80 * t92 - t79 * t93 - t111 * t567 + t375 * t10 + (t353 * ((-g(3) * t329 - t332 * t495) * t423 + (-g(3) * t332 + t329 * t495) * t414 + t434 + (-t613 + (-qJD(9) * t175 - t476 * t634 - t38) * t416) * pkin(2)) + (-pkin(12) * t353 + pkin(2) * cos(t515)) / pkin(12) * ((t175 * t416 - t177 * t407 + (t407 * t475 - t416 * t476) * mrSges(10,3)) * t400 * pkin(2) + t721)) * pkin(6) / t407 / pkin(2) - pkin(10) * t650 + (t115 * t375 + t229 * t563 - t59 * t79 - t60 * t80) * m(6) + (-t297 * t589 - (t244 * t296 + t246 * t295) * qJD(2)) * m(5) + t734 * t529 + t732 * (pkin(11) + t655) + t739 * (pkin(5) * t560 + t563) + ((-cos(qJ(3) + t512) * t650 + sin(t547) * t649) / cos(-qJ(7) - qJ(10) + t497) + (t229 * t412 + t421 * t480) * qJD(4) * m(6) + (t118 * t412 + t119 * t421 + (t244 * t421 - t246 * t412) * qJD(4)) * m(5)) * pkin(5); -t229 * (mrSges(6,1) * t158 + mrSges(6,2) * t157) + (Ifges(6,1) * t157 - t625) * t680 + t48 * t679 + (Ifges(6,5) * t157 - Ifges(6,6) * t158) * t678 - t59 * t92 + t60 * t93 + (t230 * t423 - t231 * t414) * g(3) + (t415 * t477 - t424 * t511) * g(2) + (t415 * t511 + t424 * t477) * g(1) + (t157 * t59 + t158 * t60) * mrSges(6,3) + t540 + (-Ifges(6,2) * t158 + t156 + t49) * t681 - t701;];
tau = t1(:);
