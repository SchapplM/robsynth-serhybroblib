% Calculate vector of inverse dynamics joint torques for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = palh1m1OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:37
% EndTime: 2020-04-15 19:31:10
% DurationCPUTime: 41.15s
% Computational Cost: add. (10975->947), mult. (25249->1336), div. (0->0), fcn. (20172->36), ass. (0->452)
t690 = pkin(2) * m(10);
t750 = mrSges(9,1) + t690;
t412 = sin(qJ(2));
t421 = cos(qJ(2));
t326 = -mrSges(3,1) * t412 - mrSges(3,2) * t421;
t401 = qJ(2) + qJ(7);
t384 = pkin(19) - t401;
t338 = pkin(4) * sin(t384);
t393 = t412 * pkin(1);
t386 = sin(t401);
t389 = cos(t401);
t358 = -qJ(10) + t384;
t339 = sin(t358);
t340 = cos(t358);
t511 = -t339 * mrSges(11,1) + t340 * mrSges(11,2);
t710 = mrSges(8,1) * t386 + mrSges(8,2) * t389 - t511;
t400 = qJ(2) + qJ(8);
t385 = sin(t400);
t388 = cos(t400);
t391 = qJ(9) + t400;
t361 = sin(t391);
t363 = cos(t391);
t562 = t361 * mrSges(10,1) + t363 * mrSges(10,2);
t713 = mrSges(9,2) * t388 + t750 * t385 - t562;
t749 = t713 - m(11) * (t338 - t393) - t326 + t710;
t742 = m(8) + m(4);
t748 = m(9) + m(3);
t624 = mrSges(10,2) * t361;
t747 = -mrSges(9,2) * t385 + t750 * t388 + t624;
t402 = qJ(2) + qJ(3);
t392 = qJ(4) + t402;
t362 = sin(t392);
t364 = cos(t392);
t409 = sin(qJ(5));
t626 = mrSges(6,2) * t409;
t746 = t362 * t626 + t364 * (pkin(11) * m(6) + mrSges(6,3));
t411 = sin(qJ(3));
t420 = cos(qJ(3));
t300 = t411 * t421 + t412 * t420;
t279 = t300 * qJD(1);
t303 = -t411 * t412 + t420 * t421;
t282 = t303 * qJD(1);
t410 = sin(qJ(4));
t419 = cos(qJ(4));
t508 = -t279 * t410 + t419 * t282;
t188 = Ifges(5,4) * t508;
t478 = t419 * t279 + t282 * t410;
t545 = qJD(3) + qJD(4);
t383 = qJD(2) + t545;
t612 = Ifges(5,5) * t383;
t100 = Ifges(5,1) * t478 + t188 + t612;
t559 = qJD(1) * t412;
t586 = qJD(1) * pkin(15);
t328 = pkin(1) * t559 - t586;
t228 = -pkin(5) * t282 + t328;
t635 = pkin(1) * qJD(2);
t537 = t420 * t635;
t333 = t410 * t537;
t399 = qJD(2) + qJD(3);
t539 = t411 * t635;
t468 = pkin(5) * t399 + t539;
t229 = t419 * t468 + t333;
t224 = -t383 * pkin(9) - t229;
t418 = cos(qJ(5));
t496 = mrSges(6,1) * t409 + mrSges(6,2) * t418;
t456 = t224 * t496;
t149 = t383 * t418 - t409 * t478;
t189 = qJD(5) - t508;
t150 = t383 * t409 + t418 * t478;
t606 = t150 * Ifges(6,4);
t46 = t149 * Ifges(6,2) + t189 * Ifges(6,6) + t606;
t148 = Ifges(6,4) * t149;
t47 = t150 * Ifges(6,1) + t189 * Ifges(6,5) + t148;
t592 = t418 * t47;
t604 = t229 * mrSges(5,3);
t662 = t409 / 0.2e1;
t745 = -t592 / 0.2e1 + t46 * t662 + t604 - t100 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t478 - t188 / 0.2e1 - t456 - t228 * mrSges(5,2) - t612 / 0.2e1;
t607 = Ifges(6,3) * t189;
t608 = Ifges(6,6) * t149;
t611 = Ifges(6,5) * t150;
t45 = t607 + t608 + t611;
t473 = -t410 * t468 + t419 * t537;
t603 = t473 * mrSges(5,3);
t609 = Ifges(5,6) * t383;
t617 = Ifges(5,4) * t478;
t225 = pkin(11) * t383 - t473;
t81 = -pkin(9) * t508 - pkin(11) * t478 + t228;
t58 = t225 * t418 + t409 * t81;
t723 = t58 * mrSges(6,2);
t57 = -t225 * t409 + t418 * t81;
t724 = t57 * mrSges(6,1);
t98 = Ifges(5,2) * t508 + t609 + t617;
t744 = t723 + t98 / 0.2e1 - t45 / 0.2e1 - t603 + t617 / 0.2e1 - t228 * mrSges(5,1) + t609 / 0.2e1 - t724;
t741 = t421 / 0.2e1;
t738 = mrSges(5,1) * t383 + mrSges(6,1) * t149 - mrSges(6,2) * t150 - mrSges(5,3) * t478;
t387 = sin(t402);
t390 = cos(t402);
t627 = mrSges(5,2) * t364;
t736 = mrSges(4,1) * t387 + t362 * mrSges(5,1) + mrSges(4,2) * t390 + t627;
t735 = -t364 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t362;
t734 = mrSges(11,1) * t340 + mrSges(11,2) * t339;
t631 = mrSges(6,1) * t418;
t733 = t626 - t631;
t413 = sin(qJ(1));
t422 = cos(qJ(1));
t705 = g(1) * t422 + g(2) * t413;
t407 = sin(qJ(7));
t416 = cos(qJ(7));
t299 = -t407 * t421 - t412 * t416;
t278 = t299 * qJD(1);
t558 = qJD(1) * t421;
t281 = -t407 * t559 + t416 * t558;
t403 = sin(pkin(19));
t404 = cos(pkin(19));
t158 = (-t278 * t404 - t281 * t403) * pkin(4) + t328;
t398 = qJD(2) + qJD(7);
t381 = qJD(10) + t398;
t665 = -t381 / 0.2e1;
t636 = sin(qJ(10));
t637 = cos(qJ(10));
t290 = t403 * t636 + t404 * t637;
t436 = -t403 * t637 + t404 * t636;
t127 = -t278 * t436 - t281 * t290;
t679 = -t127 / 0.2e1;
t730 = -t158 * mrSges(11,2) + Ifges(11,1) * t679 + Ifges(11,5) * t665;
t406 = sin(qJ(8));
t415 = cos(qJ(8));
t298 = -t406 * t421 - t412 * t415;
t277 = t298 * qJD(1);
t280 = -t406 * t559 + t415 * t558;
t405 = sin(qJ(9));
t414 = cos(qJ(9));
t480 = -t277 * t414 + t280 * t405;
t187 = Ifges(10,4) * t480;
t243 = -pkin(2) * t277 - t586;
t394 = qJDD(2) + qJDD(8);
t397 = qJD(2) + qJD(8);
t551 = qJD(9) * t405;
t255 = (-t394 * t414 + t397 * t551) * pkin(2);
t550 = qJD(9) * t414;
t256 = (-t394 * t405 - t397 * t550) * pkin(2);
t379 = qJDD(9) + t394;
t382 = qJD(9) + t397;
t479 = t277 * t405 + t280 * t414;
t548 = qJD(1) * qJD(2);
t312 = qJDD(1) * t421 - t412 * t548;
t313 = -qJDD(1) * t412 - t421 * t548;
t452 = t298 * qJD(8);
t159 = qJD(1) * t452 + t312 * t415 + t313 * t406;
t476 = t406 * t412 - t415 * t421;
t451 = t476 * qJD(8);
t160 = qJD(1) * t451 - t312 * t406 + t313 * t415;
t53 = qJD(9) * t480 - t159 * t414 - t160 * t405;
t54 = qJD(9) * t479 + t159 * t405 - t160 * t414;
t605 = t479 * Ifges(10,4);
t673 = -t479 / 0.2e1;
t97 = Ifges(10,2) * t480 + t382 * Ifges(10,6) - t605;
t99 = -Ifges(10,1) * t479 + t382 * Ifges(10,5) + t187;
t729 = t255 * mrSges(10,1) - t256 * mrSges(10,2) + Ifges(10,5) * t53 + Ifges(10,6) * t54 + Ifges(10,3) * t379 - (Ifges(10,5) * t480 + Ifges(10,6) * t479) * t382 / 0.2e1 - t243 * (-mrSges(10,1) * t479 + mrSges(10,2) * t480) + t97 * t673 + (Ifges(10,1) * t480 + t605) * t479 / 0.2e1 - (Ifges(10,2) * t479 + t187 + t99) * t480 / 0.2e1;
t509 = -t290 * t278 + t281 * t436;
t122 = Ifges(11,4) * t509;
t450 = t299 * qJD(7);
t161 = qJD(1) * t450 + t312 * t416 + t313 * t407;
t475 = t407 * t412 - t416 * t421;
t449 = t475 * qJD(7);
t162 = qJD(1) * t449 - t312 * t407 + t313 * t416;
t265 = t290 * qJD(10);
t266 = t436 * qJD(10);
t25 = -t161 * t290 - t162 * t436 - t265 * t278 + t266 * t281;
t26 = t161 * t436 - t162 * t290 + t265 * t281 + t266 * t278;
t395 = qJDD(2) + qJDD(7);
t372 = qJDD(10) + t395;
t540 = t407 * t635;
t652 = pkin(4) * t403;
t286 = t398 * t652 + t540;
t538 = t416 * t635;
t651 = pkin(4) * t404;
t287 = t398 * t651 + t538;
t129 = t286 * t290 + t287 * t436;
t583 = t129 * mrSges(11,3);
t634 = pkin(1) * qJD(7);
t530 = qJD(2) * t634;
t589 = pkin(1) * qJDD(2);
t291 = -t407 * t530 + t416 * t589;
t241 = t395 * t651 + t291;
t292 = t407 * t589 + t416 * t530;
t242 = t395 * t652 + t292;
t60 = -t241 * t436 - t242 * t290 - t265 * t287 + t266 * t286;
t61 = -t241 * t290 + t242 * t436 + t265 * t286 + t266 * t287;
t63 = t127 * Ifges(11,4) + Ifges(11,2) * t509 + t381 * Ifges(11,6);
t64 = t127 * Ifges(11,1) + t381 * Ifges(11,5) + t122;
t680 = -t509 / 0.2e1;
t727 = t61 * mrSges(11,1) - t60 * mrSges(11,2) + Ifges(11,5) * t25 + Ifges(11,6) * t26 + Ifges(11,3) * t372 + (t122 + t64) * t680 - (t158 * mrSges(11,1) + Ifges(11,4) * t679 + Ifges(11,2) * t680 + Ifges(11,6) * t665 + t583 - t63 / 0.2e1) * t127;
t408 = sin(qJ(6));
t663 = t408 / 0.2e1;
t661 = -t412 / 0.2e1;
t70 = -mrSges(11,1) * t509 + mrSges(11,2) * t127;
t719 = m(11) * t158 + t70;
t461 = pkin(1) * (t290 * t407 + t416 * t436);
t718 = -qJD(2) * t461 + (t265 * t403 + t266 * t404) * pkin(4);
t460 = pkin(1) * (-t290 * t416 + t407 * t436);
t717 = -qJD(2) * t460 + (-t265 * t404 + t266 * t403) * pkin(4);
t716 = -mrSges(4,1) * t282 - mrSges(8,1) * t278 + mrSges(4,2) * t279 + mrSges(8,2) * t281;
t629 = mrSges(10,1) * t363;
t316 = t413 * t629;
t715 = t413 * t747 - t316;
t319 = t422 * t629;
t714 = t422 * t747 - t319;
t565 = t734 * t413;
t625 = mrSges(8,2) * t386;
t712 = -t413 * t625 - t565;
t564 = t734 * t422;
t711 = -t422 * t625 - t564;
t561 = t364 * pkin(9) + t362 * pkin(11);
t515 = t390 * mrSges(4,1) - mrSges(4,2) * t387;
t709 = t364 * t733 + t735;
t215 = t300 * t419 + t303 * t410;
t552 = qJD(5) * t418;
t213 = t300 * t410 - t419 * t303;
t448 = t300 * qJD(3);
t222 = -qJD(2) * t300 - t448;
t447 = t303 * qJD(3);
t223 = qJD(2) * t303 + t447;
t86 = -qJD(4) * t213 + t222 * t410 + t223 * t419;
t466 = t215 * t552 + t409 * t86;
t553 = qJD(5) * t409;
t708 = -t57 * t552 - t58 * t553;
t107 = -mrSges(5,1) * t508 + mrSges(5,2) * t478;
t704 = m(5) * t228 + t107;
t499 = mrSges(3,1) * t421 - mrSges(3,2) * t412;
t618 = Ifges(3,4) * t421;
t703 = (-Ifges(3,1) * t412 - t618) * t741 - pkin(15) * t499;
t649 = pkin(5) * t387;
t656 = pkin(1) * t421;
t314 = -t649 - t656;
t531 = t362 * t631;
t630 = mrSges(8,1) * t389;
t646 = pkin(9) * t362;
t653 = pkin(4) * cos(t384);
t702 = -m(6) * (t314 - t646) + t531 + m(8) * t656 + t630 - m(11) * (-t653 - t656);
t168 = -mrSges(5,2) * t383 + mrSges(5,3) * t508;
t90 = -mrSges(6,2) * t189 + mrSges(6,3) * t149;
t91 = mrSges(6,1) * t189 - mrSges(6,3) * t150;
t700 = -t409 * t91 + t418 * t90 + t168;
t417 = cos(qJ(6));
t325 = -mrSges(7,1) * t417 + mrSges(7,2) * t408;
t445 = m(7) * pkin(14) + t325;
t657 = pkin(1) * t420;
t535 = qJD(3) * t657;
t294 = qJD(2) * t535 + t411 * t589;
t396 = qJDD(2) + qJDD(3);
t246 = pkin(5) * t396 + t294;
t556 = qJD(3) * t411;
t293 = (qJD(2) * t556 - qJDD(2) * t420) * pkin(1);
t118 = qJD(4) * t473 + t419 * t246 - t410 * t293;
t163 = -qJD(1) * t448 + t312 * t420 + t313 * t411;
t164 = qJD(1) * t447 + t312 * t411 - t313 * t420;
t56 = -qJD(4) * t478 + t163 * t419 - t164 * t410;
t380 = qJDD(4) + t396;
t55 = qJD(4) * t508 + t163 * t410 + t164 * t419;
t18 = -qJD(5) * t150 + t380 * t418 - t409 * t55;
t52 = qJDD(5) - t56;
t10 = -mrSges(6,2) * t52 + mrSges(6,3) * t18;
t17 = qJD(5) * t149 + t380 * t409 + t418 * t55;
t9 = mrSges(6,1) * t52 - mrSges(6,3) * t17;
t699 = t418 * t10 - t409 * t9 - t91 * t552 - t90 * t553;
t549 = qJDD(1) * pkin(15);
t259 = -pkin(1) * t313 - t549;
t115 = -pkin(5) * t163 + t259;
t11 = -pkin(9) * t56 - pkin(11) * t55 + t115;
t117 = qJD(4) * t229 + t410 * t246 + t419 * t293;
t111 = pkin(11) * t380 + t117;
t2 = qJD(5) * t57 + t11 * t409 + t111 * t418;
t3 = -qJD(5) * t58 + t11 * t418 - t111 * t409;
t698 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t484 = t409 * t58 + t418 * t57;
t697 = m(6) * t484 + t409 * t90 + t418 * t91;
t109 = pkin(9) * t478 - pkin(11) * t508;
t359 = pkin(5) * t390;
t696 = -m(6) * (t359 + t561) - t515 + t709;
t641 = t3 * t409;
t432 = -qJD(5) * t484 - t641;
t642 = t2 * t418;
t695 = m(6) * (t432 + t642) + t699;
t694 = mrSges(2,2) - mrSges(9,3) - mrSges(7,3) - mrSges(3,3) - mrSges(10,3) - mrSges(8,3) - mrSges(5,3) - mrSges(4,3) - mrSges(11,3);
t560 = t359 - t393;
t305 = pkin(15) + t560;
t360 = -t393 + pkin(15);
t693 = m(6) * (t305 + t561) + mrSges(2,1) + m(5) * t305 + t515 - t445 - t735 + t742 * t360 + (m(11) + m(10) + t748) * pkin(15) - t749;
t689 = m(6) * pkin(9);
t688 = t17 / 0.2e1;
t687 = t18 / 0.2e1;
t684 = t52 / 0.2e1;
t678 = -t149 / 0.2e1;
t677 = -t150 / 0.2e1;
t676 = t150 / 0.2e1;
t675 = -t189 / 0.2e1;
t669 = t279 / 0.2e1;
t668 = t280 / 0.2e1;
t667 = t281 / 0.2e1;
t660 = pkin(1) * t407;
t659 = pkin(1) * t411;
t658 = pkin(1) * t416;
t655 = pkin(2) * (-mrSges(10,1) * t480 - mrSges(10,2) * t479);
t650 = pkin(5) * t279;
t648 = pkin(5) * t410;
t647 = pkin(5) * t419;
t638 = -qJD(6) / 0.2e1;
t623 = mrSges(6,3) * t409;
t622 = mrSges(6,3) * t418;
t621 = mrSges(8,3) * t281;
t620 = mrSges(10,3) * t397;
t619 = Ifges(3,4) * t412;
t616 = Ifges(6,4) * t409;
t615 = Ifges(6,4) * t418;
t614 = Ifges(7,4) * t408;
t613 = Ifges(7,4) * t417;
t610 = Ifges(7,2) * t417;
t602 = t279 * Ifges(4,4);
t601 = t280 * Ifges(9,4);
t600 = t281 * Ifges(8,4);
t598 = t405 * (-mrSges(10,2) * t379 + mrSges(10,3) * t54);
t597 = t408 * Ifges(7,1);
t585 = t509 * mrSges(11,3);
t128 = t286 * t436 - t287 * t290;
t584 = t128 * mrSges(11,3);
t582 = t215 * t409;
t581 = t215 * t418;
t580 = (-mrSges(8,2) * t398 + mrSges(8,3) * t278) * t416;
t488 = t610 + t614;
t573 = t408 * (Ifges(7,6) * qJD(6) + qJD(1) * t488);
t572 = t409 * t413;
t571 = t409 * t422;
t570 = t410 * t420;
t569 = t413 * t418;
t371 = qJD(1) * t613;
t568 = t417 * (Ifges(7,5) * qJD(6) + qJD(1) * t597 + t371);
t567 = t418 * t422;
t566 = t419 * t420;
t557 = qJD(2) * t421;
t555 = qJD(4) * t410;
t554 = qJD(4) * t419;
t547 = qJD(1) * qJD(6);
t544 = t243 * t690;
t543 = pkin(2) * t620;
t542 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t52;
t541 = pkin(1) * t566;
t536 = pkin(1) * t556;
t375 = pkin(1) * t557;
t523 = t592 / 0.2e1;
t518 = -t553 / 0.2e1;
t202 = -pkin(5) * t222 + t375;
t517 = -t548 / 0.2e1;
t513 = t405 * t543;
t512 = t414 * t543;
t507 = mrSges(4,3) * t539;
t506 = mrSges(4,3) * t537;
t505 = mrSges(8,3) * t540;
t504 = mrSges(8,3) * t538;
t244 = -pkin(5) * t303 - t360;
t366 = pkin(5) + t659;
t269 = t410 * t366 - t541;
t501 = t746 * t413;
t500 = t746 * t422;
t495 = mrSges(7,1) * t408 + mrSges(7,2) * t417;
t492 = t421 * Ifges(3,1) - t619;
t491 = Ifges(6,1) * t418 - t616;
t490 = -t412 * Ifges(3,2) + t618;
t489 = -Ifges(6,2) * t409 + t615;
t487 = -Ifges(3,5) * t412 - Ifges(3,6) * t421;
t486 = Ifges(6,5) * t418 - Ifges(6,6) * t409;
t485 = Ifges(7,5) * t417 - Ifges(7,6) * t408;
t483 = -t409 * t57 + t418 * t58;
t212 = -t298 * t414 - t405 * t476;
t477 = t298 * t405 - t414 * t476;
t470 = m(11) * t653 + t630;
t268 = pkin(1) * t570 + t366 * t419;
t465 = t215 * t553 - t418 * t86;
t463 = pkin(14) * t495;
t459 = t149 * t489;
t458 = t150 * t491;
t457 = t189 * t486;
t455 = t408 * (Ifges(7,1) * t417 - t614);
t454 = t412 * (-Ifges(3,2) * t421 - t619);
t176 = (-t278 * t403 + t281 * t404) * pkin(4);
t89 = t109 + t650;
t431 = m(6) * (-t646 - t649) - t531;
t430 = t627 + (mrSges(5,1) + t631 + t689) * t362;
t112 = -t380 * pkin(9) - t118;
t6 = t17 * Ifges(6,4) + t18 * Ifges(6,2) + t52 * Ifges(6,6);
t7 = t17 * Ifges(6,1) + t18 * Ifges(6,4) + t52 * Ifges(6,5);
t427 = -t117 * mrSges(5,2) + t2 * t622 + t418 * t6 / 0.2e1 + t7 * t662 + t112 * t733 + t118 * mrSges(5,1) + (Ifges(6,1) * t409 + t615) * t688 + (Ifges(6,2) * t418 + t616) * t687 + (Ifges(6,5) * t409 + Ifges(6,6) * t418) * t684 + Ifges(5,3) * t380 + t46 * t518 + Ifges(5,6) * t56 + Ifges(5,5) * t55 + (t456 + t523) * qJD(5) + (t459 + t458 + t457) * qJD(5) / 0.2e1;
t169 = mrSges(10,1) * t382 + mrSges(10,3) * t479;
t181 = t277 * Ifges(9,2) + t397 * Ifges(9,6) + t601;
t262 = Ifges(9,4) * t277;
t184 = t280 * Ifges(9,1) + t397 * Ifges(9,5) + t262;
t426 = t479 * t513 + (mrSges(9,1) * t280 + mrSges(9,2) * t277) * t586 - t280 * t655 + t181 * t668 + pkin(2) * t169 * t551 - t397 * (Ifges(9,5) * t277 - Ifges(9,6) * t280) / 0.2e1 + Ifges(9,3) * t394 + Ifges(9,5) * t159 + Ifges(9,6) * t160 - t280 * (Ifges(9,1) * t277 - t601) / 0.2e1 + (-t255 * t414 - t256 * t405) * t690 - t280 * t544 - (-Ifges(9,2) * t280 + t184 + t262) * t277 / 0.2e1 + t729;
t182 = t278 * Ifges(8,2) + t398 * Ifges(8,6) + t600;
t263 = Ifges(8,4) * t278;
t185 = t281 * Ifges(8,1) + t398 * Ifges(8,5) + t263;
t425 = -t292 * mrSges(8,2) - t328 * (mrSges(8,1) * t281 + mrSges(8,2) * t278) + t182 * t667 - t281 * (Ifges(8,1) * t278 - t600) / 0.2e1 + Ifges(8,6) * t162 + Ifges(8,5) * t161 - t398 * (Ifges(8,5) * t278 - Ifges(8,6) * t281) / 0.2e1 + t278 * t504 + t291 * mrSges(8,1) + Ifges(8,3) * t395 - (-Ifges(8,2) * t281 + t185 + t263) * t278 / 0.2e1 + (t584 + t730) * t509 + t727;
t183 = t282 * Ifges(4,2) + t399 * Ifges(4,6) + t602;
t264 = Ifges(4,4) * t282;
t186 = t279 * Ifges(4,1) + t399 * Ifges(4,5) + t264;
t424 = -t3 * t623 + t183 * t669 - t328 * (mrSges(4,1) * t279 + mrSges(4,2) * t282) - t279 * t506 + t282 * t507 - t399 * (Ifges(4,5) * t282 - Ifges(4,6) * t279) / 0.2e1 + Ifges(4,3) * t396 + t294 * mrSges(4,1) - t293 * mrSges(4,2) + Ifges(4,6) * t163 + Ifges(4,5) * t164 - t279 * (Ifges(4,1) * t282 - t602) / 0.2e1 + t427 + t708 * mrSges(6,3) - (-Ifges(4,2) * t279 + t186 + t264) * t282 / 0.2e1 + (Ifges(6,5) * t677 + Ifges(6,6) * t678 + Ifges(6,3) * t675 + t744) * t478 + (t486 * t675 + t489 * t678 + t491 * t677 + t57 * t622 + t58 * t623 + t745) * t508;
t367 = -pkin(9) - t647;
t327 = t651 + t658;
t324 = t652 + t660;
t311 = qJDD(1) * t408 + t417 * t547;
t310 = qJDD(1) * t417 - t408 * t547;
t276 = Ifges(3,5) * qJD(2) + qJD(1) * t492;
t274 = Ifges(3,6) * qJD(2) + qJD(1) * t490;
t261 = t419 * t539 + t333;
t260 = (t410 * t411 - t566) * t635;
t254 = t364 * t567 + t572;
t253 = -t364 * t571 + t569;
t252 = -t364 * t569 + t571;
t251 = t364 * t572 + t567;
t235 = -mrSges(4,2) * t399 + mrSges(4,3) * t282;
t234 = mrSges(8,1) * t398 - t621;
t233 = mrSges(4,1) * t399 - mrSges(4,3) * t279;
t221 = qJD(2) * t475 + t449;
t220 = qJD(2) * t299 + t450;
t219 = qJD(2) * t476 + t451;
t218 = qJD(2) * t298 + t452;
t204 = (-t290 * t403 - t404 * t436) * pkin(4);
t203 = (-t290 * t404 + t403 * t436) * pkin(4);
t167 = -mrSges(10,2) * t382 + mrSges(10,3) * t480;
t166 = -t290 * t324 - t327 * t436;
t165 = -t290 * t327 + t324 * t436;
t139 = -pkin(2) * t160 - t549;
t114 = mrSges(11,1) * t381 - mrSges(11,3) * t127;
t113 = -mrSges(11,2) * t381 + t585;
t93 = qJD(7) * t461 + t265 * t324 + t266 * t327;
t92 = qJD(7) * t460 - t265 * t327 + t266 * t324;
t87 = qJD(4) * t215 - t419 * t222 + t223 * t410;
t85 = qJD(9) * t477 + t218 * t405 - t219 * t414;
t84 = qJD(9) * t212 - t218 * t414 - t219 * t405;
t80 = t109 * t409 + t229 * t418;
t79 = t109 * t418 - t229 * t409;
t73 = t261 * t418 + t409 * t89;
t72 = -t261 * t409 + t418 * t89;
t71 = (-t161 * t403 - t162 * t404) * pkin(4) + t259;
t43 = t220 * t436 - t221 * t290 - t265 * t475 + t266 * t299;
t42 = -t220 * t290 - t221 * t436 - t265 * t299 - t266 * t475;
t39 = -mrSges(5,2) * t380 + mrSges(5,3) * t56;
t38 = mrSges(5,1) * t380 - mrSges(5,3) * t55;
t36 = mrSges(10,1) * t379 - mrSges(10,3) * t53;
t20 = -mrSges(11,2) * t372 + mrSges(11,3) * t26;
t19 = mrSges(11,1) * t372 - mrSges(11,3) * t25;
t8 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t1 = [(-mrSges(10,1) * t139 + mrSges(10,3) * t256 + Ifges(10,4) * t53 + Ifges(10,2) * t54 + Ifges(10,6) * t379) * t212 + t87 * t724 + (m(11) * t71 - mrSges(11,1) * t26 + mrSges(11,2) * t25) * ((-t299 * t404 + t403 * t475) * pkin(4) - t360) + (mrSges(11,2) * t71 - mrSges(11,3) * t61 + Ifges(11,1) * t25 + Ifges(11,4) * t26 + Ifges(11,5) * t372) * (t290 * t475 - t299 * t436) + (-mrSges(11,1) * t71 + mrSges(11,3) * t60 + Ifges(11,4) * t25 + Ifges(11,2) * t26 + Ifges(11,6) * t372) * (-t290 * t299 - t436 * t475) - (mrSges(8,2) * t259 - mrSges(8,3) * t291 + Ifges(8,1) * t161 + Ifges(8,4) * t162 + Ifges(8,5) * t395) * t475 - (-mrSges(9,2) * t549 + Ifges(9,1) * t159 + Ifges(9,4) * t160 + Ifges(9,5) * t394) * t476 + (Ifges(3,4) * t312 + Ifges(3,2) * t313) * t661 + t509 * (Ifges(11,4) * t42 + Ifges(11,2) * t43) / 0.2e1 + t480 * (Ifges(10,4) * t84 + Ifges(10,2) * t85) / 0.2e1 + t508 * (Ifges(5,4) * t86 - Ifges(5,2) * t87) / 0.2e1 + t704 * t202 + t703 * t548 + (mrSges(4,1) * t163 + mrSges(8,1) * t162 - mrSges(4,2) * t164 - mrSges(8,2) * t161) * t360 + (Ifges(3,1) * t312 + Ifges(3,4) * t313) * t741 + (t568 / 0.2e1 + t485 * qJD(6) / 0.2e1) * qJD(6) + (-t252 * mrSges(6,1) - t251 * mrSges(6,2) + t413 * t693 + t422 * t694) * g(1) + (-t254 * mrSges(6,1) - t253 * mrSges(6,2) + t413 * t694 - t422 * t693) * g(2) + (t276 * t661 + t487 * qJD(2) / 0.2e1) * qJD(2) + t417 * (Ifges(7,4) * t311 + Ifges(7,2) * t310) / 0.2e1 + t742 * (-t259 * t360 + t328 * t375) + t311 * (t597 + t613) / 0.2e1 + (mrSges(4,2) * t259 - mrSges(4,3) * t294 + Ifges(4,1) * t164 + Ifges(4,4) * t163 + Ifges(4,5) * t396) * t300 - t43 * t583 + (-t2 * t582 - t3 * t581 + t465 * t57 - t466 * t58) * mrSges(6,3) - t6 * t582 / 0.2e1 + (-mrSges(4,1) * t259 + mrSges(4,3) * t293 + Ifges(4,4) * t164 + Ifges(4,2) * t163 + Ifges(4,6) * t396) * t303 - t274 * t557 / 0.2e1 + (pkin(15) ^ 2 * t748 + t445 * pkin(14) + Ifges(2,3)) * qJDD(1) - t42 * t584 - (-mrSges(9,1) * t219 + mrSges(9,2) * t218) * t586 - t85 * t513 + t87 * t603 + t478 * (Ifges(5,1) * t86 - Ifges(5,4) * t87) / 0.2e1 + (Ifges(8,1) * t220 + Ifges(8,4) * t221) * t667 + (Ifges(9,1) * t218 + Ifges(9,4) * t219) * t668 + (Ifges(4,1) * t223 + Ifges(4,4) * t222) * t669 + (Ifges(10,1) * t84 + Ifges(10,4) * t85) * t673 - t219 * t544 + t463 * t547 + t326 * t549 - t220 * t504 + (m(5) * t115 - mrSges(5,1) * t56 + mrSges(5,2) * t55) * t244 + (Ifges(7,1) * t311 + Ifges(7,4) * t310) * t663 + (t115 * mrSges(5,2) - t118 * mrSges(5,3) + Ifges(5,1) * t55 + Ifges(5,4) * t56 + Ifges(5,5) * t380 + t112 * t496 + t47 * t518 + t486 * t684 + t489 * t687 + t491 * t688) * t215 + t573 * t638 - t223 * t507 - t86 * t604 + t697 * (pkin(9) * t87 - pkin(11) * t86 + t202) + (t542 / 0.2e1 - Ifges(5,6) * t380 - Ifges(5,4) * t55 - Ifges(5,2) * t56 + t115 * mrSges(5,1) + Ifges(6,3) * t684 + Ifges(6,6) * t687 + Ifges(6,5) * t688 - t117 * mrSges(5,3) + t698) * t213 - t219 * t655 + (t90 * t552 + m(6) * (qJD(5) * t483 + t2 * t409 + t3 * t418) + t418 * t9 + t409 * t10 - t91 * t553) * (pkin(9) * t213 - pkin(11) * t215 + t244) - (mrSges(10,2) * t139 - mrSges(10,3) * t255 + Ifges(10,1) * t53 + Ifges(10,4) * t54 + Ifges(10,5) * t379) * t477 + (m(10) * t139 - mrSges(10,1) * t54 + mrSges(10,2) * t53) * (-pkin(2) * t298 - pkin(15)) + (t417 * (-Ifges(7,2) * t408 + t613) + t455) * t547 / 0.2e1 - t222 * t506 + t7 * t581 / 0.2e1 + (-mrSges(4,1) * t222 - mrSges(8,1) * t221 + mrSges(4,2) * t223 + mrSges(8,2) * t220) * t328 + (mrSges(9,1) * t549 + Ifges(9,4) * t159 + Ifges(9,2) * t160 + Ifges(9,6) * t394) * t298 + t221 * t505 + t84 * t512 + (-mrSges(8,1) * t259 + mrSges(8,3) * t292 + Ifges(8,4) * t161 + Ifges(8,2) * t162 + Ifges(8,6) * t395) * t299 + t397 * (Ifges(9,5) * t218 + Ifges(9,6) * t219) / 0.2e1 + t398 * (Ifges(8,5) * t220 + Ifges(8,6) * t221) / 0.2e1 + t399 * (Ifges(4,5) * t223 + Ifges(4,6) * t222) / 0.2e1 + t382 * (Ifges(10,5) * t84 + Ifges(10,6) * t85) / 0.2e1 + t383 * (Ifges(5,5) * t86 - Ifges(5,6) * t87) / 0.2e1 + t381 * (Ifges(11,5) * t42 + Ifges(11,6) * t43) / 0.2e1 + pkin(14) * (-mrSges(7,1) * t310 + mrSges(7,2) * t311) + t282 * (Ifges(4,4) * t223 + Ifges(4,2) * t222) / 0.2e1 + t277 * (Ifges(9,4) * t218 + Ifges(9,2) * t219) / 0.2e1 + t278 * (Ifges(8,4) * t220 + Ifges(8,2) * t221) / 0.2e1 + t243 * (-mrSges(10,1) * t85 + mrSges(10,2) * t84) + t228 * (mrSges(5,1) * t87 + mrSges(5,2) * t86) + t223 * t186 / 0.2e1 + t312 * t492 / 0.2e1 + t222 * t183 / 0.2e1 + t218 * t184 / 0.2e1 + t219 * t181 / 0.2e1 + t310 * t488 / 0.2e1 + t220 * t185 / 0.2e1 + t221 * t182 / 0.2e1 + t313 * t490 / 0.2e1 + t158 * (-mrSges(11,1) * t43 + mrSges(11,2) * t42) + t224 * (mrSges(6,1) * t466 - mrSges(6,2) * t465) + t149 * (-Ifges(6,4) * t465 - Ifges(6,2) * t466 + Ifges(6,6) * t87) / 0.2e1 + t189 * (-Ifges(6,5) * t465 - Ifges(6,6) * t466 + Ifges(6,3) * t87) / 0.2e1 + t127 * (Ifges(11,1) * t42 + Ifges(11,4) * t43) / 0.2e1 + (mrSges(3,1) * t313 + mrSges(9,1) * t160 - mrSges(3,2) * t312 - mrSges(9,2) * t159) * pkin(15) + t85 * t97 / 0.2e1 - t87 * t98 / 0.2e1 + t84 * t99 / 0.2e1 + t86 * t100 / 0.2e1 + t87 * t45 / 0.2e1 + t43 * t63 / 0.2e1 + t42 * t64 / 0.2e1 + (-Ifges(6,1) * t465 - Ifges(6,4) * t466 + Ifges(6,5) * t87) * t676 - t466 * t46 / 0.2e1 + t716 * t375 + t719 * (t375 + (-t220 * t403 - t221 * t404) * pkin(4)) - t87 * t723 + t454 * t517 + t86 * t523 + (Ifges(3,5) * t421 + 0.2e1 * Ifges(3,6) * t661) * qJDD(2) + (0.2e1 * Ifges(7,5) * t663 + Ifges(7,6) * t417) * qJDD(6); m(11) * (t128 * t93 - t129 * t92 + t165 * t61 + t166 * t60) + (m(5) * t229 - m(6) * t224 + t738) * (-t366 * t555 - t410 * t536 + t541 * t545) + ((-t742 * t328 - t697 - t704 - t716 - t719) * pkin(1) + t274 / 0.2e1) * t558 + (-m(5) * t473 + m(6) * t483 + t700) * (t366 * t554 + (t420 * t555 + (t411 * t419 + t570) * qJD(3)) * pkin(1)) + (t454 / 0.2e1 - t703) * qJD(1) ^ 2 + t695 * (pkin(11) + t269) + t705 * (m(4) * t656 - m(5) * t314 + t499 + t736) + (-t234 * t407 + t580) * t634 + (m(8) * (t291 * t416 + t292 * t407) + m(4) * (-t293 * t420 + t294 * t411)) * pkin(1) - t697 * t89 - t704 * t650 - t719 * t176 + t276 * t559 / 0.2e1 + t426 + m(5) * (t117 * t269 + t118 * t268) - t480 * t512 + (mrSges(8,1) * t395 - mrSges(8,3) * t161) * t658 + (mrSges(4,1) * t396 - mrSges(4,3) * t164) * t659 + (-mrSges(8,2) * t395 + mrSges(8,3) * t162) * t660 + (-t167 * t550 - t36 * t414 - t598) * pkin(2) + (m(6) * t112 + t8) * (-pkin(9) - t268) + t281 * t505 + t424 + Ifges(3,5) * t312 + Ifges(3,6) * t313 + t268 * t38 + t269 * t39 + t165 * t19 + t166 * t20 + Ifges(3,3) * qJDD(2) + t92 * t113 + t93 * t114 + (t422 * t702 - t500 + t711 + t714) * g(1) + (t413 * t702 - t501 + t712 + t715) * g(2) - (-mrSges(4,2) * t396 + mrSges(4,3) * t163) * t657 + t233 * t535 + t235 * t536 + (-m(5) * t560 + (m(6) + t742) * t393 + t696 + t749) * g(3) + t487 * t517 + t425; -t107 * t650 - t233 * t537 - t235 * t539 + t38 * t647 + t39 * t648 + t424 + t367 * t8 - g(1) * (t422 * t431 + t500) - g(2) * (t413 * t431 + t501) - t261 * t168 - t73 * t90 - t72 * t91 + (t112 * t367 + (t224 * t410 + t419 * t483) * qJD(4) * pkin(5) - t224 * t260 - t57 * t72 - t58 * t73) * m(6) + (-t228 * t650 + t229 * t260 + t473 * t261 + (t117 * t410 + t118 * t419 + (-t229 * t410 - t419 * t473) * qJD(4)) * pkin(5)) * m(5) + t700 * pkin(5) * t554 + (-m(5) * t359 + t696) * g(3) + t695 * (pkin(11) + t648) + t738 * (-pkin(5) * t555 + t260) + (m(5) * t649 + t736) * t705; (-t607 / 0.2e1 - t608 / 0.2e1 - t611 / 0.2e1 + t744) * t478 + (m(6) * (-t641 + t642 + t708) + t699) * pkin(11) + (-m(6) * t561 + t709) * g(3) - t738 * t473 + (-t457 / 0.2e1 - t459 / 0.2e1 - t458 / 0.2e1 + t484 * mrSges(6,3) + t745) * t508 - t112 * t689 - t229 * t168 - t80 * t90 - t79 * t91 - m(6) * (-t224 * t473 + t57 * t79 + t58 * t80) - pkin(9) * t8 + t427 + (t422 * t430 - t500) * g(1) + t432 * mrSges(6,3) + (t413 * t430 - t501) * g(2); -t224 * (mrSges(6,1) * t150 + mrSges(6,2) * t149) + (Ifges(6,1) * t149 - t606) * t677 + t46 * t676 + (Ifges(6,5) * t149 - Ifges(6,6) * t150) * t675 - t57 * t90 + t58 * t91 - g(1) * (mrSges(6,1) * t253 - mrSges(6,2) * t254) - g(2) * (-mrSges(6,1) * t251 + mrSges(6,2) * t252) + g(3) * t496 * t362 + (t149 * t57 + t150 * t58) * mrSges(6,3) + t542 + (-Ifges(6,2) * t150 + t148 + t47) * t678 + t698; Ifges(7,5) * t311 + Ifges(7,6) * t310 + Ifges(7,3) * qJDD(6) + g(3) * t325 + (-t568 / 0.2e1 + t573 / 0.2e1 - t417 * t371 / 0.2e1 + t485 * t638 + (-t463 - t455 / 0.2e1 + t610 * t663) * qJD(1)) * qJD(1) + t705 * t495; t718 * t114 + t717 * t113 + t710 * g(3) + (-t580 + (t234 + t621) * t407) * t635 + (t422 * t470 + t711) * g(1) + (t413 * t470 + t712) * g(2) + t203 * t19 + t204 * t20 - t176 * t70 + t425 + (-t338 * g(3) + t128 * t718 - t129 * t717 - t158 * t176 + t203 * t61 + t204 * t60) * m(11); t713 * g(3) + t715 * g(2) + t714 * g(1) + t426 + (-t598 + (-qJD(9) * t167 - t480 * t620 - t36) * t414) * pkin(2); -g(3) * t562 - g(1) * (-t422 * t624 + t319) - g(2) * (-t413 * t624 + t316) + (t167 * t414 - t169 * t405 + (t405 * t479 - t414 * t480) * mrSges(10,3)) * t397 * pkin(2) + t729; -t129 * t114 - g(3) * t511 - g(1) * t564 - g(2) * t565 + (-t113 + t585) * t128 + t730 * t509 + t727; 0; 0; 0;];
tau = t1;
