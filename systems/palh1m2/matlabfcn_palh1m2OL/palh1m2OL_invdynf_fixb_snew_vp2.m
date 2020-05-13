% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x11]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = palh1m2OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: qJD has to be [13x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [13 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: qJDD has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_invdynf_fixb_snew_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:28:51
% EndTime: 2020-05-02 21:29:17
% DurationCPUTime: 27.62s
% Computational Cost: add. (8255->759), mult. (11902->1014), div. (0->0), fcn. (8689->22), ass. (0->420)
t378 = cos(qJ(2));
t364 = sin(qJ(9));
t601 = cos(qJ(9));
t302 = -t364 * mrSges(10,1) - mrSges(10,2) * t601;
t307 = mrSges(10,1) * t601 - t364 * mrSges(10,2);
t365 = sin(qJ(8));
t373 = cos(qJ(8));
t415 = t365 * t302 + t307 * t373;
t517 = t378 * t415;
t371 = sin(qJ(2));
t416 = t302 * t373 - t365 * t307;
t655 = t416 * t371;
t108 = t517 + t655;
t516 = t378 * t416;
t109 = t415 * t371 - t516;
t358 = sin(pkin(19));
t605 = pkin(4) * m(11);
t570 = t358 * t605 - mrSges(8,2);
t359 = cos(pkin(19));
t293 = mrSges(11,1) * t359 + mrSges(11,2) * t358;
t361 = sin(qJ(10));
t551 = t293 * t361;
t292 = mrSges(11,1) * t358 - mrSges(11,2) * t359;
t362 = cos(qJ(10));
t552 = t292 * t362;
t629 = t552 - t551;
t147 = -t629 + t570;
t366 = sin(qJ(7));
t136 = t147 * t366;
t550 = t293 * t362;
t553 = t292 * t361;
t192 = t550 + t553;
t571 = t359 * t605 + mrSges(8,1);
t148 = t192 - t571;
t374 = cos(qJ(7));
t138 = t148 * t374;
t691 = t138 - t136;
t580 = t378 * t691;
t670 = t366 * t148;
t67 = t147 * t374 + t670;
t25 = t67 * t371 - t580;
t367 = sin(qJ(6));
t375 = cos(qJ(6));
t303 = t367 * mrSges(7,1) + t375 * mrSges(7,2);
t387 = qJD(9) ^ 2;
t388 = qJD(8) ^ 2;
t389 = qJD(7) ^ 2;
t390 = qJD(6) ^ 2;
t392 = qJD(4) ^ 2;
t393 = qJD(3) ^ 2;
t394 = qJD(2) ^ 2;
t395 = qJD(1) ^ 2;
t602 = cos(qJ(5));
t341 = mrSges(6,1) * t602;
t368 = sin(qJ(5));
t588 = mrSges(6,2) * t368;
t305 = -t341 + t588;
t385 = pkin(9) * m(6);
t262 = -mrSges(5,1) + t305 - t385;
t369 = sin(qJ(4));
t536 = t369 * t262;
t608 = pkin(11) * m(6);
t328 = -mrSges(5,2) + mrSges(6,3) + t608;
t376 = cos(qJ(4));
t546 = t328 * t376;
t203 = t546 + t536;
t199 = mrSges(4,2) - t203;
t377 = cos(qJ(3));
t683 = t199 * t377;
t382 = m(5) + m(6);
t346 = pkin(5) * t382;
t447 = mrSges(4,1) + t346;
t547 = t328 * t369;
t557 = t262 * t376;
t628 = t557 - t547;
t167 = t447 - t628;
t370 = sin(qJ(3));
t685 = t167 * t370;
t94 = t683 + t685;
t579 = t378 * t94;
t684 = t199 * t370;
t690 = t167 * t377 - t684;
t43 = t371 * t690 + t579;
t572 = t690 * t378;
t45 = -t94 * t371 + t572;
t540 = t366 * t629;
t97 = t192 * t374 + t540;
t577 = t378 * t97;
t539 = t366 * t192;
t660 = t629 * t374 - t539;
t46 = t371 * t660 + t577;
t669 = t370 * t628;
t671 = t203 * t377;
t111 = -t671 - t669;
t519 = t378 * t111;
t672 = t203 * t370;
t678 = t377 * t628 - t672;
t54 = t678 * t371 - t519;
t627 = -t375 * mrSges(7,1) + t367 * mrSges(7,2);
t263 = m(10) * pkin(2) + mrSges(9,1) - t307;
t543 = t365 * t263;
t290 = -mrSges(9,2) - t302;
t554 = t290 * t373;
t171 = -t554 + t543;
t155 = t171 * t378;
t169 = -t263 * t373 - t365 * t290;
t661 = t169 * t371 - t155;
t518 = t378 * t678;
t688 = t111 * t371 + t518;
t581 = t378 * t67;
t689 = t371 * t691 + t581;
t352 = m(4) + t382;
t380 = m(8) + m(11);
t653 = t380 + t352;
t14 = -t653 * pkin(1) - mrSges(3,1) + t169 + t691 - t94;
t20 = t690 + t67 - mrSges(3,2) - t171;
t702 = -t14 * t371 - t378 * t20;
t381 = m(9) + m(10);
t313 = m(3) + t381 + t653;
t4 = -pkin(15) * t313 + t702;
t705 = pkin(14) * m(7) - mrSges(2,1) + t4 + t627;
t562 = t171 * t371;
t563 = t169 * t378;
t88 = t562 + t563;
t582 = t378 * t14;
t9 = t20 * t371 - t582;
t706 = -t9 * qJDD(2) + t46 * qJDD(10) - t43 * qJDD(3) + t54 * qJDD(4) - t303 * qJDD(6) - t25 * qJDD(7) + t88 * qJDD(8) + t108 * qJDD(9) - t109 * t387 - t661 * t388 - t689 * t389 + t627 * t390 + t688 * t392 - t45 * t393 + t702 * t394 + t705 * t395;
t357 = t378 ^ 2;
t704 = (-t357 + 0.1e1) * t395;
t344 = pkin(1) * t380;
t697 = t344 - t691;
t446 = 0.2e1 * t447;
t294 = 0.2e1 * t547;
t469 = 0.2e1 * t557;
t496 = t294 - t469;
t145 = t446 + t496;
t133 = t145 * t377;
t85 = t133 - 0.2e1 * t684;
t75 = t85 * t371;
t696 = (t75 + 0.2e1 * t579) * qJD(3);
t632 = 0.2e1 * t552 - 0.2e1 * t551;
t633 = -0.2e1 * t570 + t632;
t668 = t633 * t374;
t676 = t668 - 0.2e1 * t670;
t62 = t676 * t371;
t468 = 0.2e1 * t554;
t654 = t468 - 0.2e1 * t543;
t659 = t632 * t374 - 0.2e1 * t539;
t40 = t371 * t659 + 0.2e1 * t577;
t681 = t40 * qJD(10);
t695 = (0.2e1 * t580 + t62) * qJD(7) - (t654 * t371 - 0.2e1 * t563) * qJD(8) + t681;
t306 = t368 * mrSges(6,1) + mrSges(6,2) * t602;
t488 = t306 * qJD(5);
t520 = t377 * t376;
t531 = t370 * t369;
t270 = -t520 + t531;
t535 = t369 * t377;
t271 = t370 * t376 + t535;
t630 = t270 * t378 + t371 * t271;
t438 = t630 * t488;
t121 = -0.2e1 * t438;
t611 = 0.2e1 * t378;
t183 = 0.2e1 * t672;
t465 = -0.2e1 * t547;
t197 = t469 + t465;
t651 = -t197 * t377 + t183;
t667 = t651 * t371;
t686 = t121 + (t111 * t611 + t667) * qJD(4);
t544 = t352 * t370;
t412 = pkin(1) * t544 + t167;
t578 = t378 * t660;
t47 = t97 * t371 - t578;
t666 = t346 - t628;
t680 = t666 * t370 - t671;
t372 = sin(qJ(1));
t379 = cos(qJ(1));
t421 = g(1) * t379 + t372 * g(2);
t175 = 0.2e1 * t550 + 0.2e1 * t553;
t129 = -0.2e1 * t571 + t175;
t634 = 0.2e1 * t546 + 0.2e1 * t536;
t176 = -0.2e1 * mrSges(4,2) + t634;
t677 = t129 * t366 + t176 * t370 - 0.2e1 * mrSges(3,2) + t654 - t668;
t664 = qJD(1) * t305;
t662 = t129 * t374 - 0.2e1 * t136;
t436 = qJD(5) * t664;
t427 = -0.2e1 * t436;
t355 = t376 ^ 2;
t558 = t262 * t355;
t407 = -0.2e1 * t558 + t262;
t420 = g(1) * t372 - g(2) * t379;
t311 = t328 * t355;
t441 = t376 * t536;
t652 = t441 + t311;
t90 = t175 * t374 + 0.2e1 * t540;
t610 = 2 * qJD(2);
t650 = qJD(3) * t610;
t291 = mrSges(5,3) + t306;
t649 = t291 * qJDD(1) - t421 * t382;
t609 = 0.2e1 * qJD(3);
t647 = t392 + (t609 + t610) * qJD(4);
t356 = t377 ^ 2;
t646 = 0.2e1 * t356;
t645 = 0.2e1 * t415;
t353 = t373 ^ 2;
t644 = 0.2e1 * t353 - 0.1e1;
t641 = t346 * t376;
t525 = t371 * t378;
t640 = 0.2e1 * t652;
t639 = t376 * t465 - t407;
t391 = qJD(5) ^ 2;
t497 = -t305 * qJDD(5) - t306 * t391;
t638 = -qJDD(1) * pkin(15) - t420;
t637 = -t374 * t344 + t148;
t636 = t290 * t353 - t373 * t543;
t635 = t421 * t380;
t413 = -t306 * qJDD(5) + t305 * t391;
t164 = t371 * t270 - t271 * t378;
t490 = qJD(6) * t303;
t626 = qJD(1) * (-t46 * qJD(10) + t43 * qJD(3) - qJD(4) * t54 + t25 * qJD(7) - qJD(8) * t88 - t108 * qJD(9) - t438 + t490);
t624 = -t328 + t640;
t285 = mrSges(4,3) + t291;
t623 = -t285 * qJDD(1) + t421 * t352;
t363 = mrSges(10,3) + mrSges(9,3);
t622 = -qJDD(1) * t363 + t421 * t381;
t621 = m(7) * (pkin(14) * t395 - t421) + t395 * t627;
t603 = 0.4e1 * t355 - 0.2e1;
t619 = t603 * t328 + 0.4e1 * t441;
t618 = -qJD(9) * t610 - t387;
t617 = qJD(8) * t610 + t388 + t394;
t386 = qJD(10) ^ 2;
t437 = qJD(4) * t488;
t429 = 0.2e1 * t437;
t481 = qJD(9) * qJD(8);
t448 = -0.2e1 * t481;
t470 = qJD(10) * qJD(7);
t566 = qJD(4) * t688;
t615 = t427 + (-qJD(3) * t45 - qJD(7) * t689 - qJD(8) * t661 - qJD(9) * t109 + t566) * t610 + t566 * t609 - t164 * t429 + t109 * t448 + (-qJD(10) * t610 - t386 - 0.2e1 * t470) * t47 + t706;
t614 = 0.2e1 * t346;
t354 = t374 ^ 2;
t613 = 0.2e1 * t354;
t612 = -0.2e1 * t378;
t607 = m(7) * g(3);
t335 = pkin(1) * t352;
t598 = g(3) * t352;
t597 = g(3) * t380;
t596 = g(3) * t381;
t595 = g(3) * t382;
t594 = t370 * pkin(5);
t342 = t371 * g(3);
t592 = t378 * g(3);
t317 = t371 * pkin(15) - pkin(1);
t586 = pkin(15) * t381;
t585 = qJD(1) * t9;
t360 = mrSges(11,3) + mrSges(8,3);
t439 = t164 * t488;
t569 = ((-t634 * t377 - 0.2e1 * t669) * t371 + 0.2e1 * t518) * qJD(4) - 0.2e1 * t439;
t568 = (-0.2e1 * t519 - t667) * qJD(4) + 0.2e1 * t438;
t567 = mrSges(7,3) * qJD(1);
t564 = t167 * t356;
t318 = -t371 * pkin(1) + pkin(15);
t545 = t352 * t318;
t538 = t366 * t371;
t537 = t366 * t378;
t534 = t369 * t382;
t532 = t370 * t262;
t530 = t370 * t377;
t527 = t371 * t365;
t526 = t371 * t370;
t222 = mrSges(3,3) + t285 + t360 + t363;
t217 = -mrSges(2,2) + mrSges(7,3) + t222;
t524 = t372 * t217;
t523 = t372 * t305;
t522 = t372 * t306;
t521 = t373 * t371;
t515 = t378 * t370;
t514 = t378 * t395;
t513 = t379 * t305;
t512 = t379 * t306;
t511 = t380 * t317;
t510 = t380 * t318;
t509 = t382 * (t317 * t370 - pkin(5));
t508 = t382 * (t317 - t594);
t320 = pkin(1) + 0.2e1 * t594;
t507 = t382 * t320;
t321 = t370 * pkin(1) + pkin(5);
t506 = t382 * t321;
t451 = 0.2e1 * t488;
t504 = (t197 * t370 + 0.2e1 * t671) * qJD(4) + t270 * t451;
t452 = -0.2e1 * t488;
t503 = t651 * qJD(4) + t271 * t452;
t500 = t634 * qJD(4) + t376 * t452;
t499 = t197 * qJD(4) + t369 * t451;
t348 = qJD(2) + qJD(8);
t331 = qJD(9) + t348;
t494 = qJD(1) * t331;
t350 = qJD(2) + qJD(3);
t332 = qJD(4) + t350;
t493 = qJD(1) * t332;
t349 = qJD(2) + qJD(7);
t492 = qJD(1) * t349;
t483 = qJD(1) * qJD(2);
t482 = qJD(1) * qJD(6);
t478 = t630 * qJDD(1);
t268 = t374 * t378 - t538;
t477 = t268 * qJDD(1);
t269 = t374 * t371 + t537;
t476 = t269 * qJDD(1);
t463 = qJDD(2) + qJDD(3);
t347 = qJDD(2) + qJDD(7);
t462 = qJDD(2) + qJDD(8);
t450 = t365 * t586;
t445 = t199 * t530;
t443 = t290 * t527;
t442 = t374 * t670;
t440 = t328 * t531;
t435 = t378 * t483;
t433 = 0.4e1 * t440;
t402 = t395 * pkin(15) + t421;
t214 = t378 * t402 + t342;
t430 = -0.2e1 * t437;
t426 = 0.2e1 * t436;
t424 = t603 * t532;
t213 = t402 * t371 - t592;
t274 = t371 * t377 + t515;
t241 = t274 * qJD(1);
t275 = t378 * t377 - t526;
t244 = t275 * qJD(1);
t141 = -t369 * t241 + t376 * t244;
t266 = t365 * t364 - t373 * t601;
t267 = t373 * t364 + t365 * t601;
t417 = t371 * t266 - t267 * t378;
t273 = t378 * t365 + t521;
t272 = t373 * t378 - t527;
t104 = 0.2e1 * t108;
t144 = -t446 + t197;
t132 = t144 * t377;
t410 = ((t132 + 0.2e1 * t684) * t371 - 0.2e1 * t579) * qJD(3) + t104 * qJD(9) + t568 + t695;
t408 = 0.2e1 * t636;
t404 = -t222 * qJDD(1) + t426;
t278 = pkin(5) * t534 + t328;
t194 = 0.2e1 * t416;
t299 = m(7) + m(2) + t313;
t280 = -t371 * qJDD(1) - t435;
t399 = t638 + (-t280 + t435) * pkin(1);
t396 = t410 * qJD(2) + t568 * qJD(3) - t702 * qJDD(2) - t413 * t164 - t9 * t394 + t109 * qJDD(9) + t108 * t387 - t630 * t430 + t689 * qJDD(7) - t25 * t389 - t627 * qJDD(6) - t303 * t390 + t40 * t470 + t45 * qJDD(3) - t43 * t393 + t47 * qJDD(10) + t46 * t386 - t688 * qJDD(4) + t54 * t392 + t661 * qJDD(8) + t88 * t388 + t104 * t481;
t345 = t349 ^ 2;
t330 = qJD(10) + t349;
t327 = qJDD(4) + t463;
t326 = qJDD(9) + t462;
t325 = t332 ^ 2;
t322 = qJDD(10) + t347;
t301 = -qJD(6) * mrSges(7,2) + t375 * t567;
t300 = qJD(6) * mrSges(7,1) - t367 * t567;
t289 = m(1) + t299;
t282 = t378 * qJDD(1) - t371 * t483;
t281 = t375 * qJDD(1) - t367 * t482;
t279 = t367 * qJDD(1) + t375 * t482;
t264 = -t320 * t371 + pkin(15);
t261 = t358 * t361 + t359 * t362;
t260 = t358 * t362 - t359 * t361;
t247 = -pkin(5) * t371 + t318 * t370;
t243 = t268 * qJD(1);
t242 = t272 * qJD(1);
t240 = t269 * qJD(1);
t239 = t273 * qJD(1);
t215 = -t262 + t641;
t212 = t217 * t379;
t202 = (-pkin(1) * t520 + t321 * t369) * t382 + t328;
t187 = (-t371 ^ 2 * t395 - t394) * pkin(1) + t213;
t165 = (-t371 * t514 + qJDD(2)) * pkin(1) + t214;
t160 = t266 * t378 + t371 * t267;
t159 = -0.2e1 * t169;
t158 = t269 * t492 - t477;
t157 = t268 * t492 + t476;
t156 = (pkin(1) * t535 + t321 * t376) * t382 - t262;
t154 = -t377 * t335 - t199;
t143 = t376 * t241 + t369 * t244;
t142 = t364 * t239 - t601 * t242;
t140 = t601 * t239 + t364 * t242;
t127 = qJD(5) - t141;
t126 = (t240 * t359 - t243 * t358) * pkin(4);
t125 = t366 * t344 + t147;
t118 = t379 * t439;
t117 = t372 * t439;
t116 = t331 * mrSges(10,1) - t142 * mrSges(10,3);
t115 = -t331 * mrSges(10,2) + t140 * mrSges(10,3);
t114 = t602 * t143 + t368 * t332;
t113 = -t368 * t143 + t602 * t332;
t107 = -t260 * t240 - t261 * t243;
t106 = t261 * t240 - t260 * t243;
t105 = -t371 * t645 + 0.2e1 * t516;
t91 = pkin(1) * t382 + t680;
t83 = ((t614 + t496) * t377 + t183) * t371;
t82 = t330 * mrSges(11,1) - t107 * mrSges(11,3);
t81 = -t330 * mrSges(11,2) + t106 * mrSges(11,3);
t80 = t335 + t94;
t79 = t373 * t213 + t365 * t214 + (-t239 ^ 2 - t348 ^ 2) * pkin(2);
t73 = -qJDD(1) * t164 - t630 * t493;
t72 = -t417 * qJDD(1) - t160 * t494;
t71 = t160 * qJDD(1) - t417 * t494;
t70 = -t164 * t493 + qJDD(5) + t478;
t65 = -t141 * pkin(9) - t143 * pkin(11);
t64 = -t140 * mrSges(10,1) + t142 * mrSges(10,2);
t61 = -t365 * t213 + t373 * t214 + (-t242 * t239 + t462) * pkin(2);
t60 = -t377 * t165 + t370 * t187 + (-t244 ^ 2 - t350 ^ 2) * pkin(5);
t59 = t127 * mrSges(6,1) - t114 * mrSges(6,3);
t58 = -t127 * mrSges(6,2) + t113 * mrSges(6,3);
t57 = -t113 * mrSges(6,1) + t114 * mrSges(6,2);
t56 = t370 * t165 + t377 * t187 + (t241 * t244 + t463) * pkin(5);
t52 = -t106 * mrSges(11,1) + t107 * mrSges(11,2);
t42 = -t90 * t371 + 0.2e1 * t578;
t36 = t113 * qJD(5) + t368 * t327 + t602 * t73;
t35 = -t114 * qJD(5) + t602 * t327 - t368 * t73;
t32 = -t240 * t126 + t366 * t165 + t374 * t187 + (-t345 * t359 + t347 * t358) * pkin(4);
t31 = -t243 * t126 + t374 * t165 - t366 * t187 + (t345 * t358 + t347 * t359) * pkin(4);
t28 = t106 * qJD(10) - t157 * t260 + t158 * t261;
t27 = -t107 * qJD(10) + t157 * t261 + t158 * t260;
t11 = -t325 * pkin(9) + t327 * pkin(11) + t141 * t65 + t369 * t56 + t376 * t60;
t10 = (-t141 * t332 - t73) * pkin(11) + (t478 + (-qJD(1) * t164 + t143) * t332) * pkin(9) + (-t370 * t280 - t377 * t282 + (qJD(3) + t350) * t241) * pkin(5) + t399;
t1 = [(-t513 * t630 - t522) * t391 + (t512 * t630 - t523) * qJDD(5) - t524 * t395 + (t372 * t705 + t212) * qJDD(1) - t289 * g(1) - t118 * t609 + (t372 * t585 - t118) * t610 + 0.2e1 * t372 * t626 + t615 * t379, -t413 * t630 + t217 * qJDD(1) - t47 * t386 - t421 * t299 + t569 * qJD(3) + t42 * t470 + t105 * t481 + (((-t176 * t377 + 0.2e1 * t685) * t371 - 0.2e1 * t572) * qJD(3) + (-t371 * t662 - 0.2e1 * t581) * qJD(7) + t42 * qJD(10) + (t159 * t371 + 0.2e1 * t155) * qJD(8) + t105 * qJD(9) + t569) * qJD(2) + 0.2e1 * (-qJD(4) * t164 * t306 - t664) * qJD(5) + t706, t67 * qJDD(7) + ((t144 * t370 - 0.2e1 * t683) * qJD(3) + t662 * qJD(7) + t90 * qJD(10) + 0.2e1 * t169 * qJD(8) + t645 * qJD(9) + t504) * qJD(2) - t660 * qJDD(10) + t690 * qJDD(3) + t645 * t481 + t20 * qJDD(2) - t111 * t392 - t94 * t393 + t270 * t429 + t97 * t386 + t415 * t387 + t169 * t388 + t691 * t389 + t504 * qJD(3) - t678 * qJDD(4) - t416 * qJDD(9) + t90 * t470 - t313 * t592 - t171 * qJDD(8) + t413 * t271 + (t20 * t514 + t313 * t402 + t404) * t371 + (t394 + t704) * t14, (t144 * qJD(3) + t499) * qJD(2) + t499 * qJD(3) - t274 * t598 + ((t144 * t356 + t412 + 0.2e1 * t445) * t357 + ((0.2e1 * t167 * t526 - t545) * t377 + (t646 - 0.1e1) * t371 * t199) * t378 + t564 - t445 + t317 * t544 - t167) * t395 - t167 * t393 - t199 * qJDD(3) + t369 * t429 + t628 * t392 + t203 * qJDD(4) - t412 * t394 + t154 * qJDD(2) + t413 * t376 + (t427 - t623) * t275, t164 * t595 - t156 * t394 + (((t639 - t641) * t646 + (t369 * t507 - t370 * t619) * t377 + (t294 + t506) * t376 + t407) * t357 + (-t382 * t264 * t520 + t247 * t534) * t378 + ((t294 + t346) * t376 + t407) * t356 + (t369 * t508 + t370 * t624) * t377 + (t628 + t509) * t376 + (0.4e1 * (-t311 + t608 / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + (-t557 + t346 / 0.2e1) * t369) * t356 + (t433 * t376 - t424) * t377 + t624) * t525) * t395 + t278 * qJDD(3) + t202 * qJDD(2) + t328 * qJDD(4) + t413 + (-t393 - t650) * t215 + t647 * t262 - (t427 + t649) * t630, m(6) * (t368 * t10 + t602 * t11) + t35 * mrSges(6,3) - t70 * mrSges(6,2) + t113 * t57 - t127 * t59, -qJDD(6) * mrSges(7,2) + t281 * mrSges(7,3) - qJD(6) * t300 - t367 * t607 + t375 * t621, -t268 * t597 + t175 * t470 + (t175 * qJD(10) + t129 * qJD(7)) * qJD(2) + ((-t129 * t354 + (-t366 * t633 + t344) * t374 + t148) * t357 + t510 * t537 + (t691 + t511) * t374 + (t147 * t613 - t147 + 0.2e1 * t442) * t525) * t395 + t125 * qJDD(2) + t637 * t394 - t629 * qJDD(10) + t192 * t386 + t147 * qJDD(7) + t148 * t389 - t360 * t476 + t635 * t269, -t272 * t596 - t302 * qJDD(9) + t290 * qJDD(2) + (t365 * t468 * t357 + t378 * t450 + t169 * t373 + ((-t290 + t408) * t378 + t373 * t586) * t371) * t395 + t290 * qJDD(8) + (t644 * t357 * t395 - t617) * t263 + (0.2e1 * t481 - t618) * t307 + t622 * t273, m(10) * (-t364 * t61 - t601 * t79) + t72 * mrSges(10,3) - t326 * mrSges(10,2) + t140 * t64 - t331 * t116, m(11) * (t260 * t31 - t261 * t32) + t27 * mrSges(11,3) - t322 * mrSges(11,2) + t106 * t52 - t330 * t82, 0, 0, 0, 0, 0, 0; (t522 * t630 + t513) * qJDD(5) + (-t523 * t630 + t512) * t391 + (-t379 * t705 + t524) * qJDD(1) + t212 * t395 - t289 * g(2) - t117 * t609 + (-t379 * t585 - t117) * t610 - 0.2e1 * t379 * t626 + t615 * t372, -t705 * qJDD(1) + t217 * t395 + t420 * t299 + (-0.2e1 * t490 + ((t132 - t677) * t371 + 0.2e1 * t582) * qJD(2) + t410) * qJD(1) - t497, -t691 * qJDD(7) - t97 * qJDD(10) + t313 * t342 - t194 * t481 + (-qJD(10) * t659 + t85 * qJD(3) - t676 * qJD(7) + t654 * qJD(8) - t194 * qJD(9) + t503) * qJD(2) - t415 * qJDD(9) + t111 * qJDD(4) - t14 * qJDD(2) - t678 * t392 + t690 * t393 + t20 * t394 + t271 * t430 - t660 * t386 - t416 * t387 - t171 * t388 + t67 * t389 + t503 * qJD(3) - t659 * t470 + t94 * qJDD(3) - t169 * qJDD(8) + t413 * t270 + (t313 * t421 - t4 * t395 + t404) * t378, (t176 * qJD(3) + t500) * qJD(2) + t500 * qJD(3) - t275 * t598 + t154 * t394 + t412 * qJDD(2) + ((-t176 * t356 + (t145 * t370 + t335) * t377 - t199) * t357 + t515 * t545 - (-t352 * t317 + t94) * t377 + (-t167 - 0.2e1 * t445 + 0.2e1 * t564) * t525) * t395 - t199 * t393 + t167 * qJDD(3) + t203 * t392 - t628 * qJDD(4) + t376 * t430 + t413 * t369 + (t426 + t623) * t274, t630 * t595 + (((-0.4e1 * t311 + 0.2e1 * t608 - 0.2e1 * mrSges(5,2) + 0.2e1 * mrSges(6,3) + (-0.4e1 * t557 + t614) * t369) * t356 + ((t433 + t507) * t376 - t424) * t377 - t369 * t506 + t624) * t357 + ((t247 * t376 + t264 * t535) * t382 + (-0.4e1 * (t558 + (-t547 - t346 / 0.2e1) * t376 - t588 / 0.2e1 + t341 / 0.2e1 + t385 / 0.2e1 + mrSges(5,1) / 0.2e1) * t356 + t619 * t530 + t639) * t371) * t378 + (-t278 + t640) * t356 + ((-0.2e1 * t440 + t508) * t376 + (0.2e1 * t355 - 0.1e1) * t532) * t377 - t369 * t509 - t652) * t395 + t156 * qJDD(2) + t202 * t394 + t278 * t650 + t332 * t452 + t278 * t393 + t215 * qJDD(3) - t262 * qJDD(4) - (t426 - t649) * t164 + (t395 + t647) * t328, m(6) * (t602 * t10 - t368 * t11) - t36 * mrSges(6,3) + t70 * mrSges(6,1) - t114 * t57 + t127 * t58, qJDD(6) * mrSges(7,1) - t279 * mrSges(7,3) + qJD(6) * t301 - t367 * t621 - t375 * t607, t269 * t597 - t632 * t470 + (-qJD(10) * t632 - qJD(7) * t633) * qJD(2) + t147 * t389 + ((-t633 * t354 + (0.2e1 * t138 - t344) * t366 - t147) * t357 + ((-0.2e1 * t147 * t538 + t510) * t374 + (t613 - 0.1e1) * t371 * t148) * t378 - t147 * t354 - t442 - t366 * t511 + t147) * t395 - t192 * qJDD(10) - t629 * t386 - t148 * qJDD(7) - t637 * qJDD(2) + t125 * t394 - t360 * t477 + t635 * t268, t273 * t596 + (t408 * t357 + ((-0.2e1 * t443 + t586) * t373 - t644 * t371 * t263) * t378 - t371 * t450 - t636) * t395 - t307 * qJDD(9) + t263 * qJDD(2) + t263 * qJDD(8) + (t704 + t617) * t290 + (t448 + t618) * t302 + t622 * t272, m(10) * (t364 * t79 - t601 * t61) - t71 * mrSges(10,3) + t326 * mrSges(10,1) - t142 * t64 + t331 * t115, m(11) * (-t260 * t32 - t261 * t31) - t28 * mrSges(11,3) + t322 * mrSges(11,1) - t107 * t52 + t330 * t81, 0, 0, 0, 0, 0, 0; -t289 * g(3) + t396, -t299 * g(3) + t396, t4 * qJDD(1) - t222 * t395 - t420 * t313 + (((t133 + t677) * t371 - 0.2e1 * t582) * qJD(2) + t696 + (-t111 * t612 + t667) * qJD(4) + t121 + (-t194 * t371 - 0.2e1 * t517) * qJD(9) - t695) * qJD(1) + t497, (t80 * t371 - t572) * qJDD(1) - t285 * t395 + t638 * t352 + ((t80 * t611 + t75) * qJD(2) + t696 + t686) * qJD(1) + t497, (t91 * t371 + (-t377 * t666 - t672) * t378) * qJDD(1) - t291 * t395 + t638 * t382 + ((t91 * t611 + t83) * qJD(2) + (t680 * t611 + t83) * qJD(3) + t686) * qJD(1) + t497, m(6) * (-t327 * pkin(9) - t325 * pkin(11) + t143 * t65 + t369 * t60 - t376 * t56) + t36 * mrSges(6,2) - t35 * mrSges(6,1) + t114 * t59 - t113 * t58, m(7) * (pkin(14) * qJDD(1) - t420) + t279 * mrSges(7,2) - t281 * mrSges(7,1) + (t367 * t300 - t375 * t301) * qJD(1), (t697 * t371 - t581) * qJDD(1) - t360 * t395 + t638 * t380 + ((-t697 * t612 - t62) * qJD(2) + (t612 * t691 - t62) * qJD(7) - t681) * qJD(1), (t263 * t521 + t155 + t443) * qJDD(1) - t363 * t395 + t638 * t381 + ((-t378 * t645 - 0.2e1 * t655) * qJD(9) + t348 * (t159 * t378 - 0.2e1 * t562)) * qJD(1), m(10) * ((-t373 * t280 + t365 * t282 + (qJD(8) + t348) * t242) * pkin(2) + t638) + t71 * mrSges(10,2) - t72 * mrSges(10,1) + t142 * t116 - t140 * t115, m(11) * t399 + t28 * mrSges(11,2) - t27 * mrSges(11,1) + t107 * t82 - t106 * t81 + (-(-t240 * qJD(7) + t366 * t280 + t374 * t282) * t358 - (-t243 * qJD(7) + t374 * t280 - t366 * t282) * t359 + (t240 * t358 + t243 * t359) * t349) * t605, 0, 0, 0, 0, 0, 0;];
f_new = t1;
